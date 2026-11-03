# AIM ---------------------------------------------------------------------
# Exploratory QC BEFORE deciding the final sample selection and the model.
# Starting from the table of raw counts (00_import_data.R) the script
#   1. describes the experimental design (balance, confounding between the variables)
#   2. checks the library size / number of detected genes per sample
#   3. runs sample clustering, sample distance (euclidean on vst, Poisson) and PCA with a blind transformation (no design used)
#   4. quantifies which metadata explain each PC and flags outlier samples
#   5. checks which candidate designs are estimable with the available samples
#   6. infers the sex of the samples from XIST / Y-linked genes (heatmap) and checks the agreement within each clone
#   7. shows the PCA after removing a candidate batch (clone), for visualisation only
#   8. repeats the main checks (contrasts replication, estimable designs, PCA, clustering, sex, p-value histograms)
#      on the dataset left after removing the samples in drop_final (final section, plots in plot/explore_final)
# Nothing here changes the dds used downstream: the choices made after reading these plots go into 01_import_data.R (samples + design) and 02_QC.R (filter group).
# Run from the project root (analysis_R45_pixi), paths are relative to it.

# libraries ---------------------------------------------------------------
library(tidyverse)
library(DESeq2)
library(edgeR)
library(ComplexHeatmap)
library(RNAseqQC)
library(PoiClaClu)
library(ggrepel)
library(patchwork)
library(viridis)
library(DescTools)
library(shadowtext)
library(pals)

# parameters --------------------------------------------------------------
data_dir <- "../../data"
out_dir <- "../../out"

# samples to drop from the exploration (fill after looking at the plots, e.g. c("24A8"));
# the same list has to be applied in 01_import_data.R
exclude_samples <- character(0)

# samples that would be removed for the DE: only used by the final section ("assessment after removal"),
# everything above still describes the full set of samples.
# 24A8, 6A8: half the library size / ~10.7k genes detected, extremes of PC1 (technical failures)
# 6B3, 6B7: the only Iron / Iron_myelin samples (no replicate, cannot be tested)
drop_final <- c("24A8", "6A8", "6B3", "6B7")

# metadata used to colour the plots / explain the PCs
vars_meta <- c("clone", "mut", "time", "treat")

# group used by filterByExpr: keep a gene if expressed in the samples of the smallest group
# (groups of 1 sample make the filter lenient, here mut x treat is used)
filter_group <- c("mut", "treat")

n_pc <- 6 # number of PCs tested against the metadata

# genes used to infer the sex of the samples (symbol = gene_id in the GTF)
sex_genes <- c(XIST = "female",
               TSIX = "female",
               DDX3Y = "male",
               RPS4Y1 = "male",
               USP9Y = "male",
               UTY = "male",
               EIF1AY = "male",
               ZFY = "male")

plot_dir <- file.path(out_dir, "plot", "explore")

# read in the data --------------------------------------------------------
counts <- readRDS(file.path(out_dir, "object", "raw_count_salmon_GTF.rds"))
samplesheet <- read_csv(file.path(data_dir, "samplesheet.csv"), show_col_types = FALSE)

# samples with counts and in the samplesheet, in the order of the matrix
coldata <- data.frame(sample = colnames(counts)) %>%
  left_join(samplesheet, by = c("sample" = "sample_id")) %>%
  filter(!sample %in% exclude_samples) %>%
  mutate(across(all_of(vars_meta), factor),
         group = paste(mut, treat, time, sep = "|"),
         sample_name = paste(sample, clone, mut, time, treat, sep = "_")) %>%
  column_to_rownames("sample") %>%
  mutate(sample = rownames(.))

counts <- counts[, rownames(coldata)]
storage.mode(counts) <- "integer"

message("samples with counts: ", ncol(counts), " | samples in the samplesheet: ", nrow(samplesheet))
message("samplesheet samples without counts: ", paste(setdiff(samplesheet$sample_id, colnames(counts)), collapse = ", "))

# design balance ----------------------------------------------------------
# how many samples per combination of the variables
coldata %>%
  dplyr::count(time, mut, treat) %>%
  as_tibble() %>%
  print(n = Inf)

# clone x condition: a clone seen in only one condition cannot be used as blocking factor for it
coldata %>%
  dplyr::count(clone, mut, time, treat) %>%
  pivot_wider(names_from = c(mut, time, treat), values_from = n, values_fill = 0) %>%
  as_tibble() %>%
  print(width = Inf)

# confounding between the variables (Cramer's V: 1 = one variable determines the other)
# correct = TRUE: bias-corrected V, the uncorrected one is inflated with few samples and many levels
mat_cramer <- DescTools::PairApply(coldata[vars_meta],
                                   function(x, y) DescTools::CramerV(x, y, correct = TRUE),
                                   symmetric = TRUE)
as.data.frame(mat_cramer) %>%
  rownames_to_column("v1") %>%
  pivot_longer(-v1, names_to = "v2", values_to = "cramer_v") %>%
  ggplot(aes(v1, v2, fill = cramer_v, label = round(cramer_v, 2))) +
  geom_tile() + geom_text() +
  scale_fill_viridis_c(limits = c(0, 1)) +
  labs(title = "Association between the variables (bias-corrected Cramer's V)", x = NULL, y = NULL) +
  theme_minimal()
ggsave(file.path(plot_dir, "confounding_variables.pdf"), width = 5, height = 4)

# library size and detected genes ----------------------------------------
df_lib <- coldata %>%
  mutate(tot_counts = colSums(counts)[rownames(.)],
         n_genes_detected = colSums(counts >= 10)[rownames(.)])

p_lib <- ggplot(df_lib, aes(x = reorder(sample_name, tot_counts), y = tot_counts, fill = treat)) +
  geom_col() + coord_flip() + labs(x = NULL, y = "total counts") + theme_bw()
p_det <- ggplot(df_lib, aes(x = reorder(sample_name, tot_counts), y = n_genes_detected, fill = treat)) +
  geom_col() + coord_flip() + labs(x = NULL, y = "genes with >= 10 counts") + theme_bw() +
  theme(axis.text.y = element_blank())
p_lib + p_det + plot_layout(guides = "collect")
ggsave(file.path(plot_dir, "libsize_detected_genes.pdf"), width = 11, height = 5)

df_lib %>%
  dplyr::select(sample_name, tot_counts, n_genes_detected) %>%
  arrange(tot_counts) %>%
  print()

# filter + blind transformation ------------------------------------------
# blind = TRUE: the design is not used, so the plots show the structure of the data as it is
dds_base <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ 1)
# keep the unfiltered object: the Y-linked genes are lowly expressed in females and may be removed by the filter
dds_all <- dds_base
feat_keep <- edgeR::filterByExpr(counts(dds_base), group = do.call(paste, c(coldata[filter_group], sep = "|")))
table(feat_keep)
dds_base <- dds_base[feat_keep, ]

vds <- vst(dds_base, blind = TRUE)
vds$sample_label <- vds$sample_name

# sample distance ---------------------------------------------------------
dist_vst <- dist(t(assay(vds)))
mat_dist <- as.matrix(dist_vst)

# euclidean distance on vst, clustered, annotated with all the variables
pdf(file.path(plot_dir, "heatmap_sampleCluster_vst.pdf"), width = 9, height = 6)
set.seed(1)
hm <- plot_sample_clustering(vds, anno_vars = vars_meta, distance = "euclidean")
draw(hm, heatmap_legend_side = "left", annotation_legend_side = "left",
     padding = unit(c(2, 2, 2, 30), "mm"))
dev.off()

# Poisson distance on the raw counts (uses the counts, not the vst)
poisd <- PoissonDistance(t(counts(dds_base)))
mat_pois <- as.matrix(poisd$dd)
rownames(mat_pois) <- colnames(mat_pois) <- vds$sample_label
pdf(file.path(plot_dir, "heatmap_poisson.pdf"), width = 8, height = 6)
Heatmap(mat_pois, name = "Poisson dist",
        col = circlize::colorRamp2(seq(min(mat_pois), max(mat_pois), length.out = 255), viridis(255)),
        clustering_distance_rows = poisd$dd, clustering_distance_columns = poisd$dd,
        top_annotation = HeatmapAnnotation(df = as.data.frame(colData(vds)[, vars_meta])),
        show_column_names = FALSE) %>%
  draw()
dev.off()

# outlier score: median distance of each sample to all the others (z-score of the median distance)
df_outlier <- tibble(sample = rownames(mat_dist),
                     median_dist = apply(mat_dist, 1, function(x) median(x[x > 0]))) %>%
  mutate(z = as.numeric(scale(median_dist))) %>%
  left_join(df_lib %>% dplyr::select(sample, sample_name, tot_counts, n_genes_detected), by = "sample") %>%
  arrange(desc(z))
print(df_outlier)
write_tsv(df_outlier, file.path(out_dir, "table", "explore_outlier_score.tsv"))

# PCA ---------------------------------------------------------------------
# all the genes passing the filter, and the 500 most variable (plotPCA default) for comparison
pca_all <- prcomp(t(assay(vds)))
pv <- 100 * pca_all$sdev^2 / sum(pca_all$sdev^2)

# scree
tibble(PC = factor(seq_along(pv)), var = pv) %>%
  ggplot(aes(PC, var)) + geom_col() + labs(y = "% variance") + theme_bw()
ggsave(file.path(plot_dir, "PCA_scree.pdf"), width = 5, height = 3.5)

df_pca <- cbind(as.data.frame(colData(vds)[, c("sample", "sample_name", vars_meta)]),
                pca_all$x[, seq_len(max(n_pc, 4))]) %>%
  left_join(df_lib %>% dplyr::select(sample, tot_counts, n_genes_detected), by = "sample")

# pairs of PCs coloured by treat, clone and mut: one page per covariate, the pairs of PCs wrapped on the page
pcs_plot <- paste0("PC", 1:4)
p_pca <- crossing(covariate = c("treat", "clone", "mut"),
                  PCx = pcs_plot,
                  PCy = pcs_plot) %>%
  filter(PCx < PCy) %>%
  group_split(covariate) %>%
  map(function(df_pairs) {
    df_pairs %>%
      pmap(function(covariate, PCx, PCy) {
        df_pca %>%
          ggplot(aes(x = .data[[PCx]], y = .data[[PCy]], col = .data[[covariate]], label = sample)) +
          geom_point(size = 3, alpha = 0.8) +
          ggrepel::geom_text_repel(show.legend = FALSE, size = 3) +
          labs(x = sprintf("%s: %.0f%%", PCx, pv[parse_number(PCx)]),
               y = sprintf("%s: %.0f%%", PCy, pv[parse_number(PCy)]),
               col = covariate) +
          theme_bw() + theme(strip.background = element_blank())
      }) %>%
      wrap_plots(nrow = 2, guides = "collect") +
      plot_annotation(title = paste("PCA coloured by", unique(df_pairs$covariate)))
  })

# multi-page pdf, one page per covariate
pdf(file.path(plot_dir, "PCA_byVariable.pdf"), width = 16, height = 8)
walk(p_pca, print)
dev.off()

# which metadata explain each PC -----------------------------------------
# variance of each PC explained by each variable (R2 of a linear model, p-value of the ANOVA);
# the library size (log10) and the number of detected genes are included as technical covariates
df_pca <- df_pca %>% mutate(log10_tot_counts = log10(tot_counts))
df_pc_assoc <- crossing(PC = paste0("PC", seq_len(n_pc)),
                         variable = c(vars_meta, "log10_tot_counts", "n_genes_detected")) %>%
  mutate(res = map2(PC, variable, function(pc, v) {
    fit <- lm(df_pca[[pc]] ~ df_pca[[v]])
    tibble(R2 = summary(fit)$r.squared, p = anova(fit)$`Pr(>F)`[1])
  })) %>%
  unnest(res) %>%
  mutate(PC = factor(PC, levels = paste0("PC", seq_len(n_pc))),
         label = sprintf("%.2f%s", R2, if_else(p < 0.05, "*", "")))
write_tsv(df_pc_assoc, file.path(out_dir, "table", "explore_PC_vs_metadata.tsv"))

ggplot(df_pc_assoc, aes(PC, variable, fill = R2, label = label)) +
  geom_tile() +
  geom_shadowtext(size = 3, color = "white", bg.color = "black", bg.r = 0.15) +
  scale_fill_viridis_c(limits = c(0, 1)) +
  labs(title = "Variance of each PC explained by each variable (R2, * = p < 0.05)", x = NULL, y = NULL) +
  theme_minimal()
ggsave(file.path(plot_dir, "PCA_PC_vs_metadata_R2.pdf"), width = 7, height = 4)

# PC scores by metadata
df_pca %>%
  dplyr::select(sample, all_of(vars_meta), paste0("PC", seq_len(n_pc))) %>%
  mutate(across(all_of(vars_meta), as.character)) %>%
  pivot_longer(all_of(vars_meta), names_to = "variable", values_to = "level") %>%
  pivot_longer(starts_with("PC"), names_to = "PC", values_to = "score") %>%
  ggplot(aes(level, score)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.15) +
  facet_grid(PC ~ variable, scales = "free") +
  theme_bw() + theme(strip.background = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(plot_dir, "PCA_PC_by_metadata.pdf"), width = 10, height = 12)

# sex check ---------------------------------------------------------------
# the samplesheet has no sex: infer it from XIST (female) and the Y-linked genes (male) and check that the samples of
# the same clone agree (a clone cannot change sex: discordant samples can be a swap or a contamination).
# The vst is computed on ALL the genes (dds_all): the Y-linked genes are lowly expressed in females and could be removed by the filter.
vds_all <- vst(dds_all, blind = TRUE)

annotations <- readRDS(file.path(out_dir, "object", "annotation_GTF.rds")) %>%
  dplyr::select(GENEID, SYMBOL) %>%
  distinct()

df_sex_genes <- annotations %>%
  dplyr::filter(SYMBOL %in% names(sex_genes), GENEID %in% rownames(vds_all)) %>%
  mutate(sex = unname(sex_genes[SYMBOL])) %>%
  arrange(factor(sex, levels = c("female", "male")), SYMBOL)

message("sex genes found: ", paste(df_sex_genes$SYMBOL, collapse = ", "))

mat_sex <- assay(vds_all)[df_sex_genes$GENEID, , drop = FALSE]
rownames(mat_sex) <- df_sex_genes$SYMBOL

# score per sample: mean vst of the Y genes minus the vst of XIST (> 0 male-like, < 0 female-like)
df_sex <- tibble(sample = colnames(mat_sex),
                 XIST = mat_sex["XIST", ],
                 Y_mean = colMeans(mat_sex[df_sex_genes$sex == "male", , drop = FALSE])) %>%
  mutate(sex_score = Y_mean - XIST) %>%
  left_join(coldata %>% dplyr::select(sample, sample_name, all_of(vars_meta)), by = "sample") %>%
  arrange(clone, sex_score)

# the score range within each clone: a large range = samples of the same clone with a different inferred sex
df_sex %>%
  group_by(clone) %>%
  summarise(n = n(), min_score = min(sex_score), max_score = max(sex_score),
            range = max_score - min_score, .groups = "drop") %>%
  print()

# heatmap: z-score of the vst per gene across the samples
mat_sex_z <- t(scale(t(mat_sex)))
mat_sex_z[is.na(mat_sex_z)] <- 0 # genes with no variation
cd_sex <- coldata[colnames(mat_sex_z), ]

col_clone <- setNames(pals::glasbey(nlevels(cd_sex$clone)), levels(cd_sex$clone))
col_mut   <- setNames(pals::okabe(nlevels(cd_sex$mut)), levels(cd_sex$mut))
col_time  <- setNames(pals::alphabet(nlevels(cd_sex$time)), levels(cd_sex$time))
col_treat <- setNames(pals::polychrome(nlevels(cd_sex$treat)), levels(cd_sex$treat))
column_ha <- HeatmapAnnotation(clone = cd_sex$clone,
                               mut = cd_sex$mut,
                               time = cd_sex$time,
                               treat = cd_sex$treat,
                               # sex_score = anno_barplot(df_sex$sex_score[match(colnames(mat_sex_z), df_sex$sample)]),
                               col = list(clone = col_clone,
                                          mut = col_mut,
                                          time = col_time,
                                          treat = col_treat))

colnames(mat_sex_z) <- cd_sex$sample_name
ht_sex <- Heatmap(mat_sex_z, name = "z-score vst",
                  col = circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
                  top_annotation = column_ha,
                  # row_split = df_sex_genes$sex,
                  cluster_row_slices = FALSE,
                  column_title = "sex genes (XIST: female; Y-linked: male)")
pdf(file.path(plot_dir, "heatmap_sex_genes.pdf"), width = 9, height = 6)
draw(ht_sex, heatmap_legend_side = "left", annotation_legend_side = "left")
dev.off()

# assessment after removal of the samples in drop_final --------------------
# The same checks on the dataset that would go to the DE: are the two contrasts of interest (Myelin vs BL, KO vs WT)
# replicated and estimable, what explains the variance now, are the samples consistent within each clone.
plot_dir_final <- file.path(out_dir, "plot", "explore_final")
dir.create(plot_dir_final, recursive = TRUE, showWarnings = FALSE)

coldata_f <- coldata %>%
  filter(!sample %in% drop_final) %>%
  mutate(across(where(is.factor), droplevels))
counts_f <- counts[, rownames(coldata_f)]

# variables that still vary (e.g. time is constant once the only 24h sample is removed)
vars_f <- vars_meta[map_lgl(coldata_f[vars_meta], ~ nlevels(.x) > 1)]
message("samples removed: ", paste(intersect(drop_final, colnames(counts)), collapse = ", "),
        " | samples left: ", ncol(counts_f))
message("variables still varying: ", paste(vars_f, collapse = ", "))

# replication of the contrasts --------------------------------------------
# samples per condition, and per clone: the clone is used as blocking factor, so a contrast is only
# informative within the clones that have both of its levels
coldata_f %>%
  dplyr::count(mut, treat) %>%
  as_tibble() %>%
  print(n = Inf)

coldata_f %>%
  dplyr::count(clone, mut, treat) %>%
  pivot_wider(names_from = c(mut, treat), values_from = n, values_fill = 0) %>%
  as_tibble() %>%
  print(width = Inf)

# estimable designs -------------------------------------------------------
# rank = number of columns and residual df > 0 (the interaction is expected to be not estimable: no KO Myelin sample)
cd_design <- coldata_f %>%
  mutate(treat = relevel(treat, ref = "BL"), mut = relevel(mut, ref = "WT"))

candidate_designs_f <- list(
  treat_mut = ~ treat + mut,
  clone_treat_mut = ~ clone + treat + mut,
  mut_x_treat = ~ mut * treat,
  clone_mut_x_treat = ~ clone + mut * treat
)

df_designs_f <- imap_dfr(candidate_designs_f, function(f, nm) {
  m <- model.matrix(f, cd_design)
  tibble(design = nm, formula = deparse(f), n_samples = nrow(m), n_coef = ncol(m),
         rank = qr(m)$rank, residual_df = nrow(m) - qr(m)$rank,
         full_rank = qr(m)$rank == ncol(m))
})
print(df_designs_f, width = Inf)
write_tsv(df_designs_f, file.path(out_dir, "table", "explore_candidate_designs_final.tsv"))

# filter + blind transformation -------------------------------------------
dds_f <- DESeqDataSetFromMatrix(countData = counts_f, colData = coldata_f, design = ~ 1)
dds_f_all <- dds_f # unfiltered, kept for the sex genes
feat_keep_f <- edgeR::filterByExpr(counts(dds_f), group = do.call(paste, c(coldata_f[filter_group], sep = "|")))
table(feat_keep_f)
dds_f <- dds_f[feat_keep_f, ]

vds_f <- vst(dds_f, blind = TRUE)
vds_f$sample_label <- vds_f$sample_name

# sample distance and outlier score ---------------------------------------
pdf(file.path(plot_dir_final, "heatmap_sampleCluster_vst.pdf"), width = 9, height = 6)
set.seed(1)
hm_f <- plot_sample_clustering(vds_f, anno_vars = vars_f, distance = "euclidean")
draw(hm_f, heatmap_legend_side = "left", annotation_legend_side = "left",
     padding = unit(c(2, 2, 2, 30), "mm"))
dev.off()

mat_dist_f <- as.matrix(dist(t(assay(vds_f))))
df_outlier_f <- tibble(sample = rownames(mat_dist_f),
                       median_dist = apply(mat_dist_f, 1, function(x) median(x[x > 0]))) %>%
  mutate(z = as.numeric(scale(median_dist))) %>%
  left_join(df_lib %>% dplyr::select(sample, sample_name, tot_counts, n_genes_detected), by = "sample") %>%
  arrange(desc(z))
# the z-score is relative to the samples left: someone is always close to 2, look at the median distance too
print(df_outlier_f)
write_tsv(df_outlier_f, file.path(out_dir, "table", "explore_outlier_score_final.tsv"))

# PCA ---------------------------------------------------------------------
pca_f <- prcomp(t(assay(vds_f)))
pv_f <- 100 * pca_f$sdev^2 / sum(pca_f$sdev^2)

tibble(PC = factor(seq_along(pv_f)), var = pv_f) %>%
  ggplot(aes(PC, var)) + geom_col() + labs(y = "% variance") + theme_bw()
ggsave(file.path(plot_dir_final, "PCA_scree.pdf"), width = 5, height = 3.5)

df_pca_f <- cbind(as.data.frame(colData(vds_f)[, c("sample", "sample_name", vars_f)]),
                  pca_f$x[, seq_len(n_pc)]) %>%
  left_join(df_lib %>% dplyr::select(sample, tot_counts, n_genes_detected), by = "sample") %>%
  mutate(log10_tot_counts = log10(tot_counts))

# pairs of PCs, one page per covariate
pcs_plot_f <- paste0("PC", 1:4)
p_pca_f <- crossing(covariate = intersect(c("treat", "clone", "mut"), vars_f),
                    PCx = pcs_plot_f,
                    PCy = pcs_plot_f) %>%
  filter(PCx < PCy) %>%
  group_split(covariate) %>%
  map(function(df_pairs) {
    df_pairs %>%
      pmap(function(covariate, PCx, PCy) {
        df_pca_f %>%
          ggplot(aes(x = .data[[PCx]], y = .data[[PCy]], col = .data[[covariate]], label = sample)) +
          geom_point(size = 3, alpha = 0.8) +
          ggrepel::geom_text_repel(show.legend = FALSE, size = 3) +
          labs(x = sprintf("%s: %.0f%%", PCx, pv_f[parse_number(PCx)]),
               y = sprintf("%s: %.0f%%", PCy, pv_f[parse_number(PCy)]),
               col = covariate) +
          theme_bw() + theme(strip.background = element_blank())
      }) %>%
      wrap_plots(nrow = 2, guides = "collect") +
      plot_annotation(title = paste("PCA coloured by", unique(df_pairs$covariate), "(samples after removal)"))
  })
pdf(file.path(plot_dir_final, "PCA_byVariable.pdf"), width = 16, height = 8)
walk(p_pca_f, print)
dev.off()

# which metadata explain each PC (R2 of a linear model, * = ANOVA p < 0.05)
df_pc_assoc_f <- crossing(PC = paste0("PC", seq_len(n_pc)),
                          variable = c(vars_f, "log10_tot_counts", "n_genes_detected")) %>%
  mutate(res = map2(PC, variable, function(pc, v) {
    fit <- lm(df_pca_f[[pc]] ~ df_pca_f[[v]])
    tibble(R2 = summary(fit)$r.squared, p = anova(fit)$`Pr(>F)`[1])
  })) %>%
  unnest(res) %>%
  mutate(PC = factor(PC, levels = paste0("PC", seq_len(n_pc))),
         label = sprintf("%.2f%s", R2, if_else(p < 0.05, "*", "")))
write_tsv(df_pc_assoc_f, file.path(out_dir, "table", "explore_PC_vs_metadata_final.tsv"))

ggplot(df_pc_assoc_f, aes(PC, variable, fill = R2, label = label)) +
  geom_tile() +
  geom_shadowtext(size = 3, color = "white", bg.color = "black", bg.r = 0.15) +
  scale_fill_viridis_c(limits = c(0, 1)) +
  labs(title = "Variance of each PC explained by each variable (R2, * = p < 0.05), samples after removal",
       x = NULL, y = NULL) +
  theme_minimal()
ggsave(file.path(plot_dir_final, "PCA_PC_vs_metadata_R2.pdf"), width = 7, height = 4)

# clone removed (visualisation only) --------------------------------------
# remove the clone effect keeping treat + mut: if the contrasts are real they should separate the samples
cd_vds_f <- as.data.frame(colData(vds_f)) %>%
  mutate(treat = relevel(treat, ref = "BL"), mut = relevel(mut, ref = "WT"))
mat_batch_f <- limma::removeBatchEffect(assay(vds_f), batch = cd_vds_f$clone, design = model.matrix(~ treat + mut, cd_vds_f))
  
pca_batch_f <- prcomp(t(mat_batch_f))
pv_batch_f <- 100 * pca_batch_f$sdev^2 / sum(pca_batch_f$sdev^2)
cbind(cd_vds_f[, c("sample", "clone", "mut", "treat")], pca_batch_f$x[, 1:3]) %>%
  ggplot(aes(PC1, PC2, col = treat, shape = mut, label = sample)) +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(show.legend = FALSE, size = 3) +
  labs(x = sprintf("PC1: %.0f%%", pv_batch_f[1]), y = sprintf("PC2: %.0f%%", pv_batch_f[2]),
       title = "PCA after removing the clone effect (visualisation only)") +
  theme_bw()
ggsave(file.path(plot_dir_final, "PCA_cloneRemoved.pdf"), width = 6, height = 5)

# sex check on the samples left -------------------------------------------
# the sex score and the vst are the ones computed above on the full set (the Y genes are not filtered out)
df_sex %>%
  filter(sample %in% rownames(coldata_f)) %>%
  group_by(clone) %>%
  summarise(n = n(), min_score = min(sex_score), max_score = max(sex_score),
            range = max_score - min_score, .groups = "drop") %>%
  print()

mat_sex_f <- mat_sex[, rownames(coldata_f), drop = FALSE]
mat_sex_f_z <- t(scale(t(mat_sex_f)))
mat_sex_f_z[is.na(mat_sex_f_z)] <- 0
colnames(mat_sex_f_z) <- coldata_f$sample_name
ha_sex_f <- HeatmapAnnotation(clone = coldata_f$clone,
                              mut = coldata_f$mut,
                              treat = coldata_f$treat,
                              col = list(clone = col_clone, mut = col_mut, treat = col_treat))
pdf(file.path(plot_dir_final, "heatmap_sex_genes.pdf"), width = 9, height = 6)
draw(Heatmap(mat_sex_f_z, name = "z-score vst",
             col = circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
             top_annotation = ha_sex_f,
             cluster_row_slices = FALSE,
             column_title = "sex genes (XIST: female; Y-linked: male), samples after removal"),
     heatmap_legend_side = "left", annotation_legend_side = "left")
dev.off()

# p-value histograms of the two contrasts (sanity check only) --------------
# ~ clone + treat + mut on the raw counts: a flat histogram with a peak near 0 = calibrated test with some signal;
# the DE used downstream is defined in 01_import_data.R / 03_differential_analysis.R
dds_de <- DESeqDataSetFromMatrix(countData = counts_f[feat_keep_f, ],
                                 colData = cd_design,
                                 design = ~ clone + treat + mut)
dds_de <- DESeq(dds_de, quiet = TRUE)
message("coefficients: ", paste(resultsNames(dds_de), collapse = ", "))

res_f <- list(KO_vs_WT = results(dds_de, name = "mut_KO_vs_WT"),
              Myelin_vs_BL = results(dds_de, name = "treat_Myelin_vs_BL"))

map_dfr(res_f, function(x) tibble(n_tested = sum(!is.na(x$pvalue)), n_padj_0.05 = sum(x$padj < 0.05, na.rm = TRUE)),
        .id = "contrast") %>%
  print()

map_dfr(res_f, ~ as_tibble(as.data.frame(.x)), .id = "contrast") %>%
  filter(!is.na(pvalue)) %>%
  ggplot(aes(pvalue)) +
  geom_histogram(bins = 20, boundary = 0, fill = "grey60", col = "white") +
  facet_wrap(~ contrast) +
  labs(title = "p-value distribution (~ clone + treat + mut)", x = "p-value", y = "genes") +
  theme_bw() + theme(strip.background = element_blank())
ggsave(file.path(plot_dir_final, "pvalue_hist_contrasts.pdf"), width = 7, height = 3.5)

res_f$KO_vs_WT %>%
  data.frame() %>%
  arrange(padj) %>%
  rownames_to_column("gene") %>%
  mutate(rank = rank(pvalue)) %>%
  filter(str_detect(gene,pattern = "C1Q"))

res_f$Myelin_vs_BL %>%
  data.frame() %>%
  arrange(padj) %>%
  rownames_to_column("gene") %>%
  mutate(rank = rank(pvalue))

# save --------------------------------------------------------------------
saveRDS(vds_f, file.path(out_dir, "object", "vst_blind_explore_final.rds"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "table", "sessionInfo_explore.txt"))
