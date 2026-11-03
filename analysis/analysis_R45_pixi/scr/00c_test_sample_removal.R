# AIM ---------------------------------------------------------------------
# Test the effect of removing different sets of samples BEFORE deciding the final selection.
# Each scenario (a set of samples to drop, see `scenarios`) is run through the same steps and the results are compared:
#   1. replication of the conditions / estimable designs (is the contrast still testable?)
#   2. blind vst + PCA, and which metadata explain each PC
#   3. DE with ~ clone + treat + mut (only to compare the scenarios: number of DE genes, overlap, p-value distribution)
# The DE is run on the genes that pass the filter in ALL the scenarios, so that the number of genes tested
# (and the multiple testing correction) is the same and the scenarios are comparable.
# Nothing here changes the dds used downstream.
# Run from the project root (analysis_R45_pixi), paths are relative to it.

# libraries ---------------------------------------------------------------
library(tidyverse)
library(DESeq2)
library(edgeR)
library(shadowtext)

# parameters --------------------------------------------------------------
data_dir <- "../../data"
out_dir <- "../../out"

# the scenarios to compare: name = samples to drop. Add / edit freely (e.g. drop_6A8 = "6A8").
scenarios <- list(
  all_samples = character(0),
  drop_2 = c("24A8", "6A8"),
  drop_4 = c("24A8", "6A8", "6B3", "6B7")
)

vars_meta <- c("clone", "mut", "time", "treat")
filter_group <- c("mut", "treat") # group used by filterByExpr
n_pc <- 6
padj_cutoff <- 0.05

plot_dir <- file.path(out_dir, "plot", "explore_removal")
table_dir <- file.path(out_dir, "table", "removal_test")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

# read in the data --------------------------------------------------------
counts <- readRDS(file.path(out_dir, "object", "raw_count_salmon_GTF.rds"))
samplesheet <- read_csv(file.path(data_dir, "samplesheet.csv"), show_col_types = FALSE)

coldata_all <- data.frame(sample = colnames(counts)) %>%
  left_join(samplesheet, by = c("sample" = "sample_id")) %>%
  mutate(across(all_of(vars_meta), factor),
         sample_name = paste(sample, clone, mut, time, treat, sep = "_")) %>%
  column_to_rownames("sample") %>%
  mutate(sample = rownames(.))

# check that the samples to drop exist
unknown <- setdiff(unlist(scenarios), rownames(coldata_all))
if (length(unknown) > 0) stop("samples not in the count table: ", paste(unknown, collapse = ", "))

counts <- counts[, rownames(coldata_all)]
storage.mode(counts) <- "integer"

# data of each scenario ---------------------------------------------------
# reference levels as in the DE: BL for treat, WT for mut
sc <- imap(scenarios, function(drop, nm) {
  cd <- coldata_all %>%
    filter(!sample %in% drop) %>%
    mutate(across(where(is.factor), droplevels))
  cd <- cd %>% mutate(treat = relevel(treat, ref = "BL"), mut = relevel(mut, ref = "WT"))
  cnt <- counts[, rownames(cd)]
  keep <- edgeR::filterByExpr(cnt, group = do.call(paste, c(cd[filter_group], sep = "|")))
  message(nm, ": ", ncol(cnt), " samples, ", sum(keep), " genes after the filter")
  list(coldata = cd, counts = cnt, keep = keep)
})

genes_common <- Reduce(intersect, map(sc, ~ rownames(.x$counts)[.x$keep]))
message("genes passing the filter in all the scenarios: ", length(genes_common))

# 1. replication and estimable designs -------------------------------------
# samples per condition in each scenario
df_rep <- imap_dfr(sc, ~ .x$coldata %>% dplyr::count(mut, treat), .id = "scenario")
df_rep %>%
  pivot_wider(names_from = c(mut, treat), values_from = n, values_fill = 0) %>%
  print(width = Inf)
write_tsv(df_rep, file.path(table_dir, "replication_by_scenario.tsv"))

# rank = number of columns and residual df > 0
# (the interaction is expected to be not estimable if a mut x treat combination is missing)
candidate_designs <- list(
  treat_mut = ~ treat + mut,
  clone_treat_mut = ~ clone + treat + mut,
  mut_x_treat = ~ mut * treat,
  clone_mut_x_treat = ~ clone + mut * treat
)

df_designs <- imap_dfr(sc, function(x, nm) {
  imap_dfr(candidate_designs, function(f, dn) {
    m <- model.matrix(f, x$coldata)
    tibble(design = dn, formula = deparse(f), n_samples = nrow(m), n_coef = ncol(m),
           rank = qr(m)$rank, residual_df = nrow(m) - qr(m)$rank,
           full_rank = qr(m)$rank == ncol(m))
  })
}, .id = "scenario")
print(df_designs, n = Inf, width = Inf)
write_tsv(df_designs, file.path(table_dir, "candidate_designs_by_scenario.tsv"))

# 2. vst + PCA --------------------------------------------------------------
# the vst and the PCA use the genes passing the filter of each scenario
pca_sc <- imap(sc, function(x, nm) {
  dds <- DESeqDataSetFromMatrix(countData = x$counts[x$keep, ], colData = x$coldata, design = ~ 1)
  vds <- vst(dds, blind = TRUE)
  pca <- prcomp(t(assay(vds)))
  list(pca = pca, pv = 100 * pca$sdev^2 / sum(pca$sdev^2))
})

df_pca <- imap_dfr(sc, function(x, nm) {
  pca <- pca_sc[[nm]]$pca
  cbind(x$coldata[, c("sample", "sample_name", vars_meta)], pca$x[, seq_len(n_pc)]) %>%
    mutate(tot_counts = colSums(x$counts)[sample],
           n_genes_detected = colSums(x$counts >= 10)[sample],
           log10_tot_counts = log10(tot_counts))
}, .id = "scenario") %>%
  mutate(scenario = factor(scenario, levels = names(sc)))

# PC1 vs PC2 of each scenario (the % of variance in the strip)
df_pv <- imap_dfr(pca_sc, ~ tibble(PC1 = .x$pv[1], PC2 = .x$pv[2]), .id = "scenario") %>%
  mutate(facet = sprintf("%s (PC1 %.0f%%, PC2 %.0f%%)", scenario, PC1, PC2))

df_pca %>%
  left_join(df_pv %>% dplyr::select(scenario, facet), by = "scenario") %>%
  mutate(facet = factor(facet, levels = df_pv$facet)) %>%
  ggplot(aes(PC1, PC2, col = treat, shape = mut, label = sample)) +
  geom_point(size = 3, alpha = 0.8) +
  ggrepel::geom_text_repel(show.legend = FALSE, size = 3) +
  facet_wrap(~ facet, scales = "free") +
  theme_bw() + theme(strip.background = element_blank())
ggsave(file.path(plot_dir, "PCA_scenarios.pdf"), width = 5 * length(sc), height = 5)

# same PCA coloured by clone (the clone is expected to be the main source of variation once the outliers are out)
df_pca %>%
  left_join(df_pv %>% dplyr::select(scenario, facet), by = "scenario") %>%
  mutate(facet = factor(facet, levels = df_pv$facet)) %>%
  ggplot(aes(PC1, PC2, col = clone, shape = mut, label = sample)) +
  geom_point(size = 3, alpha = 0.8) +
  ggrepel::geom_text_repel(show.legend = FALSE, size = 3) +
  facet_wrap(~ facet, scales = "free") +
  theme_bw() + theme(strip.background = element_blank())
ggsave(file.path(plot_dir, "PCA_scenarios_clone.pdf"), width = 5 * length(sc), height = 5)

# which metadata explain each PC (R2 of a linear model, * = ANOVA p < 0.05); constant variables are skipped
df_pc_assoc <- imap_dfr(sc, function(x, nm) {
  df_s <- df_pca %>% filter(scenario == nm)
  vars_s <- c(vars_meta[map_lgl(x$coldata[vars_meta], ~ nlevels(.x) > 1)], "log10_tot_counts", "n_genes_detected")
  crossing(PC = paste0("PC", seq_len(n_pc)), variable = vars_s) %>%
    mutate(res = map2(PC, variable, function(pc, v) {
      fit <- lm(df_s[[pc]] ~ df_s[[v]])
      tibble(R2 = summary(fit)$r.squared, p = anova(fit)$`Pr(>F)`[1])
    })) %>%
    unnest(res)
}, .id = "scenario") %>%
  mutate(PC = factor(PC, levels = paste0("PC", seq_len(n_pc))),
         scenario = factor(scenario, levels = names(sc)),
         label = sprintf("%.2f%s", R2, if_else(p < 0.05, "*", "")))
write_tsv(df_pc_assoc, file.path(table_dir, "PC_vs_metadata_by_scenario.tsv"))

ggplot(df_pc_assoc, aes(PC, variable, fill = R2, label = label)) +
  geom_tile() +
  geom_shadowtext(size = 3, color = "white", bg.color = "black", bg.r = 0.15) +
  scale_fill_viridis_c(limits = c(0, 1)) +
  facet_wrap(~ scenario, nrow = 1) +
  labs(title = "Variance of each PC explained by each variable (R2, * = p < 0.05)", x = NULL, y = NULL) +
  theme_minimal()
ggsave(file.path(plot_dir, "PCA_PC_vs_metadata_R2_scenarios.pdf"), width = 5 * length(sc), height = 4)

# 3. DE: ~ clone + treat + mut ----------------------------------------------
# same genes in all the scenarios; a scenario where the design is not estimable is skipped
res_sc <- imap(sc, function(x, nm) {
  tryCatch({
    dds <- DESeqDataSetFromMatrix(countData = x$counts[genes_common, ], colData = x$coldata,
                                  design = ~ clone + treat + mut)
    dds <- DESeq(dds, quiet = TRUE)
    list(KO_vs_WT = results(dds, contrast = c("mut", "KO", "WT")),
         Myelin_vs_BL = results(dds, contrast = c("treat", "Myelin", "BL")))
  }, error = function(e) {
    message(nm, ": DE skipped - ", conditionMessage(e))
    NULL
  })
}) %>%
  compact()

df_de <- imap_dfr(res_sc, function(res, nm) {
  imap_dfr(res, ~ as.data.frame(.x) %>% rownames_to_column("gene") %>% as_tibble(), .id = "contrast")
}, .id = "scenario") %>%
  mutate(scenario = factor(scenario, levels = names(sc)))

df_n_de <- df_de %>%
  group_by(scenario, contrast) %>%
  summarise(n_tested = sum(!is.na(pvalue)),
            n_padj = sum(padj < padj_cutoff, na.rm = TRUE),
            n_up = sum(padj < padj_cutoff & log2FoldChange > 0, na.rm = TRUE),
            n_down = sum(padj < padj_cutoff & log2FoldChange < 0, na.rm = TRUE),
            .groups = "drop")
print(df_n_de)

ggplot(df_n_de, aes(scenario, n_padj, fill = scenario, label = n_padj)) +
  geom_col(show.legend = FALSE) +
  geom_text(vjust = -0.3) +
  facet_wrap(~ contrast) +
  labs(y = sprintf("genes with padj < %s", padj_cutoff), x = NULL) +
  theme_bw() + theme(strip.background = element_blank())
ggsave(file.path(plot_dir, "n_DE_genes_scenarios.pdf"), width = 7, height = 4)

# distribution of the p-values: flat + peak near 0 = calibrated test with some signal
df_de %>%
  filter(!is.na(pvalue)) %>%
  ggplot(aes(pvalue)) +
  geom_histogram(bins = 20, boundary = 0, fill = "grey60", col = "white") +
  facet_grid(scenario ~ contrast) +
  labs(x = "p-value", y = "genes") +
  theme_bw() + theme(strip.background = element_blank())
ggsave(file.path(plot_dir, "pvalue_hist_scenarios.pdf"), width = 7, height = 2.5 * length(sc))

# overlap of the DE genes between the scenarios
sig_sets <- df_de %>%
  filter(!is.na(padj), padj < padj_cutoff) %>%
  group_by(contrast, scenario) %>%
  summarise(genes = list(gene), .groups = "drop")

df_overlap <- sig_sets %>%
  inner_join(sig_sets, by = "contrast", suffix = c("_a", "_b"), relationship = "many-to-many") %>%
  filter(as.integer(scenario_a) < as.integer(scenario_b)) %>%
  mutate(n_a = lengths(genes_a), n_b = lengths(genes_b),
         n_shared = map2_int(genes_a, genes_b, ~ length(intersect(.x, .y)))) %>%
  dplyr::select(contrast, scenario_a, scenario_b, n_a, n_b, n_shared)
print(df_overlap)

# log2FC of each scenario against the first one (the most complete): genes moving with the removal are off the diagonal
ref_sc <- names(res_sc)[1]
df_de %>%
  dplyr::select(scenario, contrast, gene, log2FoldChange) %>%
  pivot_wider(names_from = scenario, values_from = log2FoldChange) %>%
  pivot_longer(-c(contrast, gene, all_of(ref_sc)), names_to = "scenario", values_to = "lfc") %>%
  ggplot(aes(.data[[ref_sc]], lfc)) +
  geom_point(size = 0.4, alpha = 0.3) +
  geom_abline(col = "red") +
  coord_cartesian(xlim = c(-5, 5), ylim = c(-5, 5)) +
  facet_grid(scenario ~ contrast) +
  labs(x = paste("log2FC,", ref_sc), y = "log2FC, scenario") +
  theme_bw() + theme(strip.background = element_blank())
ggsave(file.path(plot_dir, "log2FC_vs_reference_scenario.pdf"), width = 7, height = 2.8 * length(sc))

# summary of the scenarios --------------------------------------------------
df_summary <- tibble(scenario = names(sc),
                     dropped = map_chr(scenarios, ~ if (length(.x) == 0) "-" else paste(.x, collapse = ",")),
                     n_samples = map_int(sc, ~ nrow(.x$coldata)),
                     n_genes_filter = map_int(sc, ~ sum(.x$keep))) %>%
  left_join(df_designs %>%
              filter(design == "clone_treat_mut") %>%
              dplyr::select(scenario, residual_df, full_rank), by = "scenario") %>%
  left_join(df_n_de %>%
              dplyr::select(scenario, contrast, n_padj) %>%
              mutate(scenario = as.character(scenario)) %>%
              pivot_wider(names_from = contrast, values_from = n_padj, names_prefix = "n_DE_"),
            by = "scenario") %>%
  left_join(df_pv %>% dplyr::select(scenario, PC1_var = PC1, PC2_var = PC2), by = "scenario") %>%
  left_join(df_pc_assoc %>%
              filter(PC == "PC1") %>%
              group_by(scenario) %>%
              slice_max(R2, n = 1, with_ties = FALSE) %>%
              ungroup() %>%
              transmute(scenario = as.character(scenario), PC1_main_var = variable, PC1_R2 = R2),
            by = "scenario")
print(df_summary, width = Inf)
write_tsv(df_summary, file.path(table_dir, "summary_scenarios.tsv"))
write_tsv(df_n_de, file.path(table_dir, "n_DE_genes_by_scenario.tsv"))
write_tsv(df_overlap, file.path(table_dir, "DE_overlap_scenarios.tsv"))

writeLines(capture.output(sessionInfo()), file.path(table_dir, "sessionInfo_removal_test.txt"))
