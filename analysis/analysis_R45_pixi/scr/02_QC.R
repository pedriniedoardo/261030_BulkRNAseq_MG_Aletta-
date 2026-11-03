# AIM ---------------------------------------------------------------------
# run some preprocessing and QC of the dataset to understand if there are outliers or global trends.

# libraries ---------------------------------------------------------------
library(tidyverse)
library(DESeq2)
library(ggrepel)
library(vsn)
library(hexbin)
library(viridis)
library(pheatmap)
library(PoiClaClu)
library(AnnotationDbi) 
library(AnnotationHub)
library(GGally)
library(RNAseqQC)
library(edgeR)
library(limma)
library(ComplexHeatmap)
library(pals)

# read in the data --------------------------------------------------------
ddsHTSeq <- readRDS("../../out/object/dds.rds")

# remove low expressed genes ----------------------------------------------
# To reduce the noise filter out the lowly expressed genes
# plot the total counts per sample
colSums(counts(ddsHTSeq)) %>%
  data.frame(tot_counts=.) %>%
  rownames_to_column("sample") %>% 
  # add the metadata
  left_join(colData(ddsHTSeq) %>% 
              data.frame(),by = "sample") %>% 
  ggplot(aes(x=sample_name,y = tot_counts))+geom_col()+theme_bw()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        plot.margin = margin(0.5, 0.5, 2, 2, "cm"))
ggsave("../../out/plot/barplot_tot_count.pdf",width = 7,height = 5)

# check the total size of the matrix before filtering fow low expressed genes
nrow(ddsHTSeq)

# remove potential non infirmative genes (lowly expressed genes).
# use the automatic implementation provided by edgeR
ddsHTSeq_filter <- ddsHTSeq[edgeR::filterByExpr(counts(ddsHTSeq), group = paste(colData(ddsHTSeq)$mut, colData(ddsHTSeq)$treat, sep = "|")),] # group = mut x treat, as in 00b_explore_QC.R

# this is the old implementaion 
# ddsHTSeq_filter <- ddsHTSeq[rowSums(counts(ddsHTSeq)) > 100, ]

# check the size of the dataset after filtering for lowly expressed genes
nrow(ddsHTSeq_filter)

# save the filtered object
saveRDS(ddsHTSeq_filter,file = "../../out/object/dds_filter.rds")

# scaling transformation of the data --------------------------------------
# vsd
vds_filter <- vst(ddsHTSeq_filter, blind = T)
head(assay(vds_filter), 3)

# calculate the vsd also for the whole dataset
vds_unfiltered <- vst(ddsHTSeq, blind = T)

# rlog
rld_filter <- rlog(ddsHTSeq_filter, blind = T)
head(assay(rld_filter), 3)

# diagnostic plot
meanSdPlot_vsd <- meanSdPlot(assay(vds_filter))
ggsave(plot = meanSdPlot_vsd$gg+theme_bw(),filename = "../../out/plot/meanSdPlot_vsd.pdf",width = 4,height = 4)

meanSdPlot_rlog <- meanSdPlot(assay(rld_filter))
ggsave(plot = meanSdPlot_rlog$gg+theme_bw(),filename = "../../out/plot/meanSdPlot_rlog.pdf",width = 4,height = 4)

ddsHTSeq_filter <- estimateSizeFactors(ddsHTSeq_filter)

# sample distance ---------------------------------------------------------
# using vsd scaling
sampleDists_vsd <- dist(t(assay(vds_filter)))
sampleDists_vsd

head(assay(vds_filter))

sampleDistMatrix_vsd <- as.matrix(sampleDists_vsd)

rownames(sampleDistMatrix_vsd) <- paste(vds_filter$sample_name)
colnames(sampleDistMatrix_vsd) <- NULL

map_colors<-colorRampPalette(viridis(12))(255)

hm_1 <- pheatmap(sampleDistMatrix_vsd,
                 clustering_distance_rows = sampleDists_vsd,
                 clustering_distance_cols = sampleDists_vsd,
                 col = map_colors)

pdf("../../out/plot/heatmap_vsd.pdf",width = 7,height = 5)
hm_1
dev.off()

# using rlog scaling
sampleDists_rld <- dist(t(assay(rld_filter)))
sampleDists_rld

head(assay(rld_filter))

sampleDistMatrix_rld <- as.matrix(sampleDists_rld)

rownames(sampleDistMatrix_rld) <- paste(rld_filter$sample_name)
colnames(sampleDistMatrix_rld) <- NULL

hm_1_2 <- pheatmap(sampleDistMatrix_rld,
                   clustering_distance_rows = sampleDists_rld,
                   clustering_distance_cols = sampleDists_rld,
                   col = map_colors)

pdf("../../out/plot/heatmap_rld.pdf",width = 7,height = 5)
hm_1_2
dev.off()

# using poisson distance
poisd <- PoissonDistance(t(counts(ddsHTSeq_filter,normalized = F)))

samplePoisDistMatrix <- as.matrix(poisd$dd)

rownames(samplePoisDistMatrix) <- paste(ddsHTSeq_filter$sample_name)
colnames(samplePoisDistMatrix) <- NULL

hm_p <- pheatmap(samplePoisDistMatrix,
                 clustering_distance_rows = poisd$dd,
                 clustering_distance_cols = poisd$dd,
                 col = map_colors)

pdf("../../out/plot/heatmap_poisd.pdf",width = 7,height = 5)
hm_p
dev.off()

# plot cluster alternative ------------------------------------------------
# set seed to control random annotation colors
pdf("../../out/plot/heatmap_cluster.pdf",width = 7,height = 5)
set.seed(2)
plot_sample_clustering(vds_filter,
                       anno_vars = c("clone","treat","mut","time"),
                       distance = "euclidean")
dev.off()

# PCA plot ----------------------------------------------------------------
plot_vsd <- plotPCA(vds_filter,
                    intgroup = c("clone","treat","mut","time")) +
  theme_bw()

plot_vsd$data %>%
  ggplot(aes(x=PC1,y=PC2,label = name)) +
  # ggplot(aes(x=PC1,y=PC2,col=BMP_treat,shape=condition)) +
  geom_point(aes(col=treat,shape=mut),size =3) +
  ggrepel::geom_text_repel(show.legend = F)+
  scale_x_continuous(expand = expansion(mult = 0.1))+
  theme_bw() + ylab(plot_vsd$labels[1]) + xlab(plot_vsd$labels[2])
ggsave("../../out/plot/PCA_vsd_colorTreat.pdf",width = 5,height = 4)

# pull more PC
rv <- rowVars(assay(vds_filter),useNames = TRUE)
# select the ntop genes by variance
select_var <- order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))]
# perform a PCA on the data in assay(x) for the selected genes
test <- prcomp(t(assay(vds_filter)[select_var,]))$x %>% 
  data.frame() %>% 
  rownames_to_column("sample")

# plot more PC by condition
left_join(plot_vsd$data %>% dplyr::select(-c("PC1","PC2")),test,by = c("name"="sample")) %>%
  # ggpairs(columns = 5:14,ggplot2::aes(colour=condition),upper = "blank")+
  ggpairs(columns = paste0("PC",1:10),ggplot2::aes(colour=treat),upper = "blank") +
  theme_bw()+theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))
ggsave("../../out/plot/PCA_panel_vsd_treat.pdf",width = 20,height = 20)

# explore pc score by metadata fo the samples
test2 <- left_join(plot_vsd$data %>% dplyr::select(-c("PC1","PC2")),test,by = c("name"="sample"))

test_df1 <- test2 |> 
  dplyr::select(clone,treat,mut,time,name) |> 
  # dplyr::rename(approx.time = Approx.Time.between.collection.and.processing) |> 
  pivot_longer(names_to = "var_1",values_to = "value_1",-c(name)) %>%
  filter(var_1 != "time")

test_df2 <- test2 |> 
  dplyr::select(name,PC1:PC9) |>
  pivot_longer(names_to = "var_2",values_to = "value_2",-c(name))

left_join(test_df1,test_df2,by=c("name")) |> 
  mutate(comparison = paste0(var_1,"_vs_",var_2)) |> 
  ggplot(aes(x=value_1,y=value_2)) +
  facet_wrap(~comparison,scales = "free",ncol=9) +
  # facet_grid(var_1~var_2,scales = "free_x",drop = T) +
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.2,height = 0))+
  theme_bw()+theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))
ggsave("../../out/plot/panel_metadata_PC.pdf",width = 20,height = 8)

# plot PCA also fro rlog scaled data
plot_rld <- plotPCA(rld_filter,
                    intgroup = c("clone","treat","mut","time")) +
  theme_bw()

plot_rld$data %>%
  ggplot(aes(x=PC1,y=PC2,label = name)) +
  # ggplot(aes(x=PC1,y=PC2,col=BMP_treat,shape=condition)) +
  geom_point(aes(col=treat,shape=mut),size =3) +
  ggrepel::geom_text_repel(show.legend = F)+
  scale_x_continuous(expand = expansion(mult = 0.1))+
  theme_bw() + ylab(plot_vsd$labels[1]) + xlab(plot_vsd$labels[2])
ggsave("../../out/plot/PCA_rld_colorTreat.pdf",width = 5,height = 4)

# MSD plot ----------------------------------------------------------------
mds_vsd <- as.data.frame(colData(vds_filter)) %>%
  cbind(cmdscale(sampleDistMatrix_vsd))

ggplot(mds_vsd, aes(x = `1`, y = `2`, color = treat)) +
  geom_point(size = 3)+ theme_bw()
ggsave("../../out/plot/MDS_vsd.pdf",width = 5,height = 4)

mds_rld <- as.data.frame(colData(rld_filter)) %>%
  cbind(cmdscale(sampleDistMatrix_rld))

ggplot(mds_rld, aes(x = `1`, y = `2`, color = treat)) +
  geom_point(size = 3) + theme_bw()
ggsave("../../out/plot/MDS_rld.pdf",width = 5,height = 4)

mdsPois <- as.data.frame(colData(ddsHTSeq_filter)) %>%
  cbind(cmdscale(samplePoisDistMatrix))

ggplot(mdsPois, aes(x = `1`, y = `2`, color = treat)) +
  geom_point(size = 3) + theme_bw()
ggsave("../../out/plot/PoissonDistance_scatter.pdf",width = 5,height = 4)

# sex assessment ----------------------------------------------------------
# the samplesheet has no sex: infer it from XIST (female) and the Y-linked genes (male) and check that the samples of
# the same clone agree (a clone cannot change sex: discordant samples can be a swap or a contamination).
# vst of ALL the genes (vds_unfiltered): the Y-linked genes are lowly expressed in females and could be removed by the filter.
sex_genes <- c(XIST = "female", TSIX = "female",
               DDX3Y = "male", RPS4Y1 = "male", USP9Y = "male", UTY = "male", EIF1AY = "male", ZFY = "male")

annotations <- readRDS("../../out/object/annotation_GTF.rds") %>%
  dplyr::select(GENEID, SYMBOL) %>%
  distinct()

df_sex_genes <- annotations %>%
  dplyr::filter(SYMBOL %in% names(sex_genes), GENEID %in% rownames(vds_unfiltered)) %>%
  mutate(sex = unname(sex_genes[SYMBOL])) %>%
  arrange(factor(sex, levels = c("female", "male")), SYMBOL)
message("sex genes found: ", paste(df_sex_genes$SYMBOL, collapse = ", "))

mat_sex <- assay(vds_unfiltered)[df_sex_genes$GENEID, , drop = FALSE]
rownames(mat_sex) <- df_sex_genes$SYMBOL

# # score per sample: mean vst of the Y genes minus the vst of XIST (> 0 male-like, < 0 female-like)
# df_sex <- tibble(sample = colnames(mat_sex),
#                  XIST = mat_sex["XIST", ],
#                  Y_mean = colMeans(mat_sex[df_sex_genes$sex == "male", , drop = FALSE])) %>%
#   mutate(sex_score = Y_mean - XIST) %>%
#   left_join(as.data.frame(colData(vds_unfiltered)) %>% dplyr::select(sample, sample_name, clone, mut, time, treat),
#             by = "sample") %>%
#   arrange(clone, sex_score)

# # score range within each clone: a large range = samples of the same clone with a different inferred sex
# df_sex %>%
#   group_by(clone) %>%
#   summarise(n = n(), min_score = min(sex_score), max_score = max(sex_score),
#             range = max_score - min_score, .groups = "drop") %>%
#   print()
# write_tsv(df_sex, "../../out/table/sex_score.tsv")

# heatmap: z-score of the vst per gene across the samples
mat_sex_z <- t(scale(t(mat_sex)))
mat_sex_z[is.na(mat_sex_z)] <- 0 # genes with no variation
cd_sex <- as.data.frame(colData(vds_unfiltered))[colnames(mat_sex_z), ]
colnames(mat_sex_z) <- cd_sex$sample_name

column_ha <- HeatmapAnnotation(clone = cd_sex$clone,
                               mut = cd_sex$mut,
                               treat = cd_sex$treat,
                               col = list(clone = setNames(pals::glasbey(nlevels(cd_sex$clone)), levels(cd_sex$clone)),
                                          mut = setNames(pals::okabe(nlevels(cd_sex$mut)), levels(cd_sex$mut)),
                                          treat = setNames(pals::polychrome(nlevels(cd_sex$treat)), levels(cd_sex$treat))))
pdf("../../out/plot/heatmap_sex_genes.pdf", width = 9, height = 6)
draw(Heatmap(mat_sex_z, name = "z-score vst",
             col = circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
             top_annotation = column_ha,
             cluster_row_slices = FALSE,
             column_title = "sex genes (XIST: female; Y-linked: male)"),
     heatmap_legend_side = "left", annotation_legend_side = "left")
dev.off()

# PCA after removing the clone effect -------------------------------------
# remove the clone effect keeping treat + mut: if the contrasts are real they should separate the samples (visualisation only, the DE uses clone in the design instead)
cd_vds <- as.data.frame(colData(vds_filter))
mat_batch <- limma::removeBatchEffect(assay(vds_filter),
                                      batch = cd_vds$clone,
                                      design = model.matrix(~ treat + mut, cd_vds))
pca_batch <- prcomp(t(mat_batch))
pv_batch <- 100 * pca_batch$sdev^2 / sum(pca_batch$sdev^2)

cbind(cd_vds[, c("sample", "clone", "mut", "treat")], pca_batch$x[, 1:3]) %>%
  ggplot(aes(PC1, PC2, col = treat, shape = mut, label = sample)) +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(show.legend = FALSE, size = 3) +
  labs(x = sprintf("PC1: %.0f%%", pv_batch[1]), y = sprintf("PC2: %.0f%%", pv_batch[2]),
       title = "PCA after removing the clone effect (visualisation only)") +
  theme_bw()
ggsave("../../out/plot/PCA_vsd_cloneRemoved.pdf", width = 6, height = 5)

# explore pc score by metadata after the clone removal
# same panel as above, on the PCs of the batch corrected matrix: clone should not be associated with the PCs anymore
n_pc_batch <- min(9, ncol(pca_batch$x))

test_df1_batch <- cd_vds %>%
  dplyr::select(name = sample, clone, treat, mut) %>%
  mutate(across(-name, as.character)) %>%
  pivot_longer(names_to = "var_1", values_to = "value_1", -name)

test_df2_batch <- pca_batch$x[, seq_len(n_pc_batch)] %>%
  data.frame() %>%
  rownames_to_column("name") %>%
  pivot_longer(names_to = "var_2", values_to = "value_2", -name)

left_join(test_df1_batch, test_df2_batch, by = "name") |>
  mutate(comparison = paste0(var_1, "_vs_", var_2)) |>
  ggplot(aes(x = value_1, y = value_2)) +
  facet_wrap(~comparison, scales = "free", ncol = 9) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2,height = 0)) +
  theme_bw() + theme(strip.background = element_blank(), axis.text.x = element_text(hjust = 1, angle = 45))
ggsave("../../out/plot/panel_metadata_PC_cloneRemoved.pdf", width = 20, height = 8)

# save the object of interest ---------------------------------------------
saveRDS(vds_filter,file = "../../out/object/vds_filter.rds")
saveRDS(vds_unfiltered,file = "../../out/object/vds_unfilter.rds")
saveRDS(rld_filter,file = "../../out/object/rld_filter.rds")
