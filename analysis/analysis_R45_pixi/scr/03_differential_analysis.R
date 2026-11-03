# AIM ---------------------------------------------------------------------
# run DE analysis
# model: ~ clone + treat + mut (defined in 01_import_data.R), samples in drop_final removed.
# contrasts: Myelin vs BL and KO vs WT, taken from the coefficients of the model.

# libraries ---------------------------------------------------------------
library(tidyverse)
library(DESeq2)
library(ggrepel)
library(limma)
library(ashr)
library(ComplexHeatmap)
library(pals)

# read in the data --------------------------------------------------------
LUT_gene <- readRDS("../../out/object/annotation_GTF.rds") %>%
  dplyr::select(GENEID, SYMBOL) %>%
  distinct()

ddsHTSeq_filter <- readRDS("../../out/object/dds_filter.rds")
vds_filter <- readRDS("../../out/object/vds_filter.rds")

# differential expression analyisis ---------------------------------------
ddsHTSeq_filter <- DESeq(ddsHTSeq_filter)
saveRDS(ddsHTSeq_filter, "../../out/object/dds_filter_DESeq.rds")

resultsNames(ddsHTSeq_filter)

# contrasts ---------------------------------------------------------------
# coefficients: Intercept, clone_*_vs_*, treat_Myelin_vs_BL, mut_KO_vs_WT
list_res <- list(MyelinvsBL = results(ddsHTSeq_filter, name = "treat_Myelin_vs_BL", alpha = 0.05),
                 KOvsWT = results(ddsHTSeq_filter, name = "mut_KO_vs_WT", alpha = 0.05))

lapply(list_res, summary)

# shrinkage and tables ----------------------------------------------------
list_res_shr <- lapply(list_res, function(x) lfcShrink(ddsHTSeq_filter, res = x, type = "ashr"))

# convert the results to a table (ordered by pvalue) and add the gene symbol
list_df <- lapply(list_res, function(x) {
  x %>%
    data.frame() %>%
    rownames_to_column(var = "GENEID") %>%
    arrange(pvalue) %>%
    left_join(LUT_gene, by = "GENEID")
})

list_df_shr <- lapply(list_res_shr, function(x) {
  x %>%
    data.frame() %>%
    rownames_to_column(var = "GENEID") %>%
    arrange(pvalue) %>%
    left_join(LUT_gene, by = "GENEID")
})

pwalk(list(list_df, names(list_df)), function(x, y) write_tsv(x, file = paste0("../../out/table/res_", y, ".txt")))
pwalk(list(list_df_shr, names(list_df_shr)), function(x, y) write_tsv(x, file = paste0("../../out/table/res_", y, "_shr.txt")))

# check the genes of interest
list_df$KOvsWT %>%
  mutate(rank = 1:nrow(.)) %>%
  filter(str_detect(SYMBOL,pattern = "C1QA"))

list_df_shr$KOvsWT %>%
  mutate(rank = 1:nrow(.)) %>%
  filter(str_detect(SYMBOL,pattern = "C1QA"))

# plots -------------------------------------------------------------------
padj_cut <- 0.05
# log2FC threshold for the volcano and the DEG heatmap, one per version: the shrunk log2FC is smaller than the unshrunk one
lfc_cuts <- c(unshrunk = 0.9, shrunk = 0.75)
gene_bold <- "C1QA" # gene labelled in bold
col_dir <- c(ns = "gray60", DOWN = "blue", UP = "red")

# pvalue histogram (same pvalue in the shrunk and unshrunk tables)
pwalk(list(list_df, names(list_df)), function(x, y) {
  x %>%
    dplyr::filter(baseMean > 1) %>%
    ggplot(aes(x = pvalue)) + geom_histogram(breaks = 0:20/20) + theme_bw()
  ggsave(paste0("../../out/plot/histogram_pvalue_", y, ".pdf"), width = 4, height = 3)
})

# the volcano, MA and heatmap are produced for both the unshrunk log2FC (res_<contrast>.txt, no suffix)
# and the shrunk log2FC (res_<contrast>_shr.txt, suffix _shr); padj is the same in the two versions
list_versions <- list(unshrunk = list_df, shrunk = list_df_shr)

iwalk(list_versions, function(list_x, version) {
  suffix <- if (version == "shrunk") "_shr" else ""
  lfc_cut <- lfc_cuts[[version]]

  iwalk(list_x, function(x, y) {
    # volcano: red = upregulated, blue = downregulated, grey = not significant; the gene of interest is labelled in bold
    df <- x %>%
      mutate(direction = case_when(padj < padj_cut & log2FoldChange > lfc_cut & !is.na(SYMBOL) ~ "UP",
                                   padj < padj_cut & log2FoldChange < -lfc_cut & !is.na(SYMBOL) ~ "DOWN",
                                   TRUE ~ "ns"),
             direction = factor(direction, levels = c("ns", "DOWN", "UP"))) %>%
      dplyr::filter(!is.na(padj))
    ggplot(df, aes(x = log2FoldChange, y = -log10(padj), col = direction, alpha = direction)) +
      geom_point() +
      scale_alpha_manual(values = c("ns" = 0.05, "DOWN" = 0.5, "UP" = 0.5)) +
      scale_color_manual(values = col_dir) +
      geom_vline(xintercept = c(-lfc_cut, lfc_cut), col = "gray40", linetype = "dashed") +
      geom_hline(yintercept = -log10(padj_cut), col = "gray40", linetype = "dashed") +
      # labels of the DEGs (the gene of interest is drawn separately in bold)
      ggrepel::geom_text_repel(data = df[df$direction != "ns" & !df$SYMBOL %in% gene_bold, ],
                               aes(label = SYMBOL), size = 2, show.legend = F) +
      # gene of interest: bold, always labelled and drawn on top, whatever its significance
      geom_point(data = df[df$SYMBOL %in% gene_bold, ], alpha = 1, size = 3, show.legend = F) +
      ggrepel::geom_text_repel(data = df[df$SYMBOL %in% gene_bold, ], aes(label = SYMBOL),
                               fontface = "bold", size = 4, box.padding = 1, show.legend = F) +
      labs(title = paste0(y, " (", version, " log2FC, threshold ", lfc_cut, ")")) +
      theme_bw() + theme(legend.position = "none")
    ggsave(paste0("../../out/plot/vulcano_plot_text_", y, suffix, ".pdf"), width = 12, height = 12)

    # MA plot
    x %>%
      mutate(direction = case_when(padj < padj_cut & log2FoldChange > 0 ~ "UP",
                                   padj < padj_cut & log2FoldChange < 0 ~ "DOWN",
                                   TRUE ~ "ns"),
             direction = factor(direction, levels = c("ns", "DOWN", "UP"))) %>%
      ggplot(aes(baseMean, y = log2FoldChange, col = direction)) + geom_point(alpha = 0.2) +
      scale_x_log10() + scale_color_manual(values = col_dir) +
      labs(title = paste(y, version)) + theme_bw() +
      theme(legend.position = "none")
    ggsave(paste0("../../out/plot/MA_plot_", y, suffix, ".pdf"), width = 4, height = 3)
  })
})

# heatmap of the DEGs -----------------------------------------------------
# z-score of the vst per gene across the samples, one heatmap per contrast and per version (DEGs: padj and log2FC thresholds of the volcano);
# columns annotated with clone, mut and treat (same palettes as the sex genes heatmap in 02_QC.R)
cd <- as.data.frame(colData(vds_filter))
col_ha <- list(clone = setNames(pals::glasbey(nlevels(cd$clone)), levels(cd$clone)),
               mut = setNames(pals::okabe(nlevels(cd$mut)), levels(cd$mut)),
               treat = setNames(pals::polychrome(nlevels(cd$treat)), levels(cd$treat)))
max_row_names <- 80 # row names are shown only if the heatmap has up to this number of genes

iwalk(list_versions, function(list_x, version) {
  suffix <- if (version == "shrunk") "_shr" else ""
  lfc_cut <- lfc_cuts[[version]]

  iwalk(list_x, function(x, contrast) {
    degs <- x %>%
      dplyr::filter(!is.na(padj), padj < padj_cut, abs(log2FoldChange) > lfc_cut, !is.na(SYMBOL)) %>%
      dplyr::filter(GENEID %in% rownames(vds_filter))
    message(contrast, " (", version, " log2FC): ", nrow(degs), " DEGs (padj < ", padj_cut, ", |log2FC| > ", lfc_cut, ")")
    # at least 2 genes are needed to cluster the rows
    if (nrow(degs) < 2) return(invisible(NULL))

    mat_z <- assay(vds_filter)[degs$GENEID, , drop = FALSE] %>%
      t() %>% scale() %>% t()
    mat_z[is.na(mat_z)] <- 0 # genes with no variation
    rownames(mat_z) <- degs$SYMBOL
    colnames(mat_z) <- cd[colnames(mat_z), "sample_name"]

    pdf(paste0("../../out/plot/heatmap_DEGs_", contrast, suffix, ".pdf"),
        width = 8, height = if (nrow(degs) <= max_row_names) max(5, 0.15 * nrow(degs) + 3) else 8)
    draw(Heatmap(mat_z, name = "z-score vst",
                 col = circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
                 top_annotation = HeatmapAnnotation(clone = cd$clone, mut = cd$mut, treat = cd$treat, col = col_ha),
                 show_row_names = nrow(degs) <= max_row_names,
                 row_names_gp = grid::gpar(fontsize = 7),
                 column_title = sprintf("%s (%s log2FC): %d DEGs (padj < %s, |log2FC| > %s)", contrast, version, nrow(degs), padj_cut, lfc_cut)),
         heatmap_legend_side = "left", annotation_legend_side = "left")
    dev.off()
  })
})

# gene of interest: normalized expression --------------------------------
# DESeq2 normalized counts + 0.5 pseudocount (as in plotCounts), log10 scale, samples of the same clone connected
GOI <- LUT_gene %>%
  dplyr::filter(SYMBOL %in% gene_bold) %>%
  pull(GENEID)
GOI <- intersect(GOI, rownames(ddsHTSeq_filter)) # the gene could have been removed by the expression filter

lut <- colData(ddsHTSeq_filter) %>%
  data.frame()

df_goi <- counts(ddsHTSeq_filter, normalized = TRUE) %>%
  data.frame(check.names = FALSE) %>%
  rownames_to_column("GENEID") %>%
  dplyr::filter(GENEID %in% GOI) %>%
  pivot_longer(names_to = "sample", values_to = "count", -GENEID) %>%
  left_join(LUT_gene, by = "GENEID") %>%
  left_join(lut, by = "sample") %>%
  mutate(count_norm_adj = count + 0.5)

# # treat on the x-axis, panels split by mut
# ggplot(df_goi, aes(x = treat, y = count_norm_adj)) +
#   geom_boxplot(outlier.shape = NA) +
#   geom_line(aes(group = clone, col = clone), alpha = 0.6) +
#   geom_point(aes(col = clone), size = 2, alpha = 0.8) +
#   facet_grid(rows = vars(SYMBOL), cols = vars(mut), scales = "free_y") +
#   scale_y_log10() +
#   labs(y = "normalized counts + 0.5") +
#   theme_bw() +
#   theme(strip.background = element_blank(),
#         panel.border = element_rect(colour = "black", fill = NA))
# ggsave(paste0("../../out/plot/boxplot_", gene_bold, "_treat.pdf"), width = 6, height = 4)

# mut on the x-axis, panels split by treat
ggplot(df_goi, aes(x = mut, y = count_norm_adj)) +
  geom_boxplot(outlier.shape = NA) +
  geom_line(aes(group = clone, col = clone), alpha = 0.6) +
  geom_point(aes(col = clone), size = 2, alpha = 0.8) +
  facet_grid(rows = vars(SYMBOL), cols = vars(treat), scales = "free_y") +
  scale_y_log10() +
  labs(y = "normalized counts + 0.5") +
  theme_bw() +
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))
ggsave(paste0("../../out/plot/boxplot_", gene_bold, "_mut.pdf"), width = 6, height = 4)

# dispersion
pdf("../../out/plot/ddsHTSeq_filter_dispersion.pdf", width = 5, height = 5)
plotDispEsts(ddsHTSeq_filter)
dev.off()
