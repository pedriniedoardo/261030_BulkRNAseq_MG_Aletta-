# AIM ---------------------------------------------------------------------
# run GSEA (fgsea) on the two contrasts of 03_differential_analysis.R (MyelinvsBL, KOvsWT).
# SHRUNK version: all the tables read are res_<contrast>_shr.txt (shrunk log2FC); the version on the unshrunk results is 05_GSEA.R.
# Outputs have the suffix _shr.
# with two alternative rankings of all the tested genes:
#   - logfc:       shrunk log2FC (res_<contrast>_shr.txt)
#   - signed_logp: sign(log2FC) * -log10(pvalue) (res_<contrast>_shr.txt; same p-values as the unshrunk table), see
#     https://github.com/crazyhottommy/compbio_tutorials/blob/main/scripts/07_gene_set_enrichment_RNAseq.Rmd

# libraries ---------------------------------------------------------------
library(tidyverse)
library(fgsea)
library(msigdbr)
library(ggrepel)
library(patchwork)
library(DESeq2)
library(ComplexHeatmap)
library(pals)

# parameters --------------------------------------------------------------
contrasts <- c("MyelinvsBL", "KOvsWT")
n_label <- 10 # pathways labelled in each panel of the GSEA plots (lowest padj)
species <- "Homo sapiens" # gene symbols in the GTF are human

# gene set collections: name = label used in the output files
# (use msigdbr_collections() to list the available ones)
collections <- list(Hallmark = list(collection = "H", subcollection = NULL),
                    KEGG_MEDICUS = list(collection = "C2", subcollection = "CP:KEGG_MEDICUS"),
                    Reactome = list(collection = "C2", subcollection = "CP:REACTOME"))

# ranking of the genes ----------------------------------------------------
# shrunk log2FC; mean in case of duplicated symbols
list_ranks_logfc <- map(setNames(contrasts, contrasts), function(x) {
  df <- read_tsv(paste0("../../out/table/res_", x, "_shr.txt"), show_col_types = FALSE) %>%
    # genes without symbol or without a reliable estimate of the p value are removed
    dplyr::filter(!is.na(SYMBOL), !is.na(pvalue)) %>%
    group_by(SYMBOL) %>%
    summarise(stat = mean(log2FoldChange), .groups = "drop")
  setNames(df$stat, df$SYMBOL)
})

# signed -log10(pvalue); pvalue = 0 (Inf) set to a large number
list_ranks_signed_logp <- map(setNames(contrasts, contrasts), function(x) {
  df <- read_tsv(paste0("../../out/table/res_", x, "_shr.txt"), show_col_types = FALSE) %>%
    dplyr::filter(!is.na(SYMBOL), !is.na(pvalue)) %>%
    group_by(SYMBOL) %>%
    summarise(pvalue = mean(pvalue), log2FC = mean(log2FoldChange), .groups = "drop") %>%
    mutate(neg_log10_p = -log10(pvalue),
           neg_log10_p = ifelse(is.infinite(neg_log10_p), 1000, neg_log10_p),
           stat = sign(log2FC) * neg_log10_p)
  setNames(df$stat, df$SYMBOL)
})

list_ranks_all <- list(logfc = list_ranks_logfc, signed_logp = list_ranks_signed_logp)
glimpse(list_ranks_all)

# gene sets ---------------------------------------------------------------
list_pathways <- imap(collections, function(p, nm) {
  gs <- msigdbr(species = species, collection = p$collection, subcollection = p$subcollection)
  split(x = gs$gene_symbol, f = gs$gs_name)
})

# run GSEA ----------------------------------------------------------------
# one row per ranking x collection x contrast; fgsea on each, and the redundant pathways of each run are collapsed
# (collapsePathways only needs the pathway / pval / padj columns of the fgsea output)
df_crossing <- crossing(rank_type = names(list_ranks_all),
                       collection = names(list_pathways),
                       dataset = contrasts)

# rank_type <- "logfc"
# collection <- "Hallmark"
# dataset <- "KOvsWT"
            
df_gsea <- pmap(df_crossing, function(rank_type,collection,dataset) {
  ranks <- list_ranks_all[[rank_type]][[dataset]]
  pathways <- list_pathways[[collection]]
  
  message("GSEA: ", collection, " | ", rank_type, " | ", dataset)
  # seed set in each run (not once at the top) so that the result does not depend on the runs before it; serial to make the seed effective
  set.seed(1)
  res <- fgsea(pathways, ranks, minSize = 10, maxSize = 500, BPPARAM = BiocParallel::SerialParam())
  
  res_df <- res %>%
    as_tibble() %>%
    mutate(leadingEdge = map_chr(leadingEdge, paste0, collapse = "|"))
  
  # TRUE for the pathways kept after collapsing the redundant ones
  res_full <- res_df %>%
    mutate(main_pathway = pathway %in% collapsePathways(res, pathways, ranks)$mainPathways) %>%
    mutate(rank_type_id = rank_type,
           collection_id = collection,
           dataset_id = dataset)
  
  return(res_full)
}) %>%
  bind_rows() %>%
  arrange(rank_type_id, collection_id, dataset_id, padj, -abs(NES))

# tables and plots --------------------------------------------------------
# one table + one plot (contrasts as panels) per ranking x collection, for all the pathways and for the non redundant ones
df_gsea %>%
  mutate(tag = paste(collection_id, rank_type_id, "shr", sep = "_")) %>%
  split(f = .$tag) %>%
  iwalk(function(df, tag) {
    list(all = df, nonredundant = dplyr::filter(df, main_pathway)) %>%
      iwalk(function(d, version) {
        suffix <- if (version == "all") "" else "_nonredundant"
        write_tsv(d, paste0("../../out/table/GSEA_", tag, suffix, ".tsv"))

        p <- d %>%
          # short labels (prefix removed, underscores as spaces, truncated), only for the top n_label pathways of each contrast
          mutate(pathway2 = str_remove(pathway, "^(KEGG_MEDICUS_|HALLMARK_|REACTOME_)") %>%
                   str_replace_all("_", " ") %>%
                   str_to_lower() %>%
                   str_trunc(25)) %>%
          group_by(dataset_id) %>%
          mutate(label = if_else(min_rank(padj) <= n_label, pathway2, NA_character_)) %>%
          ungroup() %>%
          ggplot(aes(y = -log10(padj), x = NES, label = label)) +
          geom_point(aes(size = size, col = padj < 0.05), alpha = 0.4) +
          scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
          facet_wrap(~ dataset_id) + theme_bw() +
          theme(strip.background = element_blank()) +
          geom_text_repel(size = 2.5, bg.color = "white", bg.r = 0.15, box.padding = 0.6, point.padding = 0.3, force = 3, segment.alpha = 0.4, max.overlaps = Inf) +
          geom_hline(yintercept = -log10(0.05), col = "gray", linetype = "dashed") +
          ggtitle(paste(tag, if (version == "all") "" else "(non redundant)"))
        ggsave(paste0("../../out/plot/GSEA_", tag, suffix, ".pdf"), plot = p, width = 12, height = 6)
      })
  })

# heatmap of the leading edge of a significant term ----------------------
# vst expression (z-score per gene) of the leading edge genes of each term, one heatmap per term x contrast where it is significant;
# columns annotated with clone, mut and treat (same palettes as the sex genes heatmap in 02_QC.R)
terms <- c("HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_COMPLEMENT")
gsea_rank <- "signed_logp" # ranking used (the term is not significant with the logfc ranking)
gsea_collection <- "Hallmark"
# DEG thresholds for the leading edge genes shown (same as the shrunk version of 03_differential_analysis.R: padj and shrunk log2FC)
padj_de <- 0.05
lfc_de <- 0.75 # shrunk log2FC

# vst of the filtered dataset (02_QC.R); the rownames are the gene symbols (gene_id of the GTF)
vds_filter <- readRDS("../../out/object/vds_filter.rds")
cd <- as.data.frame(colData(vds_filter))
col_ha <- list(clone = setNames(pals::glasbey(nlevels(cd$clone)), levels(cd$clone)),
               mut = setNames(pals::okabe(nlevels(cd$mut)), levels(cd$mut)),
               treat = setNames(pals::polychrome(nlevels(cd$treat)), levels(cd$treat)))

df_gsea %>%
  dplyr::filter(rank_type_id == gsea_rank, collection_id == gsea_collection, pathway %in% terms, padj < 0.05) %>%
  split(f = list(.$pathway, .$dataset_id), drop = TRUE) %>%
  walk(function(x) {
    term <- x$pathway
    contrast <- x$dataset_id
    # only the leading edge genes that are DEGs of the contrast (03_differential_analysis.R: shrunk log2FC)
    genes_sig <- read_tsv(paste0("../../out/table/res_", contrast, "_shr.txt"), show_col_types = FALSE) %>%
      dplyr::filter(!is.na(padj), padj < padj_de, abs(log2FoldChange) > lfc_de, !is.na(SYMBOL)) %>%
      pull(SYMBOL)
    # the genes removed by the expression filter are not in the vst
    genes_le <- str_split(x$leadingEdge, "\\|")[[1]]
    genes <- intersect(genes_le, genes_sig) %>% intersect(rownames(vds_filter))
    message(term, " | ", contrast, ": ", length(genes_le), " leading edge genes, ", length(genes), " significant in the DE")
    # nothing to show (at least 2 genes are needed to cluster the rows)
    if (length(genes) < 2) {
      message("  -> heatmap skipped")
      return(invisible(NULL))
    }

    mat_z <- assay(vds_filter)[genes, , drop = FALSE] %>%
      t() %>% scale() %>% t()
    mat_z[is.na(mat_z)] <- 0 # genes with no variation
    colnames(mat_z) <- cd[colnames(mat_z), "sample_name"]

    pdf(paste0("../../out/plot/heatmap_leadingEdge_", term, "_", contrast, "_shr.pdf"),
        width = 6, height = max(6, 0.3 * length(genes) + 3))
    draw(Heatmap(mat_z, name = "z-score vst",
                 col = circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
                 top_annotation = HeatmapAnnotation(clone = cd$clone, mut = cd$mut, treat = cd$treat, col = col_ha),
                 row_names_gp = grid::gpar(fontsize = 7),
                 column_title = sprintf("%s\nleading edge (DE padj < %s, |log2FC| > %s)\n%s (NES %.2f, padj %.2g)", term, padj_de, lfc_de, contrast, x$NES, x$padj)),
         heatmap_legend_side = "left", annotation_legend_side = "left")
    dev.off()
  })

# enrichment profile of the top pathways of one contrast / ranking --------
# example: top 9 non redundant Hallmark pathways of MyelinvsBL ranked by the logFC
# df_top <- read_tsv("../../out/table/GSEA_Hallmark_logfc_nonredundant.tsv") %>%
#   dplyr::filter(dataset_id == "MyelinvsBL") %>% dplyr::slice(1:9)
# map(df_top$pathway, ~ plotEnrichment(list_pathways$Hallmark[[.x]], list_ranks_logfc$MyelinvsBL) + labs(title = .x)) %>%
#   wrap_plots()
# ggsave("../../out/plot/GSEA_profile_Hallmark_logfc_MyelinvsBL.pdf", width = 15, height = 10)

# search for a significant term where C1QA is in the leading edges
df_gsea %>%
  filter(str_detect(leadingEdge,pattern = "C1Q")) %>%
  filter(padj<0.05)
