# AIM ---------------------------------------------------------------------
# run EnrichR on the DEGs of each contrast (03_differential_analysis.R), split by direction.
# SHRUNK version: the DEGs are defined on the shrunk log2FC (res_<contrast>_shr.txt); the version on the unshrunk results is 04_enrichr.R.
# Outputs have the suffix _shr.
# (needs internet access: the queries are sent to the Enrichr web service)

# libraries ---------------------------------------------------------------
library(tidyverse)
library(enrichR)
library(scales)
library(patchwork)

# parameters --------------------------------------------------------------
contrasts <- c("MyelinvsBL", "KOvsWT")
padj_cut <- 0.05
lfc_cut <- 0.75 # shrunk log2FC, same threshold as the volcano / DEG heatmap of the shrunk version in 03_differential_analysis.R
n_top <- 10 # terms shown per database

# DB selection ------------------------------------------------------------
# dbs <- listEnrichrDbs()
# dbs %>% filter(str_detect(libraryName, pattern = "Atlas"))
dbs_db <- c("KEGG_2021_Human", "MSigDB_Hallmark_2020", "Reactome_Pathways_2024")

# split by direction ------------------------------------------------------
# shrunk tables written by 03_differential_analysis.R (columns: GENEID, baseMean, log2FoldChange, lfcSE, pvalue, padj, SYMBOL)
df_res <- map_dfr(setNames(contrasts, contrasts),
                  ~ read_tsv(paste0("../../out/table/res_", .x, "_shr.txt"), show_col_types = FALSE),
                  .id = "comparison") %>%
  filter(!is.na(SYMBOL), !is.na(padj), padj < padj_cut, abs(log2FoldChange) > lfc_cut) %>%
  mutate(direction = if_else(log2FoldChange > 0, "UP", "DOWN"))

df_res %>% dplyr::count(comparison, direction)

# genes per comparison, one list per direction
get_genes <- function(dir) {
  df_res %>%
    filter(direction == dir) %>%
    split(f = .$comparison) %>%
    map(~ unique(.x$SYMBOL))
}
list_genes <- list(UP = get_genes("UP"), DOWN = get_genes("DOWN"))

# define the background (genes tested in the DE), uncomment to use it in enrichr()
# background <- read_tsv("../../out/table/res_MyelinvsBL.txt") %>% filter(!is.na(padj)) %>% pull(SYMBOL) %>% na.omit()

# query -------------------------------------------------------------------
run_enrichr <- function(genes) {
  out_enrich <- enrichr(genes, dbs_db)
  # drop the databases without an output
  out_enrich[map_int(out_enrich, nrow) > 0] %>%
    bind_rows(.id = "annotation")
}

list_enrichr <- map(list_genes, function(lg) {
  lg %>%
    keep(~ length(.x) >= 3) %>% # too few genes: no meaningful enrichment
    map(run_enrichr) %>%
    bind_rows(.id = "comparison")
})

iwalk(list_enrichr, function(x, dir) {
  write_tsv(x, paste0("../../out/table/enrichR_DE_", dir, "_shr.tsv"))
})

# plots -------------------------------------------------------------------
# one pdf per direction, one panel column per comparison, one row per database
iwalk(list_enrichr, function(x, dir) {
  if (nrow(x) == 0) return(invisible(NULL))
  plot_list <- x %>%
    split(f = .$comparison) %>%
    imap(function(df, comparison) {
      df %>%
        group_by(annotation) %>%
        arrange(P.value) %>%
        dplyr::slice(1:n_top) %>%
        mutate(Term = str_sub(Term, start = 1, end = 30)) %>%
        mutate(Term = fct_reorder(Term, Combined.Score, .desc = FALSE)) %>%
        ggplot(aes(y = Term, x = Combined.Score, size = Odds.Ratio, col = Adjusted.P.value)) +
        geom_point() +
        facet_wrap(~ annotation, scales = "free", ncol = 1) +
        theme_bw() +
        scale_color_gradientn(colors = c("red", "blue"),
                              values = rescale(c(0, 1)),
                              limits = c(0, 0.2)) +
        theme(strip.background = element_blank(),
              panel.border = element_rect(colour = "black", fill = NA)) +
        ggtitle(paste(comparison, dir))
    })
  p <- wrap_plots(plot_list, nrow = 1)
  ggsave(paste0("../../out/plot/enrichR_DE_", dir, "_shr.pdf"), plot = p,
         width = 9 * length(plot_list), height = 12, limitsize = FALSE)
})
