# AIM ---------------------------------------------------------------------
# same processing as 06_gene_coverage.R, but for several genes at once:
# per base coverage of each sample, exon structure, reads spanning the exon-exon junctions and multimapping fraction.
# Reads with a splice gap (N) do not contribute to the coverage of the intron.
# The outputs of each gene have the same names as in 06_gene_coverage.R (suffix _<gene>);
# in addition, the multimapping and junction tables of all the genes are written together.
#
# usage: edit `genes` below, or pass them on the command line:
#   Rscript scr/07_gene_coverage_multi.R C1QA C1QB C1QC

# libraries ---------------------------------------------------------------
library(tidyverse)
library(DESeq2)
library(GenomicAlignments)
library(Rsamtools)
library(rtracklayer)
library(patchwork)

# parameters --------------------------------------------------------------
genes <- c("C1QA", "C1QB", "C1QC","C3")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) genes <- args
bam_dir <- "/beegfs/scratch/ric.absinta/ric.absinta/pipeline/nfcore_bulkRNAseq_AlettaMG/results/star_salmon"
gtf_file <- "../../data/genes.filtered.gtf"
pad <- 500 # bases added on each side of the gene

# samples (the ones of the analysis, see 01_import_data.R) and size factors of the DE
coldata <- read_csv("../../data/LUT_samples_final.csv", show_col_types = FALSE)
dds <- readRDS("../../out/object/dds_filter_DESeq.rds")
sf <- sizeFactors(dds)
sample_meta <- coldata %>% dplyr::select(sample, mut, treat, clone)

gtf <- import(gtf_file)
gtf_exons <- gtf[gtf$type == "exon"]

# loop over the genes -----------------------------------------------------
# for each gene: coverage, junctions, multimapping, tables and plots (same steps as 06_gene_coverage.R)
# a gene without exons in the gtf is skipped (NULL)
out <- map(setNames(genes, genes), function(gene) {
  message("\n=== ", gene, " ===")

  # gene model
  exons <- gtf_exons[gtf_exons$gene_id == gene]
  if (length(exons) == 0) {
    warning(gene, ": no exons found in the gtf, skipped")
    return(NULL)
  }
  exons_df <- as.data.frame(reduce(exons)) # union of the exons of all the transcripts
  region <- GRanges(unique(as.character(seqnames(exons))), IRanges(min(start(exons)) - pad, max(end(exons)) + pad))
  chr <- as.character(seqnames(region))
  message("region: ", as.character(region), " | strand: ", unique(as.character(strand(exons))))

  # coverage and junctions
  # per base coverage normalised with the DESeq2 size factor, and junction reads (split reads) per sample
  # the NH tag (number of loci of the read) separates the uniquely mapped alignments (NH == 1) from the multimapping ones
  param <- ScanBamParam(which = region, tag = "NH")
  res <- lapply(setNames(coldata$sample, coldata$sample), function(s) {
    aln <- readGAlignments(file.path(bam_dir, paste0(s, ".markdup.sorted.bam")), param = param)
    uniq <- mcols(aln)$NH == 1
    cov <- as.numeric(coverage(aln)[[chr]][start(region):end(region)])
    cov_uniq <- as.numeric(coverage(aln[uniq])[[chr]][start(region):end(region)])
    list(cov = tibble(pos = start(region):end(region), cov = cov,
                      cov_norm = cov / sf[[s]], cov_uniq_norm = cov_uniq / sf[[s]]),
         n_reads = length(aln),
         n_multi = sum(!uniq),
         junc = summarizeJunctions(aln) %>% as.data.frame() %>% as_tibble())
  })

  df_cov <- map_dfr(res, "cov", .id = "sample") %>%
    left_join(coldata %>% dplyr::select(sample, sample_name, clone, mut, treat), by = "sample") %>%
    mutate(sample_name = factor(sample_name, levels = unique(sample_name[order(mut, treat, clone)])))

  # reads in the region and junction reads
  # multimapping alignments (NH > 1) in the region: estimate of the reads that are not specific to the gene.
  # paired-end: the two mates are counted as separate alignments; a multimapper contributes the alignments that fall in the region
  df_multi <- tibble(sample = names(res),
                     n_alignments = map_int(res, "n_reads"),
                     n_multimapping = map_int(res, "n_multi")) %>%
    mutate(frac_multimapping = n_multimapping / n_alignments) %>%
    left_join(sample_meta, by = "sample")
  print(df_multi, n = Inf)
  write_tsv(df_multi, paste0("../../out/table/multimapping_", gene, ".tsv"))

  df_junc <- map_dfr(res, "junc", .id = "sample") %>%
    dplyr::select(sample, start, end, score) %>%
    dplyr::rename(n_split_reads = score) %>%
    left_join(sample_meta, by = "sample")
  write_tsv(df_junc, paste0("../../out/table/junctions_", gene, ".tsv"))

  # plots
  # gene model (union of the exons)
  p_gene <- ggplot(exons_df) +
    geom_segment(aes(x = min(start), xend = max(end), y = 0, yend = 0)) +
    geom_rect(aes(xmin = start, xmax = end, ymin = -0.3, ymax = 0.3), fill = "black") +
    scale_x_continuous(limits = c(start(region), end(region))) +
    labs(x = paste0(gene, " (", chr, ")"), y = NULL) +
    theme_bw() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid = element_blank())

  # coverage of each sample, WT vs KO, size factor normalised
  p_cov <- ggplot(df_cov, aes(pos, cov_norm, fill = mut)) +
    geom_area() +
    facet_grid(sample_name ~ ., scales = "free_y", switch = "y") +
    scale_fill_manual(values = c(WT = "gray50", KO = "red")) +
    labs(x = NULL, y = "normalized coverage", title = gene) +
    theme_bw() +
    theme(strip.background = element_blank(), strip.text.y.left = element_text(angle = 0, size = 6),
          axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.text.x = element_blank(), panel.grid = element_blank())
  ggsave(paste0("../../out/plot/coverage_", gene, ".pdf"),
         p_cov / p_gene + plot_layout(heights = c(nrow(coldata), 1.5)),
         width = 8, height = 0.45 * nrow(coldata) + 2)

  # the same on the log scale, to see the low coverage regions
  p_cov_log <- p_cov + scale_y_continuous(trans = "log1p")
  ggsave(paste0("../../out/plot/coverage_log_", gene, ".pdf"),
         p_cov_log / p_gene + plot_layout(heights = c(nrow(coldata), 1.5)),
         width = 8, height = 0.45 * nrow(coldata) + 2)

  # overlaid traces on a common scale: one line per sample (normalized coverage), coloured by mut; thick lines = mean per mut
  df_mean <- df_cov %>%
    group_by(mut, pos) %>%
    summarise(cov_norm = mean(cov_norm), .groups = "drop")

  p_trace <- ggplot(df_cov, aes(pos, cov_norm, col = mut)) +
    geom_line(aes(group = sample), alpha = 0.4, linewidth = 0.4) +
    geom_line(data = df_mean, linewidth = 1) +
    scale_color_manual(values = c(WT = "gray40", KO = "red")) +
    labs(x = NULL, y = "normalized coverage", title = gene) +
    theme_bw() +
    theme(axis.text.x = element_blank(), panel.grid = element_blank())
  ggsave(paste0("../../out/plot/coverage_overlay_", gene, ".pdf"),
         p_trace / p_gene + plot_layout(heights = c(5, 1)), width = 8, height = 5)

  # overlaid traces: uniquely mapped alignments only (NH == 1), same y limits as the plot with all the reads
  df_mean_uniq <- df_cov %>%
    group_by(mut, pos) %>%
    summarise(cov_norm = mean(cov_uniq_norm), .groups = "drop")

  p_trace_uniq <- ggplot(df_cov, aes(pos, cov_uniq_norm, col = mut)) +
    geom_line(aes(group = sample), alpha = 0.4, linewidth = 0.4) +
    geom_line(data = df_mean_uniq, aes(y = cov_norm), linewidth = 1) +
    scale_color_manual(values = c(WT = "gray40", KO = "red")) +
    labs(x = NULL, y = "normalized coverage (unique reads)", title = gene) +
    coord_cartesian(ylim = c(0, max(df_cov$cov_norm))) +
    theme_bw() +
    theme(axis.text.x = element_blank(), strip.background = element_blank())
  ggsave(paste0("../../out/plot/coverage_overlay_unique_", gene, ".pdf"),
         p_trace_uniq / p_gene + plot_layout(heights = c(5, 1)), width = 8, height = 5)

  # fraction of multimapping alignments in the locus per sample
  p_multi <- ggplot(df_multi, aes(x = reorder(sample, frac_multimapping), y = frac_multimapping, fill = mut)) +
    geom_col() +
    scale_fill_manual(values = c(WT = "gray50", KO = "red")) +
    coord_flip() +
    labs(x = NULL, y = paste0("fraction of multimapping alignments in the ", gene, " locus")) +
    theme_bw()
  ggsave(paste0("../../out/plot/multimapping_fraction_", gene, ".pdf"), p_multi, width = 5, height = 4)

  list(multi = df_multi, junc = df_junc)
})

# tables with all the genes together --------------------------------------
skipped <- names(out)[map_lgl(out, is.null)]
if (length(skipped) > 0) message("\ngenes not processed: ", paste(skipped, collapse = ", "))

out <- purrr::compact(out)
if (length(out) > 0) {
  tag <- paste(names(out), collapse = "_")
  if (nchar(tag) > 60) tag <- paste0(length(out), "genes")
  write_tsv(map_dfr(out, "multi", .id = "gene"), paste0("../../out/table/multimapping_all_", tag, ".tsv"))
  write_tsv(map_dfr(out, "junc", .id = "gene"), paste0("../../out/table/junctions_all_", tag, ".tsv"))
}
