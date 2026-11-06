# AIM ---------------------------------------------------------------------
# coverage of a gene with many exons (default C3), restricted to the exons and binned to smooth the trend.
# The exons (union of the exons of all the transcripts) are put side by side in transcript direction (5' -> 3'),
# the introns are removed. The coverage is averaged in bins of `bin` bases inside each exon (bins never cross an exon boundary).
# Same normalisation as 06_gene_coverage.R (DESeq2 size factors; all alignments and uniquely mapped alignments NH == 1).
# Exon numbers are those of the union of the exons, not the numbering of any specific transcript.
#
# usage: edit `gene` and `bin` below, or from the command line:
#   Rscript scr/08_gene_coverage_exons.R C3 25

# libraries ---------------------------------------------------------------
library(tidyverse)
library(DESeq2)
library(GenomicAlignments)
library(Rsamtools)
library(rtracklayer)
library(patchwork)

# parameters --------------------------------------------------------------
gene <- "C3"
bin <- 50 # bases per bin (1 = no binning)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) gene <- args[1]
if (length(args) > 1) bin <- as.integer(args[2])
bam_dir <- "/beegfs/scratch/ric.absinta/ric.absinta/pipeline/nfcore_bulkRNAseq_AlettaMG/results/star_salmon"
gtf_file <- "../../data/genes.filtered.gtf"

# samples and size factors of the DE
coldata <- read_csv("../../data/LUT_samples_final.csv", show_col_types = FALSE)
dds <- readRDS("../../out/object/dds_filter_DESeq.rds")
sf <- sizeFactors(dds)

# exon model --------------------------------------------------------------
gtf <- import(gtf_file)
exons <- gtf[gtf$type == "exon" & gtf$gene_id == gene]
stopifnot("no exons found for this gene" = length(exons) > 0)
strand_gene <- unique(as.character(strand(exons)))
stopifnot("gene on more than one strand" = length(strand_gene) == 1)
chr <- unique(as.character(seqnames(exons)))
region <- GRanges(chr, IRanges(min(start(exons)), max(end(exons))))
message(gene, " | region: ", as.character(region), " | strand: ", strand_gene)

# exons in transcript order (5' -> 3'), with their offset in the concatenated exon space
exons_df <- as.data.frame(reduce(exons)) %>%
  # dplyr::select(start, end, width) %>%
  dplyr::arrange(if (strand_gene == "-") dplyr::desc(start) else start) %>%
  mutate(exon_n = row_number(), offset = cumsum(lag(width, default = 0)))
total_width <- sum(exons_df$width)
message(nrow(exons_df), " exons, ", total_width, " bases")

# one row per exonic base: genomic position, x in the concatenated space and bin
base_df <- exons_df %>%
  mutate(pos = map2(start, end, seq)) %>%
  unnest(pos) %>%
  mutate(rel = if (strand_gene == "-") end - pos else pos - start,
         x = offset + rel,
         bin_id = rel %/% bin)

# coverage ----------------------------------------------------------------
# per base coverage of the exonic bases of each sample, then mean in each bin
param <- ScanBamParam(which = region, tag = "NH")
df_base <- map_dfr(setNames(coldata$sample, coldata$sample), function(s) {
  aln <- readGAlignments(file.path(bam_dir, paste0(s, ".markdup.sorted.bam")), param = param)
  uniq <- mcols(aln)$NH == 1
  cov <- as.numeric(coverage(aln)[[chr]][start(region):end(region)])
  cov_uniq <- as.numeric(coverage(aln[uniq])[[chr]][start(region):end(region)])
  i <- base_df$pos - start(region) + 1
  base_df %>% mutate(cov_norm = cov[i] / sf[[s]], cov_uniq_norm = cov_uniq[i] / sf[[s]])
}, .id = "sample")

df_bin <- df_base %>%
  group_by(sample, exon_n, bin_id) %>%
  summarise(x = mean(x), cov_norm = mean(cov_norm), cov_uniq_norm = mean(cov_uniq_norm), .groups = "drop") %>%
  left_join(coldata %>% dplyr::select(sample, sample_name, clone, mut, treat), by = "sample") %>%
  mutate(sample_name = factor(sample_name, levels = unique(sample_name[order(mut, treat, clone)])))
write_tsv(df_bin, paste0("../../out/table/coverage_exons_binned_", gene, ".tsv"))

# mean per mut, per bin
df_mean <- df_bin %>%
  group_by(mut, exon_n, bin_id, x) %>%
  summarise(cov_norm = mean(cov_norm), cov_uniq_norm = mean(cov_uniq_norm), .groups = "drop")

# plots -------------------------------------------------------------------
# exon track: one box per exon, numbered in transcript direction
x_scale <- scale_x_continuous(limits = c(0, total_width), expand = c(0, 0))
p_gene <- ggplot(exons_df) +
  geom_rect(aes(xmin = offset, xmax = offset + width, ymin = -0.3, ymax = 0.3), fill = "black", col = "white") +
  geom_text(aes(x = offset + width / 2, y = 0.6, label = exon_n), size = 2) +
  x_scale +
  labs(x = paste0(gene, " exons, 5' -> 3' (", chr, ", ", strand_gene, " strand; bin = ", bin, " bp)"), y = NULL) +
  theme_bw() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid = element_blank())

# exon boundaries, shown in the coverage plots
bounds <- geom_vline(data = exons_df %>% filter(offset > 0), aes(xintercept = offset), col = "gray80", linewidth = 0.2)

# coverage of each sample, WT vs KO, size factor normalised
p_cov <- ggplot(df_bin, aes(x, cov_norm, fill = mut, group = exon_n)) +
  bounds +
  geom_area() +
  facet_grid(sample_name ~ ., scales = "free_y", switch = "y") +
  scale_fill_manual(values = c(WT = "gray50", KO = "red")) +
  x_scale +
  labs(x = NULL, y = "normalized coverage", title = gene) +
  theme_bw() +
  theme(strip.background = element_blank(), strip.text.y.left = element_text(angle = 0, size = 6),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_blank(), panel.grid = element_blank())
ggsave(paste0("../../out/plot/coverage_exons_", gene, ".pdf"),
       p_cov / p_gene + plot_layout(heights = c(nrow(coldata), 1.5)),
       width = 10, height = 0.45 * nrow(coldata) + 2)

# overlaid traces, all alignments: one line per sample, thick line = mean per mut
p_trace <- ggplot(df_bin, aes(x, cov_norm, col = mut)) +
  bounds +
  geom_line(aes(group = interaction(sample, exon_n)), alpha = 0.4, linewidth = 0.4) +
  geom_line(data = df_mean, aes(group = interaction(mut, exon_n)), linewidth = 1) +
  scale_color_manual(values = c(WT = "gray40", KO = "red")) +
  x_scale +
  labs(x = NULL, y = "normalized coverage", title = gene) +
  theme_bw() +
  theme(axis.text.x = element_blank(), panel.grid = element_blank())
ggsave(paste0("../../out/plot/coverage_exons_overlay_", gene, ".pdf"),
       p_trace / p_gene + plot_layout(heights = c(5, 1)), width = 10, height = 5)

# overlaid traces, uniquely mapped alignments only (NH == 1), same y limits as the plot with all the alignments
p_trace_uniq <- ggplot(df_bin, aes(x, cov_uniq_norm, col = mut)) +
  bounds +
  geom_line(aes(group = interaction(sample, exon_n)), alpha = 0.4, linewidth = 0.4) +
  geom_line(data = df_mean, aes(y = cov_uniq_norm, group = interaction(mut, exon_n)), linewidth = 1) +
  scale_color_manual(values = c(WT = "gray40", KO = "red")) +
  x_scale +
  coord_cartesian(ylim = c(0, max(df_bin$cov_norm))) +
  labs(x = NULL, y = "normalized coverage (unique reads)", title = gene) +
  theme_bw() +
  theme(axis.text.x = element_blank(), panel.grid = element_blank())
ggsave(paste0("../../out/plot/coverage_exons_overlay_unique_", gene, ".pdf"),
       p_trace_uniq / p_gene + plot_layout(heights = c(5, 1)), width = 10, height = 5)
