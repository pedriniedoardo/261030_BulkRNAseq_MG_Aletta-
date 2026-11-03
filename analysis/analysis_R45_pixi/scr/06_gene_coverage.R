# AIM ---------------------------------------------------------------------
# inspect the position of the reads along a gene (default C1QA) in the STAR BAM files of the nf-core/rnaseq run:
# per base coverage of each sample, exon structure, and reads spanning the exon-exon junctions.
# Reads with a splice gap (N) do not contribute to the coverage of the intron.

# libraries ---------------------------------------------------------------
library(tidyverse)
library(DESeq2)
library(GenomicAlignments)
library(Rsamtools)
library(rtracklayer)
library(patchwork)

# parameters --------------------------------------------------------------
gene <- "C1QA"
bam_dir <- "/beegfs/scratch/ric.absinta/ric.absinta/pipeline/nfcore_bulkRNAseq_AlettaMG/results/star_salmon"
gtf_file <- "../../data/genes.filtered.gtf"
pad <- 500 # bases added on each side of the gene

# samples (the ones of the analysis, see 01_import_data.R) and size factors of the DE
coldata <- read_csv("../../data/LUT_samples_final.csv", show_col_types = FALSE)
dds <- readRDS("../../out/object/dds_filter_DESeq.rds")
sf <- sizeFactors(dds)

# gene model --------------------------------------------------------------
gtf <- import(gtf_file)
exons <- gtf[gtf$type == "exon" & gtf$gene_id == gene]
exons_df <- as.data.frame(reduce(exons)) # union of the exons of all the transcripts
region <- GRanges(unique(as.character(seqnames(exons))), IRanges(min(start(exons)) - pad, max(end(exons)) + pad))
message("region: ", as.character(region), " | strand: ", unique(as.character(strand(exons))))

# coverage and junctions --------------------------------------------------
# per base coverage normalised with the DESeq2 size factor, and junction reads (split reads) per sample
# the NH tag (number of loci of the read) is read to separate the uniquely mapped alignments (NH == 1) from the multimapping ones
param <- ScanBamParam(which = region, tag = "NH")
chr <- as.character(seqnames(region))
cov_region <- function(a) as.numeric(coverage(a)[[chr]][start(region):end(region)])
res <- lapply(setNames(coldata$sample, coldata$sample), function(s) {
  bam <- file.path(bam_dir, paste0(s, ".markdup.sorted.bam"))
  aln <- readGAlignments(bam, param = param)
  uniq <- mcols(aln)$NH == 1
  cov <- cov_region(aln)
  cov_uniq <- cov_region(aln[uniq])
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
  left_join(coldata %>% dplyr::select(sample, mut, treat, clone), by = "sample")
print(df_multi, n = Inf)
write_tsv(df_multi, paste0("../../out/table/multimapping_", gene, ".tsv"))
df_junc <- map_dfr(res, "junc", .id = "sample") %>%
  dplyr::select(sample, start, end, score) %>%
  dplyr::rename(n_split_reads = score) %>%
  left_join(coldata %>% dplyr::select(sample, mut, treat, clone), by = "sample")
write_tsv(df_junc, paste0("../../out/table/junctions_", gene, ".tsv"))

# plots -------------------------------------------------------------------
# coverage of each sample, WT vs KO, size factor normalised
p_cov <- ggplot(df_cov, aes(pos, cov_norm, fill = mut)) +
  geom_area() +
  facet_grid(sample_name ~ ., scales = "free_y", switch = "y") +
  scale_fill_manual(values = c(WT = "gray50", KO = "red")) +
  labs(x = NULL, y = "normalized coverage") +
  theme_bw() +
  theme(strip.background = element_blank(), strip.text.y.left = element_text(angle = 0, size = 6),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_blank(), panel.grid = element_blank())

# gene model (union of the exons)
p_gene <- ggplot(exons_df) +
  geom_segment(aes(x = min(start), xend = max(end), y = 0, yend = 0)) +
  geom_rect(aes(xmin = start, xmax = end, ymin = -0.3, ymax = 0.3), fill = "black") +
  scale_x_continuous(limits = c(start(region), end(region))) +
  labs(x = paste0(gene, " (", as.character(seqnames(region)), ")"), y = NULL) +
  theme_bw() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid = element_blank())

p_cov / p_gene + plot_layout(heights = c(nrow(coldata), 1.5))
ggsave(paste0("../../out/plot/coverage_", gene, ".pdf"), width = 8, height = 0.45 * nrow(coldata) + 2)

# the same on the log scale, to see the low coverage regions
p_cov_log <- p_cov + scale_y_continuous(trans = "log1p")
p_cov_log / p_gene + plot_layout(heights = c(nrow(coldata), 1.5))
ggsave(paste0("../../out/plot/coverage_log_", gene, ".pdf"), width = 8, height = 0.45 * nrow(coldata) + 2)

# overlaid traces on a common scale ---------------------------------------
# one line per sample on the same axis (normalized coverage), coloured by mut: easier to compare the level of WT and KO.
# thick lines = mean per mut
df_mean <- df_cov %>%
  group_by(mut, pos) %>%
  summarise(cov_norm = mean(cov_norm), .groups = "drop")

p_trace <- ggplot(df_cov, aes(pos, cov_norm, col = mut)) +
  geom_line(aes(group = sample), alpha = 0.4, linewidth = 0.4) +
  geom_line(data = df_mean, linewidth = 1) +
  scale_color_manual(values = c(WT = "gray40", KO = "red")) +
  labs(x = NULL, y = "normalized coverage") +
  theme_bw() +
  theme(axis.text.x = element_blank(), panel.grid = element_blank())

p_trace / p_gene + plot_layout(heights = c(5, 1))
ggsave(paste0("../../out/plot/coverage_overlay_", gene, ".pdf"), width = 8, height = 5)

# overlaid traces: uniquely mapped reads only ----------------------------
# same plot as the overlay above (all the alignments), but using only the uniquely mapped alignments (NH == 1) in the locus
df_mean_uniq <- df_cov %>%
  group_by(mut, pos) %>%
  summarise(cov_norm = mean(cov_uniq_norm), .groups = "drop")

p_trace_uniq <- ggplot(df_cov, aes(pos, cov_uniq_norm, col = mut)) +
  geom_line(aes(group = sample), alpha = 0.4, linewidth = 0.4) +
  geom_line(data = df_mean_uniq, aes(y = cov_norm), linewidth = 1) +
  scale_color_manual(values = c(WT = "gray40", KO = "red")) +
  labs(x = NULL, y = "normalized coverage (unique reads)") +
  # same y limits as the plot with all the reads, to compare the two directly
  coord_cartesian(ylim = c(0, max(df_cov$cov_norm))) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        # panel.grid = element_blank(),
        strip.background = element_blank())

p_trace_uniq / p_gene + plot_layout(heights = c(5, 1))
ggsave(paste0("../../out/plot/coverage_overlay_unique_", gene, ".pdf"), width = 8, height = 5)

# fraction of multimapping alignments in the locus per sample
ggplot(df_multi, aes(x = reorder(sample, frac_multimapping), y = frac_multimapping, fill = mut)) +
  geom_col() +
  scale_fill_manual(values = c(WT = "gray50", KO = "red")) +
  coord_flip() +
  labs(x = NULL, y = paste0("fraction of multimapping alignments in the ", gene, " locus")) +
  theme_bw()
ggsave(paste0("../../out/plot/multimapping_fraction_", gene, ".pdf"), width = 5, height = 4)
