# AIM ---------------------------------------------------------------------
# same plot as coverage_overlay_<gene>.pdf of 06_gene_coverage.R (normalized coverage of each sample along the gene, WT vs KO),
# with the location of the guide RNA of 09_guide_region.R.
# The guide (20 bp) is very short compared with the gene: it is marked with a vertical line at the expected Cas9 cut site,
# and the guide and the PAM are drawn on the gene model.

# libraries ---------------------------------------------------------------
library(tidyverse)
library(DESeq2)
library(GenomicAlignments)
library(Rsamtools)
library(Biostrings)
library(rtracklayer)
library(patchwork)

# parameters --------------------------------------------------------------
# the sequence of the guide as given by the company and the gene that contains it (same as 09_guide_region.R)
guide <- "CTGGCATCCGGACAGGCATC"
gene <- "C1QA"
# optional command line arguments overwrite the two values above
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) guide <- toupper(args[1])
if (length(args) > 1) gene <- args[2]
spacer_strand <- NA # "+" or "-" to force the strand of the protospacer; NA = the strand with an NGG PAM
pad <- 500 # bases added on each side of the gene
bam_dir <- "/beegfs/scratch/ric.absinta/ric.absinta/pipeline/nfcore_bulkRNAseq_AlettaMG/results/star_salmon"
gtf_file <- "../../data/genes.filtered.gtf"
ref_dir <- "/beegfs/scratch/ric.absinta/ric.absinta/reference/nfcore/GRCh38"
fa <- FaFile(file.path(ref_dir, "rsem/genome.fa"), index = file.path(ref_dir, "genome.fa.fai"))

# samples (the ones of the analysis, see 01_import_data.R) and size factors of the DE
coldata <- read_csv("../../data/LUT_samples_final.csv", show_col_types = FALSE)
dds <- readRDS("../../out/object/dds_filter_DESeq.rds")
sf <- sizeFactors(dds)

# gene model --------------------------------------------------------------
gtf <- import(gtf_file)
exons <- gtf[gtf$type == "exon" & gtf$gene_id == gene]
stopifnot("no exons found for this gene" = length(exons) > 0)
exons_df <- as.data.frame(reduce(exons)) # union of the exons of all the transcripts
chr <- unique(as.character(seqnames(exons)))
gene_gr <- GRanges(chr, IRanges(min(start(exons)), max(end(exons))))
region <- GRanges(chr, IRanges(start(gene_gr) - pad, end(gene_gr) + pad))
message("region: ", as.character(region), " | strand: ", unique(as.character(strand(exons))))

# locate the guide, the PAM and the cut site (see 09_guide_region.R) ------
# the sequence is searched in the gene, as it is and as reverse complement
gene_seq <- getSeq(fa, gene_gr)[[1]]
hits <- c(start(matchPattern(DNAString(guide), gene_seq)), start(matchPattern(reverseComplement(DNAString(guide)), gene_seq)))
stopifnot("the sequence has to match exactly once in the gene" = length(hits) == 1)
g_start <- start(gene_gr) + hits - 1
g_end <- g_start + nchar(guide) - 1

# the strand of the protospacer is the one followed by an NGG PAM
pam_plus <- as.character(getSeq(fa, GRanges(chr, IRanges(g_end + 1, g_end + 3))))
pam_minus <- as.character(reverseComplement(getSeq(fa, GRanges(chr, IRanges(g_start - 3, g_start - 1)))))
if (is.na(spacer_strand)) {
  ok <- c("+" = grepl("^.GG$", pam_plus), "-" = grepl("^.GG$", pam_minus))
  stopifnot("no unique strand with an NGG PAM: set spacer_strand by hand" = sum(ok) == 1)
  spacer_strand <- names(ok)[ok]
}
flip <- spacer_strand == "-"
pam <- if (!flip) pam_plus else pam_minus
pam_start <- if (!flip) g_end + 1 else g_start - 3
pam_end <- pam_start + 2
# the cut is between cut_left and cut_left + 1 (3 bases from the PAM)
cut_left <- if (!flip) g_end - 3 else g_start + 2
cut_pos <- cut_left + 0.5
guide_label <- paste0("gRNA ", guide, " | protospacer on the ", spacer_strand, " strand, PAM ", pam,
                      "\n", chr, ":", g_start, "-", g_end, " | expected cut ", cut_left, "^", cut_left + 1)
message(guide_label)

# start and stop codon of the protein (reference) -------------------------
# they are derived from the CDS of each transcript of the annotation and checked on the genome sequence.
# The CDS bases are ordered 5' -> 3' of the transcript (also across introns, if a codon is split by an intron)
cds <- gtf[gtf$type == "CDS" & gtf$gene_id == gene]
stopifnot("no CDS found for this gene (non coding gene?)" = length(cds) > 0)
cds_strand <- unique(as.character(strand(cds)))
df_codon <- map_dfr(split(cds, as.character(cds$transcript_id)), function(tx) {
  tx <- sort(tx, ignore.strand = TRUE)
  p <- unlist(map2(start(tx), end(tx), seq)) # genomic positions of the coding bases
  if (cds_strand == "-") p <- rev(p)
  seq_at <- function(i) paste(as.character(getSeq(fa, GRanges(chr, IRanges(i, width = 1), strand = cds_strand))), collapse = "")
  start_p <- head(p, 3)
  stop_p <- tail(p, 3)
  # in some annotations the CDS does not include the stop codon: then it is the 3 bases after the CDS
  # (assumed not to be split by an intron; the sequence is checked below)
  if (!seq_at(stop_p) %in% c("TAA", "TAG", "TGA")) stop_p <- tail(p, 1) + (1:3) * if (cds_strand == "-") -1 else 1
  tibble(codon = c("start", "stop"), pos = list(start_p, stop_p), seq = c(seq_at(start_p), seq_at(stop_p)))
}, .id = "transcript_id")
stopifnot("the start codon is not ATG" = all(df_codon$seq[df_codon$codon == "start"] == "ATG"),
          "no stop codon at the end of the CDS" = all(df_codon$seq[df_codon$codon == "stop"] %in% c("TAA", "TAG", "TGA")))
# one row per distinct codon position (transcripts can share them); x = middle of the codon
df_codon <- df_codon %>%
  mutate(first = map_dbl(pos, min), last = map_dbl(pos, max), x = map_dbl(pos, median)) %>%
  distinct(codon, seq, first, last, x)
print(df_codon)
codon_label <- paste(paste0(df_codon$codon, " codon ", df_codon$seq, " ", chr, ":", df_codon$first, "-", df_codon$last), collapse = " | ")
codon_cols <- c(start = "darkgreen", stop = "purple")

# coverage ----------------------------------------------------------------
# per base coverage normalised with the DESeq2 size factor (all the alignments, as in 06_gene_coverage.R)
param <- ScanBamParam(which = region)
df_cov <- map_dfr(setNames(coldata$sample, coldata$sample), function(s) {
  aln <- readGAlignments(file.path(bam_dir, paste0(s, ".markdup.sorted.bam")), param = param)
  cov <- as.numeric(coverage(aln)[[chr]][start(region):end(region)])
  tibble(pos = start(region):end(region), cov_norm = cov / sf[[s]])
}, .id = "sample") %>%
  left_join(coldata %>% dplyr::select(sample, mut, treat, clone), by = "sample")

# thick lines = mean per mut
df_mean <- df_cov %>%
  group_by(mut, pos) %>%
  summarise(cov_norm = mean(cov_norm), .groups = "drop")

# plots -------------------------------------------------------------------
col_guide <- "steelblue"
# the cut line is the same in the two panels, so that it is aligned across them
cut_line <- geom_vline(xintercept = cut_pos, linetype = "dashed", linewidth = 0.4, colour = col_guide)
# same x limits in the two panels
x_scale <- scale_x_continuous(limits = c(start(region), end(region)))
# start (green) and stop (purple) codon: dotted lines in the two panels
codon_lines <- geom_vline(xintercept = df_codon$x, linetype = "dotted", linewidth = 0.5, colour = codon_cols[df_codon$codon])

# overlaid traces on a common scale, one line per sample
p_trace <- ggplot(df_cov, aes(pos, cov_norm, col = mut)) +
  cut_line +
  codon_lines +
  geom_line(aes(group = sample), alpha = 0.4, linewidth = 0.4) +
  geom_line(data = df_mean, linewidth = 1) +
  scale_color_manual(values = c(WT = "gray40", KO = "red")) +
  x_scale +
  labs(x = NULL, y = "normalized coverage", subtitle = paste0(guide_label, "\n", codon_label)) +
  theme_bw() +
  theme(axis.text.x = element_blank(), panel.grid = element_blank(), plot.subtitle = element_text(size = 8))

# gene model (union of the exons) with the guide and the PAM above it
# the boxes are drawn to scale (20 and 3 bases): the triangle points at the cut site, to find them on the gene
p_gene <- ggplot(exons_df) +
  cut_line +
  codon_lines +
  geom_segment(aes(x = min(start), xend = max(end), y = 0, yend = 0)) +
  geom_rect(aes(xmin = start, xmax = end, ymin = -0.3, ymax = 0.3), fill = "black") +
  annotate("rect", xmin = g_start - 0.5, xmax = g_end + 0.5, ymin = 0.4, ymax = 0.7, fill = col_guide) +
  annotate("rect", xmin = pam_start - 0.5, xmax = pam_end + 0.5, ymin = 0.4, ymax = 0.7, fill = "orange") +
  annotate("point", x = cut_pos, y = 0.9, shape = 25, size = 2, colour = col_guide, fill = col_guide) +
  annotate("text", x = cut_pos, y = 0.9, label = "  gRNA", hjust = 0, size = 3, colour = col_guide) +
  # codon labels below the gene model
  annotate("text", x = df_codon$x, y = -0.55, label = df_codon$codon, size = 3, colour = codon_cols[df_codon$codon]) +
  x_scale +
  scale_y_continuous(limits = c(-0.7, 1.1)) +
  labs(x = paste0(gene, " (", chr, ")"), y = NULL) +
  theme_bw() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid = element_blank())

p_trace / p_gene + plot_layout(heights = c(5, 1.3))
ggsave(paste0("../../out/plot/10_coverage_overlay_gRNA_", gene, ".pdf"), width = 8, height = 5.3)

# predicted protein of the reference and of the two main edited alleles ----
# The two most frequent non reference alleles of the KO samples (table written by 09_guide_region.R, run it first with the same
# guide and gene) are put in the reference mRNA, which is translated from the reference start codon to the first stop codon
# (also in the 3' UTR). This is a prediction: it assumes that the edited transcript is spliced as the reference one, uses the
# same start codon and carries no other change. It does not say if the protein is made or is stable.
n_alt <- 2 # number of edited alleles
flank_guide <- 10 # as in 09_guide_region.R: it defines the window of the allele sequences
aa_per_line <- 80

df_allele <- read_tsv(paste0("../../out/table/guide_region_alleles_", gene, ".tsv"), show_col_types = FALSE)
ko_samples <- coldata$sample[coldata$mut == "KO"]
# fraction of the reads of each KO sample (0 if the allele is not found in the sample), ranked by the mean over the samples
df_top <- df_allele %>%
  filter(sample %in% ko_samples) %>%
  dplyr::select(sample, allele, frac) %>%
  complete(sample = ko_samples, allele, fill = list(frac = 0)) %>%
  group_by(allele) %>%
  summarise(mean_frac = mean(frac), reads = paste0(sample, " ", round(100 * frac), "%", collapse = ", "), .groups = "drop") %>%
  filter(!allele %in% df_allele$allele[df_allele$is_ref]) %>%
  slice_max(mean_frac, n = n_alt, with_ties = FALSE)

# transcript used: the one with the longest CDS among those whose CDS contains the cut
cds_by_tx <- split(cds, as.character(cds$transcript_id))
cds_by_tx <- cds_by_tx[any(overlapsAny(cds_by_tx, GRanges(chr, IRanges(cut_left, cut_left + 1))))]
stopifnot("the cut is not in the coding sequence" = length(cds_by_tx) > 0)
tx_id <- names(which.max(sum(width(cds_by_tx))))
tx_exons <- sort(exons[exons$transcript_id == tx_id], ignore.strand = TRUE)
tx_exons <- tx_exons[order(start(tx_exons))]
tx_cds <- cds_by_tx[[tx_id]]
minus_gene <- cds_strand == "-"
message("protein predicted on the transcript ", tx_id)

# mRNA as genomic positions and bases, in the genomic + orientation
mrna_pos <- unlist(map2(start(tx_exons), end(tx_exons), seq))
mrna_gen <- paste(as.character(getSeq(fa, GRanges(chr, ranges(tx_exons)))), collapse = "")
# window of the allele sequences (guide +/- flank_guide, limited to the exon of the cut)
ex_cut <- tx_exons[start(tx_exons) <= cut_left & end(tx_exons) >= cut_left + 1]
a_start <- max(g_start - flank_guide, start(ex_cut))
a_end <- min(g_end + flank_guide, end(ex_cut))
i_start <- match(a_start, mrna_pos) # index of the window in the mRNA
i_end <- match(a_end, mrna_pos)
# first base of the start codon (index in the mRNA in the transcript orientation)
start_pos <- if (!minus_gene) min(start(tx_cds)) else max(end(tx_cds))
start_idx <- if (!minus_gene) match(start_pos, mrna_pos) else length(mrna_pos) - match(start_pos, mrna_pos) + 1
window_idx <- if (!minus_gene) i_start else length(mrna_pos) - i_end + 1
stopifnot("the allele window overlaps the start codon: the prediction is not valid" = window_idx > start_idx + 2)

rc <- function(x) as.character(reverseComplement(DNAString(x)))
# translation of an mRNA (given in the genomic + orientation) from the start codon to the first stop ("*", kept in the output)
translate_mrna <- function(m) {
  if (minus_gene) m <- rc(m)
  m <- substr(m, start_idx, nchar(m))
  aa <- as.character(Biostrings::translate(DNAString(substr(m, 1, nchar(m) - nchar(m) %% 3))))
  stop_at <- regexpr("*", aa, fixed = TRUE)
  if (stop_at > 0) substr(aa, 1, stop_at) else aa
}
# from an allele of the table (guide orientation, "-" = deleted base, "+ins<bases>@<left>|<right>") to the edited mRNA
edit_mrna <- function(allele) {
  window <- sub("\\+ins.*$", "", allele)
  stopifnot("allele and window have different length: is the table of this guide?" = nchar(window) == a_end - a_start + 1)
  if (flip) window <- rc(window) # genomic + orientation
  bases <- strsplit(window, "")[[1]]
  ins <- str_match(allele, "\\+ins([ACGTN]+)@(-?\\d+)\\|")
  if (!is.na(ins[1, 1])) {
    left <- as.integer(ins[1, 3]) # guide position of the base on the left (guide orientation) of the insertion
    after <- if (!flip) g_start + left - 1 else g_end - left # genomic position after which the bases are inserted
    k <- after - a_start + 1
    bases[k] <- paste0(bases[k], if (flip) rc(ins[1, 2]) else ins[1, 2])
  }
  paste0(substr(mrna_gen, 1, i_start - 1), gsub("-", "", paste(bases, collapse = "")), substr(mrna_gen, i_end + 1, nchar(mrna_gen)))
}
# check: the reference allele has to give back the reference mRNA
stopifnot(edit_mrna(df_allele$allele[df_allele$is_ref][1]) == mrna_gen)

# short description of an allele
describe <- function(allele) {
  ins <- str_match(allele, "\\+ins([ACGTN]+)@(-?\\d+)\\|(-?\\d+)")
  n_del <- str_count(sub("\\+ins.*$", "", allele), "-")
  d <- c(if (n_del > 0) paste0(n_del, " bp deletion"),
         if (!is.na(ins[1, 1])) paste0(nchar(ins[1, 2]), " bp insertion (", ins[1, 2], ")"))
  if (length(d) == 0) "substitutions" else paste(d, collapse = " + ")
}

# one row per protein: reference and edited alleles
df_prot <- bind_rows(
  tibble(name = "reference", allele = NA_character_, reads_KO = NA_character_, protein = translate_mrna(mrna_gen)),
  df_top %>% transmute(name = paste0("allele ", row_number(), ": ", map_chr(allele, describe)), allele,
                       reads_KO = reads, protein = map_chr(allele, ~ translate_mrna(edit_mrna(.x))))
)
ref_aa <- strsplit(df_prot$protein[1], "")[[1]]
stopifnot("the reference translation does not end with a stop codon" = tail(ref_aa, 1) == "*")
# per protein: length (aa), first residue different from the reference (position by position, no alignment: meant for frameshifts)
df_prot <- df_prot %>%
  mutate(length_aa = nchar(sub("\\*$", "", protein)),
         stop_found = grepl("\\*$", protein),
         first_diff = map_int(strsplit(protein, ""), function(a) {
           i <- which(a != ref_aa[seq_along(a)] | seq_along(a) > length(ref_aa))
           if (length(i) == 0) NA_integer_ else i[1]
         }),
         n_ref_residues = coalesce(first_diff - 1L, length_aa))
print(df_prot %>% dplyr::select(name, reads_KO, length_aa, stop_found, first_diff))
print(df_prot$protein)
write_tsv(df_prot, paste0("../../out/table/10_protein_alleles_gRNA_", gene, ".tsv"))

# plot: one letter per residue, wrapped on lines of `aa_per_line` residues, the proteins one below the other
# residue that contains the cut: number of coding bases before the cut (transcript orientation)
cds_pos <- unlist(map2(start(tx_cds), end(tx_cds), seq))
n_before <- if (!minus_gene) sum(cds_pos <= cut_left) else sum(cds_pos >= cut_left + 1)
cut_res <- n_before %/% 3 + 1

n_prot <- nrow(df_prot)
df_aa <- df_prot %>%
  mutate(prot_i = row_number(), aa = strsplit(protein, "")) %>%
  dplyr::select(prot_i, name, aa) %>%
  unnest(aa) %>%
  group_by(prot_i) %>% mutate(res = row_number()) %>% ungroup() %>%
  mutate(status = case_when(aa == "*" ~ "stop",
                            prot_i == 1 ~ "same as reference",
                            res <= length(ref_aa) & aa == ref_aa[pmin(res, length(ref_aa))] ~ "same as reference",
                            .default = "different"),
         line = (res - 1) %/% aa_per_line + 1,
         col = (res - 1) %% aa_per_line + 1)
line_lab <- function(l) paste0("residues ", (l - 1) * aa_per_line + 1, "-", l * aa_per_line)
df_aa <- df_aa %>% mutate(line_f = factor(line_lab(line), levels = line_lab(sort(unique(line)))))
# the cut: a dashed line inside the residue that contains it (position of the cut within the codon)
df_cut <- tibble(line_f = factor(line_lab((cut_res - 1) %/% aa_per_line + 1), levels = levels(df_aa$line_f)),
                 col = (cut_res - 1) %% aa_per_line + 0.5 + (n_before %% 3) / 3)
row_names <- c("reference", paste("allele", seq_len(n_prot - 1)))

legend_txt <- paste(paste0(df_prot$name, " | ", df_prot$length_aa, " aa",
                           if_else(is.na(df_prot$first_diff) | is.na(df_prot$allele), "",
                                   paste0(", first different residue ", df_prot$first_diff)),
                           if_else(df_prot$stop_found, "", ", no stop found"),
                           if_else(is.na(df_prot$reads_KO), "", paste0(" | reads: ", df_prot$reads_KO))),
                    collapse = "\n")

ggplot(df_aa, aes(col, prot_i)) +
  geom_tile(aes(fill = status), colour = "white", width = 1, height = 0.9) +
  geom_text(aes(label = aa, colour = status), size = 2.6, family = "mono") +
  geom_segment(data = df_cut, aes(x = col, xend = col, y = 0.4, yend = n_prot + 0.6), linetype = "dashed", colour = col_guide, linewidth = 0.6) +
  facet_wrap(~ line_f, ncol = 1) +
  scale_fill_manual(values = c("same as reference" = "gray90", different = "#D62728", stop = "black"), name = NULL) +
  scale_colour_manual(values = c("same as reference" = "gray20", different = "white", stop = "white"), guide = "none") +
  scale_x_continuous(breaks = c(1, seq(10, aa_per_line, 10)), limits = c(0.5, aa_per_line + 0.5), expand = c(0, 0)) +
  scale_y_reverse(breaks = seq_len(n_prot), labels = row_names, expand = c(0, 0.1)) +
  labs(x = "residue in the line", y = NULL,
       title = paste0(gene, " - reference protein (", tx_id, ") and predicted proteins of the edited alleles"),
       subtitle = paste0(legend_txt, "\ndashed line = expected Cas9 cut (in codon ", cut_res, "); * = stop codon")) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "bottom", strip.background = element_blank(),
        strip.text = element_text(hjust = 0, size = 8), plot.subtitle = element_text(size = 8))
ggsave(paste0("../../out/plot/10_protein_alleles_gRNA_", gene, ".pdf"), width = 11,
       height = 2 + n_distinct(df_aa$line) * (0.22 * n_prot + 0.35))
