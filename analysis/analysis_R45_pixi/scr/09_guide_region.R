# AIM ---------------------------------------------------------------------
# look at the sequence of the RNA-seq reads around the target site of a guide RNA (default: C1QA), WT vs KO.
# The region shown is the whole guide (20 bp) + the PAM + `flank_guide` bases on each side of the guide.
#
# Steps:
# 1. the given sequence is searched in the gene. It can be the guide (protospacer, 5'-spacer-3' + PAM) or its reverse complement
#    (the sequence of the other strand): both are tested and the orientation followed by an NGG PAM is used.
#    The expected Cas9 cut site is 3 bases from the PAM (between position 17 and 18 of the protospacer).
# 2. per base: depth and fraction of reads that differ from the reference (mismatch or deletion) in the region.
# 3. per sample: fraction of reads with an insertion / deletion within `near` bases of the cut site.
# 4. alleles: sequence of the reads that cover the whole region (only its exonic part), counted per sample.
# All the plots and the allele sequences are in the orientation of the guide (5' -> 3' of the protospacer, then the PAM),
# as in the Sanger / ICE figure of the company, not in the genomic + orientation.
# The reads are RNA (mature mRNA, STAR BAM of the nf-core run): only exonic sequence is seen (the intronic bases of the region are grey),
# and an edited transcript can be degraded (NMD), so the fraction of edited reads is not the fraction of edited alleles of the clone.
#
# GLOSSARY (terms used below)
# - genomic position / coordinate: the number of a base along a chromosome of the reference genome (GRCh38), counted from 1 at the
#   start of the chromosome, e.g. chr1:22638841. A region is written chr:start-end. Coordinates always follow the + strand.
# - strand (+ / -): DNA has two complementary strands. The reference sequence is the + strand; the - strand is its reverse
#   complement (A<->T, C<->G, and the order is reversed). A gene or a guide can lie on either strand.
# - protospacer / guide: the 20 bases of the target that the guide RNA matches (5' -> 3'). It must be followed by the PAM.
# - PAM: 3 bases (NGG for SpCas9) right after the protospacer on the same strand. Cas9 cuts 3 bases before the PAM,
#   i.e. between base 17 and base 18 of the protospacer: this is the cut site, where the repair makes small insertions / deletions (indels).
# - exon / intron: exons are in the mature RNA, introns are removed by splicing. RNA-seq reads come (almost) only from exons.
# - BAM: the file with the reads aligned to the genome. Each alignment has a position, a CIGAR (e.g. 40M1D60M = 40 bases matching,
#   1 deleted base, 60 matching; I = inserted bases; N = a gap = an intron; S = soft clipped bases, not aligned) and a MAPQ
#   (255 for STAR = the read maps to one place only).
# - depth: the number of reads that cover a base. pileup: the count of the bases (A/C/G/T/deletion) of all the reads at each position.

# libraries ---------------------------------------------------------------
library(tidyverse)
library(GenomicAlignments)
library(Rsamtools)
library(Biostrings)
library(rtracklayer)
library(patchwork)

# parameters --------------------------------------------------------------
# the sequence of the guide as given by the company and the gene that contains it
guide <- "CTGGCATCCGGACAGGCATC"
gene <- "C1QA"
# optional command line arguments overwrite the two values above
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) guide <- toupper(args[1])
if (length(args) > 1) gene <- args[2]
spacer_strand <- NA # "+" or "-" to force the strand of the protospacer; NA = the strand with an NGG PAM
flank_guide <- 10 # bases shown on each side of the guide
near <- 3 # an indel is counted at the cut site if it is within this distance (bases) from the cut
min_mapq <- 255 # STAR: 255 = uniquely mapped
min_depth <- 20 # minimum number of reads at a base to show the fraction of non-reference reads
n_alleles <- 10 # number of alleles (most frequent over all the samples) drawn in the aligned allele plot
# where an indel in a run of the same base is placed (all the positions in the run give the same read):
# "cut" = the position nearest the Cas9 cut site (the first one on the left if two are equally near); "left" = the leftmost position (genomic + orientation)
indel_placement <- "cut"
# BAM files of the nf-core run (one per sample), gene annotation (exon coordinates) and genome sequence used for the alignment
bam_dir <- "/beegfs/scratch/ric.absinta/ric.absinta/pipeline/nfcore_bulkRNAseq_AlettaMG/results/star_salmon"
gtf_file <- "../../data/genes.filtered.gtf"
ref_dir <- "/beegfs/scratch/ric.absinta/ric.absinta/reference/nfcore/GRCh38"
# the genome fasta (reading a region does not load the whole genome: the .fai index gives the position of each chromosome in the file)
fa <- FaFile(file.path(ref_dir, "rsem/genome.fa"), index = file.path(ref_dir, "genome.fa.fai"))

# samples of the analysis (same table of 01_import_data.R): sample id, genotype (mut), treatment, clone
coldata <- read_csv("../../data/LUT_samples_final.csv", show_col_types = FALSE)

# locate the guide, the PAM and the cut site ------------------------------
# exons of the gene from the annotation: the coordinates of every exon of the gene (and of its chromosome)
gtf <- import(gtf_file)
exons <- gtf[gtf$type == "exon" & gtf$gene_id == gene]
stopifnot("no exons found for this gene" = length(exons) > 0)
exons_red <- reduce(exons) # union of the exons of all the transcripts
chr <- unique(as.character(seqnames(exons)))
# the whole gene (first base of the first exon to last base of the last exon) and its sequence (+ strand)
gene_gr <- GRanges(chr, IRanges(min(start(exons)), max(end(exons))))
gene_seq <- getSeq(fa, gene_gr)[[1]]

# where the given sequence is found in the genome (either strand): genomic span g_start-g_end
# (it is searched as it is, and as reverse complement, because the guide can be on the - strand: then the + strand contains its reverse complement)
hits <- c(start(matchPattern(DNAString(guide), gene_seq)), start(matchPattern(reverseComplement(DNAString(guide)), gene_seq)))
stopifnot("the sequence has to match exactly once in the gene" = length(hits) == 1)
# position of the first and last base of the guide in genomic coordinates (the match is relative to the start of the gene)
g_start <- start(gene_gr) + hits - 1
g_end <- g_start + nchar(guide) - 1

# PAM if the protospacer is on the + strand (3 bases after g_end) or on the - strand (3 bases before g_start, reverse complemented)
# the right strand is the one where these 3 bases are NGG
pam_plus <- as.character(getSeq(fa, GRanges(chr, IRanges(g_end + 1, g_end + 3))))
pam_minus <- as.character(reverseComplement(getSeq(fa, GRanges(chr, IRanges(g_start - 3, g_start - 1)))))
message("PAM if the protospacer is on the + strand: ", pam_plus, " | on the - strand: ", pam_minus)
if (is.na(spacer_strand)) {
  ok <- c("+" = grepl("^.GG$", pam_plus), "-" = grepl("^.GG$", pam_minus))
  stopifnot("no unique strand with an NGG PAM: set spacer_strand by hand" = sum(ok) == 1)
  spacer_strand <- names(ok)[ok]
}
flip <- spacer_strand == "-" # the guide is on the - strand: everything is shown as reverse complement
# genomic coordinates and sequence of the PAM
pam_gr <- if (!flip) GRanges(chr, IRanges(g_end + 1, g_end + 3)) else GRanges(chr, IRanges(g_start - 3, g_start - 1))
pam <- if (!flip) pam_plus else pam_minus
# the cut is between cut_left and cut_left + 1 (genomic coordinates): between base 17 and 18 of the protospacer
# (+ strand: 3 bases before the end of the guide; - strand: the guide runs backwards, so 3 bases after its first genomic base)
cut_left <- if (!flip) g_end - 3 else g_start + 2
# the protospacer sequence as the guide RNA reads it (5' -> 3'): reverse complement of the genome if the guide is on the - strand
spacer <- as.character(getSeq(fa, GRanges(chr, IRanges(g_start, g_end), strand = if (flip) "-" else "+")))
message(gene, " | protospacer ", spacer, " on the ", spacer_strand, " strand (", chr, ":", g_start, "-", g_end, ")",
        " | PAM ", pam, " | cut between ", cut_left, " and ", cut_left + 1,
        " | cut in an exon: ", any(overlapsAny(GRanges(chr, IRanges(cut_left, cut_left + 1)), exons)))

# region: guide + PAM + flank on each side, one row per base --------------
# the region to analyse: from `flank_guide` bases before the guide to `flank_guide` bases after it (the PAM is inside the right flank)
win <- GRanges(chr, IRanges(g_start - flank_guide, g_end + flank_guide))
# one row per base of the region, with the reference base and the features of the base:
# x = position along the guide (the plots show the guide 5' -> 3' from left to right, also if it is on the - strand),
# base = reference base in that orientation, in_guide / in_pam = the base is in the guide / PAM, exonic = the base is in an exon
ref_tbl <- tibble(pos = start(win):end(win), ref = strsplit(as.character(getSeq(fa, win)[[1]]), "")[[1]]) %>%
  mutate(# position in the guide orientation: 1-20 = protospacer, 21-23 = PAM, <1 = before the guide (5'), >23 = after the PAM
         x = if (!flip) pos - g_start + 1L else g_end - pos + 1L,
         base = if (!flip) ref else chartr("ACGT", "TGCA", ref), # base in the guide orientation
         in_guide = x >= 1 & x <= nchar(guide),
         in_pam = pos %in% start(pam_gr):end(pam_gr),
         exonic = overlapsAny(GRanges(chr, IRanges(pos, pos)), exons_red))
cut_x <- ref_tbl$x[ref_tbl$pos == cut_left] + if (!flip) 0.5 else -0.5 # the cut line in the x axis of the plots
message("region (guide orientation): ", paste(ref_tbl$base[order(ref_tbl$x)], collapse = ""),
        " | intronic bases: ", sum(!ref_tbl$exonic))

# the BAM file of a sample: <sample>.markdup.sorted.bam in the nf-core results
samples <- setNames(coldata$sample, coldata$sample)
bam_of <- function(s) file.path(bam_dir, paste0(s, ".markdup.sorted.bam"))

# per base: depth and fraction of reads different from the reference ------
# for each sample, the pileup counts how many reads show A / C / G / T / a deletion at each base of the region
# (bases with quality < 20 and reads that do not map to a single place are not counted).
# frac_nonref = 1 - (reads with the reference base) / (reads covering the base): it includes mismatches and deletions
df_base <- map_dfr(samples, function(s) {
  pu <- pileup(bam_of(s), scanBamParam = ScanBamParam(which = win),
               pileupParam = PileupParam(distinguish_strands = FALSE, distinguish_nucleotides = TRUE, max_depth = 1e5,
                                         min_base_quality = 20, min_mapq = min_mapq,
                                         include_deletions = TRUE, include_insertions = FALSE))
  ref_tbl %>%
    left_join(pu %>% group_by(pos) %>%
                summarise(depth = sum(count), n_del = sum(count[nucleotide == "-"]),
                          n_ref_pu = list(setNames(count, nucleotide)), .groups = "drop"), by = "pos") %>%
    mutate(depth = replace_na(depth, 0), n_del = replace_na(n_del, 0),
           # reads with the reference base: 0 if no read covers the position (NULL) or if no read has the reference base
           n_ref = map2_dbl(n_ref_pu, ref, ~ case_when(is.null(.x) ~ 0,
                                                       .y %in% names(.x) ~ sum(.x[names(.x) == .y]),
                                                       .default = 0)),
           n_alt = depth - n_ref - n_del, # mismatches (substitutions)
           frac_nonref = if_else(depth >= min_depth, 1 - n_ref / depth, NA_real_)) %>% # not shown where the depth is too low (introns)
    dplyr::select(-n_ref_pu)
}, .id = "sample") %>%
  left_join(coldata %>% dplyr::select(sample, sample_name, mut, treat, clone), by = "sample") %>%
  mutate(sample_name = factor(sample_name, levels = unique(sample_name[order(mut, treat, clone)])))
write_tsv(df_base, paste0("../../out/table/guide_region_per_base_", gene, ".tsv"))

# per sample: indels at the cut site and alleles ---------------------------
# alleles: the exonic part of the region (reads do not cover the intron), in the exon that contains the cut
ex_cut <- exons_red[overlapsAny(exons_red, GRanges(chr, IRanges(cut_left, cut_left + 1)))]
a_start <- max(start(win), start(ex_cut))
a_end <- min(end(win), end(ex_cut))
# reads "at the cut site": they have to cover cut_left-near+1 ... cut_left+near
n_start <- cut_left - near + 1
n_end <- cut_left + near
# for each sample, the alignments of the reads in the region are read from the BAM (with their sequence and MAPQ) and
# - the CIGAR gives the deletions (D) and insertions (I) of each read, with their position on the genome
# - the sequence of the reads that cover the whole allele window is extracted, to count the different alleles
res <- lapply(samples, function(s) {
  aln <- readGAlignments(bam_of(s), param = ScanBamParam(which = win, what = c("seq", "mapq")))
  aln <- aln[!is.na(mcols(aln)$mapq) & mcols(aln)$mapq >= min_mapq] # uniquely mapped reads only
  del <- cigarRangesAlongReferenceSpace(cigar(aln), pos = start(aln), ops = "D") # genomic positions of the deleted bases of each read
  ins <- cigarRangesAlongReferenceSpace(cigar(aln), pos = start(aln), ops = "I") # zero width: position of the insertion
  # reads covering the cut site
  cover_near <- start(aln) <= n_start & end(aln) >= n_end
  # reads with a deletion / insertion within `near` bases of the cut
  has_del <- any(start(del) <= n_end & end(del) >= n_start)
  has_ins <- any(start(ins) >= n_start & start(ins) <= n_end + 1)
  # alleles: reads covering the whole window; sequence in reference space (deletions = "-", insertions removed and flagged)
  # = the read bases lined up on the genome positions of the window (genomic + orientation, reversed later)
  full <- start(aln) <= a_start & end(aln) >= a_end
  ref_seq <- sequenceLayer(mcols(aln)$seq[full], cigar(aln)[full], to = "reference", D.letter = "-", N.letter = "-")
  allele <- substr(as.character(ref_seq), a_start - start(aln)[full] + 1, a_end - start(aln)[full] + 1)
  # insertions: the inserted bases are not in the reference, so they are not in the allele sequence above.
  # For the reads with an insertion in the window, the genomic position (the insertion is just before this base) and the
  # inserted bases are taken from the CIGAR (I) and the read sequence, and added to the allele as "+ins<bases>@<position>"
  ins_seq <- extractAt(mcols(aln)$seq[full], cigarRangesAlongQuerySpace(cigar(aln)[full], ops = "I")) # inserted bases of each read
  ins_ref <- ins[full] # position of each insertion on the genome
  ins_info <- map_chr(seq_along(ins_ref), function(i) {
    p <- start(ins_ref[[i]])
    k <- which(p >= a_start & p <= a_end + 1) # insertions in the window (the first one is kept)
    if (length(k) > 0) paste0("+ins", as.character(ins_seq[[i]][k[1]]), "@", p[k[1]]) else ""
  })
  # indel: number of reads covering the cut and how many have a deletion / insertion; allele: number of reads of each sequence
  list(indel = tibble(n_cover = sum(cover_near), n_del = sum(cover_near & has_del), n_ins = sum(cover_near & has_ins),
                      n_indel = sum(cover_near & (has_del | has_ins))),
       allele = tibble(allele = paste0(allele, ins_info)) %>% dplyr::count(allele, name = "n_reads"))
})

# table of the indels at the cut: reads covering the cut, with a deletion, with an insertion, with any indel, fraction
df_indel <- map_dfr(res, "indel", .id = "sample") %>%
  mutate(frac_indel = n_indel / n_cover) %>%
  left_join(coldata %>% dplyr::select(sample, sample_name, mut, treat, clone), by = "sample")
print(df_indel, n = Inf)
write_tsv(df_indel, paste0("../../out/table/guide_region_indels_", gene, ".tsv"))

# table of the alleles: sequence of the window (guide orientation), reads and fraction of the reads of each sample
# the allele sequences are reverse complemented to the guide orientation ("-" = deleted base).
# An insertion is written after the bases of the window as "+ins<inserted bases>@<left>|<right>": the bases are inserted between the
# guide positions <left> and <right> (e.g. +insC@17|18 = a C between base 17 and 18 of the guide, which is the cut site)
ref_allele <- ref_tbl %>% filter(pos >= a_start, pos <= a_end) %>% arrange(x) %>% pull(base) %>% paste(collapse = "")
# reference bases of the allele window in the genomic + orientation (used to left-align the indels below)
ref_win_gen <- ref_tbl %>% filter(pos >= a_start, pos <= a_end) %>% arrange(pos) %>% pull(ref)
df_allele <- map_dfr(res, "allele", .id = "sample") %>%
  mutate(# split the allele in the bases of the window and the insertion (inserted bases and genomic position), if any
         ins_base = str_match(allele, "\\+ins([ACGTN]+)@")[, 2],
         ins_pos = as.integer(str_match(allele, "@(\\d+)$")[, 2]),
         allele = sub("\\+ins.*$", "", allele),
         # An insertion / deletion inside a run of the same base (e.g. the CC next to the cut) gives exactly the same read wherever
         # it is placed in the run, and the aligner chooses the placement read by read: the same event would be split in several alleles.
         # All the indels are first moved to the leftmost equivalent position (genomic + orientation), the standard left-alignment of
         # GATK / bcftools norm / vt, so that the same event is one allele. If indel_placement = "cut", the indel is then moved to the
         # equivalent position that is nearest the cut site.
         # insertion: it can move left by one base if the base before it is the same as the (last) inserted base
         ins_norm = pmap(list(ins_pos, ins_base, allele), function(p, b, a) {
           while (!is.na(b) && p - 1L >= a_start && substr(a, p - a_start, p - a_start) == ref_win_gen[p - a_start] &&
                  substr(b, nchar(b), nchar(b)) == ref_win_gen[p - a_start]) {
             b <- paste0(substr(b, nchar(b), nchar(b)), substr(b, 1, nchar(b) - 1)) # the inserted bases rotate
             p <- p - 1L
           }
           if (indel_placement == "cut" && !is.na(b)) {
             # equivalent positions to the right: the insertion can move right by one base if the base after it is the same as the
             # first inserted base. The insertion is just before position p, the cut is just before position cut_left + 1
             cand_p <- p
             cand_b <- b
             while (p <= a_end && substr(a, p - a_start + 1L, p - a_start + 1L) == ref_win_gen[p - a_start + 1L] &&
                    substr(b, 1, 1) == ref_win_gen[p - a_start + 1L]) {
               b <- paste0(substr(b, 2, nchar(b)), substr(b, 1, 1))
               p <- p + 1L
               cand_p <- c(cand_p, p)
               cand_b <- c(cand_b, b)
             }
             best <- which.min(abs(cand_p - (cut_left + 1L))) # the nearest to the cut (the leftmost if tied)
             p <- cand_p[best]
             b <- cand_b[best]
           }
           list(pos = p, base = b)
         }),
         ins_pos = map_int(ins_norm, "pos"),
         ins_base = map_chr(ins_norm, "base"),
         # deletion (a single run of "-"): it can move left by one base if the base before it is the same as the last deleted base
         allele = map_chr(allele, function(a) {
           ch <- strsplit(a, "")[[1]]
           d <- which(ch == "-")
           if (length(d) == 0 || length(d) == length(ch) || any(diff(d) != 1)) return(a) # no deletion, whole window or several runs: unchanged
           s <- d[1]
           e <- d[length(d)]
           while (s > 1 && ch[s - 1] == ref_win_gen[s - 1] && ref_win_gen[s - 1] == ref_win_gen[e]) {
             ch[s - 1] <- "-"
             ch[e] <- ref_win_gen[e]
             s <- s - 1
             e <- e - 1
           }
           if (indel_placement == "cut") {
             # equivalent positions to the right: the deletion can move right by one base if the base after it is the same as the
             # first deleted base. ic = window index of the base just left of the cut; gap = bases between the deletion and the cut
             ic <- cut_left - a_start + 1L
             chs <- list(ch)
             cand_s <- s
             cand_e <- e
             while (e < length(ch) && ch[e + 1] == ref_win_gen[e + 1] && ref_win_gen[s] == ref_win_gen[e + 1]) {
               ch[s] <- ref_win_gen[s]
               ch[e + 1] <- "-"
               s <- s + 1
               e <- e + 1
               chs[[length(chs) + 1]] <- ch
               cand_s <- c(cand_s, s)
               cand_e <- c(cand_e, e)
             }
             gap <- if_else(cand_e <= ic, ic - cand_e, if_else(cand_s >= ic + 1L, cand_s - (ic + 1L), 0L))
             ch <- chs[[which.min(gap)]] # the nearest to the cut (the leftmost if tied)
           }
           paste(ch, collapse = "")
         }),
         # reverse complement to the guide orientation (the inserted bases too)
         allele = if (flip) vapply(strsplit(chartr("ACGT", "TGCA", allele), ""), function(z) paste(rev(z), collapse = ""), "") else allele,
         ins_base = if (flip) if_else(is.na(ins_base), NA_character_,
                                      vapply(strsplit(chartr("ACGT", "TGCA", ins_base), ""), function(z) paste(rev(z), collapse = ""), "")) else ins_base,
         # the insertion is between the bases at ins_pos - 1 and ins_pos: guide position of the base on the left of the insertion
         ins_left = pmin(ref_tbl$x[match(ins_pos - 1L, ref_tbl$pos)], ref_tbl$x[match(ins_pos, ref_tbl$pos)]),
         # allele = bases of the window + "+ins<bases>@<left>|<right>" = inserted bases between the guide positions <left> and <right>
         allele = if_else(is.na(ins_base), allele, paste0(allele, "+ins", ins_base, "@", ins_left, "|", ins_left + 1))) %>%
  dplyr::select(-ins_base, -ins_pos, -ins_left) %>%
  group_by(sample, allele) %>% summarise(n_reads = sum(n_reads), .groups = "drop") %>%
  group_by(sample) %>% mutate(frac = n_reads / sum(n_reads), depth_window = sum(n_reads)) %>% ungroup() %>%
  mutate(is_ref = allele == ref_allele) %>% # is_ref = the allele is identical to the reference
  left_join(coldata %>% dplyr::select(sample, sample_name, mut, treat, clone), by = "sample")
write_tsv(df_allele, paste0("../../out/table/guide_region_alleles_", gene, ".tsv"))
message("allele window (genomic ", a_start, "-", a_end, "), reference in the guide orientation: ", ref_allele)
# the 15 most frequent alleles over all the samples
print(df_allele %>% group_by(allele, is_ref) %>% summarise(n_reads = sum(n_reads), .groups = "drop") %>% arrange(desc(n_reads)) %>% head(15))

# plots -------------------------------------------------------------------
cols <- c(WT = "gray40", KO = "red")
# elements shared by the per base plots: blue = guide, orange = PAM, dashed = cut site, grey = intronic bases (no mature RNA reads)
rect_x <- function(sel) c(min(ref_tbl$x[sel]) - 0.5, max(ref_tbl$x[sel]) + 0.5) # left and right edge of a set of bases
guide_rects <- list(
  annotate("rect", xmin = rect_x(ref_tbl$in_guide)[1], xmax = rect_x(ref_tbl$in_guide)[2], ymin = -Inf, ymax = Inf, fill = "steelblue", alpha = 0.12),
  annotate("rect", xmin = rect_x(ref_tbl$in_pam)[1], xmax = rect_x(ref_tbl$in_pam)[2], ymin = -Inf, ymax = Inf, fill = "orange", alpha = 0.25),
  geom_rect(data = ref_tbl %>% filter(!exonic), aes(xmin = x - 0.5, xmax = x + 0.5, ymin = -Inf, ymax = Inf), inherit.aes = FALSE, fill = "gray70", alpha = 0.4),
  geom_vline(xintercept = cut_x, linetype = "dashed"))
# the x axis shows the reference base at each position (guide orientation) instead of numbers
x_scale <- scale_x_continuous(breaks = ref_tbl$x, labels = ref_tbl$base, expand = expansion(add = 0.5))

# fraction of reads with a different base (mismatch) or a deleted base, per base: one line per sample
# a peak = many reads of the sample have a different base (or a deleted base) at that position.
# Insertions are not counted: an inserted base is not at any position of the reference, so the reads with an insertion look
# like the reference at every base here (they are in the allele plots). The reads are those covering the base, the alleles
# are counted on the reads covering the whole window, so the fractions of the two kinds of plot are not expected to be equal.
p_frac <- ggplot(df_base, aes(x, frac_nonref, col = mut, group = sample)) +
  guide_rects +
  geom_line(alpha = 0.5) + geom_point(size = 0.8, alpha = 0.6) +
  scale_color_manual(values = cols) + x_scale +
  labs(x = NULL, y = "fraction of reads with a different or deleted base",
       title = paste0(gene, " - guide ", spacer, " (", spacer_strand, " strand), PAM ", pam),
       subtitle = paste0("per base: mismatches + deleted bases; insertions are not counted (they are not a base of the reference). ",
                         "blue = guide, orange = PAM, dashed = cut site, grey = intron")) +
  theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank())

# depth (number of reads covering each base), to judge the previous plot: a fraction is reliable only where the depth is high
p_depth <- ggplot(df_base, aes(x, depth, col = mut, group = sample)) +
  guide_rects + geom_line(alpha = 0.5) +
  scale_color_manual(values = cols) + x_scale +
  labs(x = "reference base in the guide orientation (5' -> 3')", y = "depth (unique reads)") +
  theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank())
ggsave(paste0("../../out/plot/guide_region_per_base_", gene, ".pdf"),
       p_frac / p_depth + plot_layout(heights = c(2, 1), guides = "collect"), width = 11, height = 7)

# fraction of reads with an indel at the cut site, per sample (the number is the reads covering the site)
p_indel <- ggplot(df_indel, aes(x = reorder(sample_name, frac_indel), y = frac_indel, fill = mut)) +
  geom_col() +
  geom_text(aes(label = n_cover), hjust = -0.1, size = 2.5) +
  scale_fill_manual(values = cols) +
  coord_flip() + scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = paste0("fraction of reads with an indel within ", near, " bases of the cut (number = reads covering the site)")) +
  theme_bw()
ggsave(paste0("../../out/plot/guide_region_indels_", gene, ".pdf"), p_indel, width = 7, height = 4)

# alleles: the 6 most frequent over all the samples, the others pooled
# one bar per sample: the fraction of its reads that carry each allele (the reference allele is the unedited sequence)
top <- df_allele %>% group_by(allele) %>% summarise(n = sum(n_reads), .groups = "drop") %>% slice_max(n, n = 6, with_ties = FALSE) %>% pull(allele)
p_allele <- df_allele %>%
  mutate(allele_grp = factor(if_else(allele %in% top, allele, "other"), levels = c(top, "other"))) %>%
  group_by(sample_name, mut, allele_grp) %>% summarise(frac = sum(frac), .groups = "drop") %>%
  ggplot(aes(sample_name, frac, fill = allele_grp)) +
  geom_col() + coord_flip() +
  scale_fill_brewer(palette = "Set2", name = paste0("allele (guide orientation), reference = ", ref_allele)) +
  facet_grid(mut ~ ., scales = "free_y", space = "free_y") +
  labs(x = NULL, y = "fraction of the reads covering the window") +
  theme_bw() + theme(legend.position = "bottom", legend.direction = "vertical")
ggsave(paste0("../../out/plot/guide_region_alleles_", gene, ".pdf"), p_allele, width = 8, height = 6)

# alleles aligned to the reference ----------------------------------------
# the most frequent alleles drawn base by base, to see where they differ from the reference.
# left: one row per allele, one tile per base of the window (guide orientation, 5' -> 3'):
#   grey = same base as the reference, red = different base, black = deleted base.
#   purple = bases inserted in the allele: they are not in the reference, so they are drawn between the two bases where they are inserted.
#   Insertions and deletions in a run of the same base (e.g. the CC next to the cut) are equivalent wherever they are placed in the run:
#   they are all shown at the same position, the one nearest the cut (indel_placement = "cut") or the leftmost one of the genome ("left").
#   The guide (blue), the PAM (orange) and the expected cut site (dashed line) are marked on top.
# right: % of the reads of each sample (number of reads covering the window in brackets) that carry the allele
x_vec <- ref_tbl %>% filter(pos >= a_start, pos <= a_end) %>% arrange(x) %>% pull(x) # x of the bases of the allele window, left to right
ref_chars <- ref_tbl %>% filter(pos >= a_start, pos <= a_end) %>% arrange(x) %>% pull(base) # reference bases, same order
# the reference allele and the most frequent alleles over all the samples
top_al <- unique(c(ref_allele, df_allele %>% group_by(allele) %>% summarise(n = sum(n_reads), .groups = "drop") %>%
                     slice_max(n, n = n_alleles, with_ties = FALSE) %>% pull(allele)))
n_rows <- length(top_al) + 1 # the last row is for all the other alleles
# description of each allele: bases deleted, mismatches, insertion (the allele has the suffix "+ins" if it has an insertion)
rows_al <- tibble(allele = top_al, row = seq_along(top_al)) %>%
  mutate(has_ins = grepl("\\+ins", allele), seq = sub("\\+ins.*$", "", allele),
         ins_base = str_match(allele, "\\+ins([ACGTN]+)@")[, 2], # inserted bases
         ins_left = as.numeric(str_match(allele, "@(-?\\d+)\\|")[, 2]), # guide position of the base on the left of the insertion
         n_del = str_count(seq, "-"),
         n_mm = map_int(seq, ~ sum(strsplit(.x, "")[[1]] != ref_chars & strsplit(.x, "")[[1]] != "-")),
         desc = pmap_chr(list(n_del, n_mm, ins_base), function(d, m, i) {
           # a read that has no base in the whole window (splice gap or large deletion) is all "-"
           v <- c(if (d == length(x_vec)) "no base in the window (gap)" else if (d > 0) paste0(d, " bp deleted"),
                  if (m > 0) paste0(m, " mismatch"), if (!is.na(i)) paste0(i, " inserted"))
           if (length(v) == 0) "reference" else paste(v, collapse = " + ")
         }),
         label = paste0(row, ". ", desc))
stopifnot("allele length different from the window" = all(nchar(rows_al$seq) == length(x_vec)))
row_labels <- c(rows_al$label, "other alleles")

# one tile per allele and base: same as the reference, different base or deleted base
tiles <- rows_al %>%
  mutate(base = map(seq, ~ strsplit(.x, "")[[1]]), x = map(seq, ~ x_vec)) %>%
  dplyr::select(row, base, x) %>%
  unnest(c(base, x)) %>%
  mutate(ref_base = ref_chars[match(x, x_vec)],
         status = case_when(base == "-" ~ "deletion", base == ref_base ~ "match", .default = "mismatch"))

# guide, PAM and cut, drawn in a track above the alleles (y < 0) and as a dashed line through the alleles
# inserted bases: a purple tile between the base at ins_left and the next one, with the inserted base(s)
ins_tiles <- rows_al %>% filter(has_ins) %>%
  transmute(row, x = ins_left + 0.5, base = ins_base, status = "insertion")
p_seq <- ggplot(tiles, aes(x, row)) +
  geom_tile(aes(fill = status), colour = "white", width = 1, height = 0.9) +
  geom_text(aes(label = base, colour = status), size = 2.6) +
  geom_tile(data = tibble(x = x_vec, row = n_rows), fill = "gray97", colour = "white", width = 1, height = 0.9) + # other alleles: not drawn
  geom_tile(data = ins_tiles, aes(fill = status), colour = "white", width = 0.75, height = 0.9) +
  geom_text(data = ins_tiles, aes(label = base, colour = status), size = 2.6, fontface = "bold") +
  annotate("rect", xmin = rect_x(ref_tbl$in_guide)[1], xmax = rect_x(ref_tbl$in_guide)[2], ymin = -1.5, ymax = -0.8, fill = "steelblue", alpha = 0.5) +
  annotate("text", x = mean(rect_x(ref_tbl$in_guide)), y = -1.15, label = "guide", size = 2.8) +
  annotate("rect", xmin = rect_x(ref_tbl$in_pam)[1], xmax = rect_x(ref_tbl$in_pam)[2], ymin = -1.5, ymax = -0.8, fill = "orange", alpha = 0.7) +
  annotate("text", x = mean(rect_x(ref_tbl$in_pam)), y = -1.15, label = "PAM", size = 2.8) +
  annotate("segment", x = cut_x, xend = cut_x, y = -1.5, yend = n_rows + 0.5, linetype = "dashed", linewidth = 0.6) +
  annotate("text", x = cut_x, y = -1.9, label = "cut", size = 2.8) +
  scale_fill_manual(values = c(match = "gray90", mismatch = "#D62728", deletion = "black", insertion = "#7B3294"), name = "base vs reference",
                    labels = c(match = "same", mismatch = "different", deletion = "deleted", insertion = "inserted")) +
  scale_colour_manual(values = c(match = "gray35", mismatch = "white", deletion = "white", insertion = "white"), guide = "none") +
  scale_x_continuous(breaks = x_vec[x_vec %% 5 == 0 | x_vec == 1], expand = expansion(add = 0.5)) +
  scale_y_reverse(breaks = seq_len(n_rows), labels = row_labels, limits = c(n_rows + 0.5, -2.2), expand = c(0, 0)) +
  labs(x = "position in the guide (1-20 = guide, 21-23 = PAM), 5' -> 3'", y = NULL,
       title = paste0(gene, " - most frequent alleles aligned to the reference"),
       subtitle = paste0("reference (guide orientation): ", ref_allele)) +
  theme_bw() + theme(panel.grid = element_blank(), legend.position = "bottom")

# % of the reads of each sample with each allele (the alleles that are not shown are pooled in the last row)
sample_levels <- coldata %>% arrange(factor(mut, levels = c("WT", "KO")), treat, clone) %>% pull(sample)
df_heat <- expand_grid(sample = sample_levels, row = seq_len(n_rows)) %>%
  left_join(coldata %>% dplyr::select(sample, mut), by = "sample") %>%
  left_join(df_allele %>%
              mutate(row = match(if_else(allele %in% top_al, allele, "other"), c(top_al, "other"))) %>%
              group_by(sample, row) %>% summarise(frac = sum(frac), n_window = dplyr::first(depth_window), .groups = "drop"),
            by = c("sample", "row")) %>%
  group_by(sample) %>% mutate(n_window = max(n_window, na.rm = TRUE)) %>% ungroup() %>%
  mutate(frac = replace_na(frac, 0),
         pct = if_else(frac >= 0.005, sprintf("%.0f", 100 * frac), if_else(frac > 0, "<1", "")),
         mut = factor(mut, levels = c("WT", "KO")),
         sample_lab = factor(paste0(sample, "\n(", n_window, ")"), levels = unique(paste0(sample, "\n(", n_window, ")")[order(match(sample, sample_levels))])))
p_heat <- ggplot(df_heat, aes(sample_lab, row, fill = frac)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = pct), size = 2.4) +
  facet_grid(. ~ mut, scales = "free_x", space = "free_x") +
  scale_fill_gradient(low = "white", high = "steelblue", limits = c(0, 1), name = "fraction of reads") +
  scale_y_reverse(breaks = NULL, limits = c(n_rows + 0.5, -2.2), expand = c(0, 0)) +
  labs(x = NULL, y = NULL, title = "% of the reads of each sample") +
  theme_bw() + theme(panel.grid = element_blank(), legend.position = "bottom", axis.text.x = element_text(size = 7))
ggsave(paste0("../../out/plot/guide_region_alleles_aligned_", gene, ".pdf"),
       p_seq + p_heat + plot_layout(widths = c(3, 1.4)), width = 16, height = 0.4 * n_rows + 3.5)
