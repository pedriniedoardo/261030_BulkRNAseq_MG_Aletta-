# AIM ---------------------------------------------------------------------
# read in the raw table of counts and build the dds object

# libraries ---------------------------------------------------------------
library("tidyverse")
library("DESeq2")

# read in the data --------------------------------------------------------
test <- readRDS("../../out/object/raw_count_salmon_GTF.rds")

# read in the metadata
LUT_samples <- read_csv("../../data/samplesheet.csv")

# samples removed from the analysis (see 00b_explore_QC.R):
# 24A8, 6A8: half the library size / ~10.7k genes detected, extremes of PC1 (technical failures)
# 6B3, 6B7: the only Iron / Iron_myelin samples (no replicate, cannot be tested) also the quality is not great
drop_final <- c("24A8", "6A8", "6B3", "6B7")

# keep only the samples with counts (and in the samplesheet) and not in drop_final
LUT_samples <- LUT_samples %>%
  dplyr::filter(sample_id %in% colnames(test),
                !sample_id %in% drop_final)

mat_exp <- test %>%
  as.data.frame() %>%
  dplyr::select(LUT_samples$sample_id) %>%
  as.matrix()

# match the order of the samples in the matrix with the samplesheet
coldata <- data.frame(sample = colnames(mat_exp)) %>%
  left_join(LUT_samples, by = c("sample" = "sample_id")) %>%
  mutate(rowname = sample,
         # sample_name is used as label in the QC plots
         sample_name = paste(sample, clone, mut, time, treat, sep = "_")) %>%
  column_to_rownames("rowname")

write.csv(coldata, file = "../../data/LUT_samples_final.csv", row.names = T)

# define the model --------------------------------------------------------
# ~ clone + treat + mut; the reference levels are BL and WT, so the coefficients are
# treat_Myelin_vs_BL and mut_KO_vs_WT (clone as blocking factor)
coldata <- coldata %>%
  mutate(clone = factor(clone),
         treat = factor(treat, levels = c("BL", "Myelin")),
         mut = factor(mut, levels = c("WT", "KO")))
design <- ~ clone + treat + mut

saveRDS(design, file = "../../out/object/design.rds")

# build the object --------------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = mat_exp,
                              colData = coldata,
                              design = design)

saveRDS(dds, file = "../../out/object/dds.rds")
