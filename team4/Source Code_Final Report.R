title: "Capstone Project: R-based Mass Spectrometry Proteomics Workflow"
author: "Team 4: Hao, Tiialotta, Oluwatosin"

  ## Library
library(rpx)
library(msdata)
library(rpx)
library(Spectra)
library(tidyverse)
library(MSnbase)
library(SpectraVis)
library(tidyverse)
library(mzR)
library(PSMatch)
library(pheatmap)
library(QFeatures)
library(sageR)
library(SummarizedExperiment)
library(factoextra)
library(mzID)
library(S4Vectors)  # ensures DataFrame is available
library(cleaver)


### 1. 📁 Dataset Acquisition  

px <- PXDataset("PXD045506") 

idf <- pxget(px, grep("mzid", pxfiles(px)))   

library(R.utils)
gunzip(idf, remove=F)

### 2. PSM Object Creation & Preprocessing  

idf_unzipped <- "/Users/tiitu1/Library/Caches/org.R-project.R/R/rpx/4e51263217e8_20230908_NTUHMCLin_SAS_MGAT5_WT_KO_stepHCD.mzid"

psms <- PSM(idf_unzipped, parser = 'mzID')

dim(psms)
names(psms)

hist(psms$`byonic:score`, breaks = 100)

#There are a lot of hits with low score (between 20 and 40)
hist(psms$`byonic:score`, breaks = 100)

#Try search for decoys
idxs_psms_score_20_40 <- which(psms@listData[["byonic:score"]] < 40 & psms@listData[["byonic:score"]] > 20)

#Check accession of low score hits
psms@listData[["accession"]][idxs_psms_score_20_40]

#Check description of low score hits
psms@listData[["description"]][idxs_psms_score_20_40]

table(psms$rank)
#Seems that there are no decoys in this exported mzID file -> MaxQuant analysis


### 3. Protein & Peptide Identification  
#MaxQuant analysis of raw files

psm_file <- "C:/Users/oljabe/Documents/Courses/MSProt2025/evidence.txt"   
psms_mq <- read.delim(psm_file)
colnames(psms_mq)  

# make wider table to produce two quantity columns
psms_mq_wide <- psms_mq %>%
  pivot_wider(names_from = Raw.file, values_from = Intensity)

# make wider table to produce two quantity columns
psms_mq_wide <- psms_mq %>%
  pivot_wider(names_from = Raw.file, values_from = Intensity)

# rename quant columns
colnames(psms_mq_wide)[72:73] <- c('WT_quant', 'KO_quant')


### 4. QFeature Aggregation use QFeatures - PSMs ➡️ Peptides ➡️ Proteins

# generate QFeatures object with PSM file from MaxQuant
qf <- readQFeatures(psms_mq_wide,             
                    ecol = 72:73,
                    name = "psms",
                    fnames = 'Sequence')     

qf <- zeroIsNA(qf, i = 'psms')

# assess NA portion
nNA(qf, i = 'psms')

qf <- filterNA(qf, i = 'psms', pNA = 1/2)


### 5. Normalization & Imputation (Optional)  use QFeatures 
# imputing psms
qf <- QFeatures::impute(qf, i = 'psms', method = 'bpca', name = 'psms_imputed')

table(is.na(rowData(qf[["psms_imputed"]])$Sequence))


### 6. Protein Inference & Quantification (Optional)  
qf <- aggregateFeatures(qf,
                        i = "psms_imputed",
                        fcol = "Sequence",   #aggregate by peptide sequence
                        name = "peptides",
                        fun = colMeans)
colnames(df)

# Extract intensity data from the 'peptides' assay
intensity_data <- assay(qf[["peptides"]])

# Convert to long format for ggplot
df <- as.data.frame(intensity_data) %>%
  rownames_to_column("Peptide") %>%
  pivot_longer(-Peptide, names_to = "Sample", values_to = "Intensity")

# Plot
ggplot(df, aes(x = Intensity, color = Sample)) +
  geom_density() +
  theme_minimal() +
  labs(title = "Density Plot Before Normalization",
       x = "Intensity",
       y = "Density")

# pre-reprocessing
# log transform
qf <- logTransform(qf, base = 2, i = "peptides", name = "peptideLog")

# reverse and contaminants
qf <- filterFeatures(qf, ~ Reverse != "+", i = 'peptideLog')
qf <- filterFeatures(qf, ~ Potential.contaminant != "+", i = 'peptideLog')

# normalisation
qf <- normalize(qf, 
                i = "peptideLog", 
                name = "peptideNorm", 
                method = "center.median")

# Convert to long format for plotting
df_density <- as.data.frame(assay(qf[["peptideNorm"]])) |>
  pivot_longer(cols = everything(), names_to = "Sample", values_to = "Intensity")

ggplot(df_density, aes(x = Intensity, color = Sample)) +
  geom_density() +
  theme_minimal() +
  labs(title = "Density Plot After Normalization",
       x = "Normalized Intensity",
       y = "Density")


qf <- aggregateFeatures(qf,
                        i = "peptideNorm",
                        fcol = "Leading.razor.protein",        #aggregate to proteins or leading razor proteins
                        name = "proteins")


#psms -> peptides -> proteins


### 7. Statistical Analysis 
# create experimental design
exp_design <- data.frame(sample = c('WT_quant', 'KO_quant'), group = c('WT', 'KO'))
exp_design <- exp_design %>% column_to_rownames('sample')
colData(qf) <- exp_design
colData(qf) <- DataFrame(exp_design)  #use Dataframe not data.frame

colnames(qf)
rownames(exp_design)

# Extract intensity matrix
prot_mat <- assay(qf[["proteins"]])    #or replace with "peptideNorm" if more appropriate
prot_mat2 <- assay(qf[["peptideNorm"]])
# Check your matrix
head(prot_mat)

#log-2 FC
# Convert to data frame and add protein IDs
logFC_df <- as.data.frame(prot_mat)
logFC_df$protein <- rownames(logFC_df)

# Compute log2 fold-change: KO - WT
logFC_df$log2FC <- logFC_df$KO_quant - logFC_df$WT_quant   #KO on the right (upreg. protein in KO, downreg. in WT)
#WT on the left (upreg protein in WT, downreg in KO)

# View top changes
logFC_df |> dplyr::arrange(desc(abs(log2FC))) |> head(10)

#No p-value we used log2FC
ggplot(logFC_df, aes(x = log2FC)) +
  geom_histogram(bins = 50, fill = "steelblue") +
  labs(title = "Log2 Fold-Change (WT vs KO)",
       x = "Log2 Fold-Change", y = "Number of Proteins") +
  theme_minimal()

#alternatively
ggplot(logFC_df, aes(x = log2FC)) +
  geom_boxplot(bins = 50, fill = "steelblue") +
  labs(title = "Log2 Fold-Change (KO vs WT)",
       x = "Log2 Fold-Change", y = "Number of Proteins") +
  theme_minimal()

#alternatively
ggplot(logFC_df, aes(x = log2FC)) +
  geom_density(bins = 50, fill = "steelblue") +
  labs(title = "Log2 Fold-Change (KO vs WT)",
       x = "Log2 Fold-Change", y = "Number of Proteins") +
  theme_minimal()

nrow(qf[["peptides"]])  # Number of peptide features

#Shifting proteins
top_prots <- logFC_df |> 
  dplyr::filter(abs(log2FC) > 2)

ggplot(logFC_df, aes(x = log2FC, y = 0)) +
  geom_point(alpha = 0.6) +
  geom_text(data = top_prots, aes(label = protein), 
            vjust = -1, size = 3, check_overlap = TRUE) +
  labs(title = "Top Shifting Proteins", x = "Log2 Fold-Change", y = "") +
  theme_minimal()