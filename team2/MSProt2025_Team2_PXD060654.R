# **Team 2: Hanna, Srividhya, Jannica**
#   
#   # 🧬 Capstone Project: R-based Mass Spectrometry Proteomics Workflow
#   hello
# This document provides a step-by-step guide for your capstone project on analyzing proteomics data using R.
# 
# ---
#   
#   ## 🔁 Workflow Overview
#   
#   ### 1. 📁 Dataset Acquisition  
#   **Input:** PRIDE PXD Identifier  
# ↓  
# (*Example: `PXD0123456`, accessible via [PRIDE](https://www.ebi.ac.uk/pride/)*)
# 
# Use the `rpx` package to access metadata and download a dataset that includes `.mzID` files from a published study.

## install.package("BiocManager")

#library(BiocManager)

# Bioconductor packages
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# 
# BiocManager::install("tidyverse")
# BiocManager::install("factoextra")
# BiocManager::install("msdata")
# BiocManager::install("mzR")
# BiocManager::install("rhdf5")
# BiocManager::install("rpx")
# BiocManager::install("MsCoreUtils")
# BiocManager::install("QFeatures")
# BiocManager::install("Spectra")
# BiocManager::install("ProtGenerics")
# BiocManager::install("PSMatch")
# BiocManager::install("pheatmap")
# BiocManager::install("limma")
# BiocManager::install("MSnID")
# BiocManager::install("RforMassSpectrometry/SpectraVis")
# 
# install.packages("tidyverse")

setwd("C:/Users/janni/OneDrive/Työpöytä/HU - Proteomics 2025/Capstone_project")

# Access dataset from ProteomeXChange
library("rpx")
px <- PXDataset("PXD060654")
pxfiles(px) # file names of the object

# Download mzID file
mzID <- pxget(px, "MS_2369_mock-vs-flg-combined_031424.mzid")

---
  
#   ### 2. 🧱 PSM Object Creation & Preprocessing  
#   **Goal:** Generate a PSM (Peptide-Spectrum Match) object from `.mzID` files  
# ↓  
# - Convert `.mzID` files into `PSM` objects
# - Assess:
#   - Number of decoy hits
# - Score distributions
# - PSM rank
# - Apply filtering based on FDR or identification score
# 
  
library(mzID)
library(MSnbase)
library(MSnID)

# Read downloaded mzID file
mzid_file <- "C:/Users/janni/AppData/Local/R/cache/R/rpx/39fc2ba81653_MS_2369_mock-vs-flg-combined_031424.mzid"
mzid <- mzID::mzID(mzid_file)

# Flatten into a data frame
psms <- flatten(mzid)
head(psms)

msnid <- MSnID()
msnid <- read_mzIDs(msnid, mzid_file)

psms_df <- psms(msnid)  # this returns a data.frame of PSMs
head(psms_df)

# Access number of decoy hits
table(psms_df$isDecoy) # Returned a table with 0 rows
names(psms_df)

grep("decoy|rev|reverse|xxx", psms_df$accession, ignore.case = TRUE, value = TRUE) # No decoys found in the accession column

grep("decoy|rev|reverse|xxx", psms_df$description, ignore.case = TRUE, value = TRUE) # No decoys found in descriptions, only HXXXD

# No target-decoy FDR estimation can be performed
# Plotting the PSM scores in SEQUEST:xcorr column
library(ggplot2)

ggplot(psms_df, aes(x = `SEQUEST:xcorr`)) +
  geom_density(fill = "lightblue") +
  labs(title = "SEQUEST XCorr Score Distribution")

# Plot PSM rank distribution in a bar plot
ggplot(psms_df, aes(x = factor(rank))) +
  geom_bar(fill = "steelblue") +
  labs(
    title = "PSM Rank Distribution",
    x = "Rank",
    y = "Count"
  ) +
  theme_minimal()

# Plot Scores by Rank using boxplot
ggplot(psms_df, aes(x = factor(rank), y = `SEQUEST:xcorr`)) +
  geom_boxplot(fill = "lightgreen") +
  labs(
    title = "XCorr Score by PSM Rank",
    x = "Rank",
    y = "SEQUEST:xcorr"
  ) +
  theme_minimal()

# Apply a filtering, SEQUEST:xcorr > 2, rank = 1
psms_high_conf <- subset(psms_df, `SEQUEST:xcorr` > 2.0 & rank == 1)

  
#   ### 3. 🧬 Protein & Peptide Identification  
#   **Goal:** Determine identified peptides and proteins  
# ↓  
# - Count the number of identified peptides and proteins
# - Review peptide-to-protein mapping:
#   - Razor proteins
# - Protein groups
# 
# ---
#   
#   ### 4. 🔄 QFeature Aggregation (Optional)  
#   **Goal:** Use quantification data if available  
# ↓  
# Use `aggregateFeatures()` from the `QFeatures` package to aggregate:
#   - PSMs ➡️ Peptides ➡️ Proteins
# 
# ---
#   
#   ### 5. 🧼 Normalization & Imputation (Optional)  
#   **Goal:** Correct for technical variation and handle missing values  
# ↓  
# - Normalize using `normalize()` or `normalize_vsn()`
# - Impute missing values with `impute()`
# 
# ---
#   
#   ### 6. 🧪 Protein Inference & Quantification (Optional)  
#   **Goal:** Summarize and quantify proteins  
# ↓  
# - Aggregate peptide-level intensities into protein-level quantities
# - Optionally annotate protein IDs using external databases (e.g., UniProt)
# 
# ---
#   
#   ### 7. 📊 Statistical Analysis (Optional)  
#   **Goal:** Identify differentially abundant proteins  
# ↓  
# - Perform statistical testing (e.g., using `test_diff()`)
# - Filter by:
#   - Log2 fold-change
# - Adjusted p-value (e.g., FDR < 0.05)
# 
# ---
#   
#   ### 8. 📈 Visualization & Export (Optional)  
#   **Goal:** Visualize data for interpretation  
# ↓  
# - Create visual summaries:
#   - PCA plots
# - Volcano plots
# - Heatmaps  
# - Use packages such as `ggplot2`, `DEP`, or `limma`
# 
# ---
#   
#   ### 9. 📖 Comparison with Published Results  
#   **Goal:** Benchmark your findings  
# ↓  
# Create a comparison table including:
#   - Number of spectra identified
# - Number of peptides and proteins
# - Key figures or results (if available)
# - Comments on reproducibility
# 
# ---
#   
#   ### 10. 📝 Final Report & Interpretation  
#   **Goal:** Reflect on your analysis  
# ↓  
# Submit the following:
#   - Source code (R script or RMarkdown) in your team folder
# - A 1–2 page report that includes:
#   - Background of the dataset and study
# - Summary table of key results
# - Discussion:
#   > Why do your results match or differ from the original publication?
#   
#   ---
#   
#   ✅ **Tip:** Convert your `.Rmd` into a clean `README.md` using `knitr::knit("README.Rmd")`.

