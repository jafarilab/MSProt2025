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

install.packages("BiocManager")

library(BiocManager)

# Bioconductor packages
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
#
BiocManager::install("tidyverse")
BiocManager::install("factoextra")
BiocManager::install("msdata")
BiocManager::install("mzR")
BiocManager::install("rhdf5")
BiocManager::install("rpx")
BiocManager::install("MsCoreUtils")
BiocManager::install("QFeatures")
BiocManager::install("Spectra")
BiocManager::install("ProtGenerics")
BiocManager::install("PSMatch")
BiocManager::install("pheatmap")
BiocManager::install("limma")
BiocManager::install("MSnID")
BiocManager::install("remotes")
BiocManager::install("RforMassSpectrometry/SpectraVis")

library(rpx)
library(msdata)
library(Spectra)
library(MSnbase)
library(SpectraVis)
library(tidyverse)
library(mzR)
library(PSMatch)
library(mzID)
library(Biostrings)
library(QFeatures)
library(SummarizedExperiment)
library(limma)
library(Biostrings)

#   ### 2. 🧱 PSM Object Creation & Preprocessing  
#   **Goal:** Generate a PSM (Peptide-Spectrum Match) object from `.mzID` files  
# ↓  
# - Convert `.mzID` files into `PSM` objects
# - Assess:
#   - Number of decoy hits
# - Score distributions
# - PSM rank
# - Apply filtering based on FDR or identification score

px <- PXDataset("PXD060654")
f <-pxget(px, grep("mzid", pxfiles(px)))

mzid <- mzID("D:\\MS_course\\MS_2369_mock-vs-flg-combined_031424.mzid")

mzid_obj <- mzID("D:\\MS_course\\MS_2369_mock-vs-flg-combined_031424.mzid")
id <-"D:\\MS_course\\MS_2369_mock-vs-flg-combined_031424.mzid"
psm <- PSMatch::PSM(id)

psms <- flatten(mzid)

idtbl <- as_tibble(psms)
names(idtbl)
table(idtbl$rank, useNA = "ifany") #ranks mentioned

idtbl <- idtbl |>  filter(!isDecoy) - did not find any

idtblr1 <- idtbl |>  filter(rank == 1) # filtering only rank 1 peptides

dplyr::count(idtblr1, spectrumid)

dplyr::count(idtblr1, spectrumid) |>   filter(n >1) #remove repeated spectrums to the protein

(mltm <- dplyr::count(idtblr1, spectrumid) |>     filter(n > 1) |>     pull(spectrumid)) # creating object that pulls all the spectrum ids that match more than once to the protein

(idtbl_fil <- idtblr1 |>    filter(!spectrumid %in% mltm)) #Filter out IDs matched with more than 1 protein

length(unique(idtbl_fil$pepseq))       # Peptides
length(unique(idtbl_fil$accession))

# Read FASTA with contaminants from the PXD060654
fasta <- readAAStringSet("I:/MSProt2025/idmapping_2025_05_27.fasta/file63881f23791e.fasta")

# Extract headers
prot_names <- names(fasta)

# Optionally just IDs (e.g., split at first space)
cont_prot<- sapply(strsplit(prot_names, "\\|"), `[`, 2)

present_in_psms <- cont_prot %in% idtbl_fil$accession
nc_idtbl_fil <- idtbl_fil[!(idtbl_fil$accession %in% cont_prot), ]

length(unique(nc_idtbl_fil$pepseq))       # Peptides
length(unique(nc_idtbl_fil$accession))

# Plot SEQUEST XCorr Score Distribution in a bar plot
library(ggplot2)

ggplot(nc_idtbl_fil, aes(x = `sequest:xcorr`)) +
  geom_density(fill = "lightblue") +
  labs(title = "SEQUEST XCorr Score Distribution")

# Plot PSM rank distribution in a bar plot
ggplot(idtbl, aes(x = factor(rank))) +
  geom_bar(fill = "steelblue") +
  labs(
    title = "PSM Rank Distribution",
    x = "Rank",
    y = "Count"
  ) +
  theme_minimal()

# Plot Scores by Rank using boxplot
ggplot(idtbl, aes(x = factor(rank), y = `sequest:xcorr`)) +
  geom_boxplot(fill = "lightgreen") +
  labs(
    title = "XCorr Score by PSM Rank",
    x = "Rank",
    y = "SEQUEST:xcorr"
  ) +
  theme_minimal()

high_conf <- subset(nc_idtbl_fil, `sequest:xcorr` > 2.0 & rank == 1)

length(unique(high_conf$pepseq))       # Peptides
length(unique(high_conf$accession))

# so after everything we have 21107 unique peptides and 3288 unique proteins in the mzID file


#   ### 3. 🧬 Protein & Peptide Identification  
#   **Goal:** Determine identified peptides and proteins  
# ↓  
# - Count the number of identified peptides and proteins
# - Review peptide-to-protein mapping:
#   - Razor proteins
# - Protein groups

library(dplyr)
library(readr)

# Path to files for MaxQuant output data
evidence_path <- "I:/MSProt2025/combined (2)/combined/txt/evidence.txt"
proteinGroups_path <- "I:/MSProt2025/combined (2)/combined/txt/proteinGroups.txt"
peptides_path <- "I:/MSProt2025/combined (2)/combined/txt/peptides.txt"

# Load data
evidence_data <- read.delim(evidence_path, sep = "\t", stringsAsFactors = FALSE)
proteinGroups_data <- read.delim(proteinGroups_path, sep = "\t", stringsAsFactors = FALSE)
peptides_data <- read.delim(peptides_path, sep = "\t", stringsAsFactors = FALSE)

names(proteinGroups_data)
head(proteinGroups_data)

# Find contaminants, decoy/reverse hits
dim(evidence_data) # 128 768 rows
colnames(evidence_data)
table(evidence_data$Potential.contaminant) # < table of extent 0 >, no contaminants data
table(evidence_data$Reverse) # 36 Reverse (+)

# Filter decoys
evidence_filt <- evidence_data %>%
  filter(is.na(Reverse) | Reverse != "+")

dim(evidence_filt) # 128 656 rows

#remove contaminant
head(evidence_filt$Leading.razor.protein)
evidence_filt2 <- evidence_filt[!(evidence_filt$Leading.razor.protein %in% cont_prot), ]

dim(evidence_filt2) # 127 618 rows

# Plot Score Distribution in a bar plot
ggplot(evidence_filt2, aes(x = `Score`)) +
  geom_density(fill = "lightblue") +
  labs(title = "Score Distribution")

#number of unique peptides, 23 964
length(unique(evidence_filt2$Sequence))

#number of unique proteins, 3399
length(unique(evidence_filt2$Proteins))

#number of leading razor proteins, 2559
length(unique(evidence_filt2$Leading.razor.protein))

# Table showing how many proteins each peptide maps to
peptide_to_protein_map <- evidence_filt2 %>%
  select(Modified.sequence, Proteins) %>%
  mutate(n_proteins = lengths(strsplit(Proteins, ";"))) %>%
  group_by(n_proteins) %>%
  summarise(n_peptides = n(), .groups = "drop")

print(peptide_to_protein_map)

# Plot peptide-to-protein-mapping overview as a barplot
library(ggplot2)

ggplot(peptide_to_protein_map, aes(x = factor(n_proteins), y = n_peptides)) +
  geom_col(fill = "steelblue") +
  labs(
    x = "# of proteins matched by peptide",
    y = "# of peptides",
    title = "Peptide-to-protein mapping overview"
  ) +
  theme_minimal()

# Distribution of protein group sizes
table(proteinGroups_filt2$Number.of.proteins)

#   ### 4. 🔄 QFeature Aggregation (Was not performed due to problems with aggregation and access to data in output file)  
#   **Goal:** Use quantification data if available  
# ↓  
# Use `aggregateFeatures()` from the `QFeatures` package to aggregate:
#   - PSMs ➡️ Peptides ➡️ Proteins

#   ### 5. 🧼 Normalization & Imputation (We did not normalize or impute due to access to LFQ intensity values)  
#   **Goal:** Correct for technical variation and handle missing values  
# ↓  
# - Normalize using `normalize()` or `normalize_vsn()`
# - Impute missing values with `impute()`


#  Summarizing evidence.txt to a protein-level intensity matrix
library(dplyr)
library(tidyr)

# Select needed columns
protein_quant <- evidence_filt2 %>%
  filter(!is.na(Intensity) & Intensity > 0) %>%
  select(Leading.razor.protein, Raw.file, Intensity)
head(protein_quant)

# Summarize intensities: median
protein_summary <- protein_quant %>%
  group_by(Leading.razor.protein, Raw.file) %>%
  summarize(Intensity = median(Intensity, na.rm = TRUE), .groups = "drop")
head(protein_summary)

protein_matrix <- protein_summary %>%
  pivot_wider(names_from = Raw.file, values_from = Intensity)

# Set protein IDs as row names
protein_matrix <- as.data.frame(protein_matrix)
rownames(protein_matrix) <- protein_matrix$Leading.razor.protein
protein_matrix <- protein_matrix[, -1]

head(protein_matrix)
dim(protein_matrix)

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("DEP")

library(DEP)
library(SummarizedExperiment)
library(vsn)
library(impute)
library(dplyr)
library(tibble)
library(DEP)

# Convert protein matrix to data frame with protein IDs as a column
protein_df <- protein_matrix %>%
  as.data.frame() %>%
  rownames_to_column(var = "name")

# Create a unique identifier column manually (ID should exist or be created)
protein_df$ID <- make.unique(protein_df$name)  # create 'ID' if not already there
protein_df$name_uniq <- paste0(protein_df$name, "_", protein_df$ID)

# Identify the sample (intensity) columns — exclude ID columns
id_cols <- c("name", "name_uniq", "ID")
sample_cols <- which(!(colnames(protein_df) %in% id_cols))
sample_names <- colnames(protein_df)[sample_cols]

# Create sample annotation
sample_annotation <- data.frame(
  label = sample_names,
  condition = c("flg", "flg", "flg", "flg", "mock", "mock", "mock", "mock"),  # adjust if needed
  replicate = c(1, 2, 3, 4, 1, 2, 3, 4)
)
rownames(sample_annotation) <- sample_annotation$label

# Create SummarizedExperiment object
se <- make_se(protein_df, sample_annotation, columns = sample_cols)

# Normalize
se_norm <- normalize_vsn(se)

# Impute missing values using KNN
se_imp_knn <- impute(se_norm, fun = "knn") # 563 rows with more than 50 % entries missing; mean imputation used for these rows


#   ### 6. 🧪 Protein Inference & Quantification
#   **Goal:** Summarize and quantify proteins  
# ↓  
# - Optionally annotate protein IDs using external databases (e.g., UniProt)


#   ### 7. 📊 Statistical Analysis (Perform this using Limma)  
#   **Goal:** Identify differentially abundant proteins  
# ↓  
# - Perform statistical testing (e.g., using `test_diff()`)
# - Filter by:
#   - Log2 fold-change
# - Adjusted p-value (e.g., FDR < 0.05)

library(limma)
library(SummarizedExperiment)

# Extract normalized, imputed expression matrix
exprs_mat <- assay(se_imp_knn)

# Extract sample group information
sample_info <- colData(se_imp_knn)
group <- factor(sample_info$condition)

# Create design matrix for the comparison (e.g., flg vs mock)
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
design

# Fit the model
fit <- lmFit(exprs_mat, design)

# Create contrast (e.g., flg vs mock)
contrast_matrix <- makeContrasts(FLGvsMOCK = flg - mock, levels = design)

# Apply contrast to model
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# Get all differential proteins
results <- topTable(fit2, coef = "FLGvsMOCK", number = Inf, adjust = "fdr")

# Filter significant proteins (FDR < 0.05 and abs(log2FC) > 1)
sig_proteins <- results %>%
  dplyr::filter(adj.P.Val < 0.05 & abs(logFC) > 1)

head(sig_proteins)

#   ### 8. 📈 Visualization & Export (Optional)  
#   **Goal:** Visualize data for interpretation  
# ↓  
# - Create visual summaries:
#   - PCA plots
# - Volcano plots
# - Heatmaps  
# - Use packages such as `ggplot2`, `DEP`, or `limma`

library(ggplot2)
library(ggrepel)

# Run PCA on transposed expression matrix
pca <- prcomp(t(exprs_mat), scale. = TRUE)

# Calculate % variance explained
pca_var <- pca$sdev^2
pca_var_percent <- round(100 * pca_var / sum(pca_var), 1)

# Create PCA dataframe
pca_df <- data.frame(
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  condition = colData(se_imp_knn)$condition,
  replicate = colData(se_imp_knn)$replicate
)

# Plot with % variance on axes
ggplot(pca_df, aes(x = PC1, y = PC2, color = condition, label = replicate)) +
  geom_point(size = 4) +
  geom_text_repel(size = 3) +
  theme_minimal() +
  labs(
    title = "PCA of Protein Intensities",
    x = paste0("PC1 (", pca_var_percent[1], "%)"),
    y = paste0("PC2 (", pca_var_percent[2], "%)")
  )

library(pheatmap)

top_proteins <- results %>%
  dplyr::filter(adj.P.Val < 0.05 & abs(logFC) > 1) %>%
  dplyr::arrange(adj.P.Val) %>%
  head(30)  # adjust this number as needed

heatmap_mat <- exprs_mat[rownames(exprs_mat) %in% rownames(top_proteins), ]
heatmap_mat_scaled <- t(scale(t(heatmap_mat)))  # z-score scaling

sample_anno <- data.frame(condition = colData(se_imp_knn)$condition)
rownames(sample_anno) <- colnames(heatmap_mat_scaled)

pheatmap(heatmap_mat_scaled,
         annotation_col = sample_anno,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,  # show protein IDs
         fontsize_row = 8,      # optional: adjust font size
         main = "Top Differentially Abundant Proteins")

# # Save all and significant results
# write.csv(results, "limma_all_proteins.csv")
# write.csv(sig_proteins, "limma_significant_proteins.csv")

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
