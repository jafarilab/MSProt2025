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

# install.packages("BiocManager")
#install.packages("tidyverse")
#install.packages("factoextra")
#install.packages("pheatmap")
#install.packages("remotes")
# 
# library(BiocManager)
# 
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
############ It is very good that you left package installations here
########### However some packages are Bioconductor packages and some are from CRAN (e.g., tidyverse)
########### It is important to make this distinction here
#BiocManager::install("msdata")
#BiocManager::install("mzR")
#BiocManager::install("rhdf5")
#BiocManager::install("rpx")
#BiocManager::install("MsCoreUtils")
#BiocManager::install("QFeatures")
#BiocManager::install("Spectra")
#BiocManager::install("ProtGenerics")
#BiocManager::install("PSMatch")
#BiocManager::install("limma")
#BiocManager::install("MSnID")
#BiocManager::install("RforMassSpectrometry/SpectraVis")
#BiocManager::install("Biostrings")

#   ### 2. 🧱 PSM Object Creation & Preprocessing  
#   **Goal:** Generate a PSM (Peptide-Spectrum Match) object from `.mzID` files  
# ↓  
# - Convert `.mzID` files into `PSM` objects
# - Assess:
#   - Number of decoy hits
# - Score distributions
# - PSM rank
# - Apply filtering based on FDR or identification score

library(dplyr)
library(ggplot2)
library(PSMatch)
library(Biostrings)
library(tibble)
library(tidyr)
library(mzR)
library(Rcpp)
library(stringr)
library(rpx)
library(SummarizedExperiment)
library(vsn)
library(MSnbase)
library(impute)
library(purrr)

###########For reviewers (re-running the analysis):
  #### Open this project folder as the working directory.

  #### Place any external inputs in this working directory (e.g., files from Zenodo).

# Access dataset using ProteomeXChange accession code
px <- PXDataset("PXD060654")
f <-pxfiles(px)


path_to_file <- pxget(px , 'MS_2369_mock-vs-flg-combined_031424.mzid')

# Extract peptide-spectrum match data into object
psm <- PSMatch::PSM(path_to_file, parser='mzID')

idtbl <- as_tibble(psm)
names(idtbl)

# Distribution of ranks and include any missing values (NA)
table(idtbl$rank, useNA = "ifany")

# Filter out decoy peptides (from "isDecoy" column)
idtbl <- idtbl |>  filter(!isDecoy) # No decoys found

# Plot PSM rank distribution in a bar plot
ggplot(idtbl, aes(x = factor(rank))) +
  geom_bar(fill = "steelblue") +
  labs(
    title = "PSM Rank Distribution",
    x = "Rank",
    y = "Count"
  ) +
  theme_minimal()

# Include only rows with rank 1
idtbl_filt <- idtbl |> filter(rank == 1)

length(unique(idtbl_filt$pepseq))
length(unique(idtbl_filt$accession))

# Visualize how many times a specific spectrumid is matched to a peptide
spectrum_counts <- dplyr::count(idtbl_filt, spectrumid)
ggplot(spectrum_counts, aes(x = n)) +
  geom_histogram(binwidth = 1, fill = "lightblue", color = "black", alpha = 0.7) +
  labs(title = "Histogram of Spectrum ID Occurrences", 
       x = "Number of Occurrences", 
       y = "Frequency") +
  theme_minimal()

# Count occurrences of each spectrumid in idtbl_filt
mltm <- dplyr::count(idtbl_filt, spectrumid) |> 
  filter(n > 1) |> 
  pull(spectrumid)  # Spectrum IDs that match more than one peptide

# Filter idtblr1 to remove spectra matched to more than one peptide
idtbl_filt2 <- idtbl_filt |> 
  filter(!spectrumid %in% mltm)


fasta <- readAAStringSet("./Team2_uniprotkb_proteome_UP000006548_AND_revi_2025_05_22.fasta")

# Extract headers
prot_names <- names(fasta)

# Extract only IDs (e.g., split at first space)
cont_prot<- sapply(strsplit(prot_names, "\\|"), `[`, 2)

# Filter out rows with accession codes matching contaminants fasta
nc_idtbl_fil2 <- idtbl_filt2[!(idtbl_filt2$accession %in% cont_prot), ]

# Count unique peptides and proteins
length(unique(nc_idtbl_fil2$pepseq))
length(unique(nc_idtbl_fil2$accession))

#   ### 3. 🧬 Protein & Peptide Identification  
#   **Goal:** Determine identified peptides and proteins  
# ↓  
# - Count the number of identified peptides and proteins
# - Review peptide-to-protein mapping:
#   - Razor proteins
# - Protein groups


# Path to files for MaxQuant output data
##### If I want to execute your MaxQuant seaarch I need parameter file if settings were
##### not by default, also FASTA file you used for search. All in all, please
##### provide all files that are need to execute MQ in the same way you did it
##### I don't have these paths on my local system
evidence_path <- "./Team2_evidence.txt"
proteinGroups_path <- "./Team2_proteinGroups.txt"
msms_path <- "./Team2_msms.txt"

# Load data
evidence_data <- read.delim(evidence_path, sep = "\t", stringsAsFactors = FALSE)
proteinGroups_data <- read.delim(proteinGroups_path, sep = "\t", stringsAsFactors = FALSE)
msms_data <- read.delim(msms_path, sep = "\t", stringsAsFactors = FALSE)

colnames(evidence_data)
colnames(proteinGroups_data)
colnames(msms_data)

# Find contaminants, decoy/reverse hits
dim(evidence_data) # 128 768 rows
colnames(evidence_data)
table(evidence_data$Reverse) # 112 decoys/reverse (+)

# Filter decoys
evidence_filt <- evidence_data %>%
  filter(is.na(Reverse) | Reverse != "+")

dim(evidence_filt) # 128 656 rows

#remove contaminant
head(evidence_filt$Leading.razor.protein)
nc_evidence_filt <- evidence_filt[!(evidence_filt$Leading.razor.protein %in% cont_prot), ]

dim(nc_evidence_filt) # 127 618 rows

# Count unique peptides and proteins
length(unique(nc_evidence_filt$Sequence)) # 23 964 unique peptides
length(unique(nc_evidence_filt$Proteins)) # 3399 unique proteins
length(unique(nc_evidence_filt$Leading.razor.protein)) # 2559 leading razor proteins
length(unique(nc_evidence_filt$Protein.group.IDs)) # 3048 protein groups


# Table showing how many proteins each peptide maps to
peptide_to_protein_map <- nc_evidence_filt %>%
  select(Modified.sequence, Proteins) %>%
  mutate(n_proteins = lengths(strsplit(Proteins, ";"))) %>%
  group_by(n_proteins) %>%
  summarise(n_peptides = n(), .groups = "drop")

print(peptide_to_protein_map)

# Plot peptide-to-protein-mapping overview as a barplot
ggplot(peptide_to_protein_map, aes(x = factor(n_proteins), y = n_peptides)) +
  geom_col(fill = "steelblue") +
  labs(
    x = "# of proteins matched by peptide",
    y = "# of peptides",
    title = "Peptide-to-protein mapping overview"
  ) +
  theme_minimal()

#   ### 5. 🧼 Normalization & Imputation   
#   **Goal:** Correct for technical variation and handle missing values  
# ↓  
# - Normalize using `normalize()` or `normalize_vsn()`
# - Impute missing values with `impute()`



# Data preparation
# Filter intensity values and select relevant columns for protein quantification
protein_quant <- nc_evidence_filt %>%
  filter(!is.na(Intensity) & Intensity > 0) %>%
  select(Leading.razor.protein, Raw.file, Intensity)

# Aggregate multiple intensities per protein/sample to a single median value
protein_summary <- protein_quant %>%
  group_by(Leading.razor.protein, Raw.file) %>%
  summarize(Intensity = median(Intensity, na.rm = TRUE), .groups = "drop")

# Transform data into a wide format matrix
protein_matrix <- protein_summary %>%
  pivot_wider(names_from = Raw.file, values_from = Intensity) %>%
  as.data.frame()

# Set protein IDs as row names
rownames(protein_matrix) <- protein_matrix$Leading.razor.protein

# Move row names into a column
protein_df <- protein_matrix %>%
  rownames_to_column(var = "name")

# Generate unique names to avoid duplication
protein_df$ID <- make.unique(protein_df$name)
protein_df$name_uniq <- paste0(protein_df$name, "_", protein_df$ID)

# Isolate sample intensity columns
id_cols <- c("name", "name_uniq", "ID", "Leading.razor.protein")
sample_cols <- setdiff(colnames(protein_df), id_cols)
sample_names <- sample_cols

# Setup of metadata
# Define sample conditions and replicates
sample_annotation <- data.frame(
  label = sample_names,
  condition = c("flg", "flg", "flg", "flg", "mock", "mock", "mock", "mock"),
  replicate = c(1, 2, 3, 4, 1, 2, 3, 4)
)

# Create SummarizedExperiment
# Extract the expression matrix and assign protein names as row names
exprs_mat <- as.matrix(protein_df[, sample_cols])
rownames(exprs_mat) <- protein_df$name

# Creates row-level metadata (protein IDs)
rowData <- data.frame(Protein_ID = rownames(exprs_mat))

# Combine expression data, sample metadata, and protein metadata
se <- SummarizedExperiment(
  assays = list(counts = exprs_mat),
  rowData = rowData,
  colData = sample_annotation
)

# Normalization and imputation
# Apply VSN to stabilize variance across expression levels
vsn_fit <- vsn2(assay(se, "counts"))
exprs_norm <- predict(vsn_fit, assay(se, "counts"))

# Add normalized data as an additional assay
assay(se, "vsn") <- exprs_norm
se_norm <- se

# Wrap normalized data in MSnSet for compatibility with imputation
msnset <- MSnSet(
  exprs = exprs_norm,
  phenoData = AnnotatedDataFrame(sample_annotation)
)

# Impute missing values using KNN
msnset_imputed <- impute(msnset, method = "knn")

exprs_imputed <- exprs(msnset_imputed)

assay(se, "vsn_imputed") <- exprs_imputed
se_imp_knn <- se

# Now this works
exprs_mat <- assay(se_imp_knn, "vsn_imputed")


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

library(dplyr)
library(limma)

exprs_mat <- assay(se_imp_knn, "vsn_imputed")

# Extract sample metadata
sample_info <- colData(se_imp_knn)

# Create design matrix for the comparison (flg vs mock)
condition <- factor(sample_info$condition)
design <- model.matrix(~ 0 + condition)
colnames(design) <- levels(condition)

# Fit the model
fit <- lmFit(exprs_mat, design)

# Create contrast (e.g., flg vs mock)
contrast_matrix <- makeContrasts(FLGvsMOCK = flg - mock, levels = design)

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

library(dplyr)
library(ggplot2)
library(ggrepel)
library(pheatmap)

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

# Visualize results in heatmap
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

# Only 11 significantly different protein when performing the statistical analysis ourselves
# Filter and visualize normalized and imputed LFQ intensisities in proteinGroups.txt output file

names(proteinGroups_data)

# Remove decoys
maxquant_fil <- proteinGroups_data %>% filter(is.na(Reverse) | Reverse == "")

# Remove rows where protein IDs match contaminant
nc_maxquant_fil <- maxquant_fil[!(maxquant_fil$Majority.protein.IDs %in% cont_prot), ]

# Calculate number of unique peptides
sum(nc_maxquant_fil$Unique.peptides)

# Calculate number of razor.uniqu.peptides
sum(nc_maxquant_fil$Razor...unique.peptides)

# Generate table of number of protein grouped together due to similar peptide sequence
table(nc_maxquant_fil$Number.of.proteins)

# Plot score distribution in density plot
ggplot(nc_maxquant_fil, aes(x = `Score`)) +
  geom_density(fill = "lightblue") +
  labs(title = " Score Distribution")+
  geom_vline(xintercept = 10, 
             color = "red", linetype = "dashed", size = 0.5)


# Filter based on score >10 and Q value < 0.01
high_conf_maxquant <- subset(nc_maxquant_fil, ((`Q.value` < 0.01) & (`Score`>10)))

# Save original filtered dataset for comparison
high_conf_maxquant1 <- nc_maxquant_fil

# Extract LFQ intinsity data
lfq_cols <- grep("^LFQ\\.intensity\\.", names(high_conf_maxquant), value = TRUE)
exprs <- high_conf_maxquant[, lfq_cols]
rownames(exprs) <- high_conf_maxquant$Majority.protein.IDs

# Convert to numeric and log2 transform
enrichment <- exprs[, "LFQ.intensity.flg"] / exprs[, "LFQ.intensity.mock"]
log2FC <- log2(as.matrix(enrichment))

# Create a dataframe containing protein names, gene names and corresponding log2 FC
deg_results <- data.frame(
  Protein_names = high_conf_maxquant$Protein.names,
  Gene_names = high_conf_maxquant$Gene.names,
  log2FC = log2FC
)
# Sort by log2FC in descending order
order_results <- deg_results[order(deg_results$log2FC, decreasing = TRUE), ]
head(order_results)

# Without filtering
lfq_cols1 <- grep("^LFQ\\.intensity\\.", names(high_conf_maxquant1), value = TRUE)
exprs1 <- high_conf_maxquant1[, lfq_cols1]
rownames(exprs1) <- high_conf_maxquant1$Protein.IDs

enrichment1 <- exprs1[, "LFQ.intensity.flg"] / exprs1[, "LFQ.intensity.mock"]

log2FC1 <-   log2(as.matrix(enrichment1))

deg_results1 <- data.frame(
  Protein_names = high_conf_maxquant1$Protein.names,
  Gene_names = high_conf_maxquant1$Gene.names,
  log2FC = log2FC1
)

order_results1 <- deg_results1[order(deg_results1$log2FC, decreasing = TRUE), ]
head(order_results1)

# Extract top 30 proteins from filtered dataset (score and Q value)
top30 <- deg_results %>%  slice_max(log2FC, n = 30)

# Plot as dot plot
ggplot(top30, aes(x = log2FC, y = reorder(Gene_names, log2FC))) +
  geom_point(size = 3, color = "darkblue") +
  theme_minimal() +
  labs(title = "Top 30 Proteins by Log2 Fold Change",
       x = "Log2 Fold Change (FLG vs MOCK)",
       y = "Gene Name") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray")


# Extract top 30 proteins from unfiltered dataset
top30_1 <- deg_results1 %>%  slice_max(log2FC, n = 30)

# Plot as dot plot
ggplot(top30_1, aes(x = log2FC, y = reorder(Gene_names, log2FC))) +
  geom_point(size = 3, color = "darkblue") +
  theme_minimal() +
  labs(title = "Top 30 Proteins by Log2 Fold Change",
       x = "Log2 Fold Change (FLG vs MOCK)",
       y = "Gene Name") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray")

