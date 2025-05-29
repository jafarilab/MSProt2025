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

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

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
BiocManager::install("Biostrings")

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

# Access dataset using ProteomeXChange accession code
px <- PXDataset("PXD060654")

# Retrieve list of files in px
f <-pxfiles(px)

# Extract peptide-spectrum match data into object
psm <- PSMatch::PSM("C:/Users/janni/Downloads/39fc2ba81653_MS_2369_mock-vs-flg-combined_031424.mzid")

# Convert the PSM to a tibble (modern data frame)
psms <- flatten(mzid)
idtbl <- as_tibble(psms)
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

# Read FASTA with contaminants from the PXD060654
fasta <- readAAStringSet("C:/Users/janni/OneDrive/Työpöytä/HU - Proteomics 2025/MSProt2025_New_fork/MSProt2025_Team2/team2/file63881f23791e.fasta")

# Extract headers
prot_names <- names(fasta)

# Extract only IDs (e.g., split at first space)
cont_prot<- sapply(strsplit(prot_names, "\\|"), `[`, 2)

# Filter out rows with accession codes matching contaminants fasta
nc_idtbl_fil2 <- idtbl_filt2[!(idtbl_filt2$accession %in% cont_prot), ]

# Count unique peptides and proteins
length(unique(nc_idtbl_fil$pepseq)) # 21107 unique peptides
length(unique(nc_idtbl_fil$accession)) # 3288 unique proteins

#   ### 3. 🧬 Protein & Peptide Identification  
#   **Goal:** Determine identified peptides and proteins  
# ↓  
# - Count the number of identified peptides and proteins
# - Review peptide-to-protein mapping:
#   - Razor proteins
# - Protein groups

library(dplyr)
library(ggplot2)
library(stringr)
library(Biostrings) 

# Path to files for MaxQuant output data
evidence_path <- "C:/Users/janni/Downloads/combined (2)/combined/txt/evidence.txt"
proteinGroups_path <- "C:/Users/janni/Downloads/combined (2)/combined/txt/proteinGroups.txt"
msms_path <- "C:/Users/janni/Downloads/combined (2)/combined/txt/msms.txt"

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

library(dplyr)
library(tidyr)
library(SummarizedExperiment)
library(vsn)
library(MSnbase)
library(impute)

# Filter and select relevant columns from evidence data
protein_quant <- nc_evidence_filt %>%
  filter(!is.na(Intensity) & Intensity > 0) %>%
  select(Leading.razor.protein, Raw.file, Intensity)

# Summarize intensities by median per protein/sample
protein_summary <- protein_quant %>%
  group_by(Leading.razor.protein, Raw.file) %>%
  summarize(Intensity = median(Intensity, na.rm = TRUE), .groups = "drop")

# Pivot to wide format to create a protein-level intensity matrix
protein_matrix <- protein_summary %>%
  pivot_wider(names_from = Raw.file, values_from = Intensity) %>%
  as.data.frame()

# Set protein IDs as rownames
rownames(protein_matrix) <- protein_matrix$Leading.razor.protein

# Convert matrix to data frame with unique identifiers
protein_df <- protein_matrix %>%
  rownames_to_column(var = "name")

# Create unique identifiers if not present
protein_df$ID <- make.unique(protein_df$name)
protein_df$name_uniq <- paste0(protein_df$name, "_", protein_df$ID)

# Identify sample columns (exclude metadata/ID columns)
id_cols <- c("name", "name_uniq", "ID", "Leading.razor.protein")
sample_cols <- setdiff(colnames(protein_df), id_cols)
sample_names <- sample_cols

# Create sample annotation (adjust if needed)
sample_annotation <- data.frame(
  label = sample_names,
  condition = c("flg", "flg", "flg", "flg", "mock", "mock", "mock", "mock"),
  replicate = c(1, 2, 3, 4, 1, 2, 3, 4)
)

# Prepare the expression matrix
exprs_mat <- as.matrix(protein_df[, sample_cols])
rownames(exprs_mat) <- protein_df$name

# Create metadata for proteins (rowData)
rowData <- data.frame(Protein_ID = rownames(exprs_mat))

# Create the SummarizedExperiment object
se <- SummarizedExperiment(
  assays = list(counts = exprs_mat),
  rowData = rowData,
  colData = sample_annotation
)

# Normalize using Variance Stabilizing Normalization (VSN)
vsn_fit <- vsn2(assay(se, "counts"))
exprs_norm <- predict(vsn_fit, assay(se, "counts"))

# Add normalized data as an additional assay
assay(se, "vsn") <- exprs_norm
se_norm <- se

# Convert matrix to MSnSet
msnset <- MSnSet(
  exprs = exprs_norm,
  phenoData = AnnotatedDataFrame(sample_annotation)
)

# Impute missing values using KNN
msnset_imputed <- impute(msnset, method = "knn")


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

#   ### 9. 📖 Comparison with Published Results  
#   **Goal:** Benchmark your findings  
# ↓  
# Create a comparison table including:
#   - Number of spectra identified
# - Number of peptides and proteins
# - Key figures or results (if available)
# - Comments on reproducibility

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
