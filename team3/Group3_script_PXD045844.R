#install.packages("BiocManager")
library(BiocManager)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")


### AIM 1: Dataset Acquisition

install.packages("tidyverse")
library("rpx")

px <- PXDataset("PXD045844")

pxfiles(px, as.vector = TRUE) |> print() #TG: Dataset contains mzid, mzML, raw and msf files.

# IMPORTANT: The three-digit numbers in the raw file names represent the sample IDs (gender information, OBTAINED FROM THE AUTHORS): 
# 097 – F, 098 – F, 099 – F, 100 – F
# 094 – M, 101 – M, 102 – M, 103 – M

pxtax(px) # taxonomic name #TG: Note this study was done in mice
pxref(px) # bibliography ref
pxinstruments(px) # instrument
pxSubmissionDate(px) # submission date
pxPublicationDate(px) # publication date
pxptms(px) # PTMs


### Getting the data (PXD041638) 

pxget(px, "PatchClamp_Gender_Comparison_LFQ.mzML") # This mzML file contains only MS2, check raw files for more info! #TG: It gives you your own directory
#pxget(px, "PatchClamp_Gender_Comparison_LFQ.mzid.gz") #DONE LATER


### Exploring mzML files

library(Spectra)
library(tidyverse)
library(SpectraVis)

f1 = "/Users/tamgon/Library/Caches/org.R-project.R/R/rpx/63002267c635_PatchClamp_Gender_Comparison_LFQ.mzML" #TG directory
(sp1 <- Spectra(f1))
#print(object.size(sp1), units = "Mb")
length(sp1)
spectraVariables(sp1)
pd1 <- peaksData(sp1)
table(msLevel(sp1)) # How many MS levels and how many scans of each level

head(pd1[[1345]]) # Just to test 
head(sp1$msLevel, 300) # Just to test
msLevel(sp1)[[1234]] # Just to test
plot(pd1[[1234]], type = "h") # Just to test
head(precursorMz(sp1), 500) # Just to test

sp1_2 <- filterMsLevel(sp1, 2L) # TG: No need to filter MS1, since we don't have any in the file mzML
max(sp1$basePeakIntensity)

plot(rtime(sp1), tic(sp1), type = "l")
plot(sp1$rtime, sp1$totIonCurrent, type = "l")

spectraData(sp1) |> 
  as_tibble() |>
  filter(msLevel == 2) |>
  ggplot(aes(x = rtime,
             y = totIonCurrent)) +
  geom_line()

precScanNum(sp1_2)
table(msLevel(sp1), centroided(sp1))
anyDuplicated(precursorMz(sp1_2)) 
browseSpectra(sp1)



### AIM 2: PSM Object Creation & Preprocessing
### Goal: Generate a PSM (Peptide-Spectrum Match) object from .mzID files

### Convert .mzID files into PSM objects

library(PSMatch) 
#PSM when we can match a spectrum to a peptide using a Search engine software (third-party software) such as Mascot, PD, Sage, MaxQuant, SAGE.
idf <- pxget(PXDataset("PXD045844"), "PatchClamp_Gender_Comparison_LFQ.mzid.gz")
id <- PSM(idf, parser = "mzID") # mzID contains info for peptide and protein identification

### More info about the PSM object -Part 1

dim(id) #dimensions of the PSM object id — number of rows (PSMs) and columns (attributes such as spectrum ID, sequence, file, etc.).
names(id) #attributes

head(id$spectrumid) # first few spectrum IDs — unique identifiers for the spectra matched to peptides #TG: Make sure your file names match here
head(id$sequence) # peptide sequences that were matched to spectra
head(id$accession) # protein database accession numbers (e.g., UniProt IDs) from which the peptides were derived.


### More info -Part 2

length(unique(id$spectrumid))     # Counts the number of unique MS/MS spectra in your dataset --> 7048
length(unique(id$sequence))       # Counts the number of unique full protein sequences --> 588
length(unique(id$pepseq))         # Unique peptide sequences --> 1756
length(unique(id$accession))      # Counts the number of unique proteins (by accession ID) identified --> 588

library(tibble)
library(DT)
as.data.frame(id) |> # TABLE WITH ALL DATA!!
  as_tibble() |>
  DT::datatable()

### Assess number of decoy hits

id$isdecoy # NULL --> TG: When running str(id), variable "isdecoy" is present, but it has not been exported in the id file (?)
#table(grepl("rev|decoy", id@listData$accession, ignore.case = TRUE)) #TG: The accession column only contains real protein IDs — no typical decoy prefixes like "REV_" or "DECOY_".


### Assess score distributions

#TG: we cannot calculate decoy score (scores assigned to PSMs matched against decoy sequences (reversed proteins) because isdecoy is missing). Use another score instead.
as.data.frame(id) |>
  as_tibble() |>
  ggplot(aes(x = sequest.xcorr)) +  # TG: XCorr score measures how well a given peptide’s theoretical spectrum matches the observed MS/MS spectrum.
  geom_density(fill = "skyblue", alpha = 0.5) +
  labs(title = "Score Distribution (sequest:xcorr)",
       x = "XCorr Score",
       y = "Density") #TG: XCorr scores are skewed toward lower values with a tail at higher scores. A few matches are very strong and get high scores, but these are rarer.

table(table(id$spectrumid)) # How many scans match to 1, 2, 3... proteins --> TG: one match seems more common than multiple matches


### Assess PSM rank

library(dplyr)
library(DT)

#TG: Concatenated functions, gives a summary of specific functions - PSM objects filtered in every step
filter_with_log <- function(df) {
  initial <- nrow(df)
  message("Starting with ", initial, " PSMs")
  
  df_rank1 <- df %>% filter(rank == 1)
  message("Removed ", initial - nrow(df_rank1), " PSMs with rank > 1")  # keeps the ones with the best scores 
  
  shared <- df_rank1 %>% # identifying shared peptides
    count(spectrumid) %>%
    filter(n > 1) %>% # PSMs where one scan matches multiple proteins
    pull(spectrumid)
  
  df_final <- df_rank1 %>% filter(!spectrumid %in% shared) # filtering them out
  message("Removed ", nrow(df_rank1) - nrow(df_final), " shared peptides")
  
  message("Final count: ", nrow(df_final), " PSMs")
  return(df_final)
}

df_filtered <- filter_with_log(df) #history of specific functions


### Filtering based on FDR or identification score

# We have no decoys, so we can’t compute empirical FDR, but we can use score thresholds instead
df_filtered <- df %>%
  filter(rank == 1, sequest.xcorr >= 2.0)

cat("After rank + score filtering:", nrow(df_filtered), "PSMs retained\n")

#TG: IMPORTANT. Note that we're applying different filters here:
# df_filtered after rank + shared peptide filtering → 5,513 PSMs
# df_filtered after rank + XCorr score filtering → 8,193 PSMs




### AIM 3: Protein & Peptide Identification

psms <- PSM(idf, parser = "mzID") #renamed here

#Tyko- number of scans and PSM's
library(mzID)
id <- mzID(idf)
id
psm <- flatten(id) #this is also PSMs but a bit different (?)

n_peptides <- length(unique(psms$pepseq))
cat("Identified peptides:", n_peptides, "\n")
n_proteins <- length(unique(psms$accession))
cat("Identified proteins:", n_proteins, "\n")

#peptide to protein mapping
peptide_to_protein <- psm %>%
  select(pepseq, accession) %>%
  distinct() %>%
  group_by(pepseq) %>%
  summarise(proteins = paste(unique(accession), collapse = ";"))
view(peptide_to_protein)

#Razor proteins
protein_counts <- psm %>%
  count(accession, name = "peptide_count") %>%
  arrange(desc(peptide_count))
view(protein_counts)

razor_map <- psm %>%
  count(pepseq, accession, name = "n") %>%
  group_by(pepseq) %>%
  slice_max(order_by = n, n = 1, with_ties = FALSE) %>%  # ensure only one razor protein per peptide
  rename(razor_protein = accession)
#Protein groups

protein_to_peptides <- psm %>%
  select(accession, pepseq) %>%
  distinct() %>%
  group_by(accession) %>%
  summarise(peptides = paste(sort(unique(pepseq)), collapse = ";"), .groups = "drop")

# Group proteins with identical peptide sets
protein_groups <- protein_to_peptides %>%
  group_by(peptides) %>%
  summarise(
    protein_group = paste(accession, collapse = ";"),
    .groups = "drop"
  )

# Count number of peptides in each group
protein_groups <- protein_groups %>%
  mutate(n_peptides = lengths(strsplit(peptides, ";"))) %>%
  arrange(desc(n_peptides))  # sort by peptide count
view (protein_groups)




### AIM 4:  QFeature Aggregation (Optional). Done using a evidence.txt file, obtained from MaxQuant 
# MaxQuant: 1) Processes all the raw files in parallel, 2) Matches features across files (alignment, normalization, etc.), 3) Outputs combined quantification and identification table into a single evidence.txt file)
# Evidence.txt: Each row = one PSM (peptide-spectrum match), that has a column Raw.file to tell which sample it came from (long format)


# Libraries
library(QFeatures)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
library(SummarizedExperiment)
library(tidyverse)

#Read evidence.txt file
evidence <- read.delim("evidence.txt", stringsAsFactors = FALSE)


# To generate a QFeatures object: transform file into one row per peptide and multiple intensity columns like this: Peptide	Intensity_M1	Intensity_F1 ... (wide format)
evidence %>%
  group_by(Sequence, Raw.file) %>%
  summarise(n = n()) %>%
  filter(n > 1) #TG: It seems we have duplicates of intensities per one raw.file sample

# Take the mean intensity, for example
evidence_summary <- evidence %>%
  filter(!is.na(Intensity), Intensity > 0) %>%
  group_by(Sequence, Raw.file) %>%
  summarise(Intensity = mean(Intensity), .groups = "drop") #TG: Long data frame with one Intensity per Sequence per Raw.file, keep just those columns

# Convert to a wide format
evidence_wide <- evidence_summary %>%
  pivot_wider(names_from = Raw.file, values_from = Intensity) #TG: Convert into quantitative matrix, peptides and intensities per sample


# Build the colData dataframe with sample info
samples <- c("20230901_PatchClamp_094", "20230901_PatchClamp_097", "20230901_PatchClamp_098",
             "20230901_PatchClamp_099", "20230901_PatchClamp_100", "20230901_PatchClamp_101",
             "20230901_PatchClamp_102", "20230901_PatchClamp_103")

gender <- c("M", "F", "F", "F", "F", "M", "M", "M")

col_data <- data.frame(
  Sample = samples,
  Group = gender,
  row.names = samples
)

# Make sure columns are in the same order as col_data
assay_matrix <- evidence_wide %>%
  column_to_rownames("Sequence") %>%
  select(all_of(samples)) %>%  # select samples in right order
  as.matrix()

# Prepare rowData (feature metadata) with peptide → protein mapping
row_metadata <- evidence %>%
  distinct(Sequence, Proteins) %>%
  filter(Sequence %in% rownames(assay_matrix)) %>%
  column_to_rownames("Sequence")

# Create SummarizedExperiment and QFeatures objects
se <- SummarizedExperiment(
  assays = list(Intensity = assay_matrix),
  rowData = row_metadata[rownames(assay_matrix), , drop = FALSE],
  colData = col_data
)

qf <- QFeatures(list(psms = se)) #Now qf has: psms assay with intensities, rowData with feature metadata, colData with sample groups (gender)

# Check created object
qf
assay(qf[[1]])        # Shows the quantitative data matrix (features x samples)
rowData(qf[[1]])      # Metadata (protein names) for each peptide-spectrum match (PSM)


##Razor proteins
length(unique(evidence$`Leading.razor.protein`))
length(unique(evidence$`Leading.proteins`))
length(unique(evidence$`Protein.group.IDs`))
### AIM 5:Normalization & Imputation
library(vsn)

# Normalize intensities (e.g., Variance Stabilizing Normalization)
qf <- normalize(qf, i = "psms", method = "vsn") #TG: It corrects for systematic technical variation across samples.
meanSdPlot(assay(qf[["psms"]])) # This will plot: The mean of each feature vs. its standard deviation. Ideally showing a flat line, meaning the variance is no longer dependent on the mean — a sign of successful normalization.
# Our plot: Good VSN performance for low to mid intensities

# Impute missing values
# Options: "min", "knn", "random", "QRILC" etc. (see ?impute)
qf <- impute(qf, i = "psms", method = "min") #TG: This replaces missing values with small constants (often the lowest observed intensity).

# Check normalized and imputed assay matrix
assay(qf[["psms"]])


### AIM 6: Protein Inference & Quantification (Optional)

#rowData(qf[["psms"]]) #TG: For some reason, the column sequence was not named as such
# Add Sequence as a new column in rowData
rowData(qf[["psms"]])$Sequence <- rownames(qf[["psms"]]) #TG: Now we have 2 columns, sequence and proteins
rowData(qf[["psms"]]) 

# Use Leading.razor.protein column as ID in Proteins column
rowData(qf[["peptides"]])$Proteins <- rowData(qf[["peptides"]])$Leading.razor.protein

# See what's in Proteins column
head(rowData(qf[["peptides"]])$Proteins)

# If they are like "P12345;P67890", keep only the first
rowData(qf[["peptides"]])$Proteins <- sapply(
  strsplit(as.character(rowData(qf[["peptides"]])$Proteins), ";"),
  `[`, 1
)

qf[["proteins"]] <- NULL #If "proteins" assay already exists, remove it first
qf <- aggregateFeatures(qf, i = "peptides", fcol = "Proteins", name = "proteins", fun = colMeans, na.rm = TRUE) #aggregate

names(qf)
head(assay(qf[["proteins"]])) #inspect protein-level data

sum(is.na(assay(qf[["proteins"]]))) #1684
sum(is.nan(assay(qf[["proteins"]]))) #1684


# Impute the assay alone (a SummarizedExperiment)
imputed_proteins <- impute(qf[["proteins"]], method = "min")

# Replace the original assay with the imputed one inside the QFeatures object
qf <- replaceAssay(qf, i = "proteins", imputed_proteins)

sum(is.na(assay(qf[["proteins"]]))) #Now is 0
sum(is.nan(assay(qf[["proteins"]]))) #Now is 0

sum(is.na(assay(qf[["psms"]]))) #3476 -->  I still see NAs, even after imputation. After aggregation, the new "proteins" assay is created from the "peptides" assay, which contains NAs
qf <- impute(qf, i = "psms", method = "min", name = "psms", replace = TRUE)

#THIS NEEDS TO BE FIXED, SO FAR I'M HAVING THIS NEW ISSUE:
#If you did imputation before normalization or aggregation (e.g., peptides → proteins), it’s possible those transformations reintroduced NAs (especially with na.rm = TRUE used in aggregateFeatures()).


### I TRIED TO GENERATE SOME PLOTS ANYWAY, JUST TO CHECK WHETHER IT MAKES ANY SENSE AT ALL....

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DEP")

library(DEP)

# Convert QFeatures to SummarizedExperiment
se <- qf[["proteins"]]

# Run PCA
# Get the assay matrix
mat <- assay(se)

# Remove proteins with zero variance
mat_filtered <- mat[apply(mat, 1, function(x) sd(x) != 0), ]

# Transpose for PCA
mat_t <- t(mat_filtered)

# Run PCA
pca <- prcomp(mat_t, scale. = TRUE)

# Prepare PCA data frame
pca_df <- data.frame(pca$x, Group = colData(se)$Group)

# Plot using ggplot2
library(ggplot2)
ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 4) +
  theme_minimal() +
  labs(title = "PCA of Protein Intensities", x = "PC1", y = "PC2")

# PLOT PROTEINS IN COMMON
mat <- assay(qf[["proteins"]])
groups <- colData(qf[["proteins"]])$Group

# Logical matrix: TRUE if protein is detected
detected <- !is.na(mat)

# Get sample indices per group
group_M <- which(groups == "M")
group_F <- which(groups == "F")

# Check detection per protein across replicates
present_in_M <- rowSums(detected[, group_M]) > 0
present_in_F <- rowSums(detected[, group_F]) > 0

only_M     <- sum(present_in_M & !present_in_F)
only_F     <- sum(!present_in_M & present_in_F)
shared     <- sum(present_in_M & present_in_F)
total_proteins <- nrow(mat)

cat("Proteins only in M:", only_M, "\n") #0
cat("Proteins only in F:", only_F, "\n") #0
cat("Proteins shared:", shared, "\n") #394
cat("Total proteins:", total_proteins, "\n")

library(VennDiagram)
venn.plot <- draw.pairwise.venn(
  area1 = sum(present_in_M),
  area2 = sum(present_in_F),
  cross.area = sum(present_in_M & present_in_F),
  category = c("M", "F"),
  fill = c("skyblue", "pink"),
  lty = "blank"
)

install.packages("VennDiagram")
library(VennDiagram)

venn.plot <- draw.pairwise.venn(
  area1 = sum(present_in_M),
  area2 = sum(present_in_F),
  cross.area = sum(present_in_M & present_in_F),
  category = c("M", "F"),
  fill = c("skyblue", "pink"),
  lty = "blank"
)


