library(rpx)
library(msdata)

library(Spectra)
library(tidyverse)

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

#px <- PXDataset("PXD060654")
#f<-pxget(px, grep("mzid", pxfiles(px)))

mzid_obj <- mzID("D:\\MS_course\\MS_2369_mock-vs-flg-combined_031424.mzid")
#id <-"D:\\MS_course\\MS_2369_mock-vs-flg-combined_031424.mzid"
#psm<- PSMatch::PSM(id)

psms <- flatten(mzid_obj)

idtbl <- as_tibble(psms)
names(idtbl) # all the column names
table(idtbl$rank, useNA = "ifany") #ranks mentioned

#idtbl <- idtbl |>  filter(!isDecoy) - did not find any

idtblr1 <- idtbl |>  filter(rank == 1) # filtering only rank 1 peptides

dplyr::count(idtblr1, spectrumid)

dplyr::count(idtblr1, spectrumid) |>   filter(n >1) #remove repeated spectrums to the protein

(mltm <- dplyr::count(idtblr1, spectrumid) |>     filter(n > 1) |>     pull(spectrumid))
#bascially creating object that pulls all the spectrum ids that match more than once to the protein

(idtbl_fil <- idtblr1 |>    filter(!spectrumid %in% mltm)) #Filter out IDs matched with more than 1 protein

length(unique(idtbl_fil$pepseq))       # Peptides
length(unique(idtbl_fil$accession))


# Read FASTA
fasta <- readAAStringSet("D:\\file63881f23791e.fasta")

# Extract headers
protein_names <- names(fasta)

# Optionally just IDs (e.g., split at first space)
cont_prot<- sapply(strsplit(protein_names, "\\|"), `[`, 2)

present_in_psms <- cont_prot %in% idtbl_fil$accession
nc_idtbl_fil <- idtbl_fil[!(idtbl_fil$accession %in% cont_prot), ]

length(unique(nc_idtbl_fil$pepseq))       # Peptides
length(unique(nc_idtbl_fil$accession))


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
# so after everything we have21107 unique peptides and 3288 unique proteins

maxquant_data<-read.delim("D:\\MS_course\\combined\\txt\\proteinGroups.txt", header = TRUE, sep = "\t")
names(maxquant_data)

# Remove decoys
maxquant_fil <- maxquant_data %>% filter(is.na(Reverse) | Reverse == "")

#remove contaminant
nc_maxquant_fil <- maxquant_fil[!(maxquant_fil$Majority.protein.IDs %in% cont_prot), ]

#number of unique peptides
sum(nc_maxquant_fil$Unique.peptides)

#number of razor.uniqu.peptides
sum(nc_maxquant_fil$Razor...unique.peptides)

#To table the number of protein grouped together due to inability to separate or peptide sequence is similar
table(nc_maxquant_fil$Number.of.proteins)

#viewing score distribution
ggplot(nc_maxquant_fil, aes(x = `Score`)) +
  geom_density(fill = "lightblue") +
  labs(title = " Score Distribution")+
  geom_vline(xintercept = 10, 
             color = "red", linetype = "dashed", size = 0.5)


#filtering based on score >10
high_conf_maxquant <- subset(nc_maxquant_fil, ((`Q.value` < 0.01) & (`Score`>10)))

#without filtering
high_conf_maxquant1 <- nc_maxquant_fil

lfq_cols <- grep("^LFQ\\.intensity\\.", names(high_conf_maxquant), value = TRUE)
exprs <- high_conf_maxquant[, lfq_cols]
rownames(exprs) <- high_conf_maxquant$Majority.protein.IDs

# Convert to numeric and log2 transform
enrichment <- exprs[, "LFQ.intensity.flg"] / exprs[, "LFQ.intensity.mock"]

log2FC <-   log2(as.matrix(enrichment))

deg_results <- data.frame(
  Protein_names = high_conf_maxquant$Protein.names,
  Gene_names = high_conf_maxquant$Gene.names,
  log2FC = log2FC
)

order_results <- deg_results[order(deg_results$log2FC, decreasing = TRUE), ]
head(order_results)

#without filtering
lfq_cols1 <- grep("^LFQ\\.intensity\\.", names(high_conf_maxquant1), value = TRUE)
exprs1 <- high_conf_maxquant1[, lfq_cols1]
rownames(exprs1) <- high_conf_maxquant1$Protein.IDs

# Convert to numeric and log2 transform
enrichment1 <- exprs1[, "LFQ.intensity.flg"] / exprs1[, "LFQ.intensity.mock"]

log2FC1 <-   log2(as.matrix(enrichment1))

deg_results1 <- data.frame(
  Protein_names = high_conf_maxquant1$Protein.names,
  Gene_names = high_conf_maxquant1$Gene.names,
  log2FC = log2FC1
)

order_results1 <- deg_results1[order(deg_results$log2FC, decreasing = TRUE), ]
head(order_results1)

#with filtering based on score
top30 <- deg_results %>%  slice_max(log2FC, n = 30)

# Plot as dot plot
ggplot(top30, aes(x = log2FC, y = reorder(Gene_names, log2FC))) +
  geom_point(size = 3, color = "darkblue") +
  theme_minimal() +
  labs(title = "Top 30 Proteins by Log2 Fold Change",
       x = "Log2 Fold Change (FLG vs MOCK)",
       y = "Gene Name") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray")


#without filtering based on score
top30_1 <- deg_results1 %>%  slice_max(log2FC, n = 30)

# Plot as dot plot
ggplot(top30_1, aes(x = log2FC, y = reorder(Gene_names, log2FC))) +
  geom_point(size = 3, color = "darkblue") +
  theme_minimal() +
  labs(title = "Top 30 Proteins by Log2 Fold Change",
       x = "Log2 Fold Change (FLG vs MOCK)",
       y = "Gene Name") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray")


## Removing contaminants from the evidence file
m<-read.delim("D:\\MS_course\\combined\\txt\\evidence.txt", header = TRUE, sep = "\t")

# Remove decoys
n <- m %>% filter(is.na(Reverse) | Reverse == "")

#Remove contaminants
n <- n %>% filter(is.na(Potential.contaminant) | Potential.contaminant == "")

n_fil <- n[!(n$Leading.razor.protein %in% cont_prot), ]