library(rpx)
library(msdata)

library(Spectra)
library(tidyverse)
#library(cleaver)
library(MSnbase)
library(SpectraVis)
library(tidyverse)
library(mzR)
library(PSMatch)
library(mzID)
library(Biostrings)

#px <- PXDataset("PXD060654")
#f<-pxget(px, grep("mzid", pxfiles(px)))

mzid_obj <- mzID("D:\\MS_course\\MS_2369_mock-vs-flg-combined_031424.mzid")
#id <-"D:\\MS_course\\MS_2369_mock-vs-flg-combined_031424.mzid"
#psm<- PSMatch::PSM(id)

psms <- flatten(mzid_obj)
head(psms)
dim(psms)
table(psms$spectrumid)
idtbl <- as_tibble(psms)
names(idtbl) # all the column names
table(idtbl$rank, useNA = "ifany") #ranks mentioned

#idtbl <- idtbl |>  filter(!isDecoy) - did not find any

idtblr1 <- idtbl |>  filter(rank == 1) # filtering only rank 1 peptides

dplyr::count(idtblr1, spectrumid)

dplyr::count(idtblr1, spectrumid) |>   filter(n >1) #almost 29571 scans match more than 1 protein

(mltm <- dplyr::count(idtblr1, spectrumid) |>     filter(n > 1) |>     pull(spectrumid))
#bascially creating object that pulls all the spectrum ids that match more than once to the protein

(idtbl_fil <- idtblr1 |>    filter(!spectrumid %in% mltm)) #Filter out IDs matched with more than 1 protein

length(unique(idtbl_fil$pepseq))       # Peptides
length(unique(idtbl_fil$accession))

maxquant_data<-read.delim("D:\\MS_course\\combined\\txt\\msmsScans.txt", header = TRUE, sep = "\t")
names(maxquant_data)



# Read FASTA
fasta <- readAAStringSet("D:\\file63881f23791e.fasta")

# Extract headers
protein_names <- names(fasta)

# Optionally just IDs (e.g., split at first space)
protein_ids <- sapply(strsplit(protein_names, " "), `[`, 1)
cont_prot<- sapply(strsplit(protein_names, "\\|"), `[`, 2)

present_in_psms <- cont_prot %in% idtbl_fil$accession
nc_idtbl_fil <- idtbl_fil[!(idtbl_fil$accession %in% cont_prot), ]

length(unique(nc_idtbl_fil$pepseq))       # Peptides
length(unique(nc_idtbl_fil$accession))


mzml <- openMSfile("D:\\MS_course\\MS_2369_mock-vs-flg-combined_031424.mzML")
header(mzml)
