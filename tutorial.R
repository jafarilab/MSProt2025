#tutorial

#install package only if it is missing
# Install BiocManager if it's not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Function to install a package if it's not already installed
install_if_missing <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    BiocManager::install(package)
  }
}

# List of packages to install
packages <- c(
  "remotes",
  "tidyverse",
  "factoextra",
  "MsDataHub",
  "mzR",
  "rhdf5",
  "rpx",
  "MsCoreUtils",
  "QFeatures",
  "Spectra",
  "ProtGenerics",
  "PSMatch",
  "pheatmap",
  "limma",
  "MSnID",
  "RforMassSpectrometry/SpectraVis"
)

# Install each package if it's missing
for (pkg in packages) {
  install_if_missing(pkg)
}

########### installation done ############

#After installation, you can download some data that will be used in the latter chapter running the following:

library(rpx)
px <- PXDataset("PXD000001") ## Create a PXDataset object for the dataset with ID PXD000001
fn <- "TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01-20141210.mzML" ## Define the filename of the mzML file to be retrieved
mzf <- pxget(px, fn) ## Retrieve the mzML file specified by 'fn' from the PXDataset object 'px' and store it in 'mzf'
px <- PXDataset("PXD022816") ## Create a new PXDataset object for the dataset with ID PXD022816
pxget(px, grep("mzID", pxfiles(px))[1:3]) ## Retrieve the first three mzID files from the PXDataset object 'px'
pxget(px, grep("mzML", pxfiles(px))[1:3]) ## Retrieve the first three mzML files from the PXDataset object 'px'


#accessing data
library(rpx)
px <- PXDataset("PXD000001") 
pxfiles(px) ## Retrieve and display a list of all files associated with the PXDataset object 'px'
pxtax(px) ## Retrieve and display the taxonomy information associated with the PXDataset object 'px'
pxurl(px) ## ftp link #Retrieve and display the URL associated with the PXDataset object 'px'
pxref(px) ## publication. Retrieve and display reference information associated with the PXDataset object 'px'


fn <- "TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01-20141210.mzML"
mzf <- pxget(px, fn)

#some data may already be in MsDataHub package. 
#Data that is accessed through these hubs are cached centrally to avoid repeated downloads.
library("MsDataHub")
MsDataHub() #lists the available datasets in the above package
MsDataHub::TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.20141210.mzML.gz()

#Experiment packages
library("msdata")
#The msdata is an experiment package that directly ships raw data files relevant for both proteomics and metabolomics.
## proteomics raw data
proteomics() ## proteomics raw data
ident() ## proteomics identification data
quant() ## quantitative data

#3.1.1 The Spectra class
library(Spectra)
#Spectra is part of the R for Mass Spectrometry initiative. It defines the Spectra class that is used as a raw
#data abstraction, to manipulate MS data and metadata. The best way to learn about a data structure is to 
#create one by hand.
spd <- DataFrame(msLevel = c(1L, 2L), rtime = c(1.1, 1.2))
spd$mz <- list(c(100, 103.2, 104.3, 106.5), c(45.6, 120.4, 190.2))
spd$intensity <- list(c(200, 400, 34.2, 17), c(12.3, 15.2, 6.8))
spd

sp0 <- Spectra(spd) #convert this DataFrame into a Spectra object
sp0

#Spectra from mzML files
mzf
sp <- Spectra(mzf)
sp
length(sp) #7534 number of spectra


#adding information/variables to spectra such as retention time in minutes
sp$rtime_minute <- rtime(sp) / 60
sp$rtime_minute |> head()
table(msLevel(sp)) #many different types of mass spectrometry (MS) levels there are. how many MS1 and MS2
#in your data (sp) and counting how many scans exist for each level.
(sp2 <- filterMsLevel(sp, 2L)) #Filter for MS Level 2
sp[msLevel(sp) == 2L] #Show MS Level 2 Scans:
max(sp$basePeakIntensity)# Find Maximum Base Peak Intensity among the base peaks in your data.
rtime(sp)[1234] #retention time for the scan at index 1234.
plotSpectra(sp2[5404])
plotSpectra(sp2[1234])
plotSpectra(sp2[1230])
plotSpectra(sp2[1230:1233])
plotSpectra(sp[1])
plotSpectra(sp[1:4])
#visual plots of the spectra for specific scans, both from your filtered data (sp2) and the original data (sp). Each plot shows the mass spectrum for the specified scan(s).
  
### not limited to single file.elow we load data from two mzML files. 
#The MS data from both files in the Spectra is organized linearly 
#(first all spectra from the first file and then from the second). 
#The dataOrigin function can be used to identify spectra from the different data files. 
(fls <- dir(system.file("sciex", package = "msdata"), full.names = TRUE))

#Visualisation of raw MS data

#The chromatogram can be created by extracting the totIonCurrent and rtime variables for all MS1 spectra. Annotate the spectrum of interest.
#The chromatogram at the top displays the total ion current along the retention time. The vertical line identifies one scan in particular at retention time 1800.68 seconds (the 2807th scan).
with(spectraData(filterMsLevel(sp, 1)),
     plot(rtime, totIonCurrent, type = "l"))
abline(v = rtime(sp)[2807], col = "red") #random i guess


#The spectra on the second line represent the full MS1 spectrum marked by the red line. The vertical lines identify the 10 precursor ions that where selected for MS2 analysis. 
#The filterPrecursorScan() function can be used to retain a set parent (MS1) and children scans (MS2), as defined by an acquisition number. Use it to extract the MS1 scan of interest and all its MS2 children.
ms_2 <- filterPrecursorScan(sp, 2807)
ms_2

#Plot the MS1 spectrum of interest and highlight all the peaks that will be selected for MS2 analysis.

plotSpectra(sp[2807], xlim = c(400, 1000))
precursorMz(ms_2)
abline(v = precursorMz(ms_2)[-1], col = "grey") #[-1] excludes the first precursor m/z value (which is NA), meaning that it highlights all but the first precursor in grey.
abline(v = precursorMz(ms_2)[2], col = "red")

#Zoom in mz values 521.1 and 522.5 to reveal the isotopic envelope of that peak. The zoomed in shows one specific precursor peak of interest.
plotSpectra(sp[2807], xlim = c(521.2, 522.5), type = "l")
abline(v = precursorMz(ms_2)[2], col = "red")

#The MS2 spectra displayed along the two rows at the bottom are those resulting from the fragmentation of the 10 precursor peaks identified by the vertical bars above.
#The plotSpectra() function is used to plot all 10 MS2 spectra in one call.
plotSpectra(ms_2[-1])

## Filter MS2 level spectra and find any 2 MS2 spectra that have matching precursor peaks based on the precursor m/z values.
sp2 <- filterMsLevel(sp, 2L) # Filter for MS2 spectra (level 2) from the original mass spectrometry data
anyDuplicated(precursorMz(filterMsLevel(sp, 2))) # Check for duplicate precursor m/z values in the filtered MS2 data
# Returns the index of the first duplicate found, or 0 if none found
i <- which(precursorMz(sp2) == precursorMz(sp2)[37]) # Find the index of all MS2 spectra that have the same precursor m/z value
# as the one found at index 37
sp2i <- sp2[i] # Extract the matching MS2 spectra based on the indices found, which is 37 in this case


#Visualise the matching pair using the plotSpectraOverlay() and plotSpectraMirror() functions.
#Visualize the matching pair of MS2 spectra by overlaying them in the same plot
# Using red for the first spectrum and steelblue for the second spectrum
plotSpectraOverlay(sp2i, col = c("red", "steelblue")) 

# Optionally, visualize the matching pair of MS2 spectra in a mirrored format
# This helps in comparing the relative intensities and positions of the peaks
plotSpectraMirror(sp2i, col = c("red", "steelblue"))  #error

##############
#3.3 Raw data processing and manipulation
#example of data processing: we use below the pickPeaks() function. This function allows to convert profile mode MS data to centroid mode data (a process also referred to as centroiding).
plotSpectra(sp[2807], xlim = c(521.2, 522.5))
pickPeaks(sp[2807]) |>
  filterIntensity(1e7) |>
  plotSpectra(xlim = c(521.2, 522.5))


#####
#Chapter 4 Identification data
#Peptide identification is performed using third-party software - there is no package to run these searches directly in R.
#The example below illustrates this for 3 mzML files to be searched using MSGFplus:

(mzmls <- paste0("file_", 1:3, ".mzML"))
(mzids <- sub("mzML", "mzid", mzmls))
paste0("java -jar /path/to/MSGFPlus.jar",
       " -s ", mzmls,
       " -o ", mzids,
       " -d uniprot.fas",
       " -t 20ppm",
       " -m 0",
       " int 1")

#4.1 Identification data.frame
#Let’s use the identification from msdata:
idf <- MsDataHub::TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.20141210.mzid()
#The easiest way to read identification data in mzIdentML (often abbreviated with mzid) into R is to read it with the 
#PSM() constructor function from the PSMatch package55 Previously named PSM.. 
#The function will parse the file and return a DataFrame.
library(PSMatch)
id <- PSM(idf)
dim(id)
names(id)
#The PSM data are read as is, without any filtering. As we can see below, we still have all the hits from the forward 
#and reverse (decoy) databases.

table(id$isDecoy)
#4.2 Keeping all matches
table(table(id$spectrumID))
#The data contains also contains multiple matches for several spectra. The table below shows the number of number of spectra that have 1, 2, … up to 5 matches.

table(table(id$spectrumID))
#example Below, we can see how scan 1774 has 4 matches, all to sequence RTRYQAEVR, which itself matches to 4 different proteins:
i <- which(id$spectrumID == "controllerType=0 controllerNumber=1 scan=1774")
data.frame(id[i, ])[1:5]
#If the goal is to keep all the matches, but arranged by scan/spectrum, one can reduce the PSM object by the spectrumID variable, so that each scan correponds to a single row that still stores all values66 The rownames aren’t needed here are are removed to reduce to output in the the next code chunk display parts of id2.:
id2 <- reducePSMs(id, id$spectrumID)
id2
#The resulting object contains a single entry for scan 1774 with information for the multiple matches stored as lists within the cells.
j <- which(id2$spectrumID == "controllerType=0 controllerNumber=1 scan=1774")
id2[j, ]
id2[j, "DatabaseAccess"]

#The is the type of complete identification table that could be used to annotate an raw mass spectrometry Spectra object, as shown below.

#4.3 Filtering data
#Often, the PSM data is filtered to only retain reliable matches. The MSnID package can be used to set thresholds to attain user-defined PSM, peptide or protein-level FDRs. Here, we will simply filter out wrong identification manually.
#Here, the filter() from the dplyr package comes very handy. We will thus start by converting the DataFrame to a tibble.


library("dplyr")
id_tbl <- tidyr::as_tibble(id)
id_tbl
#Remove decoy hits
id_tbl <- id_tbl |>
  filter(!isDecoy)
id_tbl
#Keep first rank matches

id_tbl <- id_tbl |>
  filter(rank == 1)
id_tbl


#Remove shared peptides. Start by identifying scans that match different proteins. For example scan 4884 matches proteins XXX_ECA3406 and ECA3415. Scan 4099 match XXX_ECA4416_1, XXX_ECA4416_2 and XXX_ECA4416_3. Then remove the scans that match any of these proteins.

# Step 1: Group the data by spectrumID to analyze each spectrum individually
mltm <-
  id_tbl |>
  group_by(spectrumID) |>
  
  # Step 2: Create a new column 'nProts' that counts the number of unique proteins (or database accesses)
  # associated with each spectrumID
  mutate(nProts = length(unique(DatabaseAccess))) |>
  
  # Step 3: Filter to keep only those spectrumIDs that match more than one unique protein
  filter(nProts > 1) |>
  
  # Step 4: Select only the relevant columns: spectrumID and the count of unique proteins (nProts)
  select(spectrumID, nProts)

# Display the resulting data frame containing spectrumIDs with multiple associated proteins
mltm

id_tbl <-
  id_tbl |>
  filter(!spectrumID %in% mltm$spectrumID)
id_tbl
#Which leaves us with 2666 PSMs.
#his can also be achieved with the filterPSMs() function, or the individual filterPsmRank(), filterPsmDecoy and filterPsmShared() functions:

id_filtered <- filterPSMs(id) #also says how many decoy hits removed
#The describePeptides() and describeProteins() functions from the PSMatch package provide useful summaries of preptides and proteins in a PSM search result.
#describePeptides() gives the number of unique and shared peptides and for the latter, the size of their protein groups:
describePeptides(id_filtered)

#describeProteins() gives the number of proteins defined by only unique, only shared, or a mixture of unique/shared peptides:

describeProteins(id_filtered)
#The Understanding protein groups with adjacency matrices PSMatch vignette provides additional tools to explore how proteins were inferred from peptides.

#Compare the distribution of raw identification scores of the decoy and non-decoy hits. Interpret the figure.

library(ggplot2)
ggplot(id, aes(x = MS.GF.RawScore,
               colour = isDecoy)) +
  geom_density()

#The tidyverse tools are fit for data wrangling with identification data. Using the above identification dataframe, calculate the length of each peptide (you can use nchar with the peptide sequence sequence) and the number of peptides for each protein (defined as DatabaseDescription). Plot the length of the proteins against their respective number of peptides.

suppressPackageStartupMessages(library("dplyr"))
iddf <- as_tibble(id_filtered) |>
  mutate(peplen = nchar(sequence))
npeps <- iddf |>
  group_by(DatabaseAccess) |>
  tally()
iddf <- full_join(iddf, npeps)
## Joining with `by = join_by(DatabaseAccess)`
library("ggplot2")
ggplot(iddf, aes(x = n, y = DBseqLength))  + geom_point()


#4.4 Adding identification data to raw data







































