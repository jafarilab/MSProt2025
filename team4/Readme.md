**Team 4: Hao, Tiialotta, Oluwatosin**

# 🧬 Capstone Project: R-based Mass Spectrometry Proteomics Workflow

This document provides a step-by-step guide for your capstone project on analyzing proteomics data using R.

---

## 🔁 Workflow Overview

### 1. 📁 Dataset Acquisition  
**Input:** PRIDE PXD Identifier  
↓  
(*Example: `PXD0123456`, accessible via [PRIDE](https://www.ebi.ac.uk/pride/)*)

Use the `rpx` package to access metadata and download a dataset that includes `.mzID` files from a published study.

---

### 2. 🧱 PSM Object Creation & Preprocessing  
**Goal:** Generate a PSM (Peptide-Spectrum Match) object from `.mzID` files  
↓  
- Convert `.mzID` files into `PSM` objects
- Assess:
  - Number of decoy hits
  - Score distributions
  - PSM rank
- Apply filtering based on FDR or identification score

---

### 3. 🧬 Protein & Peptide Identification  
**Goal:** Determine identified peptides and proteins  
↓  
- Count the number of identified peptides and proteins
- Review peptide-to-protein mapping:
  - Razor proteins
  - Protein groups

---

### 4. 🔄 QFeature Aggregation (Optional)  
**Goal:** Use quantification data if available  
↓  
Use `aggregateFeatures()` from the `QFeatures` package to aggregate:
- PSMs ➡️ Peptides ➡️ Proteins

---

### 5. 🧼 Normalization & Imputation (Optional)  
**Goal:** Correct for technical variation and handle missing values  
↓  
- Normalize using `normalize()` or `normalize_vsn()`
- Impute missing values with `impute()`

---

### 6. 🧪 Protein Inference & Quantification (Optional)  
**Goal:** Summarize and quantify proteins  
↓  
- Aggregate peptide-level intensities into protein-level quantities
- Optionally annotate protein IDs using external databases (e.g., UniProt)

---

### 7. 📊 Statistical Analysis (Optional)  
**Goal:** Identify differentially abundant proteins  
↓  
- Perform statistical testing (e.g., using `test_diff()`)
- Filter by:
  - Log2 fold-change
  - Adjusted p-value (e.g., FDR < 0.05)

---

### 8. 📈 Visualization & Export (Optional)  
**Goal:** Visualize data for interpretation  
↓  
- Create visual summaries:
  - PCA plots
  - Volcano plots
  - Heatmaps  
- Use packages such as `ggplot2`, `DEP`, or `limma`

---

### 9. 📖 Comparison with Published Results  
**Goal:** Benchmark your findings  
↓  
Create a comparison table including:
- Number of spectra identified
- Number of peptides and proteins
- Key figures or results (if available)
- Comments on reproducibility

---

### 10. 📝 Final Report & Interpretation  
**Goal:** Reflect on your analysis  
↓  
Submit the following:
- Source code (R script or RMarkdown) in your team folder
- A 1–2 page report that includes:
  - Background of the dataset and study
  - Summary table of key results
  - Discussion:
    > Why do your results match or differ from the original publication?

---

✅ **Tip:** Convert your `.Rmd` into a clean `README.md` using `knitr::knit("README.Rmd")`.


