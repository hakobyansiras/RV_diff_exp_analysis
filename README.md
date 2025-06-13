# 🧬 Reproducible Code: *The Effects of Space Radiation on the Transcriptome of Heart Right Ventricle Tissue*

This repository contains the analysis code and necessary files to reproduce the results presented in the paper:

**"The Effects of Space Radiation on the Transcriptome of Heart Right Ventricle Tissue"**

## 📂 Contents

- **`rv_paper_results_reproduction.R`**  
  The main analysis script. Running this script in R will reproduce all figures and tables included in the manuscript.

- **`session_info.txt`**  
  Contains the R session information, including versions of all packages and their dependencies used during the analysis.

- **Custom GTF generation script**  
  Code used to generate a custom GTF file following the approach described by Sigurgeirsson et al.  
  [DOI: 10.1371/journal.pone.0091851](https://doi.org/10.1371/journal.pone.0091851)

## 📥 Reference Files

The GTF file used in the custom read-counting step:  
**`Mus_musculus.GRCm39.110.gtf`**  
You can download this file from the [Ensembl](https://www.ensembl.org) or [NCBI](https://www.ncbi.nlm.nih.gov) genome database.

## 📦 Requirements

- R (≥ 4.0.0)
- Required R packages listed in `session_info.txt`

To ensure full reproducibility, we recommend installing the exact versions of the packages as recorded in `session_info.txt`.
