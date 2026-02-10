# differential gene expression analysis

This project was designed as a hands-on practice to gain familiarity with RNA-seq data analysis in R and to develop skills in both coding and biological interpretation.

The main biological question addressed in this study is:
**How does infection with Francisella tularensis affect gene expression levels in the human monocytic cell line THP-1?**

THP-1 cells are derived from patients with acute monocytic leukemia and are commonly used as a model for studying immune responses in monocytes.

# **Dataset**
The RNA-seq dataset used in this project was obtained from the NCBI Gene Expression Omnibus (GEO):

- Accession number: GSE306199
- Organism: Homo sapiens
- Cell line: THP-1
  
# Experimental conditions:
- Infected vs Uninfected
- Treated with Dillapiole vs Vehicle control

This results in four experimental groups (n=5 per group):
- Uninfected + Vehicle
- Uninfected + Dillapiole
- Infected + Vehicle
- Infected + Dillapiole

# Methods (Overview)
A raw count matrix was generated from the dataset and used as input for differential expression analysis using the DESeq2 package in R.
Gene identifiers were retained while sample name columns were excluded to meet DESeq2 input requirements.

A metadata table was manually constructed containing experimental information for each sample.
Column names of the count matrix were matched to row names of the metadata to ensure consistency.

Reference levels were defined as:
- Uninfected for infection status
- Vehicle for treatment condition

This allowed all other conditions to be interpreted relative to the baseline (no infection, no treatment).

Lowly expressed genes were filtered out by removing genes with fewer than 10 total counts across all samples.

Differential expression analysis was performed using DESeq2, followed by variance stabilizing transformation (VST) for downstream visualization and exploratory analysis.

# Plot analysis


