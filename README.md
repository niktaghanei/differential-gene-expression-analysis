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
1. **Scatter plot** 

The dispersion plot shows the expected strong negative relationship between mean normalized counts and dispersion estimates, which is characteristic of well-behaved RNA-seq data under the negative binomial model.
Genes with low mean expression (left side of the plot) exhibit high dispersion (overdispersion), resulting in greater scatter of the gene-wise estimates (black points). As mean expression increases (moving right), dispersion decreases markedly, reflecting lower relative variability in highly expressed genes.
The fitted dispersion trend (red line) captures this relationship well, passing centrally through the cloud of points across the entire expression range.
Final shrunk dispersion estimates (blue points) closely follow the fitted trend, indicating appropriate and effective shrinkage, particularly beneficial for low-count genes where raw estimates are noisy. No severe dispersion outliers are apparent, suggesting the absence of major technical artifacts, contamination, or genes with unexplained extreme variability.
Overall, the dispersion estimation appears reliable, supporting confidence in the subsequent differential expression results (p-values and log2 fold changes). 

![Scatter plot](plots/dispersion_plot.png) 







2. **PCA plot**
PC1 and PC2 explain **53%** and **27%** of the total variance, respectively. Together, PC1 and PC2 capture **~80%** of the variation in the dataset, providing a highly representative 2D summary of the data.









3. **MA plots** 

The MA plot visualizes the relationship between mean normalized expression (baseMean, x-axis) and log2 fold change (y-axis) for the contrast Infected_Vehicle vs Uninfected_Vehicle, using shrunk log2 fold changes.

MA Plot: Infected vs Uninfected (Vehicle)
Key observations:
- Most genes cluster tightly around log2FC = 0 (gray band), as expected: the majority of genes do not show substantial expression changes under infection, reflecting stable housekeeping and baseline functions.
- A clear subset of genes exhibits strong differential expression, with log2FC values ranging up to ~+7 (upregulation) and down to ~-4 (downregulation). These extremes are biologically plausible in the context of a robust innate immune response to *Francisella tularensis* infection in THP-1 cells, where activation of inflammatory and antiviral pathways (such as NF-κB, cytokine, and interferon signaling) can drive large fold changes in key effector genes.
- Low-expressed genes (left side of the plot) show greater scatter in log2FC, which is typical due to higher relative variability (Poisson-like sampling noise + biological overdispersion) in lowly expressed transcripts.
- Highly expressed genes (right side) display much tighter fold changes, consistent with the behavior of stable housekeeping genes (such as ribosomal, actin, GAPDH) that are less responsive to infection.

The y-axis range was set to [-6, 8] to fully capture the observed fold changes without excessive empty space, while preserving detail across the expression spectrum. 
This pattern aligns with expectations for a strong biological perturbation (infection) and supports the validity of the differential expression results.  

![MA plot, Vehicle](https://github.com/niktaghanei/differential-gene-expression-analysis/blob/main/plots/ma%20plots/vehicle_MA_edited.png)















MA Plot: Infected vs Uninfected (Dillapiole)

Key observations:
Most genes cluster around log2FC = 0 (gray band), as expected for the majority of the transcriptome that remains stable. The y-axis range [-8, 10] was chosen to fully capture the observed extremes while maintaining clear visibility of all points and patterns.
- A substantial number of genes show strong differential expression, with log2FC values ranging from approximately -8 (downregulation) to +10 (upregulation). This wide spread indicates a robust transcriptional response to infection even in the presence of dillapiole.
- Low-expressed genes (left side) exhibit the greatest scatter in log2FC, typical due to higher relative variability in low-count transcripts (sampling noise + biological overdispersion).
- Highly expressed genes (right side) show much tighter fold changes, consistent with stable housekeeping genes that are generally unresponsive to infection/treatment.

Compared to the Vehicle condition, this plot shows:
- A broader cloud of points overall.
- A higher apparent number of genes with substantial |log2FC| > 1 (both up- and down-regulated).
- Reduced density of genes tightly clustered near log2FC = 0 (fewer unaffected genes).

This suggests dillapiole may lead to a more widespread or slightly amplified transcriptional response during infection, rather than suppressing changes. Biologically, this could reflect dillapiole's known role in dampening bacterial virulence gene expression, resulting in a less aggressive infection and a broader host response (activation of more diverse immune/inflammatory pathways).  

![MA plot, Dillapiole](https://github.com/niktaghanei/differential-gene-expression-analysis/blob/main/plots/ma%20plots/dillapiole_MA_edited.png)
