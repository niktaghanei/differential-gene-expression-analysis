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

**MA Plot: Infected vs Uninfected (Vehicle)**
Key observations:
- Most genes cluster tightly around log2FC = 0 (gray band), as expected: the majority of genes do not show substantial expression changes under infection, reflecting stable housekeeping and baseline functions.
- A clear subset of genes exhibits strong differential expression, with log2FC values ranging up to ~+7 (upregulation) and down to ~-4 (downregulation). These extremes are biologically plausible in the context of a robust innate immune response to *Francisella tularensis* infection in THP-1 cells, where activation of inflammatory and antiviral pathways (such as NF-κB, cytokine, and interferon signaling) can drive large fold changes in key effector genes.
- Low-expressed genes (left side of the plot) show greater scatter in log2FC, which is typical due to higher relative variability (Poisson-like sampling noise + biological overdispersion) in lowly expressed transcripts.
- Highly expressed genes (right side) display much tighter fold changes, consistent with the behavior of stable housekeeping genes (such as ribosomal, actin, GAPDH) that are less responsive to infection.

The y-axis range was set to [-6, 8] to fully capture the observed fold changes without excessive empty space, while preserving detail across the expression spectrum. 
This pattern aligns with expectations for a strong biological perturbation (infection) and supports the validity of the differential expression results.  

![MA plot, Vehicle](https://github.com/niktaghanei/differential-gene-expression-analysis/blob/main/plots/ma%20plots/vehicle_MA_edited.png)















**MA Plot: Infected vs Uninfected (Dillapiole)**

Key observations:
Most genes cluster around log2FC = 0 (gray band), as expected for the majority of the transcriptome that remains stable. The y-axis range [-8, 10] was chosen to fully capture the observed extremes while maintaining clear visibility of all points and patterns.
- A substantial number of genes show strong differential expression, with log2FC values ranging from approximately -8 (downregulation) to +10 (upregulation). This wide spread indicates a robust transcriptional response to infection even in the presence of Dillapiole.
- Low-expressed genes (left side) exhibit the greatest scatter in log2FC, typical due to higher relative variability in low-count transcripts (sampling noise + biological overdispersion).
- Highly expressed genes (right side) show much tighter fold changes, consistent with stable housekeeping genes that are generally unresponsive to infection/treatment.

Compared to the Vehicle condition, this plot shows:
- A broader cloud of points overall.
- A higher apparent number of genes with substantial |log2FC| > 1 (both up- and down-regulated).
- Reduced density of genes tightly clustered near log2FC = 0 (fewer unaffected genes).

This suggests Dillapiole may lead to a more widespread or slightly amplified transcriptional response during infection, rather than suppressing changes. Biologically, this could reflect Dillapiole's known role in dampening bacterial virulence gene expression, resulting in a less aggressive infection and a broader host response (activation of more diverse immune/inflammatory pathways).  

![MA plot, Dillapiole](https://github.com/niktaghanei/differential-gene-expression-analysis/blob/main/plots/ma%20plots/dillapiole_MA_edited.png)











**MA Plot: Dillapiole vs Vehicle (Uninfected)**

This plot shows the direct/baseline effect of Dillapiole treatment on uninfected THP-1 cells.

- Classic mean-variance trend: high scatter at low expression (left), tight fold changes at high expression (right)
- y-range [-10, 8] captures all observed changes, with some genes showing strong regulation (fold changes up to ~1000×).
- Density near log2FC = 0 is moderate. Dillapiole perturbs a noticeable fraction of the transcriptome even in the absence of infection.

**Comparison: Uninfected vs Vehicle MA Plots** 

To distinguish between dillapiole's direct effect on host cells and its potential modulation of the infection response, I compare:
- Infected + Dillapiole vs Infected + Vehicle (treatment effect during infection)
- Uninfected + Dillapiole vs Uninfected + Vehicle (baseline treatment effect without infection)

Key differences:
- **Downregulation stronger in Uninfected**: More genes and larger negative log2FC in the Uninfected plot, indicating Dillapiole alone suppresses a broader portion of the transcriptome in uninfected cells.
- **Fewer unaffected genes in Uninfected**: Reduced density near log2FC = 0 (gray band) in Uninfected compared to Vehicle, showing Dillapiole perturbs more genes even without infection.
- **Upregulation more dispersed in Vehicle**: Infected cells without treatment show wider spread in positive log2FC, consistent with robust activation of immune/inflammatory pathways.
- **Downregulation more clustered near zero in Vehicle**: Limited and weaker downregulation during pure infection, typical of innate immune activation focused on upregulation rather than suppression.

**Interpretation**:
Dillapiole exerts a predominantly suppressive effect on uninfected THP-1 cells, while infection drives widespread upregulation with minimal downregulation. These contrasting patterns suggest Dillapiole does not simply mimic or add to the infection response; its presence during infection likely modulates the host transcriptome in a way that differs from its baseline action.

**Comparison: Uninfected vs Infected MA Plots**

Key differences:
- Much higher density of unaffected genes (gray band near log2FC = 0) in the Infected plot, fewer genes overall show significant changes during infection + Dillapiole.
- Upregulation in Uninfected reaches higher extremes (up to ~+7) with greater dispersion, especially among low-expression genes (left side); in Infected, upregulation is weaker (max ~+5) and more focused on high-expression genes (right side).
- Downregulation in Uninfected is stronger among high-expression genes (right side); in Infected, it shifts toward low-expression genes (left side) with reduced overall spread.

**Interpretation**:
Dillapiole alone induces broad changes in uninfected cells, with strong upregulation in low-expression genes and downregulation in high-expression genes. During infection, its effect becomes more limited: fewer genes are perturbed, upregulation is attenuated and shifts to high-expression genes, and downregulation focuses on low-expression genes. This suggests Dillapiole modulates the host response during infection (likely by dampening bacterial virulence) leading to a less widespread and more balanced transcriptional reaction compared to its baseline action on uninfected cells. 

![MA Plot, Uninfected](https://github.com/niktaghanei/differential-gene-expression-analysis/blob/main/plots/ma%20plots/uninfected_MA_edited.png) 











**MA Plot: Infected (Infected + Vehicle vs Infected + Dillapiole)**

This MA plot shows shrunk log2 fold changes (y-axis) against mean normalized counts (baseMean, x-axis) for the contrast Infected + Vehicle vs Infected + Dillapiole. It reveals the specific effect of Dillapiole treatment on infected THP-1 cells.
Overall observations:

The y-axis range [-8, 5] fully captures the observed fold changes, with the most extreme points reaching ~−7 (downregulation) and ~+5 (upregulation).
The gray band (genes near log2FC = 0) is relatively wide and dense, indicating a large proportion of genes remain unaffected or show only small changes.
Downregulated genes (negative log2FC) show much greater dispersion and spread across the plot compared to upregulated genes.
Upregulated genes (positive log2FC) are more clustered near the gray band (closer to zero) and less dispersed overall.
Low-expression genes (left side, low baseMean) are predominantly downregulated (more blue points on the left).
High-expression genes (right side, high baseMean) are more likely to be upregulated (more blue points on the right).

This pattern suggests that Dillapiole's presence during infection leads to a more limited and asymmetric response: stronger and more widespread downregulation (especially in low-expression genes), while upregulation is weaker, more constrained, and shifted toward high-expression genes. 

**Comparison with Vehicle Plot (Vehicle + Infected vs Vehicle + Uninfected)** 

Vehicle (pure infection effect without treatment):
- Wider upregulation (up to ~+7) with greater dispersion, especially in low-expression genes.
- Downregulation limited (max ~−4) and more clustered near zero.
- More dispersed overall changes, with low-expression genes driving most upregulation.

Key differences when comparing to Infected:
Density near zero (gray band): Much higher in Infected, fewer genes are significantly perturbed when Dillapiole is present during infection.
Upregulation: Stronger and more dispersed in Vehicle (max +7, more spread in low-expression genes); weaker and more clustered in Infected (max +5, shifted to high-expression genes).
Downregulation: Limited and focused in Vehicle (max −4, clustered near zero); stronger, more dispersed, and shifted to low-expression genes in Infected (down to −7, broader spread). 

Overall effect: Pure infection (Vehicle) drives a robust, activation-dominated response (stronger upregulation, limited downregulation). Adding Dillapiole during infection (Infected) attenuates upregulation, enhances downregulation, and reduces the total number of affected genes resulting in a more suppressed and less widespread transcriptional response.

**Interpretation:**
The presence of Dillapiole during infection dampens the strong activation seen in pure infection and shifts the balance toward more pronounced downregulation, particularly in low-expression genes. This is consistent with Dillapiole reducing bacterial virulence, leading to a less aggressive host response rather than broad amplification or suppression. The reduced number of perturbed genes and constrained upregulation suggest a modulatory, protective effect on the infected cell's transcriptome.

![MA plot, Infected](https://github.com/niktaghanei/differential-gene-expression-analysis/blob/main/plots/ma%20plots/infected_MA_edited.png)



4. **Volcano plots**

Volcano Plots: Differential Expression Results 

Volcano plots are one of the most common and powerful ways to visualize RNA-seq differential expression results. They combine two key pieces of information in a single scatter plot:

- **Horizontal axis (x)**: log2 Fold Change (log2FC)  
  Shows the magnitude and direction of expression change between two conditions.  
  - Positive values (right side): gene is upregulated (higher expression in the first group)
  - Negative values (left side): gene is downregulated (lower expression in the first group)  
  - A change of log2FC = 1 means ~2-fold upregulation; log2FC = -1 means ~2-fold downregulation

- **Vertical axis (y)**: -log10(adjusted p-value) or -log10(padj)  
  Shows statistical significance (how unlikely the observed change is due to chance) 
  - Higher points = more significant (smaller padj)

**Classic thresholds** used in most RNA-seq studies:
- padj < 0.05 (or sometimes 0.01 for stricter control of false positives)
- |log2FC| > 1 (fold change > 2 or < 0.5, considered biologically meaningful)

Genes that pass both thresholds (significant + large change) are usually colored red. Genes with only large change but no significance are green, genes with only significance but small change are light blue, and non-significant/small-change genes are gray.


**Volcano Plot: Infected vs Uninfected (Vehicle)**

This volcano plot shows shrunk log2 fold changes vs -log10(padj) for the contrast Infected_Vehicle vs Uninfected_Vehicle (pure infection effect).

**Key observations**:
- Strong asymmetric response: upregulation dominates (right side), with many genes reaching high significance (red points up to -log10(padj) ~60).
- Downregulation is limited (few red points on left).
- Labeled genes include key immune/stress-related markers:  
  - **NFKBIA**: NF-κB inhibitor; feedback control of inflammation.  
  - **IL4I1**: Immune regulation and arginine metabolism; modulates inflammation.  
  - **HSPA1B**: Heat shock protein; protects cells under infection stress.  
  - **ERCC1**: DNA repair; response to infection-induced damage.   
  - **TBC1D17**: Regulates autophagy/vesicle trafficking; potential defense against intracellular pathogen.

Most significant DEGs are upregulated and enriched in immune activation, inflammation control, cellular stress response, and DNA repair pathways consistent with a robust innate immune reaction to intracellular *Francisella tularensis* infection in THP-1 cells.
Low-expression genes show higher scatter, but significant changes concentrate in medium- to high-expression genes (right-side clustering).

No technical artifacts visible; shrinkage reduces noise in low-count genes. This plot highlights infection as a strong activator of host gene expression.

![Volcano plot, Vehicle](https://github.com/niktaghanei/differential-gene-expression-analysis/blob/main/plots/volcano%20plots/vehicle_withsymbol_volcano.png)




**Volcano Plot: Infected vs Uninfected (Dillapiole)**


Key observations in this plot:
- Balanced up- and downregulation: unlike the Vehicle plot (stronger upregulation), here both directions are more symmetric, with roughly equal numbers of significant genes.
- Extremely high significance: -log10(padj) reaches ~100 (much higher than Vehicle's ~60), changes are highly significant.
- Upregulation: reaches up to ~+10, with key genes like **BEST1** (calcium-activated chloride channel, ion homeostasis), **KIF1B** (kinesin motor, vesicle/mitochondria transport), **MSC-AS1** (lncRNA, gene regulation), **ICAM4-AS1** (lncRNA, cell adhesion), **TRPA1** (cation channel, inflammation/pain sensing), **MRNIP** (DNA repair/stress response).
- Downregulation: reaches down to ~-8, with genes like **ENOSF1** (fucose metabolism), **LRRC41** (ubiquitin ligase component), **WDR5** (histone methylation), **LCORL** (transcription factor).
- Green region (large FC but non-significant) is thinner and more dispersed than in Vehicle -> fewer large non-significant changes.

**Interpretation**:
Infection in the presence of Dillapiole produces a more balanced and highly significant transcriptional response compared to pure infection (Vehicle). Upregulation is strong but less extreme, while downregulation is more prominent and dispersed. This suggests Dillapiole modulates the host response, reducing extreme activation while enhancing certain suppressive changes, likely by dampening bacterial virulence and altering the infection dynamics. 

![Volcano plot, Dillapiole](https://github.com/niktaghanei/differential-gene-expression-analysis/blob/main/plots/volcano%20plots/dillapiole_withsymbol_volcano.png)



**Volcano Plot: Uninfected (Uninfected + Dillapiole vs Uninfected + Vehicle)**

Key observations:
- Balanced up- and downregulation: roughly equal numbers of significant genes in both directions (unlike Vehicle's upregulation dominance).
- High significance: -log10(padj) reaches ~80, very strong statistical support for changes.
- Upregulation: reaches up to ~+10, with key genes like **KIF1B** (vesicle/mitochondria transport), **HDGFL3** (cell growth/repair), **TRPA1** (cation channel, inflammation sensing), **WDR13** (gene regulation/signaling).
- Downregulation: reaches down to ~-8, with top genes like **INTS8** (snRNA processing), **LRRC41** (ubiquitin ligase component), **WDR5** (histone methylation), **LCORL** (transcription factor).
- Green region (large FC but non-significant) is thinner and more dispersed than in Vehicle —> fewer large non-significant changes.

**Comparison with Vehicle (pure infection)**:
- Vehicle shows stronger upregulation (max +7) and limited downregulation (max -4), with low-expression genes driving most changes.
- Uninfected shows more balanced changes, stronger downregulation, and high-expression genes more affected.
- Fewer unaffected genes (gray band) in Uninfected, Dillapiole perturbs more of the transcriptome in healthy cells.
- Shared downregulated genes (e.g., LRRC41, WDR5, MLX, VRK2) maintain similar positions, suggesting these are direct Dillapiole effects independent of infection.

**Interpretation**:
Dillapiole alone induces a broad, balanced response in uninfected cells, with strong downregulation (especially high-expression genes) and upregulation in low-expression genes. This contrasts with pure infection's activation-dominated pattern. Dillapiole's baseline suppression may contribute to its modulatory role during infection, where the response becomes more controlled and less extreme. 

![Volcano plot, Uninfected](https://github.com/niktaghanei/differential-gene-expression-analysis/blob/main/plots/volcano%20plots/uninfected_withsymbol_volcano.png)
