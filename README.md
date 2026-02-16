# differential gene expression analysis

This project was intended to be a practical exercise to become familiar with RNA-seq data analysis in R and to develop skills in both coding and biological interpretation.

The main biological question being asked in this study is:
**What is the effect of Francisella tularensis infection on gene expression levels in the human monocytic cell line THP-1?**

The THP-1 cell line is derived from patients with acute monocytic leukemia and is widely used as a model system to study immune responses in monocytes.

# **Dataset**
The RNA-seq dataset used in this project was downloaded from the NCBI Gene Expression Omnibus (GEO) database:

- Accession number: GSE306199
- Organism: Homo sapiens
- Cell line: THP-1
  
# Experimental conditions:
- Infected vs Uninfected
- Treated with Dillapiole vs Vehicle control
  
This gives a total of four experimental groups (n=5 per group):
- Uninfected + Vehicle
- Uninfected + Dillapiole
- Infected + Vehicle
- Infected + Dillapiole

# Methods (Overview)
A raw count matrix was created from the dataset and used as input for differential expression analysis using the DESeq2 package in R.
Gene identifiers were kept while sample name columns were removed to conform to the DESeq2 input format.

A metadata table was manually constructed with experimental information for each sample.
The column names of the count matrix were matched to row names of the metadata to ensure consistency.

Reference levels were defined as:
- Uninfected for infection status
- Vehicle for treatment condition

This allowed all other conditions to be interpreted relative to the baseline (no infection, no treatment). Lowly expressed genes were filtered out by removing genes with fewer than 10 total counts across all samples. Differential expression analysis was performed using DESeq2, followed by variance stabilizing transformation (VST) for downstream visualization and exploratory analysis.

# Plot analysis
1. **Scatter plot** 

The dispersion plot shows the expected strong negative relationship between mean normalized counts and dispersion estimates, which is characteristic of well-behaved RNA-seq data under the negative binomial model.
The genes with low mean expression levels (left side of the graph) show high dispersion, leading to high scatter of gene-specific estimates (black dots). With increasing mean expression levels, the dispersion reduces significantly, indicating low relative variability for highly expressed genes. The fitted dispersion curve (red line) captures this pattern well, passing through the middle of the cloud of points across the entire expression range. The final shrunk dispersion estimates (blue dots) lie very close to the fitted curve, indicating proper and effective shrinkage. This is especially helpful for low-count genes, for which the raw estimates are noisy. There do not appear to be any serious dispersion outliers, indicating the absence of significant technical artifacts, contamination, or genes with unexplained extreme variability. The dispersion estimation procedure appears to be trustworthy, and this helps build confidence in the subsequent differential expression analysis, such as p-values and log2 fold changes.  

![Scatter plot](plots/dispersion_plot.png) 


2. **PCA plot**
PC1 and PC2 explain 53% and 27% of the total variance, respectively. PC1 and PC2 together explain ~80% of the total variance, providing a highly representative 2D summary of the data.



3. **MA plots** 

The MA plot visualizes the relationship between mean normalized expression (baseMean, x-axis) and log2 fold change (y-axis) for the contrast Infected_Vehicle vs Uninfected_Vehicle, using shrunk log2 fold changes.

**MA Plot: Infected vs Uninfected (Vehicle)**

Key observations:
Most genes are tightly grouped around log2FC = 0 (the gray area), as expected. This indicates that most genes do not display significant changes in expression levels during infection, consistent with stable housekeeping and baseline activities.

- There is a distinct group of genes that display strong differential expression, with log2FC values extending from ~+7 (upregulated) to ~-4 (downregulated). These values are biologically consistent in the context of a strong innate immune response to Francisella tularensis infection in THP-1 cells. In this scenario, the activation of inflammatory and antiviral response pathways, such as NF-κB, cytokine, and interferon responses, may mediate strong fold changes in key effector genes.
- The lowly expressed genes (left side of the plot) display more variability in log2FC, which is expected because of greater relative variability in lowly expressed genes.
- The highly expressed genes (right side of the plot) display much more compact fold changes, as expected for stable housekeeping genes such as ribosomal genes, actin, or GAPDH, which are less responsive to infection.

The y-axis was scaled to [-6, 8] to display the fold changes without too much empty space, while maintaining sufficient detail across the expression profile. This is consistent with a strong biological response and helps validate the results.

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

**Volcano Plot: Infected (Infected + Vehicle vs Infected + Dillapiole)**


Key observations:
- Relatively balanced up and downregulation, with a wide y-range (~90) indicating very high statistical significance for changes.
- Downregulation is more dispersed and extensive (down to ~-7), especially in low-expression genes (left side).
- Upregulation is weaker (max ~+5) and more clustered near zero, focusing on high-expression genes (right side).
- High density of unaffected genes (gray band) compared to other contrasts —> fewer genes are significantly perturbed.

**Comparison with Vehicle (pure infection)**:
- Vehicle shows stronger upregulation (up to ~+7) and limited downregulation (max ~-4), with changes dispersed in low-expression genes.
- Infected shows weaker upregulation (max ~+5, shifted to high-expression genes) and stronger, more dispersed downregulation (max ~-7, focused on low-expression genes).
- Much higher density near log2FC = 0 in Infected, dillapiole reduces the number of affected genes during infection.
- Downregulation is more prominent and widespread in Infected, while upregulation is attenuated and redirected.

**Interpretation**:
Dillapiole during infection attenuates the strong activation seen in pure infection (Vehicle) and enhances downregulation, particularly in low-expression genes. This results in a more limited, balanced, and suppressed transcriptional response, consistent with dillapiole dampening bacterial virulence, reducing the intensity of the host immune response, and shifting changes toward high-expression genes for a more controlled outcome. 

![Volcano plot, Infected](https://github.com/niktaghanei/differential-gene-expression-analysis/blob/main/plots/volcano%20plots/infected_withsymbol_volcano.png)


5. **Heatmaps**

Heatmaps: Expression Patterns of Top DE Genes
Heatmaps are a powerful way to visualize gene expression patterns across samples for a selected set of genes (here, the top 10 up and downregulated genes from each contrast). 
They show:
- **Rows**: Individual genes (labeled with gene symbols for readability, sourced from org.Hs.eg.db).
- **Columns**: Individual samples (grouped and annotated by infection status and treatment).
- **Color scale**: VST-transformed expression values (blue = low expression, red = high expression; white = intermediate).

Clustering is applied to both rows (genes) and columns (samples) using hierarchical clustering with Euclidean distance. This helps reveal:
- How samples group together (infected vs uninfected).
- How genes co-express across conditions.
- Patterns of up or downregulation in specific groups.

Heatmaps complement volcano and MA plots by showing **relative expression levels** across all samples, not just fold changes or significance. They are especially useful for confirming biological consistency and spotting replicate variability or outliers.

**Heatmap: Vehicle (Top DE Genes)**

Key observations
- Uninfected samples cluster tightly
- Infected samples show slightly greater spread (expected biological variability in response to intracellular infection).
- Most top genes are upregulated in infected samples (lighter colors in infected columns), including key immune/stress genes:  
  - **CCL3-AS1** (antisense to chemokine CCL3 -> immune cell recruitment)  
  - **HSPA6** (heat shock protein -> cellular stress protection)  
  - **IL4I1** (immune regulation and arginine metabolism)  
  - **BIRC3** (apoptosis inhibitor -> cell survival)  
  - **IGFBP3** (IGF binding -> growth/apoptosis control)  
  - **MAFB** (transcription factor -> macrophage differentiation)  
  - **RAB20** (vesicle trafficking -> phagosome maturation)

Downregulated genes are fewer and less noticeable with some variability (sample Infected 4 behaves slightly differently).

**Comparison with other plots**:
- Compared to Uninfected (dillapiole baseline): Vehicle shows stronger and more focused upregulation (immune activation), while Uninfected has more balanced up/down with broader suppression.
- Compared to Dillapiole (infection in presence of dillapiole): Vehicle has more extreme upregulation and less downregulation, dillapiole appears to attenuate the strong activation and shift toward more downregulation.

**Interpretation**:
Pure infection drives a strong, upregulation-dominated response in THP-1 cells, with activation of immune recruitment (CCL3-AS1), stress protection (HSPA6), and survival pathways (BIRC3, IGFBP3). This is consistent with a strong innate immune reaction to intracellular *F. tularensis*. Variability in infected samples reflects biological heterogeneity but overall reproducibility is high. 

![Heatmap, Vehicle](https://github.com/niktaghanei/differential-gene-expression-analysis/blob/main/plots/heatmaps/vehicle_withsymbol.png)


**Heatmap: Dillapiole (Top DE Genes)**

Clustering and sample patterns:
Clustering reveals considerable spread among samples with large distances between groups. This is expected given the presence of two simultaneous factors (infection and treatment). Even replicates can exhibit variable responses to the combination and with only 5 samples per condition, some dispersion is normal and biologically reasonable.

Key bservations:
An unnamed gene (identified only by ENSG ID) shows the largest shift (expression levels around 5 in uninfected + dillapiole samples rise to 7–8 in infected + dillapiole samples). This gene appears exclusively in this contrast, limiting specific functional interpretation, but the consistent upregulation across most infected samples is notable.
ATP10B: Expression remains low (4–5) in uninfected + dillapiole samples but increases to 6 in 4 out of 5 infected + dillapiole samples, a reliable and consistent pattern. ATP10B functions as a phospholipid flippase, maintaining membrane asymmetry and facilitating lipid trafficking. The upregulation in infected cells likely reflects membrane remodeling or altered lipid handling in response to intracellular infection.

**Comparison with the corresponding volcano plot (Infected vs Uninfected in Dillapiole):**

NGFR: Appears upregulated in the volcano plot (log2FC ≈ 7), but in this heatmap most infected + dillapiole samples show lower relative expression (darker blue tones). This discrepancy arises because volcano plots reflect group-level averages and statistical significance, while heatmaps display per-sample expression values. Replicate heterogeneity or threshold effects can explain the difference. NGFR (neurotrophin receptor) is involved in immune modulation and cell survival; the observed variation likely reflects sample-specific responses.
NGFR-AS1: Upregulated in the uninfected + dillapiole/vehicle volcano (log2FC ≈ 7.6), but here 4 out of 5 infected + dillapiole samples show higher expression (lighter colors), with one sample decreasing. NGFR-AS1 is an antisense lncRNA that regulates NGFR expression. The increase in infected + dillapiole samples suggests the combination of drug and infection enhances this regulatory pathway, though sample variation is evident.

**Comparison with previous plots (Uninfected and Vehicle):**

LRRC41, WDR5, MLX, and VRK2 show consistent downregulation across multiple contrasts involving dillapiole (positions and magnitude remain similar). This pattern strongly suggests these genes are direct targets of dillapiole, independent of infection status. LRRC41 participates in ubiquitin ligase complexes (protein degradation), WDR5 regulates histone methylation (gene activation), MLX is a transcription factor involved in metabolic control, and VRK2 is a stress-response kinase which their suppression may reflect dillapiole's broader modulatory or suppressive influence on host pathways.
In contrast to the Vehicle plot (pure infection), where upregulation was dominant and downregulation limited, this heatmap shows a more balanced response with stronger and more dispersed downregulation, particularly in low-expression genes. This shift indicates dillapiole attenuates infection-driven activation while enhancing suppressive effects.

**Interpretation:**
The presence of dillapiole during infection results in a relatively balanced transcriptional response, with comparable numbers of up and downregulated genes but more dispersed and prominent downregulation, especially among low-expression genes. The overall number of significantly perturbed genes is reduced compared to pure infection (Vehicle), as evidenced by greater density near baseline expression levels. Genes such as ATP10B and NGFR-AS1 show consistent increases in infected + dillapiole samples, likely related to membrane dynamics and immune regulatory pathways. Shared downregulated genes (LRRC41, WDR5, MLX, VRK2) across dillapiole-containing contrasts point to direct drug effects, while the attenuated and redirected pattern during infection supports dillapiole's role in modulating host response probably by dampening bacterial virulence and leading to a more controlled transcriptional outcome.
Reproducibility is reasonable despite visible replicate variation and is typical in infection studies with multiple factors. 

![heatmap, Dillapiole](https://github.com/niktaghanei/differential-gene-expression-analysis/blob/main/plots/heatmaps/Dillapiole_withsymbol.png)
