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

- Most genes are centered around log2FC = 0, as expected for the bulk of the transcriptome. The y-axis limit was set to [-8, 10] to encompass the observed extremes while keeping all data points visible.
- Many genes display strong differential expression, with log2FC values spread from about -8 to +10. The large spread of values indicates a strong transcriptional response to infection even in the presence of Dillapiole.
- Low-expressed genes display the largest spread of log2FC values, which is expected because of the greater relative variability of low-count transcripts.
- Highly expressed genes display much smaller fold changes, as expected for stable housekeeping genes that are largely unresponsive to infection and treatment.

In comparison to the Vehicle treatment, this plot shows:

- A larger spread of points overall.
- More genes with large |log2FC| > 1 compared to the Vehicle treatment.
- Fewer genes with log2FC close to 0.

This indicates that Dillapiole might cause a more widespread or perhaps slightly enhanced transcriptional response during infection, rather than inhibiting it. A biological explanation for this observation might be that Dillapiole suppresses the expression of bacterial virulence genes, leading to a less aggressive infection and a more widespread host response.  

![MA plot, Dillapiole](https://github.com/niktaghanei/differential-gene-expression-analysis/blob/main/plots/ma%20plots/dillapiole_MA_edited.png)



**MA Plot: Dillapiole vs Vehicle (Uninfected)**

This plot illustrates the direct baseline effect of Dillapiole treatment on uninfected THP-1 cells.

- There is a typical mean-variance relationship: high scatter at low expression and compact fold changes at high expression.
- The y-axis range [-10, 8] captures all observed changes, with some genes being highly regulated (fold changes up to ~1000×).
- The density at log2FC = 0 is moderate, indicating that Dillapiole significantly perturbs a substantial part of the transcriptome even in the absence of infection.

Comparison: Uninfected vs Vehicle MA Plots 
To distinguish between the direct effect of Dillapiole on host cells from its possible modulation of the infection response, I compared the treatment effect in infected cells to the direct treatment effect in uninfected cells.

Main differences:

- Downregulation is stronger in Uninfected: There are more genes and larger negative log2FC values in the Uninfected plot, indicating that Dillapiole alone downregulates a larger part of the transcriptome in uninfected cells.
- Fewer genes unaffected in Uninfected: There is less density at log2FC = 0 in the Uninfected plot than in the Vehicle plot, indicating that Dillapiole significantly perturbs more genes even without infection.
- Upregulation is more spread out in Vehicle: Infected cells without treatment have a more spread-out distribution of positive log2FC values, which reflects the strong induction of immune and inflammatory responses.
- Downregulation is more concentrated around zero in Vehicle: There is less and weaker downregulation in pure infection, which is expected for innate immune stimulation that primarily targets upregulation rather than downregulation.

**Interpretation:**
Dillapiole has a predominantly suppressive role on uninfected THP-1 cells, whereas infection promotes widespread upregulation and minimal downregulation. These opposing profiles indicate that the presence of Dillapiole during infection is likely to modulate the host transcriptome in a manner that differs from its suppressive role.

**Comparison: Uninfected vs Infected MA Plots**

key differences:
- There is a substantially higher density of unaffected genes in the Infected plot, indicating that fewer genes overall exhibit significant changes in the "infection + Dillapiole" scenario.
- Upregulation in Uninfected is more extreme (up to ~+7) with greater spread, particularly for low-expression genes. In the Infected group, upregulation is less extreme (up to ~+5) and more concentrated on high-expression genes.
- Downregulation in the Uninfected group is stronger in high-expression genes, whereas in the Infected group, it tilts towards low-expression genes with less overall spread.

**Interpretation:**
Uninfected cells show overall changes due to Dillapiole alone. But in the infected state, the effect is more contained: fewer genes are affected, and upregulation is lessened. This indicates that Dillapiole has a modulating effect on the host response during infection, which is likely to cause a less widespread and balanced transcriptional response compared to its original effect on uninfected cells.

![MA Plot, Uninfected](https://github.com/niktaghanei/differential-gene-expression-analysis/blob/main/plots/ma%20plots/uninfected_MA_edited.png) 


**MA Plot: Infected (Infected + Vehicle vs Infected + Dillapiole)**

This MA plot shows the direct effect of Dillapiole treatment on infected THP-1 cells.
Overall observations:

- The range of the y-axis [-8, 5] encompasses all the data points, with the most extreme points reaching ~−7 and ~+5.
- The gray shading is broad and packed, showing that most genes are not affected.
- Downregulated genes are much more spread out in the plot than upregulated genes.
- Upregulated genes are more concentrated around the origin and are less spread out.
- Low-expression genes are mostly downregulated, whereas high-expression genes are mostly upregulated.
  
This indicates that the presence of Dillapiole during infection results in a more suppressed and asymmetric response: stronger downregulation (especially in low-expression genes), while upregulation is weaker and more focused on high-expression genes.

**Comparison with Vehicle Plot (Vehicle + Infected vs Vehicle + Uninfected) in pure infection (without treatment):**
- Upregulation is wider (up to ~+7) with greater dispersion.
- Downregulation is limited (max ~−4) and more clustered around zero.

Important differences from the Treated group:
- Density around zero: This attribute is much higher in the Infected group, indicating that fewer genes are significantly affected when Dillapiole is present.
- Upregulation: This feature is stronger in the Vehicle group but weaker and more clustered when treated.
- Downregulation: This is limited in the Vehicle group, but stronger and more dispersed in the Treated group.

**Interpretation:**
Pure infection induces a strong, activation-dominant response. The presence of Dillapiole during infection weakens this upregulation, strengthens downregulation, and reduces the number of affected genes. This leads to a more suppressed and less widespread transcriptional response.

![MA plot, Infected](https://github.com/niktaghanei/differential-gene-expression-analysis/blob/main/plots/ma%20plots/infected_MA_edited.png)



4. **Volcano plots**

Volcano Plots: Differential Expression Results 

Volcano plots are one of the most popular and effective methods for visualizing RNA-seq differential expression analysis results. They provide two crucial pieces of information in a single scatter plot:

#1 **Horizontal axis (x)**: log2 Fold Change (log2FC)  
  Indicates the degree and direction of gene expression change between two conditions.  
  - Positive values (right side): gene is upregulated (higher expression in the first group)
  - Negative values (left side): gene is downregulated (lower expression in the first group)  
  - A value of log2FC = 1 indicates ~2-fold upregulation; log2FC = -1 indicates ~2-fold downregulation

#2 **Vertical axis (y)**: -log10(adjusted p-value) or -log10(padj)  
  Indicates the statistical significance (how likely the observed value is due to random chance) 
  - Higher points = more significant (smaller padj)

**Common cut-offs** used in most RNA-seq analyses:
- padj < 0.05 (or 0.01 for more stringent control of false positives)
- |log2FC| > 1 (fold change > 2 or < 0.5, biologically significant)

Genes that meet both criteria (significant and large change) are typically colored red. Genes with large change but no significance are green, genes with significance but small change are light blue, and non-significant/small change genes are gray.


**Volcano Plot: Infected vs Uninfected (Vehicle)**

This volcano plot displays shrunk log2 fold changes and -log10(padj) for the contrast Infected_Vehicle vs Uninfected_Vehicle (pure infection effect).

**Key observations**:
- Asymmetric response: upregulation is more prominent (right side), with most genes showing high significance (red dots up to -log10(padj) ~60).
- Downregulation is minimal (few red dots on the left).
- Identified genes include prominent immune/stress response markers:  
  - **NFKBIA**: NF-κB inhibitor; regulates feedback of inflammation.  
  - **IL4I1**: Immune regulation and arginine metabolism; modulates inflammation.  
  - **HSPA1B**: Heat shock protein; protects cells under infection stress.  
  - **ERCC1**: DNA repair; responds to infection-induced damage.   
  - **TBC1D17**: Regulates autophagy/vesicle trafficking; putative defense mechanism against intracellular pathogen.

Majority of highly significant DEGs are upregulated and involved in immune activation, inflammation regulation, cellular stress response, and DNA repair processes that reflect a strong innate immune response to intracellular *Francisella tularensis* infection in THP-1 cells.
Scatter is higher in low-expression genes, but significant expression is more prominent in medium to high-expression genes (right side clustering).

No technical artifact is apparent; shrinkage improves the precision of low-expression genes. This plot emphasizes the strong activation role of infection in host gene expression.

![Volcano plot, Vehicle](https://github.com/niktaghanei/differential-gene-expression-analysis/blob/main/plots/volcano%20plots/vehicle_withsymbol_volcano.png)


**Volcano Plot: Infected vs Uninfected (Dillapiole)**


Key observations:
- Symmetric up- and downregulation: unlike the Vehicle plot (more upregulation), this plot shows more balanced up- and downregulation, with similar numbers of genes.
- Highly significant: -log10(padj) ~100 (well above Vehicle ~60), highly significant changes.
- Upregulation: reaches up to ~+10, with important genes such as **BEST1** (calcium-activated chloride channel, ion homeostasis), **KIF1B** (kinesin motor, vesicle/mitochondria transport), **MSC-AS1** (lncRNA, gene regulation), **ICAM4-AS1** (lncRNA, cell adhesion), **TRPA1** (cation channel, inflammation/pain sensing), **MRNIP** (DNA repair/stress response).
- Downregulation: reaches down to ~-8, with genes such as **ENOSF1** (fucose metabolism), **LRRC41** (ubiquitin ligase component), **WDR5** (histone methylation), **LCORL** (transcription factor).
- Green area (large FC but not significant) is less dense and spread out compared to Vehicle -> fewer large non-significant changes.

**Interpretation**:
Infection with Dillapiole has a more balanced and highly significant transcriptional profile compared to pure infection (Vehicle). Upregulation is less extreme, but downregulation is more pronounced and spread out. This indicates that Dillapiole has a modulating effect on the host response, suppressing extreme upregulation and promoting more downregulation, likely by inhibiting bacterial virulence and infection dynamics.

![Volcano plot, Dillapiole](https://github.com/niktaghanei/differential-gene-expression-analysis/blob/main/plots/volcano%20plots/dillapiole_withsymbol_volcano.png)



**Volcano Plot: Uninfected (Uninfected + Dillapiole vs Uninfected + Vehicle)**

Key observations:
- More balanced up- and downregulation: similar numbers of significant genes in both categories (no longer dominated by upregulation as in Vehicle).
- Highly significant: -log10(padj) ~80, strong statistical evidence for changes.
- Upregulation: extends up to ~+10, with major genes such as **KIF1B** (vesicle/mitochondria transport), **HDGFL3** (cell growth/repair), **TRPA1** (cation channel, inflammation sensing), **WDR13** (gene regulation/signaling).
- Downregulation: extends down to ~-8, with major genes such as **INTS8** (snRNA processing), **LRRC41** (ubiquitin ligase component), **WDR5** (histone methylation), **LCORL** (transcription factor).
- Green region (large FC but non-significant) is less dense and more scattered than in Vehicle → fewer large non-significant changes.

**Comparison with Vehicle (pure infection)**:
- Vehicle has more upregulation (max +7) and less downregulation (max -4), with low-expression genes dominating most changes.
- Uninfected has more balanced changes, more downregulation, and high-expression genes are more affected.
- Fewer unaffected genes (gray band) in Uninfected, Dillapiole affects more of the transcriptome in uninfected cells.
- Common downregulated genes (e.g., LRRC41, WDR5, MLX, VRK2) retain similar positions, indicating these are direct Dillapiole effects not influenced by infection.

**Interpretation**:
Uninfected cells have a balanced response to Dillapiole alone, with strong downregulation (especially high-expression genes) and upregulation in low-expression genes. This is in contrast to the activation-biased response in pure infection. The baseline downregulation by Dillapiole could explain its modulatory action in infection, where the response is more regulated and less drastic. 

![Volcano plot, Uninfected](https://github.com/niktaghanei/differential-gene-expression-analysis/blob/main/plots/volcano%20plots/uninfected_withsymbol_volcano.png) 

**Volcano Plot: Infected (Infected + Vehicle vs Infected + Dillapiole)**

Key observations:
- Relatively balanced down- and upregulation, with large y-range (~90) reflecting very high statistical significance of the changes.
- Downregulation is more scattered and widespread (down to ~-7), particularly in low-expression genes (left side).
- Upregulation is weaker (max ~+5) and more concentrated around zero, particularly in high-expression genes (right side).
- High density of unaffected genes (gray band) compared to other comparisons --> fewer genes are significantly affected.

**Comparison with Vehicle (pure infection)**:
- Vehicle has stronger upregulation (up to ~+7) and limited downregulation (max ~-4), with scattered changes in low-expression genes.
- Infected has weaker upregulation (max ~+5, shifted to high-expression genes) and stronger, more scattered downregulation (max ~-7, in low-expression genes).
- Much higher density near log2FC = 0 in Infected, dillapiole reduces the number of affected genes during infection.
- Downregulation is more pronounced and widespread in Infected, while upregulation is weakened and shifted.

**Interpretation**:
Dillapiole during infection weakens the strong activation present in pure infection (Vehicle) and further promotes downregulation, particularly in low-expression genes. This leads to a more restricted, balanced, and suppressed transcriptional response, in agreement with dillapiole weakening bacterial virulence, reducing the strength of the host immune response, and shifting the changes to high-expression genes for a more regulated response.

![Volcano plot, Infected](https://github.com/niktaghanei/differential-gene-expression-analysis/blob/main/plots/volcano%20plots/infected_withsymbol_volcano.png)


5. **Heatmaps**

Heatmaps: Expression Patterns of Top DE Genes
Heatmaps are an excellent tool for visualizing gene expression patterns across samples for a subset of genes of interest (in this case, the top 10 up and downregulated genes for each contrast).
They display:
- **Rows**: Individual genes (displayed with gene symbols for convenience, obtained from org.Hs.eg.db).
- **Columns**: Individual samples (aggregated and labeled according to infection status and treatment).
- **Color scale**: VST-normalized expression values (blue = low expression, red = high expression; white = intermediate).

Row and column clustering (for genes and samples, respectively) are performed using hierarchical Euclidean clustering. This aids in:
- Identifying samples that group together (infected vs uninfected).
- Identifying genes that co-regulate across conditions.
- Identifying genes that are up- or downregulated in specific groups.

Heatmaps can be used in addition to volcano and MA plots to visualize **expression levels** across all samples, rather than just fold changes or significance values. They are particularly useful for checking biological consistency and identifying replicate variability or outliers.


  
**Heatmap: Vehicle (Top DE Genes)**

Key observations
- Uninfected samples cluster tightly
- - Infected samples cluster slightly more variably (reflecting expected biological variation due to intracellular infection).
- Most top genes are upregulated in infected samples (lighter colors in infected columns), including prominent immune/stress response genes:
  - **CCL3-AS1** (antisense to chemokine CCL3 -> immune cell recruitment)  
  - **HSPA6** (heat shock protein -> cellular stress protection)  
  - **IL4I1** (immune regulation and arginine metabolism)  
  - **BIRC3** (apoptosis inhibitor -> cell survival)  
  - **IGFBP3** (IGF binding -> growth/apoptosis control)  
  - **MAFB** (transcription factor -> macrophage differentiation)  
  - **RAB20** (vesicle trafficking -> phagosome maturation)
The downregulated genes are fewer and less pronounced with some variation (sample Infected 4 is slightly different).

**Comparison with other plots**:
- Relative to Uninfected (dillapiole baseline): Vehicle is more strongly and specifically upregulated (immune activation), while Uninfected has more balanced up/down regulation with wider suppression.
- Relative to Dillapiole (infection in presence of dillapiole): Vehicle is more strongly upregulated and downregulates less, dillapiole reduces strong activation and increases downregulation.

**Interpretation**:
Pure infection induces a strong, upregulation-predominant response in THP-1 cells, activating immune recruitment (CCL3-AS1), stress protection (HSPA6), and survival mechanisms (BIRC3, IGFBP3). This is as expected for a strong innate immune response to intracellular *F. tularensis*. Variability among infected samples reflects biological variability, but overall, there is high reproducibility. 

![Heatmap, Vehicle](https://github.com/niktaghanei/differential-gene-expression-analysis/blob/main/plots/heatmaps/vehicle_withsymbol.png)


**Heatmap: Dillapiole (Top DE Genes)**

This heatmap illustrates the top 10 genes that are up- and down-regulated in the comparison Infected + Dillapiole vs Uninfected + Dillapiole (the effect of infection when dillapiole is present).
Clustering and sample patterns:
The samples are relatively spread out in the clustering, with big distances between groups. This is what we should see — we have infection and drug as factors, and we're looking at them at the same time. Even with 5 replicates per group, some variation is to be expected when two different factors are being considered.
Key observations:

There is one gene (no name, only ENSG ID) that shows the most change: from about 5 in uninfected + dillapiole samples, to 7-8 in infected + dillapiole samples. This gene is only present in this comparison, so we can't say exactly what it does, but it's certainly interesting that it's upregulated in infected samples.
ATP10B: This gene is low in uninfected + dillapiole samples (4-5), but is upregulated to 6 in 4 out of 5 infected + dillapiole samples. This is a strong and consistent result. ATP10B is a phospholipid flippase that participates in membrane asymmetry and lipid transport. This upregulation is likely a result of membrane remodeling during infection.

**Comparison with the corresponding volcano plot (Infected vs Uninfected in Dillapiole):**
NGFR: Looks like it is upregulated in the volcano (log2FC ≈ 7), but in the heatmap, most of the infected + dillapiole samples have lower expression (darker blue). Note that volcano plots show group averages, while heatmaps show individual samples – this is likely a difference due to replicate variation.
NGFR-AS1: Upregulated in the uninfected + dillapiole/vehicle volcano (log2FC ≈ 7.6), but here 4 out of 5 infected + dillapiole samples have higher expression, with one sample decreasing. NGFR-AS1 is an antisense lncRNA that regulates NGFR. The fact that most infected + dillapiole samples have higher expression suggests that the drug + infection combination together increase this regulatory pathway.

Comparison with previous plots:
LRRC41, WDR5, MLX, and VRK2 are all downregulated in all contrasts involving dillapiole (similar positions and fold changes). This is a very strong indication that these genes are direct targets of dillapiole, regardless of infection status. LRRC41 is a protein degradation factor, WDR5 is a histone methylation enzyme, MLX is a transcription factor involved in metabolic regulation, and VRK2 is a stress response protein kinase. Suppression of these genes likely indicates the general regulatory role of dillapiole on host pathways.
Relative to the Vehicle plot (pure infection), this heatmap reflects a more balanced response with greater and more dispersed downregulation, particularly in low-expression genes. This suggests that dillapiole suppresses the strong activation of pure infection and favors a more suppressive

**Interpretation:**
The transcriptional response in the presence of dillapiole is a fairly balanced one, with roughly as many up-regulated and down-regulated genes, although down-regulation is more widespread, especially in low-expression genes. There are fewer genes significantly affected than in pure infection (Vehicle), as indicated by the higher density of points around the baseline expression level. Genes such as ATP10B and NGFR-AS1 are steadily up-regulated in infected + dillapiole samples, probably reflecting membrane and immune functions. The common down-regulated genes (LRRC41, WDR5, MLX, VRK2) in dillapiole-containing contrasts indicate a direct drug action, while the more carefully controlled pattern in infection suggests dillapiole’s function in modulating the host response, probably by lowering bacterial virulence.
The reproducibility is acceptable, although there is some variation between replicates, which is to be expected in infection experiments involving several factors.

![heatmap, Dillapiole](https://github.com/niktaghanei/differential-gene-expression-analysis/blob/main/plots/heatmaps/Dillapiole_withsymbol.png)

**Heatmap: Uninfected (Uninfected + Dillapiole vs Uninfected + Vehicle)**

This heatmap illustrates the VST-normalized expression of the top 10 up- and down-regulated genes from the contrast Uninfected + Dillapiole vs Uninfected + Vehicle (direct effect of dillapiole on healthy cells).
Clustering and sample distribution:
The samples are not grouped very tightly; there is some variation between the samples. This is expected because the only variable in this experiment is dillapiole treatment, and even identical samples may react in different ways to the drug. The samples treated with dillapiole are grouped a bit more tightly than the others, indicating a reaction to the drug itself.

Key gene patterns:
One unnamed gene (only ENSG ID) has the largest difference expression of about 4.5-5 in the untreated samples (all dark blue), but much higher (6-7) in the treated samples (lighter colors in all five). This gene does not come up in any other comparisons as a top hit, so it is specific to dillapiole treatment alone. The fact that it is higher in all treated samples shows that it is a direct and reliable response to the drug, probably a gene strongly induced by dillapiole in healthy cells.
- LOC389602: Very low expression in untreated samples (4-5, dark blue), but much higher (6-7) in treated samples (lighter colors in all five). This is consistent and reliable. LOC389602 is a very poorly characterized locus, but this kind of difference usually shows drug-induced transcriptional activation in non-infected cells.
- NGFR-AS1: Much higher expression in treated samples (lighter colors) than in untreated samples (darker blue, around 4). This is consistent with its upregulation in the uninfected + dillapiole/vehicle volcano (log2FC ≈ 7.6). NGFR-AS1 is an antisense lncRNA that regulates NGFR expression. The strong induction by dillapiole indicates that the drug directly enhances this regulatory pathway in healthy cells.
- LINC00543: Low in untreated samples (~4.5, dark blue), but definitely higher in treated samples (lighter colors). This is a clear pattern. LINC00543 is a long non-coding RNA with gene regulatory functions (often in stress or metabolic pathways). The elevation suggests direct dillapiole induction in healthy cells.
MIR371A, MAFB, RPSAP44, VWA5A, NGFR, CPED1: All display the same pattern: low expression in untreated samples (~4, dark blue), but definitely higher in treated samples (lighter colors, ~5-6). 
These genes have roles in the following processes:
- MIR371A: microRNA cluster - regulates pluripotency and cell fate.
- MAFB: transcription factor - regulates macrophage differentiation and inflammation.
- RPSAP44: ribosomal protein pseudogene - may have regulatory function.
- VWA5A: von Willebrand factor A domain - mediates cell adhesion.
- NGFR: neurotrophin receptor - immune regulation and cell survival.
- CPED1: cadherin-like and PC-esterase domain - mediates cell adhesion and signaling.
The same pattern of increased expression in treated samples suggests that dillapiole directly induces these genes in healthy cells, perhaps as part of a stress response.
**TP73-AS2, PLA2R1, RRM2, SPP1, LOC100128548, LOC124901073, RPL7P39:**
 All have the reverse pattern: higher expression in untreated samples (~5-6), but strongly downregulated in treated samples (dark blue, around 4). These genes are involved in:
TP73-AS2: antisense to TP73 –> regulation of apoptosis and stress response.
PLA2R1: phospholipase A2 receptor –> lipid signaling and inflammation.
RRM2: ribonucleotide reductase –> DNA synthesis and cell proliferation.
SPP1: osteopontin –> immune cell recruitment and inflammation.
LOC100128548, LOC124901073: poorly characterized loci –> likely non-coding or pseudogenes.
RPL7P39: ribosomal protein pseudogene –> possible regulatory function.
The strong downregulation in treated samples suggests that dillapiole strongly represses these genes in healthy cells possibly to inhibit proliferation, inflammation, or lipid signaling pathways.


**Comparison with previous plots:**

- Most of the downregulated genes in this plot (LRRC41, WDR5, PLA2R1, and others) also appeared downregulated in other dillapiole-containing comparisons this further supports that these are direct drug effects, and not influenced by infection.
- By contrast, in the Vehicle heatmap (pure infection), the predominant effect was upregulation, while downregulation was less prominent. However, in this heatmap, there is a more pronounced and consistent suppression of treated samples where dillapiole has a wide-ranging inhibitory effect on healthy cells.
- The trend for NGFR-AS1 and LINC00543 (upregulation with dillapiole) is consistent across uninfected contrasts and is probably direct drug induction.

**Interpretation:**
Dillapiole alone induces a strong and consistent transcriptional effect in uninfected THP-1 cells, with strong upregulation of some genes (ATP10B, NGFR-AS1, LINC00543, etc.) and strong downregulation of others (TP73-AS2, PLA2R1, RRM2, SPP1, etc.). This suggests that dillapiole has a wide-ranging suppressive effect on many pathways (proliferation, inflammation, lipid metabolism) while inducing others (possibly stress or regulatory pathways). The effect is more balanced than the Vehicle (pure infection), but suppression is more pronounced. Reproducibility is good, but sample variability is evident as expected when a single variable (the drug) is introduced.

![heatmap, uninfected](https://github.com/niktaghanei/differential-gene-expression-analysis/blob/main/plots/heatmaps/uninfected_withsymbol.png)

**Heatmap: Infected (Infected + Vehicle vs Infected + Dillapiole)**

This heatmap illustrates VST-normalized expression values for the top 10 up- and down-regulated genes from the comparison Infected + Vehicle vs Infected + Dillapiole (direct effect of dillapiole on infected cells).

**Clustering and sample patterns:**
The samples treated with dillapiole are more closely clustered than the untreated samples. This is as expected, since dillapiole is a major player in this comparison, and the treated samples are expected to be more similar to each other, even if they are infected. The untreated infected samples are more dispersed, which is also as expected, since infection alone can lead to more variability between cells (different rates of bacterial uptake, replication, or immune activation).

Key gene patterns:
- TECR: High expression in untreated infected samples (around 9), but lower in treated infected samples (around 5). This is a very clear and consistent pattern. In the Vehicle volcano (Infected vs Uninfected), TECR was upregulated (log2FC ≈ 5, adjp-value ~47). In this contrast (Infected + Vehicle vs Infected + Dillapiole), it is strongly downregulated (log2FC ≈ -6, adjp-value ~47). TECR is an enzyme involved in very long-chain fatty acid biosynthesis (trans-2,3-enoyl-CoA reductase). The high expression level in pure infection is probably a reflection of bacterial manipulation of lipid metabolism for survival or replication. The sharp decrease with dillapiole indicates that the drug inhibits this pathway, possibly by lowering bacterial virulence, so that the host does not need to activate lipid metabolism as strongly.
  
- TDRD6: Very low expression in untreated infected samples (about 6), but sharply elevated to 9-10 in treated infected samples (lighter colors across most). This is consistent and reliable. TDRD6 is a tudor domain-containing protein that is involved in piRNA biogenesis and spermatogenesis, but in immune cells, it could have roles in RNA silencing or stress response. The sharp induction with dillapiole in infected cells suggests that the drug induces this gene, perhaps as part of a protective or regulatory response when infection is present.
  
- KIT: Relatively high expression in untreated infected samples (6-7), but reduced to about 4 in treated infected samples. In the uninfected + dillapiole/vehicle volcano, KIT was downregulated (log2FC ≈ -3). KIT is a receptor tyrosine kinase that is involved in cell survival, proliferation, and differentiation (particularly in hematopoietic cells).
  
- CPED1: Low expression in untreated infected samples (about 4), but higher in treated infected samples (about 5.5). This is consistent with the Uninfected heatmap (4 in untreated, 5.5 in treated). In the uninfected volcano plot, it was upregulated (log2FC ≈ 5, adjp-value ~5.7). CPED1 has cadherin-like and PC-esterase domains, and it is involved in cell adhesion and signaling. The fact that it is consistently upregulated by dillapiole (in both uninfected and infected cells) indicates a direct drug effect, which is probably related to increased adhesion or signaling in response to the drug.

**THAP12P3, LINC00543, CCDC24, B4GALT2: These genes all show the same pattern:**
very low expression in untreated infected samples (~4), but much higher in treated infected samples (~5.5-6.6).
These genes are:
- THAP12P3: Pseudogene or non-coding RNA; role unknown, but pseudogenes can have regulatory roles.
- LINC00543: Long non-coding RNA; involved in gene regulation (often stress or metabolic responses).
- CCDC24: Coiled-coil domain protein; may be involved in protein interactions or protein structure.
- B4GALT2: Beta-1,4-galactosyltransferase; involved in glycosylation and cell surface modification.
The fact that these genes are consistently upregulated by dillapiole in infected cells indicates a direct drug effect, which is probably related to a stress or membrane modification response.

**Interpretation:**
The transcriptional profile of dillapiole-treated infected cells is a more balanced response than that of infected cells alone (Vehicle). The tighter clustering of treated samples reflects a more similar response to the drug, even in the presence of infection. The downregulation of genes such as TECR (lipid metabolism) and KIT (proliferation/survival) indicates that dillapiole inhibits processes that the bacteria may utilize during infection. The upregulation of genes such as TDRD6 (RNA silencing/stress response), CPED1 (adhesion/signaling), and others (THAP12P3, LINC00543, CCDC24, B4GALT2) indicates direct induction by the drug, likely a regulatory or protective response. In general, dillapiole seems to modulate the activation of infection while promoting a more controlled and suppressive host response.
The reproducibility is good, with expected but visible variation among replicates in this multi-factorial environment.

![heatmap, infected](https://github.com/niktaghanei/differential-gene-expression-analysis/blob/main/plots/heatmaps/infected_withsymbol.png)


## Conclusion
This study focused on the impact of Francisella tularensis infection on gene expression in THP-1 cells and whether dillapiole can modulate this process.
Pure infection (Vehicle) caused a strong, one-sided response: extremely strong induction of immune, inflammatory, and stress response genes (e.g., NFKBIA, IL4I1, HSPA6, BIRC3), with very little downregulation. This is the typical, strong innate immune response expected for an intracellular bacterial pathogen.
The effect of dillapiole alone (Uninfected) was more balanced but suppressive on normal cells, downregulating many genes (particularly highly expressed ones) and upregulating others (e.g., KIF1B, NGFR-AS1, LINC00543). This indicates that the drug also has a direct effect on the host gene expression profile, irrespective of infection.

When combined (Infected + Dillapiole), the profile was significantly different:
- The strong induction observed in the vehicle (pure infection) was clearly dampened (less strong upregulation, less extreme changes).
- Downregulation was stronger and more extensive, particularly in low-expression genes.
- Fewer genes were significantly altered, leading to a more controlled and balanced transcriptional response.

These findings indicate that dillapiole has a positive effect on F. tularensis infection. It does not only contribute its own effects but also modulates the infection by reducing the bacterial virulence, making it less aggressive and more manageable for the host. Rather than causing a massive immune response, the cells react in a more composed manner, which is consistent with the drug’s mechanism of suppressing the key virulence genes of the bacterium.
In conclusion, this analysis proves that dillapiole is a virulence modulator and not a direct killer or suppressor of the host immune system. The project was able to successfully integrate computational analysis and biological interpretation to give a clear understanding of how this natural compound affects the host-pathogen interaction at the transcriptional level.
