# load required packages
library(DESeq2)
library(biomaRt)
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)
library(here)
library(tidyverse)
library(apeglm)
library(ashr)

setwd("C:/*/*/*/checking/")
# loading raw data
raw <- read.csv("C:/*/*/*/GSE306199_raw_counts_24hr.csv.gz")
countData <- raw[, -1]
rownames(countData) <- raw[, 1]
write.csv(countData, "C:/*/*/*/checking/countData")

# sample metadata
colData <- data.frame(
  row.names = colnames(countData),
  infection = factor(rep(c('Infected', 'Uninfected'), each = 10)),
  treatment = factor(rep(c('Dillapiole', 'Vehicle'), times = 2, each = 5)))

write.csv(colData, "C:/*/*/*/checking/sample_metadata.csv")

# testing if row names and column names(samples) match
all(colnames(countData) %in% rownames(colData))
all(colnames(countData) == rownames(colData))

#reference levels
colData <- colData %>%
  mutate(infection = relevel(infection, ref = 'Uninfected')) %>%
  mutate(treatment = relevel(treatment, ref = 'Vehicle')) %>%
  mutate(group = factor(paste(infection, treatment, sep = '_'),
         levels = c('Uninfected_Vehicle',
                    'Uninfected_Dillapiole',
                    'Infected_Vehicle',
                    'Infected_Dillapiole')))




# DESeq2 object and basic filtering
dds <- DESeqDataSetFromMatrix(countData = countData, 
                       colData = colData,
                       design = ~ group)

keeps <- rowSums(counts(dds)) >= 10
dds <- dds[keeps, ]


# the actual differential gene expression analysis
dds <- DESeq(dds)
vsd <- vst(dds, blind = FALSE)

# PCA loadings
pca <- prcomp(t(assay(vsd)))

pc1_loadings <- sort(pca$rotation[,1], decreasing = TRUE)
pc2_loadings <- sort(pca$rotation[,2], decreasing = TRUE)

top_pc1_pos <- head(pc1_loadings, 20)
top_pc1_neg <- head(pc1_loadings[order(pc1_loadings)], 20)

top_pc2_pos <- head(pc2_loadings, 20)
top_pc2_neg <- head(pc2_loadings[order(pc2_loadings)], 20)

write.csv(top_pc1_pos, 'PC1_top_positive_genes.csv')
write.csv(top_pc1_neg, 'PC1_top_negative_genes.csv')
write.csv(top_pc2_pos, 'PC2_top_positive_genes.csv')
write.csv(top_pc2_neg, 'PC2_top_negative_genes.csv')


  
# PCA plots 
pcaplot <- plotPCA(vsd, intgroup = c('infection', 'treatment')) +
  theme_minimal() +
  labs(title = 'PCA plot: Infection and Treatment') +
  theme(plot.title = element_text(hjust = 0.5))

print(pcaplot)
ggsave('PCA_plot.png', pcaplot,  width = 10, height = 7, dpi = 300)

# dispersion plot
png('dispersion_plot.png', width = 800, height = 600)
plotDispEsts(dds)
dev.off()

# library sizes
lib_sizes <- colSums(counts(dds, normalized = TRUE))
cat('Library size summary:\n')
print(summary(lib_sizes))

png('sample_QC_library_sizes.png', width = 1000, height = 600, res = 120)
par(mar = c(11, 5, 4, 2) + 0.1)
barplot(lib_sizes, las = 2, col = "steelblue", border = NA,
        main = 'Normalized library sizes', ylab = 'Total normalized counts')
grid(nx = NA, ny = NULL, col = "gray80", lty = 'dotted')
dev.off()



# Raw results
res_vehicle_raw <- results(dds, contrast=c('group','Infected_Vehicle','Uninfected_Vehicle'))
res_dillapiole_raw <- results(dds, contrast=c('group','Infected_Dillapiole','Uninfected_Dillapiole'))
res_uninf_raw <- results(dds, contrast=c('group','Uninfected_Dillapiole','Uninfected_Vehicle'))
res_inf_raw <- results(dds, contrast=c('group','Infected_Dillapiole','Infected_Vehicle'))

write.csv(as.data.frame(res_vehicle_raw),'infection_vehicle_all.csv')
write.csv(as.data.frame(res_dillapiole_raw),'infection_dillapiole_all.csv')
write.csv(as.data.frame(res_uninf_raw),'treatment_uninfected_all.csv')
write.csv(as.data.frame(res_inf_raw),'treatment_infected_all.csv')


# shrink
res_vehicle_shr <- lfcShrink(dds, coef = 3, type = 'apeglm')
res_dillapiole_shr <- lfcShrink(dds, coef = 4, type = 'apeglm')
res_uninf_shr <- lfcShrink(dds, coef = 2, type = 'apeglm')
res_inf_shr <- lfcShrink(dds, contrast = c('group', 'Infected_Dillapiole', 'Infected_Vehicle'), type = 'ashr')   


# Significant (shrink)
sig_vehicle <- subset(as.data.frame(res_vehicle_shr), padj < 0.05 & abs(log2FoldChange) > 1)
sig_dillapiole <- subset(as.data.frame(res_dillapiole_shr), padj < 0.05 & abs(log2FoldChange) > 1)
sig_uninf <- subset(as.data.frame(res_uninf_shr), padj < 0.05 & abs(log2FoldChange) > 1)
sig_inf <- subset(as.data.frame(res_inf_shr), padj < 0.05 & abs(log2FoldChange) > 1)

write.csv(sig_vehicle,'infection_vehicle_sig.csv')
write.csv(sig_dillapiole,'infection_dillapiole_sig.csv')
write.csv(sig_uninf,'treatment_uninfected_sig.csv')
write.csv(sig_inf,'treatment_infected_sig.csv')


# MA plots
png('infected_vs_uninfected_vehicle_MA.png', width = 800, height = 600)
plotMA(res_vehicle_shr, main='Vehicle')
dev.off()

png('infected_vs_uninfected_dillapiole_MA.png', width = 800, height = 600)
plotMA(res_dillapiole_shr, main='Dillapiole')
dev.off()

png('dillapiole_vs_vehicle_uninfected_MA.png', width = 800, height = 600)
plotMA(res_uninf_shr, main='Uninfected')
dev.off()

png('dillapiole_vs_vehicle_infected_MA.png', width = 800, height = 600)
plotMA(res_inf_shr, main='Infected')
dev.off()


# Volcano plots
res_vehicle_shr$symbol <- mapIds(org.Hs.eg.db, keys = rownames(res_vehicle_shr), keytype = "ENSEMBL", column = "SYMBOL")
res_dillapiole_shr$symbol <- mapIds(org.Hs.eg.db, keys = rownames(res_dillapiole_shr), keytype = "ENSEMBL", column = "SYMBOL")
res_uninf_shr$symbol <- mapIds(org.Hs.eg.db, keys = rownames(res_uninf_shr), keytype = "ENSEMBL", column = "SYMBOL")
res_inf_shr$symbol <- mapIds(org.Hs.eg.db, keys = rownames(res_inf_shr), keytype = "ENSEMBL", column = "SYMBOL")

EnhancedVolcano(res_vehicle_shr,
                lab = res_vehicle_shr$symbol,
                x = 'log2FoldChange', y = 'padj',
                title = 'Infected vs Uninfected (Vehicle)')
ggsave('vehicle_withsymbol_volcano.png', width = 9, height = 7)

EnhancedVolcano(res_dillapiole_shr,
                lab = res_dillapiole_shr$symbol,
                x = 'log2FoldChange', y = 'padj',
                title = 'Infected vs Uninfected (Dillapiole)')
ggsave('dillapiole_withsymbol_volcano.png', width = 9, height = 7)

EnhancedVolcano(res_uninf_shr,
                lab = res_uninf_shr$symbol,
                x = 'log2FoldChange', y = 'padj',
                title = 'Dillapiole vs Vehicle (Uninfected)')
ggsave('uninfected_withsymbol_volcano.png', width = 9, height = 7)

EnhancedVolcano(res_inf_shr,
                lab = res_inf_shr$symbol,
                x = 'log2FoldChange', y = 'padj',
                title = 'Dillapiole vs Vehicle (Infected)')
ggsave('infected_withsymbol_volcano.png', width = 9, height = 7)


# Heatmaps 
gene_symbols <- mapIds(org.Hs.eg.db, keys = rownames(vsd), keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first")
names(gene_symbols) <- all_genes  
ann <- as.data.frame(colData(vsd)[,c('infection','treatment')])

top10_heatmap <- function(sig_res, title, filename, keep_fun) {
  up <- head(rownames(sig_res[order(sig_res$log2FoldChange, decreasing = TRUE), ]), 10)
  down <- head(rownames(sig_res[order(sig_res$log2FoldChange), ]), 10)
  genes <- unique(c(up, down))
  
  
  gene_labels <- gene_symbols[genes]
  gene_labels[is.na(gene_labels)] <- names(gene_labels)[is.na(gene_labels)]  
  
  keep <- keep_fun(colData(vsd))
  mat <- assay(vsd)[genes, keep]
  
 
  pheatmap(
    mat,
    annotation_col = ann[keep, ],
    main = title,
    filename = filename,
    labels_row = gene_labels, 
    fontsize_row = 8,          
    fontsize_col = 8,
    cluster_rows = TRUE,
    cluster_cols = TRUE
  )
}


top10_heatmap(sig_vehicle, 'Vehicle', 'vehicle_withsymbol.png', function(cd) cd$treatment == 'Vehicle')
top10_heatmap(sig_dillapiole, 'Dillapiole', 'Dillapiole_withsymbol.png', function(cd) cd$treatment == 'Dillapiole')
top10_heatmap(sig_uninf, 'Uninfected', 'uninfected_withsymbol.png', function(cd) cd$infection == 'Uninfected')
top10_heatmap(sig_inf, 'Infected', 'infected_withsymbol.png', function(cd) cd$infection == 'Infected')



# Annotation
mart <- useEnsembl(biomart = "genes", 
                   dataset = "hsapiens_gene_ensembl",
                   host = "useast.ensembl.org")
add_names <- function(df){
  ens <- rownames(df)
  info <- getBM(c("ensembl_gene_id","external_gene_name","description"), 
                "ensembl_gene_id", ens, mart)
  merge(df, info, by.x="row.names", by.y="ensembl_gene_id")
}

top_up <- head(sig_vehicle[order(sig_vehicle$log2FoldChange,decreasing=TRUE), ], 10)
top_down <- head(sig_vehicle[order(sig_vehicle$log2FoldChange), ], 10)

write.csv(add_names(top_up),'infection_vehicle_top10_up.csv')
write.csv(add_names(top_down),'infection_vehicle_top10_down.csv')

