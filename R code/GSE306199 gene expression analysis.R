# load required packages
library(DESeq2)
library(biomaRt)
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)
library(here)

  
# loading raw data
raw <- read.csv("C:/***/***/***/gene expression analysis/GSE306199_raw_counts_24hr.csv.gz")
countData <- raw[, -1]
rownames(countData) <- raw[, 1]
write.csv(countData, "C:/***/***/***/gene-expression-analysis/data/count_matrix.csv")

# sample metadata
colData <- data.frame(
  row.names = colnames(countData),
  infection = factor(rep(c('Infected', 'Uninfected'), each = 10)),
  treatment = factor(rep(c('Dillapiole', 'Vehicle'), times = 2, each = 5))
)
write.csv(colData, "C:/***/***/***/gene-expression-analysis/data/sample_metadata.csv")

# testing if row names and column names(samples) match
all(colnames(countData) %in% rownames(colData))
all(colnames(countData) == rownames(colData))

#reference levels
colData$infection <- relevel(colData$infection, ref = 'Uninfected')
colData$treatment <- relevel(colData$treatment, ref = 'Vehicle')
colData$group <- factor(paste(colData$infection, colData$treatment, sep = '_'))

  
# DESeq2 object and basic filtering
dds <- DESeqDataSetFromMatrix(countData = countData, 
                       colData = colData,
                       design = ~ group)


keeps <- rowSums(counts(dds)) >= 10
dds <- dds[keeps, ]

# the actual differential gene expression analysis
dds <- DESeq(dds)
vsd <- vst(dds, blind = FALSE)
ann <- as.data.frame(colData(vsd)[, c('infection','treatment')])



  
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



# differential expression comparisons
res_vehicle <- results(dds, contrast = c('group', 'Infected_Vehicle', 'Uninfected_Vehicle'))
res_dillapiole <- results(dds, contrast = c('group', 'Infected_Dillapiole', 'Uninfected_Dillapiole'))
res_uninf <- results(dds, contrast = c('group', 'Uninfected_Dillapiole', 'Uninfected_Vehicle'))
res_inf <- results(dds, contrast = c('group', 'Infected_Dillapiole', 'Infected_Vehicle'))

write.csv(as.data.frame(res_vehicle), 'infection_vehicle_all.csv', row.names = TRUE)
write.csv(as.data.frame(res_dillapiole), 'infection_dillapiole_all.csv', row.names = TRUE)
write.csv(as.data.frame(res_uninf), 'treatment_uninfected_all.csv', row.names = TRUE)
write.csv(as.data.frame(res_inf), 'treatment_infected_all.csv', row.names = TRUE)


# significant changes
sig_vehicle <- subset(as.data.frame(res_vehicle), padj < 0.05 & abs(log2FoldChange) > 1)
sig_dillapiole    <- subset(as.data.frame(res_dillapiole),    padj < 0.05 & abs(log2FoldChange) > 1)
sig_uninf   <- subset(as.data.frame(res_uninf),   padj < 0.05 & abs(log2FoldChange) > 1)
sig_inf     <- subset(as.data.frame(res_inf),     padj < 0.05 & abs(log2FoldChange) > 1)

write.csv(sig_vehicle, 'infection_vehicle_sig.csv', row.names = TRUE)
write.csv(sig_dillapiole, 'infection_dillapiole_sig.csv', row.names = TRUE)
write.csv(sig_uninf, 'treatment_uninfected_sig.csv', row.names = TRUE)
write.csv(sig_inf, 'treatment_infected_sig.csv' , row.names = TRUE)




# MA plots 
plotMA(res_vehicle, 'MA Plot: Infected vs Uninfected (Vehicle)', 'infected_vs_uninfected_vehicle_MA.png')
plotMA(res_dillapiole, 'MA Plot: Infected vs Uninfected (Dillapiole)', 'infected_vs_uninfected_dillapiole_MA.png')
plotMA(res_uninf, 'MA Plot: Dillapiole vs Vehicle (Uninfected)', 'dillapiole_vs_vehicle_uninfected_MA.png')
plotMA(res_inf, 'MA Plot: Dillapiole vs Vehicle (Infected)', 'dillapiole_vs_vehicle_infected_MA.png')




# volcano plots
EnhancedVolcano(res_vehicle,
                lab = rownames(res_vehicle),
                x = 'log2FoldChange',
                y = 'padj',
                title = 'Infected vs Uninfected (Vehicle)',
                pCutoff = 0.01,
                FCcutoff = 1)
ggsave('infected_vs_uninfected_vehicle_volcano.png', width = 9, height = 7)

EnhancedVolcano(res_dillapiole,
                lab = rownames(res_dillapiole),
                x = 'log2FoldChange',
                y = 'padj',
                title = 'Infected vs Uninfected (Dillapiole)',
                pCutoff = 0.01,
                FCcutoff = 1)
ggsave('infected_vs_uninfected_dillapiole_volcano.png', width = 9, height = 7)

EnhancedVolcano(res_uninf,
                lab = rownames(res_uninf),
                x = 'log2FoldChange',
                y = 'padj',
                title = 'Dillapiole vs Vehicle (Uninfected)',
                pCutoff = 0.01,
                FCcutoff = 1)
ggsave('dillapiole_vs_vehicle_uninfected_volcano.png', width = 9, height = 7)

EnhancedVolcano(res_inf,
                lab = rownames(res_inf),
                x = 'log2FoldChange',
                y = 'padj',
                title = 'Dillapiole vs Vehicle (Infected)',
                pCutoff = 0.01,
                FCcutoff = 1)
ggsave('dillapiole_vs_vehicle_infected_volcano.png', width = 9, height = 7)




# Heatmaps for top DE genes
top10_heatmap <- function(sig_res, title, filename) {
  if (is.null(sig_res) || nrow(sig_res) == 0) {
    cat("No sig genes for", title, "\n")
    return(invisible(NULL))
  }
  
  up   <- head(rownames(sig_res[order(sig_res$log2FoldChange, decreasing = TRUE), ]), 10)
  down <- head(rownames(sig_res[order(sig_res$log2FoldChange), ]), 10)
  
  genes <- c(up, down)
  
  if (length(genes) < 2) {
    cat("Too few genes for", title, "\n")
    return(invisible(NULL))
  }
  
  pheatmap(assay(vsd)[genes, ],
           annotation_col = ann,
           show_rownames = TRUE,
           fontsize_row = 8,
           main = title,
           filename = here::here("results", filename),   # یا here("results", filename) اگه library(here) لود کردی
           width = 10,
           height = 8)
}


top10_heatmap(sig_vehicle, 'Top 10 DE Genes - Inf and Uninf (Vehicle)', 'infected_vs_uninfected_vehicle_heatmap.png')
top10_heatmap(sig_dillapiole, 'Top 10 DE Genes - Inf and Uninf (Dillapiole)', 'infected_vs_uninfected_dillapiole_heatmap.png')
top10_heatmap(sig_uninf, 'Top 10 DE Genes - Treatment Effect (Uninfected)', 'dillapiole_vs_vehicle_uninfected_heatmap.png')  
top10_heatmap(sig_inf, 'Top 10 DE Genes - Treatment Effect (Infected)', 'dillapiole_vs_vehicle_infected_heatmap.png') 




# Gene Annotation
## add gene annotation 
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

add_names <- function(df) {
  if (nrow(df) == 0) return(df)
  ens <- rownames(df)
  mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  info <- getBM(attributes = c("ensembl_gene_id","external_gene_name","description"),
                filters = "ensembl_gene_id", values = ens, mart = mart)
  info$description <- gsub(" \\[Source.*\\]", "", info$description)
  merged <- merge(df, info, by.x="row.names", by.y="ensembl_gene_id", all.x=TRUE)
  merged <- merged[, c("external_gene_name","description", setdiff(names(merged), c("external_gene_name","description")))]
  rownames(merged) <- merged$Row.names
  merged$Row.names <- NULL
  return(merged)
}
  
  
top_up   <- head(sig_vehicle[order(sig_vehicle$log2FoldChange, decreasing = T),], 10)
top_down <- head(sig_vehicle[order(sig_vehicle$log2FoldChange),], 10)

write.csv(add_names(top_up), 'infection_vehicle_top10_up.csv', row.names = FALSE)
write.csv(add_names(top_down), 'infection_vehicle_top10_down.csv', row.names = FALSE)



