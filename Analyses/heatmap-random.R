topVarGenes <- head(order(rowVars(assay(mRNA_rld)), decreasing = TRUE), 20)
mat  <- assay(mRNA_rld)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(mRNA_rld)[, c("Hist","HPV")])
pheatmap(mat, annotation_col = anno)
