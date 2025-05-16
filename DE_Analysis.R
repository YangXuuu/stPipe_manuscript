################## DE Analysis ##################
FFPE_Cyt_709$dataset <- "FFPE_Cyt_709"
FFPE_Manual_709$dataset <- "FFPE_Manual_709"
OCT_Cyt_709$dataset <- "OCT_Cyt_709"
OCT_Manual_709$dataset <- "OCT_Manual_709"

library(Seurat)
library(limma)
library(ggplot2)

b_cells_173 <- subset(OCT_173_seu, idents = c(2, 1, 8))
expr_matrix_173 <- as.matrix(GetAssayData(b_cells_173, slot = "data"))
b_cells_174 <- subset(Oct_Manual_174, idents = c(2, 6, 3))
expr_matrix_174 <- as.matrix(GetAssayData(b_cells_174, slot = "data"))
b_cells_545 <- subset(OCT_WT_545_seu, idents = c(2, 14, 3))
expr_matrix_545 <- as.matrix(GetAssayData(b_cells_545, slot = "data"))
b_cells_708 <- subset(OCT_WT_708_seu, idents = c(0, 15, 4))
expr_matrix_708 <- as.matrix(GetAssayData(b_cells_708, slot = "data"))
b_cells_709 <- subset(obj, idents = c(3, 4, 10))
expr_matrix_709 <- as.matrix(GetAssayData(b_cells_709, slot = "data"))

common_genes <- Reduce(intersect, list(
  rownames(expr_matrix_173),
  rownames(expr_matrix_174),
  rownames(expr_matrix_545),
  rownames(expr_matrix_708),
  rownames(expr_matrix_709)
))

expr_matrix_173 <- expr_matrix_173[common_genes, ]
expr_matrix_174 <- expr_matrix_174[common_genes, ]
expr_matrix_545 <- expr_matrix_545[common_genes, ]
expr_matrix_708 <- expr_matrix_708[common_genes, ]
expr_matrix_709 <- expr_matrix_709[common_genes, ]

sum_counts_173 <- rowSums(expr_matrix_173)
sum_counts_174 <- rowSums(expr_matrix_174)
sum_counts_545 <- rowSums(expr_matrix_545)
sum_counts_708 <- rowSums(expr_matrix_708)
sum_counts_709 <- rowSums(expr_matrix_709)
combined_matrix <- cbind(sum_counts_173, sum_counts_174, sum_counts_545, sum_counts_708, sum_counts_709)


group <- factor(c(rep("Male", 2), rep("Female", 3)))
dge <- DGEList(counts = combined_matrix, group = group)
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
dge <- estimateDisp(dge, design)
fit <- glmQLFit(dge, design)
contrast <- makeContrasts(Male_vs_Female = Male - Female, levels = design)
qlf <- glmQLFTest(fit, contrast = contrast)

topTable <- topTags(qlf, n = Inf)
plotMD(qlf, main = "MA Plot of Male vs Female", col = c("black", "red")[as.numeric(topTable$table$logFC > 0) + 1])

logFC <- qlf$table$logFC
weights <- rep(0, length(logFC))
weights[match(XiE, rownames(qlf$table))] <- -1   
weights[match(chromo_Y, rownames(qlf$table))] <- 1  

topTable$table$significant <- ifelse(topTable$table$PValue < 0.05 & topTable$table$logFC > 1,
                                     "up",
                                     ifelse(topTable$table$PValue < 0.05 & topTable$table$logFC < -1,
                                            "down",
                                            "non-significant"))

barcodeplot(logFC,
            gene.weights = weights,
            weights.label = "logFC",
            index = significant == "up", index2 = significant == "down",
            main = "Barcode Plot of sex-specific genes in Male vs. Female B cell DE Analysis")
legend("topright", legend = c("XiE genes", "chrY genes"),
       fill = c("blue", "red"), border = "black", cex = 1.2)

chromo_Y <- c("Gm47283", "Eif2s3y", "Ddx3y", "Kdm5d", "Uty", "Gm21860", "Gm29650")

XiE <- c("Cybb", "Ddx3x", "Kdm6a", "Cfp", "Utp14a", "Firre", "Bgn", "5430427O19Rik", "Eif2s3x", "Vsig4", "Xist", "Ftx", "5530601H04Rik", "Pbdc1", "5730416F02Rik", "Kdm5c", "Tmsb4x", "Esm1")


XiE_genes <- which(rownames(topTable) %in% XiE)
XiE_results <- topTable[rownames(topTable) %in% XiE, ]


#  MA plot
plotMD(qlf, main = "MA Plot of Male vs Female", col = "black")
points(topTable[XiE_genes, "logCPM"], topTable[XiE_genes, "logFC"], col = "blue", pch = 16)

missing_genes <- XiE[!XiE %in% rownames(topTable)]
if (length(missing_genes) > 0) {
  print(paste("Missing genes in topTable:", paste(missing_genes, collapse = ", ")))
}
XiE_genes <- which(rownames(topTable) %in% XiE)
valid_XiE_genes <- XiE[XiE %in% rownames(topTable)]

XiE_genes_color <- "blue"
chromY_genes_color <- "red"

plotMD(qlf, main = "MA Plot of Male vs Female in B cell DE analysis")
valid_XiE_genes <- intersect(rownames(topTable$table), XiE)

if (length(valid_XiE_genes) > 0) {
  points(topTable$table[valid_XiE_genes, "logCPM"], topTable$table[valid_XiE_genes, "logFC"],
         col = "blue", pch = 16)
} else {
  print("No valid XiE genes found in topTable for plotting.")
}

points(topTable$table[chromo_Y, "logCPM"], topTable$table[chromo_Y, "logFC"],
       col = "red", pch = 16)

legend("topright", legend = c("Other", "XiE gene", "chrY gene"),
       pch = c(16, 16, 16), col = c("black", "blue", "red"),
       box.lwd = 1.4, cex = 1.2)

