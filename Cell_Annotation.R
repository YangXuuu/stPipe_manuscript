library(Seurat)
library(ggplot2)
library(reshape2)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(dplyr)

######################### Cell type annotation ######################################

# Extract gene IDs from the expression matrix
gene_ids <- rownames(expr_matrix_708)

# Remove version numbers from gene IDs (e.g., ENSMUSG00000121252.1 -> ENSMUSG00000121252)
gene_ids_no_version <- sub("\\..*", "", gene_ids)

# Convert Ensembl gene IDs to gene symbols using org.Mm.eg.db
gene_symbols <- mapIds(org.Mm.eg.db, keys = gene_ids_no_version,
                       column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")

# Replace row names in the expression matrix with gene symbols
expr_matrix_708 <- expr_matrix_708[!is.na(gene_symbols), ]
rownames(expr_matrix_708) <- gene_symbols[!is.na(gene_symbols)]

obj <- SCTransform(obj, assay = "Spatial", verbose = FALSE)
obj <- RunPCA(obj, assay = "SCT", verbose = FALSE)
obj <- FindNeighbors(obj, reduction = "pca", dims = 1:30)
obj <- FindClusters(obj, verbose = FALSE, resolution = 2)
obj <- RunUMAP(obj, reduction = "pca", dims = 1:30)
p1 <- DimPlot(obj, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(obj, label = TRUE, label.size = 3)
p1 + p2

# Assume we already have a Seurat object
obj <- WT_OCT_Manual_709_seu

wtf <- gsub("\\..*", "", rownames(obj))

rownames(obj) <- gsub("^\\.(\\d+)$", "ENSMUSG0000.\\1", rownames(obj))

# Define marker gene list
marker.genes <- list(
  "Macrophage" = c("Cd274", "Csf1r", "Adgre1", "Cd209b", "Cd80", "Cd68"),
  "B cell" = c("Cd19", "Cd22", "Ighd", "Cd5"),
  "Germinal centre" = c("Cxcr4", "Cd83", "Bcl6", "Rgs13", "Aicda"),
  "Neutrophil" = c("S100a9", "S100a8", "Ngp"),
  "Erythrocyte" = c("Car2", "Car1", "Klf1"),
  "Plasma cell" = c("Cd38", "Xbp1", "Irf4", "Prdm1", "Cd27"),
  "T cell" = c("Trac", "Cd3d", "Cd4", "Cd3e", "Cd8a")
)

# Original marker.genes list with Ensembl IDs
marker.genes <- list(
  "Macrophage" = c("ENSMUSG00000016496.8", "ENSMUSG00000024621.17", "ENSMUSG00000004730.16", "ENSMUSG00000065987.13" , "ENSMUSG00000075122.6", "ENSMUSG00000018774.14" ),
  "B cell" = c("ENSMUSG00000030724.8", "ENSMUSG00000030577.15",  "ENSMUSG00000104213.6",  "ENSMUSG00000024669.9" ),
  "Germinal centre" = c("ENSMUSG00000045382.7",  "ENSMUSG00000015396.5",  "ENSMUSG00000022508.6" , "ENSMUSG00000051079.9", "ENSMUSG00000040627.15" ),
  "Neutrophil" = c("ENSMUSG00000056071.13", "ENSMUSG00000056054.10" , "ENSMUSG00000032484.9" ),
  "Erythrocyte" = c("ENSMUSG00000027562.13", "ENSMUSG00000027556.16", "ENSMUSG00000054191.10" ),
  "Plasma cell" = c("ENSMUSG00000029084.6" ,"ENSMUSG00000020484.20" ,"ENSMUSG00000021356.11", "ENSMUSG00000038151.14" ,"ENSMUSG00000030336.15" ),
  "T cell" = c("ENSMUSG00000076928.6" , "ENSMUSG00000032094.9" ,"ENSMUSG00000023274.15",  "ENSMUSG00000032093.8", "ENSMUSG00000053977.14" )
)

# Remove version numbers from rownames
rownames_no_version <- sub("\\..*", "", rownames(obj))

# Define a function to match genes
match_genes <- function(genes, rownames_no_version) {
  matched_genes <- sapply(genes, function(gene) {
    match_idx <- which(rownames_no_version == gene)
    if (length(match_idx) > 0) {
      return(rownames(obj)[match_idx])
    } else {
      return(NA) # Return NA if no match found
    }
  })
  return(matched_genes)
}

# Apply matching to each cell type in marker.genes
new_marker.genes <- lapply(marker.genes, match_genes, rownames_no_version)

# Step 1: Calculate average expression of marker genes for each cluster
cluster_scores <- list()

for (cell_type in names(marker.genes)) {
  genes <- marker.genes[[cell_type]]
  
  # Filter marker genes that exist in the Seurat object
  genes_in_obj <- genes[genes %in% rownames(obj)]
  
  if (length(genes_in_obj) > 0) {
    # Calculate average expression per cluster
    cluster_avg_expr <- AverageExpression(obj, features = genes_in_obj, return.seurat = FALSE)
    cluster_scores[[cell_type]] <- cluster_avg_expr$SCT
  }
}

# Step 2: Calculate log-fold-change, excluding current cluster
log_fc_scores <- list()
# Loop through cluster expressions for each cell type
for (cell_type in names(cluster_scores)) {
  cluster_expr <- cluster_scores[[cell_type]]
  
  # Calculate log-fold-change for each cluster
  log_fc <- sapply(colnames(cluster_expr), function(curr_cluster) {
    curr_mean <- cluster_expr[, curr_cluster]  # Expression of the current cluster
    
    # Compute average of other clusters, excluding current
    other_clusters_mean <- rowMeans(cluster_expr[, colnames(cluster_expr) != curr_cluster, drop = FALSE], na.rm = TRUE)
    
    # Compute log-fold-change
    log_fc <- log2(curr_mean / other_clusters_mean)
    
    # Replace Inf values with NA
    log_fc[is.infinite(log_fc)] <- NA
    return(log_fc)
  })
  
  log_fc_scores[[cell_type]] <- log_fc
}

# Step 3: Convert results to data frame and plot heatmap
log_fc_df <- do.call(rbind, log_fc_scores)
log_fc_melted <- reshape2::melt(log_fc_df)

# Compute average value for each cell type
log_fc_melted <- log_fc_melted %>%
  group_by(Var1, Var2) %>%
  summarize(mean_value = mean(value, na.rm = TRUE))


# Map genes to cell types and calculate average values
average_log_fc <- log_fc_melted %>%
  mutate(CellType = case_when(
    Var1 %in% marker.genes$Macrophage ~ "Macrophage",
    Var1 %in% marker.genes$`B cell` ~ "B cell",
    Var1 %in% marker.genes$`Germinal centre` ~ "Germinal centre",
    Var1 %in% marker.genes$Neutrophil ~ "Neutrophil",
    Var1 %in% marker.genes$Erythrocyte ~ "Erythrocyte",
    Var1 %in% marker.genes$`Plasma cell` ~ "Plasma cell",
    Var1 %in% marker.genes$`T cell` ~ "T cell",
    TRUE ~ "Other"
  )) %>%
  group_by(CellType, Var2) %>%
  summarise(mean_value = mean(mean_value, na.rm = TRUE), .groups = 'drop')


# Step 3: Plot heatmap
ggplot(log_fc_melted, aes(x = Var2, y = Var1, fill = mean_value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limits = c(-4, 3)) +
  theme_minimal() +
  labs(title = "Cluster Scores Based on Marker Genes", x = "Cluster", y = "Cell Types") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

library(pheatmap)

# Convert data to wide format
heatmap_data <- average_log_fc %>%
  dcast(CellType ~ Var2, value.var = "mean_value")

# Remove CellType column for matrix creation
heatmap_matrix <- as.matrix(heatmap_data[, -1])

# Set row names
rownames(heatmap_matrix) <- heatmap_data$CellType

# Plot heatmap
pheatmap(heatmap_matrix,
         color = colorRampPalette(c("blue", "white", "red"))(50),
         main = "Heatmap of Average Log Fold Change",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         display_numbers = TRUE,
         fontsize_number = 10,
         fontsize_row = 10,
         fontsize_col = 10
)
dev.off()

