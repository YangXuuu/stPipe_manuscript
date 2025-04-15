######### Identifying cell type for Visium mouse spleen data ######### 

library(Seurat)

obj <- readRDS("mouse_spleen/G000218_CytAssist_FFPE_709/obj.rds")

# Normalization
obj <- SCTransform(obj, assay = "Spatial", verbose = FALSE)

# PCA 
obj <- RunPCA(obj, npcs = 30, verbose = FALSE)

# clustering
obj <- FindNeighbors(obj, dims = 1:30)
obj <- FindClusters(obj, resolution = 1.2) 

# UMAP 
obj <- RunUMAP(obj, dims = 1:30)
DimPlot(obj, reduction = "umap", label = TRUE) + NoLegend()
# find all cluster marker gene
markers <- FindAllMarkers(
  obj, 
  assay = "SCT", 
  only.pos = TRUE, # only up-regulated
  logfc.threshold = 0.25, 
  min.pct = 0.1 
)

head(markers)
# marker gene list
gene_lists <- list(
  "B cell" = c("Cd19", "Cd22", "Ighd", "Cd5"),
  "Germinal centre" = c("Cxcr4", "Cd83", "Bcl6", "Rgs13", "Aicda"),
  "Neutrophil" = c("S100a9", "S100a8", "Ngp"),
  "Plasma cell" = c("Cd38", "Xbp1", "Irf4", "Prdm1", "Cd27"),
  "T cell" = c("Trac", "Cd3d", "Cd4", "Cd3e", "Cd8a"),
  "Red pulp" = c("Ifitm3", "C1qc", "Hmox1", "Hba-a1", "Klf1")
)

markers$annotation <- NA
# assign gene list based on avg_log2FC 
for (cell_type in names(gene_lists)) {
  genes <- gene_lists[[cell_type]]
  markers$annotation[markers$gene %in% genes] <- cell_type
}
library(dplyr)
# based on cluster ID 
cluster_annotation <- markers %>%
  group_by(cluster) %>%
  summarize(
    annotation = names(sort(table(annotation), decreasing = TRUE)[1]) 
  )

obj$cluster_annotation <- cluster_annotation$annotation[match(obj$seurat_clusters, cluster_annotation$cluster)]

# Visualization
DimPlot(obj, group.by = "cluster_annotation", label = TRUE) + NoLegend()
SpatialDimPlot(obj, group.by = "cluster_annotation", pt.size.factor = 5, cols = cell_type_colors)

obj_new <- subset(obj, cells = colnames(obj)[obj@meta.data$X_coordinate != 128])
SpatialDimPlot(obj_new, group.by = "cluster_annotation", pt.size.factor = 5, cols = cell_type_colors)


library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(dplyr)


cell_type_markers <- list(
  "B cell" = c("Cd19", "Cd22", "Ighd", "Cd5"),
  "Germinal centre" = c("Cxcr4", "Cd83", "Bcl6", "Rgs13", "Aicda"),
  "Neutrophil" = c("S100a9", "S100a8", "Ngp"),
  "Plasma cell" = c("Cd38", "Xbp1", "Irf4", "Prdm1", "Cd27"),
  "T cell" = c("Trac", "Cd3d", "Cd4", "Cd3e", "Cd8a"),
  "Red pulp" = c("Ifitm3", "C1qc", "Hmox1", "Hba-a1", "Klf1")
)

fold_change_data <- data.frame()

for (cell_type in names(cell_type_markers)) {
  genes <- cell_type_markers[[cell_type]]
  avg_expression <- AverageExpression(obj, features = genes, group.by = "seurat_clusters", assays = "SCT")$SCT
  avg_fc <- rowMeans(avg_expression) 
  fold_change_data <- rbind(
    fold_change_data,
    data.frame(
      CellType = cell_type,
      Cluster = colnames(avg_expression),
      FoldChange = avg_fc
    )
  )
}

