library(Seurat)
library(Signac)
library(patchwork)
library(monocle3)
library(SeuratWrappers)

library(ggplot2)
library( dplyr)

set.seed(1234)

# read the data
ES_day8 <- readRDS("./data/fig3DEH.rds")
DefaultAssay(ES_day8) <- 'RNA'

# Fig3A
p1 <-DimPlot(ES_day8, reduction = "umap", group.by = "celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <-DimPlot(ES_day8, reduction = "umap.atac",group.by = "celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <-DimPlot(ES_day8, reduction = "umap.wnn", group.by = "celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
pdf("Fig3A.pdf", width = 15, height = 6)
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
dev.off()

# pseudotime
set.seed(22)
cds <- SeuratWrappers::as.cell_data_set(ES_day8, assay = "RNA", reduction = "umap", group.by = "celltype")
cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(ES_day8[["RNA"]])

cds <- preprocess_cds(cds, method = "PCA")
cds <- reduce_dimension(cds, preprocess_method = "PCA",umap.n_neighbors= 14L,
                        reduction_method = "UMAP")
cds <- cluster_cells(cds, reduction_method = "UMAP")
cds <- learn_graph(cds, use_partition = FALSE, close_loop = FALSE)
cell_ids <- colnames(cds)[ES_day8$seurat_clusters ==  "0"]
closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
closest_vertex <- closest_vertex[cell_ids, ]
closest_vertex <- as.numeric(names(which.max(table(closest_vertex))))
mst <- principal_graph(cds)$UMAP
root_pr_nodes <- igraph::V(mst)$name[closest_vertex]

rowData(cds)$gene_name <- rownames(cds)
rowData(cds)$gene_short_name <- rowData(cds)$gene_name
cds <- order_cells(cds, root_pr_nodes = root_pr_nodes)

# Fig3B
plot_cells(cds, color_cells_by = "pseudotime",
           label_cell_groups =T, label_leaves = F,
           label_branch_points = F,show_trajectory_graph = T,
           graph_label_size = 3,label_groups_by_cluster = T)

# Fig3C
plot_cells(cds, color_cells_by = "cluster", cell_size = 1,
           label_cell_groups = FALSE,
           show_trajectory_graph = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE)
# Fig3G
ES_day8.seur <- as.Seurat(cds, assay = NULL, clusters = "UMAP")
ES_day8.seur<-AddMetaData(ES_day8.seur,metadata= cds@principal_graph_aux$UMAP$pseudotime,
                          col.name = "monocle3_pseudotime")

FeaturePlot(ES_day8.seur,features = c("Pax7", "Myod1"),
            reduction ="UMAP",combine = T,
            blend = TRUE,
            blend.threshold = 0.0, min.cutoff = 0,
            max.cutoff = 6)

FeaturePlot(ES_day8.seur,features = c("Myog", "Myod1"),
            reduction ="UMAP",combine = T,
            blend = TRUE,
            blend.threshold = 0.0, min.cutoff = 0,
            max.cutoff = 6)
