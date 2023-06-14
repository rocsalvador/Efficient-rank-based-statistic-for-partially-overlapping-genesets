library(Seurat)
library(future)
library(gseacc)
library(dplyr)
library(patchwork)

plan("multicore", workers = 8)
options(future.globals.maxSize = 891289600000)

counts <- readCsv("data/scrna/healthy_GO_ES.csv")
counts <- counts[, sample(ncol(counts), 20000)]

full.clustering <- read.table("data/scrna/healthy_full_clustering.csv", row.names = 1, sep = ",")

healthy <- CreateSeuratObject(counts = counts, names.field = 3, names.delim = "-")
healthy <- AddMetaData(healthy, full.clustering, "full.clustering")
healthy <- ScaleData(healthy, features = rownames(healthy))
healthy <- RunPCA(healthy, features = rownames(healthy))
healthy <- FindNeighbors(healthy, dims = 1:10)
healthy <- FindClusters(healthy, resolution = 1)
healthy <- RunUMAP(healthy, dims = 1:10)

DimPlot(healthy, reduction = "umap")
DimPlot(healthy, reduction = "umap", group.by = "full.clustering")

healthy.markers <- FindAllMarkers(healthy, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
healthy.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(healthy, features = top10$gene) + NoLegend()
DoHeatmap(healthy, features = top10$gene, group.by = "full.clustering") + NoLegend()
