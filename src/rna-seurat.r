library(Seurat)
library(gseacc)
library(dplyr)
library(patchwork)

phenotypes <- readCsv("data/rna/GSE121212_samples.csv")

groups <- apply(phenotypes, 1, function(x) {
    colnames(phenotypes)[which(x == 1)]
})

counts <- readCsv("data/rna/GSE121212_psoriasis.csv", "\t")

psoriasis <- CreateSeuratObject(counts = counts)
psoriasis <- AddMetaData(psoriasis, groups, "group")
psoriasis <- NormalizeData(psoriasis)
psoriasis <- FindVariableFeatures(psoriasis, selection.method = "vst", nfeatures = 2000)
psoriasis <- ScaleData(psoriasis, features = rownames(psoriasis))
psoriasis <- RunPCA(psoriasis, features = VariableFeatures(object = psoriasis))
psoriasis <- FindNeighbors(psoriasis, dims = 1:10)
psoriasis <- FindClusters(psoriasis, resolution = 0.85)
psoriasis <- RunUMAP(psoriasis, dims = 1:10)
DimPlot(psoriasis, reduction = "umap", group.by = "group")

psoriasis.markers <- FindAllMarkers(psoriasis, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
psoriasis.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(psoriasis, features = top10$gene, group.by = "group") + NoLegend()

counts <- readCsv("data/rna/GSE121212_GO_ES.csv")

gsea <- CreateSeuratObject(counts = counts)
gsea <- AddMetaData(gsea, groups, "group")
gsea <- NormalizeData(gsea)
gsea <- FindVariableFeatures(gsea, selection.method = "vst", nfeatures = 2000)
gsea <- ScaleData(gsea, features = rownames(gsea))
gsea <- RunPCA(gsea, features = VariableFeatures(object = gsea))
gsea <- FindNeighbors(gsea, dims = 1:10)
gsea <- FindClusters(gsea, resolution = 1)
gsea <- RunUMAP(gsea, dims = 1:10)
DimPlot(gsea, reduction = "umap", group.by = "group")

gsea.markers <- FindAllMarkers(gsea, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
gsea.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(gsea, features = top10$gene, group.by = "group") + NoLegend()
