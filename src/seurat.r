library(Seurat)
library(data.table)
library(future)

plan("multicore", workers = 8)
options(future.globals.maxSize = 891289600000)

counts <- data.frame(fread("data/healthy_GO_ES.csv", sep = ",", nrows = 985), row.names = 1)

healthy <- CreateSeuratObject(counts = counts)

rm(counts)
gc()

# healthy <- LoadH5Seurat("filtered-results.h5seurat")

groups <- sample(1:100, size = ncol(healthy), replace = TRUE)
names(groups) <- colnames(healthy)
healthy <- AddMetaData(object = healthy, metadata = groups, col.name = "group")
healthy.list <- SplitObject(healthy, split.by = "group")


#healthy.list <- lapply(X = healthy.list, FUN = function(x) {
#    x <- NormalizeData(x, verbose = TRUE)
#    x <- FindVariableFeatures(x, verbose = TRUE)
#})

#features <- SelectIntegrationFeatures(object.list = healthy.list)
features <- rownames(healthy)
rm(healthy)
gc()
healthy.list <- lapply(X = healthy.list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = TRUE)
    x <- RunPCA(x, features = features, verbose = TRUE)
})

anchors <- FindIntegrationAnchors(object.list = healthy.list,
                                  anchor.features = features,
                                  scale = FALSE,
                                  reduction = "rpca",
                                  dims = 1:5,
                                  verbose = TRUE)
healthy.integrated <- IntegrateData(anchorset = anchors, dims = 1:5, verbose = TRUE)

healthy.integrated <- ScaleData(healthy.integrated, verbose = FALSE)
healthy.integrated <- RunPCA(healthy.integrated, verbose = FALSE)
healthy.integrated <- RunUMAP(healthy.integrated, dims = 1:5)
