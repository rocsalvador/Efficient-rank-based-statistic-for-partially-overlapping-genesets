library("reactome.db")
library("AnnotationDbi")
library("org.Hs.eg.db")
library("gseacc")
library("limma")

countMatrix <- as.matrix(read.csv("./data/GSE121212_psoriasis.csv",
                         sep = "\t"))
ensemblIds <- rownames(countMatrix)

xx <- as.list(org.Hs.egGO2ALLEGS)

geneSets <- lapply(xx[1:10], function(x) {
    mapIds(org.Hs.eg.db, keys = x, keytype = "ENTREZID", column = "ENSEMBL")
})

gsea(geneSets, countMatrix)
