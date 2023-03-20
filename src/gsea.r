library("reactome.db")
library("AnnotationDbi")
library("org.Hs.eg.db")
library("Rcpp")
library("gseacc")

countMatrix <- as.matrix(read.csv("./data/GSE121212_psoriasis.csv",
                            sep = "\t"))
# TODO remove rows with only 0
# TODO rpm
# TODO mean center gene rows
ensemblIds <- rownames(countMatrix)

xx <- as.list(org.Hs.egGO2ALLEGS)

# Get ENSEMBL gene IDs for all the gene sets in GO
geneSets <- lapply(xx[1:200], function(x) {
    mapIds(org.Hs.eg.db, keys = x, keytype = "ENTREZID", column = "ENSEMBL")
})

gsea(geneSets, countMatrix)
