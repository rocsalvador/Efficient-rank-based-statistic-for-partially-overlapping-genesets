library("reactome.db")
library("AnnotationDbi")
library("org.Hs.eg.db")
library("Rcpp")

rpm <- function(expr) {
    sums <- apply(expr, 2, sum)
    expr <- sweep(expr, 2, as.array(sums), "/")
    return(sweep(expr, 2, 10**6, "*"))
}

countMatrix <- as.matrix(read.csv("./data/GSE121212_psoriasis.csv",
                            sep = "\t"))
countMatrix <- rpm(countMatrix)
# TODO remove rows with only 0
# TODO rpm
# TODO mean center gene rows
ensemblIds <- rownames(countMatrix)

xx <- as.list(org.Hs.egGO2ALLEGS)

# Get ENSEMBL gene IDs for all the gene sets in GO
geneSets <- lapply(xx[1:200], function(x) {
    mapIds(org.Hs.eg.db, keys = x, keytype = "ENTREZID", column = "ENSEMBL")
})

N <- nrow(countMatrix)
enrichmentScore <- matrix(nrow = length(geneSets), ncol = ncol(countMatrix))

sourceCpp("src/rcpp/gsea-rcpp.cc")

gsea(geneSets, countMatrix)