library("reactome.db")
library("AnnotationDbi")
library("org.Hs.eg.db")
library("gseacc")
library("limma")
library("Biobase")

countMatrix <- as.matrix(read.csv("./data/GSE121212_psoriasis.csv",
                         sep = "\t"))
ensemblIds <- rownames(countMatrix)

xx <- as.list(org.Hs.egGO2ALLEGS)

geneSets <- lapply(xx[1:10], function(x) {
    mapIds(org.Hs.eg.db, keys = x, keytype = "ENTREZID", column = "ENSEMBL")
})

geneSets <- readGeneSets("data/go-gene-sets.csv", ",")

gsea(geneSets, countMatrix)

x <- read.csv("data/go-gsea.csv", sep = ",")
p <- read.csv("data/GSE121212_samples.csv", sep = ",")
eset <- ExpressionSet(assayData = as.matrix(x), phenoData = AnnotatedDataFrame(p))
design <- model.matrix(~0 + CTRL + PSOnonlesional + PSOlesional + ADnonlesional + ADlesional, data = pData(eset))
# Fit the model
fit <- lmFit(eset, design)
# Calcuate t-statistics
fit <- eBayes(fit)
contrast.matrix <- makeContrasts(CTRLvsPSOnl=CTRL-PSOnonlesional, PSOnlvsPSOl=PSOnonlesional-PSOlesional, levels=design)
results <- decideTests(fit[, "CTRL"])
summary(results)
