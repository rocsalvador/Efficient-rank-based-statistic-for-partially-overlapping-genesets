# library("reactome.db")
library("AnnotationDbi")
# library("org.Hs.eg.db")
library("gseacc")
library("limma")
library("Biobase")

countMatrix <- readCsv("./data/GSE121212_psoriasis.csv", sep = "\t")
# ensemblIds <- rownames(countMatrix

# xx <- as.list(org.Hs.egGO2ALLEGS)

# geneSets <- lapply(xx[1:10], function(x) {
#     mapIds(org.Hs.eg.db, keys = x, keytype = "ENTREZID", column = "ENSEMBL")
# })

geneSets <- readGeneSets("data/GO_gene_sets_entrez.csv")

gsea <- new(GseaRcpp, countMatrix, geneSets, 8)

gsea$normalizeExprMatrix()

gsea$run("data/results.csv", 8)

# Read GSEA results
x <- readCsv("data/GSE121212_GO_ES.csv", sep = ",")
p <- readCsv("data/GSE121212_samples.csv", sep = ",")
eset <- ExpressionSet(assayData = x, phenoData = AnnotatedDataFrame(p))
design <- model.matrix(~0 + CTRL + PSOnonlesional + PSOlesional + ADnonlesional + ADlesional, data = pData(eset))
# Fit the model
fit <- lmFit(eset, design)
# Calcuate t-statistics
fit <- eBayes(fit)
contrast.matrix <- makeContrasts(CTRLvsPSOnl=CTRL-PSOnonlesional, PSOnlvsPSOl=PSOnonlesional-PSOlesional, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
results <- decideTests(fit2)
topTable(fit2, sort.by="B", coef=1)
