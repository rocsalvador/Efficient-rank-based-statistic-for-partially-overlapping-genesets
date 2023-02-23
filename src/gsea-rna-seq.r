library("reactome.db")
library("AnnotationDbi")
library("org.Hs.eg.db")

countMatrixCsv <- as.matrix(read.csv("../data/GSE121212_psoriasis.csv", sep = "\t"))
countMatrix <- as.matrix(countMatrixCsv)
ensemblIds <- rownames(countMatrix)

xx <- as.list(org.Hs.egGO2ALLEGS)

# Get ENSEMBL gene IDs for all the gene sets in GO
geneSets <- lapply(xx[1:1], function(x) {
    mapIds(org.Hs.eg.db, keys = x, keytype = "ENTREZID", column = "ENSEMBL")
})

nGenes <- length(rownames(countMatrix))

sortedCountMatrix <- apply(countMatrix, 2, function(x) {
    sort(x, decreasing = TRUE, index.return = TRUE)
})

resultspMatrix <- matrix(ncol = length(colnames(countMatrix)),
                         nrow = length(geneSets))
resultsnMatrix <- matrix(ncol = length(colnames(countMatrix)),
                         nrow = length(geneSets))

i <- 0
for (geneSet in geneSets) {
    appearances <- lapply(rownames(countMatrix), function(x) {
        return(x %in% geneSet)
    })
    names(appearances) <- rownames(countMatrix)
    nAppearances <- length(appearances[appearances == TRUE])
    negVal <- - sqrt(nAppearances / (nGenes - nAppearances))
    posVal <- sqrt((nGenes - nAppearances) / nAppearances)

    j <- 0
    for (sample in sortedCountMatrix) {
        first <- TRUE
        for (geneId in names(sample$x)) {
            currentVal <- 0
            if (appearances[[geneId]]) currentVal <- currentVal + posVal
            else currentVal <- currentVal + negVal
            if (first) {
                first <- FALSE
                minVal <- currentVal
                maxVal <- currentVal
            }
            minVal <- min(c(minVal, currentVal))
            maxVal <- max(c(maxVal, currentVal))
        }
        resultspMatrix[i, j] <- maxVal
        resultsnMatrix[i, j] <- minVal
        j <- j + 1
        print(i)
        print(j)
    }

    i <- i + 1
}

print(resultsnMatrix)
