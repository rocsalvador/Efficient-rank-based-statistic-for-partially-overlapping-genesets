library(gseacc)
library(anndata)

raw <- read_h5ad("data/scrna/healthy_raw.h5ad")

geneSets <- readGeneSets("data/scrna/GO_gene_sets_symbol.csv")
geneIds <- raw$var_names
sampleIds <- raw$obs_names
gsea <- new(Gsea, sampleIds, geneIds, geneSets, 0)
rm(sampleIds)
rm(geneIds)

nChunks <- 100
chunkSize <- as.integer(raw$n_obs / nChunks)
offset <- raw$n_obs %% nChunks
startSample <- 0
endSample <- as.integer(chunkSize - 1)

for (i in 1:nChunks) {
    countMatrix <- raw$chunk_X(select = startSample:endSample)
    gsea$runChunked(countMatrix)
    startSample <- startSample + chunkSize
    endSample <- endSample + chunkSize
}

endSample <- startSample + offset - 1
if (startSample < endSample) {
    countMatrix <- raw$chunk_X(select = startSample:endSample)
    gsea$runChunked(countMatrix)
}

gsea$filterResults(1000, "", "filtered-results.csv")
