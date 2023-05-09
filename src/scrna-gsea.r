library(gseacc)
# library(Seurat)
library(anndata)
# library(patchwork)

raw <- read_h5ad("data/healthy-raw.h5ad")

geneIds <- raw$var_names
sampleIds <- raw$obs_names
gsea <- new(GseaRcpp, sampleIds, geneIds)
rm(sampleIds)
rm(geneIds)

nChunks <- 1000
chunkSize <- as.integer(raw$n_obs / nChunks)
startSample <- 1
endSample <- as.integer(chunkSize)

for (i in 1:nChunks) {
    countMatrix <- raw$chunk_X(select = startSample:endSample)
    gsea$runChunked(countMatrix)
    startSample <- startSample + chunkSize
    endSample <- endSample + chunkSize
}

gsea$filterResults()
