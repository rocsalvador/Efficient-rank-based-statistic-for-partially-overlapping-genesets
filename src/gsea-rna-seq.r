library("reactome.db")
library("AnnotationDbi")
library("org.Hs.eg.db")

count_matrix_csv <- read.csv("../data/GSE121212_psoriasis.csv", sep = "\t")
count_matrix <- as.matrix(count_matrix_csv)
ensembl_ids <- rownames(count_matrix)

ensembl2go <- mapIds(org.Hs.eg.db,
                     keys = ensembl_ids,
                     keytype = "ENSEMBL",
                     column = "GO")
print(sum(is.na(ensembl2go)))
print(length(ensembl2go))
ptm <- proc.time()

for (ensembl_id in ensembl_ids) {
    go_id <- ensembl2go[[ensembl_id]]
    if (!is.na(go_id)) {
        reactome_id <- reactomeGO2REACTOMEID[[go_id]]
    }
}

print(proc.time() - ptm)
