library(biomaRt)

count_matrix_csv <- read.csv("../data/GSE121212_psoriasis.csv", sep = "\t")
count_matrix <- as.matrix(count_matrix_csv)
ensembl_ids <- rownames(count_matrix)

# select mart and data set
bm <- useMart("ensembl")
bm <- useDataset("hsapiens_gene_ensembl", mart=bm)

# Get ensembl gene ids and GO terms
EG2GO <- getBM(mart = bm, attributes = c("ensembl_gene_id", "go_id"))

# examine result
head(EG2GO, 15)

# Remove blank entries
EG2GO <- EG2GO[EG2GO$go_id != '',]

# convert from table format to list format
geneID2GO <- by(EG2GO$go_id,
                EG2GO$ensembl_gene_id,
                function(x) as.character(x))

# examine result
head(geneID2GO)

for (ensembl_id in ensembl_ids) {
    print(geneID2GO[[ensembl_id]])
}
