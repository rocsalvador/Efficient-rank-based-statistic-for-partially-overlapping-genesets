#!/bin/bash

IN_FILENAME="data/GSE121212_sample_names.txt"
OUT_FILENAME="data/GSE121212_samples.csv"
group_names=("CTRL" "PSOnonlesional" "PSOlesional" "ADnonlesional" "ADlesional")
groups=("CTRL_.*" "PSO_.*_non-lesional" "PSO_.*_lesional" "AD_.*_non-lesional" "AD_.*_lesion.*")

echo -n > $OUT_FILENAME

FIRST=1
for GROUP in ${group_names[@]}; do
    echo $GROUP
    if $FIRST; then
        FIRST=0
    else
        echo -n ,$GROUP >> $OUT_FILENAME
    fi
done

echo >> $OUT_FILENAME

while IFS= read -r line; do
    curl https://www.ncbi.nlm.nih.gov/sra?term=$line > tmp.html
    echo -n $line, >> $OUT_FILENAME
    for GROUP in ${groups[@]}; do
        RET=$(grep -i "Sample: <span>$GROUP" tmp.html)
        if [ "$RET" != "" ]; then
            echo -n 1, >> $OUT_FILENAME
        else
            echo -n 0, >> $OUT_FILENAME
        fi
    done
    echo >> $OUT_FILENAME
done < $IN_FILENAME

rm tmp.html
