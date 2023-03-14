#include <Rcpp.h>
#include <iostream>
#include <vector>
#include <unordered_set>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
void gsea(List geneSets, IntegerMatrix countMatrix) {
    double nGenes = countMatrix.nrow();
    for (uint k = 0; k < geneSets.length(); ++k) {
        CharacterVector geneSetRcpp = geneSets[k];
        vector<string> geneVector = as<vector<string>> (geneSetRcpp);
        unordered_set<string> geneSet = unordered_set<string> (geneVector.begin(), geneVector.end());
        double geneSetSize = geneSetRcpp.length();
        double posScore = sqrt((nGenes - geneSetSize) / geneSetSize);
        double negScore = - sqrt((geneSetSize / (nGenes - geneSetSize)));
        for (uint j = 0; j < countMatrix.ncol(); ++j) {
            double currentValue = 0;
            double maxValue = 0;
            CharacterVector geneIdsCharVec = rownames(countMatrix);
            vector<string> geneIds = as<vector<string>> (geneIdsCharVec); 
            for (uint i = 0; i < countMatrix.nrow(); ++i) {
                // TODO Search gene ids in the countMatrix
                if (geneSet.find(geneIds[i]) != geneSet.end()) currentValue += posScore;
                else currentValue += negScore;
                if (i == 0) maxValue = currentValue;
                maxValue = max(currentValue, maxValue);
            }
        }
    }
}