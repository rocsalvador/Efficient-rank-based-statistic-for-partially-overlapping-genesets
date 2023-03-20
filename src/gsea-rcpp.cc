#include <Rcpp.h>
#include <iostream>
#include <vector>
#include <unordered_set>
#include "gsea.hh"
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
void gsea(List geneSetsRcpp, IntegerMatrix countMatrixRcpp) {
    unordered_map<string, unordered_set<string>> geneSets;
    CharacterVector geneSetsIdsRcpp = geneSetsRcpp.names();
    vector<string> geneSetsIds = as<vector<string>> (geneSetsIdsRcpp);
    for (uint i = 0; i < geneSetsRcpp.length(); ++i) {
        CharacterVector geneSetRcpp = geneSetsRcpp[i];
        vector<string> geneVector = as<vector<string>> (geneSetRcpp);
        unordered_set<string> geneSet = unordered_set<string> (geneVector.begin(), geneVector.end());
        geneSets.insert({geneSetsIds[i], geneSet});
    }

    double nGenes = countMatrixRcpp.nrow();
    double nSamples = countMatrixRcpp.ncol();
    vector<vector<GeneSample>> countMatrix (nGenes, vector<GeneSample> (nSamples));
    CharacterVector genesIdsRcpp = rownames(countMatrixRcpp);
    vector<string> genesIds = as<vector<string>> (genesIdsRcpp);
    for (int i = 0; i < nGenes; ++i) {
        for (int j = 0; j < nSamples; ++j) {
            countMatrix[i][j] = GeneSample{genesIds[i], double(countMatrixRcpp(i, j))};
        }
    }

    Gsea gsea(geneSets, countMatrix);
    gsea.run();
}