#include <Rcpp.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <unordered_set>
#include "gsea.hh"
#include "gsearcpp.hh"
using namespace Rcpp;
using namespace std;

RCPP_MODULE(GseaModule) {
    class_<GseaRcpp>("GseaRcpp")
    .constructor<CharacterVector, CharacterVector>()
    .method("runChunked", &GseaRcpp::runChunked)
    ;
}

// [[Rcpp::export]]
void gsea(List geneSetsRcpp, IntegerMatrix countMatrixRcpp)
{
    unordered_map<string, unordered_set<string>> geneSets;
    CharacterVector geneSetsIdsRcpp = geneSetsRcpp.names();
    vector<string> geneSetsIds = as<vector<string>>(geneSetsIdsRcpp);
    for (uint i = 0; i < geneSetsRcpp.length(); ++i)
    {
        CharacterVector geneSetRcpp = geneSetsRcpp[i];
        vector<string> geneVector = as<vector<string>>(geneSetRcpp);
        unordered_set<string> geneSet = unordered_set<string>(geneVector.begin(), geneVector.end());
        geneSets.insert({geneSetsIds[i], geneSet});
    }

    double nGenes = countMatrixRcpp.nrow();
    double nSamples = countMatrixRcpp.ncol();
    vector<vector<GeneSample>> countMatrix(nGenes, vector<GeneSample>(nSamples));
    CharacterVector geneIdsRcpp = rownames(countMatrixRcpp);
    vector<string> geneIds = as<vector<string>>(geneIdsRcpp);
    CharacterVector sampleIdsRcpp = colnames(countMatrixRcpp);
    vector<string> sampleIds = as<vector<string>>(sampleIdsRcpp);
    for (uint i = 0; i < nGenes; ++i)
    {
        for (uint j = 0; j < nSamples; ++j)
        {
            countMatrix[i][j] = {i, float(countMatrixRcpp(i, j))};
        }
    }

    Gsea gsea(geneSets, countMatrix, geneIds, sampleIds);
    gsea.run();
}

// [[Rcpp::export]]
void writeGeneSets(List geneSetsRcpp)
{
    ofstream file("gene-sets.csv");
    unordered_map<string, unordered_set<string>> geneSets;
    CharacterVector geneSetsIdsRcpp = geneSetsRcpp.names();
    vector<string> geneSetsIds = as<vector<string>>(geneSetsIdsRcpp);
    for (uint i = 0; i < geneSetsRcpp.length(); ++i)
    {
        file << geneSetsIds[i] << ",";
        CharacterVector geneSetRcpp = geneSetsRcpp[i];
        vector<string> geneVector = as<vector<string>>(geneSetRcpp);
        for (uint j = 0; j < geneVector.size(); ++j)
        {
            if (j != 0)
                file << ",";
            file << geneVector[j];
        }
        file << endl;
    }
}

// [[Rcpp::export]]
List readGeneSets(std::string fileName, char sep)
{
    ifstream file(fileName);
    std::string line;
    getline(file, line);
    stringstream ssLine(line);
    List geneSets;
    CharacterVector geneSetIds;
    while (getline(file, line))
    {
        CharacterVector geneSet;
        stringstream ssLine(line);
        string rowName;
        getline(ssLine, rowName, sep);
        geneSetIds.push_back(rowName);

        string valueStr;
        while (getline(ssLine, valueStr, sep))
        {
            geneSet.push_back(valueStr);
        }
        geneSets.push_back(geneSet);
    }
    geneSets.names() = geneSetIds;
    return geneSets;
}


