#include <Rcpp.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <unordered_set>
#include "gsea.hh"
#include "gsearcpp.hh"
using namespace Rcpp;
using namespace std;

RCPP_MODULE(GseaModule) {
    class_<GseaRcpp>("GseaRcpp")
    .constructor<CharacterVector, CharacterVector, List, uint>()
    .constructor<NumericMatrix, List, uint, bool>()
    .method("runChunked", &GseaRcpp::runChunked)
    .method("filterResults", &GseaRcpp::filterResults)
    .method("run", &GseaRcpp::run)
    .method("normalizeExprMatrix", &GseaRcpp::normalizeExprMatrix)
    ;
}

// [[Rcpp::export]]
void writeGeneSets(List geneSetsRcpp, String fileName)
{
    ofstream file(fileName);
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
List readGeneSets(std::string fileName)
{
    ifstream file(fileName);
    std::string line;
    getline(file, line);
    stringstream ssLine(line);
    List geneSets;
    CharacterVector geneSetIds;
    while (getline(file, line))
    {
        list<string> geneSet;
        stringstream ssLine(line);
        string rowName;
        getline(ssLine, rowName, ',');
        geneSetIds.push_back(rowName);

        string valueStr;
        while (getline(ssLine, valueStr, ','))
        {
            geneSet.push_back(valueStr);
        }
        geneSets.push_back(geneSet);
    }
    geneSets.names() = geneSetIds;
    return geneSets;
}

// [[Rcpp::export]]
NumericMatrix readCsv(String fileName, char sep) {
    ifstream file(fileName);
    std::string line;

    list<string> colNamesList;
    getline(file, line);
    stringstream ssLine(line);
    string colName;
    while (getline(ssLine, colName, sep))
    {
        colNamesList.push_back(colName);
    }

    list<string> rowNamesList;
    list<list<float>> listMatrix;
    while (getline(file, line))
    {
        stringstream ssLine(line);
        string rowName;
        getline(ssLine, rowName, sep);
        rowNamesList.push_back(rowName);

        string valueStr;
        list<float> row;
        while (getline(ssLine, valueStr, sep))
        {
            row.push_back(stof(valueStr));
        }
        listMatrix.push_back(row);
    }
    NumericMatrix matrix(listMatrix.size(), listMatrix.begin()->size());
    auto rowIt = listMatrix.begin();
    for (uint i = 0; i < matrix.nrow(); ++i, ++rowIt) {
        auto colIt = rowIt->begin();
        for (uint j = 0; j < matrix.ncol(); ++j, ++colIt) {
            matrix(i, j) = *colIt;
        }
    }
    CharacterVector rowNames(rowNamesList.size());
    auto it = rowNamesList.begin();
    for (uint i = 0; i < rowNames.length(); ++i, ++it) rowNames(i) = *it;
    CharacterVector colNames(colNamesList.size());
    it = colNamesList.begin();
    for (uint i = 0; i < colNames.length(); ++i, ++it) colNames(i) = *it;

    rownames(matrix) = rowNames;
    colnames(matrix) = colNames;
    return matrix;
}