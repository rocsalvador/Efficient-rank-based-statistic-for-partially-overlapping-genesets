#include "gsearcpp.hh"

GseaRcpp::GseaRcpp(CharacterVector sampleIdsRcpp,
                   CharacterVector geneIdsRcpp)
{
    vector<string> sampleIds = as<vector<string>> (sampleIdsRcpp);
    vector<string> geneIds = as<vector<string>> (geneIdsRcpp);
    gsea = new Gsea(sampleIds, geneIds);
}

void GseaRcpp::runChunked(const IntegerMatrix &countMatrixRcpp)
{
    double nGenes = countMatrixRcpp.ncol();
    double nSamples = countMatrixRcpp.nrow();
    vector<vector<GeneSample>> expressionMatrix(nSamples, vector<GeneSample>(nGenes));
    for (uint i = 0; i < nSamples; ++i)
    {
        for (uint j = 0; j < nGenes; ++j)
        {
            expressionMatrix[i][j] = {i, float(countMatrixRcpp(i, j))};
        }
    }

    gsea->runChunked(expressionMatrix);
}

GseaRcpp::~GseaRcpp()
{
    delete gsea;
}

