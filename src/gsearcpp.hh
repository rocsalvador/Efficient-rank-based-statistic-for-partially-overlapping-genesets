#include <Rcpp.h>
#include "gsea.hh"
using namespace Rcpp;

class GseaRcpp {
private:
    Gsea *gsea;
public:
    GseaRcpp(CharacterVector sampleIdsRcpp,
             CharacterVector geneIdsRcp);

    void runChunked(const IntegerMatrix &countMatrixRcpp);

    void filterResults();

    ~GseaRcpp();
};
