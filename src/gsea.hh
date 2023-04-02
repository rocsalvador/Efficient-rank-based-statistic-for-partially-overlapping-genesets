#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <cmath>
#include <thread>
#include <chrono>
#include <cassert>
using namespace std;
using namespace chrono;

struct GeneSample {
    string geneId;
    double count;
    string sampleId;
};

class Gsea {
private:
    string expressionMatrixFilename;
    string geneSetsFilename;
    string outputFilename;
    char expressionMatrixSep;
    char geneSetsSep;
    char outputSep;

    uint nThreads;

    bool normalizedData;

    unordered_map<string, unordered_set<string>> geneSets;

    vector<vector<GeneSample>> expressionMatrix;
    vector<string> sampleNames;

    vector<vector<double>> results;

    uint nGenes;
    uint nSamples;

    uint nGeneSets;

    static bool geneSampleComp(const GeneSample& g1, const GeneSample& g2);

    void readConfig();

    void rpm();

    void meanCenter();

    void sortColumns();

    void enrichmentScore();

    void enrichmentScoreJob(uint sampleStart, uint sampleEnd);

    void writeResults();

public:
    Gsea();
    
    Gsea(unordered_map<string, unordered_set<string>>& geneSets, vector<vector<GeneSample>>& expressionMatrix);

    void run();

    ~Gsea();
};