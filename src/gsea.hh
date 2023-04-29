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
    uint32_t geneId;
    float count;
};


class Gsea {
private:
    string expressionMatrixFilename;
    string geneSetsFilename;
    string outputFilename;
    char expressionMatrixSep;
    char geneSetsSep;
    char outputSep;
    uint ioutput;

    uint nThreads;
    uint logThread;

    bool normalizedData;
    bool scRna;

    unordered_map<string, unordered_set<string>> geneSets;

    vector<vector<GeneSample>> expressionMatrix;
    vector<string> sampleIds;
    vector<string> geneIds;

    vector<vector<float>> results;

    uint nGenes;
    uint nSamples;
    uint nGeneSets;

    system_clock::time_point startGSEATime;

    static bool geneSampleComp(const GeneSample& g1, const GeneSample& g2);

    void readRna();

    void readScRna();

    void readConfig();

    void rpm();

    void meanCenter();

    void sortColumns();

    void sortColumnsJob(uint columnStart, uint columnEnd);

    void enrichmentScore();

    void enrichmentScoreJob(uint sampleStart, uint sampleEnd);

    void writeResults();

public:
    Gsea();
    
    Gsea(unordered_map<string, unordered_set<string>>& geneSets,
         vector<vector<GeneSample>>& expressionMatrix,
         vector<string>& geneIds,
         vector<string>& sampleIds);

    void run();

    ~Gsea();
};
