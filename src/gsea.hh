#ifndef GSEA_HH
#define GSEA_HH

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

struct GeneSample
{
    uint32_t geneId;
    float count;
};

class Gsea
{
private:
    string expressionMatrixFilename;
    string geneSetsFilename;
    string outputFilename;
    char expressionMatrixSep;
    char geneSetsSep;
    char outputSep;
    uint ioutput;
    uint batchSize;

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

    static bool geneSampleComp(const GeneSample &g1, const GeneSample &g2);

    void readRna();

    void readScRna();

    void runScRna();

    void runRna();

    void readConfig();

    void rpm();

    void meanCenter();

    void sortGenes();

    void sortColumnsJob(uint columnStart, uint columnEnd);

    void enrichmentScore();

    void enrichmentScoreJob(uint sampleStart, uint sampleEnd);

    void scEnrichmentScoreJob(uint lineStart, uint lineEnd);

    void writeResults();

public:
    Gsea();

    Gsea(vector<string> &sampleIdsRcpp,
         vector<string> &geneIdsRcp);

    Gsea(unordered_map<string, unordered_set<string>> &geneSets,
         vector<vector<GeneSample>> &expressionMatrix,
         vector<string> &geneIds,
         vector<string> &sampleIds);

    void run();

    void runChunked(vector<vector<GeneSample>> &expressionMatrix);

    ~Gsea();
};

#endif