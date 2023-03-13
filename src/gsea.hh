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
using namespace std;

class Gsea {
private:
    struct GeneSample {
        string geneId;
        double count;
    };

    string expressionMatrixFilename;
    string geneSetsFilename;
    string outputFilename;
    char expressionMatrixSep;
    char geneSetsSep;
    char outputSep;

    uint nThreads;

    unordered_map<string, unordered_set<string>> geneSets;

    vector<vector<GeneSample>> expressionMatrix;
    vector<string> sampleNames;

    vector<vector<double>> results;

    uint nGenes;
    uint nSamples;

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

    void run();

    ~Gsea();
};