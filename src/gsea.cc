#include "gsea.hh"
#include <bits/chrono.h>
#include <cassert>
#include <cstddef>
#include <fstream>
#include <thread>
#include <vector>


void printTime(system_clock::time_point timePoint) {
    char timeString[9];
    time_t timePointC = system_clock::to_time_t(timePoint);
    struct tm *tm = localtime(&timePointC);
    strftime(timeString, sizeof(timeString), "%H:%M:%S", tm);
    cout << "[" << timeString << "]";
}

void Gsea::readConfig() {
    ifstream file("./gsea.config");
    if (!file.is_open()) {
        file.close();
        ofstream outFile("./gsea.config");
        string aux;
        outFile << "expression-matrix-file:     expression-matrix.csv" << endl;
        outFile << "sep:                        ," << endl;
        outFile << "gene-sets-file:             gene-sets.csv" << endl;
        outFile << "sep:                        ," << endl;
        outFile << "output-file:                results.csv" << endl;
        outFile << "sep:                        ," << endl;
        outFile << "threads-used:               0" << endl;
        outFile << "normalized-data:            0" << endl;
        outFile << "ioutput:                    100" << endl;
        outFile << "scrna:                      0" << endl;
        outFile << "batch-size:                 50" << endl;

        expressionMatrixFilename = "expression-matrix.csv";
        expressionMatrixSep = ',';
        geneSetsFilename = "gene-sets.csv";
        geneSetsSep = ',';
        outputFilename = "results.csv";
        outputSep = ',';
        nThreads = 0;
        normalizedData = false;
        ioutput = 0;
        scRna = 0;
        batchSize = 50;
        outFile.close();
    } else {
        string aux;
        file >> aux >> expressionMatrixFilename;
        file >> aux >> expressionMatrixSep;
        if (expressionMatrixSep == 't') expressionMatrixSep = '\t';
        file >> aux >> geneSetsFilename;
        file >> aux >> geneSetsSep;
        file >> aux >> outputFilename;
        file >> aux >> outputSep;
        file >> aux >> nThreads;
        file >> aux >> normalizedData;
        file >> aux >> ioutput;
        file >> aux >> scRna;
        file >> aux >> batchSize;
    }

    if (nThreads == 0) nThreads = thread::hardware_concurrency();

    cout << "[GSEA config]" << endl;
    cout << "expression-matrix-file: " << expressionMatrixFilename << endl;
    cout << "sep:                    " << expressionMatrixSep << endl;
    cout << "gene-sets-file:         " << geneSetsFilename << endl;
    cout << "sep:                    " << geneSetsSep << endl;
    cout << "output-file:            " << outputFilename << endl;
    cout << "sep:                    " << outputSep << endl;
    cout << "threads-used:           " << nThreads << endl;
    cout << "normalized-data:        " << normalizedData << endl;
    cout << "ioutput:                " << ioutput << endl;
    cout << "scrna:                  " << scRna << endl;
    cout << "batch-size:             " << batchSize << endl;
    cout << endl;

    file.close();
}

Gsea::Gsea() {
    system_clock::time_point startIOTime = system_clock::now();

    readConfig();

    if (!scRna) readRna();
    else readScRna();

    system_clock::time_point endIOTime = system_clock::now();
    cout << "IO elapsed time: " << duration_cast<milliseconds> (endIOTime - startIOTime).count() / 1000.0 << " s" << endl;
}

void Gsea::readRna() {
    ifstream file(expressionMatrixFilename);
    string line;

    // Read first row (sample ids)
    getline(file, line);
    stringstream ssLine(line);
    string colName;
    uint i = 0;
    while(getline(ssLine, colName, expressionMatrixSep)){
        sampleIds.push_back(colName);
        ++i;
    }

    i = 0;
    while(getline(file, line)) {
        stringstream ssLine(line);

        // Read first column (gene id)
        string rowName;
        getline(ssLine, rowName, expressionMatrixSep);
        geneIds.push_back(rowName);

        string valueStr;
        bool first = true;
        bool nullRow = true;
        uint j = 0;
        while(getline(ssLine, valueStr, expressionMatrixSep)){
            if (first) {
                expressionMatrix.push_back(vector<GeneSample> (sampleIds.size()));
                first = false;
            }
            float count = stof(valueStr);
            if (nullRow and count != 0) nullRow = false;
            expressionMatrix[i][j] = {i, count};
            ++j;
        }
        if (nullRow) {
            expressionMatrix.pop_back();
            geneIds.pop_back();
        } else ++i;
    }

    nGenes = expressionMatrix.size();
    if (nGenes > 0)
        nSamples = expressionMatrix[0].size();
    file.close();

    file = ifstream(geneSetsFilename);
    while(getline(file, line)) {
        stringstream ssLine(line);

        string rowName;
        getline(ssLine, rowName, geneSetsSep);

        unordered_set<string> genes;
        string valueStr;
        while(getline(ssLine, valueStr, geneSetsSep)) {
            genes.insert(valueStr);
        }
        geneSets.insert({rowName, genes});
    }
    file.close();

    results = vector<vector<float>> (geneSets.size(), vector<float> (nSamples));
}

void Gsea::readScRna() {
    // Read gene sets file
    ifstream file = ifstream(geneSetsFilename);
    string line;
    uint i = 0;
    while(getline(file, line)) {
        stringstream ssLine(line);

        string rowName;
        getline(ssLine, rowName, geneSetsSep);

        unordered_set<string> genes;
        string valueStr;
        while(getline(ssLine, valueStr, geneSetsSep)) {
            genes.insert(valueStr);
        }
        geneSets.insert({rowName, genes});
        ++i;
    }
    file.close();

    file = ifstream(expressionMatrixFilename);

    // Read first row (gene ids)
    getline(file, line);
    stringstream ssLine(line);
    string colName;
    i = 0;
    while(getline(ssLine, colName, expressionMatrixSep)){
        geneIds.push_back(colName);
        ++i;
    }
    nGenes = geneIds.size();

    nSamples = 0;
    while(getline(file, line)) ++nSamples;

    file.close();
}

void Gsea::runScRna() {
    ifstream file = ifstream(expressionMatrixFilename);
    ofstream oFile = ofstream(outputFilename);
    string line;

    // Ignore first row, already read
    getline(file, line);

    bool first = true;
    for (auto geneSet : geneSets) {
        if (first) first = false;
        else oFile << ",";
        oFile << geneSet.first;
    }
    oFile << endl;

    uint i = 0;
    uint nLines = batchSize;
    uint totalLines = nThreads * nLines;
    expressionMatrix = vector<vector<GeneSample>> (totalLines);
    vector<string> sampleNames = vector<string> (totalLines);
    results = vector<vector<float>> (totalLines, vector<float> (geneSets.size()));
    first = true;
    while(getline(file, line)) {
        stringstream ssLine(line);

        // Read first column (sample id)
        string rowName;
        getline(ssLine, rowName, expressionMatrixSep);
        sampleNames[i % totalLines] = rowName;

        string valueStr;
        uint j = 0;
        while(getline(ssLine, valueStr, expressionMatrixSep)){
            if (first) expressionMatrix[i % totalLines].push_back({j, stof(valueStr)});
            else expressionMatrix[i % totalLines][j] = {j, stof(valueStr)};
            ++j;
        }

        if (i != 0 and i % (totalLines - 1) == 0) {
            first = false;
            vector<thread> threads = vector<thread> (nThreads);
            for (uint t = 0; t < nThreads; ++t) {
                uint startLine = nLines * t;
                uint endLine = startLine + nLines;
                threads[t] = thread(&Gsea::scEnrichmentScoreJob, this, startLine, endLine);
            }

            for (thread& t : threads) t.join();

            for (uint t = 0; t < totalLines; ++t) {
                oFile << sampleNames[t];
                for (uint l = 0; l < geneSets.size(); ++l) {
                    oFile << results[t][l] << ",";
                }
                oFile << endl;
            }
            system_clock::time_point now = system_clock::now();
            printTime(now);
            cout << " Sample " << i;
            uint ETA = (nSamples - i) * duration_cast<milliseconds> (now - startGSEATime).count() / (i * 60 * 1000) ;
            cout << " ETA: " << ETA << " min" << endl;
        }
        ++i;
    }

    uint offset = i % totalLines;
    uint linesPerThread = offset / nThreads;
    uint offsetLines = linesPerThread % nThreads;
    vector<thread> threads = vector<thread> (nThreads);
    for (uint t = 0; t < nThreads; ++t) {
        uint startLine = linesPerThread * t;
        uint endLine = startLine + linesPerThread;
        if (t == nThreads - 1) endLine += offsetLines;
        threads[t] = thread(&Gsea::scEnrichmentScoreJob, this, startLine, endLine);
    }

    for (thread& t : threads) t.join();

    for (uint t = 0; t < offset; ++t) {
        oFile << sampleNames[t];
        for (uint l = 0; l < geneSets.size(); ++l) {
            oFile << results[t][l] << ",";
        }
        oFile << endl;
    }

    file.close();
}


Gsea::Gsea(unordered_map<string, unordered_set<string>>& geneSets,
           vector<vector<GeneSample>>& expressionMatrix,
           vector<string>& geneIds,
           vector<string>& sampleIds) {
    readConfig();

    this->geneSets = geneSets;
    this->expressionMatrix = expressionMatrix;
    nGenes = expressionMatrix.size();
    nSamples = sampleIds.size();
    this->sampleIds = sampleIds;
    this->geneIds = geneIds;
    results = vector<vector<float>> (geneSets.size(), vector<float> (nSamples));
}

void Gsea::rpm() {
    for (uint j = 0; j < nSamples; ++j) {
        float sum = 0;

        for (uint i = 0; i < nGenes; ++i) 
            sum += expressionMatrix[i][j].count;
        
        float multFactor = 1000000 / sum;
        for (uint i = 0; i < nGenes; ++i) 
            expressionMatrix[i][j].count *= multFactor;
    }
}

void Gsea::meanCenter() {
    for (uint i = 0; i < nGenes; ++i) {
        float mean = 0;
        for (uint j = 0; j < nSamples; ++j)
            mean += expressionMatrix[i][j].count;
        mean /= nSamples;

        for (uint j = 0; j < nSamples; ++j)
            expressionMatrix[i][j].count -= mean;
    }
}

bool Gsea::geneSampleComp(const GeneSample& g1, const GeneSample& g2) {
    return g1.count > g2.count;
}

void Gsea::sortColumnsJob(uint startSample, uint endSample) {
    for (uint j = startSample; j < endSample; ++j) {
        vector<GeneSample> column = vector<GeneSample> (nGenes);
        for (uint i = 0; i < nGenes; ++i) column[i] = expressionMatrix[i][j];
        sort(column.begin(), column.end(), &Gsea::geneSampleComp);
        for (uint i = 0; i < nGenes; ++i) expressionMatrix[i][j] = column[i];
    }
}

void Gsea::sortGenes() {
    uint samplesPerThread = nSamples / nThreads;
    uint offset = nSamples % samplesPerThread;
    vector<thread> threads = vector<thread> (nThreads);
    for (uint i = 0; i < nThreads; ++i) {
        uint startSample = i * samplesPerThread;
        uint endSample = startSample + samplesPerThread;
        if (i == nThreads - 1) endSample += offset;

        threads[i] = thread(&Gsea::sortColumnsJob, this, startSample, endSample);

        if (i == nThreads - 1) {
            auto threadId = threads[i].get_id();
            logThread = *static_cast<unsigned int*>(static_cast<void*>(&threadId));
        }
    }

    for (thread& t : threads) t.join();
}

void Gsea::enrichmentScore() {
    uint samplesPerThread = nSamples / nThreads;
    uint offset = nSamples % samplesPerThread;
    vector<thread> threads = vector<thread> (nThreads);
    for (uint i = 0; i < nThreads; ++i) {
        uint startSample = i * samplesPerThread;
        uint endSample = startSample + samplesPerThread;
        if (i == nThreads - 1) endSample += offset;
        threads[i] = thread(&Gsea::enrichmentScoreJob, this, startSample, endSample);
        if (i == nThreads - 1) {
            auto threadId = threads[i].get_id();
            logThread = *static_cast<unsigned int*>(static_cast<void*>(&threadId));
        }
    }

    for (thread& t : threads) t.join();
}

void Gsea::enrichmentScoreJob(uint startSample, uint endSample) {
    assert(endSample < nSamples);
    auto it = geneSets.begin();
    auto threadId = this_thread::get_id();
    uint id = *static_cast<unsigned int*>(static_cast<void*>(&threadId));
    for (uint k = 0; k < geneSets.size(); ++k) {
        float geneSetSize = it->second.size();
        float posScore = sqrt((nGenes - geneSetSize) / geneSetSize);
        float negScore = - sqrt((geneSetSize / (nGenes - geneSetSize)));
        for (uint j = startSample; j < endSample; ++j) {
            float currentValue = 0;
            float maxValue = 0;
            for (uint i = 0; i < nGenes; ++i) {
                if (it->second.find(geneIds[expressionMatrix[i][j].geneId]) != it->second.end()) {
                    currentValue += posScore;
                } else currentValue += negScore;

                if (i == 0) maxValue = currentValue;
                maxValue = max(currentValue, maxValue);
            }
            results[k][j] = maxValue;
        }
        if (k != 0 and id == logThread and k % ioutput == 0) {
            system_clock::time_point now = system_clock::now();
            printTime(now);
            cout << " Gene set " << k;

            uint ETA = (geneSets.size() - k) * duration_cast<milliseconds> (now - startGSEATime).count() / (k * 60 * 1000) ;
            cout << " ETA: " << ETA << " min" << endl;
        }
        it = next(it);
    }
}


void Gsea::scEnrichmentScoreJob(uint lineStart, uint lineEnd) {
    for (uint i = lineStart; i < lineEnd; ++i) {
        sort(expressionMatrix[i].begin(), expressionMatrix[i].end(), &Gsea::geneSampleComp);
    }
    auto threadId = this_thread::get_id();
    for (uint i = lineStart; i < lineEnd; ++i) {
        auto it = geneSets.begin();
        for (uint k = 0; k < geneSets.size(); ++k, ++it) {
            uint geneSetSize = it->second.size();
            float posScore = sqrt((nGenes - geneSetSize) / geneSetSize);
            float negScore = - sqrt((geneSetSize / (nGenes - geneSetSize)));
            float currentValue = 0;
            float maxValue = 0;
            for (uint j = 0; j < expressionMatrix[i].size(); ++j) {
                if (expressionMatrix[i][j].count == 0) break;

                if (it->second.find(geneIds[expressionMatrix[i][j].geneId]) != it->second.end()) {
                    currentValue += posScore;
                } else currentValue += negScore;

                if (j == 0) maxValue = currentValue;
                maxValue = max(currentValue, maxValue);
            }
            results[i][k] = maxValue;
        }
    }
}

void Gsea::writeResults() {
    ofstream file(outputFilename);
    bool first = true;
    for (string& sampleId : sampleIds) {
        if (first) first = false;
        else file << outputSep;
        file << sampleId;
    }
    file << endl;

    vector<string> geneSetNames = vector<string> (geneSets.size());
    uint i = 0;
    for (auto& it : geneSets) {
        geneSetNames[i] = it.first;
        ++i;
    } 

    i = 0;
    for (vector<float>& row : results) {
        file << geneSetNames[i];
        for (float value : row) {
            file << outputSep << value;
        }
        file << endl;
        ++i;
    }
    file.close();
}

void Gsea::runRna() {
    if (!normalizedData) {
        rpm();

        meanCenter();
    }

    sortGenes();

    enrichmentScore();

    writeResults();

    cout << "Results written in " << outputFilename << endl;
}

void Gsea::run() {
    cout << "[GSEA input size]" << endl;
    cout << "Sampled genes: " << nGenes << endl;
    cout << "Samples:       " << nSamples << endl;
    cout << "Gene sets:     " << geneSets.size() << endl;
    cout << endl;

    startGSEATime = system_clock::now();

    printTime(startGSEATime);
    cout << " Started GSEA" << endl;

    if (scRna) runScRna();
    else runRna();

    cout << "Elapsed time: " << duration_cast<minutes> (system_clock::now() - startGSEATime).count() << " min" << endl;
}

Gsea::~Gsea() {}
