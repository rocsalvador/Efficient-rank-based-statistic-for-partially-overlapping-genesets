#include "gsea.hh"

void Gsea::readConfig() {
    ifstream file("./gsea.config");
    string aux;
    file >> aux >> expressionMatrixFilename;
    file >> aux >> expressionMatrixSep;
    if (expressionMatrixSep == 't') expressionMatrixSep = '\t';
    file >> aux >> geneSetsFilename;
    file >> aux >> geneSetsSep;
    file >> aux >> outputFilename;
    file >> aux >> outputSep;
    file >> aux >> nThreads;

    if (nThreads == 0) nThreads = thread::hardware_concurrency();

    cout << "################################" << endl;
    cout << "# GSEA config                  #" << endl;
    cout << "################################" << endl;
    cout << "expression-matrix-file: " << expressionMatrixFilename << endl;
    cout << "sep:                    " << expressionMatrixSep << endl;
    cout << "gene-sets-file:         " << geneSetsFilename << endl;
    cout << "sep:                    " << geneSetsSep << endl;
    cout << "output-file:            " << outputFilename << endl;
    cout << "sep:                    " << outputSep << endl;
    cout << "threads-used:           " << nThreads << endl;
    cout << endl;

    file.close();
}

Gsea::Gsea() {
    system_clock::time_point startIOTime = system_clock::now();

    readConfig();
    ifstream file(expressionMatrixFilename);
    string line;

    getline(file, line);
    stringstream ssLine(line);
    string colName;
    uint i = 0;
    while(getline(ssLine, colName, expressionMatrixSep)){
        sampleNames.push_back(colName);
        ++i;
    }
    i = 0;
    while(getline(file, line)) {
        stringstream ssLine(line);

        string rowName;
        getline(ssLine, rowName, expressionMatrixSep);

        string valueStr;
        bool first = true;
        bool nullRow = true;
        while(getline(ssLine, valueStr, expressionMatrixSep)){
            if (first) {
                expressionMatrix.push_back(vector<GeneSample> ());
                first = false;
            }
            double count = stod(valueStr);
            if (nullRow and count != 0) nullRow = false;
            expressionMatrix[i].push_back({rowName, count});
        }
        if (nullRow) {
            expressionMatrix.pop_back();
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

    results = vector<vector<double>> (geneSets.size(), vector<double> (nSamples));

    system_clock::time_point endIOTime = system_clock::now();
    cout << "IO elapsed time: " << duration_cast<milliseconds> (endIOTime - startIOTime).count() / 1000.0 << " s" << endl;
}

Gsea::Gsea(unordered_map<string, unordered_set<string>>& geneSets, vector<vector<GeneSample>>& expressionMatrix) {
    readConfig();

    this->geneSets = geneSets;
    this->expressionMatrix = expressionMatrix;
    nGenes = expressionMatrix.size();
    if (nGenes > 0)
        nSamples = expressionMatrix[0].size();
    results = vector<vector<double>> (geneSets.size(), vector<double> (nSamples));
}

void Gsea::rpm() {
    for (uint j = 0; j < nSamples; ++j) {
        double sum = 0;

        for (uint i = 0; i < nGenes; ++i) 
            sum += expressionMatrix[i][j].count;
        
        double multFactor = 1000000 / sum;
        for (uint i = 0; i < nGenes; ++i) 
            expressionMatrix[i][j].count *= multFactor;
    }
}

void Gsea::meanCenter() {
    for (uint i = 0; i < nGenes; ++i) {
        double mean = 0;
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

void Gsea::sortColumns() {
    for (uint j = 0; j < nSamples; ++j) {
        vector<GeneSample> column = vector<GeneSample> (nGenes);
        for (uint i = 0; i < nGenes; ++i) column[i] = expressionMatrix[i][j];
        sort(column.begin(), column.end(), &Gsea::geneSampleComp);
        for (uint i = 0; i < nGenes; ++i) expressionMatrix[i][j] = column[i];
    }
}

void Gsea::enrichmentScore() {
    auto it = geneSets.begin();
    for (uint k = 0; k < geneSets.size(); ++k) {
        double geneSetSize = it->second.size();
        double posScore = sqrt((nGenes - geneSetSize) / geneSetSize);
        double negScore = - sqrt((geneSetSize / (nGenes - geneSetSize)));
        for (uint j = 0; j < nSamples; ++j) {
            double currentValue = 0;
            double maxValue = 0;
            for (uint i = 0; i < nGenes; ++i) {
                if (it->second.find(expressionMatrix[i][j].geneId) != it->second.end()) 
                    currentValue += posScore;
                else
                    currentValue += negScore;
                if (i == 0) maxValue = currentValue;
                maxValue = max(currentValue, maxValue);
            }
            results[k][j] = maxValue;
        }
        it = next(it);
    }
}

void Gsea::enrichmentScoreJob(uint startSample, uint endSample) {
    assert(endSample < nSamples);
    auto it = geneSets.begin();
    for (uint k = 0; k < geneSets.size(); ++k) {
        double geneSetSize = it->second.size();
        double posScore = sqrt((nGenes - geneSetSize) / geneSetSize);
        double negScore = - sqrt((geneSetSize / (nGenes - geneSetSize)));
        for (uint j = startSample; j < endSample; ++j) {
            double currentValue = 0;
            double maxValue = 0;
            for (uint i = 0; i < nGenes; ++i) {
                if (it->second.find(expressionMatrix[i][j].geneId) != it->second.end()) {
                    currentValue += posScore;
                } else currentValue += negScore;

                if (i == 0) maxValue = currentValue;
                maxValue = max(currentValue, maxValue);
            }
            results[k][j] = maxValue;
        }
        it = next(it);
    }
}

void Gsea::writeResults() {
    ofstream file(outputFilename);
    for (string& sampleName : sampleNames)
        file << outputSep << sampleName;
    file << endl;

    vector<string> geneSetNames = vector<string> (geneSets.size());
    uint i = 0;
    for (auto& it : geneSets) {
        geneSetNames[i] = it.first;
        ++i;
    } 

    i = 0;
    for (vector<double>& row : results) {
        file << geneSetNames[i];
        for (double value : row) {
            file << outputSep << value;
        }
        file << endl;
        ++i;
    }
    file.close();
}

void Gsea::run() {
    cout << "################################" << endl;
    cout << "# GSEA input size              #" << endl;
    cout << "################################" << endl;
    cout << "Sample genes: " << nGenes << endl;
    cout << "Samples:      " << nSamples << endl;
    cout << "Gene sets:    " << geneSets.size() << endl;
    cout << endl;

    system_clock::time_point startGSEATime = system_clock::now();

    rpm();

    meanCenter();

    sortColumns();

    uint samplesPerThread = nSamples / nThreads;
    uint offset = nSamples % samplesPerThread;
    vector<thread> threads = vector<thread> (nThreads);
    for (uint i = 0; i < nThreads; ++i) {
        uint startSample = i * samplesPerThread;
        uint endSample = startSample + samplesPerThread;
        if (i == nThreads - 1) endSample += offset;
        threads[i] = thread(&Gsea::enrichmentScoreJob, this, startSample, endSample);
    }

    for (thread& t : threads) t.join();

    system_clock::time_point endGSEATime = system_clock::now();
    cout << "Elapsed time: " << duration_cast<milliseconds> (endGSEATime - startGSEATime).count() / 1000.0 << " s" << endl;

    writeResults();

    cout << "Results written in " << outputFilename << endl;
}

Gsea::~Gsea() {}