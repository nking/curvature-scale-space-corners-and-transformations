/* 
 * File:   ParametersEncoder.cpp
 * Author: nichole
 * 
 * Created on April 27, 2014
 */

#include "ParametersEncoder.h"

namespace gev {
    
    ParametersEncoder::ParametersEncoder() {
    }
    
    ParametersEncoder::~ParametersEncoder() {
    }
  
    // TODO:  improve the format of the output file one day
    void ParametersEncoder::readFile(
        vector< unordered_set<int> >* outputEncodedVariables) {
        
        string fileName = "algorithms.curves.GEVSimilarityToolTest-output.txt";
        
        _readFile(fileName, outputEncodedVariables);
    }
    
    void ParametersEncoder::_readFile(string fileName, 
        vector< unordered_set<int> >* outputEncodedVariables) {
        
        string filePath = _getProjectBaseDirectoryPath();
        filePath = filePath.append("/");
        filePath = filePath.append(fileName);
        char buf[256];
        
        FILE *fl = NULL;
        try {
            fl = fopen(filePath.c_str(), "r");
            if (fl != NULL) {
                bool doRead = !feof(fl);
                while (doRead) {
                    
                    // values are the variable numbers
                    unordered_map<ParametersKey, int> varMap;
                                        
                    char *line = fgets(buf, 256, fl);
                    if (line == NULL) {
                        break;
                    }
                    
                    int digit0 = 0;
                    int digitn = strlen(line) - 1;
                    while ((int)line[digit0] < 33) {digit0++;}
                    while ((int)line[digitn] < 33) {digitn--;}
                    
                    int nVars = _getNVarsOfAParameter(line, digit0, digitn);
                    
                    float *k = (float *)malloc(nVars * sizeof(float));
                    float *sigma = (float *)malloc(nVars * sizeof(float));
                    float *mu = (float *)malloc(nVars * sizeof(float));
                    
                    _parseLine(line, digit0, digitn, k, sigma, mu, nVars);
                    
                    unordered_set<int> similar;
                    
                    for (int i = 0; i < nVars; i++) {
                        
                        ParametersKey *key = new ParametersKey(k[i], sigma[i], 
                            mu[i]);
                        
                        unordered_map<gev::ParametersKey, int>::const_iterator 
                            iter = varMap.find(*key);
                        
                        int varNum;
                        if (iter == varMap.end()) {                            
                            
                            varNum = varMap.empty() ? 0 : varMap.size();
                            
                            varMap.insert(
                                make_pair<gev::ParametersKey, int>(*key, varNum));                            
                        
                        } else {
                            
                            varNum = (*iter).second;
                        }
                        
                        similar.insert(varNum);
                        
                        delete key;
                    }
                    
                    outputEncodedVariables->push_back(similar);
                    
                    free(k);
                    free(sigma);
                    free(mu);
                }
            }
        } catch (std::exception e) {
            cerr << "Error: " << e.what() << endl;
            if (fl != NULL) {
                fclose(fl);
                fl = NULL;
            }
            throw;
        }
        if (fl != NULL) {
            fclose(fl);
            fl = NULL;
        }
    }
    
    void ParametersEncoder::writeFile(vector<int>* encodedCoverVariables) {
        
    }
    
    void ParametersEncoder::_writeFile(string fileName, 
        vector<int>* encodedCoverVariables) {
        
    }
    
    string ParametersEncoder::_getProjectBaseDirectoryPath() {
        
        string cwd = _getCWD();
        
        string srch = "GEVSimilarity";
        
        std::size_t found = cwd.find(srch);
        
        if (found == std::string::npos) {
            fprintf(stderr, "expecting current working directory is 'GEVSimilarity'\n");
            fflush(stderr);
            throw EINVAL;
        }
        
        string bDir = cwd.substr(0, found + srch.length());
        printf("base directory = %s\n", bDir.c_str());

        // windows OS handles forward slash?  
        
        string dir = bDir.append("/tmpdata2");

        return dir;
    }
    
    string ParametersEncoder::_getCWD() {
        
        char *wdc = getcwd(NULL, 0);
        string cwd(wdc);
        free(wdc);
        
        return cwd;
    }
    
    int ParametersEncoder::_getNVarsOfAParameter(const char *line, 
        const int digit0, const int digitn) {
        
        int nCommas = 0;
        
        bool foundK = false;
        
        for (int i = digit0; i < digitn; i++) {
            char c = line[i];
            if (foundK) {
                if (c == ',') {
                    nCommas++;
                } else if (c == ' ') {
                    break;
                }
            } else if (c == 'k') {
                foundK = true;
            }
        }
        
        if (!foundK) {
            cerr << "Error reading line: %s" << string(line).c_str() << endl;
            throw EINVAL;
        }
        
        return (nCommas + 1);
    }
    
    void ParametersEncoder::_parseLine(const char *line, const int digit0, 
        const int digitn, float *k, float *sigma, float *mu, const int nVars) {
        
        int nPos = digit0;
        
        for (int i = digit0; i < digitn; i++) {
            
            nPos = _parseLineForNextNVars(line, nPos, digitn, k, nVars);
            
            nPos = _parseLineForNextNVars(line, nPos, digitn, sigma, nVars);
            
            nPos = _parseLineForNextNVars(line, nPos, digitn, mu, nVars);            
        }        
    }
    
    int ParametersEncoder::_parseLineForNextNVars(const char *line, const int digit0, 
        const int digitn, float *a, const int nVars) {
        bool foundEq = false;
        int n = 0;
        int digitF0 = -1;
        int i = -1;
        for (i = digit0; i < digitn; i++) {
            char c = line[i];
            if (foundEq) {
                if (c == ',') {
                    // parse float from digitF0 to digitFn
                    a[n] = _convertToFloat(line, digitF0, i);
                    n++;
                    digitF0 = -1;
                } else if (c == ' ') {
                    break;
                } else if (digitF0 == -1) {
                    digitF0 = i;
                }
            } else if (c == '=') {
                foundEq = true;
            }
        }
        a[n] = _convertToFloat(line, digitF0, i);
        n++;
                
        return i;
    }
    
    float ParametersEncoder::_convertToFloat(const char *line, 
        const int digitBegin, const int digitEnd) {
        
        int len = digitEnd - digitBegin + 1;
        
        char buffer[len];
        
        for (int i = 0; i < len; i++) {
            buffer[i] = line[i + digitBegin];
        }
        
        double d = atof(buffer);
        
        return (float)d;
    }
    
} // end namespace
