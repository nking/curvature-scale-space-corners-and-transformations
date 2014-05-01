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
        vector< unordered_set<uint32_t> >* outputEncodedVariables) {
        
        string fileName = "algorithms.curves.GEVSimilarityToolTest-output.txt";
        
        _readFile(fileName, outputEncodedVariables);
    }
    
    void ParametersEncoder::_readFile(string fileName, 
        vector< unordered_set<uint32_t> >* outputEncodedVariables) {
        
        string filePath = _getProjectTmpDirectoryPath();
        filePath = filePath.append("/");
        filePath = filePath.append(fileName);
        uint32_t buffSz = 100*2048;
        char buf[buffSz];
        
        uint32_t count = 0;
        
        FILE *fl = NULL;
        try {
            fl = fopen(filePath.c_str(), "r");
            if (fl != NULL) {
                bool doRead = !feof(fl);
                while (doRead) {
                    
                    char *line = fgets(buf, buffSz, fl);
                    if (line == NULL) {
                        break;
                    }
                                        
                    int digit0 = 0;
                    int digitn = strlen(line) - 1;
                    while ((uint32_t)line[digit0] < 33) {digit0++;}
                    while ((uint32_t)line[digitn] < 33) {digitn--;}
     
                    uint32_t nVars = _getNVarsOfAParameter(line, digit0, digitn);
                     
                    float *k = (float *)malloc(nVars * sizeof(float));
                    float *sigma = (float *)malloc(nVars * sizeof(float));
                    float *mu = (float *)malloc(nVars * sizeof(float));
                    
                    _parseLine(line, digit0, digitn, k, sigma, mu, nVars);
                    
                    unordered_set<uint32_t> similar;
                    
                    for (uint32_t i = 0; i < nVars; i++) {
                        
                        ParametersKey *key = new ParametersKey(k[i], sigma[i], 
                            mu[i]);
                        
                        unordered_map<gev::ParametersKey, uint32_t>::const_iterator 
                            iter = variableMap.find(*key);
                        
                        uint32_t varNum;
                        if (iter == variableMap.end()) {                            
                            
                            varNum = variableMap.empty() ? 0 : variableMap.size();
                            
                            variableMap.insert(
                                make_pair<gev::ParametersKey, uint32_t>(*key, varNum));
   
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
                    
                    count++;
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
    
    void ParametersEncoder::writeFile(vector<uint32_t>* encodedCoverVariables) {
        
        string fileName = "gev_similarity_min_ranges-output.txt";
        
        _writeFile(fileName, encodedCoverVariables);
    }
    
    void ParametersEncoder::_writeFile(string fileName, 
        vector<uint32_t>* encodedCoverVariables) {
        
        string filePath = _getProjectTmpDirectoryPath();
        filePath = filePath.append("/");
        filePath = filePath.append(fileName);
        
        // create a reverse map to lookup vars
        unordered_map<uint32_t, gev::ParametersKey> revLookupMap;
        
        for (unordered_map<gev::ParametersKey, uint32_t>::const_iterator 
            iter = variableMap.begin(); iter != variableMap.end(); ++iter) {
             
            revLookupMap.insert( make_pair(iter->second, iter->first));
        }
        
        
        FILE *fl = NULL;
        try {
            fl = fopen(filePath.c_str(), "w");
            if (fl != NULL) {
                
                fprintf(fl, "#k        sigma        mu\n");
                fflush(fl);
                
                for (unsigned long i = 0; i < encodedCoverVariables->size(); i++) {
                
                    uint32_t varNum = (*encodedCoverVariables)[i];
                    
                    ParametersKey key = revLookupMap.find(varNum)->second;
                    
                    float k;
                    float sigma;
                    float mu;
                    
                    key.decodeToParams(&k, &sigma, &mu);
                    
                    fprintf(fl, "%e  %e  %e\n", k, sigma, mu);

                    if (i % 100 == 0) {
                        fflush(fl);
                    }
                }
                printf("wrote to outfile: %s", filePath.c_str());
            } else {
                printf("ERROR: unable to open outfile: %s", filePath.c_str());
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
    
    string ParametersEncoder::_getProjectTmpDirectoryPath() {
        
        string cwd = _getCWD();
        
        string srch = "GEVSimilarity";
        
        std::size_t found = cwd.find(srch);
        
        if (found == std::string::npos) {
            fprintf(stderr, "expecting current working directory is 'GEVSimilarity'\n");
            fflush(stderr);
            throw EINVAL;
        }
        
        string bDir = cwd.substr(0, found + srch.length());
        //printf("base directory = %s\n", bDir.c_str());

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
    
    uint32_t ParametersEncoder::_getNVarsOfAParameter(const char *line, 
        const uint32_t digit0, const uint32_t digitn) {
        
        uint32_t nCommas = 0;
        
        bool foundK = false;
        
        for (uint32_t i = digit0; i < digitn; i++) {
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
    
    void ParametersEncoder::_parseLine(const char *line, const uint32_t digit0, 
        const uint32_t digitn, float *k, float *sigma, float *mu, 
        const uint32_t nVars) {
        
        uint32_t nPos = digit0;
                    
        nPos = _parseLineForNextNVars(line, nPos, digitn, k, nVars);

        nPos = _parseLineForNextNVars(line, nPos, digitn, sigma, nVars);

        nPos = _parseLineForNextNVars(line, nPos, digitn, mu, nVars);            
    }
    
    uint32_t ParametersEncoder::_parseLineForNextNVars(const char *line, 
        const uint32_t digit0, const uint32_t digitn, float *a, 
        const uint32_t nVars) {
        
        bool foundEq = false;
        uint32_t sentinel = -1;
        uint32_t n = 0;
        uint32_t digitF0 = sentinel;
        uint32_t i = sentinel;
        for (i = digit0; i < digitn; i++) {
            char c = line[i];
            if (foundEq) {
                if (c == ',') {
                    // parse float from digitF0 to digitFn
                    a[n] = _convertToFloat(line, digitF0, i);
                    n++;
                    digitF0 = sentinel;
                } else if (c == ' ') {
                    break;
                } else if (digitF0 == sentinel) {
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
        const uint32_t digitBegin, const uint32_t digitEnd) {
        
        uint32_t len = digitEnd - digitBegin + 1;
        
        char buffer[len];
        
        for (uint32_t i = 0; i < len; i++) {
            buffer[i] = line[i + digitBegin];
        }
        
        double d = atof(buffer);
        
        return (float)d;
    }
    
} // end namespace
