/* 
 * File:   MinUnknownUniverseCover.cpp
 * Author: nichole
 * 
 * Created on April 27, 2014
 */

#include "MinUnknownUniverseCover.h"

namespace gev {
    
    MinUnknownUniverseCover::MinUnknownUniverseCover() {
    }
    
    MinUnknownUniverseCover::~MinUnknownUniverseCover() {
    }
        
    void MinUnknownUniverseCover::calculateCover(
        const vector< unordered_set<int> >* inputVariables,
        vector<int>* outputCoverVariables) {
        
        unordered_map<int, int> frequencyMap;
        
        _populateVariableFrequencyMap(inputVariables, &frequencyMap);
        
        _initializeVariableCover(&frequencyMap, outputCoverVariables);
        
        _findMinRepresentativeCover(inputVariables, outputCoverVariables);
        
    }
    
    void MinUnknownUniverseCover::_populateVariableFrequencyMap(
        const vector< unordered_set<int> >* inputVariables,
        unordered_map<int, int> *outVariableFrequencyMap) {
                
        for (unsigned long i = 0; i < inputVariables->size(); i++) {
            
            unordered_set<int> row = (*inputVariables)[i];
            
            for (unordered_set<int>::const_iterator iter = row.begin();
                iter != row.end(); ++iter) {
                
                int v = *iter;
                
                if (outVariableFrequencyMap->find(v) == outVariableFrequencyMap->end()) {
                    outVariableFrequencyMap->insert(make_pair(v, 1));
                } else {
                    (*outVariableFrequencyMap)[v]++;
                }
            }
        }
    }
    
    void MinUnknownUniverseCover::_initializeVariableCover(
        const unordered_map<int, int> *inputVariableFrequencyMap, 
        vector<int>* outputCoverVariables) {
        
        outputCoverVariables->clear();
        vector<int> counts;
        
        for (unordered_map<int, int>::const_iterator iter = 
            inputVariableFrequencyMap->begin();  
            iter != inputVariableFrequencyMap->end(); ++iter) {
            
            int var = iter->first;
            int count = iter->second;
            
            outputCoverVariables->push_back(var);
            counts.push_back(count);
        }
        
        // sort by frequency
        Sorter *sorter = new Sorter();
        sorter->sort(&counts, outputCoverVariables);
        delete sorter;
    }
    
    void MinUnknownUniverseCover::_findMinRepresentativeCover(
        const vector< unordered_set<int> >* inputVariables, 
        vector<int>* outputCoverVariables) {
        
        unsigned long len = outputCoverVariables->size();
        
        unordered_map<int, int> lookupIndexesMap;
        for (unsigned long i = 0; i < len; i++) {
            int var = (*outputCoverVariables)[i];
            lookupIndexesMap.insert(make_pair(var, i));
        }
        
        
        unordered_set<int> chose;
        
        bool *doRemove = (bool*)calloc(len, sizeof(bool));
        
        for (unsigned long i = 0; i < len; i++) {
            
            if (doRemove[i]) {
                continue;
            }
            
            int var = (*outputCoverVariables)[i];
            
            chose.insert(var);
                        
            // search each row of inputVariables, and if var is in it, mark 
            //   the remaining as doRemove
            for (unsigned long ii = 0; ii < inputVariables->size(); ii++) {
                
                unordered_set<int> row = (*inputVariables)[ii];
                
                if (row.find(var) == row.end()) {
                    continue;
                }
                
                for (unordered_set<int>::const_iterator iter = row.begin();
                iter != row.end(); ++iter) {
                
                    int v = *iter;
                                            
                    if ((v != var) && (chose.find(v) == chose.end())) { 
                        int index = lookupIndexesMap.find(v)->second;
                        doRemove[index] = true;
                    }
                }
            }
        }
        
        vector<int> keep;
        
        for (unsigned long i = 0; i < outputCoverVariables->size(); i++) {
            if (!doRemove[i]) {
                int var = (*outputCoverVariables)[i];
                keep.push_back(var);
            }
        }
        
        free(doRemove);
        
        outputCoverVariables->swap(keep);
        
        // TODO: when there are more than one variable as the possible choice to 
        //   represent a row, choose the one which is closer to answers already
        //   present in outputCoverVariables.  (closer can be roughly estimated
        //   as smaller difference between the variable numbers.)
    }
    
} // end namespace
