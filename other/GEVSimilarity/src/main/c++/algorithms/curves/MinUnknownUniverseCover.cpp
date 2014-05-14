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
        vector< unordered_set<uint32_t> >& inputVariables,
        vector<uint32_t>& outputCoverVariables) {
        
        unordered_map<uint32_t, uint32_t> frequencyMap;
        
        _populateVariableFrequencyMap(inputVariables, frequencyMap);
        
        printf("%lu unique parameter sets\n", frequencyMap.size());
        
        _initializeVariableCover(frequencyMap, outputCoverVariables);
        
        _findMinRepresentativeCover(inputVariables, outputCoverVariables);
        
    }
    
    void MinUnknownUniverseCover::_populateVariableFrequencyMap(
        vector< unordered_set<uint32_t> >& inputVariables,
        unordered_map<uint32_t, uint32_t>& outVariableFrequencyMap) {
                
        for (unsigned long i = 0; i < inputVariables.size(); i++) {
            
            unordered_set<uint32_t> row = inputVariables[i];
            
            for (unordered_set<uint32_t>::const_iterator iter = row.begin();
                iter != row.end(); ++iter) {
                
                int v = *iter;
                
                if (outVariableFrequencyMap.find(v) == outVariableFrequencyMap.end()) {
                    outVariableFrequencyMap.insert(make_pair(v, 1));
                } else {
                    outVariableFrequencyMap[v]++;
                }
            }
        }
    }
    
    void MinUnknownUniverseCover::_initializeVariableCover(
        unordered_map<uint32_t, uint32_t>& inputVariableFrequencyMap, 
        vector<uint32_t>& outputCoverVariables) {
        
        outputCoverVariables.clear();
        vector<uint32_t> counts;
        
        for (unordered_map<uint32_t, uint32_t>::const_iterator iter = 
            inputVariableFrequencyMap.begin();  
            iter != inputVariableFrequencyMap.end(); ++iter) {
            
            int var = iter->first;
            int count = iter->second;
            
            outputCoverVariables.push_back(var);
            counts.push_back(count);
        }
        
        // sort by frequency
        Sorter *sorter = new Sorter();
        sorter->sort(counts, outputCoverVariables);
        delete sorter;
    }
    
    void MinUnknownUniverseCover::_findMinRepresentativeCover(
        vector< unordered_set<uint32_t> >& inputVariables, 
        vector<uint32_t>& outputCoverVariables) {
        
        unsigned long len = outputCoverVariables.size();
        
        unordered_map<uint32_t, uint32_t> lookupIndexesMap;
        for (unsigned long i = 0; i < len; i++) {
            int var = outputCoverVariables[i];
            lookupIndexesMap.insert(make_pair(var, i));
        }
        
        
        unordered_set<uint32_t> chose; 
        
        bool *doRemove = (bool*)calloc(len, sizeof(bool));
        
        for (unsigned long i = 0; i < len; i++) {
            
            if (doRemove[i]) {
                continue;
            }
            
            int var = outputCoverVariables[i];
            
            chose.insert(var);
                        
            // search each row of inputVariables, and if var is in it, mark 
            //   the remaining as doRemove
            for (unsigned long ii = 0; ii < inputVariables.size(); ii++) {
                
                unordered_set<uint32_t> row = inputVariables[ii];
                
                if (row.find(var) == row.end()) {
                    continue;
                }
                
                for (unordered_set<uint32_t>::const_iterator iter = row.begin();
                iter != row.end(); ++iter) {
                
                    int v = *iter;
                                            
                    if ((v != var) && (chose.find(v) == chose.end())) { 
                        int index = lookupIndexesMap.find(v)->second;
                        doRemove[index] = true;
                    }
                }
            }
        }
        
        outputCoverVariables.erase( outputCoverVariables.begin(),
            outputCoverVariables.end());
        
        for (unordered_set<uint32_t>::const_iterator iter = chose.begin();
            iter != chose.end(); ++iter) {
                
            int v = *iter;
            
            outputCoverVariables.push_back(v);
        }
        
        // TODO: when there are more than one variable as the possible choice to 
        //   represent a row, choose the one which is closer to answers already
        //   present in outputCoverVariables.  (closer can be roughly estimated
        //   as smaller difference between the variable numbers.)
    }
    
} // end namespace
