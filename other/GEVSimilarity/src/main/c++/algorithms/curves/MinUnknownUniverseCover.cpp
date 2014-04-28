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
        const vector<vector<int> >* inputVariables,
        vector<int>* outputCoverVariables) {
        
        unordered_map<int, int> frequencyMap;
        
        _populateVariableFrequencyMap(inputVariables, &frequencyMap);
        
        _initializeVariableCover(&frequencyMap, outputCoverVariables);
        
        // note that frequencyMap gets modified in next method
        _findMinRepresentativeCover(inputVariables, &frequencyMap, outputCoverVariables);
        
    }
    
    void MinUnknownUniverseCover::_populateVariableFrequencyMap(
        const vector<vector<int> >* inputVariables,
        unordered_map<int, int> *outVariableFrequencyMap) {
                
        for (unsigned long i = 0; i < inputVariables->size(); i++) {
            
            vector<int> row = (*inputVariables)[i];
            
            for (unsigned long j = 0; j < row.size(); j++) {
                
                int v = row[j];
                
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
        
    }
    
    void MinUnknownUniverseCover::_findMinRepresentativeCover(
        const vector<vector<int> >* inputVariables, 
        unordered_map<int, int> *outVariableFrequencyMap,
        vector<int>* outputCoverVariables) {
        
        // TODO: when there are more than one variable as the possible choice to 
        //   represent a row, choose the one which is closer to answers already
        //   present in outputCoverVariables.  (closer can be roughly estimated
        //   as smaller difference between the variable numbers.)
    }
    
} // end namespace
