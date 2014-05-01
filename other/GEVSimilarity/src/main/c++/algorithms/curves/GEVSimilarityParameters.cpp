/* 
 * File:   GEVSimilarityParameters.cpp
 * Author: nichole
 * 
 * Created on April 25, 2014
 */

#include "GEVSimilarityParameters.h"

namespace gev {
    
    GEVSimilarityParameters::GEVSimilarityParameters() {
    }
    
    GEVSimilarityParameters::~GEVSimilarityParameters() {
    }
    
    void GEVSimilarityParameters::calculateMinSet() {
        calculateMinSet(string(""), string(""));
    }
    
    void GEVSimilarityParameters::calculateMinSet(string inputFileName, 
        string outputFileName) {        
         
        vector< unordered_set<uint32_t> > encodedVariables;
        vector<uint32_t> outputCoverVariables;
        
        ParametersEncoder *encoder = new ParametersEncoder();
    
        
        if (inputFileName.length() > 0) {
            encoder->_readFile(inputFileName, encodedVariables);
        } else {
            encoder->readFile(encodedVariables);
        }
        
        
        MinUnknownUniverseCover *coverCalculator = new MinUnknownUniverseCover();
                
        coverCalculator->calculateCover(encodedVariables, outputCoverVariables);
        
        delete coverCalculator;
        
        
        if (outputFileName.length() > 0) {
            encoder->_writeFile(outputFileName, outputCoverVariables);
        } else {
            encoder->writeFile(outputCoverVariables);
        }
        
        delete encoder;
    }
    
    
    
} // end namespace
