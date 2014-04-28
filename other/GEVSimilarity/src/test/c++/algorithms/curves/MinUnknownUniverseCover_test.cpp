/*
 * MinUnknownUniverseCover_test.cpp
 *
 *  Created on: April 27, 2014
 *      Author: nichole
 */
#ifndef ALGORITHMS_CURVES_MINUKNOWNUNIVERSECOVER_TEST_H
#define ALGORITHMS_CURVES_MINUKNOWNUNIVERSECOVER_TEST_H

// for assert
#include <assert.h>
// time is in time.h so can include that or ctime
#include <time.h>
// srandom and random are in  stdlib.h so can include that or cstdlib
#include <stdlib.h>
// for int32_t
#include <stdint.h>
// for cout
#include <iostream>
#include <vector>
#include <tr1/unordered_map>
//#include "main/c++/algorithms/curves/ParametersKey.h"
//#include "main/c++/algorithms/curves/ParametersEncoder.h"
#include "main/c++/algorithms/curves/MinUnknownUniverseCover.h"

using std::vector;
using std::tr1::unordered_map;
using gev::MinUnknownUniverseCover;
    
void addVector(vector< vector<int> >* a) {
    vector<int> b;
    a->push_back(b);
}

void test0() {

    std::cout << "test0 " ;
         
    MinUnknownUniverseCover *coverCalculator = new MinUnknownUniverseCover();
    
    vector<vector<int> > encodedVariants;
    
    /*
     0 1  2  3
       1
               4 5 6
                 5
                     7  8
     */
    addVector(&encodedVariants);
    encodedVariants[0].push_back(0);
    encodedVariants[0].push_back(1);
    encodedVariants[0].push_back(2);
    encodedVariants[0].push_back(3);
    addVector(&encodedVariants);
    encodedVariants[1].push_back(1);
    addVector(&encodedVariants);
    encodedVariants[2].push_back(4);
    encodedVariants[2].push_back(5);
    encodedVariants[2].push_back(6);
    addVector(&encodedVariants);
    encodedVariants[3].push_back(5);
    addVector(&encodedVariants);
    encodedVariants[4].push_back(7);
    encodedVariants[4].push_back(8);
    
    
    unordered_map<int, int> frequencyMap;
    
    coverCalculator->_populateVariableFrequencyMap(&encodedVariants,
        &frequencyMap);
    
    assert(frequencyMap.size() == 9);
    assert(frequencyMap[0] == 0);
    assert(frequencyMap[1] == 2);
    assert(frequencyMap[2] == 1);
    assert(frequencyMap[3] == 1);
    assert(frequencyMap[4] == 1);
    assert(frequencyMap[5] == 2);
    assert(frequencyMap[6] == 1);
    assert(frequencyMap[7] == 1);
    assert(frequencyMap[8] == 1);
    
    
    
    vector<int> outputCoverVariables;
    
    coverCalculator->_initializeVariableCover(&frequencyMap,
        &outputCoverVariables);
    
    assert(outputCoverVariables.size() == 9);
    assert(outputCoverVariables[0] == 1 || outputCoverVariables[0] == 5);
    if (outputCoverVariables[0] == 1) {
        assert(outputCoverVariables[1] == 5);
    } else {
        assert(outputCoverVariables[1] == 1);
    }
    assert(outputCoverVariables[8] == 0);
    assert(
        outputCoverVariables[2] == 2 || outputCoverVariables[2] == 3 ||
        outputCoverVariables[2] == 4 || outputCoverVariables[2] == 6 ||
        outputCoverVariables[2] == 7 || outputCoverVariables[2] == 8
    );

    
    
    coverCalculator->_findMinRepresentativeCover(&encodedVariants,
        &frequencyMap, &outputCoverVariables);
    
    // 0, 1, 5, 7 || 8
    assert(outputCoverVariables.size() == 4);
    
    bool *found = (bool*)calloc(4, sizeof(bool));
    for (unsigned long i = 0; i < outputCoverVariables.size(); i++) {
        int var = outputCoverVariables[i];
        switch(var) {
            case 0:
                found[0] = true;
                break;
            case 1:
                found[1] = true;
                break;
            case 5:
                found[2] = true;
                break;
            case 7:
            case 8:
                found[3] = true;
                break;
            default:
                break;
        }
    }
    for (unsigned long i = 0; i < 4; i++) {
        assert(found[i]);
    }
    
    free(found);
    
    delete coverCalculator;
}

int main(int argc, char** argv) {

    std::cout << "MinUnknownUniverseCover_test: ";
    
    test0();
    
    std::cout << " " << std::endl;

    return (EXIT_SUCCESS);
}

#endif
