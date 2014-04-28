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
    
void addSet(vector< unordered_set<int> >* a) {
    unordered_set<int> b;
    a->push_back(b);
}

void test0() {

    std::cout << "test0 " ;
         
    MinUnknownUniverseCover *coverCalculator = new MinUnknownUniverseCover();
    
    vector< unordered_set<int> > encodedVariants;
    
    /*
     0 1  2  3
       1
               4 5 6
                 5
                     7  8
     */
    addSet(&encodedVariants);
    encodedVariants[0].insert(0);
    encodedVariants[0].insert(1);
    encodedVariants[0].insert(2);
    encodedVariants[0].insert(3);
    addSet(&encodedVariants);
    encodedVariants[1].insert(1);
    addSet(&encodedVariants);
    encodedVariants[2].insert(4);
    encodedVariants[2].insert(5);
    encodedVariants[2].insert(6);
    addSet(&encodedVariants);
    encodedVariants[3].insert(5);
    addSet(&encodedVariants);
    encodedVariants[4].insert(7);
    encodedVariants[4].insert(8);
    
    
    unordered_map<int, int> frequencyMap;
    
    coverCalculator->_populateVariableFrequencyMap(&encodedVariants,
        &frequencyMap);
    
    assert(frequencyMap.size() == 9);
    assert(frequencyMap[0] == 1);
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
    assert(
        outputCoverVariables[2] == 0 ||
        outputCoverVariables[2] == 2 || outputCoverVariables[2] == 3 ||
        outputCoverVariables[2] == 4 || outputCoverVariables[2] == 6 ||
        outputCoverVariables[2] == 7 || outputCoverVariables[2] == 8
    );

    
    
    coverCalculator->_findMinRepresentativeCover(&encodedVariants,
        &outputCoverVariables);
        
    // 1, 5, 7 || 8
    assert(outputCoverVariables.size() == 3);
    
    bool *found = (bool*)calloc(3, sizeof(bool));
    for (unsigned long i = 0; i < outputCoverVariables.size(); i++) {
        int var = outputCoverVariables[i];
        switch(var) {
            case 1:
                found[0] = true;
                break;
            case 5:
                found[1] = true;
                break;
            case 7:
            case 8:
                found[2] = true;
                break;
            default:
                break;
        }
    }
    for (unsigned long i = 0; i < 3; i++) {
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
