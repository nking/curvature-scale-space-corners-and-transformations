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
//#include "main/c++/algorithms/curves/ParametersKey.h"
//#include "main/c++/algorithms/curves/ParametersEncoder.h"
#include "main/c++/algorithms/curves/MinUnknownUniverseCover.h"

using std::vector;

using gev::MinUnknownUniverseCover;
    
void test0() {

    std::cout << "test0 " ;
         
    MinUnknownUniverseCover *coverCalculator = new MinUnknownUniverseCover();
    
    vector<vector<int> > encodedVariants;
    
    // TODO: fill encoded variants
    
    
    vector<int> coverVariants;
    
    coverCalculator->calculateCover(&encodedVariants, &coverVariants);
    

    // TODO:  assert expected results
    
    delete coverCalculator;
}

int main(int argc, char** argv) {

    test0();
    
    std::cout << " " << std::endl;

    return (EXIT_SUCCESS);
}

#endif
