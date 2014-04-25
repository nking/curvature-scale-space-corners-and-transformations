/*
 * GEVSimilarityParameters_test.cpp
 *
 *  Created on: April 25, 2014
 *      Author: nichole
 */
#ifndef ALGORITHMS_CURVES_GEVSIMILARITYPARAMETERS_TEST_H
#define ALGORITHMS_CURVES_GEVSIMILARITYPARAMETERS_TEST_H

#include <iostream>
#include "main/c++/algorithms/curves/GEVSimilarityParameters.h"

using namespace std;

using gev::GEVSimilarityParameters;

void test0() {

    std::cout << "test0 " ;
    
    GEVSimilarityParameters *params = new GEVSimilarityParameters();
    
    delete params;
}

int main(int argc, char** argv) {

    test0();

    std::cout << " " << std::endl;

    return (EXIT_SUCCESS);
}
#endif
