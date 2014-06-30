/*
 * GEVSimilarityParameters_test.cpp
 *
 *  Created on: April 25, 2014
 *      Author: nichole
 */
#ifndef ALGORITHMS_CURVES_GEVSIMILARITYPARAMETERS_TEST_H
#define ALGORITHMS_CURVES_GEVSIMILARITYPARAMETERS_TEST_H

// for assert
#include <assert.h>
// for SYSTEM_EXIT
#include <cstdlib>
#include <iostream>
#include "main/c++/algorithms/curves/GEVSimilarityParameters.h"

using namespace std;

using gev::GEVSimilarityParameters;

void test0() {

    std::cout << "test0 " ;
    
    GEVSimilarityParameters *params = new GEVSimilarityParameters();
    
    string inFileName = "similar_curve_parameters.txt";
    string outFileName = "sim_curve_params_01.txt";
    
    params->calculateMinSet(inFileName, outFileName);
    
    delete params;
}

int main(int argc, char** argv) {

    std::cout << "GEVSimilarityParameters_test: ";
    
    test0();

    std::cout << " " << std::endl;

    return (EXIT_SUCCESS);
}
#endif
