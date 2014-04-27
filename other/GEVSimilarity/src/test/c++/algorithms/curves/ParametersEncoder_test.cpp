/*
 * ParametersEncoder_test.cpp
 *
 *  Created on: April 27, 2014
 *      Author: nichole
 */
#ifndef ALGORITHMS_CURVES_PARAMETERSENCODER_TEST_H
#define ALGORITHMS_CURVES_PARAMETERSENCODER_TEST_H

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
#include "main/c++/algorithms/curves/ParametersKey.h"
#include "main/c++/algorithms/curves/ParametersEncoder.h"

using namespace std;

using std::tr1::unordered_map;

using namespace gev;
    
void test0() {

    std::cout << "test0 " ;
         
    ParametersEncoder *encoder = new ParametersEncoder();
    
    vector<vector<int> > encodedVariants;
    
    // TODO: put test file name here
    string inFileName = "";
    
    // TODO: test readFile
    encoder->_readFile(inFileName, &encodedVariants);
    
    // TODO: assert results
    
    
    // TODO: test writeFile
    string outFileName = "";
    vector<int> encodedCover;
    
    encoder->_writeFile(outFileName, &encodedCover);
    
    // TODO: assert results
    
    delete encoder;
}

int main(int argc, char** argv) {

    test0();
    
    std::cout << " " << std::endl;

    return (EXIT_SUCCESS);
}

#endif
