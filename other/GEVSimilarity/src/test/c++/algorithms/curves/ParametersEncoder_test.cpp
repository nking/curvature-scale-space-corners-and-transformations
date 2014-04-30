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
 
void test0_0() {
    
    std::cout << "test0_0 " ;
    
    ParametersEncoder *encoder = new ParametersEncoder();
    
    float expected = -0.0011764728;
    string lineStr = string("-0.0011764728");
    const char* line = lineStr.c_str();
    
    float parsed = encoder->_convertToFloat(line, 0, lineStr.length() - 1);
    
    float diff = expected - parsed;
    if (diff < 0) {
        diff *= -1;
    }
    float eps = expected * 0.1;
    if (eps < 0) {
        eps *= -1;
    }
    
    assert(diff < eps);
    
    delete encoder;
}

void test0_1() {
    
    std::cout << "test0_1 " ;
    
    ParametersEncoder *encoder = new ParametersEncoder();
    
    float expected = -7.4505806E-9;
    string lineStr = string("-7.4505806E-9");
    const char *line = lineStr.c_str();    
    
    float parsed = encoder->_convertToFloat(line, 0, lineStr.length() - 1);
    
    float diff = expected - parsed;
    if (diff < 0) {
        diff *= -1;
    }
    float eps = expected * 0.1;
    if (eps < 0) {
        eps *= -1;
    }
    
    assert(diff < eps);
    
    delete encoder;
}

void test0_2() {
    
    std::cout << "test0_2 " ;
    
    ParametersEncoder *encoder = new ParametersEncoder();
    
    float expected = 9.000001;
    string lineStr = string("9.000001");
    const char *line = lineStr.c_str();    
    
    float parsed = encoder->_convertToFloat(line, 0, lineStr.length() - 1);
    
    float diff = expected - parsed;
    if (diff < 0) {
        diff *= -1;
    }
    float eps = expected * 0.1;
    if (eps < 0) {
        eps *= -1;
    }
    
    assert(diff < eps);
    
    delete encoder;
}

void test1_0() {
    
    std::cout << "test1_0 " ;
    
    string lineStr = "k=0.010000001,0.020000001,0.020000001,0.020000001,0.020000001";
    lineStr.append(",0.020000001,0.020000001  sigma=0.099999994,0.099999994,0.09");
    lineStr.append("9999994,0.099999994,0.099999994,0.099999994,0.099999994  ");
    lineStr.append("mu=-0.010000001,-0.008571431,-0.009090915,-0.009090915,");
    lineStr.append("-0.008571431,-0.009090915,-0.009090915");
      
    const float expected[] = {0.010000001,0.020000001,0.020000001,0.020000001,
        0.020000001,0.020000001,0.020000001};
    
    ParametersEncoder *encoder = new ParametersEncoder();
    
    const char *line = lineStr.c_str();
    
    int digit0 = 0;
    int digitn = lineStr.find("  sigma");
    int nVars = 7;
    float *array = (float*)malloc(nVars * sizeof(float));
    
    encoder->_parseLineForNextNVars(line, digit0, digitn, array, nVars);
    
    for (int i = 0; i < nVars; i++) {
        float diff = expected[i] - array[i];
        if (diff < 0) {
            diff *= -1;
        }
        float eps = expected[i] * 0.1;
        assert(diff < eps);
    }
        
    free(array);
    
    delete encoder;
}

void test1_1() {
    
    std::cout << "test1_1 " ;
    
    string lineStr = "k=0.010000001,0.020000001,0.020000001,0.020000001,0.020000001";
    lineStr.append(",0.020000001,0.020000001  sigma=0.099999994,0.099999994,0.09");
    lineStr.append("9999994,0.099999994,0.099999994,0.099999994,0.099999994  ");
    lineStr.append("mu=-0.010000001,-0.008571431,-0.009090915,-0.009090915,");
    lineStr.append("-0.008571431,-0.009090915,-0.009090915");
      
    const float expected[] = {0.099999994,0.099999994,0.099999994,
        0.099999994,0.099999994,0.099999994,0.099999994};
    
    ParametersEncoder *encoder = new ParametersEncoder();
    
    const char *line = lineStr.c_str();
    
    int digit0 = lineStr.find(" sigma");
    int digitn = lineStr.find("mu") - 2;
    
    int nVars = 7;
    float *array = (float*)malloc(nVars * sizeof(float));
    
    encoder->_parseLineForNextNVars(line, digit0, digitn, array, nVars);
    
    for (int i = 0; i < nVars; i++) {
        float diff = expected[i] - array[i];
        if (diff < 0) {
            diff *= -1;
        }
        float eps = expected[i] * 0.1;
        assert(diff < eps);
    }
        
    free(array);
    
    delete encoder;
}

void test2_0() {
    
    std::cout << "test2_0 " ;
    
    string lineStr = "k=0.010000001,0.020000001,0.020000001,0.020000001,0.020000001";
    lineStr.append(",0.020000001,0.020000001  sigma=0.099999994,0.099999994,0.09");
    lineStr.append("9999994,0.099999994,0.099999994,0.099999994,0.099999994  ");
    lineStr.append("mu=-0.010000001,-0.008571431,-0.009090915,-0.009090915,");
    lineStr.append("-0.008571431,-0.009090915,-0.009090915");
      
    const float kExpected[] = {0.010000001,0.020000001,0.020000001,0.020000001,
        0.020000001,0.020000001,0.020000001};
    const float sigmaExpected[] = {0.099999994,0.099999994,0.099999994,
        0.099999994,0.099999994,0.099999994,0.099999994};
    const float muExpected[] = {-0.010000001, -0.008571431, -0.009090915,
        -0.009090915, -0.008571431, -0.009090915, -0.009090915};
    
    ParametersEncoder *encoder = new ParametersEncoder();
    
    const char *line = lineStr.c_str();
    
    int digit0 = lineStr.find(" x-mu");
    int digitn = lineStr.length();
    
    int nVars = 7;
    float *k = (float*)malloc(nVars * sizeof(float));
    float *sigma = (float*)malloc(nVars * sizeof(float));
    float *mu = (float*)malloc(nVars * sizeof(float));
    
    int foundNVars = encoder->_getNVarsOfAParameter(line, digit0, digitn);
    assert(foundNVars == nVars);
    
    encoder->_parseLine(line, digit0, digitn, k, sigma, mu, nVars);
    
    for (int i = 0; i < nVars; i++) {
        float diff = kExpected[i] - k[i];
        if (diff < 0) {
            diff *= -1;
        }
        float eps = kExpected[i] * 0.1;
        assert(diff < eps);
    }
    for (int i = 0; i < nVars; i++) {
        float diff = sigmaExpected[i] - sigma[i];
        if (diff < 0) {
            diff *= -1;
        }
        float eps = sigmaExpected[i] * 0.1;
        assert(diff < eps);
    }
    for (int i = 0; i < nVars; i++) {
        float diff = muExpected[i] - mu[i];
        if (diff < 0) {
            diff *= -1;
        }
        float eps = muExpected[i] * 0.1;
        assert(diff < eps);
    }
    
    free(k);
    free(sigma);
    free(mu);
    
    delete encoder;
}

void test3() {
    
    std::cout << "test3 " ;
    
    ParametersEncoder *encoder = new ParametersEncoder();
    
    string cwd = encoder->_getCWD();
    assert(cwd.length() > 0);
    
    string path = encoder->_getProjectBaseDirectoryPath();
    
    assert(path.find("tmpdata2") != string::npos);
    
    delete encoder;
}

void test4() {

    std::cout << "test4 " ;
         
    ParametersEncoder *encoder = new ParametersEncoder();
    
    vector< unordered_set<int> > encodedVariables;
    
    // TODO: put test file name here
    string inFileName = "";
    
    // TODO: test readFile
    encoder->_readFile(inFileName, &encodedVariables);
    
    // TODO: assert results
    
    
    // TODO: test writeFile
    string outFileName = "";
    vector<int> encodedCover;
    
    encoder->_writeFile(outFileName, &encodedCover);
    
    // TODO: assert results
    
    delete encoder;
}

int main(int argc, char** argv) {

    std::cout << "ParametersEncoder_test: ";
    
    test0_0();
    test0_1();    
    test0_2();
    
    test1_0();
    test1_1();
    
    //test3();
    
    //test4();
    
    
    std::cout << " " << std::endl;

    return (EXIT_SUCCESS);
}

#endif
