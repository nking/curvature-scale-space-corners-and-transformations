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
// for uint32_t
#include <stdint.h>
// for remove
#include <stdio.h>
// for cout
#include <iostream>
#include <vector>
#include "main/c++/algorithms/curves/ParametersKey.h"
#include "main/c++/algorithms/curves/ParametersEncoder.h"

using namespace std;

using std::tr1::unordered_map;

using namespace gev;
 
void test0() {
    
    std::cout << "test0 " ;
    
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

void test1() {
    
    std::cout << "test1 " ;
    
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

void test2() {
    
    std::cout << "test2 " ;
    
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

void test3() {
    
    std::cout << "test3 " ;
    
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
    int digitn = lineStr.length();
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

void test6() {
    
    std::cout << "test6 " ;
    
    string lineStr =  "k=0.010000001,0.010000001,0.010000001,0.010000001,0.010000001,0.010000001,0.010000001,0.010000001,0.010000001,0.010000001,0.010000001,0.010000001,0.015000001,0.015000001,0.015000001,0.015000001,0.015000001,0.015000001,0.020000001,0.020000001,0.020000001,0.020000001,0.020000001,0.020000001,0.025000002,0.025000002,0.025000002,0.025000002,0.025000002,0.025000002,0.030000001,0.030000001,0.030000001,0.030000001,0.030000001,0.030000001,0.035000004,0.035000004,0.035000004,0.035000004,0.035000004,0.035000004,0.040000003,0.040000003,0.040000003,0.040000003,0.040000003,0.040000003,0.045,0.045,0.045,0.045,0.045,0.045,0.050000004,0.050000004,0.050000004,0.050000004,0.050000004,0.050000004,0.055000003,0.055000003,0.055000003,0.055000003,0.055000003,0.055000003,0.060000002,0.060000002,0.060000002,0.060000002,0.060000002,0.060000002,0.060000002,0.060000002,0.065000005,0.065000005,0.065000005,0.065000005,0.07000001,0.07000001,0.07000001,0.07000001,0.075,0.075,0.075,0.075,0.080000006,0.080000006,0.080000006,0.080000006,0.08500001,0.08500001,0.08500001,0.08500001,0.09,0.09,0.09,0.09,0.095000006,0.095000006,0.095000006,0.095000006,0.10000001,0.10000001,0.10000001,0.10000001,0.10000001,0.10000001,0.10000001,0.10000001,0.15,0.15,0.15,0.15,0.20000002,0.20000002,0.25000003,0.25000003,0.3,0.3,0.35000002,0.35000002  sigma=0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994,0.099999994  mu=0.010000001,0.010000001,0.015000001,0.010000001,0.010000001,0.015000001,0.010000001,0.010000001,0.015000001,0.010000001,0.010000001,0.015000001,0.010000001,0.010000001,0.015000001,0.010000001,0.010000001,0.015000001,0.010000001,0.010000001,0.015000001,0.010000001,0.010000001,0.015000001,0.010000001,0.010000001,0.015000001,0.010000001,0.010000001,0.015000001,0.010000001,0.010000001,0.015000001,0.010000001,0.010000001,0.015000001,0.010000001,0.010000001,0.015000001,0.010000001,0.010000001,0.015000001,0.010000001,0.010000001,0.015000001,0.010000001,0.010000001,0.015000001,0.010000001,0.010000001,0.015000001,0.010000001,0.010000001,0.015000001,0.010000001,0.010000001,0.015000001,0.010000001,0.010000001,0.015000001,0.010000001,0.010000001,0.015000001,0.010000001,0.010000001,0.015000001,0.010000001,0.010000001,0.015000001,0.020000001,0.010000001,0.010000001,0.015000001,0.020000001,0.015000001,0.020000001,0.015000001,0.020000001,0.015000001,0.020000001,0.015000001,0.020000001,0.015000001,0.020000001,0.015000001,0.020000001,0.015000001,0.020000001,0.015000001,0.020000001,0.015000001,0.020000001,0.015000001,0.020000001,0.015000001,0.020000001,0.015000001,0.020000001,0.015000001,0.020000001,0.015000001,0.020000001,0.015000001,0.020000001,0.015000001,0.020000001,0.015000001,0.020000001,0.015000001,0.020000001,0.020000001,0.025000002,0.020000001,0.025000002,0.030000001,0.030000001,0.035000004,0.035000004,0.040000003,0.040000003,0.045,0.045";
      
    ParametersEncoder *encoder = new ParametersEncoder();
    
    const char *line = lineStr.c_str();
    
    int digit0 = 0;
    int digitn = lineStr.length();
    
    int nVars = 122;
    
    int foundNVars = encoder->_getNVarsOfAParameter(line, digit0, digitn);
    assert(foundNVars == nVars);
   
    delete encoder;
}

void test7() {
    
    std::cout << "test7 " ;
    
    ParametersEncoder *encoder = new ParametersEncoder();
    
    string cwd = encoder->_getCWD();
    assert(cwd.length() > 0);
    
    string path = encoder->_getProjectTmpDirectoryPath();
    
    assert(path.find("tmpdata2") != string::npos);
    
    delete encoder;
}

void test8() {

    std::cout << "test8 " ;
         
    ParametersEncoder *encoder = new ParametersEncoder();
    
    vector< unordered_set<uint32_t> > encodedVariables;
    
    string inFileName = "test-data.txt";
    
    encoder->_readFile(inFileName, encodedVariables);
    
    assert(encodedVariables.size() == 1435);
        
    // spot checks for expected:
    assert(encodedVariables[0].size() <= 122);
    
    assert(encodedVariables[1434].size() <= 1248);
    
    delete encoder;
}

void test9() {

    std::cout << "test9 " ;
         
    ParametersEncoder *encoder = new ParametersEncoder();
    
    vector< unordered_set<uint32_t> > encodedVariables;
    
    string inFileName = "test-data.txt";
    
    encoder->_readFile(inFileName, encodedVariables);
    
    
    string outFileName = "delete.txt";
    string filePath = encoder->_getProjectTmpDirectoryPath();
    filePath = filePath.append("/");
    filePath = filePath.append(outFileName);
    
    std::remove(filePath.c_str());
    
    vector<uint32_t> encodedCover;
    encodedCover.push_back(1);
    encodedCover.push_back(20);
    encodedCover.push_back(30);
    
    encoder->_writeFile(outFileName, encodedCover);
    
    
    FILE *fl = NULL;
    try {
        fl = fopen(filePath.c_str(), "r");
        if (fl != NULL) {
            
            assert(!feof(fl));
            long buffSz = 256;
            char buf[buffSz];
            
            char *line = fgets(buf, buffSz, fl);
            
            int count = 0;
            
            while (line != NULL) {
                assert(strlen(line) > 3);
                line = fgets(buf, buffSz, fl);
                count++;
            }
            assert(count == 4);
        }
    } catch (std::exception e) {
        cerr << "Error: " << e.what() << endl;
        if (fl != NULL) {
            fclose(fl);
            fl = NULL;
        }
        throw;
    }
    if (fl != NULL) {
        fclose(fl);
        fl = NULL;
    }
    
    delete encoder;
}

int main(int argc, char** argv) {

    std::cout << "ParametersEncoder_test: ";
    
    test0();
    test1();
    test2();
    test3();
    test6();
    test7();
    test8();
    test9();
    
    std::cout << " " << std::endl;

    return (EXIT_SUCCESS);
}

#endif
