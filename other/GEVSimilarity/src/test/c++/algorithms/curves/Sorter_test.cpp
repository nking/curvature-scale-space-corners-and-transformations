/*
 * Sorter_test.cpp
 *
 *  Created on: April 27, 2014
 *      Author: nichole
 */
#ifndef ALGORITHMS_CURVES_SORTER_TEST_H
#define ALGORITHMS_CURVES_SORTER_TEST_H

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
#include "main/c++/algorithms/curves/Sorter.h"

using std::vector;
using gev::Sorter;

void test0() {

    std::cout << "test0 " ;
    
    Sorter *sorter = new Sorter();
    
    vector<uint32_t> a;
    vector<uint32_t> b;
    
    a.push_back(3);
    a.push_back(10);
    a.push_back(20);
    a.push_back(1);
    
    b.push_back(2);
    b.push_back(1);
    b.push_back(0);
    b.push_back(3);
        
    sorter->sort(a, b);
    
    assert(a[0] == 20);
    assert(a[1] == 10);
    assert(a[2] == 3);
    assert(a[3] == 1);
    
    assert(b[0] == 0);
    assert(b[1] == 1);
    assert(b[2] == 2);
    assert(b[3] == 3);
    
    delete sorter;
}

void test1() {

    std::cout << "test1 " ;
    
    Sorter *sorter = new Sorter();
    
    vector<uint32_t> a;
    vector<uint32_t> aIndexes;
    vector<uint32_t> aOriginal;
    
    // random tests    
    time_t seed = time(NULL);
    srandom(seed);
        
    unsigned long nIter = 200000;
    for (unsigned long i = 0; i < nIter; i++) {
        long r = random();
        uint32_t number = (uint32_t)r;
        // we're only using positive numbers the using code, but testing for +-
        a.push_back(number);
        aOriginal.push_back(number);
        aIndexes.push_back((uint32_t)i);
    }

    sorter->sort(a, aIndexes);
    
    assert(a.size() == nIter);
    assert(aIndexes.size() == nIter);
    
    for (unsigned long i = 1; i < nIter; i++) {
        assert(a[i] <= a[i-1]);
        uint32_t index = aIndexes[i];
        assert(aOriginal[index] == a[i]);
    }
    
    delete sorter;
}

void test2() {

    std::cout << "test2 " ;
    
    Sorter *sorter = new Sorter();
    
    unsigned long len = 6;
    float* a = (float *)malloc(len * sizeof(float));
    float* b = (float *)malloc(len * sizeof(float));
    float* c = (float *)malloc(len * sizeof(float));
    
    /*
     4.000000e+00  9.499999e-01  3.500000e-01
     2.000000e+00  5.500000e-01  6.500000e-01
     1.500000e-01  4.000000e-01  3.500000e-01
     3.000000e-01  4.000000e-01  3.500000e-01
     4.000000e+00  9.499999e-01  5.500000e-01
     3.000000e-01  3.000000e-01  3.500000e-01
     */
    a[0] = 4.000000e+00;
    b[0] = 9.499999e-01;
    c[0] = 3.500000e-01;
     a[1] = 2.000000e+00;
     b[1] = 5.500000e-01;
     c[1] = 6.500000e-01;
    a[2] = 1.500000e-01;
    b[2] = 4.000000e-01;
    c[2] = 3.500000e-01;
     a[3] = 3.000000e-01;
     b[3] = 4.000000e-01;
     c[3] = 3.500000e-01;
    a[4] = 4.000000e+00;
    b[4] = 9.499999e-01;
    c[4] = 5.500000e-01;
     a[5] = 3.000000e-01;
     b[5] = 3.000000e-01;
     c[5] = 3.000000e-01;
    
    sorter->sort(a, b, c, len);
    
    /*
     1.500000e-01  4.000000e-01  3.500000e-01
     3.000000e-01  3.000000e-01  3.500000e-01
     3.000000e-01  4.000000e-01  3.500000e-01
     2.000000e+00  5.500000e-01  6.500000e-01   
     4.000000e+00  9.499999e-01  3.500000e-01
     4.000000e+00  9.499999e-01  5.500000e-01     
     */
    float* ea = (float *)malloc(len * sizeof(float));
    float* eb = (float *)malloc(len * sizeof(float));
    float* ec = (float *)malloc(len * sizeof(float));
    ea[0] = 1.500000e-01;
    eb[0] = 4.000000e-01;
    ec[0] = 3.500000e-01;
     ea[1] = 3.000000e-01;
     eb[1] = 3.000000e-01;
     ec[1] = 3.500000e-01;
    ea[2] = 3.000000e-01;
    eb[2] = 4.000000e-01;
    ec[2] = 3.500000e-01;
     ea[3] = 2.000000e+00;
     eb[3] = 5.500000e-01;
     ec[3] = 6.500000e-01;
    ea[4] = 4.000000e+00;
    eb[4] = 9.499999e-01;
    ec[4] = 3.500000e-01;
     ea[5] = 4.000000e+00;
     eb[5] = 9.499999e-01;
     ec[5] = 5.500000e-01;
    //printf("\n");
    for (unsigned long i = 0; i < len; i++) {
        float diff = a[i] - ea[i];
        if (diff < 0) {
            diff *= -1;
        }
        float eps = 0.01*ea[i];
        if (eps < 0) {
            eps *= -1;
        }
        //printf("%e  %e  %e\n", a[i], b[i], c[i]);
        assert(diff < eps);
    }
    
    free(a);
    free(b);
    free(c);
    free(ea);
    free(eb);
    free(ec);
    
    delete sorter;
}
    
int main(int argc, char** argv) {

    std::cout << "Sorter_test: ";
    
    test0();
    
    test1();
    
    test2();
    
    std::cout << " " << std::endl;

    return (EXIT_SUCCESS);
}

#endif
