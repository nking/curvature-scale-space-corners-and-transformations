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
    
    vector<int> a;
    vector<int> b;
    
    a.push_back(3);
    a.push_back(10);
    a.push_back(20);
    a.push_back(1);
    
    b.push_back(2);
    b.push_back(1);
    b.push_back(0);
    b.push_back(3);
    
    sorter->sort(&a, &b);
        
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
    
    vector<int> a;
    vector<int> aIndexes;
    vector<int> aOriginal;
    
    // random tests    
    time_t seed = time(NULL);
    srandom(seed);
    
    unsigned long nIter = 200000;
    for (unsigned long i = 0; i < nIter; i++) {
        long r = random();
        int number = (int)r;
        // we're only using positive numbers the using code, but testing for +-
        a.push_back(number);
        aOriginal.push_back(number);
        aIndexes.push_back(i);
    }
    
    sorter->sort(&a, &aIndexes);
    
    assert(a.size() == nIter);
    assert(aIndexes.size() == nIter);
    
    for (unsigned long i = 1; i < nIter; i++) {
        assert(a[i] <= a[i-1]);
        int index = aIndexes[i];
        assert(aOriginal[index] == a[i]);
    }
    
    delete sorter;
}

    
int main(int argc, char** argv) {

    std::cout << "Sorter_test: ";
    
    test0();
    
    test1();
    
    std::cout << " " << std::endl;

    return (EXIT_SUCCESS);
}

#endif
