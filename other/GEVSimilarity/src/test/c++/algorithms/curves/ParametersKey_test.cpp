/*
 * ParametersKey_test.cpp
 *
 *  Created on: April 26, 2014
 *      Author: nichole
 */
#ifndef ALGORITHMS_CURVES_PARAMETERSKEY_TEST_H
#define ALGORITHMS_CURVES_PARAMETERSKEY_TEST_H

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
// for unordered_map
#include <tr1/unordered_map>
#include <vector>
#include "main/c++/algorithms/curves/ParametersKey.h"

using namespace std;

using std::tr1::unordered_map;

using gev::ParametersKey;
    
void test0() {

    std::cout << "test0 " ;
    
    // test the float to int bits methods
     
    ParametersKey *key = new ParametersKey(0.0001f, 0.001f, 0.01f);
   
    float a = 0.0001f;
    int32_t ai = key->_floatToIntBits(a);
    float aia = key->_intBitsToFloat(ai);
    assert(ai == 953267991);
    assert(a == aia);
    
    a = 0.001f;
    ai = key->_floatToIntBits(a);
    aia = key->_intBitsToFloat(ai);
    assert(ai == 981668463);
    assert(a == aia);
    
    a = 0.01f;
    ai = key->_floatToIntBits(a);
    aia = key->_intBitsToFloat(ai);
    assert(ai == 1008981770);
    assert(a == aia);
    
    // random tests    
    time_t seed = time(NULL);
    srandom(seed);
        
    int nIter = 1000;
    for (int i = 0; i < nIter; i++) {
        long r = random();
        a = key->_intBitsToFloat((uint32_t)r);
        if (a == a) { // if not NAN
            ai = key->_floatToIntBits(a);
            aia = key->_intBitsToFloat(ai);
            assert(a == aia);
        }
    }
    
    delete key;
}

void test1() {

    std::cout << "test1 " ;
    
    ParametersKey *key = new ParametersKey(0.0001f, 0.001f, 0.01f);
    
    /*
    0.0001f = 953267991  [23, 183, 209, 56]
    0.001f  = 981668463  [111, 18, 131, 58]
    0.01f   = 1008981770 [10, 215, 35, 60]
    */
    
    string keyStr = key->_getKey();
    const char *keyCStr = keyStr.c_str();
    
    assert( (int)keyCStr[0] == 23);
    
    delete key;
}

void test2() {

    std::cout << "test2 " ;
    
    // test the general use with an unordered_map
    
    float a = 1E-3;
    float b = 2E-2;
    float c = 3E-1;
    
    ParametersKey *key = new ParametersKey(a, b, c);
   
    float d = 1E-3;
    float e = 2E-2;
    float f = 3E-1;
    ParametersKey *key2 = new ParametersKey(d, e, f);
    
    assert(*key == *key2);
    
    // spot check that the hash_value and object equals are working as expected:
    unordered_map<ParametersKey, int> variantsMap;
    
    // don't put types on make_pair<gev::ParametersKey, int>  template or compiler won't find default defin
    variantsMap.insert(
        make_pair<gev::ParametersKey, int> (*key, 1));
        
    assert(variantsMap.size() == 1);
    
    assert(variantsMap.find(*key) != variantsMap.end());
    
    assert(variantsMap.find(*key2) != variantsMap.end());
    
    variantsMap.insert(
        make_pair<gev::ParametersKey, int> (*key, 2));
        
    assert(variantsMap.size() == 1);
    
    assert(variantsMap.find(*key) != variantsMap.end());    
    
    assert(variantsMap.find(*key2) != variantsMap.end());
    
    delete key;
    delete key2;
}

void test3() {
    
    std::cout << "test3 " ;
    
    // use random variables and save them to assert no collisions in map hash
    
    ParametersKey *key = new ParametersKey(0.0001f, 0.001f, 0.01f);
    
    unordered_map<ParametersKey, int> variantsMap;
    vector<float> ks;
    vector<float> sigmas;
    vector<float> mus;
    
    time_t seed = time(NULL);
    srandom(seed);
    
    float k;
    float sigma;
    float mu;
    
    int nIter = 1000;
    unsigned long nExpected = 0;
    for (int i = 0; i < nIter; i++) {
        long r = random();
        k = key->_intBitsToFloat((uint32_t)r);
        while (k != k) {
            r = random();
            k = key->_intBitsToFloat((uint32_t)r);
        }
        
        r = random();
        sigma = key->_intBitsToFloat((uint32_t)r);
        while (sigma != sigma) {
            r = random();
            sigma = key->_intBitsToFloat((uint32_t)r);
        }
        
        r = random();
        mu = key->_intBitsToFloat((uint32_t)r);
        while (mu != mu) {
            r = random();
            mu = key->_intBitsToFloat((uint32_t)r);
        }
                
        // randomly drawing 3 numbers in a row has rare chance of making same set
        //  of 3, so check before incrementing nExpected
        bool found = false;
        for (unsigned long j = 0; j < ks.size(); j++) {
            if (k == ks[j]) {
                if (sigma == sigmas[j]) {
                    if (mu == mus[j]) {
                        found = true;
                        break;
                    }
                }
            }
        }
        
        if (!found) {
            
            nExpected++;
            
            ks.push_back(k);
            sigmas.push_back(sigma);
            mus.push_back(mu);

            ParametersKey *key = new ParametersKey(k, sigma, mu);
            
            variantsMap.insert(
                make_pair<gev::ParametersKey, int> (*key, 2));
            
            delete key;
        }
    }
    assert(variantsMap.size() == nExpected);
    
    bool* allAreFound = (bool*)calloc(nExpected, sizeof(bool));
    
    const unsigned long notFound = -1;
    
    for (unordered_map<ParametersKey, int >::iterator iter = variantsMap.begin();
        iter != variantsMap.end(); ++iter ) {

        ParametersKey key = iter->first;
        
        key.decodeToParams(&k, &sigma, &mu);
                
        unsigned long idx = notFound;
        for (unsigned long j = 0; j < ks.size(); j++) {
            if (k == ks[j]) {
                if (sigma == sigmas[j]) {
                    if (mu == mus[j]) {
                        idx = j;
                        break;
                    }
                }
            }
        }
        if (idx != notFound) {
            allAreFound[idx] = true;
        }
    }
    
    for (unsigned long i = 0; i < nExpected; i++) {
        assert(allAreFound[i]);
    }
    
    free(allAreFound);
}

int main(int argc, char** argv) {

    test0();
    
    test1();
    
    test2();
    
    test3();

    std::cout << " " << std::endl;

    return (EXIT_SUCCESS);
}

#endif
