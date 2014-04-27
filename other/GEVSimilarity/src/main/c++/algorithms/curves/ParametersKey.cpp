/* 
 * File:   ParametersKey.cpp
 * Author: nichole
 * 
 * Created on April 26, 2014
 */

#include "ParametersKey.h"

namespace gev {
    
    //ParametersKey::ParametersKey() {}
    
    ParametersKey::ParametersKey(float k, float sigma, float mu) {
        kKey = _floatToIntBits(k);
        sigmaKey = _floatToIntBits(sigma);
        muKey = _floatToIntBits(mu);
        _createKey();
    }
    
    ParametersKey::ParametersKey(const ParametersKey& other)
       : kKey(other.kKey), sigmaKey(other.sigmaKey), muKey(other.muKey) {
        _createKey();
    }
    
    ParametersKey::~ParametersKey() {
    }
        
    uint32_t ParametersKey::_floatToIntBits(float f) {
        // adapted from nvdia Cg source code which is now avail under a non-restrictive, free license
        // http://www.nvidia.com/object/IO_20020719_7269.html
        if (isnan(f)) {
            return 0x7fc00000;
        } else {
            union int_float_bits u;
            u.f = f;
            return u.i;
        }
    }
    
    float ParametersKey::_intBitsToFloat(uint32_t b) {
        union int_float_bits u;
        u.i = b;
        return u.f;
        //float f;
        //memcpy(&f, &b, sizeof(f));
        //return f;
    }    
    
    bool isnan(float s) {
        // By IEEE 754 rule, NaN is not equal to NaN
        return s != s;
    }

    void ParametersKey:: decodeToParams(float* k, float* sigma, float* mu) {
        *k = _intBitsToFloat(kKey);
        *sigma = _intBitsToFloat(sigmaKey);
        *mu = _intBitsToFloat(muKey);
    }
    
    string ParametersKey::_getKey() const {
        return key;
    }
    
    void ParametersKey::_createKey() {
                    
        size_t sz = sizeof(uint32_t);

        _writeIntToCharBytes(&kKey);

        _writeIntToCharBytes(&sigmaKey);

        _writeIntToCharBytes(&muKey);

    }
    
    void ParametersKey::_writeIntToCharBytes(const uint32_t *readInt) {
        
        // assuming that string does not take more than 8 bits from a char
        for (int i = 0; i < 4; i++) {
            int shift = i * 8;
            uint32_t a = (*readInt >> shift) & 0xff;
            char b = (char) a;
            key.push_back(b);
        }
        /*
         num = 953267991  [23, 183, 209, 56]
         num = 981668463  [111, 18, 131, 58]
         num = 1008981770 [10, 215, 35, 60]
         b = []
         for i in range(0, 4) :
             shift = i * 8
             a=(num >> shift) & 0xff
             b.append(a)
         
         print b
         */
    }
    
    uint32_t ParametersKey::_getKKey() const {
        return kKey;
    }
    
    uint32_t ParametersKey::_getSigmaKey() const {
        return sigmaKey;
    }
    
    uint32_t ParametersKey::_getMuKey() const {
        return muKey;
    }
    
    bool ParametersKey::operator==(ParametersKey & other)const {
        if (kKey != other.kKey) { return false; }
        if (sigmaKey != other.sigmaKey) { return false; }
        if (muKey != other.muKey) { return false; }
        return true;
    }
    bool ParametersKey::operator!=(ParametersKey & other)const {
        if (kKey != other.kKey) { return true; }
        if (sigmaKey != other.sigmaKey) { return true; }
        if (muKey != other.muKey) { return true; }
        return false;
    }
   
    bool operator==(ParametersKey const & pk0, ParametersKey const & pk1) {
        if (pk0.kKey != pk1.kKey) { return false; }
        if (pk0.sigmaKey != pk1.sigmaKey) { return false; }
        if (pk0.muKey != pk1.muKey) { return false; }
        return true;
    }
    
    bool operator!=(ParametersKey const & pk0, ParametersKey const & pk1) {
        if (pk0.kKey != pk1.kKey) { return true; }
        if (pk0.sigmaKey != pk1.sigmaKey) { return true; }
        if (pk0.muKey != pk1.muKey) { return true; }
        return false;
    }
    
} // end namespace
