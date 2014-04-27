/* 
 * File:   ParametersKey.h
 * Author: nichole
 * 
 * class to hold the identity of a set of parameters k, sigma, and mu
 *
 * Created on April 26, 2014
 */
#ifndef ALGORITHMS_CURVES_PARAMETERSKEY_H
#define ALGORITHMS_CURVES_PARAMETERSKEY_H

#include <stdint.h>
#include <tr1/unordered_map>
#include <string>

using std::string;
using std::tr1::unordered_map;

namespace gev {
    
class ParametersKey {

public:
        
    ParametersKey(const ParametersKey&);
    
    ParametersKey(float k, float sigma, float mu);
    
    virtual ~ParametersKey();
    
    // For quickest implementation of determining the unique set of ParametersKey
    // instances, code using this class will use the boost unordered_map. 
    // The unordered_map presumably implements perfect hashing 
    // for the keys it supports: string, wstring, float, double, long double
    // when there are collisions.  the internal hash_value functions return
    // a type of size_t, so presumably the hash function can only provide
    // perfect hashing for input <= 32 bit, so the unordered_map must be using
    // methods to create internal data structures to keep track of entries
    // with different key_equals that collide.
    //
    // In order to use ParametersKey as a key, I need to either convert the
    //  identity to a string, or provide a hashing function.
    // Creating a string key is the easiest for now and allows quick testing
    // that all values even when they collide lead to successful storage and
    // retrieval from unordered_map for unique keys.
    
    // The identity for this class is k, sigma, and mu which are 32-bit floats.
    // those floats are converted to their 32-bit integer equivalents (IEEE 754).
    // Then a string is created from them as code points.
    // The char data type is usually 8 bits and string is composed of chars.
    // So the identity of this class is 3 32-bit fields written as a 12 element 
    // char array.
    
    uint32_t _floatToIntBits(float f);
    float _intBitsToFloat(uint32_t b);
    
    bool isnan(float s);

    void decodeToParams(float* k, float* sigma, float* mu);
    
    string _getKey() const;
    uint32_t _getKKey() const;
    uint32_t _getSigmaKey() const;
    uint32_t _getMuKey() const;

    bool operator==(ParametersKey & other)const;
    bool operator!=(ParametersKey & other)const;
    
    friend bool operator==(ParametersKey const & pk0, ParametersKey const & pk1);
    
    friend bool operator!=(ParametersKey const & pk0, ParametersKey const & pk1);
    
protected:
    
    union int_float_bits {
        uint32_t i;
        float f;
    };
    
    // to reduce space, could remove these and re-calculate from key upon need
    uint32_t kKey;
    uint32_t sigmaKey;
    uint32_t muKey;
    
    void _createKey();
    
private:
    
    string key;
    
    void _writeIntToCharBytes(const uint32_t *readInt);
};

} // end namespace

namespace std {
    namespace tr1 {
        template<> struct hash<gev::ParametersKey>
        : public std::unary_function<gev::ParametersKey, std::size_t> {
            std::size_t operator()(const gev::ParametersKey &pk) const {
                string tmp = pk._getKey();
                return hash<string>()(tmp);
            }
        };
    }
    
    template<> struct equal_to<gev::ParametersKey> {
        bool operator()(const gev::ParametersKey& pk0, const gev::ParametersKey& pk1) const {
            if (pk0._getKKey() != pk1._getKKey()) { return false; }
            if (pk0._getSigmaKey() != pk1._getSigmaKey()) { return false; }
            if (pk0._getMuKey() != pk1._getMuKey()) { return false; }
            return true;
        }
    };
    
}
#endif