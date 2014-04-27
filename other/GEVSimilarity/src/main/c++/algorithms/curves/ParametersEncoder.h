/* 
 * 
 * File:   ParametersEncoder.h
 * Author: nichole
 *
 * Created on April 27, 2014
 */
#ifndef ALGORITHMS_CURVES_PARAMETERSENCODER_H
#define ALGORITHMS_CURVES_PARAMETERSENCODER_H

#include "Defs.h"
#include <vector>
#include <tr1/unordered_map>
#include "ParametersKey.h"

using namespace std;
using std::tr1::unordered_map;

namespace gev {
    
class ParametersEncoder {

public:
    ParametersEncoder();
    
    virtual ~ParametersEncoder();
        
    /*
     read the parsed log file and create the internal encoding map to 
     * use to fill encoded variants arrays
     */
    void readFile(vector<vector<int> >* encodedVariants);
    
    void _readFile(string fileName, vector<vector<int> >* encodedVariants);
    
    void writeFile(vector<int>* encodedCoverVariants);
    
    void _writeFile(string fileName, vector<int>* encodedCoverVariants);

private:
    DISALLOW_COPY_AND_ASSIGN(ParametersEncoder);
        
    /* map holding the read file content w/ values = unique variant numbers
     */
    unordered_map<ParametersKey, int> variantsMap;
    
};
} // end namespace
#endif