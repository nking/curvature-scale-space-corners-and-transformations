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
//for uint32_t
#include <stdint.h>
// for EINVAL
#include <errno.h>
// for NULL
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>
#include <tr1/unordered_set>
#include <tr1/unordered_map>
#include "ParametersKey.h"

using std::string;
using std::vector;
using std::tr1::unordered_map;
using std::tr1::unordered_set;
using std::make_pair;
using std::ifstream;
using std::cout;
using std::cerr;
using std::endl;
using std::flush;

namespace gev {
    
class ParametersEncoder {

public:
    ParametersEncoder();
    
    virtual ~ParametersEncoder();
        
    /* <pre>
     * Read the parsed log file and create the internal encoding map to 
     * use to fill encoded variants arrays.
     * 
     * The input file is expected to be the logfile from the java code
     * GEVSimilarityToolTest which is placed in target/surefire-reports
     * and is called algorithms.curves.GEVSimilarityToolTest-output.txt.
     * That file should be parsed using cat and grep to keep only lines
     * containing "k=".  
     * cat algorithms.curves.GEVSimilarityToolTest-output.txt \
     *     | grep "k=" > | grep -v Uniq > out.txt
     * 
     * Each line in the file then contains without a line break:
     *    k=comma separated values <2 spaces> sigma=comma sep vals <2 spaces> \
     *      x-m=comma sep vals <2 spaces> mu=comma sep vals
     * 
     * Each line holds the k, sigma, and mu parameters that produced similar
     * curves.  Those parameters are ordered with respect to one another.
     * k=k0,k1,...kn  sigma=sigma0,sigma1,...sigman  mu=mu0,mu1,...mun
     * where the order is preserved such that parameter sets that were used
     * to generate curves were
     *    set0 is {k0, sigma0, mu0}
     *    set1 is {k1, sigma1, mu1}
     *     ...
     * 
     * Example:
     *   k=0.010000001,0.020000001,0.020000001,0.020000001,0.020000001\
     *   sigma=0.099999994,0.099999994,0.099999994,0.099999994,0.099999994\
     *   mu=-0.010000001,-0.008571431,-0.009090915,-0.009090915,-0.008571431
     * 
     * NOTE:  the code expects to find and read the file:
     *     tmpdata2/algorithms.curves.GEVSimilarityToolTest-output.txt
     * </pre>
     * 
     * @param outputEncodedVariables a vector holding each line of the file
     * as a row that represents parameters which produced curves that are
     * similar to one another.  The encoding assigns a variable number to
     * each unique combination of parameters {k, sigma, and mu} unique
     * over the entire set.
     * Note that for this file, each parameter set is already unique,
     * but the method handles the case in which they weren't too.
     */
    void readFile(vector< unordered_set<uint32_t> >* outputEncodedVariables);
    
    /*
     read log file described in readFile(vector<vector<int> >* encodedVariants)
     @param fileName file that will be found in directory tmpdata2
     @param encodedVariants 
     */
    void _readFile(string fileName, 
        vector< unordered_set<uint32_t> >* outputEncodedVariables);
    
    void writeFile(vector<uint32_t>* encodedCoverVariables);
    
    void _writeFile(string fileName, vector<uint32_t>* encodedCoverVariables);

    string _getCWD();
    
    string _getProjectTmpDirectoryPath();
    
    uint32_t _getNVarsOfAParameter(const char *line, const uint32_t digit0, 
        const uint32_t digitn);
    
    void _parseLine(const char *line, const uint32_t digit0, const uint32_t digitn, 
        float *k, float *sigma, float *mu, const uint32_t nVars);
    
    /*
     * parse the portion of the line between digit0 and digitn to fill array a
     * with the number nVars of parsed floats.
     * @param digit0 the index of the first character in line to be included in
     * the parsing.
     * @param digitn the index of the last character in line to be included in
     * the parsing.
     * @param a float array that will be filled with the parsed numbers
     * @param nVars the number of variables to be parsed
     */
    uint32_t _parseLineForNextNVars(const char *line, const uint32_t digit0, 
        const uint32_t digitn, float *a, const uint32_t nVars);
    
    /*
     * parse the float found in line from characters digitBegin to
     * digitEnd inclusive.
     * @param digitBegin the index of the first character in line to be included 
     * in the float.
     * @param digitEnd the index of the last character in line to be included in
     * the float.
     */
    float _convertToFloat(const char *line, const uint32_t digitBegin, 
        const uint32_t digitEnd);
    
private:
    DISALLOW_COPY_AND_ASSIGN(ParametersEncoder);
        
    /* map holding the read file content w/ values = unique variable numbers
     */
    unordered_map<ParametersKey, uint32_t> variableMap;
    
};
} // end namespace
#endif