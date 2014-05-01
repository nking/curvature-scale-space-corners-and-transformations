/* 
 * 
 * File:   GEVSimilarityParameters.h
 * Author: nichole
 * 
 * class to invoke all of the work in reading the default input
 * file, finding the smallest set to create all curves, and
 * writing the result to the default output file.
 * 
 * The default input file:
 *    tmpdata2/algorithms.curves.GEVSimilarityToolTest-output.txt
 * 
 * The default output file:
 *    tmpdata2/gevSimilarityTool-output.txt
 * 
 * 
 * More on the input file:
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
 * Created on April 25, 2014
 */
#ifndef ALGORITHMS_CURVES_GEVSIMILARITYPARAMETERS_H
#define ALGORITHMS_CURVES_GEVSIMILARITYPARAMETERS_H

#include "Defs.h"
// for uint32_t
#include <stdint.h>
#include <tr1/unordered_map>
#include <vector>
#include "main/c++/algorithms/curves/MinUnknownUniverseCover.h"
#include "main/c++/algorithms/curves/ParametersEncoder.h"

using std::vector;
using std::tr1::unordered_map;
using std::make_pair;

namespace gev {
    
class GEVSimilarityParameters {

public:
    GEVSimilarityParameters();
    
    virtual ~GEVSimilarityParameters();
    
    /*
     * calculate the minimum set of parameters that will generate all curves in
     * the default input file and write the results to the default output
     * file.
     * 
     * See the class level comments for a description of the input file.
     * 
     * The default input file:
     *    tmpdata2/algorithms.curves.GEVSimilarityToolTest-output.txt
     * 
     * The default output file:
     *    tmpdata2/gevSimilarityTool-output.txt
     * 
     */
    void calculateMinSet();
    
    /*
     * calculate the minimum set of parameters that will generate all curves in
     * the file inputFile.
     * 
     * @param inputFileName a file of format described in the class level 
     * comments.
     * @param outputFile the file the results are written to in directory 
     * tmpdata2 relative to the project base directory.
     */
    void calculateMinSet(string inputFileName, string outputFileName);
    
private:
    DISALLOW_COPY_AND_ASSIGN(GEVSimilarityParameters);
    
};
} // end namespace
#endif