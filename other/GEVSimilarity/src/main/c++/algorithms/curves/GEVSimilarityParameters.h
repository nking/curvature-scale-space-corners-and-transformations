/* 
 * 
 * File:   GEVSimilarityParameters.h
 * Author: nichole
 *
 * Created on April 25, 2014
 */
#ifndef ALGORITHMS_CURVES_GEVSIMILARITYPARAMETERS_H
#define ALGORITHMS_CURVES_GEVSIMILARITYPARAMETERS_H

#include "Defs.h"
#include <tr1/unordered_map>

using std::tr1::unordered_map;
using std::make_pair;

namespace gev {
    
class GEVSimilarityParameters {

public:
    GEVSimilarityParameters();
    
    virtual ~GEVSimilarityParameters();
    
private:
    DISALLOW_COPY_AND_ASSIGN(GEVSimilarityParameters);
    
};
} // end namespace
#endif // end ifdef