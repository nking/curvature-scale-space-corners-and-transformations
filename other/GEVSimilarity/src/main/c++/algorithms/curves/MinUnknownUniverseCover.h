/* 
 * File:   MinimumUnknownUniverseCover.h
 * Author: nichole
 *
 * Created on April 27, 2014
 */
#ifndef ALGORITHMS_CURVES_MINUKNOWNUNIVERSECOVER_H
#define ALGORITHMS_CURVES_MINUKNOWNUNIVERSECOVER_H

#include <vector>
#include "Defs.h"

using std::vector;

namespace gev {
class MinUnknownUniverseCover {

public:
    MinUnknownUniverseCover();
    
    virtual ~MinUnknownUniverseCover();
    
    void calculateCover(vector<vector<int> >* inputVariants,
        vector<int>* outputCoverVariants);
        
private:
    DISALLOW_COPY_AND_ASSIGN(MinUnknownUniverseCover);
        
};
} // end namespace
#endif