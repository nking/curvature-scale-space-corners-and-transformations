/* 
 * File:   MinimumUnknownUniverseCover.h
 * Author: nichole
 * 
 * <pre>
 * Given sets of similar curves which were assigned variable numbers for the 
 * unique parameters of {k, sigma, and mu}, the code finds the smallest set of 
 * the variable numbers which represent all curves.
 * 
 * For example, similar curves are listed in a row and the parameter sets 
 * are the numbers:
 * 1  2  3
 * 1
 *         4 5 6
 *           5
 * The smallest set of parameter sets representing all curves in the example
 * would be {1, 5}.
 * 
 * The "Universe" of the solution set is not known ahead of time, because
 * the goal is to only use 1 set from each row in the final "Universe" of sets.
 * 
 * runtime complexity is estimated as O(N lg_2(N)) + approx O(N^2 * m)
 * where m is the average number of variables in a row.
 * 
 * </pre>
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