/* 
 * 
 * File:   InParams.h
 * Author: nichole
 *
 * Created on April 25, 2014
 */
#ifndef ALGORITHMS_CURVES_INPARAMS_H
#define ALGORITHMS_CURVES_INPARAMS_H

#include "Defs.h"
// for NULL
#include <cstdlib>

namespace gev {
    
class InParams {

public:
    InParams();
    
    virtual ~InParams();
    
    _expandIfNeeded();
    
    addParams(float* kParams, float* sigmaParams, float* muParams);
    
protected:
    
    /*
     the number of rows in simRows
     */
    int nRows;
    /*
     pointer to an array holding the number of items in a row in simRows
     */
    int *nPerRow;
    
    /*
     pointer to pointers comprising a two-dimensional array of index numbers.
     there are nRows of pointers to the arrays per row. each array holds indexes
     to an item in the array pointed to by k.
     */
    int** simRows;
    
    /*
     pointer to an array holding values of k.  the indexes in k are parallel to
     sigma and mu, so a complete set of parameters of {k, sigma, mu} has the same
     index for all 3 arrays.
     */
    float* k;
    
    /*
     pointer to an array holding values of sigma.  the indexes in k are parallel to
     sigma and mu, so a complete set of parameters of {k, sigma, mu} has the same
     index for all 3 arrays.
     */
    float* sigma;
    
    /*
     pointer to an array holding values of mu.  the indexes in k are parallel to
     sigma and mu, so a complete set of parameters of {k, sigma, mu} has the same
     index for all 3 arrays.
     */
    float* mu;
    
private:
    DISALLOW_COPY_AND_ASSIGN(InParams);
    
};
} // end namespace
#endif // end ifndef
