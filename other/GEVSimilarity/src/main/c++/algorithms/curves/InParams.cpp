/* 
 * File:   InParams.cpp
 * Author: nichole
 * 
 * Created on April 25, 2014
 */

#include "InParams.h"

namespace gev {
    
    InParams::InParams() : nRows(), nPerRow(NULL), simRows(NULL),
        k(NULL), sigma(NULL), mu(NULL) {

    }

    InParams::~InParams() {
        
        if (nPerRow != NULL) {
            free(nPerRow);
        }
        
        if (k != NULL) {
            free(k);
        }
        
        if (sigma != NULL) {
            free(sigma);
        }
        
        if (mu != NULL) {
            free(mu);
        }
        
        if (simRows != NULL) {
            // depends upon choice of jagged or contiguous memory arrays model
        }
    }
    
    InParams::_expandIfNeeded() {
        
    }
    
    InParams::addParams(float* kParams, float* sigmaParams, float* muParams) {
        
    }
    
} // end namespace
