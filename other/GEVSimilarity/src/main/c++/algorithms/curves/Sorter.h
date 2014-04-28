/* 
 * File:   Sorter.h
 * Author: nichole
 * 
 * sort vector a in decreasing values and perform the same operations 
 * on vector b.  it uses quick sort to sort in place.
 * avg runtime is O(N*lg2(N)) w/ worse case runtime of O(N^2).
 *
 * Created on April 27, 2014
 */
#ifndef ALGORITHMS_CURVES_SORTER_H
#define ALGORITHMS_CURVES_SORTER_H

#include <vector>
#include "Defs.h"

using std::vector;

namespace gev {
class Sorter {

public:
    Sorter();
    
    virtual ~Sorter();
    
    /*
     sort vector a by decreasing value and perform same operations
     on vector b.
     @param a vector of values to be sorted in decreasing ordering
     @param b vector to receive same changes as performed on vector a.
     */
    void sort(vector<int>* a, vector<int>* b);
    
private:
    DISALLOW_COPY_AND_ASSIGN(Sorter);
        
    void _sort(vector<int>* a, vector<int>* b, int idxLo, int idxHi);
    
    int _partition(vector<int>* a, vector<int>* b, int idxLo, int idxHi);
    
};
} // end namespace
#endif