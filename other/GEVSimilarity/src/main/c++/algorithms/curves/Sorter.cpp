/* 
 * File:   Sorter.cpp
 * Author: nichole
 * 
 * Created on April 27, 2014
 */

#include "Sorter.h"

namespace gev {
    
    Sorter::Sorter() {
    }
    
    Sorter::~Sorter() {
    }
    
    // use quicksort to sort in place at avg runtime cost of O(N lg2 (N))
    void Sorter::sort(vector<uint32_t>& a, vector<uint32_t>& b) {
        if (a.empty()) {
            return;
        }
        if (a.size() != b.size()) {
            return;
        }
        _sort(a, b, 0, a.size() - 1);
    }
    
    void Sorter::_sort(vector<uint32_t>& a, vector<uint32_t>& b, long idxLo, 
        long idxHi) {
        
        if (idxLo < idxHi) {
            long idxMid = _partition(a, b, idxLo, idxHi);
            _sort(a, b, idxLo, idxMid - 1);
            _sort(a, b, idxMid + 1, idxHi);
        }
    }
    
    long Sorter::_partition(vector<uint32_t>& a, vector<uint32_t>& b, 
        long idxLo, long idxHi) {
                
        uint32_t x = a[idxHi];
        long store = idxLo - 1;
 
        for (long i = idxLo; i < idxHi; i++) {
            if (a[i] > x) {
                store++;
                if (i != store) {
                    long swap = a[i];
                    a[i] = a[store];
                    a[store] = swap;

                    swap = b[i];
                    b[i] = b[store];
                    b[store] = swap;
                }
            }
        }

        store++;
        if (store != idxHi) {
            long swap = a[store];
            a[store] = a[idxHi];
            a[idxHi] = swap;

            swap = b[store];
            b[store] = b[idxHi];
            b[idxHi] = swap;
        }
        
        return store;
    }
    
} // end namespace
