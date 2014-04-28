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
    void Sorter::sort(vector<int>* a, vector<int>* b) {
        if (a->empty()) {
            return;
        }
        if (a->size() != b->size()) {
            return;
        }
        _sort(a, b, 0, a->size() - 1);
    }
    
    void Sorter::_sort(vector<int>* a, vector<int>* b, int idxLo, int idxHi) {
        if (idxLo < idxHi) {
            int idxMid = _partition(a, b, idxLo, idxHi);
            _sort(a, b, idxLo, idxMid - 1);
            _sort(a, b, idxMid + 1, idxHi);
        }
    }
    
    int Sorter::_partition(vector<int>* a, vector<int>* b, int idxLo, int idxHi) {
        int x = (*a)[idxHi];
        int store = idxLo - 1;
        
        for (int i = idxLo; i < idxHi; i++) {
            if ((*a)[i] > x) {
                store++;
                int swap = (*a)[i];
                (*a)[i] = (*a)[store];
                (*a)[store] = swap;
                
                swap = (*b)[i];
                (*b)[i] = (*b)[store];
                (*b)[store] = swap;
            }
        }
        int swap = (*a)[store + 1];
        (*a)[store + 1] = (*a)[idxHi];
        (*a)[idxHi] = swap;
        
        swap = (*b)[store + 1];
        (*b)[store + 1] = (*b)[idxHi];
        (*b)[idxHi] = swap;
        
        return store + 1;
    }
    
} // end namespace
