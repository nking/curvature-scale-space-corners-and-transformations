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
        long swap = a[store];
        a[store] = a[idxHi];
        a[idxHi] = swap;

        swap = b[store];
        b[store] = b[idxHi];
        b[idxHi] = swap;
        
        return store;
    }
    
    //=========
    void Sorter::sort(float* a, float* b, float* c, long length) {
        if (a == NULL) {
            return;
        }
        if (length == 0) {
            return;
        }
        _sort(a, b, c, 0, length - 1);
    }
    
    void Sorter::_sort(float* a, float* b, float* c, long idxLo, long idxHi) {
        
        int q = -1;
        
        if (idxLo < idxHi) {
            q = (idxLo + idxHi)/2;
            _sort(a, b, c, idxLo, q);
            _sort(a, b, c, q + 1, idxHi);
            _merge(a, b, c, idxLo, q, idxHi);
        }
    }
    
    void Sorter::_merge(float* a, float* b, float* c, long idxLo, long idxMid,
        long idxHi) {
                
        long nLeft = idxMid + 2 - idxLo;
        long nRight =  idxHi + 2 - (idxMid + 1);
        size_t leftSz = nLeft * sizeof(float);
        size_t rightSz = nRight * sizeof(float);
        size_t leftSz0 = (nLeft - 1 ) * sizeof(float);
        size_t rightSz0 = (nRight - 1) * sizeof(float);
        
        float* leftA = (float *)malloc(leftSz);
        float* leftB = (float *)malloc(leftSz);
        float* leftC = (float *)malloc(leftSz);
        
        memcpy(leftA, a + idxLo, leftSz0);
        memcpy(leftB, b + idxLo, leftSz0);        
        memcpy(leftC, c + idxLo, leftSz0);
        
        float* rightA = (float *)malloc(rightSz);        
        float* rightB = (float *)malloc(rightSz);
        float* rightC = (float *)malloc(rightSz);
        
        memcpy(rightA, a + (idxMid + 1), rightSz0);
        memcpy(rightB, b + (idxMid + 1), rightSz0);        
        memcpy(rightC, c + (idxMid + 1), rightSz0);
        
        float sentinel = FLT_MAX;
        leftA[nLeft - 1] = sentinel;
        leftB[nLeft - 1] = sentinel;
        leftC[nLeft - 1] = sentinel;
        rightA[nRight - 1] = sentinel;
        rightB[nRight - 1] = sentinel;
        rightC[nRight - 1] = sentinel;
        
        long leftPos = 0;
        long rightPos = 0;
        
        //sort by a, then if equal b, then if equal c
        
        for (int k = idxLo; k <= idxHi; k++) {
            
            bool assignToLeft = (leftA[leftPos] < rightA[rightPos]);
            
            if (!assignToLeft && (leftA[leftPos] == rightA[rightPos])) {
                
                assignToLeft = (leftB[leftPos] < rightB[rightPos]);

                if (!assignToLeft && (leftB[leftPos] == rightB[rightPos])) {
                    
                    assignToLeft = (leftC[leftPos] <= rightC[rightPos]);
                }                
            }
            
            if (assignToLeft) {
                a[k] = leftA[leftPos];
                b[k] = leftB[leftPos];
                c[k] = leftC[leftPos];
                leftPos++;
            } else {
                a[k] = rightA[rightPos];
                b[k] = rightB[rightPos];
                c[k] = rightC[rightPos];
                rightPos++;
            }   
        }
        
        free(leftA);
        free(leftB);
        free(leftC);
        free(rightA);
        free(rightB);
        free(rightC);
    }
    
} // end namespace
