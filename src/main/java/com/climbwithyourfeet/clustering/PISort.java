package com.climbwithyourfeet.clustering;

import com.climbwithyourfeet.clustering.util.PairInt;
import java.util.Arrays;

/**
 *
 * @author nichole
 */
public class PISort {
    
    /**
     *
     * @param a
     */
    public static void mergeSortByXThenY(PairInt[] a) {
        mergeSortByXThenY(a, 0, a.length - 1);
    }
    
    /**
     *
     * @param a
     * @param idxLo
     * @param idxHi
     */
    public static void mergeSortByXThenY(PairInt[] a, int idxLo, int idxHi) {
        if (a == null) {
            throw new IllegalArgumentException("a cannot be null");
        }
        
        if (idxLo < idxHi) {

            int idxMid = (idxLo + idxHi) >> 1;

            mergeSortByXThenY(a, idxLo, idxMid);
            mergeSortByXThenY(a, idxMid + 1, idxHi);
            mergeByXThenY(a, idxLo, idxMid, idxHi);
        }
    }
    
    /**
     * for use in default, ascending sort of PairInt objects
     */
    private static PairInt maxSentinel = new PairInt(Integer.MAX_VALUE, Integer.MAX_VALUE);

    private static void mergeByXThenY(PairInt[] a, int idxLo, int idxMid, int idxHi) {
        
        PairInt[] aLeft = Arrays.copyOfRange(a, idxLo, idxMid + 2);
        PairInt[] aRight = Arrays.copyOfRange(a, idxMid + 1, idxHi + 2);
        
        aLeft[aLeft.length - 1] = maxSentinel;
        aRight[aRight.length - 1] = maxSentinel;

        int leftPos = 0;
        int rightPos = 0;

        for (int k = idxLo; k <= idxHi; k++) {
            
            int l = aLeft[leftPos].getX();
            int r = aRight[rightPos].getX();
            
            if (l == r) {
                float ly = aLeft[leftPos].getY();
                float ry = aRight[rightPos].getY();

                if (ly <= ry) {
                    a[k] = aLeft[leftPos];
                    leftPos++;
                } else {
                    a[k] = aRight[rightPos];
                    rightPos++;
                }
            } else if (l < r) {
                a[k] = aLeft[leftPos];
                leftPos++;
            } else {
                a[k] = aRight[rightPos];
                rightPos++;
            }
        }
    }
    
    /**
     *
     * @param a
     */
    public static void quickSortByXThenY(PairInt[] a) {
        quickSortByXThenY(a, 0, a.length - 1);
    }
    
    /**
     *
     * @param a
     * @param idxLo
     * @param idxHi
     */
    public static void quickSortByXThenY(PairInt[] a, int idxLo, int idxHi) {
        if (a == null) {
            throw new IllegalArgumentException("a cannot be null");
        }
        if (idxLo < idxHi) {
            int idxMid = partitionByXThenY(a, idxLo, idxHi);
            quickSortByXThenY(a, idxLo, idxMid - 1);
            quickSortByXThenY(a, idxMid + 1, idxHi);
        }
    }

    private static int partitionByXThenY(PairInt[] a, int idxLo, int idxHi) {
        
        PairInt comp = a[idxHi];
        int store = idxLo - 1;
        
        for (int i = idxLo; i < idxHi; ++i) {
            boolean doSwap = false;
            if (a[i].getX() < comp.getX()) {
                doSwap = true;
            } else if (a[i].getX() == comp.getX()) {
                if (a[i].getY() <= comp.getY()) {
                    doSwap = true;
                }
            }
            if (doSwap) {
                store++;
                PairInt swap = a[store];
                a[store] = a[i];
                a[i] = swap;
            }
        }
        
        store++;
        
        PairInt swap = a[store];
        a[store] = a[idxHi];
        a[idxHi] = swap;
        
        return store;
    }
}
