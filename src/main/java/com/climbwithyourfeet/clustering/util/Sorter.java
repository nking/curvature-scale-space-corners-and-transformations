package com.climbwithyourfeet.clustering.util;

import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author nichole
 */
public class Sorter {
    
    public static void mergeSortByXThenY(List<PairFloat> a) {
        if (a == null) {
            throw new IllegalArgumentException("a cannot be null");
        }
        if (a.size() < 2) {
            return;
        }
        mergeSortByXThenY(a, 0, a.size() - 1);
    }
    
    public static void mergeSortByXThenY(List<PairFloat> a, int idxLo, int idxHi) {
        if (a == null) {
            throw new IllegalArgumentException("a cannot be null");
        }
        
        if (idxLo < idxHi) {
            
            int idxMid = (idxLo + idxHi)/2;
            mergeSortByXThenY(a, 0, idxMid);
            mergeSortByXThenY(a, idxMid + 1, idxHi); 
            mergeByXThenY(a, idxLo, idxMid, idxHi);
        }
    }
    
    public static void mergeSortByYThenX(List<PairFloat> a) {
        if (a == null) {
            throw new IllegalArgumentException("a cannot be null");
        }
        if (a.size() < 2) {
            return;
        }
        mergeSortByYThenX(a, 0, a.size() - 1);
    }
    
    public static void mergeSortByYThenX(List<PairFloat> a, int idxLo, int idxHi) {
        if (a == null) {
            throw new IllegalArgumentException("a cannot be null");
        }
        
        if (idxLo < idxHi) {
            
            int idxMid = (idxLo + idxHi)/2;
            mergeSortByYThenX(a, 0, idxMid);
            mergeSortByYThenX(a, idxMid + 1, idxHi); 
            mergeByYThenX(a, idxLo, idxMid, idxHi);
        }
    }
    
    private final static PairFloat sentinel = new PairFloat(Float.MAX_VALUE, Float.MAX_VALUE);

    private static void mergeByXThenY(List<PairFloat> a, int idxLo, int idxMid, 
        int idxHi) {
        
        List<PairFloat> left = new ArrayList<PairFloat>(a.subList(idxLo, idxMid + 1));
        List<PairFloat> right = new ArrayList<PairFloat>(a.subList(idxMid + 1, idxHi + 1));
        left.add(sentinel);
        right.add(sentinel);
        int leftPos = 0;
        int rightPos = 0;

        for (int k = idxLo; k <= idxHi; ++k) {
            boolean setL = false;
            if (left.get(leftPos).getX() < right.get(rightPos).getX()) {
                setL = true;
            } else if (left.get(leftPos).getX() == right.get(rightPos).getX()) {
                if (left.get(leftPos).getY() < right.get(rightPos).getY()) {
                    setL = true;
                }
            }
            if (setL) {
                a.set(k, left.get(leftPos));
                leftPos++;
            } else {
                a.set(k, right.get(rightPos));
                rightPos++;
            }
        }
    }
    
    private static void mergeByYThenX(List<PairFloat> a, int idxLo, int idxMid, 
        int idxHi) {
        
        List<PairFloat> left = new ArrayList<PairFloat>(a.subList(idxLo, idxMid + 1));
        List<PairFloat> right = new ArrayList<PairFloat>(a.subList(idxMid + 1, idxHi + 1));
        left.add(sentinel);
        right.add(sentinel);
        int leftPos = 0;
        int rightPos = 0;

        for (int k = idxLo; k <= idxHi; ++k) {
            boolean setL = false;
            if (left.get(leftPos).getY() < right.get(rightPos).getY()) {
                setL = true;
            } else if (left.get(leftPos).getY() == right.get(rightPos).getY()) {
                if (left.get(leftPos).getX() < right.get(rightPos).getX()) {
                    setL = true;
                }
            }
            if (setL) {
                a.set(k, left.get(leftPos));
                leftPos++;
            } else {
                a.set(k, right.get(rightPos));
                rightPos++;
            }
        }
    }
}
