package algorithms;

import algorithms.util.PairIntWithIndex;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * @author nichole
 */
public class MergeSort {
    
    /**
     * sort by increasing value a.x then a.y
     *
     * @param a array of points to be sorted
     */
    public static void sortByXThenY(List<PairIntWithIndex> a) {
        if (a == null) {
            throw new IllegalArgumentException("a cannot be null");
        }
        
        if (a.size() < 2) {
            return;
        }
        
        sortByXThenY(a, 0, a.size() - 1);
              
    }
    /**
     * sort by increasing value a.y then a.x
     *
     * @param a array of points to be sorted
     */
    public static void sortByYThenX(List<PairIntWithIndex> a) {
        if (a == null) {
            throw new IllegalArgumentException("a cannot be null");
        }
        
        if (a.size() < 2) {
            return;
        }
        
        sortByYThenX(a, 0, a.size() - 1);
              
    }

    /**
     * sort by decreasing value a1
     *
     * @param a array of points to be sorted
     */
    public static void sortByDecr(int[] a) {
        if (a == null) {
            throw new IllegalArgumentException("a cannot be null");
        }
        
        if (a.length < 2) {
            return;
        }
        
        sortByDecr(a, 0, a.length - 1);
              
    }
    
    public static void sortByDecr(int[] a, int idxLo, int idxHi) {

        if (a.length < 2) {
            return;
        }
        
        if (idxLo < idxHi) {

            int indexMid = (idxLo + idxHi) >> 1;
            
            sortByDecr(a, idxLo, indexMid);
            
            sortByDecr(a, indexMid + 1, idxHi);
            
            mergeByDecr(a, idxLo, indexMid, idxHi);
        }
    }
    
    public static void sortByXThenY(List<PairIntWithIndex> a, int idxLo, int idxHi) {

        if (a.size() < 2) {
            return;
        }
        
        if (idxLo < idxHi) {

            int indexMid = (idxLo + idxHi) >> 1;
            
            sortByXThenY(a, idxLo, indexMid);
            
            sortByXThenY(a, indexMid + 1, idxHi);
            
            mergeByXThenY(a, idxLo, indexMid, idxHi);
        }
    }
    
    public static void sortByYThenX(List<PairIntWithIndex> a, int idxLo, int idxHi) {

        if (a.size() < 2) {
            return;
        }
        
        if (idxLo < idxHi) {

            int indexMid = (idxLo + idxHi) >> 1;
            
            sortByYThenX(a, idxLo, indexMid);
            
            sortByYThenX(a, indexMid + 1, idxHi);
            
            mergeByYThenX(a, idxLo, indexMid, idxHi);
        }
    }
    
    private static void mergeByDecr(int[] a1, int idxLo, 
        int idxMid, int idxHi) {

        int[] a1Left = Arrays.copyOfRange(a1, idxLo, idxMid + 2);
        
        int[] a1Right = Arrays.copyOfRange(a1, idxMid + 1, idxHi + 2);
        
        a1Left[a1Left.length - 1] = Integer.MIN_VALUE;
        a1Right[a1Right.length - 1] = Integer.MIN_VALUE;
        
        int leftPos = 0;
        int rightPos = 0;

        for (int k = idxLo; k <= idxHi; k++) {
            int l = a1Left[leftPos];
            int r = a1Right[rightPos];
            if (l > r) {
                a1[k] = a1Left[leftPos];
                leftPos++;
            } else {
                a1[k] = a1Right[rightPos];
                rightPos++;
            }
        }
    }
    
    private static PairIntWithIndex sentinel = new PairIntWithIndex(
        Integer.MAX_VALUE, Integer.MAX_VALUE, 0);
    
    private static void mergeByXThenY(List<PairIntWithIndex> a1, int idxLo, 
        int idxMid, int idxHi) {

        List<PairIntWithIndex> a1Left = new ArrayList<PairIntWithIndex>
            (a1.subList(idxLo, idxMid + 1));
        
        List<PairIntWithIndex> a1Right = new ArrayList<PairIntWithIndex>
            (a1.subList(idxMid + 1, idxHi + 1));
                
        a1Left.add(sentinel);
        a1Right.add(sentinel);
        
        int leftPos = 0;
        int rightPos = 0;

        for (int k = idxLo; k <= idxHi; k++) {
            PairIntWithIndex l = a1Left.get(leftPos);
            PairIntWithIndex r = a1Right.get(rightPos);
            boolean lft = false;
            if (l.getX() < r.getX()) {
                lft = true;
            } else if (l.getX() == r.getX()) {
                if (l.getY() < r.getY()) {
                    lft = true;
                }
            }
            if (lft) {
                a1.set(k, a1Left.get(leftPos));
                leftPos++;
            } else {
                a1.set(k, a1Right.get(rightPos));
                rightPos++;
            }
        }
    }
    
    private static void mergeByYThenX(List<PairIntWithIndex> a1, int idxLo, 
        int idxMid, int idxHi) {

        List<PairIntWithIndex> a1Left = new ArrayList<PairIntWithIndex>
            (a1.subList(idxLo, idxMid + 1));
        
        List<PairIntWithIndex> a1Right = new ArrayList<PairIntWithIndex>
            (a1.subList(idxMid + 1, idxHi + 1));
                
        a1Left.add(sentinel);
        a1Right.add(sentinel);
        
        int leftPos = 0;
        int rightPos = 0;

        for (int k = idxLo; k <= idxHi; k++) {
            PairIntWithIndex l = a1Left.get(leftPos);
            PairIntWithIndex r = a1Right.get(rightPos);
            boolean lft = false;
            if (l.getY() < r.getY()) {
                lft = true;
            } else if (l.getY() == r.getY()) {
                if (l.getX() < r.getX()) {
                    lft = true;
                }
            }
            if (lft) {
                a1.set(k, a1Left.get(leftPos));
                leftPos++;
            } else {
                a1.set(k, a1Right.get(rightPos));
                rightPos++;
            }
        }
    }
}
