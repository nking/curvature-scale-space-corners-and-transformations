package algorithms;

import java.util.Arrays;

/**
 * @author nichole
 */
public class MergeSort {

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
    
}
