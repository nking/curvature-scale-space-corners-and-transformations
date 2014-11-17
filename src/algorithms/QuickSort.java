package algorithms;

/**
 *
 * @author nichole
 */
public class QuickSort {
    
     /**
     * sort a from index idxLo to idxHi, inclusive
     * @param a
     */
    public static void sort(float[] a) {
        sort(a, 0, a.length - 1);
    }
    
    /**
     * sort a from index idxLo to idxHi, inclusive
     * @param a
     * @param idxLo
     * @param idxHi 
     */
    public static void sort(float[] a, int idxLo, int idxHi) {
        
        if (idxLo < idxHi) {
            int idxMid = partition(a, idxLo, idxHi);
            sort(a, idxLo, idxMid - 1);
            sort(a, idxMid + 1, idxHi);
        }
    }
    
    private static int partition(float[] a, int idxLo, int idxHi) {
        float x = a[idxHi];      // for comparison
        int store = idxLo - 1;   // store out of way to swap after pivot
        for (int i = idxLo; i < idxHi; i++) {
            if (a[i] <= x) {
                store++;
                float swap = a[store];
                a[store] = a[i];
                a[i] = swap;
            }
        }
        store++;
        float swap = a[store];
        a[store] = a[idxHi];
        a[idxHi] = swap;
        return store;
    }

}
