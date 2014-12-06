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
        if (a == null) {
            throw new IllegalArgumentException("a cannot be null");
        }
        if (a.length < 2) {
            return;
        }
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

    /**
     * sort a from index idxLo to idxHi, inclusive
     * @param a
     * @param b an array that will receive the same swap operations as are 
     * performed on a
     * @param c an array that will receive the same swap operations as are 
     * performed on a
     * @param idxLo
     * @param idxHi 
     */
    public static void sort(float[] a, int[] b, int[] c, int idxLo, int idxHi) {
        
        if (a == null) {
            throw new IllegalArgumentException("a cannot be null");
        }
        if (b == null) {
            throw new IllegalArgumentException("b cannot be null");
        }
        if (c == null) {
            throw new IllegalArgumentException("c cannot be null");
        }
        if ((a.length != b.length) || (a.length != c.length)) {
            throw new IllegalArgumentException("array lengths must be the same");
        }
        
        if (a.length < 2) {
            return;
        }
        
        if (idxLo < idxHi) {
            int idxMid = partition(a, b, c, idxLo, idxHi);
            sort(a, b, c, idxLo, idxMid - 1);
            sort(a, b, c, idxMid + 1, idxHi);
        }
    }
    
    private static int partition(float[] a, int[] b, int[] c, 
        int idxLo, int idxHi) {
        
        float x = a[idxHi];      // for comparison
        int store = idxLo - 1;   // store out of way to swap after pivot
        for (int i = idxLo; i < idxHi; i++) {
            if (a[i] <= x) {
                store++;
                float swap = a[store];
                a[store] = a[i];
                a[i] = swap;
                int swap2 = b[store];
                b[store] = b[i];
                b[i] = swap2;
                swap2 = c[store];
                c[store] = c[i];
                c[i] = swap2;
            }
        }
        store++;
        float swap = a[store];
        a[store] = a[idxHi];
        a[idxHi] = swap;
        int swap2 = b[store];
        b[store] = b[idxHi];
        b[idxHi] = swap2;
        swap2 = c[store];
        c[store] = c[idxHi];
        c[idxHi] = swap2;
        return store;
    }
}
