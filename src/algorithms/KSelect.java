package algorithms;

/**
 * identify the k-th smallest element in an unsorted array.
   Note that the arrays are sorted in place so make a copy of the array before
   using it as an argument if the array should not be modified.
   
  worse case runtime is less than O(N^2)
   avg case runtime is a little larger than O(N)
   
 * @author nichole
 */
public class KSelect {
    
    public float findMedianOfMedians(float[] a, int indexLo, int indexHi) {

        int medianIdx = findMedianOfMediansIdx(a, indexLo, indexHi);

        return a[medianIdx];
    }
    
     /**
     <pre>
     Complexity of run time:

            j = math.floor( math.log(N)/math.log(5))

            Total cost =  T(    ∑   (N/(5*i)) * const   )  + 2T(j/2) + j
                           ( ¡=1 to j                   )

            ===&gt; It's a little more than O(j) but less than O(j*lg(j)).
      </pre>
     * @param a
     * @param indexLo
     * @param indexHi
     * @return
     */
    public int findMedianOfMediansIdx(float[] a, int indexLo, int indexHi) {
        
        int nItems = (indexHi - indexLo) + 1;

        int nPerGroup = 5;

        if (nItems <= nPerGroup) {
            int i0 = indexLo;
            int i1 = i0 + nPerGroup - 1;
            if (i1 > indexHi) {
                i1 = indexHi;
            }
            float median = findMedian(a, i0, i1);      // 2T(5/2) + 5
            int medianIdx = (i1 + i0) >> 1;
            return medianIdx;
        }

        int nDiv = (int) Math.ceil((float)nItems/nPerGroup);
                                                      //    cost       times
        for (int i = 0; i < nDiv; i++) {              //                n/5
            int i0 = indexLo + i*nPerGroup;
            int i1 = i0 + nPerGroup - 1;
            if (i1 > indexHi) {
                i1 = indexHi;
            }
            float median = findMedian(a, i0, i1);      // 2T(5/2) + 5
            int medianIdx = (i1 + i0) >> 1;

            float swap = a[medianIdx];
            a[medianIdx] = a[indexLo + i];
            a[indexLo + i] = swap;
        }

        /*
        for nPerGroup = 5
        Total cost =  N/5 * cnst  +   N/(5*5) * cnst + ... while N/(5^x) < 6

            j = math.floor( math.log(N)/math.log(5))

        Total cost = summation from i = 1 to j inclusive of (N/(5*i)) * const

            so it's basically linear so far, then add the next step

        Next we sort on number of items = nDiv which is the j from
        preceding comments.  2T(j/2) + j.
             
        Adding above to below, summary:
            j = math.floor( math.log(N)/math.log(5))

            Total cost =  T(    ∑   (N/(5*i)) * const   )  + 2T(j/2) + j
                           ( ¡=1 to j                   )

            ===> It's a little more than linear on j but less than O(j*lg2(j))

            for N=1000, O(j*lg2(j)) = O(4)
        */

        float median = findMedian(a, indexLo, indexLo + nDiv);
        int medianIdx = (indexLo + nDiv + indexLo) >> 1;

        return medianIdx;
    }

    /**
     * find the median index of a, while performing the same sort (or swap)
     * operations on b and c too though not reading the later.
     <pre> 
     Complexity of run time:

            j = math.floor( math.log(N)/math.log(5))

            Total cost =  T(    ∑   (N/(5*i)) * const   )  + 2T(j/2) + j
                           ( ¡=1 to j                   )

            ===&gt; It's a little more than O(j) but less than O(j*lg(j)).
     </pre>
     * @param a
     * @param b array to perform same swap operations on that a receives
     * @param c array to perform same swap operations on that a receives
     * @param indexLo
     * @param indexHi
     * @return
     */
    public int findMedianOfMediansIdx(float[] a, int[] b, int[] c,
        int indexLo, int indexHi) {
        
        int nItems = (indexHi - indexLo) + 1;

        int nPerGroup = 5;

        if (nItems <= nPerGroup) {
            int i0 = indexLo;
            int i1 = i0 + nPerGroup - 1;
            if (i1 > indexHi) {
                i1 = indexHi;
            }
            float median = findMedian(a, b, c, i0, i1);      // 2T(5/2) + 5
            int medianIdx = (i1 + i0) >> 1;
            return medianIdx;
        }

        int nDiv = (int) Math.ceil((float)nItems/nPerGroup);
                                                      //    cost       times
        for (int i = 0; i < nDiv; i++) {              //                n/5
            int i0 = indexLo + i*nPerGroup;
            int i1 = i0 + nPerGroup - 1;
            if (i1 > indexHi) {
                i1 = indexHi;
            }
            float median = findMedian(a, b, c, i0, i1);      // 2T(5/2) + 5
            int medianIdx = (i1 + i0) >> 1;

            float swap = a[medianIdx];
            a[medianIdx] = a[indexLo + i];
            a[indexLo + i] = swap;
            
            int swap2 = b[medianIdx];
            b[medianIdx] = b[indexLo + i];
            b[indexLo + i] = swap2;
            
            swap2 = c[medianIdx];
            c[medianIdx] = c[indexLo + i];
            c[indexLo + i] = swap2;
        }

        /*
        for nPerGroup = 5
        Total cost =  N/5 * cnst  +   N/(5*5) * cnst + ... while N/(5^x) < 6

            j = math.floor( math.log(N)/math.log(5))

        Total cost = summation from i = 1 to j inclusive of (N/(5*i)) * const

            so it's basically linear so far, then add the next step

        Next we sort on number of items = nDiv which is the j from
        preceding comments.  2T(j/2) + j.
             
        Adding above to below, summary:
            j = math.floor( math.log(N)/math.log(5))

            Total cost =  T(    ∑   (N/(5*i)) * const   )  + 2T(j/2) + j
                           ( ¡=1 to j                   )

            ===> It's a little more than linear on j but less than O(j*lg2(j))

            for N=1000, O(j*lg2(j)) = O(4)
        */

        float median = findMedian(a, b, c, indexLo, indexLo + nDiv);
        int medianIdx = (indexLo + nDiv + indexLo) >> 1;

        return medianIdx;
    }

    /**
     * select the smallest kth number in an unordered list of numbers.
     *
     * the complexity of run time is greater than linear O(N) but less than
     * O(N^2) and it looks like one could show less than O(N*lg2(N)).
     *
     * It should be better than quicksort because it's using the pivot only for
     * values less than or equal to median ( so better than O(NlgN) at best and 
     * better than O(N^2) at worse).
     *
     * @param a
     * @param indexLo
     * @param indexHi
     * @param kthIndex
     * @return
     */
    protected float selectKth(float[] a, int indexLo, int indexHi, int kthIndex) {
        
        //  greater than O(j) but less than O(j*lg(j))  
        //  where j = math.floor( math.log(N)/math.log(5))
        int medianIdx = findMedianOfMediansIdx(a, indexLo, indexHi);

        float roughMedian = a[medianIdx];

        // greater than O(N) but less than O(N^2).  if medianIdx < half of 
        // sample, its closer to O(N).
        int pivotIndex = pivotPartition(a, indexLo, indexHi, medianIdx);

        int j0 = -1;
        int j1 = -1;

        // find the first and last location of median value in array a.
        // O(N)
        for (int i = indexLo; i <= indexHi; i++) {
            if (j0 == -1) {
                if (a[i] == roughMedian) {
                    j0 = i;
                    j1 = i;
                }
            } else {
                if (a[i] == roughMedian) {
                    j1 = i;
                }
            }
        }
    
        // this resembles an indexed binarySearch:
        int k = kthIndex;
        if (k < j0) {
            return selectKth(a, indexLo, j0 - 1, k);
        } else if ((k >= j0) && (k <= j1)) {
            return roughMedian; // same as a[k]
        } else if (k > j1) {
            return selectKth(a, j1 + 1, indexHi, k);
        }

        throw new IllegalArgumentException("did not find solution");
    }
    
    /**
     * given a pivotIndex re-order the array so that all values less than  a[pivotIndex]
     * have final lower indexes and all values higher than a[pivotIndex] have indexes higher.
     * returns the new location of pivotIndex.
     *
     * complexity of runtime is between O(N) best case and O(N^2) worse case
     *
     * @param a
     * @param indexLo
     * @param indexHi
     * @param pivotIndex
     * @return
     */
    public static int pivotPartition(float[] a, int indexLo, int indexHi, int pivotIndex) {
                                                          // cost  times
        float pivotValue = a[pivotIndex];                 // c00
        int i = indexLo - 1;                              // c01

        for (int j = indexLo; j < indexHi ; j++ ) {       // c1     n - 1
            if (a[j] <= pivotValue) {                     // c2     sum(part of 2->n)
                i++;                                      // c3     sum(part of 2->n)
                float swap = a[i];                        // c4     sum(part of 2->n)
                a[i] = a[j];                              // c5     sum(part of 2->n)
                a[j] = swap;                              // c6     sum(part of 2->n)
            }
        }
        
        if (i > (indexLo -1)) {                           //
            //TODO:  index might need to be [i - 1] here
            float swap = a[i];                            // c02
            a[i] = a[pivotIndex];                         // c03
            a[pivotIndex] = swap;                         // c04
        }

        /*  O(partition) = c00 + c01 + (c1*(n-1))
               + (c2*sum(part of 2->n)) + (c3*sum(part of 2->n-1)) + (c4*sum(part of 2->n-1))
               + (c5*sum(part of 2->n-1)) + (c6*sum(part of 2->n-1))
               + c02 + c03 + c04

            Best case:
     *         = (c1 + c2 + c3 + c4 + c5 + c6 )*(n-1) + c00 + c01 + c02 + c03 + c04
     *         = (c1 + c2 + c3 + c4 + c5 + c6 )*n - c1 - c2 - c3 - c4 - c5 - c6 + c00 + c01 + c02 + c03 + c04
     *         ==> LINEAR, so T(N)
     *         ==> then, order of growth results in O(N)
     *
     *     Worse case:
     *         = c00 + c01 + (c1*(n-1))
     *         + (c2*((n*(n+1)/2) - 1)) + (c3*((n*(n11)/2) - 1)) + (c4*((n*(n11)/2) - 1))
     *         + (c5*((n*(n11)/2) - 1)) + (c6*((n*(n11)/2) - 1))
     *         + c02 + c03 + c04
     *         = c00 + c01 + (c1*(n-1)) + (c2*(n^2 + n)/2) + (c3*(n^2 - n)/2) + (c4*(n^2 - n)/2)
     *         + (c5*(n^2 - n)/2) + (c6*(n^2 - n)/2) + c02 + c03 + c04
     *         = (c1 + c2/2 + c3/2 + c4/2 + c5/2 + c6/2)*n  + (c2/2 + c3/2 + c4/2 + c5/2 + c6/2)*n^2
     *            - c1 + c00 + c01 + c02 + c03 + c04 ...
     *         ==> QUADRATIC
     *         ==> then, order of growth gives O(N^2)
        */

        return i;
    }

    /**
     * find the median index for very small (indexHi - indexLo).  Note that it 
     * sorts a from indexLo to indexHi too.
     * @param a
     * @param indexLo
     * @param indexHi
     * @return
     */
    private static float findMedian(float[] a, int indexLo, int indexHi) {
               
        QuickSort.sort(a, indexLo, indexHi);

        return a[(indexHi + indexLo) >> 1];
    } 
    
    /**
     * find the median index for very small (indexHi - indexLo).  Note that it 
     * sorts a from indexLo to indexHi too.
     * @param a
     * @param indexLo
     * @param indexHi
     * @return
     */
    private static float findMedian(float[] a, int[] b, int[] c, 
        int indexLo, int indexHi) {
               
        QuickSort.sort(a, b, c, indexLo, indexHi);

        return a[(indexHi + indexLo) >> 1];
    } 
}
