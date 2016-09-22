package algorithms.misc;

import java.util.Arrays;

/**
 * class to calculate a running maximum of input in sliding window of size k.
 * 
 * It uses the same pattern as MedianSmooth, but is tailored for float data type.
 * 
 * TODO: put both MedianSmooth and this class into common class one day.
 *
 * @author nichole
 */
public class StatsInSlidingWindow {
    
    /**
     * calculate a running maximum of a window of size xWindow, yWindow.
     * runtime complexity is
     *     n_rows * ((xWindow * yWindow) + ((n_cols)*lg2(xWindow * yWindow)))
     * so is roughly O(n_pixels * lg2(window area)) where n_pixels = n_rows * n_cols
     *
     * NOTE: should only be used by a single thread.
     * 
     * NOTE: the border points outside of the window retain their 
     * initial values.
     *
     * @param input
     * @param output
     * @param xWindow
     * @param yWindow
     */
    public void calculateMaximum(float[][] input, float[][] output, int xWindow,
        int yWindow) {

        if (input == null) {
            throw new IllegalArgumentException("input cannot be null");
        }
        if (input.length < xWindow) {
            throw new IllegalArgumentException(
            "input.length must be equal to or greater than xWindow");
        }
        if (input[0].length < yWindow) {
            throw new IllegalArgumentException(
            "input[0].length must be equal to or greater than yWindow");
        }
        if (output == null) {
            throw new IllegalArgumentException("input cannot be null");
        }
        if (input.length != output.length) {
            throw new IllegalArgumentException(
            "input.length must be equal to output.length");
        }
        if (input[0].length != output[0].length) {
            throw new IllegalArgumentException(
            "input[0].length must be equal to output[0].length");
        }
        
        // tailored for maximum stats
        int xLen2 = input.length + xWindow;
        int yLen2 = input[0].length + yWindow;
        float[][] input2 = new float[xLen2][];
        for (int i = 0; i < input2.length; ++i) {
            if (i < input.length) {
                input2[i] = Arrays.copyOf(input[i], yLen2);
            } else {
                input2[i] = new float[yLen2];
            }
        }

        int nW = xWindow * yWindow;

        int xh = xWindow/2;
        int yh = yWindow/2;
        
        boolean xHEven = (xh & 1) == 0;
        boolean yHEven = (yh & 1) == 0;
        
        int h = input2[0].length;
        int w = input2.length;
        int w0 = input.length;
        int h0 = input[0].length;
        
        for (int row = 0; row <= (h - yWindow); ++row) {

            int jIdx = row;
            
            if (jIdx >= h0) {
                break;
            }
            
            SortedVector sVec = new SortedVector(nW);

            // add the first nW to the sorted vector
            for (int i = 0; i < xWindow; ++i) {
                for (int j = row; j < (row + yWindow); ++j) {
                    sVec.append(input2[i][j]);
                }
            }
            
            assert(sVec.n == sVec.a.length);
            assert(sVec.sorted);

            //O(k) + (N)*lg2(k)
            float maximum;
            
            for (int i = (xWindow - 1); i < w; ++i) {
                
                int iIdx = i - xWindow + 1;
             
                if (iIdx >= w0) {
                    break;
                }
                
                //O(1)
                maximum = sVec.getMaximum();
                
                output[iIdx][jIdx] = maximum;

                // remove each item from last column in window
                // and add each item in next column for window,

                if ((i + 1) < w) {
                    
                    for (int j = row; j < (row + yWindow); ++j) {

                        assert(sVec.n == sVec.a.length);

                        // remove : O(log_2(k))
                        sVec.remove(input2[i - xWindow + 1][j]);

                        assert(sVec.n == (sVec.a.length - 1));
                      
                        // add : O(log_2(k)) + < O(k)
                        sVec.insertIntoOpenSlot(input2[i + 1][j]);
                       
                        assert(sVec.n == sVec.a.length);
                    }       
                }
            }                        
        }
        
    }
    
    /**
     * a fixed size list that keeps the contents sorted after the capacity is
     * reached.  points are added one at a time and removed one at a time
     * and there are rules to prevent removing when list is not full or
     * adding when list is full.
     */
    static class SortedVector {
        
        protected final float[] a;

        protected int n;

        protected int availSlot;

        protected boolean sorted;

        public SortedVector(int size) {

            a = new float[size];

            n = 0;

            availSlot = -1;

            sorted = false;
        }

        /**
         * append item value onto the end of the list.  Note that if the item
         * is added to the last slot, the list is immediately sorted into
         * ascending numerical order
         * afterwards as a side effect to keep the logic in the other
         * methods consistent.
         * runtime is usually O(1), but if append is used for the last item,
         * there is a sort adding O(N*log_2(N)).
         * For best use, append(v) the first size-1 items and thereafter use
         * insertIntoOpenSlot(v).
         *
         * @param value
         */
        public void append(float value) {

            if (n == (a.length)) {
                throw new IllegalArgumentException(
                    "0) there must be an empty slot in order to append." +
                    " remove and item then try insert again or construct larger list.");
            }

            a[n] = value;

            n++;

            if (n == a.length) {

                Arrays.sort(a);

                sorted = true;
            }
        }

        /**
         * Insert the value into the list while maintaining the sorted state
         * of the list.  Note that if there is not exactly one available slot
         * in the list, an IllegalArgumentException will be thrown.
         * runtime is usually O(log_2(N)) + less than O(N), but once per class lifetime
         * the sort may occur here adding O(N*log_2(N)).
         * @param value
         */
        public void insertIntoOpenSlot(float value) {

            if (n != (a.length - 1)) {
                String err = "1) the method is meant to be used only on a full list." 
                + " a.length=" + a.length + " n=" + n;
                throw new IllegalArgumentException(err);
            }

            if (!sorted) {
                // this can happen if the user used "size - 1" append()s followed
                // by insertIntoOpenSlot.  It's only needed once for lifetime
                // of object.

                if (availSlot != -1) {
                    throw new IllegalArgumentException(
                        "Error in the algorithm... should have been sorted already");
                }

                a[n] = value;

                n++;

                Arrays.sort(a);

                sorted = true;
                
                return;
            }

            int insIdx = Arrays.binarySearch(a, value);
            if (insIdx < 0) {
                insIdx *= -1;
                insIdx--;
            }

            if (insIdx == availSlot) {

                a[availSlot] = value;

            } else if (insIdx < availSlot) {

                // move all items from insIdx to availSlot down by 1
                for (int i = (availSlot - 1); i >= insIdx; i--) {
                    a[i + 1] = a[i];
                }

                a[insIdx] = value;

            } else {

                int end = insIdx - 1;

                // move items up from availSlot +1 to insIdx - 1
                // then insert value into insIdx - 1
                for (int i = availSlot; i < end; i++) {
                    a[i] = a[i + 1];
                }

                a[insIdx - 1] = value;
            }
            n++;
            availSlot = -1;            
        }

        /**
         * remove the item from the full list of items.
         * runtime is O(log_2(N)).
         * NOTE: this could be made O(1) runtime complexity 
         * at the expense
         * of 3 * space complexity.
         * @param value
         */
        public void remove(float value) {

            if (n != a.length) {
                throw new IllegalArgumentException(
                "2) the method is meant to be used only on a full list." 
                + " a.length=" + a.length + " n=" + n);
            }

            int rmIdx = Arrays.binarySearch(a, value);

            if (rmIdx < 0) {
                throw new IllegalArgumentException("could not find item in list");
            }

            availSlot = rmIdx;

            // to keep the list in a state where the next binary search works,
            // set the empty slot value to the proceeding value or max integer.
            if (availSlot == (a.length - 1)) {
                a[availSlot] = Integer.MAX_VALUE;
            } else {
                a[availSlot] = a[availSlot + 1];
            }

            n--;            
        }

        /**
         * get maximum from the internal array.  Note that this will
         * throw an IllegalArgumentException if the list is not full.
         * runtime is O(1)
         * @return median
         */
         public float getMaximum() {

            if (n != a.length) {
                // NOTE: in the use above, this is never invoked unless the
                // list a is full so this exception should never be thrown
                throw new IllegalArgumentException(
                    "3) the method is meant to be used only on a full list." 
                    + " a.length=" + a.length + " n=" + n);
            }

            return a[n - 1];
        }
        
         /**
         * get minimum from the internal array.  Note that this will
         * throw an IllegalArgumentException if the list is not full.
         * runtime is O(1)
         * @return median
         */
         public float getMinimum() {

            if (n != a.length) {
                // NOTE: in the use above, this is never invoked unless the
                // list a is full so this exception should never be thrown
                throw new IllegalArgumentException(
                    "3) the method is meant to be used only on a full list." 
                    + " a.length=" + a.length + " n=" + n);
            }

            return a[0];
        }
         
         /**
         * get median from the internal array.  Note that this will
         * throw an IllegalArgumentException if the list is not full.
         * runtime is O(1)
         * @return median
         */
         public float getMedian() {

            if (n != a.length) {
                // NOTE: in the use above, this is never invoked unless the
                // list a is full so this exception should never be thrown
                throw new IllegalArgumentException(
                    "3) the method is meant to be used only on a full list." 
                    + " a.length=" + a.length + " n=" + n);
            }

            int midIdx = ((n & 1) == 1) ? n/2 : (n - 1)/2;

            return a[midIdx];
        }
    }
}
