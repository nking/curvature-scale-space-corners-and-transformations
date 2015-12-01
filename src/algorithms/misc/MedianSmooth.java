package algorithms.misc;

import algorithms.imageProcessing.GreyscaleImage;
import algorithms.util.PairIntArray;
import java.util.Arrays;

/**
 * class to calculate a running median of the k previous points of curveY.
 * 
 * @author nichole
 */
public class MedianSmooth {
  
    /**
     * calculate a running median of a window of size xWindow, yWindow.
     * runtime complexity is
     *     n_rows * ((xWindow * yWindow) + ((n_cols)*lg2(xWindow * yWindow)))
     * so is roughly O(n_pixels * lg2(window area)) where n_pixels = n_rows * n_cols
     * 
     * NOTE: the border points at positions half the window size are not yet
     * determined.
     * 
     * @param input
     * @param xWindow
     * @param yWindow
     * @return 
     */
    public GreyscaleImage calculate(GreyscaleImage input, int xWindow, 
        int yWindow) {
        
        if (input == null) {
            throw new IllegalArgumentException("curveY cannot be null");
        }
        if (input.getWidth() < xWindow) {
            throw new IllegalArgumentException(
            "input.getWidth() must be equal to or greater than xWindow");
        }
        if (input.getHeight() < yWindow) {
            throw new IllegalArgumentException(
            "input.getHeight() must be equal to or greater than yWindow");
        }
                       
        int nW = xWindow * yWindow;
        
        int xh = xWindow/2;
        int yh = yWindow/2;
               
        //NOTE: to use zero-padding: output = input.createWithDimensions();
        GreyscaleImage output = input.copyImage();
                
        for (int row = 0; row <= (input.getHeight() - yWindow); ++row) {
            
            SortedVector sVec = new SortedVector(nW);
            
            // add the first nW to the sorted vector
            for (int i = 0; i < xWindow; ++i) {
                for (int j = row; j < (row + yWindow); ++j) {
                    sVec.append(input.getValue(i, j));
                }
            }
                        
            //O(k) + (N)*lg2(k)
            int median;
        
            for (int i = (xWindow - 1); i < input.getWidth(); ++i) {
          
                //O(1)
                median = sVec.getMedian();
                
                output.setValue(i - xh, row + yh, median);

                // remove each item from last column in window
                // and add each item in next column for window,

                if ((i + 1) < input.getWidth()) {
                    
                    for (int j = row; j < (row + yWindow); ++j) {

                        // remove : O(log_2(k))
                        sVec.remove(input.getValue(i - xWindow + 1, j));

                        // add : O(log_2(k)) + < O(k)
                        sVec.insertIntoOpenSlot(input.getValue(i + 1, j));
                    }
                }
            }            
        }        
           
        return output;                        
    }
    
    /**
     * calculate a running median of the k previous points of curveY.
     * runtime complexity is O(k) + O(N*lg_2(k)) at most.
     * @param curveY
     * @param kPoints
     * @return 
     */
    public int[] calculate(int[] curveY, int kPoints) {
        
        if (curveY == null) {
            throw new IllegalArgumentException("curveY cannot be null");
        }
        if (curveY.length < kPoints) {
            throw new IllegalArgumentException(
            "curveY.length must be equal to or greater than kPoints");
        }
        
        int[] medians = new int[curveY.length - kPoints + 1];
        
        SortedVector sVec = new SortedVector(kPoints);
        
        // add the first k-1 to the list container
        for (int i = 0; i < (kPoints - 1); i++) {
            sVec.append(curveY[i]);
        }
        //O(k) + + (N)*lg2(k)
        int median;
        
        for (int i = (kPoints - 1); i < curveY.length; i++) {
          
            // add the kth item to the list: O(log_2(k)) + < O(k)
            // the list state is sorted at the end of the method.
            sVec.insertIntoOpenSlot(curveY[i]);
            
            //O(1)
            median = sVec.getMedian();
            
            int idx = i - kPoints + 1;
                        
            // remove the x[i - k + 1] item from sorted list : O(log_2(k))
            sVec.remove(curveY[idx]);
            
            medians[idx] = median;
        }
        
        return medians;
    }
    
    /**
     * calculate a running median of the k previous points of curve.
     * runtime complexity is O(N*k) at most.
     * @param curve
     * @param kPoints
     * @return 
     */
    public PairIntArray calculate(PairIntArray curve, int kPoints) {
        
        if (curve == null) {
            throw new IllegalArgumentException("curve cannot be null");
        }
        if (curve.getN() < kPoints) {
            throw new IllegalArgumentException(
            "curve length must be equal to or greater than kPoints");
        }
        
        PairIntArray medians = new PairIntArray(curve.getN() - kPoints + 1);
        
        SortedVector sVec = new SortedVector(kPoints);
        
        long xSum = 0;
        
        // add the first k-1 to the list container
        for (int i = 0; i < (kPoints - 1); ++i) {
            
            sVec.append(curve.getY(i));
            
            xSum += curve.getX(i);
        }
        
        int median;
        
        for (int i = (kPoints - 1); i < curve.getN(); ++i) {
          
            // add the kth item to the list: O(log_2(k)) + < O(k)
            // the list state is sorted at the end of the method.
            sVec.insertIntoOpenSlot(curve.getY(i));
            
            //O(1)
            median = sVec.getMedian();
            
            int idx = i - kPoints + 1;
                        
            // remove the x[i - k + 1] item from sorted list : O(log_2(k))
            sVec.remove(curve.getY(idx));
            
            //TODO: can be improved if know that the points are evenly sampled
            xSum += curve.getX(i);
            if (idx > 0) {
                xSum -= curve.getX(idx - 1);
            }
            
            medians.add((int)(xSum/kPoints), median);
        }
        
        return medians;
    }
    
    /**
     * a fixed size list that keeps the contents sorted after the capacity is
     * reached.  points are added one at a time and removed one at a time.
     */
    static class SortedVector {
        
        protected final int[] a;
        
        protected int n;
        
        protected int availSlot;
        
        protected boolean sorted;
        
        public SortedVector(int size) {
            
            a = new int[size];
            
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
        public void append(int value) {
            
            if (n == (a.length)) {
                throw new IllegalArgumentException(
                    "there must be an empty slot in order to append." +
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
        public void insertIntoOpenSlot(int value) {
            
            if (n != (a.length - 1)) {
                throw new IllegalArgumentException(
                "the method is meant to be used when there is only one avail slot");
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
         * NOTE: this could be made O(1) runtime complexity at the expense
         * of 3 * space complexity.
         * @param value
         */
        public void remove(int value) {
            
            if (n != a.length) {
                throw new IllegalArgumentException(
                "the method is meant to be used only on a full list");
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
         * get median from the internal array.  Note that this will
         * throw an IllegalArgumentException if the list is not full.
         * runtime is O(1)
         * @return median
         */
         public int getMedian() {
             
            if (n != a.length) {
                throw new IllegalArgumentException(
                "the method is meant to be used only on a full list");
            }
        
            int midIdx = ((n & 1) == 1) ? n/2 : (n - 1)/2;
            
            return a[midIdx];
        }
    }
}
