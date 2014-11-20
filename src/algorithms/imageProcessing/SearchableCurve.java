package algorithms.imageProcessing;

import algorithms.util.PairIntArray;
import algorithms.util.PairInt;
import java.util.Arrays;

/**
 * Stores the curve as (x, y) points and creates an index to make the ability to
 * find the closest point in a runtime complexity near that of binary search, 
 * O(lg_2(N)).
 * 
 * Note, the class has the infrastructure to track visits too
 * so that any point returned by find is removed from future 
 * search results, but that is not enabled at this time.
 * 
 * @author nichole
 */
public class SearchableCurve {
    
    private final int[] xOriginalSorted;
    private final int[] yOriginalSorted;
    private final int[] indexesToOriginalXY;
    
    /**
     * TODO: if edit to store visits: x and y will modified upon a find.
     */
    private int[] x = null;
    private int[] y = null;
    private int[] indexes = null;
    private int xLength = 0;
    
    private boolean isSorted = false;
    
    public int getN() {
        return xLength;
    }
    int[] getX() {
        return x;
    }
    int[] getY() {
        return y;
    }
        
    public SearchableCurve(PairIntArray xy) {
        
        x = Arrays.copyOf(xy.getX(), xy.getN());
        y = Arrays.copyOf(xy.getY(), xy.getN());
        xLength = xy.getN();
        
        indexes = new int[x.length];
        for (int i = 0; i < x.length; i++) {
            indexes[i] = i;
        }
        
        initialize();
        
        xOriginalSorted = Arrays.copyOf(x, x.length);
        yOriginalSorted = Arrays.copyOf(y, y.length);
        indexesToOriginalXY = Arrays.copyOf(indexes, indexes.length);
    }
    
    public void resetVisitsState() {
        x = Arrays.copyOf(xOriginalSorted, xOriginalSorted.length);
        y = Arrays.copyOf(yOriginalSorted, yOriginalSorted.length);
        indexes = Arrays.copyOf(indexesToOriginalXY, indexesToOriginalXY.length);
        xLength = x.length;
    }
    
    private void initialize() {       
        
        if (isSorted) {
            return;
        }
       
        if (xLength > 0) {
            // sort by decreasing y, then increasing x
            // (adopting same pattern as used in contours for t, sigma)
            sortBy2ndThen1st(x, y, indexes);
        }
        
        isSorted = true;
    }
    
    /**
     * use a binary search pattern to find the index of closest (x,y) point
     * w.r.t x and y.  to get the index w.r.t. the constructor's xy instance,
     * use indexesToOriginalXY[this result].
     * Note that the current implementation may not be the closest, but may
     * be the "close" match.  it looks for matching y and the closest x from
     * there.  this may change...
     * 
     * TODO: not implemented yet: Note that a side effect of this method is that the returned point is
     * marked as "visited" and will not be returned upon future visits.  resetVisitsState()
     * will reset the data so that none have been visited
     * 
     * runtime complexity is O(lg_2(N)).
     * 
     * @param xPoint
     * @param yPoint
     * @return 
     */
    protected PairInt findClosestMatch(int xPoint, int yPoint) {
        
        int idx = findClosestMatchBinarySearch(xPoint, yPoint);
        
        if (idx == -1) {
            idx = findClosestMatchForYTooHigh(xPoint, yPoint);
        }
        
        if (idx != -1) {
            //TODO: make a moveUp to delete the point in x,y,and indexes
            //   and add an xLen to not use entries past xLen
            
            PairInt pair = new PairInt(x[idx], y[idx]);
            
            return pair;
        }
        
        return null;
    }
    
    /**
     * use a binary search pattern to find the index of closest (x,y) point
     * w.r.t x and y.  to get the index w.r.t. the constructor's xy instance,
     * use indexesToOriginalXY[this result].
     * 
     * runtime complexity is O(lg_2(N)).
     * @param xPoint
     * @param yPoint
     * @return 
     */
    protected int findClosestMatchBinarySearch(int xPoint, int yPoint) {
        
        //TODO:  this needs alot more testing!  should have same bugs that 
        // contour finding had
        
        int minDiffY = Integer.MAX_VALUE;
        int minDiffX = Integer.MAX_VALUE;
        int idx = -1;
        
        int lowIdx = 0;
        int highIdx = x.length - 1;
        
        while (true) {
            
            if (lowIdx > highIdx) {
                
                break;
            
            } else {
                
                int midIdx = (lowIdx + highIdx) >> 1;
                                
                int diffY = Math.abs(y[midIdx] - yPoint);
                int diffX = Math.abs(x[midIdx] - xPoint);
                
                if ((diffY <= minDiffY) && (diffX <= minDiffX)) {
                    minDiffY = diffY;
                    minDiffX = diffX;
                    idx = midIdx;
                }
                
                if ((diffY == 0) && (diffX == 0)) {
                    
                    idx = midIdx;
                    break;
                
                } else if (diffY == 0) {
                    
                    // smaller x is at a lower index
                    if (x[midIdx] < xPoint) {
                        lowIdx = midIdx + 1;
                    } else {
                        highIdx = midIdx - 1;
                    }
                    
                } else if (y[midIdx] < yPoint) {
                    // they're at a lower index
                    highIdx = midIdx - 1;
                    
                } else {
                    // theyr'e at a higher index
                    lowIdx = midIdx + 1;
                }
            }
        }
        
        if (idx > -1) {
            return idx;
        }
        
        return -1;
    }

    /**
     * A method to be used only if findClosestMatchBinarySearch returns a -1
     * and does so because the point has more than one y which it is near,
     * but a closer x which for the lower value y.
     * The method won't usually find the closest match unless that parameter
     * space has been searched first.
     * return the index for a match w.r.t. arrays x and y, else -1 if not found.  
     * This method runtime complexity is O(N).
     * 
     * @param xPoint
     * @param yPoint
     * @return 
     */
    protected int findClosestMatchForYTooHigh(int xPoint, 
        int yPoint) {
        /*
        1, 5
        2, 4
        3, 3
        4, 3<--
        2, 2
        1, 1
        */
        int minDiffY = Integer.MAX_VALUE;
        int minDiffX = Integer.MAX_VALUE;
        int idx = -1;
        
        for (int i = 0; i < x.length; i++) {
                        
            int diffY = yPoint - y[i];
            
            int diffX = Math.abs(x[i] - xPoint);
            
            if ((diffX == 0) && (diffY == 0)) {
                return i;
            }
            
            /*
            we're looking for cases where y is too high, but the 
            x is close.
            */
            if ((diffY > 0) && (diffY <= minDiffY) 
                && (diffX <= minDiffX)) {
                
                minDiffY = diffY;
                minDiffX = diffX;
                idx = i;
            }
        }
      
        return idx;
    }
    
    /**
     * sort by decreasing a2 then increasing a1 when a2 is the same in comparison.
     * runtime complexity is O(N*lg_2(N)) + O(N/2).
     * @param a1
     * @param a2
     * @param a3
     * @param idxLo
     * @param idxHi 
     */
    private void sortBy2ndThen1st(int[] a1, int[] a2, int[] a3) {
        
        sortBy2ndThen1st(a1, a2, a3, 0, xLength - 1);
            
        // reverse the arrays
        reverse(a1, a2, a3, 0, xLength - 1);
    }
    
    /**
     * sort by decreasing a2 then increasing a1 when a2 is the same in comparison.
     * runtime complexity is O(N*lg_2(N)) + O(N/2).
     * @param a1
     * @param a2
     * @param a3
     * @param idxLo
     * @param idxHi 
     */
    private void sortBy2ndThen1st(int[] a1, int[] a2, int[] a3, int idxLo,
        int idxHi) {
        
        if (idxLo < idxHi) {
            
            int idxMid = (idxLo + idxHi) >> 1;
            
            sortBy2ndThen1st(a1, a2, a3, idxLo, idxMid);
            sortBy2ndThen1st(a1, a2, a3, idxMid + 1, idxHi);
            
            mergeBy2ndThen1st(a1, a2, a3, idxLo, idxMid, idxHi);
            
        } 
    }

    /**
     * sort by increasing a2 and when a2's are equal, sort by decreasing a1.
     * @param a1
     * @param a2
     * @param a3
     * @param idxLo
     * @param idxMid
     * @param idxHi 
     */
    private void mergeBy2ndThen1st(int[] a1, int[] a2, int[] a3, int idxLo, 
        int idxMid, int idxHi) {
        
        int[] leftX = Arrays.copyOfRange(a1, idxLo, idxMid + 2);
        int[] leftY = Arrays.copyOfRange(a2, idxLo, idxMid + 2);
        int[] leftI = Arrays.copyOfRange(a3, idxLo, idxMid + 2);
        leftX[leftX.length - 1] = Integer.MAX_VALUE;
        leftY[leftY.length - 1] = Integer.MAX_VALUE;
        leftI[leftI.length - 1] = Integer.MAX_VALUE;
        
        int[] rightX = Arrays.copyOfRange(a1, idxMid + 1, idxHi + 2);
        int[] rightY = Arrays.copyOfRange(a2, idxMid + 1, idxHi + 2);
        int[] rightI = Arrays.copyOfRange(a3, idxMid + 1, idxHi + 2);        
        rightX[rightX.length - 1] = Integer.MAX_VALUE;
        rightY[rightY.length - 1] = Integer.MAX_VALUE;
        rightI[rightI.length - 1] = Integer.MAX_VALUE;
        
        int leftPos = 0;                                               
        int rightPos = 0;
        
        for (int k = idxLo; k <= idxHi; k++) {
            int lY = leftY[leftPos];
            int rY = rightY[rightPos];
            if (lY == rY) { 
                int lX = leftX[leftPos];
                int rX = rightX[rightPos];
                if (lX > rX) {
                    a2[k] = lY;
                    a1[k] = lX;
                    a3[k] = leftI[leftPos];
                    leftPos++;
                } else {
                    a2[k] = rY;
                    a1[k] = rX;
                    a3[k] = rightI[rightPos];
                    rightPos++;
                }
            } else if (lY < rY) {
                a2[k] = lY;
                a1[k] = leftX[leftPos];
                a3[k] = leftI[leftPos];
                leftPos++;
                
            } else {
                a2[k] = rY;
                a1[k] = rightX[rightPos];
                a3[k] = rightI[rightPos];
                rightPos++;
            }
        }
    }

    /**
     * reverse the order of points in the arrays
     * @param a1
     * @param a2
     * @param a3 
     */
    private void reverse(int[] a1, int[] a2, int[] a3, int idxLo, int idxHi) {
        
        if (a1 == null) {
                throw new IllegalArgumentException("a1 cannot be null");
        }
        if (a2 == null) {
                throw new IllegalArgumentException("a2 cannot be null");
        }
        if (a3 == null) {
                throw new IllegalArgumentException("a3 cannot be null");
        }
        if ((a1.length != a2.length) || (a1.length != a3.length)) {
                throw new IllegalArgumentException("arrays must be same length");
        }
        if (a1.length < 2) {
            return;
        }
        
        int n = idxHi - idxLo + 1;
        
        int end = n >> 1;
        
        for (int i = idxLo; i < (idxLo + end); i++) {
            int idx2 = n - i - 1;
            int swap = a1[i];
            a1[i] = a1[idx2];
            a1[idx2] = swap;
            
            swap = a2[i];
            a2[i] = a2[idx2];
            a2[idx2] = swap;
            
            swap = a3[i];
            a3[i] = a3[idx2];
            a3[idx2] = swap;
        }
    }
}
