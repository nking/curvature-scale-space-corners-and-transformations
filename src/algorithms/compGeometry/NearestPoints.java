package algorithms.compGeometry;

import algorithms.sort.MultiArrayMergeSort;
import algorithms.util.PairInt;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

/**
 * Builds an unmodifiable data structure to make finding member points within
 * a distance of a location faster than O(N^2).
 * 
 * The runtime complexity is O(N*lg2(N)) for the constructor and for search is 
 * O(lg2(N)) + small number of scans surrounding the nearest x.
 * 
 * @author nichole
 */
public class NearestPoints {
    
    private final int[] x;
    private final int[] y;
    private final int[] originalIndexes;
    
    public NearestPoints(int[] xPoints, int[] yPoints) {
        
        if (xPoints == null) {
            throw new IllegalStateException("xPoints cannot be null");
        }
        if (yPoints == null) {
            throw new IllegalStateException("yPoints cannot be null");
        }
        if (xPoints.length != yPoints.length) {
            throw new IllegalStateException(
            "xPoints and yPoints must be the same length");
        }
        
        int n = xPoints.length;
        x = new int[n];
        y = new int[n];
        originalIndexes = new int[n];
        
        for (int i = 0; i < n; ++i) {
            x[i] = xPoints[i];
            y[i] = yPoints[i];
            originalIndexes[i] = i;
        }
        
        MultiArrayMergeSort.sortBy1stArgThen2nd(x, y, originalIndexes);
    }
    
    /**
     * find points within radius of (xCenter, yCenter) in the contained points.
     * @param xCenter
     * @param yCenter
     * @param radius
     * @return 
     */
    public Set<PairInt> findNeighbors(int xCenter, int yCenter, float radius) {
        
        Set<PairInt> result = new HashSet<PairInt>();
        
        if (x.length == 0) {
            return result;
        }
        
        Set<Integer> indexes = findNeighborIndexesR(xCenter, yCenter, radius);
        
        for (Integer index : indexes) {
            int i = index.intValue();
            result.add(new PairInt(x[i], y[i]));
        }
        
        return result;
    }
    
    /**
     * find points within radius of (xCenter, yCenter) in the contained points.
     * @param xCenter
     * @param yCenter
     * @param radius
     * @return 
     */
    public Set<Integer> findNeighborIndexes(int xCenter, int yCenter, float radius) {
        
        Set<Integer> result = new HashSet<Integer>();
        
        if (x.length == 0) {
            return result;
        }
        
        Set<Integer> indexes = findNeighborIndexesR(xCenter, yCenter, radius);
        
        for (Integer index : indexes) {
            int i = index.intValue();
            result.add(originalIndexes[i]);
        }
        
        return result;
    }
    
    /**
     * find points within radius of (xCenter, yCenter) in the contained points
     * and return the indexes relative to the x,y arrays
     * @param xCenter
     * @param yCenter
     * @param radius
     * @return 
     */
    private Set<Integer> findNeighborIndexesR(int xCenter, int yCenter, float radius) {
        
        Set<Integer> resultIndexes = new HashSet<Integer>();
        
        if (x.length == 0) {
            return resultIndexes;
        }
        
        // O(lg2(N))
        int idx = Arrays.binarySearch(x, xCenter);
     
        // if it's negative, (-(insertion point) - 1)
        if (idx < 0) {
            // idx = -*idx2 - 1
            idx = -1*(idx + 1);
        }
        if (idx > (x.length - 1)) {
            idx = x.length - 1;
        }
        
        double rSq = Math.sqrt(2) * radius * radius;
        
        int startIdx = idx;
        for (int i = (idx - 1); i > -1; --i) {
            int diffX = Math.abs(x[i] - xCenter);
            if (diffX > rSq) {
                break;
            }
            startIdx = i;
        }
        
        int stopIdx = idx;
        for (int i = idx; i < x.length; ++i) {
            int diffX = Math.abs(x[i] - xCenter);
            if (diffX > rSq) {
                break;
            }
            stopIdx = i;
        }
              
        // search for points within startIdx and stopIdx that are within radius
        for (int i = startIdx; i <= stopIdx; ++i) {
            int diffX = x[i] - xCenter;
            int diffY = y[i] - yCenter;
            double distSq = (diffX*diffX) + (diffY*diffY);
            if (distSq <= rSq) {
                resultIndexes.add(Integer.valueOf(i));
            }
        }
        
        return resultIndexes;
    }
    
    public PairInt getSmallestXY() {
        return new PairInt(x[0], y[0]);
    }
}
