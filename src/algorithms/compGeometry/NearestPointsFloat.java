package algorithms.compGeometry;

import algorithms.sort.MultiArrayMergeSort;
import algorithms.util.PairFloat;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

/**
 * Builds an unmodifiable data structure to make finding member points within
 * a distance of a location faster than O(N^2).
 * 
 * The runtime complexity is O(N*lg2(N)) for the constructor
 * and for search is O(lg2(N)) + small number of scans surrounding the nearest x.
 * 
 * @author nichole
 */
public class NearestPointsFloat {
    
    private final float[] x;
    private final float[] y;
    private final int[] originalIndexes;
    
    public NearestPointsFloat(float[] xPoints, float[] yPoints, int len) {
        
        if (xPoints == null) {
            throw new IllegalStateException("xPoints cannot be null");
        }
        if (yPoints == null) {
            throw new IllegalStateException("yPoints cannot be null");
        }
        if (xPoints.length < len) {
            throw new IllegalStateException("xPoints.length is less than len");
        }
        if (yPoints.length < len) {
            throw new IllegalStateException("yPoints.length is less than len");
        }
        
        int n = len;
        x = new float[n];
        y = new float[n];
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
    public Set<PairFloat> findNeighbors(float xCenter, float yCenter, float radius) {
        
        Set<PairFloat> result = new HashSet<PairFloat>();
        
        if (x.length == 0) {
            return result;
        }
        
        Set<Integer> indexes = findNeighborIndexesR(xCenter, yCenter, radius);
        
        for (Integer index : indexes) {
            int i = index.intValue();
            result.add(new PairFloat(x[i], y[i]));
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
    public Set<Integer> findNeighborIndexes(float xCenter, float yCenter, float radius) {
        
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
     * and return the indexes relative to the x,y arrays (not the original
     * array given at construct time).
     * @param xCenter
     * @param yCenter
     * @param radius
     * @return 
     */
    private Set<Integer> findNeighborIndexesR(float xCenter, float yCenter, float radius) {
        
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
        
        int startIdx = idx;
        for (int i = (idx - 1); i > -1; --i) {
            float diffX = Math.abs(x[i] - xCenter);
            if (diffX > radius) {
                break;
            }
            startIdx = i;
        }
        
        int stopIdx = idx;
        for (int i = idx; i < x.length; ++i) {
            float diffX = Math.abs(x[i] - xCenter);
            if (diffX > radius) {
                break;
            }
            stopIdx = i;
        }
       
        double rSq = Math.sqrt(2) * radius * radius;
        
        // search for points within startIdx and stopIdx that are within radius
        for (int i = startIdx; i <= stopIdx; ++i) {
            float diffX = x[i] - xCenter;
            float diffY = y[i] - yCenter;
            double distSq = (diffX*diffX) + (diffY*diffY);
            if (distSq <= rSq) {
                resultIndexes.add(Integer.valueOf(i));
            }
        }
        
        return resultIndexes;
    }
    
    /**
     * find points within radius of (xCenter, yCenter) in the contained points
     * and return the indexes relative to original given array.
     * @param xCenter
     * @param yCenter
     * @param radius
     * @return 
     */
    public Integer findClosestNeighborIndex(float xCenter, float yCenter, float radius) {
                
        if (x.length == 0) {
            return null;
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
        
        int startIdx = idx;
        for (int i = (idx - 1); i > -1; --i) {
            float diffX = Math.abs(x[i] - xCenter);
            if (diffX > radius) {
                break;
            }
            startIdx = i;
        }
        
        int stopIdx = idx;
        for (int i = idx; i < x.length; ++i) {
            float diffX = Math.abs(x[i] - xCenter);
            if (diffX > radius) {
                break;
            }
            stopIdx = i;
        }
       
        double rSq = Math.sqrt(2) * radius * radius;
        
        double minDistSq = Double.MAX_VALUE;
        int minDistIdx = -1;
        
        // search for points within startIdx and stopIdx that are within radius
        for (int i = startIdx; i <= stopIdx; ++i) {
            float diffX = x[i] - xCenter;
            float diffY = y[i] - yCenter;
            double distSq = (diffX*diffX) + (diffY*diffY);
            if (distSq <= rSq) {
                if (distSq < minDistSq) {
                    minDistSq = distSq;
                    minDistIdx = originalIndexes[i];
                }
            }
        }
        
        if (minDistIdx > -1) {
            return Integer.valueOf(minDistIdx);
        }
        
        return null;
    }
}
