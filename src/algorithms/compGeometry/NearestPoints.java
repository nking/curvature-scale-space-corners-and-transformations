package algorithms.compGeometry;

import algorithms.MultiArrayMergeSort;
import algorithms.util.PairInt;
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
public class NearestPoints {
    
    private final int[] x;
    private final int[] y;
    
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
        
        x = Arrays.copyOf(xPoints, xPoints.length);
        y = Arrays.copyOf(yPoints, yPoints.length);
        
        MultiArrayMergeSort.sortBy1stArgThen2nd(x, y);
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
        for (int i = idx; i > -1; --i) {
            int x2 = x[i];
            if (Math.abs(x2 - xCenter) > radius) {
                break;
            }
            startIdx = i;
        }
        
        int stopIdx = idx;
        for (int i = (idx + 1); i < x.length; ++i) {
            int x2 = x[i];
            if (Math.abs(x2 - xCenter) > radius) {
                break;
            }
            stopIdx = i;
        }
        
        double rSq = radius * radius;
        
        // search for points within startIdx and stopIdx that are within radius
        for (int i = startIdx; i <= stopIdx; ++i) {
            int diffX = x[i] - xCenter;
            int diffY = y[i] - yCenter;
            double distSq = (diffX*diffX) + (diffY*diffY);
            if (distSq <= rSq) {
                result.add(new PairInt(x[i], y[i]));
            }
        }
        
        return result;
    }
}
