package algorithms.compGeometry;

import algorithms.QuickSort;
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
public class NearestPoints1D {
    
    private final int[] a;
    private final int[] originalIndexes;
    
    public NearestPoints1D(int[] values) {
        
        if (values == null) {
            throw new IllegalStateException("values cannot be null");
        }
        
        int n = values.length;
        a = new int[n];
        originalIndexes = new int[n];
        
        for (int i = 0; i < n; ++i) {
            a[i] = values[i];
            originalIndexes[i] = i;
        }
        QuickSort.sortBy1stArg(a, originalIndexes);
    }
    
    public NearestPoints1D(int[] values, int length) {
        
        if (values == null) {
            throw new IllegalStateException("values cannot be null");
        }
        
        int n = length;
        a = new int[n];
        originalIndexes = new int[n];
        
        for (int i = 0; i < n; ++i) {
            a[i] = values[i];
            originalIndexes[i] = i;
        }
        QuickSort.sortBy1stArg(a, originalIndexes);
    }
    
    /**
     * find values within tolerance of value in the contained points.
     * @param value
     * @param tolerance
     * @return 
     */
    public Set<Integer> findNeighbors(int value, float tolerance) {
        
        Set<Integer> values = new HashSet<Integer>();
        
        if (a.length == 0) {
            return values;
        }
        
        Set<Integer> indexes = findNeighborIndexesR(value, tolerance);
        
        for (Integer index : indexes) {
            int i = index.intValue();
            values.add(Integer.valueOf(a[i]));
        }
        
        return values;
    }
    
    /**
     * find indexes of values within tolerance of value in the contained points.
     * @param value
     * @param tolerance
     * @return 
     */
    public Set<Integer> findNeighborIndexes(int value, float tolerance) {
        
        Set<Integer> origIndexes = new HashSet<Integer>();
        
        if (a.length == 0) {
            return origIndexes;
        }
        
        Set<Integer> indexes = findNeighborIndexesR(value, tolerance);
        
        for (Integer index : indexes) {
            int i = index.intValue();
            origIndexes.add(Integer.valueOf(originalIndexes[i]));
        }
        
        return origIndexes;
    }
    
    public int findClosestValue(int value) {
        
        // O(lg2(N))
        int idx = Arrays.binarySearch(a, value);
                    
        // if it's negative, (-(insertion point) - 1)
        if (idx < 0) {
            // idx = -*idx2 - 1
            idx = -1*(idx + 1);
        }
        if (idx > (a.length - 1)) {
            idx = a.length - 1;
        }
        
        // the answer may be in the bin below, so check both
        if (idx > 0) {
            
            int diff0 = Math.abs(a[idx] - value);
        
            int diff1 = Math.abs(a[idx - 1] - value);
            
            if (diff1 < diff0) {
                return a[idx - 1];
            }
        }
        
        return a[idx];
    }
    
    public int findClosestValueIndex(int value) {
        
        // O(lg2(N))
        int idx = Arrays.binarySearch(a, value);
                    
        // if it's negative, (-(insertion point) - 1)
        if (idx < 0) {
            // idx = -*idx2 - 1
            idx = -1*(idx + 1);
        }
        if (idx > (a.length - 1)) {
            idx = a.length - 1;
        }
        
        // the answer may be in the bin below, so check both
        if (idx > 0) {
            
            int diff0 = Math.abs(a[idx] - value);
        
            int diff1 = Math.abs(a[idx - 1] - value);
            
            if (diff1 < diff0) {
                return originalIndexes[idx - 1];
            }
        }
        
        int origIdx = originalIndexes[idx];
        
        return origIdx;
    }
    
    /**
     * find indexes of values within tolerance of value in the contained points
     * return the indexes relative to the a array
     * @param xCenter
     * @param yCenter
     * @param radius
     * @return 
     */
    private Set<Integer> findNeighborIndexesR(int value, float tolerance) {
        
        Set<Integer> relativeIndexes = new HashSet<Integer>();
        
        if (a.length == 0) {
            return relativeIndexes;
        }
        
        // O(lg2(N))
        int idx = Arrays.binarySearch(a, value);
                    
        // if it's negative, (-(insertion point) - 1)
        if (idx < 0) {
            // idx = -*idx2 - 1
            idx = -1*(idx + 1);
        }
        if (idx > (a.length - 1)) {
            idx = a.length - 1;
        }
        
        int startIdx = idx;
        for (int i = (idx - 1); i > -1; --i) {
            int diffX = Math.abs(a[i] - value);
            if (diffX > tolerance) {
                break;
            }
            startIdx = i;
        }
        
        int stopIdx = idx;
        for (int i = idx; i < a.length; ++i) {
            int diffX = Math.abs(a[i] - value);
            if (diffX > tolerance) {
                break;
            }
            stopIdx = i;
        }
               
        // search for points within startIdx and stopIdx that are within radius
        for (int i = startIdx; i <= stopIdx; ++i) {
            relativeIndexes.add(Integer.valueOf(i));
        }
        
        return relativeIndexes;
    }
    
}
