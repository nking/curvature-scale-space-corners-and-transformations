package packing;

import algorithms.QuickSort;
import algorithms.util.PixelHelper;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;

/**
 * Class to return a set of points for the intersection of 2 point sets where 
 * the set of points are separated by a given interval in x and y, that is,
 * the bins are squares.
 * The two point sets can be irregular in shape.  
 * The solution is not guaranteed to be optimal.
 * 
 * The runtime complexity is expected to be:
 *     O(min(N_points) * lg_2(N_points))
 * 
 * Considering implementing these:
 * 
 * Option 1:
 *     resembles 2D strips:
 *     -- find the intersection of the point sets by iterating over the smallest
 *        point set.
 *     -- sort the intersection data structure by x, then y
 *     -- fill the intersection space with rectangles of given x and y size.
 *        -- use IntervaleRangeSearch for collision checks
 *        
 * Option 2:
 *     find the intersection of the point sets by iterating over the smallest
 *        point set.
 *     use a medial axis built from a voronoi diagram
 *     (see algorithms.compGeometry.MedialAxis), 
 *     then fill the space using order based upon the most connected medial 
 *     axis points and then the space connected to those assigned bins, 
 *     iteratively.
 *     It also uses IntervaleRangeSearch for collision checks.
 * 
 * Option 3:
 *    find the intersection of the point sets by iterating over the smallest
 *        point set.
 *     create ordered, connected boundary points of the
 *     intersection.
 *     then walk along the border, filling in bins
 *        at cell x,y spacings and adding the connected
 *        intervals to a stack iteratively to continue filling
 *        the intersection with bins.
 
 * 
 * @author nichole
 */
public class Intersection2DPacking {
    
    /**
     * uses 2-D strip packing and cell sizes of cellSize to place points 
     * throughout the intersection of the 2-D points in the most naive greedy 
     * placement order in x, then y.  Note that the result is not guaranteed
     * to be optimal.
     * 
     * The runtime complexity is O(N * log_2(N)) 
     * where N is min(N_points1, N_points2).
     * 
     * @param points1
     * @param points2
     * @param imageWidth
     * @param cellSize
     * @return 
     */
    public TIntSet naiveStripPacking(TIntSet points1, TIntSet points2, 
        int imageWidth, int cellSize) {
        
        // O(N)
        TIntSet intersection = intersection(points1, points2);
        
        //O(N)
        int[] xs = new int[intersection.size()];
        int[] ys = new int[intersection.size()];
        _populate(intersection, imageWidth, xs, ys);
        
        //O(N*lg_2(N))
        QuickSort.sortBy1stThen2nd(ys, xs);
        
        TIntSet out = new TIntHashSet();
        
        PixelHelper ph = new PixelHelper();
        int lX = Integer.MIN_VALUE; 
        int lY = Integer.MIN_VALUE;
        int x, y;
        long pixIdx;
       
        for (int i = 0; i < xs.length; ++i) {
            x = xs[i];
            y = ys[i];
            if (x < lX) {
                lX = Integer.MIN_VALUE;
            }
            if ((x >= (lX + cellSize)) && 
                ((y == lY) || (y >= (lY + cellSize)))) {
                
                pixIdx = ph.toPixelIndex(x, y, imageWidth);
                out.add((int)pixIdx);
                lX = x;
                lY = y;
            }
        }
        return out;
    }
    
    /**
     * Find the intersection of the 2 point sets.
     * The runtime complexity is O(N) where N is min(N_points1, N_points2).
     * 
     * @param points1
     * @param points2
     * @return 
     */
    public TIntSet intersection(TIntSet points1, TIntSet points2) {
        
        TIntSet out = new TIntHashSet();
        
        TIntSet p1, p2;
        if (points1.size() <= points2.size()) {
            p1 = points1;
            p2 = points2;
        } else {
            p1 = points2;
            p2 = points1;
        }
        TIntIterator iter = p1.iterator();
        while (iter.hasNext()) {
            int pixIdx = iter.next();
            if (p2.contains(pixIdx)) {
                out.add(pixIdx);
            }
        }
        
        return out;
    }

    void _populate(TIntSet points, int imageWidth, int[] xs, int[] ys) {
    
        PixelHelper ph = new PixelHelper();
        int[] xy = new int[2];
        
        int i = 0;
        TIntIterator iter = points.iterator();
        while (iter.hasNext()) {
            int pixIdx = iter.next();
            ph.toPixelCoords(pixIdx, imageWidth, xy);
            xs[i] = xy[0];
            ys[i] = xy[1];
            ++i;
        }
    }
}
