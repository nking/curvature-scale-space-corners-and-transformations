package algorithms.compGeometry.convexHull;

import algorithms.compGeometry.convexHull.GrahamScanLong.CH;

/**
 * finds the pair of points most distant from one another a.k.a. the
 * furthest pair problem.
 * one can imagine a rough estimate by fitting an ellipse around the given points
 * and determining the major axis of the points.
 * even more precisely, one can fit a convex hull around the points and
 * determine which pair of points has the largest distance between them.
 * Note that the convex hull returns points ordered in a counter-clockwise manner
 * so is suitable for ordered traversal already.
 * An efficient way of comparing the point distances on the hull is to start
 * with the first point and determine its distance to the next point i, continuing
 * with increasing i and increasing distance.  that comparison should be
 * roughly h/2 comparisons where h is the number of points in the hull.
 * then for the i points just used in the calculation, do the same analysis
 * to find their most distant point (these points are called antipodal points).
 * the total number of "i" visited are roughly h/2 and each antipodal search is
 * h/2 = O(h^2).
 * The algorithm total runtime complexity is O(convex hull) + O(h^2).
 * 
 * If one uses Graham Scan for the convex hull algorithm, 
 * O(convex hull) = O(N*log_2(N)) where N is the number of points given.
 * If one instead uses Jarvis March for the convex hull algorithm, 
 * O(convex hull) = O(N*log_2(h)) where N is the number of points given and
 * h is the number of points in the resulting convex hull.
 * 
 * If N is greater than h, the total runtime complexity is
 * O(N*log_2(N)) if use the Graham Scan.
 * If N is equal to h, the total runtime complexity is
 * O(N^2) = O(h^2).
 * 
 * We can speculate about the size of the hull from information in Cormen
 * et al. Introduction to Algorithms, Exercise 33-5
 * for sparse-hulled distributions of a unit-radius disk, a convex polygon with k sides,
 * and a 2-D normal distribution respectively as n^(1/3), log_2(n), sqrt(log_2(n)).
 * 
 * @author nichole
 */
public class LargestDistancePair {
    
    /**
     * find the largest distance between the given pairs of points.
     * The runtime complexity is <em>O(N*log_2(N))</em>, 
     * worse case runtime complexity is if the number of hull points
     * is equivalent to the number of points, O(N^2) = O(h^2).
     * 
     * NOTE: one could make an O(N) version of this algorithm assuming h^2 is
     * less than N, by rewriting the bottleneck of the algorithm, the sorting.
     * For the polar angle counter-clockwise sorting of points, one could
     * use counting sort and integer degree bands for the angles.
     * Counting sort runtime complexity is O(max(N, range of data)), so
     * if the range of data were constrained to 0 to 359 for N >= 360,
     * the sort could be done in less than O(N) 
     * and the remaining hull and largest pair algorithm runtime complexity is still O(N).  
     * For the counting sort, if N were less than 360, and O(N) were more appealing than accuracy,
     * the 0 to 359 intervals could be made into a smaller number than 360 of evenly
     * sized bands.
     * 
     * @param x input array of x coordinates.  note that this method modifies the order of x
     * upon CCW sorting of x and y, so copy the arrays in if need to kep the
     * original order.
     * @param y x input array of x coordinates.  note that this method modifies the order of x
     * upon CCW sorting of x and y, so copy the arrays in if need to kep the
     * original order.
     * @return an array of the pair of points furthest from one another.  the
     * format in the return array is [xa, ya, yb, yb].
     * If there are more than one points with same maximum distance,
     * only one is returned, and that one is the lastest point in roughly the 1st half of the
     * convex hull and it last pairing in terms of CCW ordering of the hull points.
     * @throws algorithms.compGeometry.convexHull.GrahamScanTooFewPointsException
     * if the unique polar angles w.r.t. the smallest y point in the arrays
     * are fewer than 3 in number.
     */
    public static long[] findLargestDistancePair(long[] x, long[] y) throws GrahamScanTooFewPointsException {
        
        // x and y
        CH ch = GrahamScanLong.computeHull(x, y);
                
        long maxDist = Long.MIN_VALUE;
        int iMaxDist = 0;
        long distSq;
        long xd;
        long yd;
        
        long x0 = ch.getXH()[0];
        long y0 = ch.getYH()[0];
        
        int i0 = 0;
        
        // find furthest point from the first point.  that will be the scan range;
        int i;
        // scan to the 2nd to last point because the last point in the hull is the same as the first point
        for (i = 1; i < ch.getXH().length - 1; ++i) {
            xd = ch.getXH()[i] - x0;
            yd = ch.getYH()[i] - y0;
            distSq = xd * xd + yd * yd;
            
            if (distSq >= maxDist) {
                maxDist = distSq;
                iMaxDist = i;
            } else {
                // because the hull points are sorted in CCW order, the max dist
                // will increase or stay the same and then decrease.
                break;
            }
        }
        
        int nScan = i;
        
        int j;
        
        // try pairs i0: [2, nScan) to see if an i1 pairing has dist > maxDist
        for (i = 1; i < nScan; ++i) {
            x0 = ch.getXH()[i];
            y0 = ch.getYH()[i];
            for (j = i + 1; j < ch.getXH().length - 1; ++j) {
                xd = ch.getXH()[j] - x0;
                yd = ch.getYH()[j] - y0;
                distSq = xd * xd + yd * yd;
                if (distSq >= maxDist) {
                    i0 = i;
                    maxDist = distSq;
                    iMaxDist = j;
                } else if (j > (i + nScan)) {
                    //TODO: consider letting the scan continue further than this along the hull
                    break;
                }
            }
        }
        
        assert(maxDist > Long.MIN_VALUE);
        
        return new long[]{x[i0], y[i0], x[iMaxDist], y[iMaxDist]};
    }
}
