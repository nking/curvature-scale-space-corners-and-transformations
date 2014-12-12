package algorithms.compGeometry.convexHull;

import algorithms.compGeometry.LinesAndAngles;
import algorithms.MultiArrayMergeSort;
import java.util.Arrays;

/**
 * adapted from 
 * https://code.google.com/p/two-point-correlation/source/browse/src/test/java/algorithms/compGeometry/convexHull/
 * under MIT License (MIT), Nichole King 2013
 * 
  <pre>
  Solves the Convex Hull problem w/ a stack S of candidate points.
 
  Given a set of Q points returns the vertices of the ConvexHull(Q) in counterclockwise
  order.   a convex hull is the smallest convex polygon that will include all points in Q.
 
  Graham's Scan runs in O(n lg n).
    (in contrast to Jarvis's March which runs in O(nh) where h is the number of
    vertices in the convex hull.)
  Will adjust this after estimates...
 
  Both use a technique called 'rotational sweep' to process vertices in the order
  of the polar angles they form with a reference vertex.
 
  constructed from pseudo-code in Cormen et al. "Introduction to Algorithms
 </pre>
 
 * @author nichole
 */
public class GrahamScan {

    protected XYStack points = null;

    protected float[] xHull = null;
    protected float[] yHull = null;

	public GrahamScan() {
	}

    /**
     * find the convex hull of the given (x, y) points.  Note that the resulting
     * hull points have the same first point as last point.
     * 
     * @param x
     * @param y
     * @throws GrahamScanTooFewPointsException 
     */
    public void computeHull(float[] x, float[] y) throws GrahamScanTooFewPointsException {

        if (x == null) {
	    	throw new IllegalArgumentException("xArray cannot be null");
        }
	    if (y == null) {
	    	throw new IllegalArgumentException("yArray cannot be null");
        }
	    if (x.length != y.length) {
	    	throw new IllegalArgumentException("xArray must have the same number of items as yArray");
        }
	    if (x.length < 3) {
	        throw new IllegalArgumentException("xArray must have at least 3 items");
        }

        /*
         * Q is a stack of candidate points which have been pushed once onto the stack
         * and removed if they are not vertices of the stack.
         *
         * when complete, the stack S contains the vertices of the hull in counterclockwise order.
         *
         * Q > 3
         *
         * 1 -- let p0 be the point in Q w/ min y-coordinate, or leftmost point of a tie
         * 2 -- let <p1, p2, ... pm> be the remaining points in Q.
         *      sorted by polar angle in counter clockwise order around p0.
         *      ** if more than one point has the same angle, remove all but the one that is furthest from p0. **
         * 3 -- push p0 onto S
         * 4 -- push p1 onto S
         * 5 -- push p2 onto S
         * 6 -- for i=3 to m
         * 7 --     do while the angle formed by points NEXT-TO-TOP(S), TOP(S), and p_i makes a nonleft turn
         * 8 --         pop(S)
         * 9 --     push(S)
         * 10 -return S
         */

        // (1) let p0 be the point in Q w/ minimum yCoord,
        //     or the leftmost point if more than one w/ same minimum yCoord.
        MultiArrayMergeSort.sortBy1stArgThen2nd(y, x);
        
        points = new XYStack(x.length);

        int p0Index = 0;

        // (2) let <p1, p2, ..., pm> be the remaining points in Q, sorted
	    //     by polar angle in counterclockwise order around p0
	    //     (if more than one pt has same angle, keep only the furthest from p0)
    	int nPointsUsable = PolarAngleMergeSort.sort(x[p0Index], y[p0Index], x, y);

        if (nPointsUsable < 3) {
	        throw new GrahamScanTooFewPointsException("polar angle sorting has reduced the number of points to less than 3");
        }

        points.push( x[p0Index], y[p0Index]);
        points.push( x[1], y[1]);
        //points.push( x[2], y[2]);
        
        float topX, topY;

        // for i = 3 to m
        //    while angle between next-to-top(S), top(S) and p_i makes a nonleft turn
        //        do pop(S)
        //    push(pi, S)
        for (int i = 2; i < nPointsUsable; i++) {
    
            topX = points.peekTopX();
            topY = points.peekTopY();
            
            double direction = LinesAndAngles.direction(
                points.peekNextToTopX(), points.peekNextToTopY(), topX, topY,
                x[i], y[i]);
            
            while (!points.isEmpty() /*&& !Double.isNaN(direction)*/ && (direction > 0)) {

                points.pop();

                topX = points.peekTopX();
                topY = points.peekTopY();
                
                direction = LinesAndAngles.direction(
                    points.peekNextToTopX(), points.peekNextToTopY(), topX, topY,
                    x[i], y[i]);                
            }

            points.push(x[i], y[i]);
        }

        populateHull();
    }

    public float[] getXHull() {
        return this.xHull;
    }

    public float[] getYHull() {
        return this.yHull;
    }

    protected void populateHull() throws GrahamScanTooFewPointsException {

        if (points == null) {
            throw new GrahamScanTooFewPointsException("Points cannot be null.  Use computeHull first.");
        }
        
        int n = points.getNPoints() + 1;

        this.xHull = Arrays.copyOf(points.x, n);
        this.yHull = Arrays.copyOf(points.y, n);

        this.xHull[n-1] = points.x[0];
        this.yHull[n-1] = points.y[0];
    }

}
