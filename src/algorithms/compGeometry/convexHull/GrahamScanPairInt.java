package algorithms.compGeometry.convexHull;

import algorithms.compGeometry.LinesAndAngles;
import algorithms.QuickSort;
import algorithms.util.PairInt;
import algorithms.util.Stack;
import java.util.ArrayList;
import java.util.List;

/**
 * adapted from 
 * https://code.google.com/p/two-point-correlation/source/browse/src/test/java/algorithms/compGeometry/convexHull/
 * under MIT License (MIT), Nichole King 2013
 * 
  <pre>
  Solves the Convex Hull problem w/ a stack S of candidate points.
 
  Given a set of Q points returns the vertices of the ConvexHull(Q) in clockwise
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
public class GrahamScanPairInt<T extends PairInt> {

    protected Stack<T> points = null;

    protected List<T> hull = null;
    
	public GrahamScanPairInt() {        
	}

    /**
     * find the convex hull of the given (x, y) points.  Note that the resulting
     * hull points have the same first point as last point.
     * 
     * @param input
     * @throws GrahamScanTooFewPointsException 
     */
    @SuppressWarnings({"unchecked"})
    public void computeHull(T[] input) throws GrahamScanTooFewPointsException {

        if (input == null) {
	    	throw new IllegalArgumentException("input cannot be null");
        }
	    if (input.length < 3) {
	        throw new IllegalArgumentException("input must have at least 3 items");
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
        QuickSort.sortByYThenX(input);
        int p0Index = 0;

        // (2) let <p1, p2, ..., pm> be the remaining points in Q, sorted
	    //     by polar angle in counterclockwise order around p0
	    //     (if more than one pt has same angle, keep only the furthest from p0)
    	int nPointsUsable = PolarAngleQuickSort.sort(input[p0Index], input);

        if (nPointsUsable < 3) {
	        throw new GrahamScanTooFewPointsException(
            "polar angle sorting has reduced the number of points to less than 3");
        }
        
        points = new Stack<T>();
        
        points.push((T)(input[p0Index].copy()));
        points.push((T)input[1].copy());
        points.push((T)input[2].copy());
        
        // for i = 3 to m
        //    while angle between next-to-top(S), top(S) and p_i makes a nonleft turn
        //        do pop(S)
        //    push(pi, S)
        for (int i = 3; i < nPointsUsable; i++) {
            
            T top = points.peek();
            T nextToTop = points.peekPopNext();

            double direction = LinesAndAngles.direction(nextToTop, top, input[i]);
            
            //double direction = LinesAndAngles.direction(
            //    nextToTopX, nextToTopY, topX, topY, xi, yi);
            
            while (direction <= 0) {

                points.pop();
                
                if (points.size() < 2) {
                    break;
                }

                top = points.peek();
                nextToTop = points.peekPopNext();
                
                direction = LinesAndAngles.direction(nextToTop, top, input[i]);              
            }

            points.push((T)input[i].copy());
        }

        populateHull();
    }

    public List<T> getHull() {
        return this.hull;
    }

    protected void populateHull() throws GrahamScanTooFewPointsException {

        if (points == null) {
            throw new GrahamScanTooFewPointsException(
            "Points cannot be null.  Use computeHull first.");
        }
        
        int n = points.size() + 1;

        this.hull = new ArrayList<T>();
                
        for (int i = 0; i < (n - 1); ++i) {
            hull.add(points.pop());
        }
        
        this.hull.add((T)hull.get(0).copy());
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        if (hull != null) {
            for (int i = 0; i < hull.size(); ++i) {
                sb.append(hull.get(i));
            }
        }
        return sb.toString();
    }

}
