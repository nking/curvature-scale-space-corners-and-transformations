package algorithms.compGeometry.convexHull;

import algorithms.sort.MiscSorter;
import algorithms.util.FormatArray;
import algorithms.util.Stack;
import java.util.Arrays;

/**
 * adapted from 
 * https://code.google.com/p/two-point-correlation/source/browse/src/test/java/algorithms/compGeometry/convexHull/
 * under MIT License (MIT), Nichole King 2013
 * 
  <pre>
  Solves the Convex Hull problem w/ a stack S of candidate points.
 
  Given a set of Q points returns the vertices of the ConvexHull(Q) in counter-clockwise
  order.   a convex hull is the smallest convex polygon that will include all points in Q.
 
  This Graham's Scan runs in O(n) due to use of a linear runtime polar angle sorter
  to reduce the allowed angles to integer degrees.
 
  It uses a technique called 'rotational sweep' to process vertices in the order
  of the polar angles they form with a reference vertex.
 
  constructed from pseudo-code in Cormen et al. "Introduction to Algorithms
 </pre>
 
 * @author nichole
 */
public class GrahamScanLong {

    public static class CH {
        private long[] xH;
        private long[] yH;
        public CH(long[] x, long[] y) {
            this.xH = x;
            this.yH = y;
        }

        /**
         * @return the xH
         */
        public long[] getXH() {
            return xH;
        }

        /**
         * @return the yH
         */
        public long[] getYH() {
            return yH;
        }
        
        @Override
        public String toString() {
            if (xH == null) {
                return "[]";
            }
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < xH.length; ++i) {
                if (i > 0) {
                    sb.append(", ");
                }
                sb.append("(").append(xH[i]).append(",").append(yH[i]).append(")");
            }
            return sb.toString();
        }
    }
    
    /**
     * find the convex hull of the given (x, y) points.  Note that the resulting
     * hull points have the same first point as last point.
     * runtime complexity is O(N) because counting sort is used in the polar angle
     * sort (which is O(max(N, 360)), removing the O(N*log_2(N)) component.
     * 
     * @param x
     * @param y
     * @return convex hull of points (x, y).  note that the last point equals the first point.
     * @throws GrahamScanTooFewPointsException 
     */
    public static CH computeHull(long[] x, long[] y) throws GrahamScanTooFewPointsException {

        if (x == null) {
            throw new IllegalArgumentException("x cannot be null");
        }
        if (y == null) {
            throw new IllegalArgumentException("y cannot be null");
        }
        if (x.length != y.length) {
            throw new IllegalArgumentException("x must have the same number of items as y");
        }
        if (x.length < 3) {
            throw new IllegalArgumentException("x must have at least 3 items");
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
         * 9 --     push(pi, S)
         * 10 -return S
         */

        // (1) let p0 be the point in Q w/ minimum yCoord,
        //     or the leftmost point if more than one w/ same minimum yCoord.
        int p0Index = findIndexOfMinY(x, y);
        
        if (p0Index != 0) {
            // move the point at index [iP0] to index [0] and move the rest of the arrays as needed
            rewriteToPlaceAt0(p0Index, x, y);
            
            p0Index = 0;
        }
                
        // (2) let <p1, p2, ..., pm> be the remaining points in Q, sorted
	    //     by polar angle in counterclockwise order around p0
	    //     (if more than one pt has same angle, keep only the furthest from p0)
                
        // this step uses angles rounded to degrees between 0 and 360.
        // the runtime complexity is O( max(x.length, 360) )
    	int nPointsUsable = PolarAngleQuickSort.sortCCWBy1stPoint(x, y);

        if (nPointsUsable < 3) {
	     throw new GrahamScanTooFewPointsException("polar angle sorting has reduced the number of points to less than 3");
        }
        
        Stack<Integer> stack = new Stack<Integer>();
        
        stack.push(p0Index);
        stack.push(1);
        stack.push(2);
                
        // for i = 3 to m
        //    while angle between next-to-top(S), top(S) and p_i makes a nonleft turn
        //        do pop(S)
        //    push(pi, S)
        int top; // p1
        int nextToTop; // p2
        long dir;
        int i;
        for (i = 3; i < nPointsUsable; i++) {
            top = stack.peek();
            nextToTop = stack.peekPopNext();
            //do while the angle formed by points NEXT-TO-TOP(S), TOP(S), and p_i makes a nonleft turn
            //         pop(S)
            //push(pi, S)
            
            // takes a left turn (counter-clockwise) when dir is < 0
            // dir = (x1 - x0)(y2 - y0) - (x2 - x0)(y1 - y0) 
            // (x:top - i)(y:nextToTop - i) - (x:nextToTop - i)(y:top - i)
            dir = (x[top] - x[i]) * (y[nextToTop] - y[i]) - (x[nextToTop] - x[i]) * (y[top] - y[i]);
                        
            while (dir >= 0) {
                
                stack.pop();
                                
                if (stack.size() < 2) {
                    // cannot peak at 2 in stack
                    break;
                }

                top = stack.peek();
                nextToTop = stack.peekPopNext();
                
                dir = (x[top] - x[i]) * (y[nextToTop] - y[i]) - (x[nextToTop] - x[i]) * (y[top] - y[i]);                
            }
            
            stack.push(i);
        }
        
        int nH = stack.size();
                
        i = nH - 1;
        long[] xH = new long[nH + 1];
        long[] yH = new long[nH + 1];
        while (!stack.isEmpty()) {
            top = stack.pop();
            xH[i] = x[top];   
            yH[i] = y[top];
            --i;
        }
        xH[xH.length - 1] = x[0];
        yH[xH.length - 1] = y[0];
        CH ch = new CH(xH, yH);
        return ch;
    }

    protected static int findIndexOfMinY(long[] x, long[] y) {
        if (x.length < 1) {
            throw new IllegalArgumentException("x and y must be longer than 1");
        }
        int iMin = 0;
        long minY = y[0];
        long minX = x[0];
        int i;
        for (i = 1; i < x.length; ++i) {
            if ((y[i] < minY) || ( (y[i] == minY) && (x[i] < minX))) {
                minY = y[i];
                minX = x[i];
                iMin = i;
            }
        }
        return iMin;
    }
    
    /**
     * move the point at index [iP0] to index [0] and move the rest of the array as needed
     * @param iP0
     * @param x
     * @param y 
     */
    protected static void rewriteToPlaceAt0(final int iP0, long[] x, long[] y) {
        long x0 = x[iP0];
        long y0 = y[iP0];
        
        /*
        0---- modify
        1---- modify
        2 iP0
        3 -----  not affected by iP0 move
        4 -----
        */
        int i;
        for (i = iP0; i > 0; --i) {
            x[i] = x[i - 1];
            y[i] = y[i - 1];
        }
        x[0] = x0;
        y[0] = y0;
    }

}
