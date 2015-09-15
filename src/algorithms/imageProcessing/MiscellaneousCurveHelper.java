package algorithms.imageProcessing;

import algorithms.MultiArrayMergeSort;
import algorithms.misc.AverageUtil;
import algorithms.misc.MiscDebug;
import algorithms.misc.MiscMath;
import algorithms.util.PairIntArray;
import algorithms.util.PairIntArrayComparator;
import algorithms.util.PairFloatArray;
import algorithms.util.PairInt;
import algorithms.util.PairIntArrayWithColor;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.Stack;
import java.util.logging.Logger;
import org.ejml.simple.*;

/**
 *
 * @author nichole
 */
public class MiscellaneousCurveHelper {

    private Logger log = Logger.getLogger(this.getClass().getName());

    // choosing a minimum size empirically from looking at edges in tests
    private static int minLedgeWidth = 4;

    protected static final int[] eightNeighborsX =
        new int[]{-1, -1,  0,  1, 1, 1, 0, -1};
    protected static final int[] eightNeighborsY =
        new int[]{ 0, -1, -1, -1, 0, 1, 1,  1};

    /**
     * determine whether the closed curve points are ordered in a clockwise
     * manner.
     * calculating the cross product between adjacent edges in sequence around
     * the polygon to determine if there are fewer that are positive (CCW)
     * than negative (CW).
     * The given closedCurve cannot intersect itself or have holes in it.
     * NOTE: the answer returns true if the points are ordered in CW manner, but
     * if one needs the answer w.r.t. viewing an image which has y increasing
     * downward, need the opposite of the return here.
     *
     * @param closedCurve
     * @return
     */
    public boolean curveIsOrderedClockwise(PairIntArray closedCurve) {

        if (closedCurve.getN() < 2) {
            return false;
        }

        int nNeg = 0;
        int n = closedCurve.getN();

        for (int i = 0; i < n; i++) {

            long xm1, ym1, x, y, xp1, yp1;

            if (i == 0) {
                xm1 = closedCurve.getX(closedCurve.getN() - 1);
                ym1 = closedCurve.getY(closedCurve.getN() - 1);
                xp1 = closedCurve.getX(i + 1);
                yp1 = closedCurve.getY(i + 1);
            } else if (i == (closedCurve.getN() - 1)) {
                xm1 = closedCurve.getX(i - 1);
                ym1 = closedCurve.getY(i - 1);
                xp1 = closedCurve.getX(0);
                yp1 = closedCurve.getY(0);
            } else {
                xm1 = closedCurve.getX(i - 1);
                ym1 = closedCurve.getY(i - 1);
                xp1 = closedCurve.getX(i + 1);
                yp1 = closedCurve.getY(i + 1);
            }
            x = closedCurve.getX(i);
            y = closedCurve.getY(i);

            long dxmxm1 = (x - xm1);
            long dymym1 = (y - ym1);
            long dxp1mx = (xp1 - x);
            long dyp1my = (yp1 - y);

            //(xi - xi-1) * (yi+1 - yi) - (yi - yi-1) * (xi+1 - xi)
            long crossProduct = (dxmxm1 * dyp1my) - (dymym1 * dxp1mx);

            if (crossProduct < 0) {
                // clockwise when crossProduct is negative
                nNeg++;
            }
        }

        int nPos = n - nNeg;//n - 2 - nNeg;

        //log.info(closedCurve.toString());
        //log.info("n=" + n + " nNegative=" + nNeg + " nPositive=" + nPos);

        return ((n > 2) && (nNeg >= nPos)) || (nNeg > nPos);
     }

    public void additionalThinning45DegreeEdges(
        GreyscaleImage theta, GreyscaleImage input) {

        // thin the edges for angles 45 and -45 as suggested by
        // 1998 Mokhtarian and Suomela
        // IEEE TRANSACTIONS ON PATTERN ANALYSIS AND MACHINE INTELLIGENCE,
        //     VOL. 20, NO. 12
        //
        //compare each edge pixel which has an edge orientation of
        // 45o or -45o to one of its horizontal or vertical neighbors.
        // If the neighbor has the same orientation, the other point can be
        // removed.
        for (int i = 1; i < (input.getWidth() - 1); i++) {
            for (int j = 1; j < (input.getHeight() - 1); j++) {

                int tG = theta.getValue(i, j);

                if (((tG == 45) || (tG == -45)) && (input.getValue(i, j) > 0)) {

                    int tH0 = theta.getValue(i - 1, j);
                    int tH1 = theta.getValue(i + 1, j);
                    int tV0 = theta.getValue(i, j - 1);
                    int tV1 = theta.getValue(i, j + 1);

                    int gH0 = input.getValue(i - 1, j);
                    int gH1 = input.getValue(i + 1, j);
                    int gV0 = input.getValue(i, j - 1);
                    int gV1 = input.getValue(i, j + 1);

                    if ((tH0 == tG) && (gH0 > 0)) {
                        if ((tV0 == tG) && (gV0 > 0)) {
                            input.setValue(i, j, 0);
                        } else if ((tV1 == tG) && (gV1 > 0)) {
                            input.setValue(i, j, 0);
                        }
                    } else if ((tH1 == tG) && (gH1 > 0)) {
                        if ((tV0 == tG) && (gV0 > 0)) {
                            input.setValue(i, j, 0);
                        } else if ((tV1 == tG) && (gV1 > 0)) {
                            input.setValue(i, j, 0);
                        }
                    }
                }
            }
        }
    }

     /**
     * this is a method to combine and prune adjacent lines, but
     * it knows nothing about the overall shape. it chooses to keep the longer
     * line and append any dangling members of the adjacent line to the longest.
     * It's a time consuming method (add runtime complexity here).  It isn't
     * used anymore because the results from the CannyEdgeFilter are now
     * 1 pixel thick edges already.
     * @param edges
     * @param imageWidth the image width in pixels
     * @return
     */
    public List<PairIntArray> pruneAndIncludeAdjacentCurves(
        List<PairIntArray> edges, int imageWidth) {

        //sort to place edges with fewest points at top
        Collections.sort(edges, new PairIntArrayComparator());

        Stack<PairIntArrayWithColor> stack = new Stack<PairIntArrayWithColor>();
        List<PairIntArrayWithColor> pruneThese = new ArrayList<PairIntArrayWithColor>();
        for (int i = 0; i < edges.size(); i++) {
            // color: 0 = unvisited, 1 = processing, 2 = visited,
            //        3 = in an output edge, 4=pruned
            PairIntArrayWithColor node = new PairIntArrayWithColor(edges.get(i));
            stack.add(node);
            pruneThese.add(node);
        }

        List<PairIntArray> output = new ArrayList<PairIntArray>();

        PairIntArrayWithColor tmp = stack.peek();
        tmp.setColor(2);

        boolean foundOverlappingCurves = false;
        boolean reversedPoints = false;

        while (!stack.isEmpty()) {

            PairIntArrayWithColor uNode = stack.pop();

            if (uNode.getColor() == 4) {
                continue;
            }

            // for each neighbor v of u
            for (int i = (pruneThese.size() - 1); i > -1; i--) {

                PairIntArrayWithColor vNode = pruneThese.get(i);

                if (vNode.getColor() != 0) {
                    continue;
                }

                // because top item might be updated in processPair,
                //     place in v iter
                int uX = uNode.getX(0);
                int uY = uNode.getY(0);
                int uIdx = (uY * imageWidth) + uX;

                int vX = vNode.getX(0);
                int vY = vNode.getY(0);
                int vIdx = (vY * imageWidth) + vX;
                if (uIdx == vIdx) {
                    continue;
                }

                boolean areOverlapping = processOverlappingPair(uNode, vNode);

                if (areOverlapping) {

                    foundOverlappingCurves = true;

                    // color: 0 = unvisited, 1 = processing, 2 = visited,
                    //        3 = in an output edge, 4=pruned
                    vNode.setColor(4);

                    pruneThese.remove(vNode);

                } // else insert into stack?
            }

            if (!foundOverlappingCurves && !reversedPoints) {
                uNode.reverse();
                reversedPoints = true;
                log.fine("reversed edge to try again");
                stack.add(uNode);
                continue;
            }

            // color: 0 = unvisited, 1 = processing, 2 = visited,
            //        3 = in an output edge, 4=pruned
            uNode.setColor(3);

            output.add(uNode);

            foundOverlappingCurves = false;
            reversedPoints = false;
        }

        return output;
    }

    /**
     * given 2 edges, return true if they overlap. If they overlap
     * curve0 is given the larger curve and any outlying points.
     * @param curve0
     * @param curve1
     * @return
     */
    protected boolean processOverlappingPair(PairIntArrayWithColor curve0,
        PairIntArrayWithColor curve1) {

        boolean longerIsNode0 = (curve0.getN() >= curve1.getN());

        PairIntArray longer, shorter;
        if (longerIsNode0) {
            longer = curve0;
            shorter = curve1;
        } else {
            longer = curve1;
            shorter = curve0;
        }

        // used to return the offset w.r.t. the longest edge.
        int[] correlationOffset = new int[1];

        /*
         returns whether the curve 'check' is adjacent to the curve 'node0',
         and if so, returns the offset in the frame of the larger curve.
         the offset represents where the first point in the shorter curve
         matches in the larger curve.
         */
        boolean isAdjacent = correlation(longer, shorter,
            correlationOffset);

        if (isAdjacent) {

            /*
                -- find any points in check outside of the overlap
                   and add those to the larger node.
            */
            if (correlationOffset[0] < 0) {
                // add from the beginning of shorter if any are unmatched
                int nInsert = -1*correlationOffset[0];
                longer.insertSpaceAtTopOfArrays(nInsert);

                for (int ii = 0; ii < nInsert; ii++) {
                    longer.set(ii, shorter.getX(ii), shorter.getY(ii));
                }
            } else {
                //add from end of shorter if any are unmatched
                int n0 = longer.getN() - correlationOffset[0];
                if (n0 < shorter.getN()) {
                    for (int ii = n0; ii < shorter.getN(); ii++) {
                        longer.add(shorter.getX(ii), shorter.getY(ii));
                    }
                }
            }

            if (!longerIsNode0) {
                curve0.swapContents(curve1);
            }

            return true;
        }

        return false;
    }

    /**
     * return true if correlation shows that the 2 curves are adjacent
     * to one another.  Note that the method needs the points within the
     * curves to be ordered in a similar manner and for the endpoints of the
     * curves to be accurate.  If a point in the middle of the curve is
     * the first or last point, it may prevent comparison of it with another
     * edge's endpoints.
     *
     * @param curve0
     * @param curve1
     * @param correlationOffset offset of where the shorter curve starts
     *  with respect to the longer.  For example, an offset of -2 means that
     * the first 2 points in the shorter curve are outside of the longer curve,
     * but the next point in the longer curve is adjacent to the shorter.
     * Another example: if offset is +2, the first pixel in the shorter curve
     * is adjacent to the third pixel in the longer curve.  NOTE: the offset
     * is only useful if this method returns true;
     * @return
     */
    protected boolean correlation(PairIntArray curve0, PairIntArray curve1,
        int[] correlationOffset) {

        correlationOffset[0] = Integer.MAX_VALUE;

        //TODO: look at string matching algorithms to explore improvements here

        PairIntArray shorter, longer;
        if (curve0.getN() <= curve1.getN()) {
            shorter = curve0;
            longer = curve1;
        } else {
            shorter = curve1;
            longer = curve0;
        }

        /*
        len0 = 5; len1 = 11;
         #####
             +++++++++++
          #####
             +++++++++++
           #####
             +++++++++++
            #####
             +++++++++++
             #####
             +++++++++++

             #####
             +++++++++++

                       #####
             +++++++++++

        ccs = sqrt(sumsqdiff)/nOverlapping if nOverlapping > 0.

        if (css <= 1 pix * nOverlapping) {
            store as a possible adjacent curve
        }
        compare possible adjacent curves for the smallest css, and store that
        offset in correlationOffset and return true, else false
        */

        double cSSMin = Double.MAX_VALUE;
        int cSSMinOffset = Integer.MAX_VALUE;
        int cSSMinNOverlapping = 0;

        double sqrtTwo = Math.sqrt(2);

        for (int i = 0; i < (longer.getN() + shorter.getN() - 1); i++) {
            //siIdx is first index in shorter for comparison
            //sfIdx is last index in shorter for comparison
            //liIdx is first index of longer for comparison
            //lfIdx is last index of longer for comparison
            int siIdx, sfIdx, liIdx, lfIdx, offset;
            if (i < shorter.getN()) {
                /*
                 #####
                     +++++++++++ i=0
                  #####
                     +++++++++++
                   #####
                     +++++++++++
                    #####
                     +++++++++++
                     #####
                     +++++++++++ i=4
                */
                //sfIdx is inclusive endpoint
                sfIdx = shorter.getN() - 1;
                siIdx = sfIdx - i;
                liIdx = 0;
                //lfIdx is inclusive endpoint
                lfIdx = sfIdx - siIdx;
                offset = i - sfIdx;

            } else if (i < longer.getN() ) {

                /*
                      #####
                     +++++++++++  i=5

                       #####
                     +++++++++++

                        #####
                     +++++++++++

                         #####
                     +++++++++++

                          #####
                     +++++++++++

                           #####
                     +++++++++++ i = 10
                */
                //sfIdx is inclusive endpoint
                sfIdx = shorter.getN() - 1;
                siIdx = 0;
                liIdx = i - shorter.getN() + 1;
                //lfIdx is inclusive endpoint
                lfIdx = liIdx + shorter.getN() - 1;
                offset = i - sfIdx;

            } else {
                /*
                            #####
                     +++++++++++ i = 12
                     01234567890
                             #####
                     +++++++++++

                              #####
                     +++++++++++

                               #####
                     +++++++++++ i=15
                     01234567890
                */
                liIdx = i - shorter.getN() + 1;
                //sfIdx is inclusive endpoint
                sfIdx = longer.getN() - liIdx - 1;
                siIdx = 0;
                //lfIdx is inclusive endpoint
                lfIdx = liIdx + (sfIdx - siIdx);
                offset = liIdx;

            }

            int nOverLapping = (sfIdx - siIdx) + 1;

            if ((sfIdx - siIdx) != (lfIdx - liIdx)) {
                throw new IllegalStateException(
                    "sample ranges are not correct");
            }

            double sumSq = 0;

            int s = siIdx;
            int l = liIdx;
            while (s <= sfIdx) {
                int xs = shorter.getX(s);
                int xl = longer.getX(l);
                int dx = xs - xl;
                int ys = shorter.getY(s);
                int yl = longer.getY(l);
                int dy = ys - yl;
                sumSq += ((dx*dx) + (dy*dy));
                s++;
                l++;
            }

            double tmp = Math.sqrt(sumSq/nOverLapping);

            // assuming adjacent pixel has distance of sqrt(2) at the most
            if (tmp <= sqrtTwo) {

                if ((tmp < cSSMin) ||
                (tmp == cSSMin && (nOverLapping > cSSMinNOverlapping))
                ) {

                    cSSMin = tmp;

                    cSSMinOffset = offset;

                    cSSMinNOverlapping = nOverLapping;
                }
            }
        }

        if (cSSMin < Double.MAX_VALUE) {

            correlationOffset[0] = cSSMinOffset;

            return true;
        }

        return false;
    }

    /**
     * find the index where x is minimum value of closedCurve.  Note that when
     * there are more than one points with the same minimum x value, the point
     * with a smaller y is chosen.
     *
     * @param closedCurve
     * @return
     */
    public int findMinIdx(PairIntArray closedCurve) {

        int xMin = Integer.MAX_VALUE;
        int xMax = Integer.MIN_VALUE;

        int xMinIdx = -1;
        int xMaxIdx = -1;

        // find xMin.  when xMin==x, use yMin too.  similar pattern for maxes
        for (int i = 0; i < closedCurve.getN(); i++) {
            if (closedCurve.getX(i) < xMin) {
                xMin = closedCurve.getX(i);
                xMinIdx = i;
            } else if (closedCurve.getX(i) == xMin) {
                if (closedCurve.getY(i) < closedCurve.getY(xMinIdx)) {
                    xMin = closedCurve.getX(i);
                    xMinIdx = i;
                }
            }
            if (closedCurve.getX(i) > xMax) {
                xMax = closedCurve.getX(i);
                xMaxIdx = i;
            } else if (closedCurve.getX(i) == xMax) {
                if (closedCurve.getY(i) > closedCurve.getY(xMaxIdx)) {
                    xMax = closedCurve.getX(i);
                    xMaxIdx = i;
                }
            }
        }

        return xMinIdx;
    }

    public double[] calculateXYCentroids(PairIntArray xy, float[] weights) {

        double xc = 0;
        double yc = 0;

        for (int i = 0; i < xy.getN(); i++) {
            double x1 = xy.getX(i);
            xc += (weights[i] * x1);

            double y1 = xy.getY(i);
            yc += (weights[i] * y1);
        }

        return new double[]{xc, yc};
    }

    public double[] calculateXYCentroids(PairIntArray xy) {

        double xc = 0;
        double yc = 0;

        for (int i = 0; i < xy.getN(); i++) {

            xc += xy.getX(i);

            yc += xy.getY(i);
        }

        xc /= (double)xy.getN();

        yc /= (double)xy.getN();

        return new double[]{xc, yc};
    }

    /**
     * calculate the x and y centroids and return as
     * double[]{xCentroid, yCentroid}
     * @param xy a 3 x N matrix with column 0 being x and column 1 being y.
     * @return
     */
    public double[] calculateXYCentroids(SimpleMatrix xy) {

        double xc = 0;
        double yc = 0;

        int n = xy.numCols();

        for (int i = 0; i < n; i++) {

            xc += xy.get(0, i);

            yc += xy.get(1, i);
        }

        xc /= (double)n;

        yc /= (double)n;

        return new double[]{xc, yc};
    }

    public double[] calculateXYCentroids(PairFloatArray xy) {

        double xc = 0;
        double yc = 0;

        for (int i = 0; i < xy.getN(); i++) {

            xc += xy.getX(i);

            yc += xy.getY(i);
        }

        xc /= (double)xy.getN();

        yc /= (double)xy.getN();

        return new double[]{xc, yc};
    }

    public double[] calculateXYCentroids(List<PairIntArray> xyList) {

        double xc = 0;
        double yc = 0;

        for (PairIntArray points : xyList) {

            double[] xycen = calculateXYCentroids(points);

            xc += xycen[0];
            yc += xycen[1];
        }

        xc /= (double)xyList.size();
        yc /= (double)xyList.size();

        return new double[]{xc, yc};
    }

    public double[] calculateXYCentroids(Set<PairInt> points) {

        double xc = 0;
        double yc = 0;

        for (PairInt p : points) {

           int x = p.getX();
           int y = p.getY();

            xc += x;
            yc += y;
        }

        xc /= (double)(points.size());

        yc /= (double)(points.size());

        return new double[]{xc, yc};
    }

    public double[] calculateXYCentroids(float[] x, float[] y) {

        if (x == null) {
            throw new IllegalArgumentException("x cannot be null");
        }
        if (y == null) {
            throw new IllegalArgumentException("y cannot be null");
        }
        if (x.length != y.length) {
            throw new IllegalArgumentException("x and y must be same length");
        }

        double xc = 0;
        double yc = 0;

        for (int i = 0; i < x.length; i++) {

            xc += x[i];

            yc += y[i];
        }

        xc /= (double)(x.length);

        yc /= (double)(x.length);

        return new double[]{xc, yc};
    }

    /**
     * search for point in edge with value (x, y) within indexes lowIdx to
     * highIdx, inclusive and return true if found, else false.
     * Bounds checking is done internally, so it's safe to pass lowIdx
     * or highIdx out of range of edge.
     * @param x
     * @param i
     * @param edge
     * @param lowIdx
     * @param highIdx
     * @return
     */
    private boolean pointExists(int x, int y, PairIntArray edge, int lowIdx,
        int highIdx) {

        for (int i = lowIdx; i <= highIdx; i++) {
            if ((i < 0) || (i > (edge.getN() - 1))) {
                continue;
            }
            if ((edge.getX(i) == x) && (edge.getY(i) == y)) {
                return true;
            }
        }

        return false;
    }

    public void removeRedundantPoints(List<PairIntArray> tmpEdges) {

        //TODO: fix this method!!

        log.fine("removeRedundantPoints");

        // if there are redundant points, remove the points in between

        //TODO: replace w/ faster algorithm...

        for (int i = 0; i < tmpEdges.size(); i++) {

            List<String> points = new ArrayList<String>();

            PairIntArray edge = tmpEdges.get(i);

            removeRedundantPoints(edge);
        }
    }

    public void removeRedundantPoints(PairIntArray edge) {

        log.fine("removeRedundantPoints");

        // if there are redundant points, remove the points in between

        //TODO: replace w/ faster algorithm...

        List<String> points = new ArrayList<String>();

        for (int j = (edge.getN() - 1); j > -1; j--) {

            String str = String.format("%d:%d", edge.getX(j),
                edge.getY(j));

            int idx = points.indexOf(str);

            if (idx > -1) {

                //TODO: consider a limit for (pIdx - j)

                int pIdx = edge.getN() - idx - 1;

                edge.removeRange(j, pIdx - 1);

                // restart comparison? if we remove same region from points, we don't have to restart
                points.clear();

                j = edge.getN();

            } else {

                points.add(str);
            }
        }
    }

    public void pruneAdjacentNeighborsTo2(List<PairIntArray> tmpEdges) {

        log.fine("pruneAdjacentNeighborsTo2");

        for (int lIdx = 0; lIdx < tmpEdges.size(); lIdx++) {

            // quick check for whether an edged has 3 neighbors, then
            // compare with patterns

            PairIntArray edge = tmpEdges.get(lIdx);

            pruneAdjacentNeighborsTo2(edge);
        }
    }

    public void pruneAdjacentNeighborsTo2(PairIntArray edge) {

        log.fine("pruneAdjacentNeighborsTo2");

        // this will usually only have 2 in it, and most expected is 3
        int[] outputAdjacentNeighbors = new int[8];

        // quick check for whether an edged has 3 neighbors, then
        // compare with patterns

        for (int eIdx = 0; eIdx < edge.getN(); eIdx++) {

            int nNeighbors = getAdjacentNeighbors(edge, eIdx,
                outputAdjacentNeighbors);

            if ((nNeighbors > 2)) {

                boolean pruned = pruneAdjacentNeighborsTo2(edge, eIdx,
                    outputAdjacentNeighbors, nNeighbors);

                if (pruned) {
                    // restart iteration for easier maintenance
                    eIdx = -1;
                }
            }
        }
    }

    /**
     * does removing the point at idx create a gap between it's neighboring
     * pixels?  this uses the simplest test of only the points at idx-1
     * and idx+1.
     *
     * @param edge
     * @param idx
     * @return
     */
    public boolean doesDisconnect(PairIntArray edge, int idx) {

        // test for endpoints first
        if (idx == 0) {

            if (edge.getN() < 3) {
                return true;
            }

            // does this point currently connect to the last point?
            float diffX = edge.getX(idx) - edge.getX(edge.getN() - 1);
            if (diffX < 0) {
                diffX *= -1;
            }
            float diffY = edge.getY(idx) - edge.getY(edge.getN() - 1);
            if (diffY < 0) {
                diffY *= -1;
            }
            if (((diffX < 2) && (diffY < 2))) {
                // this is connected to the last point in the edge
                // check to see if lastPoint and idx + 1 are adjacent
                diffX = edge.getX(idx + 1) - edge.getX(edge.getN() - 1);
                if (diffX < 0) {
                    diffX *= -1;
                }
                if (diffX > 1) {
                    return true;
                }

                diffY = edge.getY(idx + 1) - edge.getY(edge.getN() - 1);
                if (diffY < 0) {
                    diffY *= -1;
                }
                if (diffY > 1) {
                    return true;
                }
            }
            return false;
        }

        if (idx == (edge.getN() - 1)) {

            if (edge.getN() < 3) {
                return true;
            }

            // does this point currently connect to the first point?
            float diffX = edge.getX(idx) - edge.getX(0);
            if (diffX < 0) {
                diffX *= -1;
            }
            float diffY = edge.getY(idx) - edge.getY(0);
            if (diffY < 0) {
                diffY *= -1;
            }
            if (((diffX < 2) && (diffY < 2))) {
                // this is connected to the first point in the edge
                // check to see if lastPoint - 1 and first point are adjacent
                diffX = edge.getX(idx - 1) - edge.getX(0);
                if (diffX < 0) {
                    diffX *= -1;
                }
                if (diffX > 1) {
                    return true;
                }

                diffY = edge.getY(idx - 1) - edge.getY(0);
                if (diffY < 0) {
                    diffY *= -1;
                }
                if (diffY > 1) {
                    return true;
                }
            }
            return false;
        }

        if ((idx + 1) < edge.getN()) {
            float diffX = edge.getX(idx - 1) - edge.getX(idx + 1);
            if (diffX < 0) {
                diffX *= -1;
            }
            if (diffX > 1) {
                return true;
            }

            float diffY = edge.getY(idx - 1) - edge.getY(idx + 1);
            if (diffY < 0) {
                diffY *= -1;
            }
            if (diffY > 1) {
                return true;
            }

            return false;
        }

        return false;
    }

    public boolean pruneAdjacentNeighborsTo2(PairIntArray edge, final int eIdx,
        int[] outputAdjacentNeighbors, int nOutputAdjacentNeighbors) {

        int h = 5;
        if (((eIdx - h) < 0) || ((eIdx + h) > (edge.getN() - 1))) {
            return false;
        }

        final int x = edge.getX(eIdx);
        final int y = edge.getY(eIdx);

        // find which point among outputAdjacentNeighbors and eIdx is furthest
        // tangentially from a line formed by the neighboring 5 points on each
        // side of (x, y) and remove that point from the edge

        int x0 = edge.getX(eIdx - h);
        int y0 = edge.getY(eIdx - h);
        int x1 = edge.getX(eIdx + h);
        int y1 = edge.getY(eIdx + h);

        // which of the 3 or more in outputAdjacentNeighbors do not disconnect
        //   the adjacent lines?
        double maxDistance = Double.MIN_VALUE;
        int maxDistanceIdx = -1;

        // if removing this point at eIdx would not disconnect the surrounding
        // points, initialize maxDistance and maxDistanceIdx with it
        if (!doesDisconnect(edge, eIdx)) {
            maxDistance = distanceFromPointToALine(x0, y0, x1, y1, x, y);
            maxDistanceIdx = eIdx;
        }

        for (int i = 0; i < nOutputAdjacentNeighbors; i++) {

            int idx = outputAdjacentNeighbors[i];

            if (!doesDisconnect(edge, idx)) {
                int xCompare = edge.getX(idx);
                int yCompare = edge.getY(idx);

                double dist = distanceFromPointToALine(x0, y0, x1, y1, xCompare,
                    yCompare);

                if (dist > maxDistance) {
                    maxDistance = dist;
                    maxDistanceIdx = idx;
                }
            }
        }

        if (maxDistanceIdx == -1) {
            return false;
        }

        log.finest("removing point (" + edge.getX(maxDistanceIdx)
            + "," + edge.getY(maxDistanceIdx) + ") " +
            "idx=" + maxDistanceIdx + " out of " + edge.getN());

        edge.removeRange(maxDistanceIdx, maxDistanceIdx);

        return true;
    }

    /**
     * looks for the immediate adjacent neighbors and return their indexes
     * in outputAdjacentNeighborIndexes and return for the method the number
     * of items to read in outputAdjacentNeighborIndexes.
     *
     * @param edge
     * @param idx
     * @param outputAdjacentNeighborIndexes
     * @return
     */
    public int getAdjacentNeighbors(PairIntArray edge, int idx,
        int[] outputAdjacentNeighborIndexes) {

        float x = edge.getX(idx);
        float y = edge.getY(idx);

        int nAdjacent = 0;
        for (int i = (idx - 2); i <= (idx + 2); i++) {
            if (i == idx) {
                continue;
            }
            if ((i < 0) || (i > (edge.getN() - 1))) {
                continue;
            }

            float diffX = edge.getX(i) - x;
            float diffY = edge.getY(i) - y;

            if ((Math.abs(diffX) < 2) && (Math.abs(diffY) < 2)) {

                outputAdjacentNeighborIndexes[nAdjacent] = i;

                nAdjacent++;
            }
        }

        return nAdjacent;
    }

    public double distanceFromPointToALine(float lineX0, float lineY0,
        float lineX1, float lineY1, float xP, float yP) {

        /*
        en.wikipedia.org/wiki/Distance_from_a_point_to_a_line

        for the edge, we have the 2 points (lineX0, lineY0) and (lineX1, lineY1)

        distance between that edge and a point (xP, yP) is

        defining diffX = lineX1 - lineX0
                 diffY = lineY1 - lineY0;

        d =
           ( diffY*xP - diffX*yP - lineX0*lineY1 + lineX1*lineY0 )
           ( --------------------------------------------------- )
           (         (diffX*diffX + diffY*diffY)^0.5             )
        )
        */

        float diffX = lineX1 - lineX0;
        float diffY = lineY1 - lineY0;

        if (diffY == 0) {
            // horizontal line
            return Math.abs(yP - lineY0);
        } else if (diffX == 0) {
            // vertical line
            return Math.abs(xP - lineX0);
        }

        double pt1 = Math.abs(diffY*xP - diffX*yP - lineX0*lineY1 + lineX1*lineY0);

        double pt2 = Math.sqrt(diffX*diffX + diffY*diffY);

        double dist = pt1/pt2;

        return dist;
    }

    public double distanceFromPointToALine(int lineX0, int lineY0,
        int lineX1, int lineY1, int xP, int yP) {

        /*
        en.wikipedia.org/wiki/Distance_from_a_point_to_a_line

        for the edge, we have the 2 points (lineX0, lineY0) and (lineX1, lineY1)

        distance between that edge and a point (xP, yP) is

        defining diffX = lineX1 - lineX0
                 diffY = lineY1 - lineY0;

        d =
           ( diffY*xP - diffX*yP - lineX0*lineY1 + lineX1*lineY0 )
           ( --------------------------------------------------- )
           (         (diffX*diffX + diffY*diffY)^0.5             )
        )
        */

        int diffX = lineX1 - lineX0;
        int diffY = lineY1 - lineY0;

        if (diffY == 0) {
            // horizontal line
            return Math.abs(yP - lineY0);
        } else if (diffX == 0) {
            // vertical line
            return Math.abs(xP - lineX0);
        }

        int pt1 = Math.abs(diffY*xP - diffX*yP - lineX0*lineY1 + lineX1*lineY0);

        double pt2 = Math.sqrt(diffX*diffX + diffY*diffY);

        double dist = pt1/pt2;

        return dist;
    }

    public void correctCheckeredSegments(List<PairIntArray> tmpEdges) {

        /*
        there are sometimes sections in the line where one pixel is displaced
               @      @ @ @@@@@@@@@
        @@@@@@@ @@@@@@ @ @

        So far, only seen in horizontal or vertical segments.
        */

        for (int lIdx = 0; lIdx < tmpEdges.size(); lIdx++) {

            PairIntArray edge = tmpEdges.get(lIdx);

            correctCheckeredSegments(edge);
        }
    }

    public void correctCheckeredSegments(PairIntArray edge) {

        /*
        there are sometimes sections in the line where one pixel is displaced
               @      @ @ @@@@@@@@@
        @@@@@@@ @@@@@@ @ @

        So far, only seen in horizontal or vertical segments.
        */

        int[] xs = new int[2];
        int[] ys = new int[2];

        for (int i = 0; i < edge.getN(); i++) {
            correctCheckeredSegments(edge, i, xs, ys);
        }
    }

    public void debugPrint(PairIntArray edge) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < edge.getN(); i++) {
             int x = edge.getX(i);
             int y = edge.getY(i);
             sb.append(String.format("%d)  (%d, %d)\n", i, x, y));
        }
        log.info(sb.toString());
    }

    public int indexOfPointsInRange(List<PairIntArray> edges,
        int xLo, int xHi, int yLo, int yHi) {

        for (int i = 0; i < edges.size(); i++) {

            PairIntArray edge = edges.get(i);

            for (int j = 0; j < edge.getN(); j++) {
                int x = edge.getX(j);
                int y = edge.getY(j);
                if ((x >= xLo) && (x <= xHi) && (y >= yLo) && (y <= yHi)) {
                    return i;
                }
            }
        }

        return -1;
    }

    public int debugFindEdgeContainingPoint(List<PairIntArray> edges, int x, int y) {
        for (int i = 0; i < edges.size(); i++) {
            PairIntArray edge = edges.get(i);
            if (debugEdgeContainsPoint(edge, x, y)) {
                return i;
            }
        }
        return -1;
    }
    public boolean debugEdgeContainsPoint(PairIntArray edge, int x, int y) {
        for (int i = 0; i < edge.getN(); i++) {
            int xc = edge.getX(i);
            if (xc != x) {
                continue;
            }
            int yc = edge.getY(i);
            if (yc != x) {
                continue;
            }
            return true;
        }
        return false;
    }

    private boolean debugIsSection1(PairIntArray edge, int idx) {
        int x = edge.getX(idx);
        int y = edge.getY(idx);
        if ((x > 92) && (x < 99) && (y > 58) && (y < 63)) {
            return true;
        }
        return false;
    }

    private boolean debugIsSection2(PairIntArray edge, int idx) {
        int x = edge.getX(idx);
        int y = edge.getY(idx);
        if ((x > 134) && (x < 147) && (y > 58) && (y < 63)) {
            return true;
        }
        return false;
    }

    private void correctCheckeredSegments(PairIntArray edge, int idx,
        int[] xs, int[] ys) {

        int h = 4;

        /*
        there are sometimes sections in the line where one pixel is displaced
               @      @ @ @@@@@@@@@
        @@@@@@@ @@@@@@ @ @

        So far, only seen in horizontal or vertical segments.
        */

        int nx = 0;
        int ny = 0;

        for (int i = (idx - h); i <= (idx + h); i++) {

            if ((i < 0) || (i > (edge.getN() - 1))) {
                continue;
            }

            int x = edge.getX(i);
            int y = edge.getY(i);

            if (nx == 0) {
                xs[nx] = x;
                nx++;
            } else if (nx == 1) {
                if (xs[0] != x) {
                    xs[nx] = x;
                    nx++;
                }
            } else if ((xs[0] != x) && (xs[1] != x)) {
                // increment to higher than "1" so we know it's not a vertical
                // checkered/fence pattern
                nx++;
            }

            if (ny == 0) {
                ys[ny] = y;
                ny++;
            } else if (ny == 1) {
                if (ys[0] != y) {
                    ys[ny] = y;
                    ny++;
                }
            } else if ((ys[0] != y) && (ys[1] != y)) {
                // increment to higher than "2" so we know it's not a horizontal
                // checkered/fence pattern
                ny++;
            }
        }

        if ((ny == 2) || (nx == 2)) {

            for (int i = (idx - h + 2); i <= (idx + h - 1); i++) {

                //TODO: condense these:
                if ((i < 0) || (i > (edge.getN() - 1))) {
                    continue;
                }
                if (((i - 1) < 0) || ((i - 1) > (edge.getN() - 1))) {
                    continue;
                }
                if (((i - 2) < 0) || ((i - 2) > (edge.getN() - 1))) {
                    continue;
                }
                if (((i + 1) < 0) || ((i + 1) > (edge.getN() - 1))) {
                    continue;
                }

                if (ny == 2) {
                    int prev2Y = edge.getY(i - 2);
                    int prev1Y = edge.getY(i - 1);

                    if (prev2Y != prev1Y) {
                        continue;
                    }

                    int prev2X = edge.getX(i - 2);
                    int prev1X = edge.getX(i - 1);
                    int x = edge.getX(i);
                    int next1X = edge.getX(i + 1);

                    int dx = prev1X - prev2X;

                    if (!(((prev1X + dx) == x) && ((x + dx) == next1X))) {
                        continue;
                    }

                    int y = edge.getY(i);
                    int next1Y = edge.getY(i + 1);
                    if ((next1Y == prev1Y) && (next1Y != y)) {
                        edge.set(i, x, next1Y);
                    }
                } else {
                    int prev2X = edge.getX(i - 2);
                    int prev1X = edge.getX(i - 1);

                    if (prev2X != prev1X) {
                        continue;
                    }

                    int prev2Y = edge.getY(i - 2);
                    int prev1Y = edge.getY(i - 1);
                    int y = edge.getY(i);
                    int next1Y = edge.getY(i + 1);

                    int dy = prev1Y - prev2Y;

                    if (!(((prev1Y + dy) == y) && ((y + dy) == next1Y))) {
                        continue;
                    }

                    int x = edge.getX(i);
                    int next1X = edge.getX(i + 1);
                    if ((next1X == prev1X) && (next1X != x)) {
                        edge.set(i, next1X, y);
                    }
                }
            }
        }
    }

    /**
     * find the jagged line segments in the curve and return the ranges
     * of the point indexes.
     * This method searches for ledges first and then within the remaining
     * space, searches for stair cases and then 45 degree lines.
     * @param curve
     * @return
     */
    public PairIntArray findJaggedLineSegments(final PairIntArray curve) {

        //TODO: use minimum curve size
        if (curve == null || curve.getN() < 5) {
            return new PairIntArray();
        }

        /*
        search for ledges first, then in the space where ledges were not
        found, search for jagged lines (these have steps of height 1 but
        varying width).
        And as a comparison, for the spaces where ledges were not found,
        merge them if they are small and close and then research the
        merged ranges for jagged lines and then the remaining space
        for ledges.
        Combine the 2 results to make the best results.
        NOTE: need to simplify and combine the logic of the 2 searches...
         */

        PairIntArray jaggedLines1 = findLedgesInCurve(curve);

        PairIntArray remainingRanges =
            writeRangesNotAlreadyIncluded(curve, jaggedLines1);

        for (int i = 0; i < remainingRanges.getN(); i++) {

            int r0 = remainingRanges.getX(i);
            int r1 = remainingRanges.getY(i);

            PairIntArray tmpStaircaseRanges =
                findJaggedLineStaircaseSegments(curve, r0, r1);

            if (tmpStaircaseRanges != null) {
                for (int j = 0; j < tmpStaircaseRanges.getN(); j++) {
                    int s0 = tmpStaircaseRanges.getX(j);
                    int s1 = tmpStaircaseRanges.getY(j);
                    jaggedLines1.add(s0, s1);
                }
            }
        }

        sortByX(jaggedLines1);

        // merge adjacent ranges
        mergeRanges(curve, jaggedLines1);

        // search for 45 degree lines
        remainingRanges =
            writeRangesNotAlreadyIncluded(curve, jaggedLines1);
        boolean changed = false;
        for (int i = 0; i < remainingRanges.getN(); i++) {
            int r0 = remainingRanges.getX(i);
            int r1 = remainingRanges.getY(i);
            PairIntArray lineRanges =
                find45DegreeSegments(curve, r0, r1);
            if (lineRanges != null) {
                for (int j = 0; j < lineRanges.getN(); j++) {
                    int s0 = lineRanges.getX(j);
                    int s1 = lineRanges.getY(j);
                    jaggedLines1.add(s0, s1);
                    changed = true;
                }
            }
        }

        if (changed) {

            sortByX(jaggedLines1);

            // merge ranges again
            mergeRanges(curve, jaggedLines1);
        }

        return jaggedLines1;
    }

    /**
     * find the jagged line segments in the curve and return the ranges
     * of the point indexes.
     * This method searches for staircases first and then within the remaining
     * space, searches for ledges and then 45 degree lines.
     * @param curve
     * @return
     */
    public PairIntArray findJaggedLineSegments2(final PairIntArray curve) {

        //TODO: use minimum curve size
        if (curve == null || curve.getN() < 5) {
            return null;
        }

        // if have a merged larger range, this suggests that starting a
        // search with staircases and following that with search for
        // ledges might result in more total accurately found jagged lines.
        PairIntArray jaggedLines2 = findJaggedLineStaircaseSegments(
            curve, 0, curve.getN() - 1);

        PairIntArray remainingRanges = writeRangesNotAlreadyIncluded(curve,
            jaggedLines2);
        for (int i = 0; i < remainingRanges.getN(); i++) {
            int r0 = remainingRanges.getX(i);
            int r1 = remainingRanges.getY(i);
            findLedgesWithinRange(curve, r0, r1, jaggedLines2);
        }

        sortByX(jaggedLines2);

        // merge adjacent ranges
        mergeRanges(curve, jaggedLines2);

        // search for 45 degree lines
        remainingRanges =
            writeRangesNotAlreadyIncluded(curve, jaggedLines2);
        for (int i = 0; i < remainingRanges.getN(); i++) {
            int r0 = remainingRanges.getX(i);
            int r1 = remainingRanges.getY(i);
            PairIntArray lineRanges =
                find45DegreeSegments(curve, r0, r1);
            if (lineRanges != null) {
                for (int j = 0; j < lineRanges.getN(); j++) {
                    int s0 = lineRanges.getX(j);
                    int s1 = lineRanges.getY(j);
                    jaggedLines2.add(s0, s1);
                }
            }
        }

        sortByX(jaggedLines2);

        // merge ranges again
        mergeRanges(curve, jaggedLines2);

        return jaggedLines2;
    }

    private void mergeRanges(PairIntArray curve, PairIntArray ranges) {

        int n = ranges.getN();
        for (int i = (n - 1); i > 0; i--) {

            int r0 = ranges.getX(i);
            int r1 = ranges.getY(i);

            // gap between end of one range and start of the next that might
            // be part of both ranges
            if ((r0 - ranges.getY(i - 1)) < 3) {

                // check slopes before merging.

                double theta10 = calcTheta(curve, r0, r1);

                int r2 = ranges.getX(i - 1);
                int r3 = ranges.getY(i - 1);
                double theta32 = calcTheta(curve, r2, r3);

                double diffTheta = Math.abs(theta10 - theta32);

                if (diffTheta > (Math.PI/4.)) {
                    continue;
                }
                if (diffTheta > 0.1) {
                    // this may be 2 regions due to 2 different methods,
                    // the ledges method and then staircase method,
                    // so retry the staircase alone for the full range to
                    // see if the entire region is found as a single region.
                    // r2 to r1
                    PairIntArray staircaseRanges =
                        findJaggedLineStaircaseSegments(curve, r2, r1);

                    if (staircaseRanges != null) {

                        for (int j = 0; j < staircaseRanges.getN(); j++) {
                            int s0 = staircaseRanges.getX(j);
                            int s1 = staircaseRanges.getY(j);
                            log.fine("       " + s0 + " : " + s1);
                        }

                        if (staircaseRanges.getN() == 0) {
                            // do not merge
                            continue;
                        } else if (staircaseRanges.getN() == 1) {
                            // does it match the whole range r2 to r1
                            // or only one of r2:r3 and r0:r1?
                            // if only matches one range, do not merge them

                            int diff10 = Math.abs((r0 - staircaseRanges.getX(0)) +
                                (r1 - staircaseRanges.getY(0)));

                            int diff32 = Math.abs((r2 - staircaseRanges.getX(0)) +
                                (r3 - staircaseRanges.getY(0)));

                            int diff21 = Math.abs((r2 - staircaseRanges.getX(0)) +
                                (r1 - staircaseRanges.getY(0)));

                            // which one is close to zero?
                            if ((diff21 < diff10) && (diff21 < diff32)) {
                                // let these merge
                            } else if ((diff32 < diff10) && (diff32 < diff21)) {
                                // matches the range r2 to r3, don't merge
                                continue;
                            } else if ((diff10 < diff32) && (diff10 < diff21)) {
                                // matches the range r0 to r1, don't merge
                                continue;
                            }

                        } else if (staircaseRanges.getN() > 1) {
                            // keep the ranges separate
                            continue;
                        }
                    } else {
                        continue;
                    }
                }

                ranges.set(i - 1, ranges.getX(i - 1), r1);

                ranges.removeRange(i, i);
            }
        }

    }

    /**
     * find the lines composed of nearly uniform stairs and return them as
     * index ranges.  For example, a jagged line that extends from point
     * 10 to point 30 inclusive is present in the returned object as a pair
     * with (x, y) = (10, 30).
     * @param curve
     * @return
     */
    private PairIntArray findJaggedLineStaircaseSegments(final PairIntArray
        curve, int startIndex, int stopIndex) {

        //TODO: use minimum curve size
        if (curve == null || (stopIndex - startIndex) < 5) {
            return new PairIntArray();
        }

        /*
        iterate over the curve to find the nearly straight line segments.
        This is useful for quickly removing false corners due to jagged
        lines.
        -- move forward and learn dx and dy.  either dx or dy must be constant
        and have value -1 or +1.  the other dimension can only change by 0
        or by the same +1 or -1 always.  the step width between the change
        must be on average a certain value and any other steps included
        can be +1 or -1 in width (for example, if step width is 2, can have
        steps with width 1 and 3 included also).
        keep a moving average and when the just stated conditions cease,
        note the endpoints.
        --> to be sure the section is a line, make an easy to remove section:
            fit the points to a line, noting the mean and stdev of the distance
            of them from the line.
            are the results consistent with a line?  mean error is?
        -- if the segment is longer than (tbd) pixels, store it
        -- repeat the above until end of curve is reached.
        */

        PairIntArray lineSegmentRanges = new PairIntArray();

        int dx = 0;
        int dy = 0;
        int i = startIndex;
        Boolean widthIsAlongX = null;
        while ((dx == 0) || (dy == 0)) {
            i++;
            if (i >= stopIndex) {
                return lineSegmentRanges;
            }
            dx = (curve.getX(i) - curve.getX(i - 1));
            dy = (curve.getY(i) - curve.getY(i - 1));
            if (dx == 0) {
                widthIsAlongX = Boolean.FALSE;
            } else if (dy == 0) {
                widthIsAlongX = Boolean.TRUE;
            }
        }
        int start = i;

        int keepDX = dx;
        int keepDY = dy;

        if (widthIsAlongX == null) {
            while ((dx != 0) && (dy != 0)) {
                i++;
                if (i >= stopIndex) {
                    return lineSegmentRanges;
                }
                dx = (curve.getX(i) - curve.getX(i - 1));
                dy = (curve.getY(i) - curve.getY(i - 1));
                if (dx == 0) {
                    widthIsAlongX = Boolean.FALSE;
                } else if (dy == 0) {
                    widthIsAlongX = Boolean.TRUE;
                }
            }
        }
        dx = keepDX;
        dy = keepDY;
        int stepStart = startIndex;
        int lineStart = startIndex;
        int nSteps = 0;
        int sumStepWidth = 0;
        float avgStepWidth = 0;
        float firstStepWidth = -1;
        float lastStepWidth = -1;

        for (i = start; i <= stopIndex; i++) {

            int x = curve.getX(i);
            int y = curve.getY(i);

            int diffX = (x - curve.getX(i - 1));

            int diffY = (y - curve.getY(i - 1));

            // if not a continuation of current step, increment step
            if (!(
                (widthIsAlongX.booleanValue() && (diffX == dx) && (diffY == 0))
                ||
                (!widthIsAlongX.booleanValue() && (diffY == dy) && (diffX == 0))
                )
                ){

                int currentStepWidth = i - stepStart;

                if (currentStepWidth > 0) {
                    if (nSteps == 0) {
                        firstStepWidth = currentStepWidth;
                    }
                    nSteps++;
                    sumStepWidth += currentStepWidth;
                    avgStepWidth = sumStepWidth/(float)nSteps;
                    stepStart = i;
                    lastStepWidth = currentStepWidth;
                }
            }

            // if an invalid dx or dy, write the lineSegment and reset the range
            if ((i == stopIndex) ||
                (widthIsAlongX.booleanValue() && (diffX != dx))
                || (!widthIsAlongX.booleanValue() && (diffY != dy))
                ||
                (widthIsAlongX.booleanValue() && (diffY != dy) && (diffY != 0))
                ||
                (!widthIsAlongX.booleanValue() && (diffX != dx) && (diffX != 0))
                ) {

                if (nSteps > 2) {

                    int avg = Math.round(avgStepWidth);

                    int[] endSegment = validateJaggedLineSegment(curve,
                        lineStart, (i - 1), avg, dx, dy,
                        widthIsAlongX);

                    // only store if has at least 3 steps (but if avg==1, 10)
                    if (
                    ((avg == 1) &&
                        ((endSegment[0] - lineStart + 1) >= 10)
                    )
                    ||
                    ((avg > 1) &&
                        ((endSegment[0] - lineStart + 1) >= 3 * avg)
                    )) {
                        if (
                            ((lastStepWidth >= 3) && ((lastStepWidth/avg) >= 2))
                            || ((lastStepWidth == 1) && (avg > 1))
                            ){

                            int endMinus = (int)(endSegment[0] - lastStepWidth);

                            if ((endMinus - lineStart + 1) > 4) {
                                lineSegmentRanges.add(lineStart, endMinus);
                            } else {
                                lineSegmentRanges.add(lineStart, endSegment[0]);
                            }
                        } else {
                            lineSegmentRanges.add(lineStart, endSegment[0]);
                        }
                    } else {
                        i = lineStart + endSegment[1];
                    }
                }

                if (i >= stopIndex) {
                    return lineSegmentRanges;
                }

                //TODO: check the stepStart index
                stepStart = i;
                lineStart = i;

                dx = 0;
                dy = 0;
                widthIsAlongX = null;
                int tmpI = i;

                while ((dx == 0) || (dy == 0)) {
                    i++;
                    if (i >= stopIndex) {
                        return lineSegmentRanges;
                    }
                    dx = (curve.getX(i) - curve.getX(i - 1));
                    dy = (curve.getY(i) - curve.getY(i - 1));
                    if (dx == 0) {
                        widthIsAlongX = Boolean.FALSE;
                    } else if (dy == 0) {
                        widthIsAlongX = Boolean.TRUE;
                    }
                }

                keepDX = dx;
                keepDY = dy;

                if (widthIsAlongX == null) {
                    while ((dx != 0) && (dy != 0)) {
                        i++;
                        if (i >= stopIndex) {
                            return lineSegmentRanges;
                        }
                        dx = (curve.getX(i) - curve.getX(i - 1));
                        dy = (curve.getY(i) - curve.getY(i - 1));
                        if (dx == 0) {
                            widthIsAlongX = Boolean.FALSE;
                        } else if (dy == 0) {
                            widthIsAlongX = Boolean.TRUE;
                        }
                    }

                } else {
                    // back track to find where the current linestart
                    // should be between tmpI and i
                    boolean iChanged = false;
                    for (int j = (i - 1); j > tmpI; j--) {
                        diffX = (curve.getX(j) - curve.getX(j - 1));
                        diffY = (curve.getY(j) - curve.getY(j - 1));
                        if (widthIsAlongX && (diffY == 0) && (diffX == dx)) {
                            i = j;
                            iChanged = true;
                        } else if (!widthIsAlongX && (diffX == 0)
                            && (diffY == dy)) {
                            i = j;
                            iChanged = true;
                        }
                    }
                    if (iChanged) {
                        i--;
                    }
                }

                dx = keepDX;
                dy = keepDY;
                sumStepWidth = 0;
                avgStepWidth = 0;
                nSteps = 0;
                lastStepWidth = -1;
                firstStepWidth = -1;

                lineStart = i;
            }
        }

        return lineSegmentRanges;

    }

    /**
     * find the lines composed of nearly uniform stairs and return them as
     * index ranges.  For example, a jagged line that extends from point
     * 10 to point 30 inclusive is present in the returned object as a pair
     * with (x, y) = (10, 30).
     * @param curve
     * @return
     */
    private PairIntArray find45DegreeSegments(final PairIntArray
        curve, int startIndex, int stopIndex) {

        //TODO: use minimum curve size
        if (curve == null || (stopIndex - startIndex) < 5) {
            return new PairIntArray();
        }

        int minNSteps = 4;

        PairIntArray lineSegmentRanges = new PairIntArray();

        int dx = 0;
        int dy = 0;
        int i = startIndex;
        while ((dx == 0) || (dy == 0)) {
            i++;
            if (i >= stopIndex) {
                return lineSegmentRanges;
            }
            dx = (curve.getX(i) - curve.getX(i - 1));
            dy = (curve.getY(i) - curve.getY(i - 1));
        }
        int start = i;

        int stepStart = startIndex;
        int lineStart = startIndex;
        int nSteps = 0;
        int sumStepWidth = 0;
        float avgStepWidth = 0;
        float firstStepWidth = -1;
        float lastStepWidth = -1;

        for (i = start; i <= stopIndex; i++) {

            int x = curve.getX(i);
            int y = curve.getY(i);

            int diffX = (x - curve.getX(i - 1));

            int diffY = (y - curve.getY(i - 1));

            int currentStepWidth = i - stepStart;

            if (currentStepWidth > 0) {
                if (nSteps == 0) {
                    firstStepWidth = currentStepWidth;
                }
                nSteps++;
                sumStepWidth += currentStepWidth;
                avgStepWidth = sumStepWidth/(float)nSteps;
                stepStart = i;
                lastStepWidth = currentStepWidth;
            }

            // if an invalid dx or dy, write the lineSegment and reset the range
            if ((i == stopIndex) || (diffX != dx) || (diffY != dy)
                ) {

                if (nSteps > minNSteps) {

                    lineSegmentRanges.add(lineStart, i - 1);

                }

                if (i >= stopIndex) {
                    return lineSegmentRanges;
                }

                //TODO: check the stepStart index
                stepStart = i;
                lineStart = i;

                dx = 0;
                dy = 0;

                while ((dx == 0) || (dy == 0)) {
                    i++;
                    if (i >= stopIndex) {
                        return lineSegmentRanges;
                    }
                    dx = (curve.getX(i) - curve.getX(i - 1));
                    dy = (curve.getY(i) - curve.getY(i - 1));
                }

                sumStepWidth = 0;
                avgStepWidth = 0;
                nSteps = 0;
                lastStepWidth = -1;
                firstStepWidth = -1;

                lineStart = i;
            }
        }

        return lineSegmentRanges;

    }

    /**
     * validate that a line segment has steps only within +- 1 of
     * step stepWidth.  returns endIndex if entire region fits those
     * characteristics, else returns the last index where it does.
     *
     * @param curve
     * @param startIndex
     * @param stopIndex last index of line segment, inclusive
     * @param stepWidth
     * @param dy
     * @param dy
     * @param widthIsAlongX
     * @return
     */
     int[] validateJaggedLineSegment(final PairIntArray curve,
        int startIndex, int stopIndex, int stepWidth, int dx, int dy,
        Boolean widthIsAlongX) {

        //TODO: use minimum curve size
        if (curve == null || curve.getN() < 5) {
            return new int[]{-1, -1};
        }

        int plusMinusWidth = 3;

        int n = curve.getN();

        int start = startIndex + 1;

        int stepStart = startIndex;

        int i;
        for (i = start; i <= stopIndex; i++) {

            int diffX = (curve.getX(i) - curve.getX(i - 1));

            int diffY = (curve.getY(i) - curve.getY(i - 1));

            if ((widthIsAlongX.booleanValue() && (diffX != dx)) ||
                (widthIsAlongX.booleanValue() && (diffY != dy)
                && (diffY != 0)) ||
                (!widthIsAlongX.booleanValue() && (diffY != dy)) ||
                (!widthIsAlongX.booleanValue() && (diffX != dx)
                && (diffX != 0)) ) {

                int currentStepWidth = i - stepStart;
                if (currentStepWidth > 0) {
                    if (Math.abs(currentStepWidth - stepWidth) > plusMinusWidth) {
                        return new int[]{(stepStart - 1), currentStepWidth};
                    }
                } else {
                    return new int[]{(stepStart - 1), currentStepWidth};
                }

                /*
                (widthIsAlongX.booleanValue() && (diffY != dy)
                && (diffY != 0))
                      --> (stepStart - 1)
                */
                return new int[]{(i - 1), 0};
            }

            // else, if just stepped up or is last index, check step size

            if ((widthIsAlongX.booleanValue() && (diffY == dy)) ||
                (!widthIsAlongX.booleanValue() && (diffX == dx)) ||
                (i == stopIndex)
            ) {

                int currentStepWidth = i - stepStart;

                if ((stepStart == 0) && (stepWidth == 1) &&
                    (currentStepWidth/stepWidth > 1)) {

                    return new int[]{0, currentStepWidth};

                } else if (currentStepWidth > 0) {

                    if (Math.abs(currentStepWidth - stepWidth) > plusMinusWidth) {
                        return new int[]{(stepStart - 1), currentStepWidth};
                    } else if (i == stopIndex) {
                        return new int[]{i, currentStepWidth};
                    }

                    stepStart = i;

                } else if (currentStepWidth == 0) {

                    return new int[]{(stepStart - 1), 0};
                }
            }
        }

        return new int[]{(i - 1), 0};
    }

     /**
      * write the set difference of the given set of ranges, indexRanges,
      * to create the set of ranges not included in indexRanges.  Note the
      * large universe that both are subsets of is curve.
      */
     private PairIntArray writeRangesNotAlreadyIncluded(PairIntArray curve,
        PairIntArray indexRanges) {

        PairIntArray output = new PairIntArray();

        int n = curve.getN();

        if (indexRanges.getN() == 0) {
            output.add(0, n - 1);
        } else {
            int idx0 = indexRanges.getX(0);
            if (idx0 > 0) {
                output.add(0, idx0 - 1);
            }
            for (int si = 1; si < indexRanges.getN(); si++) {
                output.add(indexRanges.getY(si - 1), indexRanges.getX(si));
            }
            output.add(indexRanges.getY(indexRanges.getN() - 1),
                curve.getN() - 1);
        }

        return output;
     }

     /**
      * in the curve points that are not within the staircaseSegmentRanges,
      * look for the single pixel ledge in a long stretch of a line and store
      * the entire range.  There may be more than one single pixel range
      * within a range.  a range is stored in the return array as a
      * point (x,y) = (start of range, stop of range inclusive).
      * @param curve
      * @param staircaseSegmentRanges
      * @return
      */
    PairIntArray findLedgesInCurve(PairIntArray curve) {

        /*
        looking for long stretch of line that changes by 1 pixel and then
        continues in a long line
        */

        PairIntArray allLedges = new PairIntArray();

        findLedgesWithinRange(curve, 0, curve.getN() - 1, allLedges);

        return allLedges;
    }

    /**
     * find any ledges within the range start to stop, inclusive and return
     * them as indexes of the curve.  For example, a ledge extending from
     * point 10 to point 30 inclusive is in a pair in allLedges
     * as (x,y) = (10, 30);
     * @param curve set of x,y points which comprise a curve
     * @param start first index of curve to search, inclusive
     * @param stop last index of curve to search, inclusive
     * @param allLedges the set of ranges to add the results of this too.
     * It's the output for the method.
     */
    private void findLedgesWithinRange(PairIntArray curve, int start, int stop,
        PairIntArray allLedges) {

        if ((stop - start + 1) < (2*minLedgeWidth)) {
            return;
        }

        // similar o findJaggedLineStaircaseSegments, but with a step size of
        // "1"

        int dx = 0;
        int dy = 0;
        int i = start;
        Boolean runIsAlongX = null;
        // looking for straight lines of x or y
        while (!((dx == 0) && (dy != 0)) && !((dy == 0) && (dx != 0))) {
            i++;
            if (i > (stop - 1)) {
                return;
            }
            dx = (curve.getX(i) - curve.getX(i - 1));
            dy = (curve.getY(i) - curve.getY(i - 1));
        }

        if (dx == 0) {
            runIsAlongX = Boolean.FALSE;
        } else if (dy == 0) {
            runIsAlongX = Boolean.TRUE;
        }

        int tmpI = i;

        // back track to find where the current linestart
        // should be between start and i
        for (int j = (i - 1); j >= (start + 1); j--) {
            int diffX = (curve.getX(j) - curve.getX(j - 1));
            int diffY = (curve.getY(j) - curve.getY(j - 1));
            if ((diffY != dy) || (diffX != dx)) {
                i = j;
            }
        }

        int lineStart = i;

        PairIntArray tmp = new PairIntArray();
        Boolean tmpRunIsAlongX = null;
        int tmpDX = -1;
        int tmpDY = -1;

        for (i = (lineStart + 1); i <= stop; i++) {

            int x = curve.getX(i);
            int y = curve.getY(i);

            int diffX = (x - curve.getX(i - 1));

            int diffY = (y - curve.getY(i - 1));

            /* if there's a break in the line:
                  temporarily store the section so far.

                  if the next segment is consecutive and has same runIsAlongX
                  and same diffX and diffY,
                      continue with same tmp storage,
                  else {
                     inspect storage and add to allLedges if looks like a ledge,
                     then clear the tmp storage and the last vars"
                  }
            */

            boolean runStopped = (runIsAlongX && (diffY != 0)) ||
                (!runIsAlongX && (diffX != 0));

            if ((i == stop) ||
                runStopped ||
                (runIsAlongX && (diffX != dx)) ||
                (!runIsAlongX && (diffY != dy)) ) {

                int rs = i - lineStart;

                if (rs >= minLedgeWidth) {
                    if (i == stop) {
                        if (runStopped) {
                            tmp.add(lineStart, i - 1);
                        } else {
                            tmp.add(lineStart, i);
                        }
                    } else {
                        tmp.add(lineStart, i - 1);
                    }
                    tmpRunIsAlongX = runIsAlongX;
                    tmpDX = dx;
                    tmpDY = dy;
                } else if ((i == (curve.getN() - 1)) /*&& (tmp.getN() > 0) &&
                    (Math.abs(lineStart - tmp.getY(tmp.getN() - 1)) < 2)*/) {
                    tmp.add(lineStart, i);
                } else if (tmp.getN() == 1) {
                    allLedges.add(tmp.getX(0), tmp.getY(0));
                    tmp = new PairIntArray();
                }

                if (i != stop) {
                    // find the next line segment
                    dx = 0;
                    dy = 0;
                    runIsAlongX = null;
                    // looking for straight lines of x or y
                    tmpI = i;
                    while (!((dx == 0) && (dy != 0)) &&
                        !((dy == 0) && (dx != 0))) {
                        i++;
                        if (i >= stop) {
                            break;
                        }
                        dx = (curve.getX(i) - curve.getX(i - 1));
                        dy = (curve.getY(i) - curve.getY(i - 1));
                    }

                    if (i < stop) {
                        if (dx == 0) {
                            runIsAlongX = Boolean.FALSE;
                        } else if (dy == 0) {
                            runIsAlongX = Boolean.TRUE;
                        }

                        // back track to find where the current linestart
                        // should be between tmpI and i
                        for (int j = (i - 1); j >= tmpI; j--) {
                            diffX = (curve.getX(j) - curve.getX(j - 1));
                            diffY = (curve.getY(j) - curve.getY(j - 1));
                            if ((diffY != dy) || (diffX != dx)) {
                                i = j;
                            }
                        }

                        lineStart = i;
                    }
                }

                int tmpN = tmp.getN();

                // if this is not consecutive segment,
                // decide whether to store, then reset tmp
                if ((i >= stop) || (
                    (tmp.getN() > 0) &&
                    !(
                        (runIsAlongX.compareTo(tmpRunIsAlongX) == 0)
                        && (dx == tmpDX) && (tmpDY == dy)
                        && ((lineStart - (tmp.getY(tmpN - 1)) < 3))
                    )
                    )
                    ) {

                    if (tmp.getN() > 1 || ((tmp.getN() > 0)
                        && (i == (curve.getN() - 1)))) {

                        // need to avoid removing a partial corner

                        boolean keep = true;

                        // check that the lines are not wrapping around a curve
                        if (tmp.getN() > 1) {
                            int idx0f = tmp.getX(0);
                            int idx0l = tmp.getY(0);
                            double theta0 = calcTheta(curve, idx0f, idx0l);

                            for (int j = 1; j < tmp.getN(); j++) {
                                idx0f = tmp.getX(j);
                                idx0l = tmp.getY(j);
                                double theta1 = calcTheta(curve, idx0f, idx0l);

                                // don't add corners
                                double diff = Math.abs(theta0 - theta1);
                                if (diff > Math.PI/4.) {
                                    keep = false;
                                    break;
                                }
                                theta0 = theta1;
                            }
                        }

                        if (keep) {
                            allLedges.add(tmp.getX(0), tmp.getY(tmp.getN() - 1));
                        }

                        tmp = new PairIntArray();
                    }
                }
            }
        }
    }

    public void sortByX(PairIntArray curve) {
        if (curve.getN() < 2) {
            return;
        }
        sortByX(curve, 0, curve.getN() - 1);
    }

    private void sortByX(PairIntArray curve, int idxLo, int idxHi) {
        if (idxLo < idxHi) {
            int idxMid = partitionByX(curve, idxLo, idxHi);
            sortByX(curve, idxLo, idxMid - 1);
            sortByX(curve, idxMid + 1, idxHi);
        }
    }

    private int partitionByX(PairIntArray curve, int idxLo, int idxHi) {

        int x = curve.getX(idxHi);  //for comparison
        int store = idxLo - 1;      //store to swap after pivot

        for (int i = idxLo; i < idxHi; i++) {
            if (curve.getX(i) <= x) {
                store++;
                int swapX = curve.getX(store);
                int swapY = curve.getY(store);
                curve.set(store, curve.getX(i), curve.getY(i));
                curve.set(i, swapX, swapY);
            }
        }
        store++;

        int swapX = curve.getX(store);
        int swapY = curve.getY(store);
        curve.set(store, curve.getX(idxHi), curve.getY(idxHi));
        curve.set(idxHi, swapX, swapY);

        return store;
    }

    private double calcTheta(PairIntArray curve, int idx0, int idx1) {

        int x10 = curve.getX(idx1) - curve.getX(idx0);
        int y10 = curve.getY(idx1) - curve.getY(idx0);
        double theta;
        if (x10 == 0) {
            theta = (y10 < 0) ? 1.5 * Math.PI : 0.5 * Math.PI;
        } else {
            theta = Math.atan((double) y10 / (double) x10);
        }

        return theta;
    }

    /**
     * calculate centroid and then search nearby for the first closest
     * point in points.  the result approximates a sort by x and y and
       the returned values are actually present in points which is
       not guaranteed for calculateXYCentroids alone.
       * If that fails, it sorts all of points by y then x for the same
       * y and returns the median value present.
     * @param points
     * @return
     */
    public double[] calculateXYMedian(PairIntArray points) {

        Set<PairInt> points2 = new HashSet<PairInt>();
        for (int i = 0; i < points.getN(); i++) {
            int x = points.getX(i);
            int y = points.getY(i);
            PairInt p = new PairInt(x, y);
            points2.add(p);
        }

        double[] cen = calculateXYCentroids(points);

        int xc = (int)Math.round(cen[0]);
        int yc = (int)Math.round(cen[1]);

        int xMin = MiscMath.findMin(points.getX());
        int xMax = MiscMath.findMax(points.getX());
        int dXEnd = (xMax - xc);
        if ((xc - xMin) < dXEnd) {
            dXEnd = (xc - xMin);
        }
        int yMin = MiscMath.findMin(points.getY());
        int yMax = MiscMath.findMax(points.getY());
        int dYEnd = (yMax - yc);
        if ((yc - yMin) < dYEnd) {
            dYEnd = (yc - yMin);
        }

        int[] signs = new int[]{-1, 1};

        for (int dx = 0; dx <= dXEnd; dx++) {
            for (int dy = 0; dy <= dYEnd; dy++) {
                for (int fx : signs) {
                    dx *= fx;
                    for (int fy : signs) {
                        dy *= fy;
                    }
                }

                int x = xc + dx;
                int y = yc + dx;

                PairInt p = new PairInt(x, y);

                if (points2.contains(p)) {
                    return new double[]{x, y};
                }
            }
        }

        log.warning("should not reach here.  should have found a point within "
            + " range of values");

        PairIntArray points3 = points.copy();

        MultiArrayMergeSort.sortByYThenX(points3);

        int idx = points3.getN()/2;

        return new double[]{points3.getX(idx), points3.getY(idx)};
    }

    public PairInt[] findClosestPair(Set<PairInt> set0, Set<PairInt> set1) {

        if (set0 == null) {
            throw new IllegalArgumentException("set0 cannot be null");
        }

        if (set1 == null) {
            throw new IllegalArgumentException("set1 cannot be null");
        }

        if (set0.isEmpty()) {
            throw new IllegalArgumentException("set0 cannot be empty");
        }

        if (set1.isEmpty()) {
            throw new IllegalArgumentException("set1 cannot be empty");
        }

        //TODO: consider other algorithms besides brute force

        double minDistSq = Double.MAX_VALUE;
        PairInt p0 = null;
        PairInt p1 = null;

        for (PairInt s0 : set0) {

            double x = s0.getX();
            double y = s0.getY();

            for (PairInt s1 : set1) {

                double x1 = s1.getX();
                double y1 = s1.getY();

                double diffX = x1 - x;
                double diffY = y1 - y;

                double distSq = (diffX * diffX) + (diffY * diffY);

                if (distSq < minDistSq) {
                    minDistSq = distSq;
                    p0 = s0;
                    p1 = s1;
                }
            }
        }

        return new PairInt[]{p0, p1};
    }

    public void populateGapsWithInterpolation(Set<PairInt> points) {

        // probably many ways to do this... ordered point algorithm?
        // dfs to find connected groups, then connect the closest among those?

        int[] minMaxXY = MiscMath.findMinMaxXY(points);

        DFSConnectedGroupsFinder finder = new DFSConnectedGroupsFinder();
        finder.setMinimumNumberInCluster(1);
        finder.findConnectedPointGroups(points, minMaxXY[1] + 1, minMaxXY[3] + 1);

        int nIter = 0;
        int nMaxIter = 10;

        int nGroups = finder.getNumberOfGroups();

        while ((nGroups > 1) && (nIter < nMaxIter)) {

            // find the closest pair of points between any 2 groups

            double minDistSq = Double.MAX_VALUE;
            PairInt minDistPoint0 = null;
            PairInt minDistPoint1 = null;
            int minDistGroupId0 = -1;
            int minDistGroupId1 = -1;

            for (int g0Idx = 0; g0Idx < nGroups; g0Idx++) {

                Set<PairInt> g0 = finder.getXY(g0Idx);

                for (int g1Idx = 0; g1Idx < nGroups; g1Idx++) {

                    if (g0Idx == g1Idx) {
                        continue;
                    }

                    Set<PairInt> g1 = finder.getXY(g1Idx);

                    PairInt[] closestPair = findClosestPair(g0, g1);

                    if (closestPair == null) {
                        continue;
                    }

                    double x0 = closestPair[0].getX();
                    double y0 = closestPair[0].getY();

                    double x1 = closestPair[1].getX();
                    double y1 = closestPair[1].getY();

                    double diffX = x1 - x0;
                    double diffY = y1 - y0;

                    double distSq = (diffX * diffX) + (diffY * diffY);

                    if (distSq < minDistSq) {
                        minDistSq = distSq;
                        minDistPoint0 = closestPair[0];
                        minDistPoint1 = closestPair[1];
                        minDistGroupId0 = g0Idx;
                        minDistGroupId1 = g1Idx;
                    }
                }
            }

            if (minDistPoint0 != null) {

                double x1 = minDistPoint1.getX();
                double y1 = minDistPoint1.getY();

                double dxDivDy = (minDistPoint0.getX() - x1)/
                    (minDistPoint0.getY() - y1);

                /*
                x0 - x1
                ------- = dxDivDy
                y0 - y1

                x0 - x1 = (y0 - y1) * dyDivDx;
                x0 = x1 + (y0 - y1) * dyDivDx;
                */
                int startY, stopY;
                if (minDistPoint0.getY() < minDistPoint1.getY()) {
                    startY = minDistPoint0.getY();
                    stopY = minDistPoint1.getY();
                } else {
                    startY = minDistPoint1.getY();
                    stopY = minDistPoint0.getY();
                }

                for (int y = startY; y <= stopY; y++) {

                    int x = (int)Math.round(x1 + (y - y1) * dxDivDy);

                    PairInt p = new PairInt(x, y);

                    boolean added = points.add(p);
                }

                /*
                y0 - y1
                ------- = dyDivDx
                x0 - x1

                y0 - y1 = (x0 - x1) * dxDivDy;
                y0 = y1 + (x0 - x1) * dxDivDy;
                */
                double dyDivDx = (minDistPoint0.getY() - y1)/
                    (minDistPoint0.getX() - x1);

                int startX, stopX;
                if (minDistPoint0.getX() < minDistPoint1.getX()) {
                    startX = minDistPoint0.getX();
                    stopX = minDistPoint1.getX();
                } else {
                    startX = minDistPoint1.getX();
                    stopX = minDistPoint0.getX();
                }

                for (int x = startX; x <= stopX; x++) {

                    int y = (int)Math.round(y1 + (x - x1) * dyDivDx);

                    PairInt p = new PairInt(x, y);

                    boolean added = points.add(p);
                }
            }

            finder = new DFSConnectedGroupsFinder();
            finder.setMinimumNumberInCluster(1);
            finder.findConnectedPointGroups(points, minMaxXY[1] + 1, minMaxXY[3] + 1);

            nGroups = finder.getNumberOfGroups();

            nIter++;
        }
    }

    public void pruneSpurs(List<PairIntArray> tmpEdges, int imageWidth,
        int imageHeight) {

        //TODO: improve this to follow a spur when found and remove the
        //  resulting spurs from the remove action

        int nIter = 0;
        int nMaxIter = 3;
        int nRemoved = 0;

        while (nIter < nMaxIter) {

            nRemoved = 0;

            /*
            0 0 0    0 0 0
            0 1 0    0 1 0
            0 _ _    _ _ 0
            */
            // indexes that have to be zeros that is, not within an edge
            int[] topXIdx = new int[]{-1, 0, 1, -1, 1};
            int[] topYIdx = new int[]{ 1, 1, 1,  0, 0};
            int[] leftXIdx = new int[]{-1};
            int[] leftYIdx = new int[]{-1};
            int[] rightXIdx = new int[]{1};
            int[] rightYIdx = new int[]{-1};
            // one must be a zero:
            int[] leftOrZeroXIdx = new int[]{0, 1};
            int[] leftOrZeroYIdx = new int[]{-1, -1};
            int[] rightOrZeroXIdx = new int[]{-1, 0};
            int[] rightOrZeroYIdx = new int[]{-1, -1};

            for (int r = 0; r < 4; r++) {
                if (r > 0) {
                    rotateIndexesBy90(topXIdx, topYIdx);
                    rotateIndexesBy90(leftXIdx, leftYIdx);
                    rotateIndexesBy90(rightXIdx, rightYIdx);
                    rotateIndexesBy90(leftOrZeroXIdx, leftOrZeroYIdx);
                    rotateIndexesBy90(rightOrZeroXIdx, rightOrZeroYIdx);
                }
                for (int i = 0; i < tmpEdges.size(); i++) {

                    PairIntArray edge = tmpEdges.get(i);

                    // skip the endpoints
                    for (int j = (edge.getN() - 1); j > 0; j--) {

                        if ((j < 0) || (j > (edge.getN() - 1))) {break;}

                        int x = edge.getX(j);
                        int y = edge.getY(j);

                        // TODO: consider a smaller range to search than +-10
                        int start = j - 10;
                        if (start < 0) {
                            continue;
                        }
                        int stop = j + 10;
                        if (stop > (edge.getN() - 1)) {
                            continue;
                        }
                        // "notFound" is "all are zeroes"
                        boolean notFound = notFound(edge, topXIdx, topYIdx, x, y,
                            start, stop);
                        if (notFound) {
                            notFound = notFound(edge, leftXIdx, leftYIdx, x, y,
                                start, stop);
                            if (notFound) {
                                if (hasAtLeastOneZero(edge, leftOrZeroXIdx,
                                    leftOrZeroYIdx, x, y, start, stop,
                                    imageWidth, imageHeight)) {

                                    edge.removeRange(j, j);
                                    nRemoved++;
                                }
                            } else {
                                notFound = notFound(edge, rightXIdx, rightYIdx, x, y,
                                    start, stop);
                                if (notFound) {
                                    if (hasAtLeastOneZero(edge, rightOrZeroXIdx,
                                        rightOrZeroYIdx, x, y, start, stop,
                                        imageWidth, imageHeight)) {

                                        edge.removeRange(j, j);
                                        nRemoved++;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            log.info("nRemoved=" + nRemoved);
            if (nRemoved == 0) {
                break;
            }
            nIter++;
        }
    }

    /**
     * edge contains point pairs of (x, y).  this method searches edge to assert
     * that all of the specified points are not present, else returns as soon
     * as one of the specified points is found.
     * The specified points are
     * (xIndex + xIndexOffsets[i], yIndex + yIndexOffsets[i]).
     * The search through edge is done from startIdx to stopIdx, inclusive.
     * @param edge
     * @param xIndexOffsets
     * @param yIndexOffsets
     * @param xIndex
     * @param yIndex
     * @param startIdx
     * @param stopIdx
     * @return
     */
    private boolean notFound(PairIntArray edge,
        int[] xIndexOffsets, int[] yIndexOffsets,
        int xIndex, int yIndex, int startIdx, int stopIdx) {

StringBuilder sb2 = new StringBuilder();

        for (int j = 0; j < xIndexOffsets.length; j++) {
            int xFind = xIndex + xIndexOffsets[j];
            int yFind = yIndex + yIndexOffsets[j];
sb2.append(String.format("(%d, %d)\n", xFind, yFind));
            for (int i = startIdx; i <= stopIdx; i++) {
                int x = edge.getX(i);
                int y = edge.getY(i);
                if ((x == xFind) && (y == yFind)) {
                    return false;
                }
            }
        }

StringBuilder sb = new StringBuilder();
for (int i = 0; i < edge.getN(); i++) {
     int x = edge.getX(i);
     int y = edge.getY(i);
     sb.append(String.format("%d)  (%d, %d)\n", i, x, y));
}

        return true;
    }

    /**
     * edge contains point pairs of (x, y).  this method searches edge to assert
     * that all of the specified points are not present, else returns as soon
     * as one of the specified points is found.
     * The specified points are
     * (xIndex + xIndexOffsets[i], yIndex + yIndexOffsets[i]).
     * The search through edge is done from startIdx to stopIdx, inclusive.
     * NOTE: also checks for whether nulling the pixel at (xP, yP) is
     * connected to the wall and returns false if it is.
     * @param edge
     * @param xIndexOffsets
     * @param yIndexOffsets
     * @param xIndex
     * @param yIndex
     * @param startIdx
     * @param stopIdx
     * @return
     */
    private boolean hasAtLeastOneZero(PairIntArray edge,
        int[] xIndexOffsets, int[] yIndexOffsets,
        int xP, int yP, int startIdx, int stopIdx, int imageWidth, int imageHeight) {

        // side logic of checking whether connected to wall.  return false if so
        if ((xP == 0) || (yP == 0)) {
            return false;
        }
        if ((xP == (imageWidth - 1)) || (yP == (imageHeight - 1))) {
            return false;
        }

        int nOnes = 0;
        for (int j = 0; j < xIndexOffsets.length; j++) {
            int xFind = xP + xIndexOffsets[j];
            int yFind = yP + yIndexOffsets[j];
            for (int i = startIdx; i <= stopIdx; i++) {
                int x = edge.getX(i);
                int y = edge.getY(i);
                if ((x == xFind) && (y == yFind)) {
                    nOnes++;
                }
            }
        }

        return (nOnes < (xIndexOffsets.length));
    }

    private void rotateIndexesBy90(int[] xOffsetIndexes, int[] yOffsetIndexes) {

        for (int i = 0; i < xOffsetIndexes.length; i++) {
            int xoff = xOffsetIndexes[i];
            int yoff = yOffsetIndexes[i];

            switch(xoff) {
                case -1:
                    switch (yoff) {
                        case -1:
                            xOffsetIndexes[i] = -1;
                            yOffsetIndexes[i] = 1;
                            break;
                        case 0:
                            xOffsetIndexes[i] = 0;
                            yOffsetIndexes[i] = 1;
                            break;
                        case 1:
                            xOffsetIndexes[i] = 1;
                            yOffsetIndexes[i] = 1;
                            break;
                    }
                    break;
                case 0:
                    switch (yoff) {
                        case -1:
                            xOffsetIndexes[i] = -1;
                            yOffsetIndexes[i] = 0;
                            break;
                        case 0:
                            // remains same
                            break;
                        case 1:
                            xOffsetIndexes[i] = 1;
                            yOffsetIndexes[i] = 0;
                            break;
                    }
                    break;
                case 1:
                    switch (yoff) {
                        case -1:
                            xOffsetIndexes[i] = -1;
                            yOffsetIndexes[i] = -1;
                            break;
                        case 0:
                            xOffsetIndexes[i] = 0;
                            yOffsetIndexes[i] = -1;
                            break;
                        case 1:
                            xOffsetIndexes[i] = 1;
                            yOffsetIndexes[i] = -1;
                            break;
                    }
                    break;
            }
        }
        /*
        1 2 3    7 4 1
        4 5 6      5 2
        7 _ _      6 3

        (-1,1)    (1,1)   0
        (0,1)     (1,0)   1
        (1,1)     (1,-1)  2
        (-1,0)    (0,1)   3
        (0,0)     (0,0)   4
        (1,0)     (0,-1)  5
        (-1,-1)   (-1,1)  6
        (0,-1)    (-1,0)  7
        (1,-1)    (-1,-1) 8
        */
    }

    public void straightenLines(Set<PairInt> points,
        GreyscaleImage edgeGuideImage) {

        if (edgeGuideImage == null) {
            return;
        }

        /*
        To move a pixel from one location to a better if possible
        means determining if the move does not break any line connections.
        Since the line widths have already been reduced to widths of '1',
        this should just be a matter of noting which points it is connected
        to and only choose points which are adjacent to the connected.

              V
              *                            V
           *     *  *  can be moved to  *  *  *  *

        The neighbors of p are p0 and p1 for example, so
        points which are within 1 pixel of p, p0, and p1 found as the centroid
        of them +- 1 pixel radius.

        Goal is to find if a point in points can be moved within a pixel's
        distance, to a position which is closer to the brightest pixel in the
        edgeGuideImage within range without breaking connections.

        For each point in points:
            -- find the adjacent points.
            -- determine a centroid for them and the point.
            -- iterate around the 8 neighboring pixels of point
               -- initialize maxIntensity w/ the current points's edgeGuideImage
                  intensity.
               -- if the pixel is further than 1 from the centroid, discard it,
                  else, compare the pixel's edgeGuideImage with maxIntensity and
                  keep if larger.
            -- if maxIntensityPoint is not null, move the current point to
               it (by adding point to the remove list and adding the new location
               to the add list).
        */

        int imageWidth = edgeGuideImage.getWidth();

        int imageHeight = edgeGuideImage.getHeight();

        double onePixDist = Math.sqrt(2);

        Set<PairInt> tmpPointsAdded = new HashSet<PairInt>();
        Set<PairInt> tmpPointsRemoved = new HashSet<PairInt>();

        Set<PairInt> outputNeighbors = new HashSet<PairInt>();

        for (PairInt p : points) {

            int x = p.getX();
            int y = p.getY();

            findNeighbors(x, y, outputNeighbors, points,
                tmpPointsAdded, tmpPointsRemoved, imageWidth, imageHeight);

            int nBrs = outputNeighbors.size();

            if (nBrs == 0) {
                continue;
            }

            // determine centroid
            double xc = p.getX();
            double yc = p.getY();
            for (PairInt p2 : outputNeighbors) {
                xc += p2.getX();
                yc += p2.getY();
            }
            xc /= (double)(nBrs + 1);
            yc /= (double)(nBrs + 1);

            // find highest intensity neighbor within 1 pix of centroid
            int maxIntensity = edgeGuideImage.getValue(x, y);
            PairInt maxIntensityPoint = null;

            for (int i = 0; i < eightNeighborsX.length; ++i) {
                int x2 = x + eightNeighborsX[i];
                int y2 = y + eightNeighborsY[i];
                if ((x2 < 0) || (x2 > (imageWidth - 1)) || (y2 < 0) ||
                    (y2 > (imageHeight - 1))) {
                    continue;
                }
                PairInt p2 = new PairInt(x2, y2);
                // discard if it's already a point
                if (outputNeighbors.contains(p2) || tmpPointsAdded.contains(p2)
                    || points.contains(p2)) {
                    continue;
                }

                // this is a vacant pixel.

                // check that it is within 1 pixel of (xc, yc)
                double diffX = x2 - xc;
                double diffY = y2 - yc;
                double dist = Math.sqrt((diffX * diffX) + (diffY * diffY));

                if (dist <= (onePixDist/2.)) {
                    int v = edgeGuideImage.getValue(x2, y2);
                    if (v > maxIntensity) {
                        maxIntensity = v;
                        maxIntensityPoint = p2;
                    }
                }
            }
            if (maxIntensityPoint != null) {
                // "change location" of the point.
                tmpPointsRemoved.add(p);
                tmpPointsRemoved.remove(maxIntensityPoint);
                tmpPointsAdded.add(maxIntensityPoint);
            }
        }

        int nCorrections = tmpPointsRemoved.size() + tmpPointsAdded.size();

        for (PairInt p2 : tmpPointsRemoved) {
            points.remove(p2);
        }
        for (PairInt p2 : tmpPointsAdded) {
            points.add(p2);
        }

        log.fine("method " + MiscDebug.getInvokingMethodName() + " nc=" +
            Integer.toString(nCorrections));
    }

    protected void findNeighbors(int x, int y, Set<PairInt> outputNeighbors,
        Set<PairInt> points, Set<PairInt> tmpAddedPoints,
        Set<PairInt> tmpRemovedPoints, int imageWidth, int imageHeight) {

        outputNeighbors.clear();

        for (int i = 0; i < eightNeighborsX.length; i++) {

            int x2 = x + eightNeighborsX[i];
            int y2 = y + eightNeighborsY[i];

            if ((x2 < 0) || (x2 > (imageWidth - 1)) || (y2 < 0) ||
                (y2 > (imageHeight - 1))) {
                continue;
            }

            PairInt p2 = new PairInt(x2, y2);

            if (tmpRemovedPoints.contains(p2)) {
                continue;
            }
            if (tmpAddedPoints.contains(p2) || points.contains(p2)) {
                outputNeighbors.add(p2);
            }
        }
    }

    public Set<PairInt> findNeighbors(int x, int y, Set<PairInt> points,
        int imageWidth, int imageHeight) {

        Set<PairInt> neighbors = new HashSet<PairInt>();

        for (int i = 0; i < eightNeighborsX.length; i++) {

            int x2 = x + eightNeighborsX[i];
            int y2 = y + eightNeighborsY[i];

            if ((x2 < 0) || (x2 > (imageWidth - 1)) || (y2 < 0) ||
                (y2 > (imageHeight - 1))) {
                continue;
            }

            PairInt p2 = new PairInt(x2, y2);

            if (points.contains(p2)) {
                neighbors.add(p2);
            }
        }

        return neighbors;
    }

    public void findNeighbors(int x, int y, Set<PairInt> outputNeighbors,
        Set<PairInt> points, Set<PairInt> excludePoints, int imageWidth, int imageHeight) {

        outputNeighbors.clear();

        for (int i = 0; i < eightNeighborsX.length; i++) {

            int x2 = x + eightNeighborsX[i];
            int y2 = y + eightNeighborsY[i];

            if ((x2 < 0) || (x2 > (imageWidth - 1)) || (y2 < 0) ||
                (y2 > (imageHeight - 1))) {
                continue;
            }

            PairInt p2 = new PairInt(x2, y2);

            if (excludePoints.contains(p2)) {
                continue;
            }
            if (points.contains(p2)) {
                outputNeighbors.add(p2);
            }
        }
    }

    /**
     * iterate through points, counting the number of pixels on the image
     * boundaries, and return true if the number reaches numberOfPixels.
     * @param numberOfPixels the number of pixels for which to return true
     * if they are on the image boundaries.
     * @param points
     * @param imageWidth
     * @param imageHeight
     * @return
     */
    public boolean hasNumberOfPixelsOnImageBoundaries(int numberOfPixels,
        Set<PairInt> points, int imageWidth, int imageHeight) {

        int n = 0;

        for (PairInt p : points) {

            int x = p.getX();
            int y = p.getY();

            if ((x == 0) || (y == 0) || (x == (imageWidth - 1)) ||
                (y == (imageHeight - 1))) {

                n++;

                if (n == numberOfPixels) {
                    return true;
                }
            }
        }

        return (n >= numberOfPixels);
    }

    public int countNeighbors(int x, int y, Set<PairInt> points, int imageWidth,
        int imageHeight) {

        int nn = 0;

        for (int i = 0; i < eightNeighborsX.length; i++) {

            int x2 = x + eightNeighborsX[i];
            int y2 = y + eightNeighborsY[i];

            if ((x2 < 0) || (x2 > (imageWidth - 1)) || (y2 < 0) ||
                (y2 > (imageHeight - 1))) {
                continue;
            }

            PairInt p2 = new PairInt(x2, y2);

            if (points.contains(p2)) {
                nn++;
            }
        }

        return nn;
    }

    public List<PairIntArray> smoothAndReExtractEdges(List<PairIntArray> edges,
        GreyscaleImage gradientXY, int smoothingFactor) {

        AverageUtil avgUtil = new AverageUtil();

        GreyscaleImage output = gradientXY.createWithDimensions();

        for (int i = 0; i < edges.size(); ++i) {
            PairIntArray edge = edges.get(i);
            if (edge.getN() >= smoothingFactor) {
                edge = avgUtil.calculateBoxCarAverage(edges.get(i), smoothingFactor);
                for (int j = 0; j < edge.getN(); ++j) {
                    output.setValue(edge.getX(j), edge.getY(j), 1);
                }
            }
        }

        PostLineThinnerCorrections pslt = new PostLineThinnerCorrections();
        pslt.correctForArtifacts(output);
        IEdgeExtractor edgeExtractor = new EdgeExtractorWithJunctions(output);
        edgeExtractor.removeShorterEdges(true);
        edges = edgeExtractor.findEdges();

        return edges;
    }

    public boolean hasAtLeastOneNonPointNeighbor(int x, int y,
        Set<PairInt> points, int imageWidth, int imageHeight) {

        for (int i = 0; i < eightNeighborsX.length; i++) {

            int x2 = x + eightNeighborsX[i];
            int y2 = y + eightNeighborsY[i];

            if ((x2 < 0) || (x2 > (imageWidth - 1)) || (y2 < 0) ||
                (y2 > (imageHeight - 1))) {
                continue;
            }

            if (!points.contains(new PairInt(x2, y2))) {
                return true;
            }
        }

        return false;
    }

    public void findNeighborsWithAtLeastOneNonPoint(int x, int y,
        Set<PairInt> outputNeighbors, Set<PairInt> points,
        Set<PairInt> excludePoints, int imageWidth, int imageHeight) {

        outputNeighbors.clear();

        for (int i = 0; i < eightNeighborsX.length; i++) {

            int x2 = x + eightNeighborsX[i];
            int y2 = y + eightNeighborsY[i];

            if ((x2 < 0) || (x2 > (imageWidth - 1)) || (y2 < 0) ||
                (y2 > (imageHeight - 1))) {
                continue;
            }

            boolean isPossiblyABorderPoint = hasAtLeastOneNonPointNeighbor(
                x2, y2, points, imageWidth, imageHeight);

            if (!isPossiblyABorderPoint) {
                continue;
            }

            PairInt p2 = new PairInt(x2, y2);

            if (excludePoints.contains(p2)) {
                continue;
            }
            if (points.contains(p2)) {
                outputNeighbors.add(p2);
            }
        }
    }

    public boolean isAdjacent(PairIntArray edge, int idx1, int idx2) {

        int x1 = edge.getX(idx1);
        int y1 = edge.getY(idx1);

        int x2 = edge.getX(idx2);
        int y2 = edge.getY(idx2);

        int diffX = Math.abs(x1 - x2);
        int diffY = Math.abs(y1 - y2);

        if ((diffX < 2) && (diffY < 2)) {
            return true;
        }

        return false;
    }

    /**
     * given theta and the point (xp, yp), determine which direction and hence
     * polar angle (clockwise) is perpendicular away from the centroid.
     * The reference point (xm, ym) is the point from which theta was also
     * calculated, which is probably the point for kMaxIdx.
     * @param theta
     * @param xp
     * @param yp
     * @param xm
     * @param ym
     * @param centroidXY
     * @return
     */
    public double calculatePerpendicularAngleAwayFromCentroid(
        double theta, int xp, int yp, int xm, int ym, double[] centroidXY) {

        /*
        rotate the point (xm, ym) around (xp, yp) 90 degrees and -90 degrees.
        The rotated point which is furthest from the centroid is the
        direction of the vector pointing away from the centroid.
        */

        /*
        math.cos(math.pi/2) = 0
        math.sin(math.pi/2) = 1
        math.sin(-math.pi/2) = -1

        double xr = centroidX + ((y - centroidY) * sine(angle)));
        double yr = centroidY + ((-(x - centroidX) * sine(angle)))
        */

        int xmRot90 = xp + (ym - yp);
        int ymRot90 = yp + (-(xm - xp));

        int xmRotNegative90 = xp  - (ym - yp);
        int ymRotNegative90 = yp + (xm - xp);

        double distSqRot90 = (xmRot90 - centroidXY[0]) * (xmRot90 - centroidXY[0])
            + (ymRot90 - centroidXY[1]) * (ymRot90 - centroidXY[1]);

        double distSqRotNegative90 =
            (xmRotNegative90 - centroidXY[0]) * (xmRotNegative90 - centroidXY[0])
            + (ymRotNegative90 - centroidXY[1]) * (ymRotNegative90 - centroidXY[1]);

        double perp = theta;

        if (distSqRot90 > distSqRotNegative90) {
            perp += Math.PI/2.;
        } else {
            perp -= Math.PI/2.;
        }

        if (perp >= 2*Math.PI) {
            perp = perp - 2*Math.PI;
        } else if (perp < 0) {
            perp += 2*Math.PI;
        }

        return perp;
    }

    /**
     * given theta and the point (xp, yp), determine which direction and hence
     * polar angle (clockwise) is perpendicular away from the centroid.
     * The reference point (xm, ym) is the point from which theta was also
     * calculated, which is probably the point for kMaxIdx.  The points are also
     * checked to make sure they aren't in the points set.
     *
     * @param theta
     * @param xp
     * @param yp
     * @param xm
     * @param ym
     * @param centroidXY
     * @param points
     * @return
     */
    public double calculatePerpendicularAngleAwayFromCentroid(
        double theta, int xp, int yp, int xm, int ym, double[] centroidXY,
        Set<PairInt> points) {

        /*
        rotate the point (xm, ym) around (xp, yp) 90 degrees and -90 degrees.
        The rotated point which is furthest from the centroid is the
        direction of the vector pointing away from the centroid.
        */

        /*
        math.cos(math.pi/2) = 0
        math.sin(math.pi/2) = 1
        math.sin(-math.pi/2) = -1

        double xr = centroidX + ((y - centroidY) * sine(angle)));
        double yr = centroidY + ((-(x - centroidX) * sine(angle)))
        */

        int xmRot90 = xp + (ym - yp);
        int ymRot90 = yp + (-(xm - xp));

        int xmRotNegative90 = xp  - (ym - yp);
        int ymRotNegative90 = yp + (xm - xp);

        boolean rot90IsInPoints = points.contains(
            new PairInt(Math.round(xmRot90), Math.round(ymRot90)));

        boolean rotNegative90IsInPoints = points.contains(
            new PairInt(Math.round(xmRotNegative90),
            Math.round(ymRotNegative90)));

        double distSqRot90 = (xmRot90 - centroidXY[0]) * (xmRot90 - centroidXY[0])
            + (ymRot90 - centroidXY[1]) * (ymRot90 - centroidXY[1]);

        double distSqRotNegative90 =
            (xmRotNegative90 - centroidXY[0]) * (xmRotNegative90 - centroidXY[0])
            + (ymRotNegative90 - centroidXY[1]) * (ymRotNegative90 - centroidXY[1]);

        double perp = theta;

        if (distSqRot90 > distSqRotNegative90) {
            perp += Math.PI/2.;
        } else if (distSqRot90 > distSqRotNegative90) {
            if (rot90IsInPoints && !rotNegative90IsInPoints) {
                perp -= Math.PI/2.;
            } else if (!rot90IsInPoints && rotNegative90IsInPoints) {
                perp += Math.PI/2.;
            } else {
                throw new IllegalStateException("Error in algorithm:" +
                " consider changing the test 90 and -90 points so that" +
                " one will always be in points set.");
            }
        } else {
            perp -= Math.PI/2.;
        }

        if (perp >= 2*Math.PI) {
            perp = perp - 2*Math.PI;
        } else if (perp < 0) {
            perp += 2*Math.PI;
        }

        return perp;
    }

    public int adjustEdgesTowardsBrighterPixels(PairIntArray edge,
        GreyscaleImage edgeGuideImage) {

        int nEdgeReplaced = 0;

        /*
        looking at the 8 neighbor region of each pixel for which the
        pixel's preceding and next edge pixels remain connected for
           and among those, looking for a higher intensity pixel than
           the center and if found, change coords to that.
        */
        for (int i = 1; i < (edge.getN() - 1); i++) {
            int x = edge.getX(i);
            int y = edge.getY(i);
            int prevX = edge.getX(i - 1);
            int prevY = edge.getY(i - 1);
            int nextX = edge.getX(i + 1);
            int nextY = edge.getY(i + 1);

            int maxValue = edgeGuideImage.getValue(x, y);
            int maxValueX = x;
            int maxValueY = y;
            boolean changed = false;

            for (int col = (prevX - 1); col <= (prevX + 1); col++) {

                if ((col < 0) || (col > (edgeGuideImage.getWidth() - 1))) {
                    continue;
                }

                for (int row = (prevY - 1); row <= (prevY + 1); row++) {

                    if ((row < 0) || (row > (edgeGuideImage.getHeight() - 1))) {
                        continue;
                    }
                    if ((col == prevX) && (row == prevY)) {
                        continue;
                    }
                    if ((col == nextX) && (row == nextY)) {
                        continue;
                    }

                    // skip if pixel is not next to (nextX, nextY)
                    int diffX = Math.abs(nextX - col);
                    int diffY = Math.abs(nextY - row);
                    if ((diffX > 1) || (diffY > 1)) {
                        continue;
                    }

                    if (edgeGuideImage.getValue(col, row) > maxValue) {
                        maxValue = edgeGuideImage.getValue(col, row);
                        maxValueX = col;
                        maxValueY = row;
                        changed = true;
                    }
                }
            }
            if (changed) {
                nEdgeReplaced++;
                edge.set(i, maxValueX, maxValueY);
            }
        }

        return nEdgeReplaced;
    }

    public int adjustEdgesTowardsDarkerPixels(PairIntArray edge,
        GreyscaleImage edgeGuideImage) {

        int nEdgeReplaced = 0;

        /*
        looking at the 8 neighbor region of each pixel for which the
        pixel's preceding and next edge pixels remain connected for
           and among those, looking for a lower intensity pixel than
           the center and if found, change coords to that.
        */
        for (int i = 1; i < (edge.getN() - 1); i++) {
            int x = edge.getX(i);
            int y = edge.getY(i);
            int prevX = edge.getX(i - 1);
            int prevY = edge.getY(i - 1);
            int nextX = edge.getX(i + 1);
            int nextY = edge.getY(i + 1);

            int minValue = edgeGuideImage.getValue(x, y);
            int minValueX = x;
            int minValueY = y;
            boolean changed = false;

            for (int col = (prevX - 1); col <= (prevX + 1); col++) {

                if ((col < 0) || (col > (edgeGuideImage.getWidth() - 1))) {
                    continue;
                }

                for (int row = (prevY - 1); row <= (prevY + 1); row++) {

                    if ((row < 0) || (row > (edgeGuideImage.getHeight() - 1))) {
                        continue;
                    }
                    if ((col == prevX) && (row == prevY)) {
                        continue;
                    }
                    if ((col == nextX) && (row == nextY)) {
                        continue;
                    }

                    // skip if pixel is not next to (nextX, nextY)
                    int diffX = Math.abs(nextX - col);
                    int diffY = Math.abs(nextY - row);
                    if ((diffX > 1) || (diffY > 1)) {
                        continue;
                    }

                    if (edgeGuideImage.getValue(col, row) < minValue) {
                        minValue = edgeGuideImage.getValue(col, row);
                        minValueX = col;
                        minValueY = row;
                        changed = true;
                    }
                }
            }
            if (changed) {
                nEdgeReplaced++;
                edge.set(i, minValueX, minValueY);
            }
        }

        return nEdgeReplaced;
    }

    /**
     * Find sections of the closed curve that are two pixels wide and separate
     * loops in the curve:
     * <pre>
     * for example:         #  #
     *    #  #  #         #     #
     *  #         #  #  #      #
     *   #   #  # #  #  #  #  # 
     * </pre>
     * These sections are thinned to width of '1' by the line thinner,
     * so need to be restored afterwards or prevented from being removed.
     * @param closedCurve
     * @return list of points that are part of the 2 pixel width patterns in
     * a curve where the curve closes, but is still connected.
     */
    public List<Set<PairInt>> findButterflySections(PairIntArray closedCurve) {
        
        /*
        endpoints for vert:
                       #           #
            # .      = - .      -  - .
          - - .  or  - - .  or     # .
            #          #
        
        endpoints for horiz:
                -        - -       -
              # - #    # - - #   # - #
              . .        . .       . .
        
        endpoints for diag:
        
            -  #             -  #               -  #
            -    .        -  -  .            #  -  .
            #  .   .      #  .    .          -  .    .
                 .              .                  .
        
        The sections of line which are 2 pixels wide and 1 further from the
        endpoint have 3 non-point neighbors each and 5 point set neighbors
        An area limit further constrains the geometry.
        For sections matching the patterns below, could consider storing
        the pattern for each pix as 'v', 'h', or 'd'...not an apparent use for 
        that yet though.
        
        Segment patterns between endpoints:
                       4
           -  -  -  -  3
           .  #  #  .  2
           .  #  #  .  1
           -  -  -  -  0
        0  1  2  3  4
        
                          4
           -  -  -  -  -  3
           .  #  #  #  .  2
           .  #  #  #  .  1
           -  -  -  -  -  0
        0  1  2  3  4  5
        
           # # - -   3
           - # # - - 2
           - - # #   1
           - - - - - 0
        0  1 2 3 4 5
        
        Scan the line,
           if a point fits one of the 4 segment patterns (4th is diag transformed by x=-x), 
           add it to a group and add the remaining pts fitting the pattern to a stack
           -- traverse the stack adding contiguous points to the group that fit the
              pattern.
           -- note where the first point in the pattern started, because
              when there are no more contiguous pattern matching points,
              the scan will continue at the next point after that first,
              but will skip those already added to a group.
          
        When the scan for groups has finished,
             for each group, need to apply the above endpoint patterns to see
             if the candidate segment is surrounded by 2 endpoints.
           
             test all candidate group points as adjacent to potential endpoints.
        
             When a match is found, have to exclude all of the matching pattern
             from the oppossing endpoint tests.
        
             This is the smallest pattern which will match that suggestion:
              - - - -
            # # @ # #  
          - - # @ # -
            # - - - #
             The '@'s are the candidate group points.  The #'s are points 
             matching endoint patterns.
            
             The found endpoints for one end, the left for example, would
             be excluded from a search for matching to the other endpoints.
     
        Note that this pattern and variants of it as very short sections and
        endpoints should be scanned after the above to find the shortest
        butterfly segments.
              - -
            # # @ #   
          - - # @ -
            # - - #     
        
        For each segment group which has 2 matching endpoints, those should 
        be stored as butterfly sections in a set.  Each one of those
        should be passed back in a list as the return of this method.
        
        runtime complexity is linear in the number of points in the given
        closed curve.
        */
        
        throw new UnsupportedOperationException("not yet implemented");
    }
    
    public boolean correctForXCrossings(PairIntArray closedCurve) throws Exception {
                
        // method will look for x-sections in the curve and then assert that the
        // point order in the curve does not cross the x-section, instead
        // enters on one side and later exits on the other side of the x-section.
        // points that violate that are corrected and counted.
        
        // at the end of the method, returns whether changes were made or not.
        
        /*
        examples:
            9
          0 1 8 7
              2
        
            9
            1 8 7
          0   2
        
        then pattern transformed by x and y times -1
        
        the method will make a Set of PairInts for every point.
        while nIter==0 or nChanges > 0
           scan each point in closedCurve looking for the x-crossing
           pattern and if found:
               check that points are ordered as expected
               and if not, change the point order and restart the loop with nIter++
        
        The check that points are ordered as expected needs to make corrections
        for wrap around for endpoints.
                   
        A violation of the crossing would be:
            9                   9
          0 1 2 3               1 2 3
              8               0   8
        
        The check that points are ordered as expected should probably 
        just try each of 4 correct patterns for point order and stop 
        with a "ordered=true" when found.
        */
        
        throw new UnsupportedOperationException("not yet implemented");
    }
}
