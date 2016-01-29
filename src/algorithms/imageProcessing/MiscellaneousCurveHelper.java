package algorithms.imageProcessing;

import algorithms.compGeometry.convexHull.GrahamScanPairInt;
import algorithms.compGeometry.convexHull.GrahamScanTooFewPointsException;
import algorithms.imageProcessing.features.CornerRegion;
import algorithms.imageProcessing.scaleSpace.CurvatureScaleSpaceContour;
import algorithms.imageProcessing.util.AngleUtil;
import algorithms.imageProcessing.util.PairIntWithIndex0;
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
     * determine whether the closed curve points are ordered in a counter clockwise
     * manner 
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

    /**
     * determine whether the closed curve points are ordered in a counter clockwise
     * manner by first computing the convex hull then
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
    public boolean curveIsOrderedClockwise2(PairIntArray closedCurve) {

        if (closedCurve.getN() < 2) {
            return false;
        } else if (closedCurve.getN() < 4) {
            return curveIsOrderedClockwise(closedCurve);
        }
        
        int n = closedCurve.getN();
        
        PairIntWithIndex0[] p = new PairIntWithIndex0[n];
        for (int i = 0; i < n; ++i) {
            p[i] = new PairIntWithIndex0(closedCurve.getX(i), closedCurve.getY(i),  i);
        }
        
        GrahamScanPairInt<PairIntWithIndex0> scan = new GrahamScanPairInt<PairIntWithIndex0>();
        try {
            scan.computeHull(p);
            
            // hull returns points in clockwise order
            
            n = scan.getHull().size() - 1;
            //PairIntArray hull = new PairIntArray(n);
            //List<Integer> hullCurveIndexes = new ArrayList<Integer>();
            //int[] deltaIndexes = new int[n];
            
            // nPos or nNeg might be 1 and then other n-2 if there is wrap-around
            int nNeg = 0;
            int nPos = 0;
            for (int i = 0; i < n; ++i) {
                
                PairIntWithIndex0 p0 = scan.getHull().get(i);
                
                //hull.add(Math.round(p0.getX()), Math.round(p0.getY()));
                //hullCurveIndexes.add(Integer.valueOf(p0.getPixIndex()));
                
                // for CW input, expect these to be + numbers
                int deltaIndex = scan.getHull().get(i + 1).getPixIndex() - p0.getPixIndex();
                if (deltaIndex > 0) {
                    nPos++;
                } else {
                    nNeg++;
                }
            }
            
            //boolean isCW = curveIsOrderedClockwise(hull);
            //assert(isCW);
            
            if (nPos > nNeg) {
                return true;
            }
            
            return false;
            
        } catch (GrahamScanTooFewPointsException ex) {
            return curveIsOrderedClockwise(closedCurve);
        }
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
    
    public <T extends CornerRegion> double[] calculateXYCentroids0(List<T> list) {
        
        double xc = 0;
        double yc = 0;

        for (CornerRegion cr : list) {
            double x = cr.getX()[cr.getKMaxIdx()];
            double y = cr.getY()[cr.getKMaxIdx()];
            xc += x;
            yc += y;
        }
        xc /= (double)list.size();
        yc /= (double)list.size();

        return new double[]{xc, yc};
    }
    
    public double[] calculateXYCentroids1(List<CurvatureScaleSpaceContour> list) {
        
        double xc = 0;
        double yc = 0;

        for (CurvatureScaleSpaceContour cr : list) {
            double x = cr.getPeakDetails()[0].getXCoord();
            double y = cr.getPeakDetails()[0].getYCoord();
            xc += x;
            yc += y;
        }
        xc /= (double)list.size();
        yc /= (double)list.size();

        return new double[]{xc, yc};
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

    public Set<PairInt> findNeighbors(int x, int y, Set<PairInt> points) {

        Set<PairInt> neighbors = new HashSet<PairInt>();

        for (int i = 0; i < eightNeighborsX.length; i++) {

            int x2 = x + eightNeighborsX[i];
            int y2 = y + eightNeighborsY[i];

            PairInt p2 = new PairInt(x2, y2);

            if (points.contains(p2)) {
                neighbors.add(p2);
            }
        }

        return neighbors;
    }
        
    public void findNeighbors(int x, int y, Set<PairInt> points, 
        Set<PairInt> excludePoints, int[] dxs, int[] dys, 
        Set<PairInt> outputNeighbors) {
        
        outputNeighbors.clear();
        
        for (int i = 0; i < dxs.length; i++) {
            
            int x2 = x + dxs[i];
            int y2 = y + dys[i];
            
            PairInt p2 = new PairInt(x2, y2);
            
            if (excludePoints.contains(p2)) {
                continue;
            }
            if (points.contains(p2)) {
                outputNeighbors.add(p2);
            }
        }
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
    
    public int countNeighbors(int x, int y, Set<PairInt> points) {

        int nn = 0;

        for (int i = 0; i < eightNeighborsX.length; i++) {

            int x2 = x + eightNeighborsX[i];
            int y2 = y + eightNeighborsY[i];

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
     * given 3 counter-clockwise ordered points on a curve, calculate the angle 
     * along the curve at the middle point, its direction is from p0 to p1.
     * <pre>
     * For example:
     * 
     * 135 degrees
     *       .---
     *       | .
     *           p2   
     *             p1
     *                p0
     * </pre>
     * @param x1
     * @param y1
     * @param x2
     * @param y2
     * @param x3
     * @param y3
     * @return
     */
    public double calculateAngleAtMidpoint(int x1, int y1, 
        int x2, int y2, int x3, int y3) {

        /*
        given the points (x1, y1) (x2, y2) and (x3, y3), 
        calculates the angle at the midpoint (x2, y2) for the path along
        the points.
        */
        
        double theta1 = AngleUtil.polarAngleCCW(x2 - x1, y2 - y1);
        
        double theta2 = AngleUtil.polarAngleCCW(x3 - x2, y3 - y2);
        
        double theta = AngleUtil.getAngleAverageInRadians(theta1, theta2);
                
        return theta;
    }
    
    /**
     * given 3 counter-clockwise ordered points on a curve, calculate the angle 
     * tangent to the curve at the middle point - its direction follows
     * the right hand rule.
     * <pre>
     * For example:
     *                  45 degrees
     *               __
     *               . |
     *       p2    .
     *          p1
     *             p0
     * </pre>
     * @param x1
     * @param y1
     * @param x2
     * @param y2
     * @param x3
     * @param y3
     * @return
     */
    public double calculateAngleTangentToMidpoint(int x1, int y1, 
        int x2, int y2, int x3, int y3) {

        double theta = calculateAngleAtMidpoint(x1, y1, x2, y2, x3, y3);
               
        double thetaMinus90 = theta - Math.PI/2.;
        if (thetaMinus90 < 0) {
            thetaMinus90 += (2.*Math.PI);
        }
        
        return thetaMinus90;
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

    /**
     * calculate theta in degrees as a value between -pi and pi for the given point
     * which should have a value greater than 0 at (x,y) and should be part
     * of a curve thinned to a width of 1.
     <pre>
                 90
           135    |    45
                  |
        180 ---------------  0
                  |
          -135    |   -45
                 -90
     </pre>
     * @param x
     * @param y
     * @param img
     * @return 
     */
    public int calculateThetaForPointOnEdge(int x, int y, GreyscaleImage img) {
        
        double[] gXY = calculateGradientsForPointOnEdge(x, y, img);
        
        /*
        Math.atan arc tangent, angle is in the range -pi/2 through pi/2
        Math.atan2 conversion of rectangular coordinates (x, y) 
            to polar coordinates (r, theta). 
            This method computes the phase theta by computing an arc tangent 
            of y/x in the range of -pi to pi.
        */
        double t = Math.atan2(gXY[1], gXY[0]);
      
        return (int)Math.round(t);
    }

    /**
     * calculate gradient x and gradient y for the given point
     * which should have a value greater than 0 at (x,y) and should be part
     * of a curve thinned to a width of 1.
     
     * @param x
     * @param y
     * @param img
     * @return double{gradX, gradY}
     */
    public double[] calculateGradientsForPointOnEdge(int x, int y, GreyscaleImage img) {
        
        Kernel1DHelper kernel1DHelper = new Kernel1DHelper();
        
        float sigma = 0.42466090014400953f;
        
        float[] kernel = Gaussian1D.getKernel(sigma);

        float[] kernel2 = Gaussian1D.getKernel(sigma * 1.6f);
        
        boolean calcForX = true;
               
        double convX1 = kernel1DHelper.convolvePointWithKernel(
            img, x, y, kernel, calcForX);
        
        double convX2 = kernel1DHelper.convolvePointWithKernel(
            img, x, y, kernel2, calcForX);
             
        calcForX = false;
        
        double convY1 = kernel1DHelper.convolvePointWithKernel(
            img, x, y, kernel, calcForX);
        
        double convY2 = kernel1DHelper.convolvePointWithKernel(
            img, x, y, kernel2, calcForX);
        
        double gX = convX2 - convX1;
        double gY = convY2 - convY1;
        
        return new double[]{gX, gY};
    }
    
    /**
     * calculate gradient x and gradient y for the given point
     * which should have a value greater than 0 at (x,y) and should be part
     * of a curve thinned to a width of 1.
     * Note that the magnitudes have not been calibrated because the main
     * using method uses the results to calculate the polar angle, so 
     * factor applied to both not necessary.
       Note that the magnitudes have not been calibrated because the main
     * using method uses the results to calculate the polar angle, so 
     * factor applied to both not necessary.
     * @param x
     * @param y
     * @param points
     * @return double{gradX, gradY}
     */
    public double[] calculateGradientsForPointOnEdge(int x, int y, 
        Set<PairInt> points) {
        
        Kernel1DHelper kernel1DHelper = new Kernel1DHelper();
        
        float sigma = 0.42466090014400953f;
        
        float[] kernel = Gaussian1D.getKernel(sigma);

        float[] kernel2 = Gaussian1D.getKernel(sigma * 1.6f);
        
        boolean calcForX = true;
        
        double convX1 = kernel1DHelper.convolvePointWithKernel(
            points, x, y, kernel, calcForX);
        
        double convX2 = kernel1DHelper.convolvePointWithKernel(
            points, x, y, kernel2, calcForX);
             
        calcForX = false;
        
        double convY1 = kernel1DHelper.convolvePointWithKernel(
            points, x, y, kernel, calcForX);
        
        double convY2 = kernel1DHelper.convolvePointWithKernel(
            points, x, y, kernel2, calcForX);
             
        double gX = convX2 - convX1;
        double gY = convY2 - convY1;
        
        return new double[]{gX, gY};
    }

    public double calculateArea(PairIntArray closedCurve) {
        
        int n = closedCurve.getN();
        
        double sum = 0;
        
        for (int i = 0; i < (n - 1); ++i) {
            
            double t = 0.5 * (closedCurve.getY(i + 1) + closedCurve.getY(i)) *
                (closedCurve.getX(i + 1) - closedCurve.getX(i));
            
            sum += t;
        }
        
        sum += ((closedCurve.getY(0) + closedCurve.getY(n - 1)) *
                (closedCurve.getX(0) - closedCurve.getX(n -1)));
        
        return sum;
    }

}
