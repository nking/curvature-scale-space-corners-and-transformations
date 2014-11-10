package algorithms.imageProcessing;

import Jama.Matrix;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Stack;
import java.util.logging.Logger;

/**
 * 
 * @author nichole
 */
public class MiscellaneousCurveHelper {
    
    private Logger log = Logger.getLogger(this.getClass().getName());
    
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
      
        if (closedCurve.getN() < 3) {
            return true;
        }
        
        int nNeg = 0;
        int n = closedCurve.getN();
        
        for (int i = 1; i < (n - 1); i++) {
            
            //(xi - xi-1) * (yi+1 - yi) - (yi - yi-1) * (xi+1 - xi)
            long crossProduct = ((closedCurve.getX(i) - closedCurve.getX(i - 1))
                * (closedCurve.getY(i + 1) - closedCurve.getY(i)))
                - ((closedCurve.getY(i) - closedCurve.getY(i - 1))*
                (closedCurve.getX(i + 1) - closedCurve.getX(i)));
            
            if (crossProduct < 0) {
                nNeg++;
            }
        }
        
        Logger.getLogger(this.getClass().getName()).fine(
            "nNeg=" + nNeg + " nTot=" + n);
        
        Logger.getLogger(this.getClass().getName()).fine(
            String.format("(%d,%d)", closedCurve.getX(0), closedCurve.getY(0))
            + String.format(" (%d,%d)", closedCurve.getX(1), closedCurve.getY(1))
            + String.format(" (%d,%d)", closedCurve.getX(2), closedCurve.getY(2))
            + " ... "
            + String.format(" (%d,%d)", closedCurve.getX(n - 2), closedCurve.getY(n-2))
            + String.format(" (%d,%d)", closedCurve.getX(n - 1), closedCurve.getY(n-1))
        );
        
        int nPos = n - 2 - nNeg;
        
        // note, may want to adjust this for image perspective where
        // positive y is in downward direction.
        return (nNeg >= nPos);
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
        int[] crossCorrelationOffset = new int[1];
          
        /*
         returns whether the curve 'check' is adjacent to the curve 'node0',
         and if so, returns the offset in the frame of the larger curve.
         the offset represents where the first point in the shorter curve
         matches in the larger curve.                
         */
        boolean isAdjacent = crossCorrelation(longer, shorter,
            crossCorrelationOffset);

        if (isAdjacent) {

            /*
                -- find any points in check outside of the overlap
                   and add those to the larger node.
            */
            if (crossCorrelationOffset[0] < 0) {
                // add from the beginning of shorter if any are unmatched
                int nInsert = -1*crossCorrelationOffset[0];
                longer.insertSpaceAtTopOfArrays(nInsert);

                for (int ii = 0; ii < nInsert; ii++) {
                    longer.set(ii, shorter.getX(ii), shorter.getY(ii));
                }
            } else {
                //add from end of shorter if any are unmatched
                int n0 = longer.getN() - crossCorrelationOffset[0];
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
     * return true if cross-correlation shows that the 2 curves are adjacent
     * to one another.  Note that the method needs the points within the
     * curves to be ordered in a similar manner and for the endpoints of the
     * curves to be accurate.  If a point in the middle of the curve is 
     * the first or last point, it may prevent comparison of it with another
     * edge's endpoints.
     * 
     * @param curve0
     * @param curve1
     * @param crossCorrelationOffset offset of where the shorter curve starts
     *  with respect to the longer.  For example, an offset of -2 means that
     * the first 2 points in the shorter curve are outside of the longer curve,
     * but the next point in the longer curve is adjacent to the shorter.
     * Another example: if offset is +2, the first pixel in the shorter curve
     * is adjacent to the third pixel in the longer curve.  NOTE: the offset
     * is only useful if this method returns true;
     * @return 
     */
    protected boolean crossCorrelation(PairIntArray curve0, PairIntArray curve1, 
        int[] crossCorrelationOffset) {
        
        crossCorrelationOffset[0] = Integer.MAX_VALUE;
        
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
        offset in crossCorrelationOffset and return true, else false        
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

            crossCorrelationOffset[0] = cSSMinOffset;

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
    public double[] calculateXYCentroids(Matrix xy) {
        
        double xc = 0;
        double yc = 0;
        
        int n = xy.getArray()[0].length;
        
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
    
    public double[] calculateXYCentroids(SearchableCurve xy) {
        
        double xc = 0;
        double yc = 0;
        
        for (int i = 0; i < xy.getN(); i++) {
            
            xc += xy.getX()[i];
            
            yc += xy.getY()[i];
        }
        
        xc /= (double)xy.getN();
        
        yc /= (double)xy.getN();
        
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
        
        log.info("removeRedundantPoints");
        
        // if there are redundant points, remove the points in between
        
        //TODO: replace w/ faster algorithm...
        
        for (int i = 0; i < tmpEdges.size(); i++) {
            
            List<String> points = new ArrayList<String>();
            
            PairIntArray edge = tmpEdges.get(i);
            
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
    }

    public void pruneAdjacentNeighborsTo2(List<PairIntArray> tmpEdges) {
       
        log.info("pruneAdjacentNeighborsTo2");
        
        // this will usually only have 2 in it, and most expected is 3
        int[] outputAdjacentNeighbors = new int[8];
        
        for (int lIdx = 0; lIdx < tmpEdges.size(); lIdx++) {
            
            // quick check for whether an edged has 3 neighbors, then
            // compare with patterns
            
            PairIntArray edge = tmpEdges.get(lIdx);
            
            for (int eIdx = 0; eIdx < edge.getN(); eIdx++) {
                
                int nNeighbors = getAdjacentNeighbors(edge, eIdx,
                    outputAdjacentNeighbors);
                
                if ((nNeighbors > 2)) {
                    
                    boolean pruned = pruneAdjacentNeighborsTo2(edge, eIdx, 
                        outputAdjacentNeighbors, nNeighbors);
   
                    if (pruned) {
                        // restart iteration for easier maintenance
                        eIdx = -1;
                        //237,201  edge0
                        log.info("removed a point from edge=" + lIdx);
                    }
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
        
        log.info("removing point (" + edge.getX(maxDistanceIdx) 
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
        
        int[] xs = new int[2];
        int[] ys = new int[2];
        
        for (int lIdx = 0; lIdx < tmpEdges.size(); lIdx++) {
            
            PairIntArray edge = tmpEdges.get(lIdx);
            
            for (int i = 0; i < edge.getN(); i++) {
                
                correctCheckeredSegments(edge, i, xs, ys);
            }
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
     * @param curve
     * @return 
     */
    public PairIntArray findJaggedLineSegments(final PairIntArray curve) {
        
        //TODO: use minimum curve size
        if (curve == null || curve.getN() < 5) {
            return null;
        }

        PairIntArray ledges = findLedgesInCurve(curve);
        
        PairIntArray remainingRanges = 
            writeRangesNotAlreadyIncluded(curve, ledges);
        
        for (int i = 0; i < remainingRanges.getN(); i++) {
            
            int r0 = remainingRanges.getX(i);
            int r1 = remainingRanges.getY(i);
            
            PairIntArray staircaseRanges = 
                findJaggedLineStaircaseSegments(curve, r0, r1);
            
            if (staircaseRanges != null) {
                for (int j = 0; j < staircaseRanges.getN(); j++) {
                    int s0 = staircaseRanges.getX(j);
                    int s1 = staircaseRanges.getY(j);
                    ledges.add(s0, s1);
                }
            }
        }
        
        sortByX(ledges);
        
        // merge adjacent ranges
        mergeRanges(curve, ledges);
       
        // search for 45 degree lines 
        remainingRanges = 
            writeRangesNotAlreadyIncluded(curve, ledges);
        for (int i = 0; i < remainingRanges.getN(); i++) {            
            int r0 = remainingRanges.getX(i);
            int r1 = remainingRanges.getY(i);
            PairIntArray lineRanges = 
                find45DegreeSegments(curve, r0, r1);
            if (lineRanges != null) {
                for (int j = 0; j < lineRanges.getN(); j++) {
                    int s0 = lineRanges.getX(j);
                    int s1 = lineRanges.getY(j);
                    ledges.add(s0, s1);
                }
            }
        }
        
        sortByX(ledges);
        
        // merge ranges again
        mergeRanges(curve, ledges);
        
        return ledges;
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
                            log.info("       " + s0 + " : " + s1);
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
            return null;
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
            dx = (int) (curve.getX(i) - curve.getX(i - 1));
            dy = (int) (curve.getY(i) - curve.getY(i - 1));
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
                dx = (int) (curve.getX(i) - curve.getX(i - 1));
                dy = (int) (curve.getY(i) - curve.getY(i - 1));
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
                    
            int x = (int)curve.getX(i);
            int y = (int)curve.getY(i);
            
            int diffX = (int)(x - curve.getX(i - 1));
            
            int diffY = (int)(y - curve.getY(i - 1));
            
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
                    dx = (int) (curve.getX(i) - curve.getX(i - 1));
                    dy = (int) (curve.getY(i) - curve.getY(i - 1));
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
                        dx = (int) (curve.getX(i) - curve.getX(i - 1));
                        dy = (int) (curve.getY(i) - curve.getY(i - 1));
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
                        diffX = (int)(curve.getX(j) - curve.getX(j - 1));
                        diffY = (int)(curve.getY(j) - curve.getY(j - 1));
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
            return null;
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
            dx = (int) (curve.getX(i) - curve.getX(i - 1));
            dy = (int) (curve.getY(i) - curve.getY(i - 1));
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
                    
            int x = (int)curve.getX(i);
            int y = (int)curve.getY(i);
            
            int diffX = (int)(x - curve.getX(i - 1));
            
            int diffY = (int)(y - curve.getY(i - 1));
          
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
                    dx = (int) (curve.getX(i) - curve.getX(i - 1));
                    dy = (int) (curve.getY(i) - curve.getY(i - 1));
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
       
            int diffX = (int)(curve.getX(i) - curve.getX(i - 1));
            
            int diffY = (int)(curve.getY(i) - curve.getY(i - 1));

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
        
        // choosing a minimum size empirically from looking at edges in tests
        int minLedgeWidth = 4;
        
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
            dx = (int) (curve.getX(i) - curve.getX(i - 1));
            dy = (int) (curve.getY(i) - curve.getY(i - 1));
        }

        if (dx == 0) {
            runIsAlongX = Boolean.FALSE;
        } else if (dy == 0) {
            runIsAlongX = Boolean.TRUE;
        }

        int lineStart = i - 1;
                
        PairIntArray tmp = new PairIntArray();
        Boolean tmpRunIsAlongX = null;
        int tmpDX = -1;
        int tmpDY = -1;
        
        for (i = (lineStart + 1); i <= stop; i++) {
           
            int x = curve.getX(i);
            int y = curve.getY(i);
            
            int diffX = (int) (x - curve.getX(i - 1));

            int diffY = (int) (y - curve.getY(i - 1));

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
                // choosing a minimum size of 6 from looking at edges in tests
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
                    int tmpI = i;
                    while (!((dx == 0) && (dy != 0)) &&
                        !((dy == 0) && (dx != 0))) {
                        i++;
                        if (i >= stop) {
                            break;
                        }
                        dx = (int) (curve.getX(i) - curve.getX(i - 1));
                        dy = (int) (curve.getY(i) - curve.getY(i - 1));
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
                            diffX = (int) (curve.getX(j) - curve.getX(j - 1));
                            diffY = (int) (curve.getY(j) - curve.getY(j - 1));
                            
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
    
    public boolean isWithinARange(PairIntArray lineRangeSegments, int idx) {
        
        if (lineRangeSegments == null || lineRangeSegments.getN() == 0) {
            return false;
        }
        
        return isWithinARange(lineRangeSegments, idx, 0);
    }

    /**
     * search for idx within ranges in lineRangeSegments and return the index of
     * lineRangeSegments in which it is found, else -1.  Note that lineRangeSegments
     * have to be ordered by x and unique.
     * @param lineRangeSegments
     * @param idx
     * @param minDistFromEnds
     * @return 
     */
    public boolean isWithinARange(PairIntArray lineRangeSegments, int idx,
        int minDistFromEnds) {
        
        if (lineRangeSegments == null || lineRangeSegments.getN() == 0) {
            return false;
        }
        
        // search outside of bounds first:
        if (idx < lineRangeSegments.getX(0)) {
            return false;
        } else if (idx > lineRangeSegments.getY(lineRangeSegments.getN() - 1)) {
            return false;
        }
        
        for (int i = 0; i < lineRangeSegments.getN(); i++) {
            
            int idx0 = lineRangeSegments.getX(i);
            int idx1 = lineRangeSegments.getY(i);
            
            if ((idx >= (idx0 + minDistFromEnds)) 
                && (idx <= (idx1 - minDistFromEnds))) {
                
                return true;
            }
        }
        
        return false;
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

}
