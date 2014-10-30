package algorithms.imageProcessing;

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

}
