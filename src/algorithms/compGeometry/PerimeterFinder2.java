package algorithms.compGeometry;

import algorithms.imageProcessing.MiscellaneousCurveHelper;
import algorithms.imageProcessing.util.AngleUtil;
import algorithms.misc.Misc;
import algorithms.misc.MiscMath;
import algorithms.search.NearestNeighbor2D;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.PolygonAndPointPlotter;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

/**
 * class to hold some of the newer methods
 * w.r.t. boundary point extractions.
 * 
 * @author nichole
 */
public class PerimeterFinder2 {
   
    /**
     * finds the pixels with neighbors not in given point
     * set, contiguousPoints.  Note, the contiguous points
     * fill out the shape for which the border is found.
     * 
     * @param contiguousPoints
     * @return 
     */
    public Set<PairInt> extractBorder(Set<PairInt> contiguousPoints) {

        Set<PairInt> border = new HashSet<PairInt>();
        
        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;
        
        for (PairInt p : contiguousPoints) {
            int x = p.getX();
            int y = p.getY();
            for (int i = 0; i < dxs.length; ++i) {
                int x2 = x + dxs[i];
                int y2 = y + dys[i];
                PairInt p2 = new PairInt(x2, y2);
                if (!contiguousPoints.contains(p2)) {
                    border.add(p);
                    break;
                }
            }
        }
        
        return border;
    }

    /**
     * NOT READY FOR USE...needs alot more testing.
     * 
     * given a contiguous set of points, extract the border and 
     * order the points.  NOTE: it is up to the invoker to 
     * give the method a point set which is contiguous.
     * The contiguous points
     * fill out the shape for which the border is found and they're
     * used to calculate medial axes and calculate the intersection
     * of adjacent non-shape neighbors.
     * 
     * @param contiguousPoints
     * @return 
     */
    public PairIntArray extractOrderedBorder(Set<PairInt> 
        contiguousPoints) {
        
        //O(8*N)
        Set<PairInt> boundary = extractBorder(contiguousPoints);
        
        // O(N*log_2(N))
        MedialAxis medAxis = new MedialAxis(contiguousPoints, 
            boundary);
        medAxis.fastFindMedialAxis();
        Set<PairInt> medAxisPts = medAxis.getMedialAxisPoints();
        
        return extractOrderedBorder(boundary, medAxisPts,
            contiguousPoints);
    }

    /**
     * NOT READY FOR USE...needs alot more testing.
     * 
     * given the border of a contiguous set of points and
     * given the medial axis points of that same set of
     * points (that is, the medial axis of the shape filling 
     * points within the border), return a clockwise set of 
     * ordered border points.
     * Note the borders having extended single pixel width
     * regions will have gaps in the returned ordered points.
     * 
     * @param borderPoints
     * @param medialAxisPoints
     * @param contiguousShapePoints
     * @return 
     */
    public PairIntArray extractOrderedBorder(Set<PairInt> 
        borderPoints, Set<PairInt> medialAxisPoints,
        Set<PairInt> contiguousShapePoints) {
        
        /*
        -- note points which only have 1 neighbor.  these
           are the single pixel wide spurs that result in
           gaps in the output ordered point segment
           -- the spur starts where the number of neighbors
              is 3 or more.
              that's where the curve continues after reaching
              the end of the spur.
              --> can keep track of these points and junctions
              with a linked list having key=index.
        
        starting w/ smallest x, smallest y,
        using DFS traversal through points and extracting points
        as they are added to the ordered points.
        - when there is more than one adjacent choice for next
          step, have to use the medial axis to determine
          the best choice:        
        */

        // constructed only if needed
        NearestNeighbor2D nn = null;   
        
        PairIntArray output = new PairIntArray(borderPoints.size());
        
        Set<PairInt> remaining = new HashSet<PairInt>(borderPoints);
        
        PairInt pt1 = findMinXY(borderPoints);

        output.add(pt1.getX(), pt1.getY());        
        remaining.remove(pt1);
        
        int[] neighborsX = new int[8];
        int[] neighborsY = new int[8];
        double[] angle = new double[8];
        
        // when still have points in remaining, and
        // have no neighbors in remaining for a point,
        // use this to backtrack to previous junction.
        LinkedList<Integer> junctionNodes = new LinkedList<Integer>();
        
        // because the juntionNodes may be inserted correctly, 
        // especially at the end of a complex perimeter
        // with remaining points not yet added,
        // keeping a separate list of junctionNodes
        // that are added to output while backtracking in order
        // to revise them later if needed.
        LinkedHashSet<Integer> jassocInserted = 
            new LinkedHashSet<Integer>();
        
        while (!remaining.isEmpty()) {
            int prevIdx = output.getN() - 1;
            int x = output.getX(prevIdx);
            int y = output.getY(prevIdx);
            
            // if there more than one adjacent neighbor
            // in remaining set, use the medial axis
            // and angles to determine next point
        
            int ns = findNeighbors(x, y, 
                remaining, neighborsX, neighborsY);
            
            if (ns == 1) {
                output.add(neighborsX[0], neighborsY[0]);
            } else if (ns == 0) {
                int ns2 = 0;
                while (ns2 == 0) {                    
                    // TODO: revisit this for complex shapes
                    
                    Integer index = junctionNodes.pollLast();
                    if (index == null) {
                        throw new IllegalStateException("Error in algorithm:"
                            + " no adjacent points for (" + x + " " + y
                            + ") but remaining is not empty and "
                            + " junction list is empty");
                    }
                    int x3 = output.getX(index.intValue());
                    int y3 = output.getY(index.intValue());

                    ns2 = findNeighbors(x3, y3, remaining,
                        neighborsX, neighborsY);
                
                    if (ns2 == 1) {
                        output.add(neighborsX[0], neighborsY[0]);
                    } else if (ns2 > 1) {
                        // find smallest angle subtended
                        // and add that point to output,
                        // then add output.size to junctionNodes
                        
                        if (nn == null) {
                            int[] minMaxY = MiscMath.findMinMaxXY(
                                borderPoints);
                            nn = new NearestNeighbor2D(medialAxisPoints,
                                minMaxY[1], minMaxY[3]);
                        }
     
                        int minAngleIdx = calculateMinAngles( 
                            x3, y3, ns2, neighborsX, neighborsY, nn,
                            contiguousShapePoints);
                        
                        junctionNodes.add(Integer.valueOf(output.getN()));
                        
                        output.add(neighborsX[minAngleIdx], 
                            neighborsY[minAngleIdx]);                        
                    }
                    jassocInserted.add(output.getN() - 1);
                }
            } else if (ns > 1) {
                // find smallest angle subtended
                // and add that point to output,
                // then add output.size to junctionNodes
                
                if (nn == null) {
                    int[] minMaxY = MiscMath.findMinMaxXY(borderPoints);
                    nn = new NearestNeighbor2D(medialAxisPoints, 
                        minMaxY[1], minMaxY[3]);
                }
          
                int minAngleIdx = calculateMinAngles(
                    x, y, ns, neighborsX, neighborsY, nn,
                    contiguousShapePoints);
                
                junctionNodes.add(Integer.valueOf(output.getN()));
             
//if (minAngleIdx == -1) {
    try {
        
    int[] xPolygon = null;
    int[] yPolygon = null;
    PolygonAndPointPlotter plotter = new
        PolygonAndPointPlotter();
    int[] xminmxyminmiac = MiscMath.findMinMaxXY(
         contiguousShapePoints);
    int[] xp, yp;
    int n, count;
    
    n = contiguousShapePoints.size();
    xp = new int[n];
    yp = new int[n];
    count = 0;
    for (PairInt p : contiguousShapePoints) {
        xp[count] = p.getX();
        yp[count] = p.getY();
        count++;
    }
    plotter.addPlot(xminmxyminmiac[0], xminmxyminmiac[1],
        xminmxyminmiac[2], xminmxyminmiac[3],
        xp, yp, xPolygon, yPolygon, "shape");
    plotter.writeFile2();
    
    n = medialAxisPoints.size();
    xp = new int[n];
    yp = new int[n];
    count = 0;
    for (PairInt p : medialAxisPoints) {
        xp[count] = p.getX();
        yp[count] = p.getY();
        count++;
    }    
    plotter.addPlot((float)xminmxyminmiac[0], 
        (float)xminmxyminmiac[1],
        (float)xminmxyminmiac[2], (float)xminmxyminmiac[3],
        xp, yp, xPolygon, yPolygon, "med axis");
    
    n = output.getN();
    xp = new int[n];
    yp = new int[n];
    for (int i = 0; i < output.getN(); ++i) {
        xp[i] = output.getX(i);
        yp[i] = output.getY(i);
    }
    plotter.addPlot(xminmxyminmiac[0], xminmxyminmiac[1],
        xminmxyminmiac[2], xminmxyminmiac[3],
        xp, yp, xPolygon, yPolygon, "ordered bunds");
    plotter.writeFile2();
    
    System.out.println("output=" + output.toString());
    } catch (Throwable t) {
        
    }
//}   
                
                output.add(neighborsX[minAngleIdx], neighborsY[minAngleIdx]);                
            }
            
            boolean rmvd = remaining.remove(
                new PairInt(output.getX(output.getN() - 1),
                output.getY(output.getN() - 1)));
            assert(rmvd);
        }
        
        if (!jassocInserted.isEmpty()) {
            /*System.out.println("re-check these:");
            for (Integer index : jassocInserted) {
                int idx = index.intValue();
                System.out.println("x=" + output.getX(idx) +
                    "," + output.getY(idx));
            }*/
            List<Integer> list = new ArrayList<Integer>(jassocInserted);
            Collections.sort(list);
            int[] xAdd = new int[list.size()];
            int[] yAdd = new int[list.size()];
            for (int i = (list.size() - 1); i > -1; --i) {
                int rmIdx = list.get(i).intValue();
                xAdd[i] = output.getX(rmIdx);
                yAdd[i] = output.getY(rmIdx);
                output.removeRange(rmIdx, rmIdx);
            }
            // find best place to insert them starting from end
            for (int i = (list.size() - 1); i > -1; --i) {
                int xp = xAdd[i];
                int yp = yAdd[i];
                double minDist = Double.MAX_VALUE;
                int minDistIdx = -1;
                for (int j = (output.getN() - 1); j > -1; --j) {
                    int x2 = output.getX(j);
                    int y2 = output.getY(j);
                    double d = Math.sqrt(distSq(xp, yp, x2, y2));
                    if (d < minDist) {
                        minDist = d;
                        minDistIdx = j;
                    } else if (d >= minDist && (minDist < 2)) {
                        break;
                    }
                }
                // default insert is minDistIdx + 1 unless
                // the existing distance is 1, then instead, use minDistIdx
                if (minDistIdx == (output.getN() - 1)) {
                    output.add(xp, yp);
                } else {
                    if (((minDistIdx - 1) > -1) 
                        && ((minDistIdx + 1) < output.getN())) {
                        
                        double dPrev0 = Math.sqrt(distSq(
                        output.getX(minDistIdx - 1), 
                        output.getY(minDistIdx - 1),
                        xp, yp));
                                                
                        double dNext0 = Math.sqrt(distSq(
                        output.getX(minDistIdx + 1), 
                        output.getY(minDistIdx + 1),
                        output.getX(minDistIdx), 
                        output.getY(minDistIdx)));
                        
                        double d0 = dPrev0 + minDist + dNext0;
                        
                        double dPrev1 = Math.sqrt(distSq(
                        output.getX(minDistIdx - 1), 
                        output.getY(minDistIdx - 1),
                        output.getX(minDistIdx), 
                        output.getY(minDistIdx)));
                        
                        double dNext1 = Math.sqrt(distSq(
                        output.getX(minDistIdx + 1), 
                        output.getY(minDistIdx + 1), xp, yp));
                        
                        double d1 = dPrev1 + minDist + dNext1;
                       
                        if ((d1 < d0) && (dPrev1 == minDist)) {
                            minDistIdx++;
                        }
                    }
                    output.insert(minDistIdx, xp, yp);
                }
            }
            // if any adjacent have dist > 2, see if swapping would
            // reduce that
            for (int i = 0; i < (output.getN() - 1); ++i) {
                int x = output.getX(i);
                int y = output.getY(i);
                int nextX = output.getX(i + 1);
                int nextY = output.getY(i + 1);
                double d1 = Math.sqrt(distSq(nextX, nextY, x, y));
                if (d1 >= 2 && (i > 0)) {
                    /*
                    prev  
                    x      next
                    next   x
                           nextnext
                    */
                    double dpn = Math.sqrt(distSq(nextX, nextY, 
                        output.getX(i - 1), output.getY(i - 1)));
                    if (dpn <= 2) {
                        
                    }
                }
            }
        }
        
        return output;
    }

    /**
     * find neighbors of point (x,y) present in set "remaining"
     * and place them in neighborsX and neighborsY and return
     * the number of them.
     * @param x
     * @param y
     * @param remaining
     * @param neighborsX
     * @param neighborsY
     * @return 
     */    
    private int findNeighbors(int x, int y,
        Set<PairInt> remaining, int[] neighborsX,
        int[] neighborsY) {
        
        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;
        
        int ns = 0;
        for (int k = 0; k < dxs.length; ++k) {
            int x2 = x + dxs[k];
            int y2 = y + dys[k];
            PairInt p2 = new PairInt(x2, y2);
            if (!remaining.contains(p2)) {
                continue;
            }
            neighborsX[ns] = x2;
            neighborsY[ns] = y2;
            ++ns;
        }
        
        return ns;
    }

    /**
     * find the smallest x point and the smallest y among those.
     * @param set
     * @return 
     */
    protected PairInt findMinXY(Set<PairInt> set) {

        int minX = Integer.MAX_VALUE;
        int minXsY = Integer.MAX_VALUE;        
        for (PairInt p : set) {
            if (p.getX() < minX) {
                minX = p.getX();
                minXsY = p.getY();
            } else if (p.getX() == minX) {
                if (p.getY() < minXsY) {
                    minX = p.getX();
                    minXsY = p.getY();
                }
            }
        }
        
        return new PairInt(minX, minXsY);
    }

    /**
     * calculate the minimum angle where the vertex is
     * the nearest medial axis point and the ends
     * of the segments are (x,y) and each point in 
     * neighborsX, neighborsY.
     * @param x
     * @param y
     * @param nXY
     * @param neighborsX
     * @param neighborsY
     * @param nn
     * @param medAxis0 closest medial axis point to (x, y)
     * @return index of neighborsX, neighborsY that
     * has the smallest angle subtended by the medial axis
     */
    private int calculateMinAngles(int x, int y,
        int nXY, int[] neighborsX, int[] neighborsY, 
        NearestNeighbor2D nn,
        Set<PairInt> contiguousShapePoints) {

        /*
        to determine the next best point among more than
        one choice, need to determine if the points are
        a corner for example where the diagonal and
        horizontal or vertical are both present,
        OR whether this is 2 borders adjacent to one 
        another due to the region being very narrow
        for example.  This later case requires alot more
        logic than initially present here to avoid making
        crossing perimeter paths.
        
        May need to revisit this:
           For the narrow width regions that are 2 adjacent
        boundaries, the number of remaining neighbors on the
        first pass through the region appears to be > 2,
        and for that case, cannot use the medial axis, but
        can use the number of non-shape neighbors in common.
        */

        if (nXY > 2) {
            // modifies neighborsX and neighborsY
            nXY = filterForLargestNonShapeIntersection(
                x, y, nXY, neighborsX, neighborsY,
                contiguousShapePoints);
            assert(nXY != 0);
            if (nXY == 1) {
                return 0;
            }
        }
        
        // determine whether this is a convex or concave
        // section of the boundary curve.
        
        double xCen = x;
        double yCen = y;
        for (int i = 0; i < nXY; ++i) {
            xCen += neighborsX[i];
            yCen += neighborsY[i];
        }
        xCen /= ((double)nXY + 1);
        yCen /= ((double)nXY + 1);
        
        /*
        convex:
           centroid of the subset of boundary points
           are inside the shape bounds.
           
        concave:
           centroid of the subset of boundary points 
           are outside the shape bounds.
        */
        
        int xc = (int)Math.round(xCen);
        int yc = (int)Math.round(yCen);
          
        Set<PairInt> medAxisCenClosest = nn.findClosest(xc, yc);
        assert(!medAxisCenClosest.isEmpty());
        PairInt medAxisCen = medAxisCenClosest.iterator().next();
        
        PairInt medAxis0 = nn.findClosest(x, y).iterator().next();
        double d0 = distSq(x, y, medAxisCen.getX(), medAxisCen.getY());
        double d0Cen = distSq(xCen, yCen, medAxisCen.getX(), medAxisCen.getY());
          
        boolean isOutside = (d0Cen > d0);
        
        
        System.out.println(
            String.format(
            "**med axis pt=%s (%d,%d) d0=%.4f --> (%.3f,%.3f) d0cen=%.4f",
            medAxisCen.toString(), x, y, (float) d0,
            xCen, yCen, (float) d0Cen));
        
        
        if (!isOutside) {
            for (int i = 0; i < nXY; ++i) {
                int x2 = neighborsX[i];
                int y2 = neighborsY[i];

                double d2 = distSq(x2, y2, 
                    medAxisCen.getX(), medAxisCen.getY());

                double d2Cen = distSq(xCen, yCen, 
                    medAxisCen.getX(), medAxisCen.getY());

                
                System.out.println(
                String.format(
                "  med Axis pt=%s (%d,%d) d2=%.4f --> d2cen=%.4f",
                medAxisCen.toString(), x2, y2, (float) d2,
                (float) d2Cen));
                
                
                if (d2Cen > d2) {
                    isOutside = true;
                }
            }
        }
    
        if (isOutside) {
            return calculateMinAnglesForConcave(x, y, 
                nXY, neighborsX, neighborsY, medAxis0,
                contiguousShapePoints);
        } else {
            return calculateMinAnglesForConvex(x, y, 
                nXY, neighborsX, neighborsY, medAxis0);
        }
    }
    
    /**
     * calculate the minimum angle where the vertex is
     * the nearest medial axis point and the ends
     * of the segments are (x,y) and each point in 
     * neighborsX, neighborsY for a concave section 
     * of the curve.
     * @param x
     * @param y
     * @param nXY
     * @param neighborsX
     * @param neighborsY
     * @param nn
     * @param medAxis0 closest medial axis point to (x, y)
     * @return index of neighborsX, neighborsY that
     * has the smallest angle subtended by the medial axis
     */
    private int calculateMinAnglesForConcave(int x, int y,
        int nXY, int[] neighborsX, int[] neighborsY, 
        PairInt medAxis0, Set<PairInt> contiguousShapePoints) {
        
        double maxAngle = Double.MIN_VALUE;
        int minIdx = -1;
        
        double angle0 = AngleUtil.polarAngleCW(
            x - medAxis0.getX(), y - medAxis0.getY());
        
        for (int i = 0; i < nXY; ++i) {
            int x2 = neighborsX[i];
            int y2 = neighborsY[i];
        
            double angle = AngleUtil.polarAngleCW(
                x2 - medAxis0.getX(), y2 - medAxis0.getY());
            
            
            System.out.println(
                String.format("concave: (%d,%d) a=%.4f --> (%d,%d) a=%.4f",
                    x, y, (float) angle0, 
                    x2,y2, (float) angle));
            
            
            // for convex section:
            if ((angle > angle0) && (angle > maxAngle)) {
                maxAngle = angle;
                minIdx = i;
            }
        }
        
        if (minIdx != -1) {
            return minIdx;
        }
        
        // if region is so thin that medial axis
        // points do not represent the local center
        // of the shape, might end up here when the
        // nearest medial axis point is CCW from 
        // (x,y) and neighborsX,neigbhorsY, etc.
        // - one possible way to better determine
        //   "outward" for the curve for this case
        //   is to use the centroid of the adjacent
        //   non-shape points and convex angle rules
        //   to determine best next idx
        Set<PairInt> spaceNeighbors = findNonShapeNeighbors(
            x, y, contiguousShapePoints);

        MiscellaneousCurveHelper curveHelper = new
            MiscellaneousCurveHelper();

        double[] xyCen = curveHelper.calculateXYCentroids(
            spaceNeighbors);

        PairInt xyCenP = new PairInt(
            (int)Math.round(xyCen[0]),
            (int)Math.round(xyCen[1]));

        // convex angles, w/ CCW rules for 
        // shape being on other side of xyCenP
        // so concave rules again.
                
        double minAngle = Double.MAX_VALUE;
        minIdx = -1;
        angle0 = AngleUtil.polarAngleCCW(
            x - medAxis0.getX(), y - medAxis0.getY());
        for (int i = 0; i < nXY; ++i) {
            int x2 = neighborsX[i];
            int y2 = neighborsY[i];
            double angle = AngleUtil.polarAngleCCW(
                x2 - medAxis0.getX(), y2 - medAxis0.getY());
            
            System.out.println(
            String.format(
            "convex: (%d,%d) a=%.4f --> (%d,%d) a=%.4f",
                    x, y, (float) angle0,
                    x2,y2, (float) angle));
            
            // for convex section:
            if ((angle > angle0) && (angle < minAngle)) {
                minAngle = angle;
                minIdx = i;
            }
        }
        
if (minIdx == -1) {
int z = 1;         
}
        
        return minIdx;
    }
    
    /**
     * calculate the minimum angle where the vertex is
     * the nearest medial axis point and the ends
     * of the segments are (x,y) and each point in 
     * neighborsX, neighborsY for a convex section
     * of the boundary.
     * @param x
     * @param y
     * @param nXY
     * @param neighborsX
     * @param neighborsY
     * @param nn
     * @param medAxis0 closest medial axis point to (x, y)
     * @return index of neighborsX, neighborsY that
     * has the smallest angle subtended by the medial axis
     */
    private int calculateMinAnglesForConvex(int x, int y,
        int nXY, int[] neighborsX, int[] neighborsY, 
        PairInt medAxis0) {
        
        double minAngle = Double.MAX_VALUE;
        int minIdx = -1;
        
        double angle0 = AngleUtil.polarAngleCW(
            x - medAxis0.getX(), y - medAxis0.getY());
        
        for (int i = 0; i < nXY; ++i) {
            int x2 = neighborsX[i];
            int y2 = neighborsY[i];
        
            double angle = AngleUtil.polarAngleCW(
                x2 - medAxis0.getX(), y2 - medAxis0.getY());
     
            /*
            System.out.println(
            String.format(
            "convex: (%d,%d) a=%.4f --> (%d,%d) a=%.4f",
                    x, y, (float) angle0, 
                    x2,y2, (float) angle));
            */
            
            // for convex section:
            if ((angle > angle0) && (angle < minAngle)) {
                minAngle = angle;
                minIdx = i;
            }
        }
        
        assert(minIdx != -1);
        
        return minIdx;
    }

    private double distSq(double x1, double y1, double x2, double y2) {
        double diffX = x1 - x2;
        double diffY = y1 - y2;
        return (diffX * diffX + diffY * diffY);
    }

    private int filterForLargestNonShapeIntersection(
        int x, int y, int nXY, int[] neighborsX, int[] neighborsY, 
        Set<PairInt> contiguousShapePoints) {

        int maxNS = Integer.MIN_VALUE;
        int[] maxNSIdx = new int[8];
        int n = 0;
        for (int i = 0; i < nXY; ++i) {
            int nIterNS = 
                calcNumberOfNonShapeNeigbhorsInCommon
                (x, y, neighborsX[i], neighborsY[i],
                contiguousShapePoints);
            if (nIterNS > maxNS) {
                maxNS = nIterNS;
                n = 0;
                maxNSIdx[n] = i;
                n++;
            } else if (nIterNS == maxNS) {
                maxNSIdx[n] = i;
                n++;
            }
        }

        if (n == nXY) {
            // NOTE: does this happen?
            return nXY;
        }
        
        if (n == 0) {
            return 0;
        }

        // place maxNSIdx at top of neighbors
        assert(nXY < (neighborsX.length - 1));
        
        // shift neighborsX,Y down by 1
        for (int i = (neighborsX.length - 1); i > 0; --i) {
            neighborsX[i] = neighborsX[i - 1];
            neighborsY[i] = neighborsY[i - 1];
        }
        for (int i = 0; i < n; ++i) {
            int idx = maxNSIdx[i] + 1;
            assert(idx > i);
            neighborsX[i] = neighborsX[idx];
            neighborsY[i] = neighborsY[idx];
        }

        return n;
    }

    private int calcNumberOfNonShapeNeigbhorsInCommon(
        int x1, int y1, int x2, int y2, 
        Set<PairInt> contiguousShapePoints) {

        Set<PairInt> set1 = findNonShapeNeighbors(x1, y1, 
            contiguousShapePoints);
        
        Set<PairInt> set2 = findNonShapeNeighbors(x2, y2, 
            contiguousShapePoints);
        
        // set1 = (a + intersection)
        // set2 = (b + intersection)
        // nItersection = set1.size - (set1.removeall(set2)).size
        
        int n1 = set1.size();
        set1.removeAll(set2);
        int nInterNS = n1 - set1.size();
        
        return nInterNS;
    }

    private Set<PairInt> findNonShapeNeighbors(
        int x, int y, Set<PairInt> contiguousShapePoints) {
        
        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;

        Set<PairInt> pointsNS = new HashSet<PairInt>();
        for (int k = 0; k < dxs.length; ++k) {
            int x2 = x + dxs[k];
            int y2 = y + dys[k];
            PairInt p2 = new PairInt(x2, y2);
            if (!contiguousShapePoints.contains(p2)) {
                pointsNS.add(p2);
            }
        }
        
        return pointsNS;
    }

}
