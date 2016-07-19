package algorithms.compGeometry;

import algorithms.imageProcessing.util.AngleUtil;
import algorithms.misc.Misc;
import algorithms.misc.MiscMath;
import algorithms.search.NearestNeighbor2D;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import java.util.HashSet;
import java.util.LinkedList;
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
     * set, contiguousPoints.
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
     * @param contiguousPoints
     * @return 
     */
    public PairIntArray extractOrderedBorder(Set<PairInt> 
        contiguousPoints) {
        
        //O(8*N)
        Set<PairInt> boundary = extractBorder(contiguousPoints);
        
        //dependent upon complexity of shape:
        MedialAxis1 medAxis1 = new MedialAxis1(contiguousPoints, 
            boundary);
        medAxis1.findMedialAxis();
        Set<PairInt> medAxisPts = medAxis1.getMedialAxisPoints();
        
        return extractOrderedBorder(boundary, medAxisPts);
    }
private PairIntArray debug = null;

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
     * @return 
     */
    public PairIntArray extractOrderedBorder(Set<PairInt> 
        borderPoints, Set<PairInt> medialAxisPoints) {
        
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
                            x3, y3, ns2, neighborsX, neighborsY, nn);
                        
                        junctionNodes.add(Integer.valueOf(output.getN()));
                        
                        output.add(neighborsX[minAngleIdx], 
                            neighborsY[minAngleIdx]);                        
                    }
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
 this.debug = output;               
                int minAngleIdx = calculateMinAngles(
                    x, y, ns, neighborsX, neighborsY, nn);

                junctionNodes.add(Integer.valueOf(output.getN()));
                
                output.add(neighborsX[minAngleIdx], neighborsY[minAngleIdx]);                
            }
            
            boolean rmvd = remaining.remove(
                new PairInt(output.getX(output.getN() - 1),
                output.getY(output.getN() - 1)));
            assert(rmvd);
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
        NearestNeighbor2D nn) {
     
        /*
        to determine the next best point among more than
        one choice, need to determine if the points are
        a corner for example where the diagonal and
        horizontal or vertical are both present,
        OR whether this is 2 borders adjacent to one 
        another due to the region being very narrow
        for example.  This later case require alot more
        logic than initially present here.
        
        -- might need to make a longest path algorithm
        over all border points.  the points left out
        because they are part of the concave or
        convex corners can be added back using the
        logic below afterwards.
        Need to refactor all of the dependentt code
        and this method for that change in logic.
      
        */
        
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
             
        PairInt medAxisCen = nn.findClosest(xc, yc).iterator().next();
        
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

                double d2 = distSq(x2, y2, medAxisCen.getX(), medAxisCen.getY());

                double d2Cen = distSq(xCen, yCen, medAxisCen.getX(), medAxisCen.getY());

                
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
                nXY, neighborsX, neighborsY, medAxis0);
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
        PairInt medAxis0) {
        
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
      
        assert(minIdx != -1);
        
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
    
}
