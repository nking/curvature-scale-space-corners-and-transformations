package algorithms.compGeometry;

import algorithms.imageProcessing.ContiguousGapFinder;
import algorithms.imageProcessing.PostLineThinnerCorrections;
import algorithms.imageProcessing.SpurRemover;
import algorithms.imageProcessing.ZhangSuenLineThinner;
import algorithms.misc.Misc;
import algorithms.misc.MiscMath;
import algorithms.search.NearestNeighbor2D;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.PolygonAndPointPlotter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;
import java.util.logging.Logger;

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
     * finds any gaps embedded in the contiguous points.
     * @param contiguousPoints
     * @return 
     */
    public Set<PairInt> findEmbeddedGaps(Set<PairInt> contiguousPoints) {
        
        ContiguousGapFinder finder = new ContiguousGapFinder(
            contiguousPoints);
        finder.setMinimumNumberInCluster(1);
        
        int minX = finder.getMinX();
        int maxX = finder.getMaxX();
        int minY = finder.getMinY();
        int maxY = finder.getMaxY();
                
        finder.findGaps();
        
        int nGroups = finder.getNumberOfGapGroups();
        
        Set<PairInt> embedded = new HashSet<PairInt>();
                
        for (int i = 0; i < nGroups; ++i) {
            Set<PairInt> set = finder.getXY(i);
            boolean foundEdgePoint = false;
            for (PairInt p : set) {
                int x = p.getX();
                int y = p.getY();
                if (x == minX || x == maxX || y == minY || y == maxY) {
                    foundEdgePoint = true;
                    break;
                }
            }
            if (!foundEdgePoint) {
                embedded.addAll(set);
            }
        }
        
        return embedded;
    }
    
    /**
     * finds the outer boundary points of the contiguous points.
     * any embedded holes in the contiguousPoints are not included
     * in the result boundary.
     * @param contiguousPoints
     * @return 
     */
    public Set<PairInt> extractOuterBorder(Set<PairInt> contiguousPoints) {
        
        Set<PairInt> embedded = findEmbeddedGaps(contiguousPoints);
        
        Set<PairInt> set2 = new HashSet<PairInt>(contiguousPoints);
        set2.addAll(embedded);        
        
        return extractBorder(set2);
    }
    
    public void thinTheBoundary(Set<PairInt> boundary,
        Set<PairInt> removedPoints) {
        
        Set<PairInt> b = new HashSet<PairInt>(boundary);
        
        int[] minMaxXY = MiscMath.findMinMaxXY(b); 
        
        ZhangSuenLineThinner lt = new ZhangSuenLineThinner();
        lt.applyLineThinner(b, 
            minMaxXY[0] - 1, minMaxXY[1] + 1, 
            minMaxXY[2] - 1, minMaxXY[3] + 1);
        
        PostLineThinnerCorrections pltc = new PostLineThinnerCorrections();
        pltc._correctForArtifacts(b, minMaxXY[1] + 3, 
            minMaxXY[3] + 3);
        
        SpurRemover spurRm = new SpurRemover();
        spurRm.remove(b, minMaxXY[1] + 3, 
            minMaxXY[3] + 3);
        
        // store the removed points
        removedPoints.addAll(boundary);
        removedPoints.removeAll(b);
        
        //System.out.println("line thinning removed " +
        //    removedPoints.size() + " from the boundary");
        
        boundary.clear();
        boundary.addAll(b);
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
       
        Set<PairInt> outputMedialAxisPoints = null;
        
        return extractOrderedBorder(contiguousPoints, outputMedialAxisPoints);
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
     * @param outputMedialAxisPoints - if not null, the medial axis
     * points are added to this variable
     * @return 
     */
    public PairIntArray extractOrderedBorder(Set<PairInt> 
        contiguousPoints, Set<PairInt> outputMedialAxisPoints) {
        
        if (contiguousPoints.size() < 4) {
            PairIntArray output = new PairIntArray(contiguousPoints.size());
            for (PairInt p : contiguousPoints) {
                output.add(p.getX(), p.getY());
            }
            return output;
        }
        
        //O(8*N)
        //Set<PairInt> boundary = extractOuterBorder(contiguousPoints);
        Set<PairInt> embedded = findEmbeddedGaps(contiguousPoints);
        Set<PairInt> set2 = new HashSet<PairInt>(contiguousPoints);
        set2.addAll(embedded);        
        Set<PairInt> boundary = extractBorder(set2);
        
        Set<PairInt> rmPts = new HashSet<PairInt>();
        
        //TODO: methods used for this could be improved
        thinTheBoundary(boundary, rmPts);
        
        Set<PairInt> cPts = new HashSet<PairInt>(set2);
        // NOTE: if this is small cavity, removing these points
        // may remove important shape points.
        // a more exact method should test for whether each point is
        // within boundary and if not, remove it.
        cPts.removeAll(rmPts);
        
        // O(N*log_2(N))
        MedialAxis medAxis = new MedialAxis(cPts, boundary);
        medAxis.fastFindMedialAxis();
        
        Set<PairInt> medAxisPts = medAxis.getMedialAxisPoints();
       
        if (outputMedialAxisPoints != null) {
            outputMedialAxisPoints.addAll(medAxisPts);
        }
        
        return extractOrderedBorder(boundary, medAxisPts, cPts);
        //return extractOrderedBorder00(boundary, medAxisPts);
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
        LinkedHashSet<Integer> jassocInserted = new LinkedHashSet<Integer>();
        
        boolean err = false;
        while (!remaining.isEmpty()) {
            err = false;
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
                        debug(contiguousShapePoints, output);      
                        //throw new IllegalStateException(
                        Logger.getLogger(this.getClass().getName()).warning(
                            "Error in algorithm:"
                            + " no adjacent points for (" + x + " " + y
                            + ") but remaining is not empty and "
                            + " junction list is empty");
                        err = true;
                        break;
                    }
                    int x3 = output.getX(index.intValue());
                    int y3 = output.getY(index.intValue());

                    ns2 = findNeighbors(x3, y3, remaining,
                        neighborsX, neighborsY);
                    
                    while (ns2 == 0 && !junctionNodes.isEmpty()) {
                        index = junctionNodes.pollLast();
                        x3 = output.getX(index.intValue());
                        y3 = output.getY(index.intValue());
                        ns2 = findNeighbors(x3, y3, remaining,
                             neighborsX, neighborsY);
                    }
                
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
     
                        int prevX = -1;
                        int prevY = -1;
                        if (output.getN() > 1) {
                            prevX = output.getX(output.getN() - 2);
                            prevY = output.getY(output.getN() - 2);
                        }
                        int minAngleIdx = calculateMinAngles( 
                            x3, y3, ns2, neighborsX, neighborsY, nn,
                            contiguousShapePoints, prevX, prevY);
 
                        {//DEBUG
                            if (minAngleIdx == -1) {
                                try {

                                    int[] xPolygon = null;
                                    int[] yPolygon = null;
                                    PolygonAndPointPlotter plotter = new PolygonAndPointPlotter();
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
                                        xp, yp, xPolygon, yPolygon, "shape 1");
                                    
                                    n = borderPoints.size();
                                    xp = new int[n];
                                    yp = new int[n];
                                    count = 0;
                                    for (PairInt p : borderPoints) {
                                        xp[count] = p.getX();
                                        yp[count] = p.getY();
                                        count++;
                                    }
                                    plotter.addPlot(xminmxyminmiac[0], xminmxyminmiac[1],
                                        xminmxyminmiac[2], xminmxyminmiac[3],
                                        xp, yp, xPolygon, yPolygon, "border 1");
                                    
                                    n = medialAxisPoints.size();
                                    xp = new int[n];
                                    yp = new int[n];
                                    count = 0;
                                    for (PairInt p : medialAxisPoints) {
                                        xp[count] = p.getX();
                                        yp[count] = p.getY();
                                        count++;
                                    }
                                    plotter.addPlot((float) xminmxyminmiac[0],
                                        (float) xminmxyminmiac[1],
                                        (float) xminmxyminmiac[2], (float) xminmxyminmiac[3],
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

                                    System.out.println("boundary.size=" + borderPoints.size()
                                        + " output.size=" + output.getN());

                                    System.out.println("output=" + output.toString());
                                } catch (Throwable t) {

                                }
                            }
                        }
                        
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
                
                int prevX = -1;
                int prevY = -1;
                if (output.getN() > 1) {
                    prevX = output.getX(output.getN() - 2);
                    prevY = output.getY(output.getN() - 2);
                }
          
                int minAngleIdx = calculateMinAngles(
                    x, y, ns, neighborsX, neighborsY, nn,
                    contiguousShapePoints, prevX, prevY);
                
                {//DEBUG
                    if (minAngleIdx == -1) {
                        try {

                            int[] xPolygon = null;
                            int[] yPolygon = null;
                            PolygonAndPointPlotter plotter = new PolygonAndPointPlotter();
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
                                xp, yp, xPolygon, yPolygon, "shape 2");

                            n = borderPoints.size();
                            xp = new int[n];
                            yp = new int[n];
                            count = 0;
                            for (PairInt p : borderPoints) {
                                xp[count] = p.getX();
                                yp[count] = p.getY();
                                count++;
                            }
                            plotter.addPlot(xminmxyminmiac[0], xminmxyminmiac[1],
                                xminmxyminmiac[2], xminmxyminmiac[3],
                                xp, yp, xPolygon, yPolygon, "border 2");
                                    
                            n = medialAxisPoints.size();
                            xp = new int[n];
                            yp = new int[n];
                            count = 0;
                            for (PairInt p : medialAxisPoints) {
                                xp[count] = p.getX();
                                yp[count] = p.getY();
                                count++;
                            }
                            plotter.addPlot((float) xminmxyminmiac[0],
                                (float) xminmxyminmiac[1],
                                (float) xminmxyminmiac[2], (float) xminmxyminmiac[3],
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

                            System.out.println("boundary.size=" + borderPoints.size()
                                + " output.size=" + output.getN());

                            System.out.println("output=" + output.toString());
                        } catch (Throwable t) {

                        }
                    }
                }

                if (minAngleIdx == -1) {
                    return output;
                }
                
                junctionNodes.add(Integer.valueOf(output.getN()));
     
                output.add(neighborsX[minAngleIdx], neighborsY[minAngleIdx]);                
            }
            
            if (err) {
                System.err.println("not adding back: " +
                    remaining.toString());
                remaining.clear();
            } else {           
                boolean rmvd = remaining.remove(
                    new PairInt(output.getX(output.getN() - 1),
                    output.getY(output.getN() - 1)));
                assert(rmvd);
            }
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
            /*
            // if any adjacent have dist > 2, see if swapping would
            // reduce that
            for (int i = 0; i < (output.getN() - 1); ++i) {
                int x = output.getX(i);
                int y = output.getY(i);
                int nextX = output.getX(i + 1);
                int nextY = output.getY(i + 1);
                double d1 = Math.sqrt(distSq(nextX, nextY, x, y));
                if (d1 >= 2 && (i > 0)) {
                    double dpn = Math.sqrt(distSq(nextX, nextY, 
                        output.getX(i - 1), output.getY(i - 1)));
                    if (dpn <= 2) {
                        
                    }
                }
            }
            */
        }
    
        
        {//DEBUG
            if (false)
            try {

                int[] xPolygon = null;
                int[] yPolygon = null;
                PolygonAndPointPlotter plotter = new PolygonAndPointPlotter();
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
                    xp, yp, xPolygon, yPolygon, "shape 3");

                n = medialAxisPoints.size();
                xp = new int[n];
                yp = new int[n];
                count = 0;
                for (PairInt p : medialAxisPoints) {
                    xp[count] = p.getX();
                    yp[count] = p.getY();
                    count++;
                }
                plotter.addPlot((float) xminmxyminmiac[0],
                    (float) xminmxyminmiac[1],
                    (float) xminmxyminmiac[2], (float) xminmxyminmiac[3],
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

                System.out.println("boundary.size=" + borderPoints.size()
                    + " output.size=" + output.getN());

                System.out.println("output=" + output.toString());
            } catch (Throwable t) {

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
        NearestNeighbor2D nn, Set<PairInt> contiguousShapePoints,
        int prevX, int prevY) {

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
        }
        if (nXY == 1) {
            return 0;
        }
       
        /*
             P1      PmedAxis

                P2
        */
        
        Set<PairInt> medAxisClosest = nn.findClosest(x, y);
        assert(!medAxisClosest.isEmpty());
        PairInt medAxisP = medAxisClosest.iterator().next();
       
        /*
        System.out.println(
            String.format("x,y=(%d,%d) medAxis=(%d,%d)", x, y, 
                medAxisP.getX(), medAxisP.getY()));
        */
        
        double minAngle = Double.MAX_VALUE;
        int minIdx = -1;
        
        /*
        System.out.println(
            String.format("    (%d,%d) (%d,%d) (%d,%d) nXY=%d", 
                prevX, prevY, x, y, medAxisP.getX(), medAxisP.getY(), nXY));
        */
        
        // if (x,y) == medAxisP or (x2,y2) == medAxisP
        // angle is NAN

        if (!(medAxisP.getX() == x && medAxisP.getY() == y)) {
            for (int i = 0; i < nXY; ++i) {
                int x2 = neighborsX[i];
                int y2 = neighborsY[i];

                double angle = LinesAndAngles.calcClockwiseAngle(
                    x, y, x2, y2, medAxisP.getX(), medAxisP.getY());

                /*
                System.out.println(
                    String.format("    (%d,%d) a=%.4f", x2, y2, (float) angle));
                */
                
                if (angle < minAngle) {
                    minAngle = angle;
                    minIdx = i;
                }
            }
        } else {
            // since the medial axis is same as (x,y), try using
            // the previous point in output as (x,y) and
            // current (x,y) for the medAxisP
            if (prevX == -1) {
                return -1;
            }
            for (int i = 0; i < nXY; ++i) {
                int x2 = neighborsX[i];
                int y2 = neighborsY[i];

                double angle = LinesAndAngles.calcClockwiseAngle(
                    prevX, prevY, x2, y2, x, y);

                /*
                System.out.println(
                    String.format("    *(%d,%d) a=%.4f", 
                        x2, y2, (float) angle));
                */
                
                if (angle < minAngle) {
                    minAngle = angle;
                    minIdx = i;
                }
            }
        }
        
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
    
    private void debug(Set<PairInt> contiguousShapePoints,
        PairIntArray output) {
        try {

            int[] xPolygon = null;
            int[] yPolygon = null;
            PolygonAndPointPlotter plotter = new PolygonAndPointPlotter();
            int[] xminmxyminmiac = MiscMath.findMinMaxXY(
                contiguousShapePoints);
            int[] xp, yp;
            int n, count;

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
    }

    private void debugPrint(PairIntArray ob1, PairIntArray ob2, 
        PairIntArray output) {
        
        int[] minMaxXY1 = MiscMath.findMinMaxXY(ob1);
        int[] minMaxXY2 = MiscMath.findMinMaxXY(ob2);
        int[] minMaxXY = MiscMath.findMinMaxXY(output);
        
        int maxX = Math.max(minMaxXY1[1], minMaxXY2[1]);
        maxX = Math.max(maxX, minMaxXY[1]);
        
        int maxY = Math.max(minMaxXY1[3], minMaxXY2[3]);
        maxY = Math.max(maxY, minMaxXY[3]);
        
        algorithms.imageProcessing.Image img 
            = new algorithms.imageProcessing.Image(maxX + 10, maxY + 10);
        
        long ts = algorithms.misc.MiscDebug.getCurrentTimeFormatted();
        
        algorithms.imageProcessing.ImageIOHelper.
             addCurveToImage(ob1, img, 
             2, 0, 255, 0);
        algorithms.imageProcessing.ImageIOHelper.
             addCurveToImage(ob2, img, 
             1, 0, 0, 255);
        algorithms.imageProcessing.ImageIOHelper.
             addCurveToImage(output, img, 
             0, 255, 0, 0);
        algorithms.misc.MiscDebug.writeImage(img, "_m_" + ts);
        System.out.println("wrote " + ts);
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
     * @return 
     */
    public PairIntArray extractOrderedBorder00(Set<PairInt> 
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
        
        // because the juntionNodes may be inserted correctly, 
        // especially at the end of a complex perimeter
        // with remaining points not yet added,
        // keeping a separate list of junctionNodes
        // that are added to output while backtracking in order
        // to revise them later if needed.
        LinkedHashSet<Integer> jassocInserted = new LinkedHashSet<Integer>();
        
        boolean err = false;
        while (!remaining.isEmpty()) {
            err = false;
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
                        //throw new IllegalStateException(
                        Logger.getLogger(this.getClass().getName()).warning(
                            "Error in algorithm:"
                            + " no adjacent points for (" + x + " " + y
                            + ") but remaining is not empty and "
                            + " junction list is empty");
                        err = true;
                        break;
                    }
                    int x3 = output.getX(index.intValue());
                    int y3 = output.getY(index.intValue());

                    ns2 = findNeighbors(x3, y3, remaining,
                        neighborsX, neighborsY);
                    
                    while (ns2 == 0 && !junctionNodes.isEmpty()) {
                        index = junctionNodes.pollLast();
                        x3 = output.getX(index.intValue());
                        y3 = output.getY(index.intValue());
                        ns2 = findNeighbors(x3, y3, remaining,
                             neighborsX, neighborsY);
                    }
                
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
     
                        int prevX = -1;
                        int prevY = -1;
                        if (output.getN() > 1) {
                            prevX = output.getX(output.getN() - 2);
                            prevY = output.getY(output.getN() - 2);
                        }
                        int minAngleIdx = calculateMinAngles00(
                            x3, y3, ns2, neighborsX, neighborsY, nn,
                            prevX, prevY);
                
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
          
                int prevX = -1;
                int prevY = -1;
                if (output.getN() > 1) {
                    prevX = output.getX(output.getN() - 2);
                    prevY = output.getY(output.getN() - 2);
                }
                int minAngleIdx = calculateMinAngles00(
                    x, y, ns, neighborsX, neighborsY, nn,
                    prevX, prevY);
                
                junctionNodes.add(Integer.valueOf(output.getN()));
     
                output.add(neighborsX[minAngleIdx], neighborsY[minAngleIdx]);                
            }
            
            if (err) {
                System.err.println("not adding back: " +
                    remaining.toString());
                remaining.clear();
            } else {           
                boolean rmvd = remaining.remove(
                    new PairInt(output.getX(output.getN() - 1),
                    output.getY(output.getN() - 1)));
                assert(rmvd);
            }
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
            /*
            // if any adjacent have dist > 2, see if swapping would
            // reduce that
            for (int i = 0; i < (output.getN() - 1); ++i) {
                int x = output.getX(i);
                int y = output.getY(i);
                int nextX = output.getX(i + 1);
                int nextY = output.getY(i + 1);
                double d1 = Math.sqrt(distSq(nextX, nextY, x, y));
                if (d1 >= 2 && (i > 0)) {
                    double dpn = Math.sqrt(distSq(nextX, nextY, 
                        output.getX(i - 1), output.getY(i - 1)));
                    if (dpn <= 2) {
                        
                    }
                }
            }
            */
        }
    
        return output;
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
    private int calculateMinAngles00(int x, int y,
        int nXY, int[] neighborsX, int[] neighborsY, 
        NearestNeighbor2D nn, int prevX, int prevY) {

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

        if (nXY == 1) {
            return 0;
        }
       
        /*
             P1      PmedAxis

                P2
        */
        
        Set<PairInt> medAxisClosest = nn.findClosest(x, y);
        assert(!medAxisClosest.isEmpty());
        PairInt medAxisP = medAxisClosest.iterator().next();
       
        /*
        System.out.println(
            String.format("x,y=(%d,%d) medAxis=(%d,%d)", x, y, 
                medAxisP.getX(), medAxisP.getY()));
        */
        
        double minAngle = Double.MAX_VALUE;
        int minIdx = -1;

        /*
        System.out.println(
            String.format("    (%d,%d) (%d,%d) (%d,%d) nXY=%d", 
                prevX, prevY, x, y, medAxisP.getX(), medAxisP.getY(), nXY));
        */
        
        // if (x,y) == medAxisP or (x2,y2) == medAxisP
        // angle is NAN

        if (!(medAxisP.getX() == x && medAxisP.getY() == y)) {
            for (int i = 0; i < nXY; ++i) {
                int x2 = neighborsX[i];
                int y2 = neighborsY[i];

                double angle = LinesAndAngles.calcClockwiseAngle(
                    x, y, x2, y2, medAxisP.getX(), medAxisP.getY());

                /*
                System.out.println(
                    String.format("    (%d,%d) a=%.4f", x2, y2, (float) angle));
                */
                
                if (angle < minAngle) {
                    minAngle = angle;
                    minIdx = i;
                }
            }
        } else {
            // since the medial axis is same as (x,y), try using
            // the previous point in output as (x,y) and
            // current (x,y) for the medAxisP
            if (prevX == -1) {
                return -1;
            }
            for (int i = 0; i < nXY; ++i) {
                int x2 = neighborsX[i];
                int y2 = neighborsY[i];

                double angle = LinesAndAngles.calcClockwiseAngle(
                    prevX, prevY, x2, y2, x, y);

                /*
                System.out.println(
                    String.format("    *(%d,%d) a=%.4f", 
                        x2, y2, (float) angle));
                */
                
                if (angle < minAngle) {
                    minAngle = angle;
                    minIdx = i;
                }
            }
        }
        
        return minIdx;
    }
}
