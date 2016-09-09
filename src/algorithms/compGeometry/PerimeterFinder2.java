package algorithms.compGeometry;

import algorithms.imageProcessing.PostLineThinnerCorrections;
import algorithms.imageProcessing.SpurRemover;
import algorithms.imageProcessing.ZhangSuenLineThinner;
import algorithms.misc.Misc;
import algorithms.misc.MiscMath;
import algorithms.search.NearestNeighbor2D;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.PolygonAndPointPlotter;
import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.hash.TObjectIntHashMap;
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
     * note, this intersection method assumes the boundaries are
     * both clockwise ordered points and that ob1 and ob2 are
     * adjacent and do not overlap, and that simple
     * cut and append or inserts are needed.  only the 
     * outer merged boundary is preserved.
     * 
     * @param ob1
     * @param ob2
     * @return 
     */
    public PairIntArray mergeAdjacentOrderedBorders(PairIntArray ob1,
        PairIntArray ob2) {
        
        if (ob1.getN() == 0 && ob2.getN() > 0) {
            return ob2.copy();
        } else if (ob2.getN() == 0 && ob1.getN() > 0) {
            return ob1.copy();
        } else if (ob1.getN() == 0 && ob2.getN() == 0) {
            return new PairIntArray();
        }
        
        PairIntArray a1 = ob1.copy();
        PairIntArray a2 = ob2.copy();
        rotateToMinXYAt0(a1);
        rotateToMinXYAt0(a2);
        
        if ((a2.getX(0) < a1.getX(0)) || (a2.getX(0) == a1.getX(0) && 
            a2.getY(0) < a1.getY(0))) {
            PairIntArray swap = a1;
            a1 = a2;
            a2 = swap;
        }
        
        PairIntArray output = new PairIntArray();
        
        int n1 = a1.getN();
        int n2 = a2.getN();
        
        TObjectIntMap<PairInt> pointIndexMap2 = new TObjectIntHashMap<PairInt>();
        for (int i = 0; i < n2; ++i) {
            PairInt p = new PairInt(a2.getX(i), a2.getY(i));
            pointIndexMap2.put(p, i);
        }
        
        // -- start w/ the leftmost lowermost of the 2 arrays and call that array1:
        //    until find point adjacent to array 2, continue to
        //       add points from array1 to output,
        //    then when find adacent point, mark it as the start of the
        //       intersecting region and break
        //    iterate from end of array1 to smaller indexes looking for the end 
        //       of the intersecting region.
        //    then, add points from array2 in between the 2 marks to output.
        //    then add the points from array1 after the 2 marks to output.
    
        int startIdx1 = -1;
        int startIdx2 = -1;
        
        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;
        for (int i = 0; i < n1; ++i) {
            final int x = a1.getX(i);
            final int y = a1.getY(i);
            // if found a startIdx2, keep the smallest index of adjacent
            for (int k = 0; k < dxs.length; ++k) {
                int x2 = x + dxs[k];
                int y2 = y + dys[k];
                PairInt p2 = new PairInt(x2, y2);
                if (pointIndexMap2.containsKey(p2)) {
                    int idx2 = pointIndexMap2.get(p2);
                    if (startIdx2 == -1) {
                        startIdx2 = idx2;
                    } else {
                        if (idx2 < startIdx2) {
                            startIdx2 = idx2;
                        }
                    }
                }
            }
            output.add(x, y);
            if (startIdx2 > -1) {
                startIdx1 = i;
                break;
            }
        }
        
        if (startIdx1 == -1) {
            return output;
        }
        
        int stopIdx1 = -1;
        int stopIdx2 = -1;
        
        // search from end for stopIdx1
        for (int i = (n1 - 1); i > startIdx1; --i) {
            int x = a1.getX(i);
            int y = a1.getY(i);
            // if stopIdx2 > -1, keep the smallest
            for (int k = 0; k < dxs.length; ++k) {
                int x2 = x + dxs[k];
                int y2 = y + dys[k];
                PairInt p2 = new PairInt(x2, y2);
                if (pointIndexMap2.containsKey(p2)) {
                    int idx2 = pointIndexMap2.get(p2);
                    if (stopIdx2 == -1) {
                        stopIdx2 = idx2;
                    } else {
                        if (idx2 < stopIdx2) {
                            stopIdx2 = idx2;
                        }
                    }
                }
            }
            if (stopIdx2 > -1) {
                stopIdx1 = i;
                break;
            }
        }
        
        if (stopIdx1 == -1) {
            stopIdx1 = startIdx1;
            stopIdx2 = startIdx2;
            // exception would be extreme case of single pixel wide
            //   boundaries... this may change, but for now, expecting
            //   that it would not exist.  exanple: @@@@ ####
            //   throw new IllegalStateException("Error in algorithm.  " +
            //       " did not find last adjacent point");
        }
        
        int end = (startIdx2 > stopIdx2) ? n2 - 1: stopIdx2;
        for (int i = startIdx2; i <= end; ++i) {
            output.add(a2.getX(i), a2.getY(i));
        }
        
        if (startIdx2 > stopIdx2) {
            for (int i = 0; i <= stopIdx2; ++i) {
                output.add(a2.getX(i), a2.getY(i));
            }
        }
        
        for (int i = stopIdx1; i < n1; ++i) {
            output.add(a1.getX(i), a1.getY(i));
        }
       
        /*
              startIdx1=7; stopIdx1=8  n1=11
              startIdx2=3; stopIdx2=1  n2=10
                 5  6
              @  @  @ @ 3
            @       @ # # # #
          @       @9  #     #6
          1 @   @10  0# # # #
            0 @         9 8
        
        
              startIdx1=7; stopIdx1=0  n1=11
              startIdx2=3; stopIdx2=0  n2=8
                      7
              @  @  @ @
            @       @ #3
        2 @       @ # #
            @   @ #  #
              @ # # #
              0   7
        
        
              startIdx1=7 ; stopIdx1=10  n1=10
              startIdx2=4 ; stopIdx2=1   n2=8
                      7      
              @  @  @ @
            @       @ #4 
        2 @       @ # #
            @   @ #2 #
              @   # #7
                  #
        
              startIdx1=4; stopIdx1=7  n1=11
              startIdx2=0; stopIdx2=4  n2=7
                 1  2  
              0  #  # #3
              # 6# 5# #4
              @  @  @ @7
            @    5  @ 
          @       @9  
          1 @   @10  
            0 @
        */
        
        return output;
    }
    
    /**
     * find the minx, min y point and if it's not at index 0, rotate the
     * array to place it there.
     * @param a 
     */
    protected void rotateToMinXYAt0(PairIntArray a) {
        
        if (a.getN() < 2) {
            return;
        }
        
        int minIdx = -1;
        int minX = Integer.MAX_VALUE;
        int minY = Integer.MAX_VALUE;
        for (int i = 0; i < a.getN(); ++i) {
            int x = a.getX(i);
            if (x < minX) {
                minX = x;
                minY = a.getY(i);
                minIdx = i;                
            } else if (x == minX) {
                int y = a.getY(i);
                if (y < minY) {
                    minX = x;
                    minY = y;
                    minIdx = i;
                }
            }
        }
        if (minIdx > 0) {
            a.rotateLeft(minIdx);
            assert(a.getX(0) == minX);
            assert(a.getY(0) == minY);
        }
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
        Set<PairInt> rmPts = new HashSet<PairInt>();
        
        //TODO: methods used for this could be improved
        thinTheBoundary(boundary, rmPts);
        
        Set<PairInt> cPts = new HashSet<PairInt>(contiguousPoints);
        // NOTE: if this is small cavity, removing these points
        // may remove important shape points.
        // a more exact method should test for whether each point is
        // within boundary and if not, remove it.
        cPts.removeAll(rmPts);
        
        // O(N*log_2(N))
        MedialAxis medAxis = new MedialAxis(cPts, boundary);
        medAxis.fastFindMedialAxis();
        
        Set<PairInt> medAxisPts = medAxis.getMedialAxisPoints();
       
        return extractOrderedBorder(boundary, medAxisPts, cPts);
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
     
                        int minAngleIdx = calculateMinAngles( 
                            x3, y3, ns2, neighborsX, neighborsY, nn,
                            contiguousShapePoints);
 
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
          
                int minAngleIdx = calculateMinAngles(
                    x, y, ns, neighborsX, neighborsY, nn,
                    contiguousShapePoints);
                
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
}
