package algorithms.compGeometry;

import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.SpurRemover;
import algorithms.misc.Misc;
import algorithms.misc.MiscMath;
import algorithms.search.NearestNeighbor2D;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.PixelHelper;
import algorithms.util.PolygonAndPointPlotter;
import algorithms.VeryLongBitString;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.list.TIntList;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.Stack;

/**
 * class to hold some of the newer methods
 * w.r.t. boundary point extractions.
 *
 * This class is deprecated.  Use methods in PerimeterFinder3 instead.
 * 
 * @author nichole
 */
public class PerimeterFinder2 {
   
    /**
     * This class is deprecated.  Use methods in PerimeterFinder3 instead.
     *
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
        
        int[] minmaxXY = MiscMath.findMinMaxXY(contiguousPoints);
        
        // visit the 1 pixel region surrounding the shape and
        // place the pixels in a stack.
        // then visit their neighbors that are not in contig points
        // until have reached them all
        
        int startX = minmaxXY[0] - 1;
        int startY = minmaxXY[2] - 1;
        int stopX = minmaxXY[1] + 1;
        int stopY = minmaxXY[3] + 1;
        
        Stack<PairInt> stack = new Stack<PairInt>();
        for (int i = startX; i <= stopX; ++i) {
            stack.add(new PairInt(i, startY));
            stack.add(new PairInt(i, stopY));
        }
        for (int j = startY+1; j <= stopY-1; ++j) {
            stack.add(new PairInt(startX, j));
            stack.add(new PairInt(stopX, j));
        }
        
        Set<PairInt> visited = new HashSet<PairInt>();
        
        Set<PairInt> surrounding = new HashSet<PairInt>();
        int[] dxs = Misc.dx4;
        int[] dys = Misc.dy4;
        while (!stack.isEmpty()) {
            PairInt s = stack.pop();
            if (visited.contains(s)) {
                continue;
            }
            surrounding.add(s);
            int x = s.getX();
            int y = s.getY();
            for (int k = 0; k < dxs.length; ++k) {
                int x2 = x + dxs[k];
                int y2 = y + dys[k];
                if (x2 < startX || y2 < startY || x2 > stopX ||
                    y2 > stopY) {
                    continue;
                }
                PairInt p2 = new PairInt(x2, y2);
                if (!contiguousPoints.contains(p2)) {
                    stack.add(p2);
                }
            }
            visited.add(s);
        }
       
        // visit entire region within min and max, and place
        // any point not in surrounding nor in contig into
        // embedded
        startX++;
        stopX--;
        startY++;
        stopY--;
        Set<PairInt> embedded = new HashSet<PairInt>();
        for (int i = startX; i <= stopX; ++i) {
            for (int j = startY; j <= stopY; ++j) {
                PairInt p = new PairInt(i, j);
                if (!contiguousPoints.contains(p) && !surrounding.contains(p)) {
                    embedded.add(p);
                }
            }
        }
        
        return embedded;
    }
    
    /**
     * finds any gaps embedded in the contiguous points.
     * @return 
     */
    public Set<PairInt> findEmbeddedGaps(TIntList contiguousXPoints,
        TIntList contiguousYPoints) {
        
        Set<PairInt> contiguousPoints = new HashSet<PairInt>();
        for (int i = 0; i < contiguousXPoints.size(); ++i) {
            contiguousPoints.add(new PairInt(
                contiguousXPoints.get(i), contiguousYPoints.get(i)));
        }
        
        return findEmbeddedGaps(contiguousPoints);
    }
    
    /**
     * This class is deprecated.  Use methods in PerimeterFinder3 instead.
     * finds any gaps embedded in the contiguous points.
     * @param contiguousPoints
     * @param outputEmbedded output variable
     * @param outputBoundary output variable
     */
    public void extractBorder2(Set<PairInt> contiguousPoints,
        Set<PairInt> outputEmbedded, Set<PairInt> outputBoundary) {
        
        int[] minmaxXY = MiscMath.findMinMaxXY(contiguousPoints);
        
        // visit the 1 pixel region surrounding the shape and
        // place the pixels in a stack.
        // then visit their neighbors that are not in contig points
        // until have reached them all
        
        int startX = minmaxXY[0] - 1;
        int startY = minmaxXY[2] - 1;
        int stopX = minmaxXY[1] + 1;
        int stopY = minmaxXY[3] + 1;
        
        Stack<PairInt> stack = new Stack<PairInt>();
        for (int i = startX; i <= stopX; ++i) {
            stack.add(new PairInt(i, startY));
            stack.add(new PairInt(i, stopY));
        }
        for (int j = startY+1; j <= stopY-1; ++j) {
            stack.add(new PairInt(startX, j));
            stack.add(new PairInt(stopX, j));
        }
        
        Set<PairInt> visited = new HashSet<PairInt>();
        Set<PairInt> surrounding = new HashSet<PairInt>();
        int[] dxs = Misc.dx4;
        int[] dys = Misc.dy4;
        while (!stack.isEmpty()) {
            PairInt s = stack.pop();
            if (visited.contains(s)) {
                continue;
            }
            surrounding.add(s);
            int x = s.getX();
            int y = s.getY();
            for (int k = 0; k < dxs.length; ++k) {
                int x2 = x + dxs[k];
                int y2 = y + dys[k];
                if (x2 < startX || y2 < startY || x2 > stopX || y2 > stopY) {
                    continue;
                }
                PairInt p2 = new PairInt(x2, y2);
                if (!contiguousPoints.contains(p2)) {
                    // add the spaces to the stack
                    stack.add(p2);
                } else {
                    // if not a space, it's an outer boundary point
                    outputBoundary.add(p2);
                }
            }
            visited.add(s);
        }
        
        // visit entire region within min and max, and place
        // any point not in surrounding nor in contig into
        // embedded
        startX++;
        stopX--;
        startY++;
        stopY--;
        for (int i = startX; i <= stopX; ++i) {
            for (int j = startY; j <= stopY; ++j) {
                PairInt p = new PairInt(i, j);
                if (!contiguousPoints.contains(p) && !surrounding.contains(p)) {
                    outputEmbedded.add(p);
                }
            }
        }
        
        /*try {

            int[] xPolygon = null;
            int[] yPolygon = null;
            PolygonAndPointPlotter plotter = new PolygonAndPointPlotter();
            int[] xminmxyminmiac = MiscMath.findMinMaxXY(contiguousPoints);
            int[] xp, yp;
            int n, count;

            n = contiguousPoints.size();
            xp = new int[n];
            yp = new int[n];
            count = 0;
            for (PairInt p : contiguousPoints) {
                xp[count] = p.getX();
                yp[count] = p.getY();
                count++;
            }
            plotter.addPlot(xminmxyminmiac[0], xminmxyminmiac[1],
                xminmxyminmiac[2], xminmxyminmiac[3],
                xp, yp, xPolygon, yPolygon, "shape");
            
            n = outputBoundary.size();
            xp = new int[n];
            yp = new int[n];
            count = 0;
            for (PairInt p : outputBoundary) {
                xp[count] = p.getX();
                yp[count] = p.getY();
                count++;
            }
            plotter.addPlot(xminmxyminmiac[0], xminmxyminmiac[1],
                xminmxyminmiac[2], xminmxyminmiac[3],
                xp, yp, xPolygon, yPolygon, "boundary");
          
            n = outputEmbedded.size();
            xp = new int[n];
            yp = new int[n];
            count = 0;
            for (PairInt p : outputEmbedded) {
                xp[count] = p.getX();
                yp[count] = p.getY();
                count++;
            }
            plotter.addPlot(xminmxyminmiac[0], xminmxyminmiac[1],
                xminmxyminmiac[2], xminmxyminmiac[3],
                xp, yp, xPolygon, yPolygon, "embedded");
            
            String fl = plotter.writeFile3();

            System.out.println("output=" + fl);
        } catch (Throwable t) {

        }*/
    }
    
    /**
     * This class is deprecated.  Use methods in PerimeterFinder3 instead.
     * finds any gaps embedded in the contiguous points.
     * The runtime complexity is roughly O(N_contiguousPoints).
     * 
     * @param contiguousPoints
     * @param outputEmbedded output variable
     * @param outputBoundary output variable
     */
    public void extractBorder2(TIntSet contiguousPoints,
        TIntSet outputEmbedded, TIntSet outputBoundary,
        int imgWidth) {
        
        int[] minmaxXY = MiscMath.findMinMaxXY(contiguousPoints, 
            imgWidth);
        
        // visit the 1 pixel region surrounding the shape and
        // place the pixels in a stack.
        // then visit their neighbors that are not in contig points
        // until have reached them all
        
        int startX = minmaxXY[0] - 1;
        int startY = minmaxXY[2] - 1;
        int stopX = minmaxXY[1] + 1;
        int stopY = minmaxXY[3] + 1;
        
        Stack<PairInt> stack = new Stack<PairInt>();
        for (int i = startX; i <= stopX; ++i) {
            stack.add(new PairInt(i, startY));
            stack.add(new PairInt(i, stopY));
        }
        for (int j = startY+1; j <= stopY-1; ++j) {
            stack.add(new PairInt(startX, j));
            stack.add(new PairInt(stopX, j));
        }
        
        PixelHelper ph = new PixelHelper();
        
        TIntSet visited = new TIntHashSet();
        TIntSet surrounding = new TIntHashSet();
        int[] dxs = Misc.dx4;
        int[] dys = Misc.dy4;
        while (!stack.isEmpty()) {
            PairInt s = stack.pop();
            int x = s.getX();
            int y = s.getY();
            int pixIdx = (int)ph.toPixelIndex(x, y, imgWidth);
            if (visited.contains(pixIdx)) {
                continue;
            }
            surrounding.add(pixIdx);
            
            for (int k = 0; k < dxs.length; ++k) {
                int x2 = x + dxs[k];
                int y2 = y + dys[k];
                if (x2 < startX || y2 < startY || x2 > stopX || y2 > stopY) {
                    continue;
                }
                int pixIdx2 = (y2 * imgWidth) + x2;
                if (!contiguousPoints.contains(pixIdx2)) {
                    // add the spaces to the stack
                    stack.add(new PairInt(x2, y2));
                } else {
                    // if not a space, it's an outer boundary point
                    outputBoundary.add(pixIdx2);
                }
            }
            visited.add(pixIdx);
        }
        
        // visit entire region within min and max, and place
        // any point not in surrounding nor in contig into
        // embedded
        startX++;
        stopX--;
        startY++;
        stopY--;
        for (int i = startX; i <= stopX; ++i) {
            for (int j = startY; j <= stopY; ++j) {
                int pixIdx = (int)ph.toPixelIndex(i, j, imgWidth);
                if (!contiguousPoints.contains(pixIdx) && !surrounding.contains(pixIdx)) {
                    outputEmbedded.add(pixIdx);
                }
            }
        }
        
    }

    /**
     * This class is deprecated.  Use methods in PerimeterFinder3 instead.
     *
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
        Set<PairInt> removedPoints,
        Set<PairInt> contiguousPoints) {
        
        Set<PairInt> b = new HashSet<PairInt>(boundary);
        
        int[] minMaxXY = MiscMath.findMinMaxXY(b); 
        
        ImageProcessor imp = new ImageProcessor();
        imp.applyThinning(boundary, minMaxXY[1] + 1, 
            minMaxXY[3] + 1);
        
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
     * This class is deprecated.  Use methods in PerimeterFinder3 instead.
     *
     * NOT READY FOR USE...needs alot more testing.
     * 
     * given a set of contiguous points, fills in embedded points
     * and then extracts the outer boundary points, and then
     * orders the outer boundary points.
     * For best results, the contiguous points might need to be pre-processed
     * to smooth the curve and in some cases, one might want to separate the
     * region into more than one where the region is thin enough to prevent
     * two pixel paths through it as needed by a spatially sequential closed curve.
     * 
     * NOTE: it is up to the invoker to 
     * give the method a point set which is contiguous.
     * 
     * @param contiguousPoints
     * @return 
     */
    public PairIntArray extractOrderedBorder(Set<PairInt> contiguousPoints) {
       
        if (contiguousPoints.size() < 4) {
            PairIntArray output = new PairIntArray(contiguousPoints.size());
            for (PairInt p : contiguousPoints) {
                output.add(p.getX(), p.getY());
            }
            return output;
        }
        
        //O(8*N)
        Set<PairInt> embedded = new HashSet<PairInt>();
        Set<PairInt> boundary = new HashSet<PairInt>();
        extractBorder2(contiguousPoints, embedded, boundary);
        Set<PairInt> set2 = new HashSet<PairInt>(contiguousPoints);
        set2.addAll(embedded);
        
        if (boundary == null || boundary.size() == 0) {
            return null;
        }
        
        /*
        the algorithm finds the leftmost and smallest xy point in the boundary, 
        then traverses its unadded boundary neighbors to find the neighbor 
        with the smallest clockwise angle from it and previous point.
        For the case that an immediate unadded boundary neighbor is not present
        as is the case at the end of a single pixel wide spike, for example,
        the algorithm walks back up the output list looking for a point which
        has an unadded boundary neighbor and continues from there.
        
        note also, the junctions are found ahead of time and removed
        from the junctionSet as they are added to the output array.
        corners are a subset of junctions and they are excluded by the small
        clockwise angle ordering, but if they are not added on a closed
        curve return to that region, they are added later.
        
        NOTE: for best results, the user might want to have pre-processed the
        points to remove such features.
        */
        
        //using a list so can use bit vectors for quick set intersections
        List<PairInt> boundaryList = new ArrayList<PairInt>(boundary);
        
        TObjectIntMap<PairInt> pointIndexes = new TObjectIntHashMap<PairInt>();
        for (int i = 0; i < boundaryList.size(); ++i) {
            pointIndexes.put(boundaryList.get(i), i);
        }
        
        // index = index of boundaryList, item = bitstring with neighbors set
        VeryLongBitString[] pointNeighbors = createPointNeighborArray(
            boundaryList, pointIndexes);
        
        // corner junctions
        TIntSet junctions = findJunctions(pointNeighbors);
       
        PairIntArray orderedOutput = new PairIntArray(boundary.size());
        
        Set<PairInt> remaining = new HashSet<PairInt>(boundary);
                
        PairInt curr = findSmallestXY(remaining);
        remaining.remove(curr);
        orderedOutput.add(curr.getX(), curr.getY());
        junctions.remove(pointIndexes.get(curr));
        
        // fake point for angle calc refs.  since curr is smallest x, should be
        // safe to make a point to left of it.
        PairInt prev = new PairInt(curr.getX() - 1, curr.getY());
        
        // idx w.r.t. orderedOutput.  used for walking back up array
        int outIdx = -1;
        
        while (!remaining.isEmpty()) {
       
            // among the neighbors of curr that are in remaining, chose
            // the one which has the smallest clockwise angle w/ prev and curr
            
            double minAngle = Double.MAX_VALUE;
            PairInt minP = null;
            
            int x = curr.getX();
            int y = curr.getY();
            int xPrev = prev.getX();
            int yPrev = prev.getY();
            
            int[] currNbrs = pointNeighbors[pointIndexes.get(curr)].getSetBits();
            
            for (int nbrIdx : currNbrs) {
                PairInt p2 = boundaryList.get(nbrIdx);
                if (!remaining.contains(p2)) {
                    continue;
                }
                assert(!prev.equals(p2));
                double angle = LinesAndAngles.calcClockwiseAngle(
                    xPrev, yPrev, p2.getX(), p2.getY(), x, y);
                if (Double.isNaN(angle)) {
                    // can occur is prev==p or p2==p
                    throw new IllegalStateException("error: unexpected NaN, "
                        + " due to equal points");
                } 
                    
                if (angle < minAngle) {
                    minAngle = angle;
                    minP = p2;
                }
            }
            
            if (minP == null) {
                
                // when remaining.size is >= (bounds.size - junctions.size)
                //   start checking for whether have reached the starting
                //   point, and look for whether the remaining points are
                //   junctions that need to be inserted in their 90 degree 
                //   positions.
                if (remaining.size() <= (boundaryList.size() - junctions.size())) {
                    if (isAdjacent(curr, orderedOutput.getX(0), 
                        orderedOutput.getY(0))) {
                        
                        TIntSet remIdxs = new TIntHashSet();
                        for (PairInt p : remaining) {
                            remIdxs.add(pointIndexes.get(p));
                        }
                        remIdxs.removeAll(junctions);
                        if (remIdxs.isEmpty()) {
                            break;
                        }
                    }
                }
                
                // walk back up the output to try previous points, looking for
                //   one w/ an unassigned neighbor
                outIdx--;
                
                if (outIdx < 1) {
                    
                    int nj = junctions.size();
                    int nb = boundaryList.size();
                    
                    debug(contiguousPoints, boundary, orderedOutput);
                    
                    System.out.println("output.n=" 
                        + orderedOutput.getN() 
                        + " rem.n=" + remaining.size() + " " + 
                        boundary.size());
                    
                    throw new IllegalStateException("Error in closed curve shape.");
                }
                
                curr = new PairInt(orderedOutput.getX(outIdx), 
                    orderedOutput.getY(outIdx));
                
                prev = new PairInt(orderedOutput.getX(outIdx - 1), 
                    orderedOutput.getY(outIdx - 1));
                                
                continue;
            }
            
            assert(minP != null);
            
            prev = curr;
            curr = minP;
        
            outIdx = orderedOutput.getN();
            
            remaining.remove(curr);
            orderedOutput.add(curr.getX(), curr.getY());
            junctions.remove(pointIndexes.get(curr));
        }
        
        if (!junctions.isEmpty()) {
            //insert between junction neigbhors
            //TODO: this could be improved, but for now, will use
            //  a pattern that has worse case runtime near
            //  O(junctions.size * output.size)
            TIntIterator iter = junctions.iterator();
            while (iter.hasNext()) {
                int jIdx = iter.next();
                PairInt p = boundaryList.get(jIdx);
                int[] nbrs = pointNeighbors[jIdx].getSetBits();
                Set<PairInt> nbrsSet = new HashSet<PairInt>();
                for (int nbrIdx : nbrs) {
                    nbrsSet.add(boundaryList.get(nbrIdx));
                }
                
                // add junction after first of it's neighbors it finds
                boolean found = false;
                for (int ii = 0; ii < orderedOutput.getN(); ++ii) {
                    
                    PairInt nbr = new PairInt(orderedOutput.getX(ii),
                        orderedOutput.getY(ii));
                    
                    if (nbrsSet.contains(nbr)) {
                        // this cannot be last point, because other neighbor
                        // would have been found before, so can assume
                        // the next point is one of the other neighbors
                        assert(ii < (orderedOutput.getN() - 1));
                    
                        int x2 = orderedOutput.getX(ii + 1);
                        int y2 = orderedOutput.getY(ii + 1);
                        
                        int diffX = Math.abs(x2 - nbr.getX());
                        int diffY = Math.abs(y2 - nbr.getY());
                        assert((diffX == 0 && diffY == 1) || 
                            (diffX == 1 && diffY == 0));
                        PairInt p2 = new PairInt(x2, y2);
                        
                        assert(nbrsSet.contains(p2));
                        
                        orderedOutput.insert(ii, nbr.getX(), nbr.getY());
                        
                        found = true;
                        break;
                    }
                }
                assert(found);
            }
        }

        return orderedOutput;    
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
        
        //System.out.println(
        //    String.format("x,y=(%d,%d) medAxis=(%d,%d)", x, y, 
        //        medAxisP.getX(), medAxisP.getY()));
        
        double minAngle = Double.MAX_VALUE;
        int minIdx = -1;
        
        //System.out.println(
        //    String.format("    (%d,%d) (%d,%d) (%d,%d) nXY=%d", 
        //        prevX, prevY, x, y, medAxisP.getX(), medAxisP.getY(), nXY));
        
        // if (x,y) == medAxisP or (x2,y2) == medAxisP
        // angle is NAN

        if (!(medAxisP.getX() == x && medAxisP.getY() == y)) {
            for (int i = 0; i < nXY; ++i) {
                int x2 = neighborsX[i];
                int y2 = neighborsY[i];

                double angle = LinesAndAngles.calcClockwiseAngle(
                    x, y, x2, y2, medAxisP.getX(), medAxisP.getY());

                if (Double.isNaN(angle)) {
                    if (x == x2 && x == medAxisP.getX()) {
                        if (y < y2) {
                            angle = 0;
                        } else {
                            angle = Math.PI;
                        }
                    } else if (y == y2 && y == medAxisP.getY()) {
                        if (x < x2) {
                            angle = Math.PI/2.;
                        } else {
                            angle = 3. * Math.PI/2.;
                        }
                    }
                }
               
                //System.out.println(String.format(
                //    "    (%d,%d) a=%.4f", 
                //    x2, y2, (float) angle));
                
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

                if (Double.isNaN(angle)) {
                    if (x == x2 && x == medAxisP.getX()) {
                        if (y < y2) {
                            angle = 0;
                        } else {
                            angle = Math.PI;
                        }
                    } else if (y == y2 && y == medAxisP.getY()) {
                        if (x < x2) {
                            angle = Math.PI/2.;
                        } else {
                            angle = 3. * Math.PI/2.;
                        }
                    }
                }
                
                //System.out.println(
                //    String.format("    *(%d,%d) a=%.4f", 
                //        x2, y2, (float) angle));
                
                
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
                xp, yp, xPolygon, yPolygon, "ordered bounds");

            n = contiguousShapePoints.size();
            xp = new int[n];
            yp = new int[n];
            int i = 0;
            for (PairInt p : contiguousShapePoints) {
                xp[i] = p.getX();
                yp[i] = p.getY();
                i++;
            }
            plotter.addPlot(xminmxyminmiac[0], xminmxyminmiac[1],
                xminmxyminmiac[2], xminmxyminmiac[3],
                xp, yp, xPolygon, yPolygon, "filled shape");

            plotter.writeFile2();

            System.out.println("output=" + output.toString());
        } catch (Throwable t) {

        }
    }
    
    private void debug(TIntSet contiguousShapePoints,
        PairIntArray output, int imgWidth) {
        PixelHelper ph = new PixelHelper();
        try {

            int[] xPolygon = null;
            int[] yPolygon = null;
            PolygonAndPointPlotter plotter = new PolygonAndPointPlotter();
            int[] minMaxXY = MiscMath.findMinMaxXY(
                contiguousShapePoints, imgWidth);
            int[] xp, yp;
            int n, count;

            n = output.getN();
            xp = new int[n];
            yp = new int[n];
            for (int i = 0; i < output.getN(); ++i) {
                xp[i] = output.getX(i);
                yp[i] = output.getY(i);
            }
            plotter.addPlot(minMaxXY[0], minMaxXY[1], minMaxXY[2], 
                minMaxXY[3], xp, yp, xPolygon, yPolygon, "ordered bounds");
            
            n = contiguousShapePoints.size();
            xp = new int[n];
            yp = new int[n];
            int[] xyout = new int[2];
            int i = 0;
            TIntIterator iter = contiguousShapePoints.iterator();
            while (iter.hasNext()) {
                int pixIdx = iter.next();
                ph.toPixelCoords(pixIdx, imgWidth, xyout);
                yp[i] = xyout[1];
                xp[i] = xyout[0];
                ++i;
            }
            plotter.addPlot(minMaxXY[0], minMaxXY[1], minMaxXY[2], minMaxXY[3],
                xp, yp, xPolygon, yPolygon, "filled shape");
            
            plotter.writeFile2();

            System.out.println("output=" + output.toString());
        } catch (Throwable t) {

        }
    }    

    private void debug(Set<PairInt> contiguousShapePoints,
        Set<PairInt> boundary, PairIntArray output) {
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
                xp, yp, xPolygon, yPolygon, "ordered bounds");
            
            n = boundary.size();
            xp = new int[n];
            yp = new int[n];
            int i = 0;
            for (PairInt p : boundary) {
                xp[i] = p.getX();
                yp[i] = p.getY();
                i++;
            }
            plotter.addPlot(xminmxyminmiac[0], xminmxyminmiac[1],
                xminmxyminmiac[2], xminmxyminmiac[3],
                xp, yp, xPolygon, yPolygon, "boundary");
            
            n = contiguousShapePoints.size();
            xp = new int[n];
            yp = new int[n];
            i = 0;
            for (PairInt p : contiguousShapePoints) {
                xp[i] = p.getX();
                yp[i] = p.getY();
                i++;
            }
            plotter.addPlot(xminmxyminmiac[0], xminmxyminmiac[1],
                xminmxyminmiac[2], xminmxyminmiac[3],
                xp, yp, xPolygon, yPolygon, "filled shape");
            
            plotter.writeFile2();

            System.out.println("output=" + output.toString());
        } catch (Throwable t) {

        }
    }

    private void debug(TIntSet contiguousShapePoints, TIntSet boundary, 
        PairIntArray output, int imgWidth) {
        
        try {

            int[] xPolygon = null;
            int[] yPolygon = null;
            PolygonAndPointPlotter plotter = new PolygonAndPointPlotter();
            int[] minMaxXY = MiscMath.findMinMaxXY(
                contiguousShapePoints, imgWidth);
            int[] xp, yp;
            int n, count;

            PixelHelper ph = new PixelHelper();
            
            n = output.getN();
            xp = new int[n];
            yp = new int[n];
            for (int i = 0; i < output.getN(); ++i) {
                xp[i] = output.getX(i);
                yp[i] = output.getY(i);
            }
            plotter.addPlot(minMaxXY[0], minMaxXY[1], minMaxXY[2], minMaxXY[3],
                xp, yp, xPolygon, yPolygon, "ordered bounds");
            
            n = boundary.size();
            xp = new int[n];
            yp = new int[n];
            int[] xyout = new int[2];
            int i = 0;
            TIntIterator iter = boundary.iterator();
            while (iter.hasNext()) {
                int pixIdx = iter.next();
                ph.toPixelCoords(pixIdx, imgWidth, xyout);
                yp[i] = xyout[1];
                xp[i] = xyout[0];
                
                ++i;
            }
            plotter.addPlot(minMaxXY[0], minMaxXY[1], minMaxXY[2], minMaxXY[3],
                xp, yp, xPolygon, yPolygon, "boundary");
            
            n = contiguousShapePoints.size();
            xp = new int[n];
            yp = new int[n];
            i = 0;
            iter = boundary.iterator();
            while (iter.hasNext()) {
                int pixIdx = iter.next();
                ph.toPixelCoords(pixIdx, imgWidth, xyout);
                yp[i] = xyout[1];
                xp[i] = xyout[0];
                ++i;
            }
            plotter.addPlot(minMaxXY[0], minMaxXY[1], minMaxXY[2], minMaxXY[3],
                xp, yp, xPolygon, yPolygon, "filled shape");
            
            String file = plotter.writeFile2();
            System.out.println("wrote debug file=" + file);
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

                if (Double.isNaN(angle)) {
                    if (x == x2 && x == medAxisP.getX()) {
                        if (y < y2) {
                            angle = 0;
                        } else {
                            angle = Math.PI;
                        }
                    } else if (y == y2 && y == medAxisP.getY()) {
                        if (x < x2) {
                            angle = Math.PI/2.;
                        } else {
                            angle = 3. * Math.PI/2.;
                        }
                    }
                }
                
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

                if (Double.isNaN(angle)) {
                    if (x == x2 && x == medAxisP.getX()) {
                        if (y < y2) {
                            angle = 0;
                        } else {
                            angle = Math.PI;
                        }
                    } else if (y == y2 && y == medAxisP.getY()) {
                        if (x < x2) {
                            angle = Math.PI/2.;
                        } else {
                            angle = 3. * Math.PI/2.;
                        }
                    }
                }
                
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

    private PairInt findSmallestXY(Set<PairInt> points) {

        PairInt minP = null;
        
        for (PairInt p : points) {
            if (minP == null) {
                minP = p;
            } else if (p.getX() < minP.getX()) {
                minP = p;
            } else if (p.getX() == minP.getX()) {
                if (p.getY() < minP.getY()) {
                    minP = p;
                }
            }
        }
        
        return minP;
    }
    
    private int findSmallestXY(TIntSet pixIdxs, int imgWidth) {

        if (pixIdxs.isEmpty()) {
            return -1;
        }
        
        int minPX = -1;
        int minPY = -1;
        PixelHelper ph = new PixelHelper();
        int[] xyout = new int[2];
        TIntIterator iter = pixIdxs.iterator();
        while (iter.hasNext()) {
            int pixIdx = iter.next();
            ph.toPixelCoords(pixIdx, imgWidth, xyout);
            int y = xyout[1];
            int x = xyout[0];                
            if (minPX == -1) {
                minPX = x;
                minPY = y;
            } else if (x < minPX) {
                minPX = x;
                minPY = y;
            } else if (x == minPX) {
                if (y < minPY) {
                    minPX = x;
                    minPY = y;
                }
            }
        }
        
        if (minPY == -1) {
            throw new IllegalStateException("no minimum was found");
        }
        
        int minP = (minPY * imgWidth) + minPX;
        
        return minP;
    }

    private VeryLongBitString[] createPointNeighborArray(
        List<PairInt> points, TObjectIntMap<PairInt> pointIndexes) {
        
        int n = points.size();
        
        VeryLongBitString[] out = new VeryLongBitString[n];
                
        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;
        
        for (int i = 0; i < n; ++i) {
            
            PairInt p = points.get(i);
            
            out[i] = new VeryLongBitString(n);
                        
            int x = p.getX();
            int y = p.getY();
                        
            for (int k = 0; k < dxs.length; ++k) {
                int x2 = x + dxs[k];
                int y2 = y + dys[k];
                PairInt p2 = new PairInt(x2, y2);
                if (pointIndexes.containsKey(p2)) {
                    out[i].setBit(pointIndexes.get(p2));
                }
            }
        }
        
        return out;
    }
    
    private VeryLongBitString[] createPointNeighborArray(
        int[] pixIdxs, TIntIntMap pointIndexes, int imgWidth, int imgHeight) {
        
        int n = pixIdxs.length;
        
        VeryLongBitString[] out = new VeryLongBitString[n];
                
        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;
        
        PixelHelper ph = new PixelHelper();
        int[] xyout = new int[2];
        
        for (int i = 0; i < n; ++i) {
            
            int pixIdx = pixIdxs[i];
            
            out[i] = new VeryLongBitString(n);
                  
            ph.toPixelCoords(pixIdx, imgWidth, xyout);
            int y = xyout[1];
            int x = xyout[0];
                        
            for (int k = 0; k < dxs.length; ++k) {
                int x2 = x + dxs[k];
                int y2 = y + dys[k];
                if (x2 < 0 || y2 < 0 || (x2 > (imgWidth - 1)) || 
                    (y2 > (imgHeight - 1))) {
                    continue;
                }
                int pixIdx2 = (y2 * imgWidth) + x2;
                if (pointIndexes.containsKey(pixIdx2)) {
                    out[i].setBit(pointIndexes.get(pixIdx2));
                }
            }
        }
        
        return out;
    }

    private TIntSet findJunctions(VeryLongBitString[] pointNeighbors) {

        /*
        looking for corner junctions.
        
        p4
           p3 p2
           p1
              p0
        p3 map has p1, p2, p4
        p1 map has p0, p2, p3
        the intersection of the 2 contains p2 which is the corner junction
        
        p4
           p3 p2
              p1
              p0
        p3 map has p1, p2, p4
        p1 map has p0, p2, p3
        the intersection contains p2 which is the corner junction
        
        */
        
        TIntSet out = new TIntHashSet();
        
        int n = pointNeighbors.length;
        for (int i = 0; i < n; ++i) {
            VeryLongBitString bs1 = pointNeighbors[i];
            for (int j = (i + 1); j < n; ++j) {
                VeryLongBitString bs2 = pointNeighbors[j];
                
                VeryLongBitString intersection = bs1.and(bs2);
                
                int[] setBits = intersection.getSetBits();
            
                if (setBits.length > 0) {
                    //System.out.println("i=" + i + " j=" + j +
                    //    " intersection=" + Arrays.toString(setBits));
                    out.addAll(setBits);
                }
            }
        }
        
        return out;
    }

    private boolean isAdjacent(PairInt p, int x, int y) {

        int diffX = Math.abs(p.getX() - x);
        int diffY = Math.abs(p.getY() - y);
        if ((diffX == 1 && diffY <= 1) || (diffY == 1 && diffX <= 1)) {
            return true;
        } 
        
        return false;
    }
    
    private boolean isAdjacent(int x0, int y0, int x, int y) {

        int diffX = Math.abs(x0 - x);
        int diffY = Math.abs(y0 - y);
        if ((diffX == 1 && diffY <= 1) || (diffY == 1 && diffX <= 1)) {
            return true;
        } 
        
        return false;
    }

    /**
     * NOT READY FOR USE...needs alot more testing.
     * 
     * given a set of contiguous points, fills in embedded points
     * and then extracts the outer boundary points, and then
     * orders the outer boundary points.
     * For best results, the contiguous points might need to be pre-processed
     * to smooth the curve and in some cases, one might want to separate the
     * region into more than one where the region is thin enough to prevent
     * two pixel paths through it as needed by a spatially sequential closed curve.
     * 
     * NOTE: it is up to the invoker to 
     * give the method a point set which is contiguous.
     * 
     * @param contiguousPoints
     * @return 
     */
    public PairIntArray extractOrderedBorder(TIntSet contiguousPoints, 
        int imgWidth, int imgHeight) throws Exception {
       
        PixelHelper ph = new PixelHelper();
        int[] xyout = new int[2];
        if (contiguousPoints.size() < 4) {
            PairIntArray output = new PairIntArray(contiguousPoints.size());
            TIntIterator iter = contiguousPoints.iterator();
            while (iter.hasNext()) {
                int pixIdx = iter.next();
                ph.toPixelCoords(pixIdx, imgWidth, xyout);
                int y = xyout[1];
                int x = xyout[0];
                output.add(x, y);
            }
            return output;
        }
                
        //O(8*N)
        TIntSet embedded = new TIntHashSet();
        TIntSet boundary = new TIntHashSet();
        extractBorder2(contiguousPoints, embedded, boundary, imgWidth);
        TIntSet set2 = new TIntHashSet(contiguousPoints);
        set2.addAll(embedded);
        
        if (boundary == null || boundary.size() == 0) {
            return null;
        }
        
        return orderTheBoundary(boundary, imgWidth, imgHeight);
    }
     
    public PairIntArray orderTheBoundary(TIntSet boundary, int imgWidth, 
        int imgHeight) throws Exception {
        
        PixelHelper ph = new PixelHelper();
        int[] xyout = new int[2];
        
        /*
        the algorithm finds the leftmost and smallest xy point in the boundary, 
        then traverses its unadded boundary neighbors to find the neighbor 
        with the smallest clockwise angle from it and previous point.
        For the case that an immediate unadded boundary neighbor is not present
        as is the case at the end of a single pixel wide spike, for example,
        the algorithm walks back up the output list looking for a point which
        has an unadded boundary neighbor and continues from there.
        
        note also, the junctions are found ahead of time and removed
        from the junctionSet and they are added back afterwrds to the output array.
        corners are a subset of junctions and they are excluded by the small
        clockwise angle ordering, but if they are not added on a closed
        curve return to that region, they are added later.
        
        NOTE: for best results, the user might want to have pre-processed the
        points to remove such features.
        */
        
        //using a list so can use bit vectors for quick set intersections
        int[] boundaryArray = boundary.toArray(new int[boundary.size()]);
                 
         // key = pixIdx
        TIntIntMap pointIndexes = new TIntIntHashMap();
     
        for (int i = 0; i < boundaryArray.length; ++i) {
            pointIndexes.put(boundaryArray[i], i);
        }
        
        // index = index of boundaryArray, item = bitstring with neighbors set
        //                              where bits are indexes in boundaryArray
        VeryLongBitString[] pointNeighbors = createPointNeighborArray(
            boundaryArray, pointIndexes, imgWidth, imgHeight);
        
        // corner junctions as indexes of boundaryArray
        TIntSet junctions = findJunctions(pointNeighbors);
  
        PairIntArray orderedOutput = new PairIntArray(boundary.size());
        
        TIntSet remaining = new TIntHashSet(boundary);
                
        int currPix = findSmallestXY(remaining, imgWidth);
        remaining.remove(currPix);
        ph.toPixelCoords(currPix, imgWidth, xyout);
        int currY = xyout[1];
        int currX = xyout[0];
        
        orderedOutput.add(currX, currY);
        junctions.remove(pointIndexes.get(currPix));
        
        // fake point for angle calc refs.  since curr is smallest x, should be
        // safe to make a point to left of it.
        int prevPixX = currX - 1;
        int prevPixY = currY;
        
        // idx w.r.t. orderedOutput.  used for walking back up array
        int outIdx = -1;
        
        while (!remaining.isEmpty()) {
       
            // among the neighbors of curr that are in remaining, chose
            // the one which has the smallest clockwise angle w/ prev and curr
            
            double minAngle = Double.MAX_VALUE;
            int minPX = -1;
            int minPY = -1;
            int minP = -1;
            
            int x = currX;
            int y = currY;
            int prevX = prevPixX;
            int prevY = prevPixY;
            
            // these are indexes in boundaryArray
            int[] currNbrs = pointNeighbors[pointIndexes.get(currPix)].getSetBits();
            
            for (int nbrIdx : currNbrs) {
                int pixIdx2 = boundaryArray[nbrIdx];
                
                if (!remaining.contains(pixIdx2)) {
                    continue;
                }
                ph.toPixelCoords(pixIdx2, imgWidth, xyout);
                int y2 = xyout[1];
                int x2 = xyout[0];        
                assert(!(prevX == x2 && prevY == y2));
                double angle = LinesAndAngles.calcClockwiseAngle(
                    prevX, prevY, x2, y2, x, y);
                if (Double.isNaN(angle)) {
                    // can occur is prev==p or p2==p
                    throw new IllegalStateException("error: unexpected NaN, "
                        + " due to equal points");
                }

                if (angle < minAngle) {
                    minAngle = angle;
                    minP = pixIdx2;
                    minPX = x2;
                    minPY = y2;
                }
                
                //System.out.format("(%d,%d) (%d,%d) (%d,%d) angle=%.4f  minP=%d\n",
                //    prevX, prevY, x, y, x2, y2, (float)angle, minP);
            }
            
            if (minP == -1) {
                
                // when remaining.size is >= (bounds.size - junctions.size)
                //   start checking for whether have reached the starting
                //   point, and look for whether the remaining points are
                //   junctions that need to be inserted in their 90 degree 
                //   positions.
                if (remaining.size() <= (boundaryArray.length - junctions.size())) {
                    if (isAdjacent(currX, currY, orderedOutput.getX(0), 
                        orderedOutput.getY(0))) {
                        
                        // remIdxs holds indexes of boudaryArray
                        TIntSet remIdxs = new TIntHashSet();
                        TIntIterator iter3 = remaining.iterator();
                        while (iter3.hasNext()) {
                            int pixIdx3 = iter3.next();
                            remIdxs.add(pointIndexes.get(pixIdx3));
                        }
                        remIdxs.removeAll(junctions);
                        if (remIdxs.isEmpty()) {
                            break;
                        }
                    }
                }
                
                // walk back up the output to try previous points, looking for
                //   one w/ an unassigned neighbor
                outIdx--;
                
                if (outIdx < 1) {
                    
                    int nj = junctions.size();
                    int nb = boundaryArray.length;
                    
                    debug(boundary, boundary, orderedOutput, imgWidth);
                    
                    System.out.println("output.n=" 
                        + orderedOutput.getN() 
                        + " rem.n=" + remaining.size() + " " + 
                        boundary.size());
                    
                    throw new Exception("Error in closed curve shape.");
                }
                
                currX = orderedOutput.getX(outIdx);
                currY = orderedOutput.getY(outIdx);
                currPix = (currY * imgWidth) + currX;
                
                prevPixX = orderedOutput.getX(outIdx - 1);
                prevPixY = orderedOutput.getY(outIdx - 1);
                       
                continue;
            }
            
            assert(minP != -1);
            
            prevPixX = currX;
            prevPixY = currY;
            currX = minPX;
            currY = minPY;
            currPix = (int)ph.toPixelIndex(currX, currY, imgWidth);
            
            outIdx = orderedOutput.getN();
            
            remaining.remove(currPix);
            orderedOutput.add(currX, currY);
            junctions.remove(pointIndexes.get(currPix));
        }
        
        if (!junctions.isEmpty()) {
            //insert between junction neigbhors
            //TODO: this could be improved, but for now, will use
            //  a pattern that has worse case runtime near
            //  O(junctions.size * output.size)
            TIntIterator iter = junctions.iterator();
            while (iter.hasNext()) {
                int jIdx = iter.next();
                
                int pixIdx = boundaryArray[jIdx];
                
                int[] nbrs = pointNeighbors[jIdx].getSetBits();
                TIntSet nbrsSet = new TIntHashSet();
                for (int nbrIdx : nbrs) {
                    nbrsSet.add(boundaryArray[nbrIdx]);
                }
                
                // add junction after first of it's neighbors it finds
                boolean found = false;
                for (int ii = 0; ii < orderedOutput.getN(); ++ii) {
                    
                    int nbrX = orderedOutput.getX(ii);
                    int nbrY = orderedOutput.getY(ii);
                    int nbr = (nbrY * imgWidth) + nbrX;
                    
                    if (nbrsSet.contains(nbr)) {
                        // this cannot be last point, because other neighbor
                        // would have been found before, so can assume
                        // the next point is one of the other neighbors
                        assert(ii < (orderedOutput.getN() - 1));
                    
                        int x2 = orderedOutput.getX(ii + 1);
                        int y2 = orderedOutput.getY(ii + 1);
                        
                        int diffX = Math.abs(x2 - nbrX);
                        int diffY = Math.abs(y2 - nbrY);
              
                        assert(diffX <= 1 && diffY <= 1);
                        
                        //int pixIdx2 = (int)ph.toPixelIndex(x2, y2, imgWidth);
                        //assert(nbrsSet.contains(pixIdx2));
                        
                        orderedOutput.insert(ii, nbrX, nbrY);
                        
                        found = true;
                        break;
                    }
                }
                assert(found);
            }
        }

        return orderedOutput;    
    }   
    
}
