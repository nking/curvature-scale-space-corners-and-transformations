package algorithms.imageProcessing;

import algorithms.CountingSort;
import algorithms.MultiArrayMergeSort;
import algorithms.QuickSort;
import algorithms.imageProcessing.util.MatrixUtil;
import algorithms.misc.Misc;
import algorithms.misc.MiscDebug;
import algorithms.util.PairIntArray;
import algorithms.util.PairInt;
import algorithms.util.PairIntArrayComparator;
import algorithms.util.PairIntArrayDescendingComparator;
import algorithms.util.PairIntArrayWithColor;
import algorithms.util.PointPairInt;
import java.io.IOException;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.Stack;
import java.util.logging.Level;
import thirdparty.HungarianAlgorithm;

/**
 * Edge extractor operates on an image that has already been reduced to 
 * single pixel width lines and extracts edges from it, attempting to make
 * the longest edges it can - this specific edge extractor stores and uses
 * junction information, so is not as fast as EdgeExtractor.java.
 * 
 Edge extraction
    Local Methods:
        (1) DFS walk through adjacent non-zero pixels to form a sequence of 
            pixels called an edge.
          
        (2) find join points
        
        (3) join edges using join points
        
        (4) find junction points
         
        (5) using the junction points, splice together the edges to improve
            the edges.  Currently, improve means to make the edges longer
            (which is good for an image expected to have one contour for 
             example, and doesn't harm corner finding for non-contour goals).
         
        (7) remove edges shorter than a minimum length (this may of may not
            be currently enabled)
        
  @see AbstractEdgeExtractor
  * 
 * @author nichole
 */
public class EdgeExtractorWithJunctions extends AbstractEdgeExtractor {
    
    private boolean debug = false;
    
    /**
     * map with key = center of junction pixel coordinates; 
     * value = set of adjacent pixels when there are more than the preceding 
     * and next.
     */
    private Map<Integer, Set<Integer>> junctionMap = new HashMap<Integer, Set<Integer>>();

    /**
     * map with key = pixel coordinates of all pixels involved in junctions;
     * value = PairInt holding index of edge that pixel is located in and
     * holding the index within that edge of the pixel.
     * for example, a pixel located in edges(0) at offset=100
     * would have PairInt(0, 100) as a value.
     */
    private Map<Integer, PairInt> junctionLocationMap = new HashMap<Integer, PairInt>();

    private int defaultMaxIterJunctionJoin = 1;
    
    private boolean singleClosedEdge = false;
    
    // ordered so that closest are first, then the diagonals
    private static int[] dxs8 = new int[]{-1, 0, 1,  0, -1,   1, 1, -1};
    private static int[] dys8 = new int[]{ 0, 1, 0, -1, -1,  -1, 1,  1};
    
    /**
     * NOTE:  input should have a black (empty) background and edges should
     * have values > 125 counts.  Edges should also have width of 1 and no larger.
     * 
     * @param input 
     */
    public EdgeExtractorWithJunctions(GreyscaleImage input) {
        
        super(input);
    }
    
    /**
     * NOTE:  input should have a black (empty) background and edges should
     * have values > 125 counts.  Edges should also have width of 1 and no larger.
     * The guide image is used to alter the extracted edges back towards the
     * highest intensity pixels of the guide image.  The guide image is expected
     * to be the combined X and Y gradient image from earlier processing
     * stages.
     * 
     * @param input
     * @param anEdgeGuideImage
     */
    public EdgeExtractorWithJunctions(GreyscaleImage input, 
        final GreyscaleImage anEdgeGuideImage) {
       
        super(input, anEdgeGuideImage);
    }
    
    public void overrideMaxNumberIterationsJunctionSplice(int nMaxIter) {
        if (nMaxIter < 2) {
            return;
        }
        defaultMaxIterJunctionJoin = nMaxIter;
    }
    
    /**
     * get a hashmap of the junctions within the edges, that is the regions
     * where there are more than 2 adjacent pixels to a pixel.
     * @return hashmap with key = image pixel indexes, values = set of adjacent
     * pixels as pixel indexes.
     */
    public Map<Integer, Set<Integer> > getJunctionMap() {
        return junctionMap;
    }
    
    /**
     * contains information needed to find any of the pixels that are part of 
     * a junction.
     * @return hashmap with key = image pixel indexes, values = a PairInt with
     * x holding the edge index and y holding the offset index of the point 
     * within the edge.
     */
    public Map<Integer, PairInt> getLocatorForJunctionAssociatedPoints() {
        return junctionLocationMap;
    }
    
    /**
     * find the edges and return as a list of points.  The method uses a
     * DFS search through all points in the image with values > 0 to link
     * adjacent sequential points into edges.
     * As a side effect, the method also populates
     * member variables edgeJunctionMap and outputIndexLocatorForJunctionPoints.
     * 
     * @return 
     */
    @Override
    public List<PairIntArray> findEdgesIntermediateSteps(List<PairIntArray> edges) {
        
        List<PairIntArray> output = edges;
    
        //O(N)
        Map<PairInt, PairInt> joinPoints = findJoinPoints(output);
        
        if (debug) {
            algorithms.misc.MiscDebug.printJoinPoints(joinPoints, output);
            algorithms.misc.MiscDebug.writeJoinPointsImage(joinPoints, output,
                img);
        }
        
        log.fine("edges.size()=" + output.size() + " before join-points");
        
        output = joinOnJoinPoints(joinPoints, output);
 
        log.fine("edges.size()=" + output.size() + " after join-points");
        
        int nMaxIter = defaultMaxIterJunctionJoin;
        int nIter = 0;
        int nSplices = 0;
        
        while ((nIter == 0) || ((nIter < nMaxIter) && (nSplices > 0))) {
        
            if (nSplices > 0) {
                removeEdgesShorterThan(output, 1);
            }
            
            findJunctions(output);
        
        if (output.size() > 10000) {
            return output;
        }
            
            if (debug) {
                algorithms.misc.MiscDebug.printJunctions(this.junctionMap, 
                    output, img);
            }
            
            //TODO: could improve thisfor nIter > 0 by only visiting the
            // edges which had changed (been spliced and fused on previous iter)

            nSplices = spliceEdgesAtJunctionsIfImproves(output);
            
            if (singleClosedEdge && ((nIter > (nMaxIter/2)) || (nSplices == 0))) {
                // because re-order and insertion both need corrected junction
                // maps, cannot invoke both on this iteration unless the first
                // did not alter anything.
                findJunctions(output);
                boolean ins = ((nIter % 1) == 0);
                if ((nIter % 1) == 1) {
                    int nReordered = reorderPointsInEdgesForClosedCurve(output);
                    nSplices += nReordered;
                    if (nReordered == 0) {
                        ins = true;
                    }
                }
                if (ins) {
                    int nIns = insertAdjacentForClosedCurve(output);
                    nSplices += nIns;
                }
            }
            
            ++nIter;
        }
                
        return output;
    }
    
    /**
     * Find the join points of edges which are adjacent to one another and to
     * no other edges.  This method is intended to follow the first creation
     * of edges from the line thinned edge image to find the unambiguously
     * connectable edges.  (Note that the method findJunctions(...) is later 
     * used to find where there are more than 2 adjacent pixels to any pixel.)
     * 
     * runtime complexity is O(N).
     * 
     * @param edges 
     * @return hashmap with key = edge index and index within edge of pixel;
     * value = the adjacent pixel's edge index and index within its edge
     */
    protected Map<PairInt, PairInt> findJoinPoints(List<PairIntArray> edges) {
           
        int n = edges.size();
        
        int w = img.getWidth();
        int h = img.getHeight();
        
        //with key = x, y of point
        //with value = edge index, index within edge
        Map<PairInt, PairInt> endPointMap = new HashMap<PairInt, PairInt>(2 * n);
 
        // holds x,y points of first and last points of the edges
        ArrayDeque<PairInt> endPointQueue = new ArrayDeque<PairInt>(2 * n);

        //with key = edge index, index within edge of one join point 
        //                 (the join point with the smaller edge number) 
        //       with value = edge index, index within edge of other join point
        Set<PointPairInt> theJoinPoints = new HashSet<PointPairInt>(2 * n);

        //with key = edge number
        //with value = set of the entries in theJoinPoints for which one has a 
        //first end point in this edge
        Map<Integer, Set<PointPairInt>> edgeFirstEndPointMap = 
            new HashMap<Integer, Set<PointPairInt>>();
    
        Map<Integer, Set<PointPairInt>>  edgeLastEndPointMap = 
            new HashMap<Integer, Set<PointPairInt>> ();

        // 2* O(N_edges)
        for (int edgeIdx = 0; edgeIdx < n; edgeIdx++) {
            
            PairIntArray edge = edges.get(edgeIdx);
            
            int nEdge = edge.getN();
            
            if (nEdge == 0) {
                continue;
            } 
            
            PairInt xy = new PairInt(edge.getX(0), edge.getY(0));
            PairInt loc = new PairInt(edgeIdx, 0);
            
            endPointMap.put(xy, loc);
            endPointQueue.add(xy);
            
            if (nEdge == 1) {
                continue;
            }
            
            xy = new PairInt(edge.getX(nEdge - 1), edge.getY(nEdge - 1));
            loc = new PairInt(edgeIdx, nEdge - 1);
            
            endPointMap.put(xy, loc);
            endPointQueue.add(xy);
            
            if (nEdge == 2) {
                continue;
            }
            
            // add the second and second to last to endPointMap only
            xy = new PairInt(edge.getX(1), edge.getY(1));
            loc = new PairInt(edgeIdx, 1);
            endPointMap.put(xy, loc);
            
            if (nEdge == 3) {
                continue;
            }
            
            xy = new PairInt(edge.getX(nEdge - 2), edge.getY(nEdge - 2));
            loc = new PairInt(edgeIdx, nEdge - 2);
            endPointMap.put(xy, loc);
        }
        
        assert(endPointQueue.size() >= n);
        assert(endPointMap.size() >= 2*n);
        
        while (!endPointQueue.isEmpty()) {
            
            PairInt uXY = endPointQueue.poll();
            int uX = uXY.getX();
            int uY = uXY.getY();
            
            PairInt uLoc = endPointMap.get(uXY);
            
            assert(uLoc != null);
            
            int uEdgeIdx = uLoc.getX();
            
            Set<Integer> neighborEdgeNumbers = new HashSet<Integer>();

            PairInt vClosestXY = null;
            PairInt vClosestLoc = null;
            int closestDistSq = Integer.MAX_VALUE;
            int vClosestCanBeReordered = -99;
            
            for (int nIdx = 0; nIdx < dxs8.length; nIdx++) {
                int vX = uX + dxs8[nIdx];
                int vY = uY + dys8[nIdx];
                if ((vX < 0) || (vX > (w - 1)) || (vY < 0) || (vY > (h - 1))) {
                    continue;
                }
                PairInt vXY = new PairInt(vX, vY);
                
                PairInt vLoc = endPointMap.get(vXY);
                
                if ((vLoc == null) || (vLoc.getX() == uEdgeIdx)) {
                    continue;
                }
                                   
                int vEdgeIdx = vLoc.getX();

                neighborEdgeNumbers.add(Integer.valueOf(vEdgeIdx));

                if (neighborEdgeNumbers.size() > 1) {
                    // this is a junction, so break and continue at next uXY
                    break;
                }
                    
                PairIntArray vEdge = edges.get(vEdgeIdx);

                int canBeReordered = canBeReordered(vLoc, vEdge);

                if (canBeReordered == -1) {
                    continue;
                }

                int diffX = uX - vX;
                int diffY = uY - vY;
                int distSq = (diffX * diffX) + (diffY * diffY);

                if ((distSq < closestDistSq) || ((distSq == closestDistSq) 
                    && (canBeReordered == 0))) {

                    closestDistSq = distSq;
                    vClosestLoc = vLoc;
                    vClosestXY = vXY;
                    vClosestCanBeReordered = canBeReordered;
                }
            }
           
            if ((neighborEdgeNumbers.size() != 1) || (vClosestXY == null)) {
                continue;
            }
            
            // make the smaller edge number as variables '0'
            PairInt xy0 = uXY;
            PairInt xy1 = vClosestXY;
            PairInt loc0 = uLoc;
            PairInt loc1 = vClosestLoc;
            boolean smallestEdgeIdxIsV = false;
            if (loc0.getX() > loc1.getX()) {
                PairInt swap = loc0;
                loc0 = loc1;
                loc1 = swap;
                smallestEdgeIdxIsV = true;
                swap = xy0;
                xy0 = xy1;
                xy1 = swap;
            }
            int n0 = edges.get(loc0.getX()).getN();
            int n1 = edges.get(loc1.getX()).getN();
            
            theJoinPoints.add(new PointPairInt(loc0, loc1));
            
            // add the location points to the edgeFirstEndPointMap and 
            // edgeLastEndPointMap maps
            for (int ii = 0; ii < 2; ++ii) {
                int locX = loc0.getX();
                int locY = loc0.getY();
                int nEdge = n0;
                
                if (ii == 1) {
                    locX = loc1.getX();
                    locY = loc1.getY();
                    nEdge = n1;
                }
                
                boolean addToFirst = (locY < 2);
                boolean addToLast = (locY > (nEdge - 3));
                if (nEdge == 2) {
                    if (locY == 1) {
                        addToFirst = false;
                        addToLast = true;
                    } else {
                        addToFirst = true;
                        addToLast = false;
                    }
                }
            
                if ((ii == 0) && smallestEdgeIdxIsV && (vClosestCanBeReordered == 1) && (n0 == 3)) {
                    log.fine(String.format(
                    "creating a join point for an edge with size 3: xy=(%d,%d) loc=%d:%d",
                    vClosestLoc.getX(), vClosestLoc.getY(), vClosestLoc.getX(), 
                    vClosestLoc.getY()));
                    // an entry is added to both edge maps because we don't know yet
                    // which to use.
                }

                if (addToFirst) {
                    Integer key = Integer.valueOf(locX);
                    Set<PointPairInt> joinPointSet = edgeFirstEndPointMap.get(key);
                    if (joinPointSet == null) {
                        joinPointSet = new HashSet<PointPairInt>();
                    }
                    PointPairInt joinPoint = new PointPairInt(loc0, loc1);
                    joinPointSet.add(joinPoint);
                    edgeFirstEndPointMap.put(key, joinPointSet);
                }
            
                if (addToLast) {
                    Integer key = Integer.valueOf(locX);
                    Set<PointPairInt> joinPointSet = edgeLastEndPointMap.get(key);
                    if (joinPointSet == null) {
                        joinPointSet = new HashSet<PointPairInt>();
                    }
                    PointPairInt joinPoint = new PointPairInt(loc0, loc1);
                    joinPointSet.add(joinPoint);
                    edgeLastEndPointMap.put(key, joinPointSet);
                }
            }
        }
       
        reduceMultipleEndpointsForEdge(edges, edgeFirstEndPointMap,
            edgeLastEndPointMap, endPointMap, theJoinPoints);

        if (debug) {
            // 2 of the same join points for edges of size 3 in which 
            // the join point is movable to the first or last point in edge.
            // the next block removes that possibility
            boolean skipForSize3 = true;
            algorithms.misc.MiscDebug.assertConsistentJoinPointStructures(edges, 
                edgeFirstEndPointMap, edgeLastEndPointMap, theJoinPoints, 
                skipForSize3);
        }        
        
        // edges that are size 3 and have join points that are not in the
        // first or last position, can be moved to either, so keep
        // track of them to avoid moving to the same position twice
        Map<PairInt, Set<Integer>> edgeSize3NonEndPoints = new HashMap<PairInt, 
            Set<Integer>>();
        for (PointPairInt entry : theJoinPoints) {
            PairInt loc0 = entry.getKey();
            PairInt loc1 = entry.getValue();            
            int n0 = edges.get(loc0.getX()).getN();
            int n1 = edges.get(loc1.getX()).getN();
            if ((n0 == 3) && (loc0.getY() == 1)) {
                if (!edgeSize3NonEndPoints.containsKey(loc0)) {
                    edgeSize3NonEndPoints.put(loc0, new HashSet<Integer>());
                }
            }
            if ((n1 == 3) && (loc1.getY() == 1)) {
                if (!edgeSize3NonEndPoints.containsKey(loc1)) {
                    edgeSize3NonEndPoints.put(loc1, new HashSet<Integer>());
                }
            }
        }
        
        // 2 * O(N_join_points)
        // re-order the endpoints which are not first or last position in edge        
        for (PointPairInt entry : theJoinPoints) {
            PairInt loc0 = entry.getKey();
            PairInt loc1 = entry.getValue();
            
            PairIntArray edge0 = edges.get(loc0.getX());
            PairIntArray edge1 = edges.get(loc1.getX());
            
            PairInt origLoc0 = new PairInt(loc0.getX(), loc0.getY());
            
            int swapIdx = -1; 
            if (edgeSize3NonEndPoints.containsKey(loc0)) {
                swapIdx = 0;
                Set<Integer> takenPositions = edgeSize3NonEndPoints.get(loc0);
                if (takenPositions.size() > 1) {
                    throw new IllegalStateException(
                        "then are more than 2 join points for the same edge of size 3." + 
                        " edge " + loc0.getX() + " xy=(" + edge0.getX(loc0.getY())
                        + "," + edge0.getY(loc0.getY()) +")");
                }
                if (!takenPositions.isEmpty()) {
                    swapIdx = 2;
                }
                takenPositions.add(Integer.valueOf(swapIdx));
            } else  {
                int canBeReordered = canBeReordered(loc0, edge0);
                if (canBeReordered == 1) {
                    if (loc0.getY() < 2) {
                        swapIdx = 0;
                    } else {
                        swapIdx = edge0.getN() - 1;
                    }
                }
            }
            if (swapIdx > -1) {                
                int swapX = edge0.getX(swapIdx);
                int swapY = edge0.getY(swapIdx);
                edge0.set(swapIdx, edge0.getX(1), edge0.getY(1));
                edge0.set(1, swapX, swapY);
                loc0.setY(swapIdx);
            }
            
            swapIdx = -1; 
            if (edgeSize3NonEndPoints.containsKey(loc1)) {
                swapIdx = 0;
                Set<Integer> takenPositions = edgeSize3NonEndPoints.get(loc1);
                if (takenPositions.size() > 1) {
                    throw new IllegalStateException(
                        "then are more than 2 join points for the same edge of size 3." + 
                        " edge " + loc1.getX() + " xy=(" + edge1.getX(loc1.getY())
                        + "," + edge1.getY(loc1.getY()) +")");
                }
                if (!takenPositions.isEmpty()) {
                    swapIdx = 2;
                }
                takenPositions.add(Integer.valueOf(swapIdx));
            } else  {
                int canBeReordered = canBeReordered(loc1, edge1);
                if (canBeReordered == 1) {
                    if (loc1.getY() < 2) {
                        swapIdx = 0;
                    } else {
                        swapIdx = edge1.getN() - 1;
                    }
                }
            }
            if (swapIdx > -1) {                
                int swapX = edge1.getX(swapIdx);
                int swapY = edge1.getY(swapIdx);
                edge1.set(swapIdx, edge1.getX(1), edge1.getY(1));
                edge1.set(1, swapX, swapY);
                loc1.setY(swapIdx);
            }
        }
 
        if (debug) {
            boolean skipForSize3 = false;
            algorithms.misc.MiscDebug.assertConsistentJoinPointStructures(edges, 
                edgeFirstEndPointMap, edgeLastEndPointMap, theJoinPoints, 
                skipForSize3);
        }
        
        Map<PairInt, PairInt> result = new HashMap<PairInt, PairInt>();
        for (PointPairInt entry : theJoinPoints) {
            PairInt loc0 = entry.getKey();
            PairInt loc1 = entry.getValue();
            result.put(loc0, loc1);
        }
        
        return result;       
    }
    
    /**
     * Iterate over each point looking for its neighbors and noting when
     * there are more than 2, and store the found junctions as member variables.
     * 
     * runtime complexity is O(N)
     * 
     * @param edges 
     */
    protected void findJunctions(List<PairIntArray> edges) {
        
        // key = center of junction pixel coordinates
        // value = set of adjacent pixels when there are more than the preceding
        //         and next.
        Map<Integer, Set<Integer> > theJunctionMap = new HashMap<Integer, Set<Integer>>();

        // key = pixel coordinates of all pixels involved in junctions
        // value = PairInt holding index of edge that pixel is located in and
        //         holding the index within that edge of the pixel.
        //         for example, a pixel located in edges(0) at offset=100
        //         would have PairInt(0, 100) as a value.
        Map<Integer, PairInt> theJunctionLocationMap = new HashMap<Integer, PairInt>();

        int n = edges.size();
        
        // key = image pixel index, 
        // value = pairint of edge index, and index within edge
        Map<Integer, PairInt> pointLocator = new HashMap<Integer, PairInt>();
        
        // key = center of junction pixel coordinates
        // value = number of times this point is a value in the first theJunctionMap
        Map<Integer, Integer> theJunctionFrequencyMap = new HashMap<Integer, Integer>();
        
        // O(N)
        for (int edgeIdx = 0; edgeIdx < n; edgeIdx++) {
            
            PairIntArray edge = edges.get(edgeIdx);
            
            for (int uIdx = 0; uIdx < edge.getN(); uIdx++) {
                
                int pixIdx = img.getIndex(edge.getX(uIdx), edge.getY(uIdx));
                
                pointLocator.put(Integer.valueOf(pixIdx), new PairInt(edgeIdx, 
                    uIdx));
            }
        }
        
        int w = img.getWidth();
        int h = img.getHeight();
        
        // 8 * O(N)
        for (int edgeIdx = 0; edgeIdx < n; edgeIdx++) {
            
            PairIntArray edge = edges.get(edgeIdx);
            
            for (int iEdgeIdx = 0; iEdgeIdx < edge.getN(); iEdgeIdx++) {
                
                int col = edge.getX(iEdgeIdx);
                int row = edge.getY(iEdgeIdx);
                
                int uIdx = img.getIndex(col, row);
                
                Set<PairInt> neighbors = new HashSet<PairInt>();
                                
                for (int nIdx = 0; nIdx < dxs8.length; nIdx++) {
                    
                    int x = col + dxs8[nIdx];
                    int y = row + dys8[nIdx];
                    
                    if ((x < 0) || (x > (w - 1)) || (y < 0) || (y > (h - 1))) {
                        continue;
                    }
                    
                    int vIdx = img.getIndex(x, y);
                    
                    PairInt vLoc = pointLocator.get(Integer.valueOf(vIdx));
                    
                    if (vLoc != null) {
                        neighbors.add(vLoc);
                    }
                }

                if (neighbors.size() > 2) {
                    
                    Set<Integer> indexes = new HashSet<Integer>();
                                        
                    for (PairInt p : neighbors) {
                        
                        int edge2Idx = p.getX();
                        int iEdge2Idx = p.getY();
                        
                        PairIntArray vEdge = edges.get(edge2Idx);
                        
                        int vIdx = img.getIndex(vEdge.getX(iEdge2Idx), 
                            vEdge.getY(iEdge2Idx));
                        
                        Integer key = Integer.valueOf(vIdx);
                        
                        theJunctionLocationMap.put(key, p);
                        
                        indexes.add(key);
                        
                        Integer freq = theJunctionFrequencyMap.get(key);
                        if (freq == null) {
                            theJunctionFrequencyMap.put(key, Integer.valueOf(1));
                        } else {
                            theJunctionFrequencyMap.put(key, 
                                Integer.valueOf(freq.intValue() + 1));
                        }
                    }
                                        
                    theJunctionMap.put(Integer.valueOf(uIdx), indexes);
                    
                    theJunctionLocationMap.put(Integer.valueOf(uIdx), 
                        new PairInt(edgeIdx, iEdgeIdx));
                    
                } 
                // if (neighbors.size() == 2) is handled in findJoinPoints
            }
        }
        
        // visit each junction and make sure the real center is the only one
        // stored
        Set<Integer> remove = new HashSet<Integer>();
        Set<Integer> doNotRemove = new HashSet<Integer>();
        for (Entry<Integer, Set<Integer>> entry : theJunctionMap.entrySet()) {
            
            Integer pixelIndex = entry.getKey();
            
            Set<Integer> neighborIndexes = entry.getValue();
            int nN = neighborIndexes.size();
            
            int maxN = Integer.MIN_VALUE;
            
            for (Integer pixelIndex2 : neighborIndexes) {
                
                if (remove.contains(pixelIndex2)) {
                    continue;
                }
                            
                Set<Integer> neighborIndexes2 = theJunctionMap.get(pixelIndex2);
                
                if (neighborIndexes2 == null) {
                    continue;
                }
                                
                if (neighborIndexes2.size() > nN) {
                    if (neighborIndexes2.size() > maxN) {
                        maxN = neighborIndexes2.size();
                    }
                }
            }
            
            if (maxN > Integer.MIN_VALUE) {

                if (!doNotRemove.contains(pixelIndex)) {
                    remove.add(pixelIndex);
                }
                
            } else {
                doNotRemove.add(pixelIndex);
                // remove the neighbors from the junction map
                for (Integer pixelIndex2 : neighborIndexes) {
                    int mapFreq2 = theJunctionFrequencyMap.get(pixelIndex2);
                    
                    if (!doNotRemove.contains(pixelIndex2) && (mapFreq2 == 1)) {
                        remove.add(pixelIndex2);
                    }
                }
            }
        }
        
        for (Integer pixelIndex : remove) {
            theJunctionMap.remove(pixelIndex);
        }
        
        junctionMap = theJunctionMap;
        junctionLocationMap = theJunctionLocationMap;
    }
    
    private Map<Integer, PairIntArray> createIndexedMap(List<PairIntArray> edges) {
        
        Map<Integer, PairIntArray> output = new HashMap<Integer, PairIntArray>();
        for (int i = 0; i < edges.size(); i++) {
            output.put(Integer.valueOf(i), edges.get(i));
        }
           
        return output;
    }
    
    /**
     * join edges using information in the member variable joinPoints
     * and update the junction and joinPoint information after the
     * changes.  Note that the join points must be first or last positions
     * in an edge in edges.
     * @param joinPoints map with key = the PairInt with x being the edge index
     * and y being the index of the first join point within the edge, value =
     * the PairInt with x being the edge index
     * and y being the index of the other join point within the edge
     * @param edges
     * @return 
     */
    protected List<PairIntArray> joinOnJoinPoints(Map<PairInt, PairInt> 
        joinPoints, List<PairIntArray> edges) {
        
        //order the join points by the edge number of the first join point
        
        int[] indexes = new int[joinPoints.size()];
        PairInt[][] edgeJoins = new PairInt[joinPoints.size()][2];
        int count = 0;
        for (Entry<PairInt, PairInt> entry : joinPoints.entrySet()) {
            PairInt loc0 = entry.getKey();
            PairInt loc1 = entry.getValue();
            indexes[count] = loc0.getX();
            edgeJoins[count] = new PairInt[]{loc0, loc1};
            count++;
        }
        QuickSort.sortBy1stArg(indexes, edgeJoins);
        
long nPointsBefore = countPixelsInEdges(edges);
        
        // put the edges in a map to remove and search updates faster
        Map<Integer, PairIntArray> edgesMap = createIndexedMap(edges);
                
        int n = edgeJoins.length;
        
        for (int i = (n - 1); i > -1; --i) {
           
//algorithms.misc.MiscDebug.printJoinPoints(edgeJoins, 0, i, edgesMap);
           
            PairInt[] entry = edgeJoins[i];

            PairInt loc0 = entry[0];
            
            PairInt loc1 = entry[1];
            
            // if they've already been merged
            if (loc0.getX() == loc1.getX()) {
                continue;
            }
            
            // the smaller edge index should be loc1 and the larger loc0,
            if (loc1.getX() > loc0.getX()) {
                PairInt tmp = loc0;
                loc0 = loc1;
                loc1 = tmp;
            }
            final int edge0Idx = loc0.getX();
            final int edge1Idx = loc1.getX();
            
            // edge to move
            int removedEdgeIdx = edge0Idx;
            PairIntArray edge0 = edgesMap.remove(Integer.valueOf(removedEdgeIdx));

            int n0 = edge0.getN();
            
            if ((loc0.getY() != 0) && (loc0.getY() != (n0 - 1))){
                /*
                this can happen if there is an error in the join point algorithm
                resulting in the same join point to 2 different edges.  an update
                in the location will eventually be an error in one of them.
                
                After more testing, will change this to discard the join point
                and warn of error.
                */
                throw new IllegalStateException("ERROR in the updates? " + 
                    " loc0=" + edge0Idx + "," + loc0.getY() + " n=" + n0 +
                    " i=" + i + " (nedges=" + n + ") to append edge " 
                    + edge0Idx + " to edge " + edge1Idx);
            }
            
            // join point should be at the beginning, so reverse if not
            if (loc0.getY() == (n0 - 1)) {
                
                edge0.reverse();
                loc0.setY(0);
                
                // everything with smaller index than i in edgeJoins with
                // edgeIndex==loc0.getX() needs to be updated for this reversal.
                // idx becomes n-idx-1
                for (int j = (i - 1); j > -1; --j) {
                    PairInt[] vEntry = edgeJoins[j];
                    for (int k = 0; k < 2; k++) {
                        PairInt vLoc = vEntry[k];
                        if (vLoc.getX() == edge0Idx) {
                            int idxRev = n0 - vLoc.getY() - 1;
                            vLoc.setY(idxRev);
                        }
                    }
                }
            }
            // edge to receive new edge
            PairIntArray edge1 = edgesMap.get(Integer.valueOf(loc1.getX()));
            int n1 = edge1.getN();
            // endpoint should be at the end, so reverse if not
            if (loc1.getY() == 0) {
                
                edge1.reverse();
                loc1.setY(n1 - 1);
                
                // everything with smaller index than i in edgeJoins that has
                // edgeIndex==loc1.getX() needs to be updated for this reversal.
                // idx becomes n-idx-1
                for (int j = (i - 1); j > -1; --j) {
                    PairInt[] vEntry = edgeJoins[j];
                    for (int k = 0; k < 2; k++) {
                        PairInt vLoc = vEntry[k];
                        if (vLoc.getX() == edge1Idx) {
                            int idxRev = n1 - vLoc.getY() - 1;
                            vLoc.setY(idxRev);
                        }
                    }
                }
            }
            
            final int indexWithinEdge0 = loc0.getY();
            final int indexWithinEdge1 = loc1.getY();
            
            // --- append edge0 to edge1 ----
            edge1.addAll(edge0);
            
            // for earlier items in array edgeJoins
            // need to update all edge indexes and indexes within edge.
            
            // loc0 got appended to loc1 --> [edge1][edge0]  
            // points in edge0 need n1 added to them
            
            // first, will only update the indexes within the edge for 
            // edgeIndex == loc0.getX()
            //     the new offset index = n0 + index
            for (int j = (i - 1); j > -1; --j) {
                PairInt[] vEntry = edgeJoins[j];
                for (int k = 0; k < 2; k++) {
                    PairInt vLoc = vEntry[k];
                    if (vLoc.getX() == edge0Idx) {
                        int idxEdit = vLoc.getY() + n1;
                        vLoc.setY(idxEdit);
                    }
                }
            }
            
            // any edge index > loc0 gets reduced by one, but if == it gets replaced by loc1.getX()
            for (int j = (i - 1); j > -1; --j) {
                PairInt[] vEntry = edgeJoins[j];
                for (int k = 0; k < 2; k++) {
                    PairInt vLoc = vEntry[k];
                    if (vLoc.getX() > edge0Idx) {
                        int editX = vLoc.getX() - 1;
                        vLoc.setX(editX);
                    } else if (vLoc.getX() == edge0Idx) {
                        int editX = edge1Idx;
                        vLoc.setX(editX);
                    }
                }
            }

            // the output map keeps 0 to loc1.getx(),
            // but loc0.getX() + 1 gets moved to loc0.getX()
            //    and on until have reached size() - 1
            for (int j = (removedEdgeIdx + 1); j <= edgesMap.size(); ++j) {
                PairIntArray v = edgesMap.remove(Integer.valueOf(j));
                assert(v != null);
                edgesMap.put(Integer.valueOf(j - 1), v);
            }
       
        }
        
        List<PairIntArray> output = new ArrayList<PairIntArray>();
        for (int i = 0; i < edgesMap.size(); i++) {
            Integer key = Integer.valueOf(i);
            PairIntArray v = edgesMap.get(key);
            assert(v != null);
            output.add(v);
        }
        
        return output;
    }
    
    /**
     * return whether the point can be swapped with the nearest endpoint.
     * <pre>
     * returns 
     *     0: for no need to reorder
     *     1: for it not being an endpoint but can be re-ordered
     *     -1: for cannot be re-ordered
     * </pre>
     * @param pointEdgeLocation
     * @param edge
     * @return code 
       <pre> returns 
           0: for no need to reorder
           1: for it not being an endpoint but can be re-ordered
           -1: for cannot be re-ordered,
       </pre>
     */
    private int canBeReordered(PairInt pointEdgeLocation, PairIntArray edge) {
        
        int n = edge.getN();
        
        if ((pointEdgeLocation.getY() == 0) || (pointEdgeLocation.getY() == (n - 1))) {
            return 0;
        }
        
        if ((pointEdgeLocation.getY() > 2) && (pointEdgeLocation.getY() < (n - 3))) {
            return -1;
        }
        
        int idxOrig = pointEdgeLocation.getY();
        
        int check0Idx = 0;
        int check1Idx = 2;
        
        if ((n > 3) && (idxOrig == (n - 2))) {
            check0Idx = n - 3;
            check1Idx = n - 1;
        }
        
        int diffX = edge.getX(check0Idx) - edge.getX(check1Idx);
        int diffY = edge.getY(check0Idx) - edge.getY(check1Idx);
        if ((Math.abs(diffX) > 1) || (Math.abs(diffY) > 1)) {
            return -1;
        }
        
        return 1;            
    }
    
    /**
     * If pointEdgeLocation is near the beginning or end of edge, 
     * swap it with the point that is the endpoint and update the given 
     * data structures with the updated information.
     * @param pointEdgeLocation
     * @param edge
     * @return 1 for did re-order, 0 for no need to re-order, else -1 for cannot
     * re-order
     */
    private int reorderIfNearEnd(PairInt pointEdgeLocation, 
        PairIntArray edge,
        PairInt connectingEdgeLocation, PairIntArray connectingEdge) {
        
        int n = edge.getN();
        
        if ((pointEdgeLocation.getY() == 0) || (pointEdgeLocation.getY() == (n - 1))) {
            return 0;
        }
        
        /*
        looks like max offset from end is 2
         @ . [....]

         @ .
         . [....]
        */

        if ((pointEdgeLocation.getY() > 2) && (pointEdgeLocation.getY() < (n - 3))) {
            return 0;
        }
                            
        int idxOrig = pointEdgeLocation.getY();
        int idxSwap;
        if (idxOrig == 1) {
            idxSwap = 0;
        } else if (idxOrig == (n - 2)) {
            idxSwap = n - 1;
        } else {            
            /*
            can just swap it with first or last point if the swapped point
            is adjacent to the point right before the point to be swapped,
            in other words, does not break a connection in the edge
            */
            if (idxOrig == 2) {
                idxSwap = 0;
            } else {
                idxSwap = n - 1;
            }
            int prevX = edge.getX(idxOrig - 1);
            int prevY = edge.getY(idxOrig - 1);
            int endX = edge.getX(0);
            int endY = edge.getY(0);
            if ((Math.abs(prevX - endX) > 1) || (Math.abs(prevY - endY) > 1)) {
                return -1;
            }
        }
        
        /*
        if swap position is still adjacent to the connecting point,
        can complete the change.
        */
        int connectedX = connectingEdge.getX(connectingEdgeLocation.getY());
        int connectedY = connectingEdge.getY(connectingEdgeLocation.getY());
        
        int swapX = edge.getX(idxSwap);
        int swapY = edge.getY(idxSwap);
        
        if ((Math.abs(connectedX - swapX) > 1) || (Math.abs(connectedY - swapY) > 1)) {
            return -1;
        }
 
        pointEdgeLocation.setY(idxSwap);

        return 1;            
    }

    private int spliceEdgesAtJunctionsIfImproves(List<PairIntArray> edges) {
                       
        /*        
        The main goal is to make better contours.
        
        Edges with junctions can sometimes be spliced and rejoined with another
        junction edge to make a better edge where better edge may be a longer
        edge or a closed contour useful for determining transformations
        between images containing the contour.
        
        could be used with shape templates...
        
        For now, will make a simple algorithm which tries to increase the
        length of the longest edges in a junction.
        */
        
        int nSplices = 0;
        
        // key = edge index.  
        // value = pixel indexes.
        //   the pixel indexes are used to find values in junctionLocatorMap
        //   to update it as points are moved to and from edges.
        // key = edge index.  
        // value = pixel indexes.
        //   the pixel indexes are used to find values in junctionLocatorMap
        //   to update it as points are moved to and from edges.
        Map<Integer, Set<Integer>> theEdgeToPixelIndexMap = createEdgeToPixelIndexMap();
        
        if (debug) {
            algorithms.misc.MiscDebug.assertConsistentEdgeCapacity(
                theEdgeToPixelIndexMap, junctionLocationMap, junctionMap, 
                edges);
        }
        
        /*
        TODO: consider handling splices of edges of size 1 after this
        block.  They can be inserted in a way that does not affect line
        continuity
        
        For example:
        
            0
             1 < where '<' can be inserted after 1 and before 2 or vice versa
             2
             3
        */
        
        for (Entry<Integer, Set<Integer>> entry : junctionMap.entrySet()) {
                        
            Integer centerPixelIndex = entry.getKey();
            PairInt centerLoc = junctionLocationMap.get(centerPixelIndex);
            assert(centerLoc != null);

            /*
            if (debug) {
                
                log.info("processing junction w/ center pixel index=" + centerPixelIndex
                    + " and loc=" + centerLoc.getX() + ":" + centerLoc.getY());

                String str = algorithms.misc.MiscDebug.printJunctionsToString(
                    junctionLocationMap, junctionMap, edges, img);
                log.info(str);
                
            }
            */
            
            Set<Integer> adjIndexes = entry.getValue();
            
            int[] pixIndexes = new int[adjIndexes.size() + 1];
            
            // lengths holds the edge up until the junction.  splices the edge
            // at the junction figuratively and counts number of pixels before
            // and after splice and keeps the longest.
            int[] lengths = new int[pixIndexes.length];
            
            pixIndexes[0] = centerPixelIndex.intValue();
            lengths[0] = (new Splice(edges.get(centerLoc.getX())))
                .getLengthOfLongestSide(centerLoc.getY());
            
            int maxN = lengths[0];
            int maxNPixIdx = pixIndexes[0];
            boolean foundAnotherEdge = false;
            
            int count = 1;
            
            // for junctions on the same edge, need to count the portions
            // of their edge which are on oppossite sides of the junction,
            // so need to track the maxN, then correct the lengths of splice
            // from same edge afterwards
            
            for (Integer pixIndex : adjIndexes) {
                
                pixIndexes[count] = pixIndex.intValue();
                
                PairInt loc = junctionLocationMap.get(pixIndex);
            
                if ((centerLoc.getX() == loc.getX())) {
                    count++;
                    continue;
                }
               
                foundAnotherEdge = true;
                
                Splice splice = new Splice(edges.get(loc.getX()));
                                                
                lengths[count] = splice.getLengthOfLongestSide(loc.getY());
           
                if (lengths[count] > maxN) {
                    maxN = lengths[count];
                    maxNPixIdx = pixIndex;
                }
                
                count++;
            }
            
            if (!foundAnotherEdge) {
                continue;
            }

            //-- correct the lengths for the junctions from edge maxNEdgeIdx --
            PairInt maxNLoc = junctionLocationMap.get(maxNPixIdx);
            Splice maxNSplice = new Splice(edges.get(maxNLoc.getX()));
            int maxNSpliceOtherSideN = maxNSplice.splice(new int[]{maxNLoc.getY()})[1].getN();
            // correct the splice lengths if they are from the same edge
            for (int ii = 0; ii < count; ++ii) {
                int pixIdx = pixIndexes[ii];
                PairInt loc = junctionLocationMap.get(pixIdx);
                if ((loc.getX() != maxNLoc.getX()) || (pixIdx == maxNPixIdx)) {
                    continue;
                }
                lengths[ii] = maxNSpliceOtherSideN;
            }
            
            if ((maxN > lengths.length) || (maxN > 10000000)) {
                MultiArrayMergeSort.sortByDecr(lengths, pixIndexes);
            } else {
                CountingSort.sortByDecr(lengths, pixIndexes, maxN);
            }
            
            int pixIdx0 = pixIndexes[0];
            
            int pixIdx1 = pixIndexes[1];
            
            PairInt loc0 = junctionLocationMap.get(Integer.valueOf(pixIdx0));
            final int edge0Idx = loc0.getX();
            final int indexWithinEdge0 = loc0.getY();
            
            PairInt loc1 = junctionLocationMap.get(Integer.valueOf(pixIdx1));
            final int edge1Idx = loc1.getX();
            final int indexWithinEdge1 = loc1.getY();
                      
            if (edge0Idx != edge1Idx && (lengths[0] != 0) && (lengths[1] != 0)) {
                
                Set<Integer> edgePixelIndexes = 
                    theEdgeToPixelIndexMap.get(Integer.valueOf(edge0Idx));
                assert(edgePixelIndexes != null);
                
                int[] splice0Y = new int[]{indexWithinEdge0};
                Splice splice0 = new Splice(edges.get(edge0Idx));

                Map<Integer, PairInt> smallerSpliceLocations0 = new
                    HashMap<Integer, PairInt>();
                Map<Integer, PairInt> largerSpliceLocations0 = new
                    HashMap<Integer, PairInt>();
                // splice splice0 into 2 edges and put updated i
                PairIntArray[] spliced0 = splice0.splice(splice0Y, 
                    edgePixelIndexes, junctionLocationMap, 
                    smallerSpliceLocations0, largerSpliceLocations0);
                
                int[] splice1Y = new int[]{indexWithinEdge1};
                Splice splice1 = new Splice(edges.get(edge1Idx));  
                
                edgePixelIndexes = 
                    theEdgeToPixelIndexMap.get(Integer.valueOf(edge1Idx));
                assert(edgePixelIndexes != null);
                Map<Integer, PairInt> smallerSpliceLocations1 = new
                    HashMap<Integer, PairInt>();
                Map<Integer, PairInt> largerSpliceLocations1 = new
                    HashMap<Integer, PairInt>();
                PairIntArray[] spliced1 = splice1.splice(splice1Y,
                    edgePixelIndexes, junctionLocationMap, 
                    smallerSpliceLocations1, largerSpliceLocations1);
                
                // if not spliced, continue.
                if ((spliced0[0].getN() == 0) || (spliced1[0].getN() == 0)) {
                    continue;
                }
                
                int splice01InsertedIdx = edges.size();
                
                // add the smaller part of spliced0 and spliced1 to edges
                edges.add(spliced0[1]);
                
                int splice11InsertedIdx = edges.size();
                
                edges.add(spliced1[1]);
                
                int splice00N = spliced0[0].getN();
                int splice10N = spliced1[0].getN();
                
                // if splice0Y[0] is first point, reverse the edge before append
                if (splice0Y[0] == 0) {
                    spliced0[0].reverse();
                    splice0Y[0] = splice00N - 1;
                    
                    // --- update largerSpliceLocations0 to reverse index within edge
                    for (Entry<Integer, PairInt> sEntry : largerSpliceLocations0.entrySet()) {
                        PairInt sLoc = sEntry.getValue();
                        int idxRev = splice00N - sLoc.getY() - 1;
                        
                        sLoc.setY(idxRev);
                        assert(idxRev < splice00N);                        
                    }
                }
                
                // if splice1Y[0] is not the first point, reverse the edge before append
                if (splice1Y[0] != 0) {
                    spliced1[0].reverse();
                    splice1Y[0] = 0;
                    
                    // --- update largerSpliceLocations1 to reverse index within edge
                    for (Entry<Integer, PairInt> sEntry : largerSpliceLocations1.entrySet()) {
                        PairInt sLoc = sEntry.getValue();
                        int idxRev = splice10N - sLoc.getY() - 1;
                        
                        sLoc.setY(idxRev);
                        assert(idxRev < splice10N);
                    }
                }
                
                // append splice1 to splice0                
                spliced0[0].addAll(spliced1[0]);
                
                ++nSplices;
                
                PairIntArray edge0 = edges.get(edge0Idx);
                edge0.swapContents(spliced0[0]);
                spliced0[0] = null;
                
                // --- update location map for information in
                // --- smallerSpliceLocations0, largerSpliceLocations0 and
                // --- smallerSpliceLocations1, largerSpliceLocations1
                for (Entry<Integer, PairInt> sEntry : largerSpliceLocations0.entrySet()) {
                    Integer sPixelIndex = sEntry.getKey();
                    PairInt sLoc = sEntry.getValue();
                    PairInt loc = junctionLocationMap.get(sPixelIndex);
                    assert (loc != null);
                    
                    // edge index remains same, but the index within edge may have changed
                    loc.setY(sLoc.getY());
                    assert(sLoc.getY() < edge0.getN());
                }
                for (Entry<Integer, PairInt> sEntry : smallerSpliceLocations0.entrySet()) {
                    Integer sPixelIndex = sEntry.getKey();
                    PairInt sLoc = sEntry.getValue();
                    PairInt loc = junctionLocationMap.get(sPixelIndex);
                    assert (loc != null);
                    
                    // location is new edge with edited index within edge
                    loc.setX(splice01InsertedIdx);
                    loc.setY(sLoc.getY());
                    
                    //-- move item from theEdgeToPixelIndexMap too
                    //-- remove sPixelIndex from edge loc0.getX() and add it to splice01InsertedIdx
                    Set<Integer> sPixelIndexes = theEdgeToPixelIndexMap.get(Integer.valueOf(edge0Idx));
                    boolean removed = sPixelIndexes.remove(sPixelIndex);
                    assert(removed == true);
                    
                    sPixelIndexes = theEdgeToPixelIndexMap.get(Integer.valueOf(splice01InsertedIdx));
                    if (sPixelIndexes == null) {
                        sPixelIndexes = new HashSet<Integer>();
                    }
                    sPixelIndexes.add(sPixelIndex);
                    theEdgeToPixelIndexMap.put(Integer.valueOf(splice01InsertedIdx), sPixelIndexes);
                }
                for (Entry<Integer, PairInt> sEntry : smallerSpliceLocations1.entrySet()) {
                    Integer sPixelIndex = sEntry.getKey();
                    PairInt sLoc = sEntry.getValue();
                    PairInt loc = junctionLocationMap.get(sPixelIndex);
                    assert (loc != null);
                    
                    // location is new edge with edited index within edge
                    loc.setX(splice11InsertedIdx);
                    loc.setY(sLoc.getY());
                    
                    //-- move item from theEdgeToPixelIndexMap too
                    //-- remove sPixelIndex from edge loc1.getX() and add it to splice11InsertedIdx
                    Set<Integer> sPixelIndexes = 
                        theEdgeToPixelIndexMap.get(Integer.valueOf(edge1Idx));
                    boolean removed = sPixelIndexes.remove(sPixelIndex);
                    assert(removed == true);
                    theEdgeToPixelIndexMap.put(Integer.valueOf(edge1Idx), sPixelIndexes);
                    
                    sPixelIndexes = theEdgeToPixelIndexMap.get(Integer.valueOf(splice11InsertedIdx));
                    if (sPixelIndexes == null) {
                        sPixelIndexes = new HashSet<Integer>();
                    }
                    sPixelIndexes.add(sPixelIndex);
                    theEdgeToPixelIndexMap.put(Integer.valueOf(splice11InsertedIdx), sPixelIndexes);
                }
                for (Entry<Integer, PairInt> sEntry : largerSpliceLocations1.entrySet()) {
                    Integer sPixelIndex = sEntry.getKey();
                    PairInt sLoc = sEntry.getValue();
                    PairInt loc = junctionLocationMap.get(sPixelIndex);
                    assert (loc != null);
                   
                    // location after append is loc0.getX() with offset of splice00N
                    loc.setX(edge0Idx);
                    loc.setY(sLoc.getY() + splice00N);
                    assert(loc.getY() < edge0.getN());
                    
                    //-- move item from theEdgeToPixelIndexMap too
                    //-- remove sPixelIndex from edge loc1.getX() and add it to loc0.getX()
                    Integer key = Integer.valueOf(edge1Idx);
                    Set<Integer> sPixelIndexes = theEdgeToPixelIndexMap.get(key);
                    boolean removed = sPixelIndexes.remove(sPixelIndex);
                    assert(removed == true);
                    theEdgeToPixelIndexMap.put(key, sPixelIndexes);
                    
                    key = Integer.valueOf(edge0Idx);
                    sPixelIndexes = theEdgeToPixelIndexMap.get(key);
                    if (sPixelIndexes == null) {
                        sPixelIndexes = new HashSet<Integer>();
                    }
                    sPixelIndexes.add(sPixelIndex);
                    theEdgeToPixelIndexMap.put(key, sPixelIndexes);
                }
                
                // remove edge1Idx from edges as it's now part of edge edge0Idx
                PairIntArray removedEdge = edges.remove(edge1Idx);
                assert(removedEdge != null);

                if (debug) {
                    // assert that theEdgeToPixelIndexMap does not have any
                    // entries for the edge that was just removed.
                    // all entries in theEdgeToPixelIndexMap and the junction
                    // maps should have been updated above for edge1Idx
                    Set<Integer> pixelIndexes = theEdgeToPixelIndexMap.get(edge1Idx);
                    assert(pixelIndexes == null || pixelIndexes.isEmpty());
                }

                // ---- update all entries in junctionLocationMap and pixelIndexes for the removal
                // ---- there shouldn't be any entries in for loc1.getX() at this point though
                /*
                0                0
                ...             ...
                49               49
                50 <-- rm        prev 51, now 50.
                51 
                */
                for (int eIdx = (edge1Idx + 1); eIdx <= edges.size(); ++eIdx) {
                                        
                    Integer edgeIndex = Integer.valueOf(eIdx);
                    
                    Set<Integer> pixelIndexes = theEdgeToPixelIndexMap.get(edgeIndex);
                   
                    if (pixelIndexes != null) {
                        for (Integer pixelIndex : pixelIndexes) {
                            PairInt loc = junctionLocationMap.get(pixelIndex);
                            loc.setX(loc.getX() - 1);
                        }
                        
                        theEdgeToPixelIndexMap.remove(edgeIndex);
                        
                        Integer edgeIndexEarlier = Integer.valueOf(eIdx - 1);
                    
                        theEdgeToPixelIndexMap.put(edgeIndexEarlier, pixelIndexes);
                    }                                        
                }
                
                if (debug) {
                    algorithms.misc.MiscDebug.assertConsistentEdgeCapacity(
                        theEdgeToPixelIndexMap, junctionLocationMap, 
                        junctionMap, edges);
                }
            }
        } 
        
        return nSplices;
    }

    protected void reduceMultipleEndpointsForEdge(
        List<PairIntArray> edges,
        Map<Integer, Set<PointPairInt>> edgeFirstEndPointMap,
        Map<Integer, Set<PointPairInt>> edgeLastEndPointMap, 
        Map<PairInt, PairInt> endPointMap, Set<PointPairInt> theJoinPoints) {
        
        Map<Integer, Set<PointPairInt>> tmpFirstRemoveMap = new HashMap<Integer, Set<PointPairInt>>();
        
        Map<Integer, Set<PointPairInt>> tmpLastRemoveMap = new HashMap<Integer, Set<PointPairInt>>();
        
        for (int type = 0; type < 2; ++type) {
        
            Map<Integer, Set<PointPairInt>> edgeMap = edgeFirstEndPointMap;
            if (type == 1) {
                edgeMap = edgeLastEndPointMap;
            }
            
            for (Entry<Integer, Set<PointPairInt>> entry : edgeMap.entrySet()) {

                Set<PointPairInt> joinPointSet = entry.getValue();

                if (joinPointSet.size() < 2) {
                    continue;
                }

                // has more than one joint point for this endpoint, though some
                // might have already been removed and are in tmpFirstRemoveMap
                // if so
     
                if (log.isLoggable(Level.FINE)) {
                    // construct a warning
                    StringBuilder sb = new StringBuilder("edge ");
                    sb.append(entry.getKey().toString())
                        .append(" has ").append(Integer.toString(joinPointSet.size()))
                        .append(" first endpoints, so deciding between them:");
                
                    for (PointPairInt joinPoint : entry.getValue()) {
                        PairInt loc0 = joinPoint.getKey();
                        PairInt loc1 = joinPoint.getValue();            
                        PairIntArray edge0 = edges.get(loc0.getX());
                        PairIntArray edge1 = edges.get(loc1.getX());            
                        int x0 = edge0.getX(loc0.getY());
                        int y0 = edge0.getY(loc0.getY());
                        int x1 = edge1.getX(loc1.getY());
                        int y1 = edge1.getY(loc1.getY());                
                        sb.append("\n")
                            .append(" ").append(Integer.toString(loc0.getX())).append(":")
                            .append(Integer.toString(loc0.getY()))
                            .append(" (").append(Integer.toString(x0)).append(",")
                            .append(Integer.toString(y0)).append(")<-->")
                            .append(" ").append(Integer.toString(loc1.getX())).append(":")
                            .append(Integer.toString(loc1.getY()))
                            .append(" (").append(Integer.toString(x1)).append(",")
                            .append(Integer.toString(y1)).append(")");
                    }
                    log.warning(sb.toString());
                }
                
                // choose closest join point pair, and break ties with those that
                // have higher number of members that are already endpoints.

                PointPairInt closestJoinPoint = null;
                int closestDistSq = Integer.MAX_VALUE;
                int closestCanBeReordered0 = -99;
                int closestCanBeReordered1 = -99;
                
                for (PointPairInt joinPoint : joinPointSet) {

                    PairInt loc0 = joinPoint.getKey();
                    PairIntArray edge0 = edges.get(loc0.getX());
                    int n0 = edge0.getN();
                    
                    Set<PointPairInt> removedFirstSet = 
                        tmpFirstRemoveMap.get(Integer.valueOf(loc0.getX()));
                    if (removedFirstSet == null) {
                        removedFirstSet = new HashSet<PointPairInt>();
                        tmpFirstRemoveMap.put(Integer.valueOf(loc0.getX()), removedFirstSet);
                    }
                    
                    Set<PointPairInt> removedLastSet = 
                        tmpLastRemoveMap.get(Integer.valueOf(loc0.getX()));
                    if (removedLastSet == null) {
                        removedLastSet = new HashSet<PointPairInt>();
                        tmpLastRemoveMap.put(Integer.valueOf(loc0.getX()), removedLastSet);
                    }

                    if ((n0 != 3) && (loc0.getY() < 2) && removedFirstSet.contains(joinPoint)) {
                        continue;
                    }
                    if ((n0 != 3) && (loc0.getY() > (n0 - 3)) && removedLastSet.contains(joinPoint)) {
                        continue;
                    }

                    PairInt loc1 = joinPoint.getValue();
                    PairIntArray edge1 = edges.get(loc1.getX());

                    assert(edge0 != null);
                    assert(edge1 != null);

                    int diffX = edge0.getX(loc0.getY()) - edge1.getX(loc1.getY());
                    int diffY = edge0.getY(loc0.getY()) - edge1.getY(loc1.getY());
                    int distSq = (diffX * diffX) + (diffY * diffY);
                    if (distSq < closestDistSq) {
                        closestDistSq = distSq;
                        closestJoinPoint = joinPoint;
                        closestCanBeReordered0 = canBeReordered(loc0, edge0);
                        closestCanBeReordered1 = canBeReordered(loc1, edge1);
                    } else if (distSq == closestDistSq) {
                        int cr0 = canBeReordered(loc0, edge0);
                        int cr1 = canBeReordered(loc1, edge1);;
                        if (closestCanBeReordered0 == 0 && closestCanBeReordered1 == 0) {
                            continue;
                        }
                        if (((cr0 == 0) && (cr1 == 0)) ||
                            (closestCanBeReordered0 != 0 && closestCanBeReordered1 != 0
                                && ((cr0 == 0) || cr1 == 0))
                            ) {
                            closestDistSq = distSq;
                            closestJoinPoint = joinPoint;
                            closestCanBeReordered0 = cr0;
                            closestCanBeReordered1 = cr1;
                        }
                    }
                }

                //TODO: revisit this 
                if (closestJoinPoint == null) {
                    continue;
                }                
                //assert (closestJoinPoint != null);

                for (PointPairInt joinPoint : joinPointSet) {

                    if (joinPoint.equals(closestJoinPoint)) {
                        continue;
                    }

                    PairInt loc0 = joinPoint.getKey();
                    PairIntArray edge0 = edges.get(loc0.getX());
                    int n0 = edge0.getN();

                    Set<PointPairInt> removedFirstSet0 = 
                        tmpFirstRemoveMap.get(Integer.valueOf(loc0.getX()));
                   
                    Set<PointPairInt> removedLastSet0 = 
                        tmpLastRemoveMap.get(Integer.valueOf(loc0.getX()));
                    
                    if ((n0 != 3) && (loc0.getY() < 2) && removedFirstSet0.contains(joinPoint)) {
                        continue;
                    }
                    if ((n0 != 3) && (loc0.getY() > (n0 - 3)) && removedLastSet0.contains(joinPoint)) {
                        continue;
                    }

                    PairInt loc1 = joinPoint.getValue();
                    PairIntArray edge1 = edges.get(loc1.getX());
                    int n1 = edge1.getN();

                    assert(edge0 != null);
                    assert(edge1 != null);

                    if ((n0 == 3) && (loc0.getY() == 1)) {
                        log.fine(String.format(
                            "removing point (%d,%d) but it is in the middle of an edge of size 3",
                            edge0.getX(1), edge0.getY(1)));
                    }
                    
                    Set<PointPairInt> removedFirstSet1 = 
                        tmpFirstRemoveMap.get(Integer.valueOf(loc1.getX()));
                    if (removedFirstSet1 == null) {
                        removedFirstSet1 = new HashSet<PointPairInt>();
                        tmpFirstRemoveMap.put(Integer.valueOf(loc1.getX()), removedFirstSet1);
                    }
                    
                    Set<PointPairInt> removedLastSet1 = 
                        tmpLastRemoveMap.get(Integer.valueOf(loc1.getX()));
                    if (removedLastSet1 == null) {
                        removedLastSet1 = new HashSet<PointPairInt>();
                        tmpLastRemoveMap.put(Integer.valueOf(loc1.getX()), removedLastSet1);
                    }

                    boolean isFirstEndPoint0 = false;
                    if (loc0.getY() == 0) {
                        isFirstEndPoint0 = true;
                    } else if ((loc0.getY() != (n0 - 1)) && (loc0.getY() < 2)) {
                        isFirstEndPoint0 = true;
                    }
                    boolean isFirstEndPoint1 = false;
                    if (loc1.getY() == 0) {
                        isFirstEndPoint1 = true;
                    } else if ((loc1.getY() != (n1 - 1)) && (loc1.getY() < 2)) {
                        isFirstEndPoint1 = true;
                    }

                    if (isFirstEndPoint0 && removedFirstSet0.contains(joinPoint)) {
                        continue;
                    } else if (!isFirstEndPoint0 && removedLastSet0.contains(joinPoint)) {
                        continue;
                    }

                    if (isFirstEndPoint0) {
                        removedFirstSet0.add(joinPoint);
                    } else {
                        removedLastSet0.add(joinPoint);
                    }
                    if (isFirstEndPoint1) {
                        removedFirstSet1.add(joinPoint);
                    } else {
                        removedLastSet1.add(joinPoint);
                    }

                    boolean removed = theJoinPoints.remove(joinPoint);
                    assert(removed);
                    
                    tmpFirstRemoveMap.put(Integer.valueOf(loc0.getX()), removedFirstSet0);
                    tmpLastRemoveMap.put(Integer.valueOf(loc0.getX()), removedLastSet0);
                    
                    tmpFirstRemoveMap.put(Integer.valueOf(loc1.getX()), removedFirstSet1);
                    tmpLastRemoveMap.put(Integer.valueOf(loc1.getX()), removedLastSet1);
                }                
            }
            
            // --- update the edge maps for the changes ---
            
            Map<Integer, Set<PointPairInt>> tmpMap = new HashMap<Integer, Set<PointPairInt>>();
            
            for (Entry<Integer, Set<PointPairInt>> entry : edgeFirstEndPointMap.entrySet()) {
                
                Set<PointPairInt> rmJoinPointSet = tmpFirstRemoveMap.get(entry.getKey());
                
                if (rmJoinPointSet == null || rmJoinPointSet.isEmpty()) {
                    tmpMap.put(entry.getKey(), entry.getValue());
                    continue;
                }
                
                Set<PointPairInt> tmpJoinPointSet = new HashSet<PointPairInt>();
                for (PointPairInt joinPoint : entry.getValue()) {
                    if (!rmJoinPointSet.contains(joinPoint)) {
                        tmpJoinPointSet.add(joinPoint);
                    }
                }
                tmpMap.put(entry.getKey(), tmpJoinPointSet);
            }
            edgeFirstEndPointMap.clear();
            edgeFirstEndPointMap.putAll(tmpMap);
            
            tmpMap = new HashMap<Integer, Set<PointPairInt>>();
            
            for (Entry<Integer, Set<PointPairInt>> entry : edgeLastEndPointMap.entrySet()) {
                
                Set<PointPairInt> rmJoinPointSet = tmpLastRemoveMap.get(entry.getKey());
                
                if (rmJoinPointSet == null || rmJoinPointSet.isEmpty()) {
                    tmpMap.put(entry.getKey(), entry.getValue());
                    continue;
                }
                
                Set<PointPairInt> tmpJoinPointSet = new HashSet<PointPairInt>();
                for (PointPairInt joinPoint : entry.getValue()) {
                    if (!rmJoinPointSet.contains(joinPoint)) {
                        tmpJoinPointSet.add(joinPoint);
                    }
                }
                tmpMap.put(entry.getKey(), tmpJoinPointSet);
            }
            edgeLastEndPointMap.clear();
            edgeLastEndPointMap.putAll(tmpMap);
        }
        
    }

    public PairIntArray findAsSingleClosedEdge() {
        
        this.singleClosedEdge = true;
        
        List<PairIntArray> output = connectPixelsViaDFSForBounds();
        
Image img0 = ImageIOHelper.convertImage(img);
for (int i = 0; i < output.size(); ++i) {
    PairIntArray pa = output.get(i);
    for (int j = 0; j < pa.getN(); ++j) {
        int x = pa.getX(j);
        int y = pa.getY(j);
        if (i == 0) {
            if (j == 0 || (j == (pa.getN() - 1))) {
                ImageIOHelper.addPointToImage(x, y, img0, 0, 200, 100, 0);
            } else {
                ImageIOHelper.addPointToImage(x, y, img0, 0, 255, 0, 0);
            }
        } else if (i == 1) {
            ImageIOHelper.addPointToImage(x, y, img0, 0, 0, 255, 0);
        } else {
            ImageIOHelper.addPointToImage(x, y, img0, 0, 0, 0, 255);
        }
    }
}
MiscDebug.writeImageCopy(img0, "output_before_merges_" + MiscDebug.getCurrentTimeFormatted() + ".png");       
        
        output = findEdgesIntermediateSteps(output);
                
Image img2 = ImageIOHelper.convertImage(img);
for (int i = 0; i < output.size(); ++i) {
    PairIntArray pa = output.get(i);
    for (int j = 0; j < pa.getN(); ++j) {
        int x = pa.getX(j);
        int y = pa.getY(j);
        if (i == 0) {
            if (j == 0 || (j == (pa.getN() - 1))) {
                ImageIOHelper.addPointToImage(x, y, img2, 0, 200, 100, 0);
            } else {
                ImageIOHelper.addPointToImage(x, y, img2, 0, 255, 0, 0);
            }
        } else if (i == 1) {
            ImageIOHelper.addPointToImage(x, y, img2, 0, 0, 255, 0);
        } else {
            ImageIOHelper.addPointToImage(x, y, img2, 0, 0, 0, 255);
        }
    }
}
MiscDebug.writeImageCopy(img2, "output_after_merges_" + MiscDebug.getCurrentTimeFormatted() + ".png");

        if (output.size() > 1) {
            Collections.sort(output, new PairIntArrayDescendingComparator());            
        }
        
        PairIntArray out = output.get(0);
        
        reorderEndpointsIfNeeded(out);
                
        if (endPointsAreAdjacent(out)) {
            PairIntArrayWithColor p = new PairIntArrayWithColor(out);
            p.setColor(1);
            out = p;
        }
                    
        return out;
    }
    
    private Map<Integer, Set<Integer>> createEdgeToPixelIndexMap() {
        // key = edge index.  
        // value = pixel indexes.
        //   the pixel indexes are used to find values in junctionLocatorMap
        //   to update it as points are moved to and from edges.
        Map<Integer, Set<Integer>> theEdgeToPixelIndexMap = new HashMap<Integer, Set<Integer>>();
        
        for (Entry<Integer, PairInt> entry : junctionLocationMap.entrySet()) {
            Integer pixelIndex = entry.getKey();
            PairInt loc = entry.getValue();
            Integer edgeIndex = Integer.valueOf(loc.getX());
            Set<Integer> pixelIndexes = theEdgeToPixelIndexMap.get(edgeIndex);
            if (pixelIndexes == null) {
                pixelIndexes = new HashSet<Integer>();
            }
            pixelIndexes.add(pixelIndex);
            theEdgeToPixelIndexMap.put(edgeIndex, pixelIndexes);
        }
        
        return theEdgeToPixelIndexMap;
    }

    private int reorderPointsInEdgesForClosedCurve(List<PairIntArray> output) {
        
        if (output == null || output.isEmpty()) {
            return 0;
        }
        
        int maxEdgeIdx = 0;
        int nMaxEdge = output.get(0).getN();
        for (int i = 1; i < output.size(); ++i) {
            int n = output.get(i).getN();
            if (n > nMaxEdge) {
                nMaxEdge = n;
                maxEdgeIdx = i;
            }
        }
        
        int nChanged = 0;
        
        // key = edge index.  
        // value = pixel indexes.
        //   the pixel indexes are used to find values in junctionLocatorMap
        //   to update it as points are moved to and from edges.
        Map<Integer, Set<Integer>> theEdgeToPixelIndexMap = createEdgeToPixelIndexMap();
        
        /**
        * map with key = center of junction pixel coordinates; 
        * value = set of adjacent pixels when there are more than the preceding 
        * and next.
        Map<Integer, Set<Integer>> junctionMap
        
        * map with key = pixel coordinates of all pixels involved in junctions;
        * value = PairInt holding index of edge that pixel is located in and
        * holding the index within that edge of the pixel.
        * for example, a pixel located in edges(0) at offset=100
        * would have PairInt(0, 100) as a value.
        Map<Integer, PairInt> junctionLocationMap
        */
              
Image img2 = ImageIOHelper.convertImage(img);
for (int i = 0; i < output.size(); ++i) {
    PairIntArray pa = output.get(i);
    for (int j = 0; j < pa.getN(); ++j) {
        int x = pa.getX(j);
        int y = pa.getY(j);
        if (i == 0) {
            if (j == 0 || (j == (pa.getN() - 1))) {
                ImageIOHelper.addPointToImage(x, y, img2, 0, 200, 100, 0);
            } else {
                ImageIOHelper.addPointToImage(x, y, img2, 0, 255, 0, 0);
            }
        } else if (i == 1) {
            ImageIOHelper.addPointToImage(x, y, img2, 0, 0, 255, 0);
        } else {
            ImageIOHelper.addPointToImage(x, y, img2, 0, 0, 0, 255);
        }
    }
}
MiscDebug.writeImageCopy(img2, "output_" + MiscDebug.getCurrentTimeFormatted() + ".png");
int z0 = 1;

        PairIntArray maxEdge = output.get(maxEdgeIdx);
        
        for (int i = 0; i < output.size(); ++i) {
            
            if (i == maxEdgeIdx || (output.get(i).getN() < 3)) {
                continue;
            }
            
            Map<PairInt, Set<PairInt>> edgeLocAdjToMax = new HashMap<PairInt, Set<PairInt>>();
            
            populateWithAdjacentLocations(theEdgeToPixelIndexMap, maxEdgeIdx, i, 
                edgeLocAdjToMax);
           
            // need at least 2 adjacencies
            if (edgeLocAdjToMax.size() < 2) {
                continue;
            }
            
            Set<PairInt> maxAdjToFirstPointLoc = null;
            Set<PairInt> maxAdjToLastPointLoc = null;
            
            for (Entry<PairInt, Set<PairInt>> entry : edgeLocAdjToMax.entrySet()) {
                PairInt edgeLoc = entry.getKey();
                if (edgeLoc.getY() == 0) {
                    maxAdjToFirstPointLoc = entry.getValue();
                } else if (edgeLoc.getY() == (output.get(i).getN() - 1)) {
                    maxAdjToLastPointLoc = entry.getValue();
                }
            }
            
            /*
            sometimes, this edge is prevented from merging because of the order
            of its points.
            
            For example, if the larger edge is '@''s and this edge is the 
            numbers in current order, can see that the splice algorithm wont get 
            the "potential added by merging line" correct because the line is not
            well formed yet for a merge.
            That splice algorithm will see that it could add 3 points and would
            then dead end for merging this edge.  If this edge were re-ordered, 
            the splice algorithm would see it could add 6 points.
                               @ @ @ @
                             @
                           @
                    0 1 2 @
                    5 4 3 @
                            @ @
            
            (The splice algorithm run twice should then merge on the remaining
            part of the originally longer main edge.
            
            Note that if this re-ordered line is shorter than the segment on the
            main line already part of the main line, the re-ordered line will 
            not be chosen, so an insertion algorithm is needed after all splices).            
            */
            
            // if this edge's first and last points are in the adj list,
            // this is already well ordered and should instead
            // be considered for insertion.
            if ((maxAdjToFirstPointLoc != null) && (maxAdjToLastPointLoc != null)) {
                // this can be handled in the insert algorithm
                continue;
            }
            
            // find the 2 closest pairs in edgeLocAdjToMax and those will be
            // the re-ordered first and last points.
            // need to use optimal pairing
            
            // since these are small lists, will not use a set even though search is up to O(m)
            List<PairInt> currentEdgePoints = new ArrayList<PairInt>();
            List<PairInt> maxEdgePoints = new ArrayList<PairInt>();
           
            for (Entry<PairInt, Set<PairInt>> entry : edgeLocAdjToMax.entrySet()) {
                PairInt p1 = entry.getKey();
                if (!currentEdgePoints.contains(p1)) {
                    currentEdgePoints.add(p1);
                }
                for (PairInt p2 : entry.getValue()) {
                    if (!maxEdgePoints.contains(p2)) {
                        maxEdgePoints.add(p2);
                    }
                }
            }
            
            PairIntArray currentEdge = output.get(i);
            
            float[][] distances = new float[currentEdgePoints.size()][];
            float[][] distancesCopy = new float[currentEdgePoints.size()][];
            
            for (int ii = 0; ii < currentEdgePoints.size(); ++ii) {
                PairInt pLoc = currentEdgePoints.get(ii);
                int x = currentEdge.getX(pLoc.getY());
                int y = currentEdge.getY(pLoc.getY());
                distances[ii] = new float[maxEdgePoints.size()];
                distancesCopy[ii] = new float[maxEdgePoints.size()];
                for (int jj = 0; jj < maxEdgePoints.size(); ++jj) {
                    PairInt p2Loc = maxEdgePoints.get(jj);
                    int x2 = maxEdge.getX(p2Loc.getY());
                    int y2 = maxEdge.getY(p2Loc.getY());
                    int diffX = Math.abs(x - x2);
                    int diffY = Math.abs(y - y2);
                    // if not adjacent, give a very large number
                    if ((diffX > 1) || (diffY > 1)) {
                        distances[ii][jj] = Float.MAX_VALUE;
                        distancesCopy[ii][jj] = Float.MAX_VALUE;
                    } else {
                        float dist = (float)Math.sqrt((diffX*diffX) + (diffY*diffY));
                        distances[ii][jj] = dist;
                        distancesCopy[ii][jj] = dist;
                    }
                }
            }
            boolean transposed = false;
            if (currentEdgePoints.size() > maxEdgePoints.size()) {
                distancesCopy = MatrixUtil.transpose(distancesCopy);
                transposed = true;
            }
            
            HungarianAlgorithm hAlg = new HungarianAlgorithm();
            int[][] match = hAlg.computeAssignments(distances);
            
            // idx1, idx2  dist  sort by smallest dist
            PairInt[] idxs = new PairInt[match.length];
            float[] d = new float[idxs.length];
                        
            int count = 0;
            for (int ii = 0; ii < match.length; ii++) {
                int idx1 = match[ii][0];
                int idx2 = match[ii][1];
                if (idx1 == -1 || idx2 == -1) {
                    continue;
                }
                if (idx1 == Float.MAX_VALUE || idx2 == Float.MAX_VALUE) {
                    continue;
                }
                if (transposed) {
                    int swap = idx1;
                    idx1 = idx2;
                    idx2 = swap;
                }
                float dist = distancesCopy[idx1][idx2];
                if (dist == Float.MAX_VALUE) {
                    continue;
                }
                idxs[count] = new PairInt(idx1, idx2);
                d[count] = dist;
                
                count++;
            }
            
            if (count < 2) {
                continue;
            }
            
            idxs = Arrays.copyOf(idxs, count);
            d = Arrays.copyOf(d, count);
            
            QuickSort.sortBy1stArg(d, idxs);

            PairInt cFirstLoc = currentEdgePoints.get(idxs[0].getX());
            PairInt mFirstLoc = maxEdgePoints.get(idxs[0].getY());
            
            PairInt cLastLoc = currentEdgePoints.get(idxs[1].getX());
            PairInt mLastLoc = maxEdgePoints.get(idxs[1].getY());
            
            // reorder currentEdge to be between cFirstLoc and cLastLoc
            PairIntArray reordered = reorder(currentEdge, cFirstLoc.getY(), 
                cLastLoc.getY());
          
            if (reordered == null) {
                continue;   
            }   

            //insert
            if (mFirstLoc.getY() < mLastLoc.getY()) {
                       
                maxEdge.insertAll(mFirstLoc.getY() + 1, currentEdge);
                        
            } else {

                currentEdge.reverse();

                maxEdge.insertAll(mLastLoc.getY() + 1, currentEdge);
            }
            
            nChanged++;
            
            return nChanged;
        }
        
        return nChanged;
    }
    
    private int insertAdjacentForClosedCurve(List<PairIntArray> output) {
        
        if (output == null || output.isEmpty()) {
            return 0;
        }
        
        int maxEdgeIdx = 0;
        int nMaxEdge = output.get(0).getN();
        for (int i = 1; i < output.size(); ++i) {
            int n = output.get(i).getN();
            if (n > nMaxEdge) {
                nMaxEdge = n;
                maxEdgeIdx = i;
            }
        }
        PairIntArray maxEdge = output.get(maxEdgeIdx);
        
        // key = edge index.  
        // value = pixel indexes.
        //   the pixel indexes are used to find values in junctionLocatorMap
        //   to update it as points are moved to and from edges.
        Map<Integer, Set<Integer>> theEdgeToPixelIndexMap = createEdgeToPixelIndexMap();
        
        /**
        * map with key = center of junction pixel coordinates; 
        * value = set of adjacent pixels when there are more than the preceding 
        * and next.
        Map<Integer, Set<Integer>> junctionMap
        
        * map with key = pixel coordinates of all pixels involved in junctions;
        * value = PairInt holding index of edge that pixel is located in and
        * holding the index within that edge of the pixel.
        * for example, a pixel located in edges(0) at offset=100
        * would have PairInt(0, 100) as a value.
        Map<Integer, PairInt> junctionLocationMap
        */
             
        int nChanged = 0;
        
        for (int i = 0; i < output.size(); ++i) {
            
            if (i == maxEdgeIdx || (output.get(i).getN() < 1)) {
                continue;
            }
            
            Map<PairInt, Set<PairInt>> edgeLocAdjToMax = new HashMap<PairInt, Set<PairInt>>();
            
            populateWithAdjacentLocations(theEdgeToPixelIndexMap, maxEdgeIdx, i, 
                edgeLocAdjToMax);
            
            Set<PairInt> maxAdjToFirstPointLoc = null;
            Set<PairInt> maxAdjToLastPointLoc = null;

            for (Entry<PairInt, Set<PairInt>> entry : edgeLocAdjToMax.entrySet()) {
                PairInt edgeLoc = entry.getKey();
                if (edgeLoc.getY() == 0) {
                    maxAdjToFirstPointLoc = entry.getValue();
                } else if (edgeLoc.getY() == (output.get(i).getN() - 1)) {
                    maxAdjToLastPointLoc = entry.getValue();
                }
            }
            
            if ((maxAdjToFirstPointLoc != null) && (maxAdjToLastPointLoc == null)) {
                // see if points near the end can be re-ordered
                int z = 1;
            } else if ((maxAdjToFirstPointLoc == null) && (maxAdjToLastPointLoc != null)) {
                // see if points near the beginning can be re-ordered
                PairInt closestToZero = null;
                Set<PairInt> maxAdjToFirstPointLocTmp = null;
                int minIdx = Integer.MAX_VALUE;
                for (Entry<PairInt, Set<PairInt>> entry : edgeLocAdjToMax.entrySet()) {
                    PairInt edgeLoc = entry.getKey();
                    if (edgeLoc.getY() < minIdx) {
                        minIdx = edgeLoc.getY();
                        closestToZero = edgeLoc;
                        maxAdjToFirstPointLocTmp = entry.getValue();
                    }
                }
                if (closestToZero != null) {
                    int x1 = output.get(i).getX(0);
                    int y1 = output.get(i).getY(0);
                    /*
                                  (18,21) maxEdge 0
                                  (18,22) 2  ---closestToZero
                   (17,23) 3      (18,23) 1
                                  (18,24) 0
                    
                    if point loc at index closestToZero.getY() + 1 is adjacent
                    to the point at idx=0
                    can reverse the order of points from 0 to closestToZero.getY()
                    and set maxAdjToFirstPointLoc to (closestToZero.getX(), 0)
                    */
                    int lastSwapIdx = closestToZero.getY();
                    int xNext = output.get(i).getX(lastSwapIdx + 1);
                    int yNext = output.get(i).getY(lastSwapIdx + 1);
                    if ((Math.abs(xNext - x1) > 1) || (Math.abs(yNext - y1) > 1)) {
                        continue;
                    }
                    reverse0toIdx(output.get(i), lastSwapIdx);
                    maxAdjToFirstPointLoc = maxAdjToFirstPointLocTmp;
                }
                int z = 1;
            }
            
            /*
                               @ @ @ @
                             @
                           @
                    2 1 0 @
                    3 4 5 @
                            @ @                      
            */
           
            if ((maxAdjToFirstPointLoc != null) && (maxAdjToLastPointLoc != null)) {
                
                PairIntArray currentEdge = output.get(i);
               
                /*
                can insert another line for 2 different cases:
                (1) if the insertion point is 2 adjacent points on maxEdge
                (2) OR, if maxEdge is not a closed curve and the insertion
                points are it's endpoints. (Note that the possible insertion, 
                if not at endpoint of the open curve, may indicate that the
                open curve has a segment that needs to be reversed before this
                stage.  
                */
                
                // if can insert, need to break and return because the
                // junction maps need to be updated for this change
                
                /*
                find if there is a set of adjacent points, one from 
                maxAdjToFirstPointLoc and the other from maxAdjToLastPointLoc.
                The closest of those is a potential insertion location.
                
                Since these are small sets, will use brute force comparisons
                instead of closest pair, but should consider implementing
                closest pair for n >= 3.
                */
                //p1 and pn are w.r.t. maxEdge
                PairInt p1 = null;
                PairInt pn = null;
                double minDist = Integer.MAX_VALUE;
                for (PairInt firstAdj : maxAdjToFirstPointLoc) {
                    int xf = maxEdge.getX(firstAdj.getY());
                    int yf = maxEdge.getY(firstAdj.getY());
                    for (PairInt lastAdj : maxAdjToLastPointLoc) {
                        if (firstAdj.equals(lastAdj)) {
                            continue;
                        }
                        int xl = maxEdge.getX(lastAdj.getY());
                        int yl = maxEdge.getY(lastAdj.getY());
                        int diffX = Math.abs(xf - xl);
                        int diffY = Math.abs(yf - yl);
                        if (diffX > 1 || diffY > 1) {
                            continue;
                        }
                        double dist = Math.sqrt((diffX*diffX) + (diffY*diffY));
                        if (dist < minDist) {
                            minDist = dist;
                            p1 = firstAdj;
                            pn = lastAdj;
                        }
                        // TODO: could consider == case and the order of 
                        // points for a preferred best
                    }
                }
      
                if ((p1 != null) && (pn != null)) {
                    
                    if (p1.getY() < pn.getY()) {
                        
                        /*
                                   @ @ @ @
                                 @
                               @
                        2 1 0 @
                        3 4 5 @
                                @ @                      
                        */
                        
                        maxEdge.insertAll(p1.getY() + 1, currentEdge);
                        
                    } else {
                        
                        /*
                                   @ @ @ @
                                 @5
                               @4
                        2 1 0 @3 p1
                        3 4 5 @2 pn
                                @ @ 0
                        
                        currentEdge's 0  adjacent to p1
                        reversing currentEdge makes currentEdge's 0 adjacent to pn
                        so would insert after pn
                        */
                        
                        currentEdge.reverse();
                        
                        maxEdge.insertAll(pn.getY() + 1, currentEdge);
                    }
                        
                    nChanged++;
                    
                    //TODO: consider a visit pattern where do not need to
                    //  exit method for any change
                    
                    output.remove(i);
                    
                    return nChanged;
                } 
            }            
        }
        
        return nChanged;
    }
    
    /**
     * populate the output lists with adjacent pixel locations for 
     * the edge next to the reference edge (where reference edge is usually the
     * one with largest number of points).
     * 
     * @param theEdgeToPixelIndexMap
     * @param referenceEdgeIdx
     * @param edgeIdx
     * @param outputEdgeLocAdjToReference
     */
    private void populateWithAdjacentLocations(
        Map<Integer, Set<Integer>> theEdgeToPixelIndexMap, 
        int referenceEdgeIdx, int edgeIdx, 
        Map<PairInt, Set<PairInt>> outputEdgeLocAdjToReference) {
        
        Set<Integer> maxEdgePixelIndexes = theEdgeToPixelIndexMap.get(Integer.valueOf(referenceEdgeIdx));
        
        for (Integer pixIndex : maxEdgePixelIndexes) {
            PairInt loc = junctionLocationMap.get(pixIndex);
            // see if it is adjacent to this edge by looking for it in junctionMap
            Set<Integer> adjToMax = junctionMap.get(pixIndex);
            if (adjToMax == null) {
                continue;
            }
            for (Integer adjToMaxPixIndex : adjToMax) {
                PairInt adjToMaxLoc = junctionLocationMap.get(adjToMaxPixIndex);
                if (adjToMaxLoc.getX() == edgeIdx) {
                    // adjToMaxLoc is next to loc
                    Set<PairInt> adjMax = outputEdgeLocAdjToReference.get(adjToMaxLoc);
                    if (adjMax == null) {
                        adjMax = new HashSet<PairInt>();
                        outputEdgeLocAdjToReference.put(adjToMaxLoc, adjMax);
                    }
                    adjMax.add(loc);
                }
            }
        }
        for (Integer pixIndex : theEdgeToPixelIndexMap.get(Integer.valueOf(edgeIdx))) {
            PairInt loc = junctionLocationMap.get(pixIndex);
            // see if it is adjacent to this edge by looking for it in junctionMap
            Set<Integer> adjToPix = junctionMap.get(pixIndex);
            if (adjToPix == null) {
                // search the 8 neighbors to see if any are in maxEdgePixelIndexes
                int x = img.getCol(pixIndex);
                int y = img.getRow(pixIndex);
                adjToPix = new HashSet<Integer>();
                for (int i = 0; i < dxs8.length; ++i) {
                    int x2 = x + dxs8[i];
                    int y2 = y + dys8[i];
                    if (x2 < 0 || (y2 < 0) || (x2 > (img.getWidth() - 1) || (y2 > (img.getHeight() - 1)))) {
                        continue;
                    }
                    int pixIdx2 = img.getIndex(x2, y2);
                    if (maxEdgePixelIndexes.contains(Integer.valueOf(pixIdx2))) {
                        adjToPix.add(Integer.valueOf(pixIdx2));
                    }
                }
                if (adjToPix.isEmpty()) {
                    continue;
                }
            }
            for (Integer adjToPixIndex : adjToPix) {
                PairInt adjToLoc = junctionLocationMap.get(adjToPixIndex);
                if (adjToLoc.getX() == referenceEdgeIdx) {
                    Set<PairInt> adjMax = outputEdgeLocAdjToReference.get(loc);
                    if (adjMax == null) {
                        adjMax = new HashSet<PairInt>();
                        outputEdgeLocAdjToReference.put(loc, adjMax);
                    }
                    adjMax.add(adjToLoc);
                }
            }
        }
    }

    private PairIntArray reorder(final PairIntArray edge, int firstIdx, int lastIdx) {
        
        /*
        //note, single pixel wide spurs should have been removed already
        */
        
        /*
                             @ @ @ @
                             @
                           @
                    0 1 2 @
                    5 4 3 @
                            @ @
        */
        
        Set<PairInt> remainingPoints = Misc.convert(edge);
        
        Stack<PairInt> firstStack = new Stack<PairInt>();
        Stack<PairInt> lastStack = new Stack<PairInt>();
        
        PairIntArray ordered = new PairIntArray(edge.getN());
        
        // setting the values with a DFS traversal, but walked from first to middle and
        // last to middle, so initialize ordered with fake values
        for (int i = 0; i < edge.getN(); ++i) {
            ordered.add(-1, -1);
        }
        
        int idx1 = 0;
        int idxn = edge.getN();
        int midIdx = idxn >> 1;
        
        while ((idx1 == 0) || !remainingPoints.isEmpty()) {
                                            
            if (idx1 == 0) {
                
                ordered.set(idx1, edge.getX(idx1), edge.getY(idx1));
                remainingPoints.remove(new PairInt(edge.getX(idx1), edge.getY(idx1)));

                ordered.set(idxn, edge.getX(idxn), edge.getY(idxn));
                remainingPoints.remove(new PairInt(edge.getX(idxn), edge.getY(idxn)));

                boolean added = false;
                for (int nIdx = 0; nIdx < dxs8.length; nIdx++) {
                    int vX = edge.getX(idx1) + dxs8[nIdx];
                    int vY = edge.getY(idx1) + dys8[nIdx];
                    PairInt p = new PairInt(vX, vY);
                    if (remainingPoints.contains(p)) {
                        firstStack.add(p);
                        added = true;
                        break;
                    }
                }
                if (!added) {
                    return null;
                }
                
                added = false;
                for (int nIdx = 0; nIdx < dxs8.length; nIdx++) {
                    int vX = edge.getX(idxn) + dxs8[nIdx];
                    int vY = edge.getY(idxn) + dys8[nIdx];
                    PairInt p = new PairInt(vX, vY);
                    if (remainingPoints.contains(p)) {
                        lastStack.add(p);
                        added = true;
                        break;
                    }
                }
                if (!added) {
                    return null;
                }

                idx1++;
                idxn--;
                continue;
            }
            
            if (!firstStack.isEmpty()) {
                
                PairInt xy = firstStack.pop();
                if (remainingPoints.contains(xy)) {
                    ordered.set(idx1, xy.getX(), xy.getY());
                    remainingPoints.remove(xy);
                    idx1++;
                } else {
                    // fetch the previously set value to use in finding next neighbor below
                    xy = new PairInt(ordered.getX(idx1 - 1), ordered.getY(idx1 - 1));
                }
                
                if (!lastStack.isEmpty()) {
                    PairInt xy2 = lastStack.pop();
                    if (remainingPoints.contains(xy2)) {
                        ordered.set(idxn, xy2.getX(), xy2.getY());
                        remainingPoints.remove(xy2);
                        
                        for (int nIdx = 0; nIdx < dxs8.length; nIdx++) {
                            int vX = edge.getX(idxn) + dxs8[nIdx];
                            int vY = edge.getY(idxn) + dys8[nIdx];
                            PairInt p = new PairInt(vX, vY);
                            if (remainingPoints.contains(p)) {
                                lastStack.add(p);
                                break;
                            }
                        }
                        idxn--;
                    } else {
                        PairInt xy3 = new PairInt(ordered.getX(idxn - 1), ordered.getY(idxn - 1));
                        for (int nIdx = 0; nIdx < dxs8.length; nIdx++) {
                            int vX = xy3.getX() + dxs8[nIdx];
                            int vY = xy3.getY() + dys8[nIdx];
                            PairInt p = new PairInt(vX, vY);
                            if (remainingPoints.contains(p)) {
                                lastStack.add(p);
                                break;
                            }
                        }
                    } 
                }
                
                for (int nIdx = 0; nIdx < dxs8.length; nIdx++) {
                    int vX = xy.getX() + dxs8[nIdx];
                    int vY = xy.getY() + dys8[nIdx];
                    PairInt p = new PairInt(vX, vY);
                    if (remainingPoints.contains(p)) {
                        firstStack.add(p);
                        break;
                    }
                }
                
            } else if (!lastStack.isEmpty()) {
 
                PairInt xy = lastStack.pop();
                if (remainingPoints.contains(xy)) {
                    
                    ordered.set(idxn, xy.getX(), xy.getY());
                    remainingPoints.remove(xy);
                    
                    for (int nIdx = 0; nIdx < dxs8.length; nIdx++) {
                        int vX = edge.getX(idxn) + dxs8[nIdx];
                        int vY = edge.getY(idxn) + dys8[nIdx];
                        PairInt p = new PairInt(vX, vY);
                        if (remainingPoints.contains(p)) {
                            lastStack.add(p);
                            break;
                        }
                    }
                    idxn--;
                } else {
                    PairInt xy2 = new PairInt(ordered.getX(idxn - 1), ordered.getY(idxn - 1));
                    for (int nIdx = 0; nIdx < dxs8.length; nIdx++) {
                        int vX = xy2.getX() + dxs8[nIdx];
                        int vY = xy2.getY() + dys8[nIdx];
                        PairInt p = new PairInt(vX, vY);
                        if (remainingPoints.contains(p)) {
                            lastStack.add(p);
                            break;
                        }
                    }
                }
                
            } else if (!remainingPoints.isEmpty()) {
                
                // look for an adj point to last entered for idx1 or idxn
                if (idx1 <= midIdx) {
                    PairInt xy = new PairInt(ordered.getX(idx1 - 1), ordered.getY(idx1 - 1));
                    for (int nIdx = 0; nIdx < dxs8.length; nIdx++) {
                        int vX = xy.getX() + dxs8[nIdx];
                        int vY = xy.getY() + dys8[nIdx];
                        PairInt p = new PairInt(vX, vY);
                        if (remainingPoints.contains(p)) {
                            firstStack.add(p);
                            break;
                        }
                    }
                    
                } else if (idxn >= midIdx) {
                    PairInt xy = new PairInt(ordered.getX(idxn + 1), ordered.getY(idxn + 1));
                    for (int nIdx = 0; nIdx < dxs8.length; nIdx++) {
                        int vX = xy.getX() + dxs8[nIdx];
                        int vY = xy.getY() + dys8[nIdx];
                        PairInt p = new PairInt(vX, vY);
                        if (remainingPoints.contains(p)) {
                            lastStack.add(p);
                            break;
                        }
                    }
                } else {
                    throw new IllegalStateException("Error in algorithm!  indexes w.r.t. to midIdx have passed midpoint");
                }
                
            } else {
                return null;
            }
          
        }
        
        return ordered;
    }

    private void reorderEndpointsIfNeeded(PairIntArray out) {
        
        if (out == null || out.getN() < 3) {
            return;
        }
        
        boolean endPointsAreAdjacent = endPointsAreAdjacent(out);

        if (endPointsAreAdjacent) {
            return;
        }
        
        int x1 = out.getX(0);
        int y1 = out.getY(0);
        
        int xn = out.getX(out.getN() - 1);
        int yn = out.getY(out.getN() - 1);
        
        int closestPivotIdx = -1;
        
        for (int idx = (out.getN() - 2); idx > 0; --idx) {
            
            int x = out.getX(idx);
            int y = out.getY(idx);
            
            if ((Math.abs(x - x1) < 2) && (Math.abs(y - y1) < 2)) {
                closestPivotIdx = idx;
                break;
            }
        }
        
        if (closestPivotIdx == -1) {
            // it's not a closed curve
            return;
        }
      
        if (closestPivotIdx > (out.getN() - closestPivotIdx - 1)) {
            // closest to end
            /*
             2  1
           3      0
           4     7   
            5  6  8 
                 9  
            If point before closestPivotIdx is next to (xn, yn), can just 
            reverse the points from closestPivotIdx to n-1.
            */
            int xPrev = out.getX(closestPivotIdx - 1);
            int yPrev = out.getY(closestPivotIdx - 1);
            if ((Math.abs(xPrev - xn) > 1) && (Math.abs(yPrev - yn) > 1)) {
                return;
            }
            
            int count = 0;
            int nSep = (out.getN() - closestPivotIdx) >> 1;
            for (int idx = closestPivotIdx; idx < (closestPivotIdx + nSep); ++idx) {
                int idx2 = out.getN() - count - 1;
                int swapX = out.getX(idx);
                int swapY = out.getY(idx);
                out.set(idx, out.getX(idx2), out.getY(idx2));
                out.set(idx2, swapX, swapY);
                count++;
            }
            
        } else {
            // closest to beginning
            /*
             7  8
           6      9
           5     2   
            4  3  1 
                 0  
            If point after closestPivotIdx is next to (x0, y0), can just 
            reverse the points from 0 to nTo1ClosestIdx.
            */
            
            closestPivotIdx = -1;
        
            for (int idx = 1; idx < (out.getN() - 1); ++idx) {

                int x = out.getX(idx);
                int y = out.getY(idx);

                if ((Math.abs(x - xn) < 2) && (Math.abs(y - yn) < 2)) {
                    closestPivotIdx = idx;
                    break;
                }
            }

            if (closestPivotIdx == -1) {
                // it's not a closed curve
                return;
            }
            
            int xNext = out.getX(closestPivotIdx + 1);
            int yNext = out.getY(closestPivotIdx + 1);
            if ((Math.abs(xNext - x1) > 1) || (Math.abs(yNext - y1) > 1)) {
                return;
            }
            /*45, 15 <-- 0       99---> 45, 12
              44, 15     1
              45, 14     2
              44, 13   <-- pivot
            */
            /*
             8  9
           7      10
           6      3   
            5  4  2 
               1  0  
            If point after closestPivotIdx is next to (x0, y0), can just 
            reverse the points from 0 to nTo1ClosestIdx.
            */
            reverse0toIdx(out, closestPivotIdx);
        }
    }
    
    private void reverse0toIdx(PairIntArray out, int lastSwapIdx) {
        int nSep = (lastSwapIdx + 1) >> 1;
        for (int idx = 0; idx < nSep; ++idx) {
            int idx2 = lastSwapIdx - idx;
            int swapX = out.getX(idx);
            int swapY = out.getY(idx);
            out.set(idx, out.getX(idx2), out.getY(idx2));
            out.set(idx2, swapX, swapY);
        }
    }

    private boolean endPointsAreAdjacent(PairIntArray out) {
        
        int x1 = out.getX(0);
        int y1 = out.getY(0);
        
        int xn = out.getX(out.getN() - 1);
        int yn = out.getY(out.getN() - 1);
        
        if ((Math.abs(x1 - xn) > 1) || (Math.abs(y1 - yn) > 1)) {
            return false;
        }
        
        return true;
    }
    
}
