package algorithms.imageProcessing;

import algorithms.CountingSort;
import algorithms.MultiArrayMergeSort;
import algorithms.QuickSort;
import algorithms.util.PairIntArray;
import algorithms.util.PairInt;
import java.io.IOException;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

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
    
    private boolean debug = true;
    
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
        
printJoinPoints(joinPoints, output);
writeJoinPointsImage(joinPoints, output);
   
        //output = joinOnJoinPoints(joinPoints, output);
 
        //findJunctions(output);
        
        //printJunctions(output);

        //spliceEdgesAtJunctionsIfImproves(output);
    
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
        Set<JoinPointEntry> theJoinPoints = new HashSet<JoinPointEntry>(2 * n);

        //with key = edge number
        //with value = set of the entries in theJoinPoints for which one has a first end point in this edge
        Map<Integer, Set<JoinPointEntry>> edgeFirstEndPointMap = 
            new HashMap<Integer, Set<JoinPointEntry>>();
    
        Map<Integer, Set<JoinPointEntry>>  edgeLastEndPointMap = 
            new HashMap<Integer, Set<JoinPointEntry>> ();

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
        
        int[] dxs = new int[]{-1, -1,  0,  1, 1, 1, 0, -1};
        int[] dys = new int[]{ 0, -1, -1, -1, 0, 1, 1,  1};
        
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
            
            for (int nIdx = 0; nIdx < dxs.length; nIdx++) {
                int vX = uX + dxs[nIdx];
                int vY = uY + dys[nIdx];
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
            
            theJoinPoints.add(new JoinPointEntry(loc0, loc1));
            
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
                    log.warning(String.format(
                    "creating a join point for an edge with size 3: xy=(%d,%d) loc=%d:%d",
                    vClosestLoc.getX(), vClosestLoc.getY(), vClosestLoc.getX(), 
                    vClosestLoc.getY()));
                    // an entry is added to both edge maps because we don't know yet
                    // which to use.
                }

                if (addToFirst) {
                    Integer key = Integer.valueOf(locX);
                    Set<JoinPointEntry> joinPointSet = edgeFirstEndPointMap.get(key);
                    if (joinPointSet == null) {
                        joinPointSet = new HashSet<JoinPointEntry>();
                    }
                    JoinPointEntry joinPoint = new JoinPointEntry(loc0, loc1);
                    joinPointSet.add(joinPoint);
                    edgeFirstEndPointMap.put(key, joinPointSet);
                }
            
                if (addToLast) {
                    Integer key = Integer.valueOf(locX);
                    Set<JoinPointEntry> joinPointSet = edgeLastEndPointMap.get(key);
                    if (joinPointSet == null) {
                        joinPointSet = new HashSet<JoinPointEntry>();
                    }
                    JoinPointEntry joinPoint = new JoinPointEntry(loc0, loc1);
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
            assertConsistentJoinPointStructures(edges, edgeFirstEndPointMap,
                edgeLastEndPointMap, theJoinPoints, skipForSize3);
        }        
        
        // edges that are size 3 and have join points that are not in the
        // first or last position, can be moved to either, so keep
        // track of them to avoid moving to the same position twice
        Map<PairInt, Set<Integer>> edgeSize3NonEndPoints = new HashMap<PairInt, 
            Set<Integer>>();
        for (JoinPointEntry entry : theJoinPoints) {
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
        for (JoinPointEntry entry : theJoinPoints) {
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
            assertConsistentJoinPointStructures(edges, edgeFirstEndPointMap,
                edgeLastEndPointMap, theJoinPoints, skipForSize3);
        }
        
        Map<PairInt, PairInt> result = new HashMap<PairInt, PairInt>();
        for (JoinPointEntry entry : theJoinPoints) {
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
        
        int[] dxs = new int[]{-1, -1,  0,  1, 1, 1, 0, -1};
        int[] dys = new int[]{ 0, -1, -1, -1, 0, 1, 1,  1};
                
        // 8 * O(N)
        for (int edgeIdx = 0; edgeIdx < n; edgeIdx++) {
            
            PairIntArray edge = edges.get(edgeIdx);
            
            for (int iEdgeIdx = 0; iEdgeIdx < edge.getN(); iEdgeIdx++) {
                
                int col = edge.getX(iEdgeIdx);
                int row = edge.getY(iEdgeIdx);
                
                int uIdx = img.getIndex(col, row);
                
                Set<PairInt> neighbors = new HashSet<PairInt>();
                
                for (int nIdx = 0; nIdx < dxs.length; nIdx++) {
                    
                    int x = col + dxs[nIdx];
                    int y = row + dys[nIdx];
                    
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
                        
                        theJunctionLocationMap.put(Integer.valueOf(vIdx), p);
                        
                        indexes.add(Integer.valueOf(vIdx));
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
            
            int col = img.getCol(pixelIndex.intValue());
            int row = img.getRow(pixelIndex.intValue());
            
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
                // remove this pixel
                remove.add(pixelIndex);
                assert(!doNotRemove.contains(pixelIndex));
            } else {
                doNotRemove.add(pixelIndex);
                // remove the neighbors from the junction map
                for (Integer pixelIndex2 : neighborIndexes) { 
                    if (!doNotRemove.contains(pixelIndex2)) {
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
    
    private void printEdgeEndpointMap(List<PairIntArray> edges,
        Map<Integer, Set<JoinPointEntry>> edgeEndPointMap) {
        
        StringBuilder sb = new StringBuilder("");
        
        for (Entry<Integer, Set<JoinPointEntry>> entry : edgeEndPointMap.entrySet()) {
                        
            sb.append("Edge ").append(entry.getKey().toString()).append(":\n");
                        
            for (JoinPointEntry joinPoint : entry.getValue()) {
                
                PairInt loc0 = joinPoint.getKey();
                PairInt loc1 = joinPoint.getValue();
                                        
                PairIntArray edge0 = edges.get(loc0.getX());
                int x0 = edge0.getX(loc0.getY());
                int y0 = edge0.getY(loc0.getY());

                PairIntArray edge1 = edges.get(loc1.getX());
                int x1 = edge1.getX(loc1.getY());
                int y1 = edge1.getY(loc1.getY());

                sb.append(String.format("  (%d,%d) to (%d,%d) in edges %d and %d  at positions=%d out of %d and %d out of %d\n",
                    x0, y0, x1, y1, loc0.getX(), loc1.getX(), 
                    loc0.getY(), edges.get(loc0.getX()).getN(),
                    loc1.getY(), edges.get(loc1.getX()).getN()
                ));
            }
        }
        
        log.info(sb.toString());
    }

    private void printJoinPoints(Map<PairInt, PairInt> joinPoints,
        List<PairIntArray> edges) {
        
        StringBuilder sb = new StringBuilder("join points:\n");
        
        for (Entry<PairInt, PairInt> entry : joinPoints.entrySet()) {
            
            PairInt loc0 = entry.getKey();
            
            PairInt loc1 = entry.getValue();
            
            PairIntArray edge0 = edges.get(loc0.getX());
            int x0 = edge0.getX(loc0.getY());
            int y0 = edge0.getY(loc0.getY());
            
            PairIntArray edge1 = edges.get(loc1.getX());
            int x1 = edge1.getX(loc1.getY());
            int y1 = edge1.getY(loc1.getY());
            
            sb.append(String.format("  (%d,%d) to (%d,%d) in edges %d and %d  at positions=%d out of %d and %d out of %d\n",
                x0, y0, x1, y1, loc0.getX(), loc1.getX(), 
                loc0.getY(), edges.get(loc0.getX()).getN(),
                loc1.getY(), edges.get(loc1.getX()).getN()
            ));
        }
        
        log.info(sb.toString());
    }
    
    private void printJoinPoints(Set<JoinPointEntry> joinPoints,
        List<PairIntArray> edges) {
        
        StringBuilder sb = new StringBuilder("join points:\n");
        
        for (JoinPointEntry entry : joinPoints) {
            
            PairInt loc0 = entry.getKey();
            
            PairInt loc1 = entry.getValue();
            
            PairIntArray edge0 = edges.get(loc0.getX());
            int x0 = edge0.getX(loc0.getY());
            int y0 = edge0.getY(loc0.getY());
            
            PairIntArray edge1 = edges.get(loc1.getX());
            int x1 = edge1.getX(loc1.getY());
            int y1 = edge1.getY(loc1.getY());
            
            sb.append(String.format("  (%d,%d) to (%d,%d) in edges %d and %d  at positions=%d out of %d and %d out of %d\n",
                x0, y0, x1, y1, loc0.getX(), loc1.getX(), 
                loc0.getY(), edges.get(loc0.getX()).getN(),
                loc1.getY(), edges.get(loc1.getX()).getN()
            ));
        }
        
        log.info(sb.toString());
    }
    
    private void printJoinPoints(PairInt[][] edgeJoins, int idxLo, int idxHi,
        Map<Integer, PairIntArray> edges) {
        
        // print array indexes
        int[] aIndexes = new int[edges.size()];
        int count = 0;
        for (Entry<Integer, PairIntArray> entry : edges.entrySet()) {
            Integer edgeIndex = entry.getKey();
            aIndexes[count] = edgeIndex.intValue();
            count++;
        }
        CountingSort.sort(aIndexes, 2*edges.size());
        
        StringBuilder sb = new StringBuilder("output indexes of size ");
        sb.append(Integer.toString(edges.size())).append("\n");
        for (int i = 0; i < aIndexes.length; i++) {
            sb.append(Integer.toString(aIndexes[i])).append(" ");
        }
        log.info(sb.toString());
        
        sb = new StringBuilder("join points:\n");
        
        for (int i = idxLo; i <= idxHi; i++) {
                        
            PairInt loc0 = edgeJoins[i][0];
            PairInt loc1 = edgeJoins[i][1];

            PairIntArray edge0 = edges.get(Integer.valueOf(loc0.getX()));
                  
            int x0 = edge0.getX(loc0.getY());
            int y0 = edge0.getY(loc0.getY());
            
            PairIntArray edge1 = edges.get(Integer.valueOf(loc1.getX()));
            int x1 = edge1.getX(loc1.getY());
            int y1 = edge1.getY(loc1.getY());
            
            sb.append(String.format("  (%d,%d) to (%d,%d) in edges %d and %d  at positions=%d out of %d and %d out of %d\n",
                x0, y0, x1, y1, loc0.getX(), loc1.getX(), 
                loc0.getY(), edge0.getN(),
                loc1.getY(), edge1.getN()
            ));
        }
        
        log.info(sb.toString());
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
     * changes.
     * @param joinPoints
     * @param edges
     * @return 
     */
    protected List<PairIntArray> joinOnJoinPoints(Map<PairInt, PairInt> 
        joinPoints, List<PairIntArray> edges) {
        
        //order the join points
        
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
           
//printJoinPoints(edgeJoins, 0, i, edgesMap);
           
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
    
    private void writeJoinPointsImage(Map<PairInt, PairInt> theJoinPointMap, 
        List<PairIntArray> edges) {
        
        try {
            
            Image img2 = img.copyImageToGreen();

            ImageIOHelper.addAlternatingColorCurvesToImage(edges, img2);
            int nExtraForDot = 1;
            int rClr = 255;
            int gClr = 0;
            int bClr = 100;
            for (Entry<PairInt, PairInt> entry : theJoinPointMap.entrySet()) {
                PairInt loc0 = entry.getKey();
                PairInt loc1 = entry.getValue();
                assert(!loc0.equals(loc1));

                PairIntArray edge0 = edges.get(loc0.getX());
                PairIntArray edge1 = edges.get(loc1.getX());
            
                ImageIOHelper.addPointToImage(
                    edge0.getX(loc0.getY()), edge0.getY(loc0.getY()), 
                    img2, nExtraForDot,
                    rClr, gClr, bClr);
                
                ImageIOHelper.addPointToImage(
                    edge1.getX(loc1.getY()), edge1.getY(loc1.getY()), 
                    img2, nExtraForDot,
                    rClr, gClr, bClr);
            }

            String dirPath = algorithms.util.ResourceFinder.findDirectory("bin");
            String sep = System.getProperty("file.separator");
            ImageIOHelper.writeOutputImage(dirPath + sep + "joinpoints.png", img2);
            
        } catch (IOException e) {
            
        }
    }

    private void printJunctions(List<PairIntArray> edges) {
        printJunctions(this.junctionMap, edges);
    }
    
    private void printJunctions(Map<Integer, Set<Integer>> jMap, 
        List<PairIntArray> edges) {
        
        try {
            
            Image img2 = img.copyImageToGreen();

            ImageIOHelper.addAlternatingColorCurvesToImage(edges, img2);
            int nExtraForDot = 1;
            int rClr = 255;
            int gClr = 0;
            int bClr = 100;
            for (Entry<Integer, Set<Integer>> entry : jMap.entrySet()) {
                int pixIdx = entry.getKey().intValue();
                int col = img2.getCol(pixIdx);
                int row = img2.getRow(pixIdx);
                ImageIOHelper.addPointToImage(col, row, img2, nExtraForDot,
                    rClr, gClr, bClr);
            }

            String dirPath = algorithms.util.ResourceFinder.findDirectory("bin");
            String sep = System.getProperty("file.separator");
            ImageIOHelper.writeOutputImage(dirPath + sep + "junctions.png", img2);
            
        } catch (IOException e) {
            
        }
    }
    
    private String printJunctionsToString(Map<Integer, Set<Integer>> jMap, 
        List<PairIntArray> edges) {

        StringBuilder sb = new StringBuilder("junctions:\n");
        
        for (Entry<Integer, Set<Integer>> entry : jMap.entrySet()) {
            int pixIdx = entry.getKey().intValue();
            int col = img.getCol(pixIdx);
            int row = img.getRow(pixIdx);
            sb.append(String.format("%d (%d,%d)\n", pixIdx, col, row));
        }
        
        return sb.toString();
    }
    
    private String printJunctionsToString(
        Map<Integer, PairInt> jLocationMap, Map<Integer, Set<Integer>> jMap, 
        List<PairIntArray> edges) {

        StringBuilder sb = new StringBuilder("junctions:\n");
        
        for (Entry<Integer, Set<Integer>> entry : jMap.entrySet()) {
            int pixIdx = entry.getKey().intValue();
            int col = img.getCol(pixIdx);
            int row = img.getRow(pixIdx);
            sb.append(String.format("%d (%d,%d)\n", pixIdx, col, row));
        }
        
        sb.append("junction locations:\n");
        for (Entry<Integer, PairInt> entry : jLocationMap.entrySet()) {
            Integer pixelIndex = entry.getKey();
            PairInt loc = entry.getValue();
            int edgeIdx = loc.getX();
            int indexWithinEdge = loc.getY();
            int edgeN = edges.get(edgeIdx).getN();
            int x = edges.get(edgeIdx).getX(indexWithinEdge);
            int y = edges.get(edgeIdx).getY(indexWithinEdge);
            sb.append(String.format("edge=%d idx=%d (out of %d) pixIdx=%d (%d,%d)\n", 
                edgeIdx, indexWithinEdge, edgeN, pixelIndex.intValue(), x, y));
        }//edgeIdx=49 indexWithinEdge=28
        
        return sb.toString();
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

    private void spliceEdgesAtJunctionsIfImproves(List<PairIntArray> edges) {
               
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
        
        // key = edge index.  
        // value = pixel indexes.
        //   the pixel indexes are used to find values in junctionLocatorMap
        //   to update it as points are moved to and from edges.
        Map<Integer, Set<Integer>> theEdgeToPixelIndexMap = 
            new HashMap<Integer, Set<Integer>>();
        
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
        
        if (debug) {
            assertConsistentEdgeCapacity(theEdgeToPixelIndexMap, 
                junctionLocationMap, junctionMap, edges);
        }
        
        for (Entry<Integer, Set<Integer>> entry : junctionMap.entrySet()) {
                        
            Integer centerPixelIndex = entry.getKey();
            PairInt centerLoc = junctionLocationMap.get(centerPixelIndex);
            assert(centerLoc != null);
        
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
            
            int count = 1;
            
            for (Integer pixIndex : adjIndexes) {
                
                pixIndexes[count] = pixIndex.intValue();
                
                PairInt loc = junctionLocationMap.get(pixIndex);
                
                if ((centerLoc.getX() == loc.getX())) {
                    count++;
                    continue;
                }
                                                
                lengths[count] = (new Splice(edges.get(loc.getX())))
                    .getLengthOfLongestSide(loc.getY());
           
                if (lengths[count] > maxN) {
                    maxN = lengths[count];
                }
                
                count++;
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
                    assertConsistentEdgeCapacity(theEdgeToPixelIndexMap, 
                        junctionLocationMap, junctionMap, edges);
                }
            }
        }        
    }

    private void assertConsistentEdgeCapacity(
        Map<Integer, Set<Integer>> theEdgeToPixelIndexMap, 
        Map<Integer, PairInt> junctionLocationMap, 
        Map<Integer, Set<Integer>> junctionMap, 
        List<PairIntArray> edges) {
        
        for (Entry<Integer, PairInt> entry : junctionLocationMap.entrySet()) {
            Integer pixelIndex = entry.getKey();
            assert(pixelIndex != null);
            PairInt loc = entry.getValue();
            int edgeIdx = loc.getX();
            int idx = loc.getY();
            PairIntArray edge = edges.get(edgeIdx);
            assert(edge != null);
            assert (idx < edge.getN()) :
                "idx=" + idx + " edgeN=" + edge.getN() + " edgeIndex=" + edgeIdx;
        }
        /*
        previous edge 50 had n=24
        current edge 49 has n=24.
        
        java.lang.AssertionError: idx=28 edgeN=24 edgeIndex=49
        
        
        processing junction w/ center pixel index=50427 and loc=3:15
        ...[junit] edge=49 idx=28 (out of 39) pixIdx=42810 (13,247)
        
        before splice edge 4 (137 points) to edge 3 (16 points)
        
        -----------
        splice edge 49 (29 points) to edge 50 (8 points), that is append 50 to end of 49
        50 edge size is 8
        49 edge size is 49.  spliced to 29 and 20
          'off by 1'?  spliced last point is idx=28 within edge 49 before...
        
        splice0_0: update Y for pixIdx=42810   loc 49:28 to 49:28 (edgeN=29)  (**make sure edgeIdx is same)
        
        */
        
        for (Entry<Integer, Set<Integer>> entry : junctionMap.entrySet()) {
            
            Integer pixelIndex = entry.getKey();
            Set<Integer> adjPixelIndexes = entry.getValue();
            
            PairInt loc = junctionLocationMap.get(pixelIndex);
            assert(loc != null);
            
            int edgeIdx = loc.getX();
            int idx = loc.getY();
            PairIntArray edge = edges.get(edgeIdx);
            assert(edge != null);
            assert(idx < edge.getN());
            
            for (Integer adjPixelIndex : adjPixelIndexes) {
                PairInt adjLoc = junctionLocationMap.get(adjPixelIndex);
                assert(adjLoc != null);
                edgeIdx = adjLoc.getX();
                idx = adjLoc.getY();
                edge = edges.get(edgeIdx);
                assert(edge != null);
                assert(idx < edge.getN());
            }
        }
        
        for (Entry<Integer, Set<Integer>> entry : theEdgeToPixelIndexMap.entrySet()) {
            Integer edgeIndex = entry.getKey();
            Set<Integer> pixelIndexes = entry.getValue();
            
            PairIntArray edge = edges.get(edgeIndex.intValue());
            assert(edge != null);
            
            for (Integer pixelIndex : pixelIndexes) {
                PairInt loc = junctionLocationMap.get(pixelIndex);
                assert(loc != null);

                int edgeIdx = loc.getX();
                int idx = loc.getY();
                
                assert(edgeIdx == edgeIndex.intValue());
                
                assert(idx < edge.getN());
            }
        }
    }

    private void assertConsistentJoinPointStructures(
        List<PairIntArray> edges,
        Map<Integer, Set<JoinPointEntry>> edgeFirstEndPointMap, 
        Map<Integer, Set<JoinPointEntry>> edgeLastEndPointMap, 
        Set<JoinPointEntry> theJoinPoints, boolean skipForSize3) {
        
        Set<PairInt> joinPointsSet = new HashSet<PairInt>();
        for (JoinPointEntry entry : theJoinPoints) {
            PairInt loc0 = entry.getKey();
            PairInt loc1 = entry.getValue();
            assert(!loc0.equals(loc1));
            
            PairIntArray edge0 = edges.get(loc0.getX());
            PairIntArray edge1 = edges.get(loc1.getX());
            int n0 = edge0.getN();
            int n1 = edge1.getN();
            
            assert(loc0.getY() < n0);
            assert(loc1.getY() < n1);
            
            if (!(skipForSize3 && (n0 == 3))) {
                assert(!joinPointsSet.contains(loc0));
            }
            if (!(skipForSize3 && (n1 == 3))) {
                assert(!joinPointsSet.contains(loc1));
            }
            
            joinPointsSet.add(loc0);
            joinPointsSet.add(loc1);
        }
        joinPointsSet.clear();
        for (JoinPointEntry entry : theJoinPoints) {
            PairInt loc0 = entry.getKey();
            PairInt loc1 = entry.getValue();
            assert(!loc0.equals(loc1));
            assert(!joinPointsSet.contains(loc0));
            assert(!joinPointsSet.contains(loc1));
            joinPointsSet.add(loc0);
            joinPointsSet.add(loc1);
        }

        for (Entry<Integer, Set<JoinPointEntry>> entry : edgeFirstEndPointMap.entrySet()) {
            assert(entry.getValue().size() < 2);
            for (JoinPointEntry joinPoint : entry.getValue()) {
                PairInt loc0 = joinPoint.getKey();
                PairInt loc1 = joinPoint.getValue();
                assert(loc0 != null);
                assert(loc1 != null);
                                
                boolean contains = theJoinPoints.contains(joinPoint);                   
                if (!contains) {
                    int hash = joinPoint.hashCode();
                    for (JoinPointEntry entry2 : theJoinPoints) {
                        //log.info("   hash=" + entry2.hashCode());
                        if (hash == entry2.hashCode()) {
                            contains = entry2.equals(joinPoint);
                            //log.info("Set's use of HashMap.getEntry() did not find this point.");
                        }
                    }
                }
                assert(contains);
            }
        }
        for (Entry<Integer, Set<JoinPointEntry>> entry : edgeLastEndPointMap.entrySet()) {
            assert(entry.getValue().size() < 2);
            for (JoinPointEntry joinPoint : entry.getValue()) {
                PairInt loc0 = joinPoint.getKey();
                PairInt loc1 = joinPoint.getValue();
                assert(loc0 != null);
                assert(loc1 != null);
                boolean contains = theJoinPoints.contains(joinPoint);                   
                if (!contains) {
                    int hash = joinPoint.hashCode();
                    for (JoinPointEntry entry2 : theJoinPoints) {
                        //log.info("   hash=" + entry2.hashCode());
                        if (hash == entry2.hashCode()) {
                            contains = entry2.equals(joinPoint);
                            //log.info("Set's use of HashMap.getEntry() did not find this point.");
                        }
                    }
                }
                assert(contains);
            }
        }
    }

    protected void reduceMultipleEndpointsForEdge(
        List<PairIntArray> edges,
        Map<Integer, Set<JoinPointEntry>> edgeFirstEndPointMap,
        Map<Integer, Set<JoinPointEntry>> edgeLastEndPointMap, 
        Map<PairInt, PairInt> endPointMap, Set<JoinPointEntry> theJoinPoints) {
        
        Map<Integer, Set<JoinPointEntry>> tmpFirstRemoveMap = new HashMap<Integer, Set<JoinPointEntry>>();
        
        Map<Integer, Set<JoinPointEntry>> tmpLastRemoveMap = new HashMap<Integer, Set<JoinPointEntry>>();
        
        for (int type = 0; type < 2; ++type) {
        
            Map<Integer, Set<JoinPointEntry>> edgeMap = edgeFirstEndPointMap;
            if (type == 1) {
                edgeMap = edgeLastEndPointMap;
            }
            
            for (Entry<Integer, Set<JoinPointEntry>> entry : edgeMap.entrySet()) {

                Set<JoinPointEntry> joinPointSet = entry.getValue();

                if (joinPointSet.size() < 2) {
                    continue;
                }

                // has more than one joint point for this endpoint, though some
                // might have already been removed and are in tmpFirstRemoveMap
                // if so
        
                // construct a warning
                StringBuilder sb = new StringBuilder("edge ");
                sb.append(entry.getKey().toString())
                    .append(" has ").append(Integer.toString(joinPointSet.size()))
                    .append(" first endpoints, so deciding between them:");
                for (JoinPointEntry joinPoint : entry.getValue()) {
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
            
                // choose closest join point pair, and break ties with those that
                // have higher number of members that are already endpoints.

                JoinPointEntry closestJoinPoint = null;
                int closestDistSq = Integer.MAX_VALUE;
                int closestCanBeReordered0 = -99;
                int closestCanBeReordered1 = -99;
                
                for (JoinPointEntry joinPoint : joinPointSet) {

                    PairInt loc0 = joinPoint.getKey();
                    PairIntArray edge0 = edges.get(loc0.getX());
                    int n0 = edge0.getN();
                    
                    Set<JoinPointEntry> removedFirstSet = 
                        tmpFirstRemoveMap.get(Integer.valueOf(loc0.getX()));
                    if (removedFirstSet == null) {
                        removedFirstSet = new HashSet<JoinPointEntry>();
                        tmpFirstRemoveMap.put(Integer.valueOf(loc0.getX()), removedFirstSet);
                    }
                    
                    Set<JoinPointEntry> removedLastSet = 
                        tmpLastRemoveMap.get(Integer.valueOf(loc0.getX()));
                    if (removedLastSet == null) {
                        removedLastSet = new HashSet<JoinPointEntry>();
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

                assert (closestJoinPoint != null);

                for (JoinPointEntry joinPoint : joinPointSet) {

                    if (joinPoint.equals(closestJoinPoint)) {
                        continue;
                    }

                    PairInt loc0 = joinPoint.getKey();
                    PairIntArray edge0 = edges.get(loc0.getX());
                    int n0 = edge0.getN();

                    Set<JoinPointEntry> removedFirstSet0 = 
                        tmpFirstRemoveMap.get(Integer.valueOf(loc0.getX()));
                   
                    Set<JoinPointEntry> removedLastSet0 = 
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
                        log.warning(String.format(
                            "removing point (%d,%d) but it is in the middle of an edge of size 3",
                            edge0.getX(1), edge0.getY(1)));
                    }
                    
                    Set<JoinPointEntry> removedFirstSet1 = 
                        tmpFirstRemoveMap.get(Integer.valueOf(loc1.getX()));
                    if (removedFirstSet1 == null) {
                        removedFirstSet1 = new HashSet<JoinPointEntry>();
                        tmpFirstRemoveMap.put(Integer.valueOf(loc1.getX()), removedFirstSet1);
                    }
                    
                    Set<JoinPointEntry> removedLastSet1 = 
                        tmpLastRemoveMap.get(Integer.valueOf(loc1.getX()));
                    if (removedLastSet1 == null) {
                        removedLastSet1 = new HashSet<JoinPointEntry>();
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
            
            Map<Integer, Set<JoinPointEntry>> tmpMap = new HashMap<Integer, Set<JoinPointEntry>>();
            
            for (Entry<Integer, Set<JoinPointEntry>> entry : edgeFirstEndPointMap.entrySet()) {
                
                Set<JoinPointEntry> rmJoinPointSet = tmpFirstRemoveMap.get(entry.getKey());
                
                if (rmJoinPointSet == null || rmJoinPointSet.isEmpty()) {
                    tmpMap.put(entry.getKey(), entry.getValue());
                    continue;
                }
                
                Set<JoinPointEntry> tmpJoinPointSet = new HashSet<JoinPointEntry>();
                for (JoinPointEntry joinPoint : entry.getValue()) {
                    if (!rmJoinPointSet.contains(joinPoint)) {
                        tmpJoinPointSet.add(joinPoint);
                    }
                }
                tmpMap.put(entry.getKey(), tmpJoinPointSet);
            }
            edgeFirstEndPointMap.clear();
            edgeFirstEndPointMap.putAll(tmpMap);
            
            tmpMap = new HashMap<Integer, Set<JoinPointEntry>>();
            
            for (Entry<Integer, Set<JoinPointEntry>> entry : edgeLastEndPointMap.entrySet()) {
                
                Set<JoinPointEntry> rmJoinPointSet = tmpLastRemoveMap.get(entry.getKey());
                
                if (rmJoinPointSet == null || rmJoinPointSet.isEmpty()) {
                    tmpMap.put(entry.getKey(), entry.getValue());
                    continue;
                }
                
                Set<JoinPointEntry> tmpJoinPointSet = new HashSet<JoinPointEntry>();
                for (JoinPointEntry joinPoint : entry.getValue()) {
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

    protected static class JoinPointEntry {
        
        private final PairInt key;
        private PairInt value;

        public JoinPointEntry(PairInt theKey, PairInt theValue) {
            
            if (theKey == null) {
                throw new IllegalArgumentException("theKey cannot be null");
            }
            if (theValue == null) {
                throw new IllegalArgumentException("theValue cannot be null");
            }
            
            this.key = theKey;
            this.value = theValue;
        }

        public PairInt getKey() {
            return key;
        }

        public PairInt getValue() {
            return value;
        }

        @Override
        public boolean equals(Object obj) {
            
            if (!(obj instanceof JoinPointEntry)) {
                return false;
            }

            JoinPointEntry other = (JoinPointEntry)obj;
            
            return (other.getKey().equals(key) && other.getValue().equals(value));
        }
        
        
        @Override
        public int hashCode() {

            int hash = fnvHashCode(this.key, this.value);

            return hash;
        }

        protected static int fnv321aInit = 0x811c9dc5;
        protected static int fnv32Prime = 0x01000193;

        protected int fnvHashCode(PairInt i0, PairInt i1) {

            /*
             * hash = offset_basis
             * for each octet_of_data to be hashed
             *     hash = hash xor octet_of_data
             *     hash = hash * FNV_prime
             * return hash
             *
             * Public domain:  http://www.isthe.com/chongo/src/fnv/hash_32a.c
             */
            int hash = 0;

            int sum = fnv321aInit;

            // xor the bottom with the current octet.
            sum ^= i0.getX();

            // multiply by the 32 bit FNV magic prime mod 2^32
            sum *= fnv32Prime;
            
            sum ^= i0.getY();
            sum *= fnv32Prime;

            sum ^= i1.getX();
            sum *= fnv32Prime;
            
            sum ^= i1.getY();
            sum *= fnv32Prime;

            hash = sum;

            return hash;
        }
    }
}
