package algorithms.imageProcessing;

import algorithms.CountingSort;
import algorithms.MultiArrayMergeSort;
import algorithms.QuickSort;
import algorithms.misc.MiscMath;
import algorithms.util.PairIntArray;
import algorithms.util.PairInt;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.Stack;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Edge extractor operates on an image that has already been reduced to 
 * single pixel width lines and extracts edges from it, attempting to make
 * the longest edges it can - this specific edge extractor stores and uses
 * junction information, so is not as fast as EdgeExtractor.java.
 * 
 Edge extraction
    Local Methods:
        (1) DFS walk through connected pixel to form a sequence of pixels called
            an edge.
          
        (2) find join points
        
        (3) join edges using join points
        
        (4) find junction points
         
        (5)
            
        (6) find edge endpoints which are separated from one another by a gap of
            one and fill in the gap while merging the edges.
         
        (7) remove edges shorter than a minimum length
        
  @see AbstractEdgeExtractor
  * 
 * @author nichole
 */
public class EdgeExtractorWithJunctions extends AbstractEdgeExtractor {
    
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
        
//printJoinPoints(joinPoints, output);
        
        output = joinOnJoinPoints(joinPoints, output);
  /*      
        findJunctions(output);
        
can see that can use the junction points next to make decisions
on whether to divide an edge to make longer edges
printJunctions();
        
    */    
        return output;
    }
    
    /**
     * Iterate over each point looking for its neighbors and noting when
     * it has 2 that are not in it's edge as a join point.
     * The results are stored in member variable joinPoint.
     * 
     * runtime complexity is O(N)
     * 
     * @param edges 
     * @return hashmap with key = edge index and index within edge of pixel;
     * value = the adjacent pixel's edge index and index within its edge
     */
    protected Map<PairInt, PairInt> findJoinPoints(List<PairIntArray> edges) {
        
        // join points for adjacent edge endPoints that are not junctions
        // key = edge index and index within edge of pixel
        // value = the adjacent pixel's edge index and index within its edge
        Map<PairInt, PairInt> theJoinPoints = new HashMap<PairInt, PairInt>();
        
        int n = edges.size();
        
        // key = image pixel index, 
        // value = pairint of edge index, and index within edge
        Map<Integer, PairInt> pointLocator = new HashMap<Integer, PairInt>();
        
        // O(N)
        for (int edgeIdx = 0; edgeIdx < n; edgeIdx++) {
            
            PairIntArray edge = edges.get(edgeIdx);
            
            for (int iEdgeIdx = 0; iEdgeIdx < edge.getN(); iEdgeIdx++) {
                
                int pixIdx = img.getIndex(edge.getX(iEdgeIdx), edge.getY(iEdgeIdx));
                
                pointLocator.put(Integer.valueOf(pixIdx), new PairInt(edgeIdx, 
                    iEdgeIdx));
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
                
                int pixIdx = img.getIndex(col, row);
                
                Set<PairInt> neighbors = new HashSet<PairInt>();
                int nDifferentEdgeThanUEdge = 0;
                
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
                        
                        if (vLoc.getX() != edgeIdx) {
                            nDifferentEdgeThanUEdge++;
                        }
                    }
                }
                
                // if there is a junction, but there is actually only one
                // other edge aside from edge, this is a join point w/ the closest point
                if ((neighbors.size() > 2) && (nDifferentEdgeThanUEdge == 1)) {
                    // we want a junction point between the closest members of edge
                    // to the point in the other edge.
                    // the test for this each time one of the same points is
                    // iEdgeIdx and the write of the result should be idempotent.
                    
                    PairInt nonEdge = null;
                    for (PairInt p : neighbors) {
                        if (p.getX() != edgeIdx) {
                            nonEdge = p;
                            break;
                        }
                    }
                    int nonEdgeX = edges.get(nonEdge.getX()).getX(nonEdge.getY());
                    int nonEdgeY = edges.get(nonEdge.getX()).getY(nonEdge.getY());
                    
                    // init w/ the "u" value
                    int closestIEdgeIdx = iEdgeIdx;
                    int closestDistSq = ((col - nonEdgeX)*(col - nonEdgeX)) + 
                        ((row - nonEdgeY)*(row - nonEdgeY));
                    
                    for (PairInt p : neighbors) {
                        if (p.getX() == edgeIdx) {
                            int x2 = edges.get(p.getX()).getX(p.getY());
                            int y2 = edges.get(p.getX()).getY(p.getY());
                            int diffX = x2 - nonEdgeX;
                            int diffY = y2 - nonEdgeY;
                            int distSq = (diffX * diffX) + (diffY * diffY);
                            if (distSq < closestDistSq) {
                                closestIEdgeIdx = p.getY();
                                closestDistSq = distSq;
                            }
                        }
                    }
                    
                    theJoinPoints.put(new PairInt(edgeIdx, closestIEdgeIdx), 
                        nonEdge);
                     
                } else if (neighbors.size() == 2) {
                    
                    // if they are not in the same edge already, this is a
                    // potential join point
                    
                    // key = edge index
                    // value = set of pairints of edge index and index within edge
                    Map<Integer, Set<PairInt>> edgeIndexLocatorMap 
                        = new HashMap<Integer, Set<PairInt>>();
                    
                    Set<PairInt> v = new HashSet<PairInt>();
                    v.add(new PairInt(edgeIdx, iEdgeIdx));
                    
                    edgeIndexLocatorMap.put(Integer.valueOf(edgeIdx), v);
                    
                    for (PairInt p : neighbors) {
                        
                        int eIdx = p.getX();
                        int iEIdx = p.getY();
                        
                        v = edgeIndexLocatorMap.get(Integer.valueOf(eIdx));
                        
                        if (v == null) {
                            v = new HashSet<PairInt>();
                        }
                        
                        v.add(new PairInt(eIdx, iEIdx));
                        
                        edgeIndexLocatorMap.put(Integer.valueOf(eIdx), v);
                    }
                    
                    if (edgeIndexLocatorMap.size() >= 2) {
                        
                        // store joinPoints for pairs of points in different 
                        // edges and adjacent
                        Iterator<Entry<Integer, Set<PairInt>>> iter = 
                            edgeIndexLocatorMap.entrySet().iterator();
                        
                        while (iter.hasNext()) {
                            
                            Entry<Integer, Set<PairInt>> entry = iter.next();
                            
                            int e0Idx = entry.getKey().intValue();
                            
                            for (PairInt p0 : entry.getValue()) {
                                
                                Iterator<Entry<Integer, Set<PairInt>>> iter2 = 
                                    edgeIndexLocatorMap.entrySet().iterator();
                                
                                Entry<Integer, Set<PairInt>> entry2 = iter2.next();
                                
                                if (entry.equals(entry2)) {
                                    continue;
                                }
                                
                                int e2Idx = entry2.getKey().intValue();
                                
                                if (e0Idx == e2Idx) {
                                    continue;
                                }
                                
                                //TODO: need to improve decision when not appending to first or last
                                // they might all need to be distance based, like the last
                                // else block
                                
                                int iEIdx = p0.getY();
                                PairIntArray e0 = edges.get(e0Idx);
                                int x0 = e0.getX(iEIdx);
                                int y0 = e0.getY(iEIdx);
                                int lastIdx0 = e0.getN() - 1;
                                if ((iEIdx != 0) && (iEIdx != lastIdx0)) {
                                    int dIdxStart = iEIdx;
                                    int dIdxStop = lastIdx0 - iEIdx;
                                    int reassignIdx = 0;
                                    if (dIdxStart < dIdxStop) {
                                        reassignIdx = 0;
                                    } else if (dIdxStart < dIdxStop) {
                                        reassignIdx = lastIdx0;
                                    } else {
                                        // choose by proximity
                                        int diffFirstX = e0.getX(0) - x0;
                                        int diffFirstY = e0.getY(0) - y0;
                                        int distFirstSq = (diffFirstX*diffFirstX) + (diffFirstY*diffFirstY);
                                        int diffLastX = e0.getX(lastIdx0) - x0;
                                        int diffLastY = e0.getY(lastIdx0) - y0;
                                        int distLastSq = (diffLastX*diffLastX) + (diffLastY*diffLastY);
                                        if (distFirstSq < distLastSq) {
                                            reassignIdx = 0;
                                        } else {
                                            reassignIdx = lastIdx0;
                                        }
                                    }
                                    p0.setY(reassignIdx);
                                    iEIdx = reassignIdx;
                                    x0 = e0.getX(iEIdx);
                                    y0 = e0.getY(iEIdx);
                                }
                                
                                for (PairInt p2 : entry2.getValue()) {
                                    int iE2Idx = p2.getY();
                                    PairIntArray e2 = edges.get(e2Idx);
                                    int x2 = e2.getX(iE2Idx);
                                    int y2 = e2.getY(iE2Idx);
                                    int lastIdx2 = e2.getN() - 1;                                    
                                    if ((iE2Idx != 0) && (iE2Idx != lastIdx2)) {
                                        int dIdxStart = iE2Idx;
                                        int dIdxStop = lastIdx2 - iE2Idx;
                                        int reassignIdx = 0;
                                        if (dIdxStart < dIdxStop) {
                                            reassignIdx = 0;
                                        } else if (dIdxStart < dIdxStop) {
                                            reassignIdx = lastIdx2;
                                        } else {
                                            // choose by proximity
                                            int diffFirstX = e2.getX(0) - x2;
                                            int diffFirstY = e2.getY(0) - y2;
                                            int distFirstSq = (diffFirstX*diffFirstX) + (diffFirstY*diffFirstY);
                                            int diffLastX = e2.getX(lastIdx2) - x2;
                                            int diffLastY = e2.getY(lastIdx2) - y2;
                                            int distLastSq = (diffLastX*diffLastX) + (diffLastY*diffLastY);
                                            if (distFirstSq < distLastSq) {
                                                reassignIdx = 0;
                                            } else {
                                                reassignIdx = lastIdx2;
                                            }
                                        }
                                        p2.setY(reassignIdx);
                                        iE2Idx = reassignIdx;
                                        x2 = e2.getX(iE2Idx);
                                        y2 = e2.getY(iE2Idx);
                                    }                                    
                                    
                                    int diffX = Math.abs(x2 - x0);
                                    int diffY = Math.abs(y2 - y0);
                                    if ((diffX < 2) && (diffY < 2)) {
                                        theJoinPoints.put(p0, p2);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        
        return theJoinPoints;
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
                // remove this pixel
                remove.add(pixelIndex);
            }
        }
        
        for (Integer pixelIndex : remove) {
            theJunctionMap.remove(pixelIndex);
        }
        
        junctionMap = theJunctionMap;
        junctionLocationMap = theJunctionLocationMap;
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
            
            // edge to move
            PairIntArray edge0 = edgesMap.remove(Integer.valueOf(loc0.getX()));
            int removedEdgeIdx = loc0.getX();

            int n0 = edge0.getN();
            // join point should be at the beginning, so reverse if not
            if (loc0.getY() != 0) {
                
                edge0.reverse();
                loc0.setY(0);
                
                // everything with smaller index than i in edgeJoins with
                // edgeIndex==loc0.getX() needs to be updated for this reversal.
                // idx becomes n-idx-1
                for (int j = (i - 1); j > -1; --j) {
                    PairInt[] vEntry = edgeJoins[j];
                    for (int k = 0; k < 2; k++) {
                        PairInt vLoc = vEntry[k];
                        if (vLoc.getX() == loc0.getX()) {
                            int nV = edge0.getN();
                            vLoc.setY(nV - vLoc.getY() - 1);
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
                
                // everything with smaller index than i in edgeJoins that is from
                // edgeIndex==loc1.getX() needs to be updated for this reversal.
                // idx becomes n-idx-1
                for (int j = (i - 1); j > -1; --j) {
                    PairInt[] vEntry = edgeJoins[j];
                    for (int k = 0; k < 2; k++) {
                        PairInt vLoc = vEntry[k];
                        if (vLoc.getX() == loc1.getX()) {
                            int nV = edge1.getN();
                            vLoc.setY(nV - vLoc.getY() - 1);
                        }
                    }
                }
            }

            // --- append edge0 to edge1 ----
            edge1.addAll(edge0);
            
            // for earlier items in array edgeJoins
            // need to update all edge indexes and indexes within edge.
            
            // loc0 got moved to loc1
            
            // first, will only update the indexes within the edge for 
            // edgeIndex == loc0.getX()
            //     the new offset index = n0 + index
            for (int j = (i - 1); j > -1; --j) {
                PairInt[] vEntry = edgeJoins[j];
                for (int k = 0; k < 2; k++) {
                    PairInt vLoc = vEntry[k];
                    if (vLoc.getX() == loc0.getX()) {
                        vLoc.setY(vLoc.getY() + n1);
                    }
                }
            }
            
            // any edge index > loc0 gets reduced by one, but if == it gets replaced by loc1.getX()
            for (int j = (i - 1); j > -1; --j) {
                PairInt[] vEntry = edgeJoins[j];
                for (int k = 0; k < 2; k++) {
                    PairInt vLoc = vEntry[k];
                    if (vLoc.getX() > loc0.getX()) {
                        vLoc.setX(vLoc.getX() - 1);
                    } else if (vLoc.getX() == loc0.getX()) {
                        vLoc.setX(loc1.getX());
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
 
}
