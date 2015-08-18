package algorithms.imageProcessing;

import algorithms.misc.Misc;
import algorithms.util.PairIntArray;
import algorithms.util.PairInt;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.logging.Level;

/**
 * Edge extractor operates on an image that has already been reduced to 
 * single pixel width lines and extracts edges from it, attempting to make
 * the longest edges it can. 
 * 
 * 
 * @see AbstractEdgeExtractor

Edge extraction
    Local Methods:
        (1) DFS walk through connected pixel to form a sequence of pixels called
            an edge.
          
        (2) merge adjacent edges at the endpoints.

        (3) merge curves by closest points if the outlying points can be safely
            trimmed.
            
        (4) find edge endpoints which are separated from one another by a gap of
            one and fill in the gap while merging the edges.
        
        (5) remove edges shorter than a minimum length

 * @author nichole
 */
public class EdgeExtractor extends AbstractEdgeExtractor {
            
    /**
     * NOTE:  input should have a black (empty) background and edges should
     * have values > 125 counts.  Edges should also have width of 1 and no larger.
     * 
     * @param input 
     */
    public EdgeExtractor(GreyscaleImage input) {
        
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
    public EdgeExtractor(GreyscaleImage input, 
        final GreyscaleImage anEdgeGuideImage) {
        
        super(input, anEdgeGuideImage);
    }
    
    /**
     * this method is invoked by the super class after the DFS connect of
     * points in the image to make edges.  It merges the edges that have
     * adjacent endpoints.
     * It then looks for locations within the edges where it
     * can connect to other edges by trimming a small number of pixels between
     * them.
     * 
     * @return 
     */
    @Override
    public List<PairIntArray> findEdgesIntermediateSteps(List<PairIntArray> edges) {
        
        List<PairIntArray> output = mergeAdjacentEndPoints(edges);
        
        log.log(Level.INFO, "{0} edges after merge adjacent", 
            Integer.toString(output.size()));
        
        // This helps to merge edges (that is extracted curves) at adjacent 
        // points that resemble an intersection of the lines, but it's not 
        // necessarily useful if only interested in corners and not inflection
        // points because the curvature is determined correctly 
        // whether the curves are merged or not.
        
        output = connectClosestPointsIfCanTrim(output);
        
        log.info(output.size() + " edges after connect closest");
        
        output = fillInGaps(output);
        
        log.log(Level.INFO, "{0} edges after fill in gaps", 
            new Object[]{Integer.toString(output.size())});
         
        return output;
    }
 
    /**
     * merge edges adjacent end points of given edges.
     * For best results, make sure the edges were created from an image whose
     * lines were thinned to 1 pixel widths.
     * 
     * the runtime complexity is at best O(N).
    
     * @param edges
     * @return 
     */
    protected List<PairIntArray> mergeAdjacentEndPoints(
        List<PairIntArray> edges) {
        
        List<PairIntArray> output = new ArrayList<PairIntArray>();
        
        if (edges == null || edges.isEmpty()) {
            return output;
        }
        
        // 2 * O(N)
        Map<PairInt, Integer> endPointMap = createEndPointMap(edges);
        
        // initialize the current endpoint information:
        int currentEdgeIdx = 0;
        
        PairInt currentEndPoint = findStartingPoint(
            new PairInt(edges.get(currentEdgeIdx).getX(0), 
            edges.get(currentEdgeIdx).getY(0)), endPointMap,
            edges);
        
        currentEdgeIdx = endPointMap.get(currentEndPoint).intValue();
        
        currentEndPoint = getOppositeEndPointOfEdge(currentEndPoint, 
            edges.get(currentEdgeIdx));
                
        output.add(edges.get(currentEdgeIdx));
     
        endPointMap.remove(currentEndPoint);
        endPointMap.remove(getOppositeEndPointOfEdge(currentEndPoint, 
            edges.get(currentEdgeIdx)));
        
        int[] dxs = new int[]{-1, -1,  0,  1, 1, 1, 0, -1};
        int[] dys = new int[]{ 0, -1, -1, -1, 0, 1, 1,  1};
        
        List<Integer> foundEdgesIndexes = new ArrayList<Integer>();
        List<PairInt> foundEndPoints = new ArrayList<PairInt>();
        
        // 2 * O(N) * 8
        while (!endPointMap.isEmpty()) {
            
            int maxAdjEdgesIdx = -1;
            
            int maxAdjEdgesN = Integer.MIN_VALUE;
            
            PairInt maxAdjEdgesPoint = null;
            
            foundEdgesIndexes.clear();
            foundEndPoints.clear();         
            
            for (int nIdx = 0; nIdx < dxs.length; nIdx++) {
                int x = currentEndPoint.getX() + dxs[nIdx];
                int y = currentEndPoint.getY() + dys[nIdx];
                
                PairInt p = new PairInt(x, y);
                
                Integer eIdx = endPointMap.get(p);
                if (eIdx != null) {
                    foundEdgesIndexes.add(eIdx);
                    foundEndPoints.add(p);
                    int nPoints = edges.get(eIdx.intValue()).getN();
                    if (nPoints > maxAdjEdgesN) {
                        maxAdjEdgesN = nPoints;
                        maxAdjEdgesIdx = eIdx.intValue();
                        maxAdjEdgesPoint = p;
                    }
                }
            }
            
            if (!foundEndPoints.isEmpty()) {
                if (foundEndPoints.size() > 1) {
                    // this is a junction
                    
                    // add entry too junctionMap
                    Set<PairInt> values = new HashSet<PairInt>(foundEndPoints);
                    
                    PairIntArray lastInOutput = output.get(output.size() - 1);
                    int n = lastInOutput.getN();
                    PairInt prevPoint = new PairInt(lastInOutput.getX(n - 2),
                        lastInOutput.getY(n - 2));
                    
                    values.add(prevPoint);
                    
                    currentEdgeIdx = maxAdjEdgesIdx;
                    
                    currentEndPoint = getOppositeEndPointOfEdge(
                        maxAdjEdgesPoint, edges.get(currentEdgeIdx));
                    
                } else {
                    
                    currentEdgeIdx = foundEdgesIndexes.get(0);
                    
                    currentEndPoint = getOppositeEndPointOfEdge(
                        foundEndPoints.get(0), edges.get(currentEdgeIdx));                      
                }
              
            } else {
                
                output.add(new PairIntArray());
                
                // since the last edge has ended, pick any remaining in
                // endPointMap as the start of the next.
                // then invoke a method to follow the endpoint back to find
                // the true start
                
                Entry<PairInt, Integer> tmp = 
                    endPointMap.entrySet().iterator().next();
                
                currentEndPoint = findStartingPoint(tmp.getKey(), endPointMap, 
                    edges);
                
                currentEdgeIdx = endPointMap.get(currentEndPoint).intValue();
            }
                        
            appendToOutput(output, edges.get(currentEdgeIdx));
     
            PairInt oppositeEndPoint = getOppositeEndPointOfEdge(currentEndPoint, 
                edges.get(currentEdgeIdx));
            
            Integer v0 = endPointMap.remove(currentEndPoint);
            Integer v1 = endPointMap.remove(oppositeEndPoint);
            assert(v0 != null);
            assert(v1 != null);
        }
      
        return output;
    }
    
    protected Map<PairInt, Integer> createEndPointMap(List<PairIntArray>
        edges) {
        
        Map<PairInt, Integer> endPointMap = new HashMap<PairInt, Integer>();
        
        for (int idx = 0; idx < edges.size(); idx++) {
            
            PairIntArray edge = edges.get(idx);
            
            int n = edge.getN();
            
            PairInt start = new PairInt(edge.getX(0), edge.getY(0));
            
            PairInt end = new PairInt(edge.getX(n - 1), edge.getY(n - 1));
            
            endPointMap.put(start, Integer.valueOf(idx));
            
            endPointMap.put(end, Integer.valueOf(idx));
            
        }
        
        return endPointMap;
    }
    
    /**
     * fill in gaps of '1' pixel
     * 
     * runtime complexity:
     *   2 * O(N_edges^2)
     * 
     * @param edges
     * @return 
     */
    protected List<PairIntArray> fillInGaps(List<PairIntArray>
        edges) {
            
        /*
        similar to the mergeAdjacentEndPoints:  
        compare end of uEdge to beginning of all others
        reverse uEdge and repeat
        revert reverse for next start
        */
      
        boolean[] removed = new boolean[edges.size()];
        
        List<PairIntArray> output = new ArrayList<PairIntArray>();
                
        for (int i = 0; i < edges.size(); i++) {
            
            if (removed[i]) {
                continue;
            }
            
            PairIntArray uEdge = edges.get(i);
            
            // an extra iteration for reversing uEdge
            for (int r = 0; r < 2; r++) {
                 
                // compare bottom of uEdge to top of vEdge

                for (int j = 0; j < edges.size(); j++) {

                    if (i == j) {
                        continue;
                    }
                    if (removed[j]) {
                        continue;
                    }

                    PairIntArray vEdge = edges.get(j);

                    // recalculate in case u has grown
                    int uX = uEdge.getX(uEdge.getN() - 1);
                    int uY = uEdge.getY(uEdge.getN() - 1);

                    int vX = vEdge.getX(0);
                    int vY = vEdge.getY(0);

                    int diffX = uX - vX;
                    if (diffX < 0) {
                        diffX *= -1;
                    }

                    if (diffX > 2) {
                        continue;
                    }

                    int diffY = uY - vY;
                    if (diffY < 0) {
                        diffY *= -1;
                    }

                    if (diffY > 2) {
                        continue;
                    }
                    
                    int gapX = (uX + vX)/2;
                    int gapY = (uY + vY)/2;
                                        
                    img.setValue(gapX, gapY, 255);
                    
                    uEdge.add(gapX, gapY);
                    
                    for (int k = 0; k < vEdge.getN(); k++) {
                        uEdge.add(vEdge.getX(k), vEdge.getY(k));
                    }

                    removed[j] = true;

                    // have to restart the j iteration to re-compare terms
                    j = -1;
                }
                
                if (r == 0) {
                    // just finished forward, start revers 
                    uEdge.reverse();
                } else if (r == 1) {
                    // revert the array back to other direction
                    uEdge.reverse();
                }
            }
        }
        
        for (int i = 0; i < edges.size(); i++) {
            if (!removed[i]) {
                output.add(edges.get(i));
            }
        }
        
        return output;
    }
    
    /**
     * find the closest point between the curve0 and curve1 and return the
     * points in curve0XY and curve1XY along with the method return value
     * which is the separation.
     * 
     * runtime complexity:
     *    O(N_edge1 x N_edge2)
     * 
     * @param curve0 
     * @param curve1
     * @param curve0Idx output variable to hold index to the (x, y) of the point 
     * in curve0 which is closest to curve0
     * @param curve1Idx output variable to hold index to the (x, y) of the point 
     * in curve0 which is closest to curve1
     * @return the separation between the closest pair of points in curve0 and 
     *   curve1
     */
    protected double findClosestPair(PairIntArray curve0, PairIntArray curve1, 
        int[] curve0Idx, int[] curve1Idx) {
        
        int idx0 = -1;
        int idx1 = -1;
        double min = Double.MAX_VALUE;
        
        for (int i = 0; i < curve0.getN(); i++) {
            int x0 = curve0.getX(i);
            int y0 = curve0.getY(i);
            for (int j = 0; j < curve1.getN(); j++) {
                int x1 = curve1.getX(j);
                int y1 = curve1.getY(j);
                int dx = x1 - x0;
                int dy = y1 - y0;
                double d = dx*dx + dy*dy;
                if (d < min) {
                    min = d;
                    idx0 = i;
                    idx1 = j;
                }
            }
        }
        
        curve0Idx[0] = idx0;
        
        curve1Idx[0] = idx1;
        
        return Math.sqrt(min);
    }
    
    /**
     * connect the closest points in edges if trimming the outliers does not
     * remove too many points nor add a discontinuity in either edge.
     * 
     * Runtime complexity:
     *      O(N_edge x ~N_edge x (N_edge to N_edge^2))
     * 
     * @param edges
     * @return 
     */
    protected List<PairIntArray> connectClosestPointsIfCanTrim(
        List<PairIntArray> edges) {
             
        double sqrtTwo = Math.sqrt(2) + 0.01;
        //double gapOfOne = 2*Math.sqrt(2) + 0.01;
        
        int[] edge0Idx = new int[1];
        int[] edge1Idx = new int[1];
        
        boolean[] removed = new boolean[edges.size()];
        
        List<PairIntArray> output = new ArrayList<PairIntArray>();
        
        for (int i = 0; i < edges.size(); i++) {
            
            if (removed[i]) {
                continue;
            }
            
            PairIntArray edge0 = edges.get(i);
            
            for (int j = (i + 1); j < edges.size(); j++) {
                
                if (removed[j]) {
                    continue;
                }
                
                PairIntArray edge1 = edges.get(j);
                
                // RT: O(N_edge0 x N_edge1)
                
                double sep = findClosestPair(edge0, edge1, edge0Idx, edge1Idx);

                if (sep < sqrtTwo) {
                    
                    // RT: O(N_edge0) or O(N_edge1)
                    
                    // do not merge them if the points are not near the
                    // ends of the points sets.
                    float closestFrac0 = (float)edge0Idx[0]/(float)edge0.getN();
                    
                    float closestFrac1 = (float)edge1Idx[0]/(float)edge1.getN();
                    
                    if (((closestFrac0 > 0.07) && (closestFrac0 < 0.93))) {
                        continue;
                    }
                    if (((closestFrac1 > 0.07) && (closestFrac1 < 0.93))) {
                        continue;
                    }
                    
                    boolean closest0IsNearTop = ((float)(edge0Idx[0]/
                        (edge0.getN() - edge0Idx[0]))) <= 0.5;
                    
                    if (closest0IsNearTop) {
                        // if we trim the top, is remaining bottom connected?
                        boolean isConnected = isRangeConnected(edge0, 
                            edge0Idx[0], edge0.getN() - 1);
                        if (!isConnected) {
                            continue;
                        }
                    } else {
                        // if we trim the bottom, is remaining top connected?
                        boolean isConnected = isRangeConnected(edge0, 0, 
                            edge0Idx[0]);                      
                        if (!isConnected) {
                            continue;
                        }
                    }
                                            
                    boolean closest1IsNearTop = (
                        ((float)edge1Idx[0]/
                        (float)(edge1.getN() - edge1Idx[0]))) <= 0.5;
                    
                    if (closest1IsNearTop) {
                        // if we trim the top, is remaining bottom connected?
                        boolean isConnected = isRangeConnected(edge1, 
                            edge1Idx[0], edge1.getN() - 1);
                        if (!isConnected) {
                            continue;
                        }
                    } else {
                        // if we trim the bottom, is remaining top connected?
                        boolean isConnected = isRangeConnected(edge1, 0, 
                            edge1Idx[0]);                      
                        if (!isConnected) {
                            continue;
                        }
                    }
                    
                    // if here, then can trim outside the closest points in the
                    // edges and merge the edges into edge0, and remove edge1
                    
                    if (closest0IsNearTop) {
                        if (edge0Idx[0] > 0) {
                            edge0.removeRange(0, edge0Idx[0] - 1);
                        }
                    } else {
                        if (edge0Idx[0] < (edge0.getN() - 1)) {
                            edge0.removeRange(edge0Idx[0] + 1, edge0.getN() - 1);
                        }
                    }
                    
                    //TODO:  could remove this step and adjust the add, but 
                    //       easier maintainence this way
                    if (closest1IsNearTop) {
                        if (edge1Idx[0] > 0) {
                            edge1.removeRange(0, edge1Idx[0] - 1);
                        }
                    } else {
                        if (edge1Idx[0] < (edge1.getN() - 1)) {
                            edge1.removeRange(edge1Idx[0] + 1, edge1.getN() - 1);
                        }
                    }
                    
                    if (closest0IsNearTop) {
                        // insert edge1 at top of edge0
                        edge0.insertSpaceAtTopOfArrays(edge1.getN());
                        
                        // if edge1 closest is at bottom of it's edge, just add,
                        // else reverse then add
                        if (closest1IsNearTop) {
                            edge1.reverse();
                        }

                        for (int k = 0; k < edge1.getN(); k++) {
                            edge0.set(k, edge1.getX(k), edge1.getY(k));
                        }
                    } else {
                        // append edge1 to bottom of edge0
                        
                        if (!closest1IsNearTop) {
                            edge1.reverse();
                        }

                        for (int k = 0; k < edge1.getN(); k++) {
                            edge0.add(edge1.getX(k), edge1.getY(k));
                        }
                    }
                  
                    removed[j] = true;
                    
                    // have to restart the j iteration to re-compare terms
                    j = i;
                }
            }
            output.add(edge0);
        }
        
        return output;
    }
    
    /**
     * check that points within index idxLo and idxHi, inclusive, are 
     * consecutively within 1 pixel of adjacent indexes.
     * @param edge
     * @param idxLo
     * @param idxHi
     * @return 
     */
    protected boolean isRangeConnected(PairIntArray edge, int idxLo, int idxHi) {
        
        for (int i = (idxLo + 1); i <= idxHi; i++) {
            int x0 = edge.getX(i - 1);            
            int x1 = edge.getX(i);
            int diffX = x0 - x1;
            if (diffX < 0) {
                diffX *= -1;
            }
            if (diffX > 1) {
                return false;
            }
            
            int y0 = edge.getY(i - 1);
            int y1 = edge.getY(i);
            int diffY = y0 - y1;
            if (diffY < 0) {
                diffY *= -1;
            }
            if (diffY > 1) {
                return false;
            }
        }
        
        return true;
    }

    /**
     * Search endMap "backwards" to find the earliest match to the startPoint.
     * When a junction is found, choose the longest edge.
     * <pre>
     *           -- find 
     *             \
     *              \
     *               \           start
     *                \find     *point
     *             |---/-------||------|
     *       discard
     * </pre>
     * @param startPoint
     * @param endPointMap
     * @param edges
     * @return 
     */
    protected PairInt findStartingPoint(PairInt startPoint, Map<PairInt, Integer> 
        endPointMap, List<PairIntArray> edges) {
        
        if (startPoint == null) {
            throw new IllegalArgumentException("startPoint cannot be null");
        }
        if (endPointMap == null) {
            throw new IllegalArgumentException("endPointMap cannot be null");
        }
        if (edges == null) {
            throw new IllegalArgumentException("edges cannot be null");
        }
        
        int[] dxs = new int[]{-1, -1,  0,  1, 1, 1, 0, -1};
        int[] dys = new int[]{ 0, -1, -1, -1, 0, 1, 1,  1};
        
        PairInt originalStartPoint = startPoint;
        
        PairInt lastStartPoint = startPoint;
        
        PairInt currentStartPoint = startPoint;
                
        Set<PairInt> visited = new HashSet<PairInt>();
        visited.add(currentStartPoint);
        
        int nIter = 0;
        
        while (currentStartPoint != null) {
                        
            if ((nIter > 0) && (lastStartPoint.equals(currentStartPoint))) {
                break;
            }
            
            int maxN = Integer.MIN_VALUE;
            Integer maxNIndex = null;            
            PairInt maxNPoint = null;
            
            for (int nIdx = 0; nIdx < dxs.length; nIdx++) {
                
                int x = currentStartPoint.getX() + dxs[nIdx];
                int y = currentStartPoint.getY() + dys[nIdx];
                
                PairInt p = new PairInt(x, y);
                
                Integer eIndex = endPointMap.get(p);
                
                if (eIndex != null) {
                    
                    PairIntArray pai = edges.get(eIndex.intValue());
                    
                    int n = pai.getN();
                    
                    if (n > maxN) {
                        
                        // the opposite end of this edge becomes
                        // our next potential currentStartPoint
                        PairInt reversed = getOppositeEndPointOfEdge(p, pai);
                        
                        if (visited.contains(reversed)) {
                            // if this is a closed curve, return the original
                            if (reversed.equals(originalStartPoint)) {
                                return originalStartPoint;
                            }
                            continue;
                        }
    
                        maxN = n;
                        maxNIndex = eIndex;
                        maxNPoint = reversed;
                    }
                }
            }
            
            lastStartPoint = currentStartPoint;
            
            if (maxNIndex == null) {
                currentStartPoint = null;
            } else {
                
                currentStartPoint = maxNPoint;
        
                visited.add(currentStartPoint);
                
                // if this is a closed curve, return the original value
                if (currentStartPoint.equals(originalStartPoint)) {
                    return originalStartPoint;
                }
            }
            
            nIter++;
        }
                        
        return lastStartPoint;
    }

    /**
     * append edge to output, reversing the edge if needed to minimize the 
     * distance between the last point in output and the first point 
     * appended from edge.
     * 
     * @param output
     * @param edge 
     */
    protected void appendToOutput(List<PairIntArray> output, PairIntArray edge) {
        
        PairIntArray lastOutputEdge = output.get(output.size() - 1);
        
        if (lastOutputEdge.getN() == 0) {
            lastOutputEdge.addAll(edge);
            return;
        }
        
        int lastX = lastOutputEdge.getX(lastOutputEdge.getN() - 1);
        int lastY = lastOutputEdge.getY(lastOutputEdge.getN() - 1);
        
        int x0 = edge.getX(0);
        int y0 = edge.getY(0);
        int diffX0 = x0 - lastX;
        int diffY0 = y0 - lastY;
        long dist0Sq = (diffX0 * diffX0) + (diffY0 * diffY0);
        
        int xn = edge.getX(edge.getN() - 1);
        int yn = edge.getY(edge.getN() - 1);
        
        int diffXn = xn - lastX;
        int diffYn = yn - lastY;
        long distnSq = (diffXn * diffXn) + (diffYn * diffYn);
        
        if (dist0Sq > distnSq) {
            // reverse
            PairIntArray rev = edge.copy();
            rev.reverse();
            lastOutputEdge.addAll(rev);
        } else {
            lastOutputEdge.addAll(edge);
        }
    }
    
}
