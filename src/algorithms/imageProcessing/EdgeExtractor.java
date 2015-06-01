package algorithms.imageProcessing;

import algorithms.util.PairIntArray;
import algorithms.util.PairInt;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Stack;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Edge extractor operates on an image that has already been reduced to holding
 * single pixel width lines where the original image was.  
 * The EdgeExtractor additionally accepts an argument that is an image
 * that is used as guidance to correct the extracted edges.  The guidance
 * image, for example, is expected to be the combined gradient X and Y image
 * that was used in forming the main input image.  The guidance image is used
 * to make minor corrections to the coordinates of edge pixels if the changes
 * do not disconnect the edge.  The guidance image helps corrects errors in
 * the input image due to a line thinner that doesn't make thinning decisions
 * on a slightly larger scale, for example.
 * 
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
        
        (6) if an edge guide image was provided, make adjustments to edge 
            towards highest intensity pixels in the edge guide image as long
            as the adjustment doesn't create a gap in the edge.
            This stage also reduces any redundant pixels that may be present
            in the line.

 * @author nichole
 */
public class EdgeExtractor {
    
    private final GreyscaleImage img;
    
    private GreyscaleImage edgeGuideImage = null;
    
    private long numberOfPixelsAboveThreshold = 0;
    
    private Logger log = Logger.getLogger(this.getClass().getName());
            
    /**
     * if the image is smaller than 100 on a side, this will be lowered to 5
     */
    private int edgeSizeLowerLimit = 15;
    
    private boolean repeatConnectAndTrim = false;
    
    /**
     * NOTE:  input should have a black (empty) background and edges should
     * have values > 125 counts.  Edges should also have width of 1 and no larger.
     * 
     * @param input 
     */
    public EdgeExtractor(GreyscaleImage input) {
        
        img = input;
        
        if (img.getWidth() < 100 || img.getHeight() < 100) {
            edgeSizeLowerLimit = 5;
        }
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
        
        img = input;
        
        edgeGuideImage = anEdgeGuideImage;
        
        if (img.getWidth() < 100 || img.getHeight() < 100) {
            edgeSizeLowerLimit = 5;
        }
    }
    
    public GreyscaleImage getImage() {
        return img;
    }
    
    public void overrideEdgeSizeLowerLimit(int length) {
        edgeSizeLowerLimit = length;
    }
    
    public void useRepeatConnectAndTrim() {
        repeatConnectAndTrim = true;
    }
    
    /**
     * find the edges and return as a list of points.  The method uses a
     * DFS search through all points in the image with values > 0 to link
     * adjacent sequential points into edges.
     * The method then uses method mergeAdjacentEndPoints.
     * Note that the later 2 methods are not needed if the edges will be used
     * in a corner detector only, but if the edges are to be used to
     * find inflection points in scale space maps, those methods help to
     * provide more complete shapes and better matches between the same object
     * in another scale space map.
     * 
     * @return 
     */
    public List<PairIntArray> findEdges() {
        
        int w = img.getWidth();
        int h = img.getHeight();
        
        int[] dxs = new int[]{-1, -1,  0,  1, 1, 1, 0, -1};
        int[] dys = new int[]{ 0, -1, -1, -1, 0, 1, 1,  1};
        
        // DFS search for sequential neighbors.
        
        Stack<PairInt> stack = new Stack<PairInt>();
      
        int thresh0 = 1;
        
        for (int i = 0; i < w; i++) {
            for (int j = 0; j < h; j++) {
                int v = img.getValue(i, j);
                if (v >= thresh0) {
                    stack.add(new PairInt(i, j));
                }
            }
        }
        
        numberOfPixelsAboveThreshold = stack.size();
        
        log.log(Level.FINE, "Number of pixels that exceed threshhold={0}", 
            Long.toString(numberOfPixelsAboveThreshold));
        
        List<PairIntArray> output = new ArrayList<PairIntArray>();
        int[] uNodeEdgeIdx = new int[img.getWidth() * img.getHeight()];
        Arrays.fill(uNodeEdgeIdx, -1);
           
        while (!stack.isEmpty()) {
            
            PairInt uNode = stack.pop();
            
            int uX = uNode.getX();
            int uY = uNode.getY();
            int uIdx = img.getIndex(uX, uY);
                        
            for (int nIdx = 0; nIdx < dxs.length; nIdx++) {
                
                int vX = dxs[nIdx] + uX;
                int vY = dys[nIdx] + uY;

                if ((vX < 0) || (vX > (w - 1)) || (vY < 0) || (vY > (h - 1))) {
                    continue;
                }

                int vIdx = img.getIndex(vX, vY);
                
                if (uNodeEdgeIdx[vIdx] != -1 || (uIdx == vIdx)) {
                    continue;
                }

                if (img.getValue(vX, vY) < thresh0) {
                    continue;
                }
                   
                processNeighbor(uX, uY, uIdx, vX, vY, vIdx, uNodeEdgeIdx, output);

                //TODO: do we only want 1 neighbor from the 9 as a continuation?
                // yes for now, but this requires edges be only 1 pixel wide
                // inserting back at the top of the stack assures that the 
                // search continues next from an associated point
                stack.add(new PairInt(vX, vY));
            }
        }
        
        log.fine(output.size() + " edges after DFS");
        
        int nIterMax = 100;
        int n, sz, lastSize;
        
        // count the number of points in edges
        long sum = countPixelsInEdges(output);
        
        log.log(Level.FINE, 
            "==> {0} pixels are in edges out of {1} pixels > threshhold", 
            new Object[]{Long.toString(sum), 
                Long.toString(numberOfPixelsAboveThreshold)});
        
        log.log(Level.FINE, "there are {0} edges", 
            Integer.toString(output.size()));
        
        output = mergeAdjacentEndPoints(output);
        
        log.log(Level.FINE, "{0} edges after merge adjacent", 
            Integer.toString(output.size()));
        
        // This helps to merge edges (that is extracted curves) at adjacent 
        // points that resemble an intersection of the lines, but it's not 
        // necessarily useful if only interested in corners and not inflection
        // points because the curvature is determined correctly 
        // whether the curves are merged or not.
        
        output = connectClosestPointsIfCanTrim(output);
        
        log.fine(output.size() + " edges after connect closest");
        
        
        output = fillInGaps(output);
        
        log.log(Level.FINE, "{0} edges after fill in gaps", 
            new Object[]{Integer.toString(output.size())});
         
                
        //TODO:  this may need to change
        removeEdgesShorterThan(output, edgeSizeLowerLimit);
        sz = output.size();
        
        
        log.log(Level.FINE, "{0} edges after removing those shorter", 
            new Object[]{Integer.toString(sz)});

        
        long sum2 = countPixelsInEdges(output);
               
        log.log(Level.FINE, 
            "==> {0}) pixels are in edges out of {1} pixels > threshhold", 
            new Object[]{Long.toString(sum2), 
                Long.toString(numberOfPixelsAboveThreshold)});
       
        //pruneSpurs(output);
        
        // necessary only if the LineThinner in CannyEdgeFilter was ErosionFilter
        //adjustEdgesTowardsBrightPixels(output);
        
        if (repeatConnectAndTrim) {
            
            output = connectClosestPointsIfCanTrim(output);
        }
    
        return output;
    }
    
    private long countPixelsInEdges(List<PairIntArray> edges) {
        long sum = 0;
        for (PairIntArray edge : edges) {
            sum += edge.getN();
        }
        return sum;
    }
   
    /**
     * merge edges adjacent end points of curves
     * 
     * runtime complexity:
     *   2 * O(N_edges^2)
     * 
     * @param edges
     * @return 
     */
    protected List<PairIntArray> mergeAdjacentEndPoints(List<PairIntArray> 
        edges) {
     
        /*
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

                    if (diffX > 1) {
                        continue;
                    }

                    int diffY = uY - vY;
                    if (diffY < 0) {
                        diffY *= -1;
                    }

                    if (diffY > 1) {
                        continue;
                    }

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
     * merge edges adjacent end points of given edges. 
     * 
     * the runtime complexity is at best O(N)
    
     * @param edges
     * @return 
     */
    protected List<PairIntArray> mergeAdjacentEndPoints2(
        List<PairIntArray> edges) {
        
        List<PairIntArray> output = new ArrayList<PairIntArray>();
        
        if (edges == null || edges.isEmpty()) {
            return output;
        }
        
        // key = center of junction pixel coordinates
        // value = set of adjacent pixels when there are more than the preceding
        //         and next.
        Map<PairInt, Set<PairInt> > junctionMap = new HashMap<PairInt, Set<PairInt>>();
        
        // key = pixel coordinates of all pixels involved in junctions
        // value = index in edges (NOTE: at end of method, these are converted to 
        //                        indexes of output list).
        Map<PairInt, Integer> junctionLocationMap = new HashMap<PairInt, Integer>();
        
        // key = edges index
        // value = indexes of the output edge that the pixel is in.
        Map<Integer, Integer> junctionEdgesToOutputIndexesMap = 
            new HashMap<Integer, Integer>();
        
        Map<PairInt, Integer> endPointMap = createEndPointMap(edges);
        
        // initialize the current endpoint information:
        int currentEdgeIdx = 0;
        
        
        //TODO: use findStartingIndex(currentStartPoint, endPointMap) here
        
        PairInt currentEndPoint = new PairInt(
            edges.get(currentEdgeIdx).getX(edges.get(currentEdgeIdx).getN() - 1), 
            edges.get(currentEdgeIdx).getY(edges.get(currentEdgeIdx).getN() - 1));
        
        output.add(edges.get(currentEdgeIdx));
     
        endPointMap.remove(currentEndPoint);
        endPointMap.remove(new PairInt(edges.get(currentEdgeIdx).getX(0),
            edges.get(currentEdgeIdx).getY(0)));
        
        
        int[] dxs = new int[]{-1, -1,  0,  1, 1, 1, 0, -1};
        int[] dys = new int[]{ 0, -1, -1, -1, 0, 1, 1,  1};
        
        List<Integer> foundEdgesIndexes = new ArrayList<Integer>();
        List<PairInt> foundEndPoints = new ArrayList<PairInt>();
        
        while (!endPointMap.isEmpty()) {
            
            int maxAdjEdgesIdx = -1;
            
            int maxAdjEdgesN = Integer.MIN_VALUE;
            
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
                    }
                } else if (junctionMap.containsKey(p)) {
                    eIdx = junctionLocationMap.get(p);
                    assert(eIdx != null);
                    foundEdgesIndexes.add(eIdx);
                    foundEndPoints.add(p);
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
                    
                    junctionMap.put(currentEndPoint, values);
                    
                    // add an entry to junctionLocationMap for currentEndPoint
                    junctionLocationMap.put(currentEndPoint, Integer.valueOf(currentEdgeIdx));
                    
                    // add an entry to junctionLocationMap for prevPoint
                    junctionLocationMap.put(prevPoint, Integer.valueOf(currentEdgeIdx));
                    
                    // add an entry to junctionEdgesToOutputIndexesMap for currentEndPoint
                    junctionEdgesToOutputIndexesMap.put(
                        Integer.valueOf(currentEdgeIdx), 
                        Integer.valueOf(output.size() - 1));
                    
                    // add entries to junctionLocationMap for items in
                    // foundEdgesIndexes, foundEndPoints
                    for (int i = 0; i < foundEdgesIndexes.size(); i++) {
                        junctionLocationMap.put(foundEndPoints.get(i), 
                            foundEdgesIndexes.get(i));
                    }
                      
                    currentEdgeIdx = maxAdjEdgesIdx;
                    
                } else {
                    
                    currentEdgeIdx = foundEdgesIndexes.get(0);                    
                }
                
                currentEndPoint = new PairInt(
                    edges.get(currentEdgeIdx).getX(edges.get(currentEdgeIdx).getN() - 1), 
                    edges.get(currentEdgeIdx).getY(edges.get(currentEdgeIdx).getN() - 1));
                
            } else {
                if (junctionMap.containsKey(currentEndPoint)) {
                    junctionEdgesToOutputIndexesMap.put(Integer.valueOf(currentEdgeIdx),
                        Integer.valueOf(output.size() - 1));
                }
                
                output.add(new PairIntArray());
                
                // since the last edge has ended, pick any remaining in
                // endPointMap as the start of the next.
                // then invoke a method to follow the endpoint back to find
                // the true start
                
                currentEdgeIdx = endPointMap.entrySet().iterator().next().getValue();
                
                PairInt currentStartPoint = new PairInt(
                    edges.get(currentEdgeIdx).getX(0), 
                    edges.get(currentEdgeIdx).getY(0));
               
                //TODO: currentEdgeIdx = findStartingIndex(currentStartPoint, endPointMap);
                    
                currentEndPoint = new PairInt(
                    edges.get(currentEdgeIdx).getX(edges.get(currentEdgeIdx).getN() - 1), 
                    edges.get(currentEdgeIdx).getY(edges.get(currentEdgeIdx).getN() - 1));
                
            }
            
            output.add(edges.get(currentEdgeIdx));
     
            endPointMap.remove(currentEndPoint);
            endPointMap.remove(new PairInt(edges.get(currentEdgeIdx).getX(0),
                edges.get(currentEdgeIdx).getY(0)));
        }
       
        // convert the values in junctionMap from edges indexes to output
        
        // store junctionMap and junctionLocationsMap as member variables
        
        /*
        TODO:  improve this section.  Could be done in O(N) and revised to add
        the missing junction information that subsequent operations need.
                
        (1) DONE 
           create method createEndPointMap(List<PairIntArray>) : Map<PairInt, Integer>
           need data structure to hold the searchable information of edges.  
           need to be able to search by start endpoint (x,y) and by end 
           endpoint (x,y).  the structure needs to retain reference to the item 
           in edges which it describes.
           searches are exact, so can use a HashSet for O(1) inserts and
           O(1) removals.  would like a key by start endpoint, and a key by end
           endpoint, so would use a hashmap.
           ==> Data structure needed is a HashMap:
               HashMap endPointMap
               key = PairInt of start (x,y) or end (x,y)
               value = edges index of this endpoint's edge
           ==> runtime complexity is 2 * O(N) for creating the HashMap
        (2) DONE
            retain reference to the first item in edges
            make local variables for this:
            ==> currentEdgeIdx which is the index to the edge in edges.
            ==> currentEndPoint which is a PairInt holding the last (x, y).
REVISE:  use findStartingIndex(currentStartPoint, endPointMap) to
search back until find start of that edge
        (3) DONE
            initialize local variables
            ==> create an output List<PairIntArray> 
            -- add currentEdgeIdx edge from edges to the output list
            -- remove currentEndPoint from endPointMap        
        (4) DONE
            need structures to hold junction information when more than one
            endpoint matches.  this is to be used later.  
            -- it needs to refer to the items in the output list to avoid 
            re-searching later, and it needs to retain the (x,y)'s of the junction.
            ==> Data structures needed are three HashMaps:
               HashMap = junctionMap
               key = PairInt of junction (x,y)
               value = set of PairInt coordinates that are adjacent pixels
        
               HashMap = junctionLocationMap
               key = PairInt of junction (x,y)
               value = index of the edges list of the PairIntArray holding
                       the point
                       (NOTE: before last step in method, this will be 
                       converted to index of the *output* list.)
        
               HashMap = junctionEdgesToOutputIndexesMap
               key = edges index
               value = output index
        
               where junctionMap and junctionLocationMap is added to for 
               each currentEndPoint that has more than one adjacent pixel.
               junctionEdgesToOutputIndexesMap holds lookup information for the
               currentEndPoint and for all of its adjacent pixels as they
               are added to the output list.
               junctionLocationMap holding not only the center junctions, 
               but also the adjacent points to make them findable later.
        (5) IN PROGRESS 
           search for currentEndPoint neighbors:
           -- for the 8 neighbors of currentEndPoint: 
              -- search endPointMap
              -- search junctionMap 
                 (can use junctionLocationMap if contains is true to get the
                  needed edges index)
              -- keep results in 2 parallel lists:
                 foundEdgesIndexes
                 foundEndPoints
              -- from found items from endPointMap keep the one with the
                 most edge points as maxAdjEdgesIdx and maxAdjEdgesN
           -- if any items are in foundEndPoints:
              -- if there are more than one items in foundEndPoints:
                 -- add an entry to junctionMap:
                    key = currentEndPoint PairInt(x,y)                        
                    value = foundEdgePoints.  
                            PLUS add the next to the last point in the last item 
                            of output list.
                    If the key already exists, add to the existing values.
                 -- add an entry to junctionLocationMap:
                    key = currentEndPoint PairInt(x,y)
                    value = currentEdgeIdx
                 -- add an entry to junctionLocationMap:
                    key = next to last point in output.get(output.size() - 1)
                    value = currentEdgeIdx
                 -- add to junctionEdgesToOutputIndexesMap:
                    key = currentEdgeIdx
                    value = output.size() - 1
                    Note that this makes currentEndPoint lookup of output
                    index possible and also the next to last point.
                 -- for each item in foundEdgesIndexes, foundEndPoints,
                    add and entry to junctionLocationMap
              -- set reference PairIntArray to maxPairIntArray,
                 that is currentEdgeIdx and currentEndPoint
           -- else if not found:
             -- IF currentEndPoint is found in junctionMap
                and an entry to junctionEdgesToOutputIndexesMap
                key=currentEdgeIdx
                value= size of output list - 1 (this is the index where will be added)
             -- start a new PairIntArray in output list
             -- choose a new reference PairIntArray from iter of endPointMap,
                that is currentEdgeIdx and currentEndPoint
                NOTE, this may not be the start of its true connected edge.
                can think of several ways to accomodate this.
                   a. one could create a method to search for adjacent pixels
                      to the start point, and if more than one is found, choose
                      the longest edge, and keep searching "backwards" until
                      either find no more adjacent edges or until have reached
                      the start again (as would happen for closed curves.
                      then the return of this method is used
                      as currentEdgeIdx and currentEndPoint.
                      the runtime complexity for this method is potentially as long as the
                      number of edges that compose it is so is probably the
                      fastest of the possibilities
                   b. could just use the point that is possibly not the
                      start point and make a wrapper method to invoke
                      this method with the previous output as input until
                      the output size is the same as the input size.
                      the runtime of this one could be large...
REVISE: will choose a. and implement it as
                      findStartingIndex(PairInt p, Map<PairInt, Integer> endPointMap) : int
        
           -- THEN (after else/if)
              -- append the item at edges.get(currentEdgeIdx) to the current last item in
                 output (if last is empty, can just replace it instead of append)
              -- remove currentEndPoint and it's edge's start endpoint from endPointMap.
            ==> 16*O(N)
        (6) convert the values in junctionMap from edges indexes to output
            list indexes using junctionEdgesToOutputIndexesMap.
            ==> O(N)
        (7) store junctionMap and junctionLocationsMap as member variables.
            
        
        ====> runtime complexity is O(N)
        
        
        Tests:
        
           Tests for the internal methods that mergeAdjacentEndPoints uses:
               -- for createEndPointMap(List<PairIntArray>) : Map<PairInt, Integer>
                  -- 5 or so edges, with 3 meeting at a junction and the other 
                     two connected also:
                     
                        0      4   / 3
                     |-----||-----|
                                   \
                                  1 \     2
                                     \|-------
                     
                     -- assert that the returned HashMap is size 10
                     -- assert that each point in the sketch has the expected
                        indexes in the returned HashMap
                     
                     Make several tests of the above data w/ different combinations
                     of reversed start and end points in edges 1 thru 4
        
            Tests for entire method mergeAdjacentEndPoints:
               
               -- given 5 or so edges which do not have junctions and should result
                  in a final closed curve, that is one edge.
                      
                        0      4      1      3       2
                     |-----||-----||-----||-----||------|
                                     
                  -- assert that junctionsMap and junctionLocationsMap are empty
                  -- assert the ordered points in the output final edge.

                  Make several tests of the above data w/ different combinations
                     of reversed start and end points, excepting the first
                     edge which should always be the one with the smallest
                     column and smallest row as a start endpoint.
        
               -- 5 or so edges, with 3 meeting at a junction and the other 
                     two connected also:
                     
                        0      4   / 3
                     |-----||-----|
                                   \
                                  1 \     2
                                     \|-------
        
                  -- assert that junctionMap.size() is 1 and that 
                     junctionLocationMap.size() is 4  
                  -- assert the values in junctionMap and junctionLocationMap
                  -- assert that output list contains 2 edges of expected 
                     ordered points
        
               -- 5 or so edges, with 3 meeting at a junction.
                  the final output should be 3 edges, one of which contains
                  the two non-junction containing edges plus the longest of
                  the junction containing edges.
                      
                        0         4   / 3
                     |-----|   |-----|
                                      \
                                     1 \           2
                                        \      |-------
        
                  -- assert that junctionMap.size() is 1 and that 
                     junctionLocationMap.size() is 4 
                  -- assert the values in junctionMap and junctionLocationMap
                  -- assert that output list contains 4 edges of expected 
                     points

                  Make several tests of the above data w/ different combinations
                     of reversed start and end points, excepting the first
                     edge which should always be the one with the smallest
                     column and smallest row as a start endpoint.

               -- assert that given an empty list of edges returns an empty output list
        
        */
        
        boolean[] removed = new boolean[edges.size()];

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

                    if (diffX > 1) {
                        continue;
                    }

                    int diffY = uY - vY;
                    if (diffY < 0) {
                        diffY *= -1;
                    }

                    if (diffY > 1) {
                        continue;
                    }

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

    protected void removeEdgesShorterThan(List<PairIntArray> output, 
        int minNumberOfPixelsInEdge) {
        
        for (int i = (output.size() - 1); i > -1; i--) {
            if (output.get(i).getN() < minNumberOfPixelsInEdge) {
                output.remove(i);
            }
        }
    }
    
    private void pruneSpurs(List<PairIntArray> tmpEdges) {
        
        //TODO: improve this to follow a spur when found and remove the
        //  resulting spurs from the remove action
        
        int nIter = 0;
        int nMaxIter = 3;
        int nRemoved = 0;
                
        while (nIter < nMaxIter) {
            
            nRemoved = 0;
            
            /*
            0 0 0    0 0 0
            0 1 0    0 1 0
            0 _ _    _ _ 0
            */
            // indexes that have to be zeros that is, not within an edge
            int[] topXIdx = new int[]{-1, 0, 1, -1, 1};
            int[] topYIdx = new int[]{ 1, 1, 1,  0, 0};
            int[] leftXIdx = new int[]{-1};
            int[] leftYIdx = new int[]{-1};
            int[] rightXIdx = new int[]{1};
            int[] rightYIdx = new int[]{-1};
            // one must be a zero:
            int[] leftOrZeroXIdx = new int[]{0, 1};
            int[] leftOrZeroYIdx = new int[]{-1, -1};
            int[] rightOrZeroXIdx = new int[]{-1, 0};
            int[] rightOrZeroYIdx = new int[]{-1, -1};
        
            for (int r = 0; r < 4; r++) {
                if (r > 0) {
                    rotateIndexesBy90(topXIdx, topYIdx);
                    rotateIndexesBy90(leftXIdx, leftYIdx);
                    rotateIndexesBy90(rightXIdx, rightYIdx);
                    rotateIndexesBy90(leftOrZeroXIdx, leftOrZeroYIdx);
                    rotateIndexesBy90(rightOrZeroXIdx, rightOrZeroYIdx);
                }
                for (int i = 0; i < tmpEdges.size(); i++) {
                    
                    PairIntArray edge = tmpEdges.get(i);
                    
                    // skip the endpoints
                    for (int j = (edge.getN() - 1); j > 0; j--) {
                        
                        if ((j < 0) || (j > (edge.getN() - 1))) {break;}
                        
                        int x = edge.getX(j);
                        int y = edge.getY(j);
                        
                        // TODO: consider a smaller range to search than +-10
                        int start = j - 10;
                        if (start < 0) {
                            continue;
                        }
                        int stop = j + 10;
                        if (stop > (edge.getN() - 1)) {
                            continue;
                        }
                        // "notFound" is "all are zeroes"
                        boolean notFound = notFound(edge, topXIdx, topYIdx, x, y, 
                            start, stop);
                        if (notFound) {
                            notFound = notFound(edge, leftXIdx, leftYIdx, x, y, 
                                start, stop);
                            if (notFound) {
                                if (hasAtLeastOneZero(edge, leftOrZeroXIdx,
                                    leftOrZeroYIdx, x, y, start, stop)) {
                                    
                                    edge.removeRange(j, j);
                                    nRemoved++;
                                }
                            } else {
                                notFound = notFound(edge, rightXIdx, rightYIdx, x, y, 
                                    start, stop);
                                if (notFound) {
                                    if (hasAtLeastOneZero(edge, rightOrZeroXIdx,
                                        rightOrZeroYIdx, x, y, start, stop)) {
                                        
                                        edge.removeRange(j, j);
                                        nRemoved++;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            log.info("nRemoved=" + nRemoved);
            if (nRemoved == 0) {
                break;
            }
            nIter++;
        }
    }
    
    /**
     * edge contains point pairs of (x, y).  this method searches edge to assert
     * that all of the specified points are not present, else returns as soon
     * as one of the specified points is found.  
     * The specified points are 
     * (xIndex + xIndexOffsets[i], yIndex + yIndexOffsets[i]).
     * The search through edge is done from startIdx to stopIdx, inclusive.
     * @param edge
     * @param xIndexOffsets
     * @param yIndexOffsets
     * @param xIndex
     * @param yIndex
     * @param startIdx
     * @param stopIdx
     * @return 
     */
    private boolean notFound(PairIntArray edge, 
        int[] xIndexOffsets, int[] yIndexOffsets,
        int xIndex, int yIndex, int startIdx, int stopIdx) {
        
StringBuilder sb2 = new StringBuilder();
        
        for (int j = 0; j < xIndexOffsets.length; j++) {
            int xFind = xIndex + xIndexOffsets[j];
            int yFind = yIndex + yIndexOffsets[j];
sb2.append(String.format("(%d, %d)\n", xFind, yFind));
            for (int i = startIdx; i <= stopIdx; i++) {
                int x = edge.getX(i);
                int y = edge.getY(i);
                if ((x == xFind) && (y == yFind)) {
                    return false;
                }
            }
        }
        
StringBuilder sb = new StringBuilder();
for (int i = 0; i < edge.getN(); i++) {
     int x = edge.getX(i);
     int y = edge.getY(i);
     sb.append(String.format("%d)  (%d, %d)\n", i, x, y));
}
        
        return true;
    }
    
    /**
     * edge contains point pairs of (x, y).  this method searches edge to assert
     * that all of the specified points are not present, else returns as soon
     * as one of the specified points is found.  
     * The specified points are 
     * (xIndex + xIndexOffsets[i], yIndex + yIndexOffsets[i]).
     * The search through edge is done from startIdx to stopIdx, inclusive.
     * NOTE: also checks for whether nulling the pixel at (xP, yP) is 
     * connected to the wall and returns false if it is.
     * @param edge
     * @param xIndexOffsets
     * @param yIndexOffsets
     * @param xIndex
     * @param yIndex
     * @param startIdx
     * @param stopIdx
     * @return 
     */
    private boolean hasAtLeastOneZero(PairIntArray edge, 
        int[] xIndexOffsets, int[] yIndexOffsets,
        int xP, int yP, int startIdx, int stopIdx) {
        
        // side logic of checking whether connected to wall.  return false if so
        if ((xP == 0) || (yP == 0)) {
            return false;
        }
        if ((xP == (img.getWidth() - 1)) 
            || (yP == (img.getHeight() - 1))) {
            return false;
        }
        
        int nOnes = 0;
        for (int j = 0; j < xIndexOffsets.length; j++) {
            int xFind = xP + xIndexOffsets[j];
            int yFind = yP + yIndexOffsets[j];
            for (int i = startIdx; i <= stopIdx; i++) {
                int x = edge.getX(i);
                int y = edge.getY(i);
                if ((x == xFind) && (y == yFind)) {
                    nOnes++;
                }
            }
        }        
        
        return (nOnes < (xIndexOffsets.length));
    }
    
    private void rotateIndexesBy90(int[] xOffsetIndexes, int[] yOffsetIndexes) {
        
        for (int i = 0; i < xOffsetIndexes.length; i++) {
            int xoff = xOffsetIndexes[i];
            int yoff = yOffsetIndexes[i];
           
            switch(xoff) {
                case -1:
                    switch (yoff) {
                        case -1:
                            xOffsetIndexes[i] = -1;
                            yOffsetIndexes[i] = 1;
                            break;
                        case 0:
                            xOffsetIndexes[i] = 0;
                            yOffsetIndexes[i] = 1;
                            break;
                        case 1:
                            xOffsetIndexes[i] = 1;
                            yOffsetIndexes[i] = 1;
                            break;
                    }
                    break;
                case 0:
                    switch (yoff) {
                        case -1:
                            xOffsetIndexes[i] = -1;
                            yOffsetIndexes[i] = 0;
                            break;
                        case 0:
                            // remains same
                            break;
                        case 1:
                            xOffsetIndexes[i] = 1;
                            yOffsetIndexes[i] = 0;
                            break;
                    }
                    break;
                case 1:
                    switch (yoff) {
                        case -1:
                            xOffsetIndexes[i] = -1;
                            yOffsetIndexes[i] = -1;
                            break;
                        case 0:
                            xOffsetIndexes[i] = 0;
                            yOffsetIndexes[i] = -1;
                            break;
                        case 1:
                            xOffsetIndexes[i] = 1;
                            yOffsetIndexes[i] = -1;
                            break;
                    }
                    break;
            }
        }
        /*
        1 2 3    7 4 1  
        4 5 6      5 2  
        7 _ _      6 3 
        
        (-1,1)    (1,1)   0
        (0,1)     (1,0)   1
        (1,1)     (1,-1)  2
        (-1,0)    (0,1)   3
        (0,0)     (0,0)   4
        (1,0)     (0,-1)  5
        (-1,-1)   (-1,1)  6
        (0,-1)    (-1,0)  7
        (1,-1)    (-1,-1) 8
        */
    }
    
    private void adjustEdgesTowardsBrightPixels(List<PairIntArray> tmpEdges) {
        
        if (edgeGuideImage == null) {
            return;
        }
        
        int nReplaced = 1;
        
        int nMaxIter = 100;
        int nIter = 0;
      
        while ((nIter < nMaxIter) && (nReplaced > 0)) {
           
           nReplaced = 0;
        
           for (int lIdx = 0; lIdx < tmpEdges.size(); lIdx++) {
            
                PairIntArray edge = tmpEdges.get(lIdx);

                if (edge.getN() < 3) {
                    continue;
                }

                int nEdgeReplaced = 0;
                
                /*
                looking at the 8 neighbor region of each pixel for which the
                pixel's preceding and next edge pixels remain connected for
                   and among those, looking for a higher intensity pixel than
                   the center and if found, change coords to that.
                */
                for (int i = 1; i < (edge.getN() - 1); i++) {
                    int x = edge.getX(i);                
                    int y = edge.getY(i);
                    int prevX = edge.getX(i - 1);                
                    int prevY = edge.getY(i - 1);
                    int nextX = edge.getX(i + 1);                
                    int nextY = edge.getY(i + 1);

                    int maxValue = edgeGuideImage.getValue(x, y);
                    int maxValueX = x;
                    int maxValueY = y;
                    boolean changed = false;

                    for (int col = (prevX - 1); col <= (prevX + 1); col++) {

                        if ((col < 0) || (col > (edgeGuideImage.getWidth() - 1))) {
                            continue;
                        }

                        for (int row = (prevY - 1); row <= (prevY + 1); row++) {

                            if ((row < 0) || (row > (edgeGuideImage.getHeight() - 1))) {
                                continue;
                            }
                            if ((col == prevX) && (row == prevY)) {
                                continue;
                            }
                            if ((col == nextX) && (row == nextY)) {
                                continue;
                            }

                            // skip if pixel is not next to (nextX, nextY)
                            int diffX = Math.abs(nextX - col);
                            int diffY = Math.abs(nextY - row);
                            if ((diffX > 1) || (diffY > 1)) {
                                continue;
                            }

                            if (edgeGuideImage.getValue(col, row) > maxValue) {
                                maxValue = edgeGuideImage.getValue(col, row);
                                maxValueX = col;
                                maxValueY = row;
                                changed = true;
                            }
                        }
                    }
                    if (changed) {
                        nEdgeReplaced++;
                        edge.set(i, maxValueX, maxValueY);
                    }                    
                }
                
                nReplaced += nEdgeReplaced;
            }

            log.fine("REPLACED: " + nReplaced + " nIter=" + nIter);

            nIter++;
        }
        
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        
        curveHelper.removeRedundantPoints(tmpEdges);
        
        curveHelper.pruneAdjacentNeighborsTo2(tmpEdges);
        
        curveHelper.correctCheckeredSegments(tmpEdges);
    }

    private void processNeighbor(int uX, int uY, int uIdx,
        int vX, int vY, int vIdx,
        int[] uNodeEdgeIdx, List<PairIntArray> output) {
        
        // if u is not in an edge already, create a new one
        if (uNodeEdgeIdx[uIdx] == -1) {

            PairIntArray edge = new PairIntArray();

            edge.add(uX, uY);

            uNodeEdgeIdx[uIdx] = output.size();

            output.add(edge);                        
        }

        // keep the curve points ordered

        // add v to the edge u if u is the last node in it's edge

        PairIntArray appendToNode = output.get(uNodeEdgeIdx[uIdx]);
        int aIdx = appendToNode.getN() - 1;
        if ((appendToNode.getX(aIdx) != uX) || 
            (appendToNode.getY(aIdx) != uY)) {
            return;
        }

        appendToNode.add(vX, vY);

        uNodeEdgeIdx[vIdx] = uNodeEdgeIdx[uIdx];
        
    }
}
