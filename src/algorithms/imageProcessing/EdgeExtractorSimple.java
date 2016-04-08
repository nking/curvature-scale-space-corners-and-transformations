package algorithms.imageProcessing;

import algorithms.misc.Misc;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.PairIntArrayWithColor;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Stack;

/**
 * The simplest edge extractor which does the following:
 * <pre>
 *     -- finds junctions and stores them in a set
 *     -- dfs traversal of each point w/ pixel value > 0.
 *           searches for a neighbor and if junction is within
 *           neighbors, terminates the curve and starts a new curve
 *           with the next popped point, else adds one neighbor to
 *           the curve and that to the stack.
 *     -- joins curves whose endpoints are not adjacent to junctions and whose
 *            endpoints are adjacent.
 *     -- performs a correction for an artifact pattern in the edges that tends 
 *        to be present due to the curve thinning plus concatenation operations 
 *        on the given image.
 * </pre>
 * 
 * @author nichole
 */
public class EdgeExtractorSimple {
    
    private Set<PairInt> junctions = new HashSet<PairInt>();
    
    private List<PairIntArray> theEdges = new ArrayList<PairIntArray>();
    
    private final Set<PairInt> points = new HashSet<PairInt>();
    
    private boolean finished = false;
    
    private final int n0;
    private final int n1;
    
    public EdgeExtractorSimple(int[][] edgeImage) {
        
        n0 = edgeImage.length;
        n1 = edgeImage[0].length;
        
        for (int i = 0; i < n0; ++i) {
            for (int j = 0; j < n1; ++j) {
                if (edgeImage[i][j] == 0) {
                    continue;
                }
                points.add(new PairInt(i, j));
            }
        }
        
    }
    
    public EdgeExtractorSimple(Set<PairInt> edgePoints, int imageWidth, int imageHeight) {
        
        n0 = imageWidth;
        n1 = imageHeight;
        points.addAll(edgePoints);
    }
    
    public EdgeExtractorSimple(GreyscaleImage img) {
        
        n0 = img.getWidth();
        n1 = img.getHeight();
        
        for (int i = 0; i < n0; ++i) {
            for (int j = 0; j < n1; ++j) {
                if (img.getValue(i, j) == 0) {
                    continue;
                }
                points.add(new PairInt(i, j));
            }
        }
        
    }
    
    public void extractEdges() {
        
        if (finished) {
            throw new IllegalStateException("edges have already been extracted");
        }
                
        findJunctions();
        
        findConnectedPoints();
        
        joinEdges();
                
        setEdgesAsOrdered();
        
        finished = true;
    }

    private void findJunctions() {
        
        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;
        
        // any point having more than 2 neighbors is a junctions
        
        for (PairInt p : points) {
            
            int x = p.getX();
            int y = p.getY();
            
            int count = 0;
            
            for (int k = 0; k < dxs.length; ++k) {
                int x2 = x + dxs[k];
                int y2 = y + dys[k];
                PairInt p2 = new PairInt(x2, y2);
                if (points.contains(p2)) {
                    count++;
                }
            }
            
            if (count > 2) {
                junctions.add(p);
            }
        }
    }
    
    private void findConnectedPoints() {
        
        Stack<PairInt> stack = new Stack<PairInt>();
        for (PairInt p : points) {
            if (!junctions.contains(p)) {
                stack.add(p);
            }
        }
        
        Set<PairInt> added = new HashSet<PairInt>();
        
        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;
        
        PairIntArray currentEdge = new PairIntArray();
        
        while (!stack.isEmpty()) {
            
            PairInt uPoint = stack.pop();
            
            if (added.contains(uPoint)) {
                continue;
            }
                        
            int x = uPoint.getX();
            int y = uPoint.getY();
            
            currentEdge.add(x, y);
            
            added.add(new PairInt(x, y));
            
            boolean containsJunction = false;
            
            PairInt aNeighbor = null;
                        
            for (int k = 0; k < dxs.length; ++k) {
                int x2 = x + dxs[k];
                int y2 = y + dys[k];
                PairInt p2 = new PairInt(x2, y2);
                if (!points.contains(p2) || added.contains(p2)) {
                    continue;
                }
                if (junctions.contains(p2)) {
                    containsJunction = true;
                    break;
                }
                if (aNeighbor == null) {
                    aNeighbor = p2;
                }
            }
            
            if (containsJunction || (aNeighbor == null)) {
                // store currentEdge and start a new one
                theEdges.add(currentEdge);
                currentEdge = new PairIntArray();
            } else {
                stack.add(aNeighbor);
            }            
        }
        
        if (currentEdge.getN() > 0) {
            theEdges.add(currentEdge);
        }
    }

    public Set<PairInt> getJunctions() {
        return junctions;
    }
    
    public List<PairIntArray> getEdges() {
        return theEdges;
    }

    /**
     * join edges that are adjacent to one another, but not in junctions.
     * runtime complexity is O(N_edges) + 
     */
    private void joinEdges() {
        
        Map<PairInt, Integer> endpointMap = new HashMap<PairInt, Integer>();
        
        /*
         * for each edge:
              if endpoint is not a junction, store in endpoint map
           for each edge:
              if endpoint is not a junction:
                  search for a neighbor in endpoint map and if found,
                     aggregate line with current line
                     remove both endpoints where lines are joining
                     update the index for the other line endpoint in the map
                     clear other line in edges
            vist edges in reverse order and remove the empty entries.
         */
        
        for (int edgeIdx = 0; edgeIdx < theEdges.size(); ++edgeIdx) {
            
            PairIntArray edge = theEdges.get(edgeIdx);
            
            int n = edge.getN();
            
            if (n == 0) {
                continue;
            }
                        
            PairInt p0 = new PairInt(edge.getX(0), edge.getY(0));
            
            PairInt p1 = new PairInt(edge.getX(n - 1), edge.getY(n - 1));
            
            Integer index = Integer.valueOf(edgeIdx);
            
            if (!junctions.contains(p0)) {
                endpointMap.put(p0, index);
            }
            
            if (!junctions.contains(p1)) {
                endpointMap.put(p1, index);
            }
        }
        
        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;

        PairInt[] eps = new PairInt[2];
        
        for (int edgeIdx = 0; edgeIdx < theEdges.size(); ++edgeIdx) {
            
            PairIntArray edge = theEdges.get(edgeIdx);
            
            int n = edge.getN();
            
            if (n == 0) {
                continue;
            }
            
            eps[0] = new PairInt(edge.getX(0), edge.getY(0));
            
            eps[1] = (n == 1) ? null : new PairInt(edge.getX(n - 1), edge.getY(n - 1));
            
            Integer index = Integer.valueOf(edgeIdx);
            
            for (PairInt ep : eps) {
                
                if (ep == null) {
                    continue;
                }
                                 
                for (int k = 0; k < dxs.length; ++k) {

                    int x2 = ep.getX() + dxs[k];
                    int y2 = ep.getY() + dys[k];

                    PairInt p2 = new PairInt(x2, y2);

                    Integer index2 = endpointMap.get(p2);

                    if (index2 == null || index.equals(index2)) {
                        continue;
                    }

                    PairIntArray edge2 = theEdges.get(index2.intValue());
                    assert(edge2.getN() > 0);

                    int n2 = edge2.getN();

                    // join the 2 edges
                    if ((ep.getX() == edge.getX(0)) && 
                        (ep.getY() == edge.getY(0))) {

                        if ((p2.getX() == edge2.getX(0)) && 
                            (p2.getY() == edge2.getY(0))) {

                            /*
                            e1_2
                            e1_1
                            e1_0 *
                                   * e2_0  e2_1  
                            */

                            if (n2 > 1) {
                                endpointMap.remove(p2);
                            }
                            if (n > 1) {
                                endpointMap.remove(ep);
                            }
                            PairInt otherEndpoint = new PairInt(
                                edge2.getX(n2 - 1), edge2.getY(n2 - 1));
                            edge.reverse();
                            edge.addAll(edge2);
                            theEdges.set(index2.intValue(), new PairIntArray());

                            if (!junctions.contains(otherEndpoint)) {
                                endpointMap.put(otherEndpoint, index);
                            }

                        } else {
                            // edge2 join is end
                            assert(p2.getX() == edge2.getX(n2 - 1));
                            assert(p2.getY() == edge2.getY(n2 - 1));

                            /*
                            e1_2
                            e1_1
                            e1_0 *
                                   * e2_3  e2_2  e2_1  
                            */

                            if (n2 > 1) {
                                endpointMap.remove(p2);
                            }
                            if (n > 1) {
                                endpointMap.remove(ep);
                            }
                            PairInt otherEndpoint = new PairInt(
                                edge2.getX(0), edge2.getY(0));
                            edge.reverse();
                            edge2.reverse();
                            edge.addAll(edge2);
                            theEdges.set(index2.intValue(), new PairIntArray());

                            if (!junctions.contains(otherEndpoint)) {
                                endpointMap.put(otherEndpoint, index);
                            }
                        }

                    } else {
                        // ep is at the end
                        assert(ep.getX() == edge.getX(n - 1));
                        assert(ep.getY() == edge.getY(n - 1));

                        if ((p2.getX() == edge2.getX(0)) && 
                            (p2.getY() == edge2.getY(0))) {

                            /*
                            e1_0
                            e1_1
                            e1_2 *
                                   * e2_0  e2_1  
                            */
                            if (n2 > 1) {
                                endpointMap.remove(p2);
                            }
                            if (n > 1) {
                                endpointMap.remove(ep);
                            }
                            PairInt otherEndpoint = new PairInt(
                                edge2.getX(n2 - 1), edge2.getY(n2 - 1)); 
                            edge.addAll(edge2);
                            theEdges.set(index2.intValue(), new PairIntArray());

                            if (!junctions.contains(otherEndpoint)) {
                                endpointMap.put(otherEndpoint, index);
                            }

                        } else {

                            // edge2 join is end
                            assert(p2.getX() == edge2.getX(n2 - 1));
                            assert(p2.getY() == edge2.getY(n2 - 1));

                            /*
                            e1_0
                            e1_1
                            e1_2 *
                                   * e2_3  e2_2  
                            */
                            if (n2 > 1) {
                                endpointMap.remove(p2);
                            }
                            if (n > 1) {
                                endpointMap.remove(ep);
                            }
                            PairInt otherEndpoint = new PairInt(
                                edge2.getX(0), edge2.getY(0)); 
                            edge2.reverse();
                            edge.addAll(edge2);
                            theEdges.set(index2.intValue(), new PairIntArray());

                            if (!junctions.contains(otherEndpoint)) {
                                endpointMap.put(otherEndpoint, index);
                            }
                        }
                    }
                }
            }
        }
        
        for (int i = (theEdges.size() - 1); i > -1; --i) {
            
            PairIntArray edge = theEdges.get(i);
            
            if (edge.getN() == 0) {
                theEdges.remove(i);
            }
        }
        
    }

    private void setEdgesAsOrdered() {
        
        for (int i = 0; i < theEdges.size(); ++i) {
            
            PairIntArray edge = theEdges.get(i);
            
            PairIntArrayWithColor e2 = new PairIntArrayWithColor(edge);
            
            e2.setAsOrderedCurve();
            
            theEdges.set(i, e2);
        }
        
    }
    
}
