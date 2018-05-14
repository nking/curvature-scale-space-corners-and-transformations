package algorithms.shortestPath;

import algorithms.bipartite.MinHeapForRT2012;
import algorithms.imageProcessing.HeapNode;
import algorithms.util.PairInt;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.logging.Logger;

/**
 * An extension of Dijkstra's algorithm that uses heuristics to improve the
 * performance.  runtime performance O(|E| + |V| log |V|).
 *
 * http://en.wikipedia.org/wiki/A*_search_algorithm
 *
 * "As A* traverses the graph, it follows a path of the lowest known cost,
 * keeping a sorted priority queue [minimum heap] of alternate path segments
 * along the way.   If, at any point, a segment of the path being traversed has
 * a higher cost than another encountered path segment, it abandons the
 * higher-cost path segment and traverses the lower-cost path segment instead.
 * This process continues until the goal is reached."
 *   -- maintains a minimum heap of nodes to be traversed = the open set.
 *   -- the implementation below uses a breadth first search w/ depth=1
 *
 * If heuristics aren't given to the code, it makes an assumption that all
 * nodes lie within a straight line distance to the destination node and hence
 * calculates the heuristic based upon that distance.
 *
 * * Variables:
 *     g[n] is the shortest distance traveled path from the sourceIndx to the
 *          node n
 *     h[n] is the smallest estimated cost from the node n to destinationIndx
 *     f[n] is the total search cost from sourceIndex to node n
 *          f(n) = g(n) + h(n)
 * Goal:
 *     find the path that creates the smallest f[destinationIndx]
 *
 * @author nichole
 */
public class AStar {

    //TODO: make a constructor that supplies the largest distance possible
    //    between 2 points, and replace the
    //    Heap w/ MinHeapForRT2012 which has an O(N) bucket queue internally
    //    as one implementation.
    
    protected boolean calculateHeuristics = false;
    
    private int sentinel = Integer.MAX_VALUE;
    
    private int maxDist = sentinel;

    // key is total estimate from srcIdx to destIdx for the given refIdx
    //    (that is the distance from srcIdx to refIdx + refIdx + heuristic)
    protected MinHeapForRT2012 heap = null;

    // refs to nodes internal to heap for decrease key operations
    protected HeapNode[] nodes = null;

    protected final PairInt[] points;

    protected final List<LinkedList<Integer>> adjList;

    // this is the ** g[] ** in the class comments
    protected final long[] distFromS;

    protected final HeapNode[] prevNode;

    /**
     * heuristic cost estimate is the cost left to reach the destionationIndx.
     * if heuristic isn't given to code, will use the straight line distance
     * between one node and another.
     * Note that the index is the same as for nodes.
     */
    protected final long[] heuristics;

    protected final int sourceIndx;
    protected final int destinationIndx;

    private enum State {
       INITIALIZED, COMPLETE
    }
    
    private State state = null;

    private Logger log = Logger.getLogger(this.getClass().getName());

    public AStar(final PairInt[] points, final List<LinkedList<Integer>> adjacencyList,
        final int srcIndx, final int destIndx) {

        if (points == null) {
            throw new IllegalArgumentException("points cannot be null");
        }

        int n = points.length;

        if (adjacencyList == null || (adjacencyList.size() != n)) {
            throw new IllegalArgumentException("adjacencyList cannot be null and must"
            + " be same size as points");
        }
        if ((srcIndx < 0) ||  (srcIndx > (n - 1))) {
                throw new IllegalArgumentException(
            "srcIndx is out of bounds of vertices array");
        }
        if ((destIndx < 0) ||  (destIndx > (n - 1))) {
                throw new IllegalArgumentException(
            "destIndx is out of bounds of vertices array");
        }

        this.points = points;
        this.adjList = new ArrayList<LinkedList<Integer>>();
        for (LinkedList<Integer> list : adjacencyList) {
            adjList.add(new LinkedList<Integer>(list));
        }

        sourceIndx = srcIndx;
        destinationIndx = destIndx;

        distFromS = new long[n];
        heuristics = new long[n];
        prevNode = new HeapNode[n];

        // populate w/ straight line distances as needed:
        for (int i = 0; i < points.length; ++i) {
            heuristics[i] = distBetween(i, destinationIndx);
        }

        initHeap();
        
        log.fine("src=" + points[srcIndx].toString() + " dest=" + 
            points[destIndx].toString());
    }

    private void initHeap() {
        
        int n = points.length;
        
        maxDist = calculateMaxDistance();
        
        int capacity = Math.max(n, maxDist + 1);
        
        int nBits = (int)Math.ceil(Math.log(maxDist/Math.log(2)));
        
        //int capacity, int approxN, int maxNumberOfBits
        heap = new MinHeapForRT2012(capacity, n, nBits);

        nodes = new HeapNode[n];

        // initialize all except the source node as having infinite distance
        for (int i = 0; i < points.length; ++i) {

            long dist = (i == sourceIndx) ? 0 : (maxDist + 1);

            HeapNode node = new HeapNode(dist);
            node.setData(Integer.valueOf(i));

            heap.insert(node);

            distFromS[i] = dist;

            nodes[i] = node;
        }
        
        state = State.INITIALIZED;
    }

    public int[] search() {

        int count = 0;

        HeapNode uNode = heap.extractMin();

        while (uNode != null) {

            int uIndx = ((Integer)uNode.getData()).intValue();
            
            // null the entry in nodes so it isn't used in decrease key
            nodes[uIndx] = null;
            
            LinkedList<Integer> adj = adjList.get(uIndx);

            if (adj == null) {
                continue;
            }

            Integer vIndex = adj.poll();

            while (vIndex != null) {

                int vIndx = vIndex.intValue();

                if ((distFromS[uIndx] == sentinel) || (nodes[vIndx] == null)) {
                    vIndex = adj.poll();                    
                    continue;
                }

                long distUV = distBetween(uIndx, vIndx);

                long uDistPlusCost = (distFromS[uIndx] + distUV);

                long vDist = (distFromS[vIndex] == sentinel) ?
                    sentinel : 
                    distFromS[vIndex];

                log.fine(points[uIndx].toString() + ":" + points[vIndx] + " "
                    + " dU=" + distFromS[uIndx] 
                    + " dU+UtoV=" + uDistPlusCost 
                    + " dV=" + vDist);
                    
                // relax(u,v, wght)
                if (uDistPlusCost < vDist) {

                    long vDistFromSrc = distFromS[uIndx] + distUV;

                    long vDistPlusHeuristic = vDistFromSrc + heuristics[vIndx];

                    HeapNode vNode = nodes[vIndx];

                    log.fine("dU+UtoV+H=" + vDistPlusHeuristic 
                        + " _dV=" + vNode.getKey());
                    
                    if (vDistPlusHeuristic < vNode.getKey()) {

                        //heap.decreaseKey(vNode, vDistPlusHeuristic);
                        heap.decreaseKey(vNode, vDistFromSrc);
                        
                        distFromS[vIndex] = vDistFromSrc;

                        prevNode[vIndex] = uNode;

                        count++;

                        log.fine(" ");
                    } else {

                        String str = String.format(
                        "  did not decrease v key for: u=%d v=%d uDist+distUV=%d vDist=%d currentVDistPlusHeuristic=%d rejected vDistPlusH=%d",
                            uIndx, vIndx, (int)uDistPlusCost, (int)vDist,
                            (int)vNode.getKey(), (int)vDistPlusHeuristic);
                        log.fine(str);
                    }
                }

                vIndex = adj.poll();
            }

            uNode = heap.extractMin();
        }

        if (count == 0) {
            
            state = State.COMPLETE;
            
            return null;
        }

        int[] pathIndexes = getNodeIndexesToDestination();

        state = State.COMPLETE;
        
        return pathIndexes;
    }

    /**
     * assuming that this is only invoked when adjacency has been checked already
     * @param uIdx
     * @param vIdx
     * @return
     */
    protected long distBetween(int uIdx, int vIdx) {

        int diffX = points[uIdx].getX() - points[vIdx].getX();
        int diffY = points[uIdx].getY() - points[vIdx].getY();

        long dist = (long)Math.round(Math.sqrt(diffX*diffX + diffY*diffY));

        return dist;
    }

    private int[] getNodeIndexesToDestination() {

        int[] destNodes = new int[prevNode.length];

        int count = prevNode.length - 1;
        int lastInd = destinationIndx;

        destNodes[count] = lastInd;

        count--;

        while (lastInd != this.sourceIndx) {

            HeapNode destNode = prevNode[lastInd];

            if (destNode == null) {
                return null;
            }

            lastInd = ((Integer)destNode.getData()).intValue();

            destNodes[count] = lastInd;

            count--;
        }

        int[] out = Arrays.copyOfRange(destNodes, count + 1, prevNode.length);

        return out;
    }

    /**
     * method to be used for use cases where one knows that the incomplete
     * path from search() is due to a certain type of junction 
     * (specifically, the diagonal pattern of the 
     * untraversable lobe remover class).  The returned
     * indexes are all visited nodes from the internal distance array which
     * should represent a loop for the specific junction.
     * Note that for wider use, the method could be adapted to make a 
     * partial path from the ordered distance array and then used with the
     * previous node array to complete the other half of the source to
     * source loop, possibly needing additional search for unvisited nodes.
     * The results are not currently used in an ordered manner so the
     * additional work for ordered points isn't useful currently.  The results 
     * are used to remove a sub-loop isolated by a specific type of junction 
     * in a closed curve.
     * runtime is O(|V|).
     * @return unordered indexes representing path points if the search was
     * through a specific type of junction (the diagonal pattern of the 
     * untraversable lobe remover class)
     */
    public int[] createSourceToSourcePath() {
        
        if (!state.equals(State.COMPLETE)) {
            search();
        }
        
        int[] distances = new int[distFromS.length];
        int[] indexes = new int[distances.length];
        
        int count = 0;
        
        for (int i = 0; i < distFromS.length; ++i) {
            long d = distFromS[i];
            if (d == sentinel) {
                continue;
            }
            distances[count] = (int)d;
            indexes[count] = i;
            count++;
        }
        
        if (count == 0) {
            return null;
        }
        
        //distances = Arrays.copyOf(distances, count);
        indexes = Arrays.copyOf(indexes, count);
        
        //MultiArrayMergeSort.sortByDecr(distances, indexes);
        
        return indexes;
    }

    private int calculateMaxDistance() {
        int minX = Integer.MAX_VALUE;
        int maxX = Integer.MIN_VALUE;
        int minY = Integer.MAX_VALUE;
        int maxY = Integer.MIN_VALUE;
        
        for (PairInt p : points) {
            int x = p.getX();
            int y = p.getY();
            if (x < minX) {
                minX = x;
            }
            if (x > maxX) {
                maxX = x;
            }
            if (y < minY) {
                minY = y;
            }
            if (y > maxY) {
                maxY = y;
            }
        }
        int dX = maxX - minX + 1;
        int dY = maxY - minY + 1;
        int dist = (int)Math.ceil(Math.sqrt(dX * dX + dY * dY));
        
        return dist;
    }
    
}
