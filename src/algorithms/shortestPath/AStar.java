package algorithms.shortestPath;

import algorithms.imageProcessing.DoubleLinkedCircularList;
import algorithms.imageProcessing.Heap;
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

    protected boolean calculateHeuristics = false;

    // key is total estimate from srcIdx to destIdx for the given refIdx
    //    (that is the distance from srcIdx to refIdx + refIdx + heuristic)
    protected final Heap heap = new Heap();

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
        Arrays.fill(heuristics, Long.MAX_VALUE);

        initHeap();
    }

    private void initHeap() {

        int n = points.length;

        nodes = new HeapNode[n];

        // initialize all except the source node as having infinite distance
        for (int i = 0; i < points.length; ++i) {

            long dist = (i == sourceIndx) ? 0 : Long.MAX_VALUE;

            HeapNode node = new HeapNode(dist);
            node.setData(Integer.valueOf(i));

            heap.insert(node);

            distFromS[i] = dist;

            nodes[i] = node;
        }
    }

    public int[] search() {

        int count = 0;

        HeapNode uNode = heap.extractMin();

        while (uNode != null) {

            int uIndx = ((Integer)uNode.getData()).intValue();
heap.printRootList();
            LinkedList<Integer> adj = adjList.get(uIndx);

            if (adj == null) {
                continue;
            }

            Integer vIndex = adj.poll();

            while (vIndex != null) {

                int vIndx = vIndex.intValue();

                if (distFromS[uIndx] == Long.MAX_VALUE) {
                    vIndex = adj.poll();
                    continue;
                }

                long distUV = distBetween(uIndx, vIndx);

                long uDistPlusCost = (distFromS[uIndx] + distUV);

                long vDist = (distFromS[vIndex] == Long.MAX_VALUE) ?
                    Long.MAX_VALUE : distFromS[vIndex];

                // relax(u,v, wght)
                if (uDistPlusCost < vDist) {

                    long vDistFromSrc = distFromS[uIndx] + distUV;

                    long vDistPlusHeuristic = vDistFromSrc + heuristic(vIndx);

                    HeapNode vNode = nodes[vIndx];

                    if (vDistPlusHeuristic < vNode.getKey()) {

log.info("decrease key: " + vNode.toString() + " to key=" + vDistPlusHeuristic);

                        heap.decreaseKey(vNode, vDistPlusHeuristic);
heap.printRootList();
                        distFromS[vIndex] = vDistFromSrc;

                        prevNode[vIndex] = uNode;

                        count++;

                    } else {

                        String str = String.format(
                        "did not decrease key for: u=%d v=%d uDist+distUV=%d vDist=%d currentVDistPlusHeuristic=%d rejected vDistPlusH=%d",
                            uIndx, vIndx, (int)uDistPlusCost, (int)vDist,
                            (int)vNode.getKey(), (int)vDistPlusHeuristic);
                        System.err.println(str);
                    }
                }

                vIndex = adj.poll();
            }

            uNode = heap.extractMin();
        }

        if (count == 0) {
            return null;
        }

        int[] pathIndexes = getNodeIndexesToDestination();

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

        long dist = (long)Math.round(diffX*diffX + diffY*diffY);

        return dist;
    }

    private long heuristic(int vIndx) {

        if (heuristics[vIndx] == Long.MAX_VALUE) {
            // calc straight line distance from vIndx to destination
            heuristics[vIndx] = distBetween(vIndx, destinationIndx);
        }

        return heuristics[vIndx];
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

}
