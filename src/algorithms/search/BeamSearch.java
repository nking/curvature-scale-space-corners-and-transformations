package algorithms.search;

import algorithms.util.SimpleLinkedListNode;

import algorithms.FixedSizeSortedVector; 
import algorithms.util.LinkedListCostNode;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import thirdparty.edu.princeton.cs.algs4.Queue;

        
/**
from https://www.wikiwand.com/en/Beam_search
Beam Search uses a greedy, breadth-first search to build its search tree,
but only keeps top k (beam size) nodes at each level in memory.
The next level will then be expanded from these k nodes.

Complete search variants of beam search are made by combining it with 
depth-first search resulting in beam stack search and depth-first beam search, 
or combining it with limited discrepancy search resulting in beam search 
using limited discrepancy backtracking (BULB). 
The resulting search algorithms are anytime algorithms that find good but 
likely sub-optimal solutions quickly, like beam search, then backtrack and 
continue to find improved solutions until convergence to an optimal solution.

breadth-first search:
 given a graph, it visits every node reachable from s and computes the
 distance to each of those nodes.  the predecessor tree is those nodes
 reachable from s.   the distance array can be used to find the
 shortest path.

 * @author nichole
 */
public class BeamSearch {

    protected int[] d = null;

    protected int[] predecessor = null;

    protected int[] color = null;

    protected SimpleLinkedListNode[] adjacencyList = null;

    // source index
    protected final int s;
    protected final int beamSize;

    /**
     * constructor with adjacency list, with default equal cost edges.
     * @param theAdjacencyList
     * @param sourceNodeIndex start index for search
     * @param beamSize in the breadth-first search, for each set of child nodes, 
     * only the smallest cost nodes are expanded to continue search and the
     * limited number of those smallest cost nodes is beamSize.
     */
    public BeamSearch(SimpleLinkedListNode[] theAdjacencyList, int sourceNodeIndex,
            final int beamSize) {

        this.adjacencyList = theAdjacencyList;
        this.d = new int[adjacencyList.length];
        this.predecessor = new int[adjacencyList.length];
        this.color = new int[adjacencyList.length];
        if (sourceNodeIndex < 0) {
            throw new IllegalArgumentException("sourceNodeIndex cannot be a negative number");
        }
        if (sourceNodeIndex >= adjacencyList.length) {
            throw new IllegalArgumentException("sourceNodeIndex must be an index within size limits of theAdjacencyList");
        }
        this.s = sourceNodeIndex;
        
        if (beamSize <= 0) {
            throw new IllegalArgumentException("beamSize must be a positive number greater than 0");
        }
        this.beamSize = beamSize;
    }
    
    /**
     * constructor with edge costs in adjacency list.
     * @param theAdjacencyList adjacency list that includes edge costs
     * @param sourceNodeIndex start index for search
     * @param beamSize in the breadth-first search, for each set of child nodes, 
     * only the smallest cost nodes are expanded to continue search and the
     * limited number of those smallest cost nodes is beamSize.
     */
    public BeamSearch(LinkedListCostNode[] theAdjacencyList, int sourceNodeIndex,
            final int beamSize) {

        this.adjacencyList = theAdjacencyList;
        this.d = new int[adjacencyList.length];
        this.predecessor = new int[adjacencyList.length];
        this.color = new int[adjacencyList.length];
        if (sourceNodeIndex < 0) {
            throw new IllegalArgumentException("sourceNodeIndex cannot be a negative number");
        }
        if (sourceNodeIndex >= adjacencyList.length) {
            throw new IllegalArgumentException("sourceNodeIndex must be an index within size limits of theAdjacencyList");
        }
        this.s = sourceNodeIndex;
        
        if (beamSize <= 0) {
            throw new IllegalArgumentException("beamSize must be a positive number greater than 0");
        }
        this.beamSize = beamSize;
    }

    /**
     * add source index to the bfs tree.
     * @return list of indexes searched
     */
    public TIntList search() {

        TIntList searched = new TIntArrayList();
        
        // initialize
        for (int u = 0; u < adjacencyList.length; u++) {

            SimpleLinkedListNode vNode  = adjacencyList[u];

            if (vNode != null) {

                int v = vNode.getKey();

                setColorToWhite(v);

                d[v] = Integer.MAX_VALUE;

                predecessor[v] = -1;
             }
        }

        setColorToGray(s);

        d[s] = 0;

        predecessor[s] = -1;

        Queue queue = new Queue();
        queue.enqueue(s);

        while (!queue.isEmpty()) {

            int u = (Integer)queue.dequeue();

            searched.add(u);

            SimpleLinkedListNode vNode = adjacencyList[u];

            FixedSizeSortedVector<Index> sorted = 
                new FixedSizeSortedVector(beamSize, Index.class);
        
            while (vNode != null) {

                int v = vNode.getKey();
                
                if (isColorWhite(v)) {

                    setColorToGray(v);

                    if (vNode instanceof LinkedListCostNode) {
                        d[v] = d[u] + ((LinkedListCostNode)vNode).getCost();
                    } else {
                        d[v] = d[u] + 1;
                    }

                    predecessor[v] = u;

                    sorted.add(new Index(v));
                }

                vNode = vNode.getNext();
            }
            
            Index[] sortedArray = sorted.getArray();
            for (int z = 0; z < sorted.getNumberOfItems(); z++) {
                queue.enqueue(sortedArray[z].idx);
            }

            setColorToBlack(u);
        }
        
        return searched;
    }

    protected void setColorToWhite(int nodeIndex) {
        color[nodeIndex] = 0;
        //System.out.println(nodeIndex + " = white");
    }
    protected void setColorToGray(int nodeIndex) {
        color[nodeIndex] = 1;
        //System.out.println(nodeIndex + " = gray");
    }
    protected void setColorToBlack(int nodeIndex) {
        color[nodeIndex] = 2;
        //System.out.println(nodeIndex + " = black");
    }
    protected boolean isColorWhite(int nodeIndex) {
        return (color[nodeIndex] == 0);
    }

    private class Index implements Comparable<Index> {

        private final int idx;
        
        public Index(int index) {
            this.idx = index;
        }
        
        @Override
        public int compareTo(Index c) {
            int c1 = d[idx];
            int c2 = d[c.idx];
            return Integer.compare(c1, c2);
        }
        
    }
}
