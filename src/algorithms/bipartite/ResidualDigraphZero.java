package algorithms.bipartite;

import algorithms.bipartite.MinCostUnbalancedAssignment.LeftNode;
import algorithms.bipartite.MinCostUnbalancedAssignment.PathNode;
import algorithms.bipartite.MinCostUnbalancedAssignment.RightNode;
import algorithms.bipartite.MinCostUnbalancedAssignment.SinkNode;
import algorithms.bipartite.MinCostUnbalancedAssignment.SourceNode;
import algorithms.imageProcessing.DoubleLinkedCircularList;
import algorithms.imageProcessing.HeapNode;
import algorithms.misc.Misc;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * a residual digraph for use within refine method in
 * MinCostUnbalancedAssignment.java = NOTE that unlike
 * ResidualDiGraph2, this one may have multiple mappings
 * for right to left backward links because they are
 * potential paths that have not been filtered to
 * keep only the compatible yet.
 * 
 * @author nichole
 */
public class ResidualDigraphZero {
    
    /**
     * <pre>
     * |-
     * </pre>
     */
    private final int sourceNode;
    
    /**
     * <pre>
     * -|
     * </pre>
     */
    private int sinkNode;

    /**
     * left (==X) vertices in graph
     */
    private final int nLeft;

    /**
     * right (==Y) vertices in graph
     */
    private final int nRight;

    /**
     * links X to Y (that is, left to right).
     */
    private Map<Integer, Set<Integer>> forwardLinksRM
        = new HashMap<Integer, Set<Integer>>();

    /**
     * links Y to X (that is, right to left).
     */
    private Map<Integer, Set<Integer>>backwardLinksRM
        = new HashMap<Integer, Set<Integer>>();

    /**
     * source to Left links are represented by this
     * set of Left indexes.  These are "saturated" 
     * Left vertexes at initialization from being
     * connected to a matched node.
     */
    private Set<Integer> forwardLinksSourceRM =
        new HashSet<Integer>();
    
    /**
     * Left to source links are represented by this set of
     * Left vertexes.  These are "idle" Left vertexes.
     */
    private Set<Integer> backwardLinksSourceRM =
        new HashSet<Integer>();
    
    /**
     * Right to sink links are represented by this set of
     * right vertexes.  These are "saturated" Right vertexes
     * at initialization, from being connected to a matched node.
     */
    private Set<Integer> forwardLinksSinkRM =
        new HashSet<Integer>();
    
    /**
     * sink to Right links are represented by this set of
     * Right vertexes.  These are "idle" Right vertexes.
     */
    private Set<Integer> backwardLinksSinkRM =
        new HashSet<Integer>();
    
    public ResidualDigraphZero(int nLeft, int nRight,
        int sourceNode, int sinkNode,
        DoubleLinkedCircularList[] forest) {
    
        this.nLeft = nLeft;
        this.nRight = nRight;
        this.sourceNode = sourceNode;
        this.sinkNode = sinkNode;
       
        // traverse trees in the forest
        for (int forestIdx = 0; forestIdx < forest.length; ++forestIdx) {
            DoubleLinkedCircularList tree = forest[forestIdx];
            if (tree == null) {
                continue;
            }
            long n = tree.getNumberOfNodes();
            HeapNode node = tree.getSentinel();
            int treeIdx = 0;
            while (treeIdx < n) {
                node = node.getRight();
                List<PathNode> path = MinCostUnbalancedAssignment.extractNodes(node);
                Misc.<PathNode>reverse(path);
                for (int branchIdx = 0; branchIdx < (path.size() - 1); ++branchIdx) {
                    PathNode node1 = path.get(branchIdx);
                    PathNode node2 = path.get(branchIdx + 1);
                    int l1 = (int) node1.getKey();
                    int l2 = (int) node2.getKey();
                    if (l2 != 0) {
                        continue;
                    }
                    Integer index1 = (Integer) node1.getData();
                    Integer index2 = (Integer) node2.getData();
                    if (node1 instanceof LeftNode) {
                        assert (!(node2 instanceof LeftNode));
                        if (node2 instanceof RightNode) {
                            Set<Integer> indexes = forwardLinksRM.get(index1);
                            if (indexes == null) {
                                indexes = new HashSet<Integer>();
                                forwardLinksRM.put(index1, indexes);
                            }
                            indexes.add(index2);
                        } else {
                            assert (node2 instanceof SourceNode);
                            backwardLinksSourceRM.add(index1);
                        }
                    } else if (node1 instanceof RightNode) {
                        assert (!(node2 instanceof RightNode));
                        if (node2 instanceof LeftNode) {
                            Set<Integer> indexes = backwardLinksRM.get(index1);
                            if (indexes == null) {
                                indexes = new HashSet<Integer>();
                                backwardLinksRM.put(index1, indexes);
                            }
                            indexes.add(index2);                  
                        } else {
                            assert (node2 instanceof SinkNode);
                            forwardLinksSinkRM.add(index1);
                        }
                    } else if (node1 instanceof SourceNode) {
                        assert (node2 instanceof LeftNode);
                        forwardLinksSourceRM.add(index2);
                    } else {
                        //node1 is a sink node
                        assert (node2 instanceof RightNode);
                        backwardLinksSinkRM.add(index2);
                    }
                }
                treeIdx++;
            }
        }
        
    }
  
    /**
     * @return the forwardLinksRM
     */
    public Map<Integer, Set<Integer>> getForwardLinksRM() {
        return forwardLinksRM;
    }

    /**
     * key is Right (==Y) index node, and value is
     * Left (==X) index node.
     * @return the backwardLinksRM
     */
    public Map<Integer, Set<Integer>> getBackwardLinksRM() {
        return backwardLinksRM;
    }

    /**
     * @return the sourceNode
     */
    public int getSourceNode() {
        return sourceNode;
    }

    /**
     * @return the sinkNode
     */
    public int getSinkNode() {
        return sinkNode;
    }
    
    /**
     * the arcs that are initialized from being connected
     * to matched nodes, that is Left nodes that
     * are in "saturated" arcs.
     * @return the forwardLinksSourceRM
     */
    public Set<Integer> getForwardLinksSourceRM() {
        return forwardLinksSourceRM;
    }

    /**
     * @return the backwardLinksSourceRM
     */
    public Set<Integer> getBackwardLinksSourceRM() {
        return backwardLinksSourceRM;
    }

    /**
     * the arcs that are initialized from being connected
     * to matched nodes, that is Right nodes that
     * are in "saturated" arcs.
     * @return the forwardLinksSinkRM
     */
    public Set<Integer> getForwardLinksSinkRM() {
        return forwardLinksSinkRM;
    }

    /**
     * @return the backwardLinksSinkRM
     */
    public Set<Integer> getBackwardLinksSinkRM() {
        return backwardLinksSinkRM;
    }

}
