package algorithms.bipartite;

import algorithms.imageProcessing.DoubleLinkedCircularList;
import algorithms.imageProcessing.Heap;
import algorithms.imageProcessing.HeapNode;
import algorithms.util.PairInt;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

/**
 * A solver for the min-cost, unbalanced, bipartite
 * assignment problem that also uses weight scaling
 * to solve the perfect and imperfect
 * assignment problems, but not incremental
 * with a runtime complexity of
 <pre>
 O(m * sqrt(n) * log(n * C)) 
 where m is the number of edges (a.k.a. arcs) in the graph,
 n is the maximum number of nodes in the two graphs to be matched,
 C is a constant greater than or equal to the maximum 
 edge weight and is greater than or equal to 1.

 The weight scaling is what allows the algorithm to achieve
 much better performance with respect to s
 where s is the number of matched nodes.
 
 The code below follows the paper of
 "On Minimum-Cost Assignments in Unbalanced Bipartite Graphs"
 by Ramshaw and Tarjan, 2012
 for their algorithm, FlowAssign, Refine, and other methods
 used in the paper.
 </pre>
 * @author nichole
 */
public class MinCostUnbalancedAssignment {

    public static class Graph {
                
        /**
         * left (==X) vertices in graph
         */
        Set<Integer> leftG = new HashSet<Integer>();
        
        /**
         * right (==Y) vertices in graph
         */
        Set<Integer> rightG = new HashSet<Integer>();
        
        /**
         * map of edge weights with key = pairint of left
         * index and right index and value being the edge weight.
         */
        Map<PairInt, Integer> edgeWeights = 
            new HashMap<PairInt, Integer>();
    }
    
    public static class FlowNetwork {
                
        int sourceNode = -1;
        int sinkNode = -1;
        
        /**
         * nodes that correspond to left (==X) vertices in G
         */
        Set<Integer> leftNG = new HashSet<Integer>();
        
        /**
         * nodes that correspond to right (==Y) vertices in G
         */
        Set<Integer> rightNG = new HashSet<Integer>();
        
        // forward arcs only in the flow graph
        Map<Integer, Set<Integer>> forwardArcs = 
            new HashMap<Integer, Set<Integer>>();
        
        /**
         * per unit flow costs are the edges given in G.
         * key = index in left (a.k.a. X) and index in
         * right (a.k.a. Y),
         * value = cost of the arc.
         */
        Map<PairInt, Integer> c = 
            new HashMap<PairInt, Integer>();
        
        /**
         * the flow on the arc.
         * key = index in left (a.k.a. X) and index in
         * right (a.k.a. Y),
         * value = number from 0 to 1 inclusive.
         * the value 0 is "idle" and corresponds to an
         * unmatched link in the residual graph.
         * the value 1 is "saturated" and corresponds to
         * a matched link in the residual graph.
         * 
         */
        Map<PairInt, Float> f = 
            new HashMap<PairInt, Float>(); 
        
        /*
        an alternating path is in the residual digraph, and 
           alternates between forward and backward.
        augmenting paths are in the residual digraph.
           an augmenting path starts at a maiden and ends at
           a bachelor (therefore, first and last links are
           forward links).  the links from Y to X are already
           matched (married).
       - an augmenting path takes forward steps to add a new 
         edge to the matched AND takes backward steps to 
         remove an edge from the matched.
       - an augmenting path is "tight" when all of the edges
         that underlie its links are tight, that is, have 
         zero net cost.
       - c(P) is the cost of augmenting path P.
         it's the sum of costs of the arcs of the path in 
         the residual network.
       - Let x_i_j = 1   if i is assigned to j,  else =0
         to augment along path p is to replace x_i_j by 1 for
         each arc (i,j) in P directed from X to Y, and to
         replace each x_i_j by 0 for each arc (j, i) in P
         directed from Y to X.
         - suppose x' is the matching obtained from x after
           augmenting along path P.
           the cost of x' is cx' = cx + c(P)
        
        
        Definition 2-7 on page 15 defines a proper pseudoflow
        of f and cp.
        */
    }
    
    public static class ResidualDigraph {
                
        /**
         * nodes that correspond to left (==X) vertices in G
         */
        Set<Integer> leftRM = new HashSet<Integer>();
        
        /**
         * nodes that correspond to right (==Y) vertices in G
         */
        Set<Integer> rightRM = new HashSet<Integer>();
        
        /**
         * links X to Y (that is, left to right).  
         * These are "idle" arcs, f=0, in the flow network N_G.
         * They correspond to "unmarried" in the matched
         * M graph.
         */
        Map<Integer, Set<Integer>> forwardLinksRM = 
            new HashMap<Integer, Set<Integer>>();
        
        /**
         * links Y to X (that is, right to left).  
         * These are "saturated" arcs, f=1, in the flow network N_G.
         * They correspond to "married" in the matched
         * M graph.
         */
        Map<Integer, Integer> backwardLinksRM = 
            new HashMap<Integer, Integer>();
    }
    
    /**
     * class specializing a fibonacci heap node to identify
     * a left node, a.k.a. X node
     */
    private class LeftNode extends HeapNode {
    }
    /**
     * class specializing a fibonacci heap node to identify
     * a right node, a.k.a. Y node
     */
    private class RightNode extends HeapNode {
    }
    
    /*
    input to FlowAssign is a bipartite graph G with integral 
    edge weights and a target size t. 
    
    s := min(t, ν(G))
    initialization uses Hopcroft-Karp and ignores the edge weights
    to find any matching of size t.
    If t exceeds nu(G) a warning is logged and a matching
    of size s is returned.
    (a matching of size s in the graph G corresponds to an 
    integral flow f of value |f| = s in the flow network N_G)
    
    
    input to Refine is a flow f of value s and prices p
    that together make all arcs (qε)-proper, where q is an integer 
    parameter. Refine builds a new flow f′, also of value s, 
    and prices p′ that together make all arcs ε-proper. 
    
    Note:  In FlowAssign, saturated edges will be kept proper, 
       but not necessarily tight. So, for a saturated edge 
       (x, y), we will know only that cp(x, y) ≤ 0, and we 
       will need to define lp(y ⇒ x) := −cp(x, y), to keep 
       our link lengths nonnegative. That is, for an 
       alternating path A, we will define
       lp(A):= summation over arcs X->Y (cp(X,Y)) -
                  summation over arcs Y->X (cp(X,Y))
    
    
    FlowAssign (G, t)
      (M, s) := HopcroftKarp(G, t);
      convert M into an integral flow f on N_G with |f| = s; 
      set ε := ε and, 
      for all nodes v in N_G, 
        set pd(v) := 0; 
      while ε > ε do
        ε := ε/q;
        Refine(f,p,ε); od;
      round prices to integers that make all arcs proper; 
    
    
    Refine(f,p,ε)
      S := {the s women who are matched in f};
      D := {the s men who are matched in f};
      convert the s bipartite arcs that are saturated in 
        f to idle; 
      raise the prices p, as in Figure 7.4, to make all arcs 
        ε-proper; 
      int h := s;
      while h > 0 do
        build a shortest-path forest from the current surpluses S, 
          stopping when a current deficit in D is reached;
        raise prices at forest nodes by multiples of ε, 
          shortening the discovered augmenting path to length 0;
        find a maximal set P of length-0 augmenting paths 
          that are compatible, as defined in Section 8.3;
        augment f along each of the paths in P in turn, thereby 
          reducing |S| = |D| = h by |P|;
      od;
    
    HopcroftKarp (G)
      set M to the empty matching; 
      do
        find a maximal set P of vertex-disjoint augmenting paths, 
          for M, all of the minimum possible length;
        if |P| = 0 then announce(ν(G) = |M|); 
          return; 
        fi; 
        for P in P do
          augment M along P;
          announce(M is a matching); od;
      od;    
    */
    
    protected Map<Integer, Integer> hopcroftKarp(Graph g) {
        
        Map<Integer, Integer> m = new HashMap<Integer, Integer>();
        
        ResidualDigraph rM = createResidualGraph(g, m);
        
        DoubleLinkedCircularList[] augmentingPaths =
            buildForest(rM);
        
        if (allAreEmpty(augmentingPaths)) {
            return m;
        }
    
        for (int i = 0; i < augmentingPaths.length; ++i) {
            
            DoubleLinkedCircularList path = augmentingPaths[i];
            
            if (path == null) {
                continue;
            }
            
            //augment M along path;
            
            //announce(M is a matching)
            
        }
        
        return m;
    }
    
    /*
    NOTES
    ------
    input 
       G = (X, Y, E)
       set<integer> x, set<integer> y
    output
       map<integer, integer> m0
    
    G has vertices and edges 
       (and edge is usually written (X,Y)),
    N_G, the flow network, has nodes and arcs 
       (and arc is written X → Y, and all arcs go forward.) 
    residual digraph has nodes and links, 
       (and link is written X ⇒ Y, and some links 
       go forward while others go backward.)
       - Given a matching M in G, the residual digraph 
         R_M has nodes and links that correspond precisely 
         to the vertices and edges of G, 
         but each saturated edge (X, Y) 
           becomes a backward-directed link Y ⇒ X 
           in the residual digraph R_M, 
         while each idle edge (X, Y) becomes a 
           forward-directed link X ⇒ Y. 
         That is, idle arcs (f=0) become left-to-right links in R_M, 
         while saturated arcs (f=1) become right-to-left links. 
         Note that the residual digraph R_M depends only on the 
         matching M, not on the prices p.    
    
    
       - length or distance in the residual digraph
         is a term usable when the link has length >= 0
         and the underlying edges are "tight" (have 0 net cost).
   **    - to find the shortest possible augmenting path,
           can use a variant of Dijkstra's algorithm.   
           Note that if X->Y is a forward link R_M, corresponding
           to an idle edge (x,y), 
              define the length of that link 
                to be the net cost of the edge (lp(x ⇒ y) = cp(X,Y)).
              define the length of any backward link y ⇒ x 
                to be zero: lp(y ⇒ x) := 0. 
              define the length of an alternating path to be 
                 the sum of the lengths of its links.
              Note that the term link-count is in contrast, the
                number of links regardless of edge cost.
    
    A pseudoflow is a flux in which the flow f(X,Y) along each 
      arc X → Y is nonnegative and satisfies the unit-capacity 
      constraint: 0 ≤ f (X, Y) ≤ 1. 
      - an arc is idle in f if it's a pseudoflow w/ f(X,Y) = 0 
      - an arc is saturated if f(X,Y) = 1
      - else arc is fractional flow with 0 < f(X,Y) < 1  
The value of a flow f, denoted |f|, is the total 
    flow out of the source, which is also the total 
    flow into the sink – and, for that matter, 
    the total flow over all of the bipartite arcs.
A flux f is integral when, for every arc X → Y, 
    the flow f(X,Y) over that arc is an integer. 
    We will typically be dealing with fluxes, 
    pseudoflows, and flows that are integral. 
    If a pseudoflow on some flow network is integral, 
    then every arc is either idle or saturated — 
    no arcs have fractional flow.
Given any flux f in the flow network N_G, 
    we define the cost of that flux to be the sum 
    of the costs of its arcs:
        c(f) = summation over arcs in N_G ( f(X,Y) * c(X,Y) )
Matchings in G are integral flows in N_G
     and number of matchings in M = |f|
     and c(M) = c(f)
     Thus, a min-cost matching of some size s corresponds to 
     a min-cost integral flow of value s.
-------    
    Robert B. Dial. Algorithm 360: 
       Shortest path forest with topological ordering. 
       Communications of the ACM 12 (1969) pp. 632–633.
    */
    
    /**
     * find augmenting paths of minimal length.
     * by building the shortest path forest.
     * 
     * implementing section 3.4, pg 22 of Ramshaw and Tarjan 2012.
     * 
     * TODO: edit to accept lambda, the largest link length
     * in the forest.
     * 
     * @param rM the residual digraph built from the matching 
     * graph M of graph G.
     * 
     * @return a forest of trees whose forest array indexes 
     * are the path lengths and whose items hold at the root, the
     * remaining maidens, that is unmatched left nodes (a.k.a. G's
     * X nodes).
     */
    protected DoubleLinkedCircularList[] 
        buildForest(final ResidualDigraph rM) {
        
        //TODO: revisit this.
        int lambda = Math.min(rM.leftRM.size(), 
            rM.rightRM.size())/2;
        
        DoubleLinkedCircularList[] forest 
            = new DoubleLinkedCircularList[lambda];
        
        Heap heap = new Heap();
        
        /*NOTE:
        need to run tests on this section.
        
        upon first use as hopcroft karp is initialized
        with an empty matching graph.
        there are no backward links in the residual digraph
        for an empty matching graph, and so there are no
        alternating paths, just maiden nodes from Y.
        */
        
        // init all nodes to inf length
        Map<Integer, LeftNode> leftNodes = new HashMap<Integer, LeftNode>();
        Map<Integer, RightNode> rightNodes = new HashMap<Integer, RightNode>();
        for (Integer rNode : rM.rightRM) {
            RightNode node = new RightNode();
            node.setKey(Long.MAX_VALUE);
            node.setData(rNode);
            rightNodes.put(rNode, node);
        }
        for (Integer lNode : rM.leftRM) {
            LeftNode node = new LeftNode();
            node.setKey(Long.MAX_VALUE);
            node.setData(lNode);
            leftNodes.put(lNode, node);
        }
        
        // for all maidens, that is, the keys in forwardLinksRM,
        // set key to 0 and ScanAndAdd(index)
        for (Integer lNode : rM.forwardLinksRM.keySet()) {
            LeftNode node = leftNodes.get(lNode);
            node.setKey(0);
            scanAndAdd(heap, forest, rM, rightNodes, node);
        }
        
        // at this point, the maidens are all in index 0 of
        // the forest trees.
        
        while (!heap.isEmpty()) {
                        
            // in the heap are men not in the forest who are in
            // an alternating path from a maiden.
            // the key is the length of shortest path so far
            HeapNode y = heap.extractMin();
            assert(y instanceof RightNode);
            assert(y.getData() != null);
            
            addToForest(forest, y);
        
            /*
            if y is married then
                x := wife of y;
                set l(x) := l(y) and 
                ScanAndAdd(x); 
             else
                exit(bachelor β := y reached); 
             fi;
            */
            Integer yIndex = (Integer)(y.getData());
            // the married nodes are the keys in the backward
            // links of the residual digraph
            Integer xIndex = rM.backwardLinksRM.get(yIndex);
            if (xIndex != null) {
                
                LeftNode xNode = leftNodes.get(xIndex);
                RightNode yNode = rightNodes.get(yIndex);
                
                xNode.setKey(yNode.getKey());
                // modifiy heap key if xNode is in heap?
                
                scanAndAdd(heap, forest, rM, rightNodes, xNode);
                
            } else {    
                //exit(bachelor β := y reached);
                break;
            }
        }
        
        return forest;
    }

    /**
    visit alternating paths and update the heap for the best
    path lengths in the right (matched) nodes and add the
    maiden nodes to the forest.
    
    given a maiden (left) node, search forward links to 
    right nodes in alternating paths
    to insert right nodes into the heap
    or update their heap keys for shorter paths.
    
     * @param heap
     * @param forest
     * @param rM
     * @param yNodes
     * @param xNode 
     */
    private void scanAndAdd(
        Heap heap, DoubleLinkedCircularList[] forest,
        ResidualDigraph rM, Map<Integer, RightNode> yNodes, 
        LeftNode xNode) {
        
        long lX = xNode.getKey();
        Integer x = (Integer)(xNode.getData());
        
        Set<Integer> forwardLinks = rM.forwardLinksRM.get(x);
        for (Integer y : forwardLinks) {
            
            // because the paths need to be alternating, check
            // that the y, that is the right node, is in the
            // residual graph backward links indicating it is matched.
            if (!rM.backwardLinksRM.containsKey(y)) {
                continue;
            }
            
            // by definition, the forward link in residual digraph of matches,
            // is "idle", that is f=0
            
            //link length = net cost of the edge 
            //    (lp(x ⇒ y) = cp(X,Y))
            //link length of backward link Y->X, is 0, lp(y ⇒ x) := 0.
            
            // l(x) and l(y) are the keys in the heap node
            
            RightNode yNode = yNodes.get(y);
            long lY = yNode.getKey();
            assert(((Integer)(yNode.getData())).intValue() 
                == y.intValue());
            
            //L := l(x) + lp(x ⇒ y); 
            // for the bipartite digraph logic and idle edge, lp(x ⇒ y)=0
            
            long ell = lX;// + 0
            long lOld = lY;
            if (ell < lOld) {
                lY = ell;
                if (lOld == Long.MAX_VALUE) {
                    yNode.setKey(lY);
                    heap.insert(yNode);
                } else {
                    heap.decreaseKey(yNode, lY);
                }
            }
        }
            
        //add x to the forest;
        addToForest(forest, xNode);
    }
    
    private void addToForest(DoubleLinkedCircularList[] forest,
        HeapNode node) {
        
        /*
        Section 8.1, pg 57
        We ignore any paths we find whose lengths exceed forest.length. 
        We also maintain an integer B, which stores the 
        value l(v) for the node v that was most recently 
        added to the forest. 
        We add nodes v to the forest in nondecreasing 
        order of l(v), so B never decreases. 
        To implement insert(v,k), we add v to the list Q[k].
        */
        
        if (node.getKey() < forest.length) {
            int k = (int) node.getKey();
            DoubleLinkedCircularList list = forest[k];
            if (list == null) {
                list = new DoubleLinkedCircularList();
                forest[k] = list;
            }
            list.insert(node);
        }
    }

    private boolean allAreEmpty(DoubleLinkedCircularList[] 
        augmentingPaths) {
        
        for (int i = 0; i < augmentingPaths.length; ++i) {
            DoubleLinkedCircularList path = augmentingPaths[i];
            if (path != null) {
                return false;
            }
        }
        
        return true;
    }
    
    protected ResidualDigraph createResidualGraph(Graph g, 
        Map<Integer, Integer> m) {
                
        ResidualDigraph rM = new ResidualDigraph();
        
        for (Integer x : g.leftG) {
            rM.leftRM.add(x);
        }
        
        for (Integer y : g.rightG) {
            rM.rightRM.add(y);
        }
        
        for (Entry<PairInt, Integer> entry : g.edgeWeights.entrySet()) {
            
            PairInt p = entry.getKey();
            
            // presumably, edge weight of 0 is "not connected"
            if (entry.getValue().intValue() == 0) {
                continue;
            }
            
            Integer x = Integer.valueOf(p.getX());
            Integer y = Integer.valueOf(p.getY());
            
            Set<Integer> ys = rM.forwardLinksRM.get(x);
            if (ys == null) {
                ys = new HashSet<Integer>();
                rM.forwardLinksRM.put(x, ys);
            }
            
            ys.add(y);
        }
        
        for (Entry<Integer, Integer> entry : m.entrySet()) {
            Integer x = entry.getKey();
            Integer y = entry.getValue();
            
            if (rM.forwardLinksRM.containsKey(x)) {
                rM.forwardLinksRM.get(x).remove(y);
            }
            
            rM.backwardLinksRM.put(y, x);
        }

        return rM;
    }

}
