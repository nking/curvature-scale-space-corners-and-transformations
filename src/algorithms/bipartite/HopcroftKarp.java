package algorithms.bipartite;

import gnu.trove.iterator.TIntIterator;
import gnu.trove.set.TIntSet;
import java.util.Arrays;
import java.util.logging.Logger;

/**
 * finds a maximum matching in a bipartite graph.
 * Note that the graph does not hae edge weights.
 * 
 * useful reading is Cormen et al. "Introduction to Algorithms"
 * pseudocode of Hopcroft-Karp.
 * 
 * Also helpful is 
 * http://en.wikipedia.org/wiki/Hopcroft%E2%80%93Karp_algorithm
 * and
 * http://github.com/indy256/codelibrary/blob/master/java/src/MaxMatchingHopcroftKarp.java     * @param g

 * The runtime complexity is O(sqrt(V) * E).
 * 
 * @author nichole
 */
public class HopcroftKarp {
    
    private Logger log = Logger.getLogger(this.getClass().getName());

    /*
    G is a directed bipartite graph (V, E) where V is composed
    of left L and right R.
    
    a path is a sequence of edges.
    an augmenting path in the matching M starts at an unmatched L
    and ends at an unmatched vertex in R and in between alternates
    between matched and unmatched edges, (members of E and E-M).
    an augmenting path can be composed of only two vertices and the
    edge between them.
    
    The shortest augmenting path has the fewest number of edges in it.
    
    And the symmetric difference of 2 sets is the points that are 
    not in the intersection, but are in the union of them.
    For example, sym diff of {1, 2, 3} and {3, 4} is {1, 2, 4}.
    http://en.wikipedia.org/wiki/Symmetric_difference
    
    if M is a matching within G and P is an augmenting path,
    then the summetric difference of M with P is a matching
    of size |M| + 1.
    
    input G
    
    M = 0
    repeat
        let P = {P1, P2, ...Pk} be a maximum set of vertex-disjoint
            shortest augmenting paths with respect to M
        M = the symmetric difference between M and
            (P1 union P2 union ...Pk) 
        until P = 0
    return M
    
    finding the vertex-disjoint shortest paths,
    using the pattern of single BFS followed by DFS, per L vertex
        (see wikipedia).
    
    */
    
    /**
     * implementing the version of Hopcroft-Karp that uses a
     * single round of BFS followed by DFs to find the 
     * shortest augmenting paths, all within a pattern
     * of augmenting the matching M until no new augmenting
     * paths can be found.  the matching M is returned.
     * The code below follows:
      http://github.com/indy256/codelibrary/blob/master/java/src/MaxMatchingHopcroftKarp.java     * @param g
      which uses the unlicense:
      http://github.com/indy256/codelibrary/blob/master/UNLICENSE
     
     * @return matching from perspective int[uIndex] = vIndex
     */
    public int[] hopcroftKarpV0(GraphWithoutWeights g) {
       
        int n1 = g.getNLeft();
        int n2 = g.getNRight();
        
        int[] dist = new int[n1];
        
        int[] match21 = new int[n2];
        Arrays.fill(match21, -1);
        
        // forward matching indexes, opposite mapping of match21
		int[] match12 = new int[n1];
        Arrays.fill(match12, -1);
        
        for (int res = 0; ; ) {
			
            bfs(g, match12, match21, dist);
			
            boolean[] vis = new boolean[n1];
			
            int f = 0;
			
            for (int u = 0; u < n1; ++u) {
                if ((match12[u] == -1) && 
                    dfs(g, vis, match12, match21, dist, u)) {
                    ++f;
                }
            }

            if (f == 0) {
                return match12;
            }
            res += f;
        }        
    }
    
    /**
     * note, this depends upon g
     * @param g
     * @param match12
     * @param match21
     * @param dist 
     */
    private void bfs(GraphWithoutWeights g, int[] match12,
        int[] match21, int[] dist) {
        Arrays.fill(dist, -1);
        int n1 = g.getNLeft();
        int[] Q = new int[n1];
        int sizeQ = 0;
        for (int u = 0; u < n1; ++u) {
            if (match12[u] == -1) {
                Q[sizeQ++] = u;
                dist[u] = 0;
            }
        }
        for (int i = 0; i < sizeQ; i++) {
            int u1 = Q[i];
            TIntSet neighbors = g.getAdjacencyMap().get(u1);
            if (neighbors == null) {
                continue;
            }
            TIntIterator iter = neighbors.iterator();
            while (iter.hasNext()) {
                int vIdx = iter.next();
                log.fine(String.format(
                    "bfs visiting (%d, %d)", u1, vIdx));
                int u2 = match21[vIdx];
                if (u2 > -1 && dist[u2] < 0) {
                    dist[u2] = dist[u1] + 1;
                    Q[sizeQ++] = u2;
                }
            }
        }
    }

    private boolean dfs(GraphWithoutWeights g, boolean[] vis, 
        int[] match12, int[] match21, int[] dist, int u1) {

        vis[u1] = true;
		
        TIntSet neighbors = g.getAdjacencyMap().get(u1);
        if (neighbors != null) {
            TIntIterator iter = neighbors.iterator(); 
            while (iter.hasNext()) {
                int v = iter.next();
                log.fine(String.format(
                    "DFS visiting (%d, %d)", u1, v));
                int u2 = match21[v];
                log.fine(String.format("u2=%d", u2));
                if (u2 < 0 || !vis[u2] && (dist[u2] 
                    == (dist[u1] + 1)) 
                    && dfs(g, vis, match12, match21, dist, u2)) {
                    
                    log.fine(String.format("m[%d]=%d", v, u1));
                    
                    match21[v] = u1;
                    match12[u1] = v;
                    return true;
                }
            }
        }
        
        return false;
    }

}
