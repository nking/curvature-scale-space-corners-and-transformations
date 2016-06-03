package algorithms.bipartite;

import algorithms.util.PairInt;
import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;
import java.util.HashMap;
import java.util.Map;
import java.util.logging.Logger;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class MinCostUnbalancedAssignment3Test extends TestCase {
    
    private Logger log = 
        Logger.getLogger(this.getClass().getName());
    
    public MinCostUnbalancedAssignment3Test() {
    }
    
    public void est00() {
        
        Graph g = getTestGraph00();
        
        MinCostUnbalancedAssignment bipartite = 
            new MinCostUnbalancedAssignment();
        
        Map<Integer, Integer> m = bipartite.hopcroftKarp(g, 3);
        
        assertEquals(3, m.size());
             
    }
    
    public void est0() {
        
        // test graphs on pg 49
        Graph g = getTestGraph0();
        
        MinCostUnbalancedAssignment bipartite = 
            new MinCostUnbalancedAssignment();
        
        Map<Integer, Integer> m = bipartite.flowAssign(g);
        
        assertEquals(3, m.size());
        assertTrue(m.get(Integer.valueOf(1))
            .equals(Integer.valueOf(2)));
        assertTrue(m.get(Integer.valueOf(3))
            .equals(Integer.valueOf(0)));
        assertTrue(m.get(Integer.valueOf(4))
            .equals(Integer.valueOf(3)));
    }
    
    public void est1() throws Exception {
        // size=100, scale=10 shows error in hopcroft-karp
        for (int size = 10; size <= 10; size *= 10) {
            for (int scale = 10; scale <= 100000; scale *= 10) {
                
                log.info("size=" + size + " scale=" + scale);
                
                Graph g = getTestGraph1(size, scale);
                
                MinCostUnbalancedAssignment bipartite = 
                    new MinCostUnbalancedAssignment();
                
                Map<Integer, Integer> m = bipartite.flowAssign(g);
                
                log.info("size=" + size + " scale=" + scale + 
                    " m.size=" + m.size());
                
                assertEquals(size, m.size());
                
                for (int i = 0; i < size; ++i) {
                    assertEquals(i, m.get(Integer.valueOf(i)).intValue());
                }
                
                assertNotNull(bipartite.getFinalFlowNetwork());
            }
        }
    }

    public void est2() throws Exception {
        
        int size = 10;
        int scale = 10000;//1000000;
        
        log.info("size=" + size + " scale=" + scale);

        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.nanoTime();
        sr.setSeed(seed);
        log.info("seed=" + seed);
        
        for (int nTest = 0; nTest < 10; ++nTest) {
            
            log.info("nTest=" + nTest);
            
            Graph g = getTestGraph2(size, scale, sr);

            MinCostUnbalancedAssignment bipartite = 
                new MinCostUnbalancedAssignment();

            Map<Integer, Integer> m = bipartite.flowAssign(g);

            log.info("size=" + size + " scale=" + scale + 
                " m.size=" + m.size());

            assertEquals(size, m.size());

            for (int i = 0; i < size; ++i) {
                assertEquals(i, m.get(Integer.valueOf(i)).intValue());
            }

            assertNotNull(bipartite.getFinalFlowNetwork());
        }
    }

    public void test3() throws Exception {
        
        int size = 10;
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        seed = 1464988291524L;
        sr.setSeed(seed);
        log.info("SEED=" + seed);
        
        Graph g = getTestGraph3(sr, size);

        MinCostUnbalancedAssignment bipartite = 
            new MinCostUnbalancedAssignment();

        Map<Integer, Integer> m = bipartite.flowAssign(g);

        log.info("size=" + size + " m.size=" + m.size());

        assertEquals(size, m.size());

        assertNotNull(bipartite.getFinalFlowNetwork());       
    }

    private Graph getTestGraph00() {
        
        Map<PairInt, Integer> weights 
            = new HashMap<PairInt, Integer>();
        
        /*
        0  1   2  L
        
        0  1   2  R
        */
        weights.put(new PairInt(0, 0), Integer.valueOf(1));
        weights.put(new PairInt(0, 1), Integer.valueOf(2));
        weights.put(new PairInt(1, 0), Integer.valueOf(2));
        weights.put(new PairInt(1, 1), Integer.valueOf(1));
        weights.put(new PairInt(1, 2), Integer.valueOf(1));
        weights.put(new PairInt(2, 1), Integer.valueOf(2));
        weights.put(new PairInt(2, 2), Integer.valueOf(1));
        
        Graph g = new Graph(3, 3, weights, true);

        return g;
    }
    
    private Graph getTestGraph0() {
        
        Map<PairInt, Integer> weights 
            = new HashMap<PairInt, Integer>();
        
        /*
        for best cost, w/o regard to maximizing nMatches:
        1  2
        3  0
        4  3
        
        unmatched:
        0,2,5     1,4
        ------------------
        
        for maximizing number of matches, then min cost:
        results are same as Hopcroft-Karp for this example
        which is a matching size of 5:
        5  3   or  2 3 is equivalent
        4  4
        3  0
        1  1
        0  2
        
        -----
        so far, the full algorithm progresses to this, but can find
        no further augmenting paths:
          0  1
          1  1
          3  4
          4  3  <--- this one is min-cost, optimal for the pair
                     but it prevents 5 3 
        */
                
        weights.put(new PairInt(0, 0), Integer.valueOf(2));
        weights.put(new PairInt(0, 2), Integer.valueOf(3));
        
        weights.put(new PairInt(1, 1), Integer.valueOf(2));
        weights.put(new PairInt(1, 2), Integer.valueOf(1));
        
        weights.put(new PairInt(2, 3), Integer.valueOf(2));
        
        weights.put(new PairInt(3, 0), Integer.valueOf(1));
        weights.put(new PairInt(3, 1), Integer.valueOf(2));
        weights.put(new PairInt(3, 4), Integer.valueOf(3));
    
        weights.put(new PairInt(4, 3), Integer.valueOf(1));
        weights.put(new PairInt(4, 4), Integer.valueOf(2));
        
        weights.put(new PairInt(5, 3), Integer.valueOf(2));
        
        Graph g = new Graph(6, 5, weights, true);

        return g;
    }
    
    private Graph getTestGraph1(int size, int maxCost) throws NoSuchAlgorithmException {
        
        /*
        - graph of size n for both sets
        - highly connected, that is each vertex connected
          to all right vertexes
        - one min-cost best cost for a vertex
          and all other edges have larger costs than
          those to make the graph easier to write.
        - want the cost to use a multiple so can
          test that correct solution is obtained
          no matter the magnitude of costs
          (the weight scaling needs such exploration for
          the internal values of q and eps)
        - make graphs of size of power of 10 from 10 to 
            1 million.
        */
        
        Map<PairInt, Integer> weights 
            = new HashMap<PairInt, Integer>();
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        
        int minCostUpper = maxCost/10;
        
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                int cost;
                if (j == i) {
                    cost = sr.nextInt(minCostUpper);
                } else {
                    cost = minCostUpper +
                        sr.nextInt(maxCost - minCostUpper);
                }
                weights.put(new PairInt(i, j), Integer.valueOf(cost));
            }
        }
        
        Graph g = new Graph(size, size, weights, true);

        return g;
    }
    
    private Graph getTestGraph2(int size, int maxCost, 
        SecureRandom sr) throws NoSuchAlgorithmException {
        
        /*
        - graph of size n for both sets
        - one edge only for each vertex
        - a random cost for the edges
        - graphs of size powers of 10, from 10 to 100,000?
        */
        
        Map<PairInt, Integer> weights 
            = new HashMap<PairInt, Integer>();
                        
        for (int i = 0; i < size; ++i) {
            int cost = sr.nextInt(maxCost);
            weights.put(new PairInt(i, i), Integer.valueOf(cost));
        }
        
        Graph g = new Graph(size, size, weights, true);

        return g;
    }
    
    private Graph getTestGraph3(SecureRandom sr,
        int size) throws NoSuchAlgorithmException {
        
        /*
        - graph of size n for both sets
        - random number of edges for each
        - all edges have same cost
        --> expecting same results as hopcroft-karp
        */
        
        Map<PairInt, Integer> weights 
            = new HashMap<PairInt, Integer>();
        
        int cost = sr.nextInt(100);
        
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                if (sr.nextBoolean()) {
                    continue;
                }
                weights.put(new PairInt(i, j), Integer.valueOf(cost));
            }
        }
        
        Graph g = new Graph(size, size, weights, true);

        return g;
    }
}
