package algorithms.bipartite;

import algorithms.util.PairInt;
import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.logging.Logger;
import junit.framework.TestCase;
import thirdparty.HungarianAlgorithm;

/**
 *
 * @author nichole
 */
public class MinCostUnbalancedAssignment3Test extends TestCase {
    
    private Logger log = 
        Logger.getLogger(this.getClass().getName());
    
    public MinCostUnbalancedAssignment3Test() {
    }
    
    public void est0() {
       
        log.info("test0");
        
        // test graphs on pg 49
        Graph g = getTestGraph0();
        
        /*
        for the hopcroft karp portion
        5  3   or  2 3 is equivalent
        4  4
        3  0
        1  1
        0  2
        */
                
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
    
        log.info("test1");
                
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        //seed = 1465192949978L;
        sr.setSeed(seed);
        log.info("SEED=" + seed);
        
        for (int size = 100; size <= 100; size *= 10) {
            for (int maxCost = 10; maxCost <= 10; maxCost *= 10) {
                
                log.info("size=" + size + " maxCost=" + maxCost);
                
                Graph g = getTestGraph1(size, maxCost, sr);
            
                long t0 = System.currentTimeMillis();
                
                MinCostUnbalancedAssignment bipartite = 
                    new MinCostUnbalancedAssignment();
                
                Map<Integer, Integer> m = bipartite.flowAssign(g);
                
                long t1 = System.currentTimeMillis();
                long tSec = (t1 - t0);
                System.out.println(tSec + " msec for flowAssign");

                
                log.info("size=" + size + " scale=" + maxCost + 
                    " m.size=" + m.size());
                
                assertEquals(size, m.size());
                
                for (int i = 0; i < size; ++i) {
                    assertEquals(i, m.get(Integer.valueOf(i)).intValue());
                }
                
                assertNotNull(bipartite.getFinalFlowNetwork());
            
                float[][] matrix = convert(g);
                
                long t2 = System.currentTimeMillis();
                
                HungarianAlgorithm ha = new HungarianAlgorithm();
                int[][] m2 = ha.computeAssignments(matrix);
                
                long t3 = System.currentTimeMillis();
                tSec = (t3 - t2);
                System.out.println(tSec + " msec for hungarian");
                
                log.info("size=" + size + " scale=" + maxCost + 
                    " m.size=" + m.size());
                
                for (int i = 0; i < size; ++i) {
                    int idx1 = m2[i][0];
                    int idx2 = m2[i][1];
                    assertEquals(idx1, idx2);
                }
            }
        }
    }

    public void test2() throws Exception {
        
        log.info("test2");
        
        int size = 10;
        int scale = 100;//1000000;
        
        log.info("size=" + size + " scale=" + scale);

        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        seed = 10070693215754L;
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

    public void est3() {
        
        log.info("test3");
        
        int size = 10;
        int maxCost = 10;
    
        /*
        - graph of size n for both sets
        - random number of edges for each
        - all edges have same cost
        --> expecting same results as hopcroft-karp, that
            is a maximal matching, but the random links
            may create less than maximal possibilities.
            the method is mostly to exercise the
            code to explore it further.
        */
        
        try {
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        //seed = 1464995162443L;
        sr.setSeed(seed);
        log.info("SEED=" + seed);
        
        Graph g = getTestGraph3(sr, size, maxCost);
        
        MinCostUnbalancedAssignment bipartite = 
            new MinCostUnbalancedAssignment();

        HopcroftKarp hk = new HopcroftKarp();
        int[] matched = hk.hopcroftKarpV0(new GraphWithoutWeights(g));
        int nExpected = matched.length;
        log.info("size=" + size + " hk size=" + nExpected);
        
        Map<Integer, Integer> m = bipartite.flowAssign(g);

        log.info("size=" + size + " hk size=" + nExpected 
            + " m.size=" + m.size());

        assertEquals(nExpected, m.size());

        assertNotNull(bipartite.getFinalFlowNetwork());       
    
        } catch(Throwable t) {
            //t.printStackTrace();
        }
    }
    
    public void est00() {
        
        log.info("test00");
        
        Graph g = getTestGraph00();
        
        GraphWithoutWeights g2 = new GraphWithoutWeights(g);
        
        HopcroftKarp hk = new HopcroftKarp();
        int[] matched = hk.hopcroftKarpV0(g2);
        log.info("hk matched=" + Arrays.toString(matched));

        MinCostUnbalancedAssignment bipartite = 
            new MinCostUnbalancedAssignment();
        Map<Integer, Integer> m = new HashMap<Integer, Integer>();
        ResidualDigraph rM = new ResidualDigraph(g, m);
        m = bipartite.hopcroftKarp(g, 3);
        
        
        assertEquals(3, m.size());             
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
        
        /*
        5  3   or  2 3 is equivalent
        4  4
        3  0
        1  1
        0  2
        */
        
        Graph g = new Graph(6, 5, weights, true);

        return g;
    }
    
    private Graph getTestGraph1(int size, int maxCost,
        SecureRandom sr) throws NoSuchAlgorithmException {
        
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
        
        int minCostUpper = maxCost/10;
        if (minCostUpper < 2) {
            minCostUpper = 2;
        }
        
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                int cost;
                if (j == i) {
                    cost = sr.nextInt(minCostUpper - 1) + 1;
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
    
    private float[][] convert(Graph g) {
        
        int n1 = g.getNLeft();
        int n2 = g.getNRight();
        
        float[][] cost = new float[n1][n2];
        for (int i = 0; i < n1; ++i) {
            cost[i] = new float[n2];
            Arrays.fill(cost[i], Float.MAX_VALUE);
        }
        
        for (Entry<PairInt, Integer> entry : 
            g.getEdgeWeights().entrySet()) {
        
            PairInt p = entry.getKey();
            cost[p.getX()][p.getY()] = entry.getValue().intValue();
        }
        
        return cost;
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
            int cost = sr.nextInt(maxCost - 1) + 1;
            weights.put(new PairInt(i, i), Integer.valueOf(cost));
        }
        
        Graph g = new Graph(size, size, weights, true);

        return g;
    }
    
    private Graph getTestGraph3(SecureRandom sr,
        int size, int maxCost) throws NoSuchAlgorithmException {
        
        /*
        - graph of size n for both sets
        - random number of edges for each
        - all edges have same cost
        --> expecting same results as hopcroft-karp
        */
        
        Map<PairInt, Integer> weights 
            = new HashMap<PairInt, Integer>();
        
        int cost = sr.nextInt(maxCost - 1) + 1;
        
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
