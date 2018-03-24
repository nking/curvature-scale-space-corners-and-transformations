package algorithms.bipartite;

import algorithms.util.PairInt;
import gnu.trove.iterator.TObjectIntIterator;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;
import java.util.Arrays;
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
    
    public void test00() {
        
        log.info("test00");
        
        Graph g = getTestGraph00();
        
        GraphWithoutWeights g2 = new GraphWithoutWeights(g);
        
        HopcroftKarp hk = new HopcroftKarp();
        int[] matched = hk.hopcroftKarpV0(g2);
        log.info("hk matched=" + Arrays.toString(matched));

        MinCostUnbalancedAssignment bipartite = 
            new MinCostUnbalancedAssignment();
        TIntIntMap m = new TIntIntHashMap();
        ResidualDigraph rM = new ResidualDigraph(g, m);
        m = bipartite.hopcroftKarp(g, 3);
        
        assertEquals(3, m.size());             
    }

    public void test0() {
       
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
        
       TIntIntMap m = bipartite.flowAssign(g);
        
        assertTrue(3 <= m.size());
        assertEquals(m.get(1), 2);
        assertEquals(m.get(3), 0);
        assertEquals(m.get(4), 3);
    }
    
    public void test1() throws Exception {
    
        log.info("test1");
                
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        //seed = 1466321947621L;
        sr.setSeed(seed);
        log.info("SEED=" + seed);
        
        for (int size = 10; size <= 100; size *= 10) {
            for (int maxCost = 10; maxCost <= 10000; maxCost *= 10) {
                    
                log.info("size=" + size + " maxCost=" + maxCost);
                
                Graph g = getTestGraph1(size, maxCost, sr);
       
                long t0 = System.currentTimeMillis();
                
                MinCostUnbalancedAssignment bipartite = 
                    new MinCostUnbalancedAssignment();
                
                TIntIntMap m = bipartite.flowAssign(g);
                
                long t1 = System.currentTimeMillis();
                long tSec = (t1 - t0);
                System.out.println(tSec + " msec for flowAssign");

                
                log.info("size=" + size + " scale=" + maxCost + 
                    " m.size=" + m.size());
                
                assertEquals(size, m.size());
                
                for (int i = 0; i < size; ++i) {
                    assertEquals((size - 1 - i), m.get(i));
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
                    assertEquals((size - 1 - idx1), idx2);
                }
            }
        }
    }

    public void test2()  {
        
        try {
        log.info("test2");
        
        int size = 10;
        int scale = 100;//1000000;
        
        log.info("size=" + size + " scale=" + scale);

        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        //seed = 10070693215754L;
        sr.setSeed(seed);
        log.info("seed=" + seed);
        
        for (int nTest = 0; nTest < 10; ++nTest) {
            
            log.info("nTest=" + nTest);
            
            Graph g = getTestGraph2(size, scale, sr);

            MinCostUnbalancedAssignment bipartite = 
                new MinCostUnbalancedAssignment();

            TIntIntMap m = bipartite.flowAssign(g);

            log.info("size=" + size + " scale=" + scale + 
                " m.size=" + m.size());

            assertEquals(size, m.size());

            for (int i = 0; i < size; ++i) {
                assertEquals((size - 1 - i), m.get(i));
            }

            assertNotNull(bipartite.getFinalFlowNetwork());
        }
        
        } catch(Throwable t) {
            t.printStackTrace();
            fail(t.getMessage());  
        }
    }

    public void test3() {
        
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
        
        TIntIntMap m = bipartite.flowAssign(g);

        log.info("size=" + size + " hk size=" + nExpected 
            + " m.size=" + m.size());

        assertEquals(nExpected, m.size());
        assertNotNull(bipartite.getFinalFlowNetwork());       
    
        } catch(Throwable t) {
            //t.printStackTrace();
        }
    }
    
    private Graph getTestGraph00() {
        
        TObjectIntMap<PairInt> weights 
            = new TObjectIntHashMap<PairInt>();
        
        /*
        0  1   2  L
        
        0  1   2  R
        */
        weights.put(new PairInt(0, 0), 1);
        weights.put(new PairInt(0, 1), 2);
        weights.put(new PairInt(1, 0), 2);
        weights.put(new PairInt(1, 1), 1);
        weights.put(new PairInt(1, 2), 1);
        weights.put(new PairInt(2, 1), 2);
        weights.put(new PairInt(2, 2), 1);
        
        Graph g = new Graph(3, 3, weights, true);

        return g;
    }
    
    private Graph getTestGraph0() {
        
        TObjectIntMap<PairInt> weights 
            = new TObjectIntHashMap<PairInt>();
        
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
                
        weights.put(new PairInt(0, 0), 2);
        weights.put(new PairInt(0, 2), 3);
        
        weights.put(new PairInt(1, 1), 2);
        weights.put(new PairInt(1, 2), 1);
        
        weights.put(new PairInt(2, 3), 2);
        
        weights.put(new PairInt(3, 0), 1);
        weights.put(new PairInt(3, 1), 2);
        weights.put(new PairInt(3, 4), 3);
    
        weights.put(new PairInt(4, 3), 1);
        weights.put(new PairInt(4, 4), 2);
        
        weights.put(new PairInt(5, 3), 2);
        
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
        
        TObjectIntMap<PairInt> weights 
            = new TObjectIntHashMap<PairInt>();
        
        int minCostUpper = maxCost/10;
        if (minCostUpper < 2) {
            minCostUpper = 2;
        }
        
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                int cost;
                if (j == (size - 1 - i)) {
                    cost = sr.nextInt(minCostUpper - 1) + 1;
                } else {
                    cost = minCostUpper +
                        sr.nextInt(maxCost - minCostUpper);
                }
                weights.put(new PairInt(i, j), cost);
            }
        }

        // add an extra edge with same answer as best answer for the node
        //  0  --> 3 --> 9
        //  0  --> 3 --> 8
        //  1  --> 3 --> 8
        PairInt p = new PairInt(0, size-1);
        int c = weights.get(p);
        PairInt p1 = new PairInt(1, size-2);
        int c1 = weights.get(p1);
        PairInt p2 = new PairInt(2, size-3);
        int c2 = weights.get(p2);
        if ((c <= c1) && (c <= c2)) {
            weights.put(p1, c);
            weights.put(p2, c);
        } else if ((c1 <= c) && (c1 <= c2)) {
            weights.put(p, c1);
            weights.put(p2, c1);
        } else if ((c2 <= c) && (c2 <= c1)) {
            weights.put(p1, c2);
            weights.put(p2, c2);
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

        TObjectIntIterator<PairInt> iter = g.getEdgeWeights().iterator();
        
        for (int i = g.getEdgeWeights().size(); i-- > 0;) {
            iter.advance();
            PairInt p = iter.key();
            cost[p.getX()][p.getY()] = iter.value();
        }
        
        return cost;
    }
    
    private Graph getTestGraph2(int size, int maxCost, SecureRandom sr) throws 
            NoSuchAlgorithmException {
        
        /*
        - graph of size n for both sets
        - one edge only for each vertex
        - a random cost for the edges
        - graphs of size powers of 10, from 10 to 100,000?
        */
        
        TObjectIntMap<PairInt> weights 
            = new TObjectIntHashMap<PairInt>();
                      
        for (int i = 0; i < size; ++i) {
            int cost = sr.nextInt(maxCost - 1) + 1;
            weights.put(new PairInt(i, (size - 1 - i)), cost);
        }
        
        Graph g = new Graph(size, size, weights, true);

        return g;
    }
    
    private Graph getTestGraph3(SecureRandom sr, int size, int maxCost) 
            throws NoSuchAlgorithmException {
        
        /*
        - graph of size n for both sets
        - random number of edges for each
        - all edges have same cost
        --> expecting same results as hopcroft-karp
        */
        
        TObjectIntMap<PairInt> weights 
            = new TObjectIntHashMap<PairInt>();
        
        int cost = sr.nextInt(maxCost - 1) + 1;
        
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                if (sr.nextBoolean()) {
                    continue;
                }
                weights.put(new PairInt(i, j), cost);
            }
        }
        
        Graph g = new Graph(size, size, weights, true);

        return g;
    }
}
