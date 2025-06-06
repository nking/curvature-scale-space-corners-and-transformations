package thirdparty.edu.princeton.cs.algs4;

import algorithms.misc.Misc;
import algorithms.util.PairInt;
import algorithms.util.QuadInt;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class QuadTreeInterval2DTest extends TestCase {
    
    public QuadTreeInterval2DTest() {
    }
    
    public void test0() {
        
        /*
          7
          6         & &
          5         & &
          4
          3
            -5      0 1    5
        */
        
        Interval<Double> intervalX1 = new Interval<Double>(0.0, 1.0);
        Interval<Double> intervalY1 = new Interval<Double>(5.0, 6.0);
        Interval2D<Double> box1 
            = new Interval2D<Double>(intervalX1, intervalY1);
        Interval<Double> intervalX2 = new Interval<Double>(-5.0, 5.0);
        Interval<Double> intervalY2 = new Interval<Double>(3.0, 7.0);
        Interval2D<Double> box2 = new Interval2D<Double>(
          intervalX2, intervalY2);
         
        QuadTreeInterval2D<Double, String> st2
            = new QuadTreeInterval2D<Double, String>();
         
        st2.insert(box1, "box1");
        
        List<Interval2D<Double>> results = st2.query2D(box2);
        
        assertEquals(1, results.size());
        
        assertTrue(results.get(0).equals(box1));
    
        // ---- test remove
        st2.remove(box1);
        results = st2.query2D(box2);
        assertEquals(0, results.size());
    }
    
    public void test1() {
        
        /*
          7               *
          6               *
          5               *
          4     |     |
          3     |     |
          2     |_____|      + + +
          1                  + + +
          0
            -5 -4 -3 -2 -1 0 1 2 3 4 5
        */
        Interval<Integer> x1 = new Interval<Integer>(-4, -2);
        Interval<Integer> y1 = new Interval<Integer>(2, 4);
        Interval2D<Integer> box1 = new Interval2D<Integer>(x1, y1);

        Interval<Integer> x2 = new Interval<Integer>(0, 0);
        Interval<Integer> y2 = new Interval<Integer>(5, 7);
        Interval2D<Integer> box2 = new Interval2D<Integer>(x2, y2);
        
        Interval<Integer> x3 = new Interval<Integer>(1, 3);
        Interval<Integer> y3 = new Interval<Integer>(1, 2);
        Interval2D<Integer> box3 = new Interval2D<Integer>(x3, y3);

        Interval<Integer> srchX = new Interval<Integer>(-2, 2);
        Interval<Integer> srchY = new Interval<Integer>(2, 5);
        Interval2D<Integer> srch = new Interval2D<Integer>(srchX, srchY);

        QuadTreeInterval2D<Integer, String> st2
            = new QuadTreeInterval2D<Integer, String>();
         
        st2.insert(box1, "box1");
        st2.insert(box2, "box2");
        st2.insert(box3, "box3");
        
        List<Interval2D<Integer>> results = st2.query2D(srch);
        
        assertEquals(3, results.size());
        
        Set<Interval2D<Integer>> expected = new
            HashSet<Interval2D<Integer>>();
        expected.add(box1);
        expected.add(box2);
        expected.add(box3);
        
        for (Interval2D<Integer> ii : results) {
            assertTrue(expected.remove(ii));
        }
        assertTrue(expected.isEmpty());
    
        // ---- test remove ---
        st2.remove(box2);
        results = st2.query2D(box2);
        assertEquals(0, results.size());
        
        results = st2.query2D(srch);
        assertEquals(2, results.size());
        expected = new HashSet<Interval2D<Integer>>();
        expected.add(box1);
        expected.add(box3);
        for (Interval2D<Integer> ii : results) {
            assertTrue(expected.remove(ii));
        }
        assertTrue(expected.isEmpty());
        
    }
    
    public void testInsertRemoveHierarchicalGrid0() {
        /*
        inserting subdivisions of 9 neighbors  
        from small to larger sizes.
        */
        
        QuadTreeInterval2D<Integer, PairInt> qt
            = new QuadTreeInterval2D<Integer, PairInt>();
         
        int xc = 9;
        int yc = 9;
        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;
        
        List<Interval2D<Integer>> list = 
            new ArrayList<Interval2D<Integer>>();
                 
        Interval<Integer> aX = new Interval<Integer>(9, 9);
        Interval<Integer> aY = new Interval<Integer>(9, 9);
        Interval2D<Integer> a = new Interval2D<Integer>(aX, aY);
        qt.insert(a, new PairInt(xc, yc));
        list.add(a);
        
        for (int factor = 1; factor < 3; ++factor) {
            for (int k = 0; k < dxs.length; ++k) {
                int x2 = xc + factor*dxs[k];
                int y2 = xc + factor*dys[k];
                if (xc < x2) {
                    aX = new Interval<Integer>(xc, x2);
                } else {
                    aX = new Interval<Integer>(x2, xc);
                }
                if (yc < y2) {
                    aY = new Interval<Integer>(yc, y2);
                } else {
                    aY = new Interval<Integer>(y2, yc);
                }
                a = new Interval2D<Integer>(aX, aY);
                qt.insert(a, new PairInt(x2, y2));
                list.add(a);
            }
        }
      
        List<Interval2D<Integer>> results = null;
        for (int i = 0; i < list.size(); ++i) {
            a = list.get(i);
            results = qt.query2D(a);
            boolean found = false;
            for (Interval2D<Integer> r : results) {
                if (r.equals(a)) {
                    found = true;
                    break;
                }
            }
            assertTrue(found);
        }
        
        // remove one at a time in order of insert and assert
        // remaining queries
        for (int i = 0; i < list.size(); ++i) {
            a = list.get(i);
            //System.out.println("i=" + i);
            qt.remove(a);
            results = qt.query2D(a);
            boolean found = false;
            for (Interval2D<Integer> r : results) {
                if (r.equals(a)) {
                    found = true;
                    break;
                }
            }
            assertFalse(found);
            for (int j = (i + 1); j < list.size(); ++j) {
                a = list.get(j);
                results = qt.query2D(a);
                found = false;
                for (Interval2D<Integer> r : results) {
                    if (r.equals(a)) {
                        found = true;
                        break;
                    }
                }
                assertTrue(found);
            }
        }
    }
    
    public void testInsertRemoveHierarchicalGrid1() {
        /*
        inserting subdivisions of 9 neighbors  
        from small to larger sizes.
        */
        
        QuadTreeInterval2D<Integer, PairInt> qt
            = new QuadTreeInterval2D<Integer, PairInt>();
         
        int xc = 9;
        int yc = 9;
        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;
        
        List<Interval2D<Integer>> list = 
            new ArrayList<Interval2D<Integer>>();
                 
        Interval<Integer> aX = new Interval<Integer>(9, 9);
        Interval<Integer> aY = new Interval<Integer>(9, 9);
        Interval2D<Integer> a = new Interval2D<Integer>(aX, aY);
        qt.insert(a, new PairInt(xc, yc));
        list.add(a);
        
        for (int factor = 1; factor < 3; ++factor) {
            for (int k = 0; k < dxs.length; ++k) {
                int x2 = xc + factor*dxs[k];
                int y2 = xc + factor*dys[k];
                if (xc < x2) {
                    aX = new Interval<Integer>(xc, x2);
                } else {
                    aX = new Interval<Integer>(x2, xc);
                }
                if (yc < y2) {
                    aY = new Interval<Integer>(yc, y2);
                } else {
                    aY = new Interval<Integer>(y2, yc);
                }
                a = new Interval2D<Integer>(aX, aY);
                qt.insert(a, new PairInt(x2, y2));
                list.add(a);
            }
        }
      
        List<Interval2D<Integer>> results = null;
        for (int i = 0; i < list.size(); ++i) {
            a = list.get(i);
            results = qt.query2D(a);
            boolean found = false;
            for (Interval2D<Integer> r : results) {
                if (r.equals(a)) {
                    found = true;
                    break;
                }
            }
            assertTrue(found);
        }
        
        // remove one at a time in lifo order
        for (int i = (list.size() - 1); i > -1; --i) {
            a = list.get(i);
            //System.out.println("i=" + i);
            qt.remove(a);
            results = qt.query2D(a);
            boolean found = false;
            for (Interval2D<Integer> r : results) {
                if (r.equals(a)) {
                    found = true;
                    break;
                }
            }
            assertFalse(found);
            for (int j = (i - 1); j > -1; --j) {
                a = list.get(j);
                results = qt.query2D(a);
                found = false;
                for (Interval2D<Integer> r : results) {
                    if (r.equals(a)) {
                        found = true;
                        break;
                    }
                }
                assertTrue(found);
            }
        }
    }
    
    public void testInsertRemoveHierarchicalGrid2() {
        /*
        inserting subdivisions of 9 neighbors  
        from larger to smaller sizes.
        */
        
        QuadTreeInterval2D<Integer, PairInt> qt
            = new QuadTreeInterval2D<Integer, PairInt>();
         
        int xc = 9;
        int yc = 9;
        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;
        
        List<Interval2D<Integer>> list = 
            new ArrayList<Interval2D<Integer>>();
                 
        Interval<Integer> aX = new Interval<Integer>(9, 10);
        Interval<Integer> aY = new Interval<Integer>(9, 10);
        Interval2D<Integer> a = new Interval2D<Integer>(aX, aY);
        qt.insert(a, new PairInt(xc, yc));
        list.add(a);
        
        for (int factor = 2; factor > 0; --factor) {
            for (int k = 0; k < dxs.length; ++k) {
                int x2 = xc + factor*dxs[k];
                int y2 = xc + factor*dys[k];
                if (xc < x2) {
                    aX = new Interval<Integer>(xc, x2);
                } else {
                    aX = new Interval<Integer>(x2, xc);
                }
                if (yc < y2) {
                    aY = new Interval<Integer>(yc, y2);
                } else {
                    aY = new Interval<Integer>(y2, yc);
                }
                a = new Interval2D<Integer>(aX, aY);
                qt.insert(a, new PairInt(x2, y2));
                list.add(a);
            }
        }
      
        List<Interval2D<Integer>> results = null;
        for (int i = 0; i < list.size(); ++i) {
            a = list.get(i);
            results = qt.query2D(a);
            boolean found = false;
            for (Interval2D<Integer> r : results) {
                if (r.equals(a)) {
                    found = true;
                    break;
                }
            }
            assertTrue(found);
        }
        
        // remove one at a time in order of insert and assert
        // remaining queries
        for (int i = 0; i < list.size(); ++i) {
            a = list.get(i);
            //System.out.println("i=" + i);
            qt.remove(a);
            results = qt.query2D(a);
            boolean found = false;
            for (Interval2D<Integer> r : results) {
                if (r.equals(a)) {
                    found = true;
                    break;
                }
            }
            assertFalse(found);
            for (int j = (i + 1); j < list.size(); ++j) {
                a = list.get(j);
                results = qt.query2D(a);
                found = false;
                for (Interval2D<Integer> r : results) {
                    if (r.equals(a)) {
                        found = true;
                        break;
                    }
                }
                assertTrue(found);
            }
        }
    }
    
    public void testInsertRemoveHierarchicalGrid3() {
        /*
        inserting subdivisions of 9 neighbors  
        from larger to smaller sizes.
        */
        
        QuadTreeInterval2D<Integer, PairInt> qt
            = new QuadTreeInterval2D<Integer, PairInt>();
         
        int xc = 9;
        int yc = 9;
        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;
        
        List<Interval2D<Integer>> list = 
            new ArrayList<Interval2D<Integer>>();
                 
        Interval<Integer> aX = new Interval<Integer>(9, 10);
        Interval<Integer> aY = new Interval<Integer>(9, 10);
        Interval2D<Integer> a = new Interval2D<Integer>(aX, aY);
        qt.insert(a, new PairInt(xc, yc));
        list.add(a);
        
        for (int factor = 2; factor > 0; --factor) {
            for (int k = 0; k < dxs.length; ++k) {
                int x2 = xc + factor*dxs[k];
                int y2 = xc + factor*dys[k];
                if (xc < x2) {
                    aX = new Interval<Integer>(xc, x2);
                } else {
                    aX = new Interval<Integer>(x2, xc);
                }
                if (yc < y2) {
                    aY = new Interval<Integer>(yc, y2);
                } else {
                    aY = new Interval<Integer>(y2, yc);
                }
                a = new Interval2D<Integer>(aX, aY);
                qt.insert(a, new PairInt(x2, y2));
                list.add(a);
            }
        }
      
        List<Interval2D<Integer>> results = null;
        for (int i = 0; i < list.size(); ++i) {
            a = list.get(i);
            results = qt.query2D(a);
            boolean found = false;
            for (Interval2D<Integer> r : results) {
                if (r.equals(a)) {
                    found = true;
                    break;
                }
            }
            assertTrue(found);
        }
        
        for (int i = (list.size() - 1); i > -1; --i) {
            a = list.get(i);
            //System.out.println("i=" + i);
            qt.remove(a);
            results = qt.query2D(a);
            boolean found = false;
            for (Interval2D<Integer> r : results) {
                if (r.equals(a)) {
                    found = true;
                    break;
                }
            }
            assertFalse(found);
            for (int j = (i - 1); j > -1; --j) {
                a = list.get(j);
                results = qt.query2D(a);
                found = false;
                for (Interval2D<Integer> r : results) {
                    if (r.equals(a)) {
                        found = true;
                        break;
                    }
                }
                assertTrue(found);
            }
        }
    }
    
    public void testRandomInsertRemove() throws Exception {
     
        SecureRandom r = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.nanoTime();
        //seed = 722829906140326L;
        System.out.println("SEED=" + seed);
        r.setSeed(seed);

        int xmin = 10;
        int xmax = 1000;
        int ymin = 10;
        int ymax = 1000;

        int n = 100;//1000;

        QuadTreeInterval2D<Integer, PairInt> qt
            = new QuadTreeInterval2D<Integer, PairInt>();    
        List<Interval2D<Integer>> list = 
            new ArrayList<Interval2D<Integer>>();
        Interval<Integer> aX = null;
        Interval<Integer> aY = null;
        Interval2D<Integer> a = null;
        
        Set<QuadInt> added = new HashSet<QuadInt>();
        
        int xrange = xmax - xmin;
        int yrange = ymax - ymin;
        
        // repeat a cycle of insert, query, delete some, query
        int nCycle = 10;
        for (int ii = 0; ii < nCycle; ++ii) {
        
            for (int i = 0; i < n; ++i) {
                int x1 = r.nextInt(xrange - 1);
                int xr = xrange - x1;
                x1 += xmin;
                int x2 = r.nextInt(xr) + x1;

                int y1 = r.nextInt(yrange - 1);
                int yr = yrange - y1;
                y1 += ymin;
                int y2 = r.nextInt(yr) + y1;

                QuadInt chk = new QuadInt(x1, y1, x2, y2);
                if (added.contains(chk)) {
                    continue;
                }
                added.add(chk);
                aX = new Interval<Integer>(x1, x2);
                aY = new Interval<Integer>(y1, y2);
                a = new Interval2D<Integer>(aX, aY);
                qt.insert(a, new PairInt(x1, y1));
                list.add(a);
            }
       
            // query for all
            List<Interval2D<Integer>> results = null;
            for (int i = 0; i < list.size(); ++i) {
                a = list.get(i);
                results = qt.query2D(a);
                boolean found = false;
                for (Interval2D<Integer> result : results) {
                    if (result.equals(a)) {
                        found = true;
                        break;
                    }
                }
                assertTrue(found);
            }
        
            // randomly remove several
            TIntList rmList = new TIntArrayList();
            TIntSet rmSet = new TIntHashSet();
            int nRm = list.size()/3;
            for (int i = 0; i < nRm; ++i) {
                int d = r.nextInt(n);
                if (!rmSet.contains(d)){
                    rmList.add(d);
                    rmSet.add(d);
                }
            }
            rmList.sort();
            //System.out.println("rm.sz=" + rmList.size() 
            //    + " list.sz=" + list.size());
            List<Interval2D<Integer>> treeNodes =
                new ArrayList<Interval2D<Integer>>();
            
            List<PairInt> values = new ArrayList<PairInt>();
            qt.getAllNodes(treeNodes, values);
            assertEquals(list.size(), treeNodes.size());
            //debug(treeNodes);
            
            for (int i = (rmList.size() - 1); i > -1; --i) {
                int rmIdx = rmList.get(i);
                a = list.get(rmIdx);
                assertNotNull(a);
                
                qt.remove(a);
                assertNotNull(list.remove(rmIdx)); 
                
                //treeNodes.clear(); 
                //values.clear();
                //qt.getAllNodes(treeNodes, values);
                //assertEquals(list.size(), treeNodes.size());
                //debug(treeNodes);
                
                QuadInt chk = new QuadInt(
                    a.intervalX.min(), a.intervalX.max(),
                    a.intervalY.min(), a.intervalY.max());
                added.remove(chk);
                results = qt.query2D(a);
                boolean found = false;
                for (Interval2D<Integer> result : results) {
                    if (result.equals(a)) {
                        found = true;
                        break;
                    }
                }
                assertFalse(found);
                for (int j = 0; j < list.size(); ++j) {
                    a = list.get(j);
                    //System.out.println("cycle=" + ii + " i=" + i +
                    //    " j=" + j + " list.sz=" + list.size()
                    //+ " q=" + a.toString());
                    results = qt.query2D(a);
                    found = false;
                    for (Interval2D<Integer> result : results) {
                        if (result.equals(a)) {
                            found = true;
                            break;
                        }
                    }
                    assertTrue(found);
                }
            }
        
            // query all
            for (int i = 0; i < list.size(); ++i) {
                a = list.get(i);
                //System.out.println("cycle=" + ii + " q i=" + i);
                results = qt.query2D(a);
                boolean found = false;
                for (Interval2D<Integer> result : results) {
                    if (result.equals(a)) {
                        found = true;
                        break;
                    }
                }
                assertTrue(found);
            }
        }
    }
    
    public void testRandomInsertRemove_2() throws Exception {
     
        SecureRandom r = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.nanoTime();
        //seed = 722829906140326L;
        System.out.println("SEED=" + seed);
        r.setSeed(seed);

        int xmin = 10;
        int xmax = 1000;
        int ymin = 10;
        int ymax = 1000;

        int n = 100;//1000;

        QuadTreeInterval2D<Integer, PairInt> qt
            = new QuadTreeInterval2D<Integer, PairInt>();    
        List<Interval2D<Integer>> list = 
            new ArrayList<Interval2D<Integer>>();
        Interval<Integer> aX = null;
        Interval<Integer> aY = null;
        Interval2D<Integer> a = null;
        
        Set<QuadInt> added = new HashSet<QuadInt>();
        
        int xrange = xmax - xmin;
        int yrange = ymax - ymin;
        
        // repeat a cycle of insert, query, delete some, query
        int nCycle = 10;
        for (int ii = 0; ii < nCycle; ++ii) {
        
            for (int i = 0; i < n; ++i) {
                int x1 = r.nextInt(xrange - 1);
                int xr = xrange - x1;
                x1 += xmin;
                int x2 = r.nextInt(xr) + x1;

                int y1 = r.nextInt(yrange - 1);
                int yr = yrange - y1;
                y1 += ymin;
                int y2 = r.nextInt(yr) + y1;

                QuadInt chk = new QuadInt(x1, y1, x2, y2);
                if (added.contains(chk)) {
                    continue;
                }
                added.add(chk);
                aX = new Interval<Integer>(x1, x2);
                aY = new Interval<Integer>(y1, y2);
                a = new Interval2D<Integer>(aX, aY);
                qt.insert(a, new PairInt(x1, y1));
                list.add(a);
            }
       
            // query for all
            Map<Interval2D<Integer>, PairInt> results = null;
                        
            for (int i = 0; i < list.size(); ++i) {
                a = list.get(i);
                results = qt.query2D_2(a);
                boolean found = results.containsKey(a);
                assertTrue(found);
            }
        
            // randomly remove several
            TIntList rmList = new TIntArrayList();
            TIntSet rmSet = new TIntHashSet();
            int nRm = list.size()/3;
            for (int i = 0; i < nRm; ++i) {
                int d = r.nextInt(n);
                if (!rmSet.contains(d)){
                    rmList.add(d);
                    rmSet.add(d);
                }
            }
            rmList.sort();
            //System.out.println("rm.sz=" + rmList.size() 
            //    + " list.sz=" + list.size());
            List<Interval2D<Integer>> treeNodes =
                new ArrayList<Interval2D<Integer>>();
            
            List<PairInt> values = new ArrayList<PairInt>();
            qt.getAllNodes(treeNodes, values);
            assertEquals(list.size(), treeNodes.size());
            //debug(treeNodes);
            
            for (int i = (rmList.size() - 1); i > -1; --i) {
                int rmIdx = rmList.get(i);
                a = list.get(rmIdx);
                assertNotNull(a);
                
                qt.remove(a);
                assertNotNull(list.remove(rmIdx)); 
                
                //treeNodes.clear(); 
                //values.clear();
                //qt.getAllNodes(treeNodes, values);
                //assertEquals(list.size(), treeNodes.size());
                //debug(treeNodes);
                
                QuadInt chk = new QuadInt(
                    a.intervalX.min(), a.intervalX.max(),
                    a.intervalY.min(), a.intervalY.max());
                added.remove(chk);
                results = qt.query2D_2(a);
                boolean found = results.containsKey(a);
                assertFalse(found);
                for (int j = 0; j < list.size(); ++j) {
                    a = list.get(j);
                    //System.out.println("cycle=" + ii + " i=" + i +
                    //    " j=" + j + " list.sz=" + list.size()
                    //+ " q=" + a.toString());
                    results = qt.query2D_2(a);
                    found = results.containsKey(a);
                    assertTrue(found);
                }
            }
        
            // query all
            for (int i = 0; i < list.size(); ++i) {
                a = list.get(i);
                //System.out.println("cycle=" + ii + " q i=" + i);
                results = qt.query2D_2(a);
                boolean found = false;
                results = qt.query2D_2(a);
                found = results.containsKey(a);
                assertTrue(found);
            }
        }
    }
    
    private void debug(List<Interval2D<Integer>> treeNodes) {
        StringBuilder sb = new StringBuilder();
        for (Interval2D<Integer> g : treeNodes) {
            sb.append(g.toString()).append("\n");
        }
        System.out.println(sb.toString());
    }
    
}
