package algorithms.packing;

import algorithms.util.PairInt;
import algorithms.util.PixelHelper;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.HashSet;
import java.util.Set;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class Intersection2DPackingTest extends TestCase {
    
    public Intersection2DPackingTest() {
    }
    
    public void testPopulate() {
         
        PixelHelper ph = new PixelHelper();
        int[] xy = new int[2];
        
        int w = 4;
        int h = 4;
        
        TIntSet intersection = new TIntHashSet();
        
        for (int x = 0; x < w; ++x) {
            for (int y = 0; y < h; ++y) {
                int pixIdx = (int)ph.toPixelIndex(x, y, w);
                if ((x & 1) == 1) {
                    if ((y & 1) == 1) {
                        intersection.add(pixIdx);
                    }
                }
            }
        }
        Intersection2DPacking ip = new Intersection2DPacking();
        
        int[] xs = new int[intersection.size()];
        int[] ys = new int[intersection.size()];
        
        ip._populate(intersection, w, xs, ys);
        
        for (int i = 0; i < xs.length; ++i) {
            int pixIdx = (int)ph.toPixelIndex(xs[i], ys[i], w);
            assertTrue(intersection.remove(pixIdx));
        }
        
        assertTrue(intersection.isEmpty());
    }
    
    public void testPopulate2() {
        
        int w = 4;
        int h = 4;
        
        Set<PairInt> intersection = new HashSet<PairInt>();
        
        for (int x = 0; x < w; ++x) {
            for (int y = 0; y < h; ++y) {
                PairInt p = new PairInt(x, y);
                if ((x & 1) == 1) {
                    if ((y & 1) == 1) {
                        intersection.add(p);
                    }
                }
            }
        }
        Intersection2DPacking ip = new Intersection2DPacking();
        
        int[] xs = new int[intersection.size()];
        int[] ys = new int[intersection.size()];
        
        ip._populate(intersection, xs, ys);
        
        for (int i = 0; i < xs.length; ++i) {
            PairInt p = new PairInt(xs[i], ys[i]);
            assertTrue(intersection.remove(p));
        }
        
        assertTrue(intersection.isEmpty());
    }
    
    public void testIntersection() {
        
        PixelHelper ph = new PixelHelper();
        
        int w = 4;
        int h = 4;
        
        TIntSet expected = new TIntHashSet();
        TIntSet p1 = new TIntHashSet();
        TIntSet p2 = new TIntHashSet();
        
        for (int x = 0; x < w; ++x) {
            for (int y = 0; y < h; ++y) {
                int pixIdx = (int)ph.toPixelIndex(x, y, w);
                if ((x & 1) == 1) {
                    if ((y & 1) == 1) {
                        expected.add(pixIdx);
                        p2.add(pixIdx);
                        p1.add(pixIdx);
                    } else {
                        p1.add(pixIdx);
                    }
                } else {
                    p2.add(pixIdx);
                }
                
            }
        }
        Intersection2DPacking ip = new Intersection2DPacking();
        TIntSet result = ip.intersection(p1, p2);
        
        assertEquals(expected.size(), result.size());
        
        assertTrue(expected.containsAll(result));
        
        result = ip.naiveStripPacking(expected, w, 3);
        
        assertEquals(1, result.size());
    }
    
    public void testIntersection2() {
        
        int w = 4;
        int h = 4;
        
        Set<PairInt> expected = new HashSet<PairInt>();
        Set<PairInt> p1 = new HashSet<PairInt>();
        Set<PairInt> p2 = new HashSet<PairInt>();
        
        for (int x = 0; x < w; ++x) {
            for (int y = 0; y < h; ++y) {
                PairInt p = new PairInt(x, y);
                if ((x & 1) == 1) {
                    if ((y & 1) == 1) {
                        expected.add(p);
                        p2.add(p);
                        p1.add(p);
                    } else {
                        p1.add(p);
                    }
                } else {
                    p2.add(p);
                }
                
            }
        }
        Intersection2DPacking ip = new Intersection2DPacking();
        Set<PairInt> result = ip.intersection(p1, p2);
        
        assertEquals(expected.size(), result.size());
        
        assertTrue(expected.containsAll(result));
        
        result = ip.naiveStripPacking(expected, 3);
        
        assertEquals(1, result.size());
    }
    
    public void testNaiveStripPacking() {
        
        PixelHelper ph = new PixelHelper();
        
        int w = 6;
        int h = 6;
        int cellSize = 3;
        
        /*
        5
        4
        3 *        *
        2 
        1 
        0 *        *       
          0  1  2  3  4  5  6
        */
        
        TIntSet expected = new TIntHashSet();
        TIntSet p1 = new TIntHashSet();
        TIntSet p2 = new TIntHashSet();
        
        expected.add((int)ph.toPixelIndex(0, 0, w));
        expected.add((int)ph.toPixelIndex(0, 3, w));
        expected.add((int)ph.toPixelIndex(3, 0, w));
        expected.add((int)ph.toPixelIndex(3, 3, w));
        
        for (int x = 0; x < w; ++x) {
            for (int y = 0; y < h; ++y) {
                int pixIdx = (int)ph.toPixelIndex(x, y, w);
                p2.add(pixIdx);
                p1.add(pixIdx);
            }
        }
        Intersection2DPacking ip = new Intersection2DPacking();
        TIntSet result = ip.naiveStripPacking(p1, p2, w, cellSize);
        
        assertEquals(expected.size(), result.size());
        
        assertTrue(expected.containsAll(result));
    }
    
    public void testNaiveStripPacking2() {
        
        int w = 6;
        int h = 6;
        int cellSize = 3;
        
        /*
        5
        4
        3 *        *
        2 
        1 
        0 *        *       
          0  1  2  3  4  5  6
        */
        
        Set<PairInt> expected = new HashSet<PairInt>();
        Set<PairInt> p1 = new HashSet<PairInt>();
        Set<PairInt> p2 = new HashSet<PairInt>();
        
        expected.add(new PairInt(0, 0));
        expected.add(new PairInt(0, 3));
        expected.add(new PairInt(3, 0));
        expected.add(new PairInt(3, 3));
        
        for (int x = 0; x < w; ++x) {
            for (int y = 0; y < h; ++y) {
                PairInt p = new PairInt(x, y);
                p2.add(p);
                p1.add(p);
            }
        }
        Intersection2DPacking ip = new Intersection2DPacking();
        Set<PairInt> result = ip.naiveStripPacking(p1, p2, cellSize);
        
        assertEquals(expected.size(), result.size());
        
        assertTrue(expected.containsAll(result));
    }
    
}
