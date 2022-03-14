package algorithms.compGeometry.convexHull;

import algorithms.misc.Misc;
import gnu.trove.list.TLongList;
import gnu.trove.list.array.TLongArrayList;
import gnu.trove.map.TLongObjectMap;
import gnu.trove.map.hash.TLongObjectHashMap;
import gnu.trove.set.TLongSet;
import gnu.trove.set.hash.TLongHashSet;
import java.util.Arrays;
import java.util.Random;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class LargestDistancePairTest extends TestCase {
    
    public LargestDistancePairTest() {
    }

    public void testFindLargestDistancePair() throws Exception {
        
        /*            2,6
         *
         *
         *     0,2   2,2 3,2
         *                      7,1
         *            2,0
         *         
         */
        int n = 6;
        long[] x = new long[]{0, 2, 7, 2, 2, 3};
        long[] y = new long[]{2, 2, 1, 6, 0, 2};
        
        long[] expResult = new long[]{7, 1, 2, 6};
        
        long[] result = LargestDistancePair.findLargestDistancePair(x, y);
        
        System.out.printf("result=%s\n", Arrays.toString(result));
    }
    
    public void testRandomInput() throws Exception {
        
        long seed = System.nanoTime();
        //seed = 375796260475819L;
        System.out.printf("seed=%d\n", seed);
        Random rand = Misc.getSecureRandom();
        rand.setSeed(seed);
                
        int nTests = 100;
        
        for (int i = 0; i < nTests; ++i) {
            runRandomInput(rand);
        }
    }
        
    private void runRandomInput(Random rand) throws Exception {
        // generate points in this coordinate range:
        int min = 1;
        int max = 10000;
        
        // generate this many points:
        int n = 100 + rand.nextInt(1000 - 100);
        
        //System.out.printf("n=%d  min=%d max=%d   max-min=%d\n", n, min, max, max-min);
        
        long[] x = new long[n];
        long[] y = new long[n];
        TLongObjectMap<TLongSet> points = new TLongObjectHashMap<>();
        
        int i = 0;
        TLongSet set;
        while (i < n) {
            x[i] = min + rand.nextInt(max - min);
            y[i] = min + rand.nextInt(max - min);
            set = points.get(x[i]);
            if ((set != null) && set.contains(y[i])) {
                continue;
            }
            if (set == null) {
                set = new TLongHashSet();
                points.put(x[i], set);
            }
            set.add(y[i]);
            ++i;
        }
                
        // use brute force to calculate the expected answer(s):
        TLongList expectedX0s = new TLongArrayList();
        TLongList expectedY0s = new TLongArrayList();
        TLongList expectedX1s = new TLongArrayList();
        TLongList expectedY1s = new TLongArrayList();
        
        bruteForceMaxDist(expectedX0s, expectedY0s, expectedX1s, expectedY1s, x, y);
        
        
        TLongList expectedXH0s = new TLongArrayList();
        TLongList expectedYH0s = new TLongArrayList();
        TLongList expectedXH1s = new TLongArrayList();
        TLongList expectedYH1s = new TLongArrayList();
        GrahamScanLong.CH ch = GrahamScanLong.computeHull(x, y);
        bruteForceMaxDist(expectedXH0s, expectedYH0s, expectedXH1s, expectedYH1s, 
                ch.getXH(), ch.getYH());
        
        assertTrue(expectedX0s.size() == expectedX1s.size() 
            && expectedX0s.size() == expectedY0s.size() && expectedX0s.size() == expectedY1s.size());
        
        assertFalse(expectedX0s.isEmpty());
        
        //System.out.printf("number of max dist pairs=%d\n", expectedX0s.size());
        
        // returns a pair in format [x0, y0, x1, y1]
        long[] result = LargestDistancePair.findLargestDistancePair(x, y);
        
        long xd;
        long yd;
        long distSq;
        
        long resultDistSq = ( (result[0] - result[2])*(result[0] - result[2]) +
            (result[1] - result[3])*(result[1] - result[3]));
        //System.out.printf("result=(%d,%d)(%d,%d), d=%d\n", 
        //    result[0], result[1], result[2], result[3], resultDistSq);
        
        boolean found = false;
        for (i = 0; i < expectedXH0s.size(); ++i) {
            if (expectedXH0s.contains(result[0])) {
                if (expectedYH0s.contains(result[1])) {
                    if (expectedXH1s.contains(result[2])) {
                        if (expectedYH1s.contains(result[3])) {
                            found = true;
                        }
                    }
                }
            } else if (expectedXH0s.contains(result[2])) {
                if (expectedYH0s.contains(result[3])) {
                    if (expectedXH1s.contains(result[0])) {
                        if (expectedYH1s.contains(result[1])) {
                            found = true;
                        }
                    }
                }
            }
        }
        
        
        for (i = 0; i < expectedX0s.size(); ++i) {
            xd = expectedX0s.get(i) - expectedX1s.get(i);
            yd = expectedY0s.get(i) - expectedY1s.get(i);
            distSq = xd*xd + yd*yd;
            /*System.out.printf("expected for x, y =(%d,%d)(%d,%d) distSq=%d\n", 
                expectedX0s.get(i), expectedY0s.get(i), expectedX1s.get(i), 
                expectedY1s.get(i), distSq);*/
        }
        
        distSq = 0;
        for (i = 0; i < 1/*expectedXH0s.size()*/; ++i) {
            xd = expectedXH0s.get(i) - expectedXH1s.get(i);
            yd = expectedYH0s.get(i) - expectedYH1s.get(i);
            distSq = xd*xd + yd*yd;
            /*System.out.printf("expected for xh, yh =(%d,%d)(%d,%d) distSq=%d\n", 
                expectedXH0s.get(i), expectedYH0s.get(i), expectedXH1s.get(i), 
                expectedYH1s.get(i), distSq);*/
        }
        
        assertEquals(distSq, resultDistSq);
    }

    private void bruteForceMaxDist(TLongList expectedX0s, 
        TLongList expectedY0s, TLongList expectedX1s, TLongList expectedY1s, long[] x, long[] y) {
        
        int i;
        int n = x.length;
        
        long maxDist = Long.MIN_VALUE;
        int j;
        long x0;
        long y0;
        long xd;
        long yd;
        long distSq;
        for (i = 0; i < n; ++i) {
            x0 = x[i];
            y0 = y[i];
            for (j = i + 1; j < n; ++j) {
                xd = x[j] - x0;
                yd = y[j] - y0;
                distSq = xd * xd + yd * yd;
                if (distSq > maxDist) {
                    expectedX0s.clear();
                    expectedY0s.clear();
                    expectedX1s.clear();
                    expectedY1s.clear();
                    
                    maxDist = distSq;
                    
                    expectedX0s.add(x0);
                    expectedY0s.add(y0);
                    expectedX1s.add(x[j]);
                    expectedY1s.add(y[j]);
                } else if (distSq == maxDist) {
                    
                    expectedX0s.add(x0);
                    expectedY0s.add(y0);
                    expectedX1s.add(x[j]);
                    expectedY1s.add(y[j]);;
                }
            }
        }
        
    }
}
