package algorithms.util;

import java.security.SecureRandom;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class PairIntArrayTest extends TestCase {
    
    public PairIntArrayTest(String testName) {
        super(testName);
    }

    @Override
    protected void setUp() throws Exception {
        super.setUp();
    }

    @Override
    protected void tearDown() throws Exception {
        super.tearDown();
    }

    public void testAddGet() {
        
        System.out.println("testAddGet");
        
        int xPoint = 1;
        int yPoint = 10;
        PairIntArray instance = new PairIntArray();
        instance.add(xPoint, yPoint);
        assertTrue(instance.getX(0) == xPoint);
        assertTrue(instance.getY(0) == yPoint);
        
        instance = new PairIntArray();
        SecureRandom sr = new SecureRandom();
        sr.setSeed(System.currentTimeMillis());
        int nr = 51;
        int[] x = new int[nr];
        int[] y = new int[nr];
        for (int i = 0; i < nr; i++) {
            x[i] = sr.nextInt();
            y[i] = sr.nextInt();
            instance.add(x[i], y[i]);
        }
        
        assertTrue(instance.getN() == nr);
        
        for (int i = 0; i < nr; i++) {
            assertTrue(x[i] == instance.getX(i));
            assertTrue(y[i] == instance.getY(i));
        }                
    }
    
    public void testInsertSpaceAtTopOfArrays() throws Exception {
        
        PairIntArray xy = new PairIntArray();
        for (int i = 0; i < 3; i++) {
            xy.add(i, i);
        }
        assertTrue(xy.getN() == 3);
        for (int i = 0; i < 3; i++) {
            assertTrue(xy.getX(i) == i);
            assertTrue(xy.getY(i) == i);
        }
        
        int nInsert = 10;
        
        xy.insertSpaceAtTopOfArrays(nInsert);
        
        for (int i = 0; i < 3; i++) {
            int ii = i + nInsert;
            assertTrue(xy.getX(ii) == i);
            assertTrue(xy.getY(ii) == i);
        }
        for (int i = 0; i < nInsert; i++) {
            assertTrue(xy.getX(i) == 0);
            assertTrue(xy.getY(i) == 0);
        }
                
        // test the insert space when the arrays are already large enough
        xy = new PairIntArray(20);
        for (int i = 0; i < 3; i++) {
            xy.add(i, i);
        }
        assertTrue(xy.getN() == 3);
        for (int i = 0; i < 3; i++) {
            assertTrue(xy.getX(i) == i);
            assertTrue(xy.getY(i) == i);
        }
        
        xy.insertSpaceAtTopOfArrays(nInsert);
        
        for (int i = 0; i < 3; i++) {
            assertTrue(xy.getX(i + nInsert) == i);
            assertTrue(xy.getY(i + nInsert) == i);
        }
        for (int i = 0; i < nInsert; i++) {
            assertTrue(xy.getX(i) == 0);
            assertTrue(xy.getY(i) == 0);
        }
        for (int i = 0; i < nInsert; i++) {
            xy.set(i, 1000000, 1000000);
        }
        for (int i = 0; i < nInsert; i++) {
            assertTrue(xy.getX(i) == 1000000);
            assertTrue(xy.getY(i) == 1000000);
        }
        for (int i = 0; i < 3; i++) {
            assertTrue(xy.getX(i + nInsert) == i);
            assertTrue(xy.getY(i + nInsert) == i);
        }
    }
    
    public void testReverse() throws Exception {
        
        PairIntArray xy = new PairIntArray();
        for (int i = 0; i < 5; i++) {
            xy.add(i, i);
        }
        int n = xy.getN();
        assertTrue(n == 5);
        xy.reverse();
        for (int i = 0; i < 5; i++) {
            assertTrue(xy.getX(n - i - 1) == i);
            assertTrue(xy.getY(n - i - 1) == i);
        }
        
        xy = new PairIntArray();
        for (int i = 0; i < 4; i++) {
            xy.add(i, i);
        }
        n = xy.getN();
        assertTrue(n == 4);
        xy.reverse();
        for (int i = 0; i < n; i++) {
            assertTrue(xy.getX(n - i - 1) == i);
            assertTrue(xy.getY(n - i - 1) == i);
        }
    }
    
    public void testSwap() throws Exception {
        
        PairIntArray xy0 = new PairIntArray();
        PairIntArray xy1 = new PairIntArray();
        
        int n0 = 4;
        int n1 = 14;
        
        for (int i = 0; i < n0; i++) {
            xy0.add(n0, n0);
        }
        
        for (int i = 0; i < n1; i++) {
            xy1.add(n1, n1);
        }
        
        assertTrue(xy0.getN() == n0);
        assertTrue(xy1.getN() == n1);
        
        xy0.swapContents(xy1);
        
        assertTrue(xy0.getN() == n1);
        assertTrue(xy1.getN() == n0);
        
        for (int i = 0; i < n0; i++) {
            assertTrue(xy1.getX(i) == n0);
            assertTrue(xy1.getY(i) == n0);
        }
        
        for (int i = 0; i < n1; i++) {
            assertTrue(xy0.getX(i) == n1);
            assertTrue(xy0.getY(i) == n1);
        }
    }
    
    public void testClone() throws Exception {
        
        PairIntArray xy = new PairIntArray();
        for (int i = 0; i < 5; i++) {
            xy.add(i, i);
        }
        int n = xy.getN();
        assertTrue(n == 5);
        
        PairIntArray copy = xy.copy();
        assertTrue(copy.getN() == xy.getN());
        for (int i = 0; i < 5; i++) {
           assertTrue(xy.getX(i) == copy.getX(i));
           assertTrue(xy.getY(i) == copy.getY(i));
        }
        
    }
    
    public void testRemoveRange() throws Exception {
        
        PairIntArray xy = new PairIntArray();
        for (int i = 0; i < 7; i++) {
            xy.add(i, i);
        }
        
        xy.removeRange(0, 2);
        
        assertTrue(xy.getN() == 4);
        
        assertTrue(xy.getX(0) == 3);
        assertTrue(xy.getY(0) == 3);
        assertTrue(xy.getX(1) == 4);
        assertTrue(xy.getY(1) == 4);
        assertTrue(xy.getX(2) == 5);
        assertTrue(xy.getY(2) == 5);
        assertTrue(xy.getX(3) == 6);
        assertTrue(xy.getY(3) == 6);
        
        //=======
        xy = new PairIntArray();
        for (int i = 0; i < 7; i++) {
            xy.add(i, i);
        }
        
        xy.removeRange(4, 6);
        
        assertTrue(xy.getN() == 4);
        
        assertTrue(xy.getX(0) == 0);
        assertTrue(xy.getY(0) == 0);
        assertTrue(xy.getX(1) == 1);
        assertTrue(xy.getY(1) == 1);
        assertTrue(xy.getX(2) == 2);
        assertTrue(xy.getY(2) == 2);
        assertTrue(xy.getX(3) == 3);
        assertTrue(xy.getY(3) == 3);
        
        xy = new PairIntArray();
        xy.add(0, 10);
        xy.add(12, 23);
        xy.add(23, 31);
        xy.add(34, 84);
        int i = 2;
        int r1 = 31;
        xy.set(i - 1, xy.getX(i - 1), r1);
        xy.removeRange(i, i);
        assertTrue(xy.getN() == 3);
        assertTrue((xy.getX(0) == 0) && (xy.getY(0) == 10));
        assertTrue((xy.getX(1) == 12) && (xy.getY(1) == 31));
        assertTrue((xy.getX(2) == 34) && (xy.getY(2) == 84));

    }
    
    public void testToPairFloatArray() throws Exception {
        
        PairIntArray xy = new PairIntArray();
        for (int i = 0; i < 7; i++) {
            xy.add(i, i);
        }
        
        PairFloatArray xyF = xy.toPairFloatArray();
        for (int i = 0; i < 7; i++) {
            assertTrue(Math.abs(xy.getX(i) - xyF.getX(i)) < 0.001f);
            assertTrue(Math.abs(xy.getY(i) - xyF.getY(i)) < 0.001f);
        }
    }
    
    public void testAddAll() throws Exception {
        
        int[] x0 = new int[]{0, 1, 2};
        int[] y0 = new int[]{10, 11, 12};
        
        int[] x1 = new int[]{3, 4, 5};
        int[] y1 = new int[]{13, 14, 15};
        
        PairIntArray p0 = new PairIntArray();
        for (int i = 0; i < x0.length; i++) {
            int x = x0[i];
            int y = y0[i];
            p0.add(x, y);
        }
        assertTrue(p0.getN() == x0.length);
        for (int i = 0; i < p0.getN(); i++) {
            int x = x0[i];
            int y = y0[i];
            assertTrue(p0.getX(i) == x);
            assertTrue(p0.getY(i) == y);
        }
        
        PairIntArray p1 = new PairIntArray();
        for (int i = 0; i < x1.length; i++) {
            int x = x1[i];
            int y = y1[i];
            p1.add(x, y);
        }
        assertTrue(p1.getN() == x1.length);
        for (int i = 0; i < p1.getN(); i++) {
            int x = x1[i];
            int y = y1[i];
            assertTrue(p1.getX(i) == x);
            assertTrue(p1.getY(i) == y);
        }
        
        
        PairIntArray pAppend = new PairIntArray();
        pAppend.addAll(p0);
        
        assertTrue(p0.getN() == pAppend.getN());
        for (int i = 0; i < x0.length; i++) {
            int x = x0[i];
            int y = y0[i];
            assertTrue(pAppend.getX(i) == x);
            assertTrue(pAppend.getY(i) == y);
        }
        
        pAppend.addAll(p1);
        for (int i = 0; i < x0.length; i++) {
            int x = x0[i];
            int y = y0[i];
            assertTrue(pAppend.getX(i) == x);
            assertTrue(pAppend.getY(i) == y);
        }
        for (int i = 0; i < x0.length; i++) {
            int x = x1[i];
            int y = y1[i];
            int idx = x0.length + i;
            assertTrue(pAppend.getX(idx) == x);
            assertTrue(pAppend.getY(idx) == y);
        }
        
    }
    
    public void testInsertAll() throws Exception {
        
        PairIntArray p = new PairIntArray(5);
        for (int i = 0; i < 5; ++i) {
            p.add(i, 10 + i);
        }
        
        PairIntArray p2 = new PairIntArray(5);
        for (int i = 0; i < 3; ++i) {
            p2.add(100 + i, 200 + i);
        }
        
        p.insertAll(1, p2);
        
        assertTrue(p.getN() == 8);
        
        assertTrue(p.getX(0) == 0);
        assertTrue(p.getY(0) == 10);
        
        assertTrue(p.getX(1) == (100 + 0));
        assertTrue(p.getY(1) == (200 + 0));
        assertTrue(p.getX(2) == (100 + 1));
        assertTrue(p.getY(2) == (200 + 1));
        assertTrue(p.getX(3) == (100 + 2));
        assertTrue(p.getY(3) == (200 + 2));
        
        assertTrue(p.getX(4) == (0 + 1));
        assertTrue(p.getY(4) == (10 + 1));
        assertTrue(p.getX(5) == (0 + 2));
        assertTrue(p.getY(5) == (10 + 2));
        assertTrue(p.getX(6) == (0 + 3));
        assertTrue(p.getY(6) == (10 + 3));
        assertTrue(p.getX(7) == (0 + 4));
        assertTrue(p.getY(7) == (10 + 4));
    }
    
    public void testInsertAll2() throws Exception {
        
        PairIntArray p = new PairIntArray(3);
        for (int i = 0; i < 3; ++i) {
            p.add(i, 10 + i);
        }
        
        PairIntArray p2 = new PairIntArray(1);
        for (int i = 0; i < 1; ++i) {
            p2.add(100 + i, 200 + i);
        }
        
        p.insertAll(3, p2);
        
        assertTrue(p.getN() == 4);
        
        assertTrue(p.getX(0) == 0);
        assertTrue(p.getY(0) == 10);
        assertTrue(p.getX(1) == (0 + 1));
        assertTrue(p.getY(1) == (10 + 1));
        assertTrue(p.getX(2) == (0 + 2));
        assertTrue(p.getY(2) == (10 + 2));
        
        assertTrue(p.getX(3) == (100 + 0));
        assertTrue(p.getY(3) == (200 + 0));
        
        boolean threwException = false;
        try {
            p.insertAll(-1, p2);
        } catch(Throwable t) {
            threwException = true;
        }
        assertTrue(threwException);
        
        threwException = false;
        try {
            p.insertAll(1000, p2);
        } catch(Throwable t) {
            threwException = true;
        }
        assertTrue(threwException);
    }

}
