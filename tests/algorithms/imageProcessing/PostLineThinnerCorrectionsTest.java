package algorithms.imageProcessing;

import algorithms.util.CornerArray;
import algorithms.util.PairInt;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import junit.framework.TestCase;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class PostLineThinnerCorrectionsTest extends TestCase {
    
    public PostLineThinnerCorrectionsTest() {
    }
    
    public void testCorrectForArtifacts() {
        
    }
    
    public void testReverseX() throws Exception {
        
        Set<PairInt> a = new HashSet<PairInt>();
        a.add(new PairInt(1, 1));
        a.add(new PairInt(2, 1));
        
        Set<PairInt> b = new HashSet<PairInt>();
        b.add(new PairInt(1, 1));
        b.add(new PairInt(2, 1));
        
        Set<PairInt> expected = new HashSet<PairInt>();
        expected.add(new PairInt(-1, 1));
        expected.add(new PairInt(-2, 1));
        
        PostLineThinnerCorrections pltc = new PostLineThinnerCorrections();
        pltc.reverseXs(a, b);
        for (PairInt p : a) {
            expected.remove(p);
        }
        assertTrue(expected.isEmpty());
        
        expected = new HashSet<PairInt>();
        expected.add(new PairInt(-1, 1));
        expected.add(new PairInt(-2, 1));
        
        for (PairInt p : b) {
            expected.remove(p);
        }
        assertTrue(expected.isEmpty());
        
        // ---------------------------
        a = new HashSet<PairInt>();
        a.add(new PairInt(1, 1));
        a.add(new PairInt(2, 1));
        
        b = new HashSet<PairInt>();
        b.add(new PairInt(1, 1));
        b.add(new PairInt(2, 1));
        
        Set<PairInt> c = new HashSet<PairInt>();
        c.add(new PairInt(1, 1));
        c.add(new PairInt(2, 1));
        
        Set<PairInt> d = new HashSet<PairInt>();
        d.add(new PairInt(1, 1));
        d.add(new PairInt(2, 1));
        
        pltc.reverseXs(a, b, c, d);
        
        expected = new HashSet<PairInt>();        
        expected.add(new PairInt(-1, 1));
        expected.add(new PairInt(-2, 1));
        
        for (PairInt p : a) {
            System.out.println("looking for x=" + p.getX() + " y=" + p.getY());
            assertNotNull(expected.remove(p));
        }
        assertTrue(expected.isEmpty());
        
        expected = new HashSet<PairInt>();        
        expected.add(new PairInt(-1, 1));
        expected.add(new PairInt(-2, 1));
        
        for (PairInt p : b) {
            assertNotNull(expected.remove(p));
        }
        assertTrue(expected.isEmpty());
        
        expected = new HashSet<PairInt>();        
        expected.add(new PairInt(-1, 1));
        expected.add(new PairInt(-2, 1));
        
        for (PairInt p : c) {
            expected.remove(p);
        }
        assertTrue(expected.isEmpty());
        
        expected = new HashSet<PairInt>();        
        expected.add(new PairInt(-1, 1));
        expected.add(new PairInt(-2, 1));
        
        for (PairInt p : d) {
            expected.remove(p);
        }
        assertTrue(expected.isEmpty());
    }
    
    public void testReverseY() throws Exception {
        
        Set<PairInt> a = new HashSet<PairInt>();
        a.add(new PairInt(1, 1));
        a.add(new PairInt(2, 1));
        
        Set<PairInt> b = new HashSet<PairInt>();
        b.add(new PairInt(1, 1));
        b.add(new PairInt(2, 1));
        
        Set<PairInt> expected = new HashSet<PairInt>();
        PairInt t1 = new PairInt(1, -1);
        PairInt t2 = new PairInt(2, -1);  
        expected.add(t1);
        expected.add(t2);
        
        PostLineThinnerCorrections pltc = new PostLineThinnerCorrections();
        pltc.reverseYs(a, b);
        
        expected = new HashSet<PairInt>();        
        expected.add(new PairInt(1, -1));
        expected.add(new PairInt(2, -1));
        for (PairInt p : a) {
            assertNotNull(expected.remove(p));
        }
        assertTrue(expected.isEmpty());
        
        expected = new HashSet<PairInt>();        
        expected.add(new PairInt(1, -1));
        expected.add(new PairInt(2, -1));
        for (PairInt p : b) {
            assertNotNull(expected.remove(p));
        }
        assertTrue(expected.isEmpty());
        
        // ---------------------------
        a = new HashSet<PairInt>();
        a.add(new PairInt(1, 1));
        a.add(new PairInt(2, 1));
        
        b = new HashSet<PairInt>();
        b.add(new PairInt(1, 1));
        b.add(new PairInt(2, 1));
        
        Set<PairInt> c = new HashSet<PairInt>();
        c.add(new PairInt(1, 1));
        c.add(new PairInt(2, 1));
        
        Set<PairInt> d = new HashSet<PairInt>();
        d.add(new PairInt(1, 1));
        d.add(new PairInt(2, 1));
        
        pltc.reverseYs(a, b, c, d);
        
        expected = new HashSet<PairInt>();        
        expected.add(new PairInt(1, -1));
        expected.add(new PairInt(2, -1));
        for (PairInt p : a) {
            assertNotNull(expected.remove(p));
        }
        assertTrue(expected.isEmpty());
        
        expected = new HashSet<PairInt>();        
        expected.add(new PairInt(1, -1));
        expected.add(new PairInt(2, -1));
        for (PairInt p : b) {
            assertNotNull(expected.remove(p));
        }
        assertTrue(expected.isEmpty());
        
        expected = new HashSet<PairInt>();        
        expected.add(new PairInt(1, -1));
        expected.add(new PairInt(2, -1));
        for (PairInt p : c) {
            assertNotNull(expected.remove(p));
        }
        assertTrue(expected.isEmpty());
        
        expected = new HashSet<PairInt>();        
        expected.add(new PairInt(1, -1));
        expected.add(new PairInt(2, -1));
        for (PairInt p : d) {
            assertNotNull(expected.remove(p));
        }
        assertTrue(expected.isEmpty());
    }
    
}
