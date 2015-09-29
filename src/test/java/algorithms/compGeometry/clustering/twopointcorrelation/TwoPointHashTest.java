package algorithms.compGeometry.clustering.twopointcorrelation;

import algorithms.compGeometry.clustering.twopointcorrelation.TwoPointHash.Node;
import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.List;
import static junit.framework.Assert.assertTrue;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class TwoPointHashTest extends TestCase {

    /**
     *
     * @param testName
     */
    public TwoPointHashTest(String testName) {
        super(testName);
    }

    /**
     *
     * @throws Exception
     */
    @Override
    protected void setUp() throws Exception {
        super.setUp();
    }

    /**
     *
     * @throws Exception
     */
    @Override
    protected void tearDown() throws Exception {
        super.tearDown();
    }

    /**
     * Test of storeIfDoesNotContain method, of class TwoPointIdentities.
     */
    public void testInsert() {
        int count = 0;
        ITwoPointIdentity identities = new TwoPointHash(121*121);
        for (int i = 1; i < 121; i++) {
            for (int j = (i + 1); j < 121; j++) {
                ((TwoPointHash)identities).insert(i, j);
                count++;
            }
        }
        boolean stored = false;
        for (int i = 1; i < 121; i++) {
            for (int j = (i + 1); j < 121; j++) {
                Node node = ((TwoPointHash)identities).search(i, j);
                assertTrue(node.a0 == i && node.a1 == j);
                stored = identities.storeIfDoesNotContain(i, j);
                assertFalse(stored);
            }
        }
        
        assertTrue(((TwoPointHash)identities).n == count);
        
    }

    /**
     *
     */
    public void testInsert2() {
        
        int count = 105000;
        
        int nRange = 10000;
        
        SecureRandom sr = new SecureRandom();
        long seed = System.currentTimeMillis();
        sr.setSeed(seed);
        
        List<Integer> c0 = new ArrayList<Integer>();
        List<Integer> c1 = new ArrayList<Integer>();
        
        ITwoPointIdentity identities = new TwoPointHash(count);
        for (int i = 0; i < count; i++) {
            
            int i0 = sr.nextInt(nRange);
            int i1 = sr.nextInt(nRange);
            
            while (((TwoPointHash)identities).search(i0, i1) != null) {
                i0 = sr.nextInt(nRange);
                i1 = sr.nextInt(nRange);
            }
            assertTrue(((TwoPointHash)identities).storeIfDoesNotContain(i0, i1));
            
            c0.add(Integer.valueOf(i0));
            c1.add(Integer.valueOf(i1));
        }
        
        assertTrue(((TwoPointHash)identities).n == count);
        
        for (int i = 0; i < count; i++) {
            int i0 = c0.get(i).intValue();
            int i1 = c1.get(i).intValue();
            
            assertNotNull(((TwoPointHash)identities).search(i0, i1));
        }
                
    }

    /**
     * Test of storeIfDoesNotContain method, of class TwoPointIdentities.
     */
    public void testStoreIfDoesNotContain() {

        int pointIndex0 = 0;
        int pointIndex1 = 10;

        ITwoPointIdentity identities = new TwoPointHash(121*121);
        boolean stored = identities.storeIfDoesNotContain(pointIndex0, pointIndex1);
        assertTrue(stored);

        stored = identities.storeIfDoesNotContain(pointIndex0, pointIndex1);
        assertFalse(stored);

        for (int i = 1; i < 121; i++) {
            for (int j = (i + 1); j < 121; j++) {
                stored = identities.storeIfDoesNotContain(i, j);
                assertTrue(stored);
            }
        }

        for (int i = 1; i < 121; i++) {
            for (int j = (i + 1); j < 121; j++) {
                stored = identities.storeIfDoesNotContain(i, j);
                assertFalse(stored);
            }
        }
    }

}
