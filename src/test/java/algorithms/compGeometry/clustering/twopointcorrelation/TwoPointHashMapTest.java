package algorithms.compGeometry.clustering.twopointcorrelation;

import static junit.framework.Assert.assertTrue;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class TwoPointHashMapTest extends TestCase {

    public TwoPointHashMapTest(String testName) {
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

    public void testHash() {

        int nData = 100;

        TwoPointHashMap identities = new TwoPointHashMap(100);

        int h = identities.hash(0, 0);
        assertTrue(h == 0);

        h = identities.hash(0, 1);
        assertTrue(h == 1);

        h = identities.hash(1, 0);
        assertTrue(h == nData);

        h = identities.hash(nData-1, nData-1);
        assertTrue(h == (nData*nData)-1);

        // 0, 0 should be index 0
        // 0, 1 should be index 1
        // 0, nMax-1 should be index nMax-1
        // 1, 0 should be index nMax
    }

    /**
     * Test of storeIfDoesNotContain method, of class TwoPointIdentities.
     */
    public void testStoreIfDoesNotContain() {

        int pointIndex0 = 0;
        int pointIndex1 = 10;

        int nData = 390;

        ITwoPointIdentity identities = new TwoPointHashMap(nData);
        boolean stored = identities.storeIfDoesNotContain(pointIndex0, pointIndex1);
        assertTrue(stored);

        stored = identities.storeIfDoesNotContain(pointIndex0, pointIndex1);
        assertFalse(stored);

        for (int i = 1; i < nData; i++) {
            for (int j = (i + 1); j < 11; j++) {
                stored = identities.storeIfDoesNotContain(i, j);
                assertTrue(stored);
            }
        }

        for (int i = 1; i < nData; i++) {
            for (int j = (i + 1); j < 11; j++) {
                stored = identities.storeIfDoesNotContain(i, j);
                assertFalse(stored);
            }
        }
    }
    
    public void testConstructor() {

        int nData = 46350;

        boolean threwException = false;
        
        try {
            TwoPointHashMap identities = new TwoPointHashMap(nData);
        } catch(IllegalArgumentException e) {
            threwException = true;
        }
        assertTrue(threwException);
        
    }

}
