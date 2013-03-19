package algorithms.compGeometry.clustering.twopointcorrelation;

import static junit.framework.Assert.assertTrue;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class StringArrayLiteTest extends TestCase {

    public StringArrayLiteTest(String testName) {
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

    /**
     * Test of storeIfDoesNotContain method, of class TwoPointIdentities.
     */
    public void testStoreIfDoesNotContain() {

        int pointIndex0 = 0;
        int pointIndex1 = 10;

        StringArrayLite identities = new StringArrayLite();
        boolean stored = identities.storeIfDoesNotContain(pointIndex0, pointIndex1);
        assertTrue(stored);

        stored = identities.storeIfDoesNotContain(pointIndex0, pointIndex1);
        assertFalse(stored);

        for (int i = 1; i < 120; i++) {
            stored = identities.storeIfDoesNotContain(i, pointIndex1);
            assertTrue(stored);
        }
    }

    public void testHashCode() {

        char[] c0 = "1023456789 1023456798".toCharArray();

        char[] c1 = "9876543210 8976543201".toCharArray();

        StringLite s0 = new StringLite(c0);
        StringLite s1 = new StringLite(c1);

        // same content of strings, but in different order produces different hash
        assertTrue( s0.hashCode() != s1.hashCode());

        StringLite s00 = new StringLite(c0);
        StringLite s10 = new StringLite(c1);
        // same content of strings in same order produces same hash
        assertTrue( s00.hashCode() == s0.hashCode());
        assertTrue( s10.hashCode() == s1.hashCode());


        c0 = "1023456789 10234".toCharArray();
        c1 = "9876543210 21034".toCharArray();
        s0 = new StringLite(c0);
        s1 = new StringLite(c1);
        // same content of strings, but in different order produces different hash
        assertTrue( s0.hashCode() != s1.hashCode());

        s00 = new StringLite(c0);
        s10 = new StringLite(c1);
        // same content of strings in same order produces same hash
        assertTrue( s00.hashCode() == s0.hashCode());
        assertTrue( s10.hashCode() == s1.hashCode());
    }
}
