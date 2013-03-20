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

    public void testCreateIdentity() {

        int i0 = Integer.MAX_VALUE; // 2^31 - 1 = 0111 1111  1111 1111  1111 1111  1111 1111
                                    // reversed 1111 1111  1111 1111  1111 1111 0111 1111
        byte[] expected0 = new byte[]{-1,-1,-1, 127};
        byte[] bytes0 = StringArrayLite.writeIntegerToBytes(i0, 4);

        int i1 = Integer.MIN_VALUE; // -2^31    = 1000 0000  0000 0000  0000 0000  0000 0000
                                    // reversed   0000 0000  0000 0000  0000 0000  1000 0000
        byte[] expected1 = new byte[]{0,0,0,-128};
        byte[] bytes1 = StringArrayLite.writeIntegerToBytes(i1, 4);

        assertTrue(bytes0.length == expected0.length);
        for (int i = 0; i < expected0.length; i++) {
            assertTrue(expected0[i] == bytes0[i]);
        }
        assertTrue(bytes1.length == expected1.length);
        for (int i = 0; i < expected1.length; i++) {
            assertTrue(expected1[i] == bytes1[i]);
        }

        byte[] identity = StringArrayLite.createIdentity(i0, i1);
        byte[] expected = new byte[]{-1,-1,-1,127, 0,0,0,-128};

        assertTrue(identity.length == expected.length);

        for (int i = 0; i < expected.length; i++) {
            assertTrue(expected[i] == identity[i]);
        }

        i0 = 123; // = 0111 1011  0000 0000  0000 0000  0000 0000
        expected0 = new byte[]{123,0,0,0};
        bytes0 = StringArrayLite.writeIntegerToBytes(i0, 4);
        assertTrue(bytes0.length == expected0.length);
        for (int i = 0; i < expected0.length; i++) {
            assertTrue(expected0[i] == bytes0[i]);
        }

        i1 = 3000; // = 0000 1011  1011 1000  0000 0000  0000 0000
        expected1 = new byte[]{-72,11,0,0};
        bytes1 = StringArrayLite.writeIntegerToBytes(i1, 4);
        assertTrue(bytes1.length == expected1.length);
        for (int i = 0; i < expected1.length; i++) {
            assertTrue(expected1[i] == bytes1[i]);
        }

        identity = StringArrayLite.createIdentity(i0, i1);
        expected = new byte[]{123,0,0,0, -72,11,0,0};

        assertTrue(identity.length == expected.length);

        for (int i = 0; i < expected.length; i++) {
            assertTrue(expected[i] == identity[i]);
        }
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

        byte[] c0 = StringArrayLite.createIdentity(1023456789, 1023456798);
        byte[] c1 = StringArrayLite.createIdentity(987654321, 897654321);

        StringLite s0 = new StringLite(c0);
        StringLite s1 = new StringLite(c1);

        // same content of strings, but in different order produces different hash
        assertTrue( s0.hashCode() != s1.hashCode());

        StringLite s00 = new StringLite(c0);
        StringLite s10 = new StringLite(c1);
        // same content of strings in same order produces same hash
        assertTrue( s00.hashCode() == s0.hashCode());
        assertTrue( s10.hashCode() == s1.hashCode());

        c0 = StringArrayLite.createIdentity(1023456789, 10234);
        c1 = StringArrayLite.createIdentity(987654321, 21034);
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
