package algorithms.util;

import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;
import java.util.HashSet;
import java.util.Set;
import junit.framework.TestCase;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class BitVectorRepresentationTest extends TestCase {
    
    public BitVectorRepresentationTest() {
    }

    public void testCreateBitstring() throws Exception {
        
        Set<PairInt> set = getSet1();
        
        Set<PairInt> subset0 = new HashSet<PairInt>(set);
        
        BitVectorRepresentation instance = new BitVectorRepresentation(set);
        
        Set<VeryLongBitString> bitstrings = new HashSet<VeryLongBitString>();
        
        while (!subset0.isEmpty()) {
        
            VeryLongBitString result = instance.createBitstring(subset0);

            assertFalse(bitstrings.contains(result));
        
            assertTrue(bitstrings.add(result));
            
            assertTrue(subset0.remove(subset0.iterator().next()));
        }
    }

    private Set<PairInt> getSet1() throws NoSuchAlgorithmException {
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        System.out.println("SEED=" + seed);
        sr.setSeed(seed);
        
        Set<PairInt> set = new HashSet<PairInt>();
        
        for (int i = 0; i < 1000; ++i) {
            int x = sr.nextInt((Integer.MAX_VALUE >> 1) - 1);
            int y = sr.nextInt((Integer.MAX_VALUE >> 1) - 1);;
            set.add(new PairInt(x, y));
        }
        
        return set;
    }
}
