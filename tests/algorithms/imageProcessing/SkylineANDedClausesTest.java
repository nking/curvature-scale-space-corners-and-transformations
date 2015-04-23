package algorithms.imageProcessing;

import junit.framework.TestCase;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class SkylineANDedClausesTest extends TestCase {
    
    public SkylineANDedClausesTest() {
    }

    @Test
    public void testGetForAllSkies() {
        
        SkylineANDedClauses instance = new SkylineANDedClauses();
        
        ANDedClauses[] a0 = instance.getForAllSkies();
        ANDedClauses[] a1 = instance.getForBlueSkies();
        ANDedClauses[] a2 = instance.getForRedSkies();
        
        int nTot = a0.length + a1.length + a2.length;
        
        ANDedClauses[] a = instance.getAllClauses();
        
        assertTrue(a.length == nTot);
        
        for (ANDedClauses ac : a) {
            assertNotNull(ac);
        }

        assertNotNull(instance.getFittableCoefficientsForAllSkies());
        assertNotNull(instance.getFittableCoefficientsForBlueSkies());
        assertNotNull(instance.getFittableCoefficientsForRedSkies());
    }

}
