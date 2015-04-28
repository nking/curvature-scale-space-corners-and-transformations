package algorithms.imageProcessing.optimization;

import junit.framework.TestCase;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class SkylineANDedClausesTest extends TestCase {
    
    public SkylineANDedClausesTest() {
    }

    public void testGetForAllSkies() {
        
        boolean skyIsRed = false;
        double pixContrast = 1;
        double pixBlueOrRedDiff = 2;
        double pixCIEXDiff = 3;
        double pixCIEYDiff = 4;
        double skyStDevContrast = 5;
        double skyStDevBlueOrRedDiff = 6;
        double skyStDevCIEX = 7;
        double skyStDevCIEY = 8;
        int r = 127;
        int g = 127;
        int b = 127;
        
        ColorData data = new ColorData(skyIsRed,
            pixContrast, pixBlueOrRedDiff, 
            pixCIEXDiff, pixCIEYDiff, 
            skyStDevContrast, skyStDevBlueOrRedDiff,
            skyStDevCIEX, skyStDevCIEY, r, g, b);
        
        SkylineANDedClauses instance = new SkylineANDedClauses();
        
        ANDedClauses[] a0 = instance.getForAllSkies();
        for (int i = 0; i < a0.length; i++) {
            ANDedClauses c = a0[i];
            assertNotNull(c.getSKYCONDITIONAL());
            for (int ii = 0; ii < c.n; ii++) {
                assertNotNull(c.getParams1(ii));
                assertNotNull(c.getParams2(ii));
                assertNotNull(c.getCoefficients(ii));
                assertNotNull(c.getGtOrLT(ii));
            }
            c.evaluate(data);
        }
        
        ANDedClauses[] a1 = instance.getForBlueSkies();
        for (int i = 0; i < a1.length; i++) {
            ANDedClauses c = a1[i];
            assertNotNull(c.getSKYCONDITIONAL());
            for (int ii = 0; ii < c.n; ii++) {
                assertNotNull(c.getParams1(ii));
                assertNotNull(c.getParams2(ii));
                assertNotNull(c.getCoefficients(ii));
                assertNotNull(c.getGtOrLT(ii));
            }
            c.evaluate(data);
        }
        
        ANDedClauses[] a2 = instance.getForRedSkies();
        for (int i = 0; i < a2.length; i++) {
            ANDedClauses c = a2[i];
            assertNotNull(c.getSKYCONDITIONAL());
            for (int ii = 0; ii < c.n; ii++) {
                assertNotNull(c.getParams1(ii));
                assertNotNull(c.getParams2(ii));
                assertNotNull(c.getCoefficients(ii));
                assertNotNull(c.getGtOrLT(ii));
            }
            c.evaluate(data);
        }
        
        int nTot = a0.length + a1.length + a2.length;
        
        ANDedClauses[] a = instance.getAllClauses();
        
        assertTrue(a.length == nTot);
        
        for (ANDedClauses ac : a) {
            assertNotNull(ac);
            assertNotNull(ac.getSKYCONDITIONAL());
            ac.evaluate(data);
        }

        assertNotNull(instance.getFittableCoefficientsForAllSkies());
        assertNotNull(instance.getFittableCoefficientsForBlueSkies());
        assertNotNull(instance.getFittableCoefficientsForRedSkies());
        
        boolean useBlueSkyImages = true;
        for (int idx = 0; idx < 2; idx++) {
            if (idx == 1) {
                useBlueSkyImages = false;
            }
            ANDedClauses[] clauses;
            float[][] coeffLowerLimits;
            float[][] coeffUpperLimits;
            if (useBlueSkyImages) {
                clauses = instance.getAllAndBlueClauses();
                coeffLowerLimits = instance.getAllAndBlueCoeffLowerLimits();
                coeffUpperLimits = instance.getAllAndBlueCoeffUpperLimits();
            } else {
                clauses = instance.getAllAndRedClauses();
                coeffLowerLimits = instance.getAllAndRedCoeffLowerLimits();
                coeffUpperLimits = instance.getAllAndRedCoeffUpperLimits();
            }
            assertTrue(clauses.length == coeffLowerLimits.length);
            assertTrue(clauses.length == coeffUpperLimits.length);
            for (int clauseIdx = 0; clauseIdx < clauses.length; clauseIdx++) {
                float[] c0 = clauses[clauseIdx].coefficients;
                float[] c1 = coeffLowerLimits[clauseIdx];
                float[] c2 = coeffUpperLimits[clauseIdx];//4
                assertTrue(c0.length == c1.length);
                assertTrue(c0.length == c2.length);
            }
        }
        
    }

}
