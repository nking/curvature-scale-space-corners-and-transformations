package algorithms.imageProcessing.optimization;

import junit.framework.TestCase;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class SkylineANDedClauses2Test extends TestCase {
    
    public SkylineANDedClauses2Test() {
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
        
        SkylineANDedClauses2 instance = new SkylineANDedClauses2();
        
        for (int type = 0; type < 9; type++) {
            
            ANDedClauses[] a0 = null;
            
            switch(type) {
                case 0:
                    a0 = instance.getAllClauses();
                    break;
                case 1:
                    a0 = instance.getAllClausesLowerLimits();
                    break;
                case 2:
                    a0 = instance.getAllClausesUpperLimits();
                    break;
                case 3:
                    a0 = instance.getGeneralAndBlueClauses();
                    break;
                case 4:
                    a0 = instance.getGeneralAndBlueClausesLowerLimits();
                    break;
                case 5:
                    a0 = instance.getGeneralAndBlueClausesUpperLimits();
                    break;
                case 6:
                    a0 = instance.getGeneralAndRedClauses();
                    break;
                case 7:
                    a0 = instance.getGeneralAndRedClausesLowerLimits();
                    break;
                default:
                    a0 = instance.getGeneralAndRedClausesUpperLimits();
                    break;
            }
           
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
        }
        
    }
}
