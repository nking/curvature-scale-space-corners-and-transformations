package algorithms.imageProcessing.optimization;

import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class ColorDataTest extends TestCase {
    
    public ColorDataTest() {
    }
    
    public void testGetParameter() throws Exception {
        
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
        
        assertFalse(data.skyIsRed());
        assertTrue(data.getParameter(PARAM.ABSOLUTE_CONTRAST) == pixContrast);
        assertTrue(data.getParameter(PARAM.CONTRAST) == pixContrast);
        assertTrue(data.getParameter(PARAM.ABSOLUTE_DIFF_BLUE_OR_RED) == 
            pixBlueOrRedDiff);
        assertTrue(data.getParameter(PARAM.DIFF_CIEX) == pixCIEXDiff);
        assertTrue(data.getParameter(PARAM.DIFF_CIEY) == pixCIEYDiff);
        assertTrue(data.getParameter(PARAM.STDEV_CONTRAST) == 
            skyStDevContrast);
        assertTrue(data.getParameter(PARAM.STDEV_BLUE_OR_RED) == 
            skyStDevBlueOrRedDiff);
        assertTrue(data.getParameter(PARAM.STDEV_CIEX) == skyStDevCIEX);
        assertTrue(data.getParameter(PARAM.STDEV_CIEY) == skyStDevCIEY);
        
        assertTrue(data.getParameter(PARAM.INT_ONE) == 1);
        
        assertTrue(data.getParameter(PARAM.RED) == r);
        assertTrue(data.getParameter(PARAM.GREEN) == g);
        assertTrue(data.getParameter(PARAM.BLUE) == b);
         
        assertTrue(Math.abs(data.getParameter(PARAM.R_DIV_TOT) - (1./3.)) < 0.01);
        assertTrue(Math.abs(data.getParameter(PARAM.G_DIV_TOT) - (1./3.)) < 0.01);
        assertTrue(Math.abs(data.getParameter(PARAM.B_DIV_TOT) - (1./3.)) < 0.01);
        
        assertTrue(data.getParameter(PARAM.DIFF_R_DIV_TOT_ONE_THIRD) < 0.01);
        assertTrue(data.getParameter(PARAM.DIFF_G_DIV_TOT_ONE_THIRD) < 0.01);
        assertTrue(data.getParameter(PARAM.DIFF_B_DIV_TOT_ONE_THIRD) < 0.01);
    }
}
