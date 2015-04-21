package algorithms.imageProcessing;

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
        
        ColorData data = new ColorData(skyIsRed,
            pixContrast, pixBlueOrRedDiff, 
            pixCIEXDiff, pixCIEYDiff, 
            skyStDevContrast, skyStDevBlueOrRedDiff,
            skyStDevCIEX, skyStDevCIEY);
        
        assertFalse(data.skyIsRed());
        assertTrue(data.getParameter(PARAM.ABSOLUTE_CONTRAST) == pixContrast);
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
    }
}
