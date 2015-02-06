package algorithms.imageProcessing;

import algorithms.compGeometry.PointInPolygon;
import algorithms.util.ArrayPair;
import algorithms.util.ResourceFinder;
import java.util.logging.Logger;
import junit.framework.TestCase;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class CIEChromaticityTest extends TestCase {
    
    private Logger log = Logger.getLogger(this.getClass().getName());
    
    public CIEChromaticityTest(String testName) {
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
    
    public void test() throws Exception {
    
        CIEChromaticity cieC = new CIEChromaticity();
        
        String filePath1 = ResourceFinder.findFileInTestResources(
            "sky_with_rainbow.jpg");
        Image img1 = ImageIOHelper.readImage(filePath1);
        int image1Width = img1.getWidth();
        int image1Height = img1.getHeight();

        int col = 332;
        int row = 179;
        int pixR = img1.getR(col, row);
        int pixG = img1.getG(col, row);
        int pixB = img1.getB(col, row);
        
        float[] pixCIEXY = cieC.rgbToXYChromaticity(pixR, pixG, pixB);
        
        ArrayPair yellowBounds = cieC.getYellowPolynomial();
        ArrayPair yellowOrangeBounds = cieC.getYellowishGreenThroughOrangePolynomial();
        ArrayPair redBounds = cieC.getRedPolynomial();
        
        PointInPolygon pInPoly = new PointInPolygon();
        
        ArrayPair t = yellowBounds;
        
        boolean isYellow = pInPoly.isInSimpleCurve(pixCIEXY[0], pixCIEXY[1], 
            t.getX(), t.getY(), t.getX().length);
        assertTrue(isYellow);
        
        t = yellowOrangeBounds;
        boolean isInYellowOrange = pInPoly.isInSimpleCurve(pixCIEXY[0], pixCIEXY[1], 
            t.getX(), t.getY(), t.getX().length);
        assertTrue(isInYellowOrange);
        
        t = redBounds;
        boolean isRed = pInPoly.isInSimpleCurve(pixCIEXY[0], pixCIEXY[1], 
            t.getX(), t.getY(), t.getX().length);
        assertFalse(isRed);
        
        col = 254;
        row = 170;
        pixR = img1.getR(col, row);
        pixG = img1.getG(col, row);
        pixB = img1.getB(col, row);
        pixCIEXY = cieC.rgbToXYChromaticity(pixR, pixG, pixB);
        ArrayPair purpleBounds = cieC.getPurplePolynomial();
        t = purpleBounds;
        boolean isPurple = pInPoly.isInSimpleCurve(pixCIEXY[0], pixCIEXY[1], 
            t.getX(), t.getY(), t.getX().length);
        assertTrue(isPurple);
        
        t = yellowBounds;
        isYellow = pInPoly.isInSimpleCurve(pixCIEXY[0], pixCIEXY[1], 
            t.getX(), t.getY(), t.getX().length);
        assertFalse(isYellow);
    }
    
}
