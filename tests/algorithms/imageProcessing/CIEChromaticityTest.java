package algorithms.imageProcessing;

import algorithms.compGeometry.PointInPolygon;
import algorithms.util.ArrayPair;
import algorithms.util.ResourceFinder;
import java.awt.Color;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
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
        ArrayPair yellowOrangeBounds = cieC.getYellowThroughOrangePolynomial();
        ArrayPair redBounds = cieC.getRedPolynomial();
        
        PointInPolygon pInPoly = new PointInPolygon();
        
        ArrayPair t = yellowOrangeBounds;
        boolean isInYellowOrange = pInPoly.isInSimpleCurve(pixCIEXY[0], pixCIEXY[1], 
            t.getX(), t.getY(), t.getX().length);
        //assertTrue(isInYellowOrange);
        
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
        boolean isYellow = pInPoly.isInSimpleCurve(pixCIEXY[0], pixCIEXY[1], 
            t.getX(), t.getY(), t.getX().length);
        assertFalse(isYellow);
       
    }
    
    public void test2() throws Exception {
    
        CIEChromaticity cieC = new CIEChromaticity();
        
        String filePath1 = ResourceFinder.findFileInTestResources(
            //"sky_with_rainbow.jpg"
            "arizona-sunrise-1342919937GHz.jpg"
        );
        Image img1 = ImageIOHelper.readImage(filePath1);
        int image1Width = img1.getWidth();
        int image1Height = img1.getHeight();

        int col = 83;
        int row = 89;
        int pixR = img1.getR(col, row);
        int pixG = img1.getG(col, row);
        int pixB = img1.getB(col, row);
        float[] hsb = new float[3];
        Color.RGBtoHSB(pixR, pixG, pixB, hsb);
        float[] pixCIEXY = cieC.rgbToXYChromaticity(pixR, pixG, pixB);
        log.fine(String.format(
            "(%d,%d)  rgb=(%d,%d,%d) pixCIEXY=%f,%f  hsb=(%f,%f,%f)", 
            col, row, pixR, pixG, pixB,
            pixCIEXY[0], pixCIEXY[1],
            hsb[0], hsb[1], hsb[2]));
        
        col = 88;
        row = 84;
        pixR = img1.getR(col, row);
        pixG = img1.getG(col, row);
        pixB = img1.getB(col, row);
        hsb = new float[3];
        Color.RGBtoHSB(pixR, pixG, pixB, hsb);
        pixCIEXY = cieC.rgbToXYChromaticity(pixR, pixG, pixB);
        log.fine(String.format(
            "(%d,%d)  rgb=(%d,%d,%d) pixCIEXY=%f,%f  hsb=(%f,%f,%f)", 
            col, row, pixR, pixG, pixB,
            pixCIEXY[0], pixCIEXY[1],
            hsb[0], hsb[1], hsb[2]));
        
        col = 93;
        row = 88;
        pixR = img1.getR(col, row);
        pixG = img1.getG(col, row);
        pixB = img1.getB(col, row);
        hsb = new float[3];
        Color.RGBtoHSB(pixR, pixG, pixB, hsb);
        pixCIEXY = cieC.rgbToXYChromaticity(pixR, pixG, pixB);
        log.fine(String.format(
            "(%d,%d)  rgb=(%d,%d,%d) pixCIEXY=%f,%f  hsb=(%f,%f,%f)", 
            col, row, pixR, pixG, pixB,
            pixCIEXY[0], pixCIEXY[1],
            hsb[0], hsb[1], hsb[2]));
        
        col = 525;
        row = 291;
        pixR = img1.getR(col, row);
        pixG = img1.getG(col, row);
        pixB = img1.getB(col, row);
        hsb = new float[3];
        Color.RGBtoHSB(pixR, pixG, pixB, hsb);
        pixCIEXY = cieC.rgbToXYChromaticity(pixR, pixG, pixB);
        log.fine(String.format(
            "(%d,%d)  rgb=(%d,%d,%d) pixCIEXY=%f,%f  hsb=(%f,%f,%f)", 
            col, row, pixR, pixG, pixB,
            pixCIEXY[0], pixCIEXY[1],
            hsb[0], hsb[1], hsb[2]));
        
        
        ArrayPair redBounds = cieC.getRedPolynomial();
        ArrayPair orangeBounds = cieC.getOrangePolynomial();
        ArrayPair yellowBounds = cieC.getYellowPolynomial();
        ArrayPair yellowOrangeBounds = cieC.getYellowThroughOrangePolynomial();
        ArrayPair greenishYellowOrangeBounds = cieC.getYellowThroughOrangePolynomial();
        
        PointInPolygon pInPoly = new PointInPolygon();
        
        float yellowTotRGBLimit = 0.3f * 768.f;
        
        float redTotRGBLimit = 0.2f * 768.f;
        
        float orangeTotRGBLimit = 0.2f * 768.f;
        
        float orangeSLimitForBrightSkies = 0.2f;
        
        float orangeSLimitForDarkSkies = 0.3f;
        
        float yellowSLimitForBrightSkies = 0.1f;
        
        float yellowSLimitForDarkSkies = 0.175f;
        
        float yellowOrangeSLimitForBrightSkies = 0.1f;
        
        float yellowOrangeSLimitForDarkSkies = 0.175f;
        
        float redSLimitForBrightSkies = 0.05f;
        
        float redSLimitForDarkSkies = 0.175f;
        
        float redPurpleSLimitForBrightSkies = 0.05f;
        
        float redPurpleSLimitForDarkSkies = 0.175f;
        
        for (col = 0; col < img1.getWidth(); col++) {
            for (row = 0; row < img1.getHeight(); row++) {
                pixR = img1.getR(col, row);
                pixG = img1.getG(col, row);
                pixB = img1.getB(col, row);
                pixCIEXY = cieC.rgbToXYChromaticity(pixR, pixG, pixB);
                
                int totRGB = pixR + pixG + pixB;
                
                hsb = new float[3];
                Color.RGBtoHSB(pixR, pixG, pixB, hsb);

                ArrayPair t = cieC.getRedThroughPurplishRedPolynomial();
                
                if ((pixR > 15) && (pixG > 15) && (pixB > 15)) {
                    
                    if (
                        //(hsb[1] > redPurpleSLimitForDarkSkies) 
                        (hsb[1] > redPurpleSLimitForBrightSkies)
                        && (hsb[2] > 0.25)
                        && (totRGB > redTotRGBLimit)
                        && (pInPoly.isInSimpleCurve(pixCIEXY[0], pixCIEXY[1],
                        t.getX(), t.getY(), t.getX().length))) {

                        img1.setRGB(col, row, 0, 0, 255);

                        log.fine(String.format(
                            "(%d,%d)  rgb=(%d,%d,%d) pixCIEXY=%f,%f  hsb=(%f,%f,%f)", 
                            col, row, pixR, pixG, pixB, 
                            pixCIEXY[0], pixCIEXY[1],
                            hsb[0], hsb[1], hsb[2]));
                      
                    } else if (
                        //(hsb[1] > orangeSLimitForDarkSkies) 
                        (hsb[1] > orangeSLimitForBrightSkies)
                        && (hsb[2] > 0.25)
                        && (totRGB > orangeTotRGBLimit)
                        && (pInPoly.isInSimpleCurve(pixCIEXY[0], pixCIEXY[1],
                        orangeBounds.getX(), orangeBounds.getY(), 
                        orangeBounds.getX().length))) {

                        img1.setRGB(col, row, 0, 0, 255);

                        log.fine(String.format(
                            "(%d,%d)  rgb=(%d,%d,%d) pixCIEXY=%f,%f  hsb=(%f,%f,%f)", 
                            col, row, pixR, pixG, pixB, 
                            pixCIEXY[0], pixCIEXY[1],
                            hsb[0], hsb[1], hsb[2]));
                    
                    } else if (
                        //(hsb[1] > yellowOrangeSLimitForDarkSkies) 
                        (hsb[1] > yellowOrangeSLimitForBrightSkies) 
                        && (hsb[2] > 0.25)
                        && (totRGB > yellowTotRGBLimit)
                        && (pInPoly.isInSimpleCurve(pixCIEXY[0], pixCIEXY[1],
                        greenishYellowOrangeBounds.getX(), 
                        greenishYellowOrangeBounds.getY(), 
                        greenishYellowOrangeBounds.getX().length)) ) {

                        img1.setRGB(col, row, 0, 0, 255);

                        log.fine(String.format(
                            "(%d,%d)  rgb=(%d,%d,%d) pixCIEXY=%f,%f  hsb=(%f,%f,%f)", 
                            col, row, pixR, pixG, pixB, 
                            pixCIEXY[0], pixCIEXY[1],
                            hsb[0], hsb[1], hsb[2]));
                    }
                }
            }
        }
        
        String dirPath = ResourceFinder.findDirectory("bin");
        String filePath2 = dirPath + "/tmp.png";
            
        //ImageDisplayer.displayImage("points", img1);
        ImageIOHelper.writeOutputImage(filePath2, img1);
        
    }
    
    public void testDeltaE2000() throws Exception {
        
        int r1, g1, b1, r2, g2, b3;
        float[] lab1, lab2;
        double deltaE;
        
        CIEChromaticity cieC = new CIEChromaticity();
        
        lab1 = cieC.rgbToCIELAB(0, 0, 0);
        lab2 = cieC.rgbToCIELAB(0, 0, 0);
        deltaE = cieC.calcDeltaECIE2000(lab1, lab2);
        assertEquals(0., deltaE);
        
        lab1 = cieC.rgbToCIELAB(0, 0, 0);
        lab2 = cieC.rgbToCIELAB(255, 255, 255);
        deltaE = cieC.calcDeltaECIE2000(lab1, lab2);//19.21
        
        lab1 = cieC.rgbToCIELAB(255, 255, 255);
        lab2 = cieC.rgbToCIELAB(255, 255, 255);
        deltaE = cieC.calcDeltaECIE2000(lab1, lab2);
        assertEquals(0., deltaE);        
    }
    
    public void testDeltaE2000_1() throws Exception {
        
        /*
        test data from
        http://www.ece.rochester.edu/~gsharma/ciede2000/dataNprograms/ciede2000testdata.txt
        % Each row of the file has a pair of CIELAB values and a color difference.

        For each row, columns 1-3 of the file correspond to CIE L*, a*, b* values,
        respectively for a reference color; 
        columns 4-6 correspond to L*, a*, b* values, respectively for a sample color; 
        and column 7 corresponds to the CIEDE2000 color difference between this pair of values.
        
        */
        
        List<Float> ell1 = new ArrayList<Float>();
        List<Float> a1 = new ArrayList<Float>();
        List<Float> b1 = new ArrayList<Float>();
        
        List<Float> ell2 = new ArrayList<Float>();
        List<Float> a2 = new ArrayList<Float>();
        List<Float> b2 = new ArrayList<Float>();
        
        List<Double> deltaE = new ArrayList<Double>();
        
        populateWithTestData(ell1, a1, b1, ell2, a2, b2, deltaE);
                
        CIEChromaticity cieC = new CIEChromaticity();
        
        for (int i = 0; i < ell1.size(); ++i) {
            
            double dE = cieC.calcDeltaECIE2000(ell1.get(i), a1.get(i), b1.get(i),
                ell2.get(i), a2.get(i), b2.get(i));
            
            double diffE = Math.abs(dE - deltaE.get(i).doubleValue());
        
            //log.info(i + ") dE=" + dE + " expected=" + deltaE.get(i) + " --> " + diffE);
            
            assertTrue(diffE < 1E-2);
        }
    }

    private void populateWithTestData(List<Float> ell1, List<Float> a1, List<Float> b1, 
        List<Float> ell2, List<Float> a2, List<Float> b2, List<Double> deltaE) {
        
        ell1.add(Float.valueOf(50.0000f)); a1.add(Float.valueOf(2.6772f)); b1.add(Float.valueOf(-79.7751f)); ell2.add(Float.valueOf(50.0000f)); a2.add(Float.valueOf(0.0000f)); b2.add(Float.valueOf(-82.7485f)); deltaE.add(Double.valueOf(2.0425));
        ell1.add(Float.valueOf(50.0000f)); a1.add(Float.valueOf(3.1571f)); b1.add(Float.valueOf(-77.2803f)); ell2.add(Float.valueOf(50.0000f)); a2.add(Float.valueOf(0.0000f)); b2.add(Float.valueOf(-82.7485f)); deltaE.add(Double.valueOf(2.8615));
        ell1.add(Float.valueOf(50.0000f)); a1.add(Float.valueOf(2.8361f)); b1.add(Float.valueOf(-74.0200f)); ell2.add(Float.valueOf(50.0000f)); a2.add(Float.valueOf(0.0000f)); b2.add(Float.valueOf(-82.7485f)); deltaE.add(Double.valueOf(3.4412));
        ell1.add(Float.valueOf(50.0000f)); a1.add(Float.valueOf(-1.3802f)); b1.add(Float.valueOf(-84.2814f)); ell2.add(Float.valueOf(50.0000f)); a2.add(Float.valueOf(0.0000f)); b2.add(Float.valueOf(-82.7485f)); deltaE.add(Double.valueOf(1.0000));
        ell1.add(Float.valueOf(50.0000f)); a1.add(Float.valueOf(-1.1848f)); b1.add(Float.valueOf(-84.8006f)); ell2.add(Float.valueOf(50.0000f)); a2.add(Float.valueOf(0.0000f)); b2.add(Float.valueOf(-82.7485f)); deltaE.add(Double.valueOf(1.0000));
        ell1.add(Float.valueOf(50.0000f)); a1.add(Float.valueOf(-0.9009f)); b1.add(Float.valueOf(-85.5211f)); ell2.add(Float.valueOf(50.0000f)); a2.add(Float.valueOf(0.0000f)); b2.add(Float.valueOf(-82.7485f)); deltaE.add(Double.valueOf(1.0000));
        ell1.add(Float.valueOf(50.0000f)); a1.add(Float.valueOf(0.0000f)); b1.add(Float.valueOf(0.0000f)); ell2.add(Float.valueOf(50.0000f)); a2.add(Float.valueOf(-1.0000f)); b2.add(Float.valueOf(2.0000f)); deltaE.add(Double.valueOf(2.3669));
        ell1.add(Float.valueOf(50.0000f)); a1.add(Float.valueOf(-1.0000f)); b1.add(Float.valueOf(2.0000f)); ell2.add(Float.valueOf(50.0000f)); a2.add(Float.valueOf(0.0000f)); b2.add(Float.valueOf(0.0000f)); deltaE.add(Double.valueOf(2.3669));
        ell1.add(Float.valueOf(50.0000f)); a1.add(Float.valueOf(2.4900f)); b1.add(Float.valueOf(-0.0010f)); ell2.add(Float.valueOf(50.0000f)); a2.add(Float.valueOf(-2.4900f)); b2.add(Float.valueOf(0.0009f)); deltaE.add(Double.valueOf(7.1792));
        ell1.add(Float.valueOf(50.0000f)); a1.add(Float.valueOf(2.4900f)); b1.add(Float.valueOf(-0.0010f)); ell2.add(Float.valueOf(50.0000f)); a2.add(Float.valueOf(-2.4900f)); b2.add(Float.valueOf(0.0010f)); deltaE.add(Double.valueOf(7.1792));
        ell1.add(Float.valueOf(50.0000f)); a1.add(Float.valueOf(2.4900f)); b1.add(Float.valueOf(-0.0010f)); ell2.add(Float.valueOf(50.0000f)); a2.add(Float.valueOf(-2.4900f)); b2.add(Float.valueOf(0.0011f)); deltaE.add(Double.valueOf(7.2195));
        ell1.add(Float.valueOf(50.0000f)); a1.add(Float.valueOf(2.4900f)); b1.add(Float.valueOf(-0.0010f)); ell2.add(Float.valueOf(50.0000f)); a2.add(Float.valueOf(-2.4900f)); b2.add(Float.valueOf(0.0012f)); deltaE.add(Double.valueOf(7.2195));
        ell1.add(Float.valueOf(50.0000f)); a1.add(Float.valueOf(-0.0010f)); b1.add(Float.valueOf(2.4900f)); ell2.add(Float.valueOf(50.0000f)); a2.add(Float.valueOf(0.0009f)); b2.add(Float.valueOf(-2.4900f)); deltaE.add(Double.valueOf(4.8045));
        ell1.add(Float.valueOf(50.0000f)); a1.add(Float.valueOf(-0.0010f)); b1.add(Float.valueOf(2.4900f)); ell2.add(Float.valueOf(50.0000f)); a2.add(Float.valueOf(0.0010f)); b2.add(Float.valueOf(-2.4900f)); deltaE.add(Double.valueOf(4.8045));
        ell1.add(Float.valueOf(50.0000f)); a1.add(Float.valueOf(-0.0010f)); b1.add(Float.valueOf(2.4900f)); ell2.add(Float.valueOf(50.0000f)); a2.add(Float.valueOf(0.0011f)); b2.add(Float.valueOf(-2.4900f)); deltaE.add(Double.valueOf(4.7461));
        ell1.add(Float.valueOf(50.0000f)); a1.add(Float.valueOf(2.5000f)); b1.add(Float.valueOf(0.0000f)); ell2.add(Float.valueOf(50.0000f)); a2.add(Float.valueOf(0.0000f)); b2.add(Float.valueOf(-2.5000f)); deltaE.add(Double.valueOf(4.3065));
        ell1.add(Float.valueOf(50.0000f)); a1.add(Float.valueOf(2.5000f)); b1.add(Float.valueOf(0.0000f)); ell2.add(Float.valueOf(73.0000f)); a2.add(Float.valueOf(25.0000f)); b2.add(Float.valueOf(-18.0000f)); deltaE.add(Double.valueOf(27.1492));
        ell1.add(Float.valueOf(50.0000f)); a1.add(Float.valueOf(2.5000f)); b1.add(Float.valueOf(0.0000f)); ell2.add(Float.valueOf(61.0000f)); a2.add(Float.valueOf(-5.0000f)); b2.add(Float.valueOf(29.0000f)); deltaE.add(Double.valueOf(22.8977));
        ell1.add(Float.valueOf(50.0000f)); a1.add(Float.valueOf(2.5000f)); b1.add(Float.valueOf(0.0000f)); ell2.add(Float.valueOf(56.0000f)); a2.add(Float.valueOf(-27.0000f)); b2.add(Float.valueOf(-3.0000f)); deltaE.add(Double.valueOf(31.9030));
        ell1.add(Float.valueOf(50.0000f)); a1.add(Float.valueOf(2.5000f)); b1.add(Float.valueOf(0.0000f)); ell2.add(Float.valueOf(58.0000f)); a2.add(Float.valueOf(24.0000f)); b2.add(Float.valueOf(15.0000f)); deltaE.add(Double.valueOf(19.4535));
        ell1.add(Float.valueOf(50.0000f)); a1.add(Float.valueOf(2.5000f)); b1.add(Float.valueOf(0.0000f)); ell2.add(Float.valueOf(50.0000f)); a2.add(Float.valueOf(3.1736f)); b2.add(Float.valueOf(0.5854f)); deltaE.add(Double.valueOf(1.0000));
        ell1.add(Float.valueOf(50.0000f)); a1.add(Float.valueOf(2.5000f)); b1.add(Float.valueOf(0.0000f)); ell2.add(Float.valueOf(50.0000f)); a2.add(Float.valueOf(3.2972f)); b2.add(Float.valueOf(0.0000f)); deltaE.add(Double.valueOf(1.0000));
        ell1.add(Float.valueOf(50.0000f)); a1.add(Float.valueOf(2.5000f)); b1.add(Float.valueOf(0.0000f)); ell2.add(Float.valueOf(50.0000f)); a2.add(Float.valueOf(1.8634f)); b2.add(Float.valueOf(0.5757f)); deltaE.add(Double.valueOf(1.0000));
        ell1.add(Float.valueOf(50.0000f)); a1.add(Float.valueOf(2.5000f)); b1.add(Float.valueOf(0.0000f)); ell2.add(Float.valueOf(50.0000f)); a2.add(Float.valueOf(3.2592f)); b2.add(Float.valueOf(0.3350f)); deltaE.add(Double.valueOf(1.0000));
        ell1.add(Float.valueOf(60.2574f)); a1.add(Float.valueOf(-34.0099f)); b1.add(Float.valueOf(36.2677f)); ell2.add(Float.valueOf(60.4626f)); a2.add(Float.valueOf(-34.1751f)); b2.add(Float.valueOf(39.4387f)); deltaE.add(Double.valueOf(1.2644));
        ell1.add(Float.valueOf(63.0109f)); a1.add(Float.valueOf(-31.0961f)); b1.add(Float.valueOf(-5.8663f)); ell2.add(Float.valueOf(62.8187f)); a2.add(Float.valueOf(-29.7946f)); b2.add(Float.valueOf(-4.0864f)); deltaE.add(Double.valueOf(1.2630));
        ell1.add(Float.valueOf(61.2901f)); a1.add(Float.valueOf(3.7196f)); b1.add(Float.valueOf(-5.3901f)); ell2.add(Float.valueOf(61.4292f)); a2.add(Float.valueOf(2.2480f)); b2.add(Float.valueOf(-4.9620f)); deltaE.add(Double.valueOf(1.8731));
        ell1.add(Float.valueOf(35.0831f)); a1.add(Float.valueOf(-44.1164f)); b1.add(Float.valueOf(3.7933f)); ell2.add(Float.valueOf(35.0232f)); a2.add(Float.valueOf(-40.0716f)); b2.add(Float.valueOf(1.5901f)); deltaE.add(Double.valueOf(1.8645));
        ell1.add(Float.valueOf(22.7233f)); a1.add(Float.valueOf(20.0904f)); b1.add(Float.valueOf(-46.6940f)); ell2.add(Float.valueOf(23.0331f)); a2.add(Float.valueOf(14.9730f)); b2.add(Float.valueOf(-42.5619f)); deltaE.add(Double.valueOf(2.0373));
        ell1.add(Float.valueOf(36.4612f)); a1.add(Float.valueOf(47.8580f)); b1.add(Float.valueOf(18.3852f)); ell2.add(Float.valueOf(36.2715f)); a2.add(Float.valueOf(50.5065f)); b2.add(Float.valueOf(21.2231f)); deltaE.add(Double.valueOf(1.4146));
        ell1.add(Float.valueOf(90.8027f)); a1.add(Float.valueOf(-2.0831f)); b1.add(Float.valueOf(1.4410f)); ell2.add(Float.valueOf(91.1528f)); a2.add(Float.valueOf(-1.6435f)); b2.add(Float.valueOf(0.0447f)); deltaE.add(Double.valueOf(1.4441));
        ell1.add(Float.valueOf(90.9257f)); a1.add(Float.valueOf(-0.5406f)); b1.add(Float.valueOf(-0.9208f)); ell2.add(Float.valueOf(88.6381f)); a2.add(Float.valueOf(-0.8985f)); b2.add(Float.valueOf(-0.7239f)); deltaE.add(Double.valueOf(1.5381));
        ell1.add(Float.valueOf(6.7747f)); a1.add(Float.valueOf(-0.2908f)); b1.add(Float.valueOf(-2.4247f)); ell2.add(Float.valueOf(5.8714f)); a2.add(Float.valueOf(-0.0985f)); b2.add(Float.valueOf(-2.2286f)); deltaE.add(Double.valueOf(0.6377));
        ell1.add(Float.valueOf(2.0776f)); a1.add(Float.valueOf(0.0795f)); b1.add(Float.valueOf(-1.1350f)); ell2.add(Float.valueOf(0.9033f)); a2.add(Float.valueOf(-0.0636f)); b2.add(Float.valueOf(-0.5514f)); deltaE.add(Double.valueOf(0.9082));
        
    }
    
}
