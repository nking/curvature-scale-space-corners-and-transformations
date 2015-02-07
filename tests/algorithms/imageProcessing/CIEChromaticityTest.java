package algorithms.imageProcessing;

import algorithms.compGeometry.PointInPolygon;
import algorithms.util.ArrayPair;
import algorithms.util.ResourceFinder;
import java.awt.Color;
import java.util.Arrays;
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
        
        ArrayPair t = yellowOrangeBounds;
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
        log.info(String.format(
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
        log.info(String.format(
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
        log.info(String.format(
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
        log.info(String.format(
            "(%d,%d)  rgb=(%d,%d,%d) pixCIEXY=%f,%f  hsb=(%f,%f,%f)", 
            col, row, pixR, pixG, pixB,
            pixCIEXY[0], pixCIEXY[1],
            hsb[0], hsb[1], hsb[2]));
        
        
        ArrayPair redBounds = cieC.getRedPolynomial();
        ArrayPair orangeBounds = cieC.getOrangePolynomial();
        ArrayPair yellowBounds = cieC.getYellowPolynomial();
        ArrayPair yellowOrangeBounds = cieC.getYellowThroughOrangePolynomial();
        ArrayPair greenishYellowOrangeBounds = cieC.getYellowishGreenThroughOrangePolynomial();
        
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

                        log.info(String.format(
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

                        log.info(String.format(
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

                        log.info(String.format(
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
            
        ImageDisplayer.displayImage("points", img1);
        ImageIOHelper.writeOutputImage(filePath2, img1);
        int z = 1;
    }
    
    public void test3() throws Exception {
    
        String filePath1 = ResourceFinder.findFileInTestResources(
            "sky_with_rainbow.jpg");
        Image img1 = ImageIOHelper.readImage(filePath1);        
        
        BrightSkyRainbowColors colors = new BrightSkyRainbowColors();
        
        for (int col = 0; col < img1.getWidth(); col++) {
            for (int row = 0; row < img1.getHeight(); row++) {
                
                int r = img1.getR(col, row);
                int g = img1.getG(col, row);
                int b = img1.getB(col, row);
                
                if (colors.isInRedThroughPurplishRed(r, g, b)) {
                    
                    img1.setRGB(col, row, 0, 0, 255);

                } else if (colors.isInOrangeRed(r, g, b)) {
                    
                    img1.setRGB(col, row, 0, 0, 255);
                
                } else if (colors.isInGreenishYellowOrange(r, g, b)) {
                 
                    img1.setRGB(col, row, 0, 0, 255);

                }
            }
        }
        
        String dirPath = ResourceFinder.findDirectory("bin");
        String filePath2 = dirPath + "/tmp3.png";
            
        ImageDisplayer.displayImage("points", img1);
        ImageIOHelper.writeOutputImage(filePath2, img1);
        int z = 1;
    }
    
    public void test4() throws Exception {
    
        String filePath1 = ResourceFinder.findFileInTestResources(
            "sky_with_rainbow2.jpg");
        Image img1 = ImageIOHelper.readImage(filePath1);        
        
        DarkSkyRainbowColors colors = new DarkSkyRainbowColors();
        
        for (int col = 0; col < img1.getWidth(); col++) {
            for (int row = 0; row < img1.getHeight(); row++) {
                
                int r = img1.getR(col, row);
                int g = img1.getG(col, row);
                int b = img1.getB(col, row);
                
                if (colors.isInRedThroughPurplishRed(r, g, b)) {
                    
                    img1.setRGB(col, row, 0, 0, 255);

                } else if (colors.isInOrangeRed(r, g, b)) {
                    
                    img1.setRGB(col, row, 0, 0, 255);
                
                } else if (colors.isInGreenishYellowOrange(r, g, b)) {
                 
                    img1.setRGB(col, row, 0, 0, 255);

                }
            }
        }
        
        String dirPath = ResourceFinder.findDirectory("bin");
        String filePath2 = dirPath + "/tmp4.png";
            
        ImageDisplayer.displayImage("points", img1);
        ImageIOHelper.writeOutputImage(filePath2, img1);
        int z = 1;
    }
}
