package algorithms.imageProcessing;

import algorithms.compGeometry.PointInPolygon;
import algorithms.util.ArrayPair;
import java.awt.Color;

/**
 *
 * @author nichole
 */
public abstract class AbstractSkyRainbowColors {
    
    protected CIEChromaticity cieC = new CIEChromaticity();
    
    protected ArrayPair orangeBounds = cieC.getOrangePolynomial();
    
    protected ArrayPair greenishYellowOrangeBounds = 
        cieC.getYellowishGreenThroughOrangePolynomial();
    
    protected ArrayPair redPurplishBounds = 
        cieC.getRedThroughPurplishRedPolynomial();
    
    protected float orangeSLimit = 0.2f;
    
    protected float greenishYellowOrangeSLimit = 0.1f;
        
    protected float redPurpleSLimit = 0.05f;
    
    protected float totRGBLimit = 0.2f * 768.f;
            
    public boolean isInRedThroughPurplishRed(int r, int g, int b) {
        
        float[] pixCIEXY = cieC.rgbToXYChromaticity(r, g, b);

        int totRGB = r + g + b;

        float[] hsb = new float[3];
        Color.RGBtoHSB(r, b, b, hsb);
        
        PointInPolygon pInPoly = new PointInPolygon();
         
        if ((r > 15) && (g > 15) && (b > 15)) {
                    
            if (
                (hsb[1] > redPurpleSLimit)
                && (hsb[2] > 0.25)
                && (totRGB > totRGBLimit)
                && (pInPoly.isInSimpleCurve(pixCIEXY[0], pixCIEXY[1],
                redPurplishBounds.getX(), redPurplishBounds.getY(), 
                redPurplishBounds.getX().length))) {
                
                return true;
            }
        }
        
        return false;
    }
    
    public boolean isInOrangeRed(int r, int g, int b) {
        
        float[] pixCIEXY = cieC.rgbToXYChromaticity(r, g, b);

        int totRGB = r + g + b;

        float[] hsb = new float[3];
        Color.RGBtoHSB(r, b, b, hsb);
        
        PointInPolygon pInPoly = new PointInPolygon();
         
        if ((r > 15) && (g > 15) && (b > 15)) {
                    
            if (
                (hsb[1] > orangeSLimit)
                && (hsb[2] > 0.25)
                && (totRGB > totRGBLimit)
                && (pInPoly.isInSimpleCurve(pixCIEXY[0], pixCIEXY[1],
                orangeBounds.getX(), orangeBounds.getY(), 
                orangeBounds.getX().length))) {
                
                return true;
            }
        }
        
        return false;
    }
    
    public boolean isInGreenishYellowOrange(int r, int g, int b) {
        
        float[] pixCIEXY = cieC.rgbToXYChromaticity(r, g, b);

        int totRGB = r + g + b;

        float[] hsb = new float[3];
        Color.RGBtoHSB(r, b, b, hsb);
        
        PointInPolygon pInPoly = new PointInPolygon();
         
        if ((r > 15) && (g > 15) && (b > 15)) {
                    
            if (
                (hsb[1] > greenishYellowOrangeSLimit)
                && (hsb[2] > 0.25)
                && (totRGB > totRGBLimit)
                && (pInPoly.isInSimpleCurve(pixCIEXY[0], pixCIEXY[1],
                greenishYellowOrangeBounds.getX(), 
                greenishYellowOrangeBounds.getY(), 
                greenishYellowOrangeBounds.getX().length))) {
                
                return true;
            }
        }
        
        return false;
    }
}
