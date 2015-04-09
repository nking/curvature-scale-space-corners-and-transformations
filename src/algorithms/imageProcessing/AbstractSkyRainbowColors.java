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

        float[] hsb = new float[3];
        Color.RGBtoHSB(r, g, b, hsb);
        
        return isInRedThroughPurplishRed(r, g, b, pixCIEXY[0], pixCIEXY[1], 
            hsb[1], hsb[2]);
    }
    
    public boolean isInRedThroughPurplishRed(ImageExt clrImg, int pixelIndex) {
        
        float cieX = clrImg.getCIEX(pixelIndex);
        float cieY = clrImg.getCIEY(pixelIndex);
        
        int r = clrImg.getR(pixelIndex);
        int g = clrImg.getG(pixelIndex);
        int b = clrImg.getB(pixelIndex);
        
        float saturation = clrImg.getSaturation(pixelIndex);
        float brightness = clrImg.getBrightness(pixelIndex);

        return isInRedThroughPurplishRed(r, g, b, cieX, cieY, saturation, 
            brightness);
    }
    
    public boolean isInRedThroughPurplishRed(int r, int g, int b, 
        float cieX, float cieY, float saturation, float brightness) {
        
        int totRGB = r + g + b;
        
        PointInPolygon pInPoly = new PointInPolygon();
         
        if ((r > 15) && (g > 15) && (b > 15)) {
                    
            if (
                (saturation > redPurpleSLimit)
                && (brightness > 0.25)
                && (totRGB > totRGBLimit)
                && (pInPoly.isInSimpleCurve(cieX, cieY,
                redPurplishBounds.getX(), redPurplishBounds.getY(), 
                redPurplishBounds.getX().length))) {
                
                return true;
            }
        }
        
        return false;
    }
    
    public boolean isInOrangeRed(ImageExt clrImg, int pixelIndex) {
        
        float cieX = clrImg.getCIEX(pixelIndex);
        float cieY = clrImg.getCIEY(pixelIndex);
        
        int r = clrImg.getR(pixelIndex);
        int g = clrImg.getG(pixelIndex);
        int b = clrImg.getB(pixelIndex);
        
        float saturation = clrImg.getSaturation(pixelIndex);
        float brightness = clrImg.getBrightness(pixelIndex);
        
        return isInOrangeRed(r, g, b, cieX, cieY, saturation, brightness);
    }
    
    public boolean isInOrangeRed(int r, int g, int b) {
        
        float[] pixCIEXY = cieC.rgbToXYChromaticity(r, g, b);

        float[] hsb = new float[3];
        Color.RGBtoHSB(r, g, b, hsb);
        
        return isInOrangeRed(r, g, b, pixCIEXY[0], pixCIEXY[1], hsb[1], hsb[2]);
    }
    
    public boolean isInOrangeRed(int r, int g, int b, 
        float cieX, float cieY, float saturation, float brightness) {
        
        int totRGB = r + g + b;

        PointInPolygon pInPoly = new PointInPolygon();
         
        if ((r > 15) && (g > 15) && (b > 15)) {
                    
            if (
                (saturation > orangeSLimit)
                && (brightness > 0.25)
                && (totRGB > totRGBLimit)
                && (pInPoly.isInSimpleCurve(cieX, cieY,
                orangeBounds.getX(), orangeBounds.getY(), 
                orangeBounds.getX().length))) {
                
                return true;
            }
        }
        
        return false;
    }
    
    public boolean isInGreenishYellowOrange(ImageExt clrImg, int pixelIndex) {
        
        float cieX = clrImg.getCIEX(pixelIndex);
        float cieY = clrImg.getCIEY(pixelIndex);
        
        int r = clrImg.getR(pixelIndex);
        int g = clrImg.getG(pixelIndex);
        int b = clrImg.getB(pixelIndex);
        
        float saturation = clrImg.getSaturation(pixelIndex);
        float brightness = clrImg.getBrightness(pixelIndex);
        
        return isInGreenishYellowOrange(r, g, b, cieX, cieY, saturation, 
            brightness);
    }
    
    public boolean isInGreenishYellowOrange(int r, int g, int b) {
        
        float[] pixCIEXY = cieC.rgbToXYChromaticity(r, g, b);

        float[] hsb = new float[3];
        Color.RGBtoHSB(r, g, b, hsb);
        
        return isInGreenishYellowOrange(r, g, b, pixCIEXY[0], pixCIEXY[1], 
            hsb[1], hsb[2]);
    }
    
    public boolean isInGreenishYellowOrange(int r, int g, int b, 
        float cieX, float cieY, float saturation, float brightness) {
        
        int totRGB = r + g + b;

        PointInPolygon pInPoly = new PointInPolygon();
                  
        if ((r > 15) && (g > 15) && (b > 15)) {
                    
            if (
                (saturation > greenishYellowOrangeSLimit)
                && (brightness > 0.25)
                && (totRGB > totRGBLimit)
                && (pInPoly.isInSimpleCurve(cieX, cieY,
                greenishYellowOrangeBounds.getX(), 
                greenishYellowOrangeBounds.getY(), 
                greenishYellowOrangeBounds.getX().length))) {
                
                return true;
            }
        }
        
        return false;
    }
}
