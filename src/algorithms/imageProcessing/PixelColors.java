package algorithms.imageProcessing;

import algorithms.imageProcessing.util.MatrixUtil;
import java.awt.Color;

/**
 *
 * @author nichole
 */
public class PixelColors {
    
    private final int red;
    private final int green;
    private final int blue;
    
    private float[] hsb = null;
    
    private double[] yuv = null;
    
    private float[] cieXY = null;
    
    public PixelColors(int r, int g, int b) {
        this.red = r;
        this.green = g;
        this.blue = b;
    }
   
    public int getRed() {
        return red;
    }
    
    public int getGreen() {
        return green;
    }
    
    public int getBlue() {
        return blue;
    }
    
    public float getHue() {
        
        if (hsb == null) {
            hsb = new float[3];
            Color.RGBtoHSB(red, green, blue, hsb);
        }
        return hsb[0];
    }
    
    public float[] getCIEXY() {
        
        if (cieXY == null) {
            CIEChromaticity cieC = new CIEChromaticity();
            cieXY = cieC.rgbToXYChromaticity(red, green, blue);
        }
        
        return cieXY;
    }
    
    public float getSaturation() {
        
        if (hsb == null) {
            hsb = new float[3];
            Color.RGBtoHSB(red, green, blue, hsb);
        }
        return hsb[1];
    }
    
    public float getBrightness() {
        
        if (hsb == null) {
            hsb = new float[3];
            Color.RGBtoHSB(red, green, blue, hsb);
        }
        return hsb[2];
    }
    
    public double calculateContrastToOther(int otherRed, int otherGreen, int otherBlue) {
        
        double[][] m = new double[3][];
        m[0] = new double[]{0.256, 0.504, 0.098};
        m[1] = new double[]{-0.148, -0.291, 0.439};
        m[2] = new double[]{0.439, -0.368, -0.072};
        
        if (yuv == null) {
            yuv = MatrixUtil.multiply(m, new double[]{red, green, blue});
        }
        
        double[] yuvOther = MatrixUtil.multiply(m, new double[]{otherRed, otherGreen, otherBlue});
        
        double contrast = (yuv[0] - yuvOther[0])/yuvOther[0];
        
        return contrast;
    }
    
    /**
     * calculate the difference between red, green, and blue of other and this,
     * added in quadrature.  if the scale is 0-255 inclusive, the maximum
     * possible value returned will be 441.7 and the minimum will be 0.
     * @param otherRed
     * @param otherGreen
     * @param otherBlue
     * @return 
     */
    public double calculateColorDiffererenceToOther(int otherRed, int otherGreen, int otherBlue) {
        
        double rDiff = otherRed - red;
        double gDiff = otherGreen - green;
        double bDiff = otherBlue - blue;
        
        double dist = Math.sqrt(rDiff*rDiff + gDiff*gDiff + bDiff*bDiff);
        
        return dist;
    }    
}
