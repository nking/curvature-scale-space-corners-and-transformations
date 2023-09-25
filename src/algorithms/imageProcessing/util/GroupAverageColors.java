package algorithms.imageProcessing.util;

import algorithms.imageProcessing.CIEChromaticity;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.compGeometry.MiscellaneousCurveHelper;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import java.util.Collection;

/**
 *
 * @author nichole
 */
public class GroupAverageColors {
    
    private int xCen;
    private int yCen;
    
    private int rAvg;
    private int gAvg;
    private int bAvg;
    
    private float avgCIEL;
    private float avgCIEA;
    private float avgCIEB;
    
    public GroupAverageColors(Image img,
        PairIntArray a) {
        
        CIEChromaticity cieC = new CIEChromaticity();
            
        MiscellaneousCurveHelper curveHelper =
            new MiscellaneousCurveHelper();
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        double[] xyCen = curveHelper.calculateXYCentroids(a);
                
        this.xCen = (int)Math.round(xyCen[0]);
        
        this.yCen = (int)Math.round(xyCen[1]);
                
        int[] avgRGB = imageProcessor.getAverageRGB(img, a);
        
        this.rAvg = avgRGB[0];
        this.gAvg = avgRGB[1];
        this.bAvg = avgRGB[2];
        
        float[] lab = cieC.rgbToCIELAB(
            avgRGB[0], avgRGB[1], avgRGB[2]);
                
        this.avgCIEL = lab[0];
        this.avgCIEA = lab[1];
        this.avgCIEB = lab[2];
    }
    
    public GroupAverageColors(GreyscaleImage rImg, GreyscaleImage gImg,
        GreyscaleImage bImg, Collection<PairInt> a) {
        
        CIEChromaticity cieC = new CIEChromaticity();
            
        MiscellaneousCurveHelper curveHelper =
            new MiscellaneousCurveHelper();
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        double[] xyCen = curveHelper.calculateXYCentroids(a);
                
        this.xCen = (int)Math.round(xyCen[0]);
        
        this.yCen = (int)Math.round(xyCen[1]);
                
        int[] avgRGB = imageProcessor.getAverageRGB(rImg, gImg, bImg, a);
        
        this.rAvg = avgRGB[0];
        this.gAvg = avgRGB[1];
        this.bAvg = avgRGB[2];
        
        float[] lab = cieC.rgbToCIELAB(
            avgRGB[0], avgRGB[1], avgRGB[2]);
                
        this.avgCIEL = lab[0];
        this.avgCIEA = lab[1];
        this.avgCIEB = lab[2];
    }
    
    public GroupAverageColors(Image img, Collection<PairInt> a) {
        
        CIEChromaticity cieC = new CIEChromaticity();
            
        MiscellaneousCurveHelper curveHelper =
            new MiscellaneousCurveHelper();
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        double[] xyCen = curveHelper.calculateXYCentroids(a);
                
        this.xCen = (int)Math.round(xyCen[0]);
        
        this.yCen = (int)Math.round(xyCen[1]);
                
        int[] avgRGB = imageProcessor.getAverageRGB(img, a);
        
        this.rAvg = avgRGB[0];
        this.gAvg = avgRGB[1];
        this.bAvg = avgRGB[2];
        
        float[] lab = cieC.rgbToCIELAB(
            avgRGB[0], avgRGB[1], avgRGB[2]);
                
        this.avgCIEL = lab[0];
        this.avgCIEA = lab[1];
        this.avgCIEB = lab[2];
    }
    
    public double calculateDeltaE2000(GroupAverageColors other) {
        
        CIEChromaticity cieC = new CIEChromaticity();
        
        double delta = 
            cieC.calcDeltaECIE2000(avgCIEL, avgCIEA, avgCIEB,
                other.avgCIEL, other.avgCIEA, other.avgCIEB);
        
        return delta;
    }

    /**
     * @return the xCen
     */
    public int getXCen() {
        return xCen;
    }

    /**
     * @param xCen the xCen to set
     */
    public void setXCen(int xCen) {
        this.xCen = xCen;
    }

    /**
     * @return the yCen
     */
    public int getYCen() {
        return yCen;
    }

    /**
     * @param yCen the yCen to set
     */
    public void setYCen(int yCen) {
        this.yCen = yCen;
    }

    /**
     * @return the rAvg
     */
    public int getR() {
        return rAvg;
    }

    /**
     * @return the gAvg
     */
    public int getG() {
        return gAvg;
    }

    /**
     * @return the bAvg
     */
    public int getB() {
        return bAvg;
    }

    /**
     * @return the avgCIEL
     */
    public float getCIEL() {
        return avgCIEL;
    }

    /**
     * @return the avgCIEA
     */
    public float getCIEA() {
        return avgCIEA;
    }

    /**
     * @return the avgCIEB
     */
    public float getCIEB() {
        return avgCIEB;
    }
}
