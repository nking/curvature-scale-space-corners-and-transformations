package algorithms.imageProcessing.util;

import algorithms.imageProcessing.CIEChromaticity;
import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.MiscellaneousCurveHelper;
import algorithms.util.PairIntArray;

/**
 *
 * @author nichole
 */
public class GroupAverageColors {
    
    private int xCen;
    private int yCen;
    
    protected int rAvg;
    protected int gAvg;
    protected int bAvg;
    
    protected float avgCIEL;
    protected float avgCIEA;
    protected float avgCIEB;
    
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
}
