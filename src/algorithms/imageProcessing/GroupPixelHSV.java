package algorithms.imageProcessing;

import algorithms.util.PairInt;
import java.awt.Color;
import java.util.Set;

/**
 *
 * @author nichole
 */
public class GroupPixelHSV {
    
    private float avgH;
    private float avgS;
    private float avgV;
    private int nPoints;
    
    private float stdDevH;
    private float stdDevS;
    private float stdDevV;
    
    private GroupPixelCIELAB rgbAndCIEClrs;

    public GroupPixelHSV() {
    }

    public void calculateColors(final Set<PairInt> points, ImageExt colorImage) {
        
        GroupPixelCIELAB rgb = new GroupPixelCIELAB(points, colorImage);
        rgbAndCIEClrs = rgb ;
        
        float[] hsvAvg = new float[3];
        Color.RGBtoHSB(Math.round(rgb.getAvgRed()), 
            Math.round(rgb.getAvgGreen()), 
            Math.round(rgb.getAvgBlue()), hsvAvg);
        
        this.nPoints = points.size();
        this.avgH = hsvAvg[0];
        this.avgS = hsvAvg[1];
        this.avgV = hsvAvg[2];
        
        float[] hsv = new float[3];
        double sumH = 0;
        double sumS = 0;
        double sumV = 0;
        for (PairInt p : points) {
            
            int x = p.getX();
            int y = p.getY();
            
            float diffH = colorImage.getHue(x, y) - hsvAvg[0];
            float diffS = colorImage.getSaturation(x, y) - hsvAvg[1];
            float diffV = colorImage.getBrightness(x, y) - hsvAvg[2];
            
            sumH += (diffH * diffH);
            sumS += (diffS * diffS);
            sumV += (diffV * diffV);
        }
        
        this.stdDevH = (float)Math.sqrt(sumH/(nPoints - 1.));
        this.stdDevS = (float)Math.sqrt(sumS/(nPoints - 1.));
        this.stdDevV = (float)Math.sqrt(sumV/(nPoints - 1.));
    }

    /**
     * @return the avgH
     */
    public float getAvgH() {
        return avgH;
    }

    /**
     * @return the avgS
     */
    public float getAvgS() {
        return avgS;
    }

    /**
     * @return the avgV
     */
    public float getAvgV() {
        return avgV;
    }

    /**
     * @return the nPoints
     */
    public int getNPoints() {
        return nPoints;
    }

    /**
     * @return the stdDevH
     */
    public float getStdDevH() {
        return stdDevH;
    }

    /**
     * @return the stdDevS
     */
    public float getStdDevS() {
        return stdDevS;
    }

    /**
     * @return the stdDevV
     */
    public float getStdDevV() {
        return stdDevV;
    }
   
     /**
     * @return the avgL
     */
    public float getAvgL() {
        return rgbAndCIEClrs.getAvgL();
    }

    /**
     * @return the stdDevL
     */
    public float getStdDevL() {
        return rgbAndCIEClrs.getStdDevL();
    }

    /**
     * @return the avgA
     */
    public float getAvgA() {
        return rgbAndCIEClrs.getAvgA();
    }

    /**
     * @return the stdDevA
     */
    public float getStdDevA() {
        return rgbAndCIEClrs.getStdDevA();
    }

    /**
     * @return the avgB
     */
    public float getAvgB() {
        return rgbAndCIEClrs.getAvgB();
    }

    /**
     * @return the stdDevB
     */
    public float getStdDevB() {
        return rgbAndCIEClrs.getStdDevB();
    }
    
}
