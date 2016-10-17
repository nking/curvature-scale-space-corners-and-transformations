
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

    public GroupPixelHSV() {
    }

    public void calculateColors(final Set<PairInt> points, ImageExt colorImage) {
        
        GroupPixelRGB0 rgb = new GroupPixelRGB0();
       
        rgb.calculateColors(points, colorImage, 0, 0);
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
            
            int r = colorImage.getR(x, y);
            int g = colorImage.getG(x, y);
            int b = colorImage.getB(x, y);
            
            Color.RGBtoHSB(r, g, b, hsv);
            
            float diffH = hsv[0] - hsvAvg[0];
            float diffS = hsv[1] - hsvAvg[1];
            float diffV = hsv[2] - hsvAvg[2];
            
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
}
