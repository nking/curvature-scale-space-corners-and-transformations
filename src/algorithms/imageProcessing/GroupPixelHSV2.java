package algorithms.imageProcessing;

import algorithms.util.PairInt;
import java.awt.Color;
import java.util.Set;

/**
 *
 * @author nichole
 */
public class GroupPixelHSV2 {
    
    private double sumH;
    private double sumS;
    private double sumV;
    private int nPoints;
    
    public GroupPixelHSV2() {
    }

    /**
     * resets the internal sums to 0 and sets them with calculated HSB
     * of points.
     * 
     * @param points
     * @param img 
     */
    public void calculateColors(final Set<PairInt> points, ImageExt img) {
        
        float[] hsv = new float[3];
        
        nPoints = 0;
        sumH = 0;
        sumS = 0;
        sumV = 0;
        
        for (PairInt p : points) {
            
            int pixIdx = img.getInternalIndex(p);
            
            Color.RGBtoHSB(img.getR(pixIdx), img.getG(pixIdx), 
                img.getB(pixIdx), hsv);
            
            sumH += hsv[0];
            sumS += hsv[1];
            sumV += hsv[2];
        }
                
        nPoints = points.size();
    }
    
    /**
     * adds the calculated HSB of points to the current instance sums.
     * 
     * @param points
     * @param img 
     */
    public void add(final Set<PairInt> points, ImageExt img) {
        
        float[] hsv = new float[3];
        
        for (PairInt p : points) {
            
            int pixIdx = img.getInternalIndex(p);
            
            Color.RGBtoHSB(img.getR(pixIdx), img.getG(pixIdx), 
                img.getB(pixIdx), hsv);
            
            sumH += hsv[0];
            sumS += hsv[1];
            sumV += hsv[2];
        }
                
        nPoints += points.size();
    }
    
    public void addPoint(final PairInt point, ImageExt img) {
        
        float[] hsv = new float[3];
                    
        int pixIdx = img.getInternalIndex(point);

        Color.RGBtoHSB(img.getR(pixIdx), img.getG(pixIdx), 
            img.getB(pixIdx), hsv);

        sumH += hsv[0];
        sumS += hsv[1];
        sumV += hsv[2];
        nPoints++;
    }

    /**
     * @return the avgH
     */
    public float getAvgH() {
        return (float)(sumH/(double)nPoints);
    }

    /**
     * @return the avgS
     */
    public float getAvgS() {
        return (float)(sumS/(double)nPoints);
    }

    /**
     * @return the avgV
     */
    public float getAvgV() {
        return (float)(sumV/(double)nPoints);
    }

    /**
     * @return the nPoints
     */
    public int getNPoints() {
        return nPoints;
    }
    
}
