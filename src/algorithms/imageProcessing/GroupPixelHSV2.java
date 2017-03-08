package algorithms.imageProcessing;

import algorithms.util.PairInt;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.set.TIntSet;
import java.awt.Color;
import java.util.Set;

/**
 *
 * @author nichole
 */
public class GroupPixelHSV2 {
    
    private double sumR;
    private double sumG;
    private double sumB;
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
        sumR = 0;
        sumG = 0;
        sumB = 0;
        
        for (PairInt p : points) {
            
            int pixIdx = img.getInternalIndex(p);
            
            int r = img.getR(pixIdx);
            int g = img.getG(pixIdx);
            int b = img.getB(pixIdx);
            
            Color.RGBtoHSB(r, g, b, hsv);
            
            sumH += hsv[0];
            sumS += hsv[1];
            sumV += hsv[2];
            
            sumR += r;
            sumG += g;
            sumB += b;
        }
                
        nPoints = points.size();
    }
    
    /**
     * resets the internal sums to 0 and sets them with calculated HSB
     * of points.
     * 
     * @param pIdxs pixel indexes
     * @param img 
     */
    public void calculateColors(final TIntSet pIdxs, ImageExt img) {
        
        float[] hsv = new float[3];
        
        nPoints = 0;
        sumH = 0;
        sumS = 0;
        sumV = 0;
        sumR = 0;
        sumG = 0;
        sumB = 0;
        
        TIntIterator iter = pIdxs.iterator();
        
        while (iter.hasNext()) {
            
            int pixIdx = iter.next();
            
            int r = img.getR(pixIdx);
            int g = img.getG(pixIdx);
            int b = img.getB(pixIdx);
            
            Color.RGBtoHSB(r, g, b, hsv);
            
            sumH += hsv[0];
            sumS += hsv[1];
            sumV += hsv[2];
        
            sumR += r;
            sumG += g;
            sumB += b;
        }
                
        nPoints = pIdxs.size();
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
            
            int r = img.getR(pixIdx);
            int g = img.getG(pixIdx);
            int b = img.getB(pixIdx);
            
            Color.RGBtoHSB(r, g, b, hsv);
            
            sumR += r;
            sumG += g;
            sumB += b;
                
            sumH += hsv[0];
            sumS += hsv[1];
            sumV += hsv[2];
        }
                
        nPoints += points.size();
    }
    
    /**
     * adds the calculated HSB of points to the current instance sums.
     * 
     * @param pIdxs pixel indexes
     * @param img 
     */
    public void add(final TIntSet pIdxs, ImageExt img) {
        
        float[] hsv = new float[3];
        
        TIntIterator iter = pIdxs.iterator();
        
        while (iter.hasNext()) {
            
            int pixIdx = iter.next();
            
            int r = img.getR(pixIdx);
            int g = img.getG(pixIdx);
            int b = img.getB(pixIdx);
            
            Color.RGBtoHSB(r, g, b, hsv);
            
            sumR += r;
            sumG += g;
            sumB += b;
            
            sumH += hsv[0];
            sumS += hsv[1];
            sumV += hsv[2];
        }
                
        nPoints += pIdxs.size();
    }
    
    /**
     * adds the calculated HSB of points to the current instance sums.
     * 
     * @param other 
     */
    public void add(final GroupPixelHSV2 other) {
        
        sumR += other.sumR;
        sumG += other.sumG;
        sumB += other.sumB;

        sumH += other.sumH;
        sumS += other.sumS;
        sumV += other.sumV;
                
        nPoints += other.nPoints;
    }
    
    public void addPoint(final PairInt point, ImageExt img) {
                            
        int pixIdx = img.getInternalIndex(point);

        addPoint(pixIdx, img);
    }
    
    public void addPoint(final int pixIdx, ImageExt img) {
        
        float[] hsv = new float[3];
          
        int r = img.getR(pixIdx);
        int g = img.getG(pixIdx);
        int b = img.getB(pixIdx);

        Color.RGBtoHSB(r, g, b, hsv);
        
        sumR += r;
        sumG += g;
        sumB += b;

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
     * @return the avgR
     */
    public float getAvgR() {
        return (float)(sumR/(double)nPoints);
    }
    
    /**
     * @return the avgG
     */
    public float getAvgG() {
        return (float)(sumG/(double)nPoints);
    }

    /**
     * @return the avgB
     */
    public float getAvgB() {
        return (float)(sumB/(double)nPoints);
    }
    
    /**
     * @return the nPoints
     */
    public int getNPoints() {
        return nPoints;
    }

    public float calculateDifference(GroupPixelHSV2 hsv2) {

        float sumDiff = Math.abs(getAvgH() - hsv2.getAvgH()) +
            Math.abs(getAvgS() - hsv2.getAvgS()) + 
            Math.abs(getAvgV() - hsv2.getAvgV());
    
        sumDiff /= 3.f;
        
        return sumDiff;
    }
    
    public float[] calculateDifferences(GroupPixelHSV2 hsv2) {

        float[] a = new float[]{
            Math.abs(getAvgH() - hsv2.getAvgH()),
            Math.abs(getAvgS() - hsv2.getAvgS()), 
            Math.abs(getAvgV() - hsv2.getAvgV())};
            
        return a;
    }
    
    /**
     * values are normalized from 0 to 255 to 0.f to 1.f
     * @param hsv2
     * @return 
     */
    public float calculateRGBDifference(GroupPixelHSV2 hsv2) {

        float sumDiff = Math.abs(getAvgR() - hsv2.getAvgR()) +
            Math.abs(getAvgG() - hsv2.getAvgG()) + 
            Math.abs(getAvgB() - hsv2.getAvgB());
    
        sumDiff /= (255.f * 3.f);
        
        return sumDiff;
    }
    
    /**
     * values are normalized from 0 to 255 to 0.f to 1.f
     * @param hsv2
     * @return 
     */
    public float[] calculateRGBDifferences(GroupPixelHSV2 hsv2) {

        float[] a = new float[]{
            Math.abs(getAvgR() - hsv2.getAvgR()),
            Math.abs(getAvgG() - hsv2.getAvgG()), 
            Math.abs(getAvgB() - hsv2.getAvgB())};
        
        for (int i = 0; i < a.length; ++i) {
            a[i] /= 255.f;
        }
        
        return a;
    }
    
    public boolean isGrey(int limit) {
        
        int r = Math.round(getAvgR());
        int g = Math.round(getAvgG());
        int b = Math.round(getAvgB());
        
        // looking at whether color is grey
        int avgRGB = (r + g + b)/3;
        
        /*
        System.out.format("    -> (%d,%d,%d) %d,%d,%d\n",
            r, g, b,
            (Math.abs(r - avgRGB)),
            (Math.abs(g - avgRGB)),
            (Math.abs(b - avgRGB)));
        */
        
        if ((Math.abs(r - avgRGB) < limit) &&
            (Math.abs(g - avgRGB) < limit) &&
            (Math.abs(b - avgRGB) < limit)) {
            return true;
        }
        
        return false;
    }

    public boolean isGreen() {
        
        int r = Math.round(getAvgR());
        int g = Math.round(getAvgG());
        int b = Math.round(getAvgB());
    
        //System.out.println("r=" + r + " g=" + g + " b=" + b);
        
        int limit = 8;//3% of 256
        
        return ((g - r) >= limit) && ((g - b) >= limit);
    }
}
