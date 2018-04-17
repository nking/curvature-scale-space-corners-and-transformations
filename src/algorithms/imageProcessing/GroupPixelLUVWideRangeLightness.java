package algorithms.imageProcessing;

import algorithms.util.PairInt;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.set.TIntSet;
import java.util.Set;

/**
 *
 * @author nichole
 */
public class GroupPixelLUVWideRangeLightness {
    
    private float[] avgLUV;
    private long avgR;
    private long avgG;
    private long avgB;
    private int nPoints;
    
    public GroupPixelLUVWideRangeLightness() {
    }

    /**
     * resets the internal sums to 0 and sets them with calculated HSB
     * of points.
     * 
     * @param points
     * @param img 
     */
    public void calculateColors(final Set<PairInt> points, ImageExt img) {
                
        nPoints = 0;
        avgR = 0;
        avgG = 0;
        avgB = 0;
        
        for (PairInt p : points) {
            
            int pixIdx = img.getInternalIndex(p);
            
            avgR += img.getR(pixIdx);
            avgG += img.getG(pixIdx);
            avgB += img.getB(pixIdx);
        }
                
        nPoints = points.size();
        
        avgR /= nPoints;
        avgG /= nPoints;
        avgB /= nPoints;
        
        CIEChromaticity cieC = new CIEChromaticity();
        avgLUV = cieC.rgbToCIELUV_WideRangeLightness((int)getAvgR(), (int)getAvgG(), (int)getAvgB());
    }
    
    /**
     * resets the internal sums to 0 and sets them with calculated HSB
     * of points.
     * 
     * @param pIdxs pixel indexes
     * @param img 
     */
    public void calculateColors(final TIntSet pIdxs, ImageExt img) {
        
        nPoints = 0;
        avgR = 0;
        avgG = 0;
        avgB = 0;
        
        TIntIterator iter = pIdxs.iterator();
        
        while (iter.hasNext()) {
            
            int pixIdx = iter.next();
            
            avgR += img.getR(pixIdx);
            avgG += img.getG(pixIdx);
            avgB += img.getB(pixIdx);
        }
                
        nPoints = pIdxs.size();
        
        avgR /= nPoints;
        avgG /= nPoints;
        avgB /= nPoints;
        
        CIEChromaticity cieC = new CIEChromaticity();
        avgLUV = cieC.rgbToCIELUV_WideRangeLightness((int)getAvgR(), (int)getAvgG(), (int)getAvgB());
    }
    
    /**
     * adds the calculated LUV of points to the current instance sums.
     * 
     * @param points
     * @param img 
     */
    public void add(final Set<PairInt> points, ImageExt img) {
        
        avgR *= nPoints;
        avgG *= nPoints;
        avgB *= nPoints;
        
        for (PairInt p : points) {
            
            int pixIdx = img.getInternalIndex(p);
            
            avgR += img.getR(pixIdx);
            avgG += img.getG(pixIdx);
            avgB += img.getB(pixIdx);
        }
                
        nPoints += points.size();
        
        avgR /= nPoints;
        avgG /= nPoints;
        avgB /= nPoints;
        
        CIEChromaticity cieC = new CIEChromaticity();
        avgLUV = cieC.rgbToCIELUV_WideRangeLightness((int)getAvgR(), (int)getAvgG(), (int)getAvgB());
    }
    
    /**
     * adds the calculated HSB of points to the current instance sums.
     * 
     * @param pIdxs pixel indexes
     * @param img 
     */
    public void add(final TIntSet pIdxs, ImageExt img) {
        
        avgR *= nPoints;
        avgG *= nPoints;
        avgB *= nPoints;
        
        TIntIterator iter = pIdxs.iterator();
        
        while (iter.hasNext()) {
            
            int pixIdx = iter.next();
            
            avgR += img.getR(pixIdx);
            avgG += img.getG(pixIdx);
            avgB += img.getB(pixIdx);
        }
                
        nPoints += pIdxs.size();
        
        avgR /= nPoints;
        avgG /= nPoints;
        avgB /= nPoints;
        
        CIEChromaticity cieC = new CIEChromaticity();
        avgLUV = cieC.rgbToCIELUV_WideRangeLightness((int)getAvgR(), (int)getAvgG(), (int)getAvgB());
    }
    
    /**
     * adds the calculated HSB of points to the current instance sums.
     * 
     * @param other 
     */
    public void add(final GroupPixelLUVWideRangeLightness other) {
        
        avgR *= nPoints;
        avgG *= nPoints;
        avgB *= nPoints;
        
        avgR += other.getAvgR();
        avgG += other.getAvgG();
        avgB += other.getAvgB();

        nPoints += other.nPoints;
        
        avgR /= nPoints;
        avgG /= nPoints;
        avgB /= nPoints;
        
        CIEChromaticity cieC = new CIEChromaticity();
        avgLUV = cieC.rgbToCIELUV_WideRangeLightness((int)getAvgR(), (int)getAvgG(), (int)getAvgB());
    }
    
    public void addPoint(final PairInt point, ImageExt img) {
                            
        int pixIdx = img.getInternalIndex(point);

        addPoint(pixIdx, img);
    }
    
    public void addPoint(final int pixIdx, ImageExt img) {
        
        avgR *= nPoints;
        avgG *= nPoints;
        avgB *= nPoints;
        
        avgR += img.getR(pixIdx);
        avgG += img.getG(pixIdx);
        avgB += img.getB(pixIdx);
                
        nPoints++;
        
        avgR /= nPoints;
        avgG /= nPoints;
        avgB /= nPoints;
        
        CIEChromaticity cieC = new CIEChromaticity();
        avgLUV = cieC.rgbToCIELUV_WideRangeLightness((int)getAvgR(), (int)getAvgG(), (int)getAvgB());
        
    }

    /**
     * @return the nPoints
     */
    public int getNPoints() {
        return nPoints;
    }
    
    public float[] getAvgLUV() {
        return avgLUV;
    }

    public float calculateDifference(GroupPixelLUVWideRangeLightness luv2) {

        CIEChromaticity cieC = new CIEChromaticity();
        
        float diff = cieC.calcNormalizedDifferenceLUV_WideRangeLightness(avgLUV, 
            luv2.getAvgLUV());
        
        return diff;
    }
   
    /**
     * values are normalized from 0 to 255 to 0.f to 1.f
     * @param luv2
     * @return 
     */
    public float calculateRGBDifference(GroupPixelLUVWideRangeLightness luv2) {

        float sumDiff = Math.abs(getAvgR() - luv2.getAvgR()) +
            Math.abs(getAvgG() - luv2.getAvgG()) + 
            Math.abs(getAvgB() - luv2.getAvgB());
    
        sumDiff /= (255.f * 3.f);
        
        return sumDiff;
    }
    
    /**
     * values are normalized from 0 to 255 to 0.f to 1.f
     * @param luv2
     * @return 
     */
    public float[] calculateRGBDifferences(GroupPixelLUVWideRangeLightness luv2) {

        float[] a = new float[]{
            Math.abs(getAvgR() - luv2.getAvgR()),
            Math.abs(getAvgG() - luv2.getAvgG()), 
            Math.abs(getAvgB() - luv2.getAvgB())};
        
        for (int i = 0; i < a.length; ++i) {
            a[i] /= 255.f;
        }
        
        return a;
    }

    /**
     * @return the avgR
     */
    public long getAvgR() {
        return avgR;
    }

    /**
     * @return the avgG
     */
    public long getAvgG() {
        return avgG;
    }

    /**
     * @return the avgB
     */
    public long getAvgB() {
        return avgB;
    }
    
}
