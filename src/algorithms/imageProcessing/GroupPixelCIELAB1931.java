package algorithms.imageProcessing;

import algorithms.util.PairInt;
import java.util.Set;

/**
 *
 * @author nichole
 */
public class GroupPixelCIELAB1931 extends GroupPixelRGB0 {
    
    private float avgL;
    private float stdDevL;
    private float avgA;
    private float stdDevA;
    private float avgB;
    private float stdDevB;
    
    public GroupPixelCIELAB1931(final Set<PairInt> points, ImageExt colorImage) {
               
        if (points == null) {
            throw new IllegalArgumentException("points cannot be null");
        }
        if (colorImage == null) {
            throw new IllegalArgumentException("colorImage cannot be null");
        }
        
        calculateColors(points, colorImage, 0, 0);
        
        CIEChromaticity cieC = new CIEChromaticity();
        float[] cieLABAvg = cieC.rgbToCIELAB1931(
            Math.round(getAvgRed()), 
            Math.round(getAvgGreen()), 
            Math.round(getAvgBlue()));
        
        this.nPoints = points.size();
        this.avgL = cieLABAvg[0];
        this.avgA = cieLABAvg[1];
        this.avgB = cieLABAvg[2];
        
        double sumL = 0;
        double sumA = 0;
        double sumB = 0;
        
        for (PairInt p : points) {
            
            int x = p.getX();
            int y = p.getY();
            
            float[] lab = cieC.rgbToCIELAB1931(
                colorImage.getR(x, y), colorImage.getG(x, y), 
                colorImage.getB(x, y));
            
            float diffL = lab[0] - cieLABAvg[0];
            float diffA = lab[1] - cieLABAvg[1];
            float diffB = lab[2] - cieLABAvg[2];
            
            sumL += (diffL * diffL);
            sumA += (diffA * diffA);
            sumB += (diffB * diffB);
        }
        
        this.stdDevL = (float)Math.sqrt(sumL/(nPoints - 1.));
        this.stdDevA = (float)Math.sqrt(sumA/(nPoints - 1.));
        this.stdDevB = (float)Math.sqrt(sumB/(nPoints - 1.));
    }

    /**
     * calculate the difference in L, A, B 1931 between this and other and
     * normalize the values to a sum of "1" using the range of values
     * possible from the use of a standard illuminant, D65.
     * 
     * @param other
     * @return 
     */
    public float calcNormalizedDifference(GroupPixelCIELAB1931 other) {
       
        /*
        * using the standard illuminant of daylight, D65,
        * the range of return values is
        * L*    0 to 104.5
        * a* -190 to 103
        * b* -113 to 99
        */
        
        float[] diffs = calcDifference(other);
        diffs[0] /= 104.5f;
        diffs[1] /= (103f + 190f);
        diffs[2] /= (99f + 113.f);
        
        return (diffs[0] + diffs[1] + diffs[2])/3.f;
    }
    
    /**
     * calculate the difference in L, U, V between this and other 
     * using the DeltaE2000 formula.
     * 
     * @param other
     * @return 
     */
    public float calcDeltaE2000(GroupPixelCIELAB1931 other) {
       
        CIEChromaticity cieC = new CIEChromaticity();
        
        double deltaE = cieC.calcDeltaECIE2000(
            avgL, avgA, avgB, other.avgL, other.avgA, other.avgB);
    
        return (float)deltaE;
    }
    
    /**
     * calculate the absolute difference in L, U, V between this and other
     * 
     * @param other
     * @return 
     */
    public float[] calcDifference(GroupPixelCIELAB1931 other) {
       
        float[] diffs = new float[3];
        diffs[0] = Math.abs(avgL - other.getAvgL());
        diffs[1] = Math.abs(avgA - other.getAvgA());
        diffs[2] = Math.abs(avgB - other.getAvgB());
        
        return diffs;
    }

    /**
     * @return the avgL
     */
    public float getAvgL() {
        return avgL;
    }

    /**
     * @return the stdDevL
     */
    public float getStdDevL() {
        return stdDevL;
    }

    /**
     * @return the avgA
     */
    public float getAvgA() {
        return avgA;
    }

    /**
     * @return the stdDevA
     */
    public float getStdDevA() {
        return stdDevA;
    }

    /**
     * @return the avgB
     */
    public float getAvgB() {
        return avgB;
    }

    /**
     * @return the stdDevB
     */
    public float getStdDevB() {
        return stdDevB;
    }
}
