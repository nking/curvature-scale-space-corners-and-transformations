package algorithms.imageProcessing;

import algorithms.util.PairInt;
import java.util.Set;

/**
 *
 * @author nichole
 */
public class GroupPixelCIELUV extends GroupPixelRGB0 {
    
    private float avgL;
    private float stdDevL;
    private float avgU;
    private float stdDevU;
    private float avgV;
    private float stdDevV;
    
    public GroupPixelCIELUV(final Set<PairInt> points, ImageExt colorImage) {
               
        if (points == null) {
            throw new IllegalArgumentException("points cannot be null");
        }
        if (colorImage == null) {
            throw new IllegalArgumentException("colorImage cannot be null");
        }
        
        calculateColors(points, colorImage, 0, 0);
        
        CIEChromaticity cieC = new CIEChromaticity();
        float[] cieLUVAvg = cieC.rgbToCIELUV(
            Math.round(getAvgRed()), 
            Math.round(getAvgGreen()), 
            Math.round(getAvgBlue()));
        
        this.nPoints = points.size();
        this.avgL = cieLUVAvg[0];
        this.avgU = cieLUVAvg[1];
        this.avgV = cieLUVAvg[2];
        
        double sumL = 0;
        double sumU = 0;
        double sumV = 0;
        
        for (PairInt p : points) {
            
            int x = p.getX();
            int y = p.getY();
            
            float[] lab = cieC.rgbToCIELUV(
                colorImage.getR(x, y), colorImage.getG(x, y), 
                colorImage.getB(x, y));
            
            float diffL = lab[0] - cieLUVAvg[0];
            float diffU = lab[1] - cieLUVAvg[1];
            float diffV = lab[2] - cieLUVAvg[2];
            
            sumL += (diffL * diffL);
            sumU += (diffU * diffU);
            sumV += (diffV * diffV);
        }
        
        this.stdDevL = (float)Math.sqrt(sumL/(nPoints - 1.));
        this.stdDevU = (float)Math.sqrt(sumU/(nPoints - 1.));
        this.stdDevV = (float)Math.sqrt(sumV/(nPoints - 1.));
    }
    
    /**
     * calculate the difference in L, U, V between this and other and
     * normalize the values to a sum of "1" using the range of values
     * possible from the use of a standard illumant, D65.
     * 
     * @param other
     * @return 
     */
    public float calcNormalizedDifference(GroupPixelCIELUV other) {
       
        /*
        * using the standard illuminant of daylight, D65,
        * the range of return values is
        * L       0 to 104.5
        * u   -86.9 to 183.8
        * v  -141.4 to 112.3
        */
        
        float[] diffs = calcDifference(other);
        diffs[0] /= 104.5f;
        diffs[1] /= (183.8f + 86.9f);
        diffs[2] /= (112.3f + 141.4f);
        
        return (diffs[0] + diffs[1] + diffs[2])/3.f;
    }
    
    /**
     * calculate the difference in L, U, V between this and other 
     * using the DeltaE2000 formula.
     * 
     * @param other
     * @return 
     */
    public float calcDeltaE2000(GroupPixelCIELUV other) {
       
        CIEChromaticity cieC = new CIEChromaticity();
        
        double deltaE = cieC.calcDeltaECIE2000(
            avgL, avgU, avgV, other.avgL, other.avgU, other.avgV);
    
        return (float)deltaE;
    }
    
    /**
     * calculate the absolute difference in L, U, V between this and other
     * 
     * @param other
     * @return 
     */
    public float[] calcDifference(GroupPixelCIELUV other) {
       
        float[] diffs = new float[3];
        diffs[0] = Math.abs(avgL - other.getAvgL());
        diffs[1] = Math.abs(avgU - other.getAvgU());
        diffs[2] = Math.abs(avgV - other.getAvgV());
        
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
     * @return the avgU
     */
    public float getAvgU() {
        return avgU;
    }

    /**
     * @return the stdDevU
     */
    public float getStdDevU() {
        return stdDevU;
    }

    /**
     * @return the avgV
     */
    public float getAvgV() {
        return avgV;
    }

    /**
     * @return the stdDevV
     */
    public float getStdDevV() {
        return stdDevV;
    }
}
