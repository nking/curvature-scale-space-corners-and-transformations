package algorithms.imageProcessing;

import algorithms.misc.Histogram;
import algorithms.misc.HistogramHolder;
import algorithms.misc.MiscMath;
import algorithms.util.Errors;
import algorithms.util.PairInt;
import java.util.Set;

/**
 * cylindrical LUV
 * 
 * @author nichole
 */
public class GroupPixelCIELCH extends GroupPixelRGB0 {
    
    private float avgL;
    private float stdDevL;
    private float avgC;
    private float stdDevC;
    private float avgH;
    private float stdDevH;
    
    public GroupPixelCIELCH(final Set<PairInt> points, ImageExt colorImage) {
               
        if (points == null) {
            throw new IllegalArgumentException("points cannot be null");
        }
        if (colorImage == null) {
            throw new IllegalArgumentException("colorImage cannot be null");
        }
        
        calculateColors(points, colorImage, 0, 0);
        
        CIEChromaticity cieC = new CIEChromaticity();
                
        // compare to indiv determined LCH
        float[] cieLCHAvg = cieC.rgbToCIELCH(
            Math.round(getAvgRed()), Math.round(getAvgGreen()), 
            Math.round(getAvgBlue()));
        
        this.nPoints = points.size();
        
        // the angle component, h, needs wraparound corrections
        //    so frequency has to be calculated for it first.
        
        float[] ls = new float[nPoints];
        float[] cs = new float[nPoints];
        float[] hs = new float[nPoints];
        double sumL = 0;
        double sumC = 0;
        int count = 0;
        for (PairInt p : points) {
            int x = p.getX();
            int y = p.getY();
            float[] lch = cieC.rgbToCIELCH(colorImage.getR(x, y), 
                colorImage.getG(x, y), colorImage.getB(x, y));
            ls[count] = lch[0];
            cs[count] = lch[1];
            hs[count] = lch[2];
            
            sumL += lch[0];
            sumC += lch[1];
            count++;
        }
        this.avgL = (float)(sumL/nPoints);
        this.avgC = (float)(sumC/nPoints);
        
        HistogramHolder hHist = Histogram.createSimpleHistogram(hs, 
            Errors.populateYErrorsBySqrt(hs));
      
        int wrapAround = 255;
        
        int peakIdx = MiscMath.findYMaxIndex(hHist.getYHist());
        float peakValue = hHist.getXHist()[peakIdx];
        double sumH = 0;
        for (int i = 0; i < hs.length; ++i) {
            float v = hs[i];
            float v2 = v + wrapAround; 
            if (Math.abs(v - peakValue) <= Math.abs(v2 - peakValue)) {
                sumH += v;
            } else {
                sumH += v2;
            }
        }
        this.avgH = (float)(sumH/nPoints);
        
        /*
        System.out.println("Lavg indiv = " + avgL);
        System.out.println("Cavg indiv = " + avgC);
        System.out.println("Havg indiv = " + avgH);
        
        System.out.println("Lavg group = " + cieLCHAvg[0]);
        System.out.println("Cavg group = " + cieLCHAvg[1]);
        System.out.println("Havg group = " + cieLCHAvg[2]);
        */
        
        sumL = 0;
        sumC = 0;
        sumH = 0;
        
        for (int i = 0; i < hs.length; ++i) {       
            float diffL = ls[i] - this.avgL;
            float diffC = cs[i] - this.avgC;
            sumL += (diffL * diffL);
            sumC += (diffC * diffC);
            
            float diffH;
            float v = hs[i];
            float v2 = v + wrapAround; 
            if (Math.abs(v - avgH) <= Math.abs(v2 - avgH)) {
                diffH = v - avgH;
            } else {
                diffH = v2 - avgH;
            }            
            sumH += (diffH * diffH);
        }
        
        this.stdDevL = (float)Math.sqrt(sumL/(nPoints - 1.));
        this.stdDevC = (float)Math.sqrt(sumC/(nPoints - 1.));
        this.stdDevH = (float)Math.sqrt(sumH/(nPoints - 1.));
    }
    
    /**
     * calculate the difference in C and H between this and other and
     * normalize the values to a sum of "1" using the range of values
     * possible from the use of a standard illumant, D65.
     * 
     * @param other
     * @return 
     */
    public float calcNormalizedCHDifference(GroupPixelCIELCH other) {
       
        /*
        * using the standard illuminant of daylight, D65,
        * the range of return values is
        *   magnitude, C:  0 to 139 
        *   angle,     H:  0 to 359, but image scales it to 255
        */
        
        float d1 = calcDifferenceH(other)/255.f;
        
        float d2 = calcDifferenceC(other)/139.f;
        
        return 0.5f * (d1 + d2);
    }
    
    public float calcDifferenceC(GroupPixelCIELCH other) {
                
        return Math.abs(avgC - other.avgC);
    }

    public float calcDifferenceH(GroupPixelCIELCH other) {
        
        int wrapAround = 255;
        
        float v1 = avgH;
        float v3 = other.avgH;
        
        if (v1 > v3) {
            float v4 = v3 + wrapAround;
            if (Math.abs(v1 - v3) > Math.abs(v1 - v4)) {
                return Math.abs(v1 - v4);
            } else {
                return Math.abs(v1 - v3);
            }
        } else {
            float v2 = v1 + wrapAround;
            if (Math.abs(v1 - v3) > Math.abs(v2 - v3)) {
                return Math.abs(v2 - v3);
            } else {
                return Math.abs(v1 - v3);
            }
        }
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
     * @return the avgC
     */
    public float getAvgC() {
        return avgC;
    }

    /**
     * @return the stdDevC
     */
    public float getStdDevC() {
        return stdDevC;
    }

    /**
     * @return the avgH
     */
    public float getAvgH() {
        return avgH;
    }

    /**
     * @return the stdDevH
     */
    public float getStdDevH() {
        return stdDevH;
    }
}
