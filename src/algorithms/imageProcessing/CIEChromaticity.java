package algorithms.imageProcessing;

import algorithms.util.ArrayPair;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;

/**
 *
 * see http://hyperphysics.phy-astr.gsu.edu/hbase/vision/cie.html
 * and 
 * http://en.wikipedia.org/wiki/CIE_1931_color_space#mediaviewer/File:CIE1931xy_blank.svg
 * 
 * @author nichole
 */
public class CIEChromaticity {
    
    /**
     * get the bounds of yellow in the CIE 1931 xy chromaticity diagram
     * as x and y coordinates.  Note that the last and first point are
     * the same.
     * 
     * @return 
     */
    public ArrayPair getYellowPolynomial() {
        /*
        516, 486
        488, 512  
        436, 450
        327, 340 <=== 346, 360
        466, 444
        */
        
        ArrayPair p = new ArrayPair(
            new float[]{.516f, .488f, .436f, .346f, .466f, .516f},
            new float[]{.486f, .512f, .450f, .360f, .444f, .486f}
        );
        
        return p;
    }
    
    public ArrayPair getGreenThroughYellowGreenPolynomial() {
        /*
        0.0125, 0.4827
        0.3, 0.35
        0.36,  0.35
        0.45, 0.55
        0.2, 0.8
        0.1, 0.84
        0.0, 0.8
        0.0, 0.4827
        0.0125, 0.4827
        */
        
        ArrayPair p = new ArrayPair(
            new float[]{.0125f, .3f, .36f, .45f, .2f, .1f, 0.f, 0.f, .0125f},
            new float[]{.4827f, .35f, .35f, .55f, .8f, .84f, .8f, .4827f, .4827f}
        );
        
        return p;
    }
    
    public ArrayPair getYellowThroughOrangePolynomial() {
        
        /*
        685, 312
        488, 512
        435, 454
                  343, 373
        342, 344
        420, 340
        550, 349
        */
        
        ArrayPair p = new ArrayPair(
            new float[]{.685f, .488f, .435f, .343f, .342f, .420f, .550f, .685f},
            new float[]{.312f, .512f, .454f, .373f, .344f, .340f, .349f, .312f}
        );
        
        return p;
    }
    
    public ArrayPair getYellowishGreenThroughOrangePolynomial() {
        /*
        685, 312
        456, 544
        393, 448
                  340, 373
        342, 344
        420, 340
        550, 349
        */
        
        ArrayPair p = new ArrayPair(
            new float[]{.685f, .456f, .39f, .340f, .342f, .420f, .550f, .685f},
            new float[]{.312f, .544f, .448f, .373f, .344f, .340f, .349f, .312f}
        );
        
        return p;
    }
            
    public ArrayPair getRedThroughPurplishRedPolynomial() {
        
        /*
        690, 310
        552, 346
        422, 341    
        337, 344  339, 326
        554, 182
        738, 264
        */
        
        ArrayPair p = new ArrayPair(
            new float[]{.688f, .552f, .422f, .337f, .339f, .554f, .738f, .688f},
            new float[]{.313f, .346f, .341f, .344f, .326f, .182f, .264f, .313f}
        );
        
        return p;
    }
    
    public ArrayPair getRedPolynomial() {
        
        /*
        690, 310
        552, 346
        422, 341    382, 342
        
        527, 291
        659, 229
        738, 264
        */
        
        ArrayPair p = new ArrayPair(
            new float[]{.688f, .552f, .422f, .367f, .527f, .659f, .738f, .688f},
            new float[]{.313f, .346f, .341f, .342f, .291f, .229f, .264f, .313f}
        );
        
        return p;
    }
    
    public ArrayPair getOrangePolynomial() {
        
        /*
        558, 439 <== 515, 482
        490, 416 <== 467, 442
        386, 365 <== 331,341
        420, 340
        550, 349
        687, 314
        */
        
        ArrayPair p = new ArrayPair(
            new float[]{.515f, .467f, .331f, .420f, .550f, .687f, .515f},
            new float[]{.482f, .442f, .341f, .340f, .349f, .314f, .482f}
        );
        
        return p;
    }
    
    public ArrayPair getPurplePolynomial() {
        
        /*
        
        including center is larger polygon w/ fewer points:
        168,   3
        237, 193
        300, 324
        340, 324
        552, 180
        
        */
        
        ArrayPair p = new ArrayPair(
            new float[]{.168f, .237f, .300f, .340f, .552f, .168f},
            new float[]{0.003f, .193f, .324f, .324f, .180f, .003f}
        );
        
        return p;
    }
    
    /**
     * convert from CIE XYZ 1931, to CIE XY chromaticity 1931.
     * 
     * useful for more information:
     * http://www.efg2.com/Lab/Graphics/Colors/Chromaticity.htm
     * http://en.wikipedia.org/wiki/CIE_1931_color_space
     * http://hyperphysics.phy-astr.gsu.edu/hbase/vision/cie.html
     * 
     * @param r
     * @param g
     * @param b
     * @return 
     */
    public float[] rgbToXYChromaticity(int r, int g, int b) {
        
        float[] capXYZ = rgbToCIEXYZ(r, g, b);
        
        float x = capXYZ[0]/(capXYZ[0] + capXYZ[1] + capXYZ[2]);
        
        float y = capXYZ[1]/(capXYZ[0] + capXYZ[1] + capXYZ[2]);
        
        return new float[]{x, y};
    }
    
    /**
     * convert rgb to CIE XYZ (1931).
     * 
     * uses http://en.wikipedia.org/wiki/CIE_1931_color_space#Experimental_results:_the_CIE_RGB_color_space
     * 
     * @param r
     * @param g
     * @param b
     * @return 
     */
    public float[] rgbToCIEXYZ(int r, int g, int b) {
        
        /*        
            | X |       1     | 0.49     0.31      0.20     |   | R |
            | Y | = --------- | 0.17697  0.81240   0.01063  | * | G |
            | Z |    0.17697  | 0.00     0.01      0.99     |   | B |
        */
        
        float capX = (0.49f * r +  0.31f * g + 0.20f * b)/0.17697f;
        
        float capY = (0.17697f * r +  0.81240f * g + 0.01063f * b)/0.17697f;
        
        float capZ = (0.01f * g + 0.99f * b)/0.17697f;
        
        return new float[]{capX, capY, capZ};
    }

    public List<Double> calcAvgAndStdDevXY(int[] r, int[] g, int[] b) {
        
        double xSum = 0;
        double ySum = 0;
        float[] x = new float[r.length];
        float[] y = new float[r.length];
        for (int i = 0; i < r.length; i++) {
            float[] xy = rgbToXYChromaticity(r[i], g[i], b[i]);
            x[i] = xy[0];
            y[i] = xy[1];
            xSum += xy[0];
            ySum += xy[1];
            i++;
        }
        double avgX = xSum/(double)r.length;
        double avgY = ySum/(double)r.length;
        
        xSum = 0;
        ySum = 0;
        for (int i = 0; i < r.length; i++) {
            double diffX = x[i] - avgX;
            double diffY = y[i] - avgY;
            xSum += (diffX * diffX);
            ySum += (diffY * diffY);
        }
        double stDevX = Math.sqrt(xSum/((double)r.length - 1));
        double stDevY = Math.sqrt(ySum/((double)r.length - 1));
        
        List<Double> list = new ArrayList<Double>();
        list.add(Double.valueOf(avgX));
        list.add(Double.valueOf(avgY));
        list.add(Double.valueOf(stDevX));
        list.add(Double.valueOf(stDevY));
        
        return list;
    }
}
