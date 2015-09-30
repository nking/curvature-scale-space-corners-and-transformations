package algorithms.imageProcessing;

import algorithms.misc.MiscMath;
import algorithms.util.ArrayPair;
import java.util.ArrayList;
import java.util.List;

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
     * the offset from (0.35, 0.35) considered "white"
     */
    private static final float deltaWhite = 0.0125f;
    
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
    
    public ArrayPair getYellowishGreenThroughYellowGreenPolynomial() {
        
        ArrayPair p = new ArrayPair(
            new float[]{0.18f, 0.45f, 0.45f, 0.32f, 0.31f, 0.18f},
            new float[]{0.8f,  0.8f, 0.55f, 0.32f,  0.4f,  0.8f}
        );
        
        return p;
    }
    
    public ArrayPair getGreenishYellowThroughOrangePolynomial() {
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
    
    public ArrayPair getGreenPolynomial() {
        
        ArrayPair p = new ArrayPair(
            new float[]{0.03f, 0.00f,  0.0f,  0.35f, 0.35f, 0.03f},
            new float[]{0.37f, 0.485f, 0.83f, 0.83f, 0.35f, 0.37f}
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
        
        if ((r == 0) && (g == 0) && (b == 0)) {
            // not really defined on the diagram since chromaticity is color
            // without intensity.  since all 0's is the lack of all color
            // will return 0,0, but it's N/A
            return new float[]{0, 0};
        }
        
        float[] capXYZ = rgbToCIEXYZ(r, g, b);
        
        float zz = capXYZ[0] + capXYZ[1] + capXYZ[2];
        
        float x = capXYZ[0]/zz;
        
        float y = capXYZ[1]/zz;
        
        return new float[]{x, y};
    }
    
    private float[] cieXYTmpHolder = new float[3];
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
     * @param outputcieXY the float array of length 2 that the output will be 
     * placed in.
     */
    public void rgbToXYChromaticity(int r, int g, int b, float[] outputcieXY) {
        
        if (outputcieXY == null || outputcieXY.length != 2) {
            throw new IllegalArgumentException("outputcieXY must be length 2");
        }
        
        if ((r == 0) && (g == 0) && (b == 0)) {
            // not really defined on the diagram since chromaticity is color
            // without intensity.  since all 0's is the lack of all color
            // will return 0,0, but it's N/A
            outputcieXY[0] = 0;
            outputcieXY[1] = 0;
            return;
        }
        
        rgbToCIEXYZ(r, g, b, cieXYTmpHolder);
        
        float x = cieXYTmpHolder[0]/(cieXYTmpHolder[0] + cieXYTmpHolder[1] + cieXYTmpHolder[2]);
        
        float y = cieXYTmpHolder[1]/(cieXYTmpHolder[0] + cieXYTmpHolder[1] + cieXYTmpHolder[2]);
        
        outputcieXY[0] = x;
        outputcieXY[1] = y;
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

    /**
     * convert rgb to CIE XYZ (1931).
     * 
     * uses http://en.wikipedia.org/wiki/CIE_1931_color_space#Experimental_results:_the_CIE_RGB_color_space
     * 
     * @param r
     * @param g
     * @param b
     * @param outputCIEXYZ output array of length 3 that will be populated with
     * cieX, cieY and cieZ.
     */
    public void rgbToCIEXYZ(int r, int g, int b, float[] outputCIEXYZ) {
        
        if (outputCIEXYZ == null || outputCIEXYZ.length != 3) {
            throw new IllegalArgumentException("outputCIEXYZ has to be length 3");
        }
        
        /*        
            | X |       1     | 0.49     0.31      0.20     |   | R |
            | Y | = --------- | 0.17697  0.81240   0.01063  | * | G |
            | Z |    0.17697  | 0.00     0.01      0.99     |   | B |
        */
        
        float capX = (0.49f * r +  0.31f * g + 0.20f * b)/0.17697f;
        
        float capY = (0.17697f * r +  0.81240f * g + 0.01063f * b)/0.17697f;
        
        float capZ = (0.01f * g + 0.99f * b)/0.17697f;
        
        outputCIEXYZ[0] = capX;
        outputCIEXYZ[1] = capY;
        outputCIEXYZ[2] = capZ;
    }

    /**
     * convert CIE XYZ (1931) to rgb.
     * 
     * uses http://en.wikipedia.org/wiki/CIE_1931_color_space#Experimental_results:_the_CIE_RGB_color_space
     * 
     * @param cieX
     * @param cieY
     * @param cieZ
     * @return 
     */
    public int[] cieXYZToRGB(float cieX, float cieY, float cieZ) {
        
        /*        
            | R |   |  0.41847    -0.15866   -0.082835 |   | X |
            | G | = | -0.091169    0.25243    0.015708 | * | Y |
            | B |   |  0.00092090 -0.0025498  0.17860  |   | Z |
        */
        
        float capR = (0.41847f * cieX +  -0.15866f * cieY + -0.082835f * cieZ);
        
        float capG = (-0.091169f * cieX +  0.25243f * cieY + 0.015708f * cieZ);
        
        float capB = (0.00092090f * cieX +  -0.0025498f * cieY + 0.17860f * cieZ);
        
        return new int[]{(int)capR, (int)capG, (int)capB};
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
    
    /**
     * returns roughly whether the CIE (X,Y) coordinate lands within the large
     * white zone in the center of the diagram.
     * @param cieX
     * @param cieY
     * @return 
     */
    public boolean isWhite(float cieX, float cieY) {
        
        double diffX = Math.abs(cieX - 0.35);
        
        double diffY = Math.abs(cieY - 0.35);
      
        return ((diffX <= deltaWhite) && (diffY <= deltaWhite));
    }
    
    /**
     * returns roughly whether the CIE (X,Y) coordinate lands within the large
     * white zone in the center of the diagram.
     * @param cieX
     * @param cieY
     * @return 
     */
    public boolean isWhite2(float cieX, float cieY) {
        
        double diffX = Math.abs(cieX - 0.35);
        
        double diffY = Math.abs(cieY - 0.35);
      
        double delta = 2.*(0.1/7.);
        
        return ((diffX <= delta) && (diffY <= delta));        
    }
    
    /**
     * calculate the angle in radians of the point (cieX, cieY) in the 
     * CIE chromaticity
     * diagram with respect to an origin of (0.35, 0.35).
     * Note that one should check for white before using this as the resulting
     * angle will not be a significant answer if it is.
     * The angles are
     * <pre>    90(=pi/2)
     *            |
     *            |
     *   180 ----------- 0
     *   (=pi)    |
     *            | 
     *          270(=3pi/2)
     * </pre>
     * @param cieX
     * @param cieY
     * @return the angle in radians of the point (cieX, cieY) with respect to
     * an origin of (0.35, 0.35). The invoker should check the original r,g,b
     * for black and white to avoid use here.  A value of -1 is returned for
     * white.
     */
    public double calculateXYTheta(float cieX, float cieY) {
                        
        if (cieX == 0.35) {
            if (cieY > (0.35 + deltaWhite)) {
                return Math.PI/2.;
            } else if (cieY < (0.35 - deltaWhite)) {
                return 1.5 * Math.PI;
            } else {
                // this should have been found as white to avoid use here
                return -1;
            }
        }
        
        double theta = MiscMath.calculatePolarTheta(cieX - 0.35f, cieY - 0.35f);
                    
        return theta;
    }
    
}
