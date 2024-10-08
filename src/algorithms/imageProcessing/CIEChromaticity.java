package algorithms.imageProcessing;

import algorithms.misc.MiscMath;
import java.util.ArrayList;
import java.util.List;

/**
 *  convenience methods for color space not present
    in jdk.
 * see http://hyperphysics.phy-astr.gsu.edu/hbase/vision/cie.html
 * and
 * http://en.wikipedia.org/wiki/CIE_1931_color_space#mediaviewer/File:CIE1931xy_blank.svg
 *
 * and
 *
 * @author nichole
 */
public class CIEChromaticity {

    /**
     * the offset from (0.35, 0.35) considered "white"
     */
    private static final float deltaWhite = 0.0125f;

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

        if (zz == 0.f) {
            return new float[]{0, 0};
        }

        float x = capXYZ[0]/zz;

        float y = capXYZ[1]/zz;

        return new float[]{x, y};
    }

    /**
     * convert from CIE XYZ 1931, to CIE XY chromaticity 1931.
     *
     * useful for more information:
     * http://www.efg2.com/Lab/Graphics/Colors/Chromaticity.htm
     * http://en.wikipedia.org/wiki/CIE_1931_color_space
     * http://hyperphysics.phy-astr.gsu.edu/hbase/vision/cie.html
     *
     * expects r, g, b in range 0 to 255, inclusive.
     * @param r
     * @param g
     * @param b
     * @return
     */
    public float[] _rgbToXYChromaticity(int r, int g, int b) {

        if ((r == 0) && (g == 0) && (b == 0)) {
            // not really defined on the diagram since chromaticity is color
            // without intensity.  since all 0's is the lack of all color
            // will return 0,0, but it's N/A
            return new float[]{0, 0};
        }
        
        return _rgbToXYChromaticity((float)r/255.f, (float)g/255.f, (float)b/255.f);
    }

    /**
     * convert from CIE XYZ 1931, to CIE XY chromaticity 1931.
     *
     * useful for more information:
     * http://www.efg2.com/Lab/Graphics/Colors/Chromaticity.htm
     * http://en.wikipedia.org/wiki/CIE_1931_color_space
     * http://hyperphysics.phy-astr.gsu.edu/hbase/vision/cie.html
     *
     * expects r, g, b in range 0 to 1.0, inclusive.
     * @param r
     * @param g
     * @param b
     * @return
     */
    public float[] _rgbToXYChromaticity(float r, float g, float b) {

        if ((r == 0) && (g == 0) && (b == 0)) {
            // not really defined on the diagram since chromaticity is color
            // without intensity.  since all 0's is the lack of all color
            // will return 0,0, but it's N/A
            return new float[]{0, 0};
        }

        float[] capXYZ = _rgbToCIEXYZ(r, g, b);

        float zz = capXYZ[0] + capXYZ[1] + capXYZ[2];

        if (zz == 0.f) {
            return new float[]{0, 0};
        }

        float x = capXYZ[0]/zz;

        float y = capXYZ[1]/zz;

        return new float[]{x, y};
    }

    // NOTE: not thread safe:
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
     * convert rgb to CIE LAB.
     *
     * uses http://en.wikipedia.org/wiki/CIE_1931_color_space#Experimental_results:_the_CIE_RGB_color_space
     * and http://en.wikipedia.org/wiki/Lab_color_space#Forward_transformation
     *
     * range of return values is
     * <pre>
     *    L    0 to 28.5
     *    A  -46.9  62.5
     *    B  -45.7  48.0
     * </pre>
     *
     * @param r
     * @param g
     * @param b
     * @return
     */
    public float[] rgbToCIELAB(int r, int g, int b) {

        if (r < 0 || r > 255 || g < 0 || g > 255 || b < 0 || b > 255) {
            throw new IllegalArgumentException("r, g, and b must be 0 to 255");
        }
        if ((r == 0) && (g == 0) && (b == 0)) {
            return new float[]{0, 0, 0};
        }
        
        float[] a = _rgbToCIEXYZ(r, g, b);

        a = cieXYZToCIELAB(a);
        
        return a;
    }

    /**
     * convert rgb to CIE LAB 1931.
     *
     * uses http://en.wikipedia.org/wiki/CIE_1931_color_space#Experimental_results:_the_CIE_RGB_color_space
     * and http://en.wikipedia.org/wiki/Lab_color_space#Forward_transformation
     *
     * range of return values for
     * L* 0 to 105
     * a* -190 to 103
     * b* -113 to 99
     *
     * @param r
     * @param g
     * @param b
     * @return
     */
    public float[] rgbToCIELAB1931(int r, int g, int b) {

        if ((r == 0) && (g == 0) && (b == 0)) {
            return new float[]{0, 0, 0};
        }
        
        float[] a = _rgbToCIEXYZ2((float)r/255.f, (float)g/255.f, (float)b/255.f);

        a = cieXYZToCIELAB(a);

        return a;
    }

    /**
     * convert rgb to CIE LCH, polar coordinates of CIELAB1931.
     *
     * uses http://en.wikipedia.org/wiki/CIE_1931_color_space#Experimental_results:_the_CIE_RGB_color_space
     * and http://en.wikipedia.org/wiki/Lab_color_space#Forward_transformation
     *
     * range of return values for
     * L* 0 to 105
     * c  0 to 139
     * h  0 to 359
     *
     * @param r
     * @param g
     * @param b
     * @return
     */
    public float[] rgbToCIELCH(int r, int g, int b) {

        if ((r == 0) && (g == 0) && (b == 0)) {
            return new float[]{0, 0, 0};
        }
        
        float[] xyz = _rgbToCIEXYZ2((float)r/255.f, (float)g/255.f, (float)b/255.f);

        float[] lch = cieXYZToCIELCH(xyz);

        return lch;
    }

    /**
     * convert rgb to polar coordinates of CIELUV (a.k.a. CIEL*A*B* 1976?)
     *
     * CIE 1976 (L*, u*, v*).
     * 
     * uses http://en.wikipedia.org/wiki/CIE_1931_color_space#Experimental_results:_the_CIE_RGB_color_space
     * and http://en.wikipedia.org/wiki/Lab_color_space#Forward_transformation
     *
     * <pre>
     * using the standard illuminant of daylight, D65,
     * the range of return values is
     range of CIE LUV using default standard illumination of
        D65 daylight is:
        L       0 to 104.5
        u   -86.9 to 183.8
        v  -141.4 to 112.3
        luminosity L*  0 to 104.5
        magnitude, m:  sqrt(2) * 183.8 = 260
        angle,     a:  0 to 359
     * </pre>
     * @param r in range 0 to 255, inclusive
     * @param g in range 0 to 255, inclusive
     * @param b in range 0 to 255, inclusive
     * @return luv[0], (float)m, (float)t where
     * the first item is luminosity, the
     * second is the magnitude of the color, that
     * is sqrt of square sums of U and V,
     * and the third is the polar theta angle
     */
    public float[] rgbToPolarCIELUV(int r, int g, int b) {

        if ((r == 0) && (g == 0) && (b == 0)) {
            return new float[]{0, 0, 0};
        }
        
        float[] xyz = _rgbToCIEXYZ2((float)r/255.f, (float)g/255.f, (float)b/255.f);

        float[] luv = cieXYZToCIELUV(xyz);

        double t = Math.atan2(luv[2], luv[1]);
        t *= (180. / Math.PI);
        if (t < 0) {
            t += 360;
        } else if (t > 359) {
            t -= 360;
        }

        double m = Math.sqrt(luv[1]*luv[1] + luv[2]*luv[2]);

        return new float[]{luv[0], (float)m, (float)t};
    }

    /**
     * convert rgb to polar coordinates of CIELUV.
     *  CIE 1976 (L*, u*, v*) 
     * 
     * uses http://en.wikipedia.org/wiki/CIE_1931_color_space#Experimental_results:_the_CIE_RGB_color_space
     * and http://en.wikipedia.org/wiki/Lab_color_space#Forward_transformation
     *
     * range of return values for
     * L*
     * magnitude:
     * angle:      0 to 359
     *
     * @param r
     * @param g
     * @param b
     * @param Xn standard illuminant
     * @param Yn standard illuminant
     * @param Zn standard illuminant
     * @return
     */
    public float[] rgbToPolarCIELUV(int r, int g, int b, float Xn, float Yn,
        float Zn) {

        if ((r == 0) && (g == 0) && (b == 0)) {
            return new float[]{0, 0, 0};
        }
        
        float[] xyz = _rgbToCIEXYZ2((float)r/255.f, (float)g/255.f, (float)b/255.f);

        float[] luv = cieXYZToCIELUV(xyz, Xn, Yn, Zn);

        double t = Math.atan2(luv[2], luv[1]);
        t *= (180. / Math.PI);
        if (t < 0) {
            t += 360;
        } else if (t > 359) {
            t -= 360;
        }

        double m = Math.sqrt(luv[1]*luv[1] + luv[2]*luv[2]);

        return new float[]{luv[0], (float)m, (float)t};
    }
    
    /**
     * convert rgb to CIELUV. CIE 1976 (L*, u*, v*) 
     * 
     * uses http://en.wikipedia.org/wiki/CIE_1931_color_space#Experimental_results:_the_CIE_RGB_color_space
     * and http://en.wikipedia.org/wiki/Lab_color_space#Forward_transformation
     *
     * uses standard illuminant "wide range lightness".
     * Wide-range Lightness data were generated by Fairchild et al., 
     * who conducted two different experiments to scale lightness above and 
     * below diffuse white (CIE * L = 100 ). 
     * In the Scaling Lightness Experiment 1 (SL1) they used a luminance
       range from 156 to 2 3692cd m with 2 842 Y cd m n = 
       (Yn represents the luminance of reference white) 
      whereas in the Scaling Lightness Experiment 2 (SL2) the
      luminance range was extended from 0 to 2 7432cd m with 2 997 Y cd m n = . 
     The SL2 data set was used to drive the adapted lightness ( z J ) 
     formula of the proposed color space (see later) while the SL1 data set was 
     used as a test data set. Each of the sets includes 19 samples. 

     * https://www.osapublishing.org/DirectPDFAccess/B810E9AE-C7C9-E594-C72DC7FBE1424F0A_368272/oe-25-13-15131.pdf?da=1&id=368272&seq=0&mobile=no
     * 
     * SL2 Training D65/2° (x,y,z)=968.08 997 883.51 
     * L_a=199 
     * C=0.69. N_C=1, F=1 
     * 
     * range of return values is:
     * L       0 to 104.8
     * u   -99.9 to 176.8
     * v  -150.0 to 95.0
     * 
     * @param r in range 0 to 255, inclusive
     * @param g in range 0 to 255, inclusive
     * @param b in range 0 to 255, inclusive
     * @return
     */
    public float[] rgbToCIELUV_WideRangeLightness(int r, int g, int b) {

        if ((r == 0) && (g == 0) && (b == 0)) {
            return new float[]{0, 0, 0};
        }
        
        float Xn = 96.808f;//0.96808f;
        float Yn = 99.7f;//0.997f;
        float Zn = 88.35f;//0.88351f;
        
        float[] xyz = _rgbToCIEXYZ2((float)r/255.f, (float)g/255.f, (float)b/255.f);

        float[] luv = cieXYZToCIELUV(xyz, Xn, Yn, Zn);

        return luv;
    }

    /**
     * convert rgb to CIE LUV (a.k.a. CIEL*A*B* 1976)
     * CIE 1976 (L*, u*, v*) 
     * 
     * Note that differences in CIELUV are simply the differences
     * in each component added in quadrature (no deltaE formula).
     * 
     * uses http://en.wikipedia.org/wiki/CIE_1931_color_space#Experimental_results:_the_CIE_RGB_color_space
     * and http://en.wikipedia.org/wiki/Lab_color_space#Forward_transformation
     * http://www.easyrgb.com/index.php?X=MATH&H=16#text16
     *
     * range of return values when using default standard illumination of
     * D65 daylight is:
     * L       0 to 104.5
     * u   -86.9 to 183.8
     * v  -141.4 to 112.3
     *
     * @param r
     * @param g
     * @param b
     * @return
     */
    public float[] rgbToCIELUV(int r, int g, int b) {

        if ((r == 0) && (g == 0) && (b == 0)) {
            return new float[]{0, 0, 0};
        }
        
        float[] a = _rgbToCIEXYZ2((float)r/255.f, (float)g/255.f, (float)b/255.f);

        a = cieXYZToCIELUV(a);

        return a;
    }

    /**
     * convert rgb to CIE LUV (a.k.a. CIEL*A*B* 1976?)
     * CIE 1976 (L*, u*, v*) 
     * and use the given standard illumination in tristimulus coordinates.
     * 
      <pre>
      Incandescent:
          109.850, 100, 35.585
          range of return values for incandescent
           L       0 to 104.5
           u  -156.2 to 141.6
           v  -172.0 to 45.6
      Daylight, midday (D65):
          95.047, 100, 108.883
          range of return values for D65
           L       0 to 104.5
           u   -86.9 to 183.8
           v  -141.4 to 112.3
      Fluorescent:
          99.187, 100, 67.395
          range of return values for Fluorescent
           L       0 to 104.5
           u  -113.6 to 167.5
           v  -158.5 to 75.0
      D75:
          94.972,  100, 122.638
          range of return values for D75
           L       0 to 104.5
           u   -81.7 to 186.9
           v  -136.2 to 124.5
      </pre>
     * uses http://en.wikipedia.org/wiki/CIE_1931_color_space#Experimental_results:_the_CIE_RGB_color_space
     * and http://en.wikipedia.org/wiki/Lab_color_space#Forward_transformation
     * http://www.easyrgb.com/index.php?X=MATH&H=16#text16
     
     * @param r
     * @param g
     * @param b
     * @param Xn standard illuminant X
     * @param Yn standard illuminant Y
     * #param Zn standard illuminant Z
     * @return
     */
    public float[] rgbToCIELUV(int r, int g, int b,
        float Xn, float Yn, float Zn) {

        if ((r == 0) && (g == 0) && (b == 0)) {
            return new float[]{0, 0, 0};
        }
        
        float[] a = _rgbToCIEXYZ2((float)r/255.f, (float)g/255.f, (float)b/255.f);

        a = cieXYZToCIELUV(a, Xn, Yn, Zn);

        return a;
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

        0.436  0.385  0.143
        0.222  0.7169 0.0606
        0.0139 0.097  0.7139
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
     * normalizes the r,g,b values to lie between 0 and 1.0 assuming the
     * arguments are given as a range 0 to 255.
     *
     * for r=0, g=0, b=0, CIEXY is (0, 0, 0).
     * for r=1, g=1, b=1, CIEXY is (5.65, 5.65, 5.65)
     * @param r
     * @param g
     * @param b
     * @return
     */
    public float[] _rgbToCIEXYZ(int r, int g, int b) {

        if ((r == 0) && (g == 0) && (b == 0)) {
            return new float[]{0, 0, 0};
        }
        
        return _rgbToCIEXYZ((float)r/255.f, (float)g/255.f, (float)b/255.f);
    }

    /**
     * convert rgb to CIE XYZ (1931).
     *
     * uses http://en.wikipedia.org/wiki/CIE_1931_color_space#Experimental_results:_the_CIE_RGB_color_space
     *
     * expects r,g,b values between 0 and 1, inclusive.
     *
     * for r=0, g=0, b=0, CIEXY is (0, 0, 0).
     * for r=1, g=1, b=1, CIEXY is (5.65, 5.65, 5.65)
     * @param r
     * @param g
     * @param b
     * @return
     */
    public float[] _rgbToCIEXYZ(float r, float g, float b) {

        /*
            | X |       1     | 0.49     0.31      0.20     |   | R |
            | Y | = --------- | 0.17697  0.81240   0.01063  | * | G |
            | Z |    0.17697  | 0.00     0.01      0.99     |   | B |

        D50 matrix:
        0.436  0.385  0.143
        0.222  0.7169 0.0606
        0.0139 0.097  0.7139
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
     * expects r,g,b values between 0 and 1, inclusive.
     *
     * for r=0, g=0, b=0, CIEXY is (0, 0, 0).
     * for r=1, g=1, b=1, CIEXY is (5.65, 5.65, 5.65)
     * @param r
     * @param g
     * @param b
     * @return
     */
    public float[] _rgbToCIEXYZ2(float r, float g, float b) {

        if ((r == 0) && (g == 0) && (b == 0)) {
            return new float[]{0, 0, 0};
        }
        
        //http://www.easyrgb.com/index.php?X=MATH&H=02#text2

        float fR = 0;
        float fG = 0;
        float fB = 0;
        for (int i = 0; i < 3; ++i) {
            double d;
            if (i == 0) {
                d = r;
            } else if (i == 1) {
                d = g;
            } else {
                d = b;
            }
            double f;
            if (d > 0.04045) {
                f = Math.pow(d + 0.055/1.055, 2.4);
            } else {
                f = d/12.92;
            }
            f *= 100.f;
            if (i == 0) {
                fR = (float)f;
            } else if (i == 1) {
                fG = (float)f;
            } else {
                fB = (float)f;
            }
        }

        float capX = fR * 0.4124f + fG * 0.3576f + fB * 0.1805f;
        float capY = fR * 0.2126f + fG * 0.7152f + fB * 0.0722f;
        float capZ = fR * 0.0193f + fG * 0.1192f + fB * 0.9505f;

        return new float[]{capX, capY, capZ};
    }

    /**
     * convert CIE XYZ (1931) to CIE LAB.
     *
     * uses https://en.wikipedia.org/wiki/Lab_color_space#Forward_transformation
     *
     * L is from 0 to 100
     * a is from - to +  (negative are green, pos are red)
     * b is from - to +  (negative are blue, pos are yellow)
     *
     * to change the reference point, see
     * http://www.easyrgb.com/index.php?X=MATH&H=15#text15
     *
     *
     * @param cieXYZ
     * @return
     */
    public float[] cieXYZToCIELAB(float[] cieXYZ) {

        /*
        Incandescent:
            109.850, 100, 35.585
        Daylight, midday (D65):
            95.047, 100, 108.883
        Fluorescent:
            99.187, 100, 67.395

        D75:
           94.972,  100, 122.638


        would like tristimulus colors for a range of cloud and shade
        conditions.  Mie scattering for clouds and Rayleigh scattering for air,
        then consideration for aerosols... these aren't single scattering
        conditions, so empirically gathered data might be more useful to derive
        ranges from...many interesting databases with relevant data...
        */

        float Xn = 95.047f;
        float Yn = 100.0f;
        float Zn = 108.883f;

        return cieXYZToCIELAB(cieXYZ, Xn, Yn, Zn);
    }

    /**
     * convert CIE XYZ (1931) to cylindrical LCH
     * where H is the polar angle between V and U and C is the
     * magnitude of that vector.
     *
     * uses https://en.wikipedia.org/wiki/Lab_color_space#Forward_transformation
     *
     * http://www.easyrgb.com/index.php?X=MATH&H=15#text15
     *
     * @param cieXYZ
     * @return
     */
    public float[] cieXYZToCIELCH(float[] cieXYZ) {

        /*
        Incandescent:
            109.850, 100, 35.585
        Daylight, midday (D65):
            95.047, 100, 108.883
        Fluorescent:
            99.187, 100, 67.395

        D75:
           94.972,  100, 122.638


        would like tristimulus colors for a range of cloud and shade
        conditions.  Mie scattering for clouds and Rayleigh scattering for air,
        then consideration for aerosols... these aren't single scattering
        conditions, so empirically gathered data might be more useful to derive
        ranges from...many interesting databases with relevant data...
        */

        float Xn = 95.047f;
        float Yn = 100.0f;
        float Zn = 108.883f;

        float[] lab = cieXYZToCIELAB(cieXYZ, Xn, Yn, Zn);

        double t = Math.atan2(lab[2], lab[1]);
        t *= (180. / Math.PI);
        if (t < 0) {
            t += 360;
        } else if (t > 359) {
            t -= 360;
        }

        double m = Math.sqrt(lab[1]*lab[1] + lab[2]*lab[2]);

        return new float[]{lab[0], (float)m, (float)t};
    }

    /**
     * convert CIE XYZ (1931) to CIE LUV.
     * CIE 1976 (L*, u*, v*) 
     * 
     * uses https://en.wikipedia.org/wiki/Lab_color_space#Forward_transformation
     *
     * http://www.easyrgb.com/index.php?X=MATH&H=15#text15
     *
     * @param cieXYZ
     * @return
     */
    public float[] cieXYZToCIELUV(float[] cieXYZ) {

        /*
        Incandescent:
            109.850, 100, 35.585
        Daylight, midday (D65):
            95.047, 100, 108.883
        Fluorescent:
            99.187, 100, 67.395

        D75:
           94.972,  100, 122.638


        would like tristimulus colors for a range of cloud and shade
        conditions.  Mie scattering for clouds and Rayleigh scattering for air,
        then consideration for aerosols... these aren't single scattering
        conditions, so empirically gathered data might be more useful to derive
        ranges from...many interesting databases with relevant data...
        */

        float Xn = 95.047f;
        float Yn = 100.0f;
        float Zn = 108.883f;

        return cieXYZToCIELUV(cieXYZ, Xn, Yn, Zn);
    }

    /**
     * convert CIE XYZ (1931) to CIE LUV.
     * CIE 1976 (L*, u*, v*) 
     * 
     * uses https://en.wikipedia.org/wiki/Lab_color_space#Forward_transformation
     *
     * http://www.easyrgb.com/index.php?X=MATH&H=15#text15
     *
     * @param cieXYZ
     * @param Xn standard illuminant X
     * @param Yn standard illuminant Y
     * @param Zn standard illuminant Z
     * @return
     */
    public float[] cieXYZToCIELUV(float[] cieXYZ, float Xn, float Yn, float Zn) {

        float[] lab = cieXYZToCIELAB(cieXYZ, Xn, Yn, Zn);

        if (lab[0] == 0 && lab[1] == 0 && lab[2] == 0) {
            return lab;
        }
        
        float u = (4.f * cieXYZ[0]) /
            (cieXYZ[0] + (15.f * cieXYZ[1]) + (3.f * cieXYZ[2]));

        float v = (9.f * cieXYZ[1]) /
            (cieXYZ[0] + (15.f * cieXYZ[1]) + (3.f * cieXYZ[2]));

        float y = cieXYZ[1]/100.f;

        if (y > 0.008856) {
            y = (float)Math.pow(y, 1./3.);
        } else {
            y = (7.787f * y) + (16.f / 116.f);
        }

        float Un = (4.f * Xn) / (Xn + (15.f * Yn) + (3.f * Zn));
        float Vn = (9.f * Yn) / (Xn + (15.f * Yn) + (3.f * Zn));

        float cieL = (116.f * y) - 16.f;
        float cieU = 13.f * cieL * (u - Un);
        float cieV = 13.f * cieL * (v - Vn);

        return new float[] {cieL, cieU, cieV};
    }

    /**
     * convert CIE XYZ (1931) to CIE LAB.
     *
     * uses https://en.wikipedia.org/wiki/Lab_color_space#Forward_transformation
     *
     * L is from 0 to 100
     * a is from - to +  (negative are green, pos are red)
     * b is from - to +  (negative are blue, pos are yellow)
     *
     * to change the reference point, see
     * http://www.easyrgb.com/index.php?X=MATH&H=15#text15
     *
     *
     * @param cieXYZ
     * @param Xn standard illuminant X
     * @param Yn standard illuminant Y
     * @param Zn standard illuminant Z
     * @return
     */
    public float[] cieXYZToCIELAB(float[] cieXYZ, float Xn, float Yn, float Zn) {

        float xDiv = cieXYZ[0]/Xn;
        float yDiv = cieXYZ[1]/Yn;
        float zDiv = cieXYZ[2]/Zn;

        double deltaSq = (6./29.)*(6./29.);
        double deltaCubed = deltaSq * (6./29.);
        float fX = 0;
        float fY = 0;
        float fZ = 0;
        for (int i = 0; i < 3; ++i) {
            double d;
            if (i == 0) {
                d = xDiv;
            } else if (i == 1) {
                d = yDiv;
            } else {
                d = zDiv;
            }
            double f;
            if (d > deltaCubed) {
                f = Math.pow(d, 1./3.);
            } else {
                f = (d/(3.*deltaSq)) + (4./29.);
            }
            if (i == 0) {
                fX = (float)f;
            } else if (i == 1) {
                fY = (float)f;
            } else {
                fZ = (float)f;
            }
        }

        float ell = (116.f * fY) - 16.f;
        float a = 500.f * (fX - fY);
        float b = 200.f * (fY - fZ);

        return new float[]{ell, a, b};
    }

    /**
     * calculate the CIE76 delta E for 2 sets of CIE LAB.
     *
     * uses https://en.wikipedia.org/wiki/Color_difference
     *
     * the "Just noticeable difference", JND, begins at E_ab ~ 2.3
     *
     * @param cieLAB1
     * @param cieLAB2
     * @return deltaE
     */
    public double calcDeltaECIE76(float[] cieLAB1, float[] cieLAB2) {

        double eAB = Math.sqrt(
            ((cieLAB2[0] - cieLAB1[0])*(cieLAB2[0] - cieLAB1[0])) +
            ((cieLAB2[1] - cieLAB1[1])*(cieLAB2[1] - cieLAB1[1])) +
            ((cieLAB2[2] - cieLAB1[2])*(cieLAB2[2] - cieLAB1[2])));

        return eAB;
    }

    /**
     * calculate the CIE76 delta E for 2 sets of CIE LAB.
     *
     * uses https://en.wikipedia.org/wiki/Color_difference
     *
     * the "Just noticeable difference", JND, begins at E_ab ~ 2.3
     *
     * the range of resulting values is 0 through 28.78.
     *
     * @param cieLAB1
     * @param cieLAB2
     * @return deltaE
     */
    public double calcDeltaECIE94(float[] cieLAB1, float[] cieLAB2) {

        return calcDeltaECIE94(cieLAB1[0], cieLAB1[1], cieLAB1[2], cieLAB2[0],
            cieLAB2[1], cieLAB2[2]);

    }

    /**
     * calculate the delta E 1994 for 2 sets of CIE LAB.
     *
     * uses https://en.wikipedia.org/wiki/Color_difference
     *
     * the "Just noticeable difference", JND, begins at E_ab ~ 2.3.
     *
     * the range of resulting values is 0 through 28.78.
     *
     * @return deltaE
     */
    public double calcDeltaECIE94(float ell1, float a1, float b1,
        float ell2, float a2, float b2) {

        // use graphic arts or textiles approx for K's
        boolean useGA = true;
        double kL, K1, K2;
        if (useGA) {
            kL = 1;
            K1 = 0.045;
            K2 = 0.015;
        } else {
            kL = 2;
            K1 = 0.048;
            K2 = 0.014;
        }

        float deltaEll = ell1 - ell2;
        double deltaA = a1 - a2;
        double deltaB = b1 - b2;

        double cA1 = Math.sqrt((a1*a1) + (b1*b1));
        double cA2 = Math.sqrt((a2*a2) + (b2*b2));
        double deltaCab = cA1 - cA2;

        double aa = (deltaA*deltaA) + (deltaB*deltaB) - (deltaCab*deltaCab);
        double deltaHab;
        if (aa < 1E-10) {
            deltaHab = 0;
        } else {
            deltaHab = Math.sqrt(aa);
        }

        double sL = 1;
        double sC = 1 + (K1*cA1);
        double sH = 1 + (K2*cA1);

        double t1 = deltaEll/(kL * sL);
        t1 *= t1;
        double t2 = deltaCab/sC;
        t2 *= t2;
        double t3 = deltaHab/sH;
        t3 *= t3;

        double e94 = Math.sqrt(t1 + t2 + t3);

        return e94;
    }

    /**
     * calculate the delta E 2000 for 2 sets of CIE LAB, specifically, CIEDE2000.
     *
     * uses https://en.wikipedia.org/wiki/Color_difference
     * and
     * http://www.ece.rochester.edu/~gsharma/ciede2000/ciede2000noteCRNA.pdf
     * Sharma, Wu, and Dalal 2004
     *
     * The implementation is adapted from
     * https://github.com/wuchubuzai/OpenIMAJ/blob/master/image/image-processing/src/main/java/org/openimaj/image/analysis/colour/CIEDE2000.java
     * which has copyright 2011, The University of Southampton and allows
     * redistribution of source with or without modification.
     *
     * Note that the "Just noticeable difference", JND, begins at E_ab ~ 2.3.
     *
     * The range of resulting values if CIELUV values are given
     * is 0 to about 130-ish.
     *
     * The range of resulting values if CIELAB1931 values are given
     * is 0 to about 120-ish.
     *
     * @return deltaE
     */
    public double calcDeltaECIE2000(float[] lab1, float[] lab2) {

        return calcDeltaECIE2000(lab1[0], lab1[1], lab1[2], lab2[0], lab2[1],
            lab2[2]);
    }

    /**
     * calculate the delta E 2000 for 2 sets of CIE LAB, specifically, CIEDE2000.
     * see https://en.wikipedia.org/wiki/Color_difference
     * and
     * http://www.ece.rochester.edu/~gsharma/ciede2000/ciede2000noteCRNA.pdf
     * Sharma, Wu, and Dalal 2004
     *
     * The implementation is adapted from
     * https://github.com/wuchubuzai/OpenIMAJ/blob/master/image/image-processing/src/main/java/org/openimaj/image/analysis/colour/CIEDE2000.java
     * which has copyright 2011, The University of Southampton and allows
     * redistribution of source with or without modification.
     *
     * Note that the "Just noticeable difference", JND, begins at E_ab ~ 2.3.
     *
     * The range of resulting values if CIELUV values are given
     * is 0 to about 130-ish.
     *
     * The range of resulting values if CIELAB1931 values are given
     * is 0 to about 120-ish.
     * 
     * @return deltaE
     */
    public double calcDeltaECIE2000(float L1, float a1, float b1,
        float L2, float a2, float b2) {

        double Lmean = (L1 + L2) / 2.0;
		double C1 =  Math.sqrt(a1*a1 + b1*b1);
		double C2 =  Math.sqrt(a2*a2 + b2*b2);
		double Cmean = (C1 + C2) / 2.0;

        double pow7 = Math.pow(Cmean, 7);
        double pow257 =  Math.pow(25, 7);

		double G =  (1. - Math.sqrt(pow7/(pow7 + pow257)))/2.;
		double a1prime = a1 * (1. + G);
		double a2prime = a2 * (1. + G);

		double C1prime =  Math.sqrt(a1prime*a1prime + b1*b1);
		double C2prime =  Math.sqrt(a2prime*a2prime + b2*b2);
		double Cmeanprime = (C1prime + C2prime) / 2;

		double h1prime =  Math.atan2(b1, a1prime) + 2*Math.PI * (Math.atan2(b1, a1prime)<0 ? 1 : 0);
		double h2prime =  Math.atan2(b2, a2prime) + 2*Math.PI * (Math.atan2(b2, a2prime)<0 ? 1 : 0);
		double Hmeanprime =  ((Math.abs(h1prime - h2prime) > Math.PI) ? (h1prime + h2prime + 2*Math.PI) / 2 : (h1prime + h2prime) / 2);

		double T =  1.0 - 0.17 * Math.cos(Hmeanprime - Math.PI/6.0) + 0.24 * Math.cos(2*Hmeanprime)
            + 0.32 * Math.cos(3*Hmeanprime + Math.PI/30) - 0.2 * Math.cos(4*Hmeanprime - 21*Math.PI/60);

		double deltahprime =  ((Math.abs(h1prime - h2prime) <= Math.PI) ? h2prime - h1prime :
            (h2prime <= h1prime) ? h2prime - h1prime + 2*Math.PI : h2prime - h1prime - 2*Math.PI);

		double deltaLprime = L2 - L1;
		double deltaCprime = C2prime - C1prime;
		double deltaHprime =  2.0 * Math.sqrt(C1prime*C2prime) * Math.sin(deltahprime / 2.0);
		double SL =  1.0 + ( (0.015*(Lmean - 50)*(Lmean - 50)) / (Math.sqrt( 20 + (Lmean - 50)*(Lmean - 50) )) );
		double SC =  1.0 + 0.045 * Cmeanprime;
		double SH =  1.0 + 0.015 * Cmeanprime * T;

		double deltaTheta =  (30 * Math.PI / 180) * Math.exp(-((180/Math.PI*Hmeanprime-275)/25)*((180/Math.PI*Hmeanprime-275)/25));
		double RC =  (2 * Math.sqrt(Math.pow(Cmeanprime, 7) / (Math.pow(Cmeanprime, 7) + Math.pow(25, 7))));
		double RT =  (-RC * Math.sin(2 * deltaTheta));

		double KL = 1;
		double KC = 1;
		double KH = 1;

		double deltaE = Math.sqrt(
				((deltaLprime/(KL*SL)) * (deltaLprime/(KL*SL))) +
				((deltaCprime/(KC*SC)) * (deltaCprime/(KC*SC))) +
				((deltaHprime/(KH*SH)) * (deltaHprime/(KH*SH))) +
				(RT * (deltaCprime/(KC*SC)) * (deltaHprime/(KH*SH)))
				);

		return deltaE;
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
     * Note that this returns true for grey too, so test for lower intensities
     * first
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

    private static double wd = Math.sqrt(Math.pow((0.425 - 0.275), 2) +
        Math.pow((0.3875 - 0.2625), 2));

    /**
     * Is within the bounds of large white central region in the CIE 1931 xy
     * chromaticity diagram.  Note that this region has small amount of color
     * in it too, so you may want to further process any point in rgb too.
     * Also note that grey pixels are found here too centered at (0.33, 0.33).
     *
     * @param cieX
     * @param cieY
     * @return
     */
    public boolean isInLargeWhiteCenter(float cieX, float cieY) {

        /*              (0.425, 0.3875)
                           /
                       /           within 0.075 from the line
                    /              between these two points
          (0.275, 0.2625)


        2D point (x,y) and line (a, b, c): dist=(a*x + b*y + c)/sqrt(a^2 + b^2)

        If define (x1, y1) and (x2, y2) as points on the line and
        (x0, y0) as the point that is distant from the line:

        dist = | (x2-x1)(y1-y0) - (x1-x0)(y2-y1)  |
                 --------------------------------
                  sqrt( (x2-x1)^2 + (y2-y1)^2) )
        */
        double numer = Math.abs((0.425 - 0.275) * (0.2625 - cieY)
            - (0.275 - cieX) * (0.3875 - 0.2625));

        double dist = numer/wd;

        return (dist <= 0.03333);//(dist <= 0.075);
    }

    /**
     * calculate the difference in L,A,B 1931 between between the two sets and
     * normalize the values to a sum of "1" using the range of values
     * possible from the use of a standard illuminant, D65.
     * 
     * @return 
     */
    public float calcNormalizedDifferenceLAB31(float ell1, float a1, float b1,
        float ell2, float a2, float b2) {
       
        /*
        * using the standard illuminant of daylight, D65,
        * the range of return values is
        * L*    0 to 104.5
        * a* -190 to 103
        * b* -113 to 99
        */
        
        float diff1 = Math.abs(ell1 - ell2)/104.5f;
        float diff2 = Math.abs(a1 - a2)/(103f + 190f);
        float diff3 = Math.abs(b1 - b2)/(99f + 113.f);
        
        return (diff1 + diff2 + diff3)/3.f;
    }
    
    /**
     * calculate the difference in L, U, V 1976 WideRangeLightness
     * between between the two sets and
     * normalize the values to a sum of "1" using the range of values
     * possible.
     * 
     * * uses standard illuminant "wide range lightness".
     * Wide-range Lightness data were generated by Fairchild et al., 
     * who conducted two different experiments to scale lightness above and 
     * below diffuse white (CIE * L = 100 ). 
     * In the Scaling Lightness Experiment 1 (SL1) they used a luminance
       range from 156 to 2 3692cd m with 2 842 Y cd m n = 
       (Yn represents the luminance of reference white) 
      whereas in the Scaling Lightness Experiment 2 (SL2) the
      luminance range was extended from 0 to 2 7432cd m with 2 997 Y cd m n = . 
     The SL2 data set was used to drive the adapted lightness ( z J ) 
     formula of the proposed color space (see later) while the SL1 data set was 
     used as a test data set. Each of the sets includes 19 samples. 

     * https://www.osapublishing.org/DirectPDFAccess/B810E9AE-C7C9-E594-C72DC7FBE1424F0A_368272/oe-25-13-15131.pdf?da=1&id=368272&seq=0&mobile=no
     * 
     * SL2 Training D65/2° (x,y,z)=968.08 997 883.51 
     * L_a=199 
     * C=0.69. N_C=1, F=1 
     *
     * @return 
     */
    public float calcNormalizedDifferenceLUV_WideRangeLightness(
        int r1, int g1, int b1, int r2, int g2, int b2) {
       
        float[] luv1 = rgbToCIELUV_WideRangeLightness(r1, g1, b1);
        float[] luv2 = rgbToCIELUV_WideRangeLightness(r2, g2, b2);
        
        return calcNormalizedDifferenceLUV_WideRangeLightness(luv1, luv2);
    }
    
    /**
     * calculate the difference in L, U, V 1976 WideRangeLightness
     * between between the two sets and
     * normalize the values to a sum of "1" using the range of values
     * possible.
     * 
     * * uses standard illuminant "wide range lightness".
     * Wide-range Lightness data were generated by Fairchild et al., 
     * who conducted two different experiments to scale lightness above and 
     * below diffuse white (CIE * L = 100 ). 
     * In the Scaling Lightness Experiment 1 (SL1) they used a luminance
       range from 156 to 2 3692cd m with 2 842 Y cd m n = 
       (Yn represents the luminance of reference white) 
      whereas in the Scaling Lightness Experiment 2 (SL2) the
      luminance range was extended from 0 to 2 7432cd m with 2 997 Y cd m n = . 
     The SL2 data set was used to drive the adapted lightness ( z J ) 
     formula of the proposed color space (see later) while the SL1 data set was 
     used as a test data set. Each of the sets includes 19 samples. 

     * https://www.osapublishing.org/DirectPDFAccess/B810E9AE-C7C9-E594-C72DC7FBE1424F0A_368272/oe-25-13-15131.pdf?da=1&id=368272&seq=0&mobile=no
     * 
     * SL2 Training D65/2° (x,y,z)=968.08 997 883.51 
     * L_a=199 
     * C=0.69. N_C=1, F=1 
     *
     * @return 
     */
    public float calcNormalizedDifferenceLUV_WideRangeLightness(
        float[] luv1, int r2, int g2, int b2) {
       
        float[] luv2 = rgbToCIELUV_WideRangeLightness(r2, g2, b2);
         
        return calcNormalizedDifferenceLUV_WideRangeLightness(luv1, luv2);
    }
    
    /**
     * calculate the difference in L, U, V 1976 WideRangeLightness
     * between between the two sets and
     * normalize the values to a sum of "1" using the range of values
     * possible.
     * 
     * * uses standard illuminant "wide range lightness".
     * Wide-range Lightness data were generated by Fairchild et al., 
     * who conducted two different experiments to scale lightness above and 
     * below diffuse white (CIE * L = 100 ). 
     * In the Scaling Lightness Experiment 1 (SL1) they used a luminance
       range from 156 to 2 3692cd m with 2 842 Y cd m n = 
       (Yn represents the luminance of reference white) 
      whereas in the Scaling Lightness Experiment 2 (SL2) the
      luminance range was extended from 0 to 2 7432cd m with 2 997 Y cd m n = . 
     The SL2 data set was used to drive the adapted lightness ( z J ) 
     formula of the proposed color space (see later) while the SL1 data set was 
     used as a test data set. Each of the sets includes 19 samples. 

     * https://www.osapublishing.org/DirectPDFAccess/B810E9AE-C7C9-E594-C72DC7FBE1424F0A_368272/oe-25-13-15131.pdf?da=1&id=368272&seq=0&mobile=no
     * 
     * SL2 Training D65/2° (x,y,z)=968.08 997 883.51 
     * L_a=199 
     * C=0.69. N_C=1, F=1 
     *
     * @return 
     */
    public float calcNormalizedDifferenceLUV_WideRangeLightness(
        float[] luv1, float[] luv2) {
       
        /*
        * L       0 to 104.8
        * u   -99.9 to 176.8
        * v  -150.0 to 95.0
        */        
        float diff1 = Math.abs(luv1[0] - luv2[1])/104.8f;
        float diff2 = Math.abs(luv1[1] - luv2[1])/(176.8f + 99.9f);
        float diff3 = Math.abs(luv1[2] - luv2[2])/(95.0f + 150.0f);
        
        return (diff1 + diff2 + diff3)/3.f;
    }
    
    public float calcNormalizedDifferenceLUV(float ell1, float u1, float v1,
        float ell2, float u2, float v2) {
       
        /*
        * using the standard illuminant of daylight, D65,
        * the range of return values is
        * L       0 to 104.5
        * u   -86.9 to 183.8
        * v  -141.4 to 112.3
        */
        
        float diff1 = Math.abs(ell1 - ell2)/104.5f;
        float diff2 = Math.abs(u1 - u2)/(183.8f + 86.9f);
        float diff3 = Math.abs(v1 - v2)/(112.3f + 141.4f);
        
        return (diff1 + diff2 + diff3)/3.f;
    }
    
}
