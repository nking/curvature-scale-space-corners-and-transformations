package algorithms;

import algorithms.imageProcessing.CIEChromaticity;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class TmpTest extends TestCase {
    
    public void test0() {
        CIEChromaticity cieC = new CIEChromaticity();
        
        float[] minMaxL = new float[]{Float.MAX_VALUE, Float.NEGATIVE_INFINITY};
        float[] minMaxU = new float[]{Float.MAX_VALUE, Float.NEGATIVE_INFINITY};
        float[] minMaxV = new float[]{Float.MAX_VALUE, Float.NEGATIVE_INFINITY};
        
        float[] a = null;
        for (int r = 0; r <= 255; r += 5) {
            for (int g = 0; g <= 255; g += 5) {
                for (int b = 0; b <= 255; b += 5) {
                    a = cieC.rgbToCIELUV_WideRangeLightness(r, g, b);
                    if (a[0] < minMaxL[0]) {
                        minMaxL[0] = a[0];
                    }
                    if (a[0] > minMaxL[1]) {
                        minMaxL[1] = a[0];
                    }
                    if (a[1] < minMaxU[0]) {
                        minMaxU[0] = a[1];
                    }
                    if (a[1] > minMaxU[1]) {
                        minMaxU[1] = a[1];
                    }
                    if (a[2] < minMaxV[0]) {
                        minMaxV[0] = a[2];
                    }
                    if (a[2] > minMaxV[1]) {
                        minMaxV[1] = a[2];
                    }
                }
            }
        }
        System.out.format(
            "minMax L, U, V = (%.3f, %.3f) (%.3f, %.3f) (%.3f, %.3f)", 
            minMaxL[0], minMaxL[1], minMaxU[0], minMaxU[1],
            minMaxV[0], minMaxV[1]);
    }
}
