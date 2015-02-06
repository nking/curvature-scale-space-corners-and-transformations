package algorithms.imageProcessing;

import algorithms.util.ArrayPair;

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
        greenish yellow, yellow, and yellowish orange
        455, 355
        559, 459
        493, 480
        467, 455
        437, 446
        391, 457
        345, 474
        
        including center is larger polygon w/ fewer points:
        302, 221
        559, 459
        331, 573
        */
        
        ArrayPair p = new ArrayPair(
            new float[]{.302f, .559f, .331f, .298f, .302f},
            new float[]{.221f, .459f, .573f, .449f, .221f}
        );
        
        return p;
    }
    
    public ArrayPair getYellowishGreenThroughOrangePolynomial() {
        /*
        422, 558
        550, 551
        687, 585
        457, 357
        393, 453
        438, 445
        467, 455
        490, 479
        495, 511
        438, 524
        
        including center is larger polygon w/ fewer points:
        302, 221
        615, 510
        330, 573        
        */
        
        ArrayPair p = new ArrayPair(
            new float[]{.302f, .615f, .330f, .302f},
            new float[]{.221f, .510f, .573f, .221f}
        );
        
        return p;
    }
    
    public ArrayPair getRedPolynomial() {
        
        /*
        529, 606
        660, 670
        740, 636
        688, 588
        553, 551
        
        including center is larger polygon w/ fewer points:
        330, 573
        529, 606
        660, 670
        740, 636
        688, 588
        552, 552
        */
        
        ArrayPair p = new ArrayPair(
            new float[]{.330f, .529f, .660f, .740f, .688f, .552f, .330f},
            new float[]{.573f, .606f, .670f, .636f, .588f, .552f, .573f}
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
}
