package algorithms.imageProcessing;

import algorithms.imageProcessing.transform.TransformationParameters;
import algorithms.imageProcessing.transform.Transformer;
import algorithms.misc.MiscDebug;
import algorithms.util.ResourceFinder;
import java.util.ArrayList;
import java.util.List;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class B3SplineFunctionTest extends TestCase {
    
    public B3SplineFunctionTest() {
    }

    public void testInterp1D() {
        
        int d = 10;
        
        GreyscaleImage input = new GreyscaleImage(d, d);
        for (int i = 0; i < input.getNPixels(); ++i) {
            input.setValue(i, d);
        }
        
        B3SplineFunction instance = new B3SplineFunction();
        
        int x = 5;
        int y = 5;
        int interp = instance.interpolate1DX(x, y, input);
        assertEquals(10, interp);
        
        interp = instance.interpolate1DY(x, y, input);
        assertEquals(10, interp);
        
    }
    
    public void testInterp2D() {
        
        int d = 10;
        
        GreyscaleImage input = new GreyscaleImage(d, d);
        for (int i = 0; i < input.getNPixels(); ++i) {
            input.setValue(i, 1);
        }
        
        B3SplineFunction instance = new B3SplineFunction();
        
        int x = 5;
        int y = 5;
        int interp = instance.interpolate2D(x, y, input);
        assertEquals(1, interp);
        
        input = new GreyscaleImage(d, d);
        for (int i = 0; i < input.getNPixels(); ++i) {
            input.setValue(i, 10);
        }
        
        interp = instance.interpolate2D(x, y, input);
        assertEquals(10, interp);
    }
    
    public void testCalculate() {
        
        int d = 10;
        
        int n = d*d;
        
        double[] input = new double[n];
        
        for (int x = 0; x < d; ++x) {
            for (int y = 0; y < d; ++y) {
                int pixIdx = (y * d) + x;
                input[pixIdx] = y*256;
            }
        }
        
        /*
        input:
        0  0  0  0  0  0  0  0  0  0 
        1  1  1  1  1  1  1  1  1  1
        2  2  2  2  2  2  2  2  2  2
        3  3  3  3 *3* 3  3  3  3  3
        4  4  4  4  4  4  4  4  4  4
        5  5  5  5  5  5  5  5  5  5
        ...
        
        window:    
         (1/256) *   |   1    4   3*2    4    1 |
                     |   4   16   3*8   16    4 |
                     | 3*2  3*8   9*4  3*8  3*2 |
                     |   4   16   3*8   16    4 |
                     |   1    4   3*2    4    1 |

        for x=4, y=3  window before mult is:
           1  1  1  1  1
           2  2  2  2  2
           3  3 *3* 3  3
           4  4  4  4  4
           5  5  5  5  5        
        */
        
        int expectedSum = (1 + 4 + (3*2) + 4 + 1) + 2*(4 + 16 + (3*8) + 16 + 4)
            + 3*(3*2 + 3*8 + 9*4 + 3*8 + 3*2) + 4*(4 + 16 + (3*8) + 16 + 4)
            + 5*(1 + 4 + (3*2) + 4 + 1);
        
        B3SplineFunction instance = new B3SplineFunction();
        
        double[] output = instance.calculate(input, d, d);
        
        int pixIdx = (3 * d) + 4;
        
        assertTrue(Math.abs(output[pixIdx] - expectedSum) < 0.0000001);
    }
    
    public void test1D2D() throws Exception {
        
        // demonstrates that the two 1-D operations produce similar, but not
        // identical results as the one 2-D use of spline b3 within the transform
        String fileName1 = "campus_010.jpg";
        String fileName2 = "campus_011.jpg";
        int idx = fileName1.lastIndexOf(".");
        String fileName1Root = fileName1.substring(0, idx);
        idx = fileName2.lastIndexOf(".");
        String fileName2Root = fileName2.substring(0, idx);
                
        String filePath1 = ResourceFinder.findFileInTestResources(fileName1);
        ImageExt img1 = ImageIOHelper.readImageExt(filePath1);
        String filePath2 = ResourceFinder.findFileInTestResources(fileName2);
        ImageExt img2 = ImageIOHelper.readImageExt(filePath2);
        
        ATrousWaveletTransform awt = new ATrousWaveletTransform();
         
        List<GreyscaleImage> outputTransformed = new ArrayList<GreyscaleImage>();
        List<GreyscaleImage> outputCoeff = new ArrayList<GreyscaleImage>();
        awt.calculateWithB3SplineScalingFunction(img1.copyToGreyscale(), 
            outputTransformed, outputCoeff);
        
        List<GreyscaleImage> outputTransformed2 = new ArrayList<GreyscaleImage>();
        List<GreyscaleImage> outputCoeff2 = new ArrayList<GreyscaleImage>();
        awt.calculateWithB3SplineScalingFunction2(img1.copyToGreyscale(), 
            outputTransformed2, outputCoeff2);
        
        GreyscaleImage c1F = outputCoeff.get(outputCoeff.size() - 1);
        
        GreyscaleImage c2F = outputCoeff2.get(outputCoeff2.size() - 1);
        
        boolean plot = false;
        for (int row = 0; row < c1F.getHeight(); ++row) {
            for (int col = 0; col < c1F.getWidth(); ++col) {
                int v1 = c1F.getValue(col, row);
                int v2 = c2F.getValue(col, row);
                int diffV = v1 - v2;
                if (diffV != 0) {
                    plot = true;
                }
            }
        }
        if (true) {
            MiscDebug.writeImage(c1F, "2-1D_atrous");
            MiscDebug.writeImage(c2F, "1-2D_atrous");
        }
        
        // assert that using the 2-D function on a rotated and not rotated
        // image produces same results
        
        img1 = ImageIOHelper.readImageExt(filePath1);
        TransformationParameters params90 = new TransformationParameters();
        params90.setRotationInDegrees(90);
        params90.setOriginX(0);
        params90.setOriginY(0);
        params90.setTranslationX(0);
        params90.setTranslationY(img1.getWidth() - 1);
        Transformer transformer = new Transformer();
        ImageExt img1R90 = (ImageExt) transformer.applyTransformation(img1,
            params90, img1.getHeight(), img1.getWidth());
                
        outputTransformed = new ArrayList<GreyscaleImage>();
        outputCoeff = new ArrayList<GreyscaleImage>();
        awt.calculateWithB3SplineScalingFunction2(img1R90.copyToGreyscale(), 
            outputTransformed, outputCoeff);
        
        img1 = ImageIOHelper.readImageExt(filePath1);
        outputTransformed2 = new ArrayList<GreyscaleImage>();
        outputCoeff2 = new ArrayList<GreyscaleImage>();
        awt.calculateWithB3SplineScalingFunction2(img1.copyToGreyscale(), 
            outputTransformed2, outputCoeff2);
        
        c1F = outputCoeff.get(outputCoeff.size() - 1);
        
        c2F = outputCoeff2.get(outputCoeff2.size() - 1);
        
        assertEquals(outputCoeff.size(), outputCoeff2.size());
      
    }
}
