package algorithms.imageProcessing;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class Kernel1DHelperTest {
    
    public Kernel1DHelperTest() {
    }
    
    @Before
    public void setUp() {
    }
    
    @After
    public void tearDown() {
    }

    @Test
    public void testConvolvePointWithKernel_5args() throws Exception {
        
        // create a delta dirac function and convolve it w/ first derivative
        
        GreyscaleImage img = GaussianHelperForTests.getDeltaDiracImage();
        //GreyscaleImage img = GaussianHelperForTests.getCircle();
                
        Kernel1DHelper helper = new Kernel1DHelper();
        
        float[] kernel = Gaussian1DFirstDeriv.getBinomialKernel(
            SIGMA.ZEROPOINTFIVE);
                
        boolean calculateForX = true;
        
        GreyscaleImage outputX = new GreyscaleImage(img.getWidth(), 
            img.getHeight());
        
        GreyscaleImage outputY = outputX.copyImage();
        
        for (int col = 0; col < img.getWidth(); col++) {
            for (int row = 0; row < img.getHeight(); row++) {
                
                double conv = helper.convolvePointWithKernel(img, col, row, 
                    kernel, calculateForX);
                
                outputX.setValue(col, row, (int)conv);
                
                conv = helper.convolvePointWithKernel(img, col, row, 
                    kernel, false);
                
                outputY.setValue(col, row, (int)conv);
            }
        }
        
        /*String dirPath = ResourceFinder.findDirectory("bin");
        ImageIOHelper.writeOutputImage(dirPath + "/tmpGX.png", outputX);
        ImageIOHelper.writeOutputImage(dirPath + "/tmpGY.png", outputY);
        */
        
        GreyscaleImage output2 = new GreyscaleImage(img.getWidth(), img.getHeight());
        
        calculateForX = false;
        
        for (int col = 0; col < img.getWidth(); col++) {
            for (int row = 0; row < img.getHeight(); row++) {
                
                double conv = helper.convolvePointWithKernel(outputX, col, row, 
                    kernel, calculateForX);
                
                output2.setValue(col, row, (int)conv);
            }
        }
        
        //ImageIOHelper.writeOutputImage(dirPath + "/tmpG2DDeriv.png", output2);
        int xc = output2.getWidth() >> 1;
        int yc = xc;
        
        assertTrue(output2.getValue(xc, yc) == 0);
        int v = output2.getValue(xc - 1, yc - 1);
        assertTrue(v > 0);
        assertTrue(output2.getValue(xc - 1, yc + 1) == -1*v);
        assertTrue(output2.getValue(xc + 1, yc + 1) == v);
        assertTrue(output2.getValue(xc + 1, yc - 1) == -1*v);
    }
    
}
