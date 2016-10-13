package algorithms.imageProcessing;

import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class OtsuThresholding2Test extends TestCase {
    
    public void test2D() {
        
        /*
        3  1  9  10
        3  2  8  10
        3  3  3  10 
        3  3  3  10
        */
        double[][] img = new double[4][];
        img[0] = new double[] {3,  1,  9,  10};
        img[1] = new double[] {3,  2,  8,  10};
        img[2] = new double[] {3,  3,  3,  10};
        img[3] = new double[] {3,  3,  3,  10};
        
        OtsuThresholding thrshFinder = new OtsuThresholding();
            
        double thrsh0 = thrshFinder.calculateBinaryThreshold2D(img,
                51);
           
        assert(thrsh0 >= 3 && thrsh0 < 10);
        
    }
}
