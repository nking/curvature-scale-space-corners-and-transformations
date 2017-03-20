package algorithms.imageProcessing;

import algorithms.misc.MiscDebug;
import algorithms.util.PairInt;
import algorithms.util.ResourceFinder;
import java.util.HashSet;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class NonMaximumSuppressionTest extends TestCase {
    
    public NonMaximumSuppressionTest() {
    }

    public void test0() throws Exception {
        
        //TODO: write a better test for this...
        
        String fileNamePC = "small_shapes_pc.png";
        String fileNameTheta = "small_shapes_theta.png";

        String filePath = ResourceFinder.findFileInTestResources(fileNamePC);
        
        GreyscaleImage pcImg = ImageIOHelper.readImage(filePath).copyToGreyscale(); 
        filePath = ResourceFinder.findFileInTestResources(fileNameTheta);
        GreyscaleImage thetaImg = 
            ImageIOHelper.readImage(filePath).copyToGreyscale(); 
        
        // need to rescale the theta images which should have max 179 but is
        // read as max 255.
        float f = 179.f/255.f;
        for (int i = 0; i < thetaImg.getWidth(); ++i) {
            for (int j = 0; j < thetaImg.getHeight(); ++j) {
                int v = Math.round(f * thetaImg.getValue(i, j));
                thetaImg.setValue(i, j, v);
            }
        }
        
        double radius = 1.2;
        NonMaximumSuppression nms = new NonMaximumSuppression();
        nms.nonmaxsup(pcImg, thetaImg, radius, new HashSet<PairInt>());
        
        MiscDebug.writeImage(thetaImg, "_THETA_out");
        MiscDebug.writeImage(pcImg, "_PC_out");
       
        for (int i = 17; i <= 31; ++i) {
            for (int j = 88; j <= 88; ++j) {
                assertTrue(pcImg.getValue(i, j) > 0);
            }
            for (int j = 79; j <= 79; ++j) {
                assertTrue(pcImg.getValue(i, j) > 0);
            }
        }
        
        for (int j = 81; j <= 76; ++j) {
            for (int i = 15; i <= 15; ++i) {
                assertTrue(pcImg.getValue(i, j) > 0);
            }
            for (int i = 34; i <= 34; ++i) {
                assertTrue(pcImg.getValue(i, j) > 0);
            }
        }
        
        for (int j = 52; j <= 78; ++j) {
            for (int i = 59; i <= 59; ++i) {
                assertTrue(pcImg.getValue(i, j) > 0);
            }
        }
        
        // spot checks:
        assertTrue(pcImg.getValue(16, 15) > 0);
        
        assertTrue(pcImg.getValue(21, 20) > 0);
        
        assertTrue(pcImg.getValue(11, 18) > 0);
                
        assertTrue(pcImg.getValue(6, 13) > 0);
        
        assertTrue(pcImg.getValue(48, 25) > 0);                 
    }
    
}
