package algorithms.imageProcessing.features;

import algorithms.imageProcessing.CannyEdgeFilterAdaptive;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.misc.MiscDebug;
import algorithms.misc.MiscMath;
import algorithms.util.ResourceFinder;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class PhaseCongruencyDetectorTest extends TestCase {
    
    public PhaseCongruencyDetectorTest() {
    }
    
    public void test0() throws Exception {

        String[] fileNames = new String[]{
            "blox.gif", "lab.gif", "house.gif", "seattle.jpg", "merton_college_I_001.jpg",
            "susan-in_plus.png", "lena.jpg",
            "campus_010.jpg"
        };

        float cutOff = 0.5f;//0.3f;//0.5f;
        int nScale = 5;
        int minWavelength = 3;
        float mult = 2.1f;
        float sigmaOnf = 0.55f;
        int k = 2;
        float g = 10; 
        float deviationGain = 1.5f;
        int noiseMethod = -1;
        double tLow = 0.05;
        double tHigh = 0.3;
        boolean increaseKIfNeeded = true;

        for (String fileName : fileNames) {
        
            String filePath = ResourceFinder.findFileInTestResources(fileName);
        
            GreyscaleImage img = ImageIOHelper.readImageAsGrayScale(filePath).copyToGreyscale();
        
            PhaseCongruencyDetector phaseCDetector = new PhaseCongruencyDetector();
            phaseCDetector.setToCreateCorners();                
            PhaseCongruencyDetector.PhaseCongruencyProducts products =
                phaseCDetector.phaseCongMono(img, nScale, minWavelength, mult, 
                sigmaOnf, k, increaseKIfNeeded,
                cutOff, g, deviationGain, noiseMethod, tLow, tHigh);

            assertNotNull(products);
            double[][] pc = products.getPhaseCongruency();
            double[][] orientation = products.getOrientation();
            double[][] phaseAngle = products.getPhaseAngle();
            int[][] thinned = products.getThinned();
                
            GreyscaleImage out = img.createWithDimensions();
            for (int i = 0; i < out.getWidth(); ++i) {
                for (int j = 0; j < out.getHeight(); ++j) {
                    out.setValue(i, j, thinned[j][i]);
                }
            }
            MiscDebug.writeImage(out, "_thinned_" + cutOff + "_" + fileName + "_");  
        }
    }
    
}
