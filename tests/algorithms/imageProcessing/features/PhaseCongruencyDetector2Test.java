package algorithms.imageProcessing.features;

import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.util.ResourceFinder;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class PhaseCongruencyDetector2Test extends TestCase {
    
    public PhaseCongruencyDetector2Test() {
    }
    
    public void test0() throws Exception {

        String[] fileNames = new String[]{
           // "blox.gif", "lab.gif", "house.gif", "seattle.jpg", 
           //"merton_college_I_001.jpg",
           "lena.jpg",
           // "susan-in_plus.png", 
           // "campus_010.jpg", 
            //"android_statues_01.jpg", 
            //"android_statues_02.jpg", "android_statues_03.jpg", "android_statues_04.jpg"
        };

        float cutOff = 0.5f;//0.3f;//0.5f;
        int nScale = 5;
        int minWavelength = 3;//nScale;// 3;
        float mult = 2.1f;
        float sigmaOnf = 0.55f;
        int k = 10;//5;//2;
        float g = 10; 
        int nOrient = 6;
        boolean doStoreConvolution = false;
        int noiseMethod = -1;
        double tLow = 0.0001;
        double tHigh = 0.1;
        boolean increaseKIfNeeded = false;

        ImageProcessor imageProcessor = new ImageProcessor();
        
        String label = "label=n" + nScale + "_mw" + minWavelength + "_k" + k 
            + "_t0_" + tLow + "_t1_" + tHigh;
            
        float maxDimension = 512;
        
        for (String fileName : fileNames) {
            
            System.out.println("fileName=" + fileName);
        
            String filePath = ResourceFinder.findFileInTestResources(fileName);
        
            GreyscaleImage img = ImageIOHelper.readImageAsGrayScale(filePath).copyToGreyscale();
                    
            int w1 = img.getWidth();
            int h1 = img.getHeight();

            int binFactor1 = (int) Math.ceil(Math.max(
                (float) w1 / maxDimension,
                (float) h1 / maxDimension));

            img = imageProcessor.binImage(img, binFactor1);
            
            PhaseCongruencyDetector2 phaseCDetector = new PhaseCongruencyDetector2();
            phaseCDetector.setToDebug();
            PhaseCongruencyDetector2.PhaseCongruencyProducts products =
                phaseCDetector.phaseCong(img, nScale, nOrient, minWavelength, 
                mult, sigmaOnf, k, cutOff, g, noiseMethod, tLow, tHigh, 
                doStoreConvolution);

            assertNotNull(products);    
        }
    }
    
}
