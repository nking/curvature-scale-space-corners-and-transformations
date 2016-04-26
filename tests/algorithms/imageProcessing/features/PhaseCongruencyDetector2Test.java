package algorithms.imageProcessing.features;

import algorithms.imageProcessing.CannyEdgeFilterLite;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.ImageSegmentation;
import algorithms.imageProcessing.MedianTransform;
import algorithms.imageProcessing.SIGMA;
import algorithms.imageProcessing.SegmentationMergeThreshold;
import algorithms.misc.MiscDebug;
import algorithms.misc.MiscMath;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.ResourceFinder;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
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
           // "susan-in_plus.png", "lena.jpg",
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
            
        for (String fileName : fileNames) {
            
            System.out.println("fileName=" + fileName);
        
            String filePath = ResourceFinder.findFileInTestResources(fileName);
        
            GreyscaleImage img = ImageIOHelper.readImageAsGrayScale(filePath).copyToGreyscale();
        
            PhaseCongruencyDetector2 phaseCDetector = new PhaseCongruencyDetector2();
            PhaseCongruencyDetector2.PhaseCongruencyProducts products =
                phaseCDetector.phaseCong(img, nScale, nOrient, minWavelength, mult, 
                sigmaOnf, k, increaseKIfNeeded,
                cutOff, g, noiseMethod, tLow, tHigh, doStoreConvolution);

            assertNotNull(products);
            
        }
    }
    
}
