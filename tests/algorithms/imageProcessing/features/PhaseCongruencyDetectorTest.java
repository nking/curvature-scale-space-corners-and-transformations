package algorithms.imageProcessing.features;

import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.ImageSegmentation;
import algorithms.misc.MiscDebug;
import algorithms.util.ResourceFinder;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class PhaseCongruencyDetectorTest extends TestCase {
    
    public PhaseCongruencyDetectorTest() {
    }
    
    public void est0() throws Exception {

        String[] fileNames = new String[]{
           // "blox.gif", "lab.gif", "house.gif", "seattle.jpg", "merton_college_I_001.jpg",
           // "susan-in_plus.png", "lena.jpg",
           // "campus_010.jpg", 
            "android_statues_01.jpg", 
            "android_statues_02.jpg", "android_statues_03.jpg", "android_statues_04.jpg"
        };

        float cutOff = 0.5f;//0.3f;//0.5f;
        int nScale = 5;
        int minWavelength = 3;//nScale;// 3;
        float mult = 2.1f;
        float sigmaOnf = 0.55f;
        int k = 5;//2;
        float g = 10; 
        float deviationGain = 1.5f;
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
        
            PhaseCongruencyDetector phaseCDetector = new PhaseCongruencyDetector();
            phaseCDetector.setToCreateCorners();                
            PhaseCongruencyDetector.PhaseCongruencyProducts products =
                phaseCDetector.phaseCongMono(img, nScale, minWavelength, mult, 
                sigmaOnf, k, increaseKIfNeeded,
                cutOff, g, deviationGain, noiseMethod, tLow, tHigh);

            assertNotNull(products);
            int[][] thinned = products.getThinned();
                
            GreyscaleImage pcImg = img.createWithDimensions();
            GreyscaleImage out2 = img.createWithDimensions();
            GreyscaleImage out = img.createWithDimensions();
            for (int i = 0; i < out.getWidth(); ++i) {
                for (int j = 0; j < out.getHeight(); ++j) {
                    int vPC = (int)Math.round(255. * products.getPhaseCongruency()[j][i]);
                    if (thinned[j][i] > 0) {
                        out.setValue(i, j, thinned[j][i]);
                        out2.setValue(i, j, vPC);
                    }
                    pcImg.setValue(i, j, vPC);
                }
            }
            MiscDebug.writeImage(out, "_thinned_" + label + "_" + fileName + "_"); 
            MiscDebug.writeImage(out2, "_pc_thinned_" + label + "_" + fileName + "_");
            MiscDebug.writeImage(pcImg, "_pc_" + label + "_" + fileName + "_");
            
            // ----- make O1 edges
            ImageExt imgClr = ImageIOHelper.readImageExt(filePath);
            
            GreyscaleImage o1 = imageProcessor.createO1(imgClr);
            
            products =
                phaseCDetector.phaseCongMono(o1, nScale, minWavelength, mult, 
                sigmaOnf, k, increaseKIfNeeded,
                cutOff, g, deviationGain, noiseMethod, tLow, tHigh);
            thinned = products.getThinned();
            out = img.createWithDimensions();
            pcImg = img.createWithDimensions();
            out2 = img.createWithDimensions();
            for (int i = 0; i < out.getWidth(); ++i) {
                for (int j = 0; j < out.getHeight(); ++j) {
                    int vPC = (int)Math.round(255. * products.getPhaseCongruency()[j][i]);
                    if (thinned[j][i] > 0) {
                        out.setValue(i, j, thinned[j][i]);
                        out2.setValue(i, j, vPC);
                    }
                    pcImg.setValue(i, j, vPC);
                }
            }
        
            MiscDebug.writeImage(out, "_thinned_o1_" + label + "_" + fileName);
            
            MiscDebug.writeImage(out2, "_pc_thinned_o1_" + label + "_" + fileName);
            
            MiscDebug.writeImage(pcImg, "_pc_o1_" + label + "_" + fileName);
        }
    }
    
    public void test1() throws Exception {

        String[] fileNames = new String[]{
            //"blox.gif", "lab.gif", "house.gif", "seattle.jpg", "merton_college_I_001.jpg",
            // "susan-in_plus.png", "lena.jpg",
            // "campus_010.jpg", 
            "android_statues_01.jpg", 
            "android_statues_02.jpg", 
            "android_statues_03.jpg", "android_statues_04.jpg"
        };
        
        ImageSegmentation imageSegmentation = new ImageSegmentation();
                     
        for (String fileName : fileNames) {
            
            System.out.println("fileName=" + fileName);
        
            String filePath = ResourceFinder.findFileInTestResources(fileName);
        
            ImageExt img = ImageIOHelper.readImageExt(filePath);            
            
            GreyscaleImage edgeImage = imageSegmentation.createColorEdges(img);

            /*
            ImageProcessor imageProcessor = new ImageProcessor();
            GreyscaleImage labAImg = imageProcessor.createLabAandB(img)[1];
            GreyscaleImage gsImg = img.copyToGreyscale();
            GreyscaleImage o1Img = imageProcessor.createO1(img);
            MiscDebug.writeImage(labAImg, "_laba_A_");
            MiscDebug.writeImage(gsImg, "_grey_");
            MiscDebug.writeImage(o1Img, "_o1_");
            */
        }
    }
}
