package algorithms.imageProcessing.features;

import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.HistogramEqualization;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.MorphologicalFilter;
import algorithms.imageProcessing.NonMaximumSuppression;
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
        
        //String fileName = "blox.gif";
        String fileName = "lab.gif";
        //String fileName = "house.gif";
        //String fileName = "susan-in_plus.png";
        //String fileName = "tmp.png";
        //String fileName = "lena.jpg";
        String filePath = ResourceFinder.findFileInTestResources(fileName);
        
        GreyscaleImage img = ImageIOHelper.readImageAsGrayScale(filePath).copyToGreyscale();
        
        float cutoff = 0.5f;//0.5f;
        PhaseCongruencyDetector phaseCDetector = new PhaseCongruencyDetector();
        PhaseCongruencyDetector.PhaseCongruencyProducts products =
            //phaseCDetector.phaseCongMono(img);
            phaseCDetector.phaseCongMono(img, 5, 3, 2.1f, 0.55f, 2, cutoff, 
                10, 1.5f, -1);
        assertNotNull(products);
                
        // range of data is 0 to 
        double[][] pc = products.getPhaseCongruency();
        double[][] orientation = products.getOrientation();
        double[][] phaseAngle = products.getPhaseAngle();
        
        NonMaximumSuppression ns = new NonMaximumSuppression();
        double[][] imgTh0 = ns.nonmaxsup(products.getPhaseCongruency(), 
            products.getOrientation(), 1.5);
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        double[][] thresholded = imageProcessor.threshold(imgTh0, 0.1);
        
        
        double t1 = 0.1;
        double t2 = 0.3;
        phaseCDetector.explorePhaseAngle(products, imgTh0, t1);
        
        GreyscaleImage out = img.createWithDimensions();
        for (int i = 0; i < out.getWidth(); ++i) {
            for (int j = 0; j < out.getHeight(); ++j) {
                out.setValue(i, j, (int)(255.*pc[j][i]));
            }
        }
        MiscDebug.writeImage(out, "_phase_congruency_");
        
        out = img.createWithDimensions();
        for (int i = 0; i < out.getWidth(); ++i) {
            for (int j = 0; j < out.getHeight(); ++j) {
                if (thresholded[j][i] > 0) {
                    out.setValue(i, j, 255);
                }
            }
        }
        MiscDebug.writeImage(out, "_phase_congruency__thresholded");
        
        out = img.createWithDimensions();
        double maxV = Double.MIN_VALUE;
        double minV = Double.MAX_VALUE;
        double[] scaled = new double[out.getNPixels()];
        for (int i = 0; i < out.getWidth(); ++i) {
            for (int j = 0; j < out.getHeight(); ++j) {
                double v = orientation[j][i];
                if (v < minV) {
                    minV = v;
                }
                if (v > maxV) {
                    maxV = v;
                }
                int idx = out.getInternalIndex(i, j);
                scaled[idx] = v;
            }
        }
        System.out.println("orientation minV=" + minV + " maxV=" + maxV);
        int[] scaled2 = MiscMath.rescale(scaled, 0, 255);
        for (int i = 0; i < out.getNPixels(); ++i) {
            int v = scaled2[i];
            out.setValue(i, v);                
        }
        MiscDebug.writeImage(out, "_orientation_");
        
        out = img.createWithDimensions();
        maxV = Double.MIN_VALUE;
        minV = Double.MAX_VALUE;
        scaled = new double[out.getNPixels()];
        for (int i = 0; i < out.getWidth(); ++i) {
            for (int j = 0; j < out.getHeight(); ++j) {
                double v = phaseAngle[j][i];
                if (v < minV) {
                    minV = v;
                }
                if (v > maxV) {
                    maxV = v;
                }
                int idx = out.getInternalIndex(i, j);
                scaled[idx] = v;
            }
        }
        System.out.println("phase angle minV=" + minV + " maxV=" + maxV);
        scaled2 = MiscMath.rescale(scaled, 0, 255);
        for (int i = 0; i < out.getNPixels(); ++i) {
            int v = scaled2[i];
            out.setValue(i, v);                
        }
        MiscDebug.writeImage(out, "_phase_angle_");
        
        System.out.println("threshold=" + products.getThreshold());
        
        /*
        [PC, or] = phasecongmono(imread('lena.tif'));
        nm = nonmaxsup(PC, or, 1.5);   % nonmaxima suppression
        bw = hysthresh(nm, 0.1, 0.3);  % hysteresis thresholding 0.1 - 0.3
        */
        
        
        out = img.createWithDimensions();
        for (int i = 0; i < out.getWidth(); ++i) {
            for (int j = 0; j < out.getHeight(); ++j) {
                out.setValue(i, j, (int)(255.*imgTh0[j][i]));
            }
        }
        MiscDebug.writeImage(out, "_non_max_sup_");
                        
        int[][] binaryImage = phaseCDetector.applyHysThresh(imgTh0, 0.1, 0.3);
                
        out = img.createWithDimensions();
        for (int i = 0; i < out.getWidth(); ++i) {
            for (int j = 0; j < out.getHeight(); ++j) {
                out.setValue(i, j, binaryImage[j][i]);
            }
        }
        MiscDebug.writeImage(out, "_filtered_hysteresis_");        
    }
    
}
