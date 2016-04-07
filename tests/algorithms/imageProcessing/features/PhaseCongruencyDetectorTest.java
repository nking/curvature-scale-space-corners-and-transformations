package algorithms.imageProcessing.features;

import algorithms.compGeometry.RotatedOffsets;
import algorithms.imageProcessing.CannyEdgeFilterAdaptive;
import algorithms.imageProcessing.Gaussian1DFirstDeriv;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.HistogramEqualization;
import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.Kernel1DHelper;
import algorithms.imageProcessing.MorphologicalFilter;
import algorithms.imageProcessing.NonMaximumSuppression;
import algorithms.imageProcessing.SIGMA;
import algorithms.misc.Histogram;
import algorithms.misc.HistogramHolder;
import algorithms.misc.MiscDebug;
import algorithms.misc.MiscMath;
import algorithms.util.Errors;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.ResourceFinder;
import java.util.Arrays;
import java.util.List;
import junit.framework.TestCase;
import org.ejml.simple.SimpleMatrix;

/**
 *
 * @author nichole
 */
public class PhaseCongruencyDetectorTest extends TestCase {
    
    public PhaseCongruencyDetectorTest() {
    }
    
    public void test0() throws Exception {
        
        //String fileName = "blox.gif";
        //String fileName = "lab.gif";
        String fileName = "house.gif";
        //String fileName = "susan-in_plus.png";
        //String fileName = "lena.jpg";
        String filePath = ResourceFinder.findFileInTestResources(fileName);
        
        GreyscaleImage img = ImageIOHelper.readImageAsGrayScale(filePath).copyToGreyscale();
        
        float cutOff = 0.5f;//0.3f;//0.5f;
        int nScale = 5;
        int minWavelength = 3;
        float mult = 2.1f;
        float sigmaOnf = 0.55f;
        float k = 10;//2;
        float g = 10; 
        float deviationGain = 1.5f;
        int noiseMethod = -1;
            
        PhaseCongruencyDetector phaseCDetector = new PhaseCongruencyDetector();
        
        phaseCDetector.setToCreateCorners();
        
        PhaseCongruencyDetector.PhaseCongruencyProducts products =
            phaseCDetector.phaseCongMono(img, nScale, minWavelength, mult, 
                sigmaOnf, k, cutOff, g, deviationGain, noiseMethod);
        
        assertNotNull(products);
                
        double[][] pc = products.getPhaseCongruency();
        double[][] orientation = products.getOrientation();
        double[][] phaseAngle = products.getPhaseAngle();
        
        int[][] thinned = products.getThinned();
                
        GreyscaleImage out = img.createWithDimensions();
        for (int i = 0; i < out.getWidth(); ++i) {
            for (int j = 0; j < out.getHeight(); ++j) {
                out.setValue(i, j, (int)(255.*pc[j][i]));
            }
        }
        MiscDebug.writeImage(out, "_phase_congruency_");
        
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
                
        out = img.createWithDimensions();
        for (int i = 0; i < out.getWidth(); ++i) {
            for (int j = 0; j < out.getHeight(); ++j) {
                out.setValue(i, j, thinned[j][i]);
            }
        }
        MiscDebug.writeImage(out, "_filtered_hysteresis_" + cutOff + "_");  
        
        Image out2 = img.copyToColorGreyscale();
        ImageIOHelper.addCurveToImage(products.getCorners(), out2, 
            2, 255, 0, 0);
        MiscDebug.writeImage(out2, "_corners_");
                
        // ---- compare to CannyEdgeFilter ----
        img = ImageIOHelper.readImageAsGrayScale(filePath).copyToGreyscale();
        
        CannyEdgeFilterAdaptive canny1 = new CannyEdgeFilterAdaptive();
        canny1.applyFilter(img);
        MiscDebug.writeImage(img, "_gradient_canny_adaptivw_");
        
        img = ImageIOHelper.readImageAsGrayScale(filePath).copyToGreyscale();
        canny1 = new CannyEdgeFilterAdaptive();
        canny1.setToUseAlternateNonMaximumSuppression();
        canny1.applyFilter(img);
        MiscDebug.writeImage(img, "_gradient_canny_adaptivw__alt_nms");
    }
    
    public void test1() throws Exception {
        
        String fileName = "blox.gif";
        //String fileName = "lab.gif";
        //String fileName = "house.gif";
        //String fileName = "susan-in_plus.png";
        //String fileName = "lena.jpg";
        String filePath = ResourceFinder.findFileInTestResources(fileName);
        
        GreyscaleImage img = ImageIOHelper.readImageAsGrayScale(filePath).copyToGreyscale();
       
        CannyEdgeFilterAdaptive canny1 = new CannyEdgeFilterAdaptive();        
        canny1 = new CannyEdgeFilterAdaptive();
        canny1.setToUseAlternateNonMaximumSuppression();
        canny1.setToRestoreJunctions();
        canny1.applyFilter(img);
        MiscDebug.writeImage(img, "_gradient_canny_adaptivw__alt_nms");
    }
    
}
