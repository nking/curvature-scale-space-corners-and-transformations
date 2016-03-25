package algorithms.imageProcessing.features;

import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.HistogramEqualization;
import algorithms.imageProcessing.ImageIOHelper;
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
        //String fileName = "susan-in_plus.png";
        String fileName = "tmp.png";
        String filePath = ResourceFinder.findFileInTestResources(fileName);
        
        GreyscaleImage img = ImageIOHelper.readImageAsGrayScaleAvgRGB(filePath);
        
        PhaseCongruencyDetector phaseCDetector = new PhaseCongruencyDetector();
        PhaseCongruencyDetector.PhaseCongruencyProducts products =
            phaseCDetector.phaseCongMono(img);
        
        assertNotNull(products);
        
        double[][] pc = products.getPhaseCongruency();
        double[] pc0 = new double[pc.length * pc[0].length];
        for (int i = 0; i < pc0.length; ++i) {
            int y = i/pc.length;
            int x = i - (y * pc.length);
            pc0[i] = pc[x][y];
        }
        int[] scaledPC = MiscMath.rescale(pc0, 0, 255);
        GreyscaleImage out = new GreyscaleImage(pc.length, pc[0].length);
        for (int i = 0; i < pc0.length; ++i) {
            int y = i/pc.length;
            int x = i - (y * pc.length);            
            out.setValue(x, y, scaledPC[i]);
        }
        MiscDebug.writeImage(out, "_phase_congr_0_");
        
        /*
        [PC, or] = phasecongmono(imread('lena.tif'));
        nm = nonmaxsup(PC, or, 1.5);   % nonmaxima suppression
        bw = hysthresh(nm, 0.1, 0.3);  % hysteresis thresholding 0.1 - 0.3
        */
        NonMaximumSuppression ns = new NonMaximumSuppression();
        double[][] imgTh0 = ns.nonmaxsup(products.getPhaseCongruency(), 
            products.getOrientation(), 1.5);
                
        int[][] binaryImage = phaseCDetector.applyHysThresh(imgTh0, 0.1, 0.3);
                
        out = new GreyscaleImage(binaryImage.length, binaryImage[0].length);
        for (int i = 0; i < binaryImage.length; ++i) {
            for (int j = 0; j < binaryImage[i].length; ++j) {
                if (binaryImage[i][j] > 0) {
                    out.setValue(i, j, 255);
                }
            }
        }
        
        MiscDebug.writeImage(out, "_hyst_");
        int z = 1;
    }
    
}
