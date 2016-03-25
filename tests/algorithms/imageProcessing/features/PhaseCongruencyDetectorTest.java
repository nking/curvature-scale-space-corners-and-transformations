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
        //String fileName = "tmp.png";
        String fileName = "lena.jpg";
        String filePath = ResourceFinder.findFileInTestResources(fileName);
        
        GreyscaleImage img = ImageIOHelper.readImageAsGrayScaleAvgRGB(filePath);
        
        PhaseCongruencyDetector phaseCDetector = new PhaseCongruencyDetector();
        PhaseCongruencyDetector.PhaseCongruencyProducts products =
            phaseCDetector.phaseCongMono(img);
        
        assertNotNull(products);
                
        double[][] pc = products.getPhaseCongruency();
        
        int nRows = pc.length;
        int nCols = pc[0].length;
        
        double[] pc0 = new double[nRows * nCols];
        for (int i = 0; i < pc0.length; ++i) {
            int y = i/nCols;
            int x = i - (y * nCols);
            pc0[i] = pc[y][x];
        }
        int[] scaledPC = MiscMath.rescale(pc0, 0, 255);
        GreyscaleImage out = new GreyscaleImage(nCols, nRows);
        for (int i = 0; i < pc0.length; ++i) {
            int y = i/nCols;
            int x = i - (y * nCols);           
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
                
        out = new GreyscaleImage(nCols, nRows);
        for (int row = 0; row < nRows; ++row) {
            for (int col = 0; col < nCols; ++col) {
                if (binaryImage[row][col] > 0) {
                    out.setValue(col, row, 255);
                }
            }
        }
        
        MiscDebug.writeImage(out, "_hyst_");
        int z = 1;
    }
    
}
