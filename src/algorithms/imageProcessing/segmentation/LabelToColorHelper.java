package algorithms.imageProcessing.segmentation;

import algorithms.imageProcessing.ImageExt;
import algorithms.misc.MiscMath;

/**
 *
 * @author nichole
 */
public class LabelToColorHelper {
    
    /**
     * calculate the average r,g,b of pixels grouped by their labels and
     * reassigne those pixels the average colors.
     * @param img
     * @param labels 
     */
    public static void applyLabels(ImageExt img, int[] labels) {
        
        if (img.getNPixels() != labels.length) {
            throw new IllegalArgumentException("labels.length must equal img.nPixels");
        }
        
        int maxLabel = MiscMath.findMax(labels);
        
        long[] rSum = new long[maxLabel + 1];
        long[] gSum = new long[maxLabel + 1];
        long[] bSum = new long[maxLabel + 1];
        
        int[] count = new int[rSum.length];
        
        for (int i = 0; i < labels.length; ++i) {
            int label = labels[i];
            rSum[label] += img.getR(i);
            gSum[label] += img.getG(i);
            bSum[label] += img.getB(i);
            
            count[label]++;
        }
        
        for (int i = 0; i < rSum.length; ++i) {
            if (count[i] > 0) {
                rSum[i] /= count[i];
                gSum[i] /= count[i];
                bSum[i] /= count[i];
            }
        }
        
        img.fill(0, 0, 0);
        
        for (int i = 0; i < labels.length; ++i) {
            int label = labels[i];
            img.setRGB(i, (int)rSum[label], (int)gSum[label], (int)bSum[label]);
        }
    }
}