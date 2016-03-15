package algorithms.imageProcessing;

import algorithms.misc.Histogram;
import java.util.Arrays;

/**
 * binary thresholding with Otsu's method attempts to separate the pixel values
 * into 2 classes and tries to minimize the variance within each class
 * (which subsequently maximizes the separation between classes).
 * 
 * https://en.wikipedia.org/wiki/Otsu%27s_method#cite_note-zhu2009fast-9
 * 
 * @author nichole
 */
public class OtsuThresholding {
   
    /**
     * find the binary threshold assuming the image is filled with values 
     * between 0 and 255, inclusive.
     * runtime complexity is O(N_pixels);
     * 
     * @param img
     * @return 
     */
    public int calculateBinaryThreshold256(GreyscaleImage img) {
        
        int nPix = img.getNPixels();
        
        int[] h = Histogram.createHistogram(img, 0, 255, 256);
        
        return calculateBinaryThreshold256(h, nPix);
    }
    
    /**
     * find the binary threshold assuming the image is filled with values 
     * between 0 and 255, inclusive.   runtime complexity is O(N_pixels);
     * 
     * @param h histogram
     * @param nPix the original image's number of pixels
     * @return 
     */
    private int calculateBinaryThreshold256(int[] h, int nPix) {
        
        if (h.length != 256) {
            throw new IllegalArgumentException(
                "expecting histogram length to be 256");
        }
        
        double totalCumulativeSumHist = 0;
        
        for (int i = 1; i < 256; ++i) {
            totalCumulativeSumHist += (i * h[i]);
        }
     
        double cumulativeSumHist = 0;
        double sumHist = 0;        
        double totalMinusSumHist = 0;
        double mean1 = 0;
        double mean2 = 0;
        
        double max = 0;        
        double between = 0;
        int threshold = 0;
        
        for (int i = 0; i < 256; ++i) {
            sumHist += h[i];
            if (sumHist == 0) {
                continue;
            }
            totalMinusSumHist = nPix - sumHist;
            if (totalMinusSumHist == 0) {
                break;
            }
            
            cumulativeSumHist += (i * h[i]);
            mean1 = cumulativeSumHist / sumHist;
            mean2 = (totalCumulativeSumHist - cumulativeSumHist) / totalMinusSumHist;
            
            //maximize inter-class variance
            between = sumHist * totalMinusSumHist * Math.pow(mean1 - mean2, 2);
            if (between > max) {
                threshold = i;
                max = between;
            }
        }
        
        return threshold;
    }
    
    /**
     * apply an adaptive version of binary thresholding to the image img
     * based upon nLevels of histograms created for the nLevels of
     * neighborhood intensity bins for each pixel.
     * 
     * The current window for the neighborhood region is +- 1 pixels.
     * 
     * The runtime complexity is approximately nLevels * O(N+pixels).
     * 
     * @param img
     * @param nLevels
     */
    public void applyMultiLevelThreshold256(GreyscaleImage img, int nLevels) {
        
        int nPix = img.getNPixels();
        
        GreyscaleImage avgImg = img.copyImage();
        ImageProcessor imageProcessor = new ImageProcessor();
        int halfWindow = 1;//2;
        imageProcessor.applyCenteredMean(avgImg, halfWindow);
        // subtract the center pixel
        for (int i = 0; i < nPix; ++i) {
            double v = avgImg.getValue(i) - ((double)img.getValue(i)/9.);
            if (v < 0) {
                avgImg.setValue(i, 0);
            } else {
                avgImg.setValue(i, (int)Math.round(v));
            }
        }
        
        int[][] twoDHist = new int[nLevels][];
        for (int i = 0; i < nLevels; ++i) {
            twoDHist[i] = new int[256];
        }
                
        int binWidthAvg = (255 - 0 + 1)/nLevels;
        int binWidth = (255 - 0 + 1)/256;
        
        for (int i = 0; i < img.getNPixels(); ++i) {
            int v = img.getValue(i);
            int binNumberV = (v - 0)/binWidth;
            
            int vAvg = avgImg.getValue(i);
            int binNumberVAvg = (vAvg - 0)/binWidthAvg;
            
            twoDHist[binNumberVAvg][binNumberV]++;
        }
        
        int[] thresholds = new int[nLevels];
        for (int i = 0; i < nLevels; ++i) {
            //total number of counts in a histogram has to be nPix for algorithm
            int nTot = 0;
            int[] h = twoDHist[i];
            for (int j = 0; j < h.length; ++j) {
                nTot += h[j];
            }
            if (nTot != nPix) {
                float scale = (float)nPix/(float)nTot;
                for (int j = 0; j < h.length; ++j) {
                    h[j] = (int)(h[j] * scale);
                }
            }
            thresholds[i] = calculateBinaryThreshold256(h, nPix);
        }
        
        for (int i = 0; i < img.getNPixels(); ++i) {
            int vAvg = avgImg.getValue(i);
            int binNumberVAvg = (vAvg - 0)/binWidthAvg;
            
            int v = img.getValue(i);
            int thresh = thresholds[binNumberVAvg];
            if (v > thresh) {
                img.setValue(i, 255);
            } else {
                img.setValue(i, 0);
            }
        }
    }
    
}
