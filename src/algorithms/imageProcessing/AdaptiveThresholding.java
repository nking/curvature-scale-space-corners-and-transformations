package algorithms.imageProcessing;

import algorithms.misc.MiscMath;

/**
 A class to create an adaptive threshold value image
 based upon the use of integral images
 and using techniques similar to adaptive mean,
 but including standard deviation.
 
 The algorithm follows the recipe in:
 
  "Efficient Implementation of Local Adaptive Thresholding
  Techniques Using Integral Images"
   by Shafaita, Keysersa, and Breuelb

 * @author nichole
 */
public class AdaptiveThresholding {
    
    /**
     * create an image containing a threshold for each pixel where the threshold
     * can be used for binarization.
     * 
     * runtime complexity is O(N).
     * 
     * @param img an image with all non-negative values.
     * @param windowSize the full width of a square window to be used for the
     * adaptive grid.  authors use w = 15;
     * @param k, authors use k=0.2.  increasing k makes the result more sensitive.
     * @return 
     */
    public double[][] createAdaptiveThresholdImage(double[][] img,
        int windowSize, double k) {
                
        SummedAreaTable summed = new SummedAreaTable();
        
        double[][] mTable = summed.createAbsoluteSummedAreaTable(img);
        mTable = summed.applyMeanOfWindowFromSummedAreaTable(mTable, windowSize);
        
        double min = MiscMath.findMin(img);
        double max = MiscMath.findMax(img);
        double R = (max + min)/2.;
        
        int w = img.length;
        int h = img[0].length;
        
        /*
        st dev in window is 
            sum of diff^2 where diff = mean - pix
        then sqrt(diff/(nWindow-1)) <-- will use nWindow instead of nWindow-1
        for convenience, but should consider correcting this one day.
        */
        
        double[][] sTable = new double[w][];
        for (int i = 0; i < w; ++i) {
            sTable[i] = new double[h];
            for (int j = 0; j < h; ++j) {
                sTable[i][j] = img[i][j] - mTable[i][j];
                sTable[i][j] *= sTable[i][j];
            }
        }
        sTable = summed.createAbsoluteSummedAreaTable(sTable);
        sTable = summed.applyMeanOfWindowFromSummedAreaTable(sTable, windowSize);
        
        double[][] tImg = new double[w][];
        for (int i = 0; i < w; ++i) {
            tImg[i] = new double[h];
            for (int j = 0; j < h; ++j) {
                //= m(x,y) * (1 + k*((s(x,y)/R) - 1)
                double s = Math.sqrt(sTable[i][j]);
                tImg[i][j] = mTable[i][j] * (1 + (k * ((s/R) - 1)));
            }
        }
        
        return tImg;
    }
    
    /**
     * apply the adaptive threshold to the given image.
     * 
     * runtime complexity is O(N).
     * 
     * @param img an image with all non-negative values.
     * @param windowSize the full width of a square window to be used for the
     * adaptive grid.  authors use w = 15;
     * @param k, authors use k=0.2
     * @param highValue the value to set img to when larger than threshold.
     */
    public void applyAdaptiveThresholdImage(double[][] img,
        int windowSize, double k, double highValue) {
                        
        double[][] threshs = createAdaptiveThresholdImage(img, windowSize, k);
        
        int w = img.length;
        int h = img[0].length;
        
        for (int i = 0; i < w; ++i) {
            for (int j = 0; j < h; ++j) {
                if (img[i][j] > threshs[i][j]) {
                    img[i][j] = highValue;
                } else {
                    img[i][j] = 0;
                }
            }
        }
    }
}
