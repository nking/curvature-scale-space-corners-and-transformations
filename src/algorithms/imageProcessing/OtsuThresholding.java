package algorithms.imageProcessing;

import algorithms.misc.Histogram;
import algorithms.misc.MiscMath;
import algorithms.util.PairInt;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * binary thresholding with Otsu's method attempts to separate the pixel values
 * into 2 classes and tries to minimize the variance within each class
 * (which subsequently maximizes the separation between classes).
 * 
 * @author nichole
 */
public class OtsuThresholding {
   
    
    /**
     * find the binary threshold assuming the image is filled with values 
     * between 0 and 255, inclusive.
     * 
     * The implementation follows the one dimensional code from wikipedia
     * at https://en.wikipedia.org/wiki/Otsu%27s_method#cite_note-zhu2009fast-9

     * The runtime complexity is O(N_pixels).
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
     * sparse 1-D implementation of binary thresholding.
     * only image values > 0 are included in stats because the method
     * is meant to be used on thinned images.
     * 
     * The implementation follows the one dimensional code from wikipedia
     * at https://en.wikipedia.org/wiki/Otsu%27s_method

     * The runtime complexity is O(N_pixels).
     * 
     * @param img
     * @return 
     */
    public double calculateBinaryThreshold1DSparse(double[][] img) {
        
        int w = img.length;
        int h = img[0].length;
        
        int nNZ = 0;
        double min = Double.MAX_VALUE;
        for (int i = 0; i < w; ++i) {
            for (int j = 0; j < h; ++j) {
                double v = img[i][j];
                if (v > 0) {
                    if (v < min) {
                        min = v;
                    }
                    nNZ++;
                }
            }
        }
        
        int nBins = 256;
        
        int minBin = (int)Math.round(255. * min);
        
        int[] hist = new int[nBins];
        // (255 - 0 + 1)/256
        int binWidth = (255 - minBin + 1)/nBins;
        
        for (int i = 0; i < w; ++i) {
            for (int j = 0; j < h; ++j) {
                if (img[i][j] <= 0) {
                    continue;
                }
                int v = (int)Math.round(255. * img[i][j]);
                int binNumber = (v - minBin)/binWidth;            
                hist[binNumber]++;
            }
        }
        
        int thresh = calculateBinaryThreshold256(hist, nNZ);
        
        return (double)thresh/255.;
    }
    
    /**
     * find the binary threshold for the given points assuming the image is 
     * filled with values between 0 and 255, inclusive.
     * 
     * The implementation follows the one dimensional code from wikipedia
     * at https://en.wikipedia.org/wiki/Otsu%27s_method#cite_note-zhu2009fast-9

     * The runtime complexity is O(N_pixels).
     * 
     * @param pointValues
     * @return 
     */
    public int calculateBinaryThreshold256(Map<PairInt, Integer> pointValues) {
        
        int nPix = pointValues.size();
        
        int[] h = Histogram.createHistogram(pointValues, 0, 255, 256);
        
        return calculateBinaryThreshold256(h, nPix);
    }
    
    /**
     * one-dimensional binary thresholding.
     * find the binary threshold assuming the image is filled with values 
     * between 0 and 255, inclusive.   
     * The implementation follows the one dimiensional code from wikiepedia
     * at https://en.wikipedia.org/wiki/Otsu%27s_method#cite_note-zhu2009fast-9
     * 
     * The runtime complexity is O(N_pixels).
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
            
            //maximize the difference the between the areas
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
     * The runtime complexity is approximately nLevels * O(N_pixels).
     * 
     * @param img
     * @param nLevels
     */
    public void applyMultiLevelThreshold256(GreyscaleImage img, int nLevels) {
                 
        int[] thresholds = calculateMultiLevelThreshold256(img, nLevels);
        
        for (int i = 0; i < img.getNPixels(); ++i) {
            
            int v = img.getValue(i);
            
            int thresh = thresholds[i];
            
            if (v > thresh) {
                img.setValue(i, 255);
            } else {
                img.setValue(i, 0);
            }
        }
    }
    
    /**
     * apply an adaptive version of binary thresholding to the image img
     * based upon nLevels of histograms created for the nLevels of
     * neighborhood intensity bins for each pixel.
     * 
     * The current window for the neighborhood region is +- 1 pixels.
     * 
     * The runtime complexity is approximately nLevels * O(N_pixels).
     * 
     * returns an array of the thresholds indexed by pixel number.
     * 
     * @param img
     * @param nLevels
     */
    public int[] calculateMultiLevelThreshold256(GreyscaleImage img, int nLevels) {
        
        int nPix = img.getNPixels();
        
        GreyscaleImage avgImg = img.copyImage();
        ImageProcessor imageProcessor = new ImageProcessor();
        int halfWindow = 1;//2;
        imageProcessor.applyCenteredMean(avgImg, halfWindow);
        // subtract the center pixel
        for (int i = 0; i < nPix; ++i) {
            //TODO: revisit this.  looks like it should be img - avgImg 
            // and revisit the normalization
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
            if (binNumberV > (twoDHist[0].length - 1)) {
                binNumberV =twoDHist[0].length - 1;
            }
            
            int vAvg = avgImg.getValue(i);
            int binNumberVAvg = (vAvg - 0)/binWidthAvg;
            
            if (binNumberVAvg > (twoDHist.length - 1)) {
                binNumberVAvg = twoDHist.length - 1;
            }
            
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
        
        int[] pixelThresholds = new int[nPix];
        
        for (int i = 0; i < img.getNPixels(); ++i) {
            int vAvg = avgImg.getValue(i);
            int binNumberVAvg = (vAvg - 0)/binWidthAvg;  
            if (binNumberVAvg > (thresholds.length - 1)) {
                binNumberVAvg = thresholds.length - 1;
            }
            pixelThresholds[i] = thresholds[binNumberVAvg];
        }
        
        return pixelThresholds;
    }
    
    /**
     * NOT YET TESTED.
     * 
     * 2-D fast implementation of 2D binary thresholding.
     * 
     * The runtime complexity is at best O(N_pixels) and 
     * at worst nBins^2 if that is larger than N_pixels.
     * 
     * returns an array of the thresholds indexed by pixel number.
     * 
     * The implementation uses integral images to keep the runtime complexity 
     * low.  That topic is discussed in wikipedia and in
     * https://en.wikipedia.org/wiki/Otsu%27s_method#cite_ref-zhu2009fast_9-0
     * "A Fast 2D Otsu Thresholding Algorithm Based on Improved Histogram"
     * by Zhou, Wang, Yang, and Dai
     * 
     * @param img array with values 0 through 1, inclusive
     * @param nBins the number of bins to use when calculating the inner
     * threshold.  for example, a greyscale image with range 0 to 255 and
     * nBins = 256 would be the finest grain histogram for that data.
     */
    public double calculateBinaryThreshold2D(double[][] img, int nBins) {
        
        int w = img.length;
        int h = img[0].length;
        double nPix = w * h;
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        double[][] avgImg = imageProcessor.copy(img);
        
        int halfWindow = 1;
        imageProcessor.applyCenteredMean2(avgImg, halfWindow);
        
        double min = MiscMath.findMin(img);
        double max = MiscMath.findMax(img);
        double binWidth = (max - min + 1.)/(double)nBins;

        /*
        -- histogram from a mean image made from summed table
           w/ window = 1.
        -- histogram of original pixel intensity.
        -- M is image size
        -- f_i_j is the frequency of pixels with
             i is a pixel intensity bin
             j is a mean pixel intensity bin
        -- probability of pix is P_i_j = f_i_j/M
           assert that sum of f_i_j over all i and j == M
           assert that sum of P_i_j over all i and j == 1
          - would want to keep a set of P_i_j to read the summed
            area table for it sparsely for the calcs below.
        */
        
        double[][] p = new double[nBins][];
        for (int i = 0; i < nBins; ++i) {
            p[i] = new double[nBins];
        }
        Set<PairInt> pSet = new HashSet<PairInt>();
        for (int i = 0; i < w; ++i) {
            for (int j = 0; j < h; ++j) {
                double v = img[i][j];
                int binInt = (int)((v - min)/binWidth);
                v = avgImg[i][j];
                int binAvg = (int)((v - min)/binWidth);
                PairInt pair = new PairInt(binInt, binAvg);
                pSet.add(pair);
                p[binInt][binAvg]++;
            }
        }
                
        double minP = Double.MAX_VALUE;
        // finish array of probabilities
        for (PairInt pair : pSet) {
            int i = pair.getX();
            int j = pair.getY();
            p[i][j] /= nPix;
            if (p[i][j] < minP) {
                minP = p[i][j];
            }
        }
        
        /*
       for any given point in P_i_j, noted as (s, t) where s is i and t is j,
       would want to calculate w and m for rectangle below it
       and above it separately in the summed area table to calculate the trace.
       (the region above should include the point at (s,t)?)
       these points can be found in the sparse set kept for P_i_j.
       
       the w's below are w_0 and are the sum of the P_i_j window
            from i=0 to s-1, j=0 to t-1
       the w's above are w_1
            from i=s to nBins-1, j = t to nBins-1
       the m's below are two different metrics:
             m_0_i = (sum over i of i * P_i_j)/w_0
                 from i=0 to s-1, j = 0 to t-1
             m_0_j = (sum over j of j * P_i_j)/w_0
                 from i=0 to s-1, j = 0 to t-1
       the m's above are two different metrics:
             m_1_i = (sum over i of i * P_i_j)/w_0
                 from i=s to nBins-1, j = t to nBins-1
             m_1_j = (sum over j of j * P_i_j)/w_0
                 from i=s to nBins-1, j = t to nBins-1
       the m_total_i = (sum over i of i * P_i_j)
                 from i=0 to nBins-1, j = 0 to nBins-1
       the m_total_j = (sum over j of j * P_i_j)
                 from i=0 to nBins-1, j = 0 to nBins-1
        
       S_b = inter class discrete matrix
           = sum over 0, and 1 for the classes of w
           = w_0 * ((m_0 - m_total_0)*(m_0 - m_total_0)^T)
             + w_1 * ((m_1 - m_total_1)*(m_1 - m_total_1)^T)
        
       the optimal threshold for the point (s,t) is the maximum
       of the trace.
        
       tr(S_b) = w_0 * ((m_0_i - m_total_i)^2 + (m_0_j - m_total_j)^2) 
                  + w_1 * ((m_1_i - m_total_i)^2 + (m_1_j - m_total_j)^2)
          
          which can be further simplified to 
          
          tr(S_b) = ( (m_total_i * w_0 - m_i)^2 + (m_total_j * w_0 - m_j)^2)
                  / (w_0*(1-w_0))
            where m_i is (sum over i of i * P_i_j)/w_0
                      i=0 to s-1 and j=0 to t-1
            where m_j is (sum over j of j * P_i_j)/w_0
                      i=0 to s-1 and j=0 to t-1
            where m_total_i = (sum over i of i * P_i_j)
                      i=0 to nBins-1, j = 0 to nBins-1
            where m_total_j = (sum over j of j * P_i_j)
                      i=0 to nBins-1, j = 0 to nBins-1
            where w_0 is sum of the P_i_j window
                      i=0 to s-1, j=0 to t-1
        
       ==> need summed area table for P_i_j
           need summed area table for i * P_i_j
           need summed area table for j * P_i_j
        */
        
        // summed area table for p_i_j
        SummedAreaTable summed = new SummedAreaTable();
        double[][] pTable = summed.createAbsoluteSummedAreaTable(p);
        
        assert(allAreRealNumbers(pTable));
        
        double[] binFactors = new double[nBins];
        for (int i = 0; i < nBins; ++i) {
            binFactors[i] = i * binWidth + min;
        }     
       
        //TODO: this needs corrections to avoid overrunning
        // bounds of table data type.
        // TODO: also need to verify that the simplification of the trace
        //  printed in wikipedia is correct, else need to correct that
        //  or return to the longer trace statement that includes w0 and w1
        //  terms.
        
        // the other 2 summed area tables:
        double[][] iPTable = imageProcessor.copy(p);
        double[][] jPTable = imageProcessor.copy(p);
        for (int i = 0; i < nBins; ++i) {
            double iFactor = binFactors[i];
            for (int j = 0; j < nBins; ++j) {
                double jFactor = binFactors[j];
                double iv = 0;
                double jv = 0;
                if (i > 0 && j > 0) {
                    iv = iPTable[i - 1][j] + iPTable[i][j - 1] 
                        - iPTable[i - 1][j - 1] + iPTable[i][j];
                    jv = jPTable[i - 1][j] + jPTable[i][j - 1] 
                        - jPTable[i - 1][j - 1] + jPTable[i][j];
                    if (Double.isNaN(jv)) {
                        System.out.println("jPTable[i - 1][j]=" + jPTable[i - 1][j]
                            + " jPTable[i][j - 1]=" + jPTable[i][j - 1] 
                            + " jPTable[i - 1][j - 1]=" + jPTable[i - 1][j - 1]
                            + " jPTable[i][j]=" + jPTable[i][j]
                        );
                    }
                } else if (i > 0) {
                    iv = iPTable[i - 1][j] + iPTable[i][j];
                    jv = jPTable[i - 1][j] + jPTable[i][j];
                } else if (j > 0) {
                    iv = iPTable[i][j - 1] + iPTable[i][j];
                    jv = jPTable[i][j - 1] + jPTable[i][j];
                }
                iPTable[i][j] = iFactor * iv;
                jPTable[i][j] = jFactor * jv;
                assert(!Double.isNaN(iPTable[i][j]));
                assert(!Double.isNaN(jPTable[i][j]));
                assert(Double.isFinite(iPTable[i][j]));
                assert(Double.isFinite(jPTable[i][j]));
            }
        }
        
        /*
        tr(S_b) = ( (m_total_i * w_0 - m_i)^2 + (m_total_j * w_0 - m_j)^2)
                  / (w_0*(1-w_0))
            where m_i is (sum over i of i * P_i_j)/w_0
                      i=0 to s-1 and j=0 to t-1
            where m_j is (sum over j of j * P_i_j)/w_0
                      i=0 to s-1 and j=0 to t-1
            where m_total_i = (sum over i of i * P_i_j)
                      i=0 to nBins-1, j = 0 to nBins-1
            where m_total_j = (sum over j of j * P_i_j)
                      i=0 to nBins-1, j = 0 to nBins-1
            where w_0 is sum of the P_i_j window
                      i=0 to s-1, j=0 to t-1
        */
        
        //TODO: handle multiplication by factor over window size
        
        double mTotalI = iPTable[nBins - 1][nBins - 1];
        double mTotalJ = jPTable[nBins - 1][nBins - 1];
        
        PairInt maxPair = null;
        double maxTrace = Double.MIN_VALUE;
        
        for (PairInt pair : pSet) {
            int s = pair.getX();
            int t = pair.getY();
            
            double w0, mI, mJ;
            if (s > 0 && t > 0) {
                w0 = pTable[s - 1][t - 1];
                mI = iPTable[s - 1][t - 1];
                mJ = jPTable[s - 1][t - 1];
            } else if (s > 0) {
                assert(t == 0);
                w0 = pTable[s - 1][t];
                mI = iPTable[s - 1][t];
                mJ = jPTable[s - 1][t];
            } else if (t > 0) {
                assert(s == 0);
                w0 = pTable[s][t - 1];
                mI = iPTable[s][t - 1];
                mJ = jPTable[s][t - 1];
            } else {
                // i == 0 and j == 0
                assert(s == 0);
                assert(t == 0);
                w0 = p[s][t];
                mI = iPTable[s][t];
                mJ = jPTable[s][t];
            }
            assert(w0 <= 1.);
            mI /= w0;
            mJ /= w0;
            
            double a = mTotalI * w0 - mI;
            a *= a;
            double b = mTotalJ * w0 - mJ;
            b *= b;
            assert(Double.isFinite(w0 * (1 - w0)));
            assert(!Double.isNaN(w0 * (1 - w0)));
            double trace = (a + b)/(w0 * (1 - w0));
            
            if (trace > maxTrace) {
                maxTrace = trace;
                maxPair = pair;
            }
        }
        
        double thresh = (maxPair.getX() * binWidth) + min;
        thresh += (binWidth/2);
        
        return thresh;
    }

    private boolean allAreRealNumbers(double[][] a) {
        
        int w = a.length;
        int h = a[0].length;
        for (int i = 0; i < w; ++i) {
            for (int j = 0; j < h; ++j) {
                if (Double.isNaN(a[i][j])) {
                    return false;
                }
                if (Double.isInfinite(a[i][j])) {
                    return false;
                }
            }
        }
        return true;
    }
}
