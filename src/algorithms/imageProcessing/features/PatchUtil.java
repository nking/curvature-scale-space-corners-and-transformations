package algorithms.imageProcessing.features;

import algorithms.util.PixelHelper;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;

/**
 * carries the current integrated histogram from a region of
 * HOGS, HCPT, or HGS data.
 * Note that orientation is not used.
 * Note that the user must restrict the arguments to
 * being the same origin data.
 * 
 * @author nichole
 */
public class PatchUtil {
    
    private static float eps = 0.000001f;
    
    private final long[] h;
    
    private final TIntSet pixIdxs;
    private final int imgW;
    private final int imgH;
    
    private double sumErrSq = 0;
    
    // this is sqrt(squared sum of totals/count)
    private double blockTotals = 0;
    
    public PatchUtil(int imageWidth, int imageHeight, int nBins) {
        this.h = new long[nBins];
        this.imgW = imageWidth;
        this.imgH = imageHeight;
        this.pixIdxs = new TIntHashSet();
    }
    
    public void add(TIntSet addPixelIndexes, HOGs hogs) {
        
        if (hogs.getNumberOfBins() != h.length) {
            throw new IllegalArgumentException(
               "hog number of bins differs the expected");
        }
        if (hogs.getImageWidth() != imgW) {
            throw new IllegalArgumentException(
               "hog image width differs the expected");
        }
        if (hogs.getImageHeight() != imgH) {
            throw new IllegalArgumentException(
               "hog image height differs the expected");
        }
        if (addPixelIndexes.isEmpty()) {
            return;
        }
        
        int c0 = pixIdxs.size();
        
        // to keep adding to block totals, square and factor by count again
        double tmpBlockTotals = blockTotals;
        if (blockTotals > 0) {
                   
            double norm = 255.f/(blockTotals + eps);

            div(h, norm);
        }
        
        long tmpSum = 0;
        long tmpSumErrSq = 0;
        double maxValue;
        
        long[] tmp = new long[h.length];
        
        int[] xy = new int[2];
        PixelHelper ph = new PixelHelper();
        
        //TODO: correct to use a scan by cell size pattern
        TIntIterator iter = addPixelIndexes.iterator();
        while (iter.hasNext()) {
            int pixIdx = iter.next();
            if (pixIdxs.contains(pixIdx)) {
                continue;
            }
            pixIdxs.add(pixIdx);
            
            ph.toPixelCoords(pixIdx, imgW, xy);
            
            hogs.extractBlock(xy[0], xy[1], tmp);
           
            hogs.add(h, tmp);
                   
            tmpSum = 0;
            maxValue = Double.NEGATIVE_INFINITY;
            for (int j = 0; j < tmp.length; ++j) {
                tmpSum += tmp[j];
                if (tmp[j] > maxValue) {
                    maxValue = tmp[j];
                } 
            }
            
            maxValue += eps;
            
            tmpBlockTotals += tmpSum;
            tmpSumErrSq += ((tmpSum/maxValue)*(tmpSum/maxValue));
        }
        
        int nAdded = pixIdxs.size() - c0;
        int c1 = pixIdxs.size();
        
        if (c1 > 0) {
            this.blockTotals = tmpBlockTotals;
        }  
        
        double norm = 1./(blockTotals + eps);
        float maxBlock = 255.f;
        norm *= maxBlock;
        
        mult(h, norm);
        
        this.sumErrSq *= this.sumErrSq;
        this.sumErrSq *= c0;
        this.sumErrSq += tmpSumErrSq;
        this.sumErrSq /= (double)c1;
        this.sumErrSq = Math.sqrt(sumErrSq);
    }
    
    private void mult(long[] a, double factor) {
        double t;
        for (int j = 0; j < a.length; ++j) {
            t = factor * a[j]; 
            a[j] = (long)t;
        }
    }
    private void div(long[] a, double factor) {
        if (factor == 0) {
            throw new IllegalArgumentException("factor cannot be 0");
        }
        double t;
        for (int j = 0; j < a.length; ++j) {
            t = a[j] / factor; 
            a[j] = (long)t;
        }
    }
   
    public double getAvgErr() {
        return Math.sqrt(sumErrSq);
    }
    
    public long[] getHistogram() {
        return h;
    }
}
