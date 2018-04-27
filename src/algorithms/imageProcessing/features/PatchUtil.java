package algorithms.imageProcessing.features;

import algorithms.util.PairInt;
import algorithms.util.PixelHelper;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.HashSet;
import java.util.Set;
import algorithms.packing.Intersection2DPacking;

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
    
    private double err = 0;
    
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
           
            HOGs.add(h, tmp);
                   
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
            tmpSum /= tmp.length; 
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
        
        //TODO: examine the order of divide by count and sqrt
        this.err *= this.err;
        this.err *= c0;
        this.err += tmpSumErrSq;
        this.err /= (double)c1;
        this.err = Math.sqrt(err);
    }
    
    public void add(TIntSet addPixelIndexes, HCPT hcpt) {
        
        if (hcpt.getNumberOfBins() != h.length) {
            throw new IllegalArgumentException(
               "hog number of bins differs the expected");
        }
        if (hcpt.getImageWidth() != imgW) {
            throw new IllegalArgumentException(
               "hog image width differs the expected");
        }
        if (hcpt.getImageHeight() != imgH) {
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
            
            hcpt.extractBlock(xy[0], xy[1], tmp);
           
            HOGs.add(h, tmp);
                   
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
            tmpSum /= tmp.length; 
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
        
        //TODO: examine the order of divide by count and sqrt
        this.err *= this.err;
        this.err *= c0;
        this.err += tmpSumErrSq;
        this.err /= (double)c1;
        this.err = Math.sqrt(err);
    }
    
    public void add(TIntSet addPixelIndexes, HGS hgs) {
        
        if (hgs.getNumberOfBins() != h.length) {
            throw new IllegalArgumentException(
               "hog number of bins differs the expected");
        }
        if (hgs.getImageWidth() != imgW) {
            throw new IllegalArgumentException(
               "hog image width differs the expected");
        }
        if (hgs.getImageHeight() != imgH) {
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
            
            hgs.extractBlock(xy[0], xy[1], tmp);
           
            HOGs.add(h, tmp);
                   
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
            tmpSum /= tmp.length; 
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
        
        //TODO: examine the order of divide by count and sqrt
        this.err *= this.err;
        this.err *= c0;
        this.err += tmpSumErrSq;
        this.err /= (double)c1;
        this.err = Math.sqrt(err);
    }
    
    /**
     * returns the set of points that is a subset of the intersection, scanned
     * from center to perimeter by interval of size hog cell.
     * @param other
     * @param nCellSize
     * @return 
     */
    public TIntSet calculateDetectorWindow(PatchUtil other, int nCellSize) {
        
        Intersection2DPacking ip = new Intersection2DPacking();
        
        TIntSet seeds = ip.naiveStripPacking(pixIdxs, other.pixIdxs, imgW, 
            nCellSize);
        
        return seeds;
    }
    
    /**
     * calculate the intersection of the histograms. histograms that are
     * identical have a result of 1.0 and histograms that are completely
     * different have a result of 0.
     * 
     * @param other
     * @return 
     */
    public double intersection(PatchUtil other) {
        
        if ((h.length != other.h.length)) {
            throw new IllegalArgumentException(
                "h and other.h must be same dimensions");
        }
        
        int nBins = h.length;

        double sum = 0;
        double sumA = 0;
        double sumB = 0;
        for (int j = 0; j < nBins; ++j) {

            long yA = h[j];
            long yB = other.h[j];

            sum += Math.min(yA, yB);
            sumA += yA;
            sumB += yB;

            //System.out.println(" " + yA + " -- " + yB + " sum="+sum + ", " + sumA + "," + sumB);
        }

        double d = eps + Math.min(sumA, sumB);
        double sim = sum/d;

        return sim;
    }
    
    /**
     * calculate the difference of the histograms. histograms that are
     * identical have a result of 0.0 and histograms that are completely
     * different have a result of 1.
     * 
     * @param other
     * @return 
     */
    public double[] diff(PatchUtil other) {
        
        if ((h.length != other.h.length)) {
            throw new IllegalArgumentException(
                "h and other.h must be same dimensions");
        }
        
        int nBins = h.length;

        double tmpSumDiff = 0;
        double tmpErr = 0;
        
        for (int j = 0; j < nBins; ++j) {

            long yA = h[j];
            long yB = other.h[j];

            float maxValue = Math.max(yA, yB) + eps;

            float diff = Math.abs((yA - yB)/maxValue);
            
            tmpSumDiff += diff;

            //      already squared
            tmpErr += (diff/maxValue);
        }
        
        tmpSumDiff /= (double)nBins;
                
        tmpErr /= (double)nBins;
        tmpErr = Math.sqrt(tmpErr);
        
        return new double[]{tmpSumDiff, tmpErr};
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
        return Math.sqrt(err);
    }
    
    public long[] getHistogram() {
        return h;
    }
    
    public TIntSet getPixelIndexes() {
        return pixIdxs;
    }
    
    public Set<PairInt> getPixelSet() {
        PixelHelper ph = new PixelHelper();
        int[] xy = new int[2];
        Set<PairInt> set = new HashSet<PairInt>();
        TIntIterator iter = pixIdxs.iterator();
        while (iter.hasNext()) {
            int pixIdx = iter.next();
            ph.toPixelCoords(pixIdx, imgW, xy);
            set.add(new PairInt(xy[0], xy[1]));
        }
        return set;
    }
}
