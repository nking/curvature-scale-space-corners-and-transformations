package algorithms.imageProcessing;

import algorithms.misc.HistogramHolder;
import java.util.Arrays;

/**
 *
 * @author nichole
 */
public class ImageStatistics {
    
    private double mean = Double.MIN_VALUE;
    
    private double mode = Double.MIN_VALUE;
    
    private double median = Double.MIN_VALUE;
    
    private double min = Double.MIN_VALUE;
    
    private double max = Double.MIN_VALUE;
    
    private float[] quartiles = null;
    
    /*
    if this is set a value less than Double.MAX_VALUE, its the threshold
    used on the image likely, after these stats were calculated
    */
    private double lowThresholdApplied = Double.MAX_VALUE;
    
    private HistogramHolder histogram = null;
    
    public ImageStatistics() {
        
    }

    /**
     * @return the mean
     */
    public double getMean() {
        return mean;
    }

    /**
     * @param theMean the mean to set
     */
    public void setMean(double theMean) {
        this.mean = theMean;
    }

    /**
     * @return the mode
     */
    public double getMode() {
        return mode;
    }

    /**
     * @param theMode the mode to set
     */
    public void setMode(double theMode) {
        this.mode = theMode;
    }

    /**
     * @return the median
     */
    public double getMedian() {
        return median;
    }

    /**
     * @param theMedian the median to set
     */
    public void setMedian(double theMedian) {
        this.median = theMedian;
    }

    /**
     * @return the min
     */
    public double getMin() {
        return min;
    }

    /**
     * @param theMin the min to set
     */
    public void setMin(double theMin) {
        this.min = theMin;
    }

    /**
     * @return the max
     */
    public double getMax() {
        return max;
    }

    /**
     * @param theMax the max to set
     */
    public void setMax(double theMax) {
        this.max = theMax;
    }

    /**
     * @return the histogram
     */
    public HistogramHolder getHistogram() {
        return histogram;
    }

    /**
     * @param theHistogram the histogram to set
     */
    public void setHistogram(HistogramHolder theHistogram) {
        this.histogram = theHistogram;
    }
    
    /**
     * @return the quartiles
     */
    public float[] getQuartiles() {
        return quartiles;
    }

    /**
     * @param theQuartiles the quartiles to set
     */
    public void setQuartiles(float[] theQuartiles) {
        
        this.quartiles = theQuartiles;
    }
    
    public String toString() {
        
        StringBuilder sb = new StringBuilder();
        
        sb.append("min=").append(Double.toString(min)).append("\n");
        sb.append("max=").append(Double.toString(max)).append("\n");
        sb.append("mean=").append(Double.toString(mean)).append("\n");
        sb.append("median=").append(Double.toString(median)).append("\n");
        sb.append("mode=").append(Double.toString(mode)).append("\n");
        sb.append("quartiles=");
        if (quartiles != null) {
            sb.append(Arrays.toString(quartiles));
        }
        sb.append("\n");
        sb.append("histogram=");
        if (histogram != null) {
            sb.append(histogram.toString()).append("\n");
        }
        
        if (Double.compare(lowThresholdApplied, Double.MAX_VALUE) < 0) {
            sb.append("lowThresholdApplied=").append(lowThresholdApplied);
        }
        
        return sb.toString();
    }

    public void setLowThresholdApplied(double lowThresh) {
        lowThresholdApplied = lowThresh;
    }
    
    public double getLowThresholdApplied() {
        return lowThresholdApplied;
    }

}
