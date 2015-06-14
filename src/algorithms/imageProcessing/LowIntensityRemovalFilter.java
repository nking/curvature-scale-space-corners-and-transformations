package algorithms.imageProcessing;

import algorithms.misc.HistogramHolder;
import java.util.logging.Logger;

/**
 *
 * @author nichole
 */
public class LowIntensityRemovalFilter {
    
    private Logger log = Logger.getLogger(this.getClass().getName());

    public static float defaultLowThreshFactor = 10.9f;//1.0f;
    
    private float lowThreshFactor = 0.9f;//1.0f;

    public void overrideLowThresholdFactor(float theFactor) {
        lowThreshFactor = theFactor;
    }
    
    /**
     * determine a lower threshold for removing pixels below that
     * as background and noise and return the threshold used.
     * @param input
     * @param originalImageHistogram
     * @return 
     */
    public ImageStatistics removeLowIntensityPixels(final GreyscaleImage input, 
        HistogramHolder originalImageHistogram) {
        
        boolean useSturges = false;
                
        ImageStatistics stats = ImageStatisticsHelper.examineImageBorders(
            input, useSturges);
        //ImageStatistics stats = ImageStatisticsHelper.examineImage(
        //    input, useSturges);
        
        int lowThresh = determineLowThreshold(stats, originalImageHistogram);

        stats.setLowThresholdApplied(lowThresh);
        
        log.fine(stats.toString());
        
        if (lowThresh > 0) {
            
            removeLowIntensityPixels(input, lowThresh);
        }
        
        return stats;
    }

    protected int determineLowThreshold(ImageStatistics stats, 
        HistogramHolder originalImageHistogram) {
                 
        return (int) (lowThreshFactor * stats.getMean());
     
    }
        
    public void removeLowIntensityPixels(final GreyscaleImage input, 
        int threshold) {
       
        for (int i = 0; i < input.getWidth(); i++) {
            for (int j = 0; j < input.getHeight(); j++) {
                
                int v = input.getValue(i, j);
                
                if (v < threshold) {
                    input.setValue(i, j, 0);
                }
            }
        }
    }
    
}
