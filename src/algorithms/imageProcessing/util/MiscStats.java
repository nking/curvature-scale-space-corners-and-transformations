package algorithms.imageProcessing.util;

import algorithms.imageProcessing.features.FeatureComparisonStat;
import java.util.List;

/**
 *
 * @author nichole
 */
public class MiscStats {

    public static double calculateCombinedIntensityStat(List<FeatureComparisonStat> 
        compStats) {
        
        if (compStats.isEmpty()) {
            return Double.POSITIVE_INFINITY;
        }
        
        double sum = 0;
        
        for (FeatureComparisonStat compStat : compStats) {
            sum += compStat.getSumIntensitySqDiff();
        }
        
        sum /= (double) compStats.size();
        
        return sum;
    }
   
}
