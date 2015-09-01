package algorithms.imageProcessing;

import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author nichole
 */
public class IntensityFeatureComparisonStats implements Comparable<IntensityFeatureComparisonStats> {
    
    private final List<FeatureComparisonStat> comparisonStats 
        = new ArrayList<FeatureComparisonStat>();
    
    public void addAll(List<FeatureComparisonStat> stats) {
        comparisonStats.addAll(stats);
    }
    
    public List<FeatureComparisonStat> getComparisonStats() {
        return comparisonStats;
    }

    @Override
    public int compareTo(IntensityFeatureComparisonStats other) {
        
        double otherCS = calculateCombinedIntensityStat(other.getComparisonStats());
        
        double cs = calculateCombinedIntensityStat(comparisonStats);
        
        if (otherCS == cs) {
            return 0;
        } else if (cs < otherCS) {
            // prefer the smaller
            return -1;
        }
        return 1;
    }
    
    private double calculateCombinedIntensityStat(List<FeatureComparisonStat> compStats) {

        if (compStats.isEmpty()) {
            return Double.POSITIVE_INFINITY;
        }
        
        double sum = 0;

        for (FeatureComparisonStat compStat : compStats) {
            sum += compStat.getSumIntensitySqDiff();
        }

        sum /= (double)compStats.size();

        return sum;
    }
}
