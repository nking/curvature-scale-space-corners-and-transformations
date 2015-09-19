package algorithms.imageProcessing;

import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author nichole
 */
public class IntensityFeatureComparisonStats implements Comparable<IntensityFeatureComparisonStats> {
    
    private final double cost;
    
    private double adjustedCost;
    
    private final double scale;
    
    private final int idx1;
    
    private final int idx2;
    
    private final List<FeatureComparisonStat> comparisonStats 
        = new ArrayList<FeatureComparisonStat>();
    
    public IntensityFeatureComparisonStats(final int index1, final int index2,
        double solutionCost, double solutionScale) {
        cost = solutionCost;
        scale = solutionScale;
        idx1 = index1;
        idx2 = index2;
        adjustedCost = cost;
    }
    
    public void addAll(List<FeatureComparisonStat> stats) {
        comparisonStats.addAll(stats);
    }
    
    public List<FeatureComparisonStat> getComparisonStats() {
        return comparisonStats;
    }

    @Override
    public int compareTo(IntensityFeatureComparisonStats other) {
        
        //int comb0 = compareByCombinedIntStat(other);
        
        int comb1 = compareByCost(other);
        
        return comb1;
    }
    
    public int compareByCombinedIntStat(IntensityFeatureComparisonStats other) {
        
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
    
    /**
     * note, should only compare by cost if there is an edge in common, else
     * a correction has to be made for difference in peak sigma before use
     * here.
     * @param other
     * @return 
     */
    public int compareByCost(IntensityFeatureComparisonStats other) {
        
        double otherCost = other.adjustedCost;
        
        if (otherCost == adjustedCost) {
            return 0;
        } else if (adjustedCost < otherCost) {
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

    /**
     * @return the cost
     */
    public double getCost() {
        return cost;
    }
    
    public void setAdjustedCost(double adjCost) {
        adjustedCost = adjCost;
    }
    
    public double getAdjustedCost() {
        return adjustedCost;
    }

    /**
     * @return the scale
     */
    public double getScale() {
        return scale;
    }
    
    public int getIndex1() {
        return idx1;
    }
    
    public int getIndex2() {
        return idx2;
    }
}
