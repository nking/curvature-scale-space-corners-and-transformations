package algorithms.imageProcessing.features;

import algorithms.imageProcessing.util.MiscStats;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

/**
 *
 * @author nichole
 */
public class IntensityFeatureComparisonStats implements 
    Comparable<IntensityFeatureComparisonStats> {
    
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
        
        return compareByCost(other);
    }
    
    public int compareTo2(IntensityFeatureComparisonStats other) {
        
        /*
        if the number of matches is high and the SSDs are low,
        comparison should be by combined intensity stats.
        
        cost is dependent on number of matches so a flaw in current cost is 
        that fewer matches can be lower cost but the matches are truly worse,
        so use SSD when possible.
        */
        
        int n = comparisonStats.size();
        
        int nOther = other.comparisonStats.size();
        
        if (n == nOther) {
                        
            int cc = compareByCost(other);
            int ci = compareByCombinedIntStat(other);
            if (ci != cc) {
                Logger.getLogger(this.getClass().getName()).warning(
                    "comparison by cost differs from comparison by intensity SSD"
                    + " cc=" + cc + " ci=" + ci + 
                    " this.stat=" + this.toString() +
                    " other.stat=" + other.toString()
                );
            }
            
            return ci;
        }
        
        boolean compareByCost = decideByCost(other);

        double avgSSD = MiscStats.calculateCombinedIntensityStat(comparisonStats);
        
        double avgSSDOther = MiscStats.calculateCombinedIntensityStat(
            other.comparisonStats);
        
        if ((n > nOther) && (avgSSD < avgSSDOther)) {
            
            return compareByCombinedIntStat(other);
            
        } else if ((n < nOther) && (avgSSD > avgSSDOther)) {
            
            return compareByCombinedIntStat(other);
        }
        
        if (compareByCost) {
            
            return compareByCost(other);
            
        } else {
            
            return compareByCombinedIntStat(other);
        }
    }
    
    public int compareByCombinedIntStat(IntensityFeatureComparisonStats other) {
        
        double otherCS = MiscStats.calculateCombinedIntensityStat(
            other.getComparisonStats());
        
        double cs = MiscStats.calculateCombinedIntensityStat(comparisonStats);
        
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

    private boolean decideByCost(IntensityFeatureComparisonStats other) {
        
        if (other.getComparisonStats().size() > 2 && comparisonStats.size() > 2) {
            return true;
        }
        
        return false;        
    }

    @Override
    public String toString() {
        
        StringBuilder sb = new StringBuilder();
        
        sb.append("cost=").append(Double.toString(cost)).append(" ")
            .append("adjustedCost").append(Double.toString(adjustedCost))
            .append(" ").append(" scale=").append(Double.toString(scale))
            .append(" stats=[");
        for (FeatureComparisonStat stat : comparisonStats) {
            sb.append(stat.toString()).append(", ");
        }
        sb.append("]");
            
        return sb.toString();
    }
    
}
