package algorithms.imageProcessing;

import java.util.List;

/**
 *
 * @author nichole
 */
public class MatchingSolution {
 
    private final TransformationParameters params;

    private final List<FeatureComparisonStat> comparisonStats;
    
    private final int binFactor1, binFactor2;
    
    public MatchingSolution(TransformationParameters parameters,
        List<FeatureComparisonStat> stats, int binFactor1, int binFactor2) {
        
        this.params = parameters;
        
        this.comparisonStats = stats;
        
        this.binFactor1 = binFactor1;
        
        this.binFactor2 = binFactor2;
    }
    
    /**
     * @return the params
     */
    public TransformationParameters getParams() {
        return params;
    }

    /**
     * @return the comparisonStats
     */
    public List<FeatureComparisonStat> getComparisonStats() {
        return comparisonStats;
    }

    /**
     * @return the binFactor1
     */
    public int getBinFactor1() {
        return binFactor1;
    }

    /**
     * @return the binFactor2
     */
    public int getBinFactor2() {
        return binFactor2;
    }

}
