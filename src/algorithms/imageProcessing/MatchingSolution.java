package algorithms.imageProcessing;

import java.util.List;

/**
 *
 * @author nichole
 */
public class MatchingSolution {
 
    private final TransformationParameters params;

    private final List<FeatureComparisonStat> comparisonStats;
    
    public MatchingSolution(TransformationParameters parameters,
        List<FeatureComparisonStat> stats) {
        
        this.params = parameters;
        
        this.comparisonStats = stats;
        
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

}
