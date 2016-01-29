package algorithms.imageProcessing.transform;

import algorithms.imageProcessing.features.FeatureComparisonStat;
import algorithms.imageProcessing.matching.ErrorType;
import java.util.ArrayList;
import java.util.List;
import org.ejml.simple.SimpleMatrix;

/**
 *
 * @author nichole
 */
public class EpipolarFeatureTransformationFit extends EpipolarTransformationFit {
    
    private final List<FeatureComparisonStat> fcs;
    
    public EpipolarFeatureTransformationFit(SimpleMatrix theFundamentalMatrix,
        List<Integer> theInlierIndexes, List<FeatureComparisonStat> stats,
        ErrorType theErrorType,
        List<Double> theErrors, double theTolerance) {
        
        super(theFundamentalMatrix, theInlierIndexes, theErrorType, theErrors, 
            theTolerance);
        
        this.fcs = new ArrayList<FeatureComparisonStat>(stats);
    }
    
    public List<FeatureComparisonStat> getFeatureComparisonStats() {
        return fcs;
    }
}
