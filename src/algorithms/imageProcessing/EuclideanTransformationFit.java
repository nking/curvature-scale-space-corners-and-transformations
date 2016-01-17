package algorithms.imageProcessing;

import java.util.List;

/**
 *
 * @author nichole
 */
public class EuclideanTransformationFit extends AbstractTransformationFit {
        
    private TransformationParameters parameters;

    public EuclideanTransformationFit(TransformationParameters theParams,
        List<Integer> theInlierIndexes, List<Double> theErrors, double theTolerance) {
        
        super(theInlierIndexes, theErrors, theTolerance);
        
        parameters = theParams.copy();
    }
    
    public TransformationParameters getTransformationParameters() {
        return parameters;
    }
    
    public void setTransformationParameters(TransformationParameters theParams) {
        parameters = theParams;
    }
}
