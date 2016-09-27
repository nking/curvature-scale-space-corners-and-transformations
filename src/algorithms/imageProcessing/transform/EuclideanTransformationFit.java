package algorithms.imageProcessing.transform;

import algorithms.util.PairIntArray;
import java.util.List;

/**
 *
 * @author nichole
 */
public class EuclideanTransformationFit extends AbstractTransformationFit {
        
    private TransformationParameters parameters;
    
    // NOTE: these are possibly null
    private PairIntArray left = null;
    private PairIntArray right = null;

    public EuclideanTransformationFit(TransformationParameters theParams,
        List<Integer> theInlierIndexes, List<Double> theErrors, double theTolerance) {
        
        super(theInlierIndexes, theErrors, theTolerance);
        
        parameters = theParams.copy();
    }
    
    public void setLeft(PairIntArray theLeftPoints) {
        this.left = theLeftPoints;
    }
    
    public void setRight(PairIntArray theRightPoints) {
        this.right = theRightPoints;
    }
    
    public PairIntArray getLeft() {
        return left;
    }
    
    public PairIntArray getRight() {
        return right;
    }
    
    public TransformationParameters getTransformationParameters() {
        return parameters;
    }
    
    public void setTransformationParameters(TransformationParameters theParams) {
        parameters = theParams;
    }
}
