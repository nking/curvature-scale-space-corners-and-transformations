package algorithms.imageProcessing;

import algorithms.imageProcessing.matching.ErrorType;
import java.util.List;
import org.ejml.simple.SimpleMatrix;

/**
 *
 * @author nichole
 */
public class EpipolarTransformationFit extends AbstractTransformationFit {
    
    private final ErrorType errorType;
    
    private SimpleMatrix fundamentalMatrix;

    public EpipolarTransformationFit(SimpleMatrix theFundamentalMatrix,
        List<Integer> theInlierIndexes, ErrorType theErrorType,
        List<Double> theErrors, double theTolerance) {
        
        super(theInlierIndexes, theErrors, theTolerance);
        
        fundamentalMatrix = theFundamentalMatrix.copy();
        errorType = theErrorType;
    }
    
    public SimpleMatrix getFundamentalMatrix() {
        return fundamentalMatrix;
    }

    /**
     * @return the errorType
     */
    public ErrorType getErrorType() {
        return errorType;
    }
    
    public void setFundamentalMatrix(SimpleMatrix theFM) {
        fundamentalMatrix = theFM;
    }
}
