package algorithms.imageProcessing.transform;

import algorithms.imageProcessing.matching.ErrorType;
import java.util.List;
import no.uib.cipr.matrix.DenseMatrix;

/**
 *
 * @author nichole
 */
public class EpipolarTransformationFit extends AbstractTransformationFit {
    
    private final ErrorType errorType;
    
    private DenseMatrix fundamentalMatrix;
    
    public EpipolarTransformationFit(DenseMatrix theFundamentalMatrix,
        List<Integer> theInlierIndexes, ErrorType theErrorType,
        List<Double> theErrors, double theTolerance) {
        
        super(theInlierIndexes, theErrors, theTolerance);
        
        fundamentalMatrix = theFundamentalMatrix.copy();
        errorType = theErrorType;
    }
   
    public DenseMatrix getFundamentalMatrix() {
        return fundamentalMatrix;
    }

    /**
     * @return the errorType
     */
    public ErrorType getErrorType() {
        return errorType;
    }
    
    public void setFundamentalMatrix(DenseMatrix theFM) {
        fundamentalMatrix = theFM;
    }

    @Override
    public String toString() {
        String str = super.toString();
        StringBuilder sb = new StringBuilder(str);
        sb.append(" tol=").append(tolerance);
        if (fundamentalMatrix != null) {
            sb.append(" fm=" + fundamentalMatrix.toString());
        }
        return sb.toString();
    }
    
    
}
