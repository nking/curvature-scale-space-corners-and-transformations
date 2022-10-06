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

    //TODO: change type to double[][] when have finished changes to RANSACSolver2 and ORB2 and classes that use the previous versions
    private DenseMatrix fundamentalMatrix;

    public EpipolarTransformationFit(DenseMatrix theFundamentalMatrix,
        List<Integer> theInlierIndexes, ErrorType theErrorType,
        List<Double> theErrors, double theTolerance) {

        super(theInlierIndexes, theErrors, theTolerance);

        fundamentalMatrix = theFundamentalMatrix.copy();
        errorType = theErrorType;
    }

    public EpipolarTransformationFit(double[][] theFundamentalMatrix,
                                     List<Integer> theInlierIndexes, ErrorType theErrorType,
                                     List<Double> theErrors, double theTolerance) {

        super(theInlierIndexes, theErrors, theTolerance);

        fundamentalMatrix = new DenseMatrix(theFundamentalMatrix);
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
        sb.append("\n  tol=").append(tolerance);
        if (fundamentalMatrix != null) {
            sb.append("\n  fm=\n").append(_toString(fundamentalMatrix, "%.3e"));
        }
        return sb.toString();
    }
   
}
