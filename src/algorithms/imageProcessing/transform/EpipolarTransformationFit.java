package algorithms.imageProcessing.transform;

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

    private float effectiveTolerance = -1;
    private float effectiveToleranceStdv = -1;
    
    public EpipolarTransformationFit(SimpleMatrix theFundamentalMatrix,
        List<Integer> theInlierIndexes, ErrorType theErrorType,
        List<Double> theErrors, double theTolerance) {
        
        super(theInlierIndexes, theErrors, theTolerance);
        
        fundamentalMatrix = theFundamentalMatrix.copy();
        errorType = theErrorType;
    }
    
    public float getEffectiveTolerance() {
        return effectiveTolerance;
    }
    
    public float getEffectiveToleranceStandardDeviation() {
        return effectiveToleranceStdv;
    }
    
    /**
     * set the tolerance of the resulting error of tolerance size offsets in the
     * points that the fundamental matrix is derived from, where the
     * errors were derived using Sampson's error distances from respective
     * epipolar lines.
     * Note that this is a small sampling of the projection of tolerance
     * into the epipolar frame.  If the user has the camera matrix,
     * tolerance could be used precisely with each point.
     * 
     * Without the camera matrix, the effectiveTolerance is useful in 
     * thresholds for proximity to epipolar lines.
     * 
     * @param tol
     * @param effectiveTol
     * @param effectiveTolStdev 
     */
    public void setTolerance(float effectiveTol, float effectiveTolStdev) {
        this.effectiveTolerance = effectiveTol;
        this.effectiveToleranceStdv = effectiveTolStdev;
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
