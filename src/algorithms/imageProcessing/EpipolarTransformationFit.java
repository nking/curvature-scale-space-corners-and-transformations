package algorithms.imageProcessing;

import algorithms.imageProcessing.matching.ErrorType;
import algorithms.misc.MiscMath;
import java.util.ArrayList;
import java.util.List;
import org.ejml.simple.SimpleMatrix;

/**
 *
 * @author nichole
 */
public class EpipolarTransformationFit {
    
    private final ErrorType errorType;
    
    /**
     * when value is not Long.MIN_VALUE, it holds the maximum number of points
     * that were possible to match.  This may not always be populated, so 
     * check value before using.
     */
    private long nMaxMatchable = Long.MIN_VALUE;
    
    private final double tolerance;
    
    private double meanError = Double.MAX_VALUE;
    
    private double stDevFromMean = Double.MAX_VALUE;
    
    private final SimpleMatrix fundamentalMatrix;
    
    private final List<Integer> inlierIndexes;
    
    private final List<Double> errors;

    public EpipolarTransformationFit(SimpleMatrix theFundamentalMatrix,
        List<Integer> theInlierIndexes, ErrorType theErrorType,
        List<Double> theErrors, double theTolerance) {
        
        fundamentalMatrix = theFundamentalMatrix.copy();
        tolerance = theTolerance;
        inlierIndexes = new ArrayList<Integer>(theInlierIndexes);
        errorType = theErrorType;
        errors = new ArrayList<Double>(theErrors);
    }
    
    public void calculateErrorStatistics() {
        
        if (meanError != Double.MAX_VALUE) {
            // has already been calculated
            return;
        }
        
        double[] avgAndStDv = MiscMath.getAvgAndStDev(errors);
        
        meanError = avgAndStDv[0];
        stDevFromMean = avgAndStDv[1];
    }
    
    public SimpleMatrix getFundamentalMatrix() {
        return fundamentalMatrix;
    }
    
    /**
     * @return the nMatches
     */
    public long getNMatches() {
        return inlierIndexes.size();
    }
    
    /**
     * @param theNumberOfPossibleMatches
     */
    public void setNMaxMatchable(long theNumberOfPossibleMatches) {
        nMaxMatchable = theNumberOfPossibleMatches;
    }
    
    /**
     * @return the nMaxMatchable
     */
    public long getNMaxMatchable() {
        return nMaxMatchable;
    }

    /**
     * @return the tolerance
     */
    public double getTolerance() {
        return tolerance;
    }

    public double getMeanError() {
        return meanError;
    }

    /**
     * @return the stDevFromAvg
     */
    public double getStDevFromMean() {
        return stDevFromMean;
    }
    
    /**
     * compare to other by the number of inliers, else if tie, compares
     * mean of errors, else if tie, compares mean of standard deviation of mean
     * of errors, else returns false.
     * @param other
     * @return 
     */
    public boolean isBetter(EpipolarTransformationFit other) {
        if (other == null) {
            return true;
        }
        if (this.inlierIndexes.size() > other.inlierIndexes.size()) {
            return true;
        } else if (this.inlierIndexes.size() < other.inlierIndexes.size()) {
            return false;
        }
        calculateErrorStatistics();
        other.calculateErrorStatistics();
        if (meanError < other.meanError) {
            return true;
        } else if (meanError > other.meanError) {
            return false;
        }
        if (stDevFromMean < other.stDevFromMean) {
            return true;
        } else if (stDevFromMean > other.stDevFromMean) {
            return false;
        }
        
        // these are equivalent
        return false;
    }

    public List<Integer> getInlierIndexes() {
        
        return inlierIndexes;
    }
    
    @Override
    public String toString() {
        
        StringBuilder sb = new StringBuilder();
        
        sb.append("nMatchedPoints=").append(Long.toString(inlierIndexes.size()))
            .append(" nMaxMatchable=").append(Long.toString(nMaxMatchable))
            .append(" tolerance=").append(Double.toString(tolerance))
            .append(" meanDistFromModel=").append(Double.toString(meanError))
            .append(" stDevFromMean=").append(Double.toString(stDevFromMean))
            ;
        
        return sb.toString();
    }

    /**
     * @return the errorType
     */
    public ErrorType getErrorType() {
        return errorType;
    }

    /**
     * @return the errors
     */
    public List<Double> getErrors() {
        return errors;
    }
}
