package algorithms.imageProcessing.transform;

import algorithms.misc.MiscMath;
import java.util.ArrayList;
import java.util.List;
import no.uib.cipr.matrix.DenseMatrix;

/**
 *
 * @author nichole
 */
public abstract class AbstractTransformationFit implements ITransformationFit {
    
    /**
     * when value is not Long.MIN_VALUE, it holds the maximum number of points
     * that were possible to match.  This may not always be populated, so
     * check value before using.
     */
    protected long nMaxMatchable = Long.MIN_VALUE;
    
    protected final double tolerance;
    
    protected double meanError = Double.MAX_VALUE;
    
    protected double stDevFromMean = Double.MAX_VALUE;
    
    protected final List<Integer> inlierIndexes;
    
    protected final List<Double> errors;

    public AbstractTransformationFit(List<Integer> theInlierIndexes, 
        List<Double> theErrors, double theTolerance) {
        
        tolerance = theTolerance;
        
        inlierIndexes = new ArrayList<Integer>(theInlierIndexes);
        
        errors = new ArrayList<Double>(theErrors);
    }

    @Override
    public void calculateErrorStatistics() {
        
        if (meanError != Double.MAX_VALUE) {
            // has already been calculated
            return;
        }
        
        double[] avgAndStDv = MiscMath.getAvgAndStDev(errors);
        
        meanError = avgAndStDv[0];
        
        stDevFromMean = avgAndStDv[1];
    }

    /**
     * @return the nMatches
     */
    @Override
    public long getNMatches() {
        return inlierIndexes.size();
    }

    /**
     * @param theNumberOfPossibleMatches
     */
    @Override
    public void setNMaxMatchable(long theNumberOfPossibleMatches) {
        nMaxMatchable = theNumberOfPossibleMatches;
    }

    /**
     * @return the nMaxMatchable
     */
    @Override
    public long getNMaxMatchable() {
        return nMaxMatchable;
    }

    /**
     * @return the tolerance
     */
    @Override
    public double getTolerance() {
        return tolerance;
    }

    @Override
    public double getMeanError() {
        return meanError;
    }

    /**
     * @return the stDevFromAvg
     */
    @Override
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
    @Override
    public boolean isBetter(ITransformationFit other) {
        
        if (other == null || !(other instanceof ITransformationFit)) {
            return true;
        }
        
        if (this.inlierIndexes.size() > other.getInlierIndexes().size()) {
            return true;
        } else if (this.inlierIndexes.size() < other.getInlierIndexes().size()) {
            return false;
        }
        
        calculateErrorStatistics();
        
        other.calculateErrorStatistics();
        
        if (meanError < other.getMeanError()) {
            return true;
        } else if (meanError > other.getMeanError()) {
            return false;
        }
        
        if (stDevFromMean < other.getStDevFromMean()) {
            return true;
        } else if (stDevFromMean > other.getStDevFromMean()) {
            return false;
        }
        
        // these are equivalent
        return false;
    }
    
    /*
     * compare to other by the number of inliers, else if tie, compares
     * mean of errors, else if tie, compares mean of standard deviation of mean
     * of errors, else returns false.
     * @param other
     * @return
     */
    @Override
    public boolean isBetterByCost(ITransformationFit other) {
        
        if (other == null || !(other instanceof ITransformationFit)) {
            return true;
        }
        
        if (this.inlierIndexes.size() > other.getInlierIndexes().size()) {
            return true;
        } else if (this.inlierIndexes.size() < other.getInlierIndexes().size()) {
            return false;
        }
        
        double sumError = 0;
        for (Double error : errors) {
            sumError += error.doubleValue();
        }
        
        double sumErrorOther = 0;
        for (Double error : other.getErrors()) {
            sumErrorOther += error.doubleValue();
        }
        
        if (sumError < sumErrorOther) {
            return true;
        } else if (sumError > sumErrorOther) {
            return false;
        }
        
        if (this.inlierIndexes.size() > other.getInlierIndexes().size()) {
            return true;
        } else if (this.inlierIndexes.size() < other.getInlierIndexes().size()) {
            return false;
        }
        
        calculateErrorStatistics();
        
        other.calculateErrorStatistics();
        if (stDevFromMean < other.getStDevFromMean()) {
            return true;
        } else if (stDevFromMean > other.getStDevFromMean()) {
            return false;
        }
        
        // these are equivalent
        return false;
    }

    @Override
    public List<Integer> getInlierIndexes() {
        return inlierIndexes;
    }

    public String _toString(DenseMatrix a, String decimalFormat) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < a.numRows(); ++i) {
            for (int j = 0; j < a.numColumns(); ++j) {
                sb.append(String.format(decimalFormat, a.get(i, j)));
                if (j < (a.numColumns()-1)) {
                    sb.append(", ");
                }
            }
            sb.append("\n");
        }
        return sb.toString();
    }
    
    @Override
    public String toString() {
        
        StringBuilder sb = new StringBuilder();
        
        sb.append("nMatchedPoints=").append(Long.toString(inlierIndexes.size()))
            .append(" nMaxMatchable=").append(Long.toString(nMaxMatchable))
            .append(" tolerance=").append(Double.toString(tolerance))
            .append(" meanDistFromModel=")
                .append(String.format("%.3e", meanError))
            .append(" stDevFromMean=")
                .append(String.format("%.3e", stDevFromMean));
        
        return sb.toString();
    }

    /**
     * @return the errors
     */
    @Override
    public List<Double> getErrors() {
        return errors;
    }
    
}
