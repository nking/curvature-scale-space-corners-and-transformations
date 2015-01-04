package algorithms.imageProcessing;

import java.util.ArrayList;
import java.util.List;
import org.ejml.simple.SimpleMatrix;

/**
 *
 * @author nichole
 */
public class StereoProjectionTransformerFit {
    
    /**
     * number of matches within tolerance
     */
    private final long nMatches;
    
    /**
     * when value is not Long.MIN_VALUE, it holds the maximum number of points
     * that were possible to match.  This may not always be populated, so 
     * check value before using.
     */
    private long nMaxMatchable = Long.MIN_VALUE;
    
    /**
     * tolerance in distance of a point from the epipolar line it belongs to
     */
    private final double tolerance;
    
    /**
     * the average distance of the points from their epipolar lines
     */
    private final double meanDistance;
    
    /**
     * the standard deviation of the points from the average distance
     */
    private final double stDevFromMean;
    
    private final SimpleMatrix fundamentalMatrix;
    
    private List<Integer> inlierIndexes = new ArrayList<Integer>();

    public StereoProjectionTransformerFit(SimpleMatrix theFundamentalMatrix,
        long theNumberOfMatches, 
        double theTolerance, double theAverageDistance, 
        double theStandardDeviationFromAverage) {
        
        fundamentalMatrix = theFundamentalMatrix;
        nMatches = theNumberOfMatches;
        tolerance = theTolerance;
        meanDistance = theAverageDistance;
        stDevFromMean = theStandardDeviationFromAverage;
    }
    
    public SimpleMatrix getFundamentalMatrix() {
        return fundamentalMatrix;
    }
    
    /**
     * @return the nMatches
     */
    public long getNMatches() {
        return nMatches;
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

    /**
     * @return the avgDistance
     */
    public double getMeanDistance() {
        return meanDistance;
    }

    /**
     * @return the stDevFromAvg
     */
    public double getStDevFromMean() {
        return stDevFromMean;
    }

    /**
     * compare bestFit to this instance and return true if bestFit is better.
     * First compares the number of matches, then breaks ties with the
     * average distance, then breaks ties with the standard deviation from the
     * average, else returns false.
     * 
     * @param bestFit
     * @return 
     */
    public boolean otherIsBetter(StereoProjectionTransformerFit bestFit) {
        
        //TODO: consider converting this to compare
        
        if (bestFit == null) {
            return false;
        }
        
        if (bestFit.getNMatches() != nMatches) {
            return (bestFit.getNMatches() > nMatches);
        }
        
        // else nMatches are equal
        
        if (bestFit.getMeanDistance() != meanDistance) {
            return (bestFit.getMeanDistance() > meanDistance);
        }
        
        // else avgDistances are equal
        
        if (bestFit.getStDevFromMean() != stDevFromMean) {
            return (bestFit.getStDevFromMean() > stDevFromMean);
        }
        
        // else, all params being equal, prefer current
        return false;
    }

    public void setInlierIndexes(List<Integer> theInlierIndexes) {
        
        if (theInlierIndexes == null) {
            if (!inlierIndexes.isEmpty()) {
                inlierIndexes.clear();
            }
            return;
        }
        inlierIndexes = theInlierIndexes;
    }

    public List<Integer> getInlierIndexes() {
        
        return inlierIndexes;
    }
    
    @Override
    public String toString() {
        
        StringBuilder sb = new StringBuilder();
        
        sb.append("nMatchedPoints=").append(Long.toString(nMatches))
            .append(" nMaxMatchable=").append(Long.toString(nMaxMatchable))
            .append(" tolerance=").append(Double.toString(tolerance))
            .append(" meanDistFromModel=").append(Double.toString(meanDistance))
            .append(" stDevFromMean=").append(Double.toString(stDevFromMean))
            ;
        
        return sb.toString();
    }
}
