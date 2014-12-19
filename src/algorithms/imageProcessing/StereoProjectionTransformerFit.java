package algorithms.imageProcessing;

import java.util.ArrayList;
import java.util.List;

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
    
    private List<Integer> inlierIndexes = new ArrayList<Integer>();

    public StereoProjectionTransformerFit(long theNumberOfMatches, 
        double theTolerance, double theAverageDistance, 
        double theStandardDeviationFromAverage) {
        
        nMatches = theNumberOfMatches;
        tolerance = theTolerance;
        meanDistance = theAverageDistance;
        stDevFromMean = theStandardDeviationFromAverage;
    }
    
    /**
     * @return the nMatches
     */
    public long getNMatches() {
        return nMatches;
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
            .append(" tolerance=").append(Double.toString(tolerance))
            .append(" meanDistFromModel=").append(Double.toString(meanDistance))
            .append(" stDevFromMean=").append(Double.toString(stDevFromMean))
            ;
        
        return sb.toString();
    }
}
