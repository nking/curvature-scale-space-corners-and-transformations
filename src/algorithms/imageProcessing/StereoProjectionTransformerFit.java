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
    private final double avgDistance;
    
    /**
     * the standard deviation of the points from the average distance
     */
    private final double stDevFromAvg;
    
    private List<Integer> outlierIndexes = new ArrayList<Integer>();

    public StereoProjectionTransformerFit(long theNumberOfMatches, 
        double theTolerance, double theAverageDistance, 
        double theStandardDeviationFromAverage) {
        
        nMatches = theNumberOfMatches;
        tolerance = theTolerance;
        avgDistance = theAverageDistance;
        stDevFromAvg = theStandardDeviationFromAverage;
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
    public double getAvgDistance() {
        return avgDistance;
    }

    /**
     * @return the stDevFromAvg
     */
    public double getStDevFromAvg() {
        return stDevFromAvg;
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
        
        if (bestFit.getAvgDistance() != avgDistance) {
            return (bestFit.getAvgDistance() > avgDistance);
        }
        
        // else avgDistances are equal
        
        if (bestFit.getStDevFromAvg() != stDevFromAvg) {
            return (bestFit.getStDevFromAvg() > stDevFromAvg);
        }
        
        // else, all params being equal, prefer current
        return false;
    }

    public void setOutlierIndexes(List<Integer> theOutlierIndexes) {
        
        if (theOutlierIndexes == null) {
            if (!outlierIndexes.isEmpty()) {
                outlierIndexes.clear();
            }
            return;
        }
        outlierIndexes = theOutlierIndexes;
    }

    public List<Integer> getOutlierIndexes() {
        
        return outlierIndexes;
    }
}
