package algorithms.imageProcessing;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.logging.Logger;

/**
 * wrapper for ClosedCurveCornerMatcher that invokes ClosedCurveCornerMatcher
 * and handles cases such as scales smaller than one.
 * 
 * @author nichole
 */
public class ClosedCurveCornerMatcherWrapper {
    
    private final Logger log = Logger.getLogger(this.getClass().getName());

    /**
     * the original list of corner regions given to instance, sorted if not already.
     */
    private final List<CornerRegion> corners1;

    /**
     * the original list corner regions given to instance, sorted if not already.
     */
    private final List<CornerRegion> corners2;
    
    private final IntensityFeatures features1;
    
    private final IntensityFeatures features2;

    private final int nMaxMatchable;

    private TransformationParameters solutionParameters = null;

    private double solutionCost = Double.MAX_VALUE;

    private final List<Integer> solutionMatchedCornerIndexes1;

    private final List<Integer> solutionMatchedCornerIndexes2;
    
    private final List<FeatureComparisonStat> solutionMatchedCompStats;

    private boolean solverHasFinished = false;

    private boolean debug = true;
    
    public ClosedCurveCornerMatcherWrapper(
        final IntensityFeatures features1,
        final IntensityFeatures features2,
        final List<CornerRegion> c1, final List<CornerRegion> c2,
        boolean cornersAreAlreadySorted) {

        this.corners1 = new ArrayList<CornerRegion>(c1);

        this.corners2 = new ArrayList<CornerRegion>(c2);
        
        this.features1 = features1;
        
        this.features2 = features2;

        if (!cornersAreAlreadySorted) {

            // sort by descending |k|
            Collections.sort(c1, new DescendingKComparator());

            Collections.sort(c2, new DescendingKComparator());
        }

        nMaxMatchable = Math.min(corners1.size(), corners2.size());

        solutionMatchedCornerIndexes1 = new ArrayList<Integer>();

        solutionMatchedCornerIndexes2 = new ArrayList<Integer>();
        
        solutionMatchedCompStats = new ArrayList<FeatureComparisonStat>();
    }
    
     /**
     * match corners from the first list to the second.  the best scale can
     * be retrieved with getSolutionParameters().  if a solution was found, 
     * returns true, else returns false.
     * @return
     */
    public boolean matchCorners() {

        if (solverHasFinished) {
            throw new IllegalArgumentException(
            "matchCorners cannot be invoked more than once");
        }

        assert(this.solutionMatchedCornerIndexes1.isEmpty());

        assert(this.solutionMatchedCornerIndexes2.isEmpty());

        solverHasFinished = true;

        boolean cornersAreAlreadySorted = true;
        
        ClosedCurveCornerMatcher mDefault = new ClosedCurveCornerMatcher(
            features1, features2, corners1, corners2, cornersAreAlreadySorted);

        boolean solved = mDefault.matchCorners() &&
            (mDefault.getSolvedParameters() != null);

        log.info("default order: solved=" + solved + " ambig=" + 
            mDefault.scaleIsPossiblyAmbiguous());
        
        if (solved && !mDefault.scaleIsPossiblyAmbiguous()) {

            setSolutionToDefault(mDefault);

            return true;
        }

        return false;
    }

    public int getNMaxMatchable() {
        return nMaxMatchable;
    }

    public List<Integer> getSolutionMatchedCornerIndexes1() {
        return solutionMatchedCornerIndexes1;
    }

    public List<Integer> getSolutionMatchedCornerIndexes2() {
        return solutionMatchedCornerIndexes2;
    }

    
    public TransformationParameters getSolvedParameters() {
        return solutionParameters;
    }

    /**
     * get the curvature scale space images shift between the
     * first set of corners and the second set.
     * @return
     */
    public double getSolvedCost() {
        return solutionCost;
    }

    private void setSolutionToDefault(ClosedCurveCornerMatcher matcher) {

        if (matcher == null) {
            return;
        }

        this.solutionCost = matcher.getSolvedCost();
        
        this.solutionParameters = matcher.getSolvedParameters();

        this.solutionMatchedCornerIndexes1.addAll(matcher.getSolutionMatchedCornerIndexes1());
        this.solutionMatchedCornerIndexes2.addAll(matcher.getSolutionMatchedCornerIndexes2());
        
        this.solutionMatchedCompStats.addAll(matcher.getSolutionMatchedCompStats());
    }

    /**
     * @return the solutionMatchedCompStats
     */
    public List<FeatureComparisonStat> getSolutionMatchedCompStats() {
        return solutionMatchedCompStats;
    }
}
