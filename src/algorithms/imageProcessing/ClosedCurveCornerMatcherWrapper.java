package algorithms.imageProcessing;

import java.util.ArrayList;
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

    private final List<CornerRegion> solutionMatchedCorners1;

    private final List<CornerRegion> solutionMatchedCorners2;

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

            //TODO: sort by |k|
        }

        nMaxMatchable = Math.min(corners1.size(), corners2.size());

        solutionMatchedCorners1 = new ArrayList<CornerRegion>();

        solutionMatchedCorners2 = new ArrayList<CornerRegion>();
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

        assert(this.solutionMatchedCorners1.isEmpty());

        assert(this.solutionMatchedCorners2.isEmpty());

        solverHasFinished = true;

        boolean cornersAreAlreadySorted = true;

        ClosedCurveCornerMatcher mDefault = new ClosedCurveCornerMatcher(
            features1, features2, corners1, corners2, cornersAreAlreadySorted);

        // ------- invoke reverse if needed and analyze all results -----

        boolean solved = mDefault.matchCorners() &&
            (mDefault.getSolvedParameters() != null);

        log.info("default order: solved=" + solved + " ambig=" + 
            mDefault.scaleIsPossiblyAmbiguous() + " possibly scl < 1=" + 
            mDefault.strongestPeaksImplyScaleSmallerThanOne());
        
        if (solved && !mDefault.scaleIsPossiblyAmbiguous()
            && !mDefault.strongestPeaksImplyScaleSmallerThanOne()) {

            setSolutionToDefault(mDefault);

            return true;
        }

        // ----- possibly ambiguous default solution or no default solution ----

        ClosedCurveCornerMatcher mReverse = new ClosedCurveCornerMatcher(
            features2, features1, corners2, corners1, cornersAreAlreadySorted);

        boolean solvedReverse = mReverse.matchCorners() &&
            (mReverse.getSolvedParameters() != null);

        log.info("reverse order: solved=" + solvedReverse + " ambig=" + 
            mReverse.scaleIsPossiblyAmbiguous() + " possibly scl < 1=" + 
            mReverse.strongestPeaksImplyScaleSmallerThanOne());
        
        if (!solvedReverse && solved) {

            // take the possibly ambiguous default solution
            setSolutionToDefault(mDefault);

            return true;

        } else if (!solvedReverse) {

            return false;
        }

        MatchedPointsTransformationCalculator tc = new
            MatchedPointsTransformationCalculator();
        
        TransformationParameters solutionParametersRevRev = 
            (mReverse.getSolvedParameters() == null) ? null :
            tc.swapReferenceFrames(mReverse.getSolvedParameters().copy());
            
        if (!solved) {

            // reversed is solution
            this.solutionCost = mReverse.getSolvedCost();
            this.solutionParameters = solutionParametersRevRev;

            this.solutionMatchedCorners1.addAll(mReverse.getSolutionMatchedCorners2());
            this.solutionMatchedCorners2.addAll(mReverse.getSolutionMatchedCorners1());

            return true;
        }

        // ----- compare default solution to reversed solution ---------

        //TODO: normalize by some factor to prefer more matches?

        double costDefault = mDefault.getSolvedCost();

        double costReversed = mReverse.getSolvedCost();

        if (costDefault < costReversed) {
            // default is solution

            setSolutionToDefault(mDefault);

            return true;

        } else {
            // reversed is solution

            this.solutionCost = mReverse.getSolvedCost();
            this.solutionParameters = solutionParametersRevRev;

            this.solutionMatchedCorners1.addAll(mReverse.getSolutionMatchedCorners2());
            this.solutionMatchedCorners2.addAll(mReverse.getSolutionMatchedCorners1());

            return true;
        }
    }

    public int getNMaxMatchable() {
        return nMaxMatchable;
    }

    public List<CornerRegion> getSolutionMatchedContours1() {
        return solutionMatchedCorners1;
    }

    public List<CornerRegion> getSolutionMatchedContours2() {
        return solutionMatchedCorners2;
    }

    
    public TransformationParameters getSolvedParameters() {
        return solutionParameters;
    }

    /**
     * get the curvature scale space images shift between the
     * first set of contours and the second set.
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

        this.solutionMatchedCorners1.addAll(matcher.getSolutionMatchedCorners1());
        this.solutionMatchedCorners2.addAll(matcher.getSolutionMatchedCorners2());
    }
}
