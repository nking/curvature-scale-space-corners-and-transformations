package algorithms.imageProcessing;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.logging.Logger;

/**
 * class to use CSSContourMatcher to match the given contours lists from the
 * edges of two different images
 * with logic to retry with reverse ordered lists when needed for scales < 1..
 *
 * Based upon the algorithm contained in
 * <pre>
 * IEEE 'TRANSACTIONS ON PATTERN ANALYSIS AND MACHINE INTELLIGENCE, VOL. PAMI-8,
 * NO. 1. JANUARY 1986.  "Scale-Based Description and Recognition of Planar
 * Curves and Two-Dimensional Shapes" by FARZIN MOKHTARIAN AND ALAN MACKWORTH
 * </pre>
 *
 * This code accepts the contours from an edge in image 1 and contours from
 * an edge in image 2 and finds the best match of points between them calculating
 * scale and shift and cost as state of the match.
 *
 * @author nichole
 */
public final class CSSContourMatcherWrapper {

    private final Logger log = Logger.getLogger(this.getClass().getName());

    /**
     * the original list contours1 given to instance, sorted if not already.
     */
    private final List<CurvatureScaleSpaceContour> contours1;

    /**
     * the original list contours2 given to instance, sorted if not already.
     */
    private final List<CurvatureScaleSpaceContour> contours2;

    private final int nMaxMatchable;

    private double solutionScale = Double.MAX_VALUE;

    private double solutionShift = Double.MAX_VALUE;

    private double solutionCost = Double.MAX_VALUE;

    private final List<CurvatureScaleSpaceContour> solutionMatchedContours1;

    private final List<CurvatureScaleSpaceContour> solutionMatchedContours2;

    private boolean solverHasFinished = false;

    private boolean debug = true;

    /**
     * constructor taking required contour lists as arguments.  Note that the
     * contour lists should only be from one edge each from each image, that is
     * contours1 is from a single edge of image 1 and contours1 is from a
     * single edge of image 2.
     *
     * (Note, it should be possible to find shadows in an image too using this
     * on edges in the same image).
     *
     * constructor.  the creation of internal data structures in this method
     * has runtime complexity:
     * <pre>
     *   O(N) + O(N*lg_2(N)) + O(N_curves^2) + O(N_curves^3)
     *
     *       where each curve has a number of contours.
     *
     *       N = number of contours
     *       N_curve = number of encapsulating curves.  this number is smaller
     *                 than N.
     * </pre>
     * Note that the order of items in contours1 and contours2 will be altered
     * if alreadySorted is false.
     * @param contours1
     * @param contours2
     * @param contoursAreAlreadySorted
     */
    public CSSContourMatcherWrapper(
        final List<CurvatureScaleSpaceContour> contours1,
        final List<CurvatureScaleSpaceContour> contours2,
        boolean contoursAreAlreadySorted) {

        this.contours1 = new ArrayList<CurvatureScaleSpaceContour>(contours1);

        this.contours2 = new ArrayList<CurvatureScaleSpaceContour>(contours2);

        if (!contoursAreAlreadySorted) {

            Collections.sort(this.contours1, new DescendingSigmaComparator());

            Collections.sort(this.contours2, new DescendingSigmaComparator());
        }

        nMaxMatchable = Math.min(contours1.size(), contours2.size());

        solutionMatchedContours1 = new ArrayList<CurvatureScaleSpaceContour>();

        solutionMatchedContours2 = new ArrayList<CurvatureScaleSpaceContour>();
    }

     /**
     * match contours from the first list to the second.  the best scale and
     * shifts between the contour lists can be retrieved with getSolvedScale()
     * and getSolvedShift().  if a solution was found, returns true, else
     * returns false.
     * @return
     */
    public boolean matchContours() {

        if (solverHasFinished) {
            throw new IllegalArgumentException(
            "matchContours cannot be invoked more than once");
        }

        assert(this.solutionMatchedContours1.isEmpty());

        assert(this.solutionMatchedContours2.isEmpty());

        solverHasFinished = true;

        boolean contoursAreAlreadySorted = true;

        CSSContourMatcher mDefault = new CSSContourMatcher(contours1, contours2,
            contoursAreAlreadySorted);

        // ------- invoke reverse if needed and analyze all results -----

        boolean solved = mDefault.matchContours() &&
            (mDefault.getSolvedScale() < Double.MAX_VALUE);

        log.info("default order: solved=" + solved + " ambig=" + 
            mDefault.scaleIsPossiblyAmbiguous() + " possibly scl < 1=" + 
            mDefault.strongestPeaksImplyScaleSmallerThanOne());
        
        if (solved && !mDefault.scaleIsPossiblyAmbiguous()
            && !mDefault.strongestPeaksImplyScaleSmallerThanOne()) {

            setSolutionToDefault(mDefault);

            return true;
        }

        // ----- possibly ambiguous default solution or no default solution ----

        CSSContourMatcher mReverse = new CSSContourMatcher(contours2, contours1,
            contoursAreAlreadySorted);

        boolean solvedReverse = mReverse.matchContours() &&
            (mReverse.getSolvedScale() < Double.MAX_VALUE);

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

        double solutionScaleRevRev = 1. / mReverse.getSolvedScale();
        double solutionShiftRevRev = 1. - mReverse.getSolvedShift();

        if (!solved) {

            // reversed is solution
            this.solutionCost = mReverse.getSolvedCost();
            this.solutionScale = solutionScaleRevRev;
            this.solutionShift = solutionShiftRevRev;

            this.solutionMatchedContours1.addAll(mReverse.getSolutionMatchedContours2());
            this.solutionMatchedContours2.addAll(mReverse.getSolutionMatchedContours1());

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
            this.solutionScale = solutionScaleRevRev;
            this.solutionShift = solutionShiftRevRev;

            this.solutionMatchedContours1.addAll(mReverse.getSolutionMatchedContours2());
            this.solutionMatchedContours2.addAll(mReverse.getSolutionMatchedContours1());

            return true;
        }
    }

    public int getNMaxMatchable() {
        return nMaxMatchable;
    }

    public List<CurvatureScaleSpaceContour> getSolutionMatchedContours1() {
        return solutionMatchedContours1;
    }

    public List<CurvatureScaleSpaceContour> getSolutionMatchedContours2() {
        return solutionMatchedContours2;
    }

    /**
     * get the curvature scale space images factor of scale between the
     * first set of contours and the second set.
     * @return
     */
    public double getSolvedScale() {
        return solutionScale;
    }

    /**
     * get the curvature scale space images shift between the
     * first set of contours and the second set.
     * @return
     */
    public double getSolvedShift() {
        return solutionShift;
    }

    /**
     * get the curvature scale space images shift between the
     * first set of contours and the second set.
     * @return
     */
    public double getSolvedCost() {
        return solutionCost;
    }

    private void setSolutionToDefault(CSSContourMatcher matcher) {

        if (matcher == null) {
            return;
        }

        this.solutionCost = matcher.getSolvedCost();
        this.solutionScale = matcher.getSolvedScale();
        this.solutionShift = matcher.getSolvedShift();

        this.solutionMatchedContours1.addAll(matcher.getSolutionMatchedContours1());
        this.solutionMatchedContours2.addAll(matcher.getSolutionMatchedContours2());
    }

}
