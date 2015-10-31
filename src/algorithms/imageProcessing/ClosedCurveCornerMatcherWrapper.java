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

    protected TransformationPair2 solution = null;
    
    private double solutionCost = Double.MAX_VALUE;

    private boolean solverHasFinished = false;

    private boolean debug = true;
    
    public ClosedCurveCornerMatcherWrapper() {
    }
    
    private void resetDefaults() {
        solverHasFinished = false;
        solutionCost = Double.MAX_VALUE;
        solution = null;
    }
    
     /**
     * match corners from the first list to the second.  the best scale can
     * be retrieved with getSolutionParameters().  if a solution was found, 
     * returns true, else returns false.
     * @param features1
     * @param features2
     * @param corners1
     * @param corners2
     * @param cornersAreAlreadySorted
     * @param img1 image from which to extract descriptors for features1
     * @param img2 image from which to extract descriptors for features2
     * @return
     */
    public boolean matchCorners(final IntensityFeatures features1,
        final IntensityFeatures features2, final List<CornerRegion> corners1,
        final List<CornerRegion> corners2, boolean cornersAreAlreadySorted,
        GreyscaleImage img1, GreyscaleImage img2) {

        if (solverHasFinished) {
            resetDefaults();
        }
        
        List<CornerRegion> c1 = new ArrayList<CornerRegion>(corners1);

        List<CornerRegion> c2 = new ArrayList<CornerRegion>(corners2);
    
        if (!cornersAreAlreadySorted) {
            // sort by descending |k|
            Collections.sort(c1, new DescendingKComparator());

            Collections.sort(c2, new DescendingKComparator());
            
            cornersAreAlreadySorted = true;
        }

        solverHasFinished = true;
        
        ClosedCurveCornerMatcher mDefault = new ClosedCurveCornerMatcher();

        boolean solved = mDefault.matchCorners(features1, features2, corners1, 
            corners2, cornersAreAlreadySorted, img1, img2) &&
            (mDefault.getSolution() != null);

        log.fine("default order: solved=" + solved + " ambig=" + 
            mDefault.scaleIsPossiblyAmbiguous());
        
        if (solved && !mDefault.scaleIsPossiblyAmbiguous()) {

            solutionCost = mDefault.getSolvedCost();
            
            solution = mDefault.getSolution();

            return true;
        }

        return false;
    }

    public TransformationPair2 getSolution() {
        return solution;
    }
    
    /**
     * get the curvature scale space images shift between the
     * first set of corners and the second set.
     * @return
     */
    public double getSolvedCost() {
        return solutionCost;
    }

}
