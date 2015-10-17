package algorithms.imageProcessing;

import algorithms.compGeometry.NearestPoints;
import algorithms.util.PairIntArray;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.logging.Logger;

/**
 * A quick algorithm for matching points around ideal closed curves in
 * two images.  The algorithm uses point ordering around the curve and
 * needs the same number of corners in both curves.
 * 
 * The runtime complexity is O(N_corners).
 * 
 * It can be used and followed by ClosedCurveCornerMatcher if it fails.
 * 
 * @author nichole
 */
public class ClosedCurveCornerMatcher0 {
    
    protected final List<CornerRegion> c1;

    protected final List<CornerRegion> c2;
    
    protected final NearestPoints np;
    
    protected final IntensityFeatures features1;
    
    protected final IntensityFeatures features2;

    private float kMin1;
    private float kMax1;
    private float kMin2;
    private float kMax2;

    private TransformationParameters solutionParameters = null;

    private double solutionCost = Double.MAX_VALUE;

    private List<CornerRegion> solutionMatchedCorners1 = null;

    private List<CornerRegion> solutionMatchedCorners2 = null;
    
    private List<FeatureComparisonStat> solutionMatchedCompStats = null;

    private final Logger log = Logger.getLogger(this.getClass().getName());

    private boolean solverHasFinished = false;

    private boolean hasBeenInitialized = false;

    private boolean solutionHasSomeScalesSmallerThanOne = false;
    
    private boolean strongestPeaksImplyScaleSmallerThanOne = false;
    
    //TODO: tune this
    private int tolerance = 10;

    public ClosedCurveCornerMatcher0(final IntensityFeatures features1,
        final IntensityFeatures features2, final List<CornerRegion> corners1, 
        final List<CornerRegion> corners2, boolean cornersAreAlreadySorted) {
        
        c1 = new ArrayList<CornerRegion>(corners1.size());

        c2 = new ArrayList<CornerRegion>(corners2.size());
        
        if (c1.size() != c2.size()) {
            throw new IllegalArgumentException("c1 and c2 must be same size");
        }
        
        this.features1 = features1;
        
        this.features2 = features2;
        
        initializeVariables(corners1, corners2);

        if (!cornersAreAlreadySorted) {

            // sort by descending |k|
            Collections.sort(c1, new DescendingKComparator());

            Collections.sort(c2, new DescendingKComparator());
        }
        
        int[] xC2 = new int[c2.size()];
        int[] yC2 = new int[c2.size()];
        for (int i = 0; i < c2.size(); ++i) {
            CornerRegion cr = c2.get(i);
            xC2[i] = cr.getX()[cr.getKMaxIdx()];
            yC2[i] = cr.getY()[cr.getKMaxIdx()];
        }
        np = new NearestPoints(xC2, yC2);
    }

    public boolean matchCorners() {
        
        if (solverHasFinished) {
            throw new IllegalStateException(
            "matchContours cannot be invoked more than once");
        }

        solutionParameters = null;
        solutionCost = Double.MAX_VALUE;

        boolean solved = solve();

        return solved;        
    }

    private void initializeVariables(List<CornerRegion> corners1, 
        List<CornerRegion> corners2) {
        
        if (hasBeenInitialized) {
            return;
        }

        c1.addAll(corners1);
        c2.addAll(corners2);

        float minK = Float.MAX_VALUE;
        float maxK = Float.MIN_VALUE;
        for (int i = 0; i < c1.size(); i++) {
            CornerRegion cr = c1.get(i);
            float k = Math.abs(cr.k[cr.getKMaxIdx()]);
            if (k < minK) {
                minK = k;
            }
            if (k > maxK) {
                maxK = k;
            }
        }
        kMin1 = minK;
        kMax1 = maxK;

        minK = Float.MAX_VALUE;
        maxK = Float.MIN_VALUE;
        for (int i = 0; i < c2.size(); i++) {
            CornerRegion cr = c2.get(i);
            float k = Math.abs(cr.k[cr.getKMaxIdx()]);
            if (k < minK) {
                minK = k;
            }
            if (k > maxK) {
                maxK = k;
            }
        }
        kMin2 = minK;
        kMax2 = maxK;

        hasBeenInitialized = true;
    }
     
    /**
     * find solution using corners paired by curve order.
     * The algorithm requires c1.size() == c2.size().
     * with runtime O(n).
     * 
     */
    private boolean solve() {
        
        if (c1.size() != c2.size()) {
            throw new IllegalStateException("c1 and c2 must be same size");
        }
        
        float[] outputScaleRotTransXYStDev = new float[4];
        
        MatchedPointsTransformationCalculator 
            tc = new MatchedPointsTransformationCalculator();
                
        int delta = 0;
        
        TransformationParameters bestParams = null;
        double bestCost = Double.MAX_VALUE;
        List<CornerRegion> bestCR1 = null;
        List<CornerRegion> bestCR2 = null;
        
        for (int i = 0; i < c1.size(); ++i) {
                                         
            while (delta < c1.size()) {
                           
                // populate lists by ordered with point pairs having SSD < error
                List<CornerRegion> cr1 = new ArrayList<CornerRegion>();
                List<CornerRegion> cr2 = new ArrayList<CornerRegion>();
                populateByStartingIndexes(cr1, cr2, i, i + delta);
                
                // the remaining points follow in order and are given to
                // the transformation calculator
                PairIntArray xy1 = new PairIntArray(); 
                PairIntArray xy2 = new PairIntArray();
                
                //fill in xy1 and xy2 
                populateWithCoordinates(xy1, xy2, cr1, cr2);
                                
                // consider weighing by k
                float[] weights = new float[c1.size()];
                Arrays.fill(weights, 1.f/(float)c1.size());
                
                TransformationParameters params = tc.calulateEuclidean(xy1, xy2,
                    weights, 0, 0, outputScaleRotTransXYStDev);
                
                double cost = evaluateCost(cr1, cr2, params);
                
                if (cost < bestCost) {
                    bestParams = params;
                    bestCost = cost;
                    bestCR1 = cr1;
                    bestCR2 = cr2;
                }
            }
        }
        
        solverHasFinished = true;
        
        if (bestCost == Double.MAX_VALUE) {
            return false;
        }
        
        solutionParameters = bestParams;
        
        solutionCost = bestCost;

        solutionMatchedCorners1 = bestCR1;

        solutionMatchedCorners2 = bestCR2;
        
        return true;
    }

    private double calculateCost(CornerRegion cornerRegion,
        double x2, double y2) {

        if (cornerRegion == null) {
            return tolerance;
        }
        
        double dist = distance(cornerRegion.x[cornerRegion.getKMaxIdx()], 
            cornerRegion.y[cornerRegion.getKMaxIdx()], x2, y2);
                
        return dist;
    }

    private double distance(double x1, double y1, double x2, double y2) {
        
        double diffX = x1 - x2;
        double diffY = y1 - y2;
        
        double dist = Math.sqrt(diffX*diffX + diffY*diffY);
        
        return dist;
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

    public boolean scaleIsPossiblyAmbiguous() {
        return solutionHasSomeScalesSmallerThanOne;
    }
    public boolean strongestPeaksImplyScaleSmallerThanOne() {
        return strongestPeaksImplyScaleSmallerThanOne;
    }

    public List<CornerRegion> getSolutionMatchedCorners1() {
        return solutionMatchedCorners1;
    }

    public List<CornerRegion> getSolutionMatchedCorners2() {
        return solutionMatchedCorners2;
    }
    
    public List<FeatureComparisonStat> getSolutionMatchedCompStats() {
        return solutionMatchedCompStats;
    }

    protected double[] applyTransformation(CornerRegion corner, 
        TransformationParameters params) {
        
        int x0 = corner.getX()[corner.getKMaxIdx()];
        
        int y0 = corner.getY()[corner.getKMaxIdx()];
        
        Transformer transformer = new Transformer();
        
        double[] xyT = transformer.applyTransformation(params, x0, y0);
        
        return xyT;
    }

    /**
     * populate lists by starting indexes, but only with pairs whose
     * SSD < error
     * @param cr1
     * @param cr2
     * @param i
     * @param i0 
     */
    private void populateByStartingIndexes(List<CornerRegion> cr1, 
        List<CornerRegion> cr2, int startingIndexCurve1, int startingIndexCurve2) {
        throw new UnsupportedOperationException("Not supported yet."); 
    }

    private void populateWithCoordinates(PairIntArray xy1, PairIntArray xy2, 
        List<CornerRegion> cr1, List<CornerRegion> cr2) {
        throw new UnsupportedOperationException("Not supported yet."); 
    }

    /**
     * evaluate cost, including amounts for unmatched points
     * @param cr1
     * @param cr2
     * @param params
     * @return 
     */
    private double evaluateCost(List<CornerRegion> cr1, List<CornerRegion> cr2, 
        TransformationParameters params) {
        
        double unmatchedCost = (c1.size() - cr1.size()) * tolerance;
        
        throw new UnsupportedOperationException("Not supported yet."); 
    }
    
}
