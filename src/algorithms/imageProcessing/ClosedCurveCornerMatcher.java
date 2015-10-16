package algorithms.imageProcessing;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

/**
 *
 * @author nichole
 */
public class ClosedCurveCornerMatcher {
    
    protected Heap heap = null;

    /**
     * the costs calculated here are small fractions, so they need to be
     * multiplied by a large constant for use with the Fibonacci heap
     * which uses type long for its key (key is where cost is stored).
     * using 1E12 here
     */
    protected final static long heapKeyFactor = 1000000000000l;

    protected final List<CornerRegion> c1;

    protected final List<CornerRegion> c2;
    
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

    private final Logger log = Logger.getLogger(this.getClass().getName());

    private boolean solverHasFinished = false;

    private boolean hasBeenInitialized = false;

    private boolean solutionHasSomeScalesSmallerThanOne = false;
    
    private boolean strongestPeaksImplyScaleSmallerThanOne = false;
    
    //TODO: tune this
    private double tolerance = 10;

    public ClosedCurveCornerMatcher(final IntensityFeatures features1,
        final IntensityFeatures features2, final List<CornerRegion> corners1, 
        final List<CornerRegion> corners2, boolean cornersAreAlreadySorted) {
        
        c1 = new ArrayList<CornerRegion>(corners1.size());

        c2 = new ArrayList<CornerRegion>(corners2.size());
        
        this.features1 = features1;
        
        this.features2 = features2;
        
        initializeVariables(corners1, corners2);

        if (!cornersAreAlreadySorted) {

            // sort by descending |k|
            //Collections.sort(c1, new DescendingAbsKComparator());

            //Collections.sort(c2, new DescendingAbsKComparator());
        }

        initializeHeapNodes();
    }

    public boolean matchCorners() {
        
        throw new UnsupportedOperationException("not yet implemented");
    }
    
    /**
     * for each corner in c1, find the smallest matching SSD in c2 and return
     * the indexes w.r.t. c2.
     * Note that the code does not attempt unique (bipartite matching) because
     * the results are not used as the final match, rather a part of combinations
     * tried towards a solution.
     * runtime complexity is O(N_c1 * N_c2)
     * @return
     * @throws algorithms.imageProcessing.CornerRegion.CornerRegionDegneracyException 
     */
    protected FeatureComparisonStat[] getBestSSDC1ToC2() 
        throws CornerRegion.CornerRegionDegneracyException {
        
        FeatureComparisonStat[] stats = new FeatureComparisonStat[c1.size()];
        
        int dither = 1;
        
        FeatureMatcher featureMatcher = new FeatureMatcher();
        
        for (int i = 0; i < c1.size(); ++i) {
            
            CornerRegion region1 = c1.get(i);
            
            FeatureComparisonStat best = null;
            
            for (int j = 0; j < c2.size(); ++j) {
                
                CornerRegion region2 = c2.get(j);
             
                FeatureComparisonStat compStat =
                    featureMatcher.ditherAndRotateForBestLocation(
                    features1, features2, region1, region2, dither);
                
                if (compStat == null) {
                    continue;
                }
                
                if (best == null) {
                    best = compStat;
                    continue;
                }
                
                if (compStat.getSumIntensitySqDiff() < compStat.getImg2PointIntensityErr()) {
                    if (compStat.getSumIntensitySqDiff() < best.getSumIntensitySqDiff()) {
                        best = compStat;
                    }
                }
            }
            
            stats[i] = best;
        }
        
        return stats;
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
     * initialize the starter solution nodes for the heap
     */
    private void initializeHeapNodes() {
        
        /*
        using a limit decided by runtime to decide between 
        the n*(n-1)/2  starter points and the n*n*(n-1)/2 starter points.
              for n = 5  -->  10 vs 50
              for n = 10 -->  45 vs 450
              for n = 15 -->  105 vs 1575
              for n = 25 -->  300 vs 7500
        Since the corners are already high gradient in intensities, can expect 
        that the feature descriptors are somewhat unique so an assumption that 
        the best match is correct for at least 2 out of n is not bad for n > 10.
        */
        
        int n = Math.min(c1.size(), c2.size());
        
        if (n < 10) {
            initHeapNodes2();
        } else {
            initHeapNodes1();
        }
    }
    
    /**
     * create starter solution nodes for the heap using a pattern of combinations
     * with runtime O(n*(n-1)/2).
     * 
     */
    private void initHeapNodes1() {
        
        /*
        If one knew that that best SSD match of a point in curve1 to
               point in curve2 were true for at least 2 points in the curves,
               then one could use a pattern of solution starter points of
                  pt 1 = curve1[0] w/ best SSD match in curve2
                  pt 2 = curve1[1] w/ best SSD match in curve2
                  written as (pt1, pt2) for one solution starter
                             (pt1, pt3) for another solution starter
                             (pt1, pt4)  ""
                             (pt2, pt3)  ""
                             (pt2, pt4)  ""
                             (pt3, pt4)  ""
                             which is n!/(k!(n-k)!) since k is always 2, can rewrite as n*(n-1)/2.
                                 for a curve with 4 corners, the heap would have 6 solution starter nodes
        */
    }
    
    /**
     * create starter solution nodes for the heap using a pattern of combinations
     * with runtime O(n*n*(n-1)/2).
     */
    private void initHeapNodes2() {
        /*
        An improvement in completeness to initHeapNodes1() would be to try 
        all matches just for the first point in the two points needed in the solution.
                    The number of starter solutions would be  n*n*(n-1)/2
                    The chances of finding the correct solution are much higher.
                    It requires that only one point in the corners common to both curves
                    be a true match for it's best SSD match in curve 2.
                       for a curve with 4 corners, the heap would have 24 solution starter nodes
        */
    }
   
    protected CornerRegion findBestSSDWithinTolerance(CornerRegion corner1,
        double predictedX2, double predictedY2) 
        throws CornerRegion.CornerRegionDegneracyException {
            
        int x1 = corner1.getX()[corner1.getKMaxIdx()];
        int y1 = corner1.getY()[corner1.getKMaxIdx()];
        
        int dither = 1;
        
        FeatureMatcher featureMatcher = new FeatureMatcher();
        
        FeatureComparisonStat bestStat = null;
        
        CornerRegion best = null;
            
        for (int j = 0; j < c2.size(); ++j) {
                
            CornerRegion corner2 = c2.get(j);
            
            int x2 = corner2.getX()[corner2.getKMaxIdx()];
            int y2 = corner2.getY()[corner2.getKMaxIdx()];
            
            double dist = distance(x1, y1, x2, y2);
            
            if (dist > tolerance) {
                continue;
            }

            FeatureComparisonStat compStat =
                featureMatcher.ditherAndRotateForBestLocation(
                features1, features2, corner1, corner2, dither);

            if (compStat == null) {
                continue;
            }

            if (best == null) {
                bestStat = compStat;
                best = corner2;
                continue;
            }

            if (compStat.getSumIntensitySqDiff() < compStat.getImg2PointIntensityErr()) {
                if (compStat.getSumIntensitySqDiff() < bestStat.getSumIntensitySqDiff()) {
                    bestStat = compStat;
                    best = corner2;
                }
            }
        }
        
        return best;
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

    private double distance(double x1, double y1, double x2, double y2) {
        
        double diffX = x1 - x2;
        double diffY = y1 - y2;
        
        double dist = Math.sqrt(diffX*diffX + diffY*diffY);
        
        return dist;
    }
    
}
