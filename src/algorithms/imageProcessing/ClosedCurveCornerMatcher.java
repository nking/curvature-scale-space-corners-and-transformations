package algorithms.imageProcessing;

import algorithms.compGeometry.NearestPoints;
import algorithms.imageProcessing.util.MatrixUtil;
import algorithms.util.PairInt;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;
import thirdparty.HungarianAlgorithm;

/**
 *
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
    
    protected final NearestPoints np;
    
    protected final IntensityFeatures features1;
    
    protected final IntensityFeatures features2;
    
    protected final int dither = 3;
    
    protected final int rotationTolerance = 20;

    //TODO: tune this
    private int tolerance = 2;

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
        
    public ClosedCurveCornerMatcher(final IntensityFeatures features1,
        final IntensityFeatures features2, final List<CornerRegion> corners1, 
        final List<CornerRegion> corners2, boolean cornersAreAlreadySorted) {
        
        c1 = new ArrayList<CornerRegion>(corners1.size());

        c2 = new ArrayList<CornerRegion>(corners2.size());
        
        solutionMatchedCompStats = new ArrayList<FeatureComparisonStat>();
        
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

        initializeHeapNodes();
    }

    public boolean matchCorners() {
        
        if (solverHasFinished) {
            throw new IllegalStateException(
            "matchContours cannot be invoked more than once");
        }

        if (heap.getNumberOfNodes() == 0) {
            solverHasFinished = true;
            return false;
        }

        solutionParameters = null;
        solutionCost = Double.MAX_VALUE;

        //temporary look at evaluating entire solutions:
        HeapNode minCost = heap.extractMin();//solve();
     
        if (minCost == null) {
            return false;
        }

        //TODO:  may need to revise. minCost in solve is by SSD, so looking for best dist among the top 10 SSD
        List<HeapNode> top = new ArrayList<HeapNode>();
        top.add(minCost);
        HeapNode node = heap.extractMin();
        while ((node != null) && (top.size() < 11)) {
            int n = ((TransformationPair2)node.getData()).getNextCorner()
                .getMatchedCorners1().size();
            if (n > 2) {
                top.add(node);
            }
            node = heap.extractMin();
        }
        
        // combine similar results
        // for example, in one test dataset can see the correct solution is
        // at top[1], and top[0]
        
        for (int i = 1; i < top.size(); ++i) {
            TransformationPair2 tPair = ((TransformationPair2)top.get(i).getData());
            double dist = tPair.getCostAsDistance();
            int n2 = tPair.getNextCorner().getMatchedCorners1().size();
            
        }
        
        TransformationPair2 transformationPair = (TransformationPair2)minCost.getData();

        TransformationParameters params = transformationPair.getTransformationParameters();

        solutionParameters = params;
        solutionCost = ((double)minCost.getKey()/(double)heapKeyFactor);

        NextCorner nc = transformationPair.getNextCorner();

        solutionMatchedCorners1 = nc.getMatchedCorners1();

        solutionMatchedCorners2 = nc.getMatchedCorners2();
        
        solutionMatchedCompStats = nc.getMatchedFeatureComparisonStats();

        solverHasFinished = true;

        solutionHasSomeScalesSmallerThanOne = transformationPair.scaleIsPossiblyAmbiguous();

        return true;        
    }

    private CornersAndFeatureStat findBestMatchWithinTolerance(
        CornerRegion cornerCurve1, double[] predictedXYCurve2, 
        int rotationInDegrees, int rotationTolerance) {
        
        FeatureMatcher featureMatcher = new FeatureMatcher();

        FeatureComparisonStat bestStat = null;

        int bestIdx2 = -1;

        Set<Integer> c2Indexes = np.findNeighborIndexes(
            (int) Math.round(predictedXYCurve2[0]),
            (int) Math.round(predictedXYCurve2[1]), tolerance);

        for (Integer index : c2Indexes) {

            int idx2 = index.intValue();

            CornerRegion corner2 = c2.get(idx2);

            FeatureComparisonStat compStat = null;

            try {

                compStat = featureMatcher.ditherAndRotateForBestLocation(
                    features1, features2, cornerCurve1, corner2, dither,
                    rotationInDegrees, rotationTolerance);

            } catch (CornerRegion.CornerRegionDegneracyException ex) {

                log.severe(ex.getMessage());
            }

            if (compStat == null) {
                continue;
            }

            if (bestStat == null) {
                if (compStat.getSumIntensitySqDiff() < compStat.getImg2PointIntensityErr()) {
                    bestStat = compStat;
                    bestIdx2 = idx2;
                }
                continue;
            }

            if (compStat.getSumIntensitySqDiff() < compStat.getImg2PointIntensityErr()) {
                if (compStat.getSumIntensitySqDiff() < bestStat.getSumIntensitySqDiff()) {
                    bestStat = compStat;
                    bestIdx2 = idx2;
                }
            }           
        }
        
        if (bestIdx2 > -1) {
            CornersAndFeatureStat cfs = new CornersAndFeatureStat(cornerCurve1, 
                c2.get(bestIdx2), bestStat);
            cfs.idx2 = bestIdx2;
            return cfs;
        } else {
            return null;
        } 
    }
    
    private static class CornersAndFeatureStat {
        private final CornerRegion cr1;
        private final CornerRegion cr2;
        private int idx1 = -1;
        private int idx2 = -1;
        private final FeatureComparisonStat stat;
        public CornersAndFeatureStat(CornerRegion cornerRegion1,
            CornerRegion cornerRegion2, FeatureComparisonStat compStat) {
            cr1 = cornerRegion1;
            cr2 = cornerRegion2;
            stat = compStat;
        }
    }
    
    /**
     * for each corner in c1, find the smallest matching SSD in c2 and return
     * the indexes w.r.t. c2.
     * Note that the code does not attempt unique (bipartite matching) because
     * the results are not used as the final match, rather a part of combinations
     * tried towards a solution.
     * runtime complexity is O(N_c1 * N_c2)
     * @return
     */
    protected CornersAndFeatureStat[] getBestSSDC1ToC2() {
        
        // if no match, contains a null
        CornersAndFeatureStat[] indexes2 = new CornersAndFeatureStat[c1.size()];
                
        FeatureMatcher featureMatcher = new FeatureMatcher();
        
        for (int i = 0; i < c1.size(); ++i) {
            
            CornerRegion region1 = c1.get(i);
            
            FeatureComparisonStat best = null;
            
            int bestIdx2 = -1;
            
            for (int j = 0; j < c2.size(); ++j) {
                
                CornerRegion region2 = c2.get(j);
             
                FeatureComparisonStat compStat = null;
                
                try {
                
                    compStat = featureMatcher.ditherAndRotateForBestLocation(
                        features1, features2, region1, region2, dither);
                                    
                } catch (CornerRegion.CornerRegionDegneracyException ex) {
                    log.severe("**CONSIDER using more points in corner region");
                }
                
                if (compStat == null) {
                    continue;
                }
                
                if (best == null) {
                    if (compStat.getSumIntensitySqDiff() < compStat.getImg2PointIntensityErr()) {
                        best = compStat;
                        bestIdx2 = j;
                    }
                    continue;
                }
                
                if (compStat.getSumIntensitySqDiff() < compStat.getImg2PointIntensityErr()) {
                    if (compStat.getSumIntensitySqDiff() < best.getSumIntensitySqDiff()) {
                        best = compStat;
                        bestIdx2 = j;
                    }
                }
            }
            
            if (bestIdx2 > -1) {
                indexes2[i] = new CornersAndFeatureStat(c1.get(i), 
                    c2.get(bestIdx2), best);
            }
        }
        
        return indexes2;
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
        
        heap = new Heap();
        
        int n = Math.min(c1.size(), c2.size());
        
        /*if (n < 10) {
            initHeapNodes2();
        } else {*/
            initHeapNodes1();
        //}
    }
    
    /**
     * create starter solution nodes for the heap using a pattern of combinations
     * with runtime O(n*(n-1)/2).
     * 
     */
    private void initHeapNodes1() {
        
        // if no match, contains a null
        CornersAndFeatureStat[] indexes2 = getBestSSDC1ToC2();
        
log.info("getBestSSDC1ToC2: ");
for (int i = 0; i < indexes2.length; ++i) {
    CornersAndFeatureStat cfs = indexes2[i];
    if (cfs == null) {
        continue;
    }
    log.info("pre-search: SSD match of " + cfs.stat.getImg1Point().toString() 
        + " to " + cfs.stat.getImg2Point());
}
        
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
                     which is n!/(k!(n-k)!) 
                     since k is always 2, can rewrite it as n*(n-1)/2.
        for a curve with 4 corners, the heap would have 6 solution starter nodes
        */
        
        MatchedPointsTransformationCalculator 
            tc = new MatchedPointsTransformationCalculator();

        for (int i = 0; i < indexes2.length; ++i) {
                        
            CornerRegion cr1C1 = this.c1.get(i);
            
            if (indexes2[i] == null) {
                continue;
            }
            
            CornerRegion cr1C2 = indexes2[i].cr2;
            
            for (int j = (i + 1); j < c1.size(); ++j) {
                
                CornerRegion cr2C1 = this.c1.get(j);
                
                if (indexes2[j] == null) {
                    continue;
                }
                
                CornerRegion cr2C2 = indexes2[j].cr2;
                
                if (cr2C2.equals(cr1C2)) {
                    continue;
                }
                
                // temporarily, evaluating all corners for each starter solution:
                insertNodeTMP(tc, indexes2[i], cr1C1, cr1C2, indexes2[j], cr2C1, 
                    cr2C2);  
            }
        }
    }
    
    /**
     * create starter solution nodes for the heap using a pattern of combinations
     * with runtime O(n*n*(n-1)/2).
     */
    private void initHeapNodes2() {
        
        // if no match, has a null
        CornersAndFeatureStat[] indexes2 = getBestSSDC1ToC2();
        
        /*
        An improvement in completeness to initHeapNodes1() would be to try 
        all matches just for the first point in the two points needed in the solution.
            The number of starter solutions would be  n*n*(n-1)/2
            The chances of finding the correct solution are much higher.
            It requires that only one point in the corners common to both curves
            be a true match for it's best SSD match in curve 2.
       for a curve with 4 corners, the heap would have 24 solution starter nodes
        */
        
        FeatureMatcher featureMatcher = new FeatureMatcher();
                
        MatchedPointsTransformationCalculator 
            tc = new MatchedPointsTransformationCalculator();
                
        for (int idx1C1 = 0; idx1C1 < c1.size(); ++idx1C1) {
                        
            CornerRegion cr1C1 = this.c1.get(idx1C1);
            
            for (int idx1C2 = 0; idx1C2 < c2.size(); ++idx1C2) {
                
                CornerRegion cr1C2 = this.c2.get(idx1C2);
                
                FeatureComparisonStat compStat1 = null;
                try {                
                    compStat1 = featureMatcher.ditherAndRotateForBestLocation(
                        features1, features2, cr1C1, cr1C2, dither);
                } catch (CornerRegion.CornerRegionDegneracyException ex) {                    
                    log.severe(ex.getMessage());
                }
                
                if (compStat1 == null) {
                    continue;
                }
                
                if (compStat1.getSumIntensitySqDiff() >= compStat1.getImg2PointIntensityErr()) {
                    continue;
                }
                
                CornersAndFeatureStat cfs1 = new CornersAndFeatureStat(cr1C1, 
                    cr1C2, compStat1);
                
                for (int idx2C1 = (idx1C1 + 1); idx2C1 < c1.size(); ++idx2C1) {
                
                    CornerRegion cr2C1 = this.c1.get(idx2C1);
                    
                    if (indexes2[idx2C1] == null) {
                        continue;
                    }
                
                    CornerRegion cr2C2 = indexes2[idx2C1].cr2;
                    
                    if (cr2C2.equals(cr1C2)) {
                        continue;
                    }                    
                    
                    insertNode(tc, cfs1, cr1C1, cr1C2, indexes2[idx2C1], 
                        cr2C1, cr2C2);
                }
            }
        }
    }
    
    protected CornersAndFeatureStat findBestSSDWithinTolerance(CornerRegion corner1,
        double predictedX2, double predictedY2) {
        
        int x1 = corner1.getX()[corner1.getKMaxIdx()];
        int y1 = corner1.getY()[corner1.getKMaxIdx()];
                
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

            FeatureComparisonStat compStat = null;
            
            try {
                
                compStat = featureMatcher.ditherAndRotateForBestLocation(
                    features1, features2, corner1, corner2, dither);
                
            } catch (CornerRegion.CornerRegionDegneracyException ex) {
                
                log.severe(ex.getMessage());
            }

            if (compStat == null) {
                continue;
            }

            if (bestStat == null) {
                if (compStat.getSumIntensitySqDiff() < bestStat.getImg2PointIntensityErr()) {
                    bestStat = compStat;
                    best = corner2;
                }
                continue;
            }

            if (compStat.getSumIntensitySqDiff() < compStat.getImg2PointIntensityErr()) {
                if (compStat.getSumIntensitySqDiff() < bestStat.getSumIntensitySqDiff()) {
                    bestStat = compStat;
                    best = corner2;
                }
            }
        }
        
        if (best == null) {
            return null;
        } else {
            return new CornersAndFeatureStat(corner1, best, bestStat);
        }
        
    }

    private double calculateCost(CornerRegion cornerRegion, double x2, double y2) {

        if (cornerRegion == null) {
            
            return calculateUnMatchedCost();
        }
        
        double dist = distance(cornerRegion.x[cornerRegion.getKMaxIdx()], 
            cornerRegion.y[cornerRegion.getKMaxIdx()], x2, y2);
                
        return dist;
    }
    
    private double calculateCost(PairInt cornerPoint, double x2, double y2) {

        if (cornerPoint == null) {
            
            return calculateUnMatchedCost();
        }
        
        double dist = distance(cornerPoint.getX(), cornerPoint.getY(),
            x2, y2);
                
        return dist;
    }
    
    private double calculateUnMatchedCost() {
        //TODO: refine this.  needs to be a large number but < infinity
        
        int n = Math.max(c1.size(), c2.size());
        
        return tolerance*n;
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

    public List<CornerRegion> getSolutionMatchedCorners1() {
        return solutionMatchedCorners1;
    }

    public List<CornerRegion> getSolutionMatchedCorners2() {
        return solutionMatchedCorners2;
    }
    
    public List<FeatureComparisonStat> getSolutionMatchedCompStats() {
        return solutionMatchedCompStats;
    }

    private void insertNode(MatchedPointsTransformationCalculator tc, 
        CornersAndFeatureStat stat1, CornerRegion cr1C1, CornerRegion cr1C2, 
        CornersAndFeatureStat stat2, CornerRegion cr2C1, CornerRegion cr2C2) {
       
        // use dither corrected locations:
        final int x1C1 = stat1.stat.getImg1Point().getX();
        final int y1C1 = stat1.stat.getImg1Point().getY();
        final int x1C2 = stat1.stat.getImg2Point().getX();
        final int y1C2 = stat1.stat.getImg2Point().getY();
        
        final int x2C1 = stat2.stat.getImg1Point().getX();
        final int y2C1 = stat2.stat.getImg1Point().getY();
        final int x2C2 = stat2.stat.getImg2Point().getX();
        final int y2C2 = stat2.stat.getImg2Point().getY();
                
        // with 2 points in both image, calc transformation
        TransformationParameters params = tc.calulateEuclidean(
            x1C1, y1C1, x2C1, y2C1, 
            x1C2, y1C2, x2C2, y2C2, 0, 0);

        TransformationPair2 transformationPair = 
            new TransformationPair2(cr1C1, cr1C2, cr2C1, cr2C2);
        
        transformationPair.setTransformationParameters(params);

        NextCorner nc = new NextCorner(c1);
        nc.addMatchedCorners(cr1C1, cr1C2, stat1.stat);
        nc.addMatchedCorners(cr2C1, cr2C2, stat2.stat);

        transformationPair.setNextCorner(nc);

        // choose a third point in c1
        CornerRegion c3C1 = nc.findStrongestRemainingCorner();

        // apply the transformation to it
        double[] xy3C2 = applyTransformation(c3C1, params);

        // find the best matching SSD within tolerance of predicted, and 
        // restrict to existing rotation
        CornersAndFeatureStat cfs = findBestMatchWithinTolerance(c3C1, xy3C2,
            Math.round(params.getRotationInDegrees()), rotationTolerance);
        
        double cost;
        double costSSD; 
        
        if (cfs != null) {

            nc.addMatchedCorners(c3C1, cfs.cr2, cfs.stat);
        
            cost = calculateCost(cfs.stat.getImg2Point(), xy3C2[0], xy3C2[1]);
            
            // TODO: may need to add cost as SSD in quadrature
            costSSD = cfs.stat.getSumIntensitySqDiff();
            
            transformationPair.addToCostAsDistance(cost);
            
            transformationPair.addToCostAsSSD(costSSD);
            
        } else {
                
            cost = calculateUnMatchedCost();
            
            costSSD = 10000;
        }        
             
        // put node in heap
        long costL = (long)(costSSD * heapKeyFactor);

        HeapNode node = new HeapNode(costL);

        node.setData(transformationPair);

        heap.insert(node);
    }
    
    private void insertNodeTMP(MatchedPointsTransformationCalculator tc, 
        CornersAndFeatureStat stat1, CornerRegion cr1C1, CornerRegion cr1C2, 
        CornersAndFeatureStat stat2, CornerRegion cr2C1, CornerRegion cr2C2) {
       
        // use dither corrected locations:
        final int x1C1 = stat1.stat.getImg1Point().getX();
        final int y1C1 = stat1.stat.getImg1Point().getY();
        final int x1C2 = stat1.stat.getImg2Point().getX();
        final int y1C2 = stat1.stat.getImg2Point().getY();
        
        final int x2C1 = stat2.stat.getImg1Point().getX();
        final int y2C1 = stat2.stat.getImg1Point().getY();
        final int x2C2 = stat2.stat.getImg2Point().getX();
        final int y2C2 = stat2.stat.getImg2Point().getY();
                
        // with 2 points in both image, calc transformation
        TransformationParameters params = tc.calulateEuclidean(
            x1C1, y1C1, x2C1, y2C1, 
            x1C2, y1C2, x2C2, y2C2, 0, 0);

        TransformationPair2 transformationPair = 
            new TransformationPair2(cr1C1, cr1C2, cr2C1, cr2C2);
        
        transformationPair.setTransformationParameters(params);
        
        NextCorner nc = new NextCorner(c1);
        nc.addMatchedCorners(cr1C1, cr1C2, stat1.stat);
        nc.addMatchedCorners(cr2C1, cr2C2, stat2.stat);

        transformationPair.setNextCorner(nc);
        
        Map<PairInt, CornersAndFeatureStat> indexesCFSMap = new 
            HashMap<PairInt, CornersAndFeatureStat>();
        
        Map<PairInt, Double> indexesDistMap = new HashMap<PairInt, Double>();
        
        float[][] cost = new float[c1.size()][c2.size()];
       
        // use bipartite matching on the remaining points.
                
        for (int i = 0; i < c1.size(); ++i) {
            
            cost[i] = new float[c2.size()];
            Arrays.fill(cost[i], Float.MAX_VALUE);
            
            CornerRegion c3C1 = c1.get(i);
            if (c3C1.equals(cr1C1) || c3C1.equals(cr2C1)) {
                continue;
            }
            
            double[] xy3C2 = applyTransformation(c3C1, params);
            CornersAndFeatureStat cfs = findBestMatchWithinTolerance(c3C1, xy3C2,
                Math.round(params.getRotationInDegrees()), rotationTolerance);
   //TODO: consider discarding already chosen in findBestMatchWithinTolerance      
            if ((cfs == null) || cfs.cr2.equals(cr1C2) || cfs.cr2.equals(cr2C2)) {
                continue;
            }
            cfs.idx1 = i;
            
            assert(cfs.idx2 != -1);
            
            double dist = distance(cfs.stat.getImg2Point().getX(), 
                cfs.stat.getImg2Point().getY(), xy3C2[0], xy3C2[1]);
            
            cost[i][cfs.idx2] = (float)dist;
            
            PairInt pI = new PairInt(i, cfs.idx2);
            indexesCFSMap.put(pI, cfs);
            indexesDistMap.put(pI, Double.valueOf(dist));            
        }
        
        if (indexesCFSMap.isEmpty()) {
            return;
        }
        
        if (indexesCFSMap.size() == 1) {
            // TODO: handle without Hungarian method
        }
        
        boolean transposed = false;
        if (c1.size() > c2.size()) {
            cost = MatrixUtil.transpose(cost);
            transposed = true;
        }

        HungarianAlgorithm b = new HungarianAlgorithm();
        int[][] match = b.computeAssignments(cost);
        
        int nC = 0;
        for (int i = 0; i < match.length; i++) {
            int idx1 = match[i][0];
            int idx2 = match[i][1];
            if (idx1 == -1 || idx2 == -1) {
                continue;
            }
            // points not in map were matched with another max cost, so not matched           
            if (indexesCFSMap.containsKey(new PairInt(idx1, idx2))) {
                 nC++;
            }
        }
        if (nC == 0) {
            return;
        }
        
        for (int i = 0; i < match.length; i++) {

            int idx1 = match[i][0];
            int idx2 = match[i][1];
            if (idx1 == -1 || idx2 == -1) {
                continue;
            }

            if (transposed) {
                int swap = idx1;
                idx1 = idx2;
                idx2 = swap;
            }
            
            PairInt pI = new PairInt(idx1, idx2);
            
            CornersAndFeatureStat cfs = indexesCFSMap.get(pI);
            if (cfs == null) {
                continue;
            }
            nc.addMatchedCorners(cfs.cr1, cfs.cr2, cfs.stat);
            
            double costAsDist = indexesDistMap.get(pI).doubleValue();
            
            transformationPair.addToCostAsDistance(costAsDist);
            
            transformationPair.addToCostAsSSD(cfs.stat.getSumIntensitySqDiff());
        }      
             
        double costSSD = transformationPair.getCostAsSSD();
        
        // put node in heap
        long costL = (long)(costSSD * heapKeyFactor);

        HeapNode node = new HeapNode(costL);

        node.setData(transformationPair);

        heap.insert(node);
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
     * return best SSD match to point within tolerance of predictedXYCurve2
     * as null if no match, else an instance of IndexAndFeatureStat
     * @param cornerCurve1
     * @param predictedXYCurve2
     * @return 
     */
    private CornersAndFeatureStat findBestMatchWithinTolerance(
        CornerRegion cornerCurve1, double[] predictedXYCurve2) {
        
        FeatureMatcher featureMatcher = new FeatureMatcher();

        FeatureComparisonStat bestStat = null;

        int bestIdx2 = -1;

        Set<Integer> c2Indexes = np.findNeighborIndexes(
            (int) Math.round(predictedXYCurve2[0]),
            (int) Math.round(predictedXYCurve2[1]), tolerance);

        for (Integer index : c2Indexes) {

            int idx2 = index.intValue();

            CornerRegion corner2 = c2.get(idx2);

            FeatureComparisonStat compStat = null;

            try {

                compStat = featureMatcher.ditherAndRotateForBestLocation(
                    features1, features2, cornerCurve1, corner2, dither);

            } catch (CornerRegion.CornerRegionDegneracyException ex) {

                log.severe(ex.getMessage());
            }

            if (compStat == null) {
                continue;
            }

            if (bestStat == null) {
                if (compStat.getSumIntensitySqDiff() < compStat.getImg2PointIntensityErr()) {
                    bestStat = compStat;
                    bestIdx2 = idx2;
                }
                continue;
            }

            if (compStat.getSumIntensitySqDiff() < compStat.getImg2PointIntensityErr()) {
                if (compStat.getSumIntensitySqDiff() < bestStat.getSumIntensitySqDiff()) {
                    bestStat = compStat;
                    bestIdx2 = idx2;
                }
            }           
        }
        
        if (bestIdx2 > -1) {
            return new CornersAndFeatureStat(cornerCurve1, c2.get(bestIdx2), 
                bestStat);
        } else {
            return null;
        } 
    }

    private HeapNode solve() {
        
        HeapNode u = heap.extractMin();

        while (u != null) {

            TransformationPair2 transformationPair = (TransformationPair2)u.getData();

            NextCorner nc = transformationPair.getNextCorner();

            CornerRegion corner1s = nc.findStrongestRemainingCorner();

            if (corner1s == null) {
                return u;
            }

            TransformationParameters params 
                = transformationPair.getTransformationParameters();

            // apply the transformation to it
            double[] xyC2 = applyTransformation(corner1s, params);

            // find the best matching SSD within tolerance of predicted
            CornersAndFeatureStat cfs = findBestMatchWithinTolerance(corner1s, 
                xyC2, Math.round(params.getRotationInDegrees()), 
                rotationTolerance);
                        
            double cost2;
            double cost2SSD;
            
            if (cfs != null) {
                                                                
                nc.addMatchedCorners(corner1s, cfs.cr2, cfs.stat);
                
                cost2 = calculateCost(cfs.stat.getImg2Point(), xyC2[0], xyC2[1]);
                
                cost2SSD = cfs.stat.getSumIntensitySqDiff();
                
                // TODO: revise SSD as cost... may need to be added in quadrature
                
                transformationPair.addToCostAsDistance(cost2);
                
                transformationPair.addToCostAsSSD(cost2SSD);
                
            } else {
                
                cost2 = calculateUnMatchedCost();
                
                cost2SSD = 10000;
            }            
                        
            u.setData(transformationPair);

            u.setKey(u.getKey() + (long)(cost2SSD * heapKeyFactor));

            heap.insert(u);

            u = heap.extractMin();
        }

        return u;
    }
    
}
