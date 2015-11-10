package algorithms.imageProcessing;

import algorithms.compGeometry.NearestPoints;
import algorithms.imageProcessing.util.AngleUtil;
import algorithms.imageProcessing.util.MatrixUtil;
import algorithms.util.PairInt;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.logging.Logger;
import thirdparty.HungarianAlgorithm;

/**
 *
 *
 * @author nichole
 */
public class ClosedCurveCornerMatcher {
    
    /**
     * the costs calculated here are small fractions, so they need to be
     * multiplied by a large constant for use with the Fibonacci heap
     * which uses type long for its key (key is where cost is stored).
     * using 1E12 here
     */
    protected final static long heapKeyFactor = 1000000000000l;

    protected TransformationPair2 solution = null;
    
    protected final int dither = 4;//3;

    protected final int rotationTolerance = 20;

    private int tolerance = 2;//4;

    private double solutionCost = Double.MAX_VALUE;

    private final Logger log = Logger.getLogger(this.getClass().getName());

    private boolean solverHasFinished = false;

    private boolean solutionHasSomeScalesSmallerThanOne = false;

    private Float fixedRotationInDegrees = null;
    
    public ClosedCurveCornerMatcher() {
    }
    
    public void setFixedRotation(float rD) {
        fixedRotationInDegrees = Float.valueOf(rD);
    }

    private void resetDefaults() {
        solverHasFinished = false;
        solutionHasSomeScalesSmallerThanOne = false;
        solutionCost = Double.MAX_VALUE;
        solution = null;
    }
    
    /**
     * 
     * @param <T> T extends CornerRegion
     * @param features1
     * @param features2
     * @param corners1
     * @param corners2
     * @param cornersAreAlreadySorted
     * @param img1 image from which to extract descriptors for features1
     * @param img2 image from which to extract descriptors for features2
     * @return 
     */
    public <T extends CornerRegion> boolean matchCorners(
        final IntensityFeatures features1, final IntensityFeatures features2, 
        final List<T> corners1,final List<T> corners2, 
        boolean cornersAreAlreadySorted, 
        GreyscaleImage img1, GreyscaleImage img2) {

        if (solverHasFinished) {
            resetDefaults();
        }
        
        List<CornerRegion> c1 = new ArrayList<CornerRegion>(corners1.size());
        List<CornerRegion> c2 = new ArrayList<CornerRegion>(corners2.size());
        c1.addAll(corners1);
        c2.addAll(corners2);
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
        
        NearestPoints np = new NearestPoints(xC2, yC2);

        Heap heap = new Heap();
        
        initializeHeapNodes(heap, c1, c2, np, features1, features2, img1, img2);

        if (heap.getNumberOfNodes() == 0) {
            solverHasFinished = true;
            return false;
        }
        
        //evaluating entire solutions for top 10 below instead of just best in heap
        HeapNode minCost = heap.extractMin();

        if (minCost == null) {
            solverHasFinished = true;
            return false;
        }

        // ========= begin combining results ==========
        List<TransformationParameters> topParams
            = new ArrayList<TransformationParameters>();

        List<HeapNode> top = new ArrayList<HeapNode>();
        top.add(minCost);
        topParams.add(((TransformationPair2)minCost.getData())
            .getTransformationParameters());
        HeapNode node = heap.extractMin();
        while ((node != null) && (top.size() < 11)) {
            TransformationPair2 tPair = (TransformationPair2)node.getData();
            int n = tPair.getNextCorner().getMatchedCornerIndexes1().size();
            if (n > 2) {
                top.add(node);
                topParams.add(tPair.getTransformationParameters());
            }
            node = heap.extractMin();
        }

        // combine similar results
        // for example, in one test dataset can see the correct solution is
        // at top[1] and top[2]
        MatchedPointsTransformationCalculator tc = new MatchedPointsTransformationCalculator();
        Map<Integer, Set<Integer>> similarParams = tc.findSimilarByScaleAndRotation(
            topParams);
        
        // aggregate the points and calculate the total SSD/nPoints
        List<Set<PairInt>> topMatchedCornerIndexes = new ArrayList<Set<PairInt>>();
        List<Set<Integer>> topMatchedIndexes = new ArrayList<Set<Integer>>();
        List<Double> topMatchedSSD = new ArrayList<Double>();
        List<Map<PairInt, FeatureComparisonStat>> topMatchedList 
            = new ArrayList<Map<PairInt, FeatureComparisonStat>>();
        double bestSSD = Double.MAX_VALUE;
        int bestSSDIdx = -1;
        for (Entry<Integer, Set<Integer>> entry : similarParams.entrySet()) {
            Set<Integer> index1Added = new HashSet<Integer>();
            Set<Integer> index2Added = new HashSet<Integer>();
            Set<PairInt> combined = new HashSet<PairInt>();
            Map<PairInt, FeatureComparisonStat> fcsMap = new HashMap<PairInt, FeatureComparisonStat>();
            double ssdSum = 0;
            for (Integer index : entry.getValue()) {
                TransformationPair2 tPair = (TransformationPair2)top.get(index.intValue()).getData();
                List<Integer> indexes1 = tPair.getNextCorner().getMatchedCornerIndexes1();
                List<Integer> indexes2 = tPair.getNextCorner().getMatchedCornerIndexes2();
                List<FeatureComparisonStat> stats = tPair.getNextCorner().getMatchedFeatureComparisonStats();
                for (int i = 0; i < indexes1.size(); ++i) {
                    Integer index1 = indexes1.get(i);
                    Integer index2 = indexes2.get(i);
                    PairInt p = new PairInt(index1.intValue(), index2.intValue());
                    if (!combined.contains(p)) {
                        if (!index1Added.contains(index1) && !index2Added.contains(index2)) {
                            combined.add(p);
                            index1Added.add(index1);
                            index2Added.add(index2);
                            fcsMap.put(p, stats.get(i));
                            ssdSum += stats.get(i).getSumIntensitySqDiff();
                        }
                    }
                }
            }
            topMatchedCornerIndexes.add(combined);
            topMatchedIndexes.add(entry.getValue());
            topMatchedList.add(fcsMap);
            ssdSum /= (double)combined.size();
            topMatchedSSD.add(Double.valueOf(ssdSum));
            if (ssdSum < bestSSD) {
                bestSSD = ssdSum;
                bestSSDIdx = topMatchedCornerIndexes.size() - 1;
            }
        }
        
        if (!similarParams.isEmpty() && (topMatchedCornerIndexes.get(bestSSDIdx).size() > 2)) {
        
            // TODO: consider re-calculating the TransformationParameter instead of
            // averaging knowing that the origins (x,y) are the same
            Set<Integer> minCostIndexes = topMatchedIndexes.get(bestSSDIdx);
            int[] rotationsInDegrees = new int[minCostIndexes.size()];
            int count = 0;
            for (Integer topIndex : minCostIndexes) {
                TransformationPair2 tPair = (TransformationPair2)top.get(topIndex.intValue()).getData();
                TransformationParameters params2 = tPair.getTransformationParameters();
                rotationsInDegrees[count] = Math.round(params2.getRotationInDegrees());
                count++;
            }
            float avgRotDegrees = AngleUtil.calculateAverageWithQuadrantCorrections(
                rotationsInDegrees, rotationsInDegrees.length - 1);

            float tXSum = 0;
            float tYSum = 0;
            float scaleSum = 0;
            for (Integer topIndex : minCostIndexes) {
                TransformationPair2 tPair = (TransformationPair2)top.get(topIndex.intValue()).getData();
                TransformationParameters params2 = tPair.getTransformationParameters();
                tXSum += params2.getTranslationX();
                tYSum += params2.getTranslationY();
                scaleSum += params2.getScale();
            }
            tXSum /= (float)minCostIndexes.size();
            tYSum /= (float)minCostIndexes.size();
            scaleSum /= (float)minCostIndexes.size();
            TransformationParameters params = new TransformationParameters();
            params.setOriginX(0);
            params.setOriginY(0);
            params.setRotationInDegrees(avgRotDegrees);
            params.setScale(scaleSum);
            params.setTranslationX(tXSum);
            params.setTranslationY(tYSum);

            List<Integer> combinedIndexes1 = new ArrayList<Integer>();
            List<Integer> combinedIndexes2 = new ArrayList<Integer>();
            List<FeatureComparisonStat> fcsList = new ArrayList<FeatureComparisonStat>();
            for (PairInt p : topMatchedCornerIndexes.get(bestSSDIdx)) {
                combinedIndexes1.add(Integer.valueOf(p.getX()));
                combinedIndexes2.add(Integer.valueOf(p.getY()));
                fcsList.add(topMatchedList.get(bestSSDIdx).get(p));
            }

            TransformationPair2 transformationPair = new TransformationPair2(
                c1.get(combinedIndexes1.get(0).intValue()),
                c2.get(combinedIndexes2.get(0).intValue()),
                c1.get(combinedIndexes1.get(1).intValue()),
                c2.get(combinedIndexes2.get(1).intValue()));

            transformationPair.setTransformationParameters(params);

            NextCorner nc = new NextCorner(c1);
            //for (int i = 2; i < combinedIndexes1.size(); ++i) {
            for (int i = 0; i < combinedIndexes1.size(); ++i) {
                int idx1 = combinedIndexes1.get(i).intValue();
                int idx2 = combinedIndexes2.get(i).intValue();
                FeatureComparisonStat stat = fcsList.get(i);
                nc.addMatchedCorners(c1.get(idx1), c2.get(idx2), idx1, idx2, stat);
            }
            transformationPair.setNextCorner(nc);
            
            long costL = (long)(bestSSD * heapKeyFactor);
            minCost = new HeapNode(costL);
            minCost.setData(transformationPair);
            
            // ========= end combining results ==========
            
        } else {
            // keep previous minCost node
        }
        
        TransformationPair2 transformationPair = (TransformationPair2)minCost.getData();
        
        solution = transformationPair;
        
        solutionCost = ((double)minCost.getKey()/(double)heapKeyFactor);
        
        solverHasFinished = true;

        solutionHasSomeScalesSmallerThanOne = transformationPair.scaleIsPossiblyAmbiguous();

        return true;
    }

    private CornersAndFeatureStat findBestMatchWithinTolerance(
        List<CornerRegion> c1, List<CornerRegion> c2, NearestPoints np,
        IntensityFeatures features1, IntensityFeatures features2,
        int corner1Idx, CornerRegion cornerCurve1, double[] predictedXYCurve2,
        int rotationInDegrees, int rotationTolerance,
        GreyscaleImage img1, GreyscaleImage img2) {

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
                    rotationInDegrees, rotationTolerance, img1, img2);

            } catch (CornerRegion.CornerRegionDegneracyException ex) {
                log.fine(ex.getMessage());
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
            cfs.idx1 = corner1Idx;
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
     * @param c1
     * @param c2
     * @param features1
     * @param features2
     * @param img1 image from which to extract descriptors for features1
     * @param img2 image from which to extract descriptors for features2
     * @return
     */
    protected CornersAndFeatureStat[] getBestSSDC1ToC2(
        List<CornerRegion> c1, List<CornerRegion> c2,
        IntensityFeatures features1, IntensityFeatures features2, 
        GreyscaleImage img1, GreyscaleImage img2) {

        // if no match, contains a null
        CornersAndFeatureStat[] indexes2 = new CornersAndFeatureStat[c1.size()];

        // store by index2 number so a later index1 with a worse match to index2
        // will not be assigned index2.  (note, not using bipartite matching
        // because the method would increase from approx O(N^2) to approx O(N^3)
        Map<Integer, CornersAndFeatureStat> index2Map = new HashMap<Integer,
            CornersAndFeatureStat>();
        
        FeatureMatcher featureMatcher = new FeatureMatcher();

        for (int i = 0; i < c1.size(); ++i) {

            CornerRegion region1 = c1.get(i);

            FeatureComparisonStat best = null;

            int bestIdx2 = -1;

            for (int j = 0; j < c2.size(); ++j) {

                CornerRegion region2 = c2.get(j);

                FeatureComparisonStat compStat = null;

                try {

                    if (fixedRotationInDegrees == null) {
                        compStat = featureMatcher.ditherAndRotateForBestLocation(
                            features1, features2, region1, region2, dither, 
                            img1, img2);
                    } else {
                        compStat = featureMatcher.ditherAndRotateForBestLocation(
                            features1, features2, region1, region2, dither, 
                            fixedRotationInDegrees.intValue(), 20,
                            img1, img2);
                    }

                } catch (CornerRegion.CornerRegionDegneracyException ex) {
                    log.fine("**CONSIDER using more points in corner region");
                }

                if (compStat == null) {
                    continue;
                }

                if (compStat.getSumIntensitySqDiff() >= compStat.getImg2PointIntensityErr()) {
                    continue;
                }
                
                CornersAndFeatureStat existing = index2Map.get(Integer.valueOf(j));

                if ((existing != null) && (existing.stat.getSumIntensitySqDiff() <
                    compStat.getSumIntensitySqDiff())) {
                    continue;
                }
                
                if (best == null) {
                    best = compStat;
                    bestIdx2 = j;
                } else {
                    if (compStat.getSumIntensitySqDiff() < best.getSumIntensitySqDiff()) {
                        best = compStat;
                        bestIdx2 = j;
                    }
                }
            }

            if (bestIdx2 > -1) {
                indexes2[i] = new CornersAndFeatureStat(c1.get(i),
                    c2.get(bestIdx2), best);
                indexes2[i].idx1 = i;
                indexes2[i].idx2 = bestIdx2;
                index2Map.put(Integer.valueOf(bestIdx2), indexes2[i]);
            }
        }

        return indexes2;
    }

    /**
     * initialize the starter solution nodes for the heap
     */
    private void initializeHeapNodes(Heap heap, List<CornerRegion> c1,
        List<CornerRegion> c2, NearestPoints np, IntensityFeatures features1,
        IntensityFeatures features2, GreyscaleImage img1, GreyscaleImage img2) {

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

        initHeapNodes1(heap, c1, c2, np, features1, features2, img1, img2);
    }

    /**
     * create starter solution nodes for the heap using a pattern of combinations
     * with runtime O(n*(n-1)/2).
     *
     */
    private void initHeapNodes1(Heap heap, List<CornerRegion> c1,
        List<CornerRegion> c2, NearestPoints np, IntensityFeatures features1,
        IntensityFeatures features2, GreyscaleImage img1, GreyscaleImage img2) {

        // if no match, contains a null
        CornersAndFeatureStat[] indexes2 = getBestSSDC1ToC2(c1, c2,
            features1, features2, img1, img2);

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

        for (int i1_1 = 0; i1_1 < indexes2.length; ++i1_1) {

            CornerRegion cr1C1 = c1.get(i1_1);

            if (indexes2[i1_1] == null) {
                continue;
            }

            CornerRegion cr1C2 = indexes2[i1_1].cr2;

            for (int i1_2 = (i1_1 + 1); i1_2 < c1.size(); ++i1_2) {

                CornerRegion cr2C1 = c1.get(i1_2);

                if (indexes2[i1_2] == null) {
                    continue;
                }

                CornerRegion cr2C2 = indexes2[i1_2].cr2;

                if (cr2C2.equals(cr1C2)) {
                    continue;
                }

                // temporarily, evaluating all corners for each starter solution:
                insertNodeTMP(heap, c1, c2, np, features1, features2, 
                    tc, indexes2[i1_1], cr1C1, cr1C2, indexes2[i1_2], cr2C1, cr2C2,
                    img1, img2);
               
            }
        }
    }

    private double distance(double x1, double y1, double x2, double y2) {

        double diffX = x1 - x2;
        double diffY = y1 - y2;

        double dist = Math.sqrt(diffX*diffX + diffY*diffY);

        return dist;
    }

    public TransformationPair2 getSolution() {
        return solution;
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

    private void insertNodeTMP(Heap heap, 
        List<CornerRegion> c1, List<CornerRegion> c2, NearestPoints np,
        IntensityFeatures features1, IntensityFeatures features2,
        MatchedPointsTransformationCalculator tc,
        CornersAndFeatureStat stat1, CornerRegion cr1C1, CornerRegion cr1C2,
        CornersAndFeatureStat stat2, CornerRegion cr2C1, CornerRegion cr2C2,
        GreyscaleImage img1, GreyscaleImage img2) {

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
        
        if (params == null) {
            return;
        }

        TransformationPair2 transformationPair =
            new TransformationPair2(cr1C1, cr1C2, cr2C1, cr2C2);

        transformationPair.setTransformationParameters(params);

        NextCorner nc = new NextCorner(c1);
        nc.addMatchedCorners(cr1C1, cr1C2, stat1.idx1, stat1.idx2, stat1.stat);
        nc.addMatchedCorners(cr2C1, cr2C2, stat2.idx1, stat2.idx2, stat2.stat);

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
            CornersAndFeatureStat cfs = findBestMatchWithinTolerance(c1, c2, np,
                features1, features2,
                i, c3C1, xy3C2, Math.round(params.getRotationInDegrees()), 
                rotationTolerance, img1, img2);
            
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
            nc.addMatchedCorners(cfs.cr1, cfs.cr2, Integer.valueOf(cfs.idx1),
                Integer.valueOf(cfs.idx2), cfs.stat);

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
}
