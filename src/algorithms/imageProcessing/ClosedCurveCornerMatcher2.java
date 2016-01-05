package algorithms.imageProcessing;

import algorithms.compGeometry.NearestPoints;
import algorithms.imageProcessing.util.MiscStats;
import algorithms.util.PairInt;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;

/**
 *
 *
 * @author nichole
 */
public class ClosedCurveCornerMatcher2<T extends CornerRegion> {

    // a minimum number to use when calculating the sums of SSDs for
    // normalized intensities
    private float lowerLimitSSD = 200.f;
    
    /**
     * the costs calculated here are small fractions, so they need to be
     * multiplied by a large constant for use with the Fibonacci heap
     * which uses type long for its key (key is where cost is stored).
     * using 1E12 here
    */
    protected final static long heapKeyFactor = 1000000000000l;

    protected TransformationParameters solution = null;

    private int nEval = -1;

    protected final int rotationTolerance = 20;

    private static int tolerance = 3;//4;

    private static double maxDistance = Math.sqrt(2) * tolerance;

    private final Logger log = Logger.getLogger(this.getClass().getName());

    private double solutionCost;

    private List<FeatureComparisonStat> solutionStats = null;;

    private enum State {
        INITIALIZED, FAILED, SOLVED
    }

    private State state = null;
    
    private final int dither;

    public ClosedCurveCornerMatcher2(int dither) {
        this.dither = dither;
    }

    private void resetDefaults() {
        state = null;
        solution = null;
        nEval = -1;
    }

    /**
     *
     * @param features1
     * @param features2
     * @param corners1
     * @param corners2
     * @param img1 image from which to extract descriptors for features1
     * @param img2 image from which to extract descriptors for features2
     * @return
     */
    @SuppressWarnings({"unchecked"})
    public boolean matchCorners(
        final IntensityFeatures features1, final IntensityFeatures features2,
        final List<T> corners1,final List<T> corners2,
        GreyscaleImage img1, GreyscaleImage img2, int binFactor1,
        int binFactor2) {

        if (state != null) {
            resetDefaults();
        }

        List<T> c1 = new ArrayList<T>(corners1.size());
        List<T> c2 = new ArrayList<T>(corners2.size());
        c1.addAll(corners1);
        c2.addAll(corners2);

        int[] xC2 = new int[c2.size()];
        int[] yC2 = new int[c2.size()];
        for (int i = 0; i < c2.size(); ++i) {
            T cr = c2.get(i);
            xC2[i] = cr.getX()[cr.getKMaxIdx()];
            yC2[i] = cr.getY()[cr.getKMaxIdx()];
        }

        state = State.INITIALIZED;

        NearestPoints np = new NearestPoints(xC2, yC2);

        findBestParameters(c1, c2, np, features1, features2, img1, img2,
            binFactor1, binFactor2);

        return state.equals(State.SOLVED);

    }

    private static class CornersAndFeatureStat <T extends CornerRegion> {
        private final T cr1;
        private final T cr2;
        private int idx1 = -1;
        private int idx2 = -1;
        private final FeatureComparisonStat stat;
        private final double dist;
        private final double normalizedCost;
        public CornersAndFeatureStat(T cornerRegion1,
            T cornerRegion2, FeatureComparisonStat compStat, double distance) {
            cr1 = cornerRegion1;
            cr2 = cornerRegion2;
            stat = compStat;
            dist = distance;
            normalizedCost = stat.getSumIntensitySqDiff()
                * (dist/maxDistance);
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
    @SuppressWarnings({"unchecked", "rawtypes"})
    protected CornersAndFeatureStat<T>[] getBestSSDC1ToC2(
        List<T> c1, List<T> c2,
        IntensityFeatures features1, IntensityFeatures features2,
        GreyscaleImage img1, GreyscaleImage img2) {

        // if no match, contains a null
        CornersAndFeatureStat<T>[] indexes2 = new CornersAndFeatureStat[c1.size()];

        // store by index2 number so a later index1 with a worse match to index2
        // will not be assigned index2.  (note, not using bipartite matching
        // because the method would increase from approx O(N^2) to approx O(N^3)
        Map<Integer, CornersAndFeatureStat<T>> index2Map = new HashMap<Integer,
            CornersAndFeatureStat<T>>();
/*
double[][] xy1 = new double[c1.size()][2];
for (int i = 0; i < c1.size(); ++i) {
CornerRegion cr = c1.get(i);
xy1[i] = new double[]{cr.getX()[cr.getKMaxIdx()], cr.getY()[cr.getKMaxIdx()]};
}
double[][] xy2 = new double[c2.size()][2];
for (int i = 0; i < c2.size(); ++i) {
CornerRegion cr = c2.get(i);
xy2[i] = new double[]{cr.getX()[cr.getKMaxIdx()], cr.getY()[cr.getKMaxIdx()]};
}
*/
        FeatureMatcher featureMatcher = new FeatureMatcher();

        for (int i = 0; i < c1.size(); ++i) {

            T region1 = c1.get(i);

            double bestCost = Double.MAX_VALUE;
            FeatureComparisonStat best = null;
            int bestIdx2 = -1;

            for (int j = 0; j < c2.size(); ++j) {

                T region2 = c2.get(j);

                FeatureComparisonStat compStat = 
                    featureMatcher.ditherAndRotateForBestLocation2(
                    features1, features2, region1, region2, dither,
                    img1, img2);
                
                if ((compStat == null) ||
                    (compStat.getSumIntensitySqDiff() > compStat.getImg2PointIntensityErr())
//NOTE: temporary change for lower limit filter to remove featureless patches
 || (compStat.getImg2PointIntensityErr() < lowerLimitSSD)
                    ) {
                    continue;
                }

                if ((best == null) || (compStat.getSumIntensitySqDiff() < bestCost)) {
                    best = compStat;
                    bestIdx2 = j;
                    bestCost = compStat.getSumIntensitySqDiff();
                }
            }

            if (bestIdx2 > -1) {
                indexes2[i] = new CornersAndFeatureStat<T>(region1, 
                    c2.get(bestIdx2), best, 0);
                indexes2[i].idx1 = i;
                indexes2[i].idx2 = bestIdx2;
                index2Map.put(Integer.valueOf(bestIdx2), indexes2[i]);
            }
        }

        return indexes2;
    }
    

    /**
     * return all combinations of c1 and c2 corners.  note this method should
     * only be used when sizes of c1 and c2 are very small, it produces
     * n1 * n2 combinations.
     * @param c1
     * @param c2
     * @param features1
     * @param features2
     * @param img1 image from which to extract descriptors for features1
     * @param img2 image from which to extract descriptors for features2
     * @return
     */
    @SuppressWarnings({"unchecked", "rawtypes"})
    protected CornersAndFeatureStat<T>[] getAllSSDC1ToC2(
        List<T> c1, List<T> c2,
        IntensityFeatures features1, IntensityFeatures features2,
        GreyscaleImage img1, GreyscaleImage img2) {

        int nComb = c1.size() * c2.size();

        // if no match, contains a null
        CornersAndFeatureStat<T>[] indexes2 = new CornersAndFeatureStat[nComb];

        FeatureMatcher featureMatcher = new FeatureMatcher();

        int count = 0;

        for (int i = 0; i < c1.size(); ++i) {

            T region1 = c1.get(i);

            for (int j = 0; j < c2.size(); ++j) {

                T region2 = c2.get(j);

                FeatureComparisonStat compStat =
                    featureMatcher.ditherAndRotateForBestLocation2(
                    features1, features2, region1, region2, dither,
                    img1, img2);

                if ((compStat == null) ||
                    (compStat.getSumIntensitySqDiff() > compStat.getImg2PointIntensityErr())
                    ) {
                    continue;
                }

                indexes2[count] = new CornersAndFeatureStat<T>(region1, region2,
                    compStat, 0);
                indexes2[count].idx1 = i;
                indexes2[count].idx2 = j;
                count++;
            }
        }

        return indexes2;
    }

    /**
     * find best transformation parameters using a pattern of combinations
     * with runtime O(n*(n-1)/2).
     * TODO: revise the runtime complexity comment for new changes
     */
    @SuppressWarnings({"unchecked", "rawtypes"})
    private void findBestParameters(List<T> c1,
        List<T> c2, NearestPoints np, IntensityFeatures features1,
        IntensityFeatures features2, GreyscaleImage img1, GreyscaleImage img2,
        int binFactor1, int binFactor2) {

        /*
        for each point in curve c1, the best matching point by features in c2 
        is given.  if no match, the entry contains a null.
        */
        CornersAndFeatureStat<T>[] indexes2 = null;
        
        boolean keepAllCombinations = true;
        
        if (keepAllCombinations) {
            indexes2 = getAllSSDC1ToC2(c1, c2, features1, features2, img1, img2);
        } else {
            if (c1.size() < 5 && c2.size() < 5) {
                indexes2 = getAllSSDC1ToC2(c1, c2, features1, features2, img1, img2);
            } else {
                indexes2 = getBestSSDC1ToC2(c1, c2, features1, features2, img1, img2);
            }
        }

        int nIndexes2 = 0;
        for (CornersAndFeatureStat<T> cfs : indexes2) {
            if (cfs != null) {
                nIndexes2++;
            }
        }
        
        int nTop;
        if (keepAllCombinations) {
            nTop = nIndexes2;
        } else {
            if (c1.size() < 5 && c2.size() < 5) {
                nTop = nIndexes2;
            } else {
                nTop = (nIndexes2 < 11) ? nIndexes2 : 10;
            }
        }
        
        if (nTop != nIndexes2) {
            Heap heap = new Heap();
            for (int i = 0; i < indexes2.length; ++i) {
                CornersAndFeatureStat<T> cfs = indexes2[i];
                if (cfs == null) {
                    continue;
                }
                long costL = (long)(cfs.stat.getSumIntensitySqDiff() * heapKeyFactor);
                HeapNode node = new HeapNode(costL);
                node.setData(cfs);
                heap.insert(node);
            }
            indexes2 = new CornersAndFeatureStat[nTop];
            int count = 0;
            while ((heap.getNumberOfNodes() > 0) && (count < nTop)) {
                indexes2[count] = (CornersAndFeatureStat<T>)heap.extractMin().getData();
                count++;
            }
        }

        MatchedPointsTransformationCalculator
            tc = new MatchedPointsTransformationCalculator();

        // for the true curve to curve match, there should be at least 2
        //     correct matches in indexes2
        // so, look at the parameters derived by pairing combinations of 2
        //  points in indexes2, and consolidate the answers.

        // this only needs to be a list, but keeping a map for debugging:
        Map<PairInt, TransformationParameters> paramsMap = new
            HashMap<PairInt, TransformationParameters>();

        for (int ipt1 = 0; ipt1 < indexes2.length; ++ipt1) {
            CornersAndFeatureStat<T> cfs1 = indexes2[ipt1];
            if (cfs1 == null) {
                continue;
            }
            PairInt c1Pt1 = cfs1.stat.getImg1Point();
            PairInt c2Pt1 = cfs1.stat.getImg2Point();
            for (int ipt2 = (ipt1 + 1); ipt2 < indexes2.length; ++ipt2) {
                CornersAndFeatureStat<T> cfs2 = indexes2[ipt2];
                if (cfs2 == null) {
                    continue;
                }
                if (cfs1.cr2.equals(cfs2.cr2)) {
                    continue;
                }

                PairInt c1Pt2 = cfs2.stat.getImg1Point();
                PairInt c2Pt2 = cfs2.stat.getImg2Point();

                TransformationParameters params = tc.calulateEuclidean(
                    c1Pt1.getX(), c1Pt1.getY(),
                    c1Pt2.getX(), c1Pt2.getY(),
                    c2Pt1.getX(), c2Pt1.getY(),
                    c2Pt2.getX(), c2Pt2.getY(), 0, 0);
                
                if (params == null) {
                    //only happens when both curve1 points are the same or
                    // both curve2 points are the same.
                    // this should have been caught above.
                    continue;
                }
                
                params.setNumberOfPointsUsed(2);

                paramsMap.put(new PairInt(ipt1, ipt2), params);
            }
        }

        if (paramsMap.isEmpty()) {
            state = State.FAILED;
            return;
        }

        // requiring that a params0 has at least one similar compare
        //     before adding to combinedParams in order to reduce the number
        //     of evaluations.
        List<TransformationParameters> combinedParams = null;

        if (keepAllCombinations) {
            combinedParams = MiscStats.filterToSimilarParamSets2(paramsMap,
                binFactor1, binFactor2);
        } else {
            if (c1.size() < 5 && c2.size() < 5) {
                combinedParams = MiscStats.filterToSimilarParamSets2(paramsMap,
                    binFactor1, binFactor2);
            } else {
                combinedParams = MiscStats.filterToSimilarParamSets(paramsMap,
                    binFactor1, binFactor2);
            }
        }
        
        // --- evaluate transformations on all corners, starting w/ location ----

        Transformer transformer = new Transformer();
        
        List<TransformationParameters> combinedParams2 = 
            new ArrayList<TransformationParameters>();
        
        int topK;
        
        if (keepAllCombinations) {
            combinedParams2.addAll(combinedParams);
            topK = combinedParams2.size();
        } else {
            Heap orderedParams = new Heap();
            for (TransformationParameters params : combinedParams) {
                double sumDist = 0;
                int nEval2 = 0;
                int tolTransXY = tolerance;
                if (params.getScale() < 1) {
                    tolTransXY = Math.round(tolerance * params.getScale());
                }
                if (params.getStandardDeviations() != null) {
                    tolTransXY = (int)Math.ceil(Math.max(
                        Math.abs(params.getStandardDeviations()[2]),
                        Math.abs(params.getStandardDeviations()[3])
                        ));
                }
                if (tolTransXY == 0) {
                    tolTransXY = 1;
                }

                for (int ipt1 = 0; ipt1 < c1.size(); ++ipt1) {
                    T pt1 = c1.get(ipt1);

                    double[] xyTr = transformer.applyTransformation(params,
                        pt1.getX()[pt1.getKMaxIdx()], pt1.getY()[pt1.getKMaxIdx()]);

                    Set<PairInt> candidates = np.findNeighbors(
                        (int) Math.round(xyTr[0]), (int) Math.round(xyTr[1]),
                        tolTransXY);

                    if (candidates != null && candidates.size() > 0) {
                        for (PairInt p : candidates) {
                            double diffX = xyTr[0] - p.getX();
                            double diffY = xyTr[1] - p.getY();

                            double dist = Math.sqrt(diffX*diffX + diffY*diffY);

                            nEval2++;
                            sumDist += dist;
                        }
                    }
                }
                if (nEval2 == 0) {
                    continue;
                }

                // distance needs to be adjusted by scale, else the cost prefers
                // small scale solutions
                sumDist /= params.getScale();
                sumDist /= (double)nEval2;
                // store in heap.  use nEval and dist as cost
                float cost1 = 1.f/(float)nEval2;
                float cost2 = (float)(sumDist + 1)/(float)tolTransXY;
                float normalizedCost = cost1 * cost2;
    
                long costL = (long)(normalizedCost * heapKeyFactor);
                HeapNode node = new HeapNode(costL);
                node.setData(params);
                orderedParams.insert(node);
            }
            if (c1.size() < 5 && c2.size() < 5) {
                topK = (int) orderedParams.getNumberOfNodes();
            } else {
                topK = c2.size() / 2;
            }
            int count = 0;
            while ((orderedParams.getNumberOfNodes() > 0) && (count < topK)) {
                TransformationParameters params = (TransformationParameters)
                    orderedParams.extractMin().getData();
                combinedParams2.add(params);
                count++;
            }
        }

        // --- evaluate transformations on all corners, use features ----
        //     transform c1 and count matches to c2
        List<FeatureComparisonStat> bestStats = null;
        TransformationParameters bestParams = null;
        int maxNEval = Integer.MIN_VALUE;
        double minDist = Double.MAX_VALUE;
        double minSSD = Double.MAX_VALUE;
        double minCost = Double.MAX_VALUE;

        int nMaxMatchable = Math.min(c1.size(), c2.size());
        
        FeatureMatcher featureMatcher = new FeatureMatcher();
        
//log.info("EVAL by dist and SSD:\n");

        for (TransformationParameters params : combinedParams2) {

            int tolXY = tolerance;
            if (params.getScale() < 1) {
                tolXY = Math.round(tolerance * params.getScale());
            }
            if (params.getStandardDeviations() != null) {
                tolXY = (int)Math.ceil(Math.max(
                    Math.abs(params.getStandardDeviations()[2]),
                    Math.abs(params.getStandardDeviations()[3])
                    ));
            } else {
                if (tolXY == 0) {
                    tolXY = 1;
                }
            }
            int dither2 = Math.round(dither * params.getScale());
            if (dither2 == 0) {
                dither2 = 1;
            } else if (dither2 > dither) {
                // large dither makes runtime larger
                dither2 = dither;
            }

            int rotD = Math.round(params.getRotationInDegrees());
            
            double sumSSD = 0;
            double sumDist = 0;
            List<FeatureComparisonStat> stats = new ArrayList<FeatureComparisonStat>();
            List<Double> distances = new ArrayList<Double>();

            for (int ipt1 = 0; ipt1 < c1.size(); ++ipt1) {

                T pt1 = c1.get(ipt1);

                double[] xyTr = transformer.applyTransformation(params,
                    pt1.getX()[pt1.getKMaxIdx()],
                    pt1.getY()[pt1.getKMaxIdx()]);

                Set<Integer> indexes3 = np.findNeighborIndexes(
                    (int) Math.round(xyTr[0]), (int) Math.round(xyTr[1]), tolXY);
                
                if (indexes3 == null || indexes3.isEmpty()) {
                    continue;
                }
         
                FeatureComparisonStat minStat = null;
                    
                T minStatC2 = null;
                    
                for (Integer index2 : indexes3) {
                        
                    T corner2 = c2.get(index2.intValue());
                        
                    FeatureComparisonStat compStat = 
                        featureMatcher.ditherAndRotateForBestLocation2(
                        features1, features2, pt1, corner2, dither2,
                        rotD, rotationTolerance, img1, img2);

                    if (compStat == null || 
                        (compStat.getSumIntensitySqDiff() > compStat.getImg2PointIntensityErr())
                        ) {
                        continue;
                    }
                    if ((minStat == null)
                        || (compStat.getSumIntensitySqDiff() < minStat.getSumIntensitySqDiff())) {
                        minStat = compStat;
                        minStatC2 = corner2;
                    }
                }

                if (minStat == null) {
                    continue;
                }

                double diffX = xyTr[0] - minStatC2.getX()[minStatC2.getKMaxIdx()];
                double diffY = xyTr[1] - minStatC2.getY()[minStatC2.getKMaxIdx()];
                double dist = Math.sqrt(diffX*diffX + diffY*diffY);
                
                //if (minStat.getImg2PointIntensityErr() > lowerLimitSSD) {
                    stats.add(minStat);
                    distances.add(Double.valueOf(dist));
                    sumSSD += minStat.getSumIntensitySqDiff();
                    sumDist += dist;  
                //}
            }

            if (stats.size() < 2) {
                continue;
            }
            
            List<Integer> removedIndexes = MiscStats.filterForDegeneracy(stats);
            if (removedIndexes.size() < distances.size()) {
                for (int i = (removedIndexes.size() - 1); i > -1; --i) {
                    int rmIdx = removedIndexes.get(i);
                    distances.remove(rmIdx);
                }
                sumSSD = 0;
                sumDist = 0;
                for (int i = 0; i < stats.size(); ++i) {
                    sumSSD += stats.get(i).getSumIntensitySqDiff();
                    sumDist += distances.get(i).doubleValue();
                }
            }
            int nEval2 = stats.size();
            
            // distance needs to be adjusted by scale, else the cost prefers
            // small scale solutions
            sumDist /= params.getScale();

            sumSSD  /= (double)nEval2;
            sumDist /= (double)nEval2;
            
            //float cost1 = 1.f/((float)nMaxMatchable*(float)nEval2);
            float cost1 = 1.f/(float)nEval2;
            float cost2 = (float)sumSSD + 0.01f;
            // either distance has to be removed from the normCost or only has a value when > tolXY
            float cost3 = (sumDist < tolXY) ? 0.01f : ((float)sumDist + 0.01f)/(float)tolXY;
            //float normCost = cost1 * cost2 * cost3;
            float normCost = cost1 * cost2 * cost3;
            
            if (normCost < minCost) {
                maxNEval = nEval2;
                bestParams = params;
                minDist = cost3;
                minSSD = cost2;
                bestStats = stats;
                minCost = normCost;
            }
        }

        if (bestStats == null) {
            state = State.FAILED;
            return;
        }
                
        solutionCost =  minCost;
        solutionStats = bestStats;
        state = State.SOLVED;
        solution = bestParams;
        nEval = maxNEval;
    }

    /**
     * get the transformation parameters calculated for the given datasets,
     * but in the frame of the data (if binFactor1 or binFactor2 were not
     * 1, then the this solution cannot be applied to the full frame unbinned
     * image).
     * @return
     */
    public TransformationParameters getSolution() {
        return solution;
    }
    public double getSolutionCost() {
        return solutionCost;
    }

    public int getNEval() {
        return nEval;
    }

    /**
     * @return the solutionStats
     */
    public List<FeatureComparisonStat> getSolutionStats() {
        return solutionStats;
    }

    protected double[] applyTransformation(T corner,
        TransformationParameters params) {

        int x0 = corner.getX()[corner.getKMaxIdx()];

        int y0 = corner.getY()[corner.getKMaxIdx()];

        Transformer transformer = new Transformer();

        double[] xyT = transformer.applyTransformation(params, x0, y0);

        return xyT;
    }

}
