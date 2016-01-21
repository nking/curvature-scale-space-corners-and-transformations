package algorithms.imageProcessing.features;

import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.MiscellaneousCurveHelper;
import algorithms.imageProcessing.util.MiscStats;
import algorithms.util.PairInt;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

/**
 * create lists of singly matched points between 2 images.
 * It uses the criteria that matches are discarded if a point has a second
 * best match whose SSD is within 0.8*SSD or 0.9*SSD of best match.
 * @author nichole
 */
public class CornerMatcher<T extends CornerRegion> {
        
    private final Logger log = Logger.getLogger(this.getClass().getName());

    private List<FeatureComparisonStat> solutionStats = null;
    
    private List<FeatureComparisonStat> rejectedBy2ndBest = new ArrayList<FeatureComparisonStat>();
    
    private List<PairInt> matched1 = null;
    
    private List<PairInt> matched2 = null;

    private enum State {
        INITIALIZED, FAILED, SOLVED
    }

    private State state = null;
    
    private final int dither;

    public CornerMatcher(int dither) {
        this.dither = dither;
    }
    
    private void resetDefaults() {
        state = null;
        solutionStats = null;
        matched1 = null;
        matched2 = null;
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
        final List<T> corners1,final List<T> corners2, GreyscaleImage img1, 
        GreyscaleImage img2, int binFactor1, int binFactor2) {

        if (state != null) {
            resetDefaults();
        }
        
        List<FeatureComparisonStat> stats = new ArrayList<FeatureComparisonStat>();
        
        rejectedBy2ndBest.clear();
        
        FeatureMatcher featureMatcher = new FeatureMatcher();

        for (int i = 0; i < corners1.size(); ++i) {

            T region1 = corners1.get(i);

            FeatureComparisonStat best = null;
            int bestIdx2 = -1;
            FeatureComparisonStat best2nd = null;
            int bestIdx2_2nd = -1;

            for (int j = 0; j < corners2.size(); ++j) {

                T region2 = corners2.get(j);

                FeatureComparisonStat compStat = 
                    featureMatcher.ditherAndRotateForBestLocation2(
                    features1, features2, region1, region2, dither,
                    img1, img2);
                
                if ((compStat == null) ||
                    (compStat.getSumIntensitySqDiff() > compStat.getImg2PointIntensityErr())
                    ) {
                    continue;
                }
                
                if ((best2nd != null) && (compStat.getSumIntensitySqDiff() 
                    >= best2nd.getSumIntensitySqDiff())) {
                    continue;
                }
                
                if (best == null) {
                    best = compStat;
                    bestIdx2 = j;
                } else if (best2nd == null) {
                    if (compStat.getSumIntensitySqDiff() < best.getSumIntensitySqDiff()) {
                        // first becomes second and this becomes first
                        best2nd = best;
                        bestIdx2_2nd = bestIdx2;
                        best = compStat;
                        bestIdx2 = j;
                    } else {
                        best2nd = compStat;
                        bestIdx2_2nd = j;
                    }
                } else {
                    // we know it's better than 2nd best
                    if (compStat.getSumIntensitySqDiff() < best.getSumIntensitySqDiff()) {
                        // first becomes second and this becomes first
                        best2nd = best;
                        bestIdx2_2nd = bestIdx2;
                        best = compStat;
                        bestIdx2 = j;
                    } else {
                        // replaces 2nd best
                        best2nd = compStat;
                        bestIdx2_2nd = j;
                    }
                }
            }
            
            if (best == null) {
                continue;
            }

            if (best2nd == null) {
                stats.add(best);
            } else {
                
                //TODO: the ratio threshold may need to be revised.
                // see Mikolajczyk and Schmid 2005 and the Brown & Lowe paper
                
                float ratio = best.getSumIntensitySqDiff()/best2nd.getSumIntensitySqDiff();
                
                if (ratio < 0.8) {
                    stats.add(best);
                } else {
                    rejectedBy2ndBest.add(best);
                }
            }
        }
        
        MiscStats.filterForDegeneracy(stats);
                
        assignInstanceResults(stats);
        
        return !stats.isEmpty();
    }
    
    /**
     * match corners using the color descriptors
     * @param features1
     * @param features2
     * @param corners1
     * @param corners2
     * @param redImg1
     * @param greenImg1
     * @param blueImg1
     * @param redImg2
     * @param greenImg2
     * @param blueImg2
     * @param binFactor1
     * @param binFactor2
     * @return 
     */
    @SuppressWarnings({"unchecked"})
    public boolean matchCorners(
        final IntensityClrFeatures features1, final IntensityClrFeatures features2,
        final List<T> corners1,final List<T> corners2, 
        GreyscaleImage redImg1, GreyscaleImage greenImg1, GreyscaleImage blueImg1, 
        GreyscaleImage redImg2, GreyscaleImage greenImg2, GreyscaleImage blueImg2,
        int binFactor1, int binFactor2) {

        if (state != null) {
            resetDefaults();
        }
        
        List<FeatureComparisonStat> stats = new ArrayList<FeatureComparisonStat>();
        
        rejectedBy2ndBest.clear();
        
        FeatureMatcher featureMatcher = new FeatureMatcher();

        for (int i = 0; i < corners1.size(); ++i) {

            T region1 = corners1.get(i);

            FeatureComparisonStat best = null;
            int bestIdx2 = -1;
            FeatureComparisonStat best2nd = null;
            int bestIdx2_2nd = -1;

            for (int j = 0; j < corners2.size(); ++j) {

                T region2 = corners2.get(j);

                FeatureComparisonStat compStat = 
                    featureMatcher.findBestMatch(features1, features2, 
                        region1, region2, redImg1, greenImg1, blueImg1, 
                        redImg2, greenImg2, blueImg2);
               
                if ((compStat == null) ||
                    (compStat.getSumIntensitySqDiff() > compStat.getImg2PointIntensityErr())
                    ) {
                    continue;
                }
                
                if ((best2nd != null) && (compStat.getSumIntensitySqDiff() 
                    >= best2nd.getSumIntensitySqDiff())) {
                    continue;
                }
                                
                if (best == null) {
                    best = compStat;
                    bestIdx2 = j;
                } else if (best2nd == null) {
                    if (compStat.getSumIntensitySqDiff() < best.getSumIntensitySqDiff()) {
                        // first becomes second and this becomes first
                        best2nd = best;
                        bestIdx2_2nd = bestIdx2;
                        best = compStat;
                        bestIdx2 = j;
                    } else {
                        best2nd = compStat;
                        bestIdx2_2nd = j;
                    }
                } else {
                    // we know it's better than 2nd best
                    if (compStat.getSumIntensitySqDiff() < best.getSumIntensitySqDiff()) {
                        // first becomes second and this becomes first
                        best2nd = best;
                        bestIdx2_2nd = bestIdx2;
                        best = compStat;
                        bestIdx2 = j;
                    } else {
                        // replaces 2nd best
                        best2nd = compStat;
                        bestIdx2_2nd = j;
                    }
                }
            }
            
            if (best == null) {
                continue;
            }

            if (best2nd == null) {
                stats.add(best);
                //log.info(String.format("==>(%d,%d), (%d,%d)  %.1f  (%.1f)",
                //        best.getImg1Point().getX(), best.getImg1Point().getY(),
                //        best.getImg2Point().getX(), best.getImg2Point().getY(),
                //        best.getSumIntensitySqDiff(), best.getImg2PointIntensityErr()));
            } else {
                
                //TODO: the ratio threshold may need to be revised.
                // see Mikolajczyk and Schmid 2005 and the Brown & Lowe paper
                
                float ratio = best.getSumIntensitySqDiff()/best2nd.getSumIntensitySqDiff();

                //if (ratio < 0.8) {
                if (ratio < 0.9) {
                    stats.add(best);
                    //log.info(String.format("==>(%d,%d), (%d,%d)  %.1f  (%.1f)",
                    //    best.getImg1Point().getX(), best.getImg1Point().getY(),
                    //    best.getImg2Point().getX(), best.getImg2Point().getY(),
                    //    best.getSumIntensitySqDiff(), best.getImg2PointIntensityErr()));
                } else {
                    rejectedBy2ndBest.add(best);
                }
            }
        }
        
        MiscStats.filterForDegeneracy(stats);
                
        assignInstanceResults(stats);
        
        return !stats.isEmpty();
    }
    
    /**
     * match corners using the color descriptors
     * @param features1
     * @param features2
     * @param cornersList1
     * @param cornersList2
     * @param redImg1
     * @param greenImg1
     * @param blueImg1
     * @param redImg2
     * @param greenImg2
     * @param blueImg2
     * @param binFactor1
     * @param binFactor2
     * @return 
     */
    @SuppressWarnings({"unchecked"})
    public boolean matchCornersByBlobs(
        final IntensityClrFeatures features1, final IntensityClrFeatures features2,
        final List<List<T>> cornersList1,final List<List<T>> cornersList2, 
        GreyscaleImage redImg1, GreyscaleImage greenImg1, GreyscaleImage blueImg1, 
        GreyscaleImage redImg2, GreyscaleImage greenImg2, GreyscaleImage blueImg2,
        int binFactor1, int binFactor2) {

        if (state != null) {
            resetDefaults();
        }
        
// Temporary debugging code        
MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
log.info("blobs1:");
for (int i = 0; i < cornersList1.size(); ++i) {
    double[] xyCen = curveHelper.calculateXYCentroids0(cornersList1.get(i));
    log.info(String.format("[%d] (%d,%d)", i, (int)Math.round(xyCen[0]), (int)Math.round(xyCen[1])));
}
log.info("blobs2:");
for (int i = 0; i < cornersList2.size(); ++i) {
    double[] xyCen = curveHelper.calculateXYCentroids0(cornersList2.get(i));
    log.info(String.format("[%d] (%d,%d)", i, (int)Math.round(xyCen[0]), (int)Math.round(xyCen[1])));
}
       
        List<FeatureComparisonStat> stats = new ArrayList<FeatureComparisonStat>();
        
        rejectedBy2ndBest.clear();
        
        FeatureMatcher featureMatcher = new FeatureMatcher();

        for (int i = 0; i < cornersList1.size(); ++i) {
            
            List<T> corners1 = cornersList1.get(i);
            
            double bestCost = Double.MAX_VALUE;
            double bestCostSSD = Double.MAX_VALUE;
            int bestCostIdx2 = -1;
            List<FeatureComparisonStat> bestStats = new ArrayList<FeatureComparisonStat>();
            
            for (int j = 0; j < cornersList2.size(); ++j) {

                List<T> corners2 = cornersList2.get(j);

                // find the best match for each corners1 point, 
                // then filter for degenerate
                // TODO: could consider consistent homology here, but that
                // adds alot of computations.
                
                List<FeatureComparisonStat> bestList1 = new ArrayList<FeatureComparisonStat>();
                
                for (int ii = 0; ii < corners1.size(); ++ii) {
                    
                    T region1 = corners1.get(ii);
            
                    FeatureComparisonStat best1 = null;

                    for (int jj = 0; jj < corners2.size(); ++jj) {

                        T region2 = corners2.get(jj);

                        FeatureComparisonStat compStat = 
                            featureMatcher.findBestMatch(features1, features2, 
                                region1, region2, redImg1, greenImg1, blueImg1, 
                                redImg2, greenImg2, blueImg2);
               
                        if ((compStat == null) ||
                            (compStat.getSumIntensitySqDiff() > compStat.getImg2PointIntensityErr())
                            ) {
                            continue;
                        }
                
                        if (best1 == null) {
                            best1 = compStat;
                        } else if (compStat.getSumIntensitySqDiff() < best1.getSumIntensitySqDiff()) {
                            best1 = compStat;
                        }
                    }
                    if (best1 != null) {
                        bestList1.add(best1);
                    }
                }
                
                MiscStats.filterForDegeneracy(bestList1);
                
                if (bestList1.isEmpty()) {
                    continue;
                }
                
                double cost1 = 1./(double)bestList1.size();
                double cost2 = MiscStats.calculateCombinedIntensityStat(bestList1);
                double normCost = cost1 * cost2;
                
                if (normCost < bestCost) {
                    bestCost = normCost;
                    bestCostSSD = cost2;
                    bestStats = bestList1;
                    bestCostIdx2 = j;
                }
            }
            
            if (bestStats.isEmpty()) {
                continue;
            }
            
            stats.addAll(bestStats);
        }
        
        MiscStats.filterForDegeneracy(stats);
                
        assignInstanceResults(stats);
        
        return !stats.isEmpty();
    }
    
    private void assignInstanceResults(List<FeatureComparisonStat> stats) {
        
        this.solutionStats = stats;
        
        matched1 = new ArrayList<PairInt>();
        
        matched2 = new ArrayList<PairInt>();
        
        for (FeatureComparisonStat stat : stats) {
            
            matched1.add(stat.getImg1Point().copy());
            
            matched2.add(stat.getImg2Point().copy());
        }
    }
    
    /**
     * @return the solutionStats
     */
    public List<FeatureComparisonStat> getSolutionStats() {
        return solutionStats;
    }

    /**
     * @return the matched1
     */
    public List<PairInt> getMatched1() {
        return matched1;
    }

    /**
     * @return the matched2
     */
    public List<PairInt> getMatched2() {
        return matched2;
    }

    public List<FeatureComparisonStat> getRejectedBy2ndBest() {
        return rejectedBy2ndBest;
    }
    
}
