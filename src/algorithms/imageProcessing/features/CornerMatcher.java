package algorithms.imageProcessing.features;

import algorithms.imageProcessing.CIEChromaticity;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.transform.EpipolarFeatureTransformationFit;
import algorithms.imageProcessing.util.MiscStats;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.logging.Logger;

/**
 * create lists of singly matched points in two different lists.
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
                    featureMatcher.matchDescriptors(features1, features2, 
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
     * match corners using the color descriptors and using a 2nd best filter with
     * threshold of 0.9.
     * 
     * @param features1
     * @param features2
     * @param keypoints1 
     * @param keypoints2 
     * @param img1
     * @param redImg1
     * @param greenImg1
     * @param blueImg1
     * @param img2
     * @param redImg2
     * @param greenImg2
     * @param blueImg2
     * @param binFactor1
     * @param binFactor2
     * @return 
     */
    @SuppressWarnings({"unchecked"})
    public boolean matchPoints(
        final IntensityClrFeatures features1, final IntensityClrFeatures features2,
        final Collection<PairInt> keypoints1, final Collection<PairInt> keypoints2, 
        ImageExt img1, GreyscaleImage redImg1, GreyscaleImage greenImg1, GreyscaleImage blueImg1, 
        ImageExt img2, GreyscaleImage redImg2, GreyscaleImage greenImg2, GreyscaleImage blueImg2,
        int binFactor1, int binFactor2) {

        if (state != null) {
            resetDefaults();
        }
        double deltaELimit = 20;
        
        CIEChromaticity cieC = new CIEChromaticity();
        
        List<FeatureComparisonStat> stats = new ArrayList<FeatureComparisonStat>();
        
        rejectedBy2ndBest.clear();
        
        FeatureMatcher featureMatcher = new FeatureMatcher();

        for (PairInt keypoint1 : keypoints1) {

            FeatureComparisonStat bestStat = null;
            PairInt bestSet2 = null;
            FeatureComparisonStat bestStat_2nd = null;
            PairInt bestSet2_2nd = null;
            
            float[] lab1 = img1.getCIELAB(keypoint1.getX(), keypoint1.getY());

            for (PairInt keypoint2 : keypoints2) {
                
                float[] lab2 = img2.getCIELAB(keypoint2.getX(), keypoint2.getY());
                
                double deltaE = cieC.calcDeltaECIE94(lab1, lab2);
                
                if (deltaE > deltaELimit) {
                    continue;
                }

                FeatureComparisonStat compStat = 
                    featureMatcher.matchDescriptors(features1, features2, 
                        keypoint1.getX(), keypoint1.getY(), 
                        keypoint2.getX(), keypoint2.getY(),
                        redImg1, greenImg1, blueImg1, 
                        redImg2, greenImg2, blueImg2);
               
                if ((compStat == null) ||
                    (compStat.getSumIntensitySqDiff() > compStat.getImg2PointIntensityErr())
                    ) {
                    continue;
                }
                
                if ((bestStat_2nd != null) && (compStat.getSumIntensitySqDiff() 
                    >= bestStat_2nd.getSumIntensitySqDiff())) {
                    continue;
                }
                              
                if (bestStat == null) {
                    bestStat = compStat;
                    bestSet2 = keypoint2;
                } else if (bestStat_2nd == null) {
                    if (compStat.getSumIntensitySqDiff() < bestStat.getSumIntensitySqDiff()) {
                        // first becomes second and this becomes first
                        bestStat_2nd = bestStat;
                        bestSet2_2nd = bestSet2;
                        bestStat = compStat;
                        bestSet2 = keypoint2;
                    } else {
                        bestStat_2nd = compStat;
                        bestSet2_2nd = keypoint2;
                    }
                } else {
                    // we know it's better than 2nd best
                    if (compStat.getSumIntensitySqDiff() < bestStat.getSumIntensitySqDiff()) {
                        // first becomes second and this becomes first
                        bestStat_2nd = bestStat;
                        bestSet2_2nd = bestSet2;
                        bestStat = compStat;
                        bestSet2 = keypoint2;
                    } else {
                        // replaces 2nd best
                        bestStat_2nd = compStat;
                        bestSet2_2nd = keypoint2;
                    } 
                }
            }
            
            if (bestStat == null) {
                continue;
            }

            if (bestStat_2nd == null) {
                stats.add(bestStat);
                //log.info(String.format("==>(%d,%d), (%d,%d)  %.1f  (%.1f)",
                //        best.getImg1Point().getX(), best.getImg1Point().getY(),
                //        best.getImg2Point().getX(), best.getImg2Point().getY(),
                //        best.getSumIntensitySqDiff(), best.getImg2PointIntensityErr()));
            } else {
                
                //TODO: the ratio threshold may need to be revised.
                // see Mikolajczyk and Schmid 2005 and the Brown & Lowe paper
                
                float ratio = bestStat.getSumIntensitySqDiff()/bestStat_2nd.getSumIntensitySqDiff();

                //if (ratio < 0.8) {
                if (ratio < 0.9) {
                    stats.add(bestStat);
                    //log.info(String.format("==>(%d,%d), (%d,%d)  %.1f  (%.1f)",
                    //    best.getImg1Point().getX(), best.getImg1Point().getY(),
                    //    best.getImg2Point().getX(), best.getImg2Point().getY(),
                    //    best.getSumIntensitySqDiff(), best.getImg2PointIntensityErr()));
                } else {
                    rejectedBy2ndBest.add(bestStat);
                }
            }
        }
        
        MiscStats.filterForDegeneracy(stats);
                
        assignInstanceResults(stats);
        
        return !stats.isEmpty();
    }
    
    /**
     * match corners using the color descriptors and using a 2nd best filter with
     * threshold of 0.9.
     * 
     * @param features1
     * @param features2
     * @param keyPointsAndBounds1
     * @param bmaIndex1
     * @param keyPointsAndBounds2
     * @param bmaIndex2
     * @param keypoints1 
     * @param keypoints2 
     * @param img1
     * @param redImg1
     * @param greenImg1
     * @param blueImg1
     * @param img2
     * @param redImg2
     * @param greenImg2
     * @param blueImg2
     * @param binFactor1
     * @param binFactor2
     * @param useHalfDescriptors
     * @return 
     */
    @SuppressWarnings({"unchecked"})
    public boolean matchPoints(
        final IntensityClrFeatures features1, final IntensityClrFeatures features2,
        KeyPointsAndBounds keyPointsAndBounds1, int bmaIndex1,
        KeyPointsAndBounds keyPointsAndBounds2, int bmaIndex2,
        final Collection<PairInt> keypoints1, final Collection<PairInt> keypoints2, 
        ImageExt img1, GreyscaleImage redImg1, GreyscaleImage greenImg1, GreyscaleImage blueImg1, 
        ImageExt img2, GreyscaleImage redImg2, GreyscaleImage greenImg2, GreyscaleImage blueImg2,
        int binFactor1, int binFactor2, boolean useHalfDescriptors) {

        if (state != null) {
            resetDefaults();
        }
        
        double deltaELimit = 20;
        
        CIEChromaticity cieC = new CIEChromaticity();
        
        List<FeatureComparisonStat> stats = new ArrayList<FeatureComparisonStat>();
        
        rejectedBy2ndBest.clear();
        
        FeatureMatcher featureMatcher = new FeatureMatcher();

        for (PairInt keypoint1 : keypoints1) {

            FeatureComparisonStat bestStat = null;
            PairInt bestSet2 = null;
            FeatureComparisonStat bestStat_2nd = null;
            PairInt bestSet2_2nd = null;
            
            float[] lab1 = img1.getCIELAB(keypoint1.getX(), keypoint1.getY());

            for (PairInt keypoint2 : keypoints2) {
                
                float[] lab2 = img2.getCIELAB(keypoint2.getX(), keypoint2.getY());
                
                double deltaE = cieC.calcDeltaECIE94(lab1, lab2);
                
                if (deltaE > deltaELimit) {
                    continue;
                }

                FeatureComparisonStat compStat;
                
                if (useHalfDescriptors) {
                    
                    compStat = featureMatcher.matchHalfDescriptors(
                        features1, features2, 
                        keyPointsAndBounds1, bmaIndex1,
                        keyPointsAndBounds2, bmaIndex2,
                        keypoint1.getX(), keypoint1.getY(), 
                        keypoint2.getX(), keypoint2.getY(),
                        redImg1, greenImg1, blueImg1, 
                        redImg2, greenImg2, blueImg2);
                
                } else {
                    
                    compStat = featureMatcher.matchDescriptors(features1, 
                        features2, keypoint1.getX(), keypoint1.getY(), 
                        keypoint2.getX(), keypoint2.getY(),
                        redImg1, greenImg1, blueImg1, 
                        redImg2, greenImg2, blueImg2);
                }
               
                if ((compStat == null) ||
                    (compStat.getSumIntensitySqDiff() > compStat.getImg2PointIntensityErr())
                    ) {
                    continue;
                }
                
                if ((bestStat_2nd != null) && (compStat.getSumIntensitySqDiff() 
                    >= bestStat_2nd.getSumIntensitySqDiff())) {
                    continue;
                }
                              
                if (bestStat == null) {
                    bestStat = compStat;
                    bestSet2 = keypoint2;
                } else if (bestStat_2nd == null) {
                    if (compStat.getSumIntensitySqDiff() < bestStat.getSumIntensitySqDiff()) {
                        // first becomes second and this becomes first
                        bestStat_2nd = bestStat;
                        bestSet2_2nd = bestSet2;
                        bestStat = compStat;
                        bestSet2 = keypoint2;
                    } else {
                        bestStat_2nd = compStat;
                        bestSet2_2nd = keypoint2;
                    }
                } else {
                    // we know it's better than 2nd best
                    if (compStat.getSumIntensitySqDiff() < bestStat.getSumIntensitySqDiff()) {
                        // first becomes second and this becomes first
                        bestStat_2nd = bestStat;
                        bestSet2_2nd = bestSet2;
                        bestStat = compStat;
                        bestSet2 = keypoint2;
                    } else {
                        // replaces 2nd best
                        bestStat_2nd = compStat;
                        bestSet2_2nd = keypoint2;
                    } 
                }
            }
            
            if (bestStat == null) {
                continue;
            }

            if (bestStat_2nd == null) {
                stats.add(bestStat);
                //log.info(String.format("==>(%d,%d), (%d,%d)  %.1f  (%.1f)",
                //        best.getImg1Point().getX(), best.getImg1Point().getY(),
                //        best.getImg2Point().getX(), best.getImg2Point().getY(),
                //        best.getSumIntensitySqDiff(), best.getImg2PointIntensityErr()));
            } else {
                
                //TODO: the ratio threshold may need to be revised.
                // see Mikolajczyk and Schmid 2005 and the Brown & Lowe paper
                
                float ratio = bestStat.getSumIntensitySqDiff()/bestStat_2nd.getSumIntensitySqDiff();

                //if (ratio < 0.8) {
                if (ratio < 0.9) {
                    stats.add(bestStat);
                    //log.info(String.format("==>(%d,%d), (%d,%d)  %.1f  (%.1f)",
                    //    best.getImg1Point().getX(), best.getImg1Point().getY(),
                    //    best.getImg2Point().getX(), best.getImg2Point().getY(),
                    //    best.getSumIntensitySqDiff(), best.getImg2PointIntensityErr()));
                } else {
                    rejectedBy2ndBest.add(bestStat);
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
     * @param pointList1
     * @param pointList2
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
        final List<List<PairInt>> pointList1, final List<List<PairInt>> pointList2, 
        GreyscaleImage redImg1, GreyscaleImage greenImg1, GreyscaleImage blueImg1, 
        GreyscaleImage redImg2, GreyscaleImage greenImg2, GreyscaleImage blueImg2,
        int binFactor1, int binFactor2) {

        if (state != null) {
            resetDefaults();
        }
/*        
// Temporary debugging code        
MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
log.info("blobs1:");
for (int i = 0; i < cornersList1.size(); ++i) {
    double[] xyCen = curveHelper.calculateXYCentroids0(cornersList1.get(i));
    String str = String.format("[%d] (%d,%d)", i, (int)Math.round(xyCen[0]), (int)Math.round(xyCen[1]));
    if ((Math.abs(xyCen[0] - 68) < 15) && (Math.abs(xyCen[1] - 108) < 15)) {
        str = "**" + str;
    }
    log.info(str);
}
log.info("blobs2:");
for (int i = 0; i < cornersList2.size(); ++i) {
    double[] xyCen = curveHelper.calculateXYCentroids0(cornersList2.get(i));
    String str = String.format("[%d] (%d,%d)", i, (int)Math.round(xyCen[0]), (int)Math.round(xyCen[1]));
    if ((Math.abs(xyCen[0] - 114) < 18) && (Math.abs(xyCen[1] - 118) < 15)) {
        str = "**" + str;
    }
    log.info(str);
}
*/       
        List<FeatureComparisonStat> stats = new ArrayList<FeatureComparisonStat>();
        
        rejectedBy2ndBest.clear();
        
        FeatureMatcher featureMatcher = new FeatureMatcher();

        for (int i = 0; i < pointList1.size(); ++i) {
            
            List<PairInt> points1 = pointList1.get(i);
            
            double bestCost = Double.MAX_VALUE;
            double bestCostSSD = Double.MAX_VALUE;
            int bestCostIdx2 = -1;
            List<FeatureComparisonStat> bestStats = new ArrayList<FeatureComparisonStat>();
            
            for (int j = 0; j < pointList2.size(); ++j) {

                List<PairInt> points2 = pointList2.get(j);

                // find the best match for each corners1 point, 
                // then filter for degenerate
                // TODO: could consider consistent homology here, but that
                // adds alot of computations.
                
                List<FeatureComparisonStat> bestList1 = new ArrayList<FeatureComparisonStat>();
                
                for (int ii = 0; ii < points1.size(); ++ii) {
                    
                    PairInt point1 = points1.get(ii);
            
                    FeatureComparisonStat best1 = null;

                    for (int jj = 0; jj < points2.size(); ++jj) {

                        PairInt point2 = points2.get(jj);

                        FeatureComparisonStat compStat = 
                            /*featureMatcher.findBestMatch(features1, features2, 
                                region1, region2, redImg1, greenImg1, blueImg1, 
                                redImg2, greenImg2, blueImg2, dither);*/
                            featureMatcher.matchDescriptors(features1, features2, 
                                point1.getX(), point1.getY(), point2.getX(), point2.getY(),
                                redImg1, greenImg1, blueImg1, 
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
    
    /**
     * compare each group of points to another to get the best matching
     * features w/o close 2nd bests, then feed those to RANSAC to return
     * a subset of inliers.  Each group from set1 has a best match if any
     * from set2.  If the flag useBipartiteMatching is true, and if one
     * group in set2 is matched more than once to a group in set1, only
     * the best will be kept.
     * @param features1
     * @param features2
     * @param kpab1
     * @param kpab2
     * @param img1
     * @param redImg1
     * @param greenImg1
     * @param blueImg1
     * @param img2
     * @param redImg2
     * @param greenImg2
     * @param blueImg2
     * @param binFactor1
     * @param binFactor2
     * @param useBipartiteMatching
     * @return 
     * @throws java.security.NoSuchAlgorithmException thrown when the algorithm
     * for the SecureRandom instance is not found.
     */
    @SuppressWarnings({"unchecked"})
    public List<List<FeatureComparisonStat>> matchCornersByBlobsAndRANSAC(
        final IntensityClrFeatures features1, final IntensityClrFeatures features2,
        final KeyPointsAndBounds kpab1, final KeyPointsAndBounds kpab2, 
        ImageExt img1, GreyscaleImage redImg1, GreyscaleImage greenImg1, GreyscaleImage blueImg1, 
        ImageExt img2, GreyscaleImage redImg2, GreyscaleImage greenImg2, GreyscaleImage blueImg2,
        int binFactor1, int binFactor2, boolean useBipartiteMatching) 
        throws NoSuchAlgorithmException {

        if (state != null) {
            resetDefaults();
        }
/*
int dbgIdx1 = -1;        
log.info("groups1:");
for (int i = 0; i < kpab1.getBoundingRegions().getPerimeterList().size(); ++i) {
    double[] xyCen = kpab1.getBoundingRegions().getBlobMedialAxes().getOriginalBlobXYCentroid(i);
    String str = String.format("[%d] (%d,%d)", i, (int)Math.round(xyCen[0]), (int)Math.round(xyCen[1]));
    if ((Math.abs(xyCen[0] - 267) < 15) && (Math.abs(xyCen[1] - 79) < 15)) {
        str = "**" + str;
        dbgIdx1 = i;
    }
    log.info(str);
}
int dbgIdx2 = -1;
log.info("groups2:");
for (int i = 0; i < kpab2.getBoundingRegions().getPerimeterList().size(); ++i) {
    double[] xyCen = kpab2.getBoundingRegions().getBlobMedialAxes().getOriginalBlobXYCentroid(i);
    String str = String.format("[%d] (%d,%d)", i, (int)Math.round(xyCen[0]), (int)Math.round(xyCen[1]));
    if ((Math.abs(xyCen[0] - 318) < 15) && (Math.abs(xyCen[1] - 147) < 15)) {
        str = "**" + str;
        dbgIdx2 = i;
    }
    log.info(str);
}
*/
        double deltaELimit = 20;
        
        CIEChromaticity cieC = new CIEChromaticity();
        
        List<Set<PairInt>> keypoints1 = kpab1.getKeyPointGroups();
        
        List<Set<PairInt>> keypoints2 = kpab2.getKeyPointGroups();
        
        List<List<FeatureComparisonStat>> output = new ArrayList<List<FeatureComparisonStat>>();
        Map<Integer, Integer> outputIndexes = new HashMap<Integer, Integer>();
                
        FeatureMatcher featureMatcher = new FeatureMatcher();

        for (int i = 0; i < keypoints1.size(); ++i) {
            
            Set<PairInt> points1 = keypoints1.get(i);
            
            float[] lab1 = kpab1.getBoundingRegions().getBlobMedialAxes().getLABColors(i);
            
            double bestCost = Double.MAX_VALUE;
            double bestCostSSD = Double.MAX_VALUE;
            int bestCostIdx2 = -1;
            List<FeatureComparisonStat> bestStats = null;
            
            for (int j = 0; j < keypoints2.size(); ++j) {
                
                float[] lab2 = kpab2.getBoundingRegions().getBlobMedialAxes().getLABColors(j);
                
                double deltaE = cieC.calcDeltaECIE94(lab1, lab2);
                
                if (deltaE > deltaELimit) {
                    continue;
                }
                
                Set<PairInt> points2 = keypoints2.get(j);
                
                // for reuse and quick lookups
                Map<PairInt, PairInt> inliers1To2 = new HashMap<PairInt, PairInt>();
                Set<PairInt> inliers2 = new HashSet<PairInt>();
                List<FeatureComparisonStat> latestStats = null;
                
                resetDefaults();

                boolean useHalfDescriptors = false;
                
                int nIter = 0;
                int nMaxIter = 1;
                while ((nIter == 0) || (nIter < nMaxIter)) {
                   
                    // construct matching list from best matches 
                    Set<PairInt> set1 = new HashSet<PairInt>(points1);
                    Set<PairInt> set2 = new HashSet<PairInt>(points2);
                    set1.removeAll(inliers1To2.keySet());
                    set2.removeAll(inliers2);
                    
                    boolean matched0 = matchPoints(features1, features2,
                        kpab1, i, kpab2, j,
                        set1, set2, img1, redImg1, greenImg1, blueImg1,
                        img2, redImg2, greenImg2, blueImg2,
                        binFactor1, binFactor2, useHalfDescriptors);
                        
                    if (!matched0) {
                        break;
                    }
                    
                    // input to ransac will use inliers1To2 and results from mathPoints
                    PairIntArray left = new PairIntArray(this.matched1.size() + inliers1To2.size());
                    PairIntArray right = new PairIntArray(this.matched1.size() + inliers1To2.size());
                    for (int ii = 0; ii < this.matched1.size(); ++ii) {
                        PairInt p1 = this.matched1.get(ii);
                        left.add(p1.getX(), p1.getY());
                        PairInt p2 = this.matched2.get(ii);
                        right.add(p2.getX(), p2.getY());
                    }
                    for (Entry<PairInt, PairInt> entry : inliers1To2.entrySet()) {
                        PairInt p1 = entry.getKey();
                        left.add(p1.getX(), p1.getY());
                        PairInt p2 = entry.getValue();
                        right.add(p2.getX(), p2.getY());
                    }
                    
                    if (left.getN() < 7) {
                        break;
                    }
                                        
                    PairIntArray outputLeftXY = new PairIntArray();
                    PairIntArray outputRightXY = new PairIntArray();
                    
                    RANSACEpipolarWithFeaturesSolver solver = new RANSACEpipolarWithFeaturesSolver();
                    
                    EpipolarFeatureTransformationFit fit = solver.calculateEpipolarProjection(
                        left, right, features1, features2, 
                        kpab1, i, kpab2, j, 
                        redImg1, greenImg1, blueImg1,
                        redImg2, greenImg2, blueImg2,
                        outputLeftXY, outputRightXY, useHalfDescriptors);
                    
                    if (fit == null) {
                        break;
                    }
                    
                    if (outputLeftXY.getN() <= inliers2.size()) {
                        break;
                    }
                    inliers1To2.clear();
                    inliers2.clear();

                    for (int ii = 0; ii < outputLeftXY.getN(); ++ii) {
                        PairInt p1 = new PairInt(outputLeftXY.getX(ii), outputLeftXY.getY(ii));
                        PairInt p2 = new PairInt(outputRightXY.getX(ii), outputRightXY.getY(ii));
                        inliers1To2.put(p1, p2);
                        inliers2.add(p2);
                    }
                    
                    latestStats = new ArrayList<FeatureComparisonStat>(fit.getFeatureComparisonStats());
                    
                    nIter++;
                }
                
                if (inliers1To2.isEmpty() || latestStats.isEmpty()) {
                    continue;
                }
                
                double cost1 = 1./(double)latestStats.size();
                double cost2 = MiscStats.calculateCombinedIntensityStat(latestStats);
                double normCost = cost1 * cost2;
                
                if (normCost < bestCost) {
                    bestCost = normCost;
                    bestCostSSD = cost2;
                    bestStats = latestStats;
                    bestCostIdx2 = j;
                }
            }
            
            if (bestStats == null) {
                output.add(new ArrayList<FeatureComparisonStat>());
            } else {
                output.add(bestStats);
                outputIndexes.put(Integer.valueOf(i), Integer.valueOf(bestCostIdx2));

                /*
                    MiscDebug.plotImages(bestStats, redImg1.copyImage(), 
                        redImg2.copyImage(), 2, "_" + i + "_" + bestCostIdx2 + "_");
                    int z = 1;
                */
            }
        }
        
        resetDefaults();
        
        if (useBipartiteMatching) {
            filterForDegeneracy(output, outputIndexes);
        }
        
        removeEmptyItems(output);
        
        return output;
    }
    
    private void filterForDegeneracy(List<List<FeatureComparisonStat>> stats, 
        Map<Integer, Integer> outputIndexes) {

        // reverse mapping of outputIndexes
        Map<Integer, List<Integer>> index2To1Map = new HashMap<Integer, List<Integer>>();
        
        for (Entry<Integer, Integer> entry : outputIndexes.entrySet()) {
            Integer index1 = entry.getKey();
            Integer index2 = entry.getValue();
            
            List<Integer> indexes1 = index2To1Map.get(index2);
            if (indexes1 == null) {
                indexes1 = new ArrayList<Integer>();
                index2To1Map.put(index2, indexes1);
            }
            indexes1.add(index1);
        }
        
        for (Entry<Integer, List<Integer>> entry : index2To1Map.entrySet()) {
            
            List<Integer> indexes1 = entry.getValue();
            if (indexes1.size() < 2) {
                continue;
            }
            Integer index2 = entry.getKey();
            
            double minCost = Double.MAX_VALUE;
            int minCostIdx1 = -1;
            for (Integer index1 : indexes1) {
                int idx1 = index1.intValue();
                double cost = MiscStats.calculateCombinedIntensityStat(
                    stats.get(idx1));
                if (cost < minCost) {
                    minCost = cost;
                    minCostIdx1 = idx1;
                }
            }
            
            for (Integer index1 : indexes1) {
                int idx1 = index1.intValue();
                if (idx1 == minCostIdx1) {
                    continue;
                }
                stats.get(idx1).clear();
            }
        }
        
    }
    
    private void removeEmptyItems(List<List<FeatureComparisonStat>> output) {
        
        for (int i = (output.size() - 1); i > -1; --i) {
            if (output.get(i).isEmpty()) {
                output.remove(i);
            }
        }
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
