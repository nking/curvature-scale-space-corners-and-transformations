package algorithms.imageProcessing.matching;

import algorithms.compGeometry.RotatedOffsets;
import algorithms.imageProcessing.FixedSizeSortedVector;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.features.FeatureComparisonStat;
import algorithms.imageProcessing.features.IntensityClrFeatures2;
import algorithms.imageProcessing.features.IntensityDescriptor;
import algorithms.imageProcessing.features.RANSACEpipolarWithFeaturesSolver;
import algorithms.imageProcessing.features.RANSACSolver;
import algorithms.imageProcessing.transform.EpipolarFeatureTransformationFit;
import algorithms.imageProcessing.transform.EpipolarTransformationFit;
import algorithms.imageProcessing.transform.EuclideanTransformationFit;
import algorithms.imageProcessing.transform.ITransformationFit;
import algorithms.imageProcessing.transform.MatchedPointsTransformationCalculator;
import algorithms.imageProcessing.transform.TransformationParameters;
import algorithms.imageProcessing.util.AngleUtil;
import algorithms.imageProcessing.util.MiscStats;
import algorithms.misc.Misc;
import algorithms.misc.MiscDebug;
import algorithms.misc.MiscMath;
import algorithms.search.NearestNeighbor2D;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import gnu.trove.TIntCollection;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * this class uses keypoints, their orientations,
 * and segmented cells to make a correspondence list.
 *
 * @author nichole
 */
public class SegmentedCellDescriptorMatcher {

    private final ImageExt img1;
    private final ImageExt img2;
    private final int[][] keypoints1;
    private final int[][] keypoints2;
    private final TDoubleList orientations1;
    private final TDoubleList orientations2;
    private final Set<PairInt> segmentedCells1;
    private final List<Set<PairInt>> segmentedCells2;
    private final Set<PairInt> medialAxis1;
    private final List<Set<PairInt>> medialAxis2;

    // depending upon other use, may externalize these, that is require them
    // to be given as arguments
    private final IntensityClrFeatures2 features1;
    private final IntensityClrFeatures2 features2;
    
    // key = srch image coords,
    // value = map with key = index of segmentedCells2,  value=rotation inward
    private final Map<PairInt, TIntIntMap> srchPointMap = new HashMap<PairInt, TIntIntMap>();

    // key=index of segmentedCell2, value = nearest neighbor instances for sets
    private final TIntObjectMap<NearestNeighbor2D> index2NNMap
        = new TIntObjectHashMap<NearestNeighbor2D>();

    // key = index of keypoint2,  value = set of indexes of segmentedCell2
    private final TIntObjectMap<TIntSet> index2Map 
        = new TIntObjectHashMap<TIntSet>();
    
    // key = index of segmentedCell2,  value = set of keypoint2 indexes
    private final TIntObjectMap<TIntSet> segmented2Keypoint2Map 
        = new TIntObjectHashMap<TIntSet>();
            
    private final NearestNeighbor2D medialAxisNN1;
    
    private final int nKeypointsInSegment1;
    
    private final static double twoPi = 2. * Math.PI;

    Set<PairInt> boundary1; 
    List<Set<PairInt>> boundary2;
    
    /**
     * note, RotatedOffsets caches the results of expensive transcendental
     * operations and is used in a shared manner, that is,
     * it is a singleton.  caution is needed in it's use because it
     * contains unsynchronized unguarded cache variables that could be
     * corrupted by multi-threaded use.  The design is meant to keep access
     * to it fast.
     *
     * @param img1
     * @param img2
     * @param keypoints1
     * @param keypoints2
     * @param orientations1
     * @param orientations2
     * @param segmentedCells1
     * @param segmentedCells2
     * @param medialAxis1
     * @param medialAxis2
     * @param rotatedOffsets
     */
    public SegmentedCellDescriptorMatcher(
        final ImageExt img1, final ImageExt img2,
        int[][] keypoints1, int[][] keypoints2,
        final TDoubleList orientations1, final TDoubleList orientations2,
        final Set<PairInt> segmentedCells1,
        final List<Set<PairInt>> segmentedCells2,
        PairIntArray boundaries1, List<PairIntArray> boundaries2,
        final Set<PairInt> medialAxis1, final List<Set<PairInt>> medialAxis2,
        final RotatedOffsets rotatedOffsets) {

        assert(orientations1.size() == keypoints1.length);
        assert(orientations2.size() == keypoints2.length);

        this.img1 = img1;
        this.img2 = img2;
        this.keypoints1 = keypoints1;
        this.keypoints2 = keypoints2;
        this.orientations1 = orientations1;
        this.orientations2 = orientations2;
        this.segmentedCells1 = segmentedCells1;
        this.segmentedCells2 = segmentedCells2;
        this.medialAxis1 = medialAxis1;
        this.medialAxis2 = medialAxis2;

        debugPrint();
        
        this.features1 = new IntensityClrFeatures2(img1, 5, rotatedOffsets);
        this.features2 = new IntensityClrFeatures2(img2, 5, rotatedOffsets);
        
        populateMaps();
        
        int[] minMaxXY1 = MiscMath.findMinMaxXY(segmentedCells1);
        medialAxisNN1 = new NearestNeighbor2D(this.medialAxis1, minMaxXY1[1] + 20, 
            minMaxXY1[3] + 20);
        
        // count the number of keypoints1 in segmentedCell1 for stats
        int n1 = 0;
        for (int i = 0; i < keypoints1.length; ++i) {
            int x = keypoints1[i][0];
            int y = keypoints1[i][1];
            PairInt p1 = new PairInt(x, y);
            if (segmentedCells1.contains(p1)) {
                n1++;
            }
        }
        this.nKeypointsInSegment1 = n1;
        
        segmentedCells1.addAll(Misc.convert(boundaries1));
        for (int i = 0; i < boundaries2.size(); ++i) {
            segmentedCells2.get(i).addAll(Misc.convert(boundaries2.get(i)));
        }
    }
    
    public ITransformationFit matchPointsInGroups() {
        
        /*
        ways to solve using the point location as grouping information without
        making assumptions about projection.
           -- still thinking about this ... 
        */
        
        throw new UnsupportedOperationException("not yet implemented");
    }

    public ITransformationFit matchPointsSingly() {
        
        TIntCollection kp2 = new TIntArrayList(keypoints2.length);
        for (int i2 = 0; i2 < keypoints2.length; ++i2) {
            kp2.add(i2);
        }
        
        return matchPointsSingly(kp2);
    }
    
    protected ITransformationFit matchPointsSingly(TIntCollection kp2) {

        // key = template point coordinates,  value = statistic of match
        List<FeatureComparisonStat> stats = new ArrayList<FeatureComparisonStat>();

        int n1 = keypoints1.length;
        int n2 = keypoints2.length;

        //NOTE: may need to note when a keypoint is not a boundary point
        //      to dither in rotation when comparing descriptors

        for (int i1 = 0; i1 < n1; ++i1) {
            int x = keypoints1[i1][0];
            int y = keypoints1[i1][1];
            PairInt p1 = new PairInt(x, y);
            
            if (!segmentedCells1.contains(p1)) {
                continue;
            }
            
            double d = orientations1.get(i1);
            int rotationInward1[] = calculateOrientationInward(x, y, d, 
                medialAxisNN1);
            
            IntensityDescriptor[] hsvDesc1_1 = features1.extractIntensityHSV(x, y, 
                rotationInward1[0]);
            IntensityDescriptor[] hsvDesc1_2 = features1.extractIntensityHSV(x, y, 
                rotationInward1[1]);
            
            FeatureComparisonStat best = null;
            int bestIdx2 = -1;
            FeatureComparisonStat best2nd = null;
            int bestIdx2_2nd = -1;
            
            for (int ii = 0; ii < 2; ++ii) {
                IntensityDescriptor[] hsvDesc = hsvDesc1_1;
                if (ii == 1) {
                    hsvDesc = hsvDesc1_2;
                }
                
                if (hsvDesc == null) {
                    continue;
                }
                
            // find best among list2
            TIntIterator iter2 = kp2.iterator();
            while (iter2.hasNext()) {
                int i2 = iter2.next();
                int x2 = keypoints2[i2][0];
                int y2 = keypoints2[i2][1];
                PairInt p2 = new PairInt(x2, y2);
                       
                // key = index of segmentedCells2, value = rotation angle inward
                TIntIntMap cellIndexesAndRotationMap = srchPointMap.get(p2);
                
                if (cellIndexesAndRotationMap == null) {
                    continue;
                }
                
                TIntCollection rot2Set = cellIndexesAndRotationMap.valueCollection();
                TIntIterator iter = rot2Set.iterator();
                while (iter.hasNext()) {
                    
                    int rotationInward2 = iter.next();
                    
                    IntensityDescriptor[] hsvDesc2 
                        = features2.extractIntensityHSV(x2, y2, 
                        rotationInward2);
                    
                    if (hsvDesc2 == null) {
                        continue;
                    }
                    
                    FeatureComparisonStat compStat = 
                        IntensityClrFeatures2.calculateHalfStats(
                            hsvDesc[0], hsvDesc[1], hsvDesc[2], x, y, true, 
                            hsvDesc2[0], hsvDesc2[1], hsvDesc2[2], x2, y2, true);
                    
                    if ((compStat == null) || (compStat.getSumIntensitySqDiff() 
                        > compStat.getImg2PointIntensityErr())) {
                        continue;
                    }
                                    
                    if ((best2nd != null) && (compStat.getSumIntensitySqDiff() 
                        >= best2nd.getSumIntensitySqDiff())) {
                        continue;
                    }

                    if (best == null) {
                        best = compStat;
                        bestIdx2 = i2;
                    } else if (best2nd == null) {
                        if (compStat.getSumIntensitySqDiff() < best.getSumIntensitySqDiff()) {
                            // first becomes second and this becomes first
                            best2nd = best;
                            bestIdx2_2nd = bestIdx2;
                            best = compStat;
                            bestIdx2 = i2;
                        } else {
                            best2nd = compStat;
                            bestIdx2_2nd = i2;
                        }
                    } else // we know it's better than 2nd best
                    if (compStat.getSumIntensitySqDiff() < best.getSumIntensitySqDiff()) {
                        // first becomes second and this becomes first
                        best2nd = best;
                        bestIdx2_2nd = bestIdx2;
                        best = compStat;
                        bestIdx2 = i2;
                    } else {
                        // replaces 2nd best
                        best2nd = compStat;
                        bestIdx2_2nd = i2;
                    }
                }
            } // end loop over keypoints2 
            } 
            if (best == null) {
                continue;
            }

            if (best2nd == null) {
                stats.add(best);
            } else {
                //TODO: the ratio threshold may need to be revised.
                // see Mikolajczyk and Schmid 2005 and the Brown & Lowe paper
                float ratio = best.getSumIntensitySqDiff() / best2nd.getSumIntensitySqDiff();
                if (ratio < 0.8) {
                    stats.add(best);
                }
            }
        } // end loop over keypoints1
        
        System.out.println("before filter stats.n=" + stats.size());
        
        MiscStats.filterForDegeneracy(stats);
        
        System.out.println("after filter stats.n=" + stats.size());
        
        if (stats.size() < 2) {
            return null;
        }
        
        debugPrint(stats);
        
        PairIntArray left = new PairIntArray(stats.size());
        PairIntArray right = new PairIntArray(stats.size());

        for (FeatureComparisonStat stat : stats) {
            PairInt p1 = stat.getImg1Point();
            PairInt p2 = stat.getImg2Point();
            left.add(p1.getX(), p1.getY());
            right.add(p2.getX(), p2.getY());
            System.out.println("top=" + stat + 
                " ssd/err=" + 
                stat.getSumIntensitySqDiff()/stat.getImg2PointIntensityErr());
        }
        
        if (stats.size() < 7) {
            MatchedPointsTransformationCalculator tc =
                new MatchedPointsTransformationCalculator();
            float[] weights = new float[stats.size()];
            Arrays.fill(weights, 1.f/(float)stats.size());
            TransformationParameters params = tc.calulateEuclidean(left, 
                right, weights, 0, 0, new float[4]);
            if (params == null) {
                return null;
            }
            
            //TODO: edit a method in tc to return the indexes that were used.
            
            List<Integer> theInlierIndexes = new ArrayList<Integer>();
            List<Double> theErrors = new ArrayList<Double>();
            for (int i = 0; i < stats.size(); ++i) {
                theInlierIndexes.add(Integer.valueOf(-1));
                theErrors.add(Double.valueOf(-1));
            }
            EuclideanTransformationFit fit = new EuclideanTransformationFit(
                params, theInlierIndexes, theErrors, Double.POSITIVE_INFINITY);

            return fit;
        }

        PairIntArray outputLeftXY = new PairIntArray();
        PairIntArray outputRightXY = new PairIntArray();

        RANSACSolver solver = new RANSACSolver();

        EpipolarTransformationFit fit 
            = solver.calculateEpipolarProjection(
            left, right, outputLeftXY, outputRightXY);

        return fit;
    }
    
    public ITransformationFit matchPoints(int[] segmentedCell2Indexes,
        int[] n1N2) {
        
        TIntCollection kp2 = new TIntHashSet();
        for (int idx2 : segmentedCell2Indexes) {
            TIntSet kpSet = segmented2Keypoint2Map.get(idx2);
            if (kpSet != null) {
                kp2.addAll(kpSet);
            }
        }
        
        n1N2[0] = this.nKeypointsInSegment1;
        n1N2[1] = kp2.size();
        
        if (n1N2[1] == 0) {
            return null;
        }
        
        ITransformationFit fit = matchPointsSingly(kp2);
        
        return fit;
    }

    /**
     * runtime complexity is O(keypoints2.length * segmentedCells2.size()).
     *
     */
    private void populateMaps() {

        // variables numbered 1 are those of the template
        
        int n2 = keypoints2.length;

        for (int i2 = 0; i2 < n2; ++i2) {
            int x = keypoints2[i2][0];
            int y = keypoints2[i2][1];
            PairInt p = new PairInt(x, y);

            //range -pi to pi
            double d = orientations2.get(i2);
            if (d < 0) {
                d += twoPi;
            }
            
            boolean found = false;
            
            for (int l2 = 0; l2 < segmentedCells2.size(); ++l2) {

                Set<PairInt> set2 = segmentedCells2.get(l2);

                if (set2.contains(p)) {
                    
                    found = true;

                    // store the index for later use
                    TIntSet kp2Indexes = index2Map.get(i2);
                    if (kp2Indexes == null) {
                        kp2Indexes = new TIntHashSet();
                        index2Map.put(i2, kp2Indexes);
                    }

                    NearestNeighbor2D nn = index2NNMap.get(l2);

                    if (nn == null) {
                        Set<PairInt> medialAxisPt2 = medialAxis2.get(l2);
                        int[] minMaxXY = MiscMath.findMinMaxXY(medialAxisPt2);
                        nn = new NearestNeighbor2D(medialAxisPt2,
                            minMaxXY[1] + 1, minMaxXY[3] + 1);
                        index2NNMap.put(l2, nn);
                    }

                    int rotationInward[] = calculateOrientationInward(x, y, d, nn);

                    TIntIntMap rMap = srchPointMap.get(p);
                    if (rMap == null) {
                        rMap = new TIntIntHashMap();
                        srchPointMap.put(p, rMap);
                    }
                    rMap.put(l2, rotationInward[0]);
                    rMap.put(l2, rotationInward[1]);
                    
                    TIntSet setK = segmented2Keypoint2Map.get(l2);
                    if (setK == null) {
                        setK = new TIntHashSet();
                        segmented2Keypoint2Map.put(l2, setK);
                    }
                    setK.add(i2);
                }
            }
            
            if (!found) {
                // NOTE: currently, not matching these
                int z = 1;
            }
        }
    }

    protected int[] calculateOrientationInward(int x, int y, double d,
        NearestNeighbor2D medialAxisNN) {

        float d0 = (float)(d * 180./Math.PI);
        if (d0 < 0) {
            d0 += 360.f;
        } else if (d0 >= 360.f) {
            d0 -= 360.f;
        }
        assert(d0 >= 0 && d0 < 360.);
        
        float d180 = d0 + 180.f;
        if (d0 < 0) {
            d180 += 360.f;
        } else if (d180 >= 360.f) {
            d180 -= 360.f;
        }
        assert(d180 >= 0 && d180 < 360.);

        //TODO: consider whether to use the polar angle from the (x,y) to
        //   the medial axis as the returned angle

        // d or d +- 180 closer direction to nearest medial axis point
        Set<PairInt> closest = medialAxisNN.findClosest(x, y);
        assert(!closest.isEmpty());

        double closestDiff = Double.MAX_VALUE;
        float closestDiffAngle = -1;
        for (PairInt pc : closest) {

            double angle = Math.atan2(pc.getY() - y, pc.getX() - x);

            if (angle < 0) {
                angle += twoPi;
            }

            float angleDegrees = (float)(angle * 180./Math.PI);
            assert(angleDegrees >= 0 && angleDegrees < 360.);

            double diff0 = Math.abs(AngleUtil.getAngleDifference(
                angleDegrees, d0));

            double diff180 = Math.abs(AngleUtil.getAngleDifference(
                angleDegrees, d180));

            if (diff0 < closestDiff) {
                closestDiff = diff0;
                closestDiffAngle = d0;
            }

            if (diff180 < closestDiff) {
                closestDiff = diff180;
                closestDiffAngle = d180;
            }
        }
        assert(closestDiffAngle >= 0 && closestDiffAngle < 360.);
        
        int closestD = (int)closestDiffAngle;
        int otherD = (closestD == (int)d0) ? (int)d180 : (int)d0;
        return new int[]{closestD, otherD};
    }

    private void debugPrint() {
        
        Image img1Cp = img1.copyImage();
        
        Image img2Cp = img2.copyImage();
        
        for (int i = 0; i < keypoints1.length; ++i) {
            int x = keypoints1[i][0];
            int y = keypoints1[i][1];
            PairInt p = new PairInt(x, y);
            if (!segmentedCells1.contains(p)) {
                continue;
            }
            
            ImageIOHelper.addPointToImage(x, y, img1Cp, 2, 255, 0, 0);
        }
        
        for (int i = 0; i < keypoints2.length; ++i) {
            int x = keypoints2[i][0];
            int y = keypoints2[i][1];
            ImageIOHelper.addPointToImage(x, y, img2Cp, 2, 255, 0, 0);
        }
        
        MiscDebug.writeImage(img1Cp, "_img1_keypoints_");
        MiscDebug.writeImage(img2Cp, "_img2_keypoints_");
    }
    
    private void debugPrint(List<FeatureComparisonStat> stats) {
        
        Image img1Cp = img1.copyImage();
        
        Image img2Cp = img2.copyImage();
        
        for (FeatureComparisonStat stat : stats) {
            PairInt p1 = stat.getImg1Point();
            PairInt p2 = stat.getImg2Point();
            
            ImageIOHelper.addPointToImage(p1.getX(), p1.getY(), img1Cp, 
                2, 255, 0, 0);
            
            ImageIOHelper.addPointToImage(p2.getX(), p2.getY(), img2Cp, 
                2, 255, 0, 0);
        }
        
        MiscDebug.writeImage(img1Cp, "_img1_matched_");
        MiscDebug.writeImage(img2Cp, "_img2_matched_");
    }
    
}
