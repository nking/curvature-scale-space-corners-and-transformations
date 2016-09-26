package algorithms.imageProcessing.matching;

import algorithms.compGeometry.RotatedOffsets;
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
import algorithms.imageProcessing.util.AngleUtil;
import algorithms.imageProcessing.util.MiscStats;
import algorithms.misc.MiscDebug;
import algorithms.misc.MiscMath;
import algorithms.search.NearestNeighbor2D;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import gnu.trove.TIntCollection;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.list.TDoubleList;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.ArrayList;
import java.util.HashMap;
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
    
    private final PairIntArray boundaries1;
    
    private final List<PairIntArray> boundaries2;
    
    private final static double twoPi = 2. * Math.PI;

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
        final PairIntArray boundaries1, List<PairIntArray> boundaries2,
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
        this.boundaries1 = boundaries1;
        this.boundaries2 = boundaries2;
        this.medialAxis1 = medialAxis1;
        this.medialAxis2 = medialAxis2;

        debugPrint();
        
        this.features1 = new IntensityClrFeatures2(img1, 5, rotatedOffsets);
        this.features2 = new IntensityClrFeatures2(img2, 5, rotatedOffsets);
        
    }
    
    public EpipolarTransformationFit matchPointsInGroups() {
        throw new UnsupportedOperationException("not yet implemented");
    }

    public EpipolarTransformationFit matchPointsSingly() {

        // key = template point coordinates,  value = statistic of match
        List<FeatureComparisonStat> stats = new ArrayList<FeatureComparisonStat>();

        int n1 = keypoints1.length;
        int n2 = keypoints2.length;

        //NOTE: may need to note when a keypoint is not a boundary point
        //      to dither in rotation when comparing descriptors

        // populate maps
        populateMaps(srchPointMap, index2NNMap, index2Map);

        int[] minMaxXY1 = MiscMath.findMinMaxXY(segmentedCells1);
        NearestNeighbor2D medialAxisNN1 = new NearestNeighbor2D(medialAxis1, 
            minMaxXY1[1] + 20, minMaxXY1[3] + 20);

        for (int i1 = 0; i1 < n1; ++i1) {
            int x = keypoints1[i1][0];
            int y = keypoints1[i1][1];
            PairInt p1 = new PairInt(x, y);
            
            if (!segmentedCells1.contains(p1)) {
                continue;
            }
            
            double d = orientations1.get(i1);
            int rotationInward1 = calculateOrientationInward(x, y, d, 
                medialAxisNN1);
            
            IntensityDescriptor[] hsvDesc = features1.extractIntensityHSV(x, y, 
                rotationInward1);
            
            FeatureComparisonStat best = null;
            int bestIdx2 = -1;
            FeatureComparisonStat best2nd = null;
            int bestIdx2_2nd = -1;
        
            // find best among list2
            for (int i2 = 0; i2 < n2; ++i2) {
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
        
        debugPrint(stats);
        
        if (stats.size() < 7) {
            return null;
        }

        PairIntArray left = new PairIntArray(stats.size());
        PairIntArray right = new PairIntArray(stats.size());

        for (FeatureComparisonStat stat : stats) {
            PairInt p1 = stat.getImg1Point();
            PairInt p2 = stat.getImg2Point();
            left.add(p1.getX(), p1.getY());
            right.add(p2.getX(), p2.getY());
        }

        PairIntArray outputLeftXY = new PairIntArray();
        PairIntArray outputRightXY = new PairIntArray();

        RANSACSolver solver = new RANSACSolver();

        EpipolarTransformationFit fit 
            = solver.calculateEpipolarProjection(
            left, right, outputLeftXY, outputRightXY);

        return fit;
    }

    /**
     * runtime complexity is O(keypoints2.length * segmentedCells2.size()).
     *
     * @param srchPointMap
     * @param index2NNMap
     * @param index2Map
     */
    private void populateMaps(Map<PairInt, TIntIntMap> srchPointMap,
        TIntObjectMap<NearestNeighbor2D> index2NNMap,
        TIntObjectMap<TIntSet> index2Map) {

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

                    int rotationInward = calculateOrientationInward(x, y, d, nn);

                    TIntIntMap rMap = srchPointMap.get(p);
                    if (rMap == null) {
                        rMap = new TIntIntHashMap();
                        srchPointMap.put(p, rMap);
                    }
                    rMap.put(l2, rotationInward);
                }
            }
            
            if (!found) {
                // NOTE: currently, not matching these
                int z = 1;
            }
        }
    }

    protected int calculateOrientationInward(int x, int y, double d,
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
        
        return (int)closestDiffAngle;
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
