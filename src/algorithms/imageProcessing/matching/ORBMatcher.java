package algorithms.imageProcessing.matching;

import algorithms.QuickSort;
import algorithms.compGeometry.FurthestPair;
import algorithms.compGeometry.PerimeterFinder2;
import algorithms.imageProcessing.ColorHistogram;
import algorithms.imageProcessing.FixedSizeSortedVector;
import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.MiscellaneousCurveHelper;
import algorithms.imageProcessing.SIGMA;
import algorithms.imageProcessing.VanishingPoints;
import algorithms.imageProcessing.features.CorrespondenceList;
import algorithms.imageProcessing.features.ORB;
import algorithms.imageProcessing.features.ORB.Descriptors;
import static algorithms.imageProcessing.features.ORB.convertToImage;
import algorithms.imageProcessing.features.RANSACEuclideanSolver;
import algorithms.imageProcessing.features.RANSACSolver;
import algorithms.imageProcessing.matching.PartialShapeMatcher.Result;
import algorithms.imageProcessing.matching.ShapeFinder.ShapeFinderResult;
import algorithms.imageProcessing.transform.EpipolarTransformationFit;
import algorithms.imageProcessing.transform.EpipolarTransformer;
import algorithms.imageProcessing.transform.EuclideanTransformationFit;
import algorithms.imageProcessing.transform.MatchedPointsTransformationCalculator;
import algorithms.imageProcessing.transform.TransformationParameters;
import algorithms.imageProcessing.transform.Transformer;
import algorithms.imageProcessing.util.AngleUtil;
import algorithms.misc.Misc;
import algorithms.misc.MiscDebug;
import algorithms.misc.MiscMath;
import algorithms.search.NearestNeighbor2D;
import algorithms.util.CorrespondencePlotter;
import algorithms.util.OneDIntArray;
import algorithms.util.PairFloatArray;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.QuadInt;
import algorithms.util.TwoDFloatArray;
import algorithms.util.TwoDIntArray;
import algorithms.util.VeryLongBitString;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.iterator.TObjectIntIterator;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.TFloatList;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Stack;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.ejml.simple.SimpleMatrix;

/**
 * a class to hold various methods related to matching
 * the descriptors of ORB.
 * See also ObjectMatcher.
 *
 * @see ORB
 * @see ObjectMatcher
 *
 * @author nichole
 */
public class ORBMatcher {

    // vahnishing points for dataset2
    private VanishingPoints vp2 = null;

    public void setVanishingPointsForSet2(VanishingPoints vp) {
        vp2 = vp;
    }

    // method in progress to replace match0
    public List<CorrespondenceList> match0(ORB orb1, ORB orb2,
        Set<PairInt> labeledPoints1, List<Set<PairInt>> labeledPoints2) {

        /*
        NOTE: bounds are not smoothed within this method, so if
        a partial shape matching is ever appended to it,
        one will need to make separate smoothed bounds for it.
        */
        if (!orb1.getDescrChoice().equals(orb2.getDescrChoice())) {
            throw new IllegalStateException("orbs must contain same kind of descirptors");
        }

        int distTol = 5;

        PairIntArray bounds1 = createOrderedBoundsSansSmoothing(orb1, 
            labeledPoints1);
        if (bounds1.getN() < 7) {
            throw new IllegalStateException("the boundary of object 1 "
                + " must have at least 7 points");
        }
       
        int[] minMaxXYB1 = MiscMath.findMinMaxXY(bounds1);
        NearestNeighbor2D nnb1 = new NearestNeighbor2D(Misc.convert(bounds1),
            minMaxXYB1[1] + distTol + 1,
            minMaxXYB1[3] + distTol + 1);
        TObjectIntMap<PairInt> bounds1IndexMap = new TObjectIntHashMap<PairInt>();
        for (int i = 0; i < bounds1.getN(); ++i) {
            bounds1IndexMap.put(new PairInt(bounds1.getX(i), bounds1.getY(i)), i);
        }

        int nBands = 3;
        if (orb1.getDescrChoice().equals(ORB.DescriptorChoice.HSV)) {
            if (orb1.getDescriptorsH() == null || orb2.getDescriptorsH() == null) {
                throw new IllegalStateException("hsv descriptors must be created first");
            }
        } else if (orb1.getDescrChoice().equals(ORB.DescriptorChoice.ALT)) {
            if (orb1.getDescriptorsListAlt() == null || orb2.getDescriptorsListAlt() == null) {
                throw new IllegalStateException("alt descriptors must be created first");
            }
            nBands = 1;
        } else if (orb1.getDescrChoice().equals(ORB.DescriptorChoice.GREYSCALE)) {
            if (orb1.getDescriptorsList() == null || orb2.getDescriptorsList() == null) {
                throw new IllegalStateException("descriptors must be created first");
            }
            nBands = 1;
        }

        // NOTE: keeping coords in full size reference frames

        TFloatList scales1 = extractScales(orb1.getScalesList());
        TFloatList scales2 = extractScales(orb2.getScalesList());

        if (Math.abs(scales1.get(0) - 1) > 0.01) {
            throw new IllegalArgumentException("logic depends upon first scale" + " level being '1'");
        }
        if (Math.abs(scales2.get(0) - 1) > 0.01) {
            throw new IllegalArgumentException("logic depends upon first scale" + " level being '1'");
        }

        // --- creating maps of segmented cell sizes for first octave
        TIntIntMap sizes2Maps = new TIntIntHashMap();
        for (int i = 0; i < labeledPoints2.size(); ++i) {
            Set<PairInt> set2 = labeledPoints2.get(i);
            if (set2.size() < 7) {
                continue;
            }
            int sz = calculateObjectSize(set2);
            {
                MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
                double[] xyCen = curveHelper.calculateXYCentroids(set2);
                System.out.println("set " + i +
                    " center=" + (int) xyCen[0] + ","
                    + (int) xyCen[1] + " size_full=" + sz);
            }
            sizes2Maps.put(i, sz);
        }

        List<TObjectIntMap<PairInt>> kp1IdxMapList
            = new ArrayList<TObjectIntMap<PairInt>>();
        for (int octave = 0; octave < scales1.size(); ++octave) {
            TObjectIntMap<PairInt> keypoints1IndexMap = new TObjectIntHashMap<PairInt>();
            for (int i = 0; i < orb1.getKeyPoint1List().get(octave).size(); ++i) {
                int x = orb1.getKeyPoint1List().get(octave).get(i);
                int y = orb1.getKeyPoint0List().get(octave).get(i);
                keypoints1IndexMap.put(new PairInt(x, y), i);
            }
            kp1IdxMapList.add(keypoints1IndexMap);
        }
        
        List<TObjectIntMap<PairInt>> kp2IdxMapList
            = new ArrayList<TObjectIntMap<PairInt>>();
        for (int octave = 0; octave < scales2.size(); ++octave) {
            TObjectIntMap<PairInt> keypoints2IndexMap = new TObjectIntHashMap<PairInt>();
            for (int i = 0; i < orb2.getKeyPoint1List().get(octave).size(); ++i) {
                int x = orb2.getKeyPoint1List().get(octave).get(i);
                int y = orb2.getKeyPoint0List().get(octave).get(i);
                keypoints2IndexMap.put(new PairInt(x, y), i);
            }
            kp2IdxMapList.add(keypoints2IndexMap);
        }
        
        // NOTE: some points are on the boundaries of a labeled region
        // and may actually belong to another or both regions.
        // In the search by labeled region below, the transformation is
        // found using only the labeled region keypoints, but then is
        // applied to all keypoints, so the boundary keypoints are still
        // findable at a later stage.
        //  -- for some cases, it might be necessary to consider changing
        //     the keypoint2 maps to hold more than one labeled region for them,
        //     increasing the keypoints2 in the adjacent regions when they
        //     have too few (<3) to be found.
        

        // making a lookup map for keypoint indexes in points2 labeled sets
        List<TIntObjectMap<TIntSet>> labels2KPIdxsList =
            new ArrayList<TIntObjectMap<TIntSet>>();
        for (int octave = 0; octave < scales2.size(); ++octave) {
            TIntObjectMap<TIntSet> labels2KPIdxs = new TIntObjectHashMap<TIntSet>();
            labels2KPIdxsList.add(labels2KPIdxs);
            TObjectIntMap<PairInt> keypoints2IndexMap
                = kp2IdxMapList.get(octave);
            for (int segIdx = 0; segIdx < labeledPoints2.size(); ++segIdx) {
                for (PairInt p : labeledPoints2.get(segIdx)) {
                    if (keypoints2IndexMap.containsKey(p)) {
                        int kp2Idx = keypoints2IndexMap.get(p);
                        TIntSet kpIdxs = labels2KPIdxs.get(segIdx);
                        if (kpIdxs == null) {
                            kpIdxs = new TIntHashSet();
                            labels2KPIdxs.put(segIdx, kpIdxs);
                        }
                        kpIdxs.add(kp2Idx);
                    }
                }
            }
        }

        // -- initialize bounds2MapsList and populte on demand
        TIntObjectMap<PairIntArray> bounds2Maps
            = new TIntObjectHashMap<PairIntArray>();

        List<List<QuadInt>> correspondences = new ArrayList<List<QuadInt>>();
        TDoubleList descCosts = new TDoubleArrayList();
        TIntList nDesc = new TIntArrayList();
        TDoubleList distCosts = new TDoubleArrayList();
        TIntList octs1 = new TIntArrayList();
        TIntList octs2 = new TIntArrayList();
        TIntList segIdxs = new TIntArrayList();

        for (int octave1 = 0; octave1 < scales1.size(); ++octave1) {
        //for (int octave1 = 1; octave1 < 2; ++octave1) {

            float scale1 = scales1.get(octave1);
            
            float sz1 = calculateObjectSize(labeledPoints1)/scale1;

            int nkp1 = orb1.getKeyPoint0List().get(octave1).size();

            TObjectIntMap<PairInt> leftIdxMap = new TObjectIntHashMap<PairInt>();
            PairIntArray left = new PairIntArray(nkp1);
            for (int i = 0; i < nkp1; ++i) {
                int x = orb1.getKeyPoint1List().get(octave1).get(i);
                int y = orb1.getKeyPoint0List().get(octave1).get(i);
                left.add(x, y);
                leftIdxMap.put(new PairInt(x, y), i);
            }
            assert(left.getN() == leftIdxMap.size());

            int[] minMaxXY1 = MiscMath.findMinMaxXY(left);
            NearestNeighbor2D nn1 = new NearestNeighbor2D(Misc.convert(left),
                minMaxXY1[1] + distTol + 1,
                minMaxXY1[3] + distTol + 1);

            for (int octave2 = 0; octave2 < scales2.size(); ++octave2) {
            //for (int octave2 = 0; octave2 < 1; ++octave2) {

                float scale2 = scales2.get(octave2);

                TObjectIntMap<PairInt> keypoints2IndexMap
                    = kp2IdxMapList.get(octave2);
                TIntObjectMap<TIntSet> labels2KPIdxs =
                    labels2KPIdxsList.get(octave2);

                TwoDFloatArray img2 = orb2.getPyramidImages().get(octave2);

                TIntObjectIterator<TIntSet> iter2 = labels2KPIdxs.iterator();
                for (int i2 = 0; i2 < labels2KPIdxs.size(); ++i2) {
                    iter2.advance();

                    int segIdx = iter2.key();

                    TIntSet kp2Idxs = iter2.value();

                    float sz2 = sizes2Maps.get(segIdx)/scale2;

                    if (sz2 == 0 || kp2Idxs.size() < 2) {
                        continue;
                    }

                    if ((sz1 > sz2 && Math.abs(sz1 / sz2) > 1.2) ||
                        (sz2 > sz1 && Math.abs(sz2 / sz1) > 1.2)) {
                        continue;
                    }

 System.out.println("octave1=" + octave1 + " octave2=" + octave2 +
   " sz1=" + sz1 + " sz2=" + sz2 + " segIdx=" + segIdx +
     " nKP2=" + kp2Idxs.size());

                    PairIntArray bounds2 = getOrCreateOrderedBounds(img2,
                        bounds2Maps, segIdx, labeledPoints2.get(segIdx));

                    if (bounds2 == null || bounds2.getN() < 7) {
                        continue;
                    }
                    TObjectIntMap<PairInt> bounds2IdxMap
                        = new TObjectIntHashMap<PairInt>(kp2Idxs.size());
                    for (int j = 0; j < bounds2.getN(); ++j) {
                        bounds2IdxMap.put(new PairInt(bounds2.getX(j),
                            bounds2.getY(j)), j);
                    }

                    TObjectIntMap<PairInt> rightIdxMap
                        = new TObjectIntHashMap<PairInt>(kp2Idxs.size());
                    TIntIterator iter = kp2Idxs.iterator();
                    PairIntArray right = new PairIntArray(kp2Idxs.size());
                    while (iter.hasNext()) {
                        int kpIdx2 = iter.next();
                        int x = orb2.getKeyPoint1List().get(octave2).get(kpIdx2);
                        int y = orb2.getKeyPoint0List().get(octave2).get(kpIdx2);
                        rightIdxMap.put(new PairInt(x, y), right.getN());
                        right.add(x, y);
                    }

                    ORB.Descriptors[] desc1 = getDescriptors(orb1, octave1);
                    ORB.Descriptors[] desc2 = getDescriptors(orb2, octave2);
                    int[][] costD = ORB.calcDescriptorCostMatrix(desc1, desc2);
                    
                    // sorted greedy ordered by incr descr cost
                    PairIntArray m1 = new PairIntArray();
                    PairIntArray m2 = new PairIntArray();
                    
                    // using only the keypoints contained within the current 
                    // labeled region, combinations of point pairs are used to
                    // find the best euclidean transformation.
                    // the transformation is applied to the bounds to pick
                    // up other points used later in the total cost.
                   
                    // 0 = the salukwzde distance, that is the normalized tot cost
                    // 1 = sum of normalized keypoint descriptors
                    //     (not yet normalized by total number possible to match)
                    // 2 = sum of normalized keypoint distances from transformations
                    //     (not yet normalized by total number possible to match)
                    // 3 = number of keypoint matches (not incl boundary that aren't
                    double[] normalizedCost = new double[4];
                    TransformationParameters params = matchGreedy(segIdx,
                        left, right, nBands, costD, nn1,
                        leftIdxMap, keypoints2IndexMap,
                        m1, m2,
                        orb1.getPyramidImages().get(0).a[0].length,
                        orb1.getPyramidImages().get(0).a.length,
                        distTol,
                        bounds1, bounds2, nnb1, bounds1IndexMap, bounds2IdxMap,
                        scale1, scale2, normalizedCost);
                    assert(normalizedCost[3] <= nkp1);
                  
                    //NOTE: below here, have decided not to use epipolar
                    // fit and projections to find the missed high 
                    // projection points,
                    // because the distances from epipolar lines are not
                    // precise enough in some projections to distinguish
                    // nearest neighbors...there are work arounds, but
                    // one can see from the epipolar line distance
                    // tests for android 04 to android 02 that offsets
                    // of more than 10 pixels are within tolerance, so
                    // it is not easy to distinguish between true and
                    // false matches.
                   
                    
                    System.out.println("octave1=" + octave1 + " octave2=" +
                        octave2 + " euclid m1.n=" + m1.getN()
                        + " segIdx=" + segIdx + " c=" + normalizedCost[0]);
                    if (params == null || m1.getN() < 7) {
                        continue;
                    }

                    {// DEBUG
                        String str3 = Integer.toString(segIdx);
                        while (str3.length() < 3) {
                            str3 = "0" + str3;
                        }
                        Image img1 = ORB.convertToImage(
                            orb1.getPyramidImages().get(octave1));
                        for (int ii = 0; ii < m1.getN(); ++ii) {
                            int x = Math.round((float)m1.getX(ii)/scale1);
                            int y = Math.round((float)m1.getY(ii)/scale1);
                            ImageIOHelper.addPointToImage(x, y, img1, 1, 0, 255, 0);
                        }
                        MiscDebug.writeImage(img1, "_TMP1__"
                            + octave1 + "_" + octave2 + "_" + str3 + "_" +
                            MiscDebug.getCurrentTimeFormatted());
                        //====
                        img1 = ORB.convertToImage(
                            orb2.getPyramidImages().get(octave2));
                        for (int ii = 0; ii < m2.getN(); ++ii) {
                            int x = Math.round((float)m2.getX(ii)/scale2);
                            int y = Math.round((float)m2.getY(ii)/scale2);
                            ImageIOHelper.addPointToImage(x, y, img1, 1, 0, 255, 0);
                        }

                        MiscDebug.writeImage(img1, "_TMP2__" +
                            octave1 + "_" + octave2 + "_" +
                            str3 + "_" + MiscDebug.getCurrentTimeFormatted());
                    }
                    
                    {// DEBUG, print matched in green
                        String str3 = Integer.toString(segIdx);
                        while (str3.length() < 3) {
                            str3 = "0" + str3;
                        }
                        Image img1 = ORB.convertToImage(
                            orb1.getPyramidImages().get(octave1));
                        for (int ii = 0; ii < m1.getN(); ++ii) {
                            int x = Math.round((float)m1.getX(ii)/scale1);
                            int y = Math.round((float)m1.getY(ii)/scale1);
                            ImageIOHelper.addPointToImage(x, y, img1, 1, 0, 255, 0);
                        }
                        MiscDebug.writeImage(img1, "_TMP1_"
                            + octave1 + "_" + octave2 + "_" + str3 + "_" +
                            MiscDebug.getCurrentTimeFormatted());
                        //====
                        img1 = ORB.convertToImage(
                            orb2.getPyramidImages().get(octave2));
                        for (int ii = 0; ii < m2.getN(); ++ii) {
                            int x = Math.round((float)m2.getX(ii)/scale2);
                            int y = Math.round((float)m2.getY(ii)/scale2);
                            ImageIOHelper.addPointToImage(x, y, img1, 1, 0, 255, 0);
                        }
                        MiscDebug.writeImage(img1, "_TMP2_" +
                            octave1 + "_" + octave2 + "_" +
                            str3 + "_" + MiscDebug.getCurrentTimeFormatted());
                    }
                    
                    // apply the euclidean transformation to all keypoints
                    // to be able to also include adjacent keypoints and
                    // keypoints in associated oversegmented regions.
                    //output additionalCosts:
                    // 0 = sum of normalized keypoint descriptors
                    // 1 = sum of normalized keypoint distances from transformations
                    // 2 = number of keypoint matches added
                    double[] additionalCosts = addUnmatchedKeypoints(
                        params, m1, m2,
                        nBands, costD, left, leftIdxMap,
                        keypoints2IndexMap,
                        orb1.getPyramidImages().get(0).a[0].length,
                        orb1.getPyramidImages().get(0).a.length,
                        distTol,
                        scale1, scale2
                    );
  
                    if((normalizedCost[3] + additionalCosts[2]) > nkp1) {
                        System.out.println("normalizedCost[3]=" +
                            normalizedCost[3] + " additionalCosts[2]=" +
                            additionalCosts[2] + " nkp1=" + nkp1);
                        int z = 1;
                    }
                    
                    assert(additionalCosts[2] <= nkp1);
                    assert((normalizedCost[3] + additionalCosts[2]) <= nkp1);
                    
                    // --- build combined correspondence and sums

                    int nTot = m1.getN() ;//+ addedKPIdxs.size();
                    List<QuadInt> corres = new ArrayList<QuadInt>(nTot);
                    // for any point in result that is a keypoint,
                    //    add the descriptor cost to totalDescrSum
                    for (int j = 0; j < m1.getN(); ++j) {
                        int x1 = m1.getX(j);
                        int y1 = m1.getY(j);
                        PairInt p1 = new PairInt(x1, y1);

                        int x2 = m2.getX(j);
                        int y2 = m2.getY(j);
                        PairInt p2 = new PairInt(x2, y2);

                        corres.add(new QuadInt(p1, p2));
                    }

                    assert(corres.size() == nTot);
                    
                    //normalizedCost:
                    // 0 = the salukwzde distance, that is the normalized tot cost
                    // 1 = sum of normalized keypoint descriptors
                    //     (not yet normalized by total number possible to match)
                    // 2 = sum of normalized keypoint distances from transformations
                    //     (not yet normalized by total number possible to match)
                    // 3 = number of keypoint matches (not incl boundary that aren't
                    //additionalCosts:
                    // 0 = sum of normalized keypoint descriptors
                    //     (not yet normalized by total number possible to match)
                    // 1 = sum of normalized keypoint distances from transformations
                    //     (not yet normalized by total number possible to match)
                    // 2 = number of keypoint matches added
                    
                    correspondences.add(corres);
                    descCosts.add(normalizedCost[1] + additionalCosts[0]);
                    nDesc.add((int)normalizedCost[3] + (int)additionalCosts[2]);
                    distCosts.add(normalizedCost[2] + additionalCosts[1]);
                    octs1.add(octave1);
                    octs2.add(octave2);
                    segIdxs.add(segIdx);

                }// end loop over octave2's segIdx
            }// end loop over octave2
        }  // end loop over octave1


        int nC = correspondences.size();

        //calculate "salukwzde distances" as costs to rank results
        
        // the objective function OR the logic above
        //    needs to be edited for scale effects from number of keypoints.
        //    the normalizations are different due to differing
        //    numbers of keypoints at an octave.
        //
        //    For example, one smaller set of octaves has a total
        //       lower cost even though it has fewer matched points
        //       and is the false match compared to the true match
        //       in larger set of images.  the smaller f1 is significant.
        // So need to consider a normalization change for this.
        // instead of area, might use the ratio of number of keypoints
        // that are maximally matchable compared to the full size
        // octave... a factor to be applied to increase the total
        //     cost (cost component) of the smaller image set...
        // Or, create an ambiguous set of results because of that
        // scale factor, that then need to be further searched
        // using the aggregated ShapeMatcher for example...
        
        // There is also the observation that if these results contain
        //    an object match of the same segIdx at different octaves,
        //    only the one with the largest number of matched descriptors
        //    should be kept to reduce scale effects and if they have
        //    same number, the smallest octave pair (== largest images)
        //    should be kept.
        
        // --- filter for unique segIdx among results. ---
        // key = segIdx, value = index of corres to keep
        TIntIntMap labelIdxMap = new TIntIntHashMap();
        TIntSet skipIdx = new TIntHashSet();
        for (int i = 0; i < nC; ++i) {
            int segIdx = segIdxs.get(i);
            if (labelIdxMap.containsKey(segIdx)) {
                int cIdx = labelIdxMap.get(segIdx);
                int prevND = nDesc.get(cIdx);
                int nd = nDesc.get(i);
                if (nd < prevND) {
                    skipIdx.add(i);
                    continue;
                } else if (nd == prevND) {
                    double prevDesc = descCosts.get(cIdx);
                    double desc = descCosts.get(i);
                    if (prevDesc < desc) {
                        skipIdx.add(i);
                        continue;
                    } else if (prevDesc > desc) {
                        skipIdx.add(cIdx);
                        labelIdxMap.put(segIdx, i);
                        continue;
                    }
                    // else, both nDesc and desc are same so keep largest
                    //    image set
                    int prevOct1 = octs1.get(cIdx);
                    int oct1 = octs1.get(i);
                    if (prevOct1 < oct1) {
                        skipIdx.add(i);
                        continue;
                    } else if (prevOct1 > oct1) {
                        skipIdx.add(cIdx);
                        labelIdxMap.put(segIdx, i);
                        continue;
                    }
                    int prevOct2 = octs2.get(cIdx);
                    int oct2 = octs2.get(i);
                    if (prevOct2 < oct2) {
                        skipIdx.add(i);
                        continue;
                    } else if (prevOct2 > oct2) {
                        skipIdx.add(cIdx);
                        labelIdxMap.put(segIdx, i);
                        continue;
                    }
                    throw new IllegalStateException("cannot have same "
                        + " octaves and segIdx in results");
                }
            } else {
                labelIdxMap.put(segIdx, i);
            }
        }
        
        // find the smallest octave1 in the set.  this is then the
        // used in normalizing f1
        int minOctave1 = Integer.MAX_VALUE;
        for (int i = 0; i < nC; ++i) {
            if (skipIdx.contains(i)) {
                continue;
            }   
            int octave1 = octs1.get(i);
            if (octave1 < minOctave1) {
                minOctave1 = octave1;
            }
        }
        
        TIntObjectMap<ShapeFinder2.ShapeFinderResult> shapeResults 
            = aggregatedShapeMatch(orb1, orb2,
            labeledPoints1, labeledPoints2,
            octs1, octs2, segIdxs, scales1, scales2);
        
        System.out.println("nShapes=" + shapeResults.size());
        
        float maxDesc = nBands * 256.0f;
        TIntList indexes = new TIntArrayList(nC);
        TFloatList costs = new TFloatArrayList(nC);
        for (int i = 0; i < nC; ++i) {

            if (skipIdx.contains(i)) {
                continue;
            }
            
            int octave1 = octs1.get(i);
            int octave2 = octs2.get(i);
            
            float nKP1Min = orb1.getKeyPoint0List().get(minOctave1).size();
            float nKP1 = orb1.getKeyPoint0List().get(octave1).size();
            
            float nb1 = (float)bounds1.getN();
            
            float nd = nDesc.get(i);
            
            float descCost = (float)descCosts.get(i)/nd;
            float distCost = (float)distCosts.get(i)/nd;
        
            // calculate "fraction of whole" for keypoint descriptors
            final float f1 = 1.f - (nd/nKP1);
            
            // -- a fraction of whole for boundary matching
            float f3 = 1.f - ((float)correspondences.get(i).size() / 
                (nb1 + nKP1));
            
            //calculate the cost of kp descriptors
            float d1 = 1.f - ((nBands * 256.f - descCost)/maxDesc);
            
            float d2 = distCost;
            
            //TODO: revisit these totals
            
            float sd1 = f1 * f1 + d1 * d1;
            
            float sd2 = f1 * f1 + d2 * d2;

            float tot = f1 * f1 + d1 * d1 + d2 * d2 + f3*f3;

            indexes.add(i);
            costs.add(tot);
            
            String str1 = String.format(
                "octave1=%d octave2=%d segIdx=%d nCor=%d nDesc=%d",
                octave1, octave2, segIdxs.get(i), 
                correspondences.get(i).size(), (int)nd);
            
            String str2 = String.format(
                "i=%d descCost=%.2f f1=%.2f distCost=%.2f sd1=%.2f f2=%.2f f3=%.2f tot=%f",
                i, descCost,        f1, sd1, distCost,         f1,     f3, tot);
            
            System.out.println(str1 + " " + str2);
        }

        QuickSort.sortBy1stArg(costs, indexes);

        //System.out.println("costs: " + Arrays.toString(costs));
        //System.out.println("indexes: " + Arrays.toString(costs));
        
        System.out.println("final results=" + costs.size()
            + " bestCost=" + costs.get(0) 
            + " segIdx=" + segIdxs.get(indexes.get(0)));
        
        /*
        if the best has a close 2nd best, might need to use an aggregated
        partial shape matcher, that uses the euclidean transform as a
        constraint (requires new method in PartialShapeMatcher.java).
        */
        
        List<CorrespondenceList> results = new ArrayList<CorrespondenceList>();
        for (int i = 0; i < costs.size(); ++i) {

            int idx = indexes.get(i);

            List<QuadInt> qs = correspondences.get(idx);

            // points are in full reference frame
            results.add(new CorrespondenceList(qs));
        }

        {// DEBUG   a look at the bounds and keypoints tested by octave
            for (int i = 0; i < scales2.size(); ++i) {
                float scale2 = scales2.get(i);
                Image img1 = ORB.convertToImage(orb2.getPyramidImages().get(i));
                // plot keypoints
                TIntObjectMap<TIntSet> labelKP2IndexMap =
                    labels2KPIdxsList.get(i);
                TIntObjectIterator<TIntSet> iter1 = labelKP2IndexMap.iterator();
                for (int ii = 0; ii < labelKP2IndexMap.size(); ++ii) {
                    iter1.advance();
                    int[] clr = ImageIOHelper.getNextRGB(ii);
                    int segIdx = iter1.key();
                    TIntSet kpIdxs = iter1.value();
                    TIntIterator iter2 = kpIdxs.iterator();
                    while (iter2.hasNext()) {
                        int kpIdx = iter2.next();
                        int x = Math.round((float)
                            orb2.getKeyPoint1List().get(i).get(kpIdx)/ scale2);
                        int y = Math.round((float)
                            orb2.getKeyPoint0List().get(i).get(kpIdx)/ scale2);
                        ImageIOHelper.addPointToImage(x, y, img1, 1,
                            clr[0], clr[1], clr[2]);
                    }
                }
                MiscDebug.writeImage(img1, "_TMP3__" + i + "_"
                    + MiscDebug.getCurrentTimeFormatted());

                // plot bounds
                img1 = ORB.convertToImage(orb2.getPyramidImages().get(i));
                //TIntObjectMap<PairIntArray> bounds2Maps
                TIntObjectIterator<PairIntArray> iter3 = bounds2Maps.iterator();
                for (int ii = 0; ii < bounds2Maps.size(); ++ii) {
                    iter3.advance();
                    int[] clr = ImageIOHelper.getNextRGB(ii);
                    int segIdx = iter1.key();
                    PairIntArray bounds2 = iter3.value();
                    for (int idx = 0; idx < bounds2.getN(); ++idx) {
                        int x = Math.round((float) bounds2.getX(idx)/ scale2);
                        int y = Math.round((float) bounds2.getY(idx)/ scale2);
                        ImageIOHelper.addPointToImage(x, y, img1, 0,
                            clr[0], clr[1], clr[2]);
                    }
                }
                MiscDebug.writeImage(img1, "_TMP4__"
                    + MiscDebug.getCurrentTimeFormatted());
            }
        }

        return results;
    }

    /**
     *
     * NOT READY FOR USE yet.
     *
     * needs the orbs to contain the theta pyramidal images.
     * add usage here.
     *
     * @param orb1
     * @param orb2
     * @param labeledPoints1
     * @param labeledPoints2
     * @return
     */
    public List<CorrespondenceList> matchSmall(ORB orb1, ORB orb2, Set<PairInt> labeledPoints1, List<Set<PairInt>> labeledPoints2) {

        TFloatList scales1 = extractScales(orb1.getScalesList());
        TFloatList scales2 = extractScales(orb2.getScalesList());

        SIGMA sigma = SIGMA.ZEROPOINTFIVE;

        ImageProcessor imageProcessor = new ImageProcessor();

        ColorHistogram cHist = new ColorHistogram();

        int templateSize = calculateObjectSize(labeledPoints1);

        TIntObjectMap<Set<PairInt>> labeledPoints1Lists = new TIntObjectHashMap<Set<PairInt>>();

        // key = octave number, value = histograms of cie luv
        TIntObjectMap<TwoDIntArray> ch1s = new TIntObjectHashMap<TwoDIntArray>();

        // key = octave number, value = ordered boundaries of sets
        TIntObjectMap<PairIntArray> labeledBoundaries1 = new TIntObjectHashMap<PairIntArray>();

        for (int octave1 = 0; octave1 < scales1.size(); ++octave1) {
            float scale1 = scales1.get(octave1);
            Set<PairInt> set1 = new HashSet<PairInt>();
            for (PairInt p : labeledPoints1) {
                PairInt p1 = new PairInt(Math.round((float) p.getX() / scale1), Math.round((float) p.getY() / scale1));
                set1.add(p1);
            }
            labeledPoints1Lists.put(octave1, set1);
            Image img = ORB.convertToImage(orb1.getPyramidImages().get(octave1));
            int[][] ch = cHist.histogramCIELUV(img, set1);
            ch1s.put(octave1, new TwoDIntArray(ch));
            PairIntArray bounds = imageProcessor.extractSmoothedOrderedBoundary(
                new HashSet(set1), sigma, img.getWidth(), img.getHeight());
            labeledBoundaries1.put(octave1, bounds);
        }
        int dp = 1; //2;
        float intersectionLimit = 0.5F;

        // key = octave number, value = list of labeled sets
        TIntObjectMap<List<Set<PairInt>>> labeledPoints2Lists = new TIntObjectHashMap<List<Set<PairInt>>>();

        // key = octave number, value = list of histograms of cie lab theta
        TIntObjectMap<List<TwoDIntArray>> ch2Lists
            = new TIntObjectHashMap<List<TwoDIntArray>>();

        // key = octave number, value = list of ordered points in labeled set
        TIntObjectMap<List<PairIntArray>> labeledBoundaries2Lists = new TIntObjectHashMap<List<PairIntArray>>();

        for (int k = 0; k < labeledPoints2.size(); ++k) {
            Set<PairInt> set = labeledPoints2.get(k);
            if (set.size() < 7) {
                // NOTE: this means that subsequent datasets2 will not be
                //   lists having same indexes as labeledPoints2
                continue;
            }

            assert(Math.abs(scales2.get(0) - 1) < 0.02);
            PairIntArray bounds = imageProcessor.extractSmoothedOrderedBoundary(
                new HashSet(set), sigma,
                orb2.getPyramidImages().get(0).a[0].length,
                orb2.getPyramidImages().get(0).a.length);

            for (int octave2 = 0; octave2 < scales2.size(); ++octave2) {

                float scale2 = scales2.get(octave2);

                Image img = ORB.convertToImage(
                    orb2.getPyramidImages().get(octave2));
                int w2 = img.getWidth();
                int h2 = img.getHeight();

                Set<PairInt> set2 = new HashSet<PairInt>();
                for (PairInt p : set) {
                    int x = Math.round((float) p.getX() / scale2);
                    int y = Math.round((float) p.getY() / scale2);
                    if (x == w2) {
                        x = w2 - 1;
                    }
                    if (y == h2) {
                        y = h2 - 1;
                    }
                    PairInt p2 = new PairInt(x, y);
                    set2.add(p2);
                }
                List<Set<PairInt>> list2 = labeledPoints2Lists.get(octave2);
                if (list2 == null) {
                    list2 = new ArrayList<Set<PairInt>>();
                    labeledPoints2Lists.put(octave2, list2);
                }
                list2.add(set2);

                // create histograms for later comparison w/ template at
                // different scales
                int[][] ch = cHist.histogramCIELUV(img, set2);
                List<TwoDIntArray> ch2List = ch2Lists.get(octave2);
                if (ch2List == null) {
                    ch2List = new ArrayList<TwoDIntArray>();
                    ch2Lists.put(octave2, ch2List);
                }
                ch2List.add(new TwoDIntArray(ch));

                List<PairIntArray> list3 = labeledBoundaries2Lists.get(octave2);
                if (list3 == null) {
                    list3 = new ArrayList<PairIntArray>();
                    labeledBoundaries2Lists.put(octave2, list3);
                }
                PairIntArray bounds2 = reduceBounds(bounds, scale2);
                list3.add(bounds2);

                assert(labeledBoundaries2Lists.get(octave2).size() ==
                       labeledPoints2Lists.get(octave2).size());
                assert(labeledBoundaries2Lists.get(octave2).size() ==
                       ch2Lists.get(octave2).size());
            }
        }

        // populated on demand,  key=octave, key=segmented cell, value=size
        TObjectIntMap<PairInt> size2Map = new TObjectIntHashMap<PairInt>();

        // -- compare sets over octaves, first by color histogram intersection,
        //    then by partial shape matcher
        // delaying evaluation of results until end in order to get the
        // maximum chord differerence sum, needed for Salukwzde distance.
        // for each i, list of Results, chordDiffSums, bounds1, bounds2
        //             bundling Results and bounds into an object
        TIntObjectMap<List<PObject>> resultsMap = new TIntObjectHashMap<List<PObject>>();
        TIntObjectMap<TDoubleList> chordDiffSumsMap = new TIntObjectHashMap<TDoubleList>();
        TIntObjectMap<TFloatList> intersectionsMap = new TIntObjectHashMap<TFloatList>();
        double maxDiffChordSum = Double.MIN_VALUE;
        double maxAvgDiffChord = Double.MIN_VALUE;
        double maxAvgDist = Double.MIN_VALUE;
        for (int i = 0; i < scales1.size(); ++i) {
        //for (int i = 2; i < 3; ++i) {

            float scale1 = scales1.get(i);

            int[][] ch1 = ch1s.get(i).a;
            //Set<PairInt> templateSet = labeledPoints1Lists.get(i);
            PairIntArray bounds1 = labeledBoundaries1.get(i);
            float sz1 = calculateObjectSize(bounds1);

            List<PObject> results = new ArrayList<PObject>();
            TDoubleList chordDiffSums = new TDoubleArrayList();
            TFloatList intersections = new TFloatArrayList();

            for (int j = 0; j < scales2.size(); ++j) {
            //for (int j = 0; j < 1; ++j) {

                float scale2 = scales2.get(j);

                List<TwoDIntArray> listOfCH2s = ch2Lists.get(j);
                if (listOfCH2s == null) {
                    continue;
                }
                List<Set<PairInt>> listOfSets2 = labeledPoints2Lists.get(j);
                List<PairIntArray> listOfBounds2 = labeledBoundaries2Lists.get(j);

                for (int k = 0; k < listOfCH2s.size(); ++k) {

                    PairIntArray bounds2 = listOfBounds2.get(k);

                    PairInt octLabelKey = new PairInt(j, k);
                    float sz2;
                    if (size2Map.containsKey(octLabelKey)) {
                        sz2 = size2Map.get(octLabelKey);
                    } else {
                        sz2 = calculateObjectSize(bounds2);
                    }
                    if (sz2 == 0) {
                        continue;
                    }
                    if ((sz1 > sz2 && Math.abs((float)sz1 / (float)sz2) > 1.15) ||
                        (sz2 > sz1 && Math.abs((float)sz2 / (float)sz1) > 1.15)) {
                        continue;
                    }

                    int[][] ch2 = listOfCH2s.get(k).a;
                    float intersection = cHist.intersection(ch1, ch2);
                    if (intersection < intersectionLimit) {
                        continue;
                    }

                    System.out.println("p2=" +
                        listOfSets2.get(k).iterator().next()
                        + " sz1=" + sz1 + " sz2=" + sz2
                        + " nSet=" + listOfSets2.get(k).size());

                    PartialShapeMatcher matcher = new PartialShapeMatcher();
                    matcher.overrideSamplingDistance(dp);
                    //matcher.setToDebug();
                    //matcher.setToUseSameNumberOfPoints();
                    PartialShapeMatcher.Result r = matcher.match(bounds1, bounds2);
                    if (r == null) {
                        continue;
                    }

                    //NOTE: to increase the abilit to find projected objects
                    //  that have euclidean poses and skew, might consider
                    //  fast ways to approximate an affine and evaluate it
                    //  after the euclidean solution here.
                    //  affine transformations leave parallel lines in the
                    //    transformed space so could look for that in the
                    //    unmatched portion:
                    //  for example, if half of the object is matched,
                    //  could determine the distance of the matched to the
                    //  unmatched and use that with knowledge of the
                    //  euclidean expected distance to approximate a shear.
                    //  for the evaluations to remain easy to compare results
                    //  with other results, would not want to allow too much
                    //  shear...
                    //  in order to add the chord differences, this additional
                    //  calculation needs to be handled in the
                    //  partial shape matcher (but can be left until the end)

                    double c = r.getChordDiffSum();

                    results.add(new PObject(r, bounds1, bounds2, scale1, scale2));
                    chordDiffSums.add(r.getChordDiffSum());
                    intersections.add(intersection);

                    if (r.getChordDiffSum() > maxDiffChordSum) {
                        maxDiffChordSum = r.getChordDiffSum();
                    }
                    double avgCD = r.getChordDiffSum() / (double) r.getNumberOfMatches();
                    if (avgCD > maxAvgDiffChord) {
                        maxAvgDiffChord = avgCD;
                    }
                    double avgDist = r.getDistSum() / (double) r.getNumberOfMatches();
                    if (avgDist > maxAvgDist) {
                        maxAvgDist = avgDist;
                    }

                    System.out.println(String.format(
                        "%d %d p in set=%s  shape matcher c=%.2f np=%d  inter=%.2f  dist=%.2f avgDist=%.2f",
                        i, j, listOfSets2.get(k).iterator().next().toString(),
                        (float) c, r.getNumberOfMatches(), (float) intersection,
                        (float) r.getDistSum(), (float) avgDist));
                    try {
                        CorrespondencePlotter plotter = new CorrespondencePlotter(bounds1, bounds2);
                        for (int ii = 0; ii < r.getNumberOfMatches(); ++ii) {
                            int idx1 = r.getIdx1(ii);
                            int idx2 = r.getIdx2(ii);
                            int x1 = bounds1.getX(idx1);
                            int y1 = bounds1.getY(idx1);
                            int x2 = bounds2.getX(idx2);
                            int y2 = bounds2.getY(idx2);
                            if ((ii % 4) == 0) {
                                plotter.drawLineInAlternatingColors(x1, y1, x2, y2, 0);
                            }
                        }
                        String strI = Integer.toString(i);
                        while (strI.length() < 2) {
                            strI = "0" + strI;
                        }
                        String strJ = Integer.toString(j);
                        while (strJ.length() < 2) {
                            strJ = "0" + strJ;
                        }
                        String strK = Integer.toString(k);
                        while (strK.length() < 2) {
                            strK = "0" + strK;
                        }
                        String str = strI + strJ + strK;
                        String filePath = plotter.writeImage("_andr_" + str);
                    } catch (Throwable t) {
                    }
                } //end loop over k labeled sets of dataset 2
            } // end loop over j datasets 2
            if (!results.isEmpty()) {
                resultsMap.put(i, results);
                chordDiffSumsMap.put(i, chordDiffSums);
                intersectionsMap.put(i, intersections);
            }
        } // end loop over i dataset 1
        // calculate the Salukwdze distances

        /*
        for each i, need the max chord diff sum, nPoints in bound1, and best Results
         */

        double minSD = Double.MAX_VALUE;
        int minSDI = -1;

        TIntObjectIterator<List<PObject>> iter = resultsMap.iterator();

        for (int i = 0; i < resultsMap.size(); ++i) {

            iter.advance();

            int idx = iter.key();

            //double maxDiffChordSum = chordDiffSumsMap.get(idx).max();
            double minCost = Double.MAX_VALUE;
            int minCostIdx = -1;
            List<PObject> resultsList = resultsMap.get(idx);

            for (int j = 0; j < resultsList.size(); ++j) {

                PObject obj = resultsList.get(j);

                float costIntersection = 1.0F - intersectionsMap.get(idx).get(j);

                PartialShapeMatcher.Result r = obj.r;

                int nb1 = Math.round((float) obj.bounds1.getN() / (float) dp);

                float np = r.getNumberOfMatches();

                float countComp = 1.0F - (np / (float) nb1);
                float countCompSq = countComp * countComp;
                double chordComp = ((float) r.getChordDiffSum() / np) / maxAvgDiffChord;
                double chordCompSq = chordComp * chordComp;
                double avgDist = r.getDistSum() / np;
                double distComp = avgDist / maxAvgDist;
                double distCompSq = distComp * distComp;

                // Salukwzde uses square sums
                //double sd = r.calculateSalukwdzeDistanceSquared(
                //    maxDiffChordSum, nb1);

                // TODO: consider formal analysis of dependencies and hence
                //       error terms:
                //double sd = chordCompSq*countCompSq
                //    + distCompSq*countCompSq;

                //NOTE: The coverage of the matches is currently
                // approximated as simply numberMatched/maxNumberMatchable,
                // but a term representing the spatial distribution appears
                // to be necessary also.
                // will try largestNumberGap/maxNumberMatchable.
                // TODO: need to improve this in detail later
                int lGap = maxNumberOfGaps(obj.bounds1, r)/dp;
                float gCountComp = (float)lGap/(float)nb1;

                //double sd = chordCompSq + countCompSq + distCompSq;
                double sd = chordComp + countComp + gCountComp + distComp
                    + costIntersection;

                if (sd < minCost) {
                    minCost = sd;
                    minCostIdx = j;
                }
                if (sd < minSD) {
                    minSD = sd;
                    minSDI = idx;
                }
                System.out.println("sd=" + sd + " n1="
                    + obj.bounds1.getN() + " n2=" + obj.bounds2.getN()
                    + " origN1=" + r.getOriginalN1()
                    + " nMatches=" + r.getNumberOfMatches()
                    + String.format(
                    " chord=%.2f count=%.2f spatial=%.2f dist=%.2f inter=%.2f",
                    (float)chordComp, (float)countComp,
                    (float)gCountComp, (float)distComp,
                    (float)costIntersection)
                );
            }
            assert (minCostIdx > -1);

            TDoubleList cList = chordDiffSumsMap.get(idx);

            TFloatList iList = intersectionsMap.get(idx);

            for (int j = resultsList.size() - 1; j > -1; --j) {
                if (j != minCostIdx) {
                    resultsList.remove(j);
                    cList.removeAt(j);
                    iList.removeAt(j);
                }
            }
        }
        if (resultsMap.size() > 1) {
            // TODO: build a test for this.
            //    possibly need to transform results to same reference
            //    frame to compare.
            //    using best SD for now
            TIntSet rm = new TIntHashSet();
            iter = resultsMap.iterator();
            for (int i = 0; i < resultsMap.size(); ++i) {
                iter.advance();
                int idx = iter.key();
                if (idx != minSDI) {
                    rm.add(idx);
                }
            }
            TIntIterator iter2 = rm.iterator();
            while (iter2.hasNext()) {
                int idx = iter2.next();
                resultsMap.remove(idx);
            }
        }

        List<CorrespondenceList> topResults = new ArrayList<CorrespondenceList>();

        iter = resultsMap.iterator();

        for (int i = 0; i < resultsMap.size(); ++i) {
            iter.advance();
            int idx = iter.key();
            List<PObject> resultsList = resultsMap.get(idx);
            assert (resultsList.size() == 1);
            PObject obj = resultsList.get(0);
            int n = obj.r.getNumberOfMatches();
            if (obj.r.getTransformationParameters() == null) {
                continue;
            }
            PairInt[] m1 = new PairInt[n];
            PairInt[] m2 = new PairInt[n];
            float scale1 = obj.scale1;
            float scale2 = obj.scale2;
            for (int ii = 0; ii < n; ++ii) {
                int idx1 = obj.r.getIdx1(ii);
                int idx2 = obj.r.getIdx2(ii);
                int x1 = Math.round(obj.bounds1.getX(idx1) * scale1);
                int y1 = Math.round(obj.bounds1.getY(idx1) * scale1);
                int x2 = Math.round(obj.bounds2.getX(idx2) * scale2);
                int y2 = Math.round(obj.bounds2.getY(idx2) * scale2);
                m1[ii] = new PairInt(x1, y1);
                m2[ii] = new PairInt(x2, y2);
            }
            CorrespondenceList cor
                = new CorrespondenceList(obj.r.getTransformationParameters(), m1, m2);
            topResults.add(cor);
        }
        return topResults;
    }

    /**
     *
     * NOT READY FOR USE yet.
     *
     * searchs among aggregated adjacent labeled points to find best
     * fitting shape and color object where template is
     * dataset 1 and the searchable is dataset 2.
     *
     * @param orb1
     * @param orb2
     * @param labeledPoints1
     * @param labeledPoints2
     * @return
     */
    public List<CorrespondenceList> matchAggregatedShape(
        ORB orb1, ORB orb2, Set<PairInt> labeledPoints1,
        List<Set<PairInt>> labeledPoints2) {

        TFloatList scales1 = extractScales(orb1.getScalesList());
        TFloatList scales2 = extractScales(orb2.getScalesList());

        if (Math.abs(scales1.get(0) - 1) > 0.01) {
            throw new IllegalArgumentException("logic depends upon first scale" + " level being '1'");
        }
        if (Math.abs(scales2.get(0) - 1) > 0.01) {
            throw new IllegalArgumentException("logic depends upon first scale" + " level being '1'");
        }

        SIGMA sigma = SIGMA.ZEROPOINTFIVE;

        ImageProcessor imageProcessor = new ImageProcessor();

        ColorHistogram cHist = new ColorHistogram();

        int templateSize = calculateObjectSize(labeledPoints1);

        TIntObjectMap<Set<PairInt>> labeledPoints1Lists = new TIntObjectHashMap<Set<PairInt>>();

        // key = octave number, value = histograms of cie cie luv
        TIntObjectMap<TwoDIntArray> ch1s = new TIntObjectHashMap<TwoDIntArray>();

        // key = octave number, value = ordered boundaries of sets
        TIntObjectMap<PairIntArray> labeledBoundaries1 = new TIntObjectHashMap<PairIntArray>();

        for (int octave1 = 0; octave1 < scales1.size(); ++octave1) {
            float scale1 = scales1.get(octave1);
            Set<PairInt> set1 = new HashSet<PairInt>();
            for (PairInt p : labeledPoints1) {
                PairInt p1 = new PairInt(Math.round((float) p.getX() / scale1), Math.round((float) p.getY() / scale1));
                set1.add(p1);
            }
            labeledPoints1Lists.put(octave1, set1);
            Image img = ORB.convertToImage(orb1.getPyramidImages().get(octave1));
            int[][] ch = cHist.histogramCIELUV(img, set1);
            ch1s.put(octave1, new TwoDIntArray(ch));
            PairIntArray bounds = imageProcessor.extractSmoothedOrderedBoundary(
                new HashSet(set1), sigma, img.getWidth(), img.getHeight());
            labeledBoundaries1.put(octave1, bounds);
        }
        int dp = 1; //2;
        float intersectionLimit = 0.5F;

        // key = octave number, value = list of labeled sets
        TIntObjectMap<List<Set<PairInt>>> labeledPoints2Lists = new TIntObjectHashMap<List<Set<PairInt>>>();

        // key = octave number, value = list of histograms of cie lab theta
        TIntObjectMap<List<TwoDIntArray>> ch2Lists
            = new TIntObjectHashMap<List<TwoDIntArray>>();

        for (int k = 0; k < labeledPoints2.size(); ++k) {
            Set<PairInt> set = labeledPoints2.get(k);
            if (set.size() < 7) {
                // NOTE: this means that subsequent datasets2 will not be
                //   lists having same indexes as labeledPoints2
                continue;
            }

            assert(Math.abs(scales2.get(0) - 1) < 0.02);
            PairIntArray bounds = imageProcessor.extractSmoothedOrderedBoundary(
                new HashSet(set), sigma,
                orb2.getPyramidImages().get(0).a[0].length,
                orb2.getPyramidImages().get(0).a.length);

            for (int octave2 = 0; octave2 < scales2.size(); ++octave2) {

                float scale2 = scales2.get(octave2);

                Image img = ORB.convertToImage(
                    orb2.getPyramidImages().get(octave2));
                int w2 = img.getWidth();
                int h2 = img.getHeight();

                Set<PairInt> set2 = new HashSet<PairInt>();
                for (PairInt p : set) {
                    int x = Math.round((float) p.getX() / scale2);
                    int y = Math.round((float) p.getY() / scale2);
                    if (x == w2) {
                        x = w2 - 1;
                    }
                    if (y == h2) {
                        y = h2 - 1;
                    }
                    PairInt p2 = new PairInt(x, y);
                    set2.add(p2);
                }
                List<Set<PairInt>> list2 = labeledPoints2Lists.get(octave2);
                if (list2 == null) {
                    list2 = new ArrayList<Set<PairInt>>();
                    labeledPoints2Lists.put(octave2, list2);
                }
                list2.add(set2);

                // create histograms for later comparison w/ template at
                // different scales
                int[][] ch = cHist.histogramCIELUV(img, set2);
                List<TwoDIntArray> ch2List = ch2Lists.get(octave2);
                if (ch2List == null) {
                    ch2List = new ArrayList<TwoDIntArray>();
                    ch2Lists.put(octave2, ch2List);
                }
                ch2List.add(new TwoDIntArray(ch));

                assert(labeledPoints2Lists.get(octave2).size() ==
                    ch2Lists.get(octave2).size());
            }
        }

        // populated on demand,  key=octave, key=segmented cell, value=size
        TObjectIntMap<PairInt> size2Map = new TObjectIntHashMap<PairInt>();

        // -- compare sets over octaves:
        //    aggregated search of adjacent labeled cells to compare their combined
        //       properties of color histogram and shape to the template.
        //
        // delaying evaluation of results until end in order to get the
        // maximum chord differerence sum, needed for Salukwzde distance.
        // for each i, list of Results, chordDiffSums, bounds1, bounds2
        //             bundling Results and bounds into an object
        TIntObjectMap<List<PObject>> resultsMap = new TIntObjectHashMap<List<PObject>>();
        TIntObjectMap<TDoubleList> chordDiffSumsMap = new TIntObjectHashMap<TDoubleList>();
        TIntObjectMap<TFloatList> intersectionsMap = new TIntObjectHashMap<TFloatList>();
        double maxDiffChordSum = Double.MIN_VALUE;
        double maxAvgDiffChord = Double.MIN_VALUE;
        double maxAvgDist = Double.MIN_VALUE;

        // maps to reuse the aggregated boundaries
        // list is octave2 items
        //  each map key=segmented cell label indexes,
        //           value = index to map in octave2IndexBoundsMaps
        List<Map<OneDIntArray, PairIntArray>> octave2KeyIndexMaps
            = new ArrayList<Map<OneDIntArray, PairIntArray>>();
        for (int j = 0; j < scales2.size(); ++j) {
            octave2KeyIndexMaps.add(new HashMap<OneDIntArray, PairIntArray>());
        }

        for (int i = 0; i < scales1.size(); ++i) {
        //for (int i = 0; i < 1; ++i) {

            float scale1 = scales1.get(i);

            int[][] ch1 = ch1s.get(i).a;
            //Set<PairInt> templateSet = labeledPoints1Lists.get(i);
            PairIntArray bounds1 = labeledBoundaries1.get(i);
            float sz1 = calculateObjectSize(bounds1);

            List<PObject> results = new ArrayList<PObject>();
            TDoubleList chordDiffSums = new TDoubleArrayList();
            TFloatList intersections = new TFloatArrayList();

            for (int j = 0; j < scales2.size(); ++j) {
            //for (int j = 0; j < 1; ++j) {

                float scale2 = scales2.get(j);

                List<TwoDIntArray> listOfCH2s = ch2Lists.get(j);
                if (listOfCH2s == null) {
                    continue;
                }
                List<Set<PairInt>> listOfSets2 = labeledPoints2Lists.get(j);

                Map<OneDIntArray, PairIntArray> keyBoundsMap
                    = octave2KeyIndexMaps.get(j);

                ShapeFinder shapeFinder = new ShapeFinder(
                    bounds1, ch1, scale1, sz1,
                    orb1.getPyramidImages().get(i).a[0].length -1,
                    orb1.getPyramidImages().get(i).a.length - 1,
                    listOfSets2, listOfCH2s, scale2,
                    keyBoundsMap,
                    orb2.getPyramidImages().get(j).a[0].length -1,
                    orb2.getPyramidImages().get(j).a.length - 1,
                    intersectionLimit
                );
    shapeFinder.pyr1 = orb1.getPyramidImages().get(i);
    shapeFinder.pyr2 = orb2.getPyramidImages().get(j);
    shapeFinder.lbl = Integer.toString(i) + ":" + Integer.toString(j) + "_";
    shapeFinder.oct1 = i;
    shapeFinder.oct2 = j;
    {
    //if (i==2&&j==0) {
    Image img1 = ORB.convertToImage(orb1.getPyramidImages().get(i));
    Image img2 = ORB.convertToImage(orb2.getPyramidImages().get(j));
    MiscDebug.writeImage(img2, "AAA_2_" + j);
    MiscDebug.writeImage(img1, "AAA_1_" + i);

    for (int i2 = 0; i2 < listOfSets2.size(); ++i2) {
    int clr = ImageIOHelper.getNextColorRGB(i2);
    Set<PairInt> set = listOfSets2.get(i2);
    for (PairInt p : set) {
        ImageIOHelper.addPointToImage(p.getX(), p.getY(),
        img2, 1, clr);
    }
    }
    MiscDebug.writeImage(img2,
    "_AAA_2_s_" + j);
    //}
    }
                ShapeFinderResult r = shapeFinder.findAggregated();

                if (r == null) {
                    continue;
                }

                double c = r.getChordDiffSum();

                results.add(new PObject(r, r.bounds1, r.bounds2, scale1, scale2));
                chordDiffSums.add(r.getChordDiffSum());
                intersections.add(r.intersection);

                if (r.getChordDiffSum() > maxDiffChordSum) {
                    maxDiffChordSum = r.getChordDiffSum();
                }
                double avgCD = r.getChordDiffSum() / (double) r.getNumberOfMatches();
                if (avgCD > maxAvgDiffChord) {
                    maxAvgDiffChord = avgCD;
                }
                double avgDist = r.getDistSum() / (double) r.getNumberOfMatches();
                if (avgDist > maxAvgDist) {
                    maxAvgDist = avgDist;
                }

                System.out.println(String.format(
                    "%d %d p in set=(%d,%d)  shape matcher c=%.2f np=%d  inter=%.2f  dist=%.2f avgDist=%.2f",
                    i, j, r.bounds2.getX(0), r.bounds2.getY(0),
                    (float) c, r.getNumberOfMatches(), (float) r.intersection,
                    (float) r.getDistSum(), (float) avgDist));

            } // end loop over j datasets 2
            if (!results.isEmpty()) {
                resultsMap.put(i, results);
                chordDiffSumsMap.put(i, chordDiffSums);
                intersectionsMap.put(i, intersections);
            }
        } // end loop over i dataset 1
        // calculate the Salukwdze distances

        /*
        for each i, need the max chord diff sum, nPoints in bound1, and best Results
         */

        double minSD = Double.MAX_VALUE;
        int minSDI = -1;

        TIntObjectIterator<List<PObject>> iter = resultsMap.iterator();

        for (int i = 0; i < resultsMap.size(); ++i) {

            iter.advance();

            int idx = iter.key();

            //double maxDiffChordSum = chordDiffSumsMap.get(idx).max();
            double minCost = Double.MAX_VALUE;
            int minCostIdx = -1;
            List<PObject> resultsList = resultsMap.get(idx);

            for (int j = 0; j < resultsList.size(); ++j) {

                PObject obj = resultsList.get(j);

                float costIntersection = 1.0F - intersectionsMap.get(idx).get(j);

                PartialShapeMatcher.Result r = obj.r;

                int nb1 = Math.round((float) obj.bounds1.getN() / (float) dp);

                float np = r.getNumberOfMatches();

                float countComp = 1.0F - (np / (float) nb1);
                float countCompSq = countComp * countComp;
                double chordComp = ((float) r.getChordDiffSum() / np) / maxAvgDiffChord;
                double chordCompSq = chordComp * chordComp;
                double avgDist = r.getDistSum() / np;
                double distComp = avgDist / maxAvgDist;
                double distCompSq = distComp * distComp;

                // Salukwzde uses square sums
                //double sd = r.calculateSalukwdzeDistanceSquared(
                //    maxDiffChordSum, nb1);

                // TODO: consider formal analysis of dependencies and hence
                //       error terms:
                //double sd = chordCompSq*countCompSq
                //    + distCompSq*countCompSq;

                //NOTE: The coverage of the matches is currently
                // approximated as simply numberMatched/maxNumberMatchable,
                // but a term representing the spatial distribution appears
                // to be necessary also.
                // will try largestNumberGap/maxNumberMatchable.
                // TODO: need to improve this in detail later
                int lGap = maxNumberOfGaps(obj.bounds1, r)/dp;
                float gCountComp = (float)lGap/(float)nb1;

                //double sd = chordCompSq + countCompSq + distCompSq;
                double sd = chordComp + countComp + gCountComp + distComp
                    + costIntersection;

                if (sd < minCost) {
                    minCost = sd;
                    minCostIdx = j;
                }
                if (sd < minSD) {
                    minSD = sd;
                    minSDI = idx;
                }
                System.out.println("sd=" + sd + " n1="
                    + obj.bounds1.getN() + " n2=" + obj.bounds2.getN()
                    + " origN1=" + r.getOriginalN1()
                    + " nMatches=" + r.getNumberOfMatches()
                    + String.format(
                    " chord=%.2f count=%.2f spatial=%.2f dist=%.2f inter=%.2f",
                    (float)chordComp, (float)countComp,
                    (float)gCountComp, (float)distComp,
                    (float)costIntersection)
                );
            }
            assert (minCostIdx > -1);

            TDoubleList cList = chordDiffSumsMap.get(idx);

            TFloatList iList = intersectionsMap.get(idx);

            for (int j = resultsList.size() - 1; j > -1; --j) {
                if (j != minCostIdx) {
                    resultsList.remove(j);
                    cList.removeAt(j);
                    iList.removeAt(j);
                }
            }
        }
        if (resultsMap.size() > 1) {
            // TODO: build a test for this.
            //    possibly need to transform results to same reference
            //    frame to compare.
            //    using best SD for now
            TIntSet rm = new TIntHashSet();
            iter = resultsMap.iterator();
            for (int i = 0; i < resultsMap.size(); ++i) {
                iter.advance();
                int idx = iter.key();
                if (idx != minSDI) {
                    rm.add(idx);
                }
            }
            TIntIterator iter2 = rm.iterator();
            while (iter2.hasNext()) {
                int idx = iter2.next();
                resultsMap.remove(idx);
            }
        }

        List<CorrespondenceList> topResults = new ArrayList<CorrespondenceList>();

        iter = resultsMap.iterator();

        for (int i = 0; i < resultsMap.size(); ++i) {
            iter.advance();
            int idx = iter.key();
            List<PObject> resultsList = resultsMap.get(idx);
            assert (resultsList.size() == 1);
            PObject obj = resultsList.get(0);
            int n = obj.r.getNumberOfMatches();
            PairInt[] m1 = new PairInt[n];
            PairInt[] m2 = new PairInt[n];
            float scale1 = obj.scale1;
            float scale2 = obj.scale2;
            for (int ii = 0; ii < n; ++ii) {
                int idx1 = obj.r.getIdx1(ii);
                int idx2 = obj.r.getIdx2(ii);
                int x1 = Math.round(obj.bounds1.getX(idx1) * scale1);
                int y1 = Math.round(obj.bounds1.getY(idx1) * scale1);
                int x2 = Math.round(obj.bounds2.getX(idx2) * scale2);
                int y2 = Math.round(obj.bounds2.getY(idx2) * scale2);
                m1[ii] = new PairInt(x1, y1);
                m2[ii] = new PairInt(x2, y2);
            }
            CorrespondenceList cor = new CorrespondenceList(obj.r.getTransformationParameters(), m1, m2);
            topResults.add(cor);
        }
        return topResults;
    }

    private TFloatList extractScales(List<TFloatList> scalesList) {
        TFloatList scales = new TFloatArrayList();
        for (int i = 0; i < scalesList.size(); ++i) {
            scales.add(scalesList.get(i).get(0));
        }
        return scales;
    }

    private int calculateNMaxMatchable(List<TIntList> keypointsX1, List<TIntList> keypointsX2) {
        int nMaxM = Integer.MIN_VALUE;
        for (int i = 0; i < keypointsX1.size(); ++i) {
            int n1 = keypointsX1.get(i).size();
            for (int j = 0; j < keypointsX2.size(); ++j) {
                int n2 = keypointsX2.get(j).size();
                int min = Math.min(n1, n2);
                if (min > nMaxM) {
                    nMaxM = min;
                }
            }
        }
        return nMaxM;
    }

    private int maxSize(List<TIntList> a) {
        int maxSz = Integer.MIN_VALUE;
        for (TIntList b : a) {
            int sz = b.size();
            if (sz > maxSz) {
                maxSz = sz;
            }
        }
        return maxSz;
    }

    public static double distance(int x, int y, PairInt b) {
        int diffX = x - b.getX();
        int diffY = y - b.getY();
        double dist = Math.sqrt(diffX * diffX + diffY * diffY);
        return dist;
    }

    public static int distance(PairInt p1, PairInt p2) {
        int diffX = p1.getX() - p2.getX();
        int diffY = p1.getY() - p2.getY();
        return (int) Math.sqrt(diffX * diffX + diffY * diffY);
    }

    /**
     * greedy matching of d1 to d2 by min cost, with unique mappings for
     * all indexes.
     *
     * @param d1
     * @param d2
     * @return matches - two dimensional int array of indexes in d1 and
     * d2 which are matched.
     */
    public static int[][] matchDescriptors(VeryLongBitString[] d1, VeryLongBitString[] d2, List<PairInt> keypoints1, List<PairInt> keypoints2) {
        int n1 = d1.length;
        int n2 = d2.length;
        //[n1][n2]
        int[][] cost = ORB.calcDescriptorCostMatrix(d1, d2);
        int[][] matches = greedyMatch(keypoints1, keypoints2, cost);
        // greedy or optimal match can be performed here.
        // NOTE: some matching problems might benefit from using the spatial
        //   information at the same time.  for those, will consider adding
        //   an evaluation term for these descriptors to a specialization of
        //   PartialShapeMatcher.java
        return matches;
    }

    /**
     * greedy matching of d1 to d2 by min difference, with unique mappings for
     * all indexes.
     * NOTE that if 2 descriptors match equally well, either one
     * might get the assignment.
     * Consider using instead, matchDescriptors2 which matches
     * by descriptor and relative spatial location.
     *
     * @param d1
     * @param d2
     * @param keypoints2
     * @param keypoints1
     * @return matches - two dimensional int array of indexes in d1 and
     * d2 which are matched.
     */
    public static int[][] matchDescriptors(ORB.Descriptors[] d1, ORB.Descriptors[] d2, List<PairInt> keypoints1, List<PairInt> keypoints2) {
        if (d1.length != d2.length) {
            throw new IllegalArgumentException("d1 and d2 must" + " be same length");
        }
        int n1 = d1[0].descriptors.length;
        int n2 = d2[0].descriptors.length;
        if (n1 != keypoints1.size()) {
            throw new IllegalArgumentException("number of descriptors in " + " d1 bitstrings must be same as keypoints1 length");
        }
        if (n2 != keypoints2.size()) {
            throw new IllegalArgumentException("number of descriptors in " + " d2 bitstrings must be same as keypoints2 length");
        }
        //[n1][n2]
        int[][] cost = ORB.calcDescriptorCostMatrix(d1, d2);
        int[][] matches = greedyMatch(keypoints1, keypoints2, cost);
        // greedy or optimal match can be performed here.
        // NOTE: some matching problems might benefit from using the spatial
        //   information at the same time.  for those, will consider adding
        //   an evaluation term for these descriptors to a specialization of
        //   PartialShapeMatcher.java
        return matches;
    }

    private static void debugPrint(List<QuadInt> pairs, int i, int j, TwoDFloatArray pyr1, TwoDFloatArray pyr2, float s1, float s2) {
        Image img1 = ORB.convertToImage(pyr1);
        Image img2 = ORB.convertToImage(pyr2);
        try {
            for (int ii = 0; ii < pairs.size(); ++ii) {
                QuadInt q = pairs.get(ii);
                int x1 = Math.round(q.getA() / s1);
                int y1 = Math.round(q.getB() / s1);
                int x2 = Math.round(q.getC() / s2);
                int y2 = Math.round(q.getD() / s2);
                ImageIOHelper.addPointToImage(x1, y1, img1, 1, 255, 0, 0);
                ImageIOHelper.addPointToImage(x2, y2, img2, 1, 255, 0, 0);
            }
            String strI = Integer.toString(i);
            while (strI.length() < 3) {
                strI = "0" + strI;
            }
            String strJ = Integer.toString(j);
            while (strJ.length() < 3) {
                strJ = "0" + strJ;
            }
            String str = "_pairs_" + strI + "_" + strJ + "_";
            MiscDebug.writeImage(img1, str + "_" + strI);
            MiscDebug.writeImage(img2, str + "_" + strJ);
        } catch (Exception e) {
        }
    }

    private static void debugPrint(TwoDFloatArray pyr1, TwoDFloatArray pyr2, TIntList kpX1, TIntList kpY1, TIntList kpX2, TIntList kpY2, float scale1, float scale2, int img1Idx, int img2Idx) {
        Image img1 = ORB.convertToImage(pyr1);
        Image img2 = ORB.convertToImage(pyr2);
        try {
            for (int i = 0; i < kpX1.size(); ++i) {
                int x1 = (int) (kpX1.get(i) / scale1);
                int y1 = (int) (kpY1.get(i) / scale1);
                ImageIOHelper.addPointToImage(x1, y1, img1, 1, 255, 0, 0);
            }
            for (int i = 0; i < kpX2.size(); ++i) {
                int x2 = (int) (kpX2.get(i) / scale2);
                int y2 = (int) (kpY2.get(i) / scale2);
                ImageIOHelper.addPointToImage(x2, y2, img2, 1, 255, 0, 0);
            }
            String strI = Integer.toString(img1Idx);
            while (strI.length() < 3) {
                strI = "0" + strI;
            }
            String strJ = Integer.toString(img2Idx);
            while (strJ.length() < 3) {
                strJ = "0" + strJ;
            }
            String str = "_kp_" + strI + "_" + strJ + "_";
            MiscDebug.writeImage(img1, str + "_i");
            MiscDebug.writeImage(img2, str + "_j");
        } catch (Exception e) {
        }
    }

    private static void debugPrint(TwoDFloatArray pyr1, TwoDFloatArray pyr2, TIntList kpX1, TIntList kpY1, TIntList kpX2, TIntList kpY2, int img1Idx, int img2Idx) {
        Image img1 = ORB.convertToImage(pyr1);
        Image img2 = ORB.convertToImage(pyr2);
        try {
            for (int i = 0; i < kpX1.size(); ++i) {
                int x1 = kpX1.get(i);
                int y1 = kpY1.get(i);
                ImageIOHelper.addPointToImage(x1, y1, img1, 1, 255, 0, 0);
            }
            for (int i = 0; i < kpX2.size(); ++i) {
                int x2 = kpX2.get(i);
                int y2 = kpY2.get(i);
                ImageIOHelper.addPointToImage(x2, y2, img2, 1, 255, 0, 0);
            }
            String strI = Integer.toString(img1Idx);
            while (strI.length() < 3) {
                strI = "0" + strI;
            }
            String strJ = Integer.toString(img2Idx);
            while (strJ.length() < 3) {
                strJ = "0" + strJ;
            }
            String str = "_kp_" + strI + "_" + strJ + "_";
            MiscDebug.writeImage(img1, str + "_i");
            MiscDebug.writeImage(img2, str + "_j");
        } catch (Exception e) {
        }
    }

    private static Set<PairInt> makeSet(TIntList kpX1, TIntList kpY1) {
        Set<PairInt> set = new HashSet<PairInt>();
        for (int i = 0; i < kpX1.size(); ++i) {
            PairInt p = new PairInt(kpX1.get(i), kpY1.get(i));
            set.add(p);
        }
        return set;
    }

    private static Descriptors[] getDescriptors(ORB orb, int i) {
        ORB.Descriptors[] d = null;
        if (orb.getDescrChoice().equals(ORB.DescriptorChoice.ALT)) {
            d = new ORB.Descriptors[]{orb.getDescriptorsListAlt().get(i)};
        } else if (orb.getDescrChoice().equals(ORB.DescriptorChoice.HSV)) {
            d = new ORB.Descriptors[]{orb.getDescriptorsH().get(i), orb.getDescriptorsS().get(i), orb.getDescriptorsV().get(i)};
        } else if (orb.getDescrChoice().equals(ORB.DescriptorChoice.GREYSCALE)) {
            d = new ORB.Descriptors[]{orb.getDescriptorsList().get(i)};
        }
        return d;
    }

    private static List<QuadInt> createPairLabelIndexes(int[][] cost, int nBands, List<TIntList> pointIndexLists1, TIntList kpX1, TIntList kpY1, List<TIntList> pointIndexLists2, TIntList kpX2, TIntList kpY2) {
        int costLimit = Math.round((float) (nBands * 256) * 0.65F);
        int minP1Diff = 3;
        Set<QuadInt> exists = new HashSet<QuadInt>();
        // pairs of idx from set1 and idx from set 2
        Set<PairInt> skip = new HashSet<PairInt>();
        List<QuadInt> pairIndexes = new ArrayList<QuadInt>();
        List<PairInt> pair2Indexes = calculatePairIndexes(pointIndexLists2, kpX2, kpY2, minP1Diff);
        for (int ii = 0; ii < pointIndexLists1.size(); ++ii) {
            TIntList kpIndexes1 = pointIndexLists1.get(ii);
            if (kpIndexes1.size() < 2) {
                continue;
            }
            // draw 2 from kpIndexes1
            for (int ii1 = 0; ii1 < kpIndexes1.size(); ++ii1) {
                int idx1 = kpIndexes1.get(ii1);
                int t1X = kpX1.get(idx1);
                int t1Y = kpY1.get(idx1);
                boolean skipIdx1 = false;
                for (int ii2 = 0; ii2 < kpIndexes1.size(); ++ii2) {
                    if (ii1 == ii2) {
                        continue;
                    }
                    int idx2 = kpIndexes1.get(ii2);
                    int t2X = kpX1.get(idx2);
                    int t2Y = kpY1.get(idx2);
                    if (t1X == t2X && t1Y == t2Y) {
                        continue;
                    }
                    int diffX = t1X - t2X;
                    int diffY = t1Y - t2Y;
                    int distSq = diffX * diffX + diffY * diffY;
                    //if (distSq > limitSq) {
                    //    continue;
                    //}
                    if (distSq < minP1Diff * minP1Diff) {
                        continue;
                    }
                    for (PairInt p2Index : pair2Indexes) {
                        int idx3 = p2Index.getX();
                        int idx4 = p2Index.getY();
                        PairInt p13 = new PairInt(idx1, idx3);
                        if (skip.contains(p13)) {
                            skipIdx1 = true;
                            break;
                        }
                        PairInt p24 = new PairInt(idx2, idx4);
                        if (skip.contains(p24)) {
                            continue;
                        }
                        int c13 = cost[idx1][idx3];
                        // if idx1 and idx3 cost is above limit, skip
                        if (c13 > costLimit) {
                            skip.add(p13);
                            skipIdx1 = true;
                            break;
                        }
                        int c24 = cost[idx2][idx4];
                        if (c24 > costLimit) {
                            skip.add(p24);
                            continue;
                        }
                        QuadInt q = new QuadInt(idx1, idx2, idx3, idx4);
                        QuadInt qChk = new QuadInt(idx2, idx1, idx4, idx3);
                        if (exists.contains(q) || exists.contains(qChk)) {
                            continue;
                        }
                        /*
                        int s1X = kpX2.get(idx3);
                        int s1Y = kpY2.get(idx3);
                        int s2X = kpX2.get(idx4);
                        int s2Y = kpY2.get(idx4);
                        int diffX2 = s1X - s2X;
                        int diffY2 = s1Y - s2Y;
                        int distSq2 = diffX2 * diffX2 + diffY2 * diffY2;
                        //if (distSq2 > limitSq) {
                        //    continue;
                        //}
                        if ((distSq2 < minP1Diff * minP1Diff)) {
                        continue;
                        }
                         */
                        pairIndexes.add(q);
                        exists.add(q);
                    }
                    if (skipIdx1) {
                        break;
                    }
                }
            }
        }
        return pairIndexes;
    }

    public static int calculateObjectSize(Set<PairInt> points) {
        // O(N*lg_2(N))
        FurthestPair furthestPair = new FurthestPair();
        PairInt[] fp = furthestPair.find(points);
        if (fp == null || fp.length < 2) {
            throw new IllegalArgumentException("did not find a furthest pair" + " in points");
        }
        double dist = ORBMatcher.distance(fp[0], fp[1]);
        return (int) Math.round(dist);
    }

    public static int calculateObjectSize(PairIntArray points) {
        return calculateObjectSize(Misc.convert(points));
    }

    private static float calculateDiagonal(List<TIntList> keypointsX1, List<TIntList> keypointsY1, int idx) {
        TIntList x1 = keypointsX1.get(idx);
        TIntList y1 = keypointsY1.get(idx);
        int maxX = x1.max();
        int maxY = y1.max();
        return (float) Math.sqrt(maxX * maxX + maxY * maxY);
    }

    /**
     * NOTE: preliminary results show that this matches the right pattern as
     * a subset of the object, but needs to be followed by a slightly larger
     * aggregated search by segmentation cells using partial shape matcher
     * for example.  This was started in ShapeFinder, but needs to be
     * adjusted for a search given seed cells and possibly improved for the
     * other TODO items).
     * @param keypoints1
     * @param keypoints2
     * @param mT
     * @param mS
     * @param nn
     * @param minMaxXY2
     * @param limit
     * @param tIndexes
     * @param idx1P2CostMap
     * @param indexes
     * @param costs
     * @return
     */
    private static List<CorrespondenceList> completeUsingCombinations(List<PairInt> keypoints1, List<PairInt> keypoints2, PairIntArray mT, PairIntArray mS, NearestNeighbor2D nn, int[] minMaxXY2, int limit, TIntList tIndexes, TIntObjectMap<TObjectIntMap<PairInt>> idx1P2CostMap, PairInt[] indexes, int[] costs, int bitTolerance, int nBands) {
        int nTop = mT.getN();
        System.out.println("have " + nTop + " sets of points for " + " n of k=2 combinations");
        // need to make pairs of combinations from mT,mS
        //  to calcuate euclidean transformations and evaluate them.
        // -- can reduce the number of combinations by imposing a
        //    distance limit on separation of feasible pairs
        int limitSq = limit * limit;
        MatchedPointsTransformationCalculator tc = new MatchedPointsTransformationCalculator();
        Transformer transformer = new Transformer();
        /* can try to use the tolerance in 2 different ways:
        (1) to adjust the weight of cost.
        weight = 1 - (tolerance/(256*nBands)).
        This approach might favor regions of many false keypoints
        such as highly textured regions.
        (2) keep the minCost solution, but also keep all solutions
        that are within cost + tolerance of it.
        This doesn't discard the true solution.
        ---> these will need additional information to distinguish
        between solutions.
        procedding with (2).
         */
        // this fixed size sorted vector is faster for shorter arrays.
        // TODO: consider ways to robustly set the size from the cost
        // statistics to ensure the vector will always contain the
        // correct solution even if not in top position.
        int nt = mT.getN();
        FixedSizeSortedVector<CObject> vec = new FixedSizeSortedVector<CObject>(nt, CObject.class);
        double minCost = Double.MAX_VALUE;
        //CorrespondenceList minCostCor = null;
        //PairIntArray minCostTrT = null;
        double[] minCostI = new double[nTop];
        double[] minDistI = new double[nTop];
        // temporary storage of corresp coords until object construction
        int[] m1x = new int[nTop];
        int[] m1y = new int[nTop];
        int[] m2x = new int[nTop];
        int[] m2y = new int[nTop];
        int mCount = 0;
        for (int i = 0; i < nTop; ++i) {
            int t1X = mT.getX(i);
            int t1Y = mT.getY(i);
            int s1X = mS.getX(i);
            int s1Y = mS.getY(i);
            // choose all combinations of 2nd point within distance
            // limit of point s1.
            for (int j = i + 1; j < mS.getN(); ++j) {
                int t2X = mT.getX(j);
                int t2Y = mT.getY(j);
                int s2X = mS.getX(j);
                int s2Y = mS.getY(j);
                if ((t1X == t2X && t1Y == t2Y) || (s1X == s2X && s1Y == s2Y)) {
                    continue;
                }
                int diffX = s1X - s2X;
                int diffY = s1Y - s2Y;
                int distSq = diffX * diffX + diffY * diffY;
                if (distSq > limitSq) {
                    continue;
                }
                // -- calculate euclid transformation
                // -- evaluate the fit
                TransformationParameters params = tc.calulateEuclidean(t1X, t1Y, t2X, t2Y, s1X, s1Y, s2X, s2Y, 0, 0);
                float scale = params.getScale();
                mCount = 0;
                // template object transformed
                PairIntArray trT = transformer.applyTransformation(params, mT);
                /*
                two components to the evaluation and both need normalizations
                so that their contributions to total result are
                equally weighted.
                (1) descriptors:
                -- score is sum of each matched (3*256 - cost)
                -- the normalization is the maximum possible score,
                so will use the number of template points.
                --> norm = nTemplate * 3 * 256
                -- normalized score = (3*256 - cost)/norm
                ==> normalized cost = 1 - ((3*256 - cost)/norm)
                (2) spatial distances from transformed points:
                -- sum of distances within limit
                and replacement of distance by limit if no matching
                nearest neighbor is found.
                -- divide each distance by the transformation scale
                to compare same values
                -- divide the total sum by the total max possible
                --> norm = nTemplate * limit / scale
                Then the total cost is (1) + (2) and the min cost
                among all of these combinations is the resulting
                correspondence list
                 */
                double maxCost = nBands * 256;
                double maxDist = limit / scale;
                double sum1 = 0;
                double sum2 = 0;
                double sum = 0;
                for (int k = 0; k < trT.getN(); ++k) {
                    int xTr = trT.getX(k);
                    int yTr = trT.getY(k);
                    int idx1 = tIndexes.get(k);
                    Set<PairInt> nearest = null;
                    if ((xTr >= 0) && (yTr >= 0) && (xTr <= (minMaxXY2[1] + limit)) && (yTr <= (minMaxXY2[3] + limit))) {
                        nearest = nn.findClosest(xTr, yTr, limit);
                    }
                    int minC = Integer.MAX_VALUE;
                    PairInt minCP2 = null;
                    if (nearest != null && !nearest.isEmpty()) {
                        TObjectIntMap<PairInt> cMap = idx1P2CostMap.get(idx1);
                        for (PairInt p2 : nearest) {
                            if (!cMap.containsKey(p2)) {
                                continue;
                            }
                            int c = cMap.get(p2);
                            if (c < minC) {
                                minC = c;
                                minCP2 = p2;
                            }
                        }
                    }
                    if (minCP2 != null) {
                        double scoreNorm = (nBands * 256 - minC) / maxCost;
                        double costNorm = 1.0 - scoreNorm;
                        sum1 += costNorm;
                        double dist = ORBMatcher.distance(xTr, yTr, minCP2);
                        double distNorm = dist / maxDist;
                        sum2 += distNorm;
                        m1x[mCount] = keypoints1.get(idx1).getX();
                        m1y[mCount] = keypoints1.get(idx1).getY();
                        m2x[mCount] = minCP2.getX();
                        m2y[mCount] = minCP2.getY();
                        minCostI[mCount] = costNorm;
                        minDistI[mCount] = distNorm;
                        mCount++;
                    } else {
                        sum1 += 1;
                        sum2 += 1;
                    }
                }
                sum = sum1 + sum2;
                if ((minCost == Double.MAX_VALUE) || (sum < (minCost + bitTolerance))) {
                    if (sum < minCost) {
                        minCost = sum;
                    }
                    List<PairInt> m1 = new ArrayList<PairInt>();
                    List<PairInt> m2 = new ArrayList<PairInt>();
                    CorrespondenceList corr = new CorrespondenceList(params.getScale(), Math.round(params.getRotationInDegrees()), Math.round(params.getTranslationX()), Math.round(params.getTranslationY()), 0, 0, 0, m1, m2);
                    for (int mi = 0; mi < mCount; ++mi) {
                        m1.add(new PairInt(m1x[mi], m1y[mi]));
                        m2.add(new PairInt(m2x[mi], m2y[mi]));
                    }
                    CObject cObj = new CObject(sum, corr, trT);
                    vec.add(cObj);
                }
            }
        }
        if (vec.getNumberOfItems() == 0) {
            return null;
        }
        List<CorrespondenceList> topResults = new ArrayList<CorrespondenceList>();
        for (int i = 0; i < vec.getNumberOfItems(); ++i) {
            CObject a = vec.getArray()[i];
            if (a.cost > (minCost + bitTolerance)) {
                break;
            }
            topResults.add(a.cCor);
        }
        return topResults;
    }

    private static List<PairInt> calculatePairIndexes(List<TIntList> pointIndexLists2, TIntList kpX2, TIntList kpY2, int minPDiff) {
        List<PairInt> pairIndexes = new ArrayList<PairInt>();
        Set<PairInt> exists = new HashSet<PairInt>();
        // draw 2 pairs from other dataset
        for (int jj = 0; jj < pointIndexLists2.size(); ++jj) {
            TIntList kpIndexes2 = pointIndexLists2.get(jj);
            if (kpIndexes2.size() < 2) {
                continue;
            }
            // draw 2 from kpIndexes2
            for (int jj1 = 0; jj1 < kpIndexes2.size(); ++jj1) {
                int idx3 = kpIndexes2.get(jj1);
                int s1X = kpX2.get(idx3);
                int s1Y = kpY2.get(idx3);
                for (int jj2 = 0; jj2 < kpIndexes2.size(); ++jj2) {
                    if (jj1 == jj2) {
                        continue;
                    }
                    int idx4 = kpIndexes2.get(jj2);
                    int s2X = kpX2.get(idx4);
                    int s2Y = kpY2.get(idx4);
                    if (s1X == s2X && s1Y == s2Y) {
                        continue;
                    }
                    PairInt q = new PairInt(idx3, idx4);
                    if (exists.contains(q)) {
                        continue;
                    }
                    int diffX2 = s1X - s2X;
                    int diffY2 = s1Y - s2Y;
                    int distSq2 = diffX2 * diffX2 + diffY2 * diffY2;
                    //if (distSq2 > limitSq) {
                    //    continue;
                    //}
                    if (distSq2 < minPDiff * minPDiff) {
                        continue;
                    }
                    pairIndexes.add(q);
                    exists.add(q);
                }
            }
        }
        return pairIndexes;
    }

    private static float calculateDiagonal2(List<TwoDFloatArray> pyramidImages, int idx) {
        int w = pyramidImages.get(idx).a.length;
        int h = pyramidImages.get(idx).a[0].length;
        double diag = Math.sqrt(w * w + h * h);
        return (float) diag;
    }

    private static int[][] greedyMatch(List<PairInt> keypoints1,
        List<PairInt> keypoints2, int[][] cost) {
        int n1 = keypoints1.size();
        int n2 = keypoints2.size();
        // for the greedy match, separating the index information from the cost
        // and then sorting by cost
        int nTot = n1 * n2;
        PairInt[] indexes = new PairInt[nTot];
        int[] costs = new int[nTot];
        int count = 0;
        for (int i = 0; i < n1; ++i) {
            for (int j = 0; j < n2; ++j) {
                indexes[count] = new PairInt(i, j);
                costs[count] = cost[i][j];
                count++;
            }
        }
        assert (count == nTot);
        QuickSort.sortBy1stArg(costs, indexes);
        Set<PairInt> set1 = new HashSet<PairInt>();
        Set<PairInt> set2 = new HashSet<PairInt>();
        List<PairInt> matches = new ArrayList<PairInt>();
        // visit lowest costs (== differences) first
        for (int i = 0; i < nTot; ++i) {
            PairInt index12 = indexes[i];
            int idx1 = index12.getX();
            int idx2 = index12.getY();
            PairInt p1 = keypoints1.get(idx1);
            PairInt p2 = keypoints2.get(idx2);
            if (set1.contains(p1) || set2.contains(p2)) {
                continue;
            }
            //System.out.println("p1=" + p1 + " " + " p2=" + p2 + " cost=" + costs[i]);
            matches.add(index12);
            set1.add(p1);
            set2.add(p2);
        }
        int[][] results = new int[matches.size()][2];
        for (int i = 0; i < matches.size(); ++i) {
            results[i][0] = matches.get(i).getX();
            results[i][1] = matches.get(i).getY();
        }
        return results;
    }

    private static double[] sumKeypointDistanceDifference(TIntList a2Indexes,
        PairIntArray tr2, TIntList kpX2, TIntList kpY2, NearestNeighbor2D nn,
        TransformationParameters params, int maxX, int maxY, int pixTolerance,
        double maxDist, int[] m1x, int[] m1y, int[] m2x, int[] m2y) {
        double sum2 = 0;
        int mCount = 0;
        for (int k = 0; k < tr2.getN(); ++k) {
            int x2Tr = tr2.getX(k);
            int y2Tr = tr2.getY(k);
            int idx2 = a2Indexes.get(k);
            Set<PairInt> nearest = null;
            if ((x2Tr >= 0) && (y2Tr >= 0) && (x2Tr <= (maxX + pixTolerance)) && (y2Tr <= (maxY + pixTolerance))) {
                nearest = nn.findClosest(x2Tr, y2Tr, pixTolerance);
            }
            double minDist = Double.MAX_VALUE;
            PairInt minDistP1 = null;
            if (nearest != null && !nearest.isEmpty()) {
                for (PairInt p11 : nearest) {
                    double dist = ORBMatcher.distance(x2Tr, y2Tr, p11);
                    if (dist < minDist) {
                        minDist = dist;
                        minDistP1 = p11;
                    }
                }
            }
            if (minDistP1 != null) {
                double dist = minDist;
                double distNorm = dist / maxDist;
                sum2 += distNorm;
                m2x[mCount] = kpX2.get(idx2);
                m2y[mCount] = kpY2.get(idx2);
                m1x[mCount] = minDistP1.getX();
                m1y[mCount] = minDistP1.getY();
                mCount++;
            } else {
                sum2 += 1;
            }
        } // end loop over trnsformed set 2
        return new double[]{sum2, mCount};
    }

    private static double[] sumKeypointDescAndDist(int[][] cost, int nBands,
        TIntList a1Indexes, PairIntArray tr1, TIntList kpX1, TIntList kpY1,
        NearestNeighbor2D nn2, TObjectIntMap<PairInt> p2KPIndexMap,
        TransformationParameters params, int maxX2, int maxY2, int pixTolerance,
        double maxDist, PairInt[] m1, PairInt[] m2) {
        double sumDesc = 0;
        double sumDist = 0;
        int count = 0;
        double maxDesc = nBands * 256.0;
        for (int k = 0; k < tr1.getN(); ++k) {
            int x1Tr = tr1.getX(k);
            int y1Tr = tr1.getY(k);
            int idx1 = a1Indexes.get(k);
            Set<PairInt> nearest = null;
            if ((x1Tr >= 0) && (y1Tr >= 0) && (x1Tr <= (maxX2 + pixTolerance)) && (y1Tr <= (maxY2 + pixTolerance))) {
                nearest = nn2.findClosest(x1Tr, y1Tr, pixTolerance);
            }
            int minC = Integer.MAX_VALUE;
            PairInt minCP2 = null;
            int minIdx2 = 0;
            if (nearest != null && !nearest.isEmpty()) {
                for (PairInt p2 : nearest) {
                    int idx2 = p2KPIndexMap.get(p2);
                    int c = cost[idx1][idx2];
                    if (c < minC) {
                        minC = c;
                        minCP2 = p2;
                        minIdx2 = idx2;
                    }
                }
            }
            if (minCP2 != null) {
                double scoreNorm = (nBands * 256 - minC) / maxDesc;
                double costNorm = 1.0 - scoreNorm;
                sumDesc += costNorm;
                double dist = ORBMatcher.distance(x1Tr, y1Tr, minCP2);
                double distNorm = dist / maxDist;
                sumDist += distNorm;
                m1[count] = new PairInt(kpX1.get(idx1), kpY1.get(idx1));
                m2[count] = minCP2;
                count++;
            } else {
                sumDesc += 1;
                sumDist += 1;
            }
        }
        return new double[]{sumDesc, sumDist, count};
    }

    private static double[] sumKeypointDescAndDist(int[][] cost, int nBands,
        TIntList a1Indexes, PairIntArray tr1, TIntList kpX1, TIntList kpY1,
        NearestNeighbor2D nn2, TObjectIntMap<PairInt> p2KPIndexMap,
        int maxX2, int maxY2, int pixTolerance, double maxDist, int[] m1, int[] m2) {
        double sumDesc = 0;
        double sumDist = 0;
        int count = 0;
        double maxDesc = nBands * 256.0;
        //best first match, after nearest neighbors
        // TODO: consider optimal bipartite matching when have an
        //       implementation of multi-level-buckets
        float[] costA = new float[tr1.getN()];
        float[] costDesc = new float[tr1.getN()];
        float[] costDist = new float[tr1.getN()];
        int[] indexes = new int[tr1.getN()];
        for (int k = 0; k < tr1.getN(); ++k) {
            int x1Tr = tr1.getX(k);
            int y1Tr = tr1.getY(k);
            int idx1 = a1Indexes.get(k);
            Set<PairInt> nearest = null;
            if ((x1Tr >= 0) && (y1Tr >= 0) && (x1Tr <= (maxX2 + pixTolerance)) && (y1Tr <= (maxY2 + pixTolerance))) {
                nearest = nn2.findClosest(x1Tr, y1Tr, pixTolerance);
            }
            int minC = Integer.MAX_VALUE;
            PairInt minCP2 = null;
            int minIdx2 = 0;
            if (nearest != null && !nearest.isEmpty()) {
                for (PairInt p2 : nearest) {
                    int idx2 = p2KPIndexMap.get(p2);
                    int c = cost[idx1][idx2];
                    if (c < minC) {
                        minC = c;
                        minCP2 = p2;
                        minIdx2 = idx2;
                    }
                }
            }
            if (minCP2 != null) {
                double scoreNorm = (nBands * 256 - minC) / maxDesc;
                double costNorm = 1.0 - scoreNorm;
                sumDesc += costNorm;
                double dist = ORBMatcher.distance(x1Tr, y1Tr, minCP2);
                double distNorm = dist / maxDist;
                sumDist += distNorm;
                m1[count] = idx1;
                m2[count] = minIdx2;
                costA[count] = (float) (costNorm + distNorm);
                costDesc[count] = (float) costNorm;
                costDist[count] = (float) distNorm;
                indexes[count] = count;
                count++;
            } else {
                sumDesc += 1;
                sumDist += 1;
            }
        }
        if (count > 1) {
            costA = Arrays.copyOf(costA, count);
            indexes = Arrays.copyOf(indexes, count);
            QuickSort.sortBy1stArg(costA, indexes);
            TIntSet set1 = new TIntHashSet();
            TIntSet set2 = new TIntHashSet();
            List<PairInt> matched = new ArrayList<PairInt>();
            TIntList idxs = new TIntArrayList();
            for (int i = 0; i < count; ++i) {
                int idx = indexes[i];
                int idx1 = m1[idx];
                int idx2 = m2[idx];
                if (set1.contains(idx1) || set2.contains(idx2)) {
                    continue;
                }
                idxs.add(idx);
                matched.add(new PairInt(idx1, idx2));
                set1.add(idx1);
                set2.add(idx2);
            }
            int nRedundant = count - matched.size();
            if (nRedundant > 0) {
                sumDesc = 0;
                sumDist = 0;
                for (int i = 0; i < matched.size(); ++i) {
                    m1[i] = matched.get(i).getX();
                    m2[i] = matched.get(i).getY();
                    int idx = idxs.get(i);
                    sumDesc += costDesc[idx];
                    sumDist += costDist[idx];
                }
                sumDesc += (tr1.getN() - matched.size());
                sumDist += (tr1.getN() - matched.size());
                count = matched.size();
            }
        }
        return new double[]{sumDesc, sumDist, count};
    }

    /**
     * match left and right using the right transformed to left and nn1
     * and return the sum of the descriptor costs, distance differences
     * and number matched.
     * NOTE that the each item in a sum has been normalized before adding,
     * so for example, one item contributes between 0 and 1 to the total sum.
     * a "+1" or +maxValue was not added for points not matched.
     * @param cost
     * @param nBands
     * @param a2
     * @param a2TrTo1
     * @param nn1
     * @param p1KPIndexMap point indexes lookup w.r.t. entire keypoints1 list
     *    which has same indexes as costD first dimension.  note that the
     *    template object dataset, a.k.a. left or a1 includes all points of
     *    keypoints1...no segmentation dividing the points into more than one
     *    array list. nn1 contains all keypoints1 points.
     * @param p2KPIndexMap point indexes lookup w.r.t. entire keypoints2 list
     *     which has same indexes as costD second dimension.  these indexes are not the same
     *     indexes as a2 indexes.
     * @param img1Width
     * @param img1Height
     * @param pixTolerance
     * @param maxDist
     * @param m1
     * @param m2
     * @return
     */
    private static double[] sumKeypointDescAndDist2To1(
        int[][] cost, int nBands,
        PairIntArray a2,
        PairIntArray a2TrTo1,
        NearestNeighbor2D nn1,
        TObjectIntMap<PairInt> p1KPIndexMap,
        TObjectIntMap<PairInt> p2KPIndexMap,
        int distTol, int[] m1, int[] m2) {

        int n2 = a2TrTo1.getN();

        float distMax = (float)Math.sqrt(2) * distTol;
        
        int count = 0;
        double maxDesc = nBands * 256.0;
        //best first match, after nearest neighbors
        // TODO: consider optimal bipartite matching when have an
        //       implementation of multi-level-buckets
        float[] costA = new float[n2];
        float[] costDesc = new float[n2];
        float[] costDist = new float[n2];
        int[] indexes = new int[n2];
        for (int k = 0; k < n2; ++k) {
            int x1Tr = a2TrTo1.getX(k);
            int y1Tr = a2TrTo1.getY(k);

            if (x1Tr < 0 || y1Tr < 0) {
                continue;
            }

            int kpIdx2 = p2KPIndexMap.get(new PairInt(a2.getX(k), a2.getY(k)));

            Set<PairInt> nearest = nn1.findClosest(x1Tr, y1Tr, distTol);

            int minC = Integer.MAX_VALUE;
            PairInt minCP1 = null;
            int minCIdx1 = 0;
            if (nearest != null && !nearest.isEmpty()) {
                for (PairInt p1 : nearest) {
                    int kpIdx1 = p1KPIndexMap.get(p1);
                    int c = cost[kpIdx1][kpIdx2];
                    if (c < minC) {
                        minC = c;
                        minCP1 = p1;
                        minCIdx1 = kpIdx1;
                    }
                }
            }
            if (minCP1 != null) {
                double scoreNorm = (nBands * 256 - minC) / maxDesc;
                double costNorm = 1.0 - scoreNorm;
                double dist = ORBMatcher.distance(x1Tr, y1Tr, minCP1);
                assert(dist <= distMax);
                double distNorm = dist / distMax;
                m1[count] = minCIdx1;
                m2[count] = k;// index of a2 array
                costA[count] = (float) (costNorm + distNorm);
                costDesc[count] = (float) costNorm;
                costDist[count] = (float) distNorm;
                indexes[count] = count;
                count++;
            }
        }
        double sumDesc = 0;
        double sumDist = 0;
        if (count > 1) {
            if (count < costA.length) {
                costA = Arrays.copyOf(costA, count);
                indexes = Arrays.copyOf(indexes, count);
                costDesc = Arrays.copyOf(costDesc, count);
                costDist = Arrays.copyOf(costDist, count);
            }
            QuickSort.sortBy1stArg(costA, indexes);
            TIntSet set1 = new TIntHashSet();
            TIntSet set2 = new TIntHashSet();
            PairIntArray matched = new PairIntArray(count);
            TIntList idxs = new TIntArrayList();
            for (int i = 0; i < count; ++i) {
                int idx = indexes[i];
                int idx1 = m1[idx];
                int idx2 = m2[idx];
                if (set1.contains(idx1) || set2.contains(idx2)) {
                    continue;
                }
                idxs.add(idx);
                matched.add(idx1, idx2);
                set1.add(idx1);
                set2.add(idx2);
                sumDesc += costDesc[idx];
                sumDist += costDist[idx];
            }
            count = matched.getN();
            for (int i = 0; i < count; ++i) {
                m1[i] = matched.getX(i);
                m2[i] = matched.getY(i);
            }
        }

        return new double[]{sumDesc, sumDist, count};
    }

    /**
     * match left and right using the right transformed to left and nn1
     * and return the sum of the descriptor costs, distance differences
     * and number matched.
     * NOTE that the each item in a sum has been normalized before adding,
     * so for example, one item contributes between 0 and 1 to the total sum.
     * a "+1" or +maxValue was not added for points not matched.
     * @param a2
     * @param a2TrTo1
     * @param nn1
     * @param p1KPIndexMap
     * @param p2KPIndexMap
     * @param img1Width
     * @param img1Height
     * @param pixTolerance
     * @param maxDist
     * @param m1
     * @param m2
     * @return
     */
    private static double[] sumKeypointDescAndDist2To1(
        PairIntArray a2,
        PairIntArray a2TrTo1,
        NearestNeighbor2D nn1,
        TObjectIntMap<PairInt> p1KPIndexMap,
        TObjectIntMap<PairInt> p2KPIndexMap,
        int img1Width, int img1Height, int distTol, int[] m1, int[] m2) {

        int n2 = a2TrTo1.getN();
        
        float distMax = (float)Math.sqrt(2) * distTol;

        int count = 0;

        float[] costDist = new float[n2];
        int[] indexes = new int[n2];
        for (int k = 0; k < n2; ++k) {
            int x1Tr = a2TrTo1.getX(k);
            int y1Tr = a2TrTo1.getY(k);

            if (x1Tr < 0 || y1Tr < 0) {
                continue;
            }

            int kpIdx2 = p2KPIndexMap.get(new PairInt(a2.getX(k), a2.getY(k)));

            Set<PairInt> nearest = nn1.findClosest(x1Tr, y1Tr, distTol);

            PairInt minP1 = null;
            int minIdx1 = 0;
            if (nearest != null && !nearest.isEmpty()) {
                minP1 = nearest.iterator().next();
                minIdx1 = p1KPIndexMap.get(minP1);
            }
            if (minP1 != null) {
                double dist = ORBMatcher.distance(x1Tr, y1Tr, minP1);
                assert(dist <= distMax);
                double distNorm = dist / distMax;
                m1[count] = minIdx1;
                m2[count] = kpIdx2;
                costDist[count] = (float) distNorm;
                indexes[count] = count;
                count++;
            }
        }
        double sumDist = 0;
        if (count > 1) {
            if (count < costDist.length) {
                indexes = Arrays.copyOf(indexes, count);
                costDist = Arrays.copyOf(costDist, count);
            }
            QuickSort.sortBy1stArg(costDist, indexes);
            TIntSet set1 = new TIntHashSet();
            TIntSet set2 = new TIntHashSet();
            PairIntArray matched = new PairIntArray(count);
            TIntList idxs = new TIntArrayList();
            for (int i = 0; i < count; ++i) {
                int idx = indexes[i];
                int idx1 = m1[idx];
                int idx2 = m2[idx];
                if (set1.contains(idx1) || set2.contains(idx2)) {
                    continue;
                }
                idxs.add(idx);
                matched.add(idx1, idx2);
                set1.add(idx1);
                set2.add(idx2);
                sumDist += costDist[idx];
            }
            count = matched.getN();
            for (int i = 0; i < count; ++i) {
                m1[i] = matched.getX(i);
                m2[i] = matched.getY(i);
            }
        }

        return new double[]{sumDist, count};
    }

    private static PairIntArray trimToImageBounds(TwoDFloatArray octaveImg, PairIntArray a) {
        int n0 = octaveImg.a.length;
        int n1 = octaveImg.a[0].length;
        return trimToImageBounds(n1, n0, a);
    }

    private static PairIntArray trimToImageBounds(
        int width, int height, PairIntArray a) {
        PairIntArray b = new PairIntArray(a.getN());
        for (int i = 0; i < a.getN(); ++i) {
            int x = a.getX(i);
            int y = a.getY(i);
            if (x < 0 || x > (width - 1)) {
                continue;
            } else if (y < 0 || y > (height - 1)) {
                continue;
            }
            b.add(x, y);
        }
        return b;
    }

    private static PairIntArray reduceBounds(PairIntArray bounds, float scale) {

        Set<PairInt> added = new HashSet<PairInt>();
        PairIntArray out = new PairIntArray(bounds.getN());
        for (int i = 0; i < bounds.getN(); ++i) {
            int x = Math.round((float)bounds.getX(i)/scale);
            int y = Math.round((float)bounds.getY(i)/scale);
            PairInt p = new PairInt(x, y);
            if (added.contains(p)) {
                continue;
            }
            out.add(x, y);
            added.add(p);
        }

        return out;
    }

    public static int maxNumberOfGaps(PairIntArray bounds,
        PartialShapeMatcher.Result r) {

        TIntSet mIdxs = new TIntHashSet(r.getNumberOfMatches());
        for (int i = 0; i < r.getNumberOfMatches(); ++i) {
            mIdxs.add(r.getIdx1(i));
        }

        int maxGapStartIdx = -1;
        int maxGap = 0;
        int cStartIdx = -1;
        int cGap = 0;

        // handling for startIdx of 0 to check for wraparound
        // of gap at end of block
        int gap0 = 0;

        for (int i = 0; i < bounds.getN(); ++i) {
            if (!mIdxs.contains(i)) {
                // is a gap
                if (cStartIdx == -1) {
                    cStartIdx = i;
                }
                cGap++;
                if (i == (bounds.getN() - 1)) {
                    if (gap0 > 0) {
                        // 0 1 2 3 4 5
                        // g g     g g
                        // gap0=2
                        // cGap=2 cStartIdx=4
                        if (cStartIdx > (gap0 - 1)) {
                            gap0 += cGap;
                        }
                    }
                    if (cGap > maxGap) {
                        maxGap = cGap;
                        maxGapStartIdx = cStartIdx;
                    }
                    if (gap0 > maxGap) {
                        maxGap = gap0;
                        maxGapStartIdx = 0;
                    }
                }
            } else {
                // is not a gap
                if (cStartIdx > -1) {
                    if (cGap > maxGap) {
                        maxGap = cGap;
                        maxGapStartIdx = cStartIdx;
                    }
                    if (cStartIdx == 0) {
                        gap0 = cGap;
                    }
                    cStartIdx = -1;
                    cGap = 0;
                }
            }
        }

        return maxGap;
    }

    private PairIntArray createOrderedBounds(ORB orb1,
        Set<PairInt> labeledPoints1, SIGMA sigma) {

        ImageProcessor imageProcessor = new ImageProcessor();

        Set<PairInt> set = new HashSet<PairInt>();
        set.addAll(labeledPoints1);

        PairIntArray bounds = imageProcessor.extractSmoothedOrderedBoundary(
            set, sigma,
            orb1.getPyramidImages().get(0).a[0].length,
            orb1.getPyramidImages().get(0).a.length);

        return bounds;
    }
    
    private PairIntArray createOrderedBoundsSansSmoothing(ORB orb1,
        Set<PairInt> labeledPoints1) {

        ImageProcessor imageProcessor = new ImageProcessor();

        Set<PairInt> set = new HashSet<PairInt>();
        set.addAll(labeledPoints1);

        PerimeterFinder2 finder = new PerimeterFinder2();
        PairIntArray ordered = finder.extractOrderedBorder(
            set);
        
        return ordered;
    }

    private PairIntArray getOrCreateOrderedBounds(TwoDFloatArray img,
        TIntObjectMap<PairIntArray> boundsMap, int segIdx,
        Set<PairInt> set, SIGMA sigma) {

        PairIntArray bounds = boundsMap.get(segIdx);
        if (bounds != null) {
            return bounds;
        }

        ImageProcessor imageProcessor = new ImageProcessor();

        bounds = imageProcessor.extractSmoothedOrderedBoundary(
            new HashSet<PairInt>(set), sigma, img.a[0].length, img.a.length);

        boundsMap.put(segIdx, bounds);

        {
            if (bounds.getN() > 1) {
                MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
                double[] xyCen = curveHelper.calculateXYCentroids(bounds);
                System.out.println("bounds center=" + (int)xyCen[0] + "," +
                    (int)xyCen[1] + " size_full=" +
                    calculateObjectSize(bounds));
            }
        }

        return bounds;
    }

    private PairIntArray getOrCreateOrderedBounds(TwoDFloatArray img,
        TIntObjectMap<PairIntArray> boundsMap, int segIdx,
        Set<PairInt> set) {

        PairIntArray bounds = boundsMap.get(segIdx);
        if (bounds != null) {
            return bounds;
        }

        Set<PairInt> set2 = new HashSet<PairInt>();
        set2.addAll(set);

        PerimeterFinder2 finder = new PerimeterFinder2();
        bounds = finder.extractOrderedBorder(set2);

        boundsMap.put(segIdx, bounds);

        {
            if (bounds.getN() > 1) {
                MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
                double[] xyCen = curveHelper.calculateXYCentroids(bounds);
                System.out.println("bounds center=" + (int)xyCen[0] + "," +
                    (int)xyCen[1] + " size_full=" +
                    calculateObjectSize(bounds));
            }
        }

        return bounds;
    }

    private double sumDistances(PairFloatArray distances) {

        double sum = 0;
        for (int i = 0; i < distances.getN(); ++i) {
            float d1 = distances.getX(i);
            float d2 = distances.getY(i);
            sum += Math.sqrt(d1 * d1 + d2 * d2);
        }

        return sum;
    }

    /**
     * from matched points in list 1 to list2, choose 2 pairs that have
     * small cost and large difference.
     * @param result
     * @param costD
     * @param xPoints1
     * @param yPoints1
     * @param xPoints2
     * @param yPoints2
     * @return
     */
    private QuadInt[] choose2ReferencePoints(Result result,
        PairIntArray bounds1, PairIntArray bounds2,
        TObjectIntMap<PairInt> point1KP1Map,
        TObjectIntMap<PairInt> point2KP1Map, int[][] costD) {

        int n = result.getNumberOfMatches();
        float[] costs = new float[n];
        PairInt[] points1 = new PairInt[n];
        int count = 0;
        for (int i = 0; i < n; ++i) {
            int idx1 = result.idx1s.get(i);
            int idx2 = result.idx2s.get(i);
            int x1 = bounds1.getX(idx1);
            int y1 = bounds1.getY(idx1);
            PairInt p1 = new PairInt(x1, y1);

            int x2 = bounds2.getX(idx2);
            int y2 = bounds2.getY(idx2);
            PairInt p2 = new PairInt(x2, y2);
            if (point1KP1Map.containsKey(p1) && point2KP1Map.containsKey(p2)) {
                int kpIdx1 = point1KP1Map.get(p1);
                int kpIdx2 = point2KP1Map.get(p2);
                costs[count] = costD[kpIdx1][kpIdx2];
                // these points are in full size reference frae
                points1[count] = p1;
                count++;
            }
        }

        if (count > 1) {
            if (count < n) {
                costs = Arrays.copyOf(costs, count);
                points1 = Arrays.copyOf(points1, count);
            }

            QuickSort.sortBy1stArg(costs, points1);

            int end = (int)(0.2 * n);
            if (end < 10) {
                end = n;
            }

            Set<PairInt> points = new HashSet<PairInt>();
            for (int i = 0; i < end; i++) {
                PairInt p = points1[i];
                points.add(p);
            }
            FurthestPair fp = new FurthestPair();
            PairInt[] furthest = fp.find(points);
            assert(furthest != null);
            assert(furthest.length == 2);
            PairInt[] furthest2 = new PairInt[2];

            for (int i = 0; i < n; ++i) {
                int idx1 = result.idx1s.get(i);
                int idx2 = result.idx2s.get(i);
                PairInt p1 = new PairInt(bounds1.getX(idx1),
                    bounds1.getY(idx1));
                if (furthest2[0] == null) {
                    if (furthest[0].equals(p1)) {
                        furthest2[0] = new PairInt(bounds2.getX(idx2),
                            bounds2.getY(idx2));
                    }
                }
                if (furthest2[1] == null) {
                    if (furthest[1].equals(p1)) {
                        furthest2[1] = new PairInt(bounds2.getX(idx2),
                            bounds2.getY(idx2));
                    }
                }
                if (furthest2[0] != null && furthest2[1] != null) {
                    break;
                }
            }
            if (furthest2 != null && furthest2.length == 2) {
                QuadInt[] refs = new QuadInt[2];
                refs[0] = new QuadInt(furthest[0], furthest2[0]);
                refs[1] = new QuadInt(furthest[1], furthest2[1]);
                return refs;
            }
        }

        // re-do the calculation w/o trying to use descr cost.
        Set<PairInt> points = new HashSet<PairInt>(n);
        for (int i = 0; i < n; ++i) {
            int idx1 = result.idx1s.get(i);
            PairInt p1 = new PairInt(bounds1.getX(idx1),
                bounds1.getY(idx1));
            points.add(p1);
        }
        FurthestPair fp = new FurthestPair();
        PairInt[] furthest = fp.find(points);
        assert(furthest != null);
        assert(furthest.length == 2);
        PairInt[] furthest2 = new PairInt[2];

        for (int i = 0; i < n; ++i) {
            int idx1 = result.idx1s.get(i);
            int idx2 = result.idx2s.get(i);
            PairInt p1 = new PairInt(bounds1.getX(idx1),
                bounds1.getY(idx1));
            if (furthest2[0] == null) {
                if (furthest[0].equals(p1)) {
                    furthest2[0] = new PairInt(bounds2.getX(idx2),
                        bounds2.getY(idx2));
                }
            }
            if (furthest2[1] == null) {
                if (furthest[1].equals(p1)) {
                    furthest2[1] = new PairInt(bounds2.getX(idx2),
                        bounds2.getY(idx2));
                }
            }
            if (furthest2[0] != null && furthest2[1] != null) {
                break;
            }
        }

        assert (furthest2 != null && furthest2.length == 2);

        QuadInt[] refs = new QuadInt[2];
        refs[0] = new QuadInt(furthest[0], furthest2[0]);
        refs[1] = new QuadInt(furthest[1], furthest2[1]);
        return refs;
    }

    private List<PairInt> matchUsingFM(ORB orb1, ORB orb2, int[][] costD,
        int octave1, int octave2,
        TObjectIntMap<PairInt> keypoints1IndexMap,
        TObjectIntMap<PairInt> keypoints2IndexMap,
        SimpleMatrix fm,
        PairIntArray unmatchedKP1, PairIntArray unmatchedKP2,
        TObjectIntMap<PairInt> unmatchedKP1Idxs,
        TObjectIntMap<PairInt> unmatchedKP2Idxs,
        int nBands, float distTol, double[] output) {

        float maxDesc = nBands * 256.0f;
        
        int distTolMax = (int)Math.round(Math.sqrt(2) * distTol);
        
        // output variable to hold sums and count
        // 0 = totalDistance
        // 1 = max avg total dist
        // 2 = totalDescrSum
        // 3 = nDescr

        double maxAvgDist = Double.MIN_VALUE;

        List<PairInt> addedKPIdxs = new ArrayList<PairInt>();

        EpipolarTransformer eTransformer = new EpipolarTransformer();

        SimpleMatrix unmatchedLeft =
            eTransformer.rewriteInto3ColumnMatrix(unmatchedKP1);

        SimpleMatrix unmatchedRight =
            eTransformer.rewriteInto3ColumnMatrix(unmatchedKP2);

        SimpleMatrix rightEpipolarLines = fm.mult(unmatchedLeft);
        SimpleMatrix leftEpipolarLines = fm.transpose().mult(unmatchedRight);

        float[] outputDist = new float[2];
        int nLeftUnmatched = unmatchedLeft.numCols();
        int nRightUnmatched = unmatchedRight.numCols();
        float dist, descCost;
        double d;

        TFloatList totalCost = new TFloatArrayList();
        TIntList indexes = new TIntArrayList();
        TFloatList eDist = new TFloatArrayList();
        TFloatList dCost = new TFloatArrayList();
        TIntList idx1s = new TIntArrayList();
        TIntList idx2s = new TIntArrayList();
        for (int i = 0; i < nLeftUnmatched; ++i) {

            PairInt p1 = new PairInt(unmatchedKP1.getX(i),
                unmatchedKP1.getY(i));
            int kp1Idx = unmatchedKP1Idxs.get(p1);

            for (int j = 0; j < nRightUnmatched; ++j) {

                PairInt p2 = new PairInt(unmatchedKP2.getX(j),
                    unmatchedKP2.getY(j));
                int kp2Idx = unmatchedKP2Idxs.get(p2);

                eTransformer.calculatePerpDistFromLines(unmatchedLeft,
                    unmatchedRight, rightEpipolarLines,
                    leftEpipolarLines, i, j, outputDist);
                
                if (outputDist[0] <= distTol && outputDist[1] <= distTol) {

                    d = Math.sqrt(outputDist[0] * outputDist[0] +
                        outputDist[1] * outputDist[1]);

                    dist = (float)d/distTolMax;

                    // normalized descriptor cost
                    descCost = 1.f - ((nBands * 256.f -
                        costD[kp1Idx][kp2Idx])
                        /maxDesc);

                    eDist.add(dist);
                    dCost.add(descCost);
                    totalCost.add((float)Math.sqrt(dist*dist + descCost*descCost));
                    indexes.add(indexes.size());
                    idx1s.add(kp1Idx);
                    idx2s.add(kp2Idx);
                }
            }
        }

        QuickSort.sortBy1stArg(totalCost, indexes);

        // new matched keypoint indexes
        TIntSet added1 = new TIntHashSet();
        TIntSet added2 = new TIntHashSet();
        for (int j = 0; j < totalCost.size(); ++j) {
            int idx = indexes.get(j);
            int kpIdx1 = idx1s.get(idx);
            int kpIdx2 = idx2s.get(idx);
            if (added1.contains(kpIdx1) || added2.contains(kpIdx2)) {
                continue;
            }
            // output variable to hold sums and count
            // 0 = totalDistance
            // 1 = max avg total dist
            // 2 = totalDescrSum
            // 3 = nDescr
            added1.add(kpIdx1);
            added2.add(kpIdx2);
            addedKPIdxs.add(new PairInt(kpIdx1, kpIdx2));
            output[0] += eDist.get(idx);

            d = eDist.get(idx);
            if (d > maxAvgDist) {
                maxAvgDist = d;
            }

            output[2] += dCost.get(idx);
            output[3]++;
        }

        output[1] = maxAvgDist;

        return addedKPIdxs;
    }

    // sum, avg, maxAvg
    private double[] sumAndMaxEPDist(SimpleMatrix fm, PairIntArray m1,
        PairIntArray m2) {

        EpipolarTransformer eTransformer = new EpipolarTransformer();

        SimpleMatrix matchedLeft =
            eTransformer.rewriteInto3ColumnMatrix(m1);
        SimpleMatrix matchedRight =
            eTransformer.rewriteInto3ColumnMatrix(m2);

        // this goes into total epipolar distances at end of block
        // NOTE: this method needs testing.
        //      currently no normalization is used internally
        PairFloatArray distances = eTransformer
            .calculateDistancesFromEpipolar(fm,
            matchedLeft, matchedRight);

        double max = Double.MIN_VALUE;
        double sum = 0;
        double d;
        for (int i = 0; i < distances.getN(); ++i) {
            float d1 = distances.getX(i);
            float d2 = distances.getY(i);
            d = Math.sqrt(d1*d1 + d2*d2);
            if (d > max) {
                max = d;
            }
            sum += d;
        }

        double avg = sum/(double)distances.getN();

        return new double[]{sum, avg, max};
    }

    /**
     assumptions such as transformation being near the
     expected scale of scale1/scale2 are made
     * @param segIdx
     * @param left
     * @param right
     * @param nBands
     * @param costD
     * @param nn1
     * @param keypoints1IndexMap point indexes lookup w.r.t. entire keypoints1 list
     *    which has same indexes as costD first dimension.  note that the
     *    template object dataset, a.k.a. left or a1 includes all points of
     *    keypoints1...no segmentation dividing the points into more than one
     *    array list. nn1 contains all keypoints1 points.
     * @param keypoints2IndexMap point indexes lookup w.r.t. entire keypoints2 list
     *     which has same indexes as costD second dimension.  these indexes are not the same
     *     indexes as a2 indexes.
     * @param outLeft
     * @param outRight
     * @param img1Width
     * @param img1Height
     * @param distTol
     * @param extr1
     * @param extr2
     * @param nnExtr1
     * @param extr1IndexMap
     * @param extr2IndexMap
     * @param scale1
     * @param scale2
     * @param outputNormalizedCost
     * @return 
     */
    private TransformationParameters matchGreedy(int segIdx,
        PairIntArray left, PairIntArray right,
        int nBands, int[][] costD, NearestNeighbor2D nn1,
        TObjectIntMap<PairInt> keypoints1IndexMap,
        TObjectIntMap<PairInt> keypoints2IndexMap,
        PairIntArray outLeft, PairIntArray outRight,
        int img1Width, int img1Height,
        int distTol,
        PairIntArray extr1, PairIntArray extr2, NearestNeighbor2D nnExtr1,
        TObjectIntMap<PairInt> extr1IndexMap,
        TObjectIntMap<PairInt> extr2IndexMap, float scale1, float scale2,
        double[] outputNormalizedCost) {

        // NOTE that all coordinates are in the full reference frame
        //   and remain that way through this method

        float expectedScale = scale1/scale2;
        
        float maxDesc = nBands * 256.0f;
                
        float distMax = (float)Math.sqrt(2) * distTol;
        /*
        -- finds best 20 matches of descriptors
        -- from the best 20,
           -- forms combinations of pairs of points
           -- for each pair,
              -- calc euclid transformation
                 evaluate it on all points
        result is best starting match of points.
        */

        int n1 = left.getN();
        int n2 = right.getN();

        /*
        // for each left, find the best matching right
        int[] b = new int[n1];
        int[] indexes = new int[n1];
        float[] costs = new float[n1];

        double c;
        for (int i = 0; i < n1; ++i) {

            double bestCost = Double.MAX_VALUE;
            int bestIdx2 = -1;

            PairInt p1 = new PairInt(left.getX(i), left.getY(i));
            int kpIdx1 = keypoints1IndexMap.get(p1);
            assert(keypoints1IndexMap.containsKey(p1));
        
            for (int j = 0; j < n2; ++j) {
                PairInt p2 = new PairInt(right.getX(j), right.getY(j));
                int kpIdx2 = keypoints2IndexMap.get(p2);
                assert(keypoints2IndexMap.containsKey(p2));
                c = costD[kpIdx1][kpIdx2];
                if (c < bestCost) {
                    bestCost = c;
                    bestIdx2 = j;
                }
            }
            b[i] = bestIdx2; // index w.r.t. right array
            indexes[i] = i;
            costs[i] = (float)bestCost;
        }
        */
        
        //int n12 = (n1 * n2) - ((n1 - 1) * Math.abs(n1 - n2));
        int n12 = n1 * n2;
        int[] idxs1 = new int[n12];
        int[] idxs2 = new int[n12];
        int[] indexes = new int[n12];
        float[] costs = new float[n12];

        int count = 0;
        double c;
        for (int i = 0; i < n1; ++i) {

            PairInt p1 = new PairInt(left.getX(i), left.getY(i));
            int kpIdx1 = keypoints1IndexMap.get(p1);
            assert(keypoints1IndexMap.containsKey(p1));
            
            for (int j = 0; j < n2; ++j) {
                PairInt p2 = new PairInt(right.getX(j), right.getY(j));
                int kpIdx2 = keypoints2IndexMap.get(p2);
                assert(keypoints2IndexMap.containsKey(p2));
                c = costD[kpIdx1][kpIdx2];
                idxs1[count] = i;
                idxs2[count] = j;
                indexes[count] = count;
                costs[count] = (float)c;
                count++;
            }
        }
        if (count < n12) {
            idxs1 = Arrays.copyOf(idxs1, count);
            idxs2 = Arrays.copyOf(idxs2, count);
            indexes = Arrays.copyOf(indexes, count);
            costs = Arrays.copyOf(costs, count);
        }
        
        QuickSort.sortBy1stArg(costs, indexes);

        int topK = 20;
        //if (topK > n1) {
            topK = count;
        //}

        MatchedPointsTransformationCalculator tc = new MatchedPointsTransformationCalculator();
        Transformer transformer = new Transformer();

        double bestCost = Double.MAX_VALUE;
        double bestCost2 = -1;
        int[] bestMIdx1s = null;
        int[] bestMIdx2s = null;
        int[] bestMIdxExtr1s = null;
        int[] bestMIdxExtr2s = null;
        int bestN = -1;
        int bestExtrN = -1;
        TransformationParameters bestParams = null;

        for (int i = 0; i < topK; ++i) {
            int idx = indexes[i];
            int idxI1 = idxs1[idx];
            int idxJ1 = idxs2[idx];
            int leftX1 = left.getX(idxI1);
            int leftY1 = left.getY(idxI1);
            int rightX1 = right.getX(idxJ1);
            int rightY1 = right.getY(idxJ1);

            for (int j = (i + 1); j < topK; ++j) {
                int idx2 = indexes[j];
                int idxI2 = idxs1[idx2];
                int idxJ2 = idxs2[idx2];
                int leftX2 = left.getX(idxI2);
                int leftY2 = left.getY(idxI2);
                int rightX2 = right.getX(idxJ2);
                int rightY2 = right.getY(idxJ2);

                assert(!(leftX1!=leftX2 && leftY1!= leftY1));
                assert(!(rightX1!=rightX2 && rightY1!= rightY1));
                
                // transform dataset 2 into frame 1
                // (direction is 2 to 1 to be able to reuse nearest
                //  neighbors containing dataset1 keypoints

                TransformationParameters params = tc.calulateEuclidean(
                    rightX1, rightY1, rightX2, rightY2,
                    leftX1, leftY1, leftX2, leftY2,
                    0, 0);

                if (params == null) {
                    continue;
                }

                float scale = params.getScale();

                if (((scale >= expectedScale)
                    && (scale/expectedScale) > 1.2) ||
                    ((scale < expectedScale)
                    && (expectedScale/scale) > 1.2)) {
                    continue;
                }

                PairIntArray rightTr = transformer.applyTransformation(
                    params, right);

                // filled with indexes w.r.t left and right arrays
                int[] mIdx1s = new int[n2];
                int[] mIdx2s = new int[n2];

                //sums is []{sumDesc, sumDist, count}
                // where sumDesc is the sum of normalized descriptors
                //       each with a value of 0 to 1 and summed over 
                //       count number of them. (so max is <= count)
                // same for sumDist
                double[] sums = sumKeypointDescAndDist2To1(costD, nBands,
                    right, rightTr, nn1,
                    keypoints1IndexMap, keypoints2IndexMap,
                    distTol, mIdx1s, mIdx2s);

                final int nMatched = (int)sums[2];
               
                //TODO: may need to revise this.  essentially, need at least
                // one extra pair to get a distance term
                if (nMatched < 3) {
                    continue;
                }
                
                double sumDescr = sums[0];
                double sumDist = sums[1];
                
                float descCost = (float)sumDescr/(float)nMatched;
                float distCost = (float)sumDist/(float)nMatched;
                float d1 = 1.f - ((nBands * 256.f - descCost)/maxDesc);
                float d2 = distCost;
                final float f1 = 1.f - ((float)nMatched/(float)n1);

                float tot = d1*d1 + d2*d2 + 2.f * f1 * f1;
                
                // evaluate the extra points
                int ne1 = extr1.getN();
                int ne2 = extr2.getN();

                PairIntArray extr2Tr = transformer.applyTransformation(
                    params, extr2);

                // since the extra points do not have descriptors,
                // need to decide between adding a ne1 to sumDescr
                // or 0 to it (changing the importance of the orid keypoints).
                
                int[] mIdxExtr1s = new int[ne2];
                int[] mIdxExtr2s = new int[ne2];

                double d3 = 1;
                double f3 = 1;
                int nMatchedExtr = 0;
                if (extr2Tr.getN() > 0) {
          
                    //sumsExtr = []{sumDist, count}
                    double[] sumsExtr = sumKeypointDescAndDist2To1(
                        extr2, extr2Tr, nnExtr1,
                        extr1IndexMap, extr2IndexMap,
                        img1Width, img1Height,
                        distTol, mIdxExtr1s, mIdxExtr2s);

                    nMatchedExtr = (int)sumsExtr[1];
                    double distCostExtr = sumsExtr[0]/(float)nMatchedExtr;
                    d3 = distCostExtr;
                    f3 = 1.f - ((float)nMatchedExtr/(float)ne1);
                }

                if (tot < bestCost) {
                    
                    bestCost = tot;
                    bestCost2 = d3*d3 + f3*f3;
       
                    // 0 = the salukwzde distance, that is the normalized tot cost
                    // 1 = sum of normalized keypoint descriptors
                    // 2 = sum of normalized keypoint distances from transformations
                    // 3 = number of keypoint matches (not incl boundary that aren't
                    outputNormalizedCost[0] = bestCost;
                    outputNormalizedCost[1] = sumDescr;
                    outputNormalizedCost[2] = sumDist;
                    outputNormalizedCost[3] = nMatched;
                        
                    bestMIdx1s = mIdx1s;
                    bestMIdx2s = mIdx2s;
                    bestN = nMatched;

                    bestMIdxExtr1s = mIdxExtr1s;
                    bestMIdxExtr2s = mIdxExtr2s;
                    bestExtrN = nMatchedExtr;

                    bestParams = params;
                }
            }
        }

        if (bestMIdx1s != null) {
            Set<PairInt> addedLeft = new HashSet<PairInt>();
            Set<PairInt> addedRight = new HashSet<PairInt>();
            for (int i = 0; i < bestN; ++i) {
                int kpIdx1 = bestMIdx1s[i];
                int kpIdx2 = bestMIdx2s[i];
                PairInt pLeft = new PairInt(left.getX(kpIdx1), left.getY(kpIdx1));
                PairInt pRight = new PairInt(right.getX(kpIdx2), right.getY(kpIdx2));
                if (addedLeft.contains(pLeft) || addedRight.contains(pRight)) {
                    continue;
                }
                addedLeft.add(pLeft);
                addedRight.add(pRight);
                outLeft.add(pLeft.getX(), pLeft.getY());
                outRight.add(pRight.getX(), pRight.getY());
            }
            for (int i = 0; i < bestExtrN; ++i) {
                int idx1 = bestMIdxExtr1s[i];
                int idx2 = bestMIdxExtr2s[i];
                PairInt pLeft = new PairInt(extr1.getX(idx1), extr1.getY(idx1));
                PairInt pRight = new PairInt(extr2.getX(idx2), extr2.getY(idx2));
                if (addedLeft.contains(pLeft) || addedRight.contains(pRight)) {
                    continue;
                }
                addedLeft.add(pLeft);
                addedRight.add(pRight);
                outLeft.add(pLeft.getX(), pLeft.getY());
                outRight.add(pRight.getX(), pRight.getY());
            }
        }

        return bestParams;
    }

    private double[] addUnmatchedKeypoints(TransformationParameters params, 
        PairIntArray m1, PairIntArray m2, int nBands, int[][] costD, 
        PairIntArray keypoints1, TObjectIntMap<PairInt> keypoints1IndexMap, 
        TObjectIntMap<PairInt> keypoints2IndexMap, 
        int imageWidth, int imageHeight, int distTol, 
        float scale1, float scale2) {
        
        if (m1.getN() == keypoints1.getN()) {
            // all keypoints have been matched
            return new double[]{0, 0, 0};
        }
        
        float maxDesc = nBands * 256.0f;
    
        float distTolMax = (float)Math.sqrt(2) * distTol;
        
        // find all transformed keypoints2 that are unmatched, and match
        // a left unmatched point within distTol.
        // then sort those by total cost and assign uniquely.
        // then add results to m1, m2 and
        // put the costs in return array.
        
        /*
        output:
        // 0 = sum of normalized keypoint descriptors
        // 1 = sum of normalized keypoint distances from transformations
        // 2 = number of keypoint matches added
        */
    
        // find unmatched keypoints1 and put in a nearest neighbor instance
        Set<PairInt> matched1 = new HashSet<PairInt>();
        Set<PairInt> matched2 = new HashSet<PairInt>();
        for (int j = 0; j < m1.getN(); ++j) {
            // NOTE: m1,m2 are keypoints and boundary points
            int x1 = m1.getX(j);
            int y1 = m1.getY(j);
            int x2 = m2.getX(j);
            int y2 = m2.getY(j);
            matched1.add(new PairInt(x1, y1));
            matched2.add(new PairInt(x2, y2));
        }
        PairIntArray unmatched1 = new PairIntArray();
        for (int j = 0; j < keypoints1.getN(); ++j) {
            PairInt p = new PairInt(keypoints1.getX(j), keypoints1.getY(j));
            if (!matched1.contains(p)) {
                unmatched1.add(p.getX(), p.getY());
            }
        }
        
        int[] minMaxXY1 = MiscMath.findMinMaxXY(unmatched1);
        NearestNeighbor2D nn1 = new NearestNeighbor2D(Misc.convert(unmatched1),
            minMaxXY1[1] + distTol + 2,
            minMaxXY1[3] + distTol + 2);
        
        List<PairInt> mc1 = new ArrayList<PairInt>();
        List<PairInt> mc2 = new ArrayList<PairInt>();
        TFloatList descCosts = new TFloatArrayList();
        TFloatList distances = new TFloatArrayList();
        TFloatList totalCosts = new TFloatArrayList();
        TIntList indexes = new TIntArrayList();

        Transformer transformer = new Transformer();
        
        // visit all keypoint2 and if not matched, transform point
        //    and find nearest unmatched1 neighbor
        TObjectIntIterator<PairInt> iter2 = keypoints2IndexMap.iterator();
        for (int i = 0; i < keypoints2IndexMap.size(); ++i) {
            iter2.advance();
            PairInt p2 = iter2.key();
            if (matched2.contains(p2)) {
                continue;
            }
            double[] tr = transformer.applyTransformation(params, p2.getX(), 
                p2.getY());
            int x2Tr = (int) Math.round(tr[0]);
            int y2Tr = (int) Math.round(tr[1]);
            if (x2Tr < 0 || y2Tr < 0 || (x2Tr > (imageWidth - 1)) ||
                (y2Tr > (imageHeight - 1))) {
                continue;
            }
            PairInt p2Tr = new PairInt(x2Tr, y2Tr);
           
            Set<PairInt> nearest = nn1.findClosest(x2Tr, y2Tr, distTol);
            
            if (nearest == null || nearest.isEmpty()) {
                continue;
            }
            
            int kpIdx2 = keypoints2IndexMap.get(p2);
            
            PairInt p1Closest = null;
            if (nearest.size() == 1) {
                p1Closest = nearest.iterator().next();
            } else {
                // find closest in terms of descriptor
                float minCost = Float.MAX_VALUE;
                PairInt minCostP = null;
                for (PairInt p3 : nearest) {
                    int kpIdx1 = keypoints1IndexMap.get(p3);
                    float c = costD[kpIdx1][kpIdx2];
                    if (c < minCost) {
                        minCost = c;
                        minCostP = p3;
                    }
                }
                assert(minCostP != null);
                p1Closest = minCostP;
            }
            
            assert(p1Closest != null);
            
            int kpIdx1 = keypoints1IndexMap.get(p1Closest);
            float c = costD[kpIdx1][kpIdx2];
            float dist = distance(p2Tr, p1Closest);
            assert(dist <= distTolMax);
            
            float costNorm = 1.f - ((nBands * 256 - c) / maxDesc);
            float distNorm = dist / distTolMax;
        
            assert(!matched1.contains(p1Closest));
            assert(!matched2.contains(p2));
            
            mc1.add(p1Closest);
            mc2.add(p2);
            descCosts.add(costNorm);
            distances.add(distNorm);
            // could square these, but not necessary for this comparison:
            totalCosts.add(costNorm + distNorm);
            indexes.add(indexes.size());
        }
        
        QuickSort.sortBy1stArg(totalCosts, indexes);
        
        Set<PairInt> added1 = new HashSet<PairInt>();
        Set<PairInt> added2 = new HashSet<PairInt>();
        double sumDist = 0;
        double sumDesc = 0;
        int nAdded = 0;
        
        // adding directly to m1 and m2, uniquely
        for (int i = 0; i < indexes.size(); ++i) {
            int idx = indexes.get(i);
            PairInt p1 = mc1.get(idx);
            PairInt p2 = mc2.get(idx);
            if (added1.contains(p1) || added2.contains(p2)) {
                continue;
            }
            added1.add(p1);
            added2.add(p2);
            
            m1.add(p1.getX(), p1.getY());
            m2.add(p2.getX(), p2.getY());
            
            sumDesc += descCosts.get(idx);
            sumDist += distances.get(idx);
            nAdded++;
        }
        
        return new double[] {sumDesc, sumDist, nAdded};
    }

    private TIntObjectMap<ShapeFinder2.ShapeFinderResult> aggregatedShapeMatch(
        ORB orb1, ORB orb2,
        Set<PairInt> labeledPoints1, List<Set<PairInt>> labeledPoints2, 
        TIntList octs1, TIntList octs2, TIntList segIdxs, 
        TFloatList scales1, TFloatList scales2) {

        System.out.println("aggregatedShapeMatch");
        
        TIntObjectMap<ShapeFinder2.ShapeFinderResult> resultMap = new
            TIntObjectHashMap<ShapeFinder2.ShapeFinderResult>();
        
        TIntObjectMap<Set<PairInt>> octave1ScaledSets =
            createSetsForOctaves(octs1, labeledPoints1, scales1);
        
        TIntObjectMap<PairIntArray> octave1ScaledBounds = 
            createBounds(octave1ScaledSets, 
                orb1.getPyramidImages().get(0).a[0].length,
                orb1.getPyramidImages().get(0).a.length);
        
        TIntObjectMap<TIntObjectMap<Set<PairInt>>> octave2ScaledSets =
            createSetsForOctaves(octs2, labeledPoints2, scales2);
        
        TIntObjectMap<Map<OneDIntArray, PairIntArray>> keybounds2Maps
            = initializeMaps(octs2);
        
        int n = octs1.size();
        
        for (int i = 0; i < n; ++i) {
            int octave1 = octs1.get(i);
            int octave2 = octs2.get(i);
            int segIdx = segIdxs.get(i);
            float scale1 = scales1.get(octave1);
            float scale2 = scales2.get(octave2);
            
            PairIntArray bounds1 = octave1ScaledBounds.get(octave1);
            
            float sz1 = calculateObjectSize(bounds1);
            
            TIntObjectMap<Set<PairInt>> mapOfSets2 = new
                TIntObjectHashMap<Set<PairInt>>(octave2ScaledSets.get(octave2));
        
            // removing any set that, when combined with set for segIdx,
            // has dimensions larger than twice the size of bounds1
            removeSetsByCombinedSize(bounds1, segIdx, mapOfSets2);
            
            removeSetsNotConnected(segIdx, mapOfSets2);
            
            Map<OneDIntArray, PairIntArray> keyBoundsMap2 = 
                keybounds2Maps.get(octave2);
            
            int xMax1 = orb1.getPyramidImages().get(octave1).a[0].length -1;
            int yMax1 = orb1.getPyramidImages().get(octave1).a.length - 1;
            
            int xMax2 = orb2.getPyramidImages().get(octave2).a[0].length -1;
            int yMax2 = orb2.getPyramidImages().get(octave2).a.length - 1;
            
            System.out.println("start shape search for octave1=" + octave1 + 
                " octave2=" + octave2 + " segIdx=" + segIdx);
            
   //TODO: this needs improvements and the
   //    costs updated for use here
            
            ShapeFinder2 shapeFinder = new ShapeFinder2(bounds1, scale1, sz1, 
                xMax1, yMax1, mapOfSets2, scale2, 
                keyBoundsMap2, xMax2, yMax2);
            
            ShapeFinder2.ShapeFinderResult result = shapeFinder.findAggregated();
        
            if (result == null) {
                continue;
            }
            
            System.out.println("shapeFinder result for segIdx=" + 
                + segIdx + " " + result.toString());
            
            resultMap.put(i, result);
            
            {//DEBUG
                ShapeFinder2.ShapeFinderResult r = result;
                Image img1 = ORB.convertToImage(
                    orb1.getPyramidImages().get(octave1));
                Image img2 = ORB.convertToImage(
                    orb2.getPyramidImages().get(octave2));
                try {
                    CorrespondencePlotter plotter = new CorrespondencePlotter(
                        r.bounds1, r.bounds2);
                        //img1, img2);
                    for (int ii = 0; ii < r.getNumberOfMatches(); ++ii) {
                        int idx1 = r.getIdx1(ii);
                        int idx2 = r.getIdx2(ii);
                        int x1 = r.bounds1.getX(idx1);
                        int y1 = r.bounds1.getY(idx1);
                        int x2 = r.bounds2.getX(idx2);
                        int y2 = r.bounds2.getY(idx2);
                        if ((ii % 4) == 0) {
                            plotter.drawLineInAlternatingColors(x1, y1, x2, y2, 0);
                        }
                    }
                    
                    String str = octave1 + "_" + octave2 + "_" + segIdx 
                        + MiscDebug.getCurrentTimeFormatted();
                    
                    String filePath = plotter.writeImage("_shape_" + str);
                    
                } catch (Throwable t) {
                    System.err.println(t.getMessage());
                }
            }
        }
        
        return resultMap;
    }

    private TIntObjectMap<Set<PairInt>> createSetsForOctaves(TIntList octaves, 
        Set<PairInt> labeledPoints, TFloatList scales) {

        TIntObjectMap<Set<PairInt>> output = new TIntObjectHashMap<Set<PairInt>>();
    
        for (int i = 0; i < octaves.size(); ++i) {
            int octave = octaves.get(i);
            if (output.containsKey(octave)) {
                continue;
            }
            float scale = scales.get(octave);
            
            Set<PairInt> set2 = new HashSet<PairInt>();
            for (PairInt p : labeledPoints) {
                int x = Math.round((float)p.getX()/scale);
                int y = Math.round((float)p.getY()/scale);
                set2.add(new PairInt(x, y));
            }
            output.put(octave, set2);
        }
        
        return output;
    }

    private TIntObjectMap<TIntObjectMap<Set<PairInt>>> createSetsForOctaves(
        TIntList octaves, List<Set<PairInt>> labeledPoints, TFloatList scales) {

        TIntObjectMap<TIntObjectMap<Set<PairInt>>> output = new 
            TIntObjectHashMap<TIntObjectMap<Set<PairInt>>>();
        
        for (int i = 0; i < octaves.size(); ++i) {
            int octave = octaves.get(i);
            if (output.containsKey(octave)) {
                continue;
            }
            float scale = scales.get(octave);
        
            TIntObjectMap<Set<PairInt>> labeledOctaveMap = new 
                TIntObjectHashMap<Set<PairInt>>();
            
            int n = labeledPoints.size();
            for (int segIdx = 0; segIdx < n; ++segIdx) {
                Set<PairInt> set = labeledPoints.get(segIdx);
                Set<PairInt> set2 = new HashSet<PairInt>();
                for (PairInt p : set) {
                    int x = Math.round((float)p.getX()/scale);
                    int y = Math.round((float)p.getY()/scale);
                    set2.add(new PairInt(x, y));
                }
                labeledOctaveMap.put(segIdx, set2);
            }
            
            output.put(octave, labeledOctaveMap);
        }
        
        return output;
    }

    private TIntObjectMap<PairIntArray> createBounds(
        TIntObjectMap<Set<PairInt>> octaveScaledSets, int imgWidth, int imgHeight) {
        
        TIntObjectMap<PairIntArray> output = new TIntObjectHashMap<PairIntArray>();
        
        SIGMA sigma = SIGMA.ZEROPOINTFIVE;
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        TIntObjectIterator<Set<PairInt>> iter = octaveScaledSets.iterator();
        for (int i = 0; i < octaveScaledSets.size(); ++i) {
            iter.advance();
            int segIdx = iter.key();
            Set<PairInt> set = iter.value();
            
            PairIntArray bounds = imageProcessor.extractSmoothedOrderedBoundary(
                new HashSet<PairInt>(set), sigma, imgWidth, imgHeight);
        
            output.put(segIdx, bounds);
        }
    
        return output;
    }

    private TIntObjectMap<Map<OneDIntArray, PairIntArray>> initializeMaps(
        TIntList octaves) {
        
        TIntObjectMap<Map<OneDIntArray, PairIntArray>> output = new
             TIntObjectHashMap<Map<OneDIntArray, PairIntArray>>();
    
        for (int i = 0; i < octaves.size(); ++i) {
            int octave = octaves.get(i);
            output.put(octave, new HashMap<OneDIntArray, PairIntArray>());
        }
        
        return output;
    }

    private void removeSetsByCombinedSize(PairIntArray bounds1, 
        int segIdx2, TIntObjectMap<Set<PairInt>> mapOfBounds2) {

        float sz1 = calculateObjectSize(bounds1);
        
        TIntSet rm = new TIntHashSet();
        
        Set<PairInt> set2 = mapOfBounds2.get(segIdx2);
        
        TIntObjectIterator<Set<PairInt>> iter = mapOfBounds2.iterator();
        for (int i = 0; i < mapOfBounds2.size(); ++i) {
            iter.advance();
            int segIdx3 = iter.key();
            Set<PairInt> set3 = iter.value();
            Set<PairInt> combined = new HashSet<PairInt>(set2);
            combined.addAll(set3);
            
            float sz2 = calculateObjectSize(combined);
            
            if (sz2 > 2.* sz1) {
                rm.add(segIdx3);
            }
        }
        
        TIntIterator iter2 = rm.iterator();
        while (iter2.hasNext()) {
            int rmIdx = iter2.next();
            mapOfBounds2.remove(rmIdx);
        }        
    }

    private void removeSetsNotConnected(int segIdx, 
        TIntObjectMap<Set<PairInt>> mapOfSets) {
        
        TObjectIntMap<PairInt> pointIndexMap = new TObjectIntHashMap<PairInt>();
        TIntObjectIterator<Set<PairInt>> iter = mapOfSets.iterator();
        for (int i = 0; i < mapOfSets.size(); ++i) {
            iter.advance();
            int idx = iter.key();
            for (PairInt p : iter.value()) {
                pointIndexMap.put(p, idx);
            }
        }
        
        TIntSet connected = new TIntHashSet();

        Stack<Integer> stack = new Stack<Integer>();
        stack.add(Integer.valueOf(segIdx));
        
        int[] dxs = Misc.dx4;
        int[] dys = Misc.dy4;
        
        TIntSet addedToStack = new TIntHashSet();
        TIntSet visited = new TIntHashSet();
        
        while (!stack.isEmpty()) {
            
            int idx = stack.pop().intValue();
            
            if (visited.contains(idx)) {
                continue;
            }
            
            connected.add(idx);
            
            Set<PairInt> set = mapOfSets.get(idx);
            
            for (PairInt p : set) {
                int x = p.getX();
                int y = p.getY();
                for (int k = 0; k < dxs.length; ++k) {
                    int x2 = x + dxs[k];
                    int y2 = y + dys[k];
                    PairInt p2 = new PairInt(x2, y2);
                    if (pointIndexMap.containsKey(p2)) {
                        
                        int idx2 = pointIndexMap.get(p2);
                        
                        if (!addedToStack.contains(idx2) && !visited.contains(idx2)) {
                            addedToStack.add(idx2);
                            stack.add(Integer.valueOf(idx2));
                        }
                    }
                }
            }
            
            visited.add(idx);
        }
        
        /*
        // looks like this removes the item, but just in case the
        //    values are left in the heap, will remove singly
        //    until can test for it.
        int n = mapOfSets.size();
        TIntSet keys = mapOfSets.keySet();
        keys.removeAll(connected);
        int n2 = mapOfSets.size();
        */
        
        TIntSet rm = new TIntHashSet(mapOfSets.keySet());
        rm.removeAll(connected);
        
        TIntIterator iter2 = rm.iterator();
        while (iter2.hasNext()) {
            int rmIdx = iter2.next();
            mapOfSets.remove(rmIdx);
        }
    }

    private static class PObject {

        final PartialShapeMatcher.Result r;
        final PairIntArray bounds1;
        final PairIntArray bounds2;
        final float scale1;
        final float scale2;

        public PObject(PartialShapeMatcher.Result result, PairIntArray b1, PairIntArray b2,
            float s1, float s2) {
            r = result;
            bounds1 = b1;
            bounds2 = b2;
            scale1 = s1;
            scale2 = s2;
        }
    }

    private static class CObject implements Comparable<CObject> {

        final double cost;
        final CorrespondenceList cCor;
        final PairIntArray transformedTemplate;

        public CObject(double cost, CorrespondenceList cL,
            PairIntArray templTr) {
            this.cost = cost;
            this.cCor = cL;
            this.transformedTemplate = templTr;
        }

        @Override
        public int compareTo(CObject other) {
            if (cost < other.cost) {
                return -1;
            } else if (cost > other.cost) {
                return 1;
            } else {
                int n1 = cCor.getPoints1().size();
                int n2 = other.cCor.getPoints1().size();
                if (n1 > n2) {
                    return -1;
                } else if (n1 < n2) {
                    return 1;
                }
            }
            return 0;
        }
    }

    private static class CObject4 implements Comparable<CObject4> {

        final double cost;
        final TransformationParameters params;
        final QuadInt q;

        public CObject4(double sum, TransformationParameters params,
            QuadInt q) {
            this.cost = sum;
            this.q = q;
            this.params = params;
        }

        @Override
        public int compareTo(CObject4 other) {
            if (cost < other.cost) {
                return -1;
            } else if (cost > other.cost) {
                return 1;
            }
            return 0;
        }
    }

    private static class CObject3 implements Comparable<CObject3> {

        final double cost;
        final double costDesc;
        final double costDist;
        final double costCount;
        final int index;
        final PairInt[] m1;
        final PairInt[] m2;
        final double sumPatch;
        final TransformationParameters params;
        QuadInt q;
        int keypointCount;

        public CObject3(CObject2 cObject2, double sum, double sumPatch,
            TransformationParameters params) {
            this.sumPatch = sumPatch;
            this.cost = sum;
            this.costDesc = cObject2.costDesc;
            this.costDist = cObject2.costDist;
            this.costCount = cObject2.costCount;
            this.index = cObject2.index;
            this.m1 = cObject2.m1;
            this.m2 = cObject2.m2;
            this.params = params;
        }

        @Override
        public int compareTo(CObject3 other) {
            if (cost < other.cost) {
                return -1;
            } else if (cost > other.cost) {
                return 1;
            }
            return 0;
        }
    }

    private static class CObject2 implements Comparable<CObject2> {

        final double cost;
        final double costDesc;
        final double costDist;
        final double costCount;
        final int index;
        final PairInt[] m1;
        final PairInt[] m2;

        public CObject2(int index, double cost, double costDesc, double costDist,
            double costCount, PairInt[] matched1, PairInt[] matched2) {
            this.cost = cost;
            this.index = index;
            this.m1 = matched1;
            this.m2 = matched2;
            this.costDesc = costDesc;
            this.costDist = costDist;
            this.costCount = costCount;
        }

        @Override
        public int compareTo(CObject2 other) {
            if (cost < other.cost) {
                return -1;
            } else if (cost > other.cost) {
                return 1;
            }
            return 0;
        }
    }

    private static void debugPlot(int i, int j,
        FixedSizeSortedVector<CObject> vec,
        TwoDFloatArray pyr1, TwoDFloatArray pyr2, float s1, float s2) {

        Image img1 = convertToImage(pyr1);
        Image img2 = convertToImage(pyr2);

        try {
            CorrespondenceList cor = vec.getArray()[0].cCor;

            Image img1Cp = img1.copyImage();
            Image img2Cp = img2.copyImage();
            CorrespondencePlotter plotter = new CorrespondencePlotter(
                img1Cp, img2Cp);
            for (int ii = 0; ii < cor.getPoints1().size(); ++ii) {
                PairInt p1 = cor.getPoints1().get(ii);
                PairInt p2 = cor.getPoints2().get(ii);
                int x1 = Math.round(p1.getX() / s1);
                int y1 = Math.round(p1.getY() / s1);
                int x2 = Math.round(p2.getX() / s2);
                int y2 = Math.round(p2.getY() / s2);
                plotter.drawLineInAlternatingColors(x1, y1, x2, y2, 0);
            }
            String strI = Integer.toString(i);
            while (strI.length() < 3) {
                strI = "0" + strI;
            }
            String strJ = Integer.toString(j);
            while (strJ.length() < 3) {
                strJ = "0" + strJ;
            }
            String str = strI + "_" + strJ + "_";
            plotter.writeImage("_MATCH_" + str);
        } catch (Exception e) {
        }
    }

    private static void debugPlot2(int i, int j,
        FixedSizeSortedVector<CObject3> vec,
        TwoDFloatArray pyr1, TwoDFloatArray pyr2, float s1, float s2) {

        Image img1 = convertToImage(pyr1);
        Image img2 = convertToImage(pyr2);

        try {
            PairInt[] m1 = vec.getArray()[0].m1;
            PairInt[] m2 = vec.getArray()[0].m2;

            Image img1Cp = img1.copyImage();
            Image img2Cp = img2.copyImage();
            CorrespondencePlotter plotter = new CorrespondencePlotter(
                img1Cp, img2Cp);
            for (int ii = 0; ii < m1.length; ++ii) {
                PairInt p1 = m1[ii];
                PairInt p2 = m2[ii];
                int x1 = Math.round(p1.getX() / s1);
                int y1 = Math.round(p1.getY() / s1);
                int x2 = Math.round(p2.getX() / s2);
                int y2 = Math.round(p2.getY() / s2);
                plotter.drawLineInAlternatingColors(x1, y1, x2, y2, 0);
            }
            String strI = Integer.toString(i);
            while (strI.length() < 3) {
                strI = "0" + strI;
            }
            String strJ = Integer.toString(j);
            while (strJ.length() < 3) {
                strJ = "0" + strJ;
            }
            String str = strI + "_" + strJ + "_";
            plotter.writeImage("_MATCH_" + str);
        } catch (Exception e) {
        }
    }

    private static void debugPlot(int i, int j,
        FixedSizeSortedVector<CObject> vecD,
        FixedSizeSortedVector<CObject> vecF,
        TwoDFloatArray pyr1, TwoDFloatArray pyr2, float s1, float s2) {

        Image img1 = convertToImage(pyr1);
        Image img2 = convertToImage(pyr2);

        try {
            for (int i0 = 0; i0 < 2; ++i0) {
                CorrespondenceList cor = null;
                if (i0 == 0) {
                    cor = vecD.getArray()[0].cCor;
                } else {
                    cor = vecF.getArray()[0].cCor;
                }
                Image img1Cp = img1.copyImage();
                Image img2Cp = img2.copyImage();
                CorrespondencePlotter plotter = new CorrespondencePlotter(
                    img1Cp, img2Cp);
                for (int ii = 0; ii < cor.getPoints1().size(); ++ii) {
                    PairInt p1 = cor.getPoints1().get(ii);
                    PairInt p2 = cor.getPoints2().get(ii);
                    int x1 = Math.round(p1.getX()/s1);
                    int y1 = Math.round(p1.getY()/s1);
                    int x2 = Math.round(p2.getX()/s2);
                    int y2 = Math.round(p2.getY()/s2);
                    plotter.drawLineInAlternatingColors(x1, y1, x2, y2, 0);
                }
                String strI = Integer.toString(i);
                while (strI.length() < 3) {
                    strI = "0" + strI;
                }
                String strJ = Integer.toString(j);
                while (strJ.length() < 3) {
                    strJ = "0" + strJ;
                }
                String str = strI + "_" + strJ + "_";
                if (i0 == 0) {
                    str = str + "factor";
                } else {
                    str = str + "divisor";
                }
                plotter.writeImage("_MATCH_" + str);
            }
        } catch(Exception e) {}
    }
}
