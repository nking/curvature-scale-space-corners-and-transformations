package algorithms.imageProcessing.matching;

import algorithms.QuickSort;
import algorithms.compGeometry.FurthestPair;
import algorithms.compGeometry.PerimeterFinder2;
import algorithms.compGeometry.PointInPolygon;
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
import algorithms.imageProcessing.matching.PartialShapeMatcher.Result;
import algorithms.imageProcessing.matching.ShapeFinder.ShapeFinderResult;
import algorithms.imageProcessing.transform.EpipolarTransformer;
import algorithms.imageProcessing.transform.MatchedPointsTransformationCalculator;
import algorithms.imageProcessing.transform.TransformationParameters;
import algorithms.imageProcessing.transform.Transformer;
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
import algorithms.util.TrioInt;
import algorithms.util.TwoDFloatArray;
import algorithms.util.TwoDIntArray;
import algorithms.util.VeryLongBitString;
import gnu.trove.iterator.TIntFloatIterator;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.iterator.TObjectIntIterator;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.TFloatList;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TIntFloatMap;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.hash.TIntFloatHashMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.Stack;
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

        int distTol = 10;//20;//15;//5

        PairIntArray bounds1 = createOrderedBoundsSansSmoothing(
            labeledPoints1);
        if (bounds1.getN() < 7) {
            throw new IllegalStateException("the boundary of object 1 "
                + " must have at least 7 points");
        }
        
        //TODO: this needs a setting that allows a filter to be made to
        //   prefer keypoints furthest from the outer bounds (==inner points).
        //   This is helpful for matching objects which have different
        //   backgrounds or foregrounds because they have changed location
        //   or pose.
        //   -- all keypoints can be used for the geometric model,
        //      but for costs in the last results, might need to restrict
        //      the costs to inner descriptors.

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

        TObjectIntMap<PairInt> point2LabelMap = new TObjectIntHashMap<PairInt>();
        for (int i = 0; i < labeledPoints2.size(); ++i) {
            for (PairInt p : labeledPoints2.get(i)) {
                point2LabelMap.put(p, i);
            }
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
        Map<OneDIntArray, PairIntArray> keyBounds2Map =
            new HashMap<OneDIntArray, PairIntArray>();

        // index = index of labeledPoints, item = bitstring with neighbors set
        TIntObjectMap<VeryLongBitString> label2AdjacencyMap
            = createAdjacencyMap(point2LabelMap, labeledPoints2);
       
        //key = ordered pair of adjacent label2s, value=size of combined regions
        TObjectIntMap<PairInt> label2PairSizes = new TObjectIntHashMap<PairInt>();
        
        // keys for these data are dataCount
        int dataCount = 0;
        TIntObjectMap<List<QuadInt>> correspondences = new 
            TIntObjectHashMap<List<QuadInt>>();
        TIntObjectMap<TIntList> matchedLabels2 = new TIntObjectHashMap<TIntList>();
        TIntObjectMap<TransformationParameters> transformations = new 
            TIntObjectHashMap<TransformationParameters>();
        TIntObjectMap<Double> descCosts = new TIntObjectHashMap<Double>();
        TIntIntMap nDesc = new TIntIntHashMap();
        TIntObjectMap<Double> distCosts = new TIntObjectHashMap<Double>();
        TIntFloatMap areaFraction = new TIntFloatHashMap();
        TIntIntMap octs1 = new TIntIntHashMap();
        TIntIntMap octs2 = new TIntIntHashMap();
        TIntIntMap segIdxs = new TIntIntHashMap();

        for (int octave1 = 0; octave1 < scales1.size(); ++octave1) {
        //for (int octave1 = 0; octave1 < 1; ++octave1) {

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
            //for (int octave2 = 2; octave2 < 3; ++octave2) {

                float scale2 = scales2.get(octave2);

                // key = keypoint coords, value = index of keypoint within
                //     the orb keypoint list for this octave2
                TObjectIntMap<PairInt> keypoints2IndexMap
                    = kp2IdxMapList.get(octave2);
                TIntObjectMap<TIntSet> labels2KPIdxs =
                    labels2KPIdxsList.get(octave2);

                //key = keypoint coords, value = label that keypoint is within
                TObjectIntMap<PairInt> keypoints2LabelMap = new
                    TObjectIntHashMap<PairInt>();
                {
                    TIntObjectIterator<TIntSet> iter = labels2KPIdxs.iterator();
                    for (int ii = 0; ii < labels2KPIdxs.size(); ++ii) {
                        iter.advance();
                        int label = iter.key();
                        TIntSet kp2Idxs = iter.value();
                        TIntIterator iter2 = kp2Idxs.iterator();
                        while (iter2.hasNext()) {
                            int kp2Idx = iter2.next();
                            PairInt p = new PairInt(
                                orb2.getKeyPoint1List().get(octave2).get(kp2Idx),
                                orb2.getKeyPoint0List().get(octave2).get(kp2Idx)
                            );
                            keypoints2LabelMap.put(p, label);
                        }
                    }
                }

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

                    //NOTE: not removing small sz2
                    if ((sz2 > sz1 && Math.abs(sz2 / sz1) > 1.2)) {
                        continue;
                    }

System.out.println("octave1=" + octave1 + " octave2=" + octave2 +
   " sz1=" + sz1 + " sz2=" + sz2 + " segIdx=" + segIdx +
     " nKP2=" + kp2Idxs.size()); 

                    //keypoints in this specific label2, segIdx
                    TIntIterator iter = kp2Idxs.iterator();
                    PairIntArray right = new PairIntArray(kp2Idxs.size());
                    while (iter.hasNext()) {
                        int kpIdx2 = iter.next();
                        int x = orb2.getKeyPoint1List().get(octave2).get(kpIdx2);
                        int y = orb2.getKeyPoint0List().get(octave2).get(kpIdx2);
                        right.add(x, y);
                    }
                    
                    MiscellaneousCurveHelper curveHelper 
                        = new MiscellaneousCurveHelper();
        
                    double[] xyCen2 = curveHelper.calculateXYCentroids(
                        labeledPoints2.get(segIdx));

                    ORB.Descriptors[] desc1 = getDescriptors(orb1, octave1);
                    ORB.Descriptors[] desc2 = getDescriptors(orb2, octave2);
                    int[][] costD = ORB.calcDescriptorCostMatrix(desc1, desc2);

                    // sorted greedy ordered by incr descr cost
                    PairIntArray m1 = new PairIntArray();
                    PairIntArray m2 = new PairIntArray();
                    TIntSet mLabels2 = new TIntHashSet();
                    
                    // using only the keypoints contained within the current
                    // labeled region, combinations of point pairs are used to
                    // find the best euclidean transformation.
                    // the transformation is applied to the bounds to pick
                    // up other points used later in the total cost.

                    // 0 = the salukwzde distance, that is the normalized tot cost
                    // 1 = sum of normalized keypoint descriptors
                    // 2 = sum of normalized keypoint distances from transformations
                    // 3 = the "fraction of whole" for transformed area matching
                    // 4 = number of keypoint matches
                    double[] normalizedCost = new double[5];
                    TransformationParameters params = matchGreedy(
                        segIdx, left, right, nBands, costD, nn1,
                        leftIdxMap, keypoints2IndexMap, keypoints2LabelMap,
                        extractKeypoints(orb2, octave2),
                        m1, m2, mLabels2,
                        orb1.getPyramidImages().get(0).a[0].length,
                        orb1.getPyramidImages().get(0).a.length,
                        distTol,
                        bounds1, keyBounds2Map,
                        point2LabelMap,
                        scale1, scale2, xyCen2, 
                        label2PairSizes, labeledPoints1, labeledPoints2,
                        label2AdjacencyMap,
                        normalizedCost);

//assert(normalizedCost[3] <= nkp1);

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
                    //TODO: may need to revise this limit
                    if (params == null || m1.getN() < 3) {
                        continue;
                    }

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
                    // 2 = sum of normalized keypoint distances from transformations
                    // 3 = the "fraction of whole" for transformed area matching
                    // 4 = number of keypoint matches
                    
                    correspondences.put(dataCount, corres);
                    transformations.put(dataCount, params);
                    descCosts.put(dataCount, normalizedCost[1]);
                    nDesc.put(dataCount, (int)normalizedCost[4]);
                    areaFraction.put(dataCount, (float)normalizedCost[3]);
                    distCosts.put(dataCount, normalizedCost[2]);
                    octs1.put(dataCount, octave1);
                    octs2.put(dataCount, octave2);
                    segIdxs.put(dataCount, segIdx);
                    
                    TIntList ml2 = new TIntArrayList(mLabels2);
                    ml2.sort();
                    matchedLabels2.put(dataCount, ml2);
                    
                    dataCount++;
                    
                }// end loop over octave2's segIdx
            }// end loop over octave2
        }  // end loop over octave1

        assert(correspondences.size() == dataCount);

        System.out.println("dataCount=" + dataCount);

        Transformer transformer = new Transformer();
              
        int img1Width = orb1.getPyramidImages().get(0).a[0].length;
        int img1Height = orb1.getPyramidImages().get(0).a.length;
        float sz1_0 = calculateObjectSize(bounds1);
        
        // key = octave1, octave2, label2, value = OneDIntArray sortedkeys2
        Map<TrioInt, OneDIntArray> sortedKeysMap = new
            HashMap<TrioInt, OneDIntArray>();
 
        TIntIntMap nBounds1 = new TIntIntHashMap();
        TIntObjectMap<PairIntArray> bounds2s = new TIntObjectHashMap<PairIntArray>();
        TIntObjectMap<PairIntArray> boundsMatched = new TIntObjectHashMap<PairIntArray>();
        TIntObjectMap<Double> bDistCost = new TIntObjectHashMap<Double>();
        for (int i = 0; i < dataCount; ++i) {
            int octave1 = octs1.get(i);
            int octave2 = octs2.get(i);
            int segIdx = segIdxs.get(i);
            float scale1 = scales1.get(octave1);

            // NOTE: all coordinates are in reference frame of the full
            // size octave images, that is octave1=0 and octave2=0
            
            // key = keypoint coords, value = index of keypoint within
            //     the orb keypoint list for this octave2
            TObjectIntMap<PairInt> keypoints2IndexMap
                = kp2IdxMapList.get(octave2);
            TIntObjectMap<TIntSet> labels2KPIdxs =
                labels2KPIdxsList.get(octave2);

            //key = keypoint coords, value = label that keypoint is within
            TObjectIntMap<PairInt> keypoints2LabelMap = new
                TObjectIntHashMap<PairInt>();
            {
                TIntObjectIterator<TIntSet> iter = labels2KPIdxs.iterator();
                for (int ii = 0; ii < labels2KPIdxs.size(); ++ii) {
                    iter.advance();
                    int label = iter.key();
                    TIntSet kp2Idxs = iter.value();
                    TIntIterator iter2 = kp2Idxs.iterator();
                    while (iter2.hasNext()) {
                        int kp2Idx = iter2.next();
                        PairInt p = new PairInt(
                            orb2.getKeyPoint1List().get(octave2).get(kp2Idx),
                            orb2.getKeyPoint0List().get(octave2).get(kp2Idx)
                        );
                        keypoints2LabelMap.put(p, label);
                    }
                }
            }

            TIntList sorted2 = matchedLabels2.get(i);
            // trimmed to combine to a contiguous region
            OneDIntArray sortedKeys2 = new OneDIntArray(
                sorted2.toArray(new int[sorted2.size()]));

            sortedKeysMap.put(new TrioInt(octave1, octave2, segIdx), 
                sortedKeys2);
            
            PairIntArray bounds2;
            if (keyBounds2Map.containsKey(sortedKeys2)) {
                bounds2 = keyBounds2Map.get(sortedKeys2);
            } else {
                Set<PairInt> combSet = new HashSet<PairInt>();
                for (int label2 : sortedKeys2.a) {
                    combSet.addAll(labeledPoints2.get(label2));
                }
                bounds2 = createOrderedBoundsSansSmoothing(combSet);
                keyBounds2Map.put(sortedKeys2, bounds2);
            }

            TObjectIntMap<PairInt> bounds2IndexMap
                = new TObjectIntHashMap<PairInt>();
            for (int ii = 0; ii < bounds2.getN(); ++ii) {
                PairInt p = new PairInt(bounds2.getX(ii), bounds2.getY(ii));
                bounds2IndexMap.put(p, ii);
            }
            
            int ne1 = bounds1.getN();
            int ne2 = bounds2.getN();

            PairIntArray bounds2Tr = transformer.applyTransformation(
                transformations.get(i), bounds2);

            /*            
            MatchedPointsTransformationCalculator tc = 
                new MatchedPointsTransformationCalculator();

            TransformationParameters revParams = 
                tc.swapReferenceFrames(transformations.get(i));

            PairIntArray bounds1RevTr = transformer.applyTransformation(
                revParams, bounds1);

            Set<PairInt> unique = new HashSet<PairInt>();
            for (int k = 0; k < bounds1RevTr.getN(); ++k) {
                unique.add(new PairInt(bounds1RevTr.getX(k),
                   bounds1RevTr.getY(k)));
            }
   
            ne1 = unique.size();
            */
            
            nBounds1.put(i, ne1);
            
            float sz2Tr = calculateObjectSize(bounds2Tr);
            if ((sz1_0 > sz2Tr && ((sz1_0/sz2Tr) > 1.5)) || 
                (sz2Tr > sz1_0 && ((sz2Tr/sz1_0) > 1.5))) {
                System.out.print("ERROR: scale difference too large:" +
                    " oct1=" + octave1 + " oct2=" + octave2 +
                    " segIdx=" + segIdx);
                if (sz1_0 > sz2Tr) {
                    System.out.println(" ratio=" + (sz1_0/sz2Tr));
                } else {
                    System.out.println(" ratio=" + (sz2Tr/sz1_0));
                }
                {// DEBUG: temporarily writing out error images
                    float scale2 = scales2.get(octave2);
                    Image img2 = ORB.convertToImage(
                        orb2.getPyramidImages().get(octave2));
                    for (int ii = 0; ii < bounds2.getN(); ++ii) {
                        int x = Math.round((float)bounds2.getX(ii)/scale2);
                        int y = Math.round((float)bounds2.getY(ii)/scale2);
                        ImageIOHelper.addPointToImage(x, y, img2, 1, 0, 255, 0);
                    }
                    MiscDebug.writeImage(img2, "_ERR_"
                        + octave1 + "_" + octave2 + "_" + segIdx + "_" +
                        MiscDebug.getCurrentTimeFormatted());
                }
                
                correspondences.get(i).clear();
                bounds2s.put(i, bounds2);
                boundsMatched.put(i, new PairIntArray());
                bDistCost.put(i, Double.MAX_VALUE);
                continue;
            }
            
            int distTol2 = Math.round((float)distTol/scale1);
            if (distTol2 < 1) {
                distTol2 = 1;
            }

            PairIntArray boundsMatchingIndexes = new PairIntArray();

            //sumsExtr = []{sumDist, count}
            double[] sumsExtr = sumKeypointDescAndDist2To1(
                bounds2, bounds2Tr, nnb1,
                bounds1IndexMap, bounds2IndexMap,
                img1Width, img1Height,
                distTol2, boundsMatchingIndexes);

            int nBMatched = (int)sumsExtr[1];

            bounds2s.put(i, bounds2);
            boundsMatched.put(i, boundsMatchingIndexes);
            bDistCost.put(i, sumsExtr[0]);

            for (int ii = 0; ii < boundsMatchingIndexes.getN(); ++ii) {
                int idx1 = boundsMatchingIndexes.getX(ii);
                int idx2 = boundsMatchingIndexes.getY(ii);
                int x1 = bounds1.getX(idx1);
                int y1 = bounds1.getY(idx1);
                int x2 = bounds2.getX(idx2);
                int y2 = bounds2.getY(idx2);
                QuadInt q = new QuadInt(x1, y1, x2, y2);
                correspondences.get(i).add(q);
            }
        }
        
        assert(dataCount == nBounds1.size());
        assert(dataCount == bounds2s.size());
        assert(dataCount == boundsMatched.size());
        assert(dataCount == bDistCost.size());
        assert(dataCount == transformations.size());
       
        float maxDesc = nBands * 256.0f;
        TIntList indexes = new TIntArrayList(dataCount);
        TFloatList costs = new TFloatArrayList(dataCount);
        for (int i = 0; i < dataCount; ++i) {

            int octave1 = octs1.get(i);
            int octave2 = octs2.get(i);

            float nKP1 = orb1.getKeyPoint0List().get(octave1).size();

            float nb1 = (float)bounds1.getN();

            float nd = nDesc.get(i);

            float descCost = descCosts.get(i).floatValue()/nd;
            float distCost = distCosts.get(i).floatValue()/nd;

            // calculate "fraction of whole" for keypoint descriptors
            //final float f1 = 1.f - (nd/nKP1);
            final float f1 = areaFraction.get(i);

            float d1 = descCost;
            float d2 = distCost;

            float nBMatched = boundsMatched.get(i).getN();
            float boundaryDistCost = bDistCost.get(i).floatValue()/nBMatched;

            // -- fraction of whole for boundary matching
            float f3 = 1.f - (nBMatched / (float)nBounds1.get(i));
            float d3 = boundaryDistCost;

            float tot = 2.f*f1 * f1 + d1 * d1 + d2 * d2 + d3 * d3 + f3 * f3;
            tot = f1 * f1 + d1 * d1 + d3 * d3 + f3 * f3;
            //tot = f1 * f1 + d1 * d1;

            indexes.add(i);
            costs.add(tot);

            String str1 = String.format(
                "octave1=%d octave2=%d segIdx=%d (rot=%d, s=%.2f) nCor=%d nDesc=%d",
                octave1, octave2, segIdxs.get(i),
                Math.round(transformations.get(i).getRotationInDegrees()),
                transformations.get(i).getScale(),
                correspondences.get(i).size(), (int)nd);

            String str2 = String.format(
                "i=%d descCost=%.2f f1=%.2f distCost=%.2f boundDistcost=%.2f f3=%.2f tot=%f",
                i, d1,        f1, d2,           d3, f3, tot);

            System.out.println(str1 + " " + str2);

            if (!Float.isInfinite(tot)) {

                {// DEBUG, print matched in green
                    List<QuadInt> qs = correspondences.get(i);
                    String str3 = Integer.toString(segIdxs.get(i));
                    while (str3.length() < 3) {
                        str3 = "0" + str3;
                    }
                    float scale1 = scales1.get(octave1);
                    Image img1 = ORB.convertToImage(
                        orb1.getPyramidImages().get(octave1));
                    for (int ii = 0; ii < qs.size(); ++ii) {
                        int x = Math.round((float)qs.get(ii).getA()/scale1);
                        int y = Math.round((float)qs.get(ii).getB()/scale1);
                        ImageIOHelper.addPointToImage(x, y, img1, 1, 0, 255, 0);
                    }
                    MiscDebug.writeImage(img1, "_TMP1_"
                        + octave1 + "_" + octave2 + "_" + str3 + "_" +
                        MiscDebug.getCurrentTimeFormatted());
                    //====
                    float scale2 = scales2.get(octave2);
                    img1 = ORB.convertToImage(
                        orb2.getPyramidImages().get(octave2));
                    for (int ii = 0; ii < qs.size(); ++ii) {
                        int x = Math.round((float)qs.get(ii).getC()/scale2);
                        int y = Math.round((float)qs.get(ii).getD()/scale2);
                        ImageIOHelper.addPointToImage(x, y, img1, 1, 0, 255, 0);
                    }
                    MiscDebug.writeImage(img1, "_TMP2_" +
                        octave1 + "_" + octave2 + "_" +
                        str3 + "_" + MiscDebug.getCurrentTimeFormatted());
                }
            }
        }

        assert(dataCount == nBounds1.size());
        assert(dataCount == bounds2s.size());
        assert(dataCount == boundsMatched.size());
        assert(dataCount == bDistCost.size());
        assert(dataCount == indexes.size());
        assert(dataCount == costs.size());
       
        QuickSort.sortBy1stArg(costs, indexes);

        //System.out.println("costs: " + Arrays.toString(costs));
        //System.out.println("indexes: " + Arrays.toString(costs));

        // trim to the top results within 1.12 * best
        int lastIdx = -1;

        float limit = 1.12f * costs.get(0);
        for (int j = 0; j < costs.size(); ++j) {
            if (
                //costs.get(j) > limit || 
                Float.isInfinite(costs.get(j))) {
                break;
            }
            lastIdx = j;
            System.out.format(
                "final results=%d %d segIdx=%d cost=%.2f (rot=%.2f)\n",
                costs.size(), j, segIdxs.get(indexes.get(j)),
                costs.get(j), 
                transformations.get(indexes.get(j)).getRotationInDegrees());
        }
        
        if (lastIdx < (dataCount - 1)) {
            // remove keys beyond lastIdx
            for (int i = (lastIdx + 1); i < dataCount; ++i) {
                int key = indexes.get(i);
                correspondences.remove(key);
                transformations.remove(key);
                descCosts.remove(key);
                nDesc.remove(key);
                distCosts.remove(key);
                areaFraction.remove(key);
                octs1.remove(key);
                octs2.remove(key);
                segIdxs.remove(key);
                nBounds1.remove(key);
                bounds2s.remove(key);
                boundsMatched.remove(key);
                bDistCost.remove(key);
            }
            indexes = indexes.subList(0, lastIdx + 1);
            costs = costs.subList(0, lastIdx + 1);
            dataCount = lastIdx;
            assert(indexes.size() == correspondences.size());
        }
        
        System.out.println("dataCount=" + dataCount);

        // a look at the chord differences
        
        double maxChordAvg = Double.MIN_VALUE;
        TIntObjectMap<Double> chordDiffAvgs = new TIntObjectHashMap<Double>();
        // could rewrite indexes here or just increase the new list
        // to same size as others and use "idx" to set a field
        for (int i = 0; i < indexes.size(); ++i) {
            int key = indexes.get(i);
            chordDiffAvgs.put(key, Double.MAX_VALUE);
        }
       
        TIntIntMap nb1s = new TIntIntHashMap();
        TIntIntMap nb2s = new TIntIntHashMap();
        TIntObjectMap<CorrespondenceList> results = new TIntObjectHashMap<CorrespondenceList>();
        for (int i = 0; i < indexes.size(); ++i) {

            if (Float.isInfinite(costs.get(i))) {
                break;
            }
            
            int idx = indexes.get(i);

            int segIdx = segIdxs.get(idx);
            int octave1 = octs1.get(idx);
            int octave2 = octs2.get(idx);
            float scale1 = orb1.getScalesList().get(octave1).get(0);
            float scale2 = orb2.getScalesList().get(octave2).get(0);
            
            List<QuadInt> qs = correspondences.get(idx);

            // the curves have to be nearly the same scale
            PairIntArray p = reduceBounds(bounds1, scale1);
            PairIntArray q = reduceBounds(bounds2s.get(idx), scale2);
            
            PairIntArray matchedIndexes = boundsMatched.get(idx);
            PairIntArray mIdxs = new PairIntArray(matchedIndexes.getN());
            for (int j = 0; j < matchedIndexes.getN(); ++j) {
                int idx1 = Math.round((float)matchedIndexes.getX(j)/scale1);
                int idx2 = Math.round((float)matchedIndexes.getY(j)/scale2);
                mIdxs.add(idx1, idx2);
            }
            
            //large voundaries needs larger dp for better runtime
            int dp = 1;
            if (p.getN() > 500 || q.getN() > 500) {
                int dn = Math.max(p.getN(), q.getN());
                dp += Math.ceil((float)dn/500.f);
            }
            
            PartialShapeMatcher matcher = new PartialShapeMatcher();
            matcher.overrideSamplingDistance(dp);
            TDoubleList chordDiffs = matcher.calculateChordDiffs(p, q, 
                mIdxs);
            
            TIntList rm = new TIntArrayList();
            double thresh = .3;
            double chordDiffSum = 0;
            for (int j = 0; j < chordDiffs.size(); ++j) {
                // filter out points > thresh
                double d = chordDiffs.get(j);
                if (d > thresh) {
                    rm.add(j);
                } else {
                    chordDiffSum += d;
                }
            }
            
            for (int j = (rm.size() - 1); j > -1; --j) {
                
                int rmIdx = rm.get(j);
                
                int idx1 = mIdxs.getX(rmIdx);
                int idx2 = mIdxs.getY(rmIdx);
                PairInt p1 = new PairInt(p.getX(idx1), p.getY(idx1));
                PairInt p2 = new PairInt(q.getX(idx2), q.getY(idx2));
                                
                //TODO: improve this w/ index map when have finished changes
                for (int k = 0; k < qs.size(); ++k) {
                    QuadInt q0 = qs.get(k);
                    PairInt p1c = new PairInt(q0.getA(), q0.getB());
                    PairInt p2c = new PairInt(q0.getC(), q0.getD());
                    if (p1.equals(p1c) && p2.equals(p2c)) {
                        qs.remove(k);
                        break;
                    }
                }                
                matchedIndexes.removeRange(rmIdx, rmIdx);
                chordDiffs.removeAt(rmIdx);
            }
            
            // should correct the correpondence distances too
            
            double chordAvg = chordDiffSum/(double)matchedIndexes.getN();
            
            if (chordAvg > maxChordAvg) {
                maxChordAvg = chordAvg;
            }
            
            if (matchedIndexes.getN() > 0) {
                chordDiffAvgs.put(idx, chordAvg);
            }
            
            // points are in full reference frame
            results.put(idx, new CorrespondenceList(qs));
            nb1s.put(idx, p.getN());
            nb2s.put(idx, q.getN());
        }
       
        // re-calculate the costs to include the shape component
               
        TFloatList resultCosts = new TFloatArrayList();
        TFloatList d5s = new TFloatArrayList();
        TIntList dataIndexes = new TIntArrayList();
        for (int i = 0; i < indexes.size(); ++i) {

            if (Float.isInfinite(costs.get(i))) {
                break;
            }
            
            int idx = indexes.get(i);
            
            int octave1 = octs1.get(idx);
            int octave2 = octs2.get(idx);
            
            PairIntArray matchedIndexes = boundsMatched.get(idx);
            int nBMatched = matchedIndexes.getN();
            
            float nb1 = nb1s.get(idx);
            
            float d5 = (float)(chordDiffAvgs.get(idx)/maxChordAvg);
            float f5 = 1.f - ((float)nBMatched/nb1);
            
            d5s.add(d5);
            
            float nKP1 = orb1.getKeyPoint0List().get(octave1).size();
            float nd = nDesc.get(idx);
            float descCost = descCosts.get(idx).floatValue()/nd;
            float distCost = distCosts.get(idx).floatValue()/nd;

            // calculate "fraction of whole" for keypoint descriptors
            //final float f1 = 1.f - (nd/nKP1);
            final float f1 = areaFraction.get(idx);
            float d1 = descCost;
            float d2 = distCost;

            float boundaryDistCost = bDistCost.get(idx).floatValue()/nBMatched;

            // -- fraction of whole for boundary matching
            float f3 = 1.f - (nBMatched / (float)nBounds1.get(idx));
            float d3 = boundaryDistCost;

            float tot = 2.f * f1 * f1 + d1 * d1 + d2 * d2 + d3 * d3 + f3 * f3;
            tot = f1 * f1 + d1 * d1 + d2 * d2
                + d3 * d3 + f3 * f3
                + d5 * d5 + f5 * f5;
            //tot = f1 * f1 + d1 * d1;

            if (!Float.isFinite(tot)) {
                continue;
            }
            
            resultCosts.add(tot);
            dataIndexes.add(idx);
            
            String str1 = String.format(
                "** oct1=%d oct2=%d segIdx=%d (rot=%d, s=%.2f) nCor=%d nDesc=%d nB=%d",
                octave1, octave2, segIdxs.get(idx),
                Math.round(transformations.get(idx).getRotationInDegrees()),
                transformations.get(idx).getScale(),
                correspondences.get(idx).size(), (int)nd, (int)nb1);

            String str2 = String.format(
                "\n  i=%d d1(desc)=%.2f f1=%.2f d2(dist)=%.2f d3(bDist)=%.2f f3=%.2f d5=%.2f f5=%.2f tot=%f",
                i, d1,        f1, d2,           d3, f3, 
                d5, f5, tot);

            System.out.println(str1 + " " + str2);
            
        }
        
        assert(resultCosts.size() == dataIndexes.size());
      
        QuickSort.sortBy1stThen2nd(resultCosts, d5s, dataIndexes);
        
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
                int count = 0;
                for (Entry<OneDIntArray, PairIntArray> entry :
                    keyBounds2Map.entrySet()) {

                    
                    int[] clr = ImageIOHelper.getNextRGB(count);

                    PairIntArray bounds2 = entry.getValue();
                    for (int idx = 0; idx < bounds2.getN(); ++idx) {
                        int x = Math.round((float) bounds2.getX(idx)/ scale2);
                        int y = Math.round((float) bounds2.getY(idx)/ scale2);
                        ImageIOHelper.addPointToImage(x, y, img1, 0,
                            clr[0], clr[1], clr[2]);
                    }
                    
                }
               MiscDebug.writeImage(img1, "_TMP4_" +count + "_"
                    + MiscDebug.getCurrentTimeFormatted());
                    count++;
            }
        }
        
        if (results.size() < 2) {
            
            List<CorrespondenceList> out = new ArrayList<CorrespondenceList>(
                results.size());
            
            for (int i = 0; i < dataIndexes.size(); ++i) {
                out.add(results.get(dataIndexes.get(i)));
            }
            
            return out;
        }
        
        float c0 = resultCosts.get(0);
        float c1 = resultCosts.get(1);
        float f = (c1 - c0)/c0;
        
        System.out.println("best cost=" + c0 + 
            " nCorr=" + results.get(dataIndexes.get(0)).getPoints1().size()
            + " segIdx=" + segIdxs.get(dataIndexes.get(0)));
        System.out.println("2nd best cost=" + resultCosts.get(1) + 
            " nCorr=" + results.get(dataIndexes.get(1)).getPoints1().size()
            + " (diff is " + f + " frac of best)");

        float limitFactor = 0.2f;
        
        if (f < limitFactor) {
            
            System.out.println("calculating patch diff sums");
            
            //TODO:  this patch based comparison may need to be
            //      edited to use epipolar projection and to match
            //      along lines between bounds (similar to steps in
            //      registration).
            //      other edits such as a texture filter may be needed
            
            TIntFloatMap patchSSDMap = calcCostsOfPatchSums(
               results, resultCosts,
               dataIndexes, orb1, orb2, octs1, octs2, 
               segIdxs, labeledPoints1, labeledPoints2,
               bounds1, nb1s, bounds2s,
               sortedKeysMap, transformations, 
               limitFactor);
           
            TFloatList costs3 = new TFloatArrayList();
            TFloatList patchSSDs = new TFloatArrayList();
            TIntList dataIndexes4 = new TIntArrayList();
            
            TIntFloatIterator iter3 = patchSSDMap.iterator();
            for (int i = 0; i < patchSSDMap.size(); ++i) {
                iter3.advance();
                int idx = iter3.key();
                float patchSSD = iter3.value();
                
                float nd = nDesc.get(idx);
                float descCost = descCosts.get(idx).floatValue()/nd;
                float d1 = descCost;
                
                float f1 = areaFraction.get(idx);
            
                PairIntArray matchedIndexes = boundsMatched.get(idx);
                int nBMatched = matchedIndexes.getN();
                float f3 = 1.f - (nBMatched / (float)nBounds1.get(idx));

                float d5 = (float)(chordDiffAvgs.get(idx)/maxChordAvg);
           
                float tot = d1*d1 + f1*f1 + f3*f3 + d5*d5 + patchSSD*patchSSD;
                
                costs3.add(tot);
                patchSSDs.add(patchSSD);
                dataIndexes4.add(idx);
            }
            
            QuickSort.sortBy1stThen2nd(costs3, patchSSDs, dataIndexes4);
           
            List<CorrespondenceList> output = new ArrayList<CorrespondenceList>(
                dataIndexes4.size());
            
            for (int i = 0; i < dataIndexes4.size(); ++i) {
                
                int idx = dataIndexes4.get(i);
                
                output.add(results.get(idx));    
            
                {//DEBUG
                    TransformationParameters params = 
                        transformations.get(idx);
                    
                    System.out.format(
                        "segIdx=%d (rot=%d, s=%.2f) tot=%.3f\n",
                        segIdxs.get(idx),
                        Math.round(params.getRotationInDegrees()),
                        params.getScale(), costs3.get(i));
                }
            }
            
            return output;
        }
        
        /*
        if (followWithShapeSearch) {

            results = aggregatedShapeMatch(orb1, orb2, nBands,
                labeledPoints1, labeledPoints2,
                octs1, octs2, segIdxs, scales1, scales2);

            System.out.println("nShapes=" + results.size());
        }*/

        List<CorrespondenceList> output = new ArrayList<CorrespondenceList>(
            results.size());
        for (int i = 0; i < dataIndexes.size(); ++i) {
            output.add(results.get(dataIndexes.get(i)));
        }
        
        return output;
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
    public List<CorrespondenceList> matchSmall(ORB orb1, ORB orb2,
        Set<PairInt> labeledPoints1, List<Set<PairInt>> labeledPoints2) {

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
                    matcher._overrideToThreshhold(0.2f);
                    matcher.setToRemoveOutliers();
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

  //TODO: revisit and edit for changes present in match0

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
        // the differences in chords for true match shapes should not change
        // with scale (excepting due to resolution affecting the shapes),
        // so the maximum difference in the avg chord can be
        // used on all of the results when re-normalizing and re-calculating dist.
        // The number of matches possible for each scale does change with scale
        // and that is handled in the "fraction of whole term" already.
        // so a "re-normalization" to calculate new salukvazde costs just needs
        // to use the maximum avg chord sum diff for all.
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
                        double costNorm = minC / maxCost;
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
                double costNorm = minC / maxDesc;
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
                double costNorm = minC / maxDesc;
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
    private static double[] sumKeypointDescAndDist2To1(int segIdx,
        int[][] cost, int nBands,
        PairIntArray a2,
        PairIntArray a2TrTo1, TransformationParameters params,
        NearestNeighbor2D nn1,
        TObjectIntMap<PairInt> p1KPIndexMap,
        TObjectIntMap<PairInt> p2KPIndexMap, 
        TIntSet skip2,
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
            
            if (skip2.contains(k)) {
                continue;
            }
            
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
                double costNorm = minC / maxDesc;
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
     * @param m12 output matched indexes
     * @return
     */
    private static double[] sumKeypointDescAndDist2To1(
        PairIntArray a2,
        PairIntArray a2TrTo1,
        NearestNeighbor2D nn1,
        TObjectIntMap<PairInt> p1KPIndexMap,
        TObjectIntMap<PairInt> p2KPIndexMap,
        int img1Width, int img1Height, int distTol,
        PairIntArray m12) {

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
                m12.add(minIdx1, kpIdx2);
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
                int idx1 = m12.getX(idx);
                int idx2 = m12.getY(idx);
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
            m12.removeRange(0, m12.getN() - 1);
            for (int i = 0; i < count; ++i) {
                m12.add(matched.getX(i), matched.getY(i));
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

    private PairIntArray createOrderedBoundsSansSmoothing(
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
                    calculateObjectSize(bounds) + " segIdx=" + segIdx);
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
     * @param out1
     * @param out2
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
    private TransformationParameters matchGreedy(
        int segIdx, PairIntArray left, PairIntArray right,
        int nBands, int[][] costD, NearestNeighbor2D nn1,
        TObjectIntMap<PairInt> keypoints1IndexMap,
        TObjectIntMap<PairInt> keypoints2IndexMap,
        TObjectIntMap<PairInt> keypoints2LabelMap,
        PairIntArray keypoints2,
        PairIntArray out1, PairIntArray out2, TIntSet outLabels2,
        int img1Width, int img1Height,
        int distTol,
        PairIntArray bounds1,
        Map<OneDIntArray, PairIntArray> keyBounds2Map,
        TObjectIntMap<PairInt> points2LabelMap,
        float scale1, float scale2, double[] xyCenlabeled2,
        TObjectIntMap<PairInt> label2PairSizes,
        Set<PairInt> label1Set,
        List<Set<PairInt>> label2Sets,
        TIntObjectMap<VeryLongBitString> label2AdjacencyMap,
        double[] outputNormalizedCost) {

        // NOTE that all coordinates are in the full reference frame
        //   and remain that way through this method

        float expectedScale = scale1/scale2;

        float maxDesc = nBands * 256.0f;

        int n1 = left.getN();
        int n2 = right.getN();
        
        int sz1 = calculateObjectSize(left);

        // factors to correct pixel counts for image scales
        float area1 = 1.f/(scale1 * scale1);
        float area2 = 1.f/(scale2 * scale2);
        float set1Area = (float)label1Set.size() * area1 * area2;
         
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

        Set<PairInt> leftSet = Misc.convert(left);
        
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        double[] xyCen1 = curveHelper.calculateXYCentroids(leftSet);
        
        //TODO: this may be revised after more tests
        float costLimit = (nBands * 256) * 0.5f;

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
                if (c > costLimit) {
                    continue;
                }
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

        MatchedPointsTransformationCalculator tc
            = new MatchedPointsTransformationCalculator();
        Transformer transformer = new Transformer();
       
        double bestCost = Double.MAX_VALUE;
        double bestCost2 = -1;
        int[] bestMIdx1s = null;
        int[] bestMIdx2s = null;
        TIntSet bestLabels2 = null;
        int bestN = -1;
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

                assert(!(leftX1 != leftX2 && leftY1 != leftY1));
                assert(!(rightX1 != rightX2 && rightY1 != rightY1));

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

                // distTol is w.r.t. scale1=1, so if scale1 > 1, need to adjust
                // NOTE that distTol is applied to pixels in dataset1 reference frame.
                int distTol2 = Math.round((float)distTol/scale1);
                if (distTol2 < 1) {
                    distTol2 = 1;
                }
                int szSq = (sz1/2) + distTol2;
                szSq *= szSq;

                // points in set1 that match a transformed set2 point
                Set<PairInt> set1Intersection = new HashSet<PairInt>();
                
                // ----- transform all keypoints2 and match
                PairIntArray keypoints2Tr = transformer.applyTransformation(
                    params, keypoints2);
                
                // --- while adding key points that are the nearest matching  
                //     to transformed points within tolerance,
                //     need to filter out points that are further from
                //     a object distance from the center of 
                //     the current index for the fast first size filter.
                //     -- because the transformation is now known, more than
                //        the maximum dimension can be used.
                //        -- after keypoints of adjacent labeled regions are
                //           added, can look at the fraction of the region
                //           that intersects with the region for set1 (the template
                //           object) and remove those regions with small fractions
                //        -- NOTE that in the process, there is now information
                //           to make a total area fraction component for the
                //           cost
                
                // filter out keypoints further than sz1/2 or so
                // from the center of segIdx region.
                // to keep indexes correcting, making a skip set instead
                TIntSet labels2Contiguous = new TIntHashSet();
                TIntSet skip = new TIntHashSet();
                for (int jj = 0; jj < keypoints2Tr.getN(); ++jj) {
                    double diffX = keypoints2Tr.getX(jj) - xyCen1[0];
                    double diffY = keypoints2Tr.getY(jj) - xyCen1[1];
                    double d = diffX * diffX + diffY * diffY;
                    if (d > szSq) {
                        skip.add(jj);
                    } else {
                        PairInt p3 = new PairInt(keypoints2.getX(jj),
                            keypoints2.getY(jj));
                        if (!points2LabelMap.containsKey(p3)) {
                            skip.add(jj);
                        } else {
                            int label2 = points2LabelMap.get(p3);
                            labels2Contiguous.add(label2);
                        }
                    }
                }

                labels2Contiguous = new TIntHashSet(
                    extractContiguousFromLabels(
                    segIdx, label2AdjacencyMap, labels2Contiguous));
                labels2Contiguous.add(segIdx);
                
                TIntSet rmLabel2 = new TIntHashSet();
                TIntSet visited = new TIntHashSet();
                TIntIterator iter = labels2Contiguous.iterator();
                while (iter.hasNext()) {
                    int label2 = iter.next();
                    if (visited.contains(label2)) {
                        continue;
                    }
                    visited.add(label2);
                    
                    Set<PairInt> set2 = label2Sets.get(label2);
                    Set<PairInt> set2Tr = transformer.applyTransformation2(
                        params, set2);
                    
                    Set<PairInt> uniqueSet2Tr = new HashSet<PairInt>();
                    int nIn = 0;
                    for (PairInt p3 : set2Tr) {
                        if (uniqueSet2Tr.contains(p3)) {
                            continue;
                        }
                        uniqueSet2Tr.add(p3);
                        if (label1Set.contains(p3)) {
                            nIn++;
                        }
                    }
                   
                    float f = (float)nIn/(float)uniqueSet2Tr.size();
         if (segIdx == 3) {
         //    System.out.println("s1=" + scale1 + " s2=" + scale2 +
         //       " rot=" + Math.round(params.getRotationInDegrees()) 
         //       + " lbl2=" + label2 + " f=" + f 
         //       + " n=" + uniqueSet2Tr.size());
         }
                    if (f < 0.7) {
                        rmLabel2.add(label2);
                    } else {
                        // add the intersection to the total
                        for (PairInt p3 : set2Tr) {
                            if (label1Set.contains(p3)) {
                                set1Intersection.add(p3);
                            }
                        }
                    }
                }
                labels2Contiguous.removeAll(rmLabel2);
                
                // add keypoints w/ labels in rmLabels2 to skip
                for (int jj = 0; jj < keypoints2Tr.getN(); ++jj) {
                    if (skip.contains(jj)) {
                        continue;
                    }
                    PairInt p3 = new PairInt(keypoints2.getX(jj),
                        keypoints2.getY(jj));
                    int label2 = points2LabelMap.get(p3);
                    if (!labels2Contiguous.contains(label2)) {
                        skip.add(jj);
                    }
                }
                
                // transform the points of set1 which are not in the intersection
                Set<PairInt> missing1 = new HashSet<PairInt>(label1Set);
                missing1.removeAll(set1Intersection);
                
                TransformationParameters revParams
                    = tc.swapReferenceFrames(params);

                Set<PairInt> missing1Tr = transformer.applyTransformation2(
                    revParams, missing1);
                
                TIntObjectMap<Set<PairInt>> missing1Label2Map = 
                    new TIntObjectHashMap<Set<PairInt>>();
                for (PairInt pTr1 : missing1Tr) {
                    if (points2LabelMap.containsKey(pTr1)) {
                       int label2 = points2LabelMap.get(pTr1);
                       if (labels2Contiguous.contains(label2)) {
                           continue;
                       }
                       Set<PairInt> set2 = missing1Label2Map.get(label2);
                       if (set2 == null) {
                           set2 = new HashSet<PairInt>();
                           missing1Label2Map.put(label2, set2);
                       }
                       set2.add(pTr1);
                    }
                }
                
    if (segIdx == 3) {
    //    System.out.println(
    //       "(rot=" + Math.round(params.getRotationInDegrees())
    //       + ") " + " labels2 BEFORE set1 rev trans="
    //       + labels2Contiguous.toString());
    }
                
                int nAddedIntersection = 0;
                TIntObjectIterator<Set<PairInt>> iter2 = missing1Label2Map.iterator();
                for (int k = 0; k < missing1Label2Map.size(); ++k) {
                    iter2.advance();
                    int label2 = iter2.key();
                    Set<PairInt> inter2 = iter2.value();
                
                    float n0 = (float)label2Sets.get(label2).size();
                    float n = n0 * area2;
                    float f = (float)inter2.size()/n;
                                        
                    if (f >= 0.8) {
                        //NOTE: filter for contiguous again afterwards
                        labels2Contiguous.add(label2);
                        nAddedIntersection += inter2.size();
  if (segIdx == 3) {
  //     System.out.println(
  //         "(rot=" + Math.round(params.getRotationInDegrees())
  //         + ") " + " f=" + f + " l=" + label2);
   }
                    }
                }
               
                labels2Contiguous = new TIntHashSet(
                    extractContiguousFromLabels(
                    segIdx, label2AdjacencyMap, labels2Contiguous));
                
   if (segIdx == 3) {
   //    System.out.println(
   //        "(rot=" + Math.round(params.getRotationInDegrees())
   //        + ") " + " labels2 after set1 rev trans="
   //        + labels2Contiguous.toString());
   }   
   
                //TODO: for regions with a critical areaFraction,
                //      or for all,
                //      will "fill in" the labeled regions that are consistent
                //      with the reverse transformation of missing set1 points,
                //      and if those labeled regions are nearly completely
                //      found by the missing intersection,
                //      will add those labeled regions.
                //      note that those regions are missing because they
                //      don't have keypoints.
               
                int nM1 = set1Intersection.size() + nAddedIntersection;
                float areaFraction = (float)nM1/set1Area;
              
                // filled with indexes w.r.t left and right arrays
                int[] mIdx1s = new int[keypoints2LabelMap.size()];
                int[] mIdx2s = new int[keypoints2LabelMap.size()];

                //sums is []{sumDesc, sumDist, count}
                // where sumDesc is the sum of normalized descriptors
                //       each with a value of 0 to 1 and summed over
                //       count number of them. (so max is <= count)
                // same for sumDist
                double[] sums = sumKeypointDescAndDist2To1(
                    segIdx, costD, nBands,
                    keypoints2, keypoints2Tr, params, nn1,
                    keypoints1IndexMap, keypoints2IndexMap, 
                    skip,
                    distTol2, mIdx1s, mIdx2s);

                final int nMatched = (int)sums[2];
 
    if (segIdx == 3) {
    //    System.out.println("kp.n=" + keypoints2.getN() + 
    //        " skip.n=" + skip.size() 
    //        + " nAreaM=" + nM1 + " areaFraction=" 
    //        + areaFraction + " nKPMatched=" + nMatched);
    }
    
                if (nMatched < 3) {
                    continue;
                }

                double sumDescr = sums[0];
                double sumDist = sums[1];

                float descCost = (float)sumDescr/(float)nMatched;
                float distCost = (float)sumDist/(float)nMatched;
                float d1 = descCost;
                float d2 = distCost;
                //final float f1 = 1.f - ((float)nMatched/(float)n1);
                final float f1 = 1.f - areaFraction;
                
                //TODO: consider whether d2 should be included in tot
                float tot = d1*d1 + d2*d2 + f1 * f1;
                tot = d1 * d1 + f1 * f1;
                
if (segIdx==3) {
String str1 = String.format(
    "?? lbl=%d trRot=%d scl=%.2f nm=%d",
    segIdx, (int) params.getRotationInDegrees(),
    params.getScale(), nMatched);
String str2 = String.format(
    " dc=%.2f dt=%.2f f1=%.2f tot=%.2f",
    d1, d2, f1, tot);
System.out.println(str1 + str2);
}
                if (tot < bestCost) {

                    String str1 = String.format(
                        "< lbl=%d trRot=%d scl=%.2f nm=%d",
                        segIdx, (int)params.getRotationInDegrees(),
                        params.getScale(), nMatched);
                    String str2 = String.format(
                        " dc=%.2f dt=%.2f f1=%.2f tot=%.2f",
                        d1, d2, f1, tot);
                    System.out.println(str1 + str2);

                    bestCost = tot;

                    // 0 = the salukwzde distance, that is the normalized tot cost
                    // 1 = sum of normalized keypoint descriptors
                    // 2 = sum of normalized keypoint distances from transformations
                    // 3 = the "fraction of whole" for transformed area matching
                    // 4 = number of keypoint matches
                    outputNormalizedCost[0] = bestCost;
                    outputNormalizedCost[1] = sumDescr;
                    outputNormalizedCost[2] = sumDist;
                    outputNormalizedCost[3] = f1;
                    outputNormalizedCost[4] = nMatched;

                    bestMIdx1s = mIdx1s;
                    bestMIdx2s = mIdx2s;
                    bestN = nMatched;
                    bestLabels2 = labels2Contiguous;
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

                PairInt pRight = new PairInt(keypoints2.getX(kpIdx2), keypoints2.getY(kpIdx2));
                if (addedLeft.contains(pLeft) || addedRight.contains(pRight)) {
                    continue;
                }
                addedLeft.add(pLeft);
                addedRight.add(pRight);
                out1.add(pLeft.getX(), pLeft.getY());
                out2.add(pRight.getX(), pRight.getY());
            }
            outLabels2.clear();
            outLabels2.addAll(bestLabels2);
        } else {
            Arrays.fill(outputNormalizedCost, Double.MAX_VALUE);
            outputNormalizedCost[outputNormalizedCost.length - 1] = 0;
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

        //TODO:  could improve this by reducing the unmatched2 points
        //  to only those within sz1 transformed size of center of
        //  the labeled region.

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

            float costNorm = c / maxDesc;
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

    private List<CorrespondenceList> aggregatedShapeMatch(
        ORB orb1, ORB orb2, int nBands,
        Set<PairInt> labeledPoints1, List<Set<PairInt>> labeledPoints2,
        TIntList octs1, TIntList octs2, TIntList segIdxs,
        TFloatList scales1, TFloatList scales2) {

        System.out.println("aggregatedShapeMatch for " + octs1.size() + " results");

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

        // store each ShapeFinder2.ShapeFinderResult result
        // need to compare results from different octaves
        //    diff chord sum
        //    number matched
        //    max matchable for octave pair (==n1)
        //    max avg chord diff for octave pair

        // the differences in chords for true match shapes should not change
        // with scale (excepting due to resolution affecting the shapes),
        // so the maximum difference in the avg chord can be
        // used on all of the results when re-normalizing and re-calculating dist.
        // The number of matches possible for each scale does change with scale
        // and that is handled in the "fraction of whole term" already.
        // so a "re-normalization" to calculate new salukvazde costs just needs
        // to use the maximum avg chord sum diff for all.

        TIntList shapeResultIdxs = new TIntArrayList();
        List<ShapeFinder2.ShapeFinderResult> shapeResults = new
            ArrayList<ShapeFinder2.ShapeFinderResult>();

        double maxChordDiff = Double.MIN_VALUE;

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

            float factor = 2.f;
            // removing any set that, when combined with set for segIdx,
            // has dimensions larger than twice the size of bounds1
            removeSetsByCombinedSize(bounds1, segIdx, mapOfSets2, factor);

            removeSetsNotConnected(segIdx, mapOfSets2);

            Map<OneDIntArray, PairIntArray> keyBoundsMap2 =
                keybounds2Maps.get(octave2);

            int xMax1 = orb1.getPyramidImages().get(octave1).a[0].length -1;
            int yMax1 = orb1.getPyramidImages().get(octave1).a.length - 1;

            int xMax2 = orb2.getPyramidImages().get(octave2).a[0].length -1;
            int yMax2 = orb2.getPyramidImages().get(octave2).a.length - 1;

            System.out.println("start shape search for octave1=" + octave1 +
                " octave2=" + octave2 + " segIdx=" + segIdx);

            ShapeFinder2 shapeFinder = new ShapeFinder2(bounds1, scale1, sz1,
                xMax1, yMax1, mapOfSets2, scale2,
                keyBoundsMap2, xMax2, yMax2);

            ShapeFinder2.ShapeFinderResult result = shapeFinder.findAggregated(segIdx);

            if (result == null) {
                continue;
            }

            System.out.println("shapeFinder result for segIdx=" +
                + segIdx + " " + result.toString());

            shapeResultIdxs.add(i);
            shapeResults.add(result);

            double cda = result.getChordDiffSum()/(double)result.getNumberOfMatches();
            if (cda > maxChordDiff) {
                maxChordDiff = cda;
            }

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
                        + "_" + MiscDebug.getCurrentTimeFormatted();

                    String filePath = plotter.writeImage("_shape_" + str);

                } catch (Throwable t) {
                    System.err.println(t.getMessage());
                }
            }
        }

        float maxDesc = nBands * 256.0f;

        // sort the results by cost
        TFloatList shapeCosts = new TFloatArrayList(shapeResultIdxs.size());
        TIntList shapeCostIdxs = new TIntArrayList();

        for (int i = 0; i < shapeResultIdxs.size(); ++i) {

            int idx = shapeResultIdxs.get(i);

            ShapeFinder2.ShapeFinderResult shapeResult =
                shapeResults.get(i);

            int octave1 = octs1.get(idx);
            int octave2 = octs2.get(idx);
            //int segIdx = segIdxs.get(idx);
            //float scale1 = scales1.get(octave1);
            //float scale2 = scales2.get(octave2);

            float nb1 = shapeResult.bounds1.getN();

            float nMatched = shapeResult.getNumberOfMatches();

            double chordDiff = shapeResult.getChordDiffSum()/nMatched;
            double d1 = (chordDiff / maxDesc);

            double f1 = 1.f - (nMatched/nb1);

            double sd = d1 * d1 + f1 * f1;

            shapeCosts.add((float)sd);
            shapeCostIdxs.add(i);
        }

        QuickSort.sortBy1stArg(shapeCosts, shapeCostIdxs);

        List<CorrespondenceList> results = new ArrayList<CorrespondenceList>();
        for (int i = 0; i < shapeCostIdxs.size(); ++i) {

            int idx0 = shapeCostIdxs.get(i);

            int idx = shapeResultIdxs.get(idx0);

            ShapeFinder2.ShapeFinderResult shapeResult =
                shapeResults.get(idx0);

            int np = shapeResult.getNumberOfMatches();

            int octave1 = octs1.get(idx);
            int octave2 = octs2.get(idx);
            //int segIdx = segIdxs.get(idx);
            float scale1 = scales1.get(octave1);
            float scale2 = scales2.get(octave2);

            PairIntArray bounds1 = shapeResult.bounds1;
            PairIntArray bounds2 = shapeResult.bounds2;

            // shapeResult points need to be scaled back up to full reference
            // frame
            List<QuadInt> qs = new ArrayList<QuadInt>();
            for (int j = 0; j < np; ++j) {
                int idx1 = shapeResult.idx1s.get(j);
                int idx2 = shapeResult.idx2s.get(j);
                int x1 = Math.round((float)bounds1.getX(idx1) * scale1);
                int y1 = Math.round((float)bounds1.getY(idx1) * scale1);
                int x2 = Math.round((float)bounds2.getX(idx2) * scale2);
                int y2 = Math.round((float)bounds2.getY(idx2) * scale2);
                qs.add(new QuadInt(x1, y1, x2, y2));
            }

            // points are in full reference frame
            results.add(new CorrespondenceList(qs));
        }

        return results;
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
        int segIdx2, TIntObjectMap<Set<PairInt>> mapOfBounds2, float factor) {

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

            if (sz2 > factor * sz1) {
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

    private int[] extractAndSortLabels(int segIdx,
        List<QuadInt> matched12,
        TObjectIntMap<PairInt> keypoints2LabelMap,
        TIntObjectMap<VeryLongBitString> label2AdjacencyMap) {

        TIntSet labelSet = new TIntHashSet();
        for(int i = 0; i < matched12.size(); ++i) {
            QuadInt qs = matched12.get(i);
            PairInt p2 = new PairInt(qs.getC(), qs.getD());

            assert(keypoints2LabelMap.containsKey(p2));
            int label2 = keypoints2LabelMap.get(p2);
            labelSet.add(label2);
        }

        TIntList labelList = extractContiguousFromLabels(segIdx,
            label2AdjacencyMap, labelSet);

        labelList.sort();
        
        int[] labels = labelList.toArray(new int[labelList.size()]);

        return labels;
    }

    private TIntList extractContiguousFromLabels(int segIdx,
        TIntObjectMap<VeryLongBitString> label2AdjacencyMap,
        TIntSet labelSet) {

        //NOTE: this may change: for now, starting with src=segIdx, making a
        //   DFS traversal of neighbors within labelSet to keep labelList
        //   contiguous.
        //   Might consider adding the intersection of members in labelSet
        //   to include labels not in labelSet that connect a label which would
        //   otherwise be excluded from labelList...caveat is not to include
        //   external labels that would result in the wrong boundary.

        TIntSet added = new TIntHashSet();
        TIntSet visited = new TIntHashSet();

        TIntList labelList = new TIntArrayList();
        labelList.add(segIdx);

        Stack<Integer> stack = new Stack<Integer>();
        stack.add(segIdx);
        while (!stack.isEmpty()) {
            int label = stack.pop().intValue();
            if (visited.contains(label)) {
                continue;
            }

            int[] nbrLabels = label2AdjacencyMap.get(label).getSetBits();
            for (int label2 : nbrLabels) {
                if (!added.contains(label2) && labelSet.contains(label2)) {
                    labelList.add(label2);
                    stack.add(Integer.valueOf(label2));
                    added.add(label2);
                }
            }
            visited.add(label);
        }

        return labelList;
    }
    
    private PairIntArray extractKeypoints(ORB orb, int octave) {

        int n = orb.getKeyPoint0List().get(octave).size();
        PairIntArray out = new PairIntArray(n);
        for (int i = 0; i < n; ++i) {
            out.add(orb.getKeyPoint1List().get(octave).get(i),
                orb.getKeyPoint0List().get(octave).get(i));
        }

        return out;
    }

    private TIntObjectMap<VeryLongBitString> createAdjacencyMap(
        TObjectIntMap<PairInt> pointIndexMap, List<Set<PairInt>> labeledPoints) {

        int n = labeledPoints.size();

        TIntObjectMap<VeryLongBitString> output
            = new TIntObjectHashMap<VeryLongBitString>();

        int[] dxs = Misc.dx4;
        int[] dys = Misc.dy4;

        for (int label = 0; label < n; ++label) {
            Set<PairInt> set = labeledPoints.get(label);
            VeryLongBitString nbrs = new VeryLongBitString(n);
            output.put(label, nbrs);

            for (PairInt p : set) {
                int x = p.getX();
                int y = p.getY();
                for (int k = 0; k < dxs.length; ++k) {
                    int x2 = x + dxs[k];
                    int y2 = y + dys[k];
                    PairInt p2 = new PairInt(x2, y2);
                    if (pointIndexMap.containsKey(p2)) {
                        int label2 = pointIndexMap.get(p2);
                        if (label2 != label) {
                            nbrs.setBit(label2);
                        }
                    }
                }
            }
        }

        return output;
    }

    /**
     * calculate patch SSD costs, add to some subset of previous cost terms,
     * and trim and re-sort the given results.
     * 
     * @param results
     * @param resultCosts
     * @param dataIndexes
     * @param orb1
     * @param orb2
     * @param octs1
     * @param octs2
     * @param segIdxs
     * @param labeledPoints1
     * @param labeledPoints2
     * @param bounds1
     * @param bounds2s
     * @param sortedKeysMap
     * @param transformations
     * @param limitFactor
     * @return 
     */
    private TIntFloatMap calcCostsOfPatchSums(
        TIntObjectMap<CorrespondenceList> results, 
        TFloatList resultCosts, 
        TIntList sortedDataIndexes, 
        ORB orb1, ORB orb2, 
        TIntIntMap octs1, TIntIntMap octs2, 
        TIntIntMap segIdxs, 
        Set<PairInt> labeledPoints1, List<Set<PairInt>> labeledPoints2, 
        PairIntArray bounds1, TIntIntMap nb1s,
        TIntObjectMap<PairIntArray> bounds2s, 
        Map<TrioInt, OneDIntArray> sortedKeysMap, 
        TIntObjectMap<TransformationParameters> transformations, 
        float limitFactor) {

        // coordinates are in the reference frame of the full size octave 0 images
        // row major format:
        TwoDFloatArray imgH1 = orb1.getPyramidImagesH().get(0);
        TwoDFloatArray imgS1 = orb1.getPyramidImagesS().get(0);
        TwoDFloatArray imgV1 = orb1.getPyramidImagesV().get(0);

        // row major format:
        TwoDFloatArray imgH2 = orb2.getPyramidImagesH().get(0);
        TwoDFloatArray imgS2 = orb2.getPyramidImagesS().get(0);
        TwoDFloatArray imgV2 = orb2.getPyramidImagesV().get(0);

        Transformer transformer = new Transformer();
        
        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;
        
        float limitCost = (1.f + limitFactor) * resultCosts.get(0);
        
        TIntFloatMap ssdMap = new TIntFloatHashMap();
        
        for (int i = 0; i < sortedDataIndexes.size(); ++i) {
            if (resultCosts.get(i) > limitCost) {
                break;
            }
            int idx = sortedDataIndexes.get(i);
            
            int octave1 = octs1.get(idx);
            int octave2 = octs2.get(idx);
            int segIdx = segIdxs.get(idx);
            float scale1 = orb1.getScalesList().get(octave1).get(0);
            float scale2 = orb2.getScalesList().get(octave2).get(0);
   
            // factors to correct pixel counts for image scales
            //float area1 = 1.f / (scale1 * scale1);
            //float area2 = 1.f / (scale2 * scale2);
            //float set1Area = (float) labeledPoints1.size() * area1 * area2;
            
            OneDIntArray sortedKeys = sortedKeysMap.get(
                new TrioInt(octave1, octave2, segIdx));

            TransformationParameters params = transformations.get(idx);
            
            PairIntArray bounds2 = bounds2s.get(idx);
            
            Set<PairInt> combinedSet2 = new HashSet<PairInt>();
            for (int label2 : sortedKeys.a) {
                combinedSet2.addAll(labeledPoints2.get(label2));
            }
            
            // NOTE, because of possibility of projection and the
            //    resulting foreshortening along an axis, may need to consider
            //    a more time consuming epipolar solution, then
            //    transformation and nearest neighbor consistent
            //    with bounds along epipolar lines (similar to methods
            //    used in registration).
            // NOTE also: whichever method is used, it would be faster
            //   to transform the entire combinedSet once
            //   instead of several times for overlapping descrptor
            //   points below.
            
            // transform each point in label2 and sum the 8 neighbors
            // -- compare to the same within bounds in other.
            // note: nIn12, nIn2
            int n12 = 0;
            double patchSums = 0;

            float[] h1s = new float[8];
            float[] s1s = new float[8];
            float[] v1s = new float[8];
            float[] h2s = new float[8];
            float[] s2s = new float[8];
            float[] v2s = new float[8];
            
            Set<PairInt> visited = new HashSet<PairInt>();
            
            for (PairInt p2 : combinedSet2) {
                
                // transform the 8 neighbor region into reference frame
                // of image 1 and store for normalization
                
                int count = 0;
                double sum = 0;
                
                // to speed this up, could use the rotation matrix
                //   use in features package
                
                for (int k = 0; k < dxs.length; ++k) {
                    int x3 = p2.getX() + dxs[k];
                    int y3 = p2.getY() + dys[k];
                    PairInt p3 = new PairInt(x3, y3);
                    if (!combinedSet2.contains(p3)) {
                        continue;
                    }
                    
                    double[] tr = transformer.applyTransformation(params,
                        x3, y3);

                    int xTr = (int) Math.round(tr[0]);
                    int yTr = (int) Math.round(tr[1]);

                    PairInt pTr = new PairInt(xTr, yTr);
                    if (!labeledPoints1.contains(pTr) || visited.contains(pTr)) {
                        continue;
                    }
                    
                    visited.add(pTr);
                    
                    h1s[count] = imgH1.a[yTr][xTr];
                    h2s[count] = imgH2.a[y3][x3];
                    s1s[count] = imgS1.a[yTr][xTr];
                    s2s[count] = imgS2.a[y3][x3];
                    v1s[count] = imgV1.a[yTr][xTr];
                    v2s[count] = imgV2.a[y3][x3];
                    
                    count++;
                }
                
                if (count == 0) {
                    continue;
                }
                                
                // calc mean and substr it for the brighness arrays?
                calcMeanAndSubtr(v1s, count);
                calcMeanAndSubtr(v2s, count);
                
                double hSum = 0;
                double sSum = 0;
                double vSum = 0;
                for (int j = 0; j < count; ++j) {
                    float diffH = h1s[j] - h2s[j];
                    float diffS = s1s[j] - s2s[j];
                    float diffV = v1s[j] - v2s[j];
                    hSum += (diffH * diffH);
                    sSum += (diffS * diffS);
                    vSum += (diffV * diffV);
                }
                hSum /= (double)count;
                sSum /= (double)count;
                vSum /= (double)count;
            
                patchSums += ((hSum + sSum + vSum)/3.);
            
                n12++;
            }
            
            // each h,s,v image is scaled to values from 0 to 1.
            // so max difference possible is "1".
            // then the normalization is count * 1
            
            patchSums /= (double)n12;
            
            float d6 = (float)patchSums;
            
            ssdMap.put(idx, d6);
        }
        
        return ssdMap;
    }

    private void calcMeanAndSubtr(float[] a, int count) {
        
        double sum = 0;
        for (int i = 0; i < count; ++i) {
            sum += a[i];
        }
        
        float avg = (float)sum/(float)count;
        
        for (int i = 0; i < count; ++i) {
            a[i] -= avg;
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
    
    private static int[] numberInSet(Set<PairInt> shape1, 
        Set<PairInt> set2, TransformationParameters params) {
        
        Transformer transformer = new Transformer();
        
        Set<PairInt> set2Tr = transformer.applyTransformation2(params, set2);
        
        int nIn = 0;
        for (PairInt p : set2Tr) {
            if (shape1.contains(p)) {
                nIn++;
            }
        }
        
        return new int[]{nIn, set2Tr.size() - nIn};
    }

}
