package algorithms.imageProcessing.matching;

import algorithms.QuickSort;
import algorithms.compGeometry.FurthestPair;
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
    
    /**
     * match template image and shape in orb1 and labeledPoints1
     * with the same object which is somewhere in the
     * segmented labledPoints2 and orb2.
     *
     * NOTE that if the template or the true object match in dataset2
     * are smaller than 32 pixels across, the method may not find the
     * object very well so alternative methods should be used in that case
     * or pre-processing to correct that.
     *
     * NOTE also that if precise correspondence is needed, this method should
     * probably be followed by partial shape matcher to get better transformation
     * and then add transformed matching keypoints to that correspondence.
     *
     * NOT READY FOR USE yet.
     *
     * @param orb1
     * @param orb2
     * @param labeledPoints1
     * @param labeledPoints2
     * @return
     */
    public List<CorrespondenceList> match0(ORB orb1, ORB orb2, 
        Set<PairInt> labeledPoints1, List<Set<PairInt>> labeledPoints2) {
        
        /*
        uses the descriptors given and then optionally makes masks
        for them using the labeled points.
        -- visits each octave pair
        -- calculates cost of descriptors
        -- uses the segmentation to calculate every
        permutation of 2 pairs of points.
        -- filter out high cost pairs.
        -- filters out 2 pair combinations with transformation scales not near 1
        -- keeps only the top 10 percent cost of items
        from the 2 pair list.
        -- evaluates the transformation using the transformed
        keypoints cost difference, distance from nearest
        neighbor and number of matches
        -- keeps the best of each j
        -- further compares bestJs with SSDs of intersecting
        transformed point sets of the matching keypoints
        -- top of those best is the returned result
         */
        if (!orb1.getDescrChoice().equals(orb2.getDescrChoice())) {
            throw new IllegalStateException("orbs must contain same kind of descirptors");
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
        
        boolean useMasks = false;
        if (useMasks) {
            // initialize the masks, but discard the maps
            TObjectIntMap<PairInt> pointLabels1 = new TObjectIntHashMap<PairInt>();
            Set<PairInt> set = labeledPoints1;
            for (PairInt p : set) {
                pointLabels1.put(p, 0);
            }
            TObjectIntMap<PairInt> pointLabels2 = new TObjectIntHashMap<PairInt>();
            for (int i = 0; i < labeledPoints2.size(); ++i) {
                set = labeledPoints2.get(i);
                for (PairInt p : set) {
                    pointLabels2.put(p, i);
                }
            }
            orb1.createDescriptorMasks(pointLabels1);
            orb2.createDescriptorMasks(pointLabels2);
        }
        
        //TODO: may need to revise this or allow it as a method argument:
        int pixTolerance = 10;
        MatchedPointsTransformationCalculator tc = new MatchedPointsTransformationCalculator();
        Transformer transformer = new Transformer();
        
        TFloatList scales1 = extractScales(orb1.getScalesList());
        TFloatList scales2 = extractScales(orb2.getScalesList());
        
        if (Math.abs(scales1.get(0) - 1) > 0.01) {
            throw new IllegalArgumentException("logic depends upon first scale" + " level being '1'");
        }
        if (Math.abs(scales2.get(0) - 1) > 0.01) {
            throw new IllegalArgumentException("logic depends upon first scale" + " level being '1'");
        }
        
        // a rough estimate of maximum number of matchable points in any
        //     scale dataset comparison
        final int nMaxMatchable = Math.round(0.5F * calculateNMaxMatchable(orb1.getKeyPoint1List(), orb2.getKeyPoint1List()));
        //TODO: allow a factor to be passed in
        System.out.println("nMaxMatchable=" + nMaxMatchable);
        int nMax1 = maxSize(orb1.getKeyPoint1List());
        int nMax2 = maxSize(orb2.getKeyPoint1List());
        int nMax = nMax1 * nMax2;
        // --- best cost data ----
        double minCostTotal = Double.MAX_VALUE;
        double minCost1 = Double.MAX_VALUE;
        double minCost2 = Double.MAX_VALUE;
        double minCost3 = Double.MAX_VALUE;
        float minCostTScale = Float.MAX_VALUE;
        //runtime complexity of this vector depends upon the number of items
        // it is currently holding, so can set the capacity high and fill vector only
        // with items within bitTolerance of best, but too high might affect jvm
        // performance.
        // (note, can optimize this for very large results by occassionally ejecting
        // all values with cost > best + bitTolerance.)
        // TODO: a safe size is to set capacity to the number of unique
        // transformation parameter sets, but since that isn't known
        // until later without refactoring here, will make an assumption for now,
        // that size 100 is generous for number of top solutions.
        FixedSizeSortedVector<CObject3> minVec = new FixedSizeSortedVector<CObject3>(1, CObject3.class);
        int templateSize = calculateObjectSize(labeledPoints1);
        // populated on demand
        TObjectIntMap<OneDIntArray> labeledPointsSizes2 = 
            new TObjectIntHashMap<OneDIntArray>();
        for (int i = 0; i < scales1.size(); ++i) {
            //for (int i = 2; i < 3; ++i) {
            float scale1 = scales1.get(i);
            // coords are in ref frame of scale=1 of their pyramids
            TIntList kpX1 = orb1.getKeyPoint1List().get(i);
            TIntList kpY1 = orb1.getKeyPoint0List().get(i);
            int n1 = kpX1.size();
            TwoDFloatArray octaveImg1 = orb1.getPyramidImages().get(i);
            float diag1 = (float) Math.sqrt(octaveImg1.a.length * octaveImg1.a[0].length);
            final double maxDist = diag1;
            // create data structures in scaled reference frame
            TObjectIntMap<PairInt> p1KPIndexMap = new TObjectIntHashMap<PairInt>();
            TIntList kpX1_2 = new TIntArrayList(n1);
            TIntList kpY1_2 = new TIntArrayList(n1);
            for (int i3 = 0; i3 < n1; ++i3) {
                int x = Math.round((float) kpX1.get(i3) / scale1);
                int y = Math.round((float) kpY1.get(i3) / scale1);
                kpX1_2.add(x);
                kpY1_2.add(y);
                p1KPIndexMap.put(new PairInt(x, y), i3);
            }
            List<TIntList> pointIndexLists1 = new ArrayList<TIntList>();
            int ns = 1;
            for (int i3 = 0; i3 < ns; ++i3) {
                pointIndexLists1.add(new TIntArrayList());
            }
            TObjectIntMap<PairInt> pointLabels1 = new TObjectIntHashMap<PairInt>();
            Set<PairInt> set = labeledPoints1;
            Set<PairInt> setScaled = new HashSet<PairInt>();
            TIntList list = pointIndexLists1.get(0);
            assert (list != null);
            for (PairInt p : set) {
                int x = Math.round((float) p.getX() / scale1);
                int y = Math.round((float) p.getY() / scale1);
                PairInt p2 = new PairInt(x, y);
                pointLabels1.put(p2, 0);
                int idx = p1KPIndexMap.get(p2);
                list.add(idx);
                setScaled.add(p2);
            }
            Set<PairInt> shape = new HashSet<PairInt>(setScaled);
            int objDimension = (int) Math.round((float) templateSize / (float) scale1);
            int limit = Math.round(1.15F * objDimension);
            int limitSq = limit * limit;
            PairIntArray a1 = new PairIntArray(kpX1_2.size());
            TIntList a1Indexes = new TIntArrayList(kpX1_2.size());
            for (int ii = 0; ii < kpX1.size(); ++ii) {
                int x = kpX1.get(ii);
                int y = kpY1.get(ii);
                a1.add(x, y);
                a1Indexes.add(ii);
            }
            for (int j = 0; j < scales2.size(); ++j) {
                //for (int j = 0; j < 1; ++j) {
                float scale2 = scales2.get(j);
                // coords are in ref frame of scale=1 of their pyramids
                TIntList kpX2 = orb2.getKeyPoint1List().get(j);
                TIntList kpY2 = orb2.getKeyPoint0List().get(j);
                int n2 = kpX2.size();
                // create data structures in scaled reference frame
                TObjectIntMap<PairInt> p2KPIndexMap = new TObjectIntHashMap<PairInt>();
                TObjectIntMap<PairInt> p2KPIndexMap_2 = new TObjectIntHashMap<PairInt>();
                TIntList kpX2_2 = new TIntArrayList(n2);
                TIntList kpY2_2 = new TIntArrayList(n2);
                for (int j3 = 0; j3 < n2; ++j3) {
                    int x = Math.round((float) kpX2.get(j3) / scale2);
                    int y = Math.round((float) kpY2.get(j3) / scale2);
                    kpX2_2.add(x);
                    kpY2_2.add(y);
                    p2KPIndexMap_2.put(new PairInt(x, y), j3);
                    p2KPIndexMap.put(new PairInt(kpX2.get(j3), kpY2.get(j3)), j3);
                }
                List<TIntList> pointIndexLists2 = new ArrayList<TIntList>();
                int ns2 = labeledPoints2.size();
                for (int j3 = 0; j3 < ns2; ++j3) {
                    pointIndexLists2.add(new TIntArrayList());
                }
                TObjectIntMap<PairInt> pointLabels2 = new TObjectIntHashMap<PairInt>();
                for (int j3 = 0; j3 < ns2; ++j3) {
                    Set<PairInt> set2 = labeledPoints2.get(j3);
                    TIntList list2 = pointIndexLists2.get(j3);
                    assert (list2 != null);
                    for (PairInt p : set2) {
                        int x = Math.round((float) p.getX() / scale2);
                        int y = Math.round((float) p.getY() / scale2);
                        PairInt p2 = new PairInt(x, y);
                        pointLabels2.put(p2, j3);
                        int idx = p2KPIndexMap_2.get(p2);
                        list2.add(idx);
                    }
                }
                TwoDFloatArray octaveImg2 = orb2.getPyramidImages().get(j);
                debugPrint(octaveImg1, octaveImg2, kpX1_2, kpY1_2, kpX2_2, kpY2_2, i, j);
                int maxX2 = orb2.getPyramidImages().get(0).a[0].length;
                int maxY2 = orb2.getPyramidImages().get(0).a.length;
                int maxX2_2 = octaveImg2.a[0].length;
                int maxY2_2 = octaveImg2.a.length;
                NearestNeighbor2D nn2 = new NearestNeighbor2D(makeSet(kpX2, kpY2), maxX2 + limit, maxY2 + limit);
                int nTot = n1 * n2;
                //use descriptors with params here to reduce paramsList
                int[][] cost = null;
                if (useMasks) {
                    ORB.Descriptors[] desc1 = getDescriptors(orb1, i);
                    ORB.Descriptors[] desc2 = getDescriptors(orb2, j);
                    cost = ORB.calcMaskedDescriptorCostMatrixes(desc1, desc2, orb1.getDescriptorsMaskList().get(i), orb2.getDescriptorsMaskList().get(j))[1].a;
                } else {
                    ORB.Descriptors[] desc1 = getDescriptors(orb1, i);
                    ORB.Descriptors[] desc2 = getDescriptors(orb2, j);
                    cost = ORB.calcDescriptorCostMatrix(desc1, desc2);
                }
      
                //combinations of pairs with same labels
                // storing them all to reduce nesting
                // quadint is idx1, idx2, idx3, idx4
                //TODO: can use the cost to more quickly filter the
                // pairs at creation time
                List<QuadInt> pairIndexes = createPairLabelIndexes(cost, nBands, pointIndexLists1, kpX1_2, kpY1_2, pointIndexLists2, kpX2_2, kpY2_2);
                System.out.println("i=" + i + " j=" + j + " nPairs=" + pairIndexes.size());
                FixedSizeSortedVector<CObject4> vecP = new FixedSizeSortedVector<CObject4>(100, //Math.round(0.1f * pairIndexes.size()),
                //Math.round(0.01f * pairIndexes.size()),
                CObject4.class);
                for (int ipi = 0; ipi < pairIndexes.size(); ++ipi) {
                    QuadInt q = pairIndexes.get(ipi);
                    int t1X = kpX1_2.get(q.getA());
                    int t1Y = kpY1_2.get(q.getA());
                    int t2X = kpX1_2.get(q.getB());
                    int t2Y = kpY1_2.get(q.getB());
                    int s1X = kpX2_2.get(q.getC());
                    int s1Y = kpY2_2.get(q.getC());
                    int s2X = kpX2_2.get(q.getD());
                    int s2Y = kpY2_2.get(q.getD());
                    // transform dataset 1 into frame 2
                    TransformationParameters params = tc.calulateEuclidean(t1X, t1Y, t2X, t2Y, s1X, s1Y, s2X, s2Y, 0, 0);
                    float tScale = params.getScale();
                    if (Math.abs(tScale - 1.0) > 0.15) {
                        continue;
                    }
                    int idx1_1 = p1KPIndexMap.get(new PairInt(t1X, t1Y));
                    int idx1_2 = p1KPIndexMap.get(new PairInt(t2X, t2Y));
                    int idx2_1 = p2KPIndexMap_2.get(new PairInt(s1X, s1Y));
                    int idx2_2 = p2KPIndexMap_2.get(new PairInt(s2X, s2Y));
                    // a filter for objects too large to be the template object in
                    //    dataset 1.
                    // caveat is that cannot use partial shape matcher on all
                    //    results in same manner if filter this one out, but it's
                    //    the right logic if not oversegmented or blended into
                    //    other objects.
                    int label2 = pointLabels2.get(new PairInt(kpX2.get(q.getC()), kpY2.get(q.getC())));
                    if (labeledPoints2.get(label2).size() < 2) {
                        continue;
                    }
                    OneDIntArray key = new OneDIntArray(new int[]{label2});
                    if (!labeledPointsSizes2.containsKey(key)) {
                        Set<PairInt> set2 = labeledPoints2.get(label2);
                        if (set2.size() < 2) {
                            continue;
                        }
                        int sz = calculateObjectSize(set2);
                        labeledPointsSizes2.put(key, sz);
                    }
                    int regionSize = labeledPointsSizes2.get(key);
                    if (regionSize > (1.5 * templateSize)) {
                        continue;
                    }
                    int sum = cost[idx1_1][idx2_1] + cost[idx1_2][idx2_2];
                    CObject4 cObj = new CObject4(sum, params, q);
                    boolean added = vecP.add(cObj);
                }
                System.out.println("for i=" + i + " j=" + j + " filtered nPairs=" + vecP.getNumberOfItems());
                double minCostJTotal = Double.MAX_VALUE;
                double minCostJ1 = Double.MAX_VALUE;
                double minCostJ2 = Double.MAX_VALUE;
                double minCostJ3 = Double.MAX_VALUE;
                float minCostJTScale = Float.MAX_VALUE;
                FixedSizeSortedVector<CObject3> vecJ = new FixedSizeSortedVector<CObject3>(1, CObject3.class);
                // --- evaluate cost of all keypoints, transformed
                for (int ipi = 0; ipi < vecP.getNumberOfItems(); ++ipi) {
                    CObject4 c = vecP.getArray()[ipi];
                    TransformationParameters params = c.params;
                    float tScale = params.getScale();
                    QuadInt q = c.q;
                    int t1X = kpX1_2.get(q.getA());
                    int t1Y = kpY1_2.get(q.getA());
                    int t2X = kpX1_2.get(q.getB());
                    int t2Y = kpY1_2.get(q.getB());              
                    int s1X = kpX2_2.get(q.getC());
                    int s1Y = kpY2_2.get(q.getC());
                    int s2X = kpX2_2.get(q.getD());
                    int s2Y = kpY2_2.get(q.getD());                    
                    // ----- transform keypoints and sum the distance differences ----
                    PairIntArray tr1 = transformer.applyTransformation(params, a1);
                    // trim to image dimensions
                    tr1 = trimToImageBounds(octaveImg2, tr1);
                    if (tr1.getN() == 0) {
                        continue;
                    }
                    //the matched kpx1,kpy1 kpx2,kpy2 coordinate pairs
                    int[] mp1 = new int[kpX1.size()];
                    int[] mp2 = new int[kpX1.size()];
                    double[] distAndCount = sumKeypointDescAndDist(cost, 3, a1Indexes, tr1, kpX1, kpY1, nn2, p2KPIndexMap, maxX2, maxY2, pixTolerance, maxDist, mp1, mp2);
                    double sumDesc = distAndCount[0];
                    double sumDist = distAndCount[1];
                    int np = (int) distAndCount[2];
                    int count = np;
                    if (count < 2) {
                        continue;
                    }
                    if (count == 2 && (nMaxMatchable > 2 * count)) {
                        // TODO: may want to revise this while still discarding
                        //    false positives
                        continue;
                    }
                    if (np < mp1.length) {
                        mp1 = Arrays.copyOf(mp1, np);
                        mp2 = Arrays.copyOf(mp2, np);
                    }
                    if (count > nMaxMatchable) {
                        count = nMaxMatchable;
                    }
                    double cf = count;
                    if (cf > nMaxMatchable) {
                        cf = nMaxMatchable;
                    }
                    cf /= nMaxMatchable;
                    double sum3 = 1.0 - cf;
                    //sumDesc /= (double)count;
                    //sumDist /= (double)count;
                    sumDesc /= distAndCount[2];
                    sumDist /= distAndCount[2];
                    double sum = sumDesc + sumDist + sum3;
                    // if vecJ is filled and sum is not better than last item,
                    // continue
                    if (vecJ.getNumberOfItems() == vecJ.getFixedCapacity()) {
                        if (sum < vecJ.getArray()[vecJ.getNumberOfItems() - 1].cost) {
                            continue;
                        }
                    }
                    TIntSet labels2 = new TIntHashSet();
                    PairInt[] m1 = new PairInt[np];
                    PairInt[] m2 = new PairInt[mp1.length];
                    for (int j3 = 0; j3 < m1.length; ++j3) {
                        int idx1 = mp1[j3];
                        int idx2 = mp2[j3];
                        assert (idx1 < kpX1.size() && idx1 > -1);
                        assert (idx2 < kpX2.size() && idx2 > -1);
                        m1[j3] = new PairInt(kpX1.get(idx1), kpY1.get(idx1));
                        m2[j3] = new PairInt(kpX2.get(idx2), kpY2.get(idx2));
                        assert (labeledPoints1.contains(m1[j3]));
                        assert (p1KPIndexMap.get(new PairInt(kpX1_2.get(idx1), kpY1_2.get(idx1))) == idx1);
                        assert (p2KPIndexMap_2.get(new PairInt(kpX2_2.get(idx2), kpY2_2.get(idx2))) == idx2);
                        labels2.add(pointLabels2.get(m2[j3]));
                    }
                    // apply a size filter
                    OneDIntArray keys = new OneDIntArray(
                        labels2.toArray(new int[labels2.size()]));
                    Arrays.sort(keys.a);
                    if (!labeledPointsSizes2.containsKey(keys)) {
                        Set<PairInt> combined = new HashSet<PairInt>();
                        for (int k = 0; k < keys.a.length; ++k) {
                            combined.addAll(labeledPoints2.get(keys.a[k]));
                        }
                        if (combined.size() < 2) {
                            continue;
                        }
                        int sz = calculateObjectSize(combined);
                        labeledPointsSizes2.put(keys, sz);
                    }
                    int regionSize = labeledPointsSizes2.get(keys);
                    if (regionSize > (1.5 * templateSize)) {
                        continue;
                    }
                    CObject2 cObj2 = new CObject2(ipi, sum, sumDesc, sumDist, sum3, m1, m2);
                    CObject3 cObj = new CObject3(cObj2, sum, 0, params);
                    cObj.keypointCount = count;
                    boolean added = vecJ.add(cObj);
                    if (added) {
                        minCostJTotal = sum;
                        minCostJ1 = sumDesc;
                        minCostJ2 = sumDist;
                        minCostJ3 = sum3;
                        minCostJTScale = tScale;
                        System.out.println(String.format("i=%d j=%d ipi=%d ts=%.2f  c=%.2f c1=%.2f c2=%.2f c3=%.2f count=%d", i, j, ipi, tScale, (float) sum, (float) sumDesc, (float) sumDist, (float) sum3, count));
                        if (true) {
                            CorrespondencePlotter plotter = new CorrespondencePlotter(ORB.convertToImage(orb1.getPyramidImages().get(i)), ORB.convertToImage(orb2.getPyramidImages().get(j)));
                            for (int ii = 0; ii < cObj.m1.length; ++ii) {
                                PairInt p1 = cObj.m1[ii];
                                PairInt p2 = cObj.m2[ii];
                                int x1 = Math.round((float) p1.getX() / scale1);
                                int y1 = Math.round((float) p1.getY() / scale1);
                                int x2 = Math.round((float) p2.getX() / scale2);
                                int y2 = Math.round((float) p2.getY() / scale2);
                                plotter.drawLineInAlternatingColors(x1, y1, x2, y2, 0);
                            }
                            String str = Integer.toString(i);
                            while (str.length() < 3) {
                                str = "0" + str;
                            }
                            String str2 = Integer.toString(j);
                            while (str2.length() < 3) {
                                str2 = "0" + str2;
                            }
                            str = str + "_" + str2;
                            try {
                                plotter.writeImage("_indiv_masked_corres2_" + str + "_" + ipi);
                            } catch (IOException ex) {
                                Logger.getLogger(ORB.class.getName()).log(Level.SEVERE, null, ex);
                            }
                        }
                    }
                } // end loop over paramsList
                if (vecJ.getNumberOfItems() == 0) {
                    continue;
                }
                if (false) {
                    //DEBUG
                    for (int k = 0; k < vecJ.getNumberOfItems(); ++k) {
                        CObject3 cobj = vecJ.getArray()[k];
                        CorrespondencePlotter plotter = new CorrespondencePlotter(ORB.convertToImage(orb1.getPyramidImages().get(i)), ORB.convertToImage(orb2.getPyramidImages().get(j)));
                        for (int ii = 0; ii < cobj.m1.length; ++ii) {
                            PairInt p1 = cobj.m1[ii];
                            PairInt p2 = cobj.m2[ii];
                            int x1 = Math.round((float) p1.getX() / scale1);
                            int y1 = Math.round((float) p1.getY() / scale1);
                            int x2 = Math.round((float) p2.getX() / scale2);
                            int y2 = Math.round((float) p2.getY() / scale2);
                            plotter.drawLineInAlternatingColors(x1, y1, x2, y2, 0);
                        }
                        String str = Integer.toString(i);
                        while (str.length() < 3) {
                            str = "0" + str;
                        }
                        String str2 = Integer.toString(j);
                        while (str2.length() < 3) {
                            str2 = "0" + str2;
                        }
                        str = str + "_" + str2;
                        try {
                            plotter.writeImage("_mindiv_masked_corres3_" + str + "_" + MiscDebug.getCurrentTimeFormatted());
                        } catch (IOException ex) {
                            Logger.getLogger(ORB.class.getName()).log(Level.SEVERE, null, ex);
                        }
                        System.out.println(String.format("* %d %d ts=%.2f  c=%.2f c1=%.2f c2=%.2f c3=%.2f", i, j, cobj.params.getScale(), (float) cobj.cost, (float) cobj.costDesc, (float) cobj.costDist, (float) cobj.costCount));
                    }
                }
                if (vecJ.getNumberOfItems() == 0) {
                    System.out.println("no matches for i=" + i + " j=" + j);
                    continue;
                }
                // if expand capacity of minVec, add up to capacity here
                minVec.add(vecJ.getArray()[0]);
            } // end loop over image j
        }
        if (minVec.getNumberOfItems() == 0) {
            return null;
        }
        List<CorrespondenceList> topResults = new ArrayList<CorrespondenceList>();
        for (int i = 0; i < minVec.getNumberOfItems(); ++i) {
            CObject3 a = minVec.getArray()[i];
            if (a.cost > minCostTotal) {
                break;
            }
            CorrespondenceList cor = new CorrespondenceList(a.params, a.m1, a.m2);
            topResults.add(cor);
        }
        return topResults;
    }

    /**
     * match template image and shape in orb1 and labeledPoints1
     * with the same object which is somewhere in the
     * segmented labledPoints2 and orb2.
     *
     * this method matches points on a segmented cell basis to calculate
       the minimum cost correspondence with an objective function 
       consisting of cost from an outer point chord difference matrix,
       cost from hsv orb descriptors of keypoints, and an epipolar projection
       to remove outliers and find matching inner points (and subsequent
       addition of the later costs to the total).
     
     * NOT READY FOR USE yet.
     *
     * @param orb1
     * @param orb2
     * @param labeledPoints1
     * @param labeledPoints2
     * @return
     */
    public List<CorrespondenceList> match0Epipolar(ORB orb1, ORB orb2, 
        Set<PairInt> labeledPoints1, List<Set<PairInt>> labeledPoints2) {
                
        if (!orb1.getDescrChoice().equals(orb2.getDescrChoice())) {
            throw new IllegalStateException("orbs must contain same kind of descirptors");
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

        SIGMA sigma = SIGMA.ZEROPOINTFIVE;
        
        float distTol = 5;
        float distMax = (float)(Math.sqrt(2) * distTol);
        
        EpipolarTransformer eTransformer = new EpipolarTransformer();
        
        // --- create ordered bounds1 array.
        //     NOTE that all octaves use coordinates based in the
        //     full reference frame, so only one bounds1 is needed for
        //     all octaves.
        PairIntArray bounds1 = createOrderedBounds(orb1, labeledPoints1, sigma);
        if (bounds1.getN() < 7) {
            throw new IllegalStateException("the boundary of object 1 "
                + " must have at least 7 points");
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
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        // -- initialize bounds2MapsList and populte on demand
        TIntObjectMap<PairIntArray> bounds2Maps 
            = new TIntObjectHashMap<PairIntArray>();
       
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
        
        // making a lookup map for keypoint indexes in points2 labeled sets
        List<TIntObjectMap<TIntSet>> labels2KPIdxsList = 
            new ArrayList<TIntObjectMap<TIntSet>>();
        for (int octave = 0; octave < scales2.size(); ++octave) {
            float scale2 = scales2.get(octave);
            TIntObjectMap<TIntSet> labels2KPIdxs = new TIntObjectHashMap<TIntSet>();
            labels2KPIdxsList.add(labels2KPIdxs);
            TObjectIntMap<PairInt> keypoints2IndexMap = kp2IdxMapList.get(octave);
            for (int i = 0; i < labeledPoints2.size(); ++i) {
                for (PairInt p : labeledPoints2.get(i)) {
                    if (keypoints2IndexMap.containsKey(p)) {
                        int kp2Idx = keypoints2IndexMap.get(p);
                        TIntSet kpIdxs = labels2KPIdxs.get(i);
                        if (kpIdxs == null) {
                            kpIdxs = new TIntHashSet();
                            labels2KPIdxs.put(i, kpIdxs);
                        }
                        kpIdxs.add(kp2Idx);
                    }
                }
            }     
        }
        
        // a cache for partial shape matcher and results
        TIntObjectMap<List<Object>> psmMap = new
            TIntObjectHashMap<List<Object>>();
        
        double maxChordAvg = Double.MIN_VALUE;
        
        double maxAvgDist = Double.MIN_VALUE;
        
        List<List<QuadInt>> correspondences = new ArrayList<List<QuadInt>>();
        TIntList nLabelKP2s = new TIntArrayList();
        TDoubleList descCosts = new TDoubleArrayList();
        TIntList nDesc = new TIntArrayList();
        TFloatList descNormalizations = new TFloatArrayList();
        TDoubleList epCosts = new TDoubleArrayList();
        TDoubleList chordCosts = new TDoubleArrayList();
        TIntList octs1 = new TIntArrayList();
        TIntList octs2 = new TIntArrayList();
        TIntList segIdxs = new TIntArrayList();

        for (int octave1 = 0; octave1 < scales1.size(); ++octave1) {
        //for (int octave1 = 1; octave1 < 2; ++octave1) {
            
            float scale1 = scales1.get(octave1);
            
            TObjectIntMap<PairInt> keypoints1IndexMap = kp1IdxMapList.get(octave1);
           
            float sz1 = calculateObjectSize(bounds1)/scale1;
            
            int nkp1 = orb1.getKeyPoint0List().get(octave1).size();
            int nb1 = bounds1.getN();
            
            float normDesc = nkp1 * nBands * 256;
            
            //for (int octave2 = 0; octave2 < scales2.size(); ++octave2) {
            for (int octave2 = 0; octave2 < 1; ++octave2) {
           
                float scale2 = scales2.get(octave2);
            
                TObjectIntMap<PairInt> keypoints2IndexMap = kp2IdxMapList.get(octave2);
                TIntObjectMap<TIntSet> labels2KPIdxs = 
                    labels2KPIdxsList.get(octave2);
                
                TwoDFloatArray img2 = orb2.getPyramidImages().get(octave2);
                
                TIntObjectIterator<TIntSet> iter2 = labels2KPIdxs.iterator();
                for (int i2 = 0; i2 < labels2KPIdxs.size(); ++i2) {
                    iter2.advance();
                    
                    int segIdx = iter2.key();
                    TIntSet kp2Idxs = iter2.value();
                    
                    float sz2 = sizes2Maps.get(segIdx)/scale2;
                    if (sz2 == 0) {
                        continue;
                    }
 System.out.println("octave1=" + octave1 + " octave2=" + octave2 + 
   " sz1=" + sz1 + " sz2=" + sz2 + " segIdx=" + segIdx);
                    if ((sz1 > sz2 && Math.abs(sz1 / sz2) > 1.2) || 
                        (sz2 > sz1 && Math.abs(sz2 / sz1) > 1.2)) {
                        continue;
                    }
   
                    PairIntArray bounds2 = getOrCreateOrderedBounds(img2,
                        bounds2Maps, segIdx, labeledPoints2.get(segIdx), sigma);
       
                    if (bounds2 == null || bounds2.getN() < 7) {
                        continue;
                    }
             
                    List<Object> psmObj = psmMap.get(segIdx);
                    if (psmObj == null) {
                        PartialShapeMatcher matcher = new PartialShapeMatcher();
                        matcher.overrideSamplingDistance(1);
                        matcher._overrideToThreshhold(0.2f);
                        matcher.setToRemoveOutliers();
                        matcher.overrideToStoreMatrix();
                        PartialShapeMatcher.Result result = matcher.match(
                            bounds1, bounds2);
                        psmObj = new ArrayList<Object>();
                        psmObj.add(matcher);
                        if (result != null) {
                            psmObj.add(result);
                        }
                        psmMap.put(segIdx, psmObj);
                    }
                    if (psmObj.size() == 1) {
                        continue;
                    }
                    PartialShapeMatcher matcher = (PartialShapeMatcher)psmObj.get(0);
                    PartialShapeMatcher.Result result = 
                        (PartialShapeMatcher.Result)psmObj.get(1);

                    if ((matcher.getStoredEpipolarFit() == null) 
                        || (result.getNumberOfMatches() < 3)) {
                        continue;
                    }
      
                    int nr = result.getNumberOfMatches();
                                        
                    PairIntArray m1 = new PairIntArray(nr);
                    PairIntArray m2 = new PairIntArray(nr);
                    Set<PairInt> matched1 = new HashSet<PairInt>();
                    Set<PairInt> matched2 = new HashSet<PairInt>();
                    for (int j = 0; j < nr; ++j) {
                        int idx1 = result.idx1s.get(j);
                        int idx2 = result.idx2s.get(j);
                        int x1 = bounds1.getX(idx1);
                        int y1 = bounds1.getY(idx1);
                        int x2 = bounds2.getX(idx2);
                        int y2 = bounds2.getY(idx2);
                        m1.add(x1, y1);
                        m2.add(x2, y2);
                        matched1.add(new PairInt(x1, y1));
                        matched2.add(new PairInt(x2, y2));
                    }
                    
 //TODO: correct error in including points outside of segmentation region        
                                        
                    SimpleMatrix fm = matcher.getStoredEpipolarFit().getFundamentalMatrix();
                    
                    //TODO: this method needs to be tested...normalization effects...
                    // sum, avg, max
                    double[] avgAndMaxDist = sumAndMaxEPDist(fm, m1, m2);
                    if (avgAndMaxDist[2] > maxAvgDist) {
                        maxAvgDist = avgAndMaxDist[2];
                    }    
                    
                    // key=keypoint in this labeled region, value=kp2Index
                    PairIntArray unmatchedKP2 = new PairIntArray();
                    TObjectIntMap<PairInt> unmatchedKP2Idxs = 
                        new TObjectIntHashMap<PairInt>();
                    TIntIterator iter = kp2Idxs.iterator();
                    while (iter.hasNext()) {
                        int kp2Idx = iter.next();
                        int x = orb2.getKeyPoint1List().get(octave2).get(kp2Idx);
                        int y = orb2.getKeyPoint0List().get(octave2).get(kp2Idx);
                        PairInt p = new PairInt(x, y);
                        if (!matched2.contains(p)) {
                            unmatchedKP2Idxs.put(p, kp2Idx);
                            unmatchedKP2.add(x, y);
                        }
                    }
                    
                    PairIntArray unmatchedKP1 = new PairIntArray();
                    TObjectIntMap<PairInt> unmatchedKP1Idxs = 
                        new TObjectIntHashMap<PairInt>();
              
                    TObjectIntIterator<PairInt> iter1 = keypoints1IndexMap.iterator();
                    for (int j = 0; j < keypoints1IndexMap.size();
                        ++j) {
                        iter1.advance();
                        PairInt p = iter1.key();
                        int kpIdx1 = iter1.value();
                        if (!matched1.contains(p)) {
                            unmatchedKP1Idxs.put(p, kpIdx1);
                            unmatchedKP1.add(p.getX(), p.getY());
                        }
                    }
         
                    {// DEBUG, print bounds1 and unmatchedkp1
                        Image img1 = ORB.convertToImage(
                            orb1.getPyramidImages().get(octave1));
                        for (int ii = 0; ii < m1.getN(); ++ii) {
                            int x = Math.round((float)m1.getX(ii)/scale1);
                            int y = Math.round((float)m1.getY(ii)/scale1);
                            ImageIOHelper.addPointToImage(x, y, img1, 1, 0, 255, 0);
                        }
                        for (int ii = 0; ii < unmatchedKP1.getN(); ++ii) {
                            int x = Math.round((float)unmatchedKP1.getX(ii)/scale1);
                            int y = Math.round((float)unmatchedKP1.getY(ii)/scale1);
                            ImageIOHelper.addPointToImage(x, y, img1, 1, 255, 0, 0);
                        }
                        MiscDebug.writeImage(img1, "_TMP1_" + octave1 + "_" + 
                            MiscDebug.getCurrentTimeFormatted());
                        //====
                        img1 = ORB.convertToImage(
                            orb2.getPyramidImages().get(octave2));
                        for (int ii = 0; ii < m2.getN(); ++ii) {
                            int x = Math.round((float)m2.getX(ii)/scale2);
                            int y = Math.round((float)m2.getY(ii)/scale2);
                            ImageIOHelper.addPointToImage(x, y, img1, 1, 0, 255, 0);
                        }
                        for (int ii = 0; ii < unmatchedKP2.getN(); ++ii) {
                            int x = Math.round((float)unmatchedKP2.getX(ii)/scale2);
                            int y = Math.round((float)unmatchedKP2.getY(ii)/scale2);
                            ImageIOHelper.addPointToImage(x, y, img1, 1, 255, 0, 0);
                        }
                        String str3 = Integer.toString(segIdx);
                        while (str3.length() < 3) {
                            str3 = "0" + str3;
                        }
                        MiscDebug.writeImage(img1, "_TMP2_" + octave2 + "_" + 
                            str3 + "_" + MiscDebug.getCurrentTimeFormatted());
                    }
                
                    ORB.Descriptors[] desc1 = getDescriptors(orb1, octave1);
                    ORB.Descriptors[] desc2 = getDescriptors(orb2, octave2);
                    int[][] costD = ORB.calcDescriptorCostMatrix(desc1, desc2);
                    
                    // -- use epipolar fundamental matrix to add unmatched
                    //       points from the segmented cell's keypoints
                   
                    // output variable to hold sums and count
                    // 0 = totalChordDiffSum
                    // 1 = max avg chord diff
                    // 2 = totalDistance
                    // 3 = max avg total dist
                    // 4 = totalDescrSum
                    // 5 = nDescr
                    double[] output = new double[6];
                    
                    List<PairInt> addedKPIdxs = matchUsingFM(orb1, orb2, costD, 
                        octave1, octave2, bounds1, bounds2, matcher, result,
                        keypoints1IndexMap, keypoints2IndexMap,
                        fm, unmatchedKP1, unmatchedKP2,
                        unmatchedKP1Idxs, unmatchedKP2Idxs,
                        nBands, normDesc, distTol, output);
                    
                    if (output[1] > maxChordAvg) {
                        maxChordAvg = output[1];
                    }
                    if (output[3] > maxAvgDist) {
                        maxAvgDist = output[3];
                    }
                    
                    // add the boundary matching epipolar dist sum:
                    output[2] += avgAndMaxDist[0];
                    
                    System.out.println("nAdded inner points=" + addedKPIdxs.size());

                    // --- build combined correspondence and sums
                    
                    int nTot = result.getNumberOfMatches() + addedKPIdxs.size();
                    List<QuadInt> corres = new ArrayList<QuadInt>(nTot);
                    // for any point in result that is a keypoint,
                    //    add the descriptor cost to totalDescrSum
                    for (int j = 0; j < result.getNumberOfMatches(); ++j) {
                        int idx1 = result.idx1s.get(j);
                        int idx2 = result.idx2s.get(j);
                        
                        int x1 = bounds1.getX(idx1);
                        int y1 = bounds1.getY(idx1);
                        PairInt p1 = new PairInt(x1, y1);

                        int x2 = bounds2.getX(idx2);
                        int y2 = bounds2.getY(idx2);
                        PairInt p2 = new PairInt(x2, y2);

                        if (keypoints1IndexMap.containsKey(p1) && 
                            keypoints2IndexMap.containsKey(p2)) {
                            int kpIdx1 = keypoints1IndexMap.get(p1);
                            int kpIdx2 = keypoints2IndexMap.get(p2);
                            float c = costD[kpIdx1][kpIdx2];
                            output[4] += c;
                            output[5]++;
                        }
                        
                        // coords are in full reference frame
                        corres.add(new QuadInt(p1, p2)); 
                    }
                   
                    for (int j = 0; j < addedKPIdxs.size(); ++j) {
                        
                        int kpIdx1 = addedKPIdxs.get(j).getX();
                        int kpIdx2 = addedKPIdxs.get(j).getY();
                        
                        int x1 = orb1.getKeyPoint1List().get(octave1).get(kpIdx1);
                        int y1 = orb1.getKeyPoint0List().get(octave1).get(kpIdx1);
                        
                        int x2 = orb2.getKeyPoint1List().get(octave2).get(kpIdx2);
                        int y2 = orb2.getKeyPoint0List().get(octave2).get(kpIdx2);
                        
                        corres.add(new QuadInt(x1, y1, x2, y2));                        
                    }
                    assert(corres.size() == nTot);

                    // output variable to hold sums and count
                    // 0 = totalChordDiffSum
                    // 1 = max avg chord diff
                    // 2 = totalDistance
                    // 3 = max avg total dist
                    // 4 = totalDescrSum
                    // 5 = nDescr
                                       
                    correspondences.add(corres);
                    descCosts.add(output[4]);
                    nDesc.add((int)output[5]);
                    nLabelKP2s.add(kp2Idxs.size());
                    epCosts.add(output[2]);
                    chordCosts.add(output[0]);
                    octs1.add(octave1);
                    octs2.add(octave2);
                    segIdxs.add(segIdx);
                    descNormalizations.add(normDesc);
                    
                }// end loop over octave2's segIdx
            }// end loop over octave2            
        }  // end loop over octave1
       
        
        int nC = correspondences.size();
        
        /*        
        "salukwzde distance" separate for
        descriptor cost, chords cost, and the opipolar dist from
        model then add them
           -- descr component max matchable number is nLabelKP2s
           -- ep distances component max matchable number is 
                bounds1.getN() + nLabelKP2s (note: should remove overlapping)
           -- chords component max matchable number is 
                bounds1.getN() + nLabelKP2s (note: should remove overlapping)
        */
        int[] indexes = new int[nC];
        float[] costs = new float[nC];
        for (int i = 0; i < nC; ++i) {
            
            int octave1 = octs1.get(i);
            int octave2 = octs2.get(i);
            
            // calculate "fraction of whole" for hsv keypoint descriptors
            int nKP2 = nLabelKP2s.get(i);
            float f1 = 1.f - ((float)nDesc.get(i)/(float)nKP2);
            
            //calculate the cost of hsv kp descriptors
            float d1 = 1.f - ((nBands * 256.f 
                - (float)(descCosts.get(i)/(float)nDesc.get(i)))
                /descNormalizations.get(i));
            if (descNormalizations.get(i) == 0) {
                d1 = 1;
            }
            
            float sd1 = f1 * f1 + d1 * d1;
            
            // ------ chords --------
            
            float nb1 = bounds1.getN();
            float n = correspondences.get(i).size();
            float f2 = 1.f - (n/(nKP2 + nb1));
            
            float d2 = (float)((chordCosts.get(i)/n)/maxChordAvg);
            
            float sd2 = f2 * f2 + d2 * d2;
            
            // ----- epipolar line distances -----
            float d3 = (float)((epCosts.get(i)/n)/maxAvgDist);
            
            float sd3 = f2 * f2 + d3 * d3;
            
            // add in quadrature or linearly...
            double tot = sd1 + sd2 + sd3;
        
             System.out.println(String.format(
 "octave1=%d octave2=%d segIdx=%d nCor=%d ch=%.2f, normch=%.2f normep=%.2f normdesc=%.2f sd1=%.2f sd2=%.2f sd3=%.2f nd=%d nKP2=%d tot=%.2f",
                octave1, octave2, segIdxs.get(i), (int)n,
                (float)chordCosts.get(i),
                (float)d2,
                (float)d3, (float)d1, 
                sd1, sd2, sd3,
                nDesc.get(i), nKP2,
                (float)tot));
            
            indexes[i] = i;
            costs[i] = (float)tot;
        }
    
        QuickSort.sortBy1stArg(costs, indexes);
        
        List<CorrespondenceList> results = new ArrayList<CorrespondenceList>();
        for (int i = 0; i < costs.length; ++i) {
            
            int idx = indexes[i];
            
            List<QuadInt> qs = correspondences.get(idx);

            // points are in full reference frame            
            results.add(new CorrespondenceList(qs));
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
            PairIntArray bounds = imageProcessor.extractSmoothedOrderedBoundary(new HashSet(set1), sigma, img.getWidth(), img.getHeight());
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
            PairIntArray bounds = imageProcessor.extractSmoothedOrderedBoundary(new HashSet(set1), sigma, img.getWidth(), img.getHeight());
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
                        double scoreNorm = (3 * 256 - minC) / maxCost;
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

    private static int[][] greedyMatch(List<PairInt> keypoints1, List<PairInt> keypoints2, int[][] cost) {
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

    private static double[] sumKeypointDistanceDifference(TIntList a2Indexes, PairIntArray tr2, TIntList kpX2, TIntList kpY2, NearestNeighbor2D nn, TransformationParameters params, int maxX, int maxY, int pixTolerance, double maxDist, int[] m1x, int[] m1y, int[] m2x, int[] m2y) {
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

    private static double[] sumKeypointDescAndDist(int[][] cost, int nBands, TIntList a1Indexes, PairIntArray tr1, TIntList kpX1, TIntList kpY1, NearestNeighbor2D nn2, TObjectIntMap<PairInt> p2KPIndexMap, TransformationParameters params, int maxX2, int maxY2, int pixTolerance, double maxDist, PairInt[] m1, PairInt[] m2) {
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

    private static double[] sumKeypointDescAndDist(int[][] cost, int nBands, TIntList a1Indexes, PairIntArray tr1, TIntList kpX1, TIntList kpY1, NearestNeighbor2D nn2, TObjectIntMap<PairInt> p2KPIndexMap, int maxX2, int maxY2, int pixTolerance, double maxDist, int[] m1, int[] m2) {
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

    private static PairIntArray trimToImageBounds(TwoDFloatArray octaveImg, PairIntArray a) {
        int n0 = octaveImg.a.length;
        int n1 = octaveImg.a[0].length;
        PairIntArray b = new PairIntArray(a.getN());
        for (int i = 0; i < a.getN(); ++i) {
            int x = a.getX(i);
            int y = a.getY(i);
            if (x < 0 || x > (n1 - 1)) {
                continue;
            } else if (y < 0 || y > (n0 - 1)) {
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

    private PairIntArray getOrCreateOrderedBounds(TwoDFloatArray img, 
        TIntObjectMap<PairIntArray> boundsMap, int segIdx, 
        Set<PairInt> set, SIGMA sigma) {
        
        PairIntArray bounds = boundsMap.get(segIdx);
        if (bounds != null) {
            return bounds;
        }
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        bounds = imageProcessor.extractSmoothedOrderedBoundary(
            set, sigma, img.a[0].length, img.a.length);
        
        boundsMap.put(segIdx, bounds);
    
        {
            MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
            double[] xyCen = curveHelper.calculateXYCentroids(bounds);
            System.out.println("bounds center=" + (int)xyCen[0] + "," + 
                (int)xyCen[1] + " size_full=" + 
                calculateObjectSize(bounds));
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
        PairIntArray bounds1, PairIntArray bounds2,
        PartialShapeMatcher matcher,
        PartialShapeMatcher.Result result,
        TObjectIntMap<PairInt> keypoints1IndexMap, 
        TObjectIntMap<PairInt> keypoints2IndexMap,
        SimpleMatrix fm, 
        PairIntArray unmatchedKP1, PairIntArray unmatchedKP2,
        TObjectIntMap<PairInt> unmatchedKP1Idxs,
        TObjectIntMap<PairInt> unmatchedKP2Idxs,
        int nBands, float normDesc, float distTol, double[] output) {
        
        // output variable to hold sums and count
        // 0 = totalChordDiffSum
        // 1 = max avg chord diff
        // 2 = totalDistance
        // 3 = max avg total dist
        // 4 = totalDescrSum
        // 5 = nDescr
        
        double maxAvgChord = Double.MIN_VALUE;
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
                    
                    dist = (float)d;

                    // normalized descriptor cost
                    descCost = 1.f - ((nBands * 256.f - 
                        costD[kp1Idx][kp2Idx])
                        /normDesc);

                    eDist.add(dist);
                    totalCost.add(dist + descCost);
                    indexes.add(indexes.size());
                    idx1s.add(kp1Idx);
                    idx2s.add(kp2Idx);
                }
            }
        }

        QuickSort.sortBy1stArg(totalCost, indexes);

        // choose 2 reference points from result,
        //   preferably 2 that are keypoints w/ lowest descr costs 
        //   and are far from each other
        // returns results as 2 quadints of paired x1,y1,x2,y2
        QuadInt[] resultRefs = choose2ReferencePoints(result,
            bounds1, bounds2, 
            keypoints1IndexMap, keypoints2IndexMap, costD);

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
            // output:
            // 0 = totalChordDiffSum
            // 1 = max avg chord diff
            // 2 = totalDistance
            // 3 = max avg total dist
            // 4 = totalDescrSum
            // 5 = nDescr
            added1.add(kpIdx1);
            added2.add(kpIdx2);
            addedKPIdxs.add(new PairInt(kpIdx1, kpIdx2));
            output[2] += eDist.get(idx);

            d = eDist.get(idx);
            if (d > maxAvgDist) {
                maxAvgDist = d;
            }
            
            float descrCost = totalCost.get(j) - eDist.get(idx);
            output[4] += descrCost;
            output[5]++;

            int x1 = orb1.getKeyPoint1List().get(octave1).get(kpIdx1);
            int y1 = orb1.getKeyPoint0List().get(octave1).get(kpIdx1);

            int x2 = orb2.getKeyPoint1List().get(octave2).get(kpIdx2);
            int y2 = orb2.getKeyPoint0List().get(octave2).get(kpIdx2);

            // calc chord diff for the new points using 2 reference
            // points from result.
            double chordDiff = matcher.
                calculateAChordDifference(
                resultRefs[0].getA(), resultRefs[0].getB(),
                resultRefs[1].getA(), resultRefs[1].getB(),
                x1, y1,
                resultRefs[0].getC(), resultRefs[0].getD(),
                resultRefs[1].getC(), resultRefs[1].getD(), 
                x2, y2
            );
            output[0] += chordDiff;
            
            if (chordDiff > maxAvgChord) {
                maxAvgChord = chordDiff;
            }
        }
        
        output[1] = maxAvgChord;
        output[3] = maxAvgDist;
        // output:
        // 0 = totalChordDiffSum
        // 1 = max avg chord diff
        // 2 = totalDistance
        // 3 = max avg total dist
        // 4 = totalDescrSum
        // 5 = nDescr
            
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
