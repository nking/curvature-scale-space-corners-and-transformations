package algorithms.imageProcessing;

import algorithms.compGeometry.NearestPoints;
import algorithms.compGeometry.NearestPointsFloat;
import algorithms.imageProcessing.util.PairIntWithIndex;
import algorithms.misc.Histogram;
import algorithms.misc.HistogramHolder;
import algorithms.util.CornerArray;
import algorithms.util.Errors;
import algorithms.util.PairFloatArray;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.PairIntArrayWithColor;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 The code is implemented from interpreting several papers by the authors
 * Farzin Mokhtarian and Alan Mackworth.
 *
 * They prescribe a method for detecting features and corners that is scale
 * invariant, rotation and translation invariant and does not create
 * artifacts as a side effect.
 *
 * The method finds edges in an image and then calculates the curvature of
 * the edges to find "inflection" points along the curve.  Those points of
 * inflection as a function of scale parameter t are then findable in
 * another image that may have rotation or translation, for example, using
 * a search method such as A* to find the best matching features in t space.
 * The process of creating the scale based curves is repeated for increasing
 * sigma until no points are found with curvature = 0.
 *
 * The method uses the recursive and separable properties of Gaussians where
 * possible.  (NOTE, not finished implementing the recursive portion).
 * Also, the curves can be closed curves and derived from means other than
 * edge detectors.
 *
 * @author nichole
 */
public class CSSCornerMaker {

    protected boolean enableJaggedLineCorrections = true;

    protected float factorIncreaseForCurvatureMinimum = 1.f;

    protected boolean performWholeImageCorners = true;

    protected boolean doStoreCornerRegions = true;

    protected final int width;

    protected final int height;

    private Map<Integer, List<CornerRegion>> edgeCornerRegionMap = new
        HashMap<Integer, List<CornerRegion>>();

    protected Logger log = Logger.getLogger(this.getClass().getName());

    public CSSCornerMaker(int imageWidth, int imageHeight) {
        width = imageWidth;
        height = imageHeight;
    }

    public void enableJaggedLineCorrections() {
        enableJaggedLineCorrections = true;
    }
    public void disableJaggedLineCorrections() {
        enableJaggedLineCorrections = false;
    }

    public void doNotStoreCornerRegions() {
        doStoreCornerRegions = false;
    }

    public void increaseFactorForCurvatureMinimum(float factor) {
        factorIncreaseForCurvatureMinimum = factor;
    }

    public void resetFactorForCurvatureMinimum() {
        factorIncreaseForCurvatureMinimum = 1.f;
    }

    /**
     *
     * @param theEdges edges with points ordered counter-clockwise.
     * the ordering is important if the corner regions will be used.
     * @param junctions map w/ key=edge index,
     *     values = junctions on edge as (x,y) and edge curve indexes
     * @param doUseOutdoorMode
     * @param outputCorners
     * @return
     */
    public Map<PairIntArray, Map<SIGMA, ScaleSpaceCurve> >
    findCornersInScaleSpaceMaps(final List<PairIntArray> theEdges,
        Map<Integer, Set<PairIntWithIndex>> junctions,
        final boolean doUseOutdoorMode, final PairIntArray outputCorners) {

        Map<PairIntArray, Map<SIGMA, ScaleSpaceCurve> > scaleSpaceMaps
            = new HashMap<PairIntArray, Map<SIGMA, ScaleSpaceCurve> >();

        // perform the analysis on the edgesScaleSpaceMaps
        for (int i = 0; i < theEdges.size(); i++) {

            final PairIntArray edge = theEdges.get(i);

            boolean isClosedCurve = (edge instanceof PairIntArrayWithColor) &&
                (((PairIntArrayWithColor)edge).getColor() == 1);

            final Map<SIGMA, ScaleSpaceCurve> map =
                createLowUpperThresholdScaleSpaceMaps(edge);

            scaleSpaceMaps.put(edge, map);

            Set<PairIntWithIndex> edgeJunctions = junctions.get(Integer.valueOf(i));

            CornerArray edgeCorners = findCornersInScaleSpaceMap(edge, map,
                edgeJunctions, i, doUseOutdoorMode);

            // remove aliasing artifacts of a straight line with a single
            // step in it.
            // note that edgeCorners are already ordered by index.
            PostLineThinnerCorrections.removeSingleStairsAliasArtifacts(
                edgeCorners, isClosedCurve);

            log.log(Level.FINE,
                "{0}) number of corners adding ={1} for edge={2}",
                new Object[]{Integer.valueOf(i),
                    Integer.valueOf(edgeCorners.getN()),
                    Integer.valueOf(i)});

            //store xc and yc for the edge
            for (int ii = 0; ii < edgeCorners.getN(); ii++) {
                int x = Math.round(edgeCorners.getX(ii));
                int y = Math.round(edgeCorners.getY(ii));
                outputCorners.add(x, y);
            }

            if (doStoreCornerRegions) {

                int edgeNumber = i;
                
                Set<Integer> removeIndexes = new HashSet<Integer>();

                for (int ii = 0; ii < edgeCorners.getN(); ii++) {

                    int idx = edgeCorners.getInt(ii);
                    
                    if (removeIndexes.contains(Integer.valueOf(ii))) {
                        continue;
                    }
                    int x = Math.round(edgeCorners.getX(ii));
                    int y = Math.round(edgeCorners.getY(ii));

                    SIGMA sscSigma = edgeCorners.getSIGMA(ii);
                    ScaleSpaceCurve scaleSpace = map.get(sscSigma);
                    assert(scaleSpace != null);

                    assert(x == Math.round(scaleSpace.getX(idx)));
                    assert(y == Math.round(scaleSpace.getY(idx)));

                    boolean added = storeCornerRegion(edgeNumber, idx, scaleSpace);
                    
                    if (!added) {
                        continue;
                    }
                    
                    // -- merge adjacent or nearly adjacent corners --
                    
                    int ii2 = -1;
                    int idx2 = -1;
                    int dIdx = Integer.MAX_VALUE;
                    
                    if (isClosedCurve && (ii == 0)) {
                        ii2 = edgeCorners.getN() - 1;
                        idx2 = edgeCorners.getInt(ii2);
                        dIdx = edge.getN() - idx2 + idx;
                    } else if (isClosedCurve && (ii == (edgeCorners.getN() - 1))) {
                        ii2 = 0;
                        idx2 = edgeCorners.getInt(ii2);
                        dIdx = edge.getN() - idx + idx2;
                    } else if (ii < (edgeCorners.getN() - 1)) {
                        ii2 = ii + 1;
                        idx2 = edgeCorners.getInt(ii2);
                        dIdx = idx2 + idx;
                    }
                    if ((idx2 > -1) && (dIdx < 4)) {
                        if (removeIndexes.contains(Integer.valueOf(ii2))) {
                            continue;
                        }
                        ScaleSpaceCurve scaleSpace2 = map.get(edgeCorners.getSIGMA(ii2));
                        float k = scaleSpace.getK(idx);
                        float k2 = scaleSpace2.getK(idx2);
                        float kLimit = 0.1f * Math.abs((k + k2)/2.f);
                        //remove weakest
                        if (Math.abs(k - k2) > kLimit) {
                            if (Math.abs(k2) > Math.abs(k)) {
                                removeIndexes.add(Integer.valueOf(ii));
                            } else {
                                removeIndexes.add(Integer.valueOf(ii2));
                            }
                        } else {
                            //merge, replace ii, remove ii2
                            int x2 = Math.round(edgeCorners.getX(ii2));
                            int y2 = Math.round(edgeCorners.getY(ii2));
                            float xAvg = (edgeCorners.getX(ii) + edgeCorners.getX(ii2))/2.f;
                            float yAvg = (edgeCorners.getY(ii) + edgeCorners.getY(ii2))/2.f;
                            edgeCorners.set(ii, xAvg, yAvg, idx, sscSigma);
                            removeIndexes.add(Integer.valueOf(ii2));
                        }
                    }
                }
                if (removeIndexes.size() == 1) {
                    int rmIdx = removeIndexes.iterator().next().intValue();
                    edgeCorners.removeRange(rmIdx, rmIdx);
                } else if (removeIndexes.size() > 0) {
                    int[] rmIdxs = new int[removeIndexes.size()];
                    int count = 0;
                    for (Integer rmIndex : removeIndexes) {
                        rmIdxs[count] = rmIndex.intValue();
                        count++;
                    }
                    Arrays.sort(rmIdxs);
                    for (int j = (rmIdxs.length - 1); j > -1; --j) {
                        int rmIdx = rmIdxs[j];
                        edgeCorners.removeRange(rmIdx, rmIdx);
                    }
                }
            }
        }

        return scaleSpaceMaps;
    }

    /**
     * Construct scale space images (that is X(t, sigma), y(t, sigma), and
     * k(t, sigma) where t is the intervals of spacing along the curve
     * and is valued from 0 to 1 and k is the curvature).
     *
     * The range of the scale space maps is from sigma = 0.5 up to a maximum
     * value determined in getMaxSIGMAForECSS.
     *
     * The results are returned as a map keyed by sigma.
     */
    private Map<SIGMA, ScaleSpaceCurve> createLowUpperThresholdScaleSpaceMaps(
        final PairIntArray edge) {

        ScaleSpaceCurvature scaleSpaceHelper = new ScaleSpaceCurvature();

        Map<SIGMA, ScaleSpaceCurve> map = new HashMap<SIGMA, ScaleSpaceCurve>();

        SIGMA sigma = SIGMA.ZEROPOINTFIVE;

        //TODO: take a look at the highest resolution curvature array
        //   and see if can derive the same upper limit by S/N estimate
        //   and assuming want to smooth to 2 o3 3 times the noise
        SIGMA maxSIGMA = getMaxSIGMAForECSS(edge.getN());

        // this increases by a factor of sqrt(2)
        float resultingSigma = SIGMA.getValue(sigma);

        ScaleSpaceCurve lastCurve = null;

        while (sigma.compareTo(maxSIGMA) < 1) {

            ScaleSpaceCurve curve;

            if (lastCurve == null) {
                curve = scaleSpaceHelper.computeCurvature(edge, sigma,
                    resultingSigma);
            } else {
                curve = scaleSpaceHelper.computeCurvature(
                    lastCurve.getXYCurve(), sigma, resultingSigma);
            }

            map.put(sigma, curve);

            sigma = SIGMA.increaseToFactorBySQRT2(resultingSigma);

            resultingSigma *= Math.sqrt(2);

            lastCurve = curve;

        }

        if (!map.containsKey(maxSIGMA)) {

            ScaleSpaceCurve curve = scaleSpaceHelper.computeCurvature(edge,
                maxSIGMA, SIGMA.getValue(maxSIGMA));

            map.put(maxSIGMA, curve);
        }

        return map;
    }

    /**
     * find the corners in the given scale space map.
     * The corners are found using the curvature minima and maxima points in
     * the curve.  A lower threshold is determined and used during the maxima
     * finding and minima comparison.  Each corner candidate is larger than one
     * of the adjacent minima by a factor such as 2 or 3.
     *
     * The returned variable is the set of (x,y) points for the candidate
     * corners found.
     *
     * @param scaleSpace scale space map for an edge
     * @param scaleSpaceSigma
     * @param edgeNumber the edgeNumber of the scaleSpace.  it's passed for
     * debugging purposes.
     * @param correctForJaggedLines
     * @param isAClosedCurve
     * @param doUseOutdoorMode
     * @return
     */
    protected CornerArray findCornersInScaleSpaceMap(
        final ScaleSpaceCurve scaleSpace, final SIGMA scaleSpaceSigma,
        int edgeNumber, boolean correctForJaggedLines,
        final boolean isAClosedCurve, final boolean doUseOutdoorMode) {

        float[] k = Arrays.copyOf(scaleSpace.getK(), scaleSpace.getK().length);

        float[] outputLowThreshold = new float[1];

        List<Integer> minimaAndMaximaIndexes = findMinimaAndMaximaInCurvature(
            k, outputLowThreshold);

        List<Integer> maxCandidateCornerIndexes = findCandidateCornerIndexes(
            k, minimaAndMaximaIndexes, outputLowThreshold[0], doUseOutdoorMode);

        CornerArray xy = new CornerArray(maxCandidateCornerIndexes.size());

        if (maxCandidateCornerIndexes.isEmpty()) {
            return xy;
        }

        if (correctForJaggedLines && !doUseOutdoorMode) {

            PairIntArray jaggedLines = removeFalseCorners(
                scaleSpace.getXYCurve(), maxCandidateCornerIndexes,
                isAClosedCurve);
        }

        int minDistFromEnds = 5;
        int nPoints = scaleSpace.getSize();
        for (int ii = 0; ii < maxCandidateCornerIndexes.size(); ii++) {

            int idx = maxCandidateCornerIndexes.get(ii);

            if (doUseOutdoorMode && !scaleSpace.curveIsClosed()) {
                if ((idx < minDistFromEnds)
                    || (idx > (nPoints - minDistFromEnds))) {
                    continue;
                }
            }

            float x = scaleSpace.getX(idx);
            float y = scaleSpace.getY(idx);

            xy.add(x, y, idx, scaleSpaceSigma);
        }

        return xy;
    }

    /**
     * determine the corners in the highest sigma of those maps and refine the
     * corner locations in the smaller sigma maps.
     *
     * @param edge
     * @param scaleSpaceCurves
     * @param edgeNumber
     * @return
     */
    private CornerArray findCornersInScaleSpaceMap(final PairIntArray edge,
        final Map<SIGMA, ScaleSpaceCurve> scaleSpaceCurves,
        final Set<PairIntWithIndex> junctions,
        final int edgeNumber, final boolean doUseOutdoorMode) {

        boolean isAClosedCurve = (edge instanceof PairIntArrayWithColor) &&
            (((PairIntArrayWithColor)edge).getColor() == 1);

        SIGMA maxSIGMA = getMaxSIGMAForECSS(edge.getN());

        ScaleSpaceCurve maxScaleSpace = scaleSpaceCurves.get(maxSIGMA);

        if (maxScaleSpace == null) {
            throw new IllegalStateException(
                "could not find the scale space for max sigma");
        }

        if (maxScaleSpace.getK() == null) {
            return new CornerArray(0);
        }

        CornerArray candidateCornersXY =
            findCornersInScaleSpaceMap(maxScaleSpace, maxSIGMA, edgeNumber,
                enableJaggedLineCorrections, isAClosedCurve, doUseOutdoorMode);

        insertJunctions(maxScaleSpace, maxSIGMA, edgeNumber, junctions,
            candidateCornersXY);

        SIGMA sigma = SIGMA.divideBySQRT2(maxSIGMA);

        SIGMA prevSigma = maxSIGMA;

        //find the corners in the higher res scale space curves
        while (sigma != null) {

            ScaleSpaceCurve scaleSpace = scaleSpaceCurves.get(sigma);

            refinePrimaryCoordinates(candidateCornersXY,
                scaleSpace, sigma, prevSigma, edgeNumber, isAClosedCurve,
                doUseOutdoorMode);

            removeRedundantPoints(candidateCornersXY, scaleSpaceCurves);

            prevSigma = sigma;

            sigma = SIGMA.divideBySQRT2(sigma);
        }

        log.log(Level.FINE, "number of corners adding ={0}",
            Integer.valueOf(candidateCornersXY.getN()));

        CornerArray edgeCorners = new CornerArray(candidateCornersXY.getN());

        //store xc and yc for the edge
        for (int ii = 0; ii < candidateCornersXY.getN(); ii++) {

            int xte = Math.round(candidateCornersXY.getX(ii));
            int yte = Math.round(candidateCornersXY.getY(ii));

            edgeCorners.add(xte, yte, candidateCornersXY.getInt(ii),
                candidateCornersXY.getSIGMA(ii));
        }

        //insertJunctions(edgeCorners, scaleSpaceCurves, junctions);

        return edgeCorners;
    }

    /**
     * find the minima and maxima of the curvature k and return the low threshold
     * used in the in-out variable outputLowThreshold and a list of the
     * indexes for the minima as negative values and the maxima as positive
     * values.
     * @param k
     * @param outputLowThreshold
     * @return
     */
    protected List<Integer> findMinimaAndMaximaInCurvature(float[] k,
        float[] outputLowThreshold) {

        if ((k == null) || (k.length == 0)) {
            return new ArrayList<Integer>();
        }

        for (int ii = 1; ii < k.length; ii++) {
            if (k[ii] < 0.0f) {
                k[ii] *= -1.0f;
            }
        }

        if (k.length < 3) {
            return new ArrayList<Integer>();
        }

        float[] kQuartiles = ImageStatisticsHelper.getQuartiles(k);

        log.fine("quartiles=" + Arrays.toString(kQuartiles));

        //float kMax = MiscMath.findMax(k);

        // determine float lowThresh
        HistogramHolder h = Histogram.calculateSturgesHistogram(
            0, 2 * kQuartiles[2], k, Errors.populateYErrorsBySqrt(k));

        if (h.getXHist().length < 3) {
            return new ArrayList<Integer>();
        }

        /*
        try {
            h.plotHistogram("curvature", 283746);
        } catch (Exception e) {}
        */

        int[] firstPeakAndMinIdx = findFirstPeakAndMinimum(h);

        if (firstPeakAndMinIdx[1] > 3) {
            firstPeakAndMinIdx[1] = firstPeakAndMinIdx[0] + 1;
        }
        if (firstPeakAndMinIdx[1] >= (h.getYHist().length >> 1)) {
            firstPeakAndMinIdx[1] = 1;
        }

        // sum intensity <= firstMinIdx and then after to compare
        long sum0 = 0;
        long sum1 = 0;
        for (int i = 0; i < h.getXHist().length; i++) {
            if (i <= firstPeakAndMinIdx[1]) {
                sum0 += h.getYHist()[i];
            } else {
                sum1 += h.getYHist()[i];
            }
        }

        if (sum1 == 0) {
            return new ArrayList<Integer>();
        }

        float divSum = (float)sum0/(float)sum1;
        outputLowThreshold[0] = 0;
        if (divSum > 10) {
            outputLowThreshold[0] = h.getXHist()[firstPeakAndMinIdx[0]];
        } else if ((firstPeakAndMinIdx[1] > 0) && (sum1 > 0)) {
            outputLowThreshold[0] = (h.getXHist()[firstPeakAndMinIdx[1]]
                + h.getXHist()[firstPeakAndMinIdx[1] - 1])/2;
        } else if ((firstPeakAndMinIdx[1] > 0) && (sum1 == 0)) {
            firstPeakAndMinIdx[1]--;
        }

        log.fine("lowThresh=" + outputLowThreshold[0]
            + " sum0=" + sum0 + " sum1=" + sum1 + " divsum=" + divSum
            + " firstPeakAndMinIdx[0]=" + firstPeakAndMinIdx[0]
            + " firstPeakAndMinIdx[1]=" + firstPeakAndMinIdx[1]);

        /*
        storing the minima and maxima in the same array list.
        the minima have -1*index within k
        and the maxima keep their positive values of the index within k.
        */
        List<Integer> minMaxIndexes = new ArrayList<Integer>();

        float lastK = k[0];
        boolean incr = true;
        for (int ii = 1; ii < k.length; ii++) {

            float currentK = k[ii];

            if ((currentK < lastK) && incr) {
                if (k[ii - 1] > outputLowThreshold[0]) {
                    minMaxIndexes.add(Integer.valueOf(ii - 1));
                }
                incr = false;
            } else if ((currentK > lastK) && !incr) {
                // values below outputLowThreshold[0] are handled by
                // callers.  TODO: redesign the caller and this method
                // to not need to understand peculiarities of the data.
                minMaxIndexes.add(Integer.valueOf(-1*(ii - 1)));
                incr = true;
            }

            lastK = currentK;
        }

        if (incr) {
            // add the last point
             minMaxIndexes.add(Integer.valueOf(k.length - 1));
        }

        return minMaxIndexes;
    }

    /**
     * given histogram h, find the first peak and the subsequent first minima
     * after it.
     *
     * @param h
     * @return
     */
    protected int[] findFirstPeakAndMinimum(HistogramHolder h) {

        float[] xh = h.getXHist();
        float[] yh = h.getYHistFloat();

        int firstMinIdx = -1;
        int yFirstPeakIdx = -1;
        float lastY = yh[0];
        boolean incr = true;
        for (int i = 1; i < xh.length; i++) {
            if (yFirstPeakIdx > -1) {
                if ((yh[i] < lastY) && incr) {
                    incr = false;
                } else if ((yh[i] > lastY) && !incr) {
                    firstMinIdx = i - 1;
                    break;
                }
            } else {
                if ((yh[i] < lastY) && incr) {
                    incr = false;
                    yFirstPeakIdx = i - 1;
                }
            }
            lastY = yh[i];
        }
        firstMinIdx++;

        return new int[]{yFirstPeakIdx, firstMinIdx};
    }

    /**
     * remove false corners from maxCandidateCornerIndexes by determining if
     * the corner is due to a jagged line.  The flag "isAClosedCurve" is
     * used to preserve corners that are on the edges of the curve because
     * the calculation of curvature has used that property to do a more
     * accurate "wrap around" calculation, so the corners there should be real.
     *
     * @param xyCurve
     * @param maxCandidateCornerIndexes indexes which are reduced by this method
     * to a smaller number due to removing false corners
     * @param isAClosedCurve
     * @return returns the found jagged lines in the edge
     */
     protected PairIntArray removeFalseCorners(PairIntArray xyCurve,
        List<Integer> maxCandidateCornerIndexes, boolean isAClosedCurve) {

        // until the methods used by findJaggedLineSegments and
        // findJaggedLineSegments2 are simplified, use both separately
        // and keep the solution which results in fewer corners.

        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();

        PairIntArray jaggedLineSegments1 =
            curveHelper.findJaggedLineSegments(xyCurve);

        PairIntArray jaggedLineSegments2 =
            curveHelper.findJaggedLineSegments2(xyCurve);

        List<Integer> remove1 = new ArrayList<Integer>();

        List<Integer> remove2 = new ArrayList<Integer>();

        for (int ii = 0; ii < maxCandidateCornerIndexes.size(); ii++) {

            int idx = maxCandidateCornerIndexes.get(ii);

            if (isAClosedCurve && ((idx < 3) || (idx > (xyCurve.getN() - 4)))) {
                // keep

            } else if ((idx < 4) || (idx > (xyCurve.getN() - 5))) {

                remove1.add(Integer.valueOf(ii));
                remove2.add(Integer.valueOf(ii));

            } else {

                //TODO: make this result less sensitive to minDistFromBoundary
                int range1Idx = isWithinARange(jaggedLineSegments1, idx,
                    3);

                if (range1Idx > -1) {

                    int last1Idx = jaggedLineSegments1.getN() - 1;

                    if ((range1Idx > 0) && (range1Idx < last1Idx)) {
                        remove1.add(Integer.valueOf(ii));
                    }
                }
                int range2Idx = isWithinARange(jaggedLineSegments2, idx, 3);
                if (range2Idx > -1) {
                    remove2.add(Integer.valueOf(ii));
                }
            }
        }
        int choose1 = 0;
        if (remove1.size() == remove2.size()) {
            // choose the lines which cover more of the curve
            int sum1 = 0;
            if (jaggedLineSegments1 != null) {
                for (int ii = 0; ii < jaggedLineSegments1.getN(); ii++) {
                    sum1 += (jaggedLineSegments1.getY(ii)
                        - jaggedLineSegments1.getX(ii));
                }
            }
            int sum2 = 0;
            if (jaggedLineSegments2 != null) {
                for (int ii = 0; ii < jaggedLineSegments2.getN(); ii++) {
                    sum2 += (jaggedLineSegments2.getY(ii)
                        - jaggedLineSegments2.getX(ii));
                }
            }
            if (sum1 > sum2) {
                choose1 = 1;
            }
        } else if (remove1.size() > remove2.size()) {
            choose1 = 1;
        }
        if (choose1 == 0) {
            for (int ii = (remove2.size() - 1); ii > -1; ii--) {
                int idx = remove2.get(ii).intValue();
                maxCandidateCornerIndexes.remove(idx);
            }
            log.fine("NREMOVED=" + remove2.size());

            return jaggedLineSegments2;

        } else {
            for (int ii = (remove1.size() - 1); ii > -1; ii--) {
                int idx = remove1.get(ii).intValue();
                maxCandidateCornerIndexes.remove(idx);
            }
            log.fine("NREMOVED=" + remove1.size());

            return jaggedLineSegments1;
        }
    }

    /**
     * search for idx within ranges in lineRangeSegments and return the index of
     * lineRangeSegments in which it is found, else -1.  Note that
     * lineRangeSegments have to be ordered by x and unique.
     * @param lineRangeSegments
     * @param idx
     * @param minDistFromEnds
     * @return
     */
    private int isWithinARange(PairIntArray lineRangeSegments, int idx,
        int minDistFromEnds) {

        if ((lineRangeSegments == null) || (lineRangeSegments.getN() == 0)) {
            return -1;
        }

        // search outside of bounds first:
        if (idx < lineRangeSegments.getX(0)) {
            return -1;
        } else if (idx > lineRangeSegments.getY(lineRangeSegments.getN() - 1)) {
            return -1;
        }

        for (int i = 0; i < lineRangeSegments.getN(); i++) {

            int idx0 = lineRangeSegments.getX(i);
            int idx1 = lineRangeSegments.getY(i);

            if ((idx >= (idx0 + minDistFromEnds))
                && (idx <= (idx1 - minDistFromEnds))) {

                return i;
            }
        }

        return -1;
    }

     /**
     * refine the primary coordinates given in xc and yc using the corner
     * candidate results of the given scale space map.  Matches to a corner
     * are found using a rough separation limit determined by the previous
     * scale map's sigma.
     *
     * @param xc
     * @param yc
     * @param scaleSpace
     * @param previousSigma
     * @param edgeNumber included for debugging
     */
    private void refinePrimaryCoordinates(CornerArray candidateCornersXY,
        ScaleSpaceCurve scaleSpace, final SIGMA scaleSpaceSigma,
        final SIGMA previousSigma, final int edgeNumber,
        final boolean isAClosedCurve, final boolean doUseOutdoorMode) {

        if (scaleSpace == null || scaleSpace.getK() == null) {
            //TODO: follow up on NPE here
            return;
        }

        CornerArray xy2 =
            findCornersInScaleSpaceMap(scaleSpace, scaleSpaceSigma, edgeNumber,
                false, isAClosedCurve, doUseOutdoorMode);

        // roughly estimating maxSep as the ~FWZI of the gaussian
        //TODO: this may need to be altered to a smaller value
        float maxSepSq = Gaussian1D.estimateHWZI(previousSigma, 0.01f);
        maxSepSq *= maxSepSq;
        if (maxSepSq > 4) {
            maxSepSq = 4;
        }
        float maxSep = (float)Math.sqrt(maxSepSq);

        NearestPointsFloat np = new NearestPointsFloat(xy2.getX(), xy2.getY(),
            xy2.getN());

        // revise the points in {xc, yc} to the closest in {xc2, yc2}
        for (int j = 0; j < candidateCornersXY.getN(); j++) {
            float x = candidateCornersXY.getX(j);
            float y = candidateCornersXY.getY(j);

            Integer minSepIndex = np.findClosestNeighborIndex(x, y, maxSep);

            if (minSepIndex != null) {
                int minSepIdx = minSepIndex.intValue();
                float x3 = xy2.getX(minSepIdx);
                float y3 = xy2.getY(minSepIdx);
                candidateCornersXY.set(j, x3, y3,
                    xy2.getInt(minSepIdx), xy2.getSIGMA(minSepIdx));
            }
        }
    }

    private boolean storeCornerRegion(int edgeNumber, int cornerIdx,
        ScaleSpaceCurve scaleSpace) {

        final float[] k = scaleSpace.getK();

        //for 2 neighboring points on each side, min k is 0.2
        float kCenterAbs = Math.abs(k[cornerIdx]);
        if (kCenterAbs < 0.14f) {//0.18
            return false;
        }

        int n = scaleSpace.getSize();

        if (n < 3) {
            return false;
        }

        boolean isClosedCurve = scaleSpace.curveIsClosed();

        int nCR = 0;
        int kMaxIdx = -1;

        int count = 0;
        int nH = 2;
        if (scaleSpace.getSize() < 5) {
            nH = 1;
        }
        int[] pIdxs = new int[(2*nH) + 1];
        for (int pIdx = (cornerIdx - nH); pIdx <= (cornerIdx + nH); ++pIdx) {
            if ((pIdx < 0) && isClosedCurve) {
                pIdxs[count] = n + pIdx;
            } else if ((pIdx > (n - 1)) && isClosedCurve) {
                pIdxs[count] = pIdx - n;
            } else {
                pIdxs[count] = pIdx;
            }
            count++;
        }

        for (int pIdx : pIdxs) {
            if (pIdx < 0 || (pIdx > (n - 1))) {
                // it's out of bounds only if it is not a closed curve
                continue;
            }
            if (pIdx == cornerIdx) {
                kMaxIdx = nCR;
            }
            int x = Math.round(scaleSpace.getX(pIdx));
            int y = Math.round(scaleSpace.getY(pIdx));

            // discard if out of bounds
            if ((x < 0) || (y < 0) || (x > (width - 1)) || (y > (height - 1))) {
                return false;
            }
            nCR++;
        }
        if (nCR < 3) {
            return false;
        }

        if (!isClosedCurve) {
            if (kMaxIdx == 0 || kMaxIdx == (nCR - 1)) {
                return false;
            }
        }

        CornerRegion cr = new CornerRegion(edgeNumber, nCR, kMaxIdx);
        nCR = 0;
        for (int pIdx : pIdxs) {
            if (pIdx < 0 || (pIdx > (n - 1))) {
                // it's out of bounds only if it is not a closed curve
                continue;
            }

            int x = Math.round(scaleSpace.getX(pIdx));
            int y = Math.round(scaleSpace.getY(pIdx));

            cr.set(nCR, k[pIdx], x, y);

            if (kMaxIdx == nCR) {
                cr.setIndexWithinCurve(pIdx);
            }

            nCR++;
        }
        Integer key = Integer.valueOf(edgeNumber);
        List<CornerRegion> list = edgeCornerRegionMap.get(key);
        if (list == null) {
            list = new ArrayList<CornerRegion>();
            edgeCornerRegionMap.put(key, list);
        }
        list.add(cr);

        return true;
    }

    /**
     * maxSigma is defined by the ECSS algorithm in:
     * 2006, "Performance evaluation of corner detectors using consistency and
     * accuracy measures" by Farzin Mokhtarian and Farahnaz Mohanna in
     * Computer Vision and Image Understanding, vol 102, pp 81-94.
     * @param nPoints
     * @return
     */
    protected SIGMA getMaxSIGMAForECSS(int nPoints) {

        //ECSS:
        //    < 200, sigma = 2
        //    <= 600, sigma=3
        //    else    sigma=4
        if (nPoints < 200) {
            return SIGMA.TWO;
        } else if (nPoints < 601) {
            return SIGMA.THREE;
        } else {
            return SIGMA.FOUR;
        }
    }

    /**
     * given the curvature array and a list of the indexes of the
     * minima and maxima of the curvature array, find the candidate corner
     * indexes with respect to the k array.
     *
     * @param k
     * @param minMaxIndexes
     * @param lowThreshold
     * @param doUseOutdoorMode
     * @return
     */
    protected List<Integer> findCandidateCornerIndexes(float[] k,
        List<Integer> minMaxIndexes, float lowThreshold,
        final boolean doUseOutdoorMode) {

        // find peaks where k[ii] is > factorAboveMin* adjacent local minima

        float factorAboveMin = factorIncreaseForCurvatureMinimum * 2.5f;//3.5f;// 10 misses some corners

if (doUseOutdoorMode) {
    factorAboveMin = factorIncreaseForCurvatureMinimum * 3.5f;//10.f;
}
        log.fine("using factorAboveMin=" + factorAboveMin);

        //to limit k to curvature that shows a rise in 1 pixel over a run of 3,
        // use 0.2 for a lower limit.
        // TODO: it's not clear that kLowerLimit is a good idea.  the relative change
        // filter alone is good for all size scale corners, and adding this
        // limit biases the results.  may want to only use this bias if
        // the some amount of curvature points are >= 0.2
        float kLowerLimit = 0.05f;

        List<Integer> cornerCandidates = new ArrayList<Integer>();

        // choose candidates from minMaxIndexes that are
        //     >= factorAboveMin one adjacent min
        for (int ii = 0; ii < minMaxIndexes.size(); ii++) {

            int idx = minMaxIndexes.get(ii).intValue();

            if (idx > -1) {
                // this is maximum

                boolean found = false;

                // compare to preceding minimum
                for (int iii = (ii - 1); iii > -1; iii--) {
                    int idx2 = minMaxIndexes.get(iii).intValue();
                    if (idx2 < 0) {
                        float compare = k[-1*idx2];
                        if (compare < lowThreshold) {
                            // avoids divide by very small number sometimes
                            compare = lowThreshold;
                        }
                        if (k[idx] >= kLowerLimit && k[idx] >= factorAboveMin * compare) {
                            cornerCandidates.add(Integer.valueOf(idx));
                            found = true;
                        }
                        break;
                    }
                }
                if (found) {
                    continue;
                }

                //compare to proceeding minimum
                for (int iii = (ii + 1); iii < minMaxIndexes.size(); iii++) {
                    int idx2 = minMaxIndexes.get(iii).intValue();
                    if (idx2 < 0) {
                        float compare = k[-1*idx2];
                        if (compare < lowThreshold) {
                            // avoids divide by very small number sometimes
                            compare = lowThreshold;
                        }
                        if (k[idx] >= kLowerLimit && k[idx] >= factorAboveMin * compare) {
                            cornerCandidates.add(Integer.valueOf(idx));
                        }
                        break;
                    }
                }
            }
        }

        return cornerCandidates;
    }

    /**
     * @return the edgeCornerRegionMap
     */
    public Map<Integer, List<CornerRegion>> getEdgeCornerRegionMap() {
        return edgeCornerRegionMap;
    }

    private void removeRedundantPoints(CornerArray candidateCornersXY,
        Map<SIGMA, ScaleSpaceCurve> scaleSpaceCurves) {

        boolean hasRedundant = false;

        // looking for redundant points
        Map<PairInt, Set<Integer>> coordIndexes = new HashMap<PairInt, Set<Integer>>();
        for (int j = 0; j < candidateCornersXY.getN(); j++) {

            int x = Math.round(candidateCornersXY.getX(j));
            int y = Math.round(candidateCornersXY.getY(j));
            PairInt p = new PairInt(x, y);

            Set<Integer> indexes = coordIndexes.get(p);
            if (indexes == null) {
                indexes = new HashSet<Integer>();
                coordIndexes.put(p, indexes);
            } else {
                hasRedundant = true;
            }
            indexes.add(Integer.valueOf(j));
        }

        if (!hasRedundant) {
            return;
        }

        // remove the weaker 'k' of any redundant points
        Set<Integer> resolved = new HashSet<Integer>();
        List<Integer> remove = new ArrayList<Integer>();
        for (Entry<PairInt, Set<Integer>> entry : coordIndexes.entrySet()) {
            Set<Integer> set = entry.getValue();
            if (set.size() > 1) {
                float maxAbsK = Float.MIN_VALUE;
                Integer bestIndex1 = null;
                for (Integer index1 : set) {
                    if (resolved.contains(index1)) {
                        continue;
                    }
                    int idx = index1.intValue();
                    SIGMA sigma = candidateCornersXY.getSIGMA(idx);
                    ScaleSpaceCurve ssc = scaleSpaceCurves.get(sigma);
                    int sscIdx = candidateCornersXY.getInt(idx);
                    assert(ssc != null);
                    float kAbs = Math.abs(ssc.getK()[sscIdx]);
                    if (kAbs > maxAbsK) {
                        maxAbsK = kAbs;
                        bestIndex1 = index1;
                    }
                    resolved.add(index1);
                }
                if (bestIndex1 != null) {
                    for (Integer index1 : set) {
                        if (index1.equals(bestIndex1)) {
                            continue;
                        }
                        remove.add(index1);
                    }
                }
            }
        }

        Collections.sort(remove);
        for (int i = (remove.size() - 1); i > -1; --i) {
            int idx = remove.get(i);
            candidateCornersXY.removeRange(idx, idx);
        }
    }

    private void insertJunctions(ScaleSpaceCurve scaleSpace, SIGMA sigma,
        int edgeNumber, Set<PairIntWithIndex> junctions,
        CornerArray candidateCornersXY) {

        if (junctions == null) {
            return;
        }

        for (PairIntWithIndex p : junctions) {

            int idxWithinEdge = p.getPixIndex();

            // find where to insert in candidateCornersXY
            int idx = Arrays.binarySearch(candidateCornersXY.getYInt(),
                0, candidateCornersXY.getN(), idxWithinEdge);
            // if it's negative, (-(insertion point) - 1)
            if (idx < 0) {
                // idx = -*idx2 - 1
                idx = -1 * (idx + 1);
            }

            float x = scaleSpace.getX(idxWithinEdge);
            float y = scaleSpace.getY(idxWithinEdge);

            candidateCornersXY.insert(idx, x, y, idxWithinEdge, sigma);
        }
    }

}
