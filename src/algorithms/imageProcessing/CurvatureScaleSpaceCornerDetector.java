package algorithms.imageProcessing;

import algorithms.MultiArrayMergeSort;
import algorithms.compGeometry.NearestPoints;
import algorithms.util.PairIntArray;
import algorithms.util.PairFloatArray;
import algorithms.util.PairIntArrayWithColor;
import algorithms.misc.Histogram;
import algorithms.misc.HistogramHolder;
import algorithms.util.Errors;
import algorithms.util.PairInt;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.logging.Level;

/**
 * The code is implemented from interpreting several papers by the authors
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
 * possible.
 *
 * @author nichole
 */
public class CurvatureScaleSpaceCornerDetector extends
    AbstractCurvatureScaleSpaceMapper {

    /**
     * corners detected in the image.  the false corners have been removed
     */
    protected PairIntArray corners = new PairIntArray();

    /**
     * this is not ready for use yet.  when implemented it should hold
     * a sublist of corners that are better to use for matching the same
     * in other images.
     * TODO: populate this with edgeCornerRegionMap values
     */
    protected PairIntArray cornersForMatching = new PairIntArray();

    protected Map<Integer, List<CornerRegion>> edgeCornerRegionMap = new
        HashMap<Integer, List<CornerRegion>>();
    
    /**
     * corners populated if extractSkyline is true
     */
    protected PairIntArray skylineCorners = new PairIntArray();

    protected boolean enableJaggedLineCorrections = true;
    
    protected float factorIncreaseForCurvatureMinimum = 1.f;
    
    protected boolean performWholeImageCorners = true;
    
    protected boolean doStoreCornerRegions = true;
            
    public CurvatureScaleSpaceCornerDetector(final ImageExt input) {

        super(input);
    }
    
    public CurvatureScaleSpaceCornerDetector(final GreyscaleImage input) {

        super(input);
    }

    public CurvatureScaleSpaceCornerDetector(final ImageExt input,
        List<PairIntArray> theEdges) {

        super(input, theEdges);
    }
    
    /**
     * set the edge detector to create edges that are better for outdoor
     * conditions and calculate corners only for the skyline.  
     * Note that the skyline extraction is currently
     * a long running process.
     */
    void calculateSkylineCornersOnly() {
        
        useOutdoorModeAndExtractSkyline();
        
        performWholeImageCorners = false;
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
    
    public void findCorners() {

        initialize();

        if (extractSkyline && !skylineEdges.isEmpty()) {
            calculateSkylineCorners();
        }

        if (!performWholeImageCorners) {
            return;
        }
        
        // not re-using return maps for now, but they are available here
        // while refactoring the public method returns and signatures
        Map<PairIntArray, Map<SIGMA, ScaleSpaceCurve> > maps =
            findCornersInScaleSpaceMaps(edges, useOutdoorMode, corners);
        
        includeJunctionsInCorners();
    }
    
    @Override
    protected void reinitializeSpecialization() {
        corners = new PairIntArray();
    }

    /**
     * find corners iteratively until approximately the number desired are
     * found.  This method is useful for creating corners in stereo projection
     * images - it helps to adjust the image intensity levels so that
     * similar edges can be formed in both images.
     * Presumably, if calibration of the images were possible, findCorners()
     * should alone provide a stable result (where calibration of the images
     * are steps such as bias removal, flat field corrections, intensity
     * corrections using standard candles, and corrections to make the
     * point spread functions similar for the images).
     *
     * @param approxNumberOfCornersDesired
     * @param filterOnlyAboveThisNumberOfCorners if the default number of corners
     * produced is this larger or larger, the method will iteratively
     * increase the lower intensity filter until approxNumberOfCornersDesired
     * are produced, else if the default number of corners is less
     * than useOnlyAboveThisNumberOfCorners, the method will not filter
     * the image further.
     */
    public void findCornersIteratively(int approxNumberOfCornersDesired,
        int filterOnlyAboveThisNumberOfCorners) {

        //TODO: this method needs to be refactored
        
        if (!performWholeImageCorners) {
            throw new IllegalStateException(
                "performWholeImageCorners is currently set to false");
        }

        float lowerThresholdStepSize = 1.0f;
        
        int nCorners = corners.getN();

        float lowThreshold = (useOutdoorMode) ?
            CannyEdgeFilter.defaultOutdoorLowThreshold :
            CannyEdgeFilter.defaultLowThreshold;
        
        if ((nCorners > 0) && (nCorners < filterOnlyAboveThisNumberOfCorners)) {
            return;
        } else if ((nCorners > 0) && (nCorners < approxNumberOfCornersDesired)) {
            return;
        } else if (state.ordinal() < CurvatureScaleSpaceMapperState.INITIALIZED.ordinal()) {
            findCorners();
        }

        nCorners = corners.getN();

        if (nCorners < filterOnlyAboveThisNumberOfCorners) {
            return;
        }

        //TODO: this could be improved to assert a minimum presence of corners
        // at boundaries and throughout the image compared to the first round.
        // In other words, would not want to remove all corners for an important
        // part of the image intersection with another image.

        List<PairIntArray> prevEdges = copy(this.edges);
        PairIntArray prevCorners = corners.copy();

        int nIter = 0;
            
        while ((nCorners > 0) && (nCorners > approxNumberOfCornersDesired)) {

            log.info("nCorners=" + nCorners);

            prevEdges = copy(this.edges);
            prevCorners = this.corners.copy();

            //TODO: needs adjustments:
            float additionalBlurSigma = 0;
            if (nIter == 0) {
                //additionalBlurSigma = 1;
                /*if (nCorners > 1000) {
                    lowerThresholdStepSize  = 1.5f;
                } else {*/
                    lowerThresholdStepSize = 0.5f;
                //}
            }
            
            lowThreshold += lowerThresholdStepSize;
            
            reinitialize(lowThreshold, additionalBlurSigma);

            findCornersInScaleSpaceMaps(edges, useOutdoorMode, corners);

            includeJunctionsInCorners();
            
            nCorners = corners.getN();
            
            nIter++;
        }

        if (Math.abs(corners.getN() - approxNumberOfCornersDesired) >
            Math.abs(prevCorners.getN() - approxNumberOfCornersDesired)) {
            
            this.edges = prevEdges;
            this.corners = prevCorners;
        }
        
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
     * Find corners in the edges by creating scale space maps for each edge that
     * range from a maximum value determined in getMaxSIGMAForECSS() and
     * determine the corners in the highest sigma of those maps and refine the
     * corner locations in the smaller sigma maps.  The corners are found using
     * the curvature minima and maxima points in the curve.  A lower threshold
     * is determined and used during the maxima finding and minima comparison.
     * Each corner candidate is larger than one of the adjacent minima by
     * a factor such as 2 or 3.
     * The results are set in the instance variable corners.
     * The returned variable is the set of scale space maps which might be
     * useful for other purposes, but are no longer needed for the corner
     * determination.
     *
     * @param theEdges
     * @param doUseOutdoorMode
     * @param outputCorners
     * @return scale space maps for each edge
     */
    protected Map<PairIntArray, Map<SIGMA, ScaleSpaceCurve> >
    findCornersInScaleSpaceMaps(final List<PairIntArray> theEdges, final boolean
        doUseOutdoorMode, final PairIntArray outputCorners) {

        Map<PairIntArray, Map<SIGMA, ScaleSpaceCurve> > scaleSpaceMaps
            = new HashMap<PairIntArray, Map<SIGMA, ScaleSpaceCurve> >();
       
        // perform the analysis on the edgesScaleSpaceMaps
        for (int i = 0; i < theEdges.size(); i++) {

            final PairIntArray edge = theEdges.get(i);

            final Map<SIGMA, ScaleSpaceCurve> map =
                createLowUpperThresholdScaleSpaceMaps(edge);

            scaleSpaceMaps.put(edge, map);

            PairIntArray edgeCorners = findCornersInScaleSpaceMap(edge, map, i,
                doUseOutdoorMode);

            log.log(Level.FINE,
                "{0}) number of corners adding ={1} for edge={2}",
                new Object[]{Integer.valueOf(i),
                    Integer.valueOf(edgeCorners.getN()),
                    Integer.valueOf(i)});

            //store xc and yc for the edge
            for (int ii = 0; ii < edgeCorners.getN(); ii++) {
                outputCorners.add(edgeCorners.getX(ii), edgeCorners.getY(ii));
            }
        }

        return scaleSpaceMaps;
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
     * @param edgeNumber the edgeNumber of the scaleSpace.  it's passed for
     * debugging purposes.
     * @param correctForJaggedLines
     * @param isAClosedCurve
     * @param doUseOutdoorMode
     * @return
     */
    protected PairFloatArray findCornersInScaleSpaceMap(
        final ScaleSpaceCurve scaleSpace, int edgeNumber,
        boolean correctForJaggedLines, final boolean isAClosedCurve,
        final boolean doUseOutdoorMode) {

        float[] k = Arrays.copyOf(scaleSpace.getK(), scaleSpace.getK().length);

        float[] outputLowThreshold = new float[1];

        List<Integer> minimaAndMaximaIndexes = findMinimaAndMaximaInCurvature(
            k, outputLowThreshold);

        List<Integer> maxCandidateCornerIndexes = findCandidateCornerIndexes(
            k, minimaAndMaximaIndexes, outputLowThreshold[0], doUseOutdoorMode);

        PairFloatArray xy = new PairFloatArray(maxCandidateCornerIndexes.size());

        if (maxCandidateCornerIndexes.isEmpty()) {
            return xy;
        }
        
        boolean storeCornerRegion = doStoreCornerRegions && 
            !edgeCornerRegionMap.containsKey(Integer.valueOf(edgeNumber));
        
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
            
            if (storeCornerRegion) {
                storeCornerRegion(edgeNumber, idx, k, scaleSpace);
            }
            
            xy.add(x, y);
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
    private PairIntArray findCornersInScaleSpaceMap(final PairIntArray edge,
        final Map<SIGMA, ScaleSpaceCurve> scaleSpaceCurves,
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
            return new PairIntArray(0);
        }
        
        PairFloatArray candidateCornersXY =
            findCornersInScaleSpaceMap(maxScaleSpace, edgeNumber, 
                enableJaggedLineCorrections, isAClosedCurve, doUseOutdoorMode);

        SIGMA sigma = SIGMA.divideBySQRT2(maxSIGMA);

        SIGMA prevSigma = maxSIGMA;

        //find the corners in the higher res scale space curves
        while (sigma != null) {

            ScaleSpaceCurve scaleSpace = scaleSpaceCurves.get(sigma);

            refinePrimaryCoordinates(candidateCornersXY.getX(),
                candidateCornersXY.getY(), candidateCornersXY.getN(),
                scaleSpace, prevSigma, edgeNumber, isAClosedCurve,
                doUseOutdoorMode);

            prevSigma = sigma;

            sigma = SIGMA.divideBySQRT2(sigma);
        }

        log.log(Level.FINE, "number of corners adding ={0}",
            Integer.valueOf(candidateCornersXY.getN()));

        PairIntArray edgeCorners = new PairIntArray(candidateCornersXY.getN());

        //store xc and yc for the edge
        for (int ii = 0; ii < candidateCornersXY.getN(); ii++) {
            int xte = Math.round(candidateCornersXY.getX(ii));
            int yte = Math.round(candidateCornersXY.getY(ii));
            edgeCorners.add(xte, yte);
        }

        return edgeCorners;
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
    private void refinePrimaryCoordinates(float[] xc, float[] yc,
        final int xyLength,
        ScaleSpaceCurve scaleSpace, final SIGMA previousSigma,
        final int edgeNumber, final boolean isAClosedCurve,
        final boolean doUseOutdoorMode) {

        if (scaleSpace == null || scaleSpace.getK() == null) {
            //TODO: follow up on NPE here
            return;
        }

        PairFloatArray xy2 =
            findCornersInScaleSpaceMap(scaleSpace, edgeNumber, false,
            isAClosedCurve, doUseOutdoorMode);

        // roughly estimating maxSep as the ~FWZI of the gaussian
        //TODO: this may need to be altered to a smaller value
        float maxSepSq = Gaussian1D.estimateHWZI(previousSigma, 0.01f);
        maxSepSq *= maxSepSq;

        if (maxSepSq > 4) {
            maxSepSq = 4;
        }

        // revise the points in {xc, yc} to the closest in {xc2, yc2}
        boolean[] matchedNew = new boolean[xy2.getN()];
        for (int j = 0; j < xyLength; j++) {
            float x = xc[j];
            float y = yc[j];
            float minSepSq = Float.MAX_VALUE;
            int minSepIdx = -1;
            for (int jj = 0; jj < xy2.getN(); jj++) {
                if (matchedNew[jj]) {
                    continue;
                }
                float x2 = xy2.getX(jj);
                float y2 = xy2.getY(jj);
                float dx = x2 - x;
                float dy = y2 - y;
                float sepSq = (dx*dx) + (dy*dy);
                if (sepSq < minSepSq) {
                    minSepSq = sepSq;
                    minSepIdx = jj;
                }
            }

            if (minSepIdx > -1) {
                if (minSepSq < maxSepSq) {
                    matchedNew[minSepIdx] = true;
                    xc[j] = xy2.getX(minSepIdx);
                    yc[j] = xy2.getY(minSepIdx);
                }
            }
        }

    }

    public PairIntArray getCorners() {
        return corners;
    }

    public PairIntArray getCornersForMatching() {
        return cornersForMatching;
    }

    /**
     * <em>note, this is not ready for use yet.</em>
     * it's meant to be a sublist of
     * the variable "corners" selected to be better for matching the same
     * corners in other images.
     * @return
     */
    public PairIntArray getCornersForMatchingInOriginalReferenceFrame() {

        PairIntArray co = new PairIntArray();
        /*for (int i = 0; i < cornersForMatching.getN(); i++) {
            int x = cornersForMatching.getX(i);
            int y = cornersForMatching.getY(i);
            x += this.trimmedXOffset;
            y += this.trimmedYOffset;
            co.add(x, y);
        }*/

        return co;
    }

    public PairIntArray getSkylineCornersInOriginalReferenceFrame() {

        PairIntArray co = new PairIntArray();

        for (int i = 0; i < skylineCorners.getN(); i++) {
            int x = skylineCorners.getX(i);
            int y = skylineCorners.getY(i);
            x += this.trimmedXOffset;
            y += this.trimmedYOffset;
            co.add(x, y);
        }

        return co;
    }

    public PairIntArray getCornersInOriginalReferenceFrame() {
        
PairIntArray cp = corners.copy();
MultiArrayMergeSort.sortByYThenX(cp);
  
        PairIntArray co = new PairIntArray();
        for (int i = 0; i < corners.getN(); i++) {
            int x = corners.getX(i);
            int y = corners.getY(i);
            x += this.trimmedXOffset;
            y += this.trimmedYOffset;
            co.add(x, y);
        }

        return co;
    }

    public List<PairIntArray> getEdgesInOriginalReferenceFrame() {

        List<PairIntArray> output = new ArrayList<PairIntArray>();

        for (int i = 0; i < edges.size(); i++) {

            PairIntArray ce = new PairIntArray();

            PairIntArray edge = edges.get(i);

            for (int j = 0; j < edge.getN(); j++) {
                int x = edge.getX(j);
                int y = edge.getY(j);
                x += this.trimmedXOffset;
                y += this.trimmedYOffset;
                ce.add(x, y);
            }

            output.add(ce);
        }

        return output;
    }

    public List<PairIntArray> getSkylineEdgesInOriginalReferenceFrame() {

        List<PairIntArray> output = new ArrayList<PairIntArray>();

        for (int i = 0; i < skylineEdges.size(); i++) {

            PairIntArray ce = new PairIntArray();

            PairIntArray edge = skylineEdges.get(i);

            for (int j = 0; j < edge.getN(); j++) {
                int x = edge.getX(j);
                int y = edge.getY(j);
                x += this.trimmedXOffset;
                y += this.trimmedYOffset;
                ce.add(x, y);
            }

            output.add(ce);
        }

        return output;
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

    private List<PairIntArray> copy(List<PairIntArray> edges) {
        List<PairIntArray> copied = new ArrayList<PairIntArray>();
        for (PairIntArray edge : edges) {
            copied.add(edge.copy());
        }
        return copied;
    }

    private void includeJunctionsInCorners() {

        int nTot = corners.getN();

        int[] x = new int[nTot];
        int[] y = new int[nTot];

        System.arraycopy(corners.getX(), 0, x, 0, corners.getN());
        System.arraycopy(corners.getY(), 0, y, 0, corners.getN());

        Set<PairInt> corners2 = new HashSet<PairInt>();
        for (int i = 0; i < corners.getN(); i++) {
            int x2 = corners.getX(i);
            int y2 = corners.getY(i);
            corners2.add(new PairInt(x2, y2));
        }

        float sepDist = 4;

        final NearestPoints nearestPoints = new NearestPoints(x, y);

        for (Entry<Integer, Set<Integer>> entry : junctionMap.entrySet()) {

            int pixIdx = entry.getKey().intValue();

            int xP = img.getCol(pixIdx);
            int yP = img.getRow(pixIdx);

            Set<PairInt> overlapping = nearestPoints.findNeighbors(xP, yP,
                sepDist);

            for (PairInt p : overlapping) {
                corners2.remove(p);
            }

            corners2.add(new PairInt(xP, yP));
        }

        corners = new PairIntArray();

        for (PairInt p : corners2) {
            corners.add(p.getX(), p.getY());
        }

    }

    //NOT READY FOR USE YET
    private void calculateSkylineCorners() {

        if (!extractSkyline) {
            return;
        }

        PairIntArray theSkylineCorners = new PairIntArray();
    
        enableJaggedLineCorrections = false;
        
        Map<PairIntArray, Map<SIGMA, ScaleSpaceCurve> > maps =
            findCornersInScaleSpaceMaps(skylineEdges, false, theSkylineCorners);

        // TODO: improve this for resolution
        
        // smoothing does not make fewer points, so using min curvature factor
        
        float nGoalCorners = (gradientXY.getNPixels() < 500000) ? 50.f : 
            100.f;
        
        // want about 50 points
        float factor = (float)theSkylineCorners.getN()/nGoalCorners;
        
        if (factor > 1.2) {
            
            //TODO: consider saving the space curves to make recalc faster
            
            increaseFactorForCurvatureMinimum(4.0f * factor);
            
            PairIntArray theSkylineCorners2 = new PairIntArray();
            
            maps = findCornersInScaleSpaceMaps(skylineEdges, false, 
                theSkylineCorners2);
            
            log.info("before curvature factor change, number of skyline corners=" 
                + theSkylineCorners.getN());
            log.info("after=" + theSkylineCorners2.getN());
            
            theSkylineCorners = theSkylineCorners2;
            
            resetFactorForCurvatureMinimum();
        }
        
        if (theSkylineCorners.getN() > 0) {
            skylineCorners.addAll(theSkylineCorners);
        }
        
        enableJaggedLineCorrections = true;
        
        log.info("number of skyline corners=" + theSkylineCorners.getN());
    }

    private void storeCornerRegion(int edgeNumber, int cornerIdx, float[] k, 
        ScaleSpaceCurve scaleSpace) {
        
        //for 2 neighboring points on each side, min k is 0.2
        if (k[cornerIdx] < 0.2f) {
            return;
        }

        int nCR = 0;
        int kMaxIdx = -1;

        //log.info("edgeIdx=" + edgeNumber + " corner idx=" + cornerIdx + " (" + 
        //    Math.round(scaleSpace.getX(cornerIdx)) + "," +
        //    Math.round(scaleSpace.getY(cornerIdx)) +")");
        
        for (int pIdx = (cornerIdx - 2); pIdx <= (cornerIdx + 2); ++pIdx) {
            if (pIdx < 0 || (pIdx > (k.length - 1))) {
                continue;
            }
            if (pIdx == cornerIdx) {
                kMaxIdx = nCR;
            }
            nCR++;
            //log.info(String.format(
            //   "k[%d]=%.2f  for (%.0f, %.0f)", pIdx, k[pIdx],
            //   scaleSpace.getX(pIdx), scaleSpace.getY(pIdx))); 
        }
        if (nCR < 3) {
            return;
        }
        if (kMaxIdx == 0 || kMaxIdx == (nCR - 1)) {
            return;
        }
        //log.info(" ==> kept");
        
        CornerRegion cr = new CornerRegion(edgeNumber, nCR, kMaxIdx);
        nCR = 0;
        for (int pIdx = (cornerIdx - 2); pIdx <= (cornerIdx + 2); ++pIdx) {
            if (pIdx < 0 || (pIdx > (k.length - 1))) {
                continue;
            }
            cr.set(nCR, k[pIdx], 
                Math.round(scaleSpace.getX(pIdx)),
                Math.round(scaleSpace.getY(pIdx)));
            nCR++; 
        }
        Integer key = Integer.valueOf(edgeNumber);
        List<CornerRegion> list = edgeCornerRegionMap.get(key);
        if (list == null) {
            list = new ArrayList<CornerRegion>();
            edgeCornerRegionMap.put(key, list);
        }
        list.add(cr);
    }

    public Map<Integer, List<CornerRegion>> getEdgeCornerRegionMap() {
        return edgeCornerRegionMap;
    }
    
    public Set<CornerRegion> getEdgeCornerRegions() {
        
        Set<CornerRegion> set = new HashSet<CornerRegion>();
        
        for (Entry<Integer, List<CornerRegion>> entry : edgeCornerRegionMap.entrySet()) {
            set.addAll(entry.getValue());
        }
        
        return set;
    }
    
     public Set<CornerRegion> getEdgeCornerRegionsInOriginalReferenceFrame() {
        
        Set<CornerRegion> set = new HashSet<CornerRegion>();
        
        for (Entry<Integer, List<CornerRegion>> entry : edgeCornerRegionMap.entrySet()) {
            for (CornerRegion cr : entry.getValue()) {
                CornerRegion crCopy = cr.copy();
                for (int i = 0; i < cr.getX().length; ++i) {
                    int x = cr.getX()[i] + this.trimmedXOffset;
                    int y = cr.getY()[i] + this.trimmedYOffset;
                    crCopy.set(i, cr.getK()[i], x, y);
                }
                set.add(crCopy);
            }
        }
        
        return set;
    }
     
    public Set<CornerRegion> getEdgeCornerRegionsInOriginalReferenceFrame(
        boolean removeAmbiguousPeaks) {
        
        if (!removeAmbiguousPeaks) {
            return getEdgeCornerRegionsInOriginalReferenceFrame();
        }

        Set<CornerRegion> set = getEdgeCornerRegions(removeAmbiguousPeaks);
       
        Set<CornerRegion> edited = new HashSet<CornerRegion>();
        
        for (CornerRegion cr : set) {
            CornerRegion crCopy = cr.copy();
            for (int i = 0; i < crCopy.getX().length; ++i) {
                int x = crCopy.getX()[i] + this.trimmedXOffset;
                int y = crCopy.getY()[i] + this.trimmedYOffset;
                crCopy.set(i, crCopy.getK()[i], x, y);
            }
            edited.add(crCopy);
        }
        
        return edited;
    }
    
    public Set<CornerRegion> getEdgeCornerRegions(boolean removeAmbiguousPeaks) {
        
        if (!removeAmbiguousPeaks) {
            return getEdgeCornerRegions();
        }
        
        Map<Integer, Set<Integer> > theJunctionMap = new HashMap<Integer, Set<Integer>>();

        Map<Integer, PairInt> theJunctionLocationMap = new HashMap<Integer, PairInt>();
        
        EdgeExtractorWithJunctions.findJunctions(edges, 
            theJunctionMap, theJunctionLocationMap, img.getWidth(), img.getHeight());

        this.junctionLocationMap = theJunctionLocationMap;
        this.junctionMap = theJunctionMap;
        
        Set<CornerRegion> set = new HashSet<CornerRegion>();
       
        for (Entry<Integer, List<CornerRegion>> entry : edgeCornerRegionMap.entrySet()) {
            for (CornerRegion cr : entry.getValue()) {
                int kMaxIdx = cr.getKMaxIdx();
                float kMax = cr.getK()[kMaxIdx];
                boolean keep = true;               
                /*if ((kMaxIdx > 0) && (kMaxIdx < (cr.getX().length - 1))) {
                    if (Math.abs(cr.getK()[kMaxIdx - 1] - kMax) < 0.005) {
                        keep = false;
                    } else if (Math.abs(cr.getK()[kMaxIdx + 1] - kMax) < 0.005) {
                        keep = false;
                    }
                }*/
                
                if (keep) {
                    
                    int xCorner = cr.getX()[kMaxIdx];
                    int yCorner = cr.getY()[kMaxIdx];
                    
                    /*
                    check if in a junction, and if so either try to reconstruct
                    best path around portion of image that curvature maximum
                    is on, or discard this corner region as possibly ambiguous.
                    */
                    
                    boolean isInAJunction = isInAJunction(xCorner, yCorner);
                    
                    if (!isInAJunction) {
                        
                        set.add(cr);
                        
                    } else {
                        
                        List<CornerRegion> crWithinJunctions = 
                            searchForCornerRegionWithinJunctions(xCorner, 
                            yCorner, kMax);

                        if (crWithinJunctions != null) {
                           
                            set.addAll(crWithinJunctions);
                        }
                    }
                }
            }
        }
        
        return set;
    }

    /**
     * search for the point (x, y) within the junction maps, and if found,
     * construct best defining CornerRegion for the point.  The best defining
     * region for the curvature neighbors appears to be the lowest intensity 
     * pixels when the edge surrounds a void or dark patch.
     * This method may need to be revised or used with segmentation or 
     * object recognition (in which case it probably belongs in a class
     * to specialize handling of the junction points for the EdgeExtractorWithJunctions).
     * 
     * Note that the CornerRegions are not necessarily in the same frame
     * as the original image (they do not have trimmedXOffset and trimmedYOffset
     * added to them).
     * 
     * NOTE: there is a flaw here that this instance does not retain the
     * scale space curves or the curvature of all edges currently,
     * so 2 dummy values for the neighboring curvatures are returned
     * (a flag is set in the CornerRegion to indicate that).
     * 
     * @param x
     * @param y
     * @param k not used, just stored in the returned CornerRegion
     * @return 
     */
    protected List<CornerRegion> searchForCornerRegionWithinJunctions(int x, 
        int y, float k) {
        
        Integer pixIndex = Integer.valueOf(img.getIndex(x, y));

        // if the corner is a junction, these are it's neighbors:
        Set<Integer> junctionAdjPixIndexes = junctionMap.get(pixIndex);

        PairInt edgeLocation = junctionLocationMap.get(pixIndex);
                
        if ((junctionAdjPixIndexes == null || junctionAdjPixIndexes.isEmpty())
            && (edgeLocation == null)) {
            return null;
        }
        
        if (junctionAdjPixIndexes == null || junctionAdjPixIndexes.isEmpty()) {
            // because edgeLocation is not null, this pixel is in a junction.
            // iterate over each of the 8 neighbors and if it's in the
            // junctionLocationMap, add it to junctionAdjPixIndexes for
            // making CornerRegions
            if (junctionAdjPixIndexes == null) {
                junctionAdjPixIndexes = new HashSet<Integer>();
            }
            int[] dx8 = new int[]{-1, -1,  0,  1, 1, 1, 0, -1};
            int[] dy8 = new int[]{ 0, -1, -1, -1, 0, 1, 1,  1};
            for (int i = 0; i < dx8.length; ++i) {
                int x1 = x + dx8[i];
                int y1 = y + dy8[i];
                if (x1 < 0 || y1 < 0 || (x1 > (img.getWidth() - 1)) || 
                (y1 > (img.getHeight() - 1))) {
                    continue;
                }
                Integer pixIndex1 = Integer.valueOf(img.getIndex(x1, y1));
                if (junctionLocationMap.containsKey(pixIndex1)) {
                    junctionAdjPixIndexes.add(pixIndex1);
                }
            }
        }
        
        List<CornerRegion> output = new ArrayList<CornerRegion>();
        
        List<Integer> list = new ArrayList<Integer>(junctionAdjPixIndexes);
        
        for (int i = 0; i < list.size(); ++i) {
            
            Integer pixIndex1 = list.get(i);
            
            PairInt edgeLocation1 = junctionLocationMap.get(pixIndex1);
            
            PairIntArray edge1 = edges.get(edgeLocation1.getX());
            
            //TODO: there's a bug in the junctionMap state!!!!  needs to be refactored
            int x1 = edge1.getX(edgeLocation1.getY());
            int y1 = edge1.getY(edgeLocation1.getY());
            
            for (int j = (i + 1); j < list.size(); ++j) {
                
                Integer pixIndex2 = list.get(j);
                                
                PairInt edgeLocation2 = junctionLocationMap.get(pixIndex2);
                
                PairIntArray edge2 = edges.get(edgeLocation2.getX());
                
                //TODO: there's a bug in the junctionMap state!!!!
                int x2 = edge2.getX(edgeLocation2.getY());
                int y2 = edge2.getY(edgeLocation2.getY());
            
                // making only 3 pixel regions instead of 5 for now
                
                CornerRegion cr = new CornerRegion(edgeLocation.getX(), 3, 1);
                cr.setFlagThatNeighborsHoldDummyValues();

                // add the image offsets
                cr.set(0, 0.1f*k, x1, y1);
                cr.set(1, k, x, y);
                cr.set(2, 0.1f*k, x2, y2);
                
                output.add(cr);
            }
        }
        
        return output;
    }
    
    /**
     * search for the point (x, y) within the junction maps and return true
     * if found.
     * 
     * @param x
     * @param y
     * @return 
     */
    protected boolean isInAJunction(int x, int y) {
        
        Integer pixIndex = Integer.valueOf(img.getIndex(x, y));

        Set<Integer> junctionAdjPixIndexes = junctionMap.get(pixIndex);

        if (junctionAdjPixIndexes != null && !junctionAdjPixIndexes.isEmpty()) {
            return true;
        }
        
        return junctionLocationMap.containsKey(pixIndex);
    }
}
