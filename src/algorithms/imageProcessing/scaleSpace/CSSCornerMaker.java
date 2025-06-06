package algorithms.imageProcessing.scaleSpace;

import algorithms.QuickSort;
import algorithms.imageProcessing.ImageStatisticsHelper;
import algorithms.compGeometry.MiscellaneousCurveHelper;
import algorithms.imageProcessing.SIGMA;
import algorithms.imageProcessing.features.CornerRegion;
import algorithms.misc.Histogram;
import algorithms.misc.HistogramHolder;
import algorithms.misc.MinMaxPeakFinder;
import algorithms.util.CornerArray;
import algorithms.util.Errors;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
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

    protected float factorIncreaseForCurvatureMinimum = 1.f;

    protected boolean performWholeImageCorners = true;

    protected boolean doStoreCornerRegions = true;

    protected final int width;

    protected final int height;

    protected Logger log = Logger.getLogger(this.getClass().getName());

    public CSSCornerMaker(int imageWidth, int imageHeight) {
        width = imageWidth;
        height = imageHeight;
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
    
    public List<CornerArray> findCornersInScaleSpaceMaps(
        final List<PairIntArray> theEdges) {
        
        /*
        for each edge, create scale space maps for each sigma in the range.
            and calculate corners for the edges for each sigma.
        
        after all scale space curve corners are calculated for all sigma,
            select the representative sigma as the one in which a corner
            is seen at its maximum value.
        
        create an output list of the corners for each edge at that selected
            sigma.
        */
        
        // for each edge, an item is a list of the corners for a sigm odered 
        // from max to min sigma
        List<List<CornerArray>> edgeCornerLists = new ArrayList<List<CornerArray>>();
       
        for (PairIntArray edge : theEdges) {
            List<CornerArray> scaleSpaceCorners = findCornersInScaleSpaceMaps(edge);
            edgeCornerLists.add(scaleSpaceCorners);
        }
        
        // for each edge, find the maximum curvature and its sigma
        List<Float> maxKs = new ArrayList<Float>();
        List<SIGMA> maxSigmas = new ArrayList<SIGMA>();
        
        for (List<CornerArray> scaleSpaceCorners : edgeCornerLists) {
                        
            // find the maximum k for each point and the sigma for it
            Map<Integer, Float> maxKMap = new HashMap<Integer, Float>();
            Map<Integer, SIGMA> maxKSigmasMap = new HashMap<Integer, SIGMA>();
            
            for (CornerArray corners : scaleSpaceCorners) {
                for (int i = 0; i < corners.getN(); ++i) {
                    
                    float k = Math.abs(corners.getCurvature(i));
                    Integer key = Integer.valueOf(i);
                    
                    Float kMax = maxKMap.get(key);
                    if ((kMax == null) || (kMax.floatValue() < k)) {
                        maxKMap.put(key, Float.valueOf(k));
                        maxKSigmasMap.put(key, corners.getSIGMA());
                    }
                }
            }
            
            float maxK = Float.NEGATIVE_INFINITY;
            SIGMA maxKSigma = null;
            for (Entry<Integer, Float> entry : maxKMap.entrySet()) {
                Float k = entry.getValue().floatValue();
                if (k > maxK) {
                    maxK = k;
                    maxKSigma = maxKSigmasMap.get(entry.getKey());
                }
            }
            maxSigmas.add(maxKSigma);
            maxKs.add(Float.valueOf(maxK));
        }
        
        QuickSort.<SIGMA>descendingSort(maxKs, maxSigmas);
        
        // choose the most frequent sigma from the top quarter of k's
        int n = maxKs.size()/4;
        if (n == 0) {
            n = maxKs.size();
        }
        Map<SIGMA, Integer> sigmaCountMap = new HashMap<SIGMA, Integer>();
        for (int i = 0; i < n; ++i) {
            SIGMA sigma = maxSigmas.get(i);
            Integer count = sigmaCountMap.get(sigma);
            if (count == null) {
                sigmaCountMap.put(sigma, Integer.valueOf(1));
            } else {
                sigmaCountMap.put(sigma, Integer.valueOf(count.intValue() + 1));
            }
        }
        int maxCount = Integer.MIN_VALUE;
        SIGMA selectedSigma = null;
        for (Entry<SIGMA, Integer> entry : sigmaCountMap.entrySet()) {
            int c = entry.getValue().intValue();
            if (c > maxCount) {
                selectedSigma = entry.getKey();
                maxCount = c;
            }
        }
        
        List<CornerArray> output = new ArrayList<CornerArray>();
        for (List<CornerArray> scaleSpaceCorners : edgeCornerLists) {
            boolean found = false;
            for (CornerArray corners : scaleSpaceCorners) {
                if (corners.getSIGMA() == selectedSigma) {
                    output.add(corners);
                    found = true;
                    break;
                }
            }
            if (!found) {
                output.add(new CornerArray(selectedSigma));
            }
        }
            
        output = filterCornersForMinResolvable(output);
        
        output = filterCornersUsing2ndDerivatives(output);
        
        return output;
    }
    
    /**
     * calculate corners using curvature scale space 
     * @param edge list of curves to create scale space curves of and calculate
     * corners from.
     * @return list of corners where indexes are the same as for the edges list.
     */
    public List<CornerArray> findCornersInScaleSpaceMaps(final PairIntArray 
        edge) {
        
        // ordered from min sigma to max sigma:
        final List<CornerArray> scaleSpaceCurvesList  =
            createLowUpperThresholdScaleSpaceMaps2(edge);
       
        if (scaleSpaceCurvesList.isEmpty()) {
            return scaleSpaceCurvesList;
        }
                
        List<CornerArray> output = new ArrayList<CornerArray>();
        
        /*
        find the candidate corners in the largest sigma scale space
           then filter the smaller sigma scale space curves to keep
           those same points.
        */
        
        CornerArray scaleSpaceCurve = scaleSpaceCurvesList.get(
            scaleSpaceCurvesList.size() - 1);
        CornerArray maxCorners = findCornersInScaleSpaceCurve(scaleSpaceCurve);
        //maxCorners = pruneCloseCorners(maxCorners, scaleSpaceCurve.getN());
        output.add(maxCorners);
        
        for (int i = 0; i < scaleSpaceCurvesList.size() - 1; ++i) {
            
            scaleSpaceCurve = scaleSpaceCurvesList.get(i);
            
            CornerArray corners = new CornerArray(scaleSpaceCurve.getSIGMA(), 
                 maxCorners.getN());
            if (scaleSpaceCurve.isFromAClosedCurve()) {
                corners.setIsClosedCurve();
            }
            
            for (int j = 0; j < maxCorners.getN(); ++j) {
                
                int idx = maxCorners.getInt(j);
                
                corners.add(scaleSpaceCurve.getX(idx), scaleSpaceCurve.getY(idx),
                    scaleSpaceCurve.getCurvature(idx), 
                    scaleSpaceCurve.getXFirstDeriv(idx),
                    scaleSpaceCurve.getXSecondDeriv(idx),
                    scaleSpaceCurve.getYFirstDeriv(idx),
                    scaleSpaceCurve.getYSecondDeriv(idx), idx);
            }
            
            output.add(corners);
        }
                
        return output;
    }
    
    /**
     * 
     * @param theEdges
     * @param outEdgeScaleSpaceMaps
     * @return 
     * @Deprecated 
     */
    public List<CornerArray> findCornersInScaleSpaceMaps(
        final List<PairIntArray> theEdges,
        List<Map<SIGMA, ScaleSpaceCurve>> outEdgeScaleSpaceMaps) {
        
        outEdgeScaleSpaceMaps.clear();
        
        //NOTE: this method is calculating the scale space maps twice
        // in order to avoid adding more code for a deprecated method
        // to bridge deprecated methods.
        //    if the method is retained and used, it should be refactored
        //    to preserve the CorneryArrays from the first creationg
        //    of scale maps and then create the Map<SIGMA, ScaleSpaceCurve>
        //    frm those.
        
        List<CornerArray> corners = findCornersInScaleSpaceMaps(theEdges);
        
        for (PairIntArray edge : theEdges) {
            
            Map<SIGMA, ScaleSpaceCurve> sMap = 
                createLowUpperThresholdScaleSpaceMaps(edge);
            
            outEdgeScaleSpaceMaps.add(sMap);
        }
        
        return corners;
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

            ScaleSpaceCurve curve = scaleSpaceHelper.computeCurvature(edge, 
                sigma, resultingSigma);

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
     * Construct scale space images (that is X(t, sigma), y(t, sigma), and
     * k(t, sigma) where t is the intervals of spacing along the curve
     * and is valued from 0 to 1 and k is the curvature).
     *
     * The range of the scale space maps is from sigma = 0.5 up to a maximum
     * value determined in getMaxSIGMAForECSS.
     *
     * The results are returned as a list ordered from small to large sigma.
     */
    private List<CornerArray> createLowUpperThresholdScaleSpaceMaps2(
        final PairIntArray edge) {

        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        boolean isClosedCurve = curveHelper.isAdjacent(edge, 0, 
            edge.getN() - 1) && (edge.getN() > 2);
        
        ScaleSpaceCurvature scaleSpaceHelper = new ScaleSpaceCurvature();

        List<CornerArray> output = new ArrayList<CornerArray>();
        
        SIGMA sigma = SIGMA.ZEROPOINTFIVE;

        SIGMA maxSIGMA = getMaxSIGMAForECSS(edge.getN());

        // this increases by a factor of sqrt(2)
        float resultingSigma = SIGMA.getValue(sigma);

        while (sigma.compareTo(maxSIGMA) < 1) {

            CornerArray curve = scaleSpaceHelper.computeCurvature2(edge, sigma);
            
            if (isClosedCurve) {
                curve.setIsClosedCurve();
            }
            
            output.add(curve);

            sigma = SIGMA.increaseToFactorBySQRT2(resultingSigma);

            resultingSigma *= Math.sqrt(2);
        }

        return output;
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
     * @param scaleSpaceCurve scale space curve for an edge
     * @return
     */
    protected CornerArray findCornersInScaleSpaceCurve(final CornerArray 
        scaleSpaceCurve) {
        
        if (scaleSpaceCurve.getN() == 0) {
            return new CornerArray(scaleSpaceCurve.getSIGMA(), 0);
        }

        float[] k = Arrays.copyOf(scaleSpaceCurve.getCurvature(), 
            scaleSpaceCurve.getN());

        float[] outputLowThreshold = new float[1];

        int[] minimaAndMaximaIdxs = findMinimaAndMaximaInCurvature(
            k, outputLowThreshold);

        List<Integer> maxCandidateCornerIndexes = findCandidateCornerIndexes(
            k, minimaAndMaximaIdxs, outputLowThreshold[0],
            this.factorIncreaseForCurvatureMinimum);

        CornerArray xy = new CornerArray(scaleSpaceCurve.getSIGMA(), 
            maxCandidateCornerIndexes.size());
        if (scaleSpaceCurve.isFromAClosedCurve()) {
            xy.setIsClosedCurve();
        }

        if (maxCandidateCornerIndexes.isEmpty()) {
            return xy;
        }

        int minDistFromEnds = 2;//5;
        int nPoints = scaleSpaceCurve.getN();
        for (int ii = 0; ii < maxCandidateCornerIndexes.size(); ii++) {

            int idx = maxCandidateCornerIndexes.get(ii);

            if (!scaleSpaceCurve.isFromAClosedCurve()) {
                if ((idx < minDistFromEnds)
                    || (idx > (nPoints - minDistFromEnds - 1))) {
                    continue;
                }
            }
            
            xy.add(scaleSpaceCurve.getX(idx), scaleSpaceCurve.getY(idx),
                scaleSpaceCurve.getCurvature(idx), 
                scaleSpaceCurve.getXFirstDeriv(idx),
                scaleSpaceCurve.getXSecondDeriv(idx),
                scaleSpaceCurve.getYFirstDeriv(idx),
                scaleSpaceCurve.getYSecondDeriv(idx), idx);
        }

        return xy;
    }

    /**
     * find the minima and maxima of the curvature k and return the low threshold
     * used in the in-out variable outputLowThreshold and a list of the
     * indexes for the minima as negative values and the maxima as positive
     * values.
     * @param k
     * @param outputLowThreshold array of size 1 to receive the low threshold used.
     * @return list of indexes of minima and maxima (minima have negative values).
     */
    public static int[] findMinimaAndMaximaInCurvature(float[] k,
        float[] outputLowThreshold) {

        if ((k == null) || (k.length < 5)) {
            return new int[0];
        }

        for (int ii = 1; ii < k.length; ii++) {
            if (k[ii] < 0.0f) {
                k[ii] *= -1.0f;
            }
        }

        if (k.length < 3) {
            return new int[0];
        }
 
        // get quartiles of non-zero values
        int nz = 0;
        for (float ki : k) {
            if (ki > 0) {
                nz++;
            }
        }
        if (nz < 4) {
            return new int[0];
        }
        float[] kNZ = new float[nz];
        nz = 0;
        for (float ki : k) {
            if (ki > 0) {
                kNZ[nz] = ki;
                nz++;
            }
        }
        float[] kQuartiles = ImageStatisticsHelper.getQuartiles(kNZ);

        // determine lowThresh
        HistogramHolder h = Histogram.calculateSturgesHistogram(
            0, 2 * kQuartiles[2], k, Errors.populateYErrorsBySqrt(k));

        if (h.getXHist().length < 2) {
            return new int[0];
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
            return new int[0];
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

        /*log.fine("lowThresh=" + outputLowThreshold[0]
            + " sum0=" + sum0 + " sum1=" + sum1 + " divsum=" + divSum
            + " firstPeakAndMinIdx[0]=" + firstPeakAndMinIdx[0]
            + " firstPeakAndMinIdx[1]=" + firstPeakAndMinIdx[1]);*/

        /*
        storing the minima and maxima in the same array list.
        the minima have -1*index within k
        and the maxima keep their positive values of the index within k.
        */
        MinMaxPeakFinder minMaxFinder = new MinMaxPeakFinder();
        
        int[] minMaxIdxs = minMaxFinder.findMinimaMaxima(k, outputLowThreshold[0]);

        return minMaxIdxs;
    }

    /**
     * given histogram h, find the first peak and the subsequent first minima
     * after it.
     *
     * @param h
     * @return
     */
    protected static int[] findFirstPeakAndMinimum(HistogramHolder h) {

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

    public static CornerRegion createCornerRegion(int edgeNumber, int cornerIdx,
        ScaleSpaceCurve scaleSpace, int imageWidth, int imageHeight) {

        final float[] k = scaleSpace.getK();

        //for 2 neighboring points on each side, min k is 0.2
        float kCenterAbs = Math.abs(k[cornerIdx]);
        if (kCenterAbs < 0.14f) {//0.18
            return null;
        }

        int n = scaleSpace.getSize();

        if (n < 3) {
            return null;
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
            if ((x < 0) || (y < 0) || (x > (imageWidth - 1)) 
                || (y > (imageHeight - 1))) {
                return null;
            }
            nCR++;
        }
        if (nCR < 3) {
            return null;
        }

        if (!isClosedCurve) {
            if (kMaxIdx == 0 || kMaxIdx == (nCR - 1)) {
                return null;
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
        
        return cr;
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
     * @param minMaxIdxs
     * @param lowThreshold
     * @param curvatureFactor2 factor which is multiplied by 2.0 to result
     * in the factor above minimum for a value to be significant.
     * @return
     */
    public static List<Integer> findCandidateCornerIndexes(float[] k,
        int[] minMaxIdxs, float lowThreshold, float curvatureFactor2) {

        // find peaks where k[ii] is > factorAboveMin* adjacent local minima

        float factorAboveMin = curvatureFactor2 * 2.0f;//3.5f;// 10 misses some corners

        //log.fine("using factorAboveMin=" + factorAboveMin);

        //to limit k to curvature that shows a rise in 1 pixel over a run of 3,
        // use 0.2 for a lower limit.
        // TODO: it's not clear that kLowerLimit is a good idea.  the relative change
        // filter alone is good for all size scale corners, and adding this
        // limit biases the results.  may want to only use this bias if
        // the some amount of curvature points are >= 0.2
        float kLowerLimit = 0.005f;

        List<Integer> cornerCandidates = new ArrayList<Integer>();

        // choose candidates from minMaxIndexes that are
        //     >= factorAboveMin one adjacent min
        for (int ii = 0; ii < minMaxIdxs.length; ii++) {

            int idx = minMaxIdxs[ii];

            if (idx > -1) {
                // this is maximum

                boolean found = false;
                
                // compare to preceding minimum
                for (int iii = (ii - 1); iii > -1; iii--) {
                    int idx2 = minMaxIdxs[iii];
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
                
                if (found && (ii > 0)) {
                    continue;
                }

                //compare to proceeding minimum
                for (int iii = (ii + 1); iii < minMaxIdxs.length; iii++) {
                    int idx2 = minMaxIdxs[iii];
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

    private int distSq(int x, int y, PairInt p) {
        
        int diffX = x - p.getX();
        int diffY = y - p.getY();
        
        return (diffX * diffX) + (diffY * diffY);
    }

    private CornerArray pruneCloseCorners(CornerArray corners, int curveLength) {
        
        int n = corners.getN();
        
        if (n < 2) {
            return corners;
        }
        
        boolean isClosedCurve = corners.isFromAClosedCurve();
        
        CornerArray output = new CornerArray(corners.getSIGMA());
        if (isClosedCurve) {
            output.setIsClosedCurve();
        }
        
        int startIdx = 0;
        int stopIdx = n;
        if (isClosedCurve & (corners.getInt(0) == 0)) {
            if (corners.getInt(n - 1) == (curveLength - 1)) {
                float k0 = corners.getCurvature(0);
                float kn1 = corners.getCurvature(n - 1);
                if (k0 >= kn1) {
                    stopIdx--;
                } else {
                    startIdx++;
                }
            }
        }
        
        double sqrt2 = Math.sqrt(2.);
        
        for (int i = startIdx; i < stopIdx; ++i) {
            if (i == (stopIdx - 1)) {
                output.add(corners.getX(i), corners.getY(i),
                    corners.getCurvature(i),
                    corners.getXFirstDeriv(i),
                    corners.getXSecondDeriv(i),
                    corners.getYFirstDeriv(i),
                    corners.getYSecondDeriv(i),
                    corners.getInt(i));
                continue;
            }
            int idx1 = corners.getInt(i);
            int idx2 = corners.getInt(i + 1);
            if (idx2 > (idx1 + 1)) {
                // checking for adjacent diagonal that need to be averaged,
                // else can add the isolated corner
                /*int x1 = corners.getX(i);
                int y1 = corners.getY(i);
                int x2 = corners.getX(i + 1);
                int y2 = corners.getY(i + 1);
                int diffx = x2 - x1;
                int diffy = y2 - y1;
                double dist = Math.sqrt(diffx*diffx + diffy*diffy);
                if (dist > sqrt2) { */               
                    output.add(corners.getX(i), corners.getY(i),
                        corners.getCurvature(i),
                        corners.getXFirstDeriv(i),
                        corners.getXSecondDeriv(i),
                        corners.getYFirstDeriv(i),
                        corners.getYSecondDeriv(i),
                        corners.getInt(i));
            } else {
                float k1 = corners.getCurvature(i);
                float k2 = corners.getCurvature(i + 1);
                if (k1 >= k2) {
                    // store this and skip next
                    output.add(corners.getX(i), corners.getY(i),
                        corners.getCurvature(i),
                        corners.getXFirstDeriv(i),
                        corners.getXSecondDeriv(i),
                        corners.getYFirstDeriv(i),
                        corners.getYSecondDeriv(i),
                        corners.getInt(i));
                    ++i;
                } else {
                    // store next and skip next
                    ++i;
                    output.add(corners.getX(i), corners.getY(i),
                        corners.getCurvature(i),
                        corners.getXFirstDeriv(i),
                        corners.getXSecondDeriv(i),
                        corners.getYFirstDeriv(i),
                        corners.getYSecondDeriv(i),
                        corners.getInt(i));
                }
            }
        }
        return output;
    }

    private List<CornerArray> filterCornersForMinResolvable(List<CornerArray> cList) {
        
        if (cList.isEmpty()) {
            return cList;
        }
        
        /*
        filtering by the minimum detectable angle and hence curvature for it.
        
         solid angle where r = radius of curvature.  k=1/r.
                   .
                  /|\
                 / | \
                / r-h \r
               /   |   \
              .----|----.
                   -     bottom portion is a triangle           w
                                                        .----.-----.
                                                          .  |h .
                                                             .
         the curvature is too small to determine slopes from neighboring
         points when h is less than 1 pixel and w is 3 or more.
       
         limit to k for h=1.0 and w=3:
       
              r^2 = (r-h)^2 + w^2
       
              r^2 = r^2 - 2*r*h + h^2 + w^2
              2*r*h = h^2 + w^2
                r = (h^2 + w^2)/(2*h)
                r = 5  which is k = 0.2 for no bluring
        
        gaussian bluring results in FWHM of size 2.35 * sigma.
        
        so with a sigma=1
        
        h = 1.0 * 2.35 * 1 = 2.35
        w = 3.0 * 2.35 * 1 = 7.05
        r = 11.75, so min absolute curvature = 1/r = 0.085
        */
        
        List<CornerArray> output = new ArrayList<CornerArray>();
        for (CornerArray corners : cList) {
         
            if (corners.getN() == 0) {
                output.add(corners);
                continue;
            }
            
            SIGMA sigma = corners.getSIGMA();
            CornerArray corners2 = new CornerArray(sigma);
            if (corners.isFromAClosedCurve()) {
                corners2.setIsClosedCurve();
            }
            
            double h = 2.35 * SIGMA.getValue(sigma);
            double w = 3.*h;
            
            double minCurvature = (2. * h)/(h * h + w * w);
            
            for (int i = 0; i < corners.getN(); ++i) {
                double k = Math.abs(corners.getCurvature(i));
                if (k >= minCurvature) {
                    corners2.add(
                        corners.getX(i), corners.getY(i),
                        corners.getCurvature(i),
                        corners.getXFirstDeriv(i),
                        corners.getXSecondDeriv(i),
                        corners.getYFirstDeriv(i),
                        corners.getYSecondDeriv(i), corners.getInt(i));
                }
            }
            
            output.add(corners2);
        }
        
        return output;
    }

    private List<CornerArray> filterCornersUsing2ndDerivatives(List<CornerArray> cList) {
        
        if (cList.isEmpty()) {
            return cList;
        }
        
        /*
        filtering using the property derived by Harris and Stephens 
        R = det(G) − k(tr(G))^2 here k = 0.04
        
        and G is the matrix | I_x*I_x  I_x*I_y |
                            | I_x*I_y  I_y*I_y |
        where I_x and I_y are the 2nd derivtives of the intensity.
        see Eqn (11) in "Phase Congruency Detects Corners and Edges"
        by Kovesi.
        
        Note that the limit to use with this value is sensitive to image
        contrast, so one may prefer to use a related, but different
        measure of localizability that takes the measurements along the
        major and minor axes of the feature orienation.
        See FeatureHelper.filterByLocalizability(...)
        
        */
        
        List<CornerArray> output = new ArrayList<CornerArray>();
        for (CornerArray corners : cList) {
         
            if (corners.getN() == 0) {
                output.add(corners);
                continue;
            }
            
            SIGMA sigma = corners.getSIGMA();
            CornerArray corners2 = new CornerArray(sigma);
            if (corners.isFromAClosedCurve()) {
                corners2.setIsClosedCurve();
            }
            
            for (int i = 0; i < corners.getN(); ++i) {
                
                float x2d = corners.getXSecondDeriv(i);
                float y2d = corners.getYSecondDeriv(i);
                
                float dxdx = x2d*x2d;
                float dydy = y2d*y2d;
                float dxdy = x2d*x2d;
                
                float det = dxdx * dydy - dxdy * dxdy;
                
                float trace = dxdx + dydy;
                
                float r = Math.abs(det - 0.04f * (trace * trace));
                
                if (r > 2) {
                    corners2.add(
                        corners.getX(i), corners.getY(i),
                        corners.getCurvature(i),
                        corners.getXFirstDeriv(i),
                        corners.getXSecondDeriv(i),
                        corners.getYFirstDeriv(i),
                        corners.getYSecondDeriv(i), corners.getInt(i));
                }
            }
            
            output.add(corners2);
        }
        
        return output;
    }
}
