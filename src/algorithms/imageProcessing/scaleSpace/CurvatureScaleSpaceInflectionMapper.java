package algorithms.imageProcessing.scaleSpace;

import algorithms.MultiArrayMergeSort;
import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.transform.MatchedPointsTransformationCalculator;
import algorithms.imageProcessing.transform.TransformationParameters;
import algorithms.imageProcessing.util.AngleUtil;
import algorithms.misc.MiscDebug;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.logging.Logger;

/**
 * class to map contours from one image to another and return the
 * matched inflection points and a transformation matrix that can
 * be applied to the first to put it into the frame of the second.
 *
 * The algorithm used for matching scale space image contours is documented in
 * CSSContourMatcherWrapper
 * @see algorithms.imageProcessing.CSSContourMatcherWrapper
 *
 * @author nichole
 */
public final class CurvatureScaleSpaceInflectionMapper {

    protected final Logger log = Logger.getLogger(this.getClass().getName());
    protected boolean debug = false;
    protected boolean useLineDrawingMode = false;
    protected boolean doRefineTransformations = false;
    protected boolean initialized = false;

    protected CSSContourInflectionMaker cssIM1 = null;
    protected CSSContourInflectionMaker cssIM2 = null;
    
    protected Map<Integer, PairIntArray> matchedXY1ByEdgeInOrigRefFrame = null;
    protected Map<Integer, PairIntArray> matchedXY2ByEdgeInOrigRefFrame = null;
    protected Map<Integer, List<Float>> matchedXY1ByEdgeWeights = null;
    protected Map<Integer, List<Float>> matchedXY2ByEdgeWeights = null;
   
    /**
     * matched points from the contour lists of image 1 (matched to the same
     * in image 2) with coordinates being in the reference frames of the
     * original image 1 before any trimming.
     */
    protected PairIntArray matchedXY1 = null;
    
    /**
     * matched points from the contour lists of image 2 (matched to the same
     * in image 1) with coordinates being in the reference frames of the
     * original image 2 before any trimming.
     */
    protected PairIntArray matchedXY2 = null;
    
    /**
     * weights for points in matchedXY1 created from the peak strengths.
     */
    protected float[] matchedXY1Weights = null;
    
    /**
     * weights for points in matchedXY2 created from the peak strengths.
     */
    protected float[] matchedXY2Weights = null;
    
    /**
     * indexes for edges from edges1 which produced matching contours
     */
    protected int[] matchedEdge1Indexes = null;
    
    /**
     * indexes for edges from edges2 which produced matching contours
     */
    protected int[] matchedEdge2Indexes = null;

    protected TransformationParameters bestFittingParameters = null;

    /**
     * scale derived from matching contours.  it's not necessarily the same
     * as the final scale returned in transformation solutions, but it should
     * be close;
     */
    private double matchedScale = 1;

    private List<CurvatureScaleSpaceContour> matchedContours1 = new
        ArrayList<CurvatureScaleSpaceContour>();

    private List<CurvatureScaleSpaceContour> matchedContours2 = new
        ArrayList<CurvatureScaleSpaceContour>();

    private final ImageExt image1;
    private final ImageExt image2;
    
    public CurvatureScaleSpaceInflectionMapper(ImageExt image1, ImageExt image2) {

        this.image1 = image1;
        this.image2 = image2;
    }

    public void useLineDrawingLineMode() {
        this.useLineDrawingMode = true;
    }

    public void setToRefineTransformations() {
        doRefineTransformations = true;
    }

    public void useDebugMode() {
        debug = true;
    }

    protected void initialize() {

        if (initialized) {
            return;
        }

        initialized = true;

        this.cssIM1 = new CSSContourInflectionMaker(image1);
        
        this.cssIM2 = new CSSContourInflectionMaker(image2);
        
    }

    
    protected TransformationParameters createEuclideanTransformationImpl() {

        if (bestFittingParameters == null) {
            return null;
        }

        MatchedPointsTransformationCalculator tc = new
            MatchedPointsTransformationCalculator();

        if (debug) {
            tc.useDebugMode();
        }

        //TODO: temporarily disabling the refinement while fixing PointMatcher
        if (false && doRefineTransformations) {

            /*
            boolean reverseDatasetOrder = bestFittingParameters.getScale() < 1.0;

            log.info("BEFORE REFINEMENT:\n" + bestFittingParameters.toString());

            PairIntArray[] set1 = getMatchedEdges1InOriginalReferenceFrameArray();
            PairIntArray[] set2 = getMatchedEdges2InOriginalReferenceFrameArray();
            EdgeMatcher matcher = new EdgeMatcher();
            TransformationPointFit fit2 = null;
            if (reverseDatasetOrder) {
                fit2 = matcher.refineTransformation(set2, set1, bestFittingParameters);
            } else {
                fit2 = matcher.refineTransformation(set1, set2, bestFittingParameters);
            }

            if (reverseDatasetOrder) {
                bestFittingParameters = tc.swapReferenceFrames(bestFittingParameters);
            }

            if (fit2 != null) {
                log.info("FINAL:\n" + fit2.toString());
                bestFittingParameters = fit2.getParameters();
            }*/
        }

        return bestFittingParameters;
    }

    protected List<PairIntArray> getEdges(CurvatureScaleSpaceImageMaker imgMaker) {

        return imgMaker.getClosedCurves();
    }

    protected Image getImage1() {
        return image1;
    }

    protected Image getImage2() {
        return image2;
    }

    /**
     * NOT READY FOR USE
     */
    protected void createMatchedPointArraysFromContourPeaks() {

        if (matchedXY1 != null) {
            return;
        }

        final int centroidX1 = 0;
        final int centroidY1 = 0;
        final int centroidX2 = 0;
        final int centroidY2 = 0;

        Map<Integer, List<CurvatureScaleSpaceContour>> bestMatches1 = new
            HashMap<Integer, List<CurvatureScaleSpaceContour>>();

        Map<Integer, List<CurvatureScaleSpaceContour>> bestMatchesTo1 = new
            HashMap<Integer, List<CurvatureScaleSpaceContour>>();

        Map<Integer, Float> bestScales = new HashMap<Integer, Float>();

        Map<Integer, Integer> bestI2I1 = new HashMap<Integer, Integer>();

        TreeMap<Double, Set<Integer>> bestCosts = new TreeMap<Double, Set<Integer>>();

        Map<Integer, TransformationParameters> bestParams = new
            HashMap<Integer, TransformationParameters>();

        Map<Integer, PairIntArray> bestMatchesXY1 = new
            HashMap<Integer, PairIntArray>();

        Map<Integer, List<Float>> bestMatchesXYWeights1 = new
            HashMap<Integer, List<Float>>();

        Map<Integer, PairIntArray> bestMatchesXY2 = new
            HashMap<Integer, PairIntArray>();

        Map<Integer, List<Float>> bestMatchesXYWeights2 = new
            HashMap<Integer, List<Float>>();

        Map<Integer, Double> i2CostMap = new HashMap<Integer, Double>();

        boolean alreadySorted = true;
        
        List<List<CurvatureScaleSpaceContour>> contourLists1 = cssIM1.getContours();
        
        List<List<CurvatureScaleSpaceContour>> contourLists2 = cssIM2.getContours();

        List<PairIntArray> edges1 = cssIM1.getEdges();
        
        List<PairIntArray> edges2 = cssIM2.getEdges();
        
        for (int i1 = 0; i1 < contourLists1.size(); ++i1) {

            List<CurvatureScaleSpaceContour> contours1 = contourLists1.get(i1);

            double minCost = Double.MAX_VALUE;
            List<CurvatureScaleSpaceContour> bestM1 = null;
            List<CurvatureScaleSpaceContour> bestM2 = null;
            int bestI2ForThisI1 = -1;
            double bestScale = 1;
            double bestCost = Double.MAX_VALUE;

if (debug){
Image img3 = cssIM1.getImage().createWithDimensions();
for (int j = 0; j < edges1.get(i1).getN(); ++j) {
    int x = edges1.get(i1).getX(j);
    int y = edges1.get(i1).getY(j);
    if (j == 0 || (j == (edges1.get(i1).getN() - 1))) {
        ImageIOHelper.addPointToImage(x, y, img3, 0, 200, 150, 0);
    } else {
        ImageIOHelper.addPointToImage(x, y, img3, 0, 255, 0, 0);
    }
}
MiscDebug.writeImage(img3, "edge1_" + i1);
        }
            for (int i2 = 0; i2 < contourLists2.size(); ++i2) {

                List<CurvatureScaleSpaceContour> contours2 = contourLists2.get(i2);

                CSSContourMatcherWrapper matcher =
                    new CSSContourMatcherWrapper(contours1, contours2,
                    alreadySorted);
if (debug){
Image img3 = cssIM2.getImage().createWithDimensions();
for (int j = 0; j < edges2.get(i2).getN(); ++j) {
    int x = edges2.get(i2).getX(j);
    int y = edges2.get(i2).getY(j);
    
    if (j == 0 || (j == (edges2.get(i2).getN() - 1))) {
        ImageIOHelper.addPointToImage(x, y, img3, 0, 200, 150, 0);
    } else {
        ImageIOHelper.addPointToImage(x, y, img3, 0, 255, 0, 0);
    }
}
MiscDebug.writeImage(img3, "edge2_" + i2);
}

                boolean didMatch = matcher.matchContours();

                if (!didMatch) {
                    continue;
                }

                List<CurvatureScaleSpaceContour> m1 = matcher.getSolutionMatchedContours1();
                List<CurvatureScaleSpaceContour> m2 = matcher.getSolutionMatchedContours2();
                if (m1 == null || m2 == null || m1.isEmpty() || m2.isEmpty()) {
                    continue;
                }
                assert(m1.size() == m2.size());

                /*
                There may be insignificant low cost matches for very small
                curves, so will only keep a solution when there are as few
                as 2 contours in the match if there are no other matches.
                */
                if ((m1.size() == 2) && (bestM1 != null) && (bestM1.size() > 2)) {
                    continue;
                }

                double cost = matcher.getSolvedCost();
/*
try {
// plot xy of edge
// plot contour points
// plot space image
int flNumber = MiscDebug.getCurrentTimeFormatted();
int edgeIdx1 = m1.get(0).getEdgeNumber();
int edgeIdx2 = m2.get(0).getEdgeNumber();
PairIntArray txy1 = new PairIntArray(m1.size());
PairIntArray txy2 = new PairIntArray(m2.size());
List<Float> tweights1 = new ArrayList<Float>();
List<Float> tweights2 = new ArrayList<Float>();
extract(m1, txy1, tweights1);
extract(m2, txy2, tweights2);
MiscDebug.writeImage(txy1, ImageIOHelper.convertImage(debugImg1), "check_1_xy_edge_" + edgeIdx1 + "_" + flNumber);
MiscDebug.writeImage(txy2, ImageIOHelper.convertImage(debugImg2), "check_2_xy_edge_" + edgeIdx2 + "_" + flNumber);
MiscDebug.debugPlot(contours1, ImageIOHelper.convertImage(debugImg1), "1_edge_" + edgeIdx1 + "_" + String.valueOf(flNumber));
MiscDebug.debugPlot(contours2, ImageIOHelper.convertImage(debugImg2), "2_edge_" + edgeIdx1 + "_" + String.valueOf(flNumber));
int z = 1;
} catch (IOException ex) {

}
*/

                if ((cost < minCost) || ((bestM1 != null) && (m1.size() > 2) && (bestM1.size() < 3))) {

                    // if i2 is already matched to the best of i1 and the
                    // cost there is smaller, cannot set to best here

                    Double prevI2Cost = i2CostMap.get(Integer.valueOf(i2));

                    boolean assign = (prevI2Cost == null);

                    if (!assign) {
                        assign = (cost < prevI2Cost.doubleValue());
                    }

                    if (assign) {

                        bestI2ForThisI1 = i2;
                        minCost = cost;
                        bestM1 = m1;
                        bestM2 = m2;
                        bestScale = matcher.getSolvedScale();
                        bestCost = cost;

                        log.info(" best so far has cost=" + bestCost + " i1=" + i1 + " i2=" + i2);
                    }
                }
            }

            if (bestM1 != null) {

                // calculate the implied transformation from these matched points

                PairIntArray xy1 = new PairIntArray(bestM1.size());
                PairIntArray xy2 = new PairIntArray(bestM2.size());

                List<Float> weights1 = new ArrayList<Float>();
                List<Float> weights2 = new ArrayList<Float>();

                extract(bestM1, xy1, weights1);
                extract(bestM2, xy2, weights2);

                if (xy1.getN() < 3) {
                    continue;
                }
/*
try {
    Image tmp = cssIM1.getImage().copyImage();
    MiscDebug.writeImage(xy1, tmp,
        "check_1_xy_" + MiscDebug.getCurrentTimeFormatted());
} catch (IOException ex) {
    Logger.getLogger(CurvatureScaleSpaceInflectionMapper.class.getName()).log(Level.SEVERE, null, ex);
}
*/
                MatchedPointsTransformationCalculator tc = new
                    MatchedPointsTransformationCalculator();

                TransformationParameters params = null;

                // if scale < 1, we have to swap the order of datasets to avoid
                // numerical errors in some of the methods that are the result of
                // dividing by a small number
                boolean reverseDatasetOrder = bestScale < 1.0;
                if (reverseDatasetOrder) {
                    params = tc.calulateEuclideanGivenScale(1. / bestScale,
                        xy2, xy1, centroidX2, centroidY2);
                } else {
                    params = tc.calulateEuclideanGivenScale(bestScale,
                        xy1, xy2, centroidX1, centroidY1);
                }
                if (params == null) {
                    continue;
                }

                if (reverseDatasetOrder && (params != null)) {
                    params = tc.swapReferenceFrames(params);
                }
/*
try {
    int flNumber = MiscDebug.getCurrentTimeFormatted();
    Image tmp1 = cssIM1.getImage().copyImage();
    Image tmp2 = cssIM2.getImage().copyImage();
    MiscDebug.writeImage(xy1, tmp1, "check_1_xy_" + flNumber);
    MiscDebug.writeImage(xy2, tmp2, "check_2_xy_" + flNumber);
    Transformer transformer = new Transformer();
    PairIntArray xy1Tr = transformer.applyTransformation(params, xy1);
    MiscDebug.writeImage(xy1Tr, (ImageExt)image2.copyImage(),
        "check_1_xy_tr_" + flNumber);
} catch (IOException ex) {
    Logger.getLogger(CurvatureScaleSpaceInflectionMapper.class.getName()).log(Level.SEVERE, null, ex);
}
*/

                Integer index1 = Integer.valueOf(i1);

                bestMatches1.put(index1, bestM1);
                bestMatchesTo1.put(index1, bestM2);
                bestScales.put(index1, Double.valueOf(bestScale).floatValue());
                bestParams.put(index1, params);
                bestMatchesXY1.put(index1, xy1);
                bestMatchesXY2.put(index1, xy2);
                bestMatchesXYWeights1.put(index1, weights1);
                bestMatchesXYWeights2.put(index1, weights2);

                Double key2 = Double.valueOf(bestCost);
                if (!bestCosts.containsKey(key2)) {
                    bestCosts.put(key2, new HashSet<Integer>());
                }
                bestCosts.get(key2).add(Integer.valueOf(bestI2ForThisI1));

                i2CostMap.put(Integer.valueOf(bestI2ForThisI1), key2);

                bestI2I1.put(Integer.valueOf(bestI2ForThisI1), index1);
            }
        }

        /*
        the sigmas of the peaks of the contours need to be used here when
        combining or prefering solutions between edges having no common edge.
        Will use the "penalty" formula from the paper which adds the difference
        from the tallest first peak to all other tallest first peaks to the costs.
        */
        if (bestI2I1.size() > 1) {
            adjustCostToTallestContourPeak1(bestMatches1, bestI2I1, bestCosts, i2CostMap);
        }

        /*
        compare the solutions, starting with the smallest cost solution.
        */
        int nTransformations = bestParams.size();

        /* calculate the highest number of similar transformations and the
        lowest cost from those.  store nSimilar, indexes, cost for each iteration*/
        int[] nSimilarSummary = new int[nTransformations];
        Integer[][] indexesSummary = new Integer[nTransformations][];
        double[] costsSummary = new double[nTransformations];
        int[] mainIndexSummary = new int[nTransformations];

        int count = 0;

        for (Map.Entry<Double, Set<Integer>> entry : bestCosts.entrySet()) {

            Set<Integer> indexes2 = entry.getValue();

            for (Integer index2 : indexes2) {

                Integer index1 = bestI2I1.get(index2);

                Set<Integer> similarParamsIndexes1 = new HashSet<Integer>();

                TransformationParameters params = bestParams.get(index1);

                if (params == null) {
                    continue;
                }

                similarParamsIndexes1.add(index1);

                for (Map.Entry<Integer, TransformationParameters> entryP : bestParams.entrySet()) {
                    Integer index1P = entryP.getKey();
                    if (index1P.equals(index1)) {
                        continue;
                    }
                    TransformationParameters paramsP = entryP.getValue();
                    if (paramsP == null) {
                        continue;
                    }
                    if (Math.abs(params.getScale() - paramsP.getScale()) < 0.05) {
                        float angleDiff = AngleUtil.getAngleAverageInDegrees(
                            params.getRotationInDegrees(), paramsP.getRotationInDegrees());
                        if (Math.abs(angleDiff) < 10) {
                            if (Math.abs(params.getTranslationX() - paramsP.getTranslationX()) < 10) {
                                if (Math.abs(params.getTranslationY() - paramsP.getTranslationY()) < 10) {
                                    similarParamsIndexes1.add(index1P);
                                }
                            }
                        }
                    }
                }
                nSimilarSummary[count] = similarParamsIndexes1.size();
                indexesSummary[count] = similarParamsIndexes1.toArray(new Integer[similarParamsIndexes1.size()]);
                costsSummary[count] = entry.getKey();
                mainIndexSummary[count] = index1.intValue();
                count++;
            }
        }

        if (count == 0) {
            return;
        }
        
        nSimilarSummary = Arrays.copyOf(nSimilarSummary, count);
        indexesSummary = Arrays.copyOf(indexesSummary, count);
        costsSummary = Arrays.copyOf(costsSummary, count);
        mainIndexSummary = Arrays.copyOf(mainIndexSummary, count);

        //MultiArrayMergeSort.sortBy1stDescThen2ndAsc(nSimilarSummary, costsSummary, indexesSummary, mainIndexSummary);

        MultiArrayMergeSort.sortBy1stAscThen2ndDesc(costsSummary, nSimilarSummary, indexesSummary, mainIndexSummary, 0, costsSummary.length - 1);

        Integer[] indexes = indexesSummary[0];
        int mainIndex = mainIndexSummary[0];

        bestFittingParameters = bestParams.get(Integer.valueOf(mainIndex));
        matchedScale = bestFittingParameters.getScale();

        matchedXY1ByEdgeInOrigRefFrame = new HashMap<Integer, PairIntArray>();
        matchedXY2ByEdgeInOrigRefFrame = new HashMap<Integer, PairIntArray>();
        matchedXY1ByEdgeWeights = new HashMap<Integer, List<Float>>();
        matchedXY2ByEdgeWeights = new HashMap<Integer, List<Float>>();

        matchedEdge1Indexes = new int[indexes.length];
        matchedEdge2Indexes = new int[indexes.length];

        matchedXY1 = new PairIntArray();
        matchedXY2 = new PairIntArray();

        for (int i = 0; i < indexes.length; ++i) {

            Integer index1 = indexes[i];

            List<CurvatureScaleSpaceContour> m1 = bestMatches1.get(index1);
            List<CurvatureScaleSpaceContour> m2 = bestMatchesTo1.get(index1);
            matchedContours1.addAll(m1);
            matchedContours2.addAll(m2);

            Integer e1Index = null;
            Integer e2Index = null;

            for (int mIdx1 = 0; mIdx1 < 1; ++mIdx1) {
                CurvatureScaleSpaceContour c1 = m1.get(mIdx1);
                CurvatureScaleSpaceContour c2 = m2.get(mIdx1);
                e1Index = Integer.valueOf(c1.getEdgeNumber());
                e2Index = Integer.valueOf(c2.getEdgeNumber());
            }

            matchedXY1ByEdgeInOrigRefFrame.put(e1Index, bestMatchesXY1.get(index1));
            matchedXY2ByEdgeInOrigRefFrame.put(e2Index, bestMatchesXY2.get(index1));

            //These do not have the adjusted costs included:
            matchedXY1ByEdgeWeights.put(e1Index, bestMatchesXYWeights1.get(index1));
            matchedXY2ByEdgeWeights.put(e2Index, bestMatchesXYWeights2.get(index1));

            matchedXY1.addAll(bestMatchesXY1.get(index1));
            matchedXY2.addAll(bestMatchesXY2.get(index1));

            matchedEdge1Indexes[i] = e1Index;
            matchedEdge2Indexes[i] = e2Index;
        }
    }

    public TransformationParameters createEuclideanTransformation() {

        initialize();
        
        List<List<CurvatureScaleSpaceContour>> contourLists1 =
            cssIM1.getContours();
        
        List<List<CurvatureScaleSpaceContour>> contourLists2 =
            cssIM2.getContours();

        if (contourLists2.isEmpty() || contourLists1.isEmpty()) {
            return null;
        }

        createMatchedPointArraysFromContourPeaks();

        return createEuclideanTransformationImpl();
    }

    private void extract(List<CurvatureScaleSpaceContour> contours,
        PairIntArray outputXY, List<Float> outputSigmaWeights) {

        float sumSigma = 0;

        for (int i = 0; i < contours.size(); i++) {

            CurvatureScaleSpaceContour c = contours.get(i);

            for (int j = 0; j < c.getPeakDetails().length; j++) {

                CurvatureScaleSpaceImagePoint spaceImagePoint =
                    c.getPeakDetails()[j];

                int x = spaceImagePoint.getXCoord();
                int y = spaceImagePoint.getYCoord();

                outputXY.add(x, y);
                outputSigmaWeights.add(Float.valueOf(c.getPeakSigma()));

                sumSigma += c.getPeakSigma();
            }
        }

        for (int i = 0; i < outputSigmaWeights.size(); ++i) {
            float w = outputSigmaWeights.get(i)/sumSigma;
            outputSigmaWeights.set(i, w);
        }
    }

    protected void correctPeaks0(List<CurvatureScaleSpaceContour> matched1,
        List<CurvatureScaleSpaceContour> matched2) {

        if (matched1.size() != matched2.size()) {
            throw new IllegalArgumentException("lengths of matched1"
            + " and matchedContours2 must be the same");
        }

        // the contours extracted from scale space images using a factor of
        // 2^(1/8) for recursive convolution tend to not have a single
        // peak, so the correction here for the single peak case is not
        // usually needed.  for that rare case, the avg of the other peak
        // is stored instead of both points

        for (int i = 0; i < matched1.size(); i++) {

            CurvatureScaleSpaceContour c1 = matched1.get(i);
            CurvatureScaleSpaceContour c2 = matched2.get(i);

            if (c1.getPeakDetails().length != c2.getPeakDetails().length) {
                if (c1.getPeakDetails().length == 1) {
                    CurvatureScaleSpaceImagePoint p0 = c2.getPeakDetails()[0];
                    CurvatureScaleSpaceImagePoint p1 = c2.getPeakDetails()[1];
                    float t = p0.getScaleFreeLength();
                    float s = p0.getSigma();
                    int xAvg = Math.round((p0.getXCoord() + p1.getXCoord()) / 2.f);
                    int yAvg = Math.round((p0.getYCoord() + p1.getYCoord()) / 2.f);
                    CurvatureScaleSpaceImagePoint pAvg =
                        new CurvatureScaleSpaceImagePoint(s, t, xAvg, yAvg,
                        p0.getCoordIdx());
                    CurvatureScaleSpaceImagePoint[] p =
                        new CurvatureScaleSpaceImagePoint[]{pAvg};
                    c2.setPeakDetails(p);
                    matched2.set(i, c2);
                }  else if (c2.getPeakDetails().length == 1) {
                    CurvatureScaleSpaceImagePoint p0 = c1.getPeakDetails()[0];
                    CurvatureScaleSpaceImagePoint p1 = c1.getPeakDetails()[1];
                    float t = p0.getScaleFreeLength();
                    float s = p0.getSigma();
                    int xAvg = Math.round((p0.getXCoord() + p1.getXCoord()) / 2.f);
                    int yAvg = Math.round((p0.getYCoord() + p1.getYCoord()) / 2.f);
                    CurvatureScaleSpaceImagePoint pAvg =
                        new CurvatureScaleSpaceImagePoint(s, t, xAvg, yAvg,
                        p0.getCoordIdx());
                    CurvatureScaleSpaceImagePoint[] p =
                        new CurvatureScaleSpaceImagePoint[]{pAvg};
                    c1.setPeakDetails(p);
                    matched1.set(i, c1);
                }
            }
        }
    }

    public PairIntArray getMatchedXY1() {
        return matchedXY1;
    }

    public PairIntArray getMatchedXY2() {
        return matchedXY2;
    }

    public float[] getMatchedXY1Weights() {
        return matchedXY1Weights;
    }

    public float[] getMatchedXY2Weights() {
        return matchedXY2Weights;
    }

    public List<List<CurvatureScaleSpaceContour>> getContours1() {
        return cssIM1.getContours();
    }

    public List<List<CurvatureScaleSpaceContour>> getContours2() {
        return cssIM2.getContours();
    }

    protected List<PairIntArray> getEdges1() {
        return cssIM1.getEdges();
    }

    protected List<PairIntArray> getEdges2() {
        return cssIM2.getEdges();
    }

    protected PairIntArray[] getEdges1Array() {
        List<PairIntArray> edges1 = cssIM1.getEdges();
        PairIntArray[] oe = new PairIntArray[edges1.size()];
        for (int i = 0; i < edges1.size(); i++) {
            PairIntArray edge = edges1.get(i).copy();
            
            oe[i] = edge;
        }
        return oe;
    }

    protected PairIntArray[] getEdges2InArray() {
        List<PairIntArray> edges2 = cssIM2.getEdges();
        PairIntArray[] oe = new PairIntArray[edges2.size()];
        for (int i = 0; i < edges2.size(); i++) {
            PairIntArray edge = edges2.get(i).copy();
            
            oe[i] = edge;
        }
        return oe;
    }

    protected PairIntArray[] getMatchedEdges1() {
        if (matchedEdge1Indexes == null) {
            return new PairIntArray[0];
        }
        PairIntArray[] oe = new PairIntArray[matchedEdge1Indexes.length];
        List<PairIntArray> edges1 = cssIM1.getEdges();
        for (int i = 0; i < matchedEdge1Indexes.length; i++) {
            int eIdx = matchedEdge1Indexes[i];
            PairIntArray edge = edges1.get(eIdx).copy();
            
            oe[i] = edge;
        }
        return oe;
    }

    protected PairIntArray[] getMatchedEdges2() {
        if (matchedEdge2Indexes == null) {
            return new PairIntArray[0];
        }
        PairIntArray[] oe = new PairIntArray[matchedEdge2Indexes.length];
        List<PairIntArray> edges2 = cssIM2.getEdges();
        for (int i = 0; i < matchedEdge2Indexes.length; i++) {
            int eIdx = matchedEdge2Indexes[i];
            
            PairIntArray edge = edges2.get(eIdx).copy();
            
            oe[i] = edge;
        }
        return oe;
    }

    public List<CurvatureScaleSpaceContour> getMatchedContours1() {
        return matchedContours1;
    }

    public List<CurvatureScaleSpaceContour> getMatchedContours2() {
        return matchedContours2;
    }

    public double getMatchedScale() {
        return matchedScale;
    }

    public PairInt[] getMatchedEdgesIndexes() {
        List<Integer> idx1 = new ArrayList<Integer>();
        List<Integer> idx2 = new ArrayList<Integer>();
        for (int i = 0; i < this.matchedContours1.size(); i++) {
            Integer edge1Idx = Integer.valueOf(matchedContours1.get(i).getEdgeNumber());
            Integer edge2Idx = Integer.valueOf(matchedContours2.get(i).getEdgeNumber());
            if (!idx1.contains(edge1Idx)) {
                idx1.add(edge1Idx);
                idx2.add(edge2Idx);
            }
        }
        PairInt[] indexes = new PairInt[idx1.size()];
        for (int i = 0; i < idx1.size(); i++) {
            indexes[i] = new PairInt(idx1.get(i), idx2.get(i));
        }
        return indexes;
    }

    private float findMaxSigmaOfFirstPeaks(
        Map<Integer, List<CurvatureScaleSpaceContour>> bestMatches1) {

        float maxPeakSigma = Float.MIN_VALUE;

        for (Map.Entry<Integer, List<CurvatureScaleSpaceContour>> entry : bestMatches1.entrySet()) {

            List<CurvatureScaleSpaceContour> list = entry.getValue();
            if (list.isEmpty()) {
                continue;
            }
            CurvatureScaleSpaceContour contour = list.get(0);
            float peakSigma = contour.getPeakSigma();
            if (peakSigma > maxPeakSigma) {
                maxPeakSigma = peakSigma;
            }
        }

        return maxPeakSigma;
    }

    private void adjustCostToTallestContourPeak1(
        Map<Integer, List<CurvatureScaleSpaceContour>> bestMatches1,
        Map<Integer, Integer> bestI2I1, TreeMap<Double, Set<Integer>> bestCosts,
        Map<Integer, Double> i2CostMap) {

        float maxPeakSigma = findMaxSigmaOfFirstPeaks(bestMatches1);

        TreeMap<Double, Set<Integer>> bestCostsUpdated = new TreeMap<Double, Set<Integer>>();

        for (Map.Entry<Double, Set<Integer>> entry : bestCosts.entrySet()) {

            double cost = entry.getKey().doubleValue();

            for (Integer index2 : entry.getValue()) {

                Integer index1 = bestI2I1.get(index2);

                float peak = bestMatches1.get(index1).get(0).getPeakSigma();

                double penalty = maxPeakSigma - peak;

                Double updatedCost = Double.valueOf(cost + penalty);

                Set<Integer> set = bestCostsUpdated.get(entry.getKey());
                if (set == null) {
                    set = new HashSet<Integer>();
                }
                set.add(index2);
                bestCostsUpdated.put(updatedCost, set);

                i2CostMap.put(index2, updatedCost);
            }
        }

        bestCosts.clear();
        bestCosts.putAll(bestCostsUpdated);
    }
}
