package algorithms.imageProcessing;

import algorithms.MultiArrayMergeSort;
import algorithms.QuickSort;
import algorithms.misc.MiscDebug;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.logging.Logger;

/**
 *
 * @author nichole
 */
public abstract class AbstractCurvatureScaleSpaceInflectionMapper implements 
    ICurvatureScaleSpaceInflectionMapper {
    
    protected final Logger log = Logger.getLogger(this.getClass().getName());
    protected boolean debug = false;
    protected boolean useLineDrawingMode = false;
    protected boolean doRefineTransformations = false;
    protected boolean initialized = false;
    
    protected List<PairIntArray> edges1 = null;
    protected List<PairIntArray> edges2 = null;
    protected List<List<CurvatureScaleSpaceContour>> contourLists1 = new ArrayList<List<CurvatureScaleSpaceContour>>();
    protected List<List<CurvatureScaleSpaceContour>> contourLists2 = new ArrayList<List<CurvatureScaleSpaceContour>>();
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
    protected int offsetImageX1 = 0;
    protected int offsetImageY1 = 0;
    protected int offsetImageX2 = 0;
    protected int offsetImageY2 = 0;
    protected boolean useOutdoorMode = false;

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
    
    public AbstractCurvatureScaleSpaceInflectionMapper() {
    }

    @Override
    public void useOutdoorMode() {
        useOutdoorMode = true;
    }

    @Override
    public void useLineDrawingLineMode() {
        this.useLineDrawingMode = true;
    }

    @Override
    public void setToRefineTransformations() {
        doRefineTransformations = true;
    }

    @Override
    public void useDebugMode() {
        debug = true;
    }

    protected abstract void createEdges1();
    
    protected abstract void createEdges2();
    
    @Override
    public void initialize() {
        if (initialized) {
            return;
        }
        initialized = true;
        
        createEdges1();
        
        populateContours(edges1, contourLists1);
                
        /*
        note that when modifying the contour lists in any way, one has to
        maintain decreasing order by sigma and when sigma is equal, the
        order must be by increasing scale free parameter.
        two of the search methods in the matcher depend upon those properties.
         */
        
        createEdges2();
        
        populateContours(edges2, contourLists2);
        
    }

//TMP DEBUGGING
public GreyscaleImage debugImg1 = null;
public GreyscaleImage debugImg2 = null;
     
    protected void createMatchedPointArraysFromContourPeaks() {
        
        if (matchedXY1 != null) {
            return;
        }
        
        Map<Integer, List<CurvatureScaleSpaceContour>> bestMatches1 = new
            HashMap<Integer, List<CurvatureScaleSpaceContour>>();
        
        Map<Integer, List<CurvatureScaleSpaceContour>> bestMatchesTo1 = new
            HashMap<Integer, List<CurvatureScaleSpaceContour>>();
        
        Map<Integer, Float> bestScales = new HashMap<Integer, Float>();
        
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
        //i1=4 i2=3
        boolean alreadySorted = true;
        
        for (int i1 = 0; i1 < contourLists1.size(); ++i1) {
            
            List<CurvatureScaleSpaceContour> contours1 = contourLists1.get(i1);
            
            double minCost = Double.MAX_VALUE;
            List<CurvatureScaleSpaceContour> bestM1 = null;
            List<CurvatureScaleSpaceContour> bestM2 = null;
            double bestScale = 1;
            double bestCost = Double.MAX_VALUE;
if (debug && (debugImg1 != null)){           
Image img3 = new Image(debugImg1.getWidth(), debugImg1.getHeight());
for (int j = 0; j < edges1.get(i1).getN(); ++j) {
    int x = edges1.get(i1).getX(j);
    int y = edges1.get(i1).getY(j);
    /*if (i > 0) {
        x += xOffset;
        y += yOffset;
    }*/
    if (j == 0 || (j == (edges1.get(i1).getN() - 1))) {
        ImageIOHelper.addPointToImage(x, y, img3, 0, 200, 150, 0);
    } else {
        ImageIOHelper.addPointToImage(x, y, img3, 0, 255, 0, 0);
    }
}
MiscDebug.writeImageCopy(img3, "edge1_" + i1 + "_.png"); 
        }
            for (int i2 = 0; i2 < contourLists2.size(); ++i2) {
                
                List<CurvatureScaleSpaceContour> contours2 = contourLists2.get(i2);
                
                CSSContourMatcherWrapper matcher = 
                    new CSSContourMatcherWrapper(contours1, contours2, 
                    alreadySorted);
if (debug && (i1 == 0) && (debugImg2 != null)){               
Image img3 = new Image(debugImg2.getWidth(), debugImg2.getHeight());
for (int j = 0; j < edges2.get(i2).getN(); ++j) {
    int x = edges2.get(i2).getX(j);
    int y = edges2.get(i2).getY(j);
    /*if (i > 0) {
        x += xOffset;
        y += yOffset;
    }*/
    if (j == 0 || (j == (edges2.get(i2).getN() - 1))) {
        ImageIOHelper.addPointToImage(x, y, img3, 0, 200, 150, 0);
    } else {
        ImageIOHelper.addPointToImage(x, y, img3, 0, 255, 0, 0);
    }
}
MiscDebug.writeImageCopy(img3, "edge2_" + i2 + "_.png");     
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
extract(m1, txy1, tweights1, offsetImageX1, offsetImageY1);
extract(m2, txy2, tweights2, offsetImageX2, offsetImageY2);
MiscDebug.writeImage(txy1, ImageIOHelper.convertImage(debugImg1), "check_1_xy_edge_" + edgeIdx1 + "_" + flNumber);
MiscDebug.writeImage(txy2, ImageIOHelper.convertImage(debugImg2), "check_2_xy_edge_" + edgeIdx2 + "_" + flNumber);
MiscDebug.debugPlot(contours1, ImageIOHelper.convertImage(debugImg1), offsetImageX1, offsetImageY1, "1_edge_" + edgeIdx1 + "_" + String.valueOf(flNumber));
MiscDebug.debugPlot(contours2, ImageIOHelper.convertImage(debugImg2), offsetImageX2, offsetImageY2, "2_edge_" + edgeIdx1 + "_" + String.valueOf(flNumber));
int z = 1;
} catch (IOException ex) {
    
}               
*/                

                if ((cost < minCost) || ((bestM1 != null) && (m1.size() > 2) && (bestM1.size() < 3))) {
                    minCost = cost;
                    bestM1 = m1;
                    bestM2 = m2;
                    bestScale = matcher.getSolvedScale();
                    bestCost = cost;
                }
            }
            
            if (bestM1 != null) {
                
                // calculate the implied transformation from these matched points
                
                correctPeaks(bestM1, bestM2);
                
                PairIntArray xy1 = new PairIntArray(bestM1.size());
                PairIntArray xy2 = new PairIntArray(bestM2.size());
        
                List<Float> weights1 = new ArrayList<Float>();
                List<Float> weights2 = new ArrayList<Float>();
        
                //xy1 and xy2 have the image offsets added
                extract(bestM1, xy1, weights1, offsetImageX1, offsetImageY1);
                extract(bestM2, xy2, weights2, offsetImageX2, offsetImageY2);
                
                if (xy1.getN() < 3) {
                    continue;
                }
/*                
try {
    MiscDebug.writeImage(xy1, (ImageExt)image1.copyImage(),
        "check_1_xy_" + MiscDebug.getCurrentTimeFormatted());
} catch (IOException ex) {
    Logger.getLogger(CurvatureScaleSpaceInflectionMapper.class.getName()).log(Level.SEVERE, null, ex);
}
*/                
                MatchedPointsTransformationCalculator tc = new 
                    MatchedPointsTransformationCalculator();
                
                /*
                the xy1, xy2 coordinates are w.r.t. the original image coordinate
                reference frame (the offsets have been added back in).
                
                */
                
                int centroidX1 = 0;
                int centroidY1 = 0;
                int centroidX2 = 0;
                int centroidY2 = 0;
                
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
    MiscDebug.writeImage(xy1, (ImageExt)image1.copyImage(),
        "check_1_xy_" + flNumber);
    MiscDebug.writeImage(xy2, (ImageExt)image2.copyImage(),
        "check_2_xy_" + flNumber);
    Transformer transformer = new Transformer();
    PairIntArray xy1Tr = transformer.applyTransformation(params, xy1);
    MiscDebug.writeImage(xy1Tr, (ImageExt)image2.copyImage(),
        "check_1_xy_tr_" + flNumber);
} catch (IOException ex) {
    Logger.getLogger(CurvatureScaleSpaceInflectionMapper.class.getName()).log(Level.SEVERE, null, ex);
}
*/
                Integer key = Integer.valueOf(i1);
                
                bestMatches1.put(key, bestM1);
                bestMatchesTo1.put(key, bestM2);
                bestScales.put(key, Double.valueOf(bestScale).floatValue());
                bestParams.put(key, params);
                bestMatchesXY1.put(key, xy1);
                bestMatchesXY2.put(key, xy2);
                bestMatchesXYWeights1.put(key, weights1);
                bestMatchesXYWeights2.put(key, weights2);
                
                Double key2 = Double.valueOf(bestCost);
                if (!bestCosts.containsKey(key2)) {                    
                    bestCosts.put(key2, new HashSet<Integer>());
                }
                bestCosts.get(key2).add(key);
            }
        }
        
//TODO: need to discard the solutions with only 2 contours above
// or handle the comparison here.
// would prefer keep them in case that is the only solution
        
        
        /*
        compare the solutions, starting with the smallest cost solution.
        */
        int nTransformations = bestParams.size();
        
        /* calculate the highest number of similar transformations and the 
        lowest cost from those.
        store nSimilar, indexes, cost for each iteration
        */   
        int[] nSimilarSummary = new int[nTransformations];
        Integer[][] indexesSummary = new Integer[nTransformations][];
        double[] costsSummary = new double[nTransformations];
        int[] mainIndexSummary = new int[nTransformations];
        
        int count = 0;
        
        for (Map.Entry<Double, Set<Integer>> entry : bestCosts.entrySet()) {
                        
            Set<Integer> indexes = entry.getValue();
                                    
            for (Integer key : indexes) {
                
                Set<Integer> similar = new HashSet<Integer>();
                
                TransformationParameters params = bestParams.get(key);
                
                if (params == null) {
                    continue;
                }
                     
                similar.add(key);
                
                for (int j = 0; j < contourLists1.size(); ++j) {
                    if (j == key.intValue()) {
                        continue;
                    }                    
                    Integer key2 = Integer.valueOf(j);
                    TransformationParameters params2 = bestParams.get(key2);
                    if (params2 == null) {
                        continue;
                    }
                    if (Math.abs(params.getScale() - params2.getScale()) < 0.05) {
                        if (Math.abs(params.getRotationInDegrees() - params2.getRotationInDegrees()) < 10) {
                            if (Math.abs(params.getTranslationX() - params2.getTranslationX()) < 10) {
                                if (Math.abs(params.getTranslationY() - params2.getTranslationY()) < 10) {
                                    similar.add(key2);
                                }
                            }
                        }
                    }
                }
                nSimilarSummary[count] = similar.size();
                indexesSummary[count] = similar.toArray(new Integer[similar.size()]);
                costsSummary[count] = entry.getKey();
                mainIndexSummary[count] = key.intValue();
                count++;
            }
        }
        
        if (count == 0) {
            int z = 1;
            return;
        }
                
//==>TODO: change to make sure using unique matchings only in "indexes"
        
        /*
        nSimilarSummary[count]
        costsSummary[count]
        indexesSummary[count]
        mainIndexSummary[count]        
        */
        MultiArrayMergeSort.sortBy1stDescThen2ndAsc(nSimilarSummary, costsSummary,
            indexesSummary, mainIndexSummary);
       
        int nSimilar = nSimilarSummary[0];
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
                    
            Integer index = indexes[i];
                
            List<CurvatureScaleSpaceContour> m1 = bestMatches1.get(index);
            List<CurvatureScaleSpaceContour> m2 = bestMatchesTo1.get(index);
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

            matchedXY1ByEdgeInOrigRefFrame.put(e1Index,
                bestMatchesXY1.get(index));
            matchedXY2ByEdgeInOrigRefFrame.put(e2Index,
                bestMatchesXY2.get(index));
            matchedXY1ByEdgeWeights.put(e1Index,
                bestMatchesXYWeights1.get(index));
            matchedXY2ByEdgeWeights.put(e2Index,
                bestMatchesXYWeights2.get(index));
            
            matchedXY1.addAll(bestMatchesXY1.get(index));
            matchedXY2.addAll(bestMatchesXY2.get(index));
            
            matchedEdge1Indexes[i] = e1Index;
            matchedEdge2Indexes[i] = e2Index;
        }
        
    }
    
    
    public abstract TransformationParameters createEuclideanTransformationImpl();
    
    @Override
    public TransformationParameters createEuclideanTransformation() {
        
        //TODO:  no need to check and reverse contours in the init stage
        //   so remove those
        
        initialize();
        
        if (contourLists2.isEmpty() || contourLists1.isEmpty()) {
            return null;
        }
        
        createMatchedPointArraysFromContourPeaks();
        
        return createEuclideanTransformationImpl();
    }
    private void extract(List<CurvatureScaleSpaceContour> contours, 
        PairIntArray outputXY, List<Float> outputSigmaWeights,
        final int imageOffsetX, final int imageOffsetY) {
        
        float sumSigma = 0;
        
        for (int i = 0; i < contours.size(); i++) {
    
            CurvatureScaleSpaceContour c = contours.get(i);
        
            for (int j = 0; j < c.getPeakDetails().length; j++) {
                
                CurvatureScaleSpaceImagePoint spaceImagePoint = 
                    c.getPeakDetails()[j];
                
                int x = spaceImagePoint.getXCoord() + imageOffsetX;
                int y = spaceImagePoint.getYCoord() + imageOffsetY;
                
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
    
    protected void correctPeaks(List<CurvatureScaleSpaceContour> matched1, 
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

    @Override
    public PairIntArray getMatchedXY1() {
        return matchedXY1;
    }

    @Override
    public PairIntArray getMatchedXY2() {
        return matchedXY2;
    }

    @Override
    public float[] getMatchedXY1Weights() {
        return matchedXY1Weights;
    }

    @Override
    public float[] getMatchedXY2Weights() {
        return matchedXY2Weights;
    }

    @Override
    public List<List<CurvatureScaleSpaceContour>> getContours1() {
        return contourLists1;
    }

    @Override
    public List<List<CurvatureScaleSpaceContour>> getContours2() {
        return contourLists2;
    }

    protected List<PairIntArray> getEdges1() {
        return edges1;
    }

    protected List<PairIntArray> getEdges2() {
        return edges2;
    }

    protected List<PairIntArray> getEdges1InOriginalReferenceFrame() {
        
        List<PairIntArray> oe = new ArrayList<PairIntArray>();
        
        for (int i = 0; i < edges1.size(); i++) {
            PairIntArray edge = edges1.get(i).copy();
            for (int j = 0; j < edge.getN(); j++) {
                edge.set(j, edge.getX(j) + offsetImageX1, 
                    edge.getY(j) + offsetImageY1);
            }
            oe.add(edge);
        }
        return oe;
    }

    protected PairIntArray[] getEdges1InOriginalReferenceFrameArray() {
        PairIntArray[] oe = new PairIntArray[edges1.size()];
        for (int i = 0; i < edges1.size(); i++) {
            PairIntArray edge = edges1.get(i).copy();
            for (int j = 0; j < edge.getN(); j++) {
                edge.set(j, edge.getX(j) + offsetImageX1, 
                    edge.getY(j) + offsetImageY1);
            }
            oe[i] = edge;
        }
        return oe;
    }

    protected List<PairIntArray> getEdges2InOriginalReferenceFrame() {
        List<PairIntArray> oe = new ArrayList<PairIntArray>();
        for (int i = 0; i < edges2.size(); i++) {
            PairIntArray edge = edges2.get(i).copy();
            for (int j = 0; j < edge.getN(); j++) {
                edge.set(j, edge.getX(j) + offsetImageX2, 
                    edge.getY(j) + offsetImageY2);
            }
            oe.add(edge);
        }
        return oe;
    }

    protected PairIntArray[] getEdges2InOriginalReferenceFrameArray() {
        PairIntArray[] oe = new PairIntArray[edges2.size()];
        for (int i = 0; i < edges2.size(); i++) {
            PairIntArray edge = edges2.get(i).copy();
            for (int j = 0; j < edge.getN(); j++) {
                edge.set(j, edge.getX(j) + offsetImageX2, 
                    edge.getY(j) + offsetImageY2);
            }
            oe[i] = edge;
        }
        return oe;
    }

    protected PairIntArray[] getMatchedEdges1InOriginalReferenceFrameArray() {
        if (matchedEdge1Indexes == null) {
            return new PairIntArray[0];
        }
        PairIntArray[] oe = new PairIntArray[matchedEdge1Indexes.length];
        for (int i = 0; i < matchedEdge1Indexes.length; i++) {
            int eIdx = matchedEdge1Indexes[i];
            PairIntArray edge = edges1.get(eIdx).copy();
            for (int j = 0; j < edge.getN(); j++) {
                edge.set(j, edge.getX(j) + offsetImageX1, 
                    edge.getY(j) + offsetImageY1);
            }
            oe[i] = edge;
        }
        return oe;
    }

    protected PairIntArray[] getMatchedEdges2InOriginalReferenceFrameArray() {
        if (matchedEdge2Indexes == null) {
            return new PairIntArray[0];
        }
        PairIntArray[] oe = new PairIntArray[matchedEdge2Indexes.length];
        for (int i = 0; i < matchedEdge2Indexes.length; i++) {
            int eIdx = matchedEdge2Indexes[i];
            PairIntArray edge = edges2.get(eIdx).copy();
            for (int j = 0; j < edge.getN(); j++) {
                edge.set(j, edge.getX(j) + offsetImageX2, 
                    edge.getY(j) + offsetImageY2);
            }
            oe[i] = edge;
        }
        return oe;
    }
    
    int getOffsetImageX1() {
        return offsetImageX1;
    }

    int getOffsetImageY1() {
        return offsetImageY1;
    }

    int getOffsetImageX2() {
        return offsetImageX2;
    }

    int getOffsetImageY2() {
        return offsetImageY2;
    }

    protected abstract List<PairIntArray> getEdges(
        CurvatureScaleSpaceImageMaker imgMaker);

    private void populateContours(List<PairIntArray> edges, 
        List<List<CurvatureScaleSpaceContour>> contours) {
        
        CurvatureScaleSpaceCurvesMaker csscMaker = new CurvatureScaleSpaceCurvesMaker();
        
        // if use 2^(1/8) as a sigma factor should result in an error less than 10%
        // in determing the peak of a contour.  smaller factors have smaller
        // errors than that.
        float factor = (float)Math.pow(2, 1./32.);
        
        for (int i = 0; i < edges.size(); i++) {

            PairIntArray edge = edges.get(i);
      
            Map<Float, ScaleSpaceCurve> scaleSpaceMap = 
                csscMaker.createScaleSpaceMetricsForEdge(edge, factor,
                SIGMA.ONE, SIGMA.TWOHUNDREDANDFIFTYSIX);
           
            ScaleSpaceCurveImage scaleSpaceImage = 
                csscMaker.convertScaleSpaceMapToSparseImage(
                scaleSpaceMap, i, edge.getN());
            
                       
 try {
 String fileSuffix = "edge_" + i + "_" + MiscDebug.getCurrentTimeFormatted();
 MiscDebug.printScaleSpaceCurve(scaleSpaceImage, fileSuffix);
 int z = 1;
 } catch (IOException ex) {
 }
            
            ContourFinder contourFinder = new ContourFinder();

            List<CurvatureScaleSpaceContour> result = contourFinder.findContours(scaleSpaceImage, i);

            boolean reversed = contourFinder.reverseIfClockwise(result);

            if (reversed) {
                log.info("EDGES1: contour isCW=true");

                // these are extracted from contourFinder in order of decreasing
                // sigma already, so only need to be sorted if the list was
                // reversed
                Collections.sort(result, new DescendingSigmaComparator());
            }

            //MiscDebug.debugPlot(result, (ImageExt) image1.copyImage(), offsetImageX1, offsetImageY1,
            //    "_1_" + MiscDebug.getCurrentTimeFormatted());

            contours.add(result);
        }

        if (contours.isEmpty()) {
            log.info("no contours found in image 1");
            return;
        }
        if (edges.size() > 1) {
            Collections.sort(contours, new DescendingSigmaComparator2());
        }

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

}
