package algorithms.imageProcessing;

import algorithms.util.PairIntArray;
import algorithms.util.PairInt;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

/**
 * class to map contours from one image to another and return the
 * matched inflection points and a transformation matrix that can
 * be applied to the first to put it into the frame of the second.
 * 
 * The algorithm used for matching scale space image contours is documented in 
 * CurvatureScaleSpaceContourMatcher
 * @see algorithms.imageProcessing.CurvatureScaleSpaceContourMatcher
 * 
 * @author nichole
 */
public final class CurvatureScaleSpaceInflectionMapper extends 
    AbstractCurvatureScaleSpaceInflectionMapper {
        
    private List<CurvatureScaleSpaceContour> matchedContours1 = new 
        ArrayList<CurvatureScaleSpaceContour>();
    
    private List<CurvatureScaleSpaceContour> matchedContours2 = new 
        ArrayList<CurvatureScaleSpaceContour>();
    
    /**
     * scale derived from matching contours.  it's not necessarily the same
     * as the final scale returned in transformation solutions, but it should
     * be close;
     */
    private double matchedScale = 1;
      
    public CurvatureScaleSpaceInflectionMapper(ImageExt image1, ImageExt image2) {
        
        super(image1, image2);
        
    }
    
    protected void createMatchedPointArraysFromContourPeaks() {
        
        if (matchedXY1 != null) {
            return;
        }
        
        /**
         * TODO:
         * change to match contours from one edge against contours of another
         * edge rather than all contours at once.
         *
         */
        CurvatureScaleSpaceContourMatcher matcher = new CurvatureScaleSpaceContourMatcher();
        matcher.matchContours(contours1, contours2);
        List<CurvatureScaleSpaceContour> transAppliedTo1 = matcher.getSolutionMatchedContours1();
        List<CurvatureScaleSpaceContour> transAppliedTo2 = matcher.getSolutionMatchedContours2();
        if (transAppliedTo1 == null || transAppliedTo2 == null) {
            return;
        }
        if (transAppliedTo1.size() != transAppliedTo2.size()) {
            throw new IllegalStateException("contour matcher should have same number of contours in both lists");
        }
        matchedContours1.addAll(transAppliedTo1);
        matchedContours2.addAll(transAppliedTo2);
        matchedScale = matcher.getSolvedScale();
        log.info("Contour matcher solution scale=" + matcher.getSolvedScale());
        log.info("Contour matcher solution shift=" + matcher.getSolvedShift());
        log.info("Contour matcher solution cost=" + matcher.getSolvedCost());
        PairIntArray xy1 = new PairIntArray(transAppliedTo1.size());
        PairIntArray xy2 = new PairIntArray(transAppliedTo1.size());
        List<Float> weights1 = new ArrayList<Float>();
        List<Float> weights2 = new ArrayList<Float>();
        double sumS1 = 0;
        double sumS2 = 0;
        List<Integer> matchedE1Idxs = new ArrayList<Integer>();
        List<Integer> matchedE2Idxs = new ArrayList<Integer>();
        matchedXY1ByEdgeInOrigRefFrame = new HashMap<Integer, PairIntArray>();
        matchedXY2ByEdgeInOrigRefFrame = new HashMap<Integer, PairIntArray>();
        matchedXY1ByEdgeWeights = new HashMap<Integer, List<Float>>();
        matchedXY2ByEdgeWeights = new HashMap<Integer, List<Float>>();
        for (int i = 0; i < transAppliedTo1.size(); i++) {
            CurvatureScaleSpaceContour c1 = transAppliedTo1.get(i);
            CurvatureScaleSpaceContour c2 = transAppliedTo2.get(i);
            Integer e1Index = Integer.valueOf(c1.getEdgeNumber());
            Integer e2Index = Integer.valueOf(c2.getEdgeNumber());
            if (matchedE1Idxs.contains(e1Index)) {
                if (!matchedE2Idxs.contains(e2Index)) {
                    throw new IllegalStateException("inconsistency in matched edges for matched contours");
                }
            } else {
                matchedE1Idxs.add(e1Index);
                matchedE2Idxs.add(e2Index);
                matchedXY1ByEdgeInOrigRefFrame.put(e1Index, new PairIntArray());
                matchedXY2ByEdgeInOrigRefFrame.put(e2Index, new PairIntArray());
                matchedXY1ByEdgeWeights.put(e1Index, new ArrayList<Float>());
                matchedXY2ByEdgeWeights.put(e2Index, new ArrayList<Float>());
            }
            float sigma1 = c1.getPeakSigma();
            float sigma2 = c2.getPeakSigma();
            StringBuilder s1 = new StringBuilder();
            StringBuilder s2 = new StringBuilder();
            if (debug) {
                s1.append(String.format("CONTOUR PEAK1: (%f, %f)", c1.getPeakSigma(), c1.getPeakScaleFreeLength()));
                s2.append(String.format("CONTOUR PEAK2: (%f, %f)", c2.getPeakSigma(), c2.getPeakScaleFreeLength()));
            }
            // the contours extracted from scale space images using a factor of
            // 2^(1/8) for recursive convolution tend to not have a single
            // peak, so the correction here for the single peak case is not
            // usually needed.  for that rare case, the avg of the other peak
            // is stored instead of both points
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
                    transAppliedTo2.set(i, c2);
                } else if (c2.getPeakDetails().length == 1) {
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
                    transAppliedTo1.set(i, c1);
                }
            }
            for (int j = 0; j < c1.getPeakDetails().length; j++) {
                CurvatureScaleSpaceImagePoint spaceImagePoint = 
                    c1.getPeakDetails()[j];
                int x = spaceImagePoint.getXCoord() + offsetImageX1;
                int y = spaceImagePoint.getYCoord() + offsetImageY1;
                xy1.add(x, y);
                weights1.add(Float.valueOf(sigma1));
                sumS1 += sigma1;
                matchedXY1ByEdgeInOrigRefFrame.get(e1Index)
                    .add(x + offsetImageX1, y + offsetImageY1);
                matchedXY1ByEdgeWeights.get(e1Index).add(Float.valueOf(sigma1));
                if (debug) {
                    s1.append(String.format(" (%d, %d)", x, y));
                }
                spaceImagePoint = c2.getPeakDetails()[j];
                x = spaceImagePoint.getXCoord() + offsetImageX2;
                y = spaceImagePoint.getYCoord() + offsetImageY2;
                xy2.add(x, y);
                weights2.add(Float.valueOf(sigma2));
                sumS2 += sigma2;
                matchedXY2ByEdgeInOrigRefFrame.get(e2Index)
                    .add(x + offsetImageX2, y + offsetImageX2);
                matchedXY2ByEdgeWeights.get(e2Index).add(Float.valueOf(sigma2));
                if (debug) {
                    s2.append(String.format(" (%d, %d)", x, y));
                }
            }
            if (debug) {
                log.info(s1.toString());
                log.info(s2.toString());
            }
        }
        if (xy1.getN() < 3) {
            throw new IllegalStateException("need at least 3 points");
        }
        if (debug) {
            log.info("offsetImgX1=" + offsetImageX1 + " offsetImgY1=" 
                + offsetImageY1 + "\noffsetImgX2=" + offsetImageX2 
                + " offsetImgY2=" + offsetImageY2);
        }
        matchedEdge1Indexes = new int[matchedE1Idxs.size()];
        matchedEdge2Indexes = new int[matchedE2Idxs.size()];
        for (int i = 0; i < matchedE1Idxs.size(); i++) {
            int e1Idx = matchedE1Idxs.get(i).intValue();
            int e2Idx = matchedE2Idxs.get(i).intValue();
            matchedEdge1Indexes[i] = e1Idx;
            matchedEdge2Indexes[i] = e2Idx;
        }
        matchedXY1 = xy1;
        matchedXY2 = xy2;
        matchedXY1Weights = new float[weights1.size()];
        matchedXY2Weights = new float[weights2.size()];
        for (int i = 0; i < weights1.size(); i++) {
            double tmp = weights1.get(i).floatValue() / sumS1;
            matchedXY1Weights[i] = Float.valueOf((float) tmp);
        }
        for (int i = 0; i < weights2.size(); i++) {
            double tmp = weights2.get(i).floatValue() / sumS2;
            matchedXY2Weights[i] = Float.valueOf((float) tmp);
        }
    }
    
    public TransformationParameters createEuclideanTransformationImpl() {
        
        if (matchedXY1.getN() < 3) {
            throw new IllegalStateException("need at least 3 points");
        }
        
        MatchedPointsTransformationCalculator tc = new 
            MatchedPointsTransformationCalculator();
        
        if (debug) {
            tc.useDebugMode();
        }
        
        int centroidX1 = image1OriginalWidth >> 1;
        int centroidY1 = image1OriginalHeight >> 1;
        int centroidX2 = image2OriginalWidth >> 1;
        int centroidY2 = image2OriginalHeight >> 1;
        TransformationParameters params = null;
        // if scale < 1, we have to swap the order of datasets to avoid
        // numerical errors in some of the methods that are the result of
        // dividing by a small number
        boolean reverseDatasetOrder = matchedScale < 1.0;
        if (reverseDatasetOrder) {
            params = tc.calulateEuclideanGivenScale(1. / matchedScale, 
                matchedXY2, matchedXY2Weights, matchedXY1, matchedXY1Weights, 
                centroidX2, centroidY2);
        } else {
            params = tc.calulateEuclideanGivenScale(matchedScale, matchedXY1, 
                matchedXY1Weights, matchedXY2, matchedXY2Weights, centroidX1, 
                centroidY1);
        }
        if (params == null) {
            return null;
        }
        
        //TODO: temporarily disabling the refinement while fixing PointMatcher
        if (doRefineTransformations) {
            //TODO: this needs to add scale changes to the refinement method
            PairIntArray[] set1 = getMatchedEdges1InOriginalReferenceFrameArray();
            PairIntArray[] set2 = getMatchedEdges2InOriginalReferenceFrameArray();
            EdgeMatcher matcher = new EdgeMatcher();
            TransformationPointFit fit2 = null;
            if (reverseDatasetOrder) {
                fit2 = matcher.refineTransformation(set2, set1, params);
            } else {
                fit2 = matcher.refineTransformation(set1, set2, params);
            }
            
            if (fit2 != null) {
                log.info("FINAL:\n" + fit2.toString());
                params = fit2.getParameters();
            }
        }
        
        if (reverseDatasetOrder) {
            params = tc.swapReferenceFrames(params);            
        }
        return params;
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
    
    @Override
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

    @Override
    protected List<PairIntArray> getEdges(CurvatureScaleSpaceImageMaker imgMaker) {
        
        return imgMaker.getClosedCurves();
    }

}
