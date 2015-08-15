package algorithms.imageProcessing;

import algorithms.MultiArrayMergeSort;
import algorithms.util.PairIntArray;
import algorithms.util.PairInt;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.TreeMap;

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
    
    protected TransformationParameters bestFittingParameters = null;
      
    public CurvatureScaleSpaceInflectionMapper(ImageExt image1, ImageExt image2) {
        
        super(image1, image2);
        
    }
    
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
        
        boolean alreadySorted = true;
        
        for (int i1 = 0; i1 < contours1.size(); ++i1) {
            
            List<CurvatureScaleSpaceContour> list1 = contours1.get(i1);
            
            double minCost = Double.MAX_VALUE;
            List<CurvatureScaleSpaceContour> bestM1 = null;
            List<CurvatureScaleSpaceContour> bestM2 = null;
            double bestScale = 1;
            double bestCost = Double.MAX_VALUE;
            
            for (int i2 = 0; i2 < contours2.size(); ++i2) {
                
                List<CurvatureScaleSpaceContour> list2 = contours2.get(i2);
                
                CurvatureScaleSpaceContourMatcher matcher = 
                    new CurvatureScaleSpaceContourMatcher(list1, list2, 
                    alreadySorted);
                
                matcher.matchContours();
                
                List<CurvatureScaleSpaceContour> m1 = matcher.getSolutionMatchedContours1();
                List<CurvatureScaleSpaceContour> m2 = matcher.getSolutionMatchedContours2();
                if (m1 == null || m2 == null) {
                    continue;
                }
                assert(m1.size() == m2.size());
                
                double cost = matcher.getSolvedCost();
                
                if (cost < minCost) {
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
                extract(bestM1, xy1, weights1);
                extract(bestM2, xy2, weights2);
                
                if (xy1.getN() < 3) {
                    continue;
                }
                
                MatchedPointsTransformationCalculator tc = new 
                    MatchedPointsTransformationCalculator();
                
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
        
        for (Entry<Double, Set<Integer>> entry : bestCosts.entrySet()) {
                        
            Set<Integer> indexes = entry.getValue();
                                    
            for (Integer key : indexes) {
                
                Set<Integer> similar = new HashSet<Integer>();
                
                TransformationParameters params = bestParams.get(key);
                
                if (params == null) {
                    continue;
                }
                     
                similar.add(key);
                
                for (int j = 0; j < contours1.size(); ++j) {
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
    
    public TransformationParameters createEuclideanTransformationImpl() {
        
        if (bestFittingParameters == null) {
            return null;
        }
        
        MatchedPointsTransformationCalculator tc = new 
            MatchedPointsTransformationCalculator();
        
        if (debug) {
            tc.useDebugMode();
        }
                
        //TODO: temporarily disabling the refinement while fixing PointMatcher
        if (doRefineTransformations) {
            
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
            }
        }
        
        return bestFittingParameters;
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

    private void extract(List<CurvatureScaleSpaceContour> contours, 
        PairIntArray outputXY, List<Float> outputSigmaWeights) {
        
        float sumSigma = 0;
        
        for (int i = 0; i < contours.size(); i++) {
    
            CurvatureScaleSpaceContour c = contours.get(i);
        
            for (int j = 0; j < c.getPeakDetails().length; j++) {
                
                CurvatureScaleSpaceImagePoint spaceImagePoint = 
                    c.getPeakDetails()[j];
                
                int x = spaceImagePoint.getXCoord() + offsetImageX1;
                int y = spaceImagePoint.getYCoord() + offsetImageY1;
                
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

}
