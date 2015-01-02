package algorithms.imageProcessing;

import algorithms.MultiArrayMergeSort;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import java.util.List;

/**
 * a temporary version to explore matching inflection points for edges
 * that have been constructed with "outdoor mode" and are not closed
 * curves.
 * 
 * @author nichole
 */
public final class CurvatureScaleSpaceInflectionMapperForOpenCurves 
extends AbstractCurvatureScaleSpaceInflectionMapper {
        
    public CurvatureScaleSpaceInflectionMapperForOpenCurves(GreyscaleImage 
        image1, GreyscaleImage image2) {
        
        super(image1, image2);
        
        useOutdoorMode = true;
    }
    
    PairIntArray[] createUnmatchedXYFromContourPeaks() {
            
        initialize();
        
        if (matchedXY1 != null) {
            return null;
        }
        
        /*
        because the curves are not closed, the scale space images are not easy
        to match, so the peaks are extracted as "points of interest"
        and the PointMatcher class is used to match the points.
        */
        
        //TODO: find robust way to adjust this:        
        int sigmaLimit = 10;
        
        PairIntArray xyPeaks1 = new PairIntArray();
        
        PairIntArray xyPeaks2 = new PairIntArray();

        for (int i = 0; i < contours1.size(); i++) {
            
            CurvatureScaleSpaceContour c1 = contours1.get(i);
            
            if (c1.getPeakSigma() < sigmaLimit) {
                continue;
            }
            
            int edgeNumber = c1.getEdgeNumber();
            
            CurvatureScaleSpaceImagePoint[] peakDetails = c1.getPeakDetails();
            
            for (int j = 0; j < peakDetails.length; j++) {
                
                CurvatureScaleSpaceImagePoint p = peakDetails[j];
                
                log.info(String.format(
                    "1: (%d, %d) sigma=%f edgeLength=%d edgeNumber=%d", 
                        p.getXCoord(), p.getYCoord(), p.getSigma(),
                        edges1.get(edgeNumber).getN(), edgeNumber));
                
                xyPeaks1.add(p.getXCoord(), p.getYCoord());
            }
        }
        
        for (int i = 0; i < contours2.size(); i++) {
            
            CurvatureScaleSpaceContour c2 = contours2.get(i);
            
            if (c2.getPeakSigma() < sigmaLimit) {
                continue;
            }
            
            int edgeNumber = c2.getEdgeNumber();
            
            CurvatureScaleSpaceImagePoint[] peakDetails = c2.getPeakDetails();
            
            for (int j = 0; j < peakDetails.length; j++) {
                
                CurvatureScaleSpaceImagePoint p = peakDetails[j];
                
                log.info(String.format(
                    "2: (%d, %d) sigma=%f edgeLength=%d edgeNumber=%d", 
                        p.getXCoord(), p.getYCoord(), p.getSigma(),
                        edges2.get(edgeNumber).getN(), edgeNumber));
                
                xyPeaks2.add(p.getXCoord(), p.getYCoord());
            }
        }
        
        return new PairIntArray[]{xyPeaks1, xyPeaks2};
    }
        
    protected void createMatchedPointArraysFromContourPeaks() {
                
        PairIntArray[] xyPeaks = createUnmatchedXYFromContourPeaks();
        
        if (xyPeaks != null) {
            return;
        }
        
        PairIntArray xyPeaks1 = xyPeaks[0];
        
        PairIntArray xyPeaks2 = xyPeaks[1];
        
        PointMatcher pointMatcher = new PointMatcher();
        
        TransformationPointFit fit = 
            pointMatcher.calculateRoughTransformationForUnmatched(
                xyPeaks1, xyPeaks2, 
                (image1OriginalWidth >> 1), 
                (image1OriginalHeight >> 1));
        
        log.info("FIT: " + fit.toString());
        
        // ===== store the matched points and associated data ======
        /*
        PairIntArray xy1 = new PairIntArray();
        PairIntArray xy2 = new PairIntArray();
        List<Float> weights1 = new ArrayList<Float>();
        List<Float> weights2 = new ArrayList<Float>();
     
        double sumS1 = 0;
        double sumS2 = 0;
        
        List<Integer> matchedE1Idxs = new ArrayList<Integer>();
        List<Integer> matchedE2Idxs = new ArrayList<Integer>();
        matchedXY1ByEdgeInOrigRefFrame = new HashMap<Integer, PairIntArray>();
        matchedXY2ByEdgeInOrigRefFrame = new HashMap<Integer, PairIntArray>();
        matchedXY1ByEdgeWeights = new HashMap<Integer, List<Float> >();
        matchedXY2ByEdgeWeights = new HashMap<Integer, List<Float> >();
        
        for (int i = 0; i < transAppliedTo1.size(); i++) {
                        
            CurvatureScaleSpaceContour c1 = transAppliedTo1.get(i);
            CurvatureScaleSpaceContour c2 = transAppliedTo2.get(i);
            
            Integer e1Index = Integer.valueOf(c1.getEdgeNumber());
            Integer e2Index = Integer.valueOf(c2.getEdgeNumber());
            if (matchedE1Idxs.contains(e1Index)) {
                if (!matchedE2Idxs.contains(e2Index)) {
                    throw new IllegalStateException(
                    "inconsistency in matched edges for matched contours");
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
                s1.append(String.format("CONTOUR PEAK1: (%f, %f)", 
                    c1.getPeakSigma(), c1.getPeakScaleFreeLength()));
                s2.append(String.format("CONTOUR PEAK2: (%f, %f)", 
                    c2.getPeakSigma(), c2.getPeakScaleFreeLength()));
            }
        }

        if (xy1.getN() < 3) {
            throw new IllegalStateException("need at least 3 points");
        }
        
        if (debug) {
            log.info("offsetImgX1=" + offsetImageX1 
                + " offsetImgY1=" + offsetImageY1
                + "\noffsetImgX2=" + offsetImageX2 
                + " offsetImgY2=" + offsetImageY2
            );
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
            double tmp = weights1.get(i).floatValue()/sumS1;
            matchedXY1Weights[i] = Float.valueOf((float)tmp);
        }
        for (int i = 0; i < weights2.size(); i++) {
            double tmp = weights2.get(i).floatValue()/sumS2;
            matchedXY2Weights[i] = Float.valueOf((float)tmp);
        }
        */
    }

    @Override
    protected List<PairIntArray> getEdges(CurvatureScaleSpaceImageMaker imgMaker) {
        return imgMaker.getEdges();
    }
    
    @Override
    public TransformationParameters createEuclideanTransformationImpl() {
        
         throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public PairInt[] getMatchedEdgesIndexes() {
        throw new UnsupportedOperationException("Not supported yet."); 
    }
    
}
