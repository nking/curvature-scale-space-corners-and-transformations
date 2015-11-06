package algorithms.imageProcessing;

import algorithms.QuickSort;
import algorithms.misc.MiscDebug;
import algorithms.misc.MiscMath;
import algorithms.util.PairInt;
import java.io.IOException;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.logging.Logger;

/**
 * class encapsulating the steps from scale calculation to matching corners
 * to make correspondence lists.
 * 
 * @author nichole
 */
public class FeatureMatcherWrapper {
    
    private final ImageExt img1;
    private final ImageExt img2;
    
    private GreyscaleImage gsImg1 = null;
    private GreyscaleImage gsImg2 = null;
    
    private GreyscaleImage gXY1 = null;
    private GreyscaleImage gXY2 = null;
    
    private Set<CornerRegion> cornerRegions1 = null;
    private Set<CornerRegion> cornerRegions2 = null;
    
    private GreyscaleImage theta1 = null;
    private GreyscaleImage theta2 = null;

    private enum State {
        DID_APPLY_HIST_EQ, COULD_NOT_DETERMINE_SCALE
    }
    private Set<State> stateSet = new HashSet<State>();
    
    private final boolean doDetermineScale;
    
    private final boolean debug;
    
    private final String debugTagPrefix;
    
    private TransformationParameters params = null;
    
    private float scaleSetByUser = Float.MIN_VALUE;
        
    private float scaleTol = 0.2f;
    
    private float rotationInRadiansTol = (float)(20. * Math.PI/180.);
    
    private int transXYTol = 70;
    
    private Logger log = Logger.getLogger(this.getClass().getName());
    
    public FeatureMatcherWrapper(ImageExt image1, ImageExt image2) {
        img1 = image1;
        img2 = image2;
        doDetermineScale = true;
        debug = false;
        debugTagPrefix = "";
    }
    
    public FeatureMatcherWrapper(ImageExt image1, ImageExt image2, 
        String debugTagPrefix) {
        img1 = image1;
        img2 = image2;
        doDetermineScale = true;
        debug = true;
        this.debugTagPrefix = debugTagPrefix;
    }
    
    public FeatureMatcherWrapper(ImageExt image1, ImageExt image2, float scale) {
        img1 = image1;
        img2 = image2;
        doDetermineScale = false;
        scaleSetByUser = scale;
        debug = false;
        debugTagPrefix = "";
    }
    
    public FeatureMatcherWrapper(ImageExt image1, ImageExt image2, float scale,
        String debugTagPrefix) {
        img1 = image1;
        img2 = image2;
        doDetermineScale = false;
        scaleSetByUser = scale;
        debug = true;
        this.debugTagPrefix = debugTagPrefix;
    }
    
    public FeatureMatcherWrapper(ImageExt image1, ImageExt image2, 
        TransformationParameters parameters) {
        img1 = image1;
        img2 = image2;
        doDetermineScale = false;
        params = parameters;
        debug = false;
        this.debugTagPrefix = "";
    }
    
    public FeatureMatcherWrapper(ImageExt image1, ImageExt image2, 
        TransformationParameters parameters, String debugTagPrefix) {
        img1 = image1;
        img2 = image2;
        doDetermineScale = false;
        params = parameters;
        debug = true;
        this.debugTagPrefix = debugTagPrefix;
    }
    
    public CorrespondenceList matchFeatures() throws IOException, NoSuchAlgorithmException {
        
        /*
        options:
            (1) determine scale 
                (a) match remaining points derived in scale calc.
                    if resulting set spans the intersection,
                    make and return the correspondence list
                    else, follow (2)
            (2) given scale
                (b) extract corner regions from greyscale image
                (3) use feature matcher w/ scale to make the correspondence list
        */
        
        CorrespondenceList cl = null;
        
        if (doDetermineScale) {
            cl = solveForScale();
            if (cl != null) {
                return cl;
            }
        }
        
        if (stateSet.contains(State.COULD_NOT_DETERMINE_SCALE)) {
            return null;
        }
                
        if (cl == null) {
            applyHistEqIfNeeded();
        }
        
        extractCornerRegions();
                
        if (params != null) {
            cl = findCorrespondence(params);
        } else {
            cl = findCorrespondence(this.scaleSetByUser);
        }
        
        if (debug) {
            Collection<PairInt> m1 = cl.getPoints1();
            Collection<PairInt> m2 = cl.getPoints2();

            MiscDebug.plotCorners(gsImg1.copyImage(), m1, debugTagPrefix + "_1_matched", 2);
            MiscDebug.plotCorners(gsImg2.copyImage(), m2, debugTagPrefix + "_2_matched", 2);
        }
        
        return cl;
    }
    
    private CorrespondenceList solveForScale() throws IOException, 
        NoSuchAlgorithmException {
        
        BlobScaleFinderWrapper scaleFinder = null;
            
        if (debug) {
            scaleFinder = new BlobScaleFinderWrapper(img1, img2, debugTagPrefix);
        } else {
            scaleFinder = new BlobScaleFinderWrapper(img1, img2);
        }
        
        params = scaleFinder.calculateScale();
        
        if (params == null) {
            stateSet.add(State.COULD_NOT_DETERMINE_SCALE);
            return null;
        }
         
        boolean didApplyHist = scaleFinder.img1Helper.didApplyHistEq();
        this.gsImg1 = scaleFinder.img1Helper.getGreyscaleImage().copyImage();
        this.gsImg2 = scaleFinder.img2Helper.getGreyscaleImage().copyImage();
        
        if (didApplyHist) {
            stateSet.add(State.DID_APPLY_HIST_EQ);
        }
        
        List<FeatureComparisonStat> stats = null;
        
        CorrespondenceList cl = null;
        
        // try to match the remaining points created in the scale finder
        if (scaleFinder.getAllCornerRegions1OfSolution() != null) {
            
            List<List<CornerRegion>> transformedFilteredC1 
                = new ArrayList<List<CornerRegion>>();
            List<List<CornerRegion>> filteredC1 = new ArrayList<List<CornerRegion>>();
            List<List<CornerRegion>> filteredC2 = new ArrayList<List<CornerRegion>>();
            
            FeatureMatcher.filterForIntersection(params, transXYTol,
                scaleFinder.getAllCornerRegions1OfSolution(),
                scaleFinder.getAllCornerRegions2OfSolution(),
                transformedFilteredC1, filteredC1, filteredC2);
            
            stats = matchRemainingBlobCornerPoints(scaleFinder, 
                transformedFilteredC1, filteredC1, filteredC2);
            
            if (statsCoverIntersection(stats, filteredC2)) {
                
                List<PairInt> matched1 = new ArrayList<PairInt>();
                List<PairInt> matched2 = new ArrayList<PairInt>();
                populateLists(stats, matched1, matched2);
                
                cl = new CorrespondenceList(params.getScale(), 
                    Math.round(params.getRotationInDegrees()), 
                    Math.round(params.getTranslationX()),
                    Math.round(params.getTranslationY()), 
                    Math.round(params.getStandardDeviations()[0]),
                    Math.round(params.getStandardDeviations()[2]),
                    Math.round(params.getStandardDeviations()[3]),
                    matched1, matched2);
                
                return cl;
            }
            
            cl = extractAndMatch(params, stats);
            
        } else {
            
            // solve for contours
            List<List<CurvatureScaleSpaceContour>> transformedFilteredC1 
                = new ArrayList<List<CurvatureScaleSpaceContour>>();
            List<List<CurvatureScaleSpaceContour>> filteredC1 = 
                new ArrayList<List<CurvatureScaleSpaceContour>>();
            List<List<CurvatureScaleSpaceContour>> filteredC2 = 
                new ArrayList<List<CurvatureScaleSpaceContour>>();
            
            FeatureMatcher.filterForIntersection2(params, transXYTol,
                scaleFinder.getAllContours1OfSolution(),
                scaleFinder.getAllContours2OfSolution(),
                transformedFilteredC1, filteredC1, filteredC2);
            
            stats = matchRemainingBlobContourPoints(scaleFinder, 
                transformedFilteredC1, filteredC1, filteredC2);
            
            if (statsCoverIntersection2(stats, filteredC2)) {
                
                List<PairInt> matched1 = new ArrayList<PairInt>();
                List<PairInt> matched2 = new ArrayList<PairInt>();
                populateLists(stats, matched1, matched2);
                
                cl = new CorrespondenceList(params.getScale(), 
                    Math.round(params.getRotationInDegrees()), 
                    Math.round(params.getTranslationX()),
                    Math.round(params.getTranslationY()), 
                    Math.round(params.getStandardDeviations()[0]),
                    Math.round(params.getStandardDeviations()[2]),
                    Math.round(params.getStandardDeviations()[3]),
                    matched1, matched2);
                
                return cl;
            }
            
            cl = extractAndMatch(params, stats);
        }
        
        if (cl != null) {
            return cl;
        }
               
        return cl;
    }
    
    private void applyHistEqIfNeeded() {
        
        if (stateSet.contains(State.DID_APPLY_HIST_EQ)) {
            return;
        }
        
        if (gsImg1 != null) {
            // gs images were set during scale calculation
            return;
        }
        
        this.gsImg1 = img1.copyToGreyscale();
        this.gsImg2 = img2.copyToGreyscale();
        
        ImageStatistics stats1 = ImageStatisticsHelper.examineImage(gsImg1, true);
        ImageStatistics stats2 = ImageStatisticsHelper.examineImage(gsImg2, true);
        
        boolean performHistEq = false;
        double median1DivMedian2 = stats1.getMedian()/stats2.getMedian();
        double meanDivMedian1 = stats1.getMean()/stats1.getMedian();
        double meanDivMedian2 = stats2.getMean()/stats2.getMedian();
        if (
            ((median1DivMedian2 > 1) && ((median1DivMedian2 - 1) > 0.2)) ||
            ((median1DivMedian2 < 1) && (median1DivMedian2 < 0.8))) {
            performHistEq = true;
        } else if (
            ((meanDivMedian1 > 1) && ((meanDivMedian1 - 1) > 0.2)) ||
            ((meanDivMedian1 < 1) && (meanDivMedian1 < 0.8))) {
            performHistEq = true;
        } else if (
            ((meanDivMedian2 > 1) && ((meanDivMedian2 - 1) > 0.2)) ||
            ((meanDivMedian2 < 1) && (meanDivMedian2 < 0.8))) {
            performHistEq = true;
        }
        if (performHistEq) {
            log.info("use histogram equalization on the greyscale images");
            HistogramEqualization hEq = new HistogramEqualization(gsImg1);
            hEq.applyFilter();
            hEq = new HistogramEqualization(gsImg2);
            hEq.applyFilter();
            stateSet.add(State.DID_APPLY_HIST_EQ);
        }

    }

    private void extractCornerRegions() {
        
        ImageProcessor imageProcessor = new ImageProcessor();
        imageProcessor.blur(gsImg1, SIGMA.ONE);
        
        CurvatureScaleSpaceCornerDetector detector = new
            CurvatureScaleSpaceCornerDetector(gsImg1);
        detector.doNotPerformHistogramEqualization();
        detector.findCorners();
        cornerRegions1 = detector.getEdgeCornerRegions(true);
        //cornerRegions1 = detector.getEdgeCornerRegionsInOriginalReferenceFrame(true);
        gXY1 = detector.getEdgeFilterProducts().getGradientXY();
        //GreyscaleImage img1Grey = gsImg1.copyImage();
        //imageProcessor.shrinkImage(img1Grey, 
        //    new int[]{gXY1.getXRelativeOffset(), gXY1.getYRelativeOffset(),
        //        gXY1.getWidth(), gXY1.getHeight()
        //    });
        
        theta1 = imageProcessor.computeTheta360(
            detector.getEdgeFilterProducts().getGradientX(), 
            detector.getEdgeFilterProducts().getGradientY());
        
        //-------
        
        imageProcessor.blur(gsImg2, SIGMA.ONE);
        
        detector = new
            CurvatureScaleSpaceCornerDetector(gsImg2);
        detector.doNotPerformHistogramEqualization();
        detector.findCorners();
        cornerRegions2 = detector.getEdgeCornerRegions(true);
        //cornerRegions2 = detector.getEdgeCornerRegionsInOriginalReferenceFrame(true);
        gXY2 = detector.getEdgeFilterProducts().getGradientXY();
        //GreyscaleImage img2Grey = gsImg2.copyImage();
        //imageProcessor.shrinkImage(img1Grey, 
        //    new int[]{gXY1.getXRelativeOffset(), gXY1.getYRelativeOffset(),
        //        gXY1.getWidth(), gXY1.getHeight()
        //    });
        
        theta2 = imageProcessor.computeTheta360(
            detector.getEdgeFilterProducts().getGradientX(), 
            detector.getEdgeFilterProducts().getGradientY());
    }

    private CorrespondenceList findCorrespondence(TransformationParameters parameters) {
        
        FeatureMatcher matcher = new FeatureMatcher();
        
        CorrespondenceList cl = matcher.findSimilarFeatures(gsImg1, gXY1, theta1,
            cornerRegions1.toArray(new CornerRegion[cornerRegions1.size()]),
            gsImg2, gXY2, theta2,
            cornerRegions2.toArray(new CornerRegion[cornerRegions2.size()]), 
            parameters, scaleTol, rotationInRadiansTol, transXYTol);

        return cl;
    }
    
    private CorrespondenceList findCorrespondence(float scale) {
        
        FeatureMatcher matcher = new FeatureMatcher();
        
        CorrespondenceList cl = matcher.findSimilarFeatures(gsImg1, gXY1, theta1,
            cornerRegions1.toArray(new CornerRegion[cornerRegions1.size()]),
            gsImg2, gXY2, theta2,
            cornerRegions2.toArray(new CornerRegion[cornerRegions2.size()]), 
            scale);

        return cl;
    }
    
    private List<FeatureComparisonStat> matchRemainingBlobCornerPoints(
        BlobScaleFinderWrapper scaleFinder, 
        List<List<CornerRegion>> filteredC1Transformed,
        List<List<CornerRegion>> filteredC1, List<List<CornerRegion>> filteredC2) {
        
        if (filteredC1Transformed.size() != filteredC1.size()) {
            throw new IllegalArgumentException("filteredC1Transformed and "
            + "filteredC1 are expected to be same size");
        }
        
        // use the association w/ tranformed blobs to make the matching faster

        List<FeatureComparisonStat> compStats = 
            new ArrayList<FeatureComparisonStat>(
                scaleFinder.getSolution().getComparisonStats());
        
        /*
        choose the best for each '1' and if a high quality exists, store
        it for further quality check (theta and intensity) then add to compStats
        
        for transformedblob1
            init storage for best match to blob1
            for blob2
                if centroid within tolerance,
                    use features to match untransformed blob1 corners to blob2 corners
                    (this is the curve matcher within corner matcher for combinations?)
                    if results are high quality and better than best,
                        assign it as best
             store best in map for blob1
        remove outliers by theta and by ssd
        return combined results
        */
        
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        
        Map<Integer, IntensityFeatureComparisonStats> index1Map 
            = new HashMap<Integer, IntensityFeatureComparisonStats>();
        
        int tolXY = this.transXYTol;
        
        for (int i1 = 0; i1 < filteredC1Transformed.size(); ++i1) {
            List<CornerRegion> trC1List = filteredC1Transformed.get(i1);
            if (trC1List.isEmpty()) {
                continue;
            }
            List<CornerRegion> c1List = filteredC1.get(i1);
            
            double[] xyCen1 = curveHelper.calculateXYCentroids0(trC1List);
            
            IntensityFeatureComparisonStats best = null;
            
            for (int i2 = 0; i2 < filteredC2.size(); ++i2) {
                List<CornerRegion> c2List = filteredC2.get(i2);
                if (c2List.isEmpty()) {
                    continue;
                }
                double[] xyCen2 = curveHelper.calculateXYCentroids0(c2List);
                
                double diffX = Math.abs(xyCen1[0] - xyCen2[0]);
                double diffY = Math.abs(xyCen1[1] - xyCen2[1]);
                if ((diffX > tolXY) || (diffY > tolXY)) {
                    continue;
                }
                
                ClosedCurveCornerMatcherWrapper mapper =
                    new ClosedCurveCornerMatcherWrapper();
                
                boolean matched = mapper.matchCorners(
                    scaleFinder.getSolutionFeatures1(), 
                    scaleFinder.getSolutionFeatures2(),
                    c1List, c2List, true, gsImg1, gsImg2);
                
                if (!matched) {
                    continue;
                }
                
                TransformationPair2 transformationPair = mapper.getSolution();
                transformationPair.setCornerListIndex1(i1);
                transformationPair.setCornerListIndex2(i2);
                
                TransformationParameters params2 = 
                    transformationPair.getTransformationParameters();
                
                if (params2 == null) {
                    continue;
                }
                
                List<FeatureComparisonStat> compStats2 = 
                    transformationPair.getNextCorner().getMatchedFeatureComparisonStats();
               
                FeatureMatcher.removeDiscrepantThetaDiff(compStats2, 
                    params.getRotationInDegrees());
                
                if (compStats2.isEmpty()) {
                    continue;
                }
                
                //FeatureMatcher.removeIntensityOutliers(compStats2);
                
                IntensityFeatureComparisonStats stats2 = new 
                    IntensityFeatureComparisonStats(i1, i2,
                    mapper.getSolvedCost(), params.getScale());
                stats2.addAll(compStats2);
                
int nm = stats2.getComparisonStats().size();
int x1 = stats2.getComparisonStats().get(0).getImg1Point().getX();
int y1 = stats2.getComparisonStats().get(0).getImg1Point().getY();
int x2 = stats2.getComparisonStats().get(0).getImg2Point().getX();
int y2 = stats2.getComparisonStats().get(0).getImg2Point().getY();

                int comp = -1;
                if (best != null) {
                    comp = stats2.compareTo(best);
                }
                if (comp == -1) {
                    best = stats2;
                }
            }
            
            if (best != null) {
                index1Map.put(Integer.valueOf(i1), best);
            }
        }
//[0],      
//wrong: [1],
                
        List<FeatureComparisonStat> add = new ArrayList<FeatureComparisonStat>();
        for (Entry<Integer, IntensityFeatureComparisonStats> entry : index1Map.entrySet()) {
            // make sure not already in compStats
            for (FeatureComparisonStat stat : entry.getValue().getComparisonStats()) {
                PairInt p1 = stat.getImg1Point();
                PairInt p2 = stat.getImg2Point();
                for (FeatureComparisonStat cStat : compStats) {
                    boolean found = false;
                    PairInt p1c = cStat.getImg1Point();
                    PairInt p2c = cStat.getImg2Point();
                    if (p1c.equals(p1) || p2c.equals(p2)) {
                        found = true;
                        break;
                    }
                    int diffX1 = Math.abs(p1c.getX() - p1.getX());
                    int diffY1 = Math.abs(p1c.getY() - p1.getY());
                    if ((diffX1 < 5) && (diffY1 < 5)) {
                        found = true;
                        break;
                    }
                    int diffX2 = Math.abs(p2c.getX() - p2.getX());
                    int diffY2 = Math.abs(p2c.getY() - p2.getY());
                    if ((diffX2 < 5) && (diffY2 < 5)) {
                        found = true;
                        break;
                    }
                    if (!found) {
                        add.add(stat);
                    }
                }
            }
        }
        
        compStats.addAll(add);
            
        return compStats;
    }
    
    private List<FeatureComparisonStat> matchRemainingBlobContourPoints(
        BlobScaleFinderWrapper scaleFinder, 
        List<List<CurvatureScaleSpaceContour>> transformedFilteredC1, 
        List<List<CurvatureScaleSpaceContour>> filteredC1, 
        List<List<CurvatureScaleSpaceContour>> filteredC2) {
        
        throw new UnsupportedOperationException("Not supported yet."); 
    }

    private boolean statsCoverIntersection2(List<FeatureComparisonStat> stats, 
        List<List<CurvatureScaleSpaceContour>> filteredC2) {
                
        /*
        dividing the range in filteredC2 by 2 in x and 2 in y and returning
        true if at least one point2 in stats is found in each division.             
        */
        int n = 0;
        for (List<CurvatureScaleSpaceContour> list : filteredC2) {
            n += list.size();
        }
        float[] xPoints = new float[n];
        float[] yPoints = new float[n];
        
        n = 0;
        for (List<CurvatureScaleSpaceContour> list : filteredC2) {
            for (CurvatureScaleSpaceContour cr : list) {
                float x = cr.getPeakDetails()[0].getXCoord();
                float y = cr.getPeakDetails()[0].getYCoord();
                xPoints[n] = x;
                yPoints[n] = y;
                ++n;
            }
        }
        
        return statsCoverIntersection(stats, xPoints, yPoints);
    }
    
    private boolean statsCoverIntersection(List<FeatureComparisonStat> stats, 
        List<List<CornerRegion>> filteredC2) {
                
        /*
        dividing the range in filteredC2 by 2 in x and 2 in y and returning
        true if at least one point2 in stats is found in each division.             
        */
        int n = 0;
        for (List<CornerRegion> list : filteredC2) {
            n += list.size();
        }
        float[] xPoints = new float[n];
        float[] yPoints = new float[n];
        
        n = 0;
        for (List<CornerRegion> list : filteredC2) {
            for (CornerRegion cr : list) {
                float x = cr.getX()[cr.getKMaxIdx()];
                float y = cr.getY()[cr.getKMaxIdx()];
                xPoints[n] = x;
                yPoints[n] = y;
                ++n;
            }
        }
        
        return statsCoverIntersection(stats, xPoints, yPoints);
    }
    
    private boolean statsCoverIntersection(List<FeatureComparisonStat> stats, 
        final float[] xPoints, final float[] yPoints) {
               
        int n = xPoints.length;
        
        QuickSort.sortBy1stThen2nd(xPoints, yPoints);
        
        float minX = xPoints[0];
        float maxX = xPoints[n - 1];
        float divX = (maxX + minX)/2.f;
       
        /*
        Finding y min and max within each of these division.  The reason
        for doing this separately from y min max over all of filteredC2 is that 
        the geometry of matchable points might not be rectangular.
        
            |          |          |
            |          |          |
            |          |          |
                 0           1
        */
        
        float[] yMin = new float[2];
        Arrays.fill(yMin, Float.MAX_VALUE);
        float[] yMax = new float[2];
        Arrays.fill(yMax, Float.MIN_VALUE);
        
        for (int i = 0; i < xPoints.length; ++i) {
            float x = xPoints[i];
            float y = yPoints[i];
            int cIdx = 1;
            if (x < divX) {
                cIdx = 0;
            }
            if (y < yMin[cIdx]) {
                yMin[cIdx] = y;
            }
            if (y > yMax[cIdx]) {
                yMax[cIdx] = y;
            }
        }
        
        float yDiv12 = (yMax[0] + yMin[0])/2.f;
        float yDiv03 = (yMax[1] + yMin[1])/2.f;
        
        /*
             2       3
        |        |        |
             1       0
        */
        
        int[] counts = new int[4];
        for (FeatureComparisonStat stat : stats) {
            PairInt p2 = stat.getImg2Point();
            int x = p2.getX();
            int y = p2.getY();
            if (x < divX) {
                if (y < yDiv12) {
                    counts[1]++;
                } else {
                    counts[2]++;
                }
            } else {
                if (y < yDiv03) {
                    counts[0]++;
                } else {
                    counts[3]++;
                }
            }
        }
        
        //TODO: revise while testing
        // check that there is at least 1 in each quadrant
        if ((counts[0] > 0) && (counts[1] > 0) && (counts[2] > 0) && 
            (counts[3] > 0)) {
            
            return true;
        }
        
        return false;
    }
    
    private void populateLists(List<FeatureComparisonStat> stats, 
        List<PairInt> matched1, List<PairInt> matched2) {
        
        for (FeatureComparisonStat stat : stats) {
            
            int x1 = stat.getImg1Point().getX();
            int y1 = stat.getImg1Point().getY();
            int x2 = stat.getImg2Point().getX();
            int y2 = stat.getImg2Point().getY();
            
            matched1.add(new PairInt(x1, y1));
            
            matched2.add(new PairInt(x2, y2));
        }
    }
    
    private CorrespondenceList extractAndMatch(
        TransformationParameters parameters, List<FeatureComparisonStat> stats) {
        
        extractCornerRegions();

        CorrespondenceList cl = findCorrespondence(parameters);
        
        if (cl != null) {
            // add stats in if not already present
            int z = 1;
        }
        
        return cl;
    }
    
}
