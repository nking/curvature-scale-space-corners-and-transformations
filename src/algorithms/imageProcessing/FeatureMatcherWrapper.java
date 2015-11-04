package algorithms.imageProcessing;

import algorithms.misc.MiscDebug;
import algorithms.util.PairInt;
import java.io.IOException;
import java.security.NoSuchAlgorithmException;
import java.util.Collection;
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

    private boolean didApplyHistEq = false;
    
    private final boolean doDetermineScale;
    
    private final boolean debug;
    
    private final String debugTagPrefix;
    
    private TransformationParameters params = null;
    
    private float scaleSetByUser = Float.MIN_VALUE;
    
    private float rotationInRadianSetByUser = Float.MIN_VALUE;
    
    private float scaleTol = 0.2f;
    
    private float rotationInRadiansTol = (float)(20. * Math.PI/180.);
    
    private int transXYTol = 30;
    
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
        float rotationInRadians) {
        img1 = image1;
        img2 = image2;
        doDetermineScale = false;
        scaleSetByUser = scale;
        rotationInRadianSetByUser = rotationInRadians;
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
    
    public FeatureMatcherWrapper(ImageExt image1, ImageExt image2, float scale,
        float rotationInRadians, String debugTagPrefix) {
        img1 = image1;
        img2 = image2;
        doDetermineScale = false;
        scaleSetByUser = scale;
        rotationInRadianSetByUser = rotationInRadians;
        debug = true;
        this.debugTagPrefix = debugTagPrefix;
    }
    
    private TransformationParameters solveForScale() throws IOException, 
        NoSuchAlgorithmException {
        
        BlobScaleFinderWrapper scaleFinder = null;
            
        if (debug) {
            scaleFinder = new BlobScaleFinderWrapper(img1, img2, debugTagPrefix);
        } else {
            scaleFinder = new BlobScaleFinderWrapper(img1, img2);
        }
        
        TransformationParameters params = scaleFinder.calculateScale();
        
        if (scaleFinder.img1Helper.didApplyHistEq()) {
            this.didApplyHistEq = didApplyHistEq;
        }
        this.gsImg1 = scaleFinder.img1Helper.getGreyscaleImage();
        this.gsImg2 = scaleFinder.img2Helper.getGreyscaleImage();
        
        return params;
    }
    
    public CorrespondenceList matchFeatures() throws IOException, NoSuchAlgorithmException {
        
        if (doDetermineScale) {
            params = solveForScale();
            if (params == null) {
                //TODO: consider whether to make an assumption that scale=1
                return null;
            }
        } else {
            applyHistEqIfNeeded();
        }
        
        extractCornerRegions();
        
        CorrespondenceList cl = null;
        
        if (params != null) {
            cl = findCorrespondence(params);
        } else if (this.rotationInRadianSetByUser > Float.MIN_VALUE) {
            cl = findCorrespondence(this.scaleSetByUser, this.rotationInRadianSetByUser);
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

    private void applyHistEqIfNeeded() {
        
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
            didApplyHistEq = true;
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

    private CorrespondenceList findCorrespondence(TransformationParameters params) {
        
        FeatureMatcher matcher = new FeatureMatcher();
        
        CorrespondenceList cl = matcher.findSimilarFeatures(gsImg1, gXY1, theta1,
            cornerRegions1.toArray(new CornerRegion[cornerRegions1.size()]),
            gsImg2, gXY2, theta2,
            cornerRegions2.toArray(new CornerRegion[cornerRegions2.size()]), 
            params, scaleTol, rotationInRadiansTol, transXYTol);

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
    
    private CorrespondenceList findCorrespondence(float scale, 
        float rotationInRadians) {
        
        FeatureMatcher matcher = new FeatureMatcher();
        
        CorrespondenceList cl = matcher.findSimilarFeatures(gsImg1, gXY1, theta1,
            cornerRegions1.toArray(new CornerRegion[cornerRegions1.size()]),
            gsImg2, gXY2, theta2,
            cornerRegions2.toArray(new CornerRegion[cornerRegions2.size()]), 
            scale, rotationInRadians, scaleTol, rotationInRadiansTol);

        return cl;
    }
    
}
