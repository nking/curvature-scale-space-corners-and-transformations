package algorithms.imageProcessing;

import algorithms.imageProcessing.util.MiscStats;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 class whose goal is to find best single euclidean transformation for 
 image1 to image2.  It uses segmentation to create blobs.
 It solves for transformation of each blob in one image against the other.
 It keeps the best transformation solution of each blob in image1.
 Then evaluates all blob1 best transformations to find the best overall
 for the image corners.
 -- pros: for panorama and stereo images it does lead to a solution and it
          solves for rotation too.
 -- cons: has long runtime due to curve to curve transformation comparisons.
          also, it sometimes discards very good matching points
          because the curve needs at least 3 in order to
          estimate a euclidean transformation so any curves w/ fewer than 
          2 points are skipped.
 
 Prefer to use EuclideanSegmentFeatureMatcher2 which is faster, excepting
 for images like the checkerboard tests.  The checkerboard tests do not
 have unique features so a solution is better found by 
 EuclideanSegmentFeatureMatcher which tries many combinations from detailed
 to larger level groupings.
 
 * @author nichole
 */
public class EuclideanSegmentFeatureMatcher extends AbstractFeatureMatcher { 

    protected TransformationParameters solutionTransformation = null;
    
    public EuclideanSegmentFeatureMatcher(ImageExt image1, ImageExt image2,
        FeatureMatcherSettings settings) {
        
        super(image1, image2, settings);
    }
    
    @Override
    protected boolean match(SegmentationType type, boolean useBinned) {
        
        int binFactor1, binFactor2;
        GreyscaleImage img1, img2;
        IntensityFeatures f1, f2;
        
        if (useBinned) {
            binFactor1 = img1Helper.getBinFactor(useBinned);
            binFactor2 = img2Helper.getBinFactor(useBinned);
            img1 = img1Helper.getGreyscaleImageBinned();
            img2 = img2Helper.getGreyscaleImageBinned();
            f1 = featuresBinned1;
            f2 = featuresBinned2;
        } else {
            binFactor1 = 1;
            binFactor2 = 1;
            img1 = img1Helper.getGreyscaleImage();
            img2 = img2Helper.getGreyscaleImage();
            f1 = features1;
            f2 = features2;
            if (!f1.gradientWasCreated()) {
                f1.calculateGradientWithGreyscale(
                    img1Helper.getGreyscaleImage().copyImage());
            }
            if (!f2.gradientWasCreated()) {
                f2.calculateGradientWithGreyscale(
                    img2Helper.getGreyscaleImage().copyImage());
            }
        }
        
        int dither2;
        if (settings.doUse2ndDerivCorners()) {
            dither2 = 1;
        } else {
            dither2 = DitherDefault.dither;
        }

        //long t0 = System.currentTimeMillis();
        
        BlobCornersEuclideanCalculator bsFinder = new BlobCornersEuclideanCalculator();

        // the solution for the binFactor modified images:
        MatchingSolution soln = bsFinder.solveTransformation(img1Helper, f1,
            type, useBinned, img2Helper, f2, type, useBinned, dither2);

        //int n1 = img1Helper.sumPointsOfInterest(type, useBinned);
        //int n2 = img2Helper.sumPointsOfInterest(type, useBinned);
        // long t1 = System.currentTimeMillis();
        //long t1Sec = (t1 - t0)/1000;
        //Logger.getLogger(this.getClass().getName()).info("matching(sec)=" + t1Sec);
              
        if (soln == null || soln.getComparisonStats() == null || 
            soln.getComparisonStats().isEmpty()) {
            return false;
        }
        
        // transform images to full size
        soln = transformSolutionToFullFrames(soln, img1Helper, img2Helper, 
            binFactor1, binFactor2);
            
        TransformationParameters params = soln.getParams();
                
        assert(params.getStandardDeviations() != null);

        log.info("params for type"
            + " (" + type.name() + ", binned=" + useBinned + ")"
            + " (" + type.name() + ", binned=" + useBinned + ")"
            + " : " + params.toString());

        log.info(String.format(
            "stDev scale=%.1f  stDev rot=%.0f  stDev tX=%.0f  stDev tY=%.0f",
            params.getStandardDeviations()[0], 
            params.getStandardDeviations()[1],
            params.getStandardDeviations()[2], 
            params.getStandardDeviations()[3]));
       
        List<FeatureComparisonStat> solnStats = soln.getComparisonStats();
        
        if (type.equals(SegmentationType.GREYSCALE_CANNY)) {
            //extractMoreCorners uses canny, so exit here
            copyToInstanceVars(solnStats);
            
            solutionTransformation = soln.getParams().copy();
                
            return true;
        }
        
        // look at intersection of solution to see if need more corners

        int tolXY;
        if (params.getStandardDeviations() != null) {
            tolXY = Math.round(Math.max(params.getStandardDeviations()[2],
                params.getStandardDeviations()[3]));
            if (tolXY < 3) {
                tolXY = 3;
            }
        } else {
            tolXY = 10;
        }

        int nLimit = 16;
        
        boolean covers = MiscStats.statsCoverIntersection(solnStats, params,
            img1.getWidth(), img1.getHeight(), img2.getWidth(), img2.getHeight()
        );
        
        boolean extractMoreCorners = (solnStats.size() < nLimit) || !covers;
        
        if (!extractMoreCorners) {

            copyToInstanceVars(solnStats);
            
            solutionTransformation = soln.getParams().copy();
                
            return true;
        }
            
        Set<CornerRegion> cr1 = new HashSet<CornerRegion>(); 
        Set<CornerRegion> cr2 = new HashSet<CornerRegion>(); 
        extractCannyCornerRegions(img1, img2, cr1, cr2);

        float scaleTol = 0.2f;
        float rotationInRadiansTol = (float)(20. * Math.PI/180.);
        List<FeatureComparisonStat> extraStats = findCorrespondence(img1, img2, 
            cr1, cr2, params, f1.getRotatedOffsets(), dither2, tolXY, scaleTol, 
            rotationInRadiansTol);
        
        if (extraStats == null) {
            
            log.severe("error finding more points to add to solution");
            
            copyToInstanceVars(solnStats);
            
            solutionTransformation = soln.getParams().copy();
                
            return true;
        }
        
        solnStats.addAll(extraStats);
        
        MiscStats.filterForDegeneracy(solnStats);
                    
        copyToInstanceVars(solnStats);
            
        solutionTransformation = soln.getParams().copy();
        
        return true;
    }
    
    public TransformationParameters getSolutionTransformation() {
        return solutionTransformation;
    }
    
}
