package algorithms.imageProcessing;

import java.io.IOException;
import java.security.NoSuchAlgorithmException;
import java.util.logging.Logger;

/**
 * determine scale between 2 images using blob contours.
 * NOT READY FOR USE YET.
 *
 * @author nichole
 */
public class BlobScaleFinderWrapper {

    protected Logger log = Logger.getLogger(this.getClass().getName());

    protected boolean debug = false;

    protected final int smallestGroupLimit = 100;
    protected final int largestGroupLimit = 5000;
        
    public void setToDebug() {
        debug = true;
    }

    /**
     * NOT READY FOR USE YET.
     * From the given images, determine the scale between them and roughly
     * estimate the rotation and translation too.  Note that image processing
     * such as sky masks should be applied before using this method.
     * Also note that it is expected that it will be followed by a more rigorous
     * solver for Euclidean transforms such as the FeatureMatcher and then
     * the epipolar projection solver.
     *
     * @param img1 the first image holding objects for which a Euclidean
     * transformation is found that can be applied to the image to put it in
     * the same scale reference frame as image2.
     * @param img2 the second image representing the reference frame that
     * image1 is transformed to using the resulting parameters,
     * @return Euclidean scale to be applied to image1 to place it in the same
     * scale reference frame as image2.  Rotation and transformation are also
     * roughly solved for.
     * @throws java.io.IOException
     * @throws java.security.NoSuchAlgorithmException
     */
    public TransformationParameters calculateScale(ImageExt img1,
        ImageExt img2) throws IOException, NoSuchAlgorithmException {

        /*
        tries to solve using a black and white binned segmentation.
        if that fails, tries to solve using a full size black and white segmentation.
        if that fails, tries to solve using a more rigorous color segmentation
        on binned images. 
            tries k=2 and k=3
        if that fails, tries to solve using a more rigorous color segmentation
        on full size images,
            starting with the smallest k and increasing until a k limit or
            a solution.
        */
        
        TransformationParameters params;
        
        // try the B&W solution first
        params = calculateScaleForBW(img1, img2);
        
        if (params != null) {
            return params;
        }
        
        
        // try the Clr solution
        int k = 2;
        params =  calculateScaleForClr(img1, img2, k);
        if (params != null) {
            return params;
        }
        k = 4;
        params =  calculateScaleForClr(img1, img2, k);
        if (params != null) {
            return params;
        }
        
        return params;
    }
    
    public TransformationParameters calculateScaleForBW(ImageExt img1,
        ImageExt img2) throws IOException, NoSuchAlgorithmException {

        /*
        tries to solve using a black and white binned segmentation.
        if that fails, tries to solve using a full size black and white segmentation.
        */
        
        GreyscaleImage img1Grey = img1.copyToGreyscale();
        GreyscaleImage img2Grey = img2.copyToGreyscale();

        //TODO: this decision might belong to a preprocessor:
        boolean didApplyHistEq = applyHistogramEqualizationIfNeeded(img1Grey,
            img2Grey);
        
        final int k = 2;
        
        //---------------------------------------------
        
        float minDimension = 300.f;//200.f
        int binFactor = (int) Math.ceil(
            Math.max((float)img1.getWidth()/minDimension,
            (float)img2.getHeight()/minDimension));
        int smallestGroupLimitBinned = smallestGroupLimit/(binFactor*binFactor);
        int largestGroupLimitBinned = largestGroupLimit/(binFactor*binFactor);

        log.info("binFactor=" + binFactor);

        // prevent from being smaller than needed for a convex hull
        if (smallestGroupLimitBinned < 4) {
            smallestGroupLimitBinned = 4;
        }

        ImageProcessor imageProcessor = new ImageProcessor();

        GreyscaleImage img1GreyBinned = imageProcessor.binImage(img1Grey,
            binFactor);
        GreyscaleImage img2GreyBinned = imageProcessor.binImage(img2Grey,
            binFactor);
        
        BlobScaleFinderBW scaleFinder = new BlobScaleFinderBW();
        
        if (debug) {
            scaleFinder.setToDebug();
        }
        
        float[] outputScaleRotTransXYStDevBinned = new float[4];
        
        TransformationParameters paramsBinned = scaleFinder.calculateScale(
            img1GreyBinned, img2GreyBinned, k, 
            smallestGroupLimitBinned, largestGroupLimitBinned,
            outputScaleRotTransXYStDevBinned);
        
        if (paramsBinned != null) {
            
            // put back into unbinned image reference frame
            float transX = paramsBinned.getTranslationX();
            float transY = paramsBinned.getTranslationY();

            transX *= binFactor;
            transY *= binFactor;

            outputScaleRotTransXYStDevBinned[2] *= binFactor;
            outputScaleRotTransXYStDevBinned[3] *= binFactor;

            paramsBinned.setTranslationX(transX);
            paramsBinned.setTranslationY(transY);
        
            log.info("binned params: " + paramsBinned.toString());

            log.info(String.format(
                "stDev scale=%.1f  stDev rot=%.0f  stDev tX=%.0f  stDev tY=%.0f",
                outputScaleRotTransXYStDevBinned[0], outputScaleRotTransXYStDevBinned[1],
                outputScaleRotTransXYStDevBinned[2], outputScaleRotTransXYStDevBinned[3]));

            //TODO: review this limit
            // sometimes a single pair solution is correct even though it has
            //   stdev of NaN due to no other pairs
            if (
                (outputScaleRotTransXYStDevBinned[0]/paramsBinned.getScale()) < 0.2) {
                return paramsBinned;
            }
        }

        //---------------------------------------------
        
        float[] outputScaleRotTransXYStDev = new float[4];

        scaleFinder = new BlobScaleFinderBW();
        
        if (debug) {
            scaleFinder.setToDebug();
        }
                
        TransformationParameters params = scaleFinder.calculateScale(
            img1Grey, img2Grey, k, smallestGroupLimit, largestGroupLimit,
            outputScaleRotTransXYStDev);
        
        if (params != null) {
            
            log.info("params: " + params.toString());

            log.info(String.format(
                "stDev scale=%.1f  stDev rot=%.0f  stDev tX=%.0f  stDev tY=%.0f",
                outputScaleRotTransXYStDev[0], outputScaleRotTransXYStDev[1],
                outputScaleRotTransXYStDev[2], outputScaleRotTransXYStDev[3]));

            if (paramsBinned != null) {

                float stat0 = (outputScaleRotTransXYStDevBinned[0]/paramsBinned.getScale());
                float stat1 = (outputScaleRotTransXYStDev[0]/params.getScale());

                log.info("comparing to paramsBinned.  stat0=" + stat0 
                    + " stat1=" + stat1);
                
                if (stat1 < stat0) {
                    return params;
                } else {
                    return paramsBinned;
                }
            }

            // TODO: consider keeping even if it is only one pair because
            // the contour matching followed by features already filters for
            // relative location and descriptor... but there may be textures...
            if (!Float.isNaN(outputScaleRotTransXYStDev[0])) {
                return params;
            }
        }

        // this is either null or possibly a single pair solution
        //return paramsBinned;
        
        return null;
    }

    protected boolean applyHistogramEqualizationIfNeeded(GreyscaleImage image1,
        GreyscaleImage image2) {

        // doing this automatically for now, but can use the stats insead to
        // decide
        boolean performHistEq = true;

        boolean useStatsToDecide = false;

        if (useStatsToDecide) {
            ImageStatistics stats1 = ImageStatisticsHelper.examineImage(image1,
                true);
            ImageStatistics stats2 = ImageStatisticsHelper.examineImage(image2,
                true);
            double median1DivMedian2 = stats1.getMedian() / stats2.getMedian();
            double meanDivMedian1 = stats1.getMean() / stats1.getMedian();
            double meanDivMedian2 = stats2.getMean() / stats2.getMedian();
            if (((median1DivMedian2 > 1) && ((median1DivMedian2 - 1) > 0.2))
                || ((median1DivMedian2 < 1) && (median1DivMedian2 < 0.8))) {
                performHistEq = true;
            } else if (((meanDivMedian1 > 1) && ((meanDivMedian1 - 1) > 0.2))
                || ((meanDivMedian1 < 1) && (meanDivMedian1 < 0.8))) {
                performHistEq = true;
            } else if (((meanDivMedian2 > 1) && ((meanDivMedian2 - 1) > 0.2))
                || ((meanDivMedian2 < 1) && (meanDivMedian2 < 0.8))) {
                performHistEq = true;
            }
        }

        if (performHistEq) {
            HistogramEqualization hEq = new HistogramEqualization(image1);
            hEq.applyFilter();
            hEq = new HistogramEqualization(image2);
            hEq.applyFilter();
        }

        return performHistEq;
    }

    public TransformationParameters calculateScaleForClr(ImageExt img1,
        ImageExt img2, final int k) throws IOException, NoSuchAlgorithmException {

        /*
        tries to solve using a color binned segmentation.
        if that fails, tries to solve using a full size colorsegmentation.
        */
        
        //GreyscaleImage img1Grey = img1.copyToGreyscale();
        //GreyscaleImage img2Grey = img2.copyToGreyscale();

        /*
        consider if histogram equalization is needed
        */
                
        //---------------------------------------------
        
        float minDimension = 300.f;//200.f
        int binFactor = (int) Math.ceil(
            Math.max((float)img1.getWidth()/minDimension,
            (float)img2.getHeight()/minDimension));
        int smallestGroupLimitBinned = smallestGroupLimit/(binFactor*binFactor);
        int largestGroupLimitBinned = largestGroupLimit/(binFactor*binFactor);

        log.info("binFactor=" + binFactor);

        // prevent from being smaller than needed for a convex hull
        if (smallestGroupLimitBinned < 4) {
            smallestGroupLimitBinned = 4;
        }

        ImageProcessor imageProcessor = new ImageProcessor();

        ImageExt img1Binned = imageProcessor.binImage(img1, binFactor);
        ImageExt img2Binned = imageProcessor.binImage(img2, binFactor);
        
        BlobScaleFinderClr scaleFinder = new BlobScaleFinderClr();
        
        if (debug) {
            scaleFinder.setToDebug();
        }
        
        float[] outputScaleRotTransXYStDevBinned = new float[4];
        
        TransformationParameters paramsBinned = scaleFinder.calculateScale(
            img1Binned, img2Binned, k, 
            smallestGroupLimitBinned, largestGroupLimitBinned,
            outputScaleRotTransXYStDevBinned);
        
        if (paramsBinned != null) {
            
            // put back into unbinned image reference frame
            float transX = paramsBinned.getTranslationX();
            float transY = paramsBinned.getTranslationY();

            transX *= binFactor;
            transY *= binFactor;

            outputScaleRotTransXYStDevBinned[2] *= binFactor;
            outputScaleRotTransXYStDevBinned[3] *= binFactor;

            paramsBinned.setTranslationX(transX);
            paramsBinned.setTranslationY(transY);
        
            log.info("binned params: " + paramsBinned.toString());

            log.info(String.format(
                "stDev scale=%.1f  stDev rot=%.0f  stDev tX=%.0f  stDev tY=%.0f",
                outputScaleRotTransXYStDevBinned[0], outputScaleRotTransXYStDevBinned[1],
                outputScaleRotTransXYStDevBinned[2], outputScaleRotTransXYStDevBinned[3]));

            //TODO: review this limit
            // sometimes a single pair solution is correct even though it has
            //   stdev of NaN due to no other pairs
            if (
                (outputScaleRotTransXYStDevBinned[0]/paramsBinned.getScale()) < 0.2) {
                return paramsBinned;
            }
        }

        //---------------------------------------------
        
        float[] outputScaleRotTransXYStDev = new float[4];

        scaleFinder = new BlobScaleFinderClr();
        
        if (debug) {
            scaleFinder.setToDebug();
        }
                
        TransformationParameters params = scaleFinder.calculateScale(
            img1, img2, k, smallestGroupLimit, largestGroupLimit,
            outputScaleRotTransXYStDev);
        
        if (params != null) {
            
            log.info("params: " + params.toString());

            log.info(String.format(
                "stDev scale=%.1f  stDev rot=%.0f  stDev tX=%.0f  stDev tY=%.0f",
                outputScaleRotTransXYStDev[0], outputScaleRotTransXYStDev[1],
                outputScaleRotTransXYStDev[2], outputScaleRotTransXYStDev[3]));

            if (paramsBinned != null) {

                float stat0 = (outputScaleRotTransXYStDevBinned[0]/paramsBinned.getScale());
                float stat1 = (outputScaleRotTransXYStDev[0]/params.getScale());

                log.info("comparing to paramsBinned.  stat0=" + stat0 
                    + " stat1=" + stat1);
                
                if (stat1 < stat0) {
                    return params;
                } else {
                    return paramsBinned;
                }
            }

            return params;
        }

        // this is either null or possibly a single pair solution
        //return paramsBinned;
       
        return null;
    }

}
