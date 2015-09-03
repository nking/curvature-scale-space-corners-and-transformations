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
    
    protected final ImageExt img1;
    protected final ImageExt img2;
    
    protected final GreyscaleImage img1Grey;
    protected final GreyscaleImage img2Grey;
    
    /**
     * 
     * @param img1 the first image holding objects for which a Euclidean
     * transformation is found that can be applied to put it in
     * the same scale reference frame as image2.
     * @param img2 the second image representing the reference frame that
     * image1 is transformed to using the resulting parameters,
     */
    public BlobScaleFinderWrapper(ImageExt img1, ImageExt img2) {
        this.img1 = img1;
        this.img2 = img2;
        
        img1Grey = img1.copyToGreyscale();
        img2Grey = img2.copyToGreyscale();

        //TODO: this decision might belong to a preprocessor:
        boolean didApplyHistEq = applyHistogramEqualizationIfNeeded(img1Grey,
            img2Grey);
    }
    
    public void setToDebug() {
        debug = true;
    }

    /**
     * NOT READY FOR USE YET.
     * From the given images, determine the scale between them and roughly
     * estimate the rotation and translation too.  Note that image processing
     * such as sky masks should be applied before using this method.
     * Also note that it is expected that it will be followed by a more rigorous
     * solver such as the FeatureMatcher for a correspondence list 
     * (and a better Euclidean transform) to be used in
    
     * @return Euclidean scale to be applied to image1 to place it in the same
     * scale reference frame as image2.  Rotation and transformation are also
     * roughly solved for.
     * @throws java.io.IOException
     * @throws java.security.NoSuchAlgorithmException
     */
    public TransformationParameters calculateScale() throws IOException, 
        NoSuchAlgorithmException {

        /*
        Tries to solve using various combinations of binning and segmentation,
        starting with binned images because the solution if possible with them
        is faster.
        (1) Greyscale binned, k=2
        (2) Clr binned, k=2
            -- possibly watershed on blobs if no solution
        (3) Clr binned, k=3 or k=4
            -- possibly watershed on blobs if no solution
        (4) Greyscale not binned, k=2
        (5) Clr not binned, k=2
        */
        
        TransformationParameters params;
        
        // (0) try the greyscale solution first, binned=true ---------------------
        float[] scaleRotTransXYStDev0 = new float[4];
        params = calculateScaleForGSUseBinning(scaleRotTransXYStDev0);
   
        if (params != null) {
            return params;
        }
        
        // (1) try the Clr solution, binned=true ---------------------------
        int k = 2;
        float[] scaleRotTransXYStDev1 = new float[4];
        params = calculateScaleForClrUseBinning(img1, img2, k, 
            scaleRotTransXYStDev1);
        if (params != null) {
            return params;
        }
        
        // (2) try the greyscale solution, binned=false --------------------------
        float[] scaleRotTransXYStDev2 = new float[4];
        params = calculateScaleForGS(img1, img2, scaleRotTransXYStDev2);
        if (params != null) {
            return params;
        }
        
        // (3) try the Clr solution, binned=false ---------------------------
        k = 2;
        float[] scaleRotTransXYStDev3 = new float[4];
        params = calculateScaleForClr(img1, img2, k,  scaleRotTransXYStDev3);
        if (params != null) {
            return params;
        }
        
        k = 4;
        float[] scaleRotTransXYStDev4 = new float[4];
        params = calculateScaleForClr(img1, img2, k,  scaleRotTransXYStDev4);
        if (params != null) {
            return params;
        }
 
        return null;
    }
    
    /**
     * calculate scale using the full size image and B&W segmentation algorithm.
     * @param img1
     * @param img2
     * @param outputScaleRotTransXYStDev
     * @return
     * @throws IOException
     * @throws NoSuchAlgorithmException 
     */
    public TransformationParameters calculateScaleForGS(ImageExt img1,
        ImageExt img2, final float[] outputScaleRotTransXYStDev) 
        throws IOException, NoSuchAlgorithmException {

        final int k = 2;
        
        BlobScaleFinderGS scaleFinder = new BlobScaleFinderGS();
        
        if (debug) {
            scaleFinder.setToDebug();
        }
                
        scaleFinder = new BlobScaleFinderGS();
        
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
            
            // TODO: consider keeping even if it is only one pair because
            // the contour matching followed by features already filters for
            // relative location and descriptor... but there may be textures...
            // see above comments about NaN
            if (
                ((outputScaleRotTransXYStDev[0]/params.getScale()) < 0.2)
                ) {
                return params;
            }
        }

        // this is either null or possibly a single pair solution
        //return paramsBinned;
        
        return null;
    }
    
    TransformationParameters calculateScaleForGSUseBinning(
        final float[] outputScaleRotTransXYStDevBinned) 
        throws IOException, NoSuchAlgorithmException {

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
        
        BlobScaleFinderGS scaleFinder = new BlobScaleFinderGS();
        
        if (debug) {
            scaleFinder.setToDebug();
        }
                
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
            //   stdev of NaN due to no other pairs.
            // might assume that only high quality matches have made it
            // here at this stage...
            if (
                ((outputScaleRotTransXYStDevBinned[0]/paramsBinned.getScale()) < 0.2)
                ) {
                return paramsBinned;
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

    public TransformationParameters calculateScaleForClrUseBinning(ImageExt img1,
        ImageExt img2, final int k, 
        final float[] outputScaleRotTransXYStDevBinned) throws IOException, 
        NoSuchAlgorithmException {
          
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
            //   stdev of NaN due to no other pairs.
            // might assume that only high quality matches have made it
            // here at this stage...
            if (
                ((outputScaleRotTransXYStDevBinned[0]/paramsBinned.getScale()) < 0.2)
                || Float.isNaN(outputScaleRotTransXYStDevBinned[0])) {
                return paramsBinned;
            }
        }

        // this is either null or possibly a single pair solution
        //return paramsBinned;
        return null;
    }
    
    public TransformationParameters calculateScaleForClr(ImageExt img1,
        ImageExt img2, final int k, final float[] outputScaleRotTransXYStDev) 
        throws IOException, NoSuchAlgorithmException {

        BlobScaleFinderClr scaleFinder = new BlobScaleFinderClr();
        
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

            // TODO: consider keeping even if it is only one pair because
            // the contour matching followed by features already filters for
            // relative location and descriptor... but there may be textures...
            // see above comments about NaN
            if (
                ((outputScaleRotTransXYStDev[0]/params.getScale()) < 0.2)
                || Float.isNaN(outputScaleRotTransXYStDev[0])) {
                return params;
            }
        }

        // this is either null or possibly a single pair solution
        //return paramsBinned;
        return null;
    }

}
