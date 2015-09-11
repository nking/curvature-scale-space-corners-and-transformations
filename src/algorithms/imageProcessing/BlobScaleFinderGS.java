package algorithms.imageProcessing;

import algorithms.misc.Histogram;
import algorithms.misc.Misc;
import algorithms.misc.MiscDebug;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import java.io.IOException;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * determine scale between 2 images using blob contours and methods for
 * segmentation on greyscale images (more specifically, using KMPP for
 * finding the k bins of intensity).
 *
 * NOT READY FOR USE YET.
 *
 * @author nichole
 */
public class BlobScaleFinderGS extends AbstractBlobScaleFinder {

    /**
     * NOT READY FOR USE YET.
     * From the given images, determine the scale between them and roughly
     * estimate the rotation and translation too.  Note that image processing
     * such as sky masks should be applied before using this method.
     * Also note that it is expected that it will be followed by a more rigorous
     * solver for Euclidean transforms such as the FeatureMatcher and then
     * the epipolar projection solver.
     *
     * applies KMeansPlusPlus algorithm to the values in each input image to 
     * create kBands segmented image.
     * 
     * Note, images should be pre-processed to have similar comparable qualities,
     * such as histogram equalization, or resolution changes.
     * 
     * @param img1 the first image holding objects for which a Euclidean
     * transformation is found that can be applied to the image to put it in
     * the same scale reference frame as image2.
     * @param img2 the second image representing the reference frame that
     * image1 is transformed to using the resulting parameters,
     * @param segType 
     * <pre>
     *     0 applies KMeansPlusPlus to obtain k segmented images;
     *     1 applies KMeansPlusPlus to obtain a binary segmented image and 
     *       then uses adaptive mean thresholding.
     * </pre>
     * @return Euclidean scale to be applied to image1 to place it in the same
     * scale reference frame as image2.  Rotation and transformation are also
     * roughly solved for.
     * @throws java.io.IOException
     * @throws java.security.NoSuchAlgorithmException
     */
    public TransformationParameters calculateScale(final GreyscaleImage img1,
        final GreyscaleImage img2, final int k,
        final int smallestGroupLimit, final int largestGroupLimit,
        final float[] outputScaleRotTransXYStDev, final int segType)
        throws IOException, NoSuchAlgorithmException {
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        GreyscaleImage imgS1 = null; 
        GreyscaleImage imgS2 = null;
        
        switch(segType) {
            case 1: {
                int scl1 = (int)Math.ceil(Math.max((float)img1.getWidth()/300.f,
                    (float)img1.getHeight()/300.f));
                scl1 *= 10;
                int scl2 = (int)Math.ceil(Math.max((float)img2.getWidth()/300.f,
                    (float)img2.getHeight()/300.f));
                scl2 *= 10;
                imgS1 = img1.copyImage();
                imgS2 = img2.copyImage();
                imageProcessor.applyImageSegmentation(imgS1, 2);
                imageProcessor.applyImageSegmentation(imgS2, 2);
                imageProcessor.applyAdaptiveMeanThresholding(imgS1, scl1);
                imageProcessor.applyAdaptiveMeanThresholding(imgS2, scl2);
                break;
            }
            default: {
                // case 0:
                imgS1 = img1.copyImage();
                imageProcessor.applyImageSegmentation(imgS1, k);
                imgS2 = img2.copyImage();
                imageProcessor.applyImageSegmentation(imgS2, k);
            }
        }
        
        TransformationParameters params = calculateScaleImpl(img1, img2, 
            imgS1, imgS2, k, smallestGroupLimit, largestGroupLimit, 
            outputScaleRotTransXYStDev);

        return params;
    }
    
    /**
       <pre>
        (3) extract the top 10 blobs and their contours from k=2 segmented image:
            (3a) perform segmentation for k=2
            (3b) find blobs w/ DFSContiguousValueFinder for each intensity level
            (3c) use EdgeExtractorForBlobBorder to extract closed contour for
                 each blob
            (3d) return the top 10 longest contours and the blobs
       </pre>
     * @param k
     * @param img
     * @param segImg segmented image img
     * @param outputBlobs
     * @param outputBounds
     * @param smallestGroupLimit
     * @param largestGroupLimit
     */
    protected void extractBlobsFromSegmentedImage(final int k,
        GreyscaleImage img, GreyscaleImage segImg,
        List<Set<PairInt>> outputBlobs, List<PairIntArray> outputBounds,
        int smallestGroupLimit, int largestGroupLimit) {

        extractBlobsFromSegmentedImage(k, segImg, outputBlobs, smallestGroupLimit,
            largestGroupLimit);

        if (outputBlobs.isEmpty()) {
            return;
        }

        boolean discardWhenCavityIsSmallerThanBorder = true;

        extractBoundsOfBlobs(img, outputBlobs, outputBounds, img.getWidth(),
            img.getHeight(), discardWhenCavityIsSmallerThanBorder);

    }

    /**
     * @param k
     * @param segImg segmented image
     * @param outputBlobs
     * @param smallestGroupLimit
     * @param largestGroupLimit
     */
    protected void extractBlobsFromSegmentedImage(int k,
        GreyscaleImage segImg,
        List<Set<PairInt>> outputBlobs, int smallestGroupLimit,
        int largestGroupLimit) {

        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();

        Map<Integer, Integer> freqMap = Histogram.createAFrequencyMap(segImg);

        for (Map.Entry<Integer, Integer> entry : freqMap.entrySet()) {

            Integer pixValue = entry.getKey();

            DFSContiguousValueFinder finder = new DFSContiguousValueFinder(segImg);
            finder.setMinimumNumberInCluster(smallestGroupLimit);
            finder.findGroups(pixValue.intValue());

            int nGroups = finder.getNumberOfGroups();

            for (int i = 0; i < nGroups; ++i) {

                PairIntArray xy = finder.getXY(i);

                if (xy.getN() < largestGroupLimit) {

                    Set<PairInt> points = Misc.convert(xy);

                    // skip blobs that are on the image boundaries because they
                    // are incomplete
                    if (!curveHelper.hasNumberOfPixelsOnImageBoundaries(3,
                        points, segImg.getWidth(), segImg.getHeight())) {

                        outputBlobs.add(points);
                    }
                }
            }
        }

if (debug) {
Image img0 = ImageIOHelper.convertImage(segImg);
int c = 0;
for (int i = 0; i < outputBlobs.size(); ++i) {
    Set<PairInt> blobSet = outputBlobs.get(i);
    int clr = ImageIOHelper.getNextColorRGB(c);
    for (PairInt p : blobSet) {
        int x = p.getX();
        int y = p.getY();
        ImageIOHelper.addPointToImage(x, y, img0, 0, clr);
    }
    c++;
}
MiscDebug.writeImageCopy(img0, "blobs_" + MiscDebug.getCurrentTimeFormatted() + ".png");
}

    }

    /**
     * 
     * @param img1
     * @param img2
     * @param imgS1 segmented image 1
     * @param imgS2 segmented image 2
     * @param k
     * @param smallestGroupLimit
     * @param largestGroupLimit
     * @param outputScaleRotTransXYStDev
     * @return
     */
    protected TransformationParameters calculateScaleImpl(GreyscaleImage img1,
        GreyscaleImage img2, GreyscaleImage imgS1, GreyscaleImage imgS2,
        int k, int smallestGroupLimit,
        int largestGroupLimit, float[] outputScaleRotTransXYStDev) {

        /*
        extract the top 10 or so blobs and their contours from k=2 segmented image:
        -- perform segmentation for k
        -- find blobs w/ DFSContiguousValueFinder for each intensity level
        -- use EdgeExtractorForBlobBorder to extract closed contour for each
           blob
        -- return the top 10 longest contours and the blobs
        */
        List<Set<PairInt>> blobs1 = new ArrayList<Set<PairInt>>();
        List<Set<PairInt>> blobs2 = new ArrayList<Set<PairInt>>();
        List<PairIntArray> bounds1 = new ArrayList<PairIntArray>();
        List<PairIntArray> bounds2 = new ArrayList<PairIntArray>();

        log.info("image1:");
        extractBlobsFromSegmentedImage(k, img1, imgS1, blobs1, bounds1,
            smallestGroupLimit, largestGroupLimit);
        log.info("image2:");
        extractBlobsFromSegmentedImage(k, img2, imgS2, blobs2, bounds2,
            smallestGroupLimit, largestGroupLimit);

        log.info("nBounds1=" + bounds1.size());

        log.info("nBounds2=" + bounds2.size());

        /*
        filter out dissimilar pairings:
        -- given blobs from image1 and blobs from image 2, use feature
           matching to rule out possible pairings, resulting in possible
           matches for each.
        */
        //ImageProcessor imageProcessor = new ImageProcessor();

//TODO: put debug sections in AOP for special build after replace aspectj
if (debug) {
Image img0 = ImageIOHelper.convertImage(img1);
for (int i = 0; i < bounds1.size(); ++i) {
    PairIntArray pa = bounds1.get(i);
    for (int j = 0; j < pa.getN(); ++j) {
        int x = pa.getX(j);
        int y = pa.getY(j);
        if (i == 0) {
            if (j == 0 || (j == (pa.getN() - 1))) {
                ImageIOHelper.addPointToImage(x, y, img0, 0, 200, 100, 0);
            } else {
                ImageIOHelper.addPointToImage(x, y, img0, 0, 255, 0, 0);
            }
        } else if (i == 1) {
            ImageIOHelper.addPointToImage(x, y, img0, 0, 0, 255, 0);
        } else {
            ImageIOHelper.addPointToImage(x, y, img0, 0, 0, 0, 255);
        }
    }
}
MiscDebug.writeImageCopy(img0, "blob_contours_1_" + MiscDebug.getCurrentTimeFormatted() + ".png");
img0 = ImageIOHelper.convertImage(img2);
for (int i = 0; i < bounds2.size(); ++i) {
    PairIntArray pa = bounds2.get(i);
    for (int j = 0; j < pa.getN(); ++j) {
        int x = pa.getX(j);
        int y = pa.getY(j);
        if (i == 0) {
            if (j == 0 || (j == (pa.getN() - 1))) {
                ImageIOHelper.addPointToImage(x, y, img0, 0, 200, 100, 0);
            } else {
                ImageIOHelper.addPointToImage(x, y, img0, 0, 255, 0, 0);
            }
        } else if (i == 1) {
            ImageIOHelper.addPointToImage(x, y, img0, 0, 0, 255, 0);
        } else {
            ImageIOHelper.addPointToImage(x, y, img0, 0, 0, 0, 255);
        }
    }
}
MiscDebug.writeImageCopy(img0, "blob_contours_2_" + MiscDebug.getCurrentTimeFormatted() + ".png");
int z = 1;
}

        /*
        solve for scale:
        -- use ContourMather to get scale solutions for each pairing then
           statistical basis of combining the results and removing
           outliers.
        */

        TransformationParameters params = solveForScale(img1, img2, blobs1,
            blobs2, bounds1, bounds2, outputScaleRotTransXYStDev);

        return params;
    }

}
