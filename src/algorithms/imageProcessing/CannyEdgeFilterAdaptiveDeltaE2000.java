package algorithms.imageProcessing;

import algorithms.compGeometry.HoughTransform;
import algorithms.misc.Misc;
import algorithms.misc.MiscDebug;
import java.util.logging.Logger;
import java.util.HashSet;
import algorithms.util.PairInt;
import java.util.Map;
import java.util.Set;

/**
 * The CannyEdge filter is an algorithm to operate on an image to
 * replace objects with their edges.
 *
 * The class began by following the general advice given in
 * "Performance Analysis of Adaptive Canny Edge Detector
 * Using Bilateral Filter" by Rashmi, Kumar, Jaiswal, and Saxena,
 * but made modifications afterwards.
 * Their paper has the following qualities:
 *<pre>
 * -- instead of a Gaussian filter, uses a bilateral filter
 * -- Otsu's method is used for adaptive threshold levels
 * -- they note that their algorithm performs better on high quality images but
 *    does not find edges well in low signal-to-noise images.
 * </pre>
 * This class uses 2 one dimensional binomial filters for smoothing.
 *
 *
 * Not ready for use yet...
 *
 * @author nichole
 */
public class CannyEdgeFilterAdaptiveDeltaE2000 {

    /** the factor that the low threshold is below the high threshold in the
    2 layer filter.
    */
    protected float factorBelowHighThreshold = 2.f;

    private EdgeFilterProducts filterProducts = null;

    private boolean performNonMaxSuppr = true;

    private boolean debug = false;

    private boolean restoreJunctions = true;

    private boolean useHigherThresholdIfNeeded = false;

    /**
     * the sigma from the blur combined with the sigma present in the gradient
     * are present in this variable by the end of processing.
     * The variable is used to interpret resolution of theta angles, for example.
     * The sigmas are combined using: sigma^2 = sigma_1^2 + sigma_2^2.
     * The FWHM of a gaussian is approx 2.35 * sigma.
     * (HWZI is about 5 * sigma, by the way).
     * So, for the default use of the filter, a sigma of 1 combined with sqrt(1)/2
     * results in a minimum resolution of 3 pixels, hence about 19 degrees.
     */
    private double approxProcessedSigma = 0;

    /**
     * default number of levels used in the histogram used to find the
     * high threshold in the 2 - layer thresholding.  it's the number
     *  of levels given to Otsu's multi-level thresholding.
     * for a value of 1, the histogram is performed over the entire image.
     * for a value of 2, histograms are calculated for 2 sets of calculations
     * over the image, one in which the average intensity of neighbors is
     * in bin1 and the other in bin2, etc... the histograms are indexed
     * by the average values of the neighbors.
     */
    private int numberOfLevelsForHistogram = 1;

    private boolean useAdaptive2Layer = true;

    private float otsuScaleFactor = 0.75f;//0.65f;

    protected Logger log = Logger.getLogger(this.getClass().getName());

    protected boolean useZhangSuen = true;

    public CannyEdgeFilterAdaptiveDeltaE2000() {
    }

    public void setToDebug() {
        debug = true;
    }

    public void setToNotUseZhangSuen() {
        useZhangSuen = false;
    }

    /**
     * to enable more complete contours, use this to restore pixels that were
     * removed during non-maximum suppression that disconnected edges and
     * have values above the low threshold of the 2 layer adaptive filter.
     */
    public void setToNotRestoreJunctions() {
        restoreJunctions = false;
    }

    /**
     * by default this is 0.45.
     * @param factor
     */
    public void setOtsuScaleFactor(float factor) {
        otsuScaleFactor = factor;
    }

    public void setToNotUseNonMaximumSuppression() {
        performNonMaxSuppr = false;
    }

    /**
     * use this setting to check the number of thinned points and if the number
     * is very high, resets the parameters to nLevels=1
     * and otsuFactor=2.5 and performs the line thinning again.   Note that
     * if the image has alot of small textures, using PhaseCongruencyDetector
     * in default mode and k=1 w nsCale = 6 or 3, etc may produce better results.
     */
    public void setToUseHigherThresholdIfNeeded() {
        useHigherThresholdIfNeeded = true;
    }

    /**
     * override the default number of levels used in the histogram used to find the
     * high threshold in the 2 layer thresholding, The default is number of
     * levels is 1.  It's the number of levels given to Otsu's multi-level
     * thresholding.  For a value of 1, the histogram is performed over the
     * entire image.  For a value of 2, histograms are calculated for 2 sets of
     * calculations over the image, one in which the average intensity of
     * neighbors is in bin1 (values 0 to 126)and the other in
     * bin2 (values 127 to 255), etc... The histograms are indexed
     * by the average values of the neighbors.
     * @param nLevels
     */
    public void overrideDefaultNumberOfLevels(int nLevels) {
        if (nLevels < 1 || nLevels > 255) {
            throw new IllegalArgumentException(
                "nLevels must be between 0 and 255, inclusive");
        }
        this.numberOfLevelsForHistogram = nLevels;
    }

    public void setToUseSingleThresholdIn2LayerFilter() {
        useAdaptive2Layer = false;
    }

    /**
     * override the default factor of low threshold below high threshold, which
     * is 2.
     * @param factor
     */
    public void override2LayerFactorBelowHighThreshold(float factor) {
        factorBelowHighThreshold = factor;
    }

    public void applyFilter(final ImageExt input) {

        if (input.getWidth() < 3 || input.getHeight() < 3) {
            throw new IllegalArgumentException("images should be >= 3x3 in size");
        }
        
        SIGMA sigma = SIGMA.ZEROPOINTFIVE;
        ImageProcessor imageProcessor = new ImageProcessor();
        imageProcessor.blur(input, sigma, 0, 255);
        approxProcessedSigma = SIGMA.getValue(sigma);

        //(2) create gradient
        // uses a binomial filter for a first derivative gradient, sobel.
        filterProducts = createGradient(input);

        GreyscaleImage gradientCopyBeforeThinning = filterProducts.getGradientXY().copyImage();

        approxProcessedSigma = Math.sqrt(
            approxProcessedSigma*approxProcessedSigma + (1./4.));

        Set<PairInt> removedDisconnecting = new HashSet<PairInt>();

        if (debug) {
            MiscDebug.writeImage(filterProducts.getGradientXY(), "_before_nms_");
        }

        //(3) non-maximum suppression
        if (performNonMaxSuppr) {
            applyNonMaximumSuppression(filterProducts, removedDisconnecting);
        }

        if (debug) {
            MiscDebug.writeImage(filterProducts.getGradientXY(), "_after_nms_");
        }

        //(4) adaptive 2 layer filter
        apply2LayerFilter(filterProducts.getGradientXY(), removedDisconnecting,
            gradientCopyBeforeThinning);

        if (useHigherThresholdIfNeeded) {
            int c = countAboveZero(filterProducts.getGradientXY());
            float cFraction = (float)c/(float)input.getNPixels();
            if (cFraction > 0.1) {
                this.useAdaptive2Layer = false;
                this.numberOfLevelsForHistogram = 1;
                this.otsuScaleFactor = 2.5f;

                filterProducts.getGradientX().resetTo(gradientCopyBeforeThinning);
                if (performNonMaxSuppr) {
                    applyNonMaximumSuppression(filterProducts, removedDisconnecting);
                }
                apply2LayerFilter(filterProducts.getGradientXY(), removedDisconnecting,
                    gradientCopyBeforeThinning);
            }
        }

        if (restoreJunctions) {
            int minResolution = (int)Math.ceil(2.35 * approxProcessedSigma);
            int minResolvableAngle = (int)Math.ceil(Math.atan2(1, minResolution) * 180./Math.PI);
            if (minResolvableAngle < 0) {
                minResolvableAngle *= -1;
            }
            MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
            curveHelper.additionalThinning45DegreeEdges2(
                filterProducts.getTheta(), filterProducts.getGradientXY(),
                minResolvableAngle);
        }

        // is this necessary?
        for (int i = 0; i < filterProducts.getGradientXY().getNPixels(); ++i) {
            int v = filterProducts.getGradientXY().getValue(i);
            if (v < 0) {
                filterProducts.getGradientXY().setValue(i, 0);
            }
        }

        applyHoughBasedLineThinning(filterProducts);
    }

    /*
    Tracing edges through the image and hysteresis thresholding

       It uses 2 "lower" threshold limits, an upper and lower edge limit.

        points below low threshhold are rejected.

        points above the highest threshhold are "sure-edge" pixels.

        for points between highest and low threshold:
            if connected to "sure-edge"
                pixels, are part of an edge, else discarded.

        Canny recommends:
            high limit
           -----------   should be 2 to 3 based on SNR
            low limit

       Note that pixels that are attached to the image boundaries should not
       be nulled in the 2nd stage of the layer if that would disconnect a line.

       The threshold is determined adaptively using Otsu's multilevel threshold
        calculation.
    */
    protected void apply2LayerFilter(final GreyscaleImage gradientXY,
        Set<PairInt> removedDisconnecting,
        GreyscaleImage gradientCopyBeforeThinning) {

        int n = gradientXY.getNPixels();
        int w = gradientXY.getWidth();
        int h = gradientXY.getHeight();

        if (w < 3 || h < 3) {
            throw new IllegalArgumentException("images should be >= 3x3 in size");
        }

        int nLevels = numberOfLevelsForHistogram;

        OtsuThresholding ot = new OtsuThresholding();

        int[] pixelThresholds = null;
        if (useAdaptive2Layer && (nLevels > 1)) {
            pixelThresholds = ot.calculateMultiLevelThreshold256(gradientXY, nLevels);
        }

        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;

        /*
        the paper suggests extending the association of a bright "sure edge"
        pixel to a radius of 2 neighbors (here neighbors have intensity > t_low).

        The traversal used in other methods is to visit each pixel and if
           its intensity > thigh, add the pixel and then visit it's neighbors
           to see if their intensities are > tLow.
           (there are redundant visits, but the results are idempotent).

        The extended association is better handled from the perspective of
        only adding the pixel currently visited using rules:
            if v < t_low
                continue
            else if v > t_high
                set value in output
            else
                if any neighbor has v > t_high
                    set value of current pixel in output
            else
                if any neighbor has value > t_low
                   search the 5x5 region and if any neighbor has v > t_high,
                      set value in current pixel in output

        Note:
           If using nLevels > 1, should consider the different results
        of using the threshold determined by the current pixel versus the
        thresholds extracted from the perspective of the "sure edge" pixels
        only, and then a pattern of searching around it.
        TODO: write a second method for the later to try that too...
        */

        float tHigh = 0;
        float tLow = 0;
        if (!useAdaptive2Layer || (nLevels == 1)) {
            tHigh = otsuScaleFactor * ot.calculateBinaryThreshold256(gradientXY);
            tLow = tHigh/factorBelowHighThreshold;
        }

        GreyscaleImage img2 = gradientXY.createWithDimensions();

        for (int i = 0; i < gradientXY.getNPixels(); ++i) {

            if (useAdaptive2Layer && (nLevels > 1)) {

                tHigh = otsuScaleFactor * pixelThresholds[i];

                tLow = tHigh/factorBelowHighThreshold;
            }

            int v = gradientXY.getValue(i);

            if (v < tLow) {
                continue;
            } else if (v > tHigh) {
                img2.setValue(i, v);
                continue;
            }

            int x = gradientXY.getCol(i);
            int y = gradientXY.getRow(i);

            boolean foundHigh = false;
            boolean foundMid = false;

            for (int k = 0; k < dxs.length; ++k) {
                int x2 = x + dxs[k];
                int y2 = y + dys[k];
                if ((x2 < 0) || (y2 < 0) || (x2 > (w - 1)) || (y2 > (h - 1))) {
                    continue;
                }
                int v2 = gradientXY.getValue(x2, y2);
                if (v2 > tHigh) {
                    foundHigh = true;
                    break;
                } else if (v2 > tLow) {
                    foundMid = true;
                }
            }
            if (foundHigh) {
                img2.setValue(i, v);
                continue;
            }
            if (!foundMid) {
                continue;
            }
            // search the 5 by 5 region for a "sure edge" pixel
            for (int dx = -2; dx <= 2; ++dx) {
                int x2 = x + dx;
                if ((x2 < 0) || (x2 > (w - 1))) {
                    continue;
                }
                for (int dy = -2; dy <= 2; ++dy) {
                    int y2 = y + dy;
                    if ((y2 < 0) || (y2 > (h - 1))) {
                        continue;
                    }
                    if (x2 == x && y2 == y) {
                        continue;
                    }
                    int v2 = gradientXY.getValue(x2, y2);
                    if (v2 > tHigh) {
                        img2.setValue(i, v);
                        foundHigh = true;
                        break;
                    }
                }
                if (foundHigh) {
                    break;
                }
            }
        }

        gradientXY.resetTo(img2);

        if (debug) {
            MiscDebug.writeImage(gradientXY, "_after_2layer_before_pltc_");
        }

        applyPostLineThinningCorrections(gradientXY,
            gradientCopyBeforeThinning);

        if (debug) {
            MiscDebug.writeImage(gradientXY, "_after_linethinning_1_");
        }

        if (restoreJunctions) {

            // while have information about the thresholds, will use them to
            // make decisions about restoring pixels that disconnected edges.
            PairInt[][] neighborCoordOffsets
                = AbstractLineThinner.createCoordinatePointsForEightNeighbors(
                    0, 0);

            for (PairInt p : removedDisconnecting) {

                if (ImageSegmentation.doesDisconnect(gradientXY,
                    neighborCoordOffsets, p.getX(), p.getY())) {

                    int x = p.getX();
                    int y = p.getY();
                    int i = gradientXY.getIndex(x, y);

                    if (useAdaptive2Layer && (nLevels > 1)) {
                        tHigh = otsuScaleFactor * pixelThresholds[i];
                        tLow = tHigh/factorBelowHighThreshold;
                    }

                    int v = gradientCopyBeforeThinning.getValue(p);

                    if (v > tLow) {
                        if (isAdjacentToAHorizOrVertLine(gradientXY, x, y, 3)) {
                            gradientXY.setValue(x, y, 255);
                        }
                    }
                }
            }

            if (debug) {
                MiscDebug.writeImage(gradientXY, "_after_restore_junctions_");
            }

            applyPostLineThinningCorrections(gradientXY,
                gradientCopyBeforeThinning);

            if (debug) {
                MiscDebug.writeImage(gradientXY, "_after_linethinning_2_");
            }
        }
    }

    /**
     * construct the gradient in X, gradient in Y, their combined average and
     * theta image
     * using two 1-D passes of a Sobel 1D kernel which is the same as a
     * Gaussian first derivative with sigma = sqrt(1)/2 where FWHM=2.355*sigma.
     * The theta image has range 0 t 360.
     *
     * @param img
     * @return
     */
    protected EdgeFilterProducts createGradient(final ImageExt img) {

        ImageProcessor imageProcessor = new ImageProcessor();

        int n = img.getNPixels();
        int w = img.getWidth();
        int h = img.getHeight();

        float jnd = 2.3f;
        
        // using deltaE between pixels with coords -1 and +1 of 
        // center pixels.

        // the maximum of any point should be math.sqrt(2)*19.22     
        float[][] gradients = imageProcessor
            .calculateGradientUsingDeltaE2000(img);
          
        //TODO: consider what a good scale range would be for these
        
        //float scale1 = 19.22f;
        //float scale = 255.f/(float)(19.22f * Math.sqrt(2));
        float scale1 = 1;
        float scale = 1;
        
        GreyscaleImage gX, gY, g, theta;
        
        gX = new GreyscaleImage(w, h, GreyscaleImage.Type.Bits32FullRangeInt);
        gY = new GreyscaleImage(w, h, GreyscaleImage.Type.Bits32FullRangeInt);
        g = new GreyscaleImage(w, h, GreyscaleImage.Type.Bits32FullRangeInt);
        
        for (int i = 0; i < n; ++i) {
            float vX = gradients[0][i] * scale1;
            float vY = gradients[1][i] * scale1;
            float vXY = gradients[2][i] * scale;
            
            gX.setValue(i, Math.round(vX));
            gY.setValue(i, Math.round(vY));
            g.setValue(i, Math.round(vXY));
        }
            
        // the theta is in range 0 to 180
        theta = imageProcessor.computeTheta180(gX, gY);

        if (debug) {
            // when replace the aspect library, put these renders in the
            //   equivalent replacement
            long ts = MiscDebug.getCurrentTimeFormatted();
            MiscDebug.writeImage(theta, "_theta_" + ts);
            MiscDebug.writeImage(gX, "_gX_" + ts);
            MiscDebug.writeImage(gY, "_gY_" + ts);
            /*
            int x = 37; int y = 163;
            System.out.println("(" + x + ", " + y + ") math.atan2(" + gY.getValue(x, y)
                + "," + gX.getValue(x, y) + ")*180./math.pi=" +
                theta.getValue(x, y));*/
        }

        EdgeFilterProducts efp = new EdgeFilterProducts();
        efp.setGradientX(gX);
        efp.setGradientY(gY);
        efp.setGradientXY(g);
        efp.setTheta(theta);

        return efp;
    }

    private void applyNonMaximumSuppression(EdgeFilterProducts filterProducts,
        Set<PairInt> disconnectingRemovals) {

         NonMaximumSuppression nms = new NonMaximumSuppression();

         //TODO: radius can be adjusted from 1 to higher
         nms.nonmaxsup(filterProducts.getGradientXY(),
            filterProducts.getTheta(), 1.2, disconnectingRemovals);
    }

    private boolean isAdjacentToAHorizOrVertLine(GreyscaleImage gradientXY,
        int x, int y, int minLineSize) {

        /*
           1  1  1
           1  0  1
           1  1  1

        if find a neighbor, check it for being part of _ or | of size minLineSize
        */
        int w = gradientXY.getWidth();
        int h = gradientXY.getHeight();
        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;

        for (int k = 0; k < dxs.length; ++k) {
            int x2 = x + dxs[k];
            int y2 = y + dys[k];
            if ((x2 < 0) || (y2 < 0) || (x2 > (w - 1)) || (y2 > (h - 1))) {
                continue;
            }
            if (gradientXY.getValue(x2, y2) == 0) {
                continue;
            }

            /*
            if (x2, y2) is part of a vertical line, return true

              (x2, y2)
                 *
                 *
            */
            int c = 0;
            for (int dy2 = ((-minLineSize) + 1); dy2 < 0; ++dy2) {
                int y3 = y2 + dy2;
                if ((y3 < 0) || (y3 > (h - 1)) || (gradientXY.getValue(x2, y3) == 0)) {
                    continue;
                }
                c++;
            }
            if (c == (minLineSize - 1)) {
                return true;
            }
            /*
                   *
                (x2, y2)
                   *
            */
            c = 0;
            for (int dy2 = (-minLineSize/2); dy2 <= (minLineSize/2); ++dy2) {
                int y3 = y2 + dy2;
                if ((y3 < 0) || (y3 > (h - 1)) || (gradientXY.getValue(x2, y3) == 0)
                    || (y3 == y2)) {
                    continue;
                }
                c++;
            }
            if (c == (minLineSize - 1)) {
                return true;
            }
            /*
                   *
                   *
               (x2, y2)
            */
            c = 0;
            for (int dy2 = 1; dy2 < minLineSize; ++dy2) {
                int y3 = y2 + dy2;
                if ((y3 < 0) || (y3 > (h - 1)) || (gradientXY.getValue(x2, y3) == 0)) {
                    continue;
                }
                c++;
            }
            if (c == (minLineSize - 1)) {
                return true;
            }

            /*
            if (x2, y2) is part of a horizontal line, return true

              *   *   (x2, y2)
            */
            c = 0;
            for (int dx2 = ((-minLineSize) + 1); dx2 < 0; ++dx2) {
                int x3 = x2 + dx2;
                if ((x3 < 0) || (x3 > (w - 1)) || (gradientXY.getValue(x3, y2) == 0)) {
                    continue;
                }
                c++;
            }
            if (c == (minLineSize - 1)) {
                return true;
            }

            /*
                  *  (x2, y2)   *
            */
            c = 0;
            for (int dx2 = (-minLineSize/2); dx2 <= (minLineSize/2); ++dx2) {
                int x3 = x2 + dx2;
                if ((x3 < 0) || (x3 > (w - 1)) || (gradientXY.getValue(x3, y2) == 0)
                    || (x3 == x2)) {
                    continue;
                }
                c++;
            }
            if (c == (minLineSize - 1)) {
                return true;
            }

            /*
                  (x2, y2)   *    *
            */
            c = 0;
            for (int dx2 = 1; dx2 < minLineSize; ++dx2) {
                int x3 = x2 + dx2;
                if ((x3 < 0) || (x3 > (w - 1)) || (gradientXY.getValue(x3, y2) == 0)) {
                    continue;
                }
                c++;
            }
            if (c == (minLineSize - 1)) {
                return true;
            }
        }

        return false;
    }

    /**
     * get the filter products for gradient and orientation, and hough
     * lines.  note that the orientation image has values between 0 and 180.
     * @return the filterProducts
     */
    public EdgeFilterProducts getFilterProducts() {
        return filterProducts;
    }

    private void applyPostLineThinningCorrections(GreyscaleImage gradientXY,
        GreyscaleImage valuesBeforeThinning) {

        Set<PairInt> correctedPoints = new HashSet<PairInt>();

        for (int i = 0; i < gradientXY.getWidth(); ++i) {
            for (int j = 0; j < gradientXY.getHeight(); ++j) {
                if (gradientXY.getValue(i, j) > 0) {
                    correctedPoints.add(new PairInt(i, j));
                }
            }
        }

        int n0 = gradientXY.getWidth();
        int n1 = gradientXY.getHeight();

        // if there are too many points in correctedPoints, then
        // the image is probably still filled with pixels noise or textures
        // so do not perform line thinning in that case.
        int nP = correctedPoints.size();
        float frac = (float)nP/(float)(n0 * n1);
        if (nP < 30000) {
            PostLineThinnerCorrections pltc = new PostLineThinnerCorrections();
 //           pltc.correctForHolePattern100(correctedPoints, n0, n1);
            pltc.correctForLineHatHoriz(correctedPoints, n0, n1);
            pltc.correctForLineHatVert(correctedPoints, n0, n1);
            pltc.correctForLineSpurHoriz(correctedPoints, n0, n1);
            pltc.correctForLineSpurVert(correctedPoints, n0, n1);
            pltc.correctForIsolatedPixels(correctedPoints, n0, n1);
        }

        GreyscaleImage out = gradientXY.createWithDimensions();
        int[][] morphInput = new int[out.getWidth()][out.getHeight()];
        for (int i = 0; i < morphInput.length; ++i) {
            morphInput[i] = new int[out.getHeight()];
        }
        for (PairInt p : correctedPoints) {
            int v = valuesBeforeThinning.getValue(p);
            out.setValue(p.getX(), p.getY(), v);
            morphInput[p.getX()][p.getY()] = 1;
        }

        MorphologicalFilter mFilter = new MorphologicalFilter();
        int[][] skel = mFilter.bwMorphThin(morphInput, Integer.MAX_VALUE);
        for (int i = 0; i < n0; ++i) {
            for (int j = 0; j < n1; ++j) {
                int m = skel[i][j];
                if (m == 0) {
                    out.setValue(i, j, 0);
                }
            }
        }
        if (useZhangSuen) {
            ZhangSuenLineThinner lt = new ZhangSuenLineThinner();
            lt.applyFilter(out);
        }

        gradientXY.resetTo(out);
    }

    private int countAboveZero(GreyscaleImage img) {

        int n = 0;
        for (int i = 0; i < img.getNPixels(); ++i) {
            if (img.getValue(i) > 0) {
                n++;
            }
        }

        return n;
    }

    private void applyHoughBasedLineThinning(EdgeFilterProducts products) {

        GreyscaleImage gradientXY = products.getGradientXY();

        int n = gradientXY.getNPixels();
        int w = gradientXY.getWidth();
        int h = gradientXY.getHeight();

        Set<PairInt> points = new HashSet<PairInt>();

        for (int i = 0; i < n; ++i) {
            if (gradientXY.getValue(i) > 0) {
                points.add(new PairInt(gradientXY.getCol(i), gradientXY.getRow(i)));
            }
        }
        Set<PairInt> pointsCp = new HashSet<PairInt>(points);

        HoughTransform ht = new HoughTransform();
        Map<Set<PairInt>, PairInt> lines = ht.findContiguousLines(points, 3);

        products.setHoughLines(lines);

        PostLineThinnerCorrections pltc = new PostLineThinnerCorrections();
        pltc.thinLineStaircases(lines, points, w, h);
        pltc.correctForLineSpurHoriz(points, w, h);
        pltc.correctForLineSpurVert(points, w, h);
        pltc.correctForLine2SpurVert(points, w, h);
        pltc.correctForLine2SpurHoriz(points, w, h);
        pltc.correctForIsolatedPixels(points, w, h);

        pointsCp.removeAll(points);

        for (PairInt p : pointsCp) {
            gradientXY.setValue(p.getX(), p.getY(), 0);

            Set<PairInt> line = null;
            for (Set<PairInt> hLine : lines.keySet()) {
                if (hLine.contains(p)) {
                    line = hLine;
                    break;
                }
            }
            if (line != null) {
                line.remove(p);
            }
        }

    }

}
