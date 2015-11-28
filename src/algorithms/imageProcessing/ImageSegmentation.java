package algorithms.imageProcessing;

import algorithms.CountingSort;
import algorithms.compGeometry.NearestPoints1D;
import algorithms.compGeometry.clustering.KMeansPlusPlus;
import algorithms.compGeometry.clustering.KMeansPlusPlusFloat;
import algorithms.imageProcessing.util.PairIntWithIndex;
import algorithms.misc.Histogram;
import algorithms.misc.HistogramHolder;
import algorithms.misc.Misc;
import algorithms.misc.MiscDebug;
import algorithms.misc.MiscMath;
import algorithms.util.Errors;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.ResourceFinder;
import com.climbwithyourfeet.clustering.DTClusterFinder;
import java.io.IOException;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.Stack;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * class holding several different image segmentation methods.  Note that
 * some other techniques involving contrast for example, are elsewhere.
 *
 * @author nichole
 */
public class ImageSegmentation {

    private Logger log = Logger.getLogger(this.getClass().getName());

    /**
     * applies KMeansPlusPlus algorithm to the values in input
     * (greyscale intensities) to create kBands of clustered pixels
     * (operates on input).
     * @param input
     * @param kBands
     * @throws IOException
     * @throws NoSuchAlgorithmException
     */
    public void applyUsingKMPP(GreyscaleImage input, int kBands)
        throws IOException, NoSuchAlgorithmException {

        KMeansPlusPlus instance = new KMeansPlusPlus();
        instance.computeMeans(kBands, input);

        assignToNearestCluster(input, instance.getCenters());
    }

    /**
     * applies binary algorithm (simple thresholding) to the values in input
     * (greyscale intensities) to create pixels given above or below highest
     * frequency value.  Note that some images may need to be pre-processed
     * in order to use this one (for example, correct for illumination and
     * remove items like sky if main objects are not sky).
     * (operates on input).
     * @param input
     */
    public void applyBinaryUsingFrequency(GreyscaleImage input) {

        PairIntArray valueCounts = Histogram.createADescendingSortbyFrequencyArray(
            input);

        if (valueCounts == null || valueCounts.getN() == 0) {
            return;
        }

        //96, 30
        int divider = valueCounts.getX(0);
        int v0 = 255/4;
        int v1 = 3*v0;

        for (int i = 0; i < input.getNPixels(); ++i) {
            if (input.getValue(i) < divider) {
                input.setValue(i, v0);
            } else {
                input.setValue(i, v1);
            }
        }

    }

    /**
     * applies DTClustering algorithm to the values in input
     * (greyscale intensities) to create kBands of clustered pixels
     * (operates on input).
     * (This one is competitive with applyUsingPolarCIEXYAndFrequency
     * with lowFreqLimit 0.1f)
     * @param input
     * @param kBands
     * @throws IOException
     * @throws NoSuchAlgorithmException
     */
    public void applyUsingDTClustering(GreyscaleImage input, int kBands)
        throws IOException, NoSuchAlgorithmException {

        PairIntArray valueCounts = Histogram.createADescendingSortbyFrequencyArray(
            input);

        // first, trying clustering by value and frequency,
        // second, compare that to clustering just by giving points as value,value
        //    which is effectivly a 1D clustering
        //Set<PairInt> points = new HashSet<PairInt>();
        Set<com.climbwithyourfeet.clustering.util.PairInt> points =
            new HashSet<com.climbwithyourfeet.clustering.util.PairInt>();

        for (int i = 0; i < valueCounts.getN(); ++i) {
            com.climbwithyourfeet.clustering.util.PairInt p = new
                com.climbwithyourfeet.clustering.util.PairInt(
                    valueCounts.getX(i), valueCounts.getY(i));

            points.add(p);
        }

        DTClusterFinder<com.climbwithyourfeet.clustering.util.PairInt> cFinder
            = new DTClusterFinder<com.climbwithyourfeet.clustering.util.PairInt>(
                points, input.getWidth(), input.getHeight());

        cFinder.calculateCriticalDensity();

        cFinder.findClusters();

        int n = cFinder.getNumberOfClusters();

        int[] centers = new int[n];
        for (int i = 0; i < n; ++i) {

            Set<com.climbwithyourfeet.clustering.util.PairInt> set = cFinder.getCluster(i);

            // find centeroid for x
            double xc = 0;
            for (com.climbwithyourfeet.clustering.util.PairInt p : set) {
                double x1 = p.getX();
                xc += x1;
            }

            centers[i] = (int)Math.round(xc/(double)set.size());
        }

        if (n > kBands) {
            n = kBands;
        }

        Arrays.sort(centers);
        int[] kCenters = new int[n];
        int count = 0;
        for (int i = (centers.length - 1); i > (centers.length - 1 - n); --i) {
            kCenters[count] = centers[i];
            count++;
        }

        assignToNearestCluster(input, kCenters);
    }

    /**
     * places points by their proximity to cluster centers
     * @param input
     * @param binCenters
     */
    public void assignToNearestCluster(GreyscaleImage input, int[] binCenters) {
        
        NearestPoints1D np = new NearestPoints1D(binCenters);
                
        for (int col = 0; col < input.getWidth(); col++) {

            for (int row = 0; row < input.getHeight(); row++) {

                int v = input.getValue(col, row);
                
                int vc = np.findClosestValue(v);
                
                input.setValue(col, row, vc);
            }
        }
    }

    public List<Set<PairInt>> assignToNearestPolarCIECluster(
        Map<PairInt, Integer> polarCIEXYMap, int[] binCenters) {

        Arrays.sort(binCenters);

        int nc = binCenters.length;

        List<Set<PairInt>> groups = new ArrayList<Set<PairInt>>();
        for (int i = 0; i < nc; ++i) {
            groups.add(new HashSet<PairInt>());
        }

        for (Entry<PairInt, Integer> entry : polarCIEXYMap.entrySet()) {

            int theta = entry.getValue().intValue();

            int idx = Arrays.binarySearch(binCenters, theta);
            // if it's negative, (-(insertion point) - 1)
            if (idx < 0) {
                // idx = -*idx2 - 1
                idx = -1*(idx + 1);
            }
            if (idx > (nc - 1)) {
                idx = nc - 1;
            }

            int vc = binCenters[idx];

            if (idx == 0) {

                int bisectorBelowHalfLength = (360 - binCenters[nc - 1] + vc)/2;
                if ((vc - theta) > bisectorBelowHalfLength) {
                    idx = nc - 1;
                }

            } else {

                int bisectorBelow = ((binCenters[idx - 1] + vc) / 2);

                if (theta < bisectorBelow) {
                    idx = idx - 1;
                }
            }

            Set<PairInt> set = groups.get(Integer.valueOf(idx));
            set.add(entry.getKey());
        }

        for (int i = (groups.size() - 1); i > -1; --i) {
            if (groups.get(i).isEmpty()) {
                groups.remove(i);
            }
        }

        return groups;
    }

    /**
     * applies a blur of sigma=1 to image,
     * converts each pixel color to the polar angle of CIE XY Lab color space
     * with an origin of (0.35, 0.35) and uses a histogram binning of kColors=8,
     * then maps those bins to 0 to 255,
     * then replaces a pixel if 5,6 or its neighbors have the same color,
     * then applies histogram equalization to stretch the values to range
     * 0 to 255.
     * @param input
     * @return
     */
    public GreyscaleImage applyUsingCIEXYPolarThetaThenHistEq(ImageExt input) {

        int kColors = 8;

        return applyUsingCIEXYPolarThetaThenHistEq(input, kColors, true);
    }

    /**
     * converts each pixel color to the polar angle of CIE XY Lab color space
     * with an origin of (0.35, 0.35) and uses a histogram binning of kColors,
     * then maps those bins to 0 to 255,
     * then replaces a pixel if 5,6 or its neighbors have the same color,
     * then applies histogram equalization to stretch the values to range
     * 0 to 255.
     * @param input
     * @param kColors the number of colors to bin the image by.  max allowed value
     * is 253.
     * @param useBlur if true, a blur of sigma=1 is applied to the image before
     * processing.
     * @return
     */
    public GreyscaleImage applyUsingCIEXYPolarThetaThenHistEq(ImageExt input,
        int kColors, boolean useBlur) {

        if (kColors > 253) {
            throw new IllegalArgumentException("kColors must be <= 253");
        }
        if (kColors < 2) {
            throw new IllegalArgumentException("kColors must be >= 2");
        }

        ImageProcessor imageProcessor = new ImageProcessor();

        GreyscaleImage img;

        int minNeighborLimit;

        if (useBlur) {

            Image input2 = input.copyImage();

            imageProcessor.blur(input2, 1/*(float)Math.sqrt(2)/2.f*/);

            img = applyUsingCIEXYPolarThetaThenHistogram(input2, kColors);

            minNeighborLimit = 6;

        } else {

            img = applyUsingCIEXYPolarThetaThenHistogram(input, kColors);

            minNeighborLimit = 5;
        }

        int w = img.getWidth();
        int h = img.getHeight();

        // ----replace pixel, if 5,6 or more neighbors have same color -----
        int[] dxs = new int[]{-1, -1,  0,  1, 1, 1, 0, -1};
        int[] dys = new int[]{ 0, -1, -1, -1, 0, 1, 1,  1};

        Map<Integer, Integer> freqMap = new HashMap<Integer, Integer>();

        int nChanged = 1;
        int nIterMax = 100;
        int nIter = 0;

        while (!useBlur && (nIter < nIterMax) && (nChanged > 0)) {

            log.fine("nIter=" + nIter + " nChanged=" + nChanged);

            nChanged = 0;

            for (int col = 0; col < w; col++) {
                for (int row = 0; row < h; row++) {

                    freqMap.clear();

                    Integer maxCountValue = null;
                    int maxCount = Integer.MIN_VALUE;

                    for (int nIdx = 0; nIdx < dxs.length; nIdx++) {
                        int x = dxs[nIdx] + col;
                        int y = dys[nIdx] + row;

                        if ((x < 0) || (x > (w - 1)) || (y < 0) || (y > (h - 1))) {
                            break;
                        }

                        Integer v = Integer.valueOf(img.getValue(x, y));

                        Integer c = freqMap.get(v);
                        if (c == null) {
                            c = Integer.valueOf(1);
                        } else {
                            c = Integer.valueOf(c.intValue() + 1);
                        }
                        freqMap.put(v, c);

                        if (c.intValue() > maxCount) {
                            maxCount = c.intValue();
                            maxCountValue = v;
                        }
                    }

                    if ((maxCount >= minNeighborLimit) &&
                        (img.getValue(col, row) != maxCountValue.intValue())) {

                        img.setValue(col, row, maxCountValue.intValue());
                        nChanged++;
                    }
                }
            }

            nIter++;
        }

        // rescale the image
        HistogramEqualization hEq = new HistogramEqualization(img);
        hEq.applyFilter();

        return img;
    }

    /**
     * applies a blur of sigma=1 to image,
     * converts each pixel's color into CIE XY polar theta, then uses KMeansPlusPlus
     * to create kColors=8 bins points then remaps the points to use the
     * range 0 to 255, then replaces a pixel if it has 7 neighbors of same
     * color, then applies histogram equalization to rescale the range to be
     * between 0 and 255.
     * @param input
     * @return
     */
    public GreyscaleImage applyUsingCIEXYPolarThetaThenKMPPThenHistEq(ImageExt
        input) {

        int kColors = 8;

        return applyUsingCIEXYPolarThetaThenKMPPThenHistEq(input, kColors, true);
    }

    /**
     * converts each pixel's color into CIE XY polar theta, then uses KMeansPlusPlus
     * to create kColors bins points then remaps the points to use the
     * range 0 to 255, then replaces a pixel if it has 5,6 neighbors of same
     * color, then applies histogram equalization to rescale the range to be
     * between 0 and 255.
     *
     * @param input
     * @param kColors the number of colors to bin the image by.  max allowed value
     * is 253.
     * @param useBlur
     * @return
     */
    public GreyscaleImage applyUsingCIEXYPolarThetaThenKMPPThenHistEq(ImageExt input,
        int kColors, boolean useBlur) {

        if (kColors > 253) {
            throw new IllegalArgumentException("kColors must be <= 253");
        }
        if (kColors < 2) {
            throw new IllegalArgumentException("kColors must be >= 2");
        }

        ImageProcessor imageProcessor = new ImageProcessor();

        GreyscaleImage img;

        int minNeighborLimit;

        if (useBlur) {

            Image input2 = input.copyImage();

            imageProcessor.blur(input2, 1/*(float)Math.sqrt(2)/2.f*/);

            img = applyUsingCIEXYPolarThetaThenKMPP(input2, kColors);

            minNeighborLimit = 6;

        } else {

            img = applyUsingCIEXYPolarThetaThenKMPP(input, kColors);

            minNeighborLimit = 5;
        }

        int w = img.getWidth();
        int h = img.getHeight();

        // ----replace pixel, if 5,6 or more neighbors have same color -----
        int[] dxs = new int[]{-1, -1,  0,  1, 1, 1, 0, -1};
        int[] dys = new int[]{ 0, -1, -1, -1, 0, 1, 1,  1};

        Map<Integer, Integer> freqMap = new HashMap<Integer, Integer>();

        int nChanged = 1;
        int nIterMax = 100;
        int nIter = 0;

        while (!useBlur && (nIter < nIterMax) && (nChanged > 0)) {

            log.fine("nIter=" + nIter + " nChanged=" + nChanged);

            nChanged = 0;

            for (int col = 0; col < w; col++) {
                for (int row = 0; row < h; row++) {

                    freqMap.clear();

                    Integer maxCountValue = null;
                    int maxCount = Integer.MIN_VALUE;

                    for (int nIdx = 0; nIdx < dxs.length; nIdx++) {
                        int x = dxs[nIdx] + col;
                        int y = dys[nIdx] + row;

                        if ((x < 0) || (x > (w - 1)) || (y < 0) || (y > (h - 1))) {
                            break;
                        }

                        Integer v = Integer.valueOf(img.getValue(x, y));

                        Integer c = freqMap.get(v);
                        if (c == null) {
                            c = Integer.valueOf(1);
                        } else {
                            c = Integer.valueOf(c.intValue() + 1);
                        }
                        freqMap.put(v, c);

                        if (c.intValue() > maxCount) {
                            maxCount = c.intValue();
                            maxCountValue = v;
                        }
                    }

                    if ((maxCount >= minNeighborLimit) &&
                        (img.getValue(col, row) != maxCountValue.intValue())) {

                        img.setValue(col, row, maxCountValue.intValue());
                        nChanged++;
                    }
                }
            }

            nIter++;
        }

        // rescale the image
        HistogramEqualization hEq = new HistogramEqualization(img);
        hEq.applyFilter();

        return img;
    }

    /**
     * converts each pixel's color into CIE XY polar theta,
     * then applies histogram mapping of kColors to remap the pixels to
     * values between 0 and 255.
     *
     * @param input
     * @return
     */
    public GreyscaleImage applyUsingCIEXYPolarThetaThenHistogram(Image input) {

        return applyUsingCIEXYPolarThetaThenHistogram(input, 253);
    }

    /**
     * converts each pixel's color into CIE XY polar theta,
     * then applies histogram mapping of kColors to remap the pixels to
     * values between 0 and 255.
     *
     * @param input
     * @param kColors the number of color bins to use for the image segmentation.
     * The minimum allowed value is 2 and the maximum allowed value is 253.
     * @return
     */
    public GreyscaleImage applyUsingCIEXYPolarThetaThenHistogram(Image input,
        int kColors) {

        if (kColors > 253) {
            throw new IllegalArgumentException("kColors must be <= 253");
        }
        if (kColors < 2) {
            throw new IllegalArgumentException("kColors must be >= 2");
        }

        /*
        TODO: needs some improvements in color mapping.
           -- for regions that are very spotty, might consider using the
              intensity image to help define the region and take the highest
              number density color in that region and assign it to all.
        */

        int w = input.getWidth();
        int h = input.getHeight();

        float[] tmpColorBuffer = new float[2];

        GreyscaleImage output = new GreyscaleImage(w, h);

        Map<PairInt, Float> pixThetaMap = new HashMap<PairInt, Float>();

        CIEChromaticity cieC = new CIEChromaticity();

        float[] thetaValues = new float[input.getNPixels()];
        int thetaCount = 0;

        for (int col = 0; col < w; col++) {
            for (int row = 0; row < h; row++) {

                PairInt p = new PairInt(col, row);

                int r = input.getR(col, row);
                int g = input.getG(col, row);
                int b = input.getB(col, row);

                if ((r < 25) && (g < 25) && (b < 25)) {
                    continue;
                } else if ((r > 230) && (g > 230) && (b > 230)) {//might need to use 195 as lower limit
                    output.setValue(col, row, 255);
                    continue;
                }

                float[] cieXY = tmpColorBuffer;
                cieC.rgbToXYChromaticity(r, g, b, cieXY);

                if (cieC.isWhite(cieXY[0], cieXY[1])) {
                    output.setValue(col, row, 255);
                } else {

                    double thetaRadians = cieC.calculateXYTheta(cieXY[0], cieXY[1]);

                    thetaValues[thetaCount] = (float)thetaRadians;

                    pixThetaMap.put(p, Float.valueOf((float)thetaRadians));

                    thetaCount++;
                }
            }
        }

        thetaValues = Arrays.copyOf(thetaValues, thetaCount);

        createAndApplyHistMapping(output, pixThetaMap, thetaValues, kColors);

        return output;
    }

    /**
     * converts each pixel's color into CIE XY polar theta, then uses KMeansPlusPlus
     * to create kColors bins points then remaps the points to use the
     * range 0 to 255.
     *
     * runtime complexity is O(N) + O(N*lg_2(N))
     *
     * @param input
     * @param kColors the number of color bins to use for the image segmentation.
     * The minimum allowed value is 2 and the maximum allowed value is 253.
     * @return
     */
    public GreyscaleImage applyUsingCIEXYPolarThetaThenKMPP(Image input, int kColors) {

        if (kColors > 253) {
            throw new IllegalArgumentException("kColors must be <= 253");
        }
        if (kColors < 2) {
            throw new IllegalArgumentException("kColors must be >= 2");
        }

        /*
        TODO: needs some improvements in color mapping.
           -- for regions that are very spotty, might consider using the
              intensity image to help define the region and take the highest
              number density color in that region and assign it to all.
        */

        int w = input.getWidth();
        int h = input.getHeight();

        float[] tmpColorBuffer = new float[2];

        GreyscaleImage output = new GreyscaleImage(w, h);

        Map<PairInt, Float> pixThetaMap = new HashMap<PairInt, Float>();

        CIEChromaticity cieC = new CIEChromaticity();

        float[] thetaValues = new float[input.getNPixels()];
        int thetaCount = 0;

        for (int col = 0; col < w; col++) {
            for (int row = 0; row < h; row++) {

                PairInt p = new PairInt(col, row);

                int r = input.getR(col, row);
                int g = input.getG(col, row);
                int b = input.getB(col, row);

                if ((r < 25) && (g < 25) && (b < 25)) {
                    continue;
                } else if ((r > 230) && (g > 230) && (b > 230)) {//might need to use 195 as lower limit
                    output.setValue(col, row, 255);
                    continue;
                }

                float[] cieXY = tmpColorBuffer;
                cieC.rgbToXYChromaticity(r, g, b, cieXY);

                if (cieC.isWhite(cieXY[0], cieXY[1])) {
                    output.setValue(col, row, 255);
                } else {

                    double thetaRadians = cieC.calculateXYTheta(cieXY[0], cieXY[1]);

                    thetaValues[thetaCount] = (float)thetaRadians;

                    pixThetaMap.put(p, Float.valueOf((float)thetaRadians));

                    thetaCount++;
                }
            }
        }

        thetaValues = Arrays.copyOf(thetaValues, thetaCount);

        createAndApplyKMPPMapping(output, pixThetaMap, thetaValues, kColors);

        return output;
    }

    private void createAndApplyKMPPMapping(GreyscaleImage output,
        Map<PairInt, Float> pixThetaMap, float[] thetaValues,
        final int kColors) {

        //TODO: assert kColors.  The invoker is reserving 2 bands for
        // B & W, so nBins should probably be (kColors - 2)...
        // correct this for the invoker when testing
        int nBins = kColors;

        KMeansPlusPlusFloat kmpp = new KMeansPlusPlusFloat();
        kmpp.computeMeans(nBins, thetaValues);

        float minValue = kmpp.getMinValue();
        float maxValue = kmpp.getMaxValue();

        float[] binCenters = kmpp.getCenters();

        Iterator<Map.Entry<PairInt, Float> > iter = pixThetaMap.entrySet().iterator();

        while (iter.hasNext()) {

            Map.Entry<PairInt, Float> entry = iter.next();

            PairInt p = entry.getKey();

            float theta = entry.getValue().floatValue();

            for (int i = 0; i < binCenters.length; i++) {

                float vc = binCenters[i];

                float bisectorBelow = ((i - 1) > -1) ?
                    ((binCenters[i - 1] + vc) / 2) : minValue;

                float bisectorAbove = ((i + 1) > (binCenters.length - 1)) ?
                    maxValue : ((binCenters[i + 1] + vc) / 2);

                if ((theta >= bisectorBelow) && (theta <= bisectorAbove)) {

                    //TODO: check this
                    int mappedValue = 255 - nBins + i;

                    output.setValue(p.getX(), p.getY(), mappedValue);

                    break;
                }
            }

            /*
            // if binCenters is ordered, use binary search for faster results
            int idx = Arrays.binarySearch(startBins, theta);

            // if it's negative, (-(insertion point) - 1)
            if (idx < 0) {
                // idx = -*idx2 - 1
                idx = -1*(idx + 1);
            }
            int mappedValue = 255 - startBins.length + idx;

            output.setValue(p.getX(), p.getY(), mappedValue);
            */
        }
    }

    private void createAndApplyHistMapping(GreyscaleImage output,
        Map<PairInt, Float> pixThetaMap, float[] thetaValues,
        final int kColors) {

        float minValue = MiscMath.findMin(thetaValues);
        float maxValue = MiscMath.findMax(thetaValues);

        log.fine("minTheta=" + (minValue * 180./Math.PI) +
            " maxTheta=" + (maxValue * 180./Math.PI));

        int nReserved = 254 - kColors;

        HistogramHolder hist = Histogram.createSimpleHistogram(minValue,
            maxValue, (256 - nReserved - 1), thetaValues,
            Errors.populateYErrorsBySqrt(thetaValues));

        try {
            hist.plotHistogram("cie XY theta histogram", "cieXY_hist_"
                + MiscDebug.getCurrentTimeFormatted());
        } catch (Exception e) {}

        int nonZeroCount = 0;
        for (int i = 0; i < hist.getXHist().length; i++) {
            int c = hist.getYHist()[i];
            if (c > 0) {
                nonZeroCount++;
            }
        }

        float[] startBins = new float[nonZeroCount];

        float halfBinWidth = (hist.getXHist()[1] - hist.getXHist()[0])/2.f;

        nonZeroCount = 0;
        for (int i = 0; i < hist.getXHist().length; i++) {
            int c = hist.getYHist()[i];
            if (c > 0) {
                startBins[nonZeroCount] = hist.getXHist()[i] - halfBinWidth;
                nonZeroCount++;
            }
        }

        Iterator<Map.Entry<PairInt, Float> > iter = pixThetaMap.entrySet().iterator();

        // O(N * lg_2(N))
        while (iter.hasNext()) {

            Map.Entry<PairInt, Float> entry = iter.next();

            PairInt p = entry.getKey();

            float theta = entry.getValue().floatValue();

            int idx = Arrays.binarySearch(startBins, theta);

            // if it's negative, (-(insertion point) - 1)
            if (idx < 0) {
                // idx = -*idx2 - 1
                idx = -1*(idx + 1);
            }

            int mappedValue = 255 - startBins.length + idx;

            output.setValue(p.getX(), p.getY(), mappedValue);
        }
    }

    /**
     * NOT READY FOR USE.  STILL EXPERIMENTING.
     *
     * Calculates lists of black pixels, white pixels, grey pixels, and assigns
     * the remaining to CIEXY Lab color space, then creates a map of the
     * CIEX and CIEY points and uses density based clustering
     * (http://nking.github.io/two-point-correlation/)
     * to find clusters of points in CIE X, CIEY space,
     * then merges pixels in the grey list with adjacent clusters if
     * similar, and the final result is a list of pixel clusters, including the
     * black and white and remaining grey.
     *
     * Note that the color black is not defined in CIE XY color space and that
     * the color white is at the center of the space as a large circle so they
     * are not included in the density based clustering.
     * Note also that the remaining grey pixels which did not merge with the
     * cie xy pixel clusters, may be in separate groups already due to a
     * frequency based grouping for them.
     *
     * @param input
     * @param useBlur
     * @return
     */
    public List<Set<PairInt>> calculateUsingCIEXYAndClustering(ImageExt input,
        boolean useBlur) {

        //TODO: improve the clustering results in two ways:
        // (1) for smaller ciexy clusters, merge with adjacent clusters if
        //     similar color
        // (2) any pixel with 7 neighbors of same color should be that color too

        if (useBlur) {
            ImageProcessor imageProcessor = new ImageProcessor();
            imageProcessor.blur(input, 1.0f);
        }

        //TODO: consider making a segmentation method using CIEXY theta
        // for x and the frequency for y, both scaled to numerically
        // resolvable range < max of 5000.
        // this would be good to compare to the method here which
        // uses CIE XY Theta followed by histograms or KMeans++.
        // No need to specify the number of bins before use for suggested
        // version.

        //NOTE: the method needs to have gaps in the data given to it
        //    that is a lack of points for some region between the
        //    min and max of x and y data in integer space

        // max = 6250 unless reduce space complexity
        float factor = 2000;// learn this from numerical resolution

        // then subtract the minima in both cieX and cieY

        int minCIEX = Integer.MAX_VALUE;
        int minCIEY = Integer.MAX_VALUE;
        int maxCIEX = Integer.MIN_VALUE;
        int maxCIEY = Integer.MIN_VALUE;

        Set<PairInt> blackPixels = new HashSet<PairInt>();

        Map<Integer, Collection<PairInt>> greyPixelMap = new HashMap<Integer, Collection<PairInt>>();

        Set<PairInt> whitePixels = new HashSet<PairInt>();

        Set<PairInt> points0 = new HashSet<PairInt>();

        populatePixelLists(input, points0, blackPixels, whitePixels, greyPixelMap);

        // -------- debug -------
        int nGrey = 0;
        for (Map.Entry<Integer, Collection<PairInt>> entry : greyPixelMap.entrySet()) {
            nGrey += entry.getValue().size();
        }
        assert((points0.size() + blackPixels.size() + nGrey +
            whitePixels.size()) == input.getNPixels());
        // -------- end debug -------

        List<Set<PairInt>> greyPixelGroups = groupByPeaks(greyPixelMap);

        // ------- debug -------
        int nGrey2 = 0;
        for (Set<PairInt> set : greyPixelGroups) {
            nGrey2 += set.size();
        }
        assert(nGrey == nGrey2);
        // ------- end debug =====

        Map<PairIntWithIndex, List<PairIntWithIndex>> pointsMap0 =
            new HashMap<PairIntWithIndex, List<PairIntWithIndex>>();

        for (PairInt p : points0) {
            int idx = input.getInternalIndex(p.getX(), p.getY());
            float cx = input.getCIEX(idx);
            float cy = input.getCIEY(idx);

            int cieXInt = Math.round(factor * cx);
            int cieYInt = Math.round(factor * cy);

            PairIntWithIndex p0 = new PairIntWithIndex(cieXInt, cieYInt, idx);
            List<PairIntWithIndex> list = pointsMap0.get(p0);
            if (list == null) {
                list = new ArrayList<PairIntWithIndex>();
                pointsMap0.put(p0, list);
            }
            list.add(p0);

            if (cieXInt < minCIEX) {
                minCIEX = cieXInt;
            }
            if (cieYInt < minCIEY) {
                minCIEY = cieYInt;
            }
            if (cieXInt > maxCIEX) {
                maxCIEX = cieXInt;
            }
            if (cieYInt > maxCIEY) {
                maxCIEY = cieYInt;
            }
        }

        Map<PairIntWithIndex, List<PairIntWithIndex>> pointsMap =
            new HashMap<PairIntWithIndex, List<PairIntWithIndex>>();

        // subtract minima from the points
        for (PairIntWithIndex p : pointsMap0.keySet()) {

            int x = p.getX() - minCIEX;
            int y = p.getY() - minCIEY;

            PairIntWithIndex p2 = new PairIntWithIndex(x, y, p.getPixIndex());
            List<PairIntWithIndex> list2 = pointsMap.get(p2);
            if (list2 == null) {
                list2 = new ArrayList<PairIntWithIndex>();
                pointsMap.put(p2, list2);
            }
            // because this is a list, this will eventually be present twice:
            //list2.add(p2);

            for (PairIntWithIndex p0 : pointsMap0.get(p)) {
                PairIntWithIndex p3 = new PairIntWithIndex(
                    p0.getX() - minCIEX, p0.getY() - minCIEY, p0.getPixIndex());
                list2.add(p3);
            }
        }
        maxCIEX -= minCIEX;
        maxCIEY -= minCIEY;

        // frequency of colors:
        Map<PairIntWithIndex, Integer> freqMap = new
            HashMap<PairIntWithIndex, Integer>();
        for (Map.Entry<PairIntWithIndex, List<PairIntWithIndex>> entry :
            pointsMap.entrySet()) {
            int c = entry.getValue().size();
            freqMap.put(entry.getKey(), Integer.valueOf(c));
        }

        // ----- debug ---
        int nGreyBW = nGrey + blackPixels.size() + whitePixels.size();
        int nTot = 0;
        //Map<PairIntWithIndex, List<PairIntWithIndex>> pointsMap
        for (Entry<PairIntWithIndex, List<PairIntWithIndex>> entry : pointsMap.entrySet()) {
            nTot += entry.getValue().size();
        }
        nTot += nGreyBW;
        log.info("nTot=" + nTot + " nPixels=" + input.getNPixels());
        assert(nTot == input.getNPixels());

        // plot the points as an image to see the data first
        GreyscaleImage img = new GreyscaleImage(maxCIEX + 1, maxCIEY + 1);
        for (com.climbwithyourfeet.clustering.util.PairInt p : pointsMap.keySet()) {
            img.setValue(p.getX(), p.getY(), 255);
        }
        try {
            ImageIOHelper.writeOutputImage(
                ResourceFinder.findDirectory("bin") + "/dt_input.png", img);
        } catch (IOException ex) {
            Logger.getLogger(ImageProcessor.class.getName()).log(Level.SEVERE,
                null, ex);
        }
        // --- end debug

        //Map<PairIntWithIndex, List<PairIntWithIndex>> pointsMap

        DTClusterFinder<PairIntWithIndex> clusterFinder
            = new DTClusterFinder<PairIntWithIndex>(pointsMap.keySet(),
            maxCIEX + 1, maxCIEY + 1);

        clusterFinder.setToDebug();

        // to recover every point, set limit to 1
        clusterFinder.setMinimumNumberInCluster(1);

        clusterFinder.calculateCriticalDensity();

        clusterFinder.findClusters();

        log.info("clustering critical density=" + clusterFinder.getCriticalDensity());

        int nGroups = clusterFinder.getNumberOfClusters();

        List<Set<PairInt>> groupList = new ArrayList<Set<PairInt>>();

        for (int k = 0; k < nGroups; ++k) {

            Set<PairIntWithIndex> group = clusterFinder.getCluster(k);

            Set<PairInt> coordPoints = new HashSet<PairInt>();

            for (PairIntWithIndex p : group) {

                int idx = p.getPixIndex();
                int xCoord = input.getCol(idx);
                int yCoord = input.getRow(idx);

                PairInt pCoord = new PairInt(xCoord, yCoord);
                coordPoints.add(pCoord);

                // include the other points of/ same color
                List<PairIntWithIndex> list = pointsMap.get(p);
                assert(list != null);
                for (PairIntWithIndex p3 : list) {
                    int idx3 = p3.getPixIndex();
                    int xCoord3 = input.getCol(idx3);
                    int yCoord3 = input.getRow(idx3);
                    pCoord = new PairInt(xCoord3, yCoord3);
                    coordPoints.add(pCoord);
                }
            }

            groupList.add(coordPoints);
        }

        // ------ debug ---------
        nTot = 0;
        for (Set<PairInt> set : groupList) {
            nTot += set.size();
        }
        nTot += nGreyBW;
        log.info("nTot=" + nTot + " nPixels=" + input.getNPixels());
        assert(nTot == input.getNPixels());
        // ------ end debug -----

        mergeOrAppendGreyWithOthers(input, greyPixelGroups, groupList,
            blackPixels, whitePixels);

        // add back in blackPixels and whitePixels
        groupList.add(blackPixels);
        groupList.add(whitePixels);

        int nTot2 = 0;
        for (Set<PairInt> groups : groupList) {
            nTot2 += groups.size();
        }

        log.info("img nPix=" + input.getNPixels() + " nTot2=" + nTot2);
        assert(nTot2 == input.getNPixels());

        return groupList;
    }

    /**
     * Calculates lists of black pixels, white pixels, grey pixels, and assigns
     * the remaining to the polar angle of CIEXY Lab color space, then creates a map of the
     * polar angles and the number of pixels with those angles (=frequency)
     * and uses density based clustering
     * (http://nking.github.io/two-point-correlation/)
     * to find clusters in that space (polar CIEXY vs frequency),
     * then merges pixels in the grey list with adjacent clusters if
     * similar, and the final result is a list of pixel clusters, including the
     * black and white and remaining grey.
     * The changeable parameters are the scaling of the range of polar angles
     * and the range of frequency in order to place the data between values of 0
     * and a number.
     * <code>The scale factors are double xFactor = 2000. and int yFactor = 2000
     * for example.  Smaller scale factors result in faster runtimes,
     * but must be balanced by a large enough factor to have cluster resulution.
     * </code>
     * Note that the polar angle vs frequency maps are actually partitioned into
     * 4 maps to do density based cluster finding separately.
     * partitionFreqFracs = new float[]{0.03f, 0.15f, 0.25f} is the fraction of
     * the maximum frequency defining a partition.
     *
     * Note that the color black is not defined in CIE XY color space and that
     * the color white is at the center of the space as a large circle so they
     * are not included in the density based clustering.
     * Note also that the remaining grey pixels which did not merge with the
     * CIE xy pixel clusters, may be in separate groups already due to a
     * frequency based grouping for them.
     *
     * @param input
     * @param useBlur apply a gaussian blur of sigma=1 before the method logic
     * @return
     */
    public List<Set<PairInt>> calculateUsingPolarCIEXYAndClustering(ImageExt input,
        boolean useBlur) {

        if (useBlur) {
            ImageProcessor imageProcessor = new ImageProcessor();
            imageProcessor.blur(input, 1.0f);
        }

        //making a segmentation method using CIEXY theta
        // for x and the frequency for y, both scaled to numerically
        // resolvable range < max of 5000.
        // this would be good to compare to the method here which
        // uses CIE XY Theta followed by histograms or KMeans++.
        // No need to specify the number of bins before use for suggested
        // version.

        int w = input.getWidth();
        int h = input.getHeight();

        Set<PairInt> blackPixels = new HashSet<PairInt>();

        Map<Integer, Collection<PairInt>> greyPixelMap = new HashMap<Integer, Collection<PairInt>>();

        Set<PairInt> whitePixels = new HashSet<PairInt>();

        Set<PairInt> points0 = new HashSet<PairInt>();

        populatePixelLists(input, points0, blackPixels, whitePixels, greyPixelMap);

        double[] minMaxTheta0 = findMinMaxTheta(input, points0);

        log.info("for all non-white and non-black, minTheta=" + minMaxTheta0[0]
            + " maxTheta=" + minMaxTheta0[1]);

        // -------- debug -------
        int nGrey = 0;
        for (Map.Entry<Integer, Collection<PairInt>> entry : greyPixelMap.entrySet()) {
            nGrey += entry.getValue().size();
        }
        int nGreyBW = + blackPixels.size() + nGrey + whitePixels.size();
        int nTot = (points0.size() + nGreyBW);
        log.info("img nPix=" + input.getNPixels() + " nTot=" + nTot);
        assert(nTot == input.getNPixels());
        // -------- end debug -------

        List<Set<PairInt>> greyPixelGroups = groupByPeaks(greyPixelMap);

        List<Set<PairInt>> groupList = new ArrayList<Set<PairInt>>(greyPixelGroups.size());

        if ((minMaxTheta0[1] - minMaxTheta0[0]) == 0) {
             if (points0.isEmpty()) {
                 groupList.add(points0);
             }
             if (!blackPixels.isEmpty()) {
                 groupList.add(blackPixels);
             }
             for (Set<PairInt> set : greyPixelGroups) {
                groupList.add(set);
             }
             if (!whitePixels.isEmpty()) {
                 groupList.add(whitePixels);
             }
             return groupList;
        }

        double xFactor = 2000.;
        int yFactor = 2000;

        double[] minMaxTheta = new double[2];
        int[] minMaxFreq = new int[2];
        double thetaFactor0 = xFactor/(minMaxTheta0[1] - minMaxTheta0[0]);
        Map<Integer, List<PairInt>> thetaPointMap = createThetaCIEXYMap(points0,
            input, minMaxTheta0[0], thetaFactor0, minMaxTheta, minMaxFreq);

        // ---- debug ------
        nTot = nGreyBW;
        for (Map.Entry<Integer, List<PairInt>> entry : thetaPointMap.entrySet()) {
            nTot += entry.getValue().size();
        }
        log.info("img nPix=" + input.getNPixels() + " nTot=" + nTot);
        assert(nTot == input.getNPixels());
        // ----- end debug ------

        /* ---- create frequency maps partitioned by given fractions ----
        starting w/ partitions at 3 percent (maybe discard below),
            15 percent, and 25 percent resulting in 4 maps

        For each map:
            key is pairint w/ x=theta, y=freq,
            value is all pixels having same key
        */

        final float[] partitionFreqFracs = new float[]{0.03f, 0.15f, 0.25f};

        List<Map<com.climbwithyourfeet.clustering.util.PairInt,
            List<PairIntWithIndex>>> thetaFreqMaps =
            partitionIntoFrequencyMaps(input, thetaPointMap,
                partitionFreqFracs, minMaxFreq[1]);

        //---- debug, assert number of pixels ----
        nTot = nGreyBW;
        for (Map<com.climbwithyourfeet.clustering.util.PairInt,
        List<PairIntWithIndex>> map : thetaFreqMaps) {
            for (Map.Entry<com.climbwithyourfeet.clustering.util.PairInt,
                List<PairIntWithIndex>> entry : map.entrySet()) {
                nTot += entry.getValue().size();
            }
        }
        log.info("img nPix=" + input.getNPixels() + " nTot=" + nTot);
        assert(nTot == input.getNPixels());

        // TODO: handle wrap around values!
        //    if there are points at 0 and 360, and a gap elsewhere, can
        //    shift the values so the gap is at 360 instead.

        // ------ TODO: rescale each map by frequencies to span ~1000 pixels -----
        rescaleKeys(thetaFreqMaps, yFactor);

        int nTot2 = 0;

        for (int i = 1; i < thetaFreqMaps.size(); ++i) {

            Map<com.climbwithyourfeet.clustering.util.PairInt,
                List<PairIntWithIndex>> thetaFreqMapI = thetaFreqMaps.get(i);

            int[] maxXY = findMaxXY(thetaFreqMapI.keySet());

            if (maxXY[0] < 0 || maxXY[1] <= 0) {
                continue;
            }

            // ----- debug ---
            // plot the points as an image to see the data first
            GreyscaleImage img = new GreyscaleImage(maxXY[0] + 1, maxXY[1] + 1);
            for (com.climbwithyourfeet.clustering.util.PairInt p : thetaFreqMapI.keySet()) {
                img.setValue(p.getX(), p.getY(), 255);
            }
            try {
                ImageIOHelper.writeOutputImage(
                    ResourceFinder.findDirectory("bin") + "/dt2_input_" + i
                        + "_.png", img);
            } catch (IOException ex) {
                Logger.getLogger(ImageProcessor.class.getName()).log(Level.SEVERE,
                    null, ex);
            }
            // --- end debug

            // Map<com.climbwithyourfeet.clustering.util.PairInt,
            //     List<PairIntWithIndex>> thetaFreqMapI

            // map w/ key=(theta, freq) value=collection of coords
            DTClusterFinder<com.climbwithyourfeet.clustering.util.PairInt>
                clusterFinder
                = new DTClusterFinder<com.climbwithyourfeet.clustering.util.PairInt>(
                    thetaFreqMapI.keySet(), maxXY[0] + 1, maxXY[1] + 1);

            clusterFinder.setToDebug();

            // to recover every point, set limit to 1
            clusterFinder.setMinimumNumberInCluster(1);

            clusterFinder.calculateCriticalDensity();

            clusterFinder.findClusters();

            int nGroups = clusterFinder.getNumberOfClusters();

            for (int k = 0; k < nGroups; ++k) {

                Set<com.climbwithyourfeet.clustering.util.PairInt> group
                    = clusterFinder.getCluster(k);

                Set<PairInt> coordPoints = new HashSet<PairInt>();

                for (com.climbwithyourfeet.clustering.util.PairInt pThetaFreq : group) {

                    // include the other points of same ciexy theta, freq
                    List<PairIntWithIndex> list = thetaFreqMapI.get(pThetaFreq);
                    assert (list != null);
                    for (PairIntWithIndex p3 : list) {
                        int idx3 = p3.getPixIndex();
                        int xCoord3 = input.getCol(idx3);
                        int yCoord3 = input.getRow(idx3);
                        PairInt pCoord = new PairInt(xCoord3, yCoord3);
                        coordPoints.add(pCoord);
                    }
                }

                nTot2 += coordPoints.size();

                groupList.add(coordPoints);
            }
        }

        mergeOrAppendGreyWithOthers(input, greyPixelGroups, groupList,
            blackPixels, whitePixels);

        // add back in blackPixels and whitePixels
        groupList.add(blackPixels);
        groupList.add(whitePixels);

        //TODO: assert npixels

        return groupList;
    }

    /**
     * Calculates lists of black pixels, white pixels, grey pixels, and assigns
     * the remaining to the polar angle of CIEXY Lab color space, then creates a map of the
     * polar angles and the number of pixels with those angles (=frequency),
     * then finds peaks in theta above a fraction =0.03 of max limit then groups all
     * pixels in the map by proximity to the peak,
     * then merges pixels in the grey list with adjacent clusters if
     * similar, and the final result is a list of pixel clusters, including the
     * black and white and remaining grey.  The list is sorted by size and
     * set to values from 255 to 0.  If there are more than 255 clusters, the
     * remaining (smaller) clusters are pixels with value 0.
     *
     * Note that the color black is not defined in CIE XY color space and that
     * the color white is at the center of the space as a large circle so they
     * are not included in the density based clustering.
     * Note also that the remaining grey pixels which did not merge with the
     * CIE xy pixel clusters, may be in separate groups already due to a
     * frequency based grouping for them.
     * @param input
     * @param useBlur apply a gaussian blur of sigma=1 before the method logic
     * @return
     */
    public GreyscaleImage applyUsingPolarCIEXYAndFrequency(ImageExt input,
        boolean useBlur) {

        float fracFreqLimit = 0.03f;

        return applyUsingPolarCIEXYAndFrequency(input, fracFreqLimit, useBlur);
    }

    /**
     * Calculates lists of black pixels, white pixels, grey pixels, and assigns
     * the remaining to the polar angle of CIEXY Lab color space, then creates a map of the
     * polar angles and the number of pixels with those angles (=frequency),
     * then finds peaks in theta above a fraction =0.03 of max limit then groups all
     * pixels in the map by proximity to the peak,
     * then merges pixels in the grey list with adjacent clusters if
     * similar, and the final result is a list of pixel clusters, including the
     * black and white and remaining grey.  The list is sorted by size and
     * set to values from 255 to 0.  If there are more than 255 clusters, the
     * remaining (smaller) clusters are pixels with value 0.
     *
     * Note that the color black is not defined in CIE XY color space and that
     * the color white is at the center of the space as a large circle so they
     * are not included in the density based clustering.
     * Note also that the remaining grey pixels which did not merge with the
     * CIE xy pixel clusters, may be in separate groups already due to a
     * frequency based grouping for them.
     * @param input
     * @param fracFreqLimit
     * @param useBlur apply a gaussian blur of sigma=1 before the method logic
     * @return
     */
    public GreyscaleImage applyUsingPolarCIEXYAndFrequency(ImageExt input,
        final float fracFreqLimit, boolean useBlur) {

        int w = input.getWidth();
        int h = input.getHeight();

        List<Set<PairInt>> clusters = calculateUsingPolarCIEXYAndFrequency(
            input, fracFreqLimit, useBlur);

        int n = clusters.size();
        //assert(n < 256);

        // sort indexes by set size
        int maxSize = Integer.MIN_VALUE;
        int[] indexes = new int[n];
        int[] sizes = new int[n];
        for (int i = 0; i < n; ++i) {
            indexes[i] = i;
            sizes[i] = clusters.get(i).size();
            if (sizes[i] > maxSize) {
                maxSize = sizes[i];
            }
        }
        CountingSort.sort(sizes, indexes, maxSize);

        int delta = 256/clusters.size();
        if (delta == 0) {
            delta = 1;
        }

        GreyscaleImage img2 = new GreyscaleImage(w, h);
        for (int k = 0; k < n; ++k) {

            int idx = indexes[n - k - 1];

            Set<PairInt> set = clusters.get(idx);

            int v = 255 - delta*k;
            if (v < 1) {
                continue;
            }

            for (PairInt p : set) {
                img2.setValue(p.getX(), p.getY(), v);
            }
        }

        return img2;
    }

    /**
     * Calculates lists of black pixels, white pixels, grey pixels, and assigns
     * the remaining to the polar angle of CIEXY Lab color space, then creates a map of the
     * polar angles and the number of pixels with those angles (=frequency),
     * then finds peaks in theta above a fraction =0.03 of max limit then groups all
     * pixels in the map by proximity to the peak,
     * then merges pixels in the grey list with adjacent clusters if
     * similar, and the final result is a list of pixel clusters, including the
     * black and white and remaining grey.
     *
     * Note that the color black is not defined in CIE XY color space and that
     * the color white is at the center of the space as a large circle so they
     * are not included in the density based clustering.
     * Note also that the remaining grey pixels which did not merge with the
     * CIE xy pixel clusters, may be in separate groups already due to a
     * frequency based grouping for them.
     * @param input
     * @param useBlur apply a gaussian blur of sigma=1 before the method logic
     * @return
     */
    public List<Set<PairInt>> calculateUsingPolarCIEXYAndFrequency(ImageExt input,
        boolean useBlur) {

        float fracFreqLimit = 0.03f;

        return calculateUsingPolarCIEXYAndFrequency(input, fracFreqLimit, useBlur);
    }

    /**
     * Calculates lists of black pixels, white pixels, grey pixels, and assigns
     * the remaining to the polar angle of CIEXY Lab color space, then creates a map of the
     * polar angles and the number of pixels with those angles (=frequency),
     * then finds peaks in theta above a fraction of max limit then groups all
     * pixels in the map by proximity to the peak,
     * then merges pixels in the grey list with adjacent clusters if
     * similar, and the final result is a list of pixel clusters, including the
     * black and white and remaining grey.
     * The changeable parameter is the fracFreqLimit.  Larger values exclude
     * smaller frequency peaks.
     *
     * Note that the color black is not defined in CIE XY color space and that
     * the color white is at the center of the space as a large circle so they
     * are not included in the density based clustering.
     * Note also that the remaining grey pixels which did not merge with the
     * CIE xy pixel clusters, may be in separate groups already due to a
     * frequency based grouping for them.
     *
     * @param input image to find color clusters within
     * @param fracFreqLimit fraction of the maximum above which peaks will be found
     * @param useBlur if true, performs a gaussian blur of sigma=1 before finding
     * clusters.
     * @return
     */
    public List<Set<PairInt>> calculateUsingPolarCIEXYAndFrequency(ImageExt input,
        float fracFreqLimit, boolean useBlur) {

        if (useBlur) {
            input = input.copyToImageExt();
            ImageProcessor imageProcessor = new ImageProcessor();
            imageProcessor.blur(input, 1.0f);
        }

        //making a segmentation method using CIEXY polar theta
        // and the number of points with those colors.
        // choosing the peaks to be the cluster centers, then
        // gathering the pixels by proximity to the theta peaks
        // and when equidistant, chooses the largest peak.

        Set<PairInt> blackPixels = new HashSet<PairInt>();

        Map<Integer, Collection<PairInt>> greyPixelMap = new HashMap<Integer, Collection<PairInt>>();

        Set<PairInt> whitePixels = new HashSet<PairInt>();

        Set<PairInt> points0 = new HashSet<PairInt>();

        populatePixelLists(input, points0, blackPixels, whitePixels, greyPixelMap);

        double[] minMaxTheta0 = findMinMaxTheta(input, points0);

        log.info("for all non-white and non-black, minTheta=" + minMaxTheta0[0]
            + " maxTheta=" + minMaxTheta0[1]);

        // -------- debug -------
        int nGrey = 0;
        for (Map.Entry<Integer, Collection<PairInt>> entry : greyPixelMap.entrySet()) {
            nGrey += entry.getValue().size();
        }
        assert((points0.size() + blackPixels.size() + nGrey +
            whitePixels.size()) == input.getNPixels());
        // -------- end debug -------

        List<Set<PairInt>> greyPixelGroups = groupByPeaks(greyPixelMap);

        List<Set<PairInt>> groupList = new ArrayList<Set<PairInt>>(greyPixelGroups.size());

        if ((minMaxTheta0[1] - minMaxTheta0[0]) == 0) {
             if (points0.isEmpty()) {
                 groupList.add(points0);
             }
             if (!blackPixels.isEmpty()) {
                 groupList.add(blackPixels);
             }
             for (Set<PairInt> set : greyPixelGroups) {
                groupList.add(set);
             }
             if (!whitePixels.isEmpty()) {
                 groupList.add(whitePixels);
             }
             return groupList;
        }

        /* ----- create a map of theta and frequency ----
        need to find the peaks in frequency for frequencies larger than about
        3 percent of max frequency
        but don't want to use a spline3 to smooth, so will average every
        few pixels.
        */

        int binWidth = 3;
        Map<Integer, Collection<PairInt>> thetaPointMap = createThetaCIEXYMap(
            points0, input, binWidth);

        int n = (360/binWidth) + 1;

        int[] orderedThetaKeys = new int[n];
        for (int i = 0; i < n; ++i) {
            orderedThetaKeys[i] = i;
        }
        int maxFreq = Integer.MIN_VALUE;
        int nTot = 0;
        for (Map.Entry<Integer, Collection<PairInt>> entry : thetaPointMap.entrySet()) {
            int count = entry.getValue().size();
            if (count > maxFreq) {
                maxFreq = count;
            }
            nTot += entry.getValue().size();
        }
        nTot += (blackPixels.size() + nGrey + whitePixels.size());
        assert(nTot == input.getNPixels());

        /*
        TODO: this is where the DTClusterFinder would be good to use to find
        the peaks.
        */

        PairIntArray peaks = findPeaksInThetaPointMap(orderedThetaKeys,
            thetaPointMap,
            Math.round(fracFreqLimit * maxFreq));

        /*
        // ----- debug ---
        // plot the points as an image to see the data first
        int[] minMaxXY = MiscMath.findMinMaxXY(peaks);
        int nPoints = 0;
        int maxX = Integer.MIN_VALUE;
        int maxY = Integer.MIN_VALUE;
        for (int i : orderedThetaKeys) {
            Integer key = Integer.valueOf(i);
            Collection<PairInt> list = thetaPointMap.get(key);
            if (list == null) {
                continue;
            }
            int y = list.size();
            nPoints++;
            if (key.intValue() > maxX) {
                maxX = key.intValue();
            }
            if (y > maxY) {
                maxY = y;
            }
        }
        float[] xPoints = new float[nPoints];
        float[] yPoints = new float[nPoints];
        int count = 0;
        for (int i : orderedThetaKeys) {
            Integer key = Integer.valueOf(i);
            Collection<PairInt> list = thetaPointMap.get(key);
            if (list == null) {
                continue;
            }
            int y = list.size();
            xPoints[count] = key.intValue();
            yPoints[count] = y;
            count++;
        }
        try {
            //maxY=2000;
            PolygonAndPointPlotter plotter = new PolygonAndPointPlotter();
            plotter.addPlot(0, maxX, 0, maxY, xPoints, yPoints, xPoints, yPoints, "cieXY theta vs freq");
            plotter.writeFile("_segmentation3_");
        } catch (IOException ex) {
            Logger.getLogger(ImageProcessor.class.getName()).log(Level.SEVERE,
                null, ex);
        }
        // --- end debug
        */

        if (peaks.getN() == 0) {
            for (Map.Entry<Integer, Collection<PairInt>> entry : thetaPointMap.entrySet()) {
                groupList.add(new HashSet<PairInt>(entry.getValue()));
            }
            groupList.add(blackPixels);
            groupList.add(whitePixels);
            return groupList;
        }
        for (int i = 0; i < peaks.getN(); ++i) {
            groupList.add(new HashSet<PairInt>());
        }

        /* traverse in ordered manner thetaPointMap to compare to current theta position
           w.r.t. peaks
           then place it in the groupsList.
           points before the first peak are compared with last peak too for wrap around.
        */
        int currentPeakIdx = -1;
        for (int i : orderedThetaKeys) {
            Integer key = Integer.valueOf(i);
            Collection<PairInt> list = thetaPointMap.get(key);
            if (list == null) {
                continue;
            }
            int idx = -1;
            if ((currentPeakIdx == -1) || (currentPeakIdx == (peaks.getN() - 1))) {
                int diffL, diffF;
                if (currentPeakIdx == -1) {
                    diffL = key.intValue() + 360 - peaks.getX(peaks.getN() - 1);
                    diffF = peaks.getX(0) - key.intValue();
                    if (diffF == 0) {
                        currentPeakIdx = 0;
                    }
                } else {
                    diffL = key.intValue() - peaks.getX(currentPeakIdx);
                    diffF = peaks.getX(0) + 360 - key.intValue();
                }
                if (diffL < diffF) {
                    idx = peaks.getN() - 1;
                } else if (diffL == diffF) {
                    int freqL = peaks.getY(peaks.getN() - 1);
                    int freqF = peaks.getY(0);
                    if (freqL < freqF) {
                        idx = peaks.getN() - 1;
                    } else {
                        idx = 0;
                    }
                } else {
                    idx = 0;
                }
            } else {
                // this has to update currentPeakIdx
                int diffP = key.intValue() - peaks.getX(currentPeakIdx);
                int diffN = peaks.getX(currentPeakIdx + 1) - key.intValue();
                if (diffN == 0) {
                    currentPeakIdx++;
                    idx = currentPeakIdx;
                } else {
                    if (diffP < diffN) {
                        idx = currentPeakIdx;
                    } else if (diffP == diffN) {
                        int freqP = peaks.getY(currentPeakIdx);
                        int freqN = peaks.getY(currentPeakIdx + 1);
                        if (freqP < freqN) {
                            idx = currentPeakIdx;
                        } else {
                            idx = currentPeakIdx + 1;
                        }
                    } else {
                        idx = currentPeakIdx + 1;
                    }
                }
            }
            assert(idx != -1);
            groupList.get(idx).addAll(list);
        }

        mergeOrAppendGreyWithOthers(input, greyPixelGroups, groupList,
            blackPixels, whitePixels);

        // add back in blackPixels and whitePixels
        /*if (!blackPixels.isEmpty()) {
            groupList.add(blackPixels);
        }
        if (!whitePixels.isEmpty()) {
            groupList.add(whitePixels);
        }

        int nTot2 = 0;
        for (Set<PairInt> groups : groupList) {
            nTot2 += groups.size();
        }

        log.info("img nPix=" + input.getNPixels() + " nTot2=" + nTot2);
        assert(nTot2 == input.getNPixels());
        */
        return groupList;
    }

    private Map<Integer, List<PairInt>> createThetaCIEXYMap(Set<PairInt>
        points0, ImageExt input, double minTheta, double thetaFactor,
        double[] outputMinMaxTheta, int[] outputMinMaxFreq) {

        CIEChromaticity cieC = new CIEChromaticity();

        // key = theta, value = pixels having that key
        Map<Integer, List<PairInt>> thetaPointMap = new HashMap<Integer, List<PairInt>>();

        double minTheta0 = minTheta;
        outputMinMaxTheta[0] = Double.MIN_VALUE;
        outputMinMaxTheta[1] = Double.MAX_VALUE;

        for (PairInt p : points0) {

            int idx = input.getInternalIndex(p.getX(), p.getY());

            float cx = input.getCIEX(idx);
            float cy = input.getCIEY(idx);

            double theta = thetaFactor * (
                (cieC.calculateXYTheta(cx, cy)*180./Math.PI) - minTheta0);

            Integer thetaCIEXY = Integer.valueOf((int)Math.round(theta));

            List<PairInt> list = thetaPointMap.get(thetaCIEXY);
            if (list == null) {
                list = new ArrayList<PairInt>();
                thetaPointMap.put(thetaCIEXY, list);
            }
            list.add(p);

            if (theta < outputMinMaxTheta[0]) {
                outputMinMaxTheta[0] = theta;
            }
            if (theta > outputMinMaxTheta[1]) {
                outputMinMaxTheta[1] = theta;
            }
        }

        outputMinMaxFreq[0] = Integer.MAX_VALUE;
        outputMinMaxFreq[1] = Integer.MIN_VALUE;
        for (Map.Entry<Integer, List<PairInt>> entry : thetaPointMap.entrySet()) {
            int count = entry.getValue().size();
            if (count < outputMinMaxFreq[0]) {
                outputMinMaxFreq[0] = count;
            }
            if (count > outputMinMaxFreq[1]) {
                outputMinMaxFreq[1] = count;
            }
        }

        return thetaPointMap;
    }

    private List<Map<com.climbwithyourfeet.clustering.util.PairInt,
    List<PairIntWithIndex>>> partitionIntoFrequencyMaps(
    ImageExt input, Map<Integer, List<PairInt>> thetaPointMap,
    float[] partitionFreqFracs, int maxFreq) {

        List<Map<com.climbwithyourfeet.clustering.util.PairInt,
            List<PairIntWithIndex>>>
            mapsList =
                new ArrayList<Map<com.climbwithyourfeet.clustering.util.PairInt,
                List<PairIntWithIndex>>>();

        int nMaps = partitionFreqFracs.length + 1;

        for (int i = 0; i < nMaps; ++i) {
            mapsList.add(
                new HashMap<com.climbwithyourfeet.clustering.util.PairInt,
                List<PairIntWithIndex>>());
        }

        final int[] partitionFreqs = new int[partitionFreqFracs.length];
        for (int i = 0; i < partitionFreqs.length; ++i) {
            partitionFreqs[i] = Math.round(partitionFreqFracs[i]*maxFreq);
        }

        int[] maxXY = new int[2];
        Arrays.fill(maxXY, Integer.MIN_VALUE);

        for (Map.Entry<Integer, List<PairInt>> entry : thetaPointMap.entrySet()) {

            Integer theta = entry.getKey();

            int count = entry.getValue().size();

            List<PairIntWithIndex> list = new ArrayList<PairIntWithIndex>();
            for (PairInt p : entry.getValue()) {
                int pixIdx = input.getInternalIndex(p.getX(), p.getY());
                PairIntWithIndex p2 = new PairIntWithIndex(theta.intValue(),
                    count, pixIdx);
                list.add(p2);
            }

            int n = partitionFreqs.length;
            int idx = 0;
            for (int i = 0; i < n; ++i) {
                if (i == 0) {
                    if (count < partitionFreqs[0]) {
                        idx = 0;
                        break;
                    }
                } else if ((i == (n - 1)) && count >= partitionFreqs[i]) {
                    idx = n;// one past partitions is last list bin
                    break;
                } else {
                    if ((count >= partitionFreqs[i - 1]) && (count < partitionFreqs[i])) {
                        idx = i;
                        break;
                    }
                }
            }

            // this is unique to all maps, not stomping on existing key
            mapsList.get(idx).put(
                new com.climbwithyourfeet.clustering.util.PairInt(
                    theta.intValue(), count), list);

            if (theta.intValue() > maxXY[0]) {
                maxXY[0] = theta.intValue();
            }
            if (count > maxXY[1]) {
                maxXY[1] = count;
            }
        }

        // ----- temporary print of all pixels -------
        // ----- debug ---
        GreyscaleImage img = new GreyscaleImage(maxXY[0] + 1, maxXY[1] + 1);
        for (int i = 0; i < mapsList.size(); ++i) {
            for (com.climbwithyourfeet.clustering.util.PairInt p : mapsList.get(i).keySet()) {
                img.setValue(p.getX(), p.getY(), 255);
            }
        }
        try {
            ImageIOHelper.writeOutputImage(
                ResourceFinder.findDirectory("bin") + "/dt2_input.png", img);
        } catch (IOException ex) {
            Logger.getLogger(ImageProcessor.class.getName()).log(Level.SEVERE,
                null, ex);
        }
        // --- end debug

        return mapsList;
    }

    private void rescaleKeys(
        List<Map<com.climbwithyourfeet.clustering.util.PairInt,
            List<PairIntWithIndex>>> thetaFreqMaps, int scaleTo) {

        // --- can remove the count after debugging ----
        int nTotBefore = 0;
        for (int i = 0; i < thetaFreqMaps.size(); ++i) {
            for (Map.Entry<com.climbwithyourfeet.clustering.util.PairInt,
                List<PairIntWithIndex>> entry : thetaFreqMaps.get(i).entrySet()) {
                nTotBefore += entry.getValue().size();
            }
        }

        for (int i = 0; i < thetaFreqMaps.size(); ++i) {

            Map<com.climbwithyourfeet.clustering.util.PairInt,
                List<PairIntWithIndex>> thetaFreqMap = thetaFreqMaps.get(i);

            int[] minMax = findMinMaxOfKeyYs(thetaFreqMap.keySet());

            Map<com.climbwithyourfeet.clustering.util.PairInt,
                List<PairIntWithIndex>> thetaFreqMap2 = rescaleKeyYs(
                    thetaFreqMap, scaleTo, minMax);

            thetaFreqMaps.set(i, thetaFreqMap2);
        }

        // --- can remove the count after debugging ----
        int nTotAfter = 0;
        for (int i = 0; i < thetaFreqMaps.size(); ++i) {
            for (Map.Entry<com.climbwithyourfeet.clustering.util.PairInt,
                List<PairIntWithIndex>> entry : thetaFreqMaps.get(i).entrySet()) {
                nTotAfter += entry.getValue().size();
            }
        }
        assert(nTotBefore == nTotAfter);
    }

    private int[] findMinMaxOfKeyYs(
        Set<com.climbwithyourfeet.clustering.util.PairInt> keySet) {

        int min = Integer.MAX_VALUE;
        int max = Integer.MIN_VALUE;

        for (com.climbwithyourfeet.clustering.util.PairInt p : keySet) {
            int y = p.getY();
            if (y < min) {
                min = y;
            }
            if (y > max) {
                max = y;
            }
        }
        return new int[]{min, max};
    }

    private Map<com.climbwithyourfeet.clustering.util.PairInt,
    List<PairIntWithIndex>> rescaleKeyYs(
    Map<com.climbwithyourfeet.clustering.util.PairInt,
    List<PairIntWithIndex>> thetaFreqMap, final int scaleTo, final int[] minMaxY) {

        Map<com.climbwithyourfeet.clustering.util.PairInt,
            List<PairIntWithIndex>> scaledMap
            = new HashMap<com.climbwithyourfeet.clustering.util.PairInt,
                List<PairIntWithIndex>>();

        if ((minMaxY[1] - minMaxY[0]) == 0) {
            return thetaFreqMap;
        }

        float factor = scaleTo/(minMaxY[1] - minMaxY[0]);

        for (Map.Entry<com.climbwithyourfeet.clustering.util.PairInt,
            List<PairIntWithIndex>> entry : thetaFreqMap.entrySet()) {

            com.climbwithyourfeet.clustering.util.PairInt p = entry.getKey();

            int y = Math.round(factor * (p.getY() - minMaxY[0]));

            com.climbwithyourfeet.clustering.util.PairInt p2 = new
                com.climbwithyourfeet.clustering.util.PairInt(p.getX(), y);

            scaledMap.put(p2, entry.getValue());
        }

        return scaledMap;
    }

    private int[] findMaxXY(Set<com.climbwithyourfeet.clustering.util.PairInt>
        keySet) {

        int maxX = Integer.MIN_VALUE;
        int maxY = Integer.MIN_VALUE;

        for (com.climbwithyourfeet.clustering.util.PairInt p : keySet) {
            int x = p.getX();
            int y = p.getY();
            if (x > maxX) {
                maxX = x;
            }
            if (y > maxY) {
                maxY = y;
            }
        }
        return new int[]{maxX, maxY};
    }

    private Map<Integer, Collection<PairInt>> createThetaCIEXYMap(Set<PairInt>
        points, ImageExt input, int binWidth) {

        CIEChromaticity cieC = new CIEChromaticity();

        // key = theta, value = pixels having that key
        Map<Integer, Collection<PairInt>> thetaPointMap = new HashMap<Integer, Collection<PairInt>>();

        for (PairInt p : points) {

            int idx = input.getInternalIndex(p.getX(), p.getY());

            float cx = input.getCIEX(idx);
            float cy = input.getCIEY(idx);

            double thetaRadians = cieC.calculateXYTheta(cx, cy);
            double theta = thetaRadians * 180./Math.PI;

            int thetaCIEXY = (int)Math.round(theta);

            Integer binKey = Integer.valueOf(thetaCIEXY/binWidth);

            Collection<PairInt> list = thetaPointMap.get(binKey);
            if (list == null) {
                list = new ArrayList<PairInt>();
                thetaPointMap.put(binKey, list);
            }
            list.add(p);
        }

        return thetaPointMap;
    }

    /**
     * find peaks in the theta point map above lower limit.
     * @param orderedThetaKeys
     * @param thetaPointMap
     * @param limit
     * @return
     */
    protected PairIntArray findPeaksInThetaPointMap(final int[] orderedThetaKeys,
        final Map<Integer, Collection<PairInt>> thetaPointMap, final int limit) {

        int lastKey = -1;
        int lastValue = -1;
        boolean isIncr = false;
        PairIntArray peaks = new PairIntArray();
        int nInMap = 0;
        for (int i : orderedThetaKeys) {
            Integer key = Integer.valueOf(i);
            Collection<PairInt> list = thetaPointMap.get(key);
            if (list == null) {
                if ((nInMap > 0) && isIncr && (lastValue > limit)) {
                    peaks.add(lastKey, lastValue);
                }
                lastKey = key.intValue();
                lastValue = 0;
                isIncr = false;
                continue;
            }
            int count = list.size();
            if (nInMap == 1) {
                if (count > lastValue) {
                    isIncr = true;
                } else {
                    if (lastValue > limit) {
                        peaks.add(lastKey, lastValue);
                    }
                    isIncr = false;
                }
            } else if (nInMap != 0) {
                if (isIncr) {
                    if (count < lastValue) {
                        if (lastValue > limit) {
                            peaks.add(lastKey, lastValue);
                        }
                        isIncr = false;
                    }
                } else {
                    if (count > lastValue) {
                        isIncr = true;
                    }
                }
            }
            lastValue = count;
            lastKey = key.intValue();
            nInMap++;
        }
        if ((nInMap > 0) && isIncr && (lastValue > limit)) {
            //checking value at theta=0 to make sure this is a peak
            Integer key = orderedThetaKeys[0];
            if (key.intValue() == 0) {
                Collection<PairInt> list = thetaPointMap.get(key);
                if (list == null) {
                    peaks.add(lastKey, lastValue);
                } else if (list.size() < lastValue) {
                    peaks.add(lastKey, lastValue);
                }
            } else {
                peaks.add(lastKey, lastValue);
            }
        }

        return peaks;
    }

    private List<Set<PairInt>> groupByPeaks(
        Map<Integer, Collection<PairInt>> greyPixelMap) {

        int nTot = 0;
        int minKey = Integer.MAX_VALUE;
        int maxKey = Integer.MIN_VALUE;
        for (Map.Entry<Integer, Collection<PairInt>> entry : greyPixelMap.entrySet()) {
            int key = entry.getKey().intValue();
            if (key < minKey) {
                minKey = key;
            }
            if (key > maxKey) {
                maxKey = key;
            }
            nTot += entry.getValue().size();
        }

        int binWidth = 8;
        greyPixelMap = binByKeys(greyPixelMap, minKey, maxKey, binWidth);

        int nTot2 = 0;
        minKey = Integer.MAX_VALUE;
        maxKey = Integer.MIN_VALUE;
        int maxFreq = Integer.MIN_VALUE;
        for (Map.Entry<Integer, Collection<PairInt>> entry : greyPixelMap.entrySet()) {
            int key = entry.getKey().intValue();
            if (key < minKey) {
                minKey = key;
            }
            if (key > maxKey) {
                maxKey = key;
            }
            int y = entry.getValue().size();
            if (y > maxFreq) {
                maxFreq = y;
            }

            nTot2 += y;
        }
        assert(nTot == nTot2);

        int count = 0;

        /*
        // --- debug
        float[] xPoints = new float[greyPixelMap.size()];
        float[] yPoints = new float[greyPixelMap.size()];
        for (int i = minKey; i <= maxKey; ++i) {
            Integer key = Integer.valueOf(i);
            Collection<PairInt> set = greyPixelMap.get(key);
            if (set == null) {
                continue;
            }
            int y = set.size();
            xPoints[count] = key.intValue();
            yPoints[count] = y;
            count++;
        }
        try {
            //maxY=2000;
            PolygonAndPointPlotter plotter = new PolygonAndPointPlotter();
            plotter.addPlot(0, maxKey, 0, maxFreq, xPoints, yPoints, xPoints, yPoints, "grey avgRGB vs freq");
            plotter.writeFile("_segmentation3_grey_");
        } catch (IOException ex) {
            Logger.getLogger(ImageProcessor.class.getName()).log(Level.SEVERE,
                null, ex);
        }
        // --- end debug
        */

        final int[] orderedKeys = new int[maxKey - minKey + 1];
        count = 0;
        for (int i = minKey; i <= maxKey; ++i) {
            orderedKeys[count] = i;
            count++;
        }
        int limit = (int)(0.03 * maxFreq);
        PairIntArray peaks = findPeaksInThetaPointMap(orderedKeys, greyPixelMap, limit);

        // if there are several peaks within small range of keys, that's noise,
        // so removing them
        filterPeaksIfNoisey(peaks);

        // ---- gather points in greyPixelMap into groups around the peaks ----
        List<Set<PairInt>> groupList = new ArrayList<Set<PairInt>>(orderedKeys.length + 1);
        for (int i = 0; i < peaks.getN(); ++i) {
            groupList.add(new HashSet<PairInt>());
        }

        if (peaks.getN() == 0) {
            return groupList;
        } else if (peaks.getN() == 1) {
            for (Map.Entry<Integer, Collection<PairInt>> entry : greyPixelMap.entrySet()) {
                Collection<PairInt> set = entry.getValue();
                groupList.get(0).addAll(set);
            }
            return groupList;
        }

        int currentPeakIdx = -1;
        for (int i : orderedKeys) {
            Integer key = Integer.valueOf(i);
            Collection<PairInt> set = greyPixelMap.get(key);
            if (set == null) {
                continue;
            }
            int idx = -1;
            if (currentPeakIdx == -1) {
                currentPeakIdx = 0;
                idx = 0;
            } else if (currentPeakIdx == (peaks.getN() - 1)) {
                idx = currentPeakIdx;
            } else {
                // this has to update currentPeakIdx
                int diffP = key.intValue() - peaks.getX(currentPeakIdx);
                int diffN = peaks.getX(currentPeakIdx + 1) - key.intValue();
                if (diffN == 0) {
                    currentPeakIdx++;
                    idx = currentPeakIdx;
                } else {
                    if (diffP < diffN) {
                        idx = currentPeakIdx;
                    } else if (diffP == diffN) {
                        int freqP = peaks.getY(currentPeakIdx);
                        int freqN = peaks.getY(currentPeakIdx + 1);
                        if (freqP < freqN) {
                            idx = currentPeakIdx;
                        } else {
                            idx = currentPeakIdx + 1;
                        }
                    } else {
                        idx = currentPeakIdx + 1;
                    }
                }
            }
            assert(idx != -1);
            groupList.get(idx).addAll(set);
        }

        // ----- debug ----
        int nTot3 = 0;
        for (Collection<PairInt> set : groupList) {
            nTot3 += set.size();
        }
        assert(nTot == nTot3);
        // ----- end debug -----

        return groupList;
    }

    private Map<Integer, Collection<PairInt>> binByKeys(
        Map<Integer, Collection<PairInt>> greyPixelMap,
        int minKey, int maxKey, int binWidth) {

        Map<Integer, Collection<PairInt>> map2
            = new HashMap<Integer, Collection<PairInt>>();

        for (int i = minKey; i <= maxKey; ++i) {

            Integer key = Integer.valueOf(i);

            Collection<PairInt> c = greyPixelMap.get(key);

            if (c == null) {
                continue;
            }

            Integer binKey = Integer.valueOf(i/binWidth);

            Collection<PairInt> c2 = map2.get(binKey);
            if (c2 == null) {
                c2 = new HashSet<PairInt>();
                map2.put(binKey, c2);
            }
            c2.addAll(c);
        }

        return map2;
    }

    private void filterPeaksIfNoisey(PairIntArray peaks) {

        if (peaks.getN() == 0) {
            return;
        }

        int sumDeltaX = 0;
        for (int i = (peaks.getN() - 1); i > 0; --i) {
            sumDeltaX += (peaks.getX(i) - peaks.getX(i - 1));
        }
        //TODO: this may need to be revised:
        // if there are more than 1 peaks per delta x of 5 or so, re-bin by 4
        float deltaX = (float)sumDeltaX/((float)peaks.getN() - 1);
        if (deltaX < 5) {
            // re-bin by 4
            PairIntArray peaks2 = new PairIntArray();
            for (int i = 0; i < peaks.getN(); i += 4) {
                int sumX = 0;
                int sumY = 0;
                int count = 0;
                for (int j = i; j < (i + 5); ++j) {
                    sumX += peaks.getX(i);
                    sumY += peaks.getY(i);
                    count++;
                }
                sumX = Math.round((float)sumX/(float)count);
                sumY = Math.round((float)sumY/(float)count);
                peaks2.add(sumX, sumY);
            }

            peaks.removeRange(0, peaks.getN() - 1);
            peaks.addAll(peaks2);
        }

    }

    private void mergeOrAppendGreyWithOthers(ImageExt input,
        List<Set<PairInt>> greyPixelGroups, List<Set<PairInt>> groupList,
        Set<PairInt> blackPixels, Set<PairInt> whitePixels) {

        Map<PairInt, Integer> colorPixGroupMap = new HashMap<PairInt, Integer>();
        for (int i = 0; i < groupList.size(); ++i) {
            Set<PairInt> set = groupList.get(i);
            Integer groupKey = Integer.valueOf(i);
            for (PairInt p : set) {
                colorPixGroupMap.put(p, groupKey);
            }
        }

        // similarity limit for a grey pixel to join adjacent color pixel's cluster
        int limit = 40;

        int w = input.getWidth();
        int h = input.getHeight();

        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;
        for (Set<PairInt> greyGroup : greyPixelGroups) {

            Set<PairInt> remove = new HashSet<PairInt>();

            for (PairInt greyP : greyGroup) {
                int x = greyP.getX();
                int y = greyP.getY();

                int idx = input.getInternalIndex(x, y);
                int r = input.getR(idx);
                int g = input.getG(idx);
                int b = input.getB(idx);

                // ---- check for color similarity ------
                int minDiffRGB = Integer.MAX_VALUE;
                int colorClusterIdx = -1;
                int minDiffBlack = Integer.MAX_VALUE;
                int minDiffWhite = Integer.MAX_VALUE;

                for (int i = 0; i < dxs.length; ++i) {
                    int x2 = x + dxs[i];
                    int y2 = y + dys[i];
                    if ((x2 < 0) || (y2 < 0) || (x2 > (w - 1)) || (y2 > (h - 1))) {
                        continue;
                    }
                    PairInt p2 = new PairInt(x2, y2);

                    boolean adjIsBlack = blackPixels.contains(p2);
                    boolean adjIsWhite = whitePixels.contains(p2);

                    Integer colorClusterIndex = colorPixGroupMap.get(p2);
                    if ((colorClusterIndex == null) && !adjIsBlack && !adjIsWhite) {
                        continue;
                    }

                    int idx2 = input.getInternalIndex(x2, y2);
                    int r2 = input.getR(idx2);
                    int g2 = input.getG(idx2);
                    int b2 = input.getB(idx2);

                    int diffR = Math.abs(r2 - r);
                    int diffG = Math.abs(g2 - g);
                    int diffB = Math.abs(b2 - b);

                    int diffRGB = diffR + diffG + diffB;

                    if (adjIsBlack) {
                        if (diffR == diffG && diffR == diffB) {
                            minDiffBlack = diffR;
                        }
                    } else if (adjIsWhite) {
                        if (diffR == diffG && diffR == diffB) {
                            minDiffWhite = diffR;
                        }
                    } else {
                        if ((diffRGB < minDiffRGB) && (diffRGB < limit)) {
                            minDiffRGB = diffRGB;
                            colorClusterIdx = (colorClusterIndex == null) ? -1 :
                                colorClusterIndex.intValue();
                        }
                    }
                }
                if (minDiffBlack < 75) {
                    blackPixels.add(greyP);
                    remove.add(greyP);
                } else if (minDiffWhite < 75) {
                    whitePixels.add(greyP);
                    remove.add(greyP);
                } else {
                    if (colorClusterIdx != -1) {
                        //add to color cluster and remove from grey list
                        groupList.get(colorClusterIdx).add(greyP);
                        remove.add(greyP);
                    }
                }
            }
            for (PairInt rm : remove) {
                greyGroup.remove(rm);
            }
        }

        // any remaining points in the grey list should be added as sets to
        // the color pixels list now
        for (int i = 0; i < greyPixelGroups.size(); ++i) {
            Set<PairInt> greyGroup = greyPixelGroups.get(i);
            groupList.add(greyGroup);
        }
        greyPixelGroups.clear();
    }

    protected double[] findMinMaxTheta(ImageExt input, Set<PairInt> points0) {

        CIEChromaticity cieC = new CIEChromaticity();

        double[] minMaxTheta = new double[2];

        for (PairInt p : points0) {
            int x = p.getX();
            int y = p.getY();

            double thetaRadians = cieC.calculateXYTheta(x, y);
            double theta = thetaRadians * 180. / Math.PI;

            if (theta < minMaxTheta[0]) {
                minMaxTheta[0] = theta;
            }
            if (theta > minMaxTheta[1]) {
                minMaxTheta[1] = theta;
            }
        }

        return minMaxTheta;
    }

    /**
     * Find the 4 categories of point color as black, white, grey, and color.
     * CIE XY color space is used to place a pixel in the category, with the
     * additional determination of black, grey, and white intensity limits
     * using histograms within expected intensity limits.
     * 
     * @param input
     * @param points0
     * @param blackPixels
     * @param whitePixels
     * @param greyPixelMap 
     */
    private void populatePixelLists(ImageExt input, Set<PairInt> points0,
        Set<PairInt> blackPixels, Set<PairInt> whitePixels,
        Map<Integer, Collection<PairInt>> greyPixelMap) {

        int w = input.getWidth();
        int h = input.getHeight();

        CIEChromaticity cieC = new CIEChromaticity();

        // looking for limits in peaks of (r,g,b) <= (45,45,45) and > (191,191,191)
        int[] whiteBlackLimits = findByHistogramLimitsForBlackAndWhite(input);
        // overriding:
        int whiteLimit = 245;
        
        for (int i = 0; i < w; ++i) {
            for (int j = 0; j < h; ++j) {

                int idx = input.getInternalIndex(i, j);

                int r = input.getR(idx);
                int g = input.getG(idx);
                int b = input.getB(idx);

                int avg = (r + g + b)/3;
                if (avg <= whiteBlackLimits[0]) {
                    blackPixels.add(new PairInt(i, j));
                    continue;
                }

                float rgbTot = r + g + b;
                float bDivTot = (float)b/rgbTot;
                float rDivTot = (float)r/rgbTot;
                float gDivTot = (float)g/rgbTot;
                        
                float cx = input.getCIEX(idx);
                float cy = input.getCIEY(idx);

                if (cieC.isWhite2(cx, cy) && 
                    (Math.abs(0.333f - bDivTot) < 0.02f) && 
                    (Math.abs(0.333f - rDivTot) < 0.02f) &&
                    (Math.abs(0.333f - gDivTot) < 0.02f)) {
                    
                    //if (avg > whiteBlackLimits[1]) {
                    if ((r > whiteLimit) && (g > whiteLimit) && (b > whiteLimit)) {  
                        whitePixels.add(new PairInt(i, j));
                    } else {
                        Integer avgRGB = Integer.valueOf(avg);
                        Collection<PairInt> set = greyPixelMap.get(avgRGB);
                        if (set == null) {
                            set = new HashSet<PairInt>();
                            greyPixelMap.put(avgRGB, set);
                        }
                        set.add(new PairInt(i, j));
                    }
                } else {
                    points0.add(new PairInt(i, j));
                }
            }
        }
    }
    
    private void mergeNoise(ImageExt input, List<Set<PairInt>> groupList) {

        Map<PairInt, Integer> pixelToGroupMap = new HashMap<PairInt, Integer>();
        for (int i = 0; i < groupList.size(); ++i) {
            Set<PairInt> set = groupList.get(i);
            for (PairInt p : set) {
                pixelToGroupMap.put(p, Integer.valueOf(i));
            }
        }

        final int w = input.getWidth();
        final int h = input.getHeight();

        int[] dxs8 = Misc.dx8;
        int[] dys8 = Misc.dy8;

        float diffLimit = 0.01f;

        for (int i = 0; i < groupList.size(); ++i) {
            Set<PairInt> set = groupList.get(i);
            Set<PairInt> remove = new HashSet<PairInt>();
            for (PairInt p : set) {
                int x = p.getX();
                int y = p.getY();
                int idx = input.getInternalIndex(x, y);

                float cieX = input.getCIEX(idx);
                float cieY = input.getCIEY(idx);

                Integer groupIndex = pixelToGroupMap.get(p);

                // key=groupIndex, value=number of pixels similar
                Map<Integer, Integer> groupSimilarCount = new HashMap<Integer, Integer>();

                // use cieXY, polar theta of cieXY, or rgb?
                for (int nIdx = 0; nIdx < dxs8.length; ++nIdx) {
                    int x2 = x + dxs8[nIdx];
                    int y2 = y + dys8[nIdx];
                    if ((x2 < 0) || (y2 < 0) || (x2 > (w - 1)) || (y2 > (h - 1))) {
                        continue;
                    }
                    int idx2 = input.getInternalIndex(x2, y2);

                    Integer groupIndex2 = pixelToGroupMap.get(new PairInt(x2, y2));

                    if (groupIndex.intValue() == groupIndex2.intValue()) {
                        continue;
                    }

                    float cieX2 = input.getCIEX(idx2);
                    float cieY2 = input.getCIEY(idx2);

                    float diffCIEX = Math.abs(cieX2 - cieX);
                    float diffCIEY = Math.abs(cieY2 - cieY);

                    if ((diffCIEX > diffLimit) || (diffCIEY > diffLimit)) {
                        continue;
                    }

                    Integer count = groupSimilarCount.get(groupIndex2);
                    if (count == null) {
                        groupSimilarCount.put(groupIndex2, Integer.valueOf(1));
                    } else {
                        groupSimilarCount.put(groupIndex2, Integer.valueOf(count.intValue() + 11));
                    }

                }
                if (groupSimilarCount.isEmpty()) {
                    continue;
                }
                for (Entry<Integer, Integer> entry : groupSimilarCount.entrySet()) {
                    if (entry.getValue() >= 6) {
                        // assign this group to pixel p
                        pixelToGroupMap.put(p, entry.getKey());
                        remove.add(p);
                        groupList.get(entry.getKey().intValue()).add(p);
                        break;
                    }
                }
            }
            for (PairInt rm : remove) {
                set.remove(rm);
            }
        }
    }

    private int[] findByHistogramLimitsForBlackAndWhite(ImageExt input) {

        //looking for limits of (r,g,b) <= (45,45,45) and > (180,180,180)
        int l0 = 45;
        int l0B = 70;
        int h0 = 245;

        List<Integer> avgL = new ArrayList<Integer>();

        List<Integer> avgH = new ArrayList<Integer>();

        for (int i = 0; i < input.getNPixels(); ++i) {
            int r = input.getR(i);
            int g = input.getG(i);
            int b = input.getB(i);
            if ((r <= l0) && (g <= l0) && (b <= l0B)) {
                int avg = (r + g + b)/3;
                avgL.add(Integer.valueOf(avg));
            } else if ((r > h0) && (g > h0) && (b > h0)) {
                int avg = (r + g + b)/3;
                avgH.add(Integer.valueOf(avg));
            }
        }
        
        int[] limits = new int[2];

        if (avgL.size() > 30) {
            HistogramHolder hist = Histogram.createSimpleHistogram(avgL);
            if (hist == null) {
                limits[0] = l0 - 1;
            } else {
                /*
                try {
                    hist.plotHistogram("black pixels", "black_" + MiscDebug.getCurrentTimeFormatted());
                } catch (IOException ex) {
                    Logger.getLogger(ImageSegmentation.class.getName()).log(Level.SEVERE, null, ex);
                }
                */
                int yMaxIdx = MiscMath.findYMaxIndex(hist.getYHist());
                int lastZeroIdx = MiscMath.findLastZeroIndex(hist);
                if ((lastZeroIdx != -1) && (lastZeroIdx < hist.getXHist().length)) {
                    limits[0] = Math.round(hist.getXHist()[lastZeroIdx]);
                } else if (yMaxIdx == -1) {
                    limits[0] = l0 - 1;
                } else {
                    //limits[0] = Math.round(hist.getXHist()[indexes.get(0).intValue()]);
                    limits[0] = Math.round(hist.getXHist()[yMaxIdx]);
                }
            }
        } else {
            limits[0] = l0 - 1;
        }

        if (avgH.size() > 30) {
            
            int[] q = ImageStatisticsHelper.getQuartiles(avgH);
            /*
            HistogramHolder hist = Histogram.createSimpleHistogram(avgH);
            if (hist == null) {
                limits[1] = h0;
            } else {
                List<Integer> indexes = MiscMath.findStrongPeakIndexes(hist, 0.1f);
                try {
                    hist.plotHistogram("hite pixels", "white_" + MiscDebug.getCurrentTimeFormatted());
                } catch (IOException ex) {
                    Logger.getLogger(ImageSegmentation.class.getName()).log(Level.SEVERE, null, ex);
                }
                if (indexes == null || indexes.isEmpty()) {
                    limits[1] = h0;
                } else {
                    limits[1] = Math.round(hist.getXHist()[indexes.get(0).intValue()]);
                }
            }
            */
            limits[1] = q[3];
        } else {
            limits[1] = h0;
        }

        return limits;
    }

    public Map<PairInt, Integer> calculatePolarCIEXY(ImageExt input, Set<PairInt> points) {

        Map<PairInt, Integer> map = new HashMap<PairInt, Integer>();

        CIEChromaticity cieC = new CIEChromaticity();

        for (PairInt p : points) {

            float cieX = input.getCIEX(p.getX(), p.getY());
            float cieY = input.getCIEY(p.getX(), p.getY());

            double thetaRadians = cieC.calculateXYTheta(cieX, cieY);

            int thetaDegrees = (int)Math.round(thetaRadians * 180./Math.PI);

            map.put(p, Integer.valueOf(thetaDegrees));
        }

        return map;
    }

    public float[] getValues(Map<PairInt, Integer> map) {
        float[] values = new float[map.size()];
        int count = 0;
        for (Entry<PairInt, Integer> entry : map.entrySet()) {
            values[count] = entry.getValue().floatValue();
            count++;
        }
        return values;
    }

    /**
     * NOT READY FOR USE.  STILL EXPERIMENTING.
     *
     * Calculates lists of black pixels, white pixels, grey pixels, and assigns
     * the remaining to CIEXY Lab color space, then creates a map of the
     * CIEX and CIEY points and uses density based clustering
     * (http://nking.github.io/two-point-correlation/)
     * to find clusters of points in CIE X, CIEY space,
     * then merges pixels in the grey list with adjacent clusters if
     * similar, and the final result is a list of pixel clusters, including the
     * black and white and remaining grey.
     *
     * Note that the color black is not defined in CIE XY color space and that
     * the color white is at the center of the space as a large circle so they
     * are not included in the density based clustering.
     * Note also that the remaining grey pixels which did not merge with the
     * cie xy pixel clusters, may be in separate groups already due to a
     * frequency based grouping for them.
     *
     * @param input
     * @param useBlur
     */
    public void applyPolarCIEXY(ImageExt input, boolean useBlur) {

        //TODO: improve the clustering results in two ways:
        // (1) for smaller ciexy clusters, merge with adjacent clusters if
        //     similar color
        // (2) any pixel with 7 neighbors of same color should be that color too

        if (useBlur) {
            ImageProcessor imageProcessor = new ImageProcessor();
            imageProcessor.blur(input, 1.0f);
        }

        Set<PairInt> blackPixels = new HashSet<PairInt>();

        Map<Integer, Collection<PairInt>> greyPixelMap = new HashMap<Integer, Collection<PairInt>>();

        Set<PairInt> whitePixels = new HashSet<PairInt>();

        Set<PairInt> points0 = new HashSet<PairInt>();

        populatePixelLists(input, points0, blackPixels, whitePixels, greyPixelMap);

        Map<PairInt, Integer> clrPolarCIEXYMap = calculatePolarCIEXY(input, points0);
        float[] clrPolarCIEXY = getValues(clrPolarCIEXYMap);
        float binWidth = 1;
        HistogramHolder hist = Histogram.createSimpleHistogram(//binWidth,
            clrPolarCIEXY, Errors.populateYErrorsBySqrt(clrPolarCIEXY));
        List<Integer> indexes = MiscMath.findStrongPeakIndexesDescSort(hist,
            0.1f);
        int[] binCenters = createBinCenters360(hist, indexes);

        List<Set<PairInt>> colorPixelGroups = assignToNearestPolarCIECluster(
            clrPolarCIEXYMap, binCenters);

        List<Set<PairInt>> greyPixelGroups = groupByPeaksForGrey(input, greyPixelMap);

        ImageExt imgExt = input.copyToImageExt();

        for (PairInt p : blackPixels) {
           imgExt.setRGB(p.getX(), p.getY(), 0, 0, 0);
        }

        for (PairInt p : whitePixels) {
           imgExt.setRGB(p.getX(), p.getY(), 255, 255, 255);
        }

        int gClr = 255;
        int s = 127/colorPixelGroups.size();
        for (Set<PairInt> set : colorPixelGroups) {
            for (PairInt p : set) {
               imgExt.setRGB(p.getX(), p.getY(), 0, gClr, 0);
            }
            gClr -= s;
        }

        /*
        int rClr = 255;
        s = 127/greyPixelMap.size();
        for (Entry<Integer, Collection<PairInt>> entry : greyPixelMap.entrySet()) {
            for (PairInt p : entry.getValue()) {
                imgExt.setRGB(p.getX(), p.getY(), rClr, 0, 0);
            }
            rClr -= s;
        }
        */
        int rClr = 255;
        s = 127/greyPixelGroups.size();
        for (Set<PairInt> set : greyPixelGroups) {
            for (PairInt p : set) {
               imgExt.setRGB(p.getX(), p.getY(), rClr, 0, 0);
            }
            rClr -= s;
        }

        input.resetTo(imgExt);
    }

    private int[] createBinCenters360(HistogramHolder hist, List<Integer> indexes) {

        if (indexes.isEmpty()) {
            return new int[0];
        }

        int n = indexes.size();

        if (n == 1) {

            int[] binCenters = new int[n + 1];

            /*
            examples for n=1
                --  90
            180
                -- 270
            ----------------
                --  90-176=  360-(176-90) = 274
            4
                -- 270-176=  94
            ----------------
                --  90+(330-180)=  240
            330
                -- 270+(330-180)=  270 + 150 - 360 = 60
            */

            int vc = Math.round(hist.getXHist()[indexes.get(0).intValue()]);

            if (vc == 180) {
                binCenters[0] = 90;
                binCenters[1] = 270;
            } else if (vc < 180) {
                int delta = 180 - vc;
                binCenters[0] = 360 - (delta - 90);
                binCenters[1] = 270 - delta;
            } else {
                int delta = vc - 180;
                binCenters[0] = 90 + delta;
                binCenters[1] = 270 + delta - 360;
            }

            return binCenters;
        }

        int[] binCenters = new int[n];

        for (int i = 0; i < n; ++i) {

            int idx = indexes.get(i);

            binCenters[i] = Math.round(hist.getXHist()[idx]);
        }

        return binCenters;
    }

    private int[] createBinCenters255(HistogramHolder hist, List<Integer> indexes) {

        if (indexes.isEmpty()) {
            return new int[0];
        }

        int n = indexes.size();

        if (n == 1) {

            int[] binCenters = new int[n + 1];

            int vc = Math.round(hist.getXHist()[indexes.get(0).intValue()]);

            binCenters[0] = vc/2;
            binCenters[1] = (255 + vc)/2;

            return binCenters;
        }

        int[] binCenters = new int[n];

        for (int i = 0; i < n; ++i) {

            int idx = indexes.get(i);

            binCenters[i] = Math.round(hist.getXHist()[idx]);
        }

        return binCenters;
    }

    private List<Set<PairInt>> groupByPeaksForGrey(ImageExt input,
        Map<Integer, Collection<PairInt>> greyPixelMap) {

        Set<PairInt> points = new HashSet<PairInt>();
        for (Entry<Integer, Collection<PairInt>> entry : greyPixelMap.entrySet()) {
            points.addAll(entry.getValue());
        }

        Map<PairInt, Integer> polarCIEXYMap = calculatePolarCIEXY(input, points);
        float[] polarCIEXY = getValues(polarCIEXYMap);
        float binWidth = 1;
        HistogramHolder hist = Histogram.createSimpleHistogram(binWidth,
            polarCIEXY, Errors.populateYErrorsBySqrt(polarCIEXY));
        List<Integer> indexes = MiscMath.findStrongPeakIndexesDescSort(hist,
            0.1f);
        int[] binCenters = createBinCenters360(hist, indexes);

        List<Set<PairInt>> pixelGroups = assignToNearestPolarCIECluster(
            polarCIEXYMap, binCenters);

        return pixelGroups;

    }

    /**
     * a greyscale segmentation algorithm similar to the KMPP, but does not use
     * a random number generator to find seeds of intensity bins.
     * @param input 
     */
    public void applyGreyscaleHistogram(GreyscaleImage input) {

        float[] values = new float[input.getNPixels()];
        for (int i = 0; i < input.getNPixels(); ++i) {
            int v = input.getValue(i);
            values[i] = v;
        }
        
        int[] binCenters = determineGreyscaleBinCenters(values);
        
        assignToNearestCluster(input, binCenters);
    }
    
    /**
     * a greyscale segmentation algorithm similar to the KMPP, but does not use
     * a random number generator to find seeds of intensity bins.
     * 
     */
    public int[] determineGreyscaleBinCenters(float[] values) {

        /*
        Arrays.sort(values);
        
        float median = values[values.length/2];
        
        float iqr = values[(int)(0.75f*values.length)] - values[(int)(0.25f*values.length)];
        
        float low = values[(int)0.25f*values.length] - 1.5f*iqr;
        float high = values[(int)0.75f*values.length] + 1.5f*iqr;
        int nTot = values.length;
        int nLow = 0;
        int nHigh = 0;
        for (float v : values) {
            if (v < low) {
                nLow++;
            } else if (v > high) {
                nHigh++;
            }
        }
        float nFracOut = (float)(nLow + nHigh)/(float)nTot;
        
        float[] avgAndStDev = MiscMath.getAvgAndStDev(values);

        //3*(MEAN – MEDIAN)/STANDARD DEVIATION
        float skew0 = (3.f*(avgAndStDev[0] - median))/avgAndStDev[1];
        
        log.info(String.format(
            "mean=%f median=%f stdev=%f\niqr=%f skew0=%f nFracOut=%f",
            avgAndStDev[0], median, avgAndStDev[1], iqr, skew0, nFracOut));
        */
        
        float binWidth = 3;//10;
        HistogramHolder hist = Histogram.createSimpleHistogram(0, 255, binWidth,
            values, Errors.populateYErrorsBySqrt(values));

        List<Integer> indexes = MiscMath.findStrongPeakIndexesDescSort(hist,
            0.05f);
        /*
        try {
            hist.plotHistogram("greyscale", "1");
        } catch (IOException ex) {
            Logger.getLogger(ImageSegmentation.class.getName()).log(Level.SEVERE, null, ex);
        }*/
        
        int k = 3;
        int idx = (indexes.size() < k) ? indexes.size() : k; 
        indexes = indexes.subList(0, idx);
      
        int[] binCenters = createBinCenters255(hist, indexes);

        return binCenters;
    }

    /**
     * segmentation algorithm using cieXY and rgb to make pixel categories,
     * and grow lists, then histograms to rescale.
     * 
     * @param input
     * @return 
     */
    public GreyscaleImage createGreyscale3(ImageExt input) {

        Set<PairInt> blackPixels = new HashSet<PairInt>();

        Map<PairInt, Float> colorPixelMap = new HashMap<PairInt, Float>();
        
        Set<PairInt> unassignedPixels = new HashSet<PairInt>();

        populatePixelLists2(input, blackPixels, colorPixelMap, unassignedPixels);
       
        int nTotC = blackPixels.size() + colorPixelMap.size() + unassignedPixels.size();
        assert(nTotC == input.getNPixels());
                
        // grow black pixels from unassigned pixels if within a tolerance of rgb
        growPixelsWithinRGBTolerance(input, blackPixels, unassignedPixels, 5);
        
        // then grow colorPixelMap from unassignedPixels if within a tolerance of rgb
        Set<PairInt> addToColor = findPixelsWithinRGBTolerance(input, 
            colorPixelMap.keySet(), unassignedPixels, 5);
        
        // reassign colorPixelMap to averaged color for adjacent thetas within tolerance
        // (1) try just colorPixelMap for this step
        // (2) try adding addToColor to colorPixelMap for this step
        //    Map<PairInt, Float> colorPixelMap2 = createPolarCIEXYMap(input, addToColor);
        //    colorPixelMap.putAll(colorPixelMap2);
               
        int toleranceInValue = determineToleranceForIllumCorr(colorPixelMap, nTotC);

        if (toleranceInValue > -1 && !colorPixelMap.isEmpty()) {
            correctForIllumination(input, toleranceInValue, colorPixelMap);
        }

ImageExt tmpInput = input.copyToImageExt();
for (PairInt p : unassignedPixels) {
tmpInput.setRGB(p.getX(), p.getY(), 255, 0, 0);
}        
MiscDebug.writeImage(tmpInput, "_after_illumc0_unassigned_pix" + MiscDebug.getCurrentTimeFormatted());
MiscDebug.writeImage(input, "_after_illumc0_" + MiscDebug.getCurrentTimeFormatted());

        Map<PairInt, Float> colorPixelMap2 = createPolarCIEXYMap(input, addToColor);
        if (toleranceInValue > -1 && !colorPixelMap2.isEmpty()) {
            correctForIllumination(input, toleranceInValue, colorPixelMap2);
        }

tmpInput = input.copyToImageExt();
for (PairInt p : colorPixelMap2.keySet()) {
tmpInput.setRGB(p.getX(), p.getY(), 255, 0, 0);
}        
MiscDebug.writeImage(tmpInput, "_after2_illumc0_pix" + MiscDebug.getCurrentTimeFormatted());   
MiscDebug.writeImage(input, "_after2_illumc0_" + MiscDebug.getCurrentTimeFormatted());

        Map<PairInt, Float> colorPixelMap3 = createPolarCIEXYMap(input, unassignedPixels);
        if (toleranceInValue > -1 && !colorPixelMap3.isEmpty()) {
            correctForIllumination(input, toleranceInValue, colorPixelMap3);
        }
tmpInput = input.copyToImageExt();
for (PairInt p : colorPixelMap3.keySet()) {
tmpInput.setRGB(p.getX(), p.getY(), 255, 0, 0);
}        
MiscDebug.writeImage(tmpInput, "_after3_illumc0_pix" + MiscDebug.getCurrentTimeFormatted());   
MiscDebug.writeImage(input, "_after3_illumc0_" + MiscDebug.getCurrentTimeFormatted());

        if (nTotC != blackPixels.size()) {
            //TODO: not sure will use this... might consider same method w/ rgb
            toleranceInValue = 2;
            correctForIllumination(input, blackPixels, toleranceInValue);
        }

        GreyscaleImage img = input.copyToGreyscale();
        
        MiscDebug.writeImage(img, "_before_greyscale_bins_" + MiscDebug.getCurrentTimeFormatted());

        float[] cValues = new float[input.getNPixels()];
        for (int i = 0; i < input.getNPixels(); ++i) {
            int v = img.getValue(i);
            cValues[i] = v;
        }

        int[] binCenters = determineGreyscaleBinCenters(cValues);

        //ImageExt img2 = input.copyToImageExt();

        assignToNearestCluster(img, binCenters);
        MiscDebug.writeImage(img, "_after_greyscale_bins_" + MiscDebug.getCurrentTimeFormatted());
        
        /*    
        try {
            applyUsingKMPP(img, 3);
        } catch (IOException ex) {
            Logger.getLogger(ImageSegmentation.class.getName()).log(Level.SEVERE, null, ex);
        } catch (NoSuchAlgorithmException ex) {
            Logger.getLogger(ImageSegmentation.class.getName()).log(Level.SEVERE, null, ex);
        }*/
        
        for (PairInt p : blackPixels) {
            img.setValue(p.getX(), p.getY(), 0);
            //img2.setRGB(p.getX(), p.getY(), 255, 0, 0);
        }

MiscDebug.writeImage(img, "_end_seg_" + MiscDebug.getCurrentTimeFormatted());
                      
        return img;
    }
    
    protected Map<PairInt, Float> createPolarCIEXYMap(ImageExt input,
        Set<PairInt> points) {

        Map<PairInt, Float> thetaMap = new HashMap<PairInt, Float>();

        //TODO: remove imgCp when finished debugging
        //ImageExt imgCp = input.copyToImageExt();
        int w = input.getWidth();
        int h = input.getHeight();

        CIEChromaticity cieC = new CIEChromaticity();

        double thetaFactor = (double) 255 / (double) 360;

        // to set the non member colors, need to traverse all pixels.
        for (int col = 0; col < w; ++col) {
            for (int row = 0; row < h; row++) {
                PairInt p = new PairInt(col, row);
                if (points.contains(p)) {
                    float cieX = input.getCIEX(col, row);
                    float cieY = input.getCIEY(col, row);
                    double theta = thetaFactor * ((cieC.calculateXYTheta(cieX, cieY) * 180. / Math.PI));
          //          int thetaCIEXY = (int)Math.round(theta);
                    //          imgCp.setRGB(col, row, thetaCIEXY, thetaCIEXY, thetaCIEXY);

                    thetaMap.put(p, Float.valueOf((float) theta));
                } else {
                    //         imgCp.setRGB(col, row, 0, 255, 0);
                }
            }
        }

        //MiscDebug.writeImage(imgCp, "polarciexy_" + MiscDebug.getCurrentTimeFormatted());

        return thetaMap;
    }

    private void correctForIllumination(ImageExt input, Set<PairInt> points, 
        int toleranceInValue) {
        
        if (points.isEmpty() || (toleranceInValue == 0)) {
            return;
        }
        
        Map<PairInt, Float> thetaMap = createPolarCIEXYMap(input, points);
       
        correctForIllumination(input, toleranceInValue, thetaMap);
    }
    
    private void correctForIllumination(ImageExt input,
        int toleranceInValue, Map<PairInt, Float> thetaMap) {
        
        /*
        in CIE XY color space, reassign colors:
        
        for points:
        -- convert to polar theta in CIE XY space and plot in image
           as greyscale without changes.
        
        -- DFS traversal through the theta values of pixels in points
           to make contigous groups of pixels that are within a tolerance of
           theta value of one another.
        
        -- reassign the average rgb color to those pixels in a group
        */
               
        int w = input.getWidth();
        int h = input.getHeight();
        
        // find pixels similar in color that are contiguous
        
        DFSConnectedGroupsFinder2 groupFinder = new DFSConnectedGroupsFinder2();
        groupFinder.findConnectedPointGroups(thetaMap, 360, toleranceInValue,
            w, h, false);
        
        int nGroups = groupFinder.getNumberOfGroups();
        
        // calc avg rgb of each contiguous group and reset rgb of all in group 
        // to the avg
        
        for (int i = 0; i < nGroups; ++i) {
            
            Set<PairInt> group = groupFinder.getXY(i);
        
            if (group.size() == 1) {
                continue;
            }
            
            float sumR = 0;
            float sumG = 0;
            float sumB = 0;
            for (PairInt p : group) {
                int pixIdx = input.getInternalIndex(p.getX(), p.getY());
                sumR += input.getR(pixIdx);
                sumG += input.getG(pixIdx);
                sumB += input.getB(pixIdx);
            }
            int avgR = Math.round(sumR/(float)group.size());
            int avgG = Math.round(sumG/(float)group.size());
            int avgB = Math.round(sumB/(float)group.size());

            for (PairInt p : group) {
                input.setRGB(p.getX(), p.getY(), avgR, avgG, avgB);
            }
        }        
    }

    private int determineToleranceForIllumCorr(Map<PairInt, Float> thetaMap,
        int nTotCategoryMembers) {
        
        if (thetaMap.isEmpty()) {
            return -1;
        }
        
        // using a histogram of color in
        // determining whether toleranceInValue should be 1 or 4.
        // 1 is better for low contrast regions like the brown & lowe 2003 mtns
        // but 4 is better for many colors, large contrast, and large dynamic
        // range such as seen in the merton college test images
        
        float[] values = new float[thetaMap.size()];
        int count = 0;
        for (Entry<PairInt, Float> entry : thetaMap.entrySet()) {
            values[count] = entry.getValue().floatValue();
            count++;
        }
        
        int tolerance = 1;
        
        HistogramHolder hist = Histogram.createSimpleHistogram(values, 
            Errors.populateYErrorsBySqrt(values));
        try {
            hist.plotHistogram("cieTheta", "_cietheta_");
        } catch (IOException ex) {
            Logger.getLogger(ImageSegmentation.class.getName()).log(Level.SEVERE, null, ex);
        }
        
        if (hist != null) {
            List<Integer> indexes = MiscMath.findStrongPeakIndexes(hist, 0.05f);
            
            //TODO: for rare case that image has only one color, and the pixels
            // are all in thetaMap, we want the invoker not to remove illumination
            // because that may be the only contrast in the image.
            // an example of that happening wouldbe a blue sky and blue mountains
            // very similar in color.
            if ((indexes.size() == 1) && (nTotCategoryMembers == thetaMap.size())) {
                return -1;
            }
            
            // from lowest to highest peak, are there continuous values >= 1/4 lowest peak?
            if (indexes.size() > 1) {
                
                Collections.sort(indexes);
                
                int minPeak = Integer.MAX_VALUE;
                for (Integer index : indexes) {
                    int v = hist.getYHist()[index.intValue()];
                    if (v < minPeak) {
                        minPeak = v;
                    }
                }
                
                int nBelow = 0;
                int n = 0;
                int limit = minPeak/4;
                for (int i = indexes.get(0).intValue(); 
                    i < indexes.get(indexes.size() - 1).intValue(); ++i) {
                    int v = hist.getYHist()[i];
                    if (v < limit) {
                        nBelow++;
                    }
                    n++;
                }
                
                if (((float)nBelow/(float)n) < 0.5) {
                    tolerance = 3;
                }
            }
        }
        
        return tolerance;
    }

    private int correctPointsBoundaries(ImageExt input, 
        Set<PairInt> points, Set<PairInt> blackPixels, 
        Set<PairInt> whitePixels, Set<PairInt> greyPixels) {
        
         // boundary corrections:
        // iterate over points0 to find boundary points and their neighbors.  
        // a boundary point has at least one non-points0 point neighbor.
        // for the boundary points and each of their neigbhors,
        // determine whether should move the point to another group
        // (or add it to its own).
        // the local rgb or ciexy can be computed and kept in a hash to avoid
        // re-calculating.
        // the local color will be the 9 point region members only regions
        // 
        
        int imageWidth = input.getWidth(); 
        int imageHeight = input.getHeight();
        
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        
        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;
        
        Set<PairInt> addP = new HashSet<PairInt>();
        Set<PairInt> rmP = new HashSet<PairInt>();
        Set<PairInt> moved = new HashSet<PairInt>();
        
        for (PairInt p : points) {
            int x = p.getX();
            int y = p.getY();
            
            int nNonNeighbors = 0;            
            for (int i = 0; i < dxs.length; ++i) {
                int x2 = x + dxs[i];
                int y2 = y + dys[i];
                if ((x2 < 0) || (x2 > (imageWidth - 1)) || (y2 < 0) ||
                    (y2 > (imageHeight - 1))) {
                    continue;
                }
                PairInt p2 = new PairInt(x2, y2);                
                if (!points.contains(p2)) {
                   nNonNeighbors++;
                }
            }
            
            if (nNonNeighbors == 0) {
                continue;
            }
            
            Set<PairInt> neighbors = new HashSet<PairInt>();
            Set<PairInt> nonNB = new HashSet<PairInt>();
            Set<PairInt> nonNW = new HashSet<PairInt>();
            Set<PairInt> nonNG = new HashSet<PairInt>();
            
            for (int i = 0; i < dxs.length; ++i) {
                int x2 = x + dxs[i];
                int y2 = y + dys[i];                
                if ((x2 < 0) || (x2 > (imageWidth - 1)) || (y2 < 0) ||
                    (y2 > (imageHeight - 1))) {
                    continue;
                }
                PairInt p2 = new PairInt(x2, y2);
                if (!moved.contains(p2)) {
                    if (points.contains(p2)) {
                        neighbors.add(p2);
                    } else if (blackPixels.contains(p2)) {
                        nonNB.add(p2);
                    } else if (whitePixels.contains(p2)) {
                        nonNW.add(p2);
                    } else if (greyPixels.contains(p2)) {
                        nonNG.add(p2);
                    }
                }
            }
            
            // calc the avg cieX, cieY of neighbors
            // calc the avg cieX, cieY of non-neighbors
            // if p's cieXY is closer to neighbors, do not move it, but do
            // look at each of the nonNeighbors to see whether to move it...
            //float[] avgCIEXYNeighbors = ImageStatisticsHelper.calculateAvgCIEXY(
            //    input, neighbors);
            float[] avgRGBNeighbors = ImageStatisticsHelper.calculateAvgRGB(
                input, neighbors);
         
            //float pCIEX = input.getCIEX(x, y);
            //float pCIEY = input.getCIEY(x, y);
            int pR = input.getR(x, y);
            int pG = input.getG(x, y);
            int pB = input.getB(x, y);
            
            double distSq0 = Math.pow((avgRGBNeighbors[0] - pR), 2) + 
                Math.pow((avgRGBNeighbors[1] - pG), 2) +
                Math.pow((avgRGBNeighbors[2] - pB), 2);
            
            double minDistSq = distSq0;
            
            // 0 = neighbors, 1=blackPixels, 2=whitePixels, 3=greyPixels
            int minDistGroup = 0;
           
            if (!nonNB.isEmpty()) {
                float[] avgRGB = ImageStatisticsHelper.calculateAvgRGB(
                    input, nonNB);
                double distSq = Math.pow((avgRGB[0] - pR), 2) + 
                    Math.pow((avgRGB[1] - pG), 2) +
                    Math.pow((avgRGB[2] - pB), 2);
                if (distSq < minDistSq) {
                    minDistSq = distSq;
                    minDistGroup = 1;
                }
            }
            if (!nonNW.isEmpty()) {
                float[] avgRGB = ImageStatisticsHelper.calculateAvgRGB(
                    input, nonNW);
                double distSq = Math.pow((avgRGB[0] - pR), 2) + 
                    Math.pow((avgRGB[1] - pG), 2) +
                    Math.pow((avgRGB[2] - pB), 2);
                if (distSq < minDistSq) {
                    minDistSq = distSq;
                    minDistGroup = 2;
                }
            }
            if (!nonNG.isEmpty()) {
                float[] avgRGB = ImageStatisticsHelper.calculateAvgRGB(
                    input, nonNG);
                double distSq = Math.pow((avgRGB[0] - pR), 2) + 
                    Math.pow((avgRGB[1] - pG), 2) +
                    Math.pow((avgRGB[2] - pB), 2);
                if (distSq < minDistSq) {
                    minDistSq = distSq;
                    minDistGroup = 3;
                }
            }
            
            // 0 = neighbors, 1=blackPixels, 2=whitePixels, 3=greyPixels
            if (minDistGroup > 0) {
                if (minDistGroup == 1) {
                    blackPixels.add(p);
                } else if (minDistGroup == 2) {
                    whitePixels.add(p);
                } else if (minDistGroup == 3) {
                    greyPixels.add(p);
                }
                rmP.add(p);
                moved.add(p);
            }
            
            // determine if the non neighbors should be moved into point group
            if (!nonNB.isEmpty()) {
                for (PairInt p2 : nonNB) {
                    if (moved.contains(p2)) {
                        continue;
                    }
                    int x2 = p2.getX();
                    int y2 = p2.getY();
                    int p2R = input.getR(x2, y2);
                    int p2G = input.getG(x2, y2);
                    int p2B = input.getB(x2, y2);
                    
                    Set<PairInt> nbN = curveHelper.findNeighbors(x2, y2, blackPixels);
                    float[] avgRGB = ImageStatisticsHelper.calculateAvgRGB(input, nbN);
                    
                    double distSqToNeighbors = Math.pow((avgRGB[0] - p2R), 2) + 
                        Math.pow((avgRGB[1] - p2G), 2) + Math.pow((avgRGB[2] - p2B), 2);
                    
                    double distSqToP = Math.pow((pR - p2R), 2) + 
                        Math.pow((pG - p2G), 2) + Math.pow((pB - p2B), 2);

                    if (distSqToP < distSqToNeighbors) {
                        blackPixels.remove(p2);
                        addP.add(p2);
                        moved.add(p2);
                    }
                }
            }
            if (!nonNW.isEmpty()) {
                for (PairInt p2 : nonNW) {
                    if (moved.contains(p2)) {
                        continue;
                    }
                    int x2 = p2.getX();
                    int y2 = p2.getY();
                    int p2R = input.getR(x2, y2);
                    int p2G = input.getG(x2, y2);
                    int p2B = input.getB(x2, y2);
                    
                    Set<PairInt> nbN = curveHelper.findNeighbors(x2, y2, whitePixels);
                    float[] avgRGB = ImageStatisticsHelper.calculateAvgRGB(input, nonNW);
                    
                    double distSqToNeighbors = Math.pow((avgRGB[0] - p2R), 2) + 
                        Math.pow((avgRGB[1] - p2G), 2) + Math.pow((avgRGB[2] - p2B), 2);
                    
                    double distSqToP = Math.pow((pR - p2R), 2) + 
                        Math.pow((pG - p2G), 2) + Math.pow((pB - p2B), 2);
                    
                    if (distSqToP < distSqToNeighbors) {
                        whitePixels.remove(p2);
                        addP.add(p2);
                        moved.add(p2);
                    }
                }
            }
            if (!nonNG.isEmpty()) {
                for (PairInt p2 : nonNG) {
                    if (moved.contains(p2)) {
                        continue;
                    }
                    int x2 = p2.getX();
                    int y2 = p2.getY();
                    int p2R = input.getR(x2, y2);
                    int p2G = input.getG(x2, y2);
                    int p2B = input.getB(x2, y2);
                    
                    Set<PairInt> nbN = curveHelper.findNeighbors(x2, y2, greyPixels);
                    float[] avgRGB = ImageStatisticsHelper.calculateAvgRGB(input, nonNG);
                    
                    double distSqToNeighbors = Math.pow((avgRGB[0] - p2R), 2) + 
                        Math.pow((avgRGB[1] - p2G), 2) + Math.pow((avgRGB[2] - p2B), 2);
                    
                    double distSqToP = Math.pow((pR - p2R), 2) + 
                        Math.pow((pG - p2G), 2) + Math.pow((pB - p2B), 2);
                    
                    if (distSqToP < distSqToNeighbors) {
                        greyPixels.remove(p2);
                        addP.add(p2);
                        moved.add(p2);
                    }
                }
            }
        }
        
        for (PairInt p : addP) {
            points.add(p);
        }
        for (PairInt p : rmP) {
            points.remove(p);
        }
        
        return moved.size();
    }

    /**
     * an algorithm to traverse the image and store black pixels in a list
     * then use cieXY to place pixels in colorPixelMap if they are not in the
     * large white central space, else place them in unassignedPixels.
     * @param input
     * @param blackPixels
     * @param colorPixelMap
     * @param unassignedPixels 
     */
    private void populatePixelLists2(ImageExt input, Set<PairInt> blackPixels, 
        Map<PairInt, Float> colorPixelMap, Set<PairInt> unassignedPixels) {
                
        int w = input.getWidth();
        int h = input.getHeight();

        CIEChromaticity cieC = new CIEChromaticity();

        // looking for limits in peaks of (r,g,b) <= (45,45,45) and > (191,191,191)
        int[] whiteBlackLimits = findByHistogramLimitsForBlackAndWhite(input);
        // overriding:
        
        for (int i = 0; i < w; ++i) {
            for (int j = 0; j < h; ++j) {

                int idx = input.getInternalIndex(i, j);

                int r = input.getR(idx);
                int g = input.getG(idx);
                int b = input.getB(idx);

                int avg = (r + g + b)/3;
                if (avg <= whiteBlackLimits[0]) {
                    blackPixels.add(new PairInt(i, j));
                    continue;
                }
                
                float cieX = input.getCIEX(idx);
                float cieY = input.getCIEY(idx);
                /*    
                float c = r + g + b;
                float rd = Math.abs(0.3333f - ((float)r/c));
                float bd = Math.abs(0.3333f - ((float)b/c));
                */
                if (cieC.isInLargeWhiteCenter(cieX, cieY)) {
                    //log.info(String.format(
                    //"rgb=(%3d, %3d, %3d)  rdiv,bdiv=(%.2f, %.2f) cie=(%.3f, %.3f) xy=(%3d, %3d)", 
                    //r, g, b, rd, bd, cieX, cieY, i, j));
                    unassignedPixels.add(new PairInt(i, j));
                    continue;
                }
                
                double thetaRadians = cieC.calculateXYTheta(cieX, cieY);

                colorPixelMap.put(new PairInt(i, j), 
                    Float.valueOf((float)thetaRadians));
            }
        }
    }

    private void growPixelsWithinRGBTolerance(Image img, 
        Set<PairInt> assignedPixels, Set<PairInt> unassignedPixels, int tolRGB) {
        
        Set<PairInt> addTo = findPixelsWithinRGBTolerance(img, 
            assignedPixels, unassignedPixels, tolRGB);
        
        assignedPixels.addAll(addTo);
    }

    /**
     * traverses assignedPixels and finds adjacent unassignedPixels that are 
     * within tolRGB of r, g, b colors and removes those points from 
     * unassignedPixels and returns them as a set.
     * @param input
     * @param assignedPixels
     * @param unassignedPixels
     * @param tolRGB
     * @return 
     */
    private Set<PairInt> findPixelsWithinRGBTolerance(Image img, 
        Set<PairInt> assignedPixels, Set<PairInt> unassignedPixels, int tolRGB) {
        
        Set<PairInt> addToAssigned = new HashSet<PairInt>();
        
        int w = img.getWidth();
        int h = img.getHeight();
        
        int[] dxs8 = Misc.dx8;
        int[] dys8 = Misc.dy8;
        
        Set<PairInt> visited = new HashSet<PairInt>();
        
        Stack<PairInt> stack = new Stack<PairInt>();
        stack.addAll(assignedPixels);
        
        while (!stack.isEmpty()) {
            PairInt uPoint = stack.pop(); 
            if (visited.contains(uPoint)) {
                continue;
            }
            int x = uPoint.getX();
            int y = uPoint.getY();
            
            int r = img.getR(x, y);
            int g = img.getG(x, y);
            int b = img.getB(x, y);
            
            for (int i = 0; i < dxs8.length; ++i) {
                int x2 = x + dxs8[i];
                int y2 = y + dys8[i];
                if ((x2 < 0) || (x2 > (w - 1)) || (y2 < 0) || (y2 > (h - 1))) {
                    continue;
                }        
                PairInt vPoint = new PairInt(x2, y2);
                if (!unassignedPixels.contains(vPoint)) {
                    continue;
                }
                int r2 = img.getR(x2, y2);
                int g2 = img.getG(x2, y2);
                int b2 = img.getB(x2, y2);
                int diffR = Math.abs(r2 - r);
                int diffG = Math.abs(g2 - g);
                int diffB = Math.abs(b2 - b);
                
                if ((diffR < tolRGB) && (diffG < tolRGB) && (diffB < tolRGB)) {
                    unassignedPixels.remove(vPoint);
                    addToAssigned.add(vPoint);
                    stack.add(vPoint);
                }
            }
            
            visited.add(uPoint);
        }
        
        return addToAssigned;
    }
}
