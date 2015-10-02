package algorithms.imageProcessing;

import algorithms.compGeometry.clustering.KMeansPlusPlus;
import algorithms.compGeometry.clustering.KMeansPlusPlusFloat;
import algorithms.misc.Histogram;
import algorithms.misc.HistogramHolder;
import algorithms.misc.Misc;
import algorithms.misc.MiscDebug;
import algorithms.misc.MiscMath;
import algorithms.util.Errors;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.PolygonAndPointPlotter;
import algorithms.util.ResourceFinder;
import com.climbwithyourfeet.clustering.DTClusterFinder;
import java.io.IOException;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
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

        int[] binCenters = instance.getCenters();

        for (int col = 0; col < input.getWidth(); col++) {

            for (int row = 0; row < input.getHeight(); row++) {

                int v = input.getValue(col, row);

                for (int i = 0; i < binCenters.length; i++) {

                    int vc = binCenters[i];

                    int bisectorBelow = ((i - 1) > -1) ?
                        ((binCenters[i - 1] + vc) / 2) : 0;

                    int bisectorAbove = ((i + 1) > (binCenters.length - 1)) ?
                        255 : ((binCenters[i + 1] + vc) / 2);

                    if ((v >= bisectorBelow) && (v <= bisectorAbove)) {

                        input.setValue(col, row, vc);

                        break;
                    }
                }
            }
        }
    }
    
    /**
     * applies a blue of sigma=1 to image,
     * converts each pixel color to the polar angle of CIE XY Lab color space
     * with an origin of (0.35, 0.35) and uses a histogram binning of kColors=8,
     * then maps those bins to 0 to 255,
     * then replaces a pixel if 7 or its neighbors have the same color,
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
     * then replaces a pixel if 7 or its neighbors have the same color,
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

        // ----replace pixel, if 7 or more neighbors have same color -----
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
     * range 0 to 255, then replaces a pixel if it has 7 neighbors of same 
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

        // ----replace pixel, if 7 or more neighbors have same color -----
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

        return applyUsingCIEXYPolarThetaThenHistogram(input, 254);
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
        // method name may change to apply.  might average the cluster color and apply it to points

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

            PairIntWithIndex p2 = new PairIntWithIndex(x, y, p.pixIdx);
            List<PairIntWithIndex> list2 = pointsMap.get(p2);
            if (list2 == null) {
                list2 = new ArrayList<PairIntWithIndex>();
                pointsMap.put(p2, list2);
            }
            // because this is a list, this will eventually be present twice:
            //list2.add(p2);

            for (PairIntWithIndex p0 : pointsMap0.get(p)) {
                PairIntWithIndex p3 = new PairIntWithIndex(
                    p0.getX() - minCIEX, p0.getY() - minCIEY, p0.pixIdx);
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

        int nGroups = clusterFinder.getNumberOfClusters();

        List<Set<PairInt>> groupList = new ArrayList<Set<PairInt>>();

        for (int k = 0; k < nGroups; ++k) {
            
            Set<PairIntWithIndex> group = clusterFinder.getCluster(k);

            Set<PairInt> coordPoints = new HashSet<PairInt>();

            for (PairIntWithIndex p : group) {

                int idx = p.pixIdx;
                int xCoord = input.getCol(idx);
                int yCoord = input.getRow(idx);

                PairInt pCoord = new PairInt(xCoord, yCoord);
                coordPoints.add(pCoord);

                // include the other points of/ same color
                List<PairIntWithIndex> list = pointsMap.get(p);
                assert(list != null);
                for (PairIntWithIndex p3 : list) {
                    int idx3 = p3.pixIdx;
                    int xCoord3 = input.getCol(idx3);
                    int yCoord3 = input.getRow(idx3);
                    pCoord = new PairInt(xCoord3, yCoord3);
                    coordPoints.add(pCoord);
                }
            }

            groupList.add(coordPoints);
        }
        
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
                        int idx3 = p3.pixIdx;
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
        Map<Integer, Collection<PairInt>> thetaPointMap = createThetaCIEXYMap(points0,
            input, binWidth);
       
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
        
        PairIntArray peaks = findPeaksInThetaPointMap(orderedThetaKeys, 
            thetaPointMap, 
            Math.round(fracFreqLimit * maxFreq));
        
        int[] minMaxXY = MiscMath.findMinMaxXY(peaks);
        
        // ----- debug ---
        // plot the points as an image to see the data first
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
        
        // --- debug
        float[] xPoints = new float[greyPixelMap.size()];
        float[] yPoints = new float[greyPixelMap.size()];
        int count = 0;
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
        List<Set<PairInt>> greyPixelGroups, 
        List<Set<PairInt>> groupList, Set<PairInt> blackPixels,
        Set<PairInt> whitePixels) {
                
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
                            colorClusterIdx = colorClusterIndex.intValue();
                            diffR = Math.abs(r2 - r); 
                            diffG = Math.abs(g2 - g);
                            diffB = Math.abs(b2 - b);
                        }
                    }
                }
                if (minDiffBlack < 75) {
                    blackPixels.add(greyP);
                    remove.add(greyP);
                    continue;
                } else if (minDiffWhite < 75) {
                    whitePixels.add(greyP);
                    remove.add(greyP);
                    continue;
                } else {
                    if (colorClusterIdx != -1) {
                        //add to color cluster and remove from grey list
                        groupList.get(colorClusterIdx).add(greyP);
                        remove.add(greyP);
                        continue;
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
            greyPixelGroups.get(i).clear();
        }
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

    private void populatePixelLists(ImageExt input, Set<PairInt> points0, 
        Set<PairInt> blackPixels, Set<PairInt> whitePixels, 
        Map<Integer, Collection<PairInt>> greyPixelMap) {
        
        int w = input.getWidth();
        int h = input.getHeight();
        
        CIEChromaticity cieC = new CIEChromaticity();
        
        for (int i = 0; i < w; ++i) {
            for (int j = 0; j < h; ++j) {

                int idx = input.getInternalIndex(i, j);

                int r = input.getR(idx);
                int g = input.getG(idx);
                int b = input.getB(idx);

                //dark grey, such as r,g,b=105,105,105?
                if ((r <= 32) && (g <= 32) && (b <= 32)) {
                    blackPixels.add(new PairInt(i, j));
                    continue;
                }

                float cx = input.getCIEX(idx);
                float cy = input.getCIEY(idx);

                if (cieC.isWhite2(cx, cy)) {
                    //grey will be binned into clusters by avgRGB and peak frequency
                    if ((r <= 191) && (g <= 191) && (b <= 191)) {
                        Integer avgRGB = Integer.valueOf(Math.round((r + g + b)/3.f));
                        Collection<PairInt> set = greyPixelMap.get(avgRGB);
                        if (set == null) {
                            set = new HashSet<PairInt>();
                            greyPixelMap.put(avgRGB, set);
                        }
                        set.add(new PairInt(i, j));
                    } else {
                        whitePixels.add(new PairInt(i, j));
                    }
                    continue;
                }

                points0.add(new PairInt(i, j));                
            }
        }
    }

    public class PairIntWithIndex extends com.climbwithyourfeet.clustering.util.PairInt {

        int pixIdx;

        public PairIntWithIndex(int xPoint, int yPoint, int thePixIndex) {
            super(xPoint, yPoint);
            pixIdx = thePixIndex;
        }

        @Override
        public boolean equals(Object obj) {

            if (!(obj instanceof com.climbwithyourfeet.clustering.util.PairInt)) {
                return false;
            }

            com.climbwithyourfeet.clustering.util.PairInt other
                = (com.climbwithyourfeet.clustering.util.PairInt) obj;

            return (x == other.getX()) && (y == other.getY());
        }

        @Override
        public int hashCode() {

            int hash = fnvHashCode(this.x, this.y);

            return hash;
        }

        @Override
        public com.climbwithyourfeet.clustering.util.PairInt copy() {
             return new PairIntWithIndex(x, y, pixIdx);
        }

        @Override
        public String toString() {

            StringBuilder sb = new StringBuilder(super.toString());
            sb.append(" pixIdx=").append(Integer.toString(pixIdx));

            return sb.toString();
        }

    }
    
}
