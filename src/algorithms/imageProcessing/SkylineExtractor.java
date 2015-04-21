package algorithms.imageProcessing;

import algorithms.CountingSort;
import algorithms.MultiArrayMergeSort;
import algorithms.compGeometry.EllipseHelper;
import algorithms.compGeometry.PerimeterFinder;
import algorithms.compGeometry.PointInPolygon;
import algorithms.imageProcessing.util.MatrixUtil;
import algorithms.misc.Histogram;
import algorithms.misc.HistogramHolder;
import algorithms.misc.MiscMath;
import algorithms.util.ArrayPair;
import algorithms.util.Errors;
import algorithms.util.LinearRegression;
import algorithms.util.PairFloat;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.PolynomialFitter;
import algorithms.util.ResourceFinder;
import algorithms.util.ScatterPointPlotterPNG;
import java.awt.Color;
import java.io.IOException;
import java.security.NoSuchAlgorithmException;
import java.util.*;
import java.util.logging.Logger;

public class SkylineExtractor {

    private Logger log = Logger.getLogger(this.getClass().getName());
     /**
     * using the gradient's theta image, find the sky as the largest set of
     * contiguous 0 values and apply the edge filter to it to reduce the
     * boundary to a single pixel curve.  
     * 
     * NOTE that the theta image has a boundary that has been increased by 
     * the original image blur and then the difference of gaussians to make 
     * the gradient, so the distance of the skyline from the real image horizon 
     * is several pixels.
     * For example, the canny edge filter used in "outdoorMode" results in
     * gaussian kernels applied twice to give an effective sigma of 
     * sqrt(2*2 + 0.5*0.5) = 2.1.  The FWHM of such a spread is then 
     * 2.355*2.1 = 6 pixels.   The theta image skyline is probably blurred to a 
     * width larger than the combined FWHM, however, making it 7 or 8 pixels.  
     * Therefore, it's recommended that the image returned from this be followed 
     * with: edge extraction; then fit the edges to the intermediate canny edge 
     * filter product (the output of the 2 layer filter) by making a translation
     * of the extracted skyline edge until the correlation with the filter2
     * image is highest (should be within 10 pixel shift).  Because this
     * method does not assume orientation of the image, the invoker needs
     * to also retrieve the centroid of the sky, so that is also returned 
     * in an output variable given in the arguments.
     * 
     * @param theta
     * @param gradientXY
     * @param originalImage
     * @param outputSkyCentroid container to hold the output centroid of 
     * the sky.
     * @param edgeSettings
     * @return
     * @throws IOException
     * @throws NoSuchAlgorithmException 
     */
    public GreyscaleImage createSkyline(GreyscaleImage theta, 
        GreyscaleImage gradientXY, ImageExt originalImage,
        CannyEdgeFilterSettings edgeSettings, PairIntArray outputSkyCentroid) 
        throws IOException, NoSuchAlgorithmException {        
      
        GreyscaleImage mask = createBestSkyMask(theta, gradientXY, originalImage, 
            edgeSettings, outputSkyCentroid);
        
        if (mask != null) {
            
            ImageProcessor imageProcessor = new ImageProcessor();
            
            imageProcessor.multiply(mask, 255);
            
            CannyEdgeFilter filter = new CannyEdgeFilter();
            
            filter.setFilterImageTrim(theta.getXRelativeOffset(), 
                theta.getYRelativeOffset(), theta.getWidth(), 
                theta.getHeight());
            
            filter.applyFilter(mask);
            
            return mask;
        }
        
        return null;
    }
    
    protected int determineBinFactorForSkyMask(int numberOfThetaPixels) {
        
        //TODO: this can be adjusted by the jvm settings for stack size
        int defaultLimit = 87000;
        
        if (numberOfThetaPixels <= defaultLimit) {
            return 1;
        }
                    
        double a = (double)numberOfThetaPixels/87000.;
        // rounds down
        int f2 = (int)a/2;
        int binFactor = f2 * 2;
        if ((a - binFactor) > 0) {
            binFactor += 2;
        }
            
        return binFactor;
    }
    
    /**
     * NOT READY FOR USE
     * 
     * create a mask for what is interpreted as sky in the image and return
     * a mask with 0's for sky and 1's for non-sky.
     * 
     * Internally, the method looks at contiguous regions of zero value pixels 
     * in the theta image and it looks at color in the original image.  
     * The camera image plane can have a rotation such that the 
     * horizon might not be along rows in the image, that is the method
     * looks for sky that is not necessarily at the top of the image, but 
     * should be on the boundary of the image
     * (note that reflection of sky can be found the same way, but this
     * method does not try to find sky that is not on the boundary of the
     * image).
     * (the "boundary" logic may change, in progress...)
     * 
     * Note that if the image contains a silhouette of featureless
     * foreground and a sky full of clouds, the method will interpret the
     * foreground as sky so it is up to the invoker to invert the mask.
     * 
     * NOTE: the cloud finding logic will currently fail if originalColorImage
     * is black and white.
     * 
     * @param theta
     * @param gradientXY
     * @param originalColorImage
     * @param edgeSettings
     * @param outputSkyCentroid container to hold the output centroid of 
     * the sky.
     * @return 
     * @throws java.io.IOException 
     * @throws java.security.NoSuchAlgorithmException 
     */
    public GreyscaleImage createBestSkyMask(final GreyscaleImage theta,
        GreyscaleImage gradientXY, ImageExt originalColorImage, 
        CannyEdgeFilterSettings edgeSettings, PairIntArray outputSkyCentroid) 
        throws IOException, NoSuchAlgorithmException {
        
        if (theta == null) {
            throw new IllegalArgumentException("theta cannot be null");
        }
        
        int binFactor = determineBinFactorForSkyMask(theta.getNPixels());

        log.info("binFactor=" + binFactor);
        
        Set<PairInt> points = new HashSet<PairInt>();
        
        /*TODO:  
        adjust this algorithm to allow alternate ways of determining the seed 
        sky points (and excluded points that should not be re-added).  
        Example case where this is needed here:
            an image with a foreground large smooth snow field and dark mountain 
            ranges under a cloudy structured sky.
            The edge detector w/o changes would probably find the mountain 
            ranges best in this case, but one would still not know "sky" without
            GPS or external sensors or assumption of horizontal.
            **The sun’s position however is learnable from some images and 
            that would help determine the location of seed sky points correctly.
            **The scattered and reflected light from the sun and water might be 
            hard to distinguish. 
            Will look for a pattern in solar light as a function of distance
            along the increasing brightness axis.
            Note that when the sun is near the horizon, the light reaching the
            camera has been scattered less and so one might be able to make 
            a simple model of scattering of solar light off of
            optically thick water w/ consideration for mie and rayleigh 
            scattering in the atmosphere altering the received source light.
        */
        
        RemovedSets removedSets = new RemovedSets();
        
        GreyscaleImage threshholdedGXY = filterAndExtractSkyFromGradient(
            originalColorImage, theta, gradientXY, binFactor, points,
            removedSets);
        
        if (removedSets == null) {
            
            GreyscaleImage mask = theta.createWithDimensions();
               
            // return an image of all 1's
            mask.fill(1);
            
            return mask;
        }
        
        //now the coordinates in zeroPointLists are w.r.t. thetaImg

        PerimeterFinder perimeterFinder = new PerimeterFinder();
        int[] skyRowMinMax = new int[2];
        Set<PairInt> outputEmbeddedGapPoints = new HashSet<PairInt>();
        Map<Integer, List<PairInt>> skyRowColRange = perimeterFinder.find(points, 
            skyRowMinMax, originalColorImage.getWidth(), outputEmbeddedGapPoints);
        
        rightAndLowerDownSizingSkyPointCorrections(points, binFactor, 
            skyRowColRange, skyRowMinMax, originalColorImage,
            theta.getWidth(), theta.getHeight(),
            theta.getXRelativeOffset(), theta.getYRelativeOffset());

debugPlot(points, originalColorImage, theta.getXRelativeOffset(), theta.getYRelativeOffset(), 
"after_downsize_corrections_2");

        GreyscaleImage mask = gradientXY.createWithDimensions();
        mask.fill(1);
        for (PairInt p : points) {
            int x = p.getX();
            int y = p.getY(); 
            mask.setValue(x, y, 0);
        }
                
        int nSkyPointsBeforeFindClouds = points.size();
        
        populatePixelExtColors(points, originalColorImage, mask);
            
        findClouds(points, new HashSet<PairInt>(), originalColorImage, mask);

try {
    String dirPath = ResourceFinder.findDirectory("bin");
    ImageExt clr = (ImageExt)originalColorImage.copyImage();
    ImageIOHelper.addToImage(points, mask.getXRelativeOffset(), mask.getYRelativeOffset(), clr);
    ImageIOHelper.writeOutputImage(
        dirPath + "/sky_after_find_clouds_" + outImgNum + ".png", clr);
    outImgNum++;
} catch (IOException e) {
    log.severe("ERROR: " + e.getMessage());
}        
debugPlot(points, originalColorImage, mask.getXRelativeOffset(), mask.getYRelativeOffset(), 
"before_add_embedded");

        GroupPixelColors allSkyColor = new GroupPixelColors(points,
            originalColorImage, theta.getXRelativeOffset(), 
            theta.getYRelativeOffset());
        
        boolean skyIsDarkGrey = skyIsDarkGrey(allSkyColor);

        Set<PairInt> sunPoints = new HashSet<PairInt>();
        double[] ellipFitParams = findSunConnectedToSkyPoints(points, 
            removedSets.getReflectedSunRemoved(), originalColorImage, 
            theta.getXRelativeOffset(), 
            theta.getYRelativeOffset(), skyIsDarkGrey, sunPoints);

        //TODO: adjust this:
        if ((ellipFitParams != null) && 
            (
            ((ellipFitParams[2]/ellipFitParams[3])> 7)
            //TODO: reconsider this rule for sun on edge of image
            || ((ellipFitParams[0] < 0) || (ellipFitParams[1] < 0))
            )
            ) {
            sunPoints.clear();
        }

        // should not see sun and rainbow in same image
        Set<PairInt> rainbowPoints = sunPoints.isEmpty() ?
            findRainbowPoints(points, 
                removedSets.getReflectedSunRemoved(), 
                originalColorImage, 
                theta.getXRelativeOffset(), theta.getYRelativeOffset(),
                skyIsDarkGrey) :
            new HashSet<PairInt>();
        
        // find remaining embedded points
        skyRowMinMax = new int[2];
        Set<PairInt> embeddedPoints = new HashSet<PairInt>();
        skyRowColRange = perimeterFinder.find(points, skyRowMinMax, 
            originalColorImage.getWidth(), embeddedPoints);
        
        if (!embeddedPoints.isEmpty()) {
debugPlot(points, originalColorImage, mask.getXRelativeOffset(), mask.getYRelativeOffset(), 
"before_add_embedded");

            HistogramHolder[] brightnessHistogram = new HistogramHolder[1];

            // brightest sky is in bin [2], and dimmest is in [0]
            GroupPixelColors[] skyPartitionedByBrightness = 
                partitionInto3ByBrightness(points, originalColorImage, 
                theta.getXRelativeOffset(), theta.getYRelativeOffset(), 
                brightnessHistogram);

            // no need to update rowColRanges as this was not an external "grow"
            //add embedded pixels if they're near existing sky colors
            addIfSimilarToSky(embeddedPoints, points, 
                removedSets.getHighContrastRemoved(),
                originalColorImage, mask,
                brightnessHistogram, skyPartitionedByBrightness);        

debugPlot(points, originalColorImage, mask.getXRelativeOffset(), mask.getYRelativeOffset(), 
"after_add_embedded");

        }

        Set<PairInt> exclude = new HashSet<PairInt>();
        exclude.addAll(removedSets.getHighContrastRemoved());
        exclude.addAll(rainbowPoints);
        exclude.addAll(sunPoints);

        log.info("number of sunPoints=" + sunPoints.size() + " "
        + "reflectedSunRemoved.size()=" + removedSets.getReflectedSunRemoved().size());

        if (sunPoints.isEmpty()) {
        
            growForLowContrastLimits(points, exclude, skyRowColRange,
                skyRowMinMax, originalColorImage, mask);

debugPlot(points, originalColorImage, mask.getXRelativeOffset(), mask.getYRelativeOffset(), 
"after_low_contrast_grow");

        }

        if (!sunPoints.isEmpty()) {
            correctSkylineForSun(sunPoints, points, originalColorImage, mask, gradientXY);
        }/* else if (!rainbowPoints.isEmpty()) {
            correctSkylineForRainbow(rainbowPoints, points, colorImg, mask);
        }*/

debugPlot(points, originalColorImage, mask.getXRelativeOffset(), 
    mask.getYRelativeOffset(), "final");
        
        for (PairInt p : sunPoints) {
            int x = p.getX();
            int y = p.getY();            
            mask.setValue(x, y, 0);
        }
        
        points.addAll(sunPoints);
        
        for (PairInt p : rainbowPoints) {
            int x = p.getX();
            int y = p.getY();            
            mask.setValue(x, y, 0);
        }
        points.addAll(rainbowPoints);
        
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        
        double[] xycen = curveHelper.calculateXYCentroids(points);
 
        outputSkyCentroid.add((int)Math.round(xycen[0]), (int)Math.round(xycen[1]));

        ImageProcessor imageProcessor = new ImageProcessor();
        imageProcessor.removeSpurs(mask);
   
        return mask;
    }
    
    public Set<PairInt> combine(List<PairIntArray> points) {
        Set<PairInt> set = new HashSet<PairInt>();
        combine(points, set);
        return set;
    }
    public void combine(List<PairIntArray> points, Set<PairInt> outputCombined) {
        for (PairIntArray p : points) {
            for (int i = 0; i < p.getN(); i++) {
                int x = p.getX(i);
                int y = p.getY(i);
                PairInt pi = new PairInt(x, y);
                outputCombined.add(pi);
            }
        }
    }
    
    public List<PairIntArray> getLargestSortedContiguousZeros(GreyscaleImage theta) {
                 
        return getSortedContiguousValues(theta, 0, false, true);
    }
    
    public List<PairIntArray> getSortedContiguousZeros(GreyscaleImage theta) {
                 
        return getSortedContiguousValues(theta, 0, false, false);
    }
    
    public List<PairIntArray> getLargestSortedContiguousNonZeros(GreyscaleImage theta) {
                 
        return getSortedContiguousValues(theta, 0, true, true);
    }
    
    private List<PairIntArray> getSortedContiguousValues(GreyscaleImage theta,
        int value, boolean excludeValue, boolean limitToLargest) {
        
        DFSContiguousValueFinder zerosFinder = new DFSContiguousValueFinder(theta);
        
        if (excludeValue) {
            zerosFinder.findGroupsNotThisValue(value);
        } else {
            zerosFinder.findGroups(value);
        }
        
        int nGroups = zerosFinder.getNumberOfGroups();
        
        if (nGroups == 0) {
            return new ArrayList<PairIntArray>();
        }
        // ====== find the group(s) with the largest number of zero pixels =====
        
        int nMaxGroupN = Integer.MIN_VALUE;
        int[] groupIndexes = new int[nGroups];
        int[] groupN = new int[nGroups];
        for (int gId = 0; gId < nGroups; gId++) {
            int n = zerosFinder.getNumberofGroupMembers(gId);
            groupIndexes[gId] = gId;
            groupN[gId] = n;
            if (n > nMaxGroupN) {
                nMaxGroupN = n;
            }
        }
        
        int maxValue = MiscMath.findMax(groupN);
        if ((maxValue > groupN.length) || (nMaxGroupN > 10000000)) {
            MultiArrayMergeSort.sortByDecr(groupN, groupIndexes);
        } else {
            CountingSort.sortByDecr(groupN, groupIndexes, maxValue);
        }
        
        List<Integer> groupIds = new ArrayList<Integer>();
        groupIds.add(Integer.valueOf(groupIndexes[0]));
        
        if (nGroups > 1) {
                        
            float n0 = (float)groupN[0];
            
            for (int i = 1; i < groupN.length; i++) {
                
                if (limitToLargest) {
                    float number = groupN[i];

                    float frac = number/n0;
                    //TODO: this should be adjusted by some metric.
                    //      a histogram?
                    // since most images should have been binned to <= 300 x 300 pix,
                    // making an assumption about a group >= 100 pixels 
                    if ((1 - frac) < 0.4) {
                    //if (number > 100) {
                        groupIds.add(Integer.valueOf(groupIndexes[i]));
                    } else {
                        break;
                    }
                } else {
                    groupIds.add(Integer.valueOf(groupIndexes[i]));
                }
            }
        }

        List<PairIntArray> list = new ArrayList<PairIntArray>();
        
        for (Integer gIndex : groupIds) {
            
            int gIdx = gIndex.intValue();
            
            PairIntArray points = zerosFinder.getXY(gIdx);
            
            list.add(points);
        }
        
        return list;
    }
    
    private void reduceToLargest(List<PairIntArray> zeroPointLists) {
        
        int rmIdx = -1;
        
        if (zeroPointLists.size() > 1) {
                        
            float n0 = (float)zeroPointLists.get(0).getN();
            
            for (int i = 1; i < zeroPointLists.size(); i++) {
                
                float number = zeroPointLists.get(i).getN();

                float frac = number/n0;
                //TODO: this should be adjusted by some metric.
                //      a histogram?
                // since most images should have been binned to <= 300 x 300 pix,
                // making an assumption about a group >= 100 pixels 
                if (frac < 0.1) {
                    rmIdx = i;
                    break;
                }
            }
        }
        
        if (rmIdx > -1) {
            List<PairIntArray> out = new ArrayList<PairIntArray>();
            for (int i = 0; i < rmIdx; i++) {
                out.add(zeroPointLists.get(i));
            }
            zeroPointLists.clear();
            zeroPointLists.addAll(out);
        }
    }
    
    /**
     * remove high contrast points from the sky points.  this helps to remove
     * points that are present due to "blind spots" in gradientXY on the scale
     * of the combined convolution of gaussians that created the gradientXY.
     * For example, repetitive structure like skyscraper windows are objects
     * in the color image which may be missing an outline in the gradientXY.
     * These features have higher contrast in the color image than the normal
     * sky because they are objects, not sky, so this method tries to find
     * those and remove them from the sky points.  Note that the method prefers
     * to err on the side of over subtracting because later steps can find 
     * and re-include any connected sky as long as the majority of sky points 
     * remain at the end of this method.
     * 
     * @param zeroPointLists
     * @param originalColorImage
     * @param theta
     * @param avgY
     * @param addAlongX
     * @param addAmount 
     * @param outputRemovedPoints is populated with removed points if 
     * the object is not null
     */
    private void removeHighContrastPoints(List<PairIntArray> 
        zeroPointLists, Image originalColorImage, GreyscaleImage theta,
        double avgY, Set<PairInt> outputRemovedPoints) {
        
        int xOffset = theta.getXRelativeOffset();
        int yOffset = theta.getYRelativeOffset();
        
        // remove points that have contrast larger than tail of histogram
        HistogramHolder h = createContrastHistogram(avgY, zeroPointLists, 
            originalColorImage, xOffset, yOffset);
        
        if (h == null) {
            return;
        }

try {
    h.plotHistogram("contrast", 1);
    // printed as bin/classes/points_and_polygon1.html
} catch (IOException e) {
    log.severe(e.getMessage());
}
        List<Integer> strongPeaks = MiscMath.findStrongPeakIndexes(h, 0.1f);
         
        if (strongPeaks == null || strongPeaks.isEmpty()) {
            return;
        }
       
        int lastPeakIdx = strongPeaks.get(strongPeaks.size() - 1).intValue();
        
        //TODO: this is sensitive to the histogram formation so needs a wide
        // variety of data for testing to make sure it finds the right characteristic
        
        int yPeakIdx = lastPeakIdx;
        int tailXIdx = h.getXHist().length - 1;
        if (tailXIdx > yPeakIdx) {
            float yPeak =  h.getYHist()[yPeakIdx];
            float crit = 0.03f;
            float dy = Float.MIN_VALUE;
            for (int i = (yPeakIdx + 1); i < h.getYHist().length; i++) {
                
                float f = (float)h.getYHist()[i]/yPeak;
                dy = Math.abs(h.getYHist()[i] - h.getYHist()[i - 1]);
                
                //System.out.println("x=" + h.getXHist()[i] + " f=" + f + " dy=" + dy);
                
                if (f < crit) {
                    tailXIdx = i;
                    break;
                }
            }
        }
        
        double[][] m = new double[3][];
        m[0] = new double[]{0.256, 0.504, 0.098};
        m[1] = new double[]{-0.148, -0.291, 0.439};
        m[2] = new double[]{0.439, -0.368, -0.072};
        
        //remove points w/ contrast higher than the tail of the histogram
        double critContrast = h.getXHist()[tailXIdx];
             
        for (int gId = 0; gId < zeroPointLists.size(); gId++) {
                        
            PairIntArray points  = zeroPointLists.get(gId);
                        
            Set<PairInt> pointsSet = new HashSet<PairInt>();
            for (int i = 0; i < points.getN(); i++) {
                int x = points.getX(i);
                int y = points.getY(i);
                PairInt pi = new PairInt(x, y);
                pointsSet.add(pi);
            }

            for (int i = 0; i < points.getN(); i++) {
                
                int x = points.getX(i);
                int y = points.getY(i);
                
                int ox = x + xOffset;
                int oy = y + yOffset;

                if ((ox < 0) || (ox > (originalColorImage.getWidth() - 1))) {
                    continue;
                }
                if ((oy < 0) || (oy > (originalColorImage.getHeight() - 1))) {
                    continue;
                }
        
                int r = originalColorImage.getR(ox, oy);
                int g = originalColorImage.getG(ox, oy);
                int b = originalColorImage.getB(ox, oy);
                
                double[] rgb = new double[]{r, g, b};
                        
                double[] yuv = MatrixUtil.multiply(m, rgb);
                yuv = MatrixUtil.add(yuv, new double[]{16, 128, 128});

                float contrast = (float)((avgY - yuv[0]) / yuv[0]);
                
                if (contrast > critContrast) {
                    
                    PairInt pi0 = new PairInt(x, y);
                    
                    pointsSet.remove(pi0);
                    
                    if (outputRemovedPoints != null) {
                        outputRemovedPoints.add(pi0);
                    }

                    for (int xx = (x - 1); xx <= (x + 1); xx++) {
                        if ((xx < 0) || (xx > (theta.getWidth() - 1))) {
                            continue;
                        }
                        for (int yy = (y - 1); yy <= (y + 1); yy++) {
                            if ((yy < 0) || (yy > (theta.getHeight() - 1))) {
                                continue;
                            }
                            if ((xx == x) && (yy == y)) {
                                continue;
                            }
                            
                            PairInt pi1 = new PairInt(xx, yy);
                            
                            if (pointsSet.contains(pi1)) {
                                pointsSet.remove(pi1);
                                if (outputRemovedPoints != null) {
                                    outputRemovedPoints.add(pi1);
                                }
                            }
                        }
                    }
                }
            }

            if (pointsSet.size() != points.getN()) {
                PairIntArray points2 = new PairIntArray();
                for (PairInt pi : pointsSet) {
                    points2.add(pi.getX(), pi.getY());
                }
                points.swapContents(points2);
            }
        }
        
        // remove empty sets
        for (int i = (zeroPointLists.size() - 1); i > -1; i--) {
            PairIntArray point = zeroPointLists.get(i);
            if (point.getN() == 0) {
                zeroPointLists.remove(i);
            }
        }
        
debugPlot(outputRemovedPoints, originalColorImage, xOffset, yOffset, "filtered_out_high_contrast");

    }

     
    /**
     * if the range of contrast is large, return a contrast histogram,
     * else return null.  
     * TODO: refactor to move the logic to return null to the invoker. 
     * For now, the only use of this method is simpler if it does not
     * return a histogram when it won't be needed.  definitely should
     * be refactored...
     * 
     * @param avgY
     * @param zeroPointLists
     * @param originalColorImage
     * @param xOffset
     * @param yOffset
     * @return 
     */
    private HistogramHolder createContrastHistogram(double avgY,
        List<PairIntArray> zeroPointLists, Image originalColorImage,
        int xOffset, int yOffset) {
        
        if (zeroPointLists.isEmpty()) {
            return null;
        }
        
        double[][] m = new double[3][];
        m[0] = new double[]{0.256, 0.504, 0.098};
        m[1] = new double[]{-0.148, -0.291, 0.439};
        m[2] = new double[]{0.439, -0.368, -0.072};
                
        int nPoints = 0;
        for (int gId = 0; gId < zeroPointLists.size(); gId++) {
            nPoints += zeroPointLists.get(gId).getN();
        }
        
        float[] yValues = new float[nPoints];
        
        int count = 0;
        for (int gId = 0; gId < zeroPointLists.size(); gId++) {
            
            PairIntArray points  = zeroPointLists.get(gId);
                        
            for (int i = 0; i < points.getN(); i++) {
                
                int x = points.getX(i);
                int y = points.getY(i);
         
                int ox = x + xOffset;
                int oy = y + yOffset;

                if ((ox < 0) || (ox > (originalColorImage.getWidth() - 1))) {
                    continue;
                }
                if ((oy < 0) || (oy > (originalColorImage.getHeight() - 1))) {
                    continue;
                }
                
                int r = originalColorImage.getR(ox, oy);
                int g = originalColorImage.getG(ox, oy);
                int b = originalColorImage.getB(ox, oy);
                
                double[] rgb = new double[]{r, g, b};
                        
                double[] yuv = MatrixUtil.multiply(m, rgb);
                yuv = MatrixUtil.add(yuv, new double[]{16, 128, 128});

                float contrastValue = (float)((avgY - yuv[0]) / yuv[0]);
                
                yValues[count] = contrastValue;
                
                count++;
            }
        }
        
        float[] yErr = Errors.populateYErrorsBySqrt(yValues);
        HistogramHolder h = Histogram.createSimpleHistogram(yValues, 
            yErr);

        if (h != null) {
            
            int lastZeroIdx = MiscMath.findLastZeroIndex(h);
            
            int nBins = h.getXHist().length;
            
            int nLastZeros = nBins - lastZeroIdx;
            
            if ((nLastZeros > 4) && (lastZeroIdx > -1)) {
                
                float halfBinWidth = (h.getXHist()[1] - h.getXHist()[1]) / 2.f;

                nBins = lastZeroIdx - 1;
                float xMin = h.getXHist()[0] - halfBinWidth;
                float xMax = h.getXHist()[lastZeroIdx - 1] + halfBinWidth;

                float contrastRange = h.getXHist()[lastZeroIdx - 1]
                    - h.getXHist()[0];

                if (contrastRange > 1.5) {

                    h = Histogram.calculateSturgesHistogramRemoveZeroTail(yValues, yErr);
                    
                    if ((h != null) && (h.getXHist().length == 1)) {
                        h = Histogram.calculateSturgesHistogram(xMin, xMax, yValues, yErr);
                    }
                    
                    return h;

                }
            }
        }
        
        return null;
    }
    
    private void transformPointsToOriginalReferenceFrame(Set<PairInt> points,
        GreyscaleImage theta, boolean makeCorrectionsAlongX, int addAmount) {
        
         // transform points to original color image frame
        int totalXOffset = theta.getXRelativeOffset();
        int totalYOffset = theta.getYRelativeOffset();

        if (makeCorrectionsAlongX) {
            totalXOffset += addAmount;
        } else {
            totalYOffset += addAmount;
        }
        
        for (PairInt p : points) {
            int x = p.getX();
            int y = p.getY();
            x += totalXOffset;
            y += totalYOffset;
            p.setX(x);
            p.setY(y);
        }
    }

   /**
     * using adaptive "thresholding" to subtract intensity levels from
     * gradientXY, find the contiguous zero values connected to skyPoints
     * and add them to skyPoints.
     * 
     * @param gradientXY
     * @param skyPoints
     * @param excludeThesePoints
     * @return the threshold subtracted gradient image
     */
    public GreyscaleImage extractSkyFromGradientXY(GreyscaleImage gradientXY,
        Set<PairInt> skyPoints, Set<PairInt> excludeThesePoints) {
                
        GreyscaleImage gXY2 = gradientXY.copyImage();
                
        // x is pixelValue , y is number of pixels holding that value
        // last in array is for the smallest pixelValue
        PairIntArray gXYValues = Histogram.createADescendingSortByKeyArray(gXY2);

        float sumF = 0;
        int valueFor0Point9 = 0;
        for (int i = (gXYValues.getN() - 1); i > -1; i--) {
            float f = (float)gXYValues.getY(i)/(float)gradientXY.getNPixels();
            sumF += f;
            if (sumF > 0.9) {
                break;
            }
            valueFor0Point9 = gXYValues.getX(i);
        }
/*        
float sumFrac = 0;
StringBuilder sb = new StringBuilder("gXY:\n");
for (int i = (gXYValues.getN() - 1); i > -1; i--) {
int sumToHighValues = 0;
for (int ii = i; ii > -1; ii--) {
    sumToHighValues += gXYValues.getY(ii);
}
float frac = (float)gXYValues.getY(i)/(float)gradientXY.getNPixels();
sumFrac += frac;
sb.append(String.format(" value=%d count=%d  f=%f  sumToEnd=%d", 
gXYValues.getX(i), gXYValues.getY(i), frac, sumToHighValues));
sb.append("sumF=").append(Float.toString(sumFrac));
sb.append("\n");
}
log.info(sb.toString());
*/
        int subtract = 0;
        int lastHistIdx = gXYValues.getN();
        
        float c0 = (float)gXYValues.getY(gXYValues.getN() - 1)/(float)gradientXY.getNPixels();
        float c1 = (float)gXYValues.getY(gXYValues.getN() - 2)/(float)gradientXY.getNPixels();
        //float c01 = c1 + c0;
        float cm01 = c0 - c1;
        log.info("valueFor0Point9=" + valueFor0Point9 + " cm01=" + cm01);
        if ((cm01 < 0.0) && (valueFor0Point9 <= 2)) {
            subtract = 0;
            lastHistIdx = gXYValues.getN();
        } else if (cm01 < 0.3) {
            if (valueFor0Point9 <= 2) {
                subtract = 2;
                lastHistIdx = gXYValues.getN() - 2;
            } /*else {
                subtract = 1;
                lastHistIdx = gXYValues.getN() - 1;
            }*/
        } else if ((cm01 >= 0.3) && (cm01 <= 0.45) && (valueFor0Point9 <= 2)) {
            subtract = 1;
            lastHistIdx = gXYValues.getN() - 1;
        }
        
        int nIter = 0;
        
        float originalMaxValue = gXYValues.getY(gXYValues.getN() - 1);
          
        //TODO: the logic changed over time here to revisit this code
        // and remove iteration.
        
        while ((subtract < originalMaxValue) && (nIter == 0)) {
            
            if (nIter > 0) {
                gXY2 = gradientXY.copyImage();
            }
            
            lastHistIdx--;
            if (lastHistIdx < 1) {
                break;
            }
            subtract = gXYValues.getX(lastHistIdx);

            subtractWithCorrectForNegative(gXY2, subtract);
           
            // ==== find contiguous zeros =====  
            
            growZeroValuePoints(skyPoints, excludeThesePoints, gXY2);
                                        
            log.info("nIter=" + nIter + ")" 
                + " out of " + skyPoints.size()
                + " (level=" + ((float)gXYValues.getY(lastHistIdx)/originalMaxValue) 
                + " subtract=" + subtract + " out of max=" + gXYValues.getX(0)
                + ")"
            );
                  
            nIter++;
        }
        
/*
try {
    Image img1 = gXY2.copyImageToGreen();
    ImageIOHelper.addToImage(skyPoints, 0, 0, img1);
    ImageDisplayer.displayImage("sky points subtract=" + subtract, img1);
} catch (IOException ex) {
    log.severe(ex.getMessage());
}
*/         
        return gXY2;
    }

    private void growZeroValuePoints(Set<PairInt> points, 
        Set<PairInt> excludeThesePoints, GreyscaleImage gradientXY) {
  
        if (points.isEmpty()) {
            return;
        }
        
        java.util.Stack<PairInt> stack = new java.util.Stack<PairInt>();
        
        //O(N_sky)
        for (PairInt p : points) {
            if (!excludeThesePoints.contains(p)) {
                stack.add(p);
            }
        }
        
        // null = unvisited, presence = visited
        Set<PairInt> visited = new HashSet<PairInt>();
        visited.add(stack.peek());
        
        int width = gradientXY.getWidth();
        int height = gradientXY.getHeight();
        
        while (!stack.isEmpty()) {

            PairInt uPoint = stack.pop();
            
            int uX = uPoint.getX();
            int uY = uPoint.getY();

            //(1 + frac)*O(N) where frac is the fraction added back to stack
            for (int vX = (uX - 1); vX <= (uX + 1); vX++) {
                
                if ((vX < 0) || (vX > (width - 1))) {
                    continue;
                }
                
                for (int vY = (uY - 1); vY <= (uY + 1); vY++) {
                    
                    if ((vY < 0) || (vY > (height - 1))) {
                        continue;
                    }
                    
                    PairInt vPoint = new PairInt(vX, vY);
                    
                    if (vPoint.equals(uPoint)) {
                        continue;
                    }
                
                    if (excludeThesePoints.contains(vPoint)) {
                        continue;
                    }
                    
                    if (visited.contains(vPoint)) {
                        continue;
                    }
                                      
                    if ((vX < 0) || (vX > (width - 1))) {
                        continue;
                    }
                    if ((vY < 0) || (vY > (height - 1))) {
                        continue;
                    }
                    
                    visited.add(vPoint);
                    
                    int v = gradientXY.getValue(vX, vY);
                                        
                    if (v == 0) {
                        
                        stack.add(vPoint);
                        
                        if (!points.contains(vPoint)) {
                            points.add(vPoint);
                        }
                    }
                }
            }
        }
    }
    
    boolean isPerimeterUnbound(Map<Integer, PairInt> gRowColRange, 
        int[] gRowMinMax, Map<Integer, PairInt> boundingRowColRange, 
        int[] boundingRowMinMax,
        int xMinImage, int xMaxImage, int yMinImage, int yMaxImage) {
        
        // check top and bottom rows of group are within bounds
        if (gRowMinMax[0] == boundingRowMinMax[0]) {
            if (gRowMinMax[0] != yMinImage) {
                return true;
            }
        } else if (gRowMinMax[0] < boundingRowMinMax[0]) {
            return true;
        }
        if (gRowMinMax[1] == boundingRowMinMax[1]) {
            if (gRowMinMax[1] != yMaxImage) {
                return true;
            }
        } else if (gRowMinMax[1] > boundingRowMinMax[1]) {
            return true;
        }
        
        for (int r = gRowMinMax[0]; r <= gRowMinMax[1]; r++) {
                    
            PairInt cRange = gRowColRange.get(Integer.valueOf(r));
            
            // see if each point in cRange is on the boundary of rowColRange
            //   or within it
           
            // check left half            
            int x = cRange.getX();
            int y = r;
            
            PairInt cRange2 = boundingRowColRange.get(Integer.valueOf(y));
            if (cRange2 == null) {
                return true;
            }
            int x2 = cRange2.getX();
            
            if (x == x2) {
                if (x > xMinImage) {
                    return true;
                }
            } else if (x < x2) {
                // x is outside of the larger region defined by gRowColRange
                return true;
            }
            
            // check right half
            x2 = cRange2.getY();
            
            if (x == x2) {
                if (x < xMaxImage) {
                    return true;
                }
            } else if (x > x2) {
                // x is outside of the larger region defined by gRowColRange
                return true;
            }
        }
        
        return false;
    }
    
    private boolean isPerimeterUnbound(Map<Integer, PairInt> gRowColRange, 
        int[] gRowMinMax, Set<PairInt> skyPoints, double[] groupXYCen,
        int imageWidth, int imageHeight) {
        
        boolean unbounded = false;
        
        for (int r = gRowMinMax[0]; r <= gRowMinMax[1]; r++) {
                    
            PairInt cRange = gRowColRange.get(Integer.valueOf(r));
            
            for (int k = 0; k < 2; k++) {
                int c;
                switch(k) {
                    case 0:
                        c = cRange.getX();
                        break;
                    default:
                        c = cRange.getY();
                        break;
                }
                
                if (c < groupXYCen[0]) {
                
                    // look for points to left
                    int xt = c - 1;
                    if (xt < 0) {
                        // bounded by edge of image
                        continue;
                    }
                
                    if (r < groupXYCen[1]) {
                
                        //look for points to left and top (=lower y)                            
                        int yt = r;
                        PairInt p = new PairInt(xt, yt);
                        if (!skyPoints.contains(p)) {
                            // not bounded on left
                            unbounded = true;
                            break;
                        }
                        //found a sky point to the left
                        yt--;
                        if (yt < 0) {
                            // bounded by edge of image
                            continue;
                        } else {
                            p = new PairInt(xt, yt);
                            if (!skyPoints.contains(p)) {
                                // not bounded on left
                                unbounded = true;
                                break;
                            }
                        }
                        
                    } else {
                        
                        //look for bounding points to left, bottom (=higher y)
                        int yt = r;
                        PairInt p = new PairInt(xt, yt);
                        if (!skyPoints.contains(p)) {
                            // not bounded on left
                            unbounded = true;
                            break;
                        }
                        yt++;
                        if (yt > (imageHeight - 1)) {
                            // bounded by edge of image
                            continue;
                        } else {
                            p = new PairInt(xt, yt);
                            if (!skyPoints.contains(p)) {
                                // not bounded on left
                                unbounded = true;
                                break;
                            }
                        }
                    }
                
                } else {

                    // look for points to the right
                    int xt = c + 1;
                    if (xt > (imageWidth - 1)) {
                        // bounded by edge of image
                        continue;
                    }
                
                    if (r < groupXYCen[1]) {

                        //look for bounding points to right, top (=lower y),

                        int yt = r;
                        PairInt p = new PairInt(xt, yt);
                        if (!skyPoints.contains(p)) {
                            // not bounded on left
                            unbounded = true;
                            break;
                        }
                        yt--;
                        if (yt < 0) {
                            // bounded by edge of image
                            continue;
                        } else {
                            p = new PairInt(xt, yt);
                            if (!skyPoints.contains(p)) {
                                // not bounded on left
                                unbounded = true;
                                break;
                            }
                        }

                    } else {
                    
                        //look for bounding points to right, bottom (=higher y)

                        int yt = r;
                        PairInt p = new PairInt(xt, yt);
                        if (!skyPoints.contains(p)) {
                            // not bounded on left
                            unbounded = true;
                            break;
                        }
                        yt++;
                        if (yt > (imageHeight - 1)) {
                            // bounded by edge of image
                            continue;
                        } else {
                            p = new PairInt(xt, yt);
                            if (!skyPoints.contains(p)) {
                                // not bounded on left
                                unbounded = true;
                                break;
                            }
                        }
                    }
                }
            }
            
            if (unbounded) {
                break;
            }
        }
        
        return unbounded;
    }

    /**
     * attempt to find within pixels connected to skyPoints, pixels that
     * look like sun pixels by color (hsb) and whose x,y distribution
     * resemble and ellipse (circle w/ possible occlusion).  
     * those sun points are then added to the skyPoints.
     * Note that if the sun is present in sky and in reflection, such as
     * water, their location in x, y must be fittable by an ellipse, else they 
     * may not be found as sun points.
     * @param skyPoints
     * @param reflectedSunRemoved previously removed points assumed to be
     * reflected light
     * @param clr
     * @param xOffset
     * @param yOffset 
     * @param skyIsDarkGrey 
     * @param outputYellowPoints 
     * @return ellipse fitting parameters if there are sun points
     *  [xc, yc, a, b, alpha]
     */
    protected double[] findSunConnectedToSkyPoints(Set<PairInt> skyPoints, 
        Set<PairInt> reflectedSunRemoved,
        ImageExt clr, int xOffset, int yOffset, boolean skyIsDarkGrey,
        Set<PairInt> outputYellowPoints) {
        
        Set<PairInt> yellowPoints = new HashSet<PairInt>();

        java.util.Stack<PairInt> yellowStack = new java.util.Stack<PairInt>();
     
        SunColors sunColors = new SunColors();
        
        if (skyIsDarkGrey) {
        
            sunColors.useDarkSkiesLogic();
        }
                
        int width = clr.getWidth();
        int height = clr.getHeight();
        
        int[] dxs = new int[]{-1, -1,  0,  1, 1, 1, 0, -1};
        int[] dys = new int[]{ 0, -1, -1, -1, 0, 1, 1,  1};
       
        for (PairInt p : skyPoints) {

            int x = p.getX() + xOffset;
            int y = p.getY() + yOffset;
            
            for (int k = 0; k < dxs.length; k++) {
            
                // xx,yy are w.r.t full color img
                int xx = x + dxs[k];
                int yy = y + dys[k];
                
                if ((xx < 0) || (xx > (width - 1)) || (yy < 0) || 
                    (yy > (height - 1))) {
                    continue;
                }
            
                PairInt p2 = new PairInt(xx - xOffset, yy - yOffset);

                if (yellowPoints.contains(p2) || reflectedSunRemoved.contains(p2)) {
                    continue;
                }
                
                int vIdx = clr.getInternalIndex(xx, yy);
                
                if (sunColors.isSunCenterColor(clr, vIdx)) {
                    
                    yellowPoints.add(p2);

                    yellowStack.add(p2);
                    
                }
            }
        }
        
        log.info("found " + yellowPoints.size() + " yellow points in connected search");
        
        if (yellowStack.size() < 3) {
            return null;
        }
        
        Set<PairInt> visited = new HashSet<PairInt>();
        visited.add(yellowStack.peek());
       
        while (!yellowStack.isEmpty()) {

            PairInt uPoint = yellowStack.pop();
            
            int uX = uPoint.getX() + xOffset;
            int uY = uPoint.getY() + yOffset;

            //(1 + frac)*O(N) where frac is the fraction added back to stack
            
            for (int k = 0; k < dxs.length; k++) {
                    
                int dx = dxs[k];
                int dy = dys[k];
                
                int vX = uX + dx;
                int vY = uY + dy;
                
                if ((vX < 0) || (vX > (width - 1)) || (vY < 0) || 
                    (vY > (height - 1))) {
                    continue;
                }
            
                PairInt vPoint = new PairInt(vX - xOffset, vY - yOffset);

                //TODO: consider whether sun would already be present in skyPoints 
                // and if so, a correction may be needed in invoking code
                
                // skypoints has already been check for same color criteria
                // so no need to check again here
                if (visited.contains(vPoint) || 
                    reflectedSunRemoved.contains(vPoint)
                    || yellowPoints.contains(vPoint)) {
                    
                    continue;
                }

                visited.add(vPoint);

                int vIdx = clr.getInternalIndex(vX, vY);

                if (sunColors.isSunCenterColor(clr, vIdx)) {
                    
                    yellowPoints.add(vPoint);

                    yellowStack.add(vPoint);
                    
                }
            }
        }
        
        log.info("found " + yellowPoints.size() + " points");

        if (yellowPoints.size() < 6) {
            return null;
        }
        
        //fit ellipse to yellowPoints.  ellipse because of possible occlusion.
        EllipseHelper ellipseHelper = new EllipseHelper();
        double[] params = ellipseHelper.fitEllipseToPoints(yellowPoints);
        
        if (params == null) {
            // not close to an ellipse
            return null;
        }
        
        float xc = (float)params[0];
        float yc = (float)params[1];
        float a = (float)params[2];
        float b = (float)params[3];
        float alpha = (float)params[4];
        
        log.info("elliptical fit to yellow points: " + Arrays.toString(params));
       
        //TODO:  this may need adjustements.  top of sun rising over mountains..
        /*if ((a/b) > 6) {
            return new HashSet<PairInt>();
        }*/
        
        outputYellowPoints.addAll(yellowPoints);
        
        return params;
    }
    
    /**
     * search for rainbow colored points over the entire image then fits an
     * ellipse to them and asserts that the points have certain colors in
     * them.  If the original fit to an ellipse is not good, the
     * method divides the rainbow points by contiguous subsamples to find best
     * and similar fits.  The last step of color requirement helps to rule
     * out half ellipse patterns in rocks for instance that have only rock
     * colors in them. 
     * 
     * @param skyPoints
     * @param reflectedSunRemoved
     * @param colorImg
     * @param xRelativeOffset
     * @param yRelativeOffset
     * @return 
     */
    Set<PairInt> findRainbowPoints(Set<PairInt> skyPoints, 
        Set<PairInt> reflectedSunRemoved,
        ImageExt colorImg, int xOffset, int yOffset, boolean skyIsDarkGrey) {

        Set<PairInt> rainbowPoints = findRainbowColoredPoints(colorImg, 
            reflectedSunRemoved, xOffset, yOffset, skyIsDarkGrey);
        
        if (rainbowPoints.isEmpty()) {
            return rainbowPoints;
        }
        
        if (rainbowPoints.size() < 12) {
            return new HashSet<PairInt>();
        }
        
        // fit a polynomial to rainbow points.  
        // would prefer a circle, but the optical depth of the dispersers and the
        // orientation of groups of them is not always a slab perpendicular to 
        // the camera

        int[] minMaxXY = MiscMath.findMinMaxXY(rainbowPoints);
        log.info("rainbow range in x: " + minMaxXY[0] + " to " + minMaxXY[1]);
        
        PolynomialFitter polyFitter = new PolynomialFitter();
        //y = c0*1 + c1*x[i] + c2*x[i]*x[i]
        float[] coef = polyFitter.solveAfterRandomSampling(rainbowPoints);
        
        if (coef == null) {
            return new HashSet<PairInt>();
        }
        
        log.info("rainbow polynomial coefficients = " + Arrays.toString(coef));
        log.info("image dimensions are " + colorImg.getWidth() + " X " + 
            colorImg.getHeight() + " pixels^2");
         
        polyFitter.plotFit(coef, rainbowPoints, colorImg.getWidth(),
            colorImg.getHeight(), 23, "rainbow points");
        
        double resid = polyFitter.calcResiduals(coef, rainbowPoints);

        //TODO: determine this more accurately:
        if (resid > (rainbowPoints.size() * 5)) {
            
            Set<PairInt> bestFittingPoints = new HashSet<PairInt>();
            
            coef = polyFitter.solveForBestFittingContiguousSubSamples(
                rainbowPoints, bestFittingPoints, colorImg.getWidth(), 
                colorImg.getHeight());
            
            if (coef == null) {
                return new HashSet<PairInt>();
            }
            
            int w = colorImg.getWidth() - xOffset;
            int h = colorImg.getHeight() - yOffset;
            if (bestFittingPoints.size() > (0.3f * (float)(w*h))) {
                return new HashSet<PairInt>();
            }
 
            // assert that colors are as expected, that is, that we don't
            // have only green and blue from rocks
            
            // filter to keep only those with a significant fraction with 
            //    cieX > 0.4 and cieY < 0.3
            // filter to keep only those with red
            
            int nGTX = 0;
            int nLTY = 0;
            int nPurpRed = 0;
            int nOranRed = 0;
            int nYellow = 0;
            CIEChromaticity cieC = new CIEChromaticity();
            ArrayPair cPurpRed = cieC.getRedThroughPurplishRedPolynomial();
            ArrayPair cOranRed = cieC.getOrangePolynomial();
            ArrayPair cYellow = cieC.getYellowishGreenThroughOrangePolynomial();
            PointInPolygon pInPoly = new PointInPolygon();
            
            float minCIEX = Float.MAX_VALUE;
            float maxCIEX = Float.MIN_VALUE;
            float minCIEY = Float.MAX_VALUE;
            float maxCIEY = Float.MIN_VALUE;
            for (PairInt p : bestFittingPoints) {
                int x = p.getX();
                int y = p.getY();
                int r = colorImg.getR(x + xOffset, y + yOffset);
                int g = colorImg.getG(x + xOffset, y + yOffset);
                int b = colorImg.getB(x + xOffset, y + yOffset);
                float[] hsb = new float[3];
                Color.RGBtoHSB(r, g, b, hsb);
        
                float[] cie = cieC.rgbToXYChromaticity(r, g, b);
                if (cie[0] >= 0.39) {
                    nGTX++;
                }
                if (cie[1] <= 0.3) {
                    nLTY++;
                }
                if (cie[0] < minCIEX) {
                    minCIEX = cie[0];
                }
                if (cie[0] > maxCIEX) {
                    maxCIEX = cie[0];
                }
                if (cie[1] < minCIEY) {
                    minCIEY = cie[1];
                }
                if (cie[1] > maxCIEY) {
                    maxCIEY = cie[1];
                }
                
                if (pInPoly.isInSimpleCurve(cie[0], cie[1], cPurpRed.getX(), 
                    cPurpRed.getY(), cPurpRed.getX().length)) {
                    nPurpRed++;
                }
                if (pInPoly.isInSimpleCurve(cie[0], cie[1], cOranRed.getX(), 
                    cOranRed.getY(), cOranRed.getX().length)) {
                    nOranRed++;
                }
                if (pInPoly.isInSimpleCurve(cie[0], cie[1], cYellow.getX(), 
                    cYellow.getY(), cYellow.getX().length)) {
                    nYellow++;
                }
            }
            
            float cieXRange = maxCIEX - minCIEX;
            float cieYRange = maxCIEY - minCIEY;
    
            log.info("nGTX=" + nGTX + " nLTY=" + nLTY + " n=" 
                + bestFittingPoints.size() + " "
                + " CIE: minx=" + minCIEX + " maxx=" + maxCIEX
                + " miny=" + minCIEY + " maxy=" + maxCIEY
                + " range=(" + cieXRange + "," + cieYRange + ")"
                + " nPurpRed=" + nPurpRed + " nOranRed=" + nOranRed
                + " nYellow=" + nYellow
            );

            rainbowPoints.clear();
            
            if ((nPurpRed < 10) || 
                (((float)nPurpRed/(float)bestFittingPoints.size()) < 0.05)) {
                
                return rainbowPoints;
                
            } else if ((cieXRange < 0.08) && (cieYRange < 0.08)) {
                
                return rainbowPoints;
            }
            
            if ((nGTX > 10) || (nLTY > 10)) {
                
                float frac = (float)(nGTX + nLTY)/(float)bestFittingPoints.size();
                if (frac > 0.002) {
                    rainbowPoints.addAll(bestFittingPoints);
                }
            }
        }
        
        polyFitter.plotFit(coef, rainbowPoints, colorImg.getWidth(),
            colorImg.getHeight(), 234, "rainbow points");
        
        return rainbowPoints;
    }
   
    private Map<PairInt, PairFloat> calculateContrastAndBOrR(
        Set<PairInt> points, boolean useBlue, Image originalColorImage, 
        double[] avgYRGB, int totalXOffset, int totalYOffset) {
        
        double[][] m = new double[3][];
        m[0] = new double[]{0.256, 0.504, 0.098};
        m[1] = new double[]{-0.148, -0.291, 0.439};
        m[2] = new double[]{0.439, -0.368, -0.072};
        
        double yColor = avgYRGB[0];
        
        Map<PairInt, PairFloat> map = new HashMap<PairInt, PairFloat>();
        
        for (PairInt p : points) {
            
            int x = p.getX();
            int y = p.getY();
            
            x += totalXOffset;
            y += totalYOffset;
            
            int r = originalColorImage.getR(x, y);
            int g = originalColorImage.getG(x, y);
            int b = originalColorImage.getB(x, y);

            double[] rgb = new double[]{r, g, b};
            double[] yuv = MatrixUtil.multiply(m, rgb);
            yuv = MatrixUtil.add(yuv, new double[]{16, 128, 128});

            float contrast= (float) ((yColor - yuv[0]) / yuv[0]);
            
            PairFloat crb = new PairFloat();
            crb.setX(contrast);
            if (useBlue) {
                crb.setY(b);
            } else {
                crb.setY(r);
            }
            
            map.put(p, crb);
        }
        
        return map;
    }

    private void subtractWithCorrectForNegative(GreyscaleImage gXY2, int subtract) {

        int nz = 0;
        
        if (subtract > 0) {
            
            for (int i = 0; i < gXY2.getNPixels(); i++) {
                int v = gXY2.getValue(i);
                v -= subtract;
                if (v < 0) {
                    v = 0;
                }
                gXY2.setValue(i, v);
                if (v == 0) {
                    nz++;
                }
            }
        }

        log.info("number of set 0's=" + nz);

    }

    /**
     * make downsizing corrections if any to fill in the rightmost sky points 
     * that would be missing due to down sizing
     * then the same for the lowest sky points in the image (highest y).
     * 
     * @param skyPoints
     * @param binFactor the size of the former down sizing to find sky points.
     * @param skyRowColRange the minimum row number of sky points with respect
     * to the canny edge filter product images (theta, gradientXY).
     * @param skyRowMinMax the maximum row number of sky points with respect
     * to the canny edge filter product images (theta, gradientXY).
     * @param originalColorImage the original color image
     * @param xRelativeOffset the offset in x of the canny edge filter intermediate
     * product images from the reference frame of the originalColorImage.
     * @param yRelativeOffset the offset in y of the canny edge filter intermediate
     * product images from the reference frame of the originalColorImage.
     */
    void rightAndLowerDownSizingSkyPointCorrections(Set<PairInt> skyPoints, 
        int binFactor, Map<Integer, List<PairInt>> skyRowColRange, int[] skyRowMinMax,
        Image originalColorImage, int imageWidth, int imageHeight,
        int xRelativeOffset, int yRelativeOffset) {
        
        if (binFactor == 1) {
            return;
        }
        
        // right boundary
        
        int lastCol = (originalColorImage.getWidth() - 1) - xRelativeOffset;
        
        for (int r = skyRowMinMax[0]; r <= skyRowMinMax[1]; r++) {
            
            final int row = r;
            final Integer rowIndex = Integer.valueOf(row);
            
            List<PairInt> cRanges = skyRowColRange.get(rowIndex);
                        
            for (PairInt cRange : cRanges) {
            
                int rightCol = cRange.getY();

                if (rightCol < (lastCol - binFactor + 1)) {
                    continue;
                }

                for (int i = 1; i <= binFactor; i++) {

                    int x = rightCol + i;

                    if ((x + xRelativeOffset) > (originalColorImage.getWidth() - 1)) {
                        break;
                    }

                    PairInt rightPoint = new PairInt(x, row);

                    skyPoints.add(rightPoint);

                    cRange.setY(x);
                }
            }            
        }

        // lower boundary
        int diff = imageHeight - skyRowMinMax[1];
        if ((diff <= binFactor) && (diff > 0)) {
            
            Integer rowIndex = Integer.valueOf(skyRowMinMax[1]);
            
            List<PairInt> cRanges = skyRowColRange.get(rowIndex);
                        
            for (PairInt cRange : cRanges) {
            
                for (int row = (skyRowMinMax[1] + 1); row < imageHeight; row++) {
                    
                    PairInt colRange = new PairInt(cRange.getX(), cRange.getY());
                    
                    Integer key = Integer.valueOf(row);
                    
                    // should be null:
                    List<PairInt> cLRanges = skyRowColRange.get(key);
                    if (cLRanges == null) {
                        cLRanges = new ArrayList<PairInt>();
                        skyRowColRange.put(key, cLRanges);
                    }
                    cLRanges.add(colRange);
                    
                    for (int col = colRange.getX(); col <= colRange.getY(); col++) {
                        skyPoints.add(new PairInt(col, row));
                    }
                }
            }
            skyRowMinMax[1] += diff;
        }
    }
static int outImgNum=0;
    /**
     * given seed skyPoints to start from, use conservative limits on contrast
     * and color difference to add neighbors to skyPoints.  The conservative
     * limits are meant to help avoid overrunning the skyline for low
     * contrast such as a hazy sky and snow covered peaks, for example.
     * The conservative limits do not necessarily find all sky points.
     * 
     * @param skyPoints
     * @param originalColorImage
     * @param mask
     */
    void findClouds(Set<PairInt> skyPoints, Set<PairInt> excludePoints,
        ImageExt originalColorImage, GreyscaleImage mask) {
        
        int maskWidth = mask.getWidth();
        int maskHeight = mask.getHeight();
        
        ArrayDeque<PairInt> cloudQueue = new ArrayDeque<PairInt>(skyPoints.size());
        for (PairInt skyPoint : skyPoints) {
            cloudQueue.add(skyPoint);
        }
        
        Set<PairInt> visited = new HashSet<PairInt>();
        visited.add(cloudQueue.peek());
              
        int xOffset = mask.getXRelativeOffset();
        int yOffset = mask.getYRelativeOffset();
        
        GroupPixelColors allSkyColor = new GroupPixelColors(skyPoints,
            originalColorImage, xOffset, yOffset);
        
        double rDivB = allSkyColor.getAvgRed() / allSkyColor.getAvgBlue();
        boolean skyIsRed = (rDivB > 1);
       
        log.info("==> r/b=" + rDivB
            + " redStdev=" + allSkyColor.getStdDevRed()
            + " blueStDev=" + allSkyColor.getStdDevBlue());
        
        CIEChromaticity cieC = new CIEChromaticity();
        PointInPolygon pInPoly = new PointInPolygon();
        
        int[] dxs = new int[]{-1,  0, 1, 0};
        int[] dys = new int[]{ 0, -1, 0, 1};
        
        //int[] dxs = new int[]{-1, -1,-1,  0, 0, 1,  1, 1};
        //int[] dys = new int[]{ -1, 0, 1, -1, 1, -1, 0, 1};
        
        while (!cloudQueue.isEmpty()) {

            PairInt uPoint = cloudQueue.poll();
                        
            int uX = uPoint.getX();
            int uY = uPoint.getY();
        
            //(1 + frac)*O(N) where frac is the fraction added back to stack
            for (int k = 0; k < dxs.length; k++) {
                            
                int vX = uX + dxs[k];
                int vY = uY + dys[k];
                
                if ((vX < 0) || (vX > (maskWidth - 1)) || (vY < 0) || 
                    (vY > (maskHeight - 1))) {
                    continue;
                }
            
                PairInt vPoint = new PairInt(vX, vY);
                
                if (visited.contains(vPoint) || skyPoints.contains(vPoint) ||
                    excludePoints.contains(vPoint)) {
                    continue;
                }

                visited.add(vPoint);
                
                Set<PairInt> neighbors = getThe8NeighborPixelsWithin(
                    uPoint, skyPoints, maskWidth, maskHeight);             
                
                neighbors.add(uPoint);
                
                GroupPixelColors localSky = new GroupPixelColors(neighbors, 
                    originalColorImage, xOffset, yOffset);
                
                int vIdx = originalColorImage.getInternalIndex(vX + xOffset, 
                    vY + yOffset);

                int rV = originalColorImage.getR(vIdx);
                int gV = originalColorImage.getG(vIdx);
                int bV = originalColorImage.getB(vIdx);

                float totalRGBV = rV + gV + bV;
                
                float localSkyLuma = localSky.getAverageLuma();
                
                float lumaV = originalColorImage.getLuma(vIdx);
        
                double contrastV = (localSkyLuma - lumaV)/lumaV;

                double colorDiffV = localSky.calcColorDiffToOther(rV, gV, bV);

                double skyStDevContrast = localSky.getStdDevContrast();

                double skyStDevColorDiff = localSky.getStdDevColorDiff();

                boolean doNotAddToStack = false;
                 
 //TODO: this needs adjustments...
                float rPercentV = (float)rV/totalRGBV;
                float gPercentV = (float)gV/totalRGBV;
                float bPercentV = (float)bV/totalRGBV;
                float gDivR = (float)gV/(float)rV;
                
                float cieX = originalColorImage.getCIEX(vIdx);
                float cieY = originalColorImage.getCIEY(vIdx);
                double diffCIEX = Math.abs(cieX - localSky.getAverageCIEX());
                double diffCIEY = Math.abs(cieY - localSky.getAverageCIEY());
                
                float saturation = originalColorImage.getSaturation(vIdx);
                
                boolean isBrown = (Math.abs(rPercentV - 0.5) < 0.4)
                    && (Math.abs(gPercentV - 0.32) < 0.1)
                    && (Math.abs(bPercentV - 0.17) < 0.1);

if (check(vX, xOffset, vY, yOffset)) {
log.info("contrastV=" + contrastV + " div=" + (contrastV/skyStDevContrast)
+ " colordiff=" + colorDiffV + " div=" + (colorDiffV/skyStDevColorDiff)
+ " diffciex=" + diffCIEX + " div=" + (diffCIEX/localSky.getStdDevCIEX())
+ " diffciey=" + diffCIEY + " div=" + (diffCIEY/localSky.getStdDevCIEY())
+ " r=" + rV + " g=" + gV + " b=" + bV + " bPercentV=" + bPercentV
);
     
}

                if (isBrown) {
                    
                    // trying to skip over foreground such as land or sunset + water
                    
                    if ((colorDiffV > 15*skyStDevColorDiff) && (saturation < 0.5)) {
                        
                        if ((saturation <= 0.4) && (colorDiffV > 50*skyStDevColorDiff) 
                            ) {
log.fine("FILTER 00");
                            continue;
                        }
                         if ((saturation <= 0.4) && 
                            (Math.abs(contrastV) > 10.*Math.abs(skyStDevContrast))
                            ) {
log.fine("FILTER 01");                                                        
                            continue;
                        }
                    }
                }
                
                if ( // no contrast or color change, should be sky
                    (Math.abs(contrastV) < 0.015)
                    && (colorDiffV < 10)
                    && (diffCIEX < 0.009) && (diffCIEY < 0.009)) {

log.fine("FILTER 02");
                } else if (
                    (skyStDevContrast != 0.)
                    && (Math.abs(contrastV) > 30.*skyStDevContrast)
                    ) {
 
log.fine("FILTER 03");
                    continue;
                } else if (
                    (skyStDevContrast != 0.)
                    && ((Math.abs(contrastV) > 0.5) && (Math.abs(contrastV) > 1.5*skyStDevContrast))
                    && (Math.abs(colorDiffV) > 2.5*skyStDevColorDiff)
                    ) {

log.fine("FILTER 04");
                    continue; 
                } else if (
                    (skyStDevContrast != 0.)
                    && ((Math.abs(contrastV) > 0.4) && (Math.abs(contrastV) > 2.0*skyStDevContrast))
                    && (Math.abs(colorDiffV) > 5.*skyStDevColorDiff)
                    && (!((diffCIEX < 0.02) && (diffCIEY < 0.02)))
                    ) {

log.fine("FILTER 05");
                    continue;
                    
                } else if (
                    (skyStDevContrast != 0.)
                    && ((Math.abs(contrastV) > 0.25) && (Math.abs(contrastV) > 2.0*skyStDevContrast))
                    && (Math.abs(colorDiffV) > 3.5*skyStDevColorDiff)
                    ) {

log.fine("FILTER 06");
                    continue;
                    
                } else if (
                    (skyStDevContrast != 0.)
                    && (Math.abs(contrastV) > 0.19) 
                    && (Math.abs(contrastV) > 3.5*skyStDevContrast)
                    && (Math.abs(colorDiffV) > 3.5*skyStDevColorDiff)
                    && (diffCIEX > 2.*localSky.getStdDevCIEX())
                    ) {
log.fine("FILTER 07");
                    continue;
                } else if (
                    (skyStDevContrast != 0.)
                    && (Math.abs(contrastV) > 0.19) 
                    && (Math.abs(contrastV) > 3.5*skyStDevContrast)
                    && (Math.abs(colorDiffV) > 3.5*skyStDevColorDiff)
                    && (diffCIEY > 2.*localSky.getStdDevCIEY())
                    ) {

log.fine("FILTER 08");
                    continue;

                } else if (
                    (skyStDevContrast != 0.)
                    && ((Math.abs(contrastV) > 0.1) 
                    && (Math.abs(contrastV) > 2.5*skyStDevContrast))
                    && (!((diffCIEX < 0.02) && (diffCIEY < 0.02)))
                    ) {

                    //sometimes, sky pixels are found here

log.fine("FILTER 09");
                    continue;
                       
                } else if (
                    (skyStDevContrast != 0.)
                    && (Math.abs(contrastV) > 5.*skyStDevContrast)
                    && (Math.abs(colorDiffV) > 5.*skyStDevColorDiff)
                    && 
                    // if cieXY diffs are zero and within stdev, these are sky,
                    // so test for opposite for boundary pixel
                    (!((diffCIEX < 0.001) && (diffCIEX < 1.5*localSky.getStdDevCIEX())
                    && (diffCIEY < 0.001) && (diffCIEY < 1.5*localSky.getStdDevCIEY())))
                    
                    && (skyStDevContrast > 0.005)
                    && (skyStDevColorDiff > 1.)
                    ) {

log.fine("FILTER 10");
                    continue;
                    
                } else if (skyIsRed) {
                    
                    if (
                        // contrast is defined by luma, so might be weak near
                        // skyline near sun for example
                        (skyStDevContrast != 0.)
                        && (Math.abs(colorDiffV) > 3.*skyStDevColorDiff)
                        && (diffCIEX > 0.03) && (diffCIEX > 3.*localSky.getStdDevCIEX())
                        ) {

log.fine("FILTER 11");
                        continue;
                    } else if (
                        // contrast is defined by luma, so might be weak near
                        // skyline near sun for example
                        (skyStDevContrast != 0.)
                        && (Math.abs(colorDiffV) > 3.*skyStDevColorDiff)
                        && (diffCIEY > 0.03) && (diffCIEY > 3.*localSky.getStdDevCIEY())
                        ) {

log.fine("FILTER 12");
                        continue;
                    } else if (
                        // contrast is defined by luma, so might be weak near
                        // skyline near sun for example
                        (skyStDevContrast != 0.)
                        && (Math.abs(colorDiffV) > 1.1*skyStDevColorDiff)
                        && ((diffCIEX > 0.065) 
                        && (diffCIEX > 1.1*localSky.getStdDevCIEX()))
                        ) {

log.fine("FILTER 13");
                        continue;
                    } else if (
                        // contrast is defined by luma, so might be weak near
                        // skyline near sun for example
                        (skyStDevContrast != 0.)
                        && (Math.abs(colorDiffV) > 1.1*skyStDevColorDiff)
                        && ((diffCIEY > 0.03) 
                        && (diffCIEY > 3.*localSky.getStdDevCIEY()))
                        ) {

log.fine("FILTER 14");
                        continue;
                        
                    } else if (
                        // contrast is defined by luma, so might be weak near
                        // skyline near sun for example
                        (skyStDevContrast != 0.)
                        && (Math.abs(colorDiffV) > 10.*skyStDevColorDiff)
                        && ((diffCIEX > 0.02) 
                        && (diffCIEX > 1.5*localSky.getStdDevCIEX()))
                        ) {

log.fine("FILTER 15");
                        continue;
                    } else if (
                        // contrast is defined by luma, so might be weak near
                        // skyline near sun for example
                        (skyStDevContrast != 0.)
                        && (Math.abs(colorDiffV) > 10.*skyStDevColorDiff)
                        && ((diffCIEY > 0.02) 
                        && (diffCIEY > 1.5*localSky.getStdDevCIEY()))
                        ) {

log.fine("FILTER 16");
                        continue;
                 
                    } else if (
                        (skyStDevContrast != 0.)
                        && (contrastV > 0.01) 
                        && (colorDiffV > 15*skyStDevColorDiff)
                        && ((diffCIEX > 0.03) || (diffCIEY > 0.03))
                        ) {
                        
                        ArrayPair yellowGreenOrange = cieC.getYellowishGreenThroughOrangePolynomial();
                        if (pInPoly.isInSimpleCurve(cieX, cieY,
                            yellowGreenOrange.getX(), yellowGreenOrange.getY(),yellowGreenOrange.getX().length)) {
log.fine("FILTER 17");
                        } else {

log.fine("FILTER 18");
                            continue;
                        }
                    
                    } else if (skyStDevContrast == 0.) {
                        if (contrastV >= 0.) {
log.fine("FILTER 19");
                            doNotAddToStack = true;
                        }
                        
                    } else {
                        //TODO:  if there are sun points, need a zone of
                        // avoidance to not erode the foreground 

log.fine("FILTER 20");
                    }

                } else {
                    //blue filters
                    if (
                        (skyStDevContrast > 0.0)
                        && (contrastV < 0.05) 
                        && ((Math.abs(contrastV)/skyStDevContrast) > 1.5)
                        && ((Math.abs(colorDiffV)/skyStDevColorDiff) > 3.0)
                        && (diffCIEX < 0.01) && (diffCIEY < 0.01)
                        && (!(bPercentV > 0.37) && (bV > 199) && (gV > 199))
                        ) {

log.fine("FILTER 21");
                        continue;
                    } else if (
                        (skyStDevContrast > 0.0)
                        && (contrastV < 0.05) 
                        && ((Math.abs(colorDiffV)/skyStDevColorDiff) > 1.1)
                        && (diffCIEX > 2.5*localSky.getStdDevCIEX())
                        && (diffCIEY > 2.5*localSky.getStdDevCIEY())
                        
                        && (skyStDevContrast > 0.005)
                        && (skyStDevColorDiff > 1.)
                    ){

log.fine("FILTER 22");
                        continue;
                    } else if (
                        (skyStDevContrast > 0.0)
                        && (Math.abs(contrastV) < 0.05)
                        && ((Math.abs(colorDiffV)/skyStDevColorDiff) > 3.0)
                        && (diffCIEX > 1.5*localSky.getStdDevCIEX())
                        && (diffCIEY > 1.5*localSky.getStdDevCIEY())                        
                    ) {

log.fine("FILTER 23");
                        continue;
                    } else if (
                        (contrastV < 0.05) 
                        && (colorDiffV < 16)
                        && ((Math.abs(colorDiffV)/skyStDevColorDiff) < 2.5)
                        ) {

log.fine("FILTER 24");
                    } else if (
                        (skyStDevContrast > 0.0)
                        && (contrastV > 0.0) 
                        && ((Math.abs(contrastV)/skyStDevContrast) < 2.5)
                        && (diffCIEX < 0.01) && (diffCIEY < 0.01)
                        ) {

log.fine("FILTER 25");
                   } else if (
                        (skyStDevContrast == 0.0) && (skyStDevColorDiff == 0.0)
                        && (contrastV < 0.01) && (colorDiffV < 16)
                        ) {

log.fine("FILTER 26");
                    } else if (
                    // bright grey clouds
                    (Math.abs(contrastV) < 0.04)
                    && (diffCIEX < 0.002) && (diffCIEY < 0.003) 
                    &&
                        ((Math.abs(0.30 - rPercentV) < 0.07) 
                        && (Math.abs(0.33 - gPercentV) < 0.03) 
                        && (Math.abs(0.36 - bPercentV) < 0.06)
                        && (gV > 170) && (bV > 170))
                    ) {

log.fine("FILTER 27");
                    } else if (
                    (contrastV < 0.05)
                    && (diffCIEX < 0.005) && (diffCIEY < 0.005)
                    &&
                        !((Math.abs(0.33 - rPercentV) < 0.08) 
                        && (Math.abs(0.33 - gPercentV) < 0.03) 
                        && (Math.abs(0.33 - bPercentV) < 0.06)
                        && (gV > 199) && (bV > 199))
                    ) {

log.fine("FILTER 28");
                        continue; 
                    } else {

log.fine("FILTER 29");
                        continue;

                    }
                }
                
                skyPoints.add(vPoint);

                if (!doNotAddToStack) {
                    cloudQueue.add(vPoint);
                }
            }
        }
                
        for (PairInt p : skyPoints) {
            int x = p.getX();
            int y = p.getY();
            mask.setValue(x, y, 0);
        }
    }

private boolean check(int vX, int xOffset, int vY, int yOffset) {
    /*if (((vX + xOffset)>=97) && ((vX + xOffset)<=97) && ((vY + yOffset)==172)) {
        return true;
    }*/
    return false;
}

    public enum PARAMS {
        ABSOLUTE_CONTRAST,
        ABSOLUTE_DIFF_BLUE_OR_RED,
        STDEV_CONTRAST,
        STDEV_BLUE_OR_RED,
        DIFF_CIEX,
        DIFF_CIEY,
        STDEV_CIEX,
        STDEV_CIEY,
        INT_ONE,
        
        /*ABSOLUTE_BLUE_OR_RED,
        ABSOLUTE_DIFF_HUE,
        ABSOLUTE_DIFF_CIEX,
        ABSOLUTE_DIFF_CIEY,
        CONTRAST,
        BLUE_OR_RED,
        HUE,
        CIEX,
        CIEY,
        DIFF_BLUE_OR_RED,
        DIFF_HUE,
        STDEV_HUE,*/
    }
    
    public enum COMPARISON {
        LESS_THAN,
        GREATER_THAN
    }
    
    public enum SKYCONDITIONAL {
        ALL, BLUE, RED
    }
    
    /**
     * given seed skyPoints to start from, use conservative limits on contrast
     * and color difference to add neighbors to skyPoints.  The conservative
     * limits are meant to help avoid overrunning the skyline for low
     * contrast such as a hazy sky and snow covered peaks, for example.
     * The conservative limits do not necessarily find all sky points.
     * 
     * <pre>
     * The contrast and color filters are supplied by params1, params2, gtOrLT,
     * and coefficients of format:
     *   (params1[i]/params2[i]) gtOrLT[i] coefficients[i]
     * </pre>
     * 
     * @param skyPoints
     * @param originalColorImage
     * @param mask
     */
    void findClouds2(Set<PairInt> skyPoints, Set<PairInt> excludePoints,
        ImageExt originalColorImage, GreyscaleImage mask, 
        PARAMS[] params1, PARAMS[] params2, COMPARISON[] gtOrLT, 
        float[] coefficients, SKYCONDITIONAL[] skyConditional) {
        
        int maskWidth = mask.getWidth();
        int maskHeight = mask.getHeight();
        
        ArrayDeque<PairInt> cloudQueue = new ArrayDeque<PairInt>(skyPoints.size());
        for (PairInt skyPoint : skyPoints) {
            cloudQueue.add(skyPoint);
        }
        
        Set<PairInt> visited = new HashSet<PairInt>();
        visited.add(cloudQueue.peek());
              
        int xOffset = mask.getXRelativeOffset();
        int yOffset = mask.getYRelativeOffset();
        
        GroupPixelColors allSkyColor = new GroupPixelColors(skyPoints,
            originalColorImage, xOffset, yOffset);
        
        double rDivB = allSkyColor.getAvgRed() / allSkyColor.getAvgBlue();
        boolean skyIsRed = (rDivB > 1);
       
        log.info("==> r/b=" + rDivB
            + " redStdev=" + allSkyColor.getStdDevRed()
            + " blueStDev=" + allSkyColor.getStdDevBlue());
        
        CIEChromaticity cieC = new CIEChromaticity();
        PointInPolygon pInPoly = new PointInPolygon();
        
        int[] dxs = new int[]{-1,  0, 1, 0};
        int[] dys = new int[]{ 0, -1, 0, 1};
                
        while (!cloudQueue.isEmpty()) {

            PairInt uPoint = cloudQueue.poll();
                        
            int uX = uPoint.getX();
            int uY = uPoint.getY();
        
            //(1 + frac)*O(N) where frac is the fraction added back to stack
            for (int k = 0; k < dxs.length; k++) {
                            
                int vX = uX + dxs[k];
                int vY = uY + dys[k];
                
                if ((vX < 0) || (vX > (maskWidth - 1)) || (vY < 0) || 
                    (vY > (maskHeight - 1))) {
                    continue;
                }
            
                PairInt vPoint = new PairInt(vX, vY);
                
                if (visited.contains(vPoint) || skyPoints.contains(vPoint) ||
                    excludePoints.contains(vPoint)) {
                    continue;
                }

                visited.add(vPoint);
                
                Set<PairInt> neighbors = getThe8NeighborPixelsWithin(
                    uPoint, skyPoints, maskWidth, maskHeight);             
                
                neighbors.add(uPoint);
                
                GroupPixelColors localSky = new GroupPixelColors(neighbors, 
                    originalColorImage, xOffset, yOffset);
                
                int vIdx = originalColorImage.getInternalIndex(vX + xOffset, 
                    vY + yOffset);

                int rV = originalColorImage.getR(vIdx);
                int gV = originalColorImage.getG(vIdx);
                int bV = originalColorImage.getB(vIdx);

                float totalRGBV = rV + gV + bV;
                
                float localSkyLuma = localSky.getAverageLuma();
                
                float lumaV = originalColorImage.getLuma(vIdx);
        
                double contrastV = (localSkyLuma - lumaV)/lumaV;

                double colorDiffV = localSky.calcColorDiffToOther(rV, gV, bV);

                double skyStDevContrast = localSky.getStdDevContrast();

                double skyStDevColorDiff = localSky.getStdDevColorDiff();
                 
 //TODO: this needs adjustments...
                float rPercentV = (float)rV/totalRGBV;
                float gPercentV = (float)gV/totalRGBV;
                float bPercentV = (float)bV/totalRGBV;
                
                float cieX = originalColorImage.getCIEX(vIdx);
                float cieY = originalColorImage.getCIEY(vIdx);
                double diffCIEX = Math.abs(cieX - localSky.getAverageCIEX());
                double diffCIEY = Math.abs(cieY - localSky.getAverageCIEY());
                
                float saturation = originalColorImage.getSaturation(vIdx);
                
                boolean isBrown = (Math.abs(rPercentV - 0.5) < 0.4)
                    && (Math.abs(gPercentV - 0.32) < 0.1)
                    && (Math.abs(bPercentV - 0.17) < 0.1);

if (check(vX, xOffset, vY, yOffset)) {
log.info("contrastV=" + contrastV + " div=" + (contrastV/skyStDevContrast)
+ " colordiff=" + colorDiffV + " div=" + (colorDiffV/skyStDevColorDiff)
+ " diffciex=" + diffCIEX + " div=" + (diffCIEX/localSky.getStdDevCIEX())
+ " diffciey=" + diffCIEY + " div=" + (diffCIEY/localSky.getStdDevCIEY())
+ " r=" + rV + " g=" + gV + " b=" + bV + " bPercentV=" + bPercentV
);
     
}

                if (isBrown) {
                    
                    // trying to skip over foreground such as land or sunset + water
                    
                    if ((colorDiffV > 15*skyStDevColorDiff) && (saturation < 0.5)) {
                        
                        if ((saturation <= 0.4) && (colorDiffV > 50*skyStDevColorDiff) 
                            ) {
log.fine("FILTER 00");
                            continue;
                        }
                         if ((saturation <= 0.4) && 
                            (Math.abs(contrastV) > 10.*Math.abs(skyStDevContrast))
                            ) {
log.fine("FILTER 01");                                                        
                            continue;
                        }
                    }
                }
                
                boolean isSolarYellow = false;
                if (skyIsRed) {
                    if ((skyStDevContrast != 0.) && (contrastV > 0.01) 
                        && (colorDiffV > 15*skyStDevColorDiff)
                        && ((diffCIEX > 0.03) || (diffCIEY > 0.03))
                        ) {
                        ArrayPair yellowGreenOrange = cieC.getYellowishGreenThroughOrangePolynomial();
                        if (pInPoly.isInSimpleCurve(cieX, cieY,
                            yellowGreenOrange.getX(), yellowGreenOrange.getY(),yellowGreenOrange.getX().length)) {
                            isSolarYellow = true;
                        }
                    }
                }
                
                if ( // no contrast or color change, should be sky
                    (Math.abs(contrastV) < 0.015)
                    && (colorDiffV < 10)
                    && (diffCIEX < 0.009) && (diffCIEY < 0.009)) {
                    // this is a sky point
                    
log.fine("FILTER 02");
                } else if (isSolarYellow) {
                    // this is a sky point
log.fine("FILTER 03");
                } else {
                    
                    boolean isNotSky = false;
                    
                    for (int i = 0; i < params1.length; i++) {
                        
                        float param1;
                        
                        switch(params1[i]) {
                            case ABSOLUTE_CONTRAST:
                                param1 = (float)Math.abs(contrastV);
                                break;
                            case ABSOLUTE_DIFF_BLUE_OR_RED:
                                param1 = (float)Math.abs(colorDiffV);
                                break;
                            case STDEV_CONTRAST:
                                param1 = (float)skyStDevContrast;
                                break;
                            case STDEV_BLUE_OR_RED:
                                param1 = (float)skyStDevColorDiff;
                                break;
                            case DIFF_CIEX:
                                param1 = (float)diffCIEX;
                                break;
                            case DIFF_CIEY:
                                param1 = (float)diffCIEY;
                                break;
                            case STDEV_CIEX:
                                param1 = localSky.stdDevCIEX;
                                break;
                            case STDEV_CIEY:
                                param1 = localSky.stdDevCIEY;
                                break;
                            case INT_ONE:
                                param1 = 1.0f;
                                break;
                            default:
                                param1 = 1.0f;
                        }
                        
                        float param2;
                        
                        switch(params2[i]) {
                            case ABSOLUTE_CONTRAST:
                                param2 = (float)Math.abs(contrastV);
                                break;
                            case ABSOLUTE_DIFF_BLUE_OR_RED:
                                param2 = (float)Math.abs(colorDiffV);
                                break;
                            case STDEV_CONTRAST:
                                param2 = (float)skyStDevContrast;
                                break;
                            case STDEV_BLUE_OR_RED:
                                param2 = (float)skyStDevColorDiff;
                                break;
                            case DIFF_CIEX:
                                param2 = (float)diffCIEX;
                                break;
                            case DIFF_CIEY:
                                param2 = (float)diffCIEY;
                                break;
                            case STDEV_CIEX:
                                param2 = localSky.stdDevCIEX;
                                break;
                            case STDEV_CIEY:
                                param2 = localSky.stdDevCIEY;
                                break;
                            case INT_ONE:
                                param2 = 1.0f;
                                break;
                            default:
                                param2 = 1.0f;
                        }
                               
                        float coeff = coefficients[i];
                        
                        if (skyConditional[i].ordinal() == SKYCONDITIONAL.ALL.ordinal()) {
                            if (gtOrLT[i].ordinal() == COMPARISON.GREATER_THAN.ordinal()) {
                                if ((param1 / param2) > coeff) {
                                    isNotSky = true;
                                    break;
                                }
                            } else {
                                if ((param1 / param2) < coeff) {
                                    isNotSky = true;
                                    break;
                                }
                            }
                        } else if ((skyConditional[i].ordinal() == 
                            SKYCONDITIONAL.RED.ordinal()) && skyIsRed) {
                            if (gtOrLT[i].ordinal() == COMPARISON.GREATER_THAN.ordinal()) {
                                if ((param1 / param2) > coeff) {
                                    isNotSky = true;
                                    break;
                                }
                            } else {
                                if ((param1 / param2) < coeff) {
                                    isNotSky = true;
                                    break;
                                }
                            }
                        } else if ((skyConditional[i].ordinal() == 
                            SKYCONDITIONAL.BLUE.ordinal()) && !skyIsRed) {
                            if (gtOrLT[i].ordinal() == COMPARISON.GREATER_THAN.ordinal()) {
                                if ((param1 / param2) > coeff) {
                                    isNotSky = true;
                                    break;
                                }
                            } else {
                                if ((param1 / param2) < coeff) {
                                    isNotSky = true;
                                    break;
                                }
                            }
                        }                        
                    }
                    
                    if (isNotSky) {
                        continue;
                    }
                }
                
                skyPoints.add(vPoint);

                cloudQueue.add(vPoint);
            }
        }
 
        for (PairInt p : skyPoints) {
            int x = p.getX();
            int y = p.getY();
            mask.setValue(x, y, 0);
        }
    }
   
    private GroupPixelColors[] partitionInto3ByColorDifference(Set<PairInt> skyPoints, 
        ImageExt originalColorImage, int xOffset, int yOffset, int[] outputCounts) {
         
        GroupPixelColors allSkyColor = new GroupPixelColors(skyPoints,
            originalColorImage, xOffset, yOffset);
        
        int i = 0;
        float[] colorDiffs = new float[skyPoints.size()];
        PairInt[] ps = new PairInt[colorDiffs.length];
        for (PairInt p : skyPoints) {
            int x = p.getX() + xOffset;
            int y = p.getY() + yOffset;
            int r = originalColorImage.getR(x, y);
            int g = originalColorImage.getG(x, y);
            int b = originalColorImage.getB(x, y);
            colorDiffs[i] = allSkyColor.calcColorDiffToOther(r, g, b);
            ps[i] = p;
            i++;
        }
        
        float minColorDiff = MiscMath.findMin(colorDiffs);
        float maxColorDiff = MiscMath.findMax(colorDiffs);
        float[] yErr = Errors.populateYErrorsBySqrt(colorDiffs);
        HistogramHolder h = Histogram.createSimpleHistogram(minColorDiff,
            maxColorDiff, 3, colorDiffs, yErr);
        
        for (i = 0; i < 3; i++) {
            outputCounts[i] = h.getYHist()[i];
        }

        Set<PairInt> set0 = new HashSet<PairInt>();
        Set<PairInt> set1 = new HashSet<PairInt>();
        Set<PairInt> set2 = new HashSet<PairInt>();
        
        float binSize = h.getXHist()[1] - h.getXHist()[0];
        
        for (i = 0; i < colorDiffs.length; i++) {
            float cd = colorDiffs[i];
            int binN = (int)((cd - minColorDiff)/binSize);
            switch(binN) {
                case 0:
                    set0.add(ps[i]);
                    break;
                case 1:
                    set1.add(ps[i]);
                    break;
                default:
                    set2.add(ps[i]);
                    break;
            }
        }
        
        GroupPixelColors[] sets = new GroupPixelColors[3];
        sets[0] = new GroupPixelColors(set0, originalColorImage, xOffset, yOffset);
        sets[1] = new GroupPixelColors(set1, originalColorImage, xOffset, yOffset);
        sets[2] = new GroupPixelColors(set2, originalColorImage, xOffset, yOffset);
        
        return sets;
    }
    
    protected GroupPixelColors[] partitionInto3ByBrightness(
        Set<PairInt> skyPoints, 
        ImageExt originalColorImage, int xOffset, int yOffset, HistogramHolder[] h) {
         
        int i = 0;
        float[] brightness = new float[skyPoints.size()];
        PairInt[] ps = new PairInt[brightness.length];
        for (PairInt p : skyPoints) {
            int x = p.getX() + xOffset;
            int y = p.getY() + yOffset;
            int idx = originalColorImage.getInternalIndex(x, y);
            
            brightness[i] = originalColorImage.getBrightness(idx);
            ps[i] = p;
            i++;
        }
        
        float min = MiscMath.findMin(brightness);
        float max = MiscMath.findMax(brightness);
        float[] yErr = Errors.populateYErrorsBySqrt(brightness);
        h[0] = Histogram.createSimpleHistogram(min, max, 3, 
            brightness, yErr);
        
try {
    h[0].plotHistogram("sky brightness", 235);
} catch(Exception e) {
    
}     

        Set<PairInt> set0 = new HashSet<PairInt>();
        Set<PairInt> set1 = new HashSet<PairInt>();
        Set<PairInt> set2 = new HashSet<PairInt>();
        
        float binSize = h[0].getXHist()[1] - h[0].getXHist()[0];
        
        for (i = 0; i < brightness.length; i++) {
            float cd = brightness[i];
            int binN = (int)((cd - min)/binSize);
            switch(binN) {
                case 0:
                    set0.add(ps[i]);
                    break;
                case 1:
                    set1.add(ps[i]);
                    break;
                default:
                    set2.add(ps[i]);
                    break;
            }
        }
        
        GroupPixelColors[] sets = new GroupPixelColors[3];
        sets[0] = new GroupPixelColors(set0, originalColorImage, xOffset, yOffset);
        sets[1] = new GroupPixelColors(set1, originalColorImage, xOffset, yOffset);
        sets[2] = new GroupPixelColors(set2, originalColorImage, xOffset, yOffset);
        
        return sets;
    }
  
        /*
        scattering:
            for atmospheric particles, their size is much smaller than
            lambda, the wavelength of optical light:  
                Rayleigh scattering is prop to lambda^-4.
                this is responsible for blue or red skies compared to 
                the yellow color of the sun 
                (photosphere peak is near 5500 Angstroms).
        
            for clouds, the water droplets are the scatters and their
            size is comparable to lambda for optical light:
                Mie scattering affects optical colors roughly equally,
                so leads to less light without a color change.
        
        One could use the position of the sun to approximate the airmass for
        different depths of features and calculate the transmission of the solar 
        spectrum through atmosphere with absorption by water, scattering and 
        absorption by aerosols and Rayleigh scattering in addition to Mie 
        scattering by water droplets in the clouds
        combined with a model for water vapor in the atmosphere and a range 
        of cloud optical depths and locations and a camera response function 
        to estimate the colors expected in images...
        
        The source function upon the clouds could be modeled with
        http://rredc.nrel.gov/solar/pubs/spectral/model/section5.html
        
        Since that isn't feasible for now, could take skyPoints colors and colors
        from the outer patch of sun and look for some variations of those in the
        image close to, but unconnected to existing skyPoints.
        */


    private int count(List<PairIntArray> zeroPointLists) {
        
        int n = 0;
        for (PairIntArray p : zeroPointLists) {
            n += p.getN();
        }
        
        return n;
    }

    void debugPlot(Set<PairInt> extSkyPoints, Image originalColorImage, 
        int xOffset, int yOffset, String outputPrefixForFileName) {
        
        //plot is made in aspects
        
    }
    
    void debugPlot(Set<PairInt> r0, Set<PairInt> r1, Set<PairInt> r2, 
        Set<PairInt> r3, Image originalColorImage, int xOffset, int yOffset, 
        String filtered_out_of_clouds) {
        
        //plot is made in aspects
    }

    private void filterToAddCloudPixel(Set<PairInt> excludeFromThesePoints, 
        Set<PairInt> skyPoints, ImageExt origColorImg, 
        int xOffset, int yOffset, PairInt p,
        HistogramHolder[] brightnessHistogram, 
        GroupPixelColors[] skyBinsByBrightness,
        boolean skyIsRed, boolean skyIsPurple, boolean hasDarkGreyClouds,
        boolean useKEqualsZero, boolean pointIsEmbeddedInSky,
        Set<PairInt> outputCloudPoints, Stack<PairInt> outputStack) {
        
        if (skyPoints.contains(p) || excludeFromThesePoints.contains(p)) {
            return;
        }
        
        int col = p.getX() + xOffset;
        int row = p.getY() + yOffset;
        
        if ((col < 0) || (col > (origColorImg.getWidth() - 1)) ||
            (row < 0) || (row > (origColorImg.getHeight() - 1))) {
            return;
        }

        int r = origColorImg.getR(col, row);
        int g = origColorImg.getG(col, row);
        int b = origColorImg.getB(col, row);
        
        int idx = origColorImg.getInternalIndex(col, row);
        
        float luma = origColorImg.getLuma(idx);
        
        for (int k = 2; k >= 0; k--) {                    

            if ((k == 0) && !useKEqualsZero) {
                continue;
            }             
            
            float skyLuma = skyBinsByBrightness[k].getAverageLuma();
            
            float contrast = (skyLuma - luma)/luma;

            double diffContrast = contrast - skyLuma;
            
            if (skyIsRed && (diffContrast > 0) && !pointIsEmbeddedInSky) {
                continue;
            }
            
            double diffR = r - skyBinsByBrightness[k].getAvgRed();
            if (diffR < 0) {
                diffR *= -1;
            }
            double diffG = g - skyBinsByBrightness[k].getAvgGreen();
            if (diffG < 0) {
                diffG *= -1;
            }
            double diffB = b - skyBinsByBrightness[k].getAvgBlue();
            if (diffB < 0) {
                diffB *= -1;
            }

            if (
                (pointIsEmbeddedInSky) ||
                (
                (
                    (diffR <= skyBinsByBrightness[k].getStdDevRed())
                    ||
                    (skyIsRed && (r > skyBinsByBrightness[k].getAvgRed()))
                    ||
                    (skyIsRed && (r > 155) && 
                        (diffR <= 3.5*skyBinsByBrightness[k].getStdDevRed())
                    )
                    ||
                    (!skyIsRed && /*hasDarkGreyClouds &&*/ 
                        (diffR <= 2.0*skyBinsByBrightness[k].getStdDevRed())
                    )
                )
                && (
                    (diffG <= 1.5*skyBinsByBrightness[k].getStdDevGreen())
                    ||
                    //consider (contrast < 0.) && (contrast > -0.05)
                    (skyIsRed && (contrast < 0.) && (g > 130) && 
                        (diffG <= 2.0*skyBinsByBrightness[k].getStdDevGreen())
                    )
                    ||
                    ((skyIsPurple || skyIsRed) && (contrast > 0.) && (g < 130) && 
                        (diffG <= 2.0*skyBinsByBrightness[k].getStdDevGreen())
                    )
                    ||
                    (!skyIsRed && hasDarkGreyClouds && 
                        (diffG <= 2.0*skyBinsByBrightness[k].getStdDevGreen())
                    )
                )
                && (
                    (diffB <= 1.2*skyBinsByBrightness[k].getStdDevBlue())
                    ||
                    (skyIsPurple && (contrast > 0.) && (b < 130) && 
                        (diffB <= 2.5*skyBinsByBrightness[k].getStdDevBlue())
                    )
                    ||
                    (!skyIsRed && /*hasDarkGreyClouds &&*/
                        (diffB <= 1.5*skyBinsByBrightness[k].getStdDevBlue())
                    )
                )
                )
                ) {

                outputCloudPoints.add(p);

                outputStack.add(p);

                break;
            }
        }
    }

    void populatePixelExtColors(Set<PairInt> points, ImageExt originalColorImage, 
        GreyscaleImage mask) {
        
        int xOffset = mask.getXRelativeOffset(); 
        int yOffset = mask.getYRelativeOffset();
        
        for (PairInt p : points) {
            
            int x = p.getX() + xOffset;
            int y = p.getY() + yOffset;
            int index = mask.getIndex(x, y);
            
            originalColorImage.calculateColorIncludingNeighbors(index, 0);
        }
    }

    private Set<PairInt> removeReflectedSun(List<PairIntArray> zeroPointLists, 
        ImageExt colorImg, GreyscaleImage thetaImg) {
        
        List<Integer> remove = new ArrayList<Integer>();
        
        for (int i = 0; i < zeroPointLists.size(); i++) {
            
            PairIntArray zeroPointList = zeroPointLists.get(i);
            
            Set<PairInt> pointSet = new HashSet<PairInt>();
            
            for (int ii = 0; ii < zeroPointList.getN(); ii++) {
                int x = zeroPointList.getX(ii);
                int y = zeroPointList.getY(ii);
                PairInt p = new PairInt(x, y);
                pointSet.add(p);
            }
            
            Set<PairInt> sunPoints = new HashSet<PairInt>();
            
            double[] ellipFitParams = findSunConnectedToSkyPoints(pointSet,
                new HashSet<PairInt>(),
                colorImg, thetaImg.getXRelativeOffset(),
                thetaImg.getYRelativeOffset(), false, sunPoints);
            
            float fracSun = (float)sunPoints.size()/(float)pointSet.size();
            
            if (fracSun >= 0.8) {
                remove.add(Integer.valueOf(i));
            }
        }
                
        if (!remove.isEmpty()) {
            
            Set<PairInt> reflectedSun = new HashSet<PairInt>();
            
            for (int i = (remove.size() - 1); i > -1; i--) {
                
                int idx = remove.get(i).intValue();
                
                PairIntArray pai = zeroPointLists.remove(idx);
                
                for (int ii = 0; ii < pai.getN(); ii++) {
                    int x = pai.getX(ii);
                    int y = pai.getY(ii);
                    
                    reflectedSun.add(new PairInt(x, y));
                }
            }
            
            log.info("removed " + remove.size() + 
                " points that resembled reflected sun");
            
            return reflectedSun;
        }
        
        return new HashSet<PairInt>();
    }

    private void growForLowContrastLimits(Set<PairInt> points, 
        Set<PairInt> excludePoints, Map<Integer, List<PairInt>> skyRowColRange,
        int[] skyRowMinMax, ImageExt colorImg, GreyscaleImage mask) {

        // it tries to avoid adding snowy mountain tops to hazy sky pixels,
        // for example.

        int cWidth = colorImg.getWidth();
        int cHeight = colorImg.getHeight();
        int mWidth = mask.getWidth();
        int mHeight = mask.getHeight();
        int xOffset = mask.getXRelativeOffset();
        int yOffset = mask.getYRelativeOffset();
        int imageMaxColumn = cWidth;
        int imageMaxRow = cHeight;

        GroupPixelColors allSkyColor = new GroupPixelColors(points, colorImg, 
            xOffset, yOffset);

        //TODO: change to use border pixels and their 24 neighbors?
        PixelColors foregroundColor = getAveragedColorOfDifference(points,
            colorImg, xOffset, yOffset);

        double rDivB = allSkyColor.getAvgRed()/allSkyColor.getAvgBlue();

        boolean skyIsRed = (rDivB > 1);

        PerimeterFinder perimeterFinder = new PerimeterFinder();

        int[] dxs = new int[]{-1, -1,  0,  1, 1, 1, 0, -1};
        int[] dys = new int[]{ 0, -1, -1, -1, 0, 1, 1,  1};

        int nAdded = 0;
        int nIter = 0;
        int nMaxIter = 1;

        Set<PairInt> borderPixels = perimeterFinder.getBorderPixels(
            skyRowColRange, skyRowMinMax, imageMaxColumn, imageMaxRow);
        
        Set<PairInt> borderPixelsAndNeighbors = new HashSet<PairInt>(borderPixels);

        while ((nIter == 0) || ((nIter < nMaxIter) && (nAdded > 0))) {
            
            nAdded = 0;            

            for (PairInt uPoint : borderPixels) {

                int uX = uPoint.getX() + xOffset;
                int uY = uPoint.getY() + yOffset;
                
                // for each u, aggregate it's neighbors colors if they are 
                // in the "points" set

                for (int k = 0; k < dxs.length; k++) {

                    int vX = uX + dxs[k];
                    int vY = uY + dys[k];

                    if ((vX < 0) || (vX > (cWidth - 1)) || (vY < 0) || 
                        (vY > (cHeight - 1))) {
                        continue;
                    }

                    PairInt vPoint = new PairInt(vX - xOffset, vY - yOffset);
                    
                    if (uPoint.equals(vPoint) || !points.contains(vPoint) ||
                        excludePoints.contains(vPoint)) {
                        continue;
                    }
                                        
                    Set<PairInt> neighbors = getThe8NeighborPixelsWithin(
                        uPoint, points, mWidth, mHeight);
                        
                    borderPixelsAndNeighbors.addAll(neighbors);     
                }
            }
            
            Set<PairInt> addedPixels = growForLowContrastLimits(points, 
                excludePoints, colorImg, mask, skyIsRed, borderPixels, 
                borderPixelsAndNeighbors, foregroundColor);
            
            borderPixels = addedPixels;

            nAdded = addedPixels.size();
            nIter++;
        }
    }
        
    private Set<PairInt> growForLowContrastLimits(Set<PairInt> points, 
        Set<PairInt> excludePoints, ImageExt colorImg, GreyscaleImage mask,
        boolean skyIsRed, Set<PairInt> borderPixels,
        Set<PairInt> borderPixelsAndNeighbors, PixelColors foregroundColor) {
        
        if (borderPixels.isEmpty()) {
            return new HashSet<PairInt>();
        }
        
        int xOffset = mask.getXRelativeOffset();
        int yOffset = mask.getYRelativeOffset();
        int width = mask.getWidth();
        int height = mask.getHeight();
        
        GroupPixelColors gpc = new GroupPixelColors(borderPixelsAndNeighbors,
            colorImg, xOffset, yOffset);
        
        //TODO: ***change to BFS w/ insert at process pair insted of DFS***
        
        Stack<PairInt> stack = new Stack<PairInt>();
        stack.addAll(borderPixels);
        
        Set<PairInt> visited = new HashSet<PairInt>();
        visited.add(stack.peek());
        
        //int[] dxs = new int[]{-1,  0, 1, 0};
        //int[] dys = new int[]{ 0, -1, 0, 1};
        int[] dxs = new int[]{-1, -1,  0,  1, 1, 1, 0, -1};
        int[] dys = new int[]{ 0, -1, -1, -1, 0, 1, 1,  1};
        
        double skyMinusForegroundColor = 
            foregroundColor.calculateColorDiffererenceToOther(
            Math.round(gpc.getAvgRed()), 
            Math.round(gpc.getAvgGreen()), Math.round(gpc.getAvgBlue()));
        
        log.info("skyMinusForegroundColor=" + skyMinusForegroundColor + 
            " skyIsRed=" + skyIsRed);
        
        //TODO:  improve this
        int firstSubtrTol = 7;
        int secondSubtrTol = 2;
        if (skyMinusForegroundColor >= 200) {
            firstSubtrTol = 9;
            secondSubtrTol = 8;
        }
        
        int nAdded = 0;
        
        //Stack<PairInt> addedSky = new Stack<PairInt>();
        Set<PairInt> addedSky = new HashSet<PairInt>();
        
        Stack<PixelColors> addedSkyColorToAdd = new Stack<PixelColors>();
        
        while (!stack.isEmpty()) {

            PairInt uPoint = stack.pop();

            int uX = uPoint.getX() + xOffset;
            int uY = uPoint.getY() + yOffset;

            GroupPixelColors uLocalColors = null;
            
            for (int k = 0; k < dxs.length; k++) {

                int vX = uX + dxs[k];
                int vY = uY + dys[k];

                if ((vX < 0) || (vX > (width - 1)) || (vY < 0) || 
                    (vY > (height - 1))) {
                    continue;
                }

                PairInt vPoint = new PairInt(vX - xOffset, vY - yOffset);

                if (uPoint.equals(vPoint) || visited.contains(vPoint) ||
                    excludePoints.contains(vPoint) || points.contains(vPoint)
                    || borderPixels.contains(vPoint)) {
                    continue;
                }

                visited.add(vPoint);
                
                int r = colorImg.getR(vX, vY);
                int g = colorImg.getG(vX, vY);
                int b = colorImg.getB(vX, vY);
                
                if (uLocalColors == null) {
                    
                    Set<PairInt> neighbors = getThe8NeighborPixelsWithin(
                        uPoint, points, width, height);
                    neighbors.add(uPoint); 
                        
                    uLocalColors = new GroupPixelColors(neighbors, colorImg, 
                        xOffset, yOffset);                    
                }
                
                double contrast = uLocalColors.calcContrastToOther(r, g, b);
                
                double diffColor = b - gpc.getAvgBlue();

                double colorChange = diffColor/gpc.getStdDevBlue();

                if (skyIsRed) {
                    diffColor = r - gpc.getAvgRed();
                    colorChange = diffColor/gpc.getStdDevRed();
                }

                double diffContrast = contrast - gpc.getAvgContrast();

                double contrastChange = diffContrast/gpc.getStdDevContrast();

                boolean isBorder = false;
                
                float bFactor = 0.25f;
                float cFactor = 0.01f;
        
                if ((Math.abs(contrastChange) > cFactor) && (Math.abs(colorChange) > bFactor)) {
                    isBorder = true;
                }

                if (!isBorder || (
                    ((Math.abs(r - gpc.getAvgRed()) < firstSubtrTol) 
                    && (Math.abs(g - gpc.getAvgGreen()) < firstSubtrTol) &&
                    (Math.abs(b - gpc.getAvgBlue()) < firstSubtrTol)))) {
                    
                    nAdded++;

                    stack.add(vPoint);
                    points.add(vPoint);

                    mask.setValue(vX - xOffset, vY - yOffset, 0);

                    // store points which were just added to try an
                    //   add of all similar color neighbors
                    //   after 

                    addedSky.add(vPoint);
                    addedSkyColorToAdd.add(new PixelColors(r, g, b));
                }
               
            }
        }
        
        return addedSky;

    }

    /**
     * find sun colored points by searching entire image.
     * Note, that if the points returned are well fit by a circle
     * or a fraction of a circle that is calculable, one has a limit
     * on scale in the image (knowledge of the camera and lens are needed,
     * else such objects and points of interest can help provide them).
     * The sun's visible radius is the photosphere and that size is known.
     * @param colorImg
     * @param xOffset
     * @param yOffset
     * @param skyIsDarkGrey
     * @return 
     */
    public Set<PairInt> findSunPoints(ImageExt colorImg, int xOffset, 
        int yOffset, boolean skyIsDarkGrey) {
        
        Set<PairInt> set = new HashSet<PairInt>();
        
        SunColors sunColors = new SunColors();
        
        if (skyIsDarkGrey) {
            sunColors.useDarkSkiesLogic();
        }
        
        for (int col = 0; col < colorImg.getWidth(); col++) {
            for (int row = 0; row < colorImg.getHeight(); row++) {
                
                int idx = colorImg.getInternalIndex(col, row);

                if (sunColors.isSunCenterColor(colorImg, idx)) {
                    set.add(new PairInt(col - xOffset, row - yOffset));
                }
            }
        }
        
debugPlot(set, colorImg, xOffset, yOffset, 
"sunpoints");

        return set;
    }

    /**
     * search over the entire image for pixels that are rainbow colored.
     * 
     * @param colorImg
     * @param reflectedSunRemoved
     * @param xOffset
     * @param yOffset
     * @param skyIsDarkGrey
     * @return rainbowPoints pixels from the entire image containing rainbow 
     * colors.
     */
    Set<PairInt> findRainbowColoredPoints(ImageExt colorImg, 
        Set<PairInt> reflectedSunRemoved,
        int xOffset, int yOffset, boolean skyIsDarkGrey) {
        
        AbstractSkyRainbowColors colors = skyIsDarkGrey ?
            new DarkSkyRainbowColors() : new BrightSkyRainbowColors();
        
        Set<PairInt> set = new HashSet<PairInt>();
        
        for (int col = 0; col < colorImg.getWidth(); col++) {
            for (int row = 0; row < colorImg.getHeight(); row++) {
                
                PairInt p = new PairInt(col - xOffset, row - yOffset);
                
                if (reflectedSunRemoved.contains(p)) {
                    continue;
                }
                
                int idx = colorImg.getInternalIndex(col, row);
                
                if (colors.isInRedThroughPurplishRed(colorImg, idx)) {
                    
                    set.add(p);

                } else if (colors.isInOrangeRed(colorImg, idx)) {
                    
                    set.add(p);
                
                } else if (colors.isInGreenishYellowOrange(colorImg, idx)) {
                 
                    set.add(p);
                }
            }
        }
        
        return set;
    }
    
    Set<PairInt> findRainbowColoredPointsConnectedToSkyPoints(
        Set<PairInt> skyPoints, Set<PairInt> reflectedSunRemoved,
        ImageExt colorImg, int xOffset, int yOffset, boolean skyIsDarkGrey) {
        
        Set<PairInt> rainbowPoints = new HashSet<PairInt>();
        
        AbstractSkyRainbowColors rainbowColors = skyIsDarkGrey ?
            new DarkSkyRainbowColors() : new BrightSkyRainbowColors();
        
        int width = colorImg.getWidth();
        int height = colorImg.getHeight();
        
        int[] dxs = new int[]{-1, -1,  0,  1, 1, 1, 0, -1};
        int[] dys = new int[]{ 0, -1, -1, -1, 0, 1, 1,  1};
        
        java.util.Stack<PairInt> stack = new java.util.Stack<PairInt>();

        for (PairInt p : skyPoints) {
            
            int x = p.getX() + xOffset;
            int y = p.getY() + yOffset;
            
            int idx = colorImg.getInternalIndex(x, y);

            if (rainbowColors.isInGreenishYellowOrange(colorImg, idx)
                || rainbowColors.isInOrangeRed(colorImg, idx)
                || rainbowColors.isInRedThroughPurplishRed(colorImg, idx)) {
                
                stack.add(p);
                
                rainbowPoints.add(p);
            }
        }
        
        Set<PairInt> visited = new HashSet<PairInt>();
        visited.add(stack.peek());

        while (!stack.isEmpty()) {

            PairInt uPoint = stack.pop();

            int uX = uPoint.getX() + xOffset;
            int uY = uPoint.getY() + yOffset;

            //(1 + frac)*O(N) where frac is the fraction added back to stack

            for (int k = 0; k < dxs.length; k++) {

                int dx = dxs[k];
                int dy = dys[k];

                int vX = uX + dx;
                int vY = uY + dy;

                if ((vX < 0) || (vX > (width - 1)) || (vY < 0) || 
                    (vY > (height - 1))) {
                    continue;
                }

                PairInt vPoint = new PairInt(vX - xOffset, vY - yOffset);

                if (visited.contains(vPoint) || rainbowPoints.contains(vPoint)
                    || reflectedSunRemoved.contains(vPoint)) {
                    continue;
                }

                visited.add(vPoint);

                int idx = colorImg.getInternalIndex(vX, vY);
                        
                if (rainbowColors.isInGreenishYellowOrange(colorImg, idx) ||
                    rainbowColors.isInOrangeRed(colorImg, idx) ||
                    rainbowColors.isInRedThroughPurplishRed(colorImg, idx)) {
                    stack.add(vPoint);
                    rainbowPoints.add(vPoint);
                }
            }
        }
        
        if (rainbowPoints.isEmpty()) {
            return rainbowPoints;
        }
        
        if (rainbowPoints.size() < 12) {
            return new HashSet<PairInt>();
        }
        
        return rainbowPoints;
    }

    boolean skyIsDarkGrey(GroupPixelColors allSkyColor) {
        
        float rDiffG = Math.abs(allSkyColor.getAvgRed() -
            allSkyColor.getAvgGreen());

        float gDiffB = Math.abs(allSkyColor.getAvgGreen() -
            allSkyColor.getAvgBlue());

        float rDivB = Math.abs(allSkyColor.getAvgRed()/
            allSkyColor.getAvgBlue());

        if ((rDivB >= 0.85) && (rDivB <= 1.15) && (rDiffG < 20) &&
            (gDiffB < 20) && (allSkyColor.getAvgGreen() < 125)) {

            return true;
        }
        
        return false;
    }

    private Set<PairInt> getThe8NeighborPixelsWithin(PairInt uPoint, 
        Set<PairInt> points, int width, int height) {
        
        return getTheNeighborPixels(uPoint, points, width, height, 1);
    }
    
    private Set<PairInt> getThe24NeighborPixelsWithin(PairInt uPoint, 
        Set<PairInt> points, int width, int height) {
        
        return getTheNeighborPixels(uPoint, points, width, height, 2);
    }

    Set<PairInt> getTheNeighborPixels(PairInt uPoint, 
        Set<PairInt> points, int width, int height, int radius) {
        
        int uX = uPoint.getX();
        int uY = uPoint.getY();
        
        Set<PairInt> set = new HashSet<PairInt>();
        
        for (int vX = (uX - radius); vX <= (uX + radius); vX++) {

            if ((vX < 0) || (vX > (width - 1))) {
                continue;
            }

            for (int vY = (uY - radius); vY <= (uY + radius); vY++) {
                
                if ((vY < 0) || (vY > (height - 1))) {
                    continue;
                }
                
                PairInt vPoint = new PairInt(vX, vY);
                
                if (uPoint.equals(vPoint) || !points.contains(vPoint)) {
                    continue;
                }
                
                set.add(vPoint);
            }
        }
        
        return set;
    }

    private Set<PairInt> getThe8NeighborPixelsWithin(PairInt uPoint, 
        Set<PairInt> points, Set<PairInt> excludePoints, int width, int height) {
        
        int uX = uPoint.getX();
        int uY = uPoint.getY();
        
        Set<PairInt> set = new HashSet<PairInt>();
        
        for (int vX = (uX - 1); vX <= (uX + 1); vX++) {

            if ((vX < 0) || (vX > (width - 1))) {
                continue;
            }

            for (int vY = (uY - 1); vY <= (uY + 1); vY++) {
                
                if ((vY < 0) || (vY > (height - 1))) {
                    continue;
                }
                
                PairInt vPoint = new PairInt(vX, vY);
                
                if (uPoint.equals(vPoint) || !points.contains(vPoint)
                    || !excludePoints.contains(vPoint)) {
                    continue;
                }
                
                set.add(vPoint);
            }
        }
        
        return set;
    }

    /**
     * a debug method to print scatter diagrams and histograms if needed
     * of the colors in search of best indicator that blue or red skies
     * have many clouds.
     * @param points
     * @param colorImg
     * @param mask 
     */
    private void plotSkyColor(Set<PairInt> points, ImageExt colorImg, 
        GreyscaleImage mask) {

        int n = points.size();
        
        int xOffset = mask.getXRelativeOffset();
        int yOffset = mask.getYRelativeOffset();
        
        CIEChromaticity cieC = new CIEChromaticity();
        
        float[] cieX = new float[n];
        float[] cieY = new float[n];
        float[] hue = new float[n];
        float[] saturation = new float[n];
        float[] brightness = new float[n];
        float[] r = new float[n];
        float[] g = new float[n];
        float[] b = new float[n];
        
        int i = 0;
        for (PairInt p : points) {
            int x = p.getX() + xOffset;
            int y = p.getY() + yOffset;
            int rr = colorImg.getR(x, y);
            int gg = colorImg.getG(x, y);
            int bb = colorImg.getB(x, y);
            
            r[i] = rr;
            g[i] = gg;
            b[i] = bb;
            
            float[] cieXY = cieC.rgbToXYChromaticity(rr, gg, bb);
            cieX[i] = cieXY[0];
            cieY[i] = cieXY[1];
            
            float[] hsb = new float[3];
            Color.RGBtoHSB(rr, gg, bb, hsb);
            hue[i] = hsb[0];
            saturation[i] = hsb[1];
            brightness[i] = hsb[2];
            
            i++;
        }
        
        double t0 = System.currentTimeMillis();
        double t = t0 - ((int)(t0/1.E9)) * 1E9;
        int plotNumber = (int)t;
        
        try {
            ScatterPointPlotterPNG plotter = new ScatterPointPlotterPNG();
            
            plotter.plot(0.0f, 0.8f, 0.0f, 0.9f, 
                cieX, cieY, "CIE X vs Y", "CIEX", "CIEY");
            plotter.writeFile(Integer.valueOf(plotNumber));
            
            plotter = new ScatterPointPlotterPNG();
            plotter.plot(0.0f, 256f, 0.0f, 1.1f, 
                b, brightness, "B vs brightness", "B", "brightness");
            plotter.writeFile(Integer.valueOf(plotNumber + 1));
            
            plotter = new ScatterPointPlotterPNG();
            plotter.plot(0.0f, 256f, 0.0f, 1.1f, 
                r, brightness, "R vs brightness", "R", "brightness");
            plotter.writeFile(Integer.valueOf(plotNumber + 2));
            
            plotter = new ScatterPointPlotterPNG();
            plotter.plot(0.0f, 256f, 0.0f, 256f, 
                b, r, "B vs R", "B", "R");
            plotter.writeFile(Integer.valueOf(plotNumber + 3));
            
            plotter = new ScatterPointPlotterPNG();
            plotter.plot(0.0f, 1.1f, 0.0f, 1.1f, 
                hue, saturation, "hue vs saturation", "hue", "saturation");
            plotter.writeFile(Integer.valueOf(plotNumber + 4));
            
            plotter = new ScatterPointPlotterPNG();
            plotter.plot(0.0f, 0.8f, 0.0f, 1.1f, 
                cieX, brightness, "CIE X vs brightness", "CIEX", "brightness");
            plotter.writeFile(Integer.valueOf(plotNumber + 5));
            
            plotter = new ScatterPointPlotterPNG();
            plotter.plot(0.0f, 0.9f, 0.0f, 1.1f, 
                cieY, brightness, "CIE Y vs brightness", "CIEY", "brightness");
            plotter.writeFile(Integer.valueOf(plotNumber + 6));
                        
            // for hue histograms, because it's space is formed in 0 to 360 degrees,
            // need to wrap around the bins.
            // for that reason, will look for an empty bin and cut and merge
            // the wrap
            
            HistogramHolder hueHist = Histogram.createSimpleHistogram(
                0.f, 1.0f, 10, hue, Errors.populateYErrorsBySqrt(hue));
                        
            int[] yh = Arrays.copyOf(hueHist.getYHist(), hueHist.getYHist().length);
            int zeroBinIdx = -1;
            for (int ii = 0; ii < yh.length; ii++) {
                if (yh[ii] == 0) {
                    zeroBinIdx = ii;
                    break;
                }
            }
            if (zeroBinIdx > -1) {
                int[] yShifted = new int[yh.length];
                float[] yShiftedF = new float[yh.length];
                int count = 0;
                for (int ii = (zeroBinIdx + 1); ii < yh.length; ii++) {
                    yShifted[count] = yh[ii];
                    yShiftedF[count] = yh[ii];
                    count++;
                }
                for (int ii = 0; ii <= zeroBinIdx; ii++) {
                    yShifted[count] = yh[ii];
                    yShiftedF[count] = yh[ii];
                    count++;
                }
                hueHist.setYHist(yShifted);
                hueHist.setYHistFloat(yShiftedF);
            }
            
            float fwhmHue = Histogram.measureFWHMOfStrongestPeak(hueHist);
            
            hueHist.plotHistogram("hue shifted by " + (zeroBinIdx + 0) + " bins",
                plotNumber);
            
            log.info("fwhm hue=" + fwhmHue);

            HistogramHolder saturationHist = Histogram.createSimpleHistogram(
                0.f, 1.0f, 10, 
                saturation, Errors.populateYErrorsBySqrt(saturation));
            
            saturationHist.plotHistogram("saturation", plotNumber + 1);
            
            float[] fwhmSaturation = Histogram.measureFWHMOfAllPeaks(
                saturationHist, 0.1f);
            
            log.info("fwhm saturation=" + Arrays.toString(fwhmSaturation));
            
        } catch (IOException e) {
            
            log.severe(e.getMessage());
        }
    }

    private PixelColors getAveragedColorOfDifference(Set<PairInt> pointsToExclude, 
        Image colorImg, int xOffset, int yOffset) {
        
        long rSum = 0;
        long gSum = 0;
        long bSum = 0;
        int n = 0;
        
        for (int col = 0; col < colorImg.getWidth(); col++) {
            for (int row = 0; row < colorImg.getHeight(); row++) {
                
                PairInt p = new PairInt(col - xOffset, row - yOffset);
                
                if (pointsToExclude.contains(p)) {
                    continue;
                }
                
                rSum += colorImg.getR(col, row);
                gSum += colorImg.getG(col, row);
                bSum += colorImg.getB(col, row);
                n++;
            }
        }
        
        int rAvg = (int)Math.round((double)rSum/(double)n);
        int gAvg = (int)Math.round((double)gSum/(double)n);
        int bAvg = (int)Math.round((double)bSum/(double)n);
        
        PixelColors avgColor = new PixelColors(rAvg, gAvg, bAvg);
        
        return avgColor;
    }

    private void removeBinFactorArtifacts(Set<PairInt> points, int binFactor, 
        int pointsImageWidth, int pointsImageHeight,
        int pointsImageXOffset, int pointsImageYOffset) {
        
        if (binFactor == 0) {
            return;
        }
        
        int yLimit = (pointsImageHeight + pointsImageYOffset - 1);
        
        List<PairInt> remove = new ArrayList<PairInt>();
        
        // making vertical corrections only
        for (PairInt p : points) {
            
            int x = p.getX();
            int y = p.getY();
            
            boolean isProtruding = true;
            
            for (int i = 1; i < 4; i++) {
                                
                int y2 = y + i;
                
                if (y2 > yLimit) {
                    isProtruding = false;
                    break;
                }
                
                PairInt p2 = new PairInt(x, y2);
                if (!points.contains(p2)) {
                    isProtruding = false;
                    break;
                }
            }
            
            if (isProtruding) {
                
                // check that there is no point at i=5, enclosing the protrusion
                PairInt p2 = new PairInt(x, y + 5);
                if ((y < yLimit) && points.contains(p2)) {
                    continue;
                }
                for (int i = 0; i < 4; i++) {
                    int y2 = y + i;
                    p2 = new PairInt(x, y2);
                    remove.add(p2);
                }
            }
        }
        
        for (PairInt p : remove) {
            points.remove(p);
        }
    }
                    
    private Set<PairInt> removeSetsWithNonCloudColors(List<PairIntArray> 
        zeroPointLists, Image originalColorImage, GreyscaleImage theta,
        boolean addAlongX, int addAmount) {
        
        int colorLimit = 100;
        
        List<Integer> remove = new ArrayList<Integer>();
        
        int xOffset = theta.getXRelativeOffset();
        int yOffset = theta.getYRelativeOffset();
        
        for (int gId = 0; gId < zeroPointLists.size(); gId++) {
            
            PairIntArray points  = zeroPointLists.get(gId);
            
            int nBelowLimit = 0;
            int nBelowLimit2 = 0;
            int nBelowLimit3 = 0;
            
            for (int i = 0; i < points.getN(); i++) {
                
                int x = points.getX(i);
                int y = points.getY(i);
                
                int ox = x + xOffset;
                int oy = y + yOffset;
                
                if (addAlongX) {
                    ox += addAmount;
                } else {
                    oy += addAmount;
                }

                if ((ox < 0) || (ox > (originalColorImage.getWidth() - 1))) {
                    continue;
                }
                if ((oy < 0) || (oy > (originalColorImage.getHeight() - 1))) {
                    continue;
                }
                
                int r = originalColorImage.getR(ox, oy);
                int g = originalColorImage.getG(ox, oy);
                int b = originalColorImage.getB(ox, oy);
                
                if ((r < colorLimit) && (b < colorLimit) && (g < colorLimit)) {
                    nBelowLimit++;
                }
                // tracking color of sand
                if ((Math.abs((r/255.) - 0.35) < .10) 
                    && (Math.abs((g/255.) - 0.35) < .10) 
                    && (Math.abs((b/255.) - 0.35) < .10)
                    && (b < g)) {
                    nBelowLimit2++;
                }
                // TODO: there has to be a better way to filter green and cyan
                // green colors
                if ((Math.abs((r/255.) - 0.55) < .05) 
                    && (Math.abs((g/255.) - 0.62) < .05) 
                    && (Math.abs((b/255.) - 0.25) < .05)
                    && (b < g)) {
                    nBelowLimit3++;
                } else if (((g - r) >= 20) && ((g - b) >= 20)) {
                    nBelowLimit3++;
                } else if ((Math.abs((r/255.) - 0.25) < .03) 
                    && (Math.abs((g/255.) - 0.38) < .05) 
                    && (Math.abs((b/255.) - 0.40) < .05)
                    && (g > 80)) {
                    nBelowLimit3++;
                } else if ((Math.abs((r/255.) - 0.41) < .03) 
                    && (Math.abs((g/255.) - 0.48) < .03) 
                    && (Math.abs((b/255.) - 0.49) < .03)
                    && (g > 80)) {
                    nBelowLimit3++;
                }
            }
            
            log.fine(gId + ") nBelowLimit=" + nBelowLimit
                + " (" + ((double)nBelowLimit/(double)points.getN()) + ")");
         
            if (((double)nBelowLimit/(double)points.getN()) > 0.5) {
                
                remove.add(Integer.valueOf(gId));
                
            } else if (((double)nBelowLimit2/(double)points.getN()) > 0.5) {
                
                remove.add(Integer.valueOf(gId));
                
            } else if (((double)nBelowLimit3/(double)points.getN()) > 0.5) {
                
                remove.add(Integer.valueOf(gId));
            }
        }
        
        Set<PairInt> removed = new HashSet<PairInt>();
        
        if (!remove.isEmpty()) {
            
            for (int i = (remove.size() - 1); i > -1; i--) {
                
                PairIntArray r = zeroPointLists.remove(remove.get(i).intValue());
                
                for (int ii = 0; ii < r.getN(); ii++) {
                    
                    int x = r.getX(ii);
                    int y = r.getY(ii);
                    PairInt p = new PairInt(x, y);

                    removed.add(p);
                }
            }
        }
        
        return removed;
    }
    
    Set<PairInt> findEmbeddedNonPoints(Set<PairInt> points,
        Set<PairInt> exclude0, Set<PairInt> exclude1, int imageMaxColumn) {
        
        Set<PairInt> embeddedPoints = new HashSet<PairInt>();
        
        PerimeterFinder finder = new PerimeterFinder();
        int[] rowMinMax = new int[2];
        Map<Integer, List<PairInt>> rowColRange = finder.find(points, 
            rowMinMax, imageMaxColumn, embeddedPoints);
        
        for (PairInt exclude : exclude0) {
            embeddedPoints.remove(exclude);
        }
        
        for (PairInt exclude : exclude1) {
            embeddedPoints.remove(exclude);
        }
        
        return embeddedPoints;
    }
    
    private void correctSkylineForSun(Set<PairInt> sunPoints, 
        Set<PairInt> skyPoints, Image colorImg, GreyscaleImage mask, 
        GreyscaleImage gradientXY) {
        
        MiscellaneousCurveHelper ch = new MiscellaneousCurveHelper();
        double[] xySunCen = ch.calculateXYCentroids(sunPoints);
        
        int h = mask.getHeight();
        int w = mask.getWidth();
        
        int srchRadius = 125;
        
        int maxValue = Integer.MIN_VALUE;
        
        Set<PairInt> pointsNearSun = new HashSet<PairInt>();
                
        for (int col = ((int)xySunCen[0] - srchRadius); 
            col < ((int)xySunCen[0] + srchRadius); col++) {
        
            for (int row = ((int)xySunCen[1] - srchRadius); 
                row < ((int)xySunCen[1] + srchRadius); row++) {
                
                if ((col < 0) || (col > (w - 1)) || (row < 0) || (row > (h - 1))) {
                    continue;
                }
                
                int v = gradientXY.getValue(col, row);
                
                if (v > maxValue) {
                    maxValue = v;
                }
                
                if (v > 0) {
                    pointsNearSun.add(new PairInt(col, row));
                }
            }
        }
        
        PairIntArray valueCounts = 
            Histogram.createADescendingSortByKeyArray(pointsNearSun, gradientXY);
        int countIsDoubledValue = maxValue;
        int lastCount = Integer.MAX_VALUE;
        for (int i = 1; i < valueCounts.getN(); i++) {
            int count = valueCounts.getY(i);
            if ((count/lastCount) >= 2) {
                countIsDoubledValue = valueCounts.getX(i);
                break;
            }
            lastCount = count;
        }
        
        assert(maxValue == valueCounts.getX(0));
        
        int valueTolerance = 10;
        if ((maxValue - countIsDoubledValue) > valueTolerance) {
            valueTolerance = maxValue - countIsDoubledValue;
        }
        
        int minX = Integer.MAX_VALUE;
        int minY = Integer.MAX_VALUE;
        int maxX = Integer.MIN_VALUE;
        int maxY = Integer.MIN_VALUE;
        
        Set<PairInt> set = new HashSet<PairInt>();
        
        for (int col = ((int)xySunCen[0] - srchRadius); 
            col < ((int)xySunCen[0] + srchRadius); col++) {
        
            for (int row = ((int)xySunCen[1] - srchRadius); 
                row < ((int)xySunCen[1] + srchRadius); row++) {
                
                if ((col < 0) || (col > (w - 1)) || (row < 0) || (row > (h - 1))) {
                    continue;
                }
                
                int v = gradientXY.getValue(col, row);
                
                if (Math.abs(v - maxValue) < valueTolerance) {
                    
                    set.add(new PairInt(col, row));
                    
                    if (col < minX) {
                        minX = col;
                    }
                    if (col > maxX) {
                        maxX = col;
                    }
                    if (row < minY) {
                        minY = row;
                    }
                    if (row > maxY) {
                        maxY = row;
                    }
                }
            }
        }

debugPlot(set, colorImg, mask.getXRelativeOffset(), mask.getYRelativeOffset(), 
"horizon_near_sun_before_thinning");
        
        // points in set now represent the skyline near the sun.
        // the points are a line widened by convolution so need to be thinned
        // to the line centroid.
        // ErosionFilter isn't currently able to provide a line that is centered,
        // it may be a line closer to one side than the other of the original
        // thick band of points.
        
        // so trying Zhang-Suen: 
        
        ZhangSuenLineThinner zsLT = new ZhangSuenLineThinner();
        zsLT.applyLineThinner(set, minX, maxX, minY, maxY);
        
        ch.populateGapsWithInterpolation(set);

debugPlot(set, colorImg, mask.getXRelativeOffset(), mask.getYRelativeOffset(), 
"horizon_near_sun");

        // if know skyline direction and any points are "below" set towards
        // the skyline, those should be removed from skyPoints
        
        removePointsUnderSkyline(set, sunPoints, skyPoints, mask.getWidth(), 
            mask.getHeight());

        skyPoints.addAll(set);
        skyPoints.addAll(sunPoints);
        
        // fill in the gaps
        //Set<PairInt> embeddedPoints = findEmbeddedNonPoints(skyPoints);
        //skyPoints.addAll(embeddedPoints);
        
        mask.fill(1);
        for (PairInt p : skyPoints) {
            int x = p.getX();
            int y = p.getY(); 
            mask.setValue(x, y, 0);
        }
        
debugPlot(skyPoints, colorImg, mask.getXRelativeOffset(), mask.getYRelativeOffset(), 
"horizon_near_sun_final");

    }
    
    private int[] determineChangeTowardsSkyline(Set<PairInt> skyline, 
        Set<PairInt> sunPoints, Set<PairInt> skyPoints) {
        
        //TODO:  needs alot of testing with images that are tilted and
        // have sun in them.
        
        // need to fit a line to skyline to get the slope and hence
        // know what is perpendicular to the skyline near the sun.
        
        int[] xsl = new int[skyline.size()];
        int[] ysl = new int[xsl.length];
        int i = 0;
        for (PairInt p : skyline) {
            xsl[i] = p.getX();
            ysl[i] = p.getY();
            i++;
        }
        
        MiscellaneousCurveHelper ch = new MiscellaneousCurveHelper();
            
        double[] skyXYCen = ch.calculateXYCentroids(skyPoints);
            
        double[] sunXYCen = ch.calculateXYCentroids(sunPoints);
        
        LinearRegression lr = new LinearRegression();
        float[] yInterceptAndSlope = lr.calculateTheilSenEstimatorParams(
            xsl, ysl);
        
        int dSkylineX = 0;
        int dSkylineY = 0;
        
        /*
        inf
           1.7
         |    1.0
         |   /    0.57
         | /     
        @ ----- 0        
        */
        
        if (Math.abs(yInterceptAndSlope[1]) < 0.57) {
            // ~ horizontal line, less than 30 degrees from horiz
            if (skyXYCen[1] > sunXYCen[1]) {
                dSkylineY = -1;
            } else {
                dSkylineY = 1;
            }
        } else if (Math.abs(yInterceptAndSlope[1]) < 1.73) {
            // ~ diagonal line, between 30 degrees and 60 degrees from horiz
            if (skyXYCen[0] > sunXYCen[0]) {
                dSkylineX = -1;
            } else {
                dSkylineX = 1;
            }
            if (skyXYCen[1] > sunXYCen[1]) {
                dSkylineY = -1;
            } else {
                dSkylineY = 1;
            }
        } else  {
            // approx w/ a vertical line
            if (skyXYCen[0] > sunXYCen[0]) {
                dSkylineX = -1;
            } else {
                dSkylineX = 1;
            }
        }
        
        return new int[]{dSkylineX, dSkylineY};
    }

    private void removePointsUnderSkyline(Set<PairInt> skyline, 
        Set<PairInt> sunPoints, Set<PairInt> skyPoints, 
        int width, int height) {
                
        //TODO:  needs alot of testing with images that are tilted and
        // have sun in them.
        
        int[] dSkylineXY = determineChangeTowardsSkyline(skyline, sunPoints,
            skyPoints);
        
        int dSkylineX = Math.round(dSkylineXY[0]);
        int dSkylineY = Math.round(dSkylineXY[1]);
        
        if ((dSkylineX == 0) && (dSkylineY == 0)) {
            // do nothing
            
        } else if ((dSkylineX >= 0) && (dSkylineY >= 0)) {
            
            for (PairInt s : skyline) {
                int vX = s.getX() + dSkylineX;
                int vY = s.getY() + dSkylineY;
                while ((vX < (width - 1)) && (vY < (height - 1))) {
                    PairInt p = new PairInt(vX, vY);
                    skyPoints.remove(p);
                    vX += dSkylineX;
                    vY += dSkylineY;
                }
            }
            
        } else if ((dSkylineX == -1) && (dSkylineY >= 0)) {
            
            for (PairInt s : skyline) {
                int vX = s.getX() + dSkylineX;
                int vY = s.getY() + dSkylineY;
                while ((vX > -1) && (vY < (height - 1))) {
                    PairInt p = new PairInt(vX, vY);
                    skyPoints.remove(p);
                    vX += dSkylineX;
                    vY += dSkylineY;
                }
            }
        
        } else if ((dSkylineX >= 0) && (dSkylineY == -1)) {
            
            for (PairInt s : skyline) {
                int vX = s.getX() + dSkylineX;
                int vY = s.getY() + dSkylineY;
                while ((vX < (width - 1)) && (vY > -1)) {
                    PairInt p = new PairInt(vX, vY);
                    skyPoints.remove(p);
                    vX += dSkylineX;
                    vY += dSkylineY;
                }
            }
            
        } else if ((dSkylineX == -1) && (dSkylineY == -1)) {
            
            for (PairInt s : skyline) {
                int vX = s.getX() + dSkylineX;
                int vY = s.getY() + dSkylineY;
                while ((vX > -1) && (vY > -1)) {
                    PairInt p = new PairInt(vX, vY);
                    skyPoints.remove(p);
                    vX += dSkylineX;
                    vY += dSkylineY;
                }
            }
        
        }
    }
    
    GreyscaleImage filterAndExtractSkyFromGradient(ImageExt colorImg, 
        GreyscaleImage thetaImg, GreyscaleImage gXYImg, int binFactor,
        Set<PairInt> outputSkyPoints, RemovedSets removedSets) {
        
        ImageProcessor imageProcessor = new ImageProcessor();

        GreyscaleImage originalThetaImg = thetaImg.copyImage();
        
        if (binFactor > 1) {
            thetaImg = imageProcessor.binImage(thetaImg, binFactor);
            colorImg = imageProcessor.binImage(colorImg, binFactor);
            gXYImg = imageProcessor.binImage(gXYImg, binFactor);
        }

        List<PairIntArray> zeroPointLists = getSortedContiguousZeros(thetaImg);
        
        if (zeroPointLists.isEmpty()) {
            return null;
        }
      
        Set<PairInt> removedNonCloudColors = removeSetsWithNonCloudColors(
            zeroPointLists, colorImg, thetaImg, true, 0);
        
        if (zeroPointLists.isEmpty()) {
            return null;
        }
        
        //TODO:  assumes that largest smooth component of image is sky.
        // if sky is small and a foreground object is large and featureless
        // and not found as dark, this will fail. 
        // will adjust for that one day, possibly with color validation
        reduceToLargest(zeroPointLists);
 
        if (zeroPointLists.isEmpty()) {
            return null;
        }
        
        double[] avgYRGB = imageProcessor.calculateYRGB(zeroPointLists.get(0),
            colorImg,
            thetaImg.getXRelativeOffset(), thetaImg.getYRelativeOffset(),
            true, 0);
        
        int nBeforeHighContrastRemoval = count(zeroPointLists);
        
        Set<PairInt> highContrastRemoved = new HashSet<PairInt>();
        removeHighContrastPoints(zeroPointLists, colorImg, 
            thetaImg, avgYRGB[0], highContrastRemoved);
  
        if (zeroPointLists.isEmpty()) {
            return null;
        }

        int nAfterHighContrastRemoval = count(zeroPointLists);

        Set<PairInt> reflectedSunRemoved = removeReflectedSun(zeroPointLists, 
            colorImg, thetaImg);
 
        if (zeroPointLists.isEmpty()) {
            return null;
        }
        
        combine(zeroPointLists, outputSkyPoints);
       
        GreyscaleImage thresholdedGXY = extractSkyFromGradientXY(gXYImg, 
            outputSkyPoints, highContrastRemoved);
        
        removedSets.setHighContrastRemoved(highContrastRemoved);
        removedSets.setReflectedSunRemoved(reflectedSunRemoved);
        removedSets.setRemovedNonCloudColors(removedNonCloudColors);
        removedSets.setBeforeHighContrastRemoval(nBeforeHighContrastRemoval);
        removedSets.setAfterHighContrastRemoval(nAfterHighContrastRemoval);
        
        if (binFactor > 1) {
            
            Set<PairInt> tmp = imageProcessor.unbinZeroPointLists(
                outputSkyPoints, binFactor);
            outputSkyPoints.clear();
            outputSkyPoints.addAll(tmp);
            
            tmp = imageProcessor.unbinZeroPointLists(
                removedNonCloudColors, binFactor);
            removedSets.setRemovedNonCloudColors(tmp);
            
            tmp = imageProcessor.unbinZeroPointLists(
                highContrastRemoved, binFactor);
            removedSets.setHighContrastRemoved(tmp);
            
            tmp = imageProcessor.unbinZeroPointLists(
                reflectedSunRemoved, binFactor);
            removedSets.setReflectedSunRemoved(tmp);
            
            thresholdedGXY = imageProcessor.unbinMask(gXYImg, 
                binFactor, originalThetaImg);
            
            if (!highContrastRemoved.isEmpty()) {
                                
                removeBinFactorArtifacts(outputSkyPoints, binFactor,
                    originalThetaImg.getWidth(), originalThetaImg.getHeight(),
                    originalThetaImg.getXRelativeOffset(), 
                    originalThetaImg.getYRelativeOffset());
            }
        }
        
        return thresholdedGXY;
    }

    protected void addIfSimilarToSky(Set<PairInt> embeddedPoints, 
        Set<PairInt> skyPoints, Set<PairInt> highContrastRemoved, 
        ImageExt originalColorImage, GreyscaleImage mask, 
        HistogramHolder[] brightnessHistogram, 
        GroupPixelColors[] skyPartitionedByBrightness) {
        
        Set<PairInt> excludePoints = highContrastRemoved;
        int xOffset = mask.getXRelativeOffset();
        int yOffset = mask.getYRelativeOffset();
        
        boolean pointsAreEmbeddedInSky = true;
        
        double rDivB = skyPartitionedByBrightness[2].getAvgRed() /
            skyPartitionedByBrightness[2].getAvgBlue();

        boolean skyIsRed = (rDivB > 1);

        boolean skyIsPurple = skyIsRed && (rDivB < 1.5) &&
            (Math.abs(
                skyPartitionedByBrightness[2].getAvgGreen() -
                skyPartitionedByBrightness[2].getAvgBlue())
                <
                0.1 * skyPartitionedByBrightness[2].getAvgBlue()
            );
        boolean hasDarkGreyClouds = false;
//only perform k=0 if the sky has a narrow dark section
        boolean useKEqualsZero =
            (skyPartitionedByBrightness[0].getStdDevGreen() < 10) &&
            (skyPartitionedByBrightness[0].getStdDevBlue() < 10) &&
            (skyPartitionedByBrightness[0].getStdDevContrast() < 10)
            && (skyPartitionedByBrightness[0].getStdDevColorDiff()
            < 10);
        if (useKEqualsZero) {
            if (skyIsRed) {
                useKEqualsZero = useKEqualsZero &&
                (skyPartitionedByBrightness[0].getStdDevRed() < 20);
            } else {
                useKEqualsZero = useKEqualsZero &&
                (skyPartitionedByBrightness[0].getStdDevRed() < 10);
            }
        } else {
            if (!skyIsRed) {
                // if the sky is grey, may have dark clouds
                if (skyIsDarkGrey(skyPartitionedByBrightness[0])) {
                    useKEqualsZero = true;
                    hasDarkGreyClouds = true;
                }
            }
        }
        
        Set<PairInt> outputCloudPoints = new HashSet<PairInt>();
        Stack<PairInt> outputStack = new java.util.Stack<PairInt>();
        
        for (PairInt embeddedPoint : embeddedPoints) {

            if (excludePoints.contains(embeddedPoint)) {

                continue;
            }

            filterToAddCloudPixel(excludePoints, skyPoints, 
                originalColorImage, xOffset, yOffset, embeddedPoint,
                brightnessHistogram, skyPartitionedByBrightness,
                skyIsRed, skyIsPurple, hasDarkGreyClouds, useKEqualsZero,
                pointsAreEmbeddedInSky,
                outputCloudPoints, outputStack
            );
        }
        skyPoints.addAll(outputCloudPoints);
    }

    public class RemovedSets {
        private Set<PairInt> removedNonCloudColors = null;
        private Set<PairInt> highContrastRemoved = null;
        private Set<PairInt> reflectedSunRemoved = null;
        private int nBeforeHighContrastRemoval = Integer.MIN_VALUE;
        private int nAfterHighContrastRemoval = Integer.MIN_VALUE;
        
        public RemovedSets(){};
        public void setRemovedNonCloudColors(Set<PairInt> points) {
            this.removedNonCloudColors = points;
        }
        public void setHighContrastRemoved(Set<PairInt> points) {
            this.highContrastRemoved = points;
        }
        public void setReflectedSunRemoved(Set<PairInt> points) {
            this.reflectedSunRemoved = points;
        }
        public Set<PairInt> getRemovedNonCloudColors() {
            return removedNonCloudColors;
        }
        public Set<PairInt> getHighContrastRemoved() {
            return highContrastRemoved;
        }
        public Set<PairInt> getReflectedSunRemoved() {
            return reflectedSunRemoved;
        }

        private void setBeforeHighContrastRemoval(int nBeforeHighContrastRemoval) {
            this.nBeforeHighContrastRemoval = nBeforeHighContrastRemoval;
        }

        private void setAfterHighContrastRemoval(int nAfterHighContrastRemoval) {
            this.nAfterHighContrastRemoval = nAfterHighContrastRemoval;
        }
        
        private int getNBeforeHighContrastRemoval() {
            return nBeforeHighContrastRemoval;
        }

        private int getNAfterHighContrastRemoval() {
            return nAfterHighContrastRemoval;
        }
    }
}
