package algorithms.imageProcessing;

import algorithms.CountingSort;
import algorithms.MultiArrayMergeSort;
import algorithms.compGeometry.EllipseHelper;
import algorithms.compGeometry.PerimeterFinder;
import algorithms.compGeometry.PointInPolygon;
import algorithms.imageProcessing.optimization.ANDedClauses;
import algorithms.imageProcessing.optimization.ColorData;
import algorithms.imageProcessing.optimization.SKYCONDITIONAL;
import algorithms.imageProcessing.optimization.SkylineANDedClauses;
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
import algorithms.util.PolygonAndPointPlotter;
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
    
    public static String debugName = "";
    public static void setDebugName(String name) {
        debugName = name;
    }
    
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
       
        RemovedSets removedSets = new RemovedSets();
        
        Set<PairInt> points = createBestSkyMaskPt1(theta, gradientXY, 
            originalColorImage, edgeSettings, outputSkyCentroid,
            removedSets);
       
        if (points.isEmpty()) {
            
            GreyscaleImage mask = theta.createWithDimensions();
               
            // return an image of all 1's
            mask.fill(1);
            
            return mask;
        }
        
        
        GroupPixelColors allSkyColor = new GroupPixelColors(points,
            originalColorImage, theta.getXRelativeOffset(), 
            theta.getYRelativeOffset());
        
        boolean skyIsDarkGrey = skyIsDarkGrey(allSkyColor);
        
        Set<PairInt> sunPoints = new HashSet<PairInt>();
        double[] ellipFitParams = findSunPoints(
            originalColorImage, 
            theta.getXRelativeOffset(), theta.getYRelativeOffset(),
            skyIsDarkGrey, sunPoints);
        
        double density = (ellipFitParams != null) ?
            (sunPoints.size()/(2*Math.PI*ellipFitParams[2]*ellipFitParams[3]))
            : 0;
                
        if (density < 0.5) {
            sunPoints.clear();
        } else if (ellipFitParams == null) {
            sunPoints.clear();
        } else if ((ellipFitParams[2]/ellipFitParams[3])> 6) {
            sunPoints.clear();
        }

if (!sunPoints.isEmpty()) {
System.out.println(" density=" + density +
" ellipFitParams[2]/ellipFitParams[3]=" 
+ ellipFitParams[2]/ellipFitParams[3] + " " + debugName);
debugPlot(sunPoints, originalColorImage, 
theta.getXRelativeOffset(), theta.getYRelativeOffset(), 
"sun_points_" + debugName);
}

        // should not see sun and rainbow in same image
        Set<PairInt> rainbowPoints = new HashSet<PairInt>();
        float[] rainbowCoeff = null;
        if (sunPoints.isEmpty()) {
            rainbowCoeff = findRainbowPoints(points, 
                removedSets.getReflectedSunRemoved(), 
                originalColorImage, 
                theta.getXRelativeOffset(), theta.getYRelativeOffset(),
                skyIsDarkGrey, rainbowPoints);
            
            if (rainbowPoints.size() < 10) {
                rainbowPoints.clear();
                rainbowCoeff = null;
            }
        }

        Set<PairInt> excludeRainbow = new HashSet<PairInt>();
        Hull rainbowHull = null;
        if (rainbowCoeff != null) {
            
            // note that this adds the outer points to the sky too
            rainbowHull = createRainbowHull(points, rainbowCoeff, rainbowPoints,
                originalColorImage, 
                theta.getXRelativeOffset(), theta.getYRelativeOffset());
            
            if (rainbowHull != null) {
                
                int minXHull = (int)MiscMath.findMin(rainbowHull.xHull);
                int maxXHull = (int)Math.ceil(MiscMath.findMax(rainbowHull.xHull));
                int minYHull = (int)MiscMath.findMin(rainbowHull.yHull);
                int maxYHull = (int)Math.ceil(MiscMath.findMax(rainbowHull.yHull));
                PointInPolygon p = new PointInPolygon();
                for (int col = minXHull; col <= maxXHull; col++) {
                    for (int row = minYHull; row <= maxYHull; row++) {
                        boolean in = p.isInSimpleCurve(col, row, 
                            rainbowHull.xHull, rainbowHull.yHull, 
                            rainbowHull.yHull.length);
                        if (in) {
                            excludeRainbow.add(new PairInt(col, row));
                        }
                    }
                }
                
                if (!excludeRainbow.isEmpty()) {
                    // addRainbow to Hull, but only if there are sky points adjacent to hull
                    addRainbowToPoints(points, excludeRainbow, rainbowHull,
                        theta.getWidth() - 1, theta.getHeight() - 1);
                }
            }
        }

        int nSkyPointsBeforeFindClouds = points.size();
        
        findClouds(points, excludeRainbow, originalColorImage, theta);  

        if (!excludeRainbow.isEmpty()) {
            // addRainbow to Hull, but only if there are sky points adjacent to hull
            addRainbowToPoints(points, excludeRainbow, rainbowHull,
                theta.getWidth() - 1, theta.getHeight() - 1);
        }
        
        GreyscaleImage mask = gradientXY.createWithDimensions();
        mask.fill(1);
        for (PairInt p : points) {
            int x = p.getX();
            int y = p.getY(); 
            mask.setValue(x, y, 0);
        }
      
debugPlot(points, originalColorImage, mask.getXRelativeOffset(), mask.getYRelativeOffset(), 
"before_add_embedded");
             
        addEmbeddedIfSimilarToSky(points, originalColorImage, mask,
            removedSets);

debugPlot(points, originalColorImage, mask.getXRelativeOffset(), mask.getYRelativeOffset(), 
"after_add_embedded");

        Set<PairInt> exclude = new HashSet<PairInt>();
        exclude.addAll(removedSets.getHighContrastRemoved());
        //exclude.addAll(rainbowPoints);
        exclude.addAll(sunPoints);

        log.info("number of sunPoints=" + sunPoints.size() + " "
        + "reflectedSunRemoved.size()=" + removedSets.getReflectedSunRemoved().size());

        boolean addEmbedded = false;
        
        //if (sunPoints.isEmpty()) {
        
            growForLowContrastLimits(points, exclude, originalColorImage, mask,
                determineBinFactorForSkyMask(theta.getNPixels()));

            addEmbedded = true;
debugPlot(points, originalColorImage, mask.getXRelativeOffset(), mask.getYRelativeOffset(), 
"after_low_contrast_grow");

        //}

        if (!sunPoints.isEmpty()) {
            
            correctSkylineForSun(sunPoints, points, originalColorImage, mask, 
                gradientXY);
            
            addEmbedded = true;
        }
        
        if (addEmbedded) {
            addEmbeddedIfSimilarToSky(points, originalColorImage, mask,
                removedSets);
        }
        
debugPlot(points, originalColorImage, mask.getXRelativeOffset(), 
mask.getYRelativeOffset(), "final");
        
        for (PairInt p : sunPoints) {
            int x = p.getX();
            int y = p.getY();            
            mask.setValue(x, y, 0);
        }
        
        points.addAll(sunPoints);
        
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();        
        double[] xycen = curveHelper.calculateXYCentroids(points);
        outputSkyCentroid.add((int)Math.round(xycen[0]), (int)Math.round(xycen[1]));

        ImageProcessor imageProcessor = new ImageProcessor();
        imageProcessor.removeSpurs(mask);
   
        return mask;
    }
    
    public Set<PairInt> createBestSkyMaskPt1(final GreyscaleImage theta,
        GreyscaleImage gradientXY, ImageExt originalColorImage, 
        CannyEdgeFilterSettings edgeSettings, PairIntArray outputSkyCentroid,
        RemovedSets removedSets) {
        
        if (theta == null) {
            throw new IllegalArgumentException("theta cannot be null");
        }
        
        int binFactor = determineBinFactorForSkyMask(theta.getNPixels());

        log.info("binFactor=" + binFactor);
        
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
       
        Set<PairInt> points = new HashSet<PairInt>();
        
        GreyscaleImage threshholdedGXY = filterAndExtractSkyFromGradient(
            originalColorImage, theta, gradientXY, binFactor, points,
            removedSets);
        
        if (threshholdedGXY == null) {
            
            return new HashSet<PairInt>();
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
           
        populatePixelExtColors(points, originalColorImage, theta);
        
        return points;
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
                
 System.out.println("x=" + h.getXHist()[i] + " f=" + f + " dy=" + dy);
                
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
     * If the range of contrast is large, return a contrast histogram,
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
       
        //TODO:  this may need adjustments.  top of sun rising over mountains..
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
     * @param outputRainbowPoints output variable to return found rainbow
     * points if any
     * @return polynomial fit coefficients to 
     * y = c0*1 + c1*x[i] + c2*x[i]*x[i].  this may be null if a fit wasn't
     * possible.
     */
    float[] findRainbowPoints(Set<PairInt> skyPoints, 
        Set<PairInt> reflectedSunRemoved,
        ImageExt colorImg, int xOffset, int yOffset, boolean skyIsDarkGrey,
        Set<PairInt> outputRainbowPoints) {

        Set<PairInt> rainbowPoints = findRainbowColoredPoints(colorImg, 
            reflectedSunRemoved, xOffset, yOffset, skyIsDarkGrey);
        
        if (rainbowPoints.isEmpty()) {
            return null;
        }
        
        if (rainbowPoints.size() < 12) {
            return null;
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
            return null;
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
                return null;
            }
            
            int w = colorImg.getWidth() - xOffset;
            int h = colorImg.getHeight() - yOffset;
            if (bestFittingPoints.size() > (0.3f * (float)(w*h))) {
                return null;
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
    
            log.fine("nGTX=" + nGTX + " nLTY=" + nLTY + " n=" 
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
                
                return null;
                
            } else if ((cieXRange < 0.08) && (cieYRange < 0.08)) {
                
                return null;
            }
            
            if ((nGTX > 10) || (nLTY > 10)) {
                
                float frac = (float)(nGTX + nLTY)/(float)bestFittingPoints.size();
                if (frac > 0.002) {
                    rainbowPoints.addAll(bestFittingPoints);
                }
            }
        }
        
        outputRainbowPoints.addAll(rainbowPoints);
        
        polyFitter.plotFit(coef, outputRainbowPoints, colorImg.getWidth(),
            colorImg.getHeight(), 234, "rainbow points");
        
System.out.println("rainbow points size=" + outputRainbowPoints.size());
 
        return coef;
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
        ImageExt originalColorImage, GreyscaleImage thetaImg) {
        
        if (true) {
            
            SkylineANDedClauses skylineANDedClauses = new SkylineANDedClauses();
            
            findClouds(skyPoints, excludePoints, originalColorImage, 
                thetaImg, skylineANDedClauses.getAllClauses());
            
            return;
        }
        
        int maskWidth = thetaImg.getWidth();
        int maskHeight = thetaImg.getHeight();
        
        ArrayDeque<PairInt> cloudQueue = new ArrayDeque<PairInt>(skyPoints.size());
        for (PairInt skyPoint : skyPoints) {
            cloudQueue.add(skyPoint);
        }
        
        Set<PairInt> visited = new HashSet<PairInt>();
        visited.add(cloudQueue.peek());
              
        int xOffset = thetaImg.getXRelativeOffset();
        int yOffset = thetaImg.getYRelativeOffset();
        
        GroupPixelColors allSkyColor = new GroupPixelColors(skyPoints,
            originalColorImage, xOffset, yOffset);
        
        double rDivB = allSkyColor.getAvgRed() / allSkyColor.getAvgBlue();
        boolean skyIsRed = (rDivB > 1);
       
        log.fine("==> r/b=" + rDivB
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
System.out.println("(" + (vX + xOffset) + "," + (vY + yOffset) + ") "
+   "contrastV=" + contrastV + " div=" + (contrastV/skyStDevContrast)
+ " colordiff=" + colorDiffV + " div=" + (colorDiffV/skyStDevColorDiff)
+ " diffciex=" + diffCIEX + " div=" + (diffCIEX/localSky.getStdDevCIEX())
+ " diffciey=" + diffCIEY + " div=" + (diffCIEY/localSky.getStdDevCIEY())
+ " r=" + rV + " g=" + gV + " b=" + bV + " bPercentV=" + bPercentV);
}

                if (isBrown) {
                    
                    // trying to skip over foreground such as land or sunset + water
                    
                    if (((colorDiffV/skyStDevColorDiff) > 15) && (saturation < 0.5)) {
                        
                        if ((saturation <= 0.4) && ((colorDiffV/skyStDevColorDiff) > 50) 
                            ) {
log.fine("FILTER 00");
                            continue;
                        }
                         if ((saturation <= 0.4) && 
                            ((Math.abs(contrastV)/Math.abs(skyStDevContrast)) > 10.)
                            ) {
log.fine("FILTER 01");                                                        
                            continue;
                        }
                    }
                }
                
                if ( // no contrast or color change, should be sky
                    (Math.abs(contrastV) < 0.01)
                    && (colorDiffV < 10)
                    && (diffCIEX < 0.009) && (diffCIEY < 0.009)) {
log.fine("FILTER 02");
                } else if (
                    (skyStDevContrast != 0.)
                    && ((Math.abs(contrastV)/skyStDevContrast) > 10.)
                    ) {
log.fine("FILTER 03");
                    continue;
                } else if (
                    (skyStDevContrast != 0.)
                    && ((Math.abs(contrastV) > 0.1) 
                    && ((Math.abs(contrastV)/skyStDevContrast) > (1.5 + (Math.abs(contrastV)-0.5)*(-2.0))))
                    && ((Math.abs(colorDiffV)/skyStDevColorDiff) > 2.5)
                    ) {
log.fine("FILTER 04");
                    continue;                     
                } else if (
                    (skyStDevContrast > 0.005)
                    && ((Math.abs(contrastV)/skyStDevContrast) > 5.)
                    && ((Math.abs(colorDiffV)/skyStDevColorDiff) > 5.)
                    && 
                    // if cieXY diffs are zero and within stdev, these are sky,
                    // so test for opposite for boundary pixel
                    (((diffCIEX > 0.001) 
                        || ((diffCIEX/localSky.getStdDevCIEX()) > 1.5)
                        || (diffCIEY > 0.001) 
                        || ((diffCIEY/localSky.getStdDevCIEY()) > 1.5)))                    
                    && (skyStDevColorDiff > 1.)
                    ) {
log.fine("FILTER 05");
                    continue;
                } else if (skyIsRed) {
                    
                    /* TODO:
                    The original filters have been replaced by a simpler
                    relationship which may need to be returned to this
                    after adding more tests.
                    if ((skyStDevContrast > 0.)
                        && ((Math.abs(colorDiffV)/skyStDevColorDiff) > 3.)
                        && (diffCIEX > 0.03) 
                        && ((diffCIEX/localSky.getStdDevCIEX()) > 3.)
                        ) {
                        continue;
                    } else if ((skyStDevContrast > 0.)
                        && ((Math.abs(colorDiffV)/skyStDevColorDiff) > 1.1)
                        && ((diffCIEX > 0.065) 
                        && ((diffCIEX/localSky.getStdDevCIEX()) > 1.1))
                        ) {
                        continue;
                     } else if ((skyStDevContrast > 0.)
                        && ((Math.abs(colorDiffV)/skyStDevColorDiff) > 10.)
                        && ((diffCIEX > 0.02) 
                        && ((diffCIEX/localSky.getStdDevCIEX()) > 1.5))
                        ) {
                        continue;
                    } else if ((skyStDevContrast > 0.)
                        && ((Math.abs(colorDiffV)/skyStDevColorDiff) > 3.)
                        && (diffCIEY > 0.03) 
                        && ((diffCIEY/localSky.getStdDevCIEY()) > 3.)
                        ) {
                        continue;
                    } else if ((skyStDevContrast > 0.)
                        && ((Math.abs(colorDiffV)/skyStDevColorDiff) > 1.1)
                        && ((diffCIEY > 0.03) 
                        && ((diffCIEY/localSky.getStdDevCIEY()) > 3.))
                        ) {
                        continue;                        
                    } else if ((skyStDevContrast > 0.)
                        && ((Math.abs(colorDiffV)/skyStDevColorDiff) > 10.)
                        && ((diffCIEY > 0.02) 
                        && ((diffCIEY/localSky.getStdDevCIEY()) > 1.5))
                        ) {
                        continue;
                     ...
                    */
                    
                    
                    if (
                        // contrast is defined by luma, so might be weak near
                        // skyline near sun for example
                        (skyStDevContrast != 0.)
                        && ((Math.abs(colorDiffV)/skyStDevColorDiff) > 15.*diffCIEX)
                        && (diffCIEX > 0.03) 
                        && ((diffCIEX/localSky.getStdDevCIEX()) > 15.*diffCIEX)
                        ) {
log.fine("FILTER 06");
                        continue;
                    } else if (
                        // contrast is defined by luma, so might be weak near
                        // skyline near sun for example
                        (skyStDevContrast != 0.)
                        && ((Math.abs(colorDiffV)/skyStDevColorDiff) > 15.*diffCIEY)
                        && (diffCIEY > 0.03) 
                        && ((diffCIEY/localSky.getStdDevCIEY()) > 15.*diffCIEY)
                        ) {
log.fine("FILTER 07");
                        continue;
                    } else if (skyStDevContrast == 0.) {
                        if (contrastV >= 0.) {
log.fine("FILTER 08");
                            doNotAddToStack = true;
                        }
                        
                    } else {
                        //TODO:  if there are sun points, need a zone of
                        // avoidance to not erode the foreground 
log.fine("FILTER 09");
                    }

                } else {
                    //blue filters
                    if (
                        (skyStDevContrast > 0.0)
                        && (contrastV < 0.05) 
                        && ((Math.abs(contrastV)/skyStDevContrast) > 1.5)
                        && ((Math.abs(colorDiffV)/skyStDevColorDiff) > 3.0)
                        && ((bPercentV < 0.37) && (bV > 199) && (gV > 199))
                        ) {
log.fine("FILTER 10");
                        continue;
                    } else if (
                        (skyStDevContrast > 0.0)
                        && (Math.abs(contrastV) > 0.05)
                        && ((Math.abs(colorDiffV)/skyStDevColorDiff) > 1.5)
                        && ((diffCIEX/localSky.getStdDevCIEX()) > 0.9) // ?
                        && ((diffCIEY/localSky.getStdDevCIEY()) > 0.9) // ?                       
                        && (skyStDevColorDiff > 0.)
                    ){
log.fine("FILTER 11");
                        continue;
                    } else if (
                    //(contrastV < 0.05)
                    (Math.abs(contrastV) < 0.05)
                    && (diffCIEX < 0.005) && (diffCIEY < 0.005)
                    &&
                        ((Math.abs(0.33 - rPercentV) > 0.08)
                        || (Math.abs(0.33 - gPercentV) > 0.03) 
                        || (Math.abs(0.33 - bPercentV) > 0.03)
                        || (gV > 199) || (bV > 199))
                    ) {
log.fine("FILTER 12");
                        continue;
                    } else {
log.fine("FILTER 13");
                        continue;
                    }
                }
                
                skyPoints.add(vPoint);

                if (!doNotAddToStack) {
                    cloudQueue.add(vPoint);
                }
            }
        }
    }

private boolean check(int vX, int xOffset, int vY, int yOffset) {
    //if (((vX + xOffset)>=281) && ((vX + xOffset)<=285) && ((vY + yOffset)==73)) {
    //    return true;
    //}
    return false;
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
     * and a set of those ANDed is present in each ANDedClauses[] clauses.
     * For example, 
     *    a clause to determine that a pixel is a border pixel because
     *    it has high contrast and high color difference might be:
     *       (((Math.abs(contrast)/1.0) > coeff0) && 
     *        ((Math.abs(contrast)/skyStDevContrast) > 1.5*coeff0*coeff1) &&
     *        ((Math.abs(colorDiff)/1.0) > coeff3) && 
     *        ((Math.abs(colorDiff)/skyStDevcolorDiff) > 1.5*coeff3*coeff4))
     *    (Note though that the use of 'coeff0*coeff1' has to be handled by the
     *    invoking code.  The ANDedClauses would actually just receive a
     *    coeff2, unaware that is it '1.5*coeff0*coeff1' to the invoker).
     * 
     * Then each ANDedClauses is 'OR'-ed with one another, trying to quickly
     * act if an ANDedClauses[i] evaluates to 'T'.  This is disjunctive 
     * normal form.
     * 
     *  if (is border pix)
     *      break
     *  else if (is border pix)
     *      break
     *  else if (is border pix)
     *      break
     *   
     *  a pixel making it to here looks like a sky pixel.
     * 
     * Note that internally, if the sky is red, the blue rules are removed
     * from the clauses and vice versa for blue skies.
     * </pre>
     * 
     * @param skyPoints
     * @param originalColorImage
     */
    public void findClouds(Set<PairInt> skyPoints, Set<PairInt> excludePoints,
        ImageExt originalColorImage, GreyscaleImage thetaImage, 
        ANDedClauses[] clauses) {
        
        int maskWidth = thetaImage.getWidth();
        int maskHeight = thetaImage.getHeight();
        
        ArrayDeque<PairInt> cloudQueue = new ArrayDeque<PairInt>(skyPoints.size());
        for (PairInt skyPoint : skyPoints) {
            cloudQueue.add(skyPoint);
        }

        Set<PairInt> visited = new HashSet<PairInt>();
        visited.add(cloudQueue.peek());
              
        int xOffset = thetaImage.getXRelativeOffset();
        int yOffset = thetaImage.getYRelativeOffset();
        
        GroupPixelColors allSkyColor = new GroupPixelColors(skyPoints,
            originalColorImage, xOffset, yOffset);
        
        double rDivB = allSkyColor.getAvgRed() / allSkyColor.getAvgBlue();
        boolean skyIsRed = (rDivB > 1);
        
        clauses = filterForSkyColor(clauses, skyIsRed);
       
        log.fine("==> r/b=" + rDivB
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
                
                if (visited.contains(vPoint) || /*skyPoints.contains(vPoint) ||*/
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
System.out.println("(" + (vX + xOffset) + "," + (vY + yOffset) + ") "
+   "contrastV=" + contrastV + " div=" + (contrastV/skyStDevContrast)
+ " colordiff=" + colorDiffV + " div=" + (colorDiffV/skyStDevColorDiff)
+ " diffciex=" + diffCIEX + " div=" + (diffCIEX/localSky.getStdDevCIEX())
+ " diffciey=" + diffCIEY + " div=" + (diffCIEY/localSky.getStdDevCIEY())
+ " r=" + rV + " g=" + gV + " b=" + bV + " bPercentV=" + bPercentV 
+ " isBrown=" + isBrown + " skyIsRes=" + skyIsRed);
int z = 1;
}

                if (isBrown) {
                    // trying to skip over foreground such as land or sunset + water                    
                    if ((colorDiffV > 15*skyStDevColorDiff) && (saturation < 0.5)) {
                        if ((saturation <= 0.4) && (colorDiffV > 50*skyStDevColorDiff) 
                            ) {
log.fine("FILTER 00");
if (check(vX, xOffset, vY, yOffset)) {
System.out.println("(" + (vX + xOffset) + "," + (vY + yOffset) + ") rejected by filter00");
}
                            continue;
                        }
                         if ((saturation <= 0.4) && 
                            (Math.abs(contrastV) > 10.*Math.abs(skyStDevContrast))
                            ) {
log.fine("FILTER 01"); 
if (check(vX, xOffset, vY, yOffset)) {
System.out.println("(" + (vX + xOffset) + "," + (vY + yOffset) + ") rejected by filter01 ");
}
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
                        ArrayPair yellowGreenOrange = 
                            cieC.getYellowishGreenThroughOrangePolynomial();
                        if (pInPoly.isInSimpleCurve(cieX, cieY,
                            yellowGreenOrange.getX(), yellowGreenOrange.getY(),
                            yellowGreenOrange.getX().length)) {
                            isSolarYellow = true;
                        }
                    }
                }
                
                ColorData data = new ColorData(skyIsRed,
                    contrastV, colorDiffV,
                    diffCIEX, diffCIEY,
                    skyStDevContrast, skyStDevColorDiff,
                    localSky.getStdDevCIEX(), localSky.getStdDevCIEY(), 
                    rV, gV, bV
                );
                
                if (// no contrast or color change, should be sky
                    (Math.abs(contrastV) < 0.015)
                    && (Math.abs(colorDiffV) < 10)
                    && (diffCIEX < 0.009) && (diffCIEY < 0.009)) {
                    // this is a sky point
log.fine("FILTER 02");

                //} else if (isSolarYellow) {
                    // this is a sky point
//log.fine("FILTER 03");
                } else {
                    
                    // evaluate clauses that evaluate to 'T' when the pixel 
                    // looks like a border (non-sky) pixel
                    boolean isNotSky = false;

                    for (ANDedClauses clause : clauses) {
                        if (clause.evaluate(data)) {
                            isNotSky = true;
                            break;
                        }
                    }
            
                    if (isNotSky) {
                        continue;
                    }
                }
                
                skyPoints.add(vPoint);

                boolean doNotAddToStack = false;
                if (skyIsRed && (skyStDevContrast == 0.) && (contrastV >= 0.)) {
                    doNotAddToStack = true;
                }
                if (!doNotAddToStack) {
                    cloudQueue.add(vPoint);
                }
            }
        }

    }
   
    public static void getEmbeddedAndBorderPoints(Set<PairInt> skyPoints,
        int width, int height, Set<PairInt> outputEmbeddedGapPoints,
        Set<PairInt> outputBorderPoints) {
        
        PerimeterFinder perimeterFinder = new PerimeterFinder();

        int imageMaxColumn = width - 1;
        int imageMaxRow = height - 1;
       
        int[] skyRowMinMax = new int[2];
        
        Map<Integer, List<PairInt>> skyRowColRanges = perimeterFinder.find(
            skyPoints, skyRowMinMax, imageMaxColumn, 
            outputEmbeddedGapPoints);

        // update the perimeter for "filling in" embedded points
        perimeterFinder.updateRowColRangesForAddedPoints(skyRowColRanges, 
            skyRowMinMax, imageMaxColumn, outputEmbeddedGapPoints);
        
        Set<PairInt> borderPixels = perimeterFinder.getBorderPixels(
            skyRowColRanges, skyRowMinMax, imageMaxColumn, imageMaxRow);
        
        outputBorderPoints.addAll(borderPixels);
        
    }
    
    public static void getEmbeddedPoints(Set<PairInt> skyPoints,
        int width, int height, Set<PairInt> outputEmbeddedPoints) {
        
        PerimeterFinder perimeterFinder = new PerimeterFinder();

        int imageMaxColumn = width - 1;
        int imageMaxRow = height - 1;
       
        int[] skyRowMinMax = new int[2];
        
        Map<Integer, List<PairInt>> skyRowColRanges = perimeterFinder.find(
            skyPoints, skyRowMinMax, imageMaxColumn, 
            outputEmbeddedPoints);
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
        Set<PairInt> excludePoints, ImageExt colorImg, GreyscaleImage mask,
        int binFactor) {

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

        double rDivB = allSkyColor.getAvgRed()/allSkyColor.getAvgBlue();

        boolean skyIsRed = (rDivB > 1);

        PerimeterFinder perimeterFinder = new PerimeterFinder();

        int[] dxs = new int[]{-1, -1,  0,  1, 1, 1, 0, -1};
        int[] dys = new int[]{ 0, -1, -1, -1, 0, 1, 1,  1};

        int nAdded = 0;
        int nIter = 0;
        int nMaxIter = binFactor;

        int[] skyRowMinMax = new int[2];
        Set<PairInt> embeddedPoints = new HashSet<PairInt>();
        Map<Integer, List<PairInt>> skyRowColRange = perimeterFinder.find(points, skyRowMinMax, 
            colorImg.getWidth(), embeddedPoints);
        Set<PairInt> borderPixels = null;
        
        while ((nIter == 0) || ((nIter < nMaxIter) && (nAdded > 0))) {
            
            if (borderPixels != null) {
                skyRowColRange = perimeterFinder.find(
                    points, skyRowMinMax, imageMaxColumn, 
                    embeddedPoints);
            }
            borderPixels = perimeterFinder.getBorderPixels(
                skyRowColRange, skyRowMinMax, imageMaxColumn, imageMaxRow);
            
            nAdded = 0;          

            for (PairInt uPoint : borderPixels) {

                int uX = uPoint.getX() + xOffset;
                int uY = uPoint.getY() + yOffset;
                
                Set<PairInt> neighborsInSky = new HashSet<PairInt>();
                Set<PairInt> neighborsNotInSky = new HashSet<PairInt>();
                                        
                for (int k = 0; k < dxs.length; k++) {

                    int vX = uX + dxs[k];
                    int vY = uY + dys[k];

                    if ((vX < 0) || (vX > (cWidth - 1)) || (vY < 0) || 
                        (vY > (cHeight - 1))) {
                        continue;
                    }

                    PairInt vPoint = new PairInt(vX - xOffset, vY - yOffset);
                    
                    if (uPoint.equals(vPoint) || excludePoints.contains(vPoint)) {
                        continue;
                    }
                    
                    if (points.contains(vPoint)) {
                        neighborsInSky.add(vPoint);
                    } else {
                        neighborsNotInSky.add(vPoint);
                    }
                }
                
                if (!neighborsNotInSky.isEmpty()) {
                    
                    GroupPixelColors gpc = new GroupPixelColors(neighborsInSky,
                        colorImg, xOffset, yOffset);
                    
                    for (PairInt p : neighborsNotInSky) {
                        
                        int x = p.getX();
                        int y = p.getY();
                        
                        int r = colorImg.getR(x, y);
                        int g = colorImg.getG(x, y);
                        int b = colorImg.getB(x, y);
                        
                        float contrast = gpc.calcContrastToOther(r, g, b);
                        double contrastDivStDev = contrast/gpc.getStdDevContrast();
                        float colorDiff = gpc.calcColorDiffToOther(r, g, b);
                        double colorDiffForSkyColor = skyIsRed ? 
                            (r - gpc.getAvgRed()) : (b - gpc.getAvgBlue());
                        double colorDiffForSkyColorDivStDev = skyIsRed ?
                            colorDiffForSkyColor/gpc.getStdDevRed() :
                            colorDiffForSkyColor/gpc.getStdDevBlue();
                        double diffCIEX = colorImg.getCIEX(x, y) - gpc.getAverageCIEX();
                        double diffCIEY = colorImg.getCIEY(x, y) - gpc.getAverageCIEY();
                        double diffCIEXDivStDev = Math.abs(diffCIEX)/gpc.getStdDevCIEX();
                        double diffCIEYDivStDev = Math.abs(diffCIEY)/gpc.getStdDevCIEY();
                           
                        if (
                            (diffCIEX < 0.01) && (diffCIEY < 0.01)
                            && (Math.abs(contrast) < 0.01)  //0.045 for less conservative
                            ) {
                            points.add(new PairInt(x - xOffset, y - yOffset));
                            nAdded++;
                        } else if (
                            (diffCIEX < 0.005) && (diffCIEY < 0.005)
                            && (Math.abs(contrast) < 0.005)
                            ) {
                            points.add(new PairInt(x - xOffset, y - yOffset));
                            nAdded++;
                        }
                    }
                }
            }
            
            System.out.println("nAdded=" + nAdded);
            
            nIter++;
        }
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
    public double[] findSunPoints(ImageExt colorImg, int xOffset, 
        int yOffset, boolean skyIsDarkGrey, Set<PairInt> outputSunPoints) {
                
        SunColors sunColors = new SunColors();
        
        if (skyIsDarkGrey) {
            sunColors.useDarkSkiesLogic();
        }
        
        for (int col = 0; col < colorImg.getWidth(); col++) {
            for (int row = 0; row < colorImg.getHeight(); row++) {
                
                int idx = colorImg.getInternalIndex(col, row);

                if (sunColors.isSunCenterColor(colorImg, idx)) {
                    outputSunPoints.add(new PairInt(col - xOffset, row - yOffset));
                }
            }
        }
        
        log.info("found " + outputSunPoints.size() 
            + " sun points from image of width " + colorImg.getWidth());

        if (outputSunPoints.size() < 6) {
            return null;
        }
        
        //fit ellipse to yellowPoints.  ellipse because of possible occlusion.
        EllipseHelper ellipseHelper = new EllipseHelper();
        double[] params = ellipseHelper.fitEllipseToPoints(outputSunPoints);
        
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
       
        //TODO:  this may need adjustments.  top of sun rising over mountains..
        /*if ((a/b) > 6) {
            return new HashSet<PairInt>();
        }*/
        
        return params;
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
    
    public Set<PairInt> findEmbeddedNonPoints(Set<PairInt> points, 
        int imageMaxColumn) {
        
        PerimeterFinder perimeterFinder = new PerimeterFinder();
        
        int[] skyRowMinMax = new int[2];
        
        Set<PairInt> outputEmbeddedGapPoints = new HashSet<PairInt>();
        
        Map<Integer, List<PairInt>> skyRowColRange = 
            perimeterFinder.find(points, 
            skyRowMinMax, imageMaxColumn, outputEmbeddedGapPoints);
        
        return outputEmbeddedGapPoints;
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
        
        /*
        usually, the largest smooth component of the image is sky and for that
        case, one just needs to retain the largest contiguous list in
        zeroPointLists as the seed of the sky points from which to find
        the remaining sky points.
        
        For a case such as a large smooth foreground field that has colors 
        similar to sky, such as a featureless snow field, one needs to 
        use other information to distinguish which points are sky points.
        --> this can be seen when all lists in zeroPointLists, that is all points
        in zeroPointLists have cie X and Y white.
        */
        // if all are white, reduce to the bluest lists
        boolean didReduceToBluest = reduceToBluestIfAllAreWhite(zeroPointLists, 
            colorImg, thetaImg);
        
        //TODO:  assumes that largest smooth component of image is sky.
        // if sky is small and a foreground object is large and featureless
        // and not found as dark, this will fail. 
        // will adjust for that one day, possibly with color validation
        //if (!didReduceToBluest) {
            reduceToLargest(zeroPointLists);
        //}
        
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
    
    private ANDedClauses[] filterForSkyColor(ANDedClauses[] clauses, 
        boolean skyIsRed) {
        
        int n = 0;
        for (int i = 0; i < clauses.length; i++) {
            
            ANDedClauses c = clauses[i];
            
            if (c.getSKYCONDITIONAL().ordinal() == SKYCONDITIONAL.ALL.ordinal()) {
                n++;
            } else if (skyIsRed && (c.getSKYCONDITIONAL().ordinal() == 
                SKYCONDITIONAL.RED.ordinal())) {
                n++;
            } else if (!skyIsRed && (c.getSKYCONDITIONAL().ordinal() == 
                SKYCONDITIONAL.BLUE.ordinal())) {
                n++;
            }
        }
        
        ANDedClauses[] filtered = new ANDedClauses[n];
        int ii = 0;
        for (int i = 0; i < clauses.length; i++) {
            
            ANDedClauses c = clauses[i];
            
            if (c.getSKYCONDITIONAL().ordinal() == SKYCONDITIONAL.ALL.ordinal()) {
                filtered[ii] = c;
                ii++;
            } else if (skyIsRed && (c.getSKYCONDITIONAL().ordinal() == 
                SKYCONDITIONAL.RED.ordinal())) {
                filtered[ii] = c;
                ii++;
            } else if (!skyIsRed && (c.getSKYCONDITIONAL().ordinal() == 
                SKYCONDITIONAL.BLUE.ordinal())) {
                filtered[ii] = c;
                ii++;
            }
        }
        
        return filtered;
    }

    protected List<PairInt> orderByProximity(Set<PairInt> points) {
        
        /*
        starting from (startX, startY), order the points into a list of 
        sequentially closest points.
        for points with more than one adjacent point, the order should
        be such that horiz adjacent or vert adjacent is preferred over diag,
        and if there is more than one horiz or vert adjacent,
        it shouldn't matter which is chosen.
        */
        
        // O( N * lg_2(N) ) + O(N)
        DoubleLinkedCircularList tmp = new DoubleLinkedCircularList();
        Map<Integer, HeapNode> firstXIndexes = 
            new HashMap<Integer, HeapNode>();
        orderByIncreasingXThenY(points, tmp, firstXIndexes);
                
        List<PairInt> ordered = new ArrayList<PairInt>();
        
        int lastX = -1;
        int lastY = -1;
        
        /*
        Because this is ordered already, the runtime complexity should be 
        close to O(N) when there are no gaps which is usually the case for
        the context this method was created for, else the runtime complexity is 
        close to O(N!) when the entire remaining list has to be compared.
        TODO: replace w/ a data structure that has indexed y ordering too.
        */
        while (tmp.getSentinel().getLeft() != tmp.getSentinel()) {
                            
            //if have adjacent points, stop search in tmp when dx > 1
            boolean atLeastOneIsAdjacent = false;
            
            int minDistSq = Integer.MAX_VALUE;
            HeapNode closest = null;
            
            HeapNode currentNode = tmp.getSentinel().getLeft();
            while (currentNode != tmp.getSentinel()) {
                
                PairInt p = (PairInt)currentNode.getData();

                if (lastX == -1) {
                    closest = currentNode;
                    break;
                } else {
                    int dx = Math.abs(p.getX() - lastX);
                    int dy = Math.abs(p.getY() - lastY);
                    int distSq = (dx * dx) + (dy * dy);
                    if ((dx <= 1) && (dy <= 1)) {
                        atLeastOneIsAdjacent = true;
                        if (distSq < minDistSq) {
                            closest = currentNode;
                            minDistSq = distSq;
                            if (distSq == 1) {
                                // this is horiz or vertical adjacent which are preferred
                                break;
                            }
                        }
                    } else if (!atLeastOneIsAdjacent) {
                        // there may not be an adjacent, so store the next closest
                        if (distSq < minDistSq) {
                            closest = currentNode;
                            minDistSq = distSq;
                        }
                    } else if (atLeastOneIsAdjacent && (dx > 1)) {
                        // we've gone past the nearest, so can break
                        break;
                    }
                }
                
                currentNode = currentNode.getLeft();
            }
            
            //O(1)
            ordered.add((PairInt)closest.getData());
            
            //O(1)
            tmp.remove(closest);
            
            lastX = ((PairInt)closest.getData()).getX();
            
            lastY = ((PairInt)closest.getData()).getY();
        }
        
        return ordered;
    }

    protected void orderByIncreasingXThenY(Set<PairInt> points,
        DoubleLinkedCircularList outputDLCList, 
        Map<Integer, HeapNode> outputFirstXLocator) {

        int[] x = new int[points.size()];
        int[] y = new int[x.length];
        int i = 0;
        for (PairInt p : points) {
            x[i] = p.getX();
            y[i] = p.getY();
            i++;
        }
        //O(N*lg_2(N))
        MultiArrayMergeSort.sortBy1stArgThen2nd(x, y);
        
        //O(N)
        for (i = 0; i < x.length; ++i) {
            
            PairInt p = new PairInt(x[i], y[i]);
            
            HeapNode node = new HeapNode(i);
            node.setData(p);
            outputDLCList.insert(node);
            
            Integer xCoord = Integer.valueOf(p.getX());
            if (!outputFirstXLocator.containsKey(xCoord)) {
                outputFirstXLocator.put(xCoord, node);
            }
        }
    }

    private void addRainbowToPoints(Set<PairInt> skyPoints, 
        Set<PairInt> rainbowHullPoints, Hull rainbowHull,
        int lastImgCol, int lastImgRow) {
        
        //TODO: Note, it may be necessary to build a hull from a spine of
        // more than 10 points for complex images w/ rainbow intersecting
        // multiple times with structured foreground and sky
        
        /* n=21
         0,20    1    2    3    4    5    6    7    8    9
                                                  
         19   18   17   16   15   14   13   12   11   10
        
        Points 0 and 19 are opposite points on the hull, so are 1 and 18, etc.
        So can do quick check that for each segment such as 0->1->18->19->0
        that there are skyPoints surrounding them.
        If not, and if the coverage is partial, will reduce the segment polygon
        size and only add points to skypoints within the segment polygon.
        */
        
        int nHull = rainbowHull.xHull.length;
        int nHalf = nHull >> 1;
        for (int c = 0; c < nHalf; c++) {
            
            int count0 = c;
            int count1 = nHull - 2 - c;
            
            double dy0 = rainbowHull.yHull[count0 + 1] - rainbowHull.yHull[count0];
            double dx0 = rainbowHull.xHull[count0 + 1] - rainbowHull.xHull[count0];
            int dist0 = (int)Math.sqrt(dx0*dx0 + dy0*dy0);

            double dy1 = rainbowHull.yHull[count1 - 1] - rainbowHull.yHull[count1];
            double dx1 = rainbowHull.xHull[count1 - 1] - rainbowHull.xHull[count1];
            int dist1 = (int)Math.sqrt(dx1*dx1 + dy1*dy1);

            boolean removeSection = false;
            
            int dist = (dist0 < dist1) ? dist0 : dist1;
            for (int i = 0; i < dist; i++) {
                int x0 = (int)(rainbowHull.xHull[count0] + i*dx0);
                int y0 = (int)(rainbowHull.yHull[count0] + i*dy0);

                int x1 = (int)(rainbowHull.xHull[count1] + i*dx1);
                int y1 = (int)(rainbowHull.yHull[count1] + i*dy1);

                int n0Sky = 0;
                int n0SkyPossible = 0;
                int n1Sky = 0;
                int n1SkyPossible = 0;
                for (int type = 0; type < 2; type++) {
                    int x = (type == 0) ? x0 : x1;
                    int y = (type == 0) ? y0 : y1;
                    int n = 0;
                    int nPossible = 0;
                    for (int col = (x - 1); col <= (x + 1); col++) {
                        if ((col < 0) || (col > lastImgCol)) {
                            continue;
                        }
                        for (int row = (y - 1); row <= (y + 1); row++) {
                            if ((row < 0) || (row > lastImgRow)) {
                                continue;
                            }
                            PairInt p = new PairInt(col, row);
                            if (!rainbowHullPoints.contains(p)) {
                                nPossible++;
                                if (skyPoints.contains(p)) {
                                    n++;
                                }
                            }
                        }
                    }
                    if (type == 0) {
                        n0Sky = n;
                        n0SkyPossible = nPossible;
                    } else {
                        n1Sky = n;
                        n1SkyPossible = nPossible;
                    }
                }
                // evaluate the n's and shorten rainbowHull values for count0 and count1 if needed
                float n0Div = (float)n0Sky/(float)n0SkyPossible;
                float n1Div = (float)n1Sky/(float)n1SkyPossible;
                if ((n0Div < 0.5) || (n1Div < 0.5)) {
                    removeSection = true;
                    break;
                }
            } // end for i
            if (removeSection) {
                
                // remove points from rainbowHull
                float[] xh = new float[]{
                    rainbowHull.xHull[count0], rainbowHull.xHull[count0 + 1],
                    rainbowHull.xHull[count1 - 1], rainbowHull.xHull[count1],
                    rainbowHull.xHull[count0]};
                float[] yh = new float[]{
                    rainbowHull.yHull[count0], rainbowHull.yHull[count0 + 1],
                    rainbowHull.yHull[count1 - 1], rainbowHull.yHull[count1],
                    rainbowHull.yHull[count0]};
                
                int minXHull = (int)MiscMath.findMin(xh);
                int maxXHull = (int)Math.ceil(MiscMath.findMax(xh));
                int minYHull = (int)MiscMath.findMin(yh);
                int maxYHull = (int)Math.ceil(MiscMath.findMax(yh));
                PointInPolygon p = new PointInPolygon();
                for (int col = minXHull; col <= maxXHull; col++) {
                    for (int row = minYHull; row <= maxYHull; row++) {
              
                        boolean in = p.isInSimpleCurve(col, row, xh, yh, 
                            xh.length);
                        if (in) {
                            rainbowHullPoints.remove(new PairInt(col, row));
                        }
                    }
                }
            }
        }
        
        skyPoints.addAll(rainbowHullPoints);
        
    }

    private void addEmbeddedIfSimilarToSky(Set<PairInt> points, 
        ImageExt originalColorImage, GreyscaleImage mask, 
        RemovedSets removedSets) {
        
        Set<PairInt> embeddedPoints = new HashSet<PairInt>();
        getEmbeddedPoints(points, originalColorImage.getWidth(), 
            originalColorImage.getHeight(),
            embeddedPoints);
        
        if (embeddedPoints.isEmpty()) {
            return;
        }
        
        HistogramHolder[] brightnessHistogram = new HistogramHolder[1];

        // brightest sky is in bin [2], and dimmest is in [0]
        GroupPixelColors[] skyPartitionedByBrightness = 
            partitionInto3ByBrightness(points, originalColorImage, 
            mask.getXRelativeOffset(), mask.getYRelativeOffset(), 
            brightnessHistogram);

        // no need to update rowColRanges as this was not an external "grow"
        //add embedded pixels if they're near existing sky colors
        addIfSimilarToSky(embeddedPoints, points, 
            removedSets.getHighContrastRemoved(),
            originalColorImage, mask,
            brightnessHistogram, skyPartitionedByBrightness);        

    }

    private void createDebugPlots(List<PairIntArray> zeroPointLists, 
        ImageExt colorImg, GreyscaleImage thetaImg) {
        
        // scatter diagrams of cie X vs cie Y,  a plot of their (x,y) locations
        // for the image dimensions, and scatter diagram of red/total vs blue/total.
        
        PolygonAndPointPlotter plotter = null;
        try {
            plotter = new PolygonAndPointPlotter();
        } catch (IOException e) {
        }
        
        int xOffset = thetaImg.getXRelativeOffset();
        int yOffset = thetaImg.getYRelativeOffset();
        int count = 0;
        for (PairIntArray pai : zeroPointLists) {
            int n = pai.getN();
            float[] xP = new float[n];
            float[] yP = new float[n];
            float[] rDivTotal = new float[n];
            float[] bDivTotal = new float[n];
            CIEChromaticity cieC = new CIEChromaticity();
            float[] cieX = new float[n];
            float[] cieY = new float[n];
            for (int i = 0; i < n; i++) {
                int x = pai.getX(i) + xOffset;
                int y = pai.getY(i) + yOffset;
                
                xP[i] = x;
                yP[i] = y;
                
                int r = colorImg.getR(x, y);
                int g = colorImg.getG(x, y);
                int b = colorImg.getB(x, y);
                int rgbTotal = r + g + b;
                rDivTotal[i] = (float)r/(float)rgbTotal;
                bDivTotal[i] = (float)b/(float)rgbTotal;
                
                float[] cie = cieC.rgbToXYChromaticity(r, g, b);
                cieX[i] = cie[0];
                cieY[i] = cie[1];                
            }
        
            // scatter diagrams of cie X vs cie Y,  a plot of their (x,y) locations
            // for the image dimensions, and scatter diagram of red/total vs blue/total.

            plotter.addPlot(0.f, 1.f, 0.f, 1.f, cieX, cieY, 
                null, null, "cie X vs cie Y (" + count + ")");

            plotter.addPlot(0.f, 0.4f, 0.f, 0.4f, rDivTotal, bDivTotal, 
                null, null, "R/RGB vs B/RGB (" + count + ")");

            plotter.addPlot(0.f, colorImg.getWidth(), 0, colorImg.getHeight(), 
                xP, yP, null, null, "x vs y (" + count + ")");

            count++;
        }
        
        try {
            String fileName = plotter.writeFile(Integer.valueOf(700 + count));

            System.out.println("fileName=" + fileName);
        
        } catch (IOException e) {
            Logger.getLogger(this.getClass().getName()).severe(e.getMessage());
        }    
    }
    
    private void calculateAverageColors(List<PairIntArray> zeroPointLists, 
        ImageExt colorImg, GreyscaleImage thetaImg, 
        float[] outputCIEX, float[] outputCIEY, 
        float[] outputRDIVRGB, float[] outputBDIVRGB) {
        
        int xOffset = thetaImg.getXRelativeOffset();
        int yOffset = thetaImg.getYRelativeOffset();
        int count = 0;
        CIEChromaticity cieC = new CIEChromaticity();
        for (PairIntArray pai : zeroPointLists) {
                        
            double sumCIEX = 0;
            double sumCIEY = 0;
            double sumRDIVGRB = 0;
            double sumBDIVGRB = 0;
            
            int n = pai.getN();
            
            for (int i = 0; i < n; i++) {
                int x = pai.getX(i) + xOffset;
                int y = pai.getY(i) + yOffset;
                                
                int r = colorImg.getR(x, y);
                int g = colorImg.getG(x, y);
                int b = colorImg.getB(x, y);
                
                sumRDIVGRB += ((double)r/(r + g + b));
                sumBDIVGRB += ((double)b/(r + g + b));
                
                float[] cie = cieC.rgbToXYChromaticity(r, g, b);
                if (cie[0] < 0.2 || cie[1] < 0.2) {
                    int z = 1;
                }
                sumCIEX += cie[0];
                sumCIEY += cie[1];
            }
            outputCIEX[count] = (float)(sumCIEX/(float)n);
            outputCIEY[count] = (float)(sumCIEY/(float)n);
            
            outputRDIVRGB[count] = (float)(sumRDIVGRB/(float)n);
            outputBDIVRGB[count] = (float)(sumBDIVGRB/(float)n);
            
            count++;
        }
    }

    private boolean reduceToBluestIfAllAreWhite(List<PairIntArray> zeroPointLists, 
        ImageExt colorImg, GreyscaleImage thetaImg) {
        
        //createDebugPlots(zeroPointLists, colorImg, thetaImg);
        
        float[] cieX = new float[zeroPointLists.size()];
        float[] cieY = new float[zeroPointLists.size()];
        float[] rDIVRGB = new float[zeroPointLists.size()];
        float[] bDIVRGB = new float[zeroPointLists.size()];
                
        calculateAverageColors(zeroPointLists, colorImg, thetaImg, 
            cieX, cieY, rDIVRGB, bDIVRGB);        
        
        //TODO: may need to revise the ranges here
        boolean allAreWhite = true;
        for (int i = 0; i < cieX.length; i++) {
            float x = cieX[i];
            float y = cieY[i];
            float xDivY = x/y;
            if (!((x >= 0.23) && (x <= 0.43) && (Math.abs(xDivY - 1) < 0.2))){
                allAreWhite = false;
                break;
            }
            if ((x > 0.33) && (y < 0.28)) {
                allAreWhite = false;
                break;
            }
            if ((x > 0.25) & (y < 0.2)) {
                allAreWhite = false;
                break;
            }
        }
        if (!allAreWhite) {
System.out.println("cieX=" + Arrays.toString(cieX));
System.out.println("cieY=" + Arrays.toString(cieY));
            return false;
        }
        
        // Looking for differences between a snow field and clouds.
        // The Warren 1982 paper suggests that snow will have a redder slope
        // over the optical bandpass than optically thick clouds.  
        // The optically thick clouds will be grey if not slightly blue.
        
        HistogramHolder bHist = Histogram.createSimpleHistogram(3,
            bDIVRGB, Errors.populateYErrorsBySqrt(bDIVRGB));
        System.out.println("bHist y=" + Arrays.toString(bHist.getYHist()) 
            + " x=" + Arrays.toString(bHist.getXHist()));
        
        /*try {
            bHist.plotHistogram("blueHist", 8001);
        } catch (IOException e) {
        }*/
        
        List<PairIntArray> keep = new ArrayList<PairIntArray>();
        float limit = 0.5f*(bHist.getXHist()[1] + bHist.getXHist()[0]);
        for (int i = 0; i < bDIVRGB.length; i++) {
            if (bDIVRGB[i] >= limit) {
                keep.add(zeroPointLists.get(i));
            }
        }
        
        zeroPointLists.clear();
        zeroPointLists.addAll(keep);
        
        return true;
    }
    
    private class Hull {
        float[] xHull;
        float[] yHull;
    }
    
    /**
     * 
     * @param skyPoints
     * @param rainbowCoeff coefficients for a 2nd order polynomial
     * @param rainbowPoints  rainbow colored points in the image that fit a polynomial.
     * there should be 10 or more points at least.
     * @param originalColorImage
     * @param xOffset
     * @param yOffset
     * @return 
     */
    private Hull createRainbowHull(Set<PairInt> skyPoints, 
        float[] rainbowCoeff, Set<PairInt> rainbowPoints, 
        ImageExt originalColorImage, int xOffset, int yOffset) {
        
        /*
        need to know the furthest closest distanc to the polynomial, that is the
        furthest perpendicular point to the polynomial in order to expand the 
        region around the polynomial to become a hull that encloses all rainbow 
        points.
        
        (1) Could determine for every point, the distance to the polynomial:
        Robust and Efficient Computation of the Closest Point on a Spline Curve" 
        Hongling Wang, Joseph Kearney, and Kendall Atkinson
        http://homepage.cs.uiowa.edu/~kearney/pubs/CurvesAndSufacesClosestPoint.pdf
        implemented by:
        http://www.ogre3d.org/tikiwiki/tiki-index.php?page=Nearest+point+on+a+Spline
        
        -- The method for nearest point on a spline requires creating a
        spline out of the polynomial.  can assume that something like 1000
        points would be necessary, though maybe 100 would be fine.
        -- Then the method requires a good starting guess for 3 points on the
        polynomial that would be near the true perpendicular point.
        One could guess the first 3 splines by making the polynomial roughly
        10 splines separately and doing a distance test to each spline
        for each point (with shortcuts for not needing to check).
        -- Then about 10 iterations or less to find the best answer.
        
        The runtime might look roughly like
        N_poly +  N_points * 10 distance tests + N_points * 10 iterations of the 
        closestToSpline algorithm.
        
        (2) Or could instead, characterize the hull by 10 points and calculate 
        the coordinates of the points projected perpendicular to them above and 
        below at distances that are the hull of the rainbow.
        Then evaluation of the size would be "point in polygon" tests for
        each point.
        Improved estimates of the hull size can use a binomial search pattern
        to expand or shrink the distance estimate used for furthest extension
        of the hull.
        The first estimate of the hull size can use the rough distance point test
        just as above, by comparing each point to the polynomial as 10 segments.
        
        runtime is roughly:
            N_points * 10 distance tests +
            N_iter * N_points *
               "point in polygon" (which is roughly O(N_poly) where N_poly would be ~20)
            = N_points * 10 + N_iter * N_points * 20 for each iteration
        
        expect that N_iter is probably <= 5 because it's close already.
        
        So the 2nd method will be fine for the purposes here, but the first
        method is interesting for use cases which need precision.
        */
        
        int width = originalColorImage.getWidth() - xOffset;
        int height = originalColorImage.getHeight() - yOffset;
        
        float[] xc = new float[10];
        float[] yc = new float[xc.length];
        generatePolynomialPoints(rainbowPoints, rainbowCoeff, xc, yc);
       
        float maxOfPointMinDistances = maxOfPointMinDistances(rainbowPoints,
            xc, yc);
        
        float high = 2 * maxOfPointMinDistances;
        float low = maxOfPointMinDistances / 2;
        int nMatched = 0;
        
        /* 
        n=21
         0,20    1    2    3    4    5    6    7    8    9
                                                  
         19   18   17   16   15   14   13   12   11   10
        */
        
        float[] xPoly = new float[2 * xc.length + 1];
        float[] yPoly = new float[xPoly.length];
        int nMaxIter = 5;
        int nIter = 0;
        
        int eps = (int)(0.01f * rainbowPoints.size());
        
        while ((low < high) && (nIter < nMaxIter)) {
            
            float mid = (high + low)/2.f;
            
            populatePolygon(xc, yc, mid, xPoly, yPoly, rainbowCoeff, width, 
                height);
            
            nMatched = nPointsInPolygon(rainbowPoints, xPoly, yPoly);
            
System.out.println("low=" + low + " high=" + high + " mid=" + mid + " nMatched=" + nMatched);

            if (Math.abs(nMatched - rainbowPoints.size()) < eps) {
                if (low < mid) {
                    // decrease high so next mid is lower
                    high = mid;
                } else {
                    break;
                }
            } else {
                // nMatched < rainbowPoints.size()
                // increase low so next mid is higher
                low = mid;
            }
            
            nIter++;
        }
        
debugPlot(rainbowPoints, xPoly, yPoly, originalColorImage.getWidth(),
originalColorImage.getHeight(), "hull around rainbow");
int z = 1;   
        Hull hull = new Hull();
        hull.xHull = xPoly;
        hull.yHull = yPoly;

        return hull;
    }
    
    protected void generatePolynomialPoints(Set<PairInt> points, 
        float[] polyCoeff, float[] outputX, float[] outputY) {
        
        /*
        determine good endpoints for the polynomial solution.
        it depends upon the orientation of the rainbow polynomial.
        
        should be able to determine the average to the 100 or so median
        values as the middle of the rainbow,
        then would find the smallest residual points from the model
        polynomial that are located furthest from the median location.
        */
        
        // sort points by x then y
        int[] x = new int[points.size()];
        int[] y = new int[x.length];
        int i = 0;
        for (PairInt p : points) {
            x[i] = p.getX();
            y[i] = p.getY();
            i++;
        }
        //O(N*lg_2(N))
        MultiArrayMergeSort.sortBy1stArgThen2nd(x, y);
        
        // average of central 10 or so median
        int mid = points.size() >> 1;
        float medianX;
        float medianY;
        if (points.size() < 10) {
            medianX = x[points.size()/2];
            medianY = y[points.size()/2];
        } else {
            double sumX = 0;
            double sumY = 0;
            int nMed = 5;
            for (i = (mid - nMed); i <= (mid + nMed); i++) {
                sumX += x[i];
                sumY += y[i];
            }
            medianX = (float)sumX/(float)(2*nMed);
            medianY = (float)sumY/(float)(2*nMed);
        }
        
        // find the furthest points that have the smallest residuals on each side
        int minX = -1;
        int yForMinX = -1;
        int maxX = -1;
        int yForMaxX = -1;
        for (int half = 0; half < 2; half++) {
            double minResid = Double.MAX_VALUE;
            int minResidIdx = -1;
            double distFromMedianSq = Double.MIN_VALUE;
            if (half == 0) {
                for (i = 0; i < mid; ++i) {
                    float yV = polyCoeff[0] * (polyCoeff[1] * x[i]) +
                        (polyCoeff[2] * x[i] * x[i]);
                    double resid = Math.abs(y[i] - yV);

                    double distX = x[i] - medianX;
                    double distY = y[i] - medianY;

                    double distSq = (distX * distX) + (distY * distY);

                    if ((resid < minResid) && (distSq >= distFromMedianSq)) {
                        minResid = resid;
                        minResidIdx = i;
                        distFromMedianSq = distSq;
                    }
                }
                minX = x[minResidIdx];
                yForMinX = y[minResidIdx];
            } else {
                for (i = (points.size() - 1); i > mid; --i) {
                    float yV = polyCoeff[0] * (polyCoeff[1] * x[i]) +
                        (polyCoeff[2] * x[i] * x[i]);
                    double resid = Math.abs(y[i] - yV);

                    double distX = x[i] - medianX;
                    double distY = y[i] - medianY;

                    double distSq = (distX * distX) + (distY * distY);

                    if ((resid < minResid) && (distSq >= distFromMedianSq)) {
                        minResid = resid;
                        minResidIdx = i;
                        distFromMedianSq = distSq;
                    }
                }
                maxX = x[minResidIdx];
                yForMaxX = y[minResidIdx];
            }
        } 
        
//[204.83781, 0.40015212, -3.880514E-4]
System.out.println("polyCoeff=" + Arrays.toString(polyCoeff) 
+ " endpoints=(" + minX + "," + yForMinX + ") (" + maxX + "," + yForMaxX + ")");
       
        int n = outputX.length;
        
        // max-min divided by 9 gives 8 more points
        float deltaX = (maxX - minX)/(float)(n - 1);
        
        outputX[0] = minX;
        outputY[0] = yForMinX;
        for (i = 1; i < (n - 1); i++) {
            outputX[i] = outputX[i - 1] + deltaX;
            outputY[i] = polyCoeff[0] + polyCoeff[1] * outputX[i] 
                + polyCoeff[2] * outputX[i] * outputX[i];
        }
        outputX[n - 1] = maxX;
        outputY[n - 1] = yForMaxX;
        
    }

    private float maxOfPointMinDistances(Set<PairInt> rainbowPoints, float[] xc, 
        float[] yc) {
        
        double maxDistSq = Double.MIN_VALUE;
        
        for (PairInt p : rainbowPoints) {
            int x = p.getX();
            int y = p.getY();
            double minDistSq = Double.MAX_VALUE;
            for (int i = 0; i < xc.length; i++) {
                float diffX = xc[i] - x;
                float diffY = yc[i] - y;
                float dist = (diffX * diffX) + (diffY * diffY);
                if (dist < minDistSq) {
                    minDistSq = dist;
                }
            }
            if (minDistSq > maxDistSq) {
                maxDistSq = minDistSq;
            }
        }
        
        return (float)Math.sqrt(maxDistSq);
    }

    /**
    populate outputXPoly and outputYPoly with points perpendicular to x and y
    * at distances dist.  note that the lengths of outputXPoly and outputYPoly
    * should be 2*x.length+1
    */
    protected void populatePolygon(float[] x, float[] y, float dist, 
        float[] outputXPoly, float[] outputYPoly, float[] polynomialCoeff,
        int imgWidth, int imgHeight) {
        
        if (x == null || y == null) {
            throw new IllegalArgumentException("neither x nor y can be null");
        }
        if (x.length != y.length) {
            throw new IllegalArgumentException("x and y must be the same length");
        }
        if (outputXPoly == null || outputYPoly == null) {
            throw new IllegalArgumentException(
                "neither outputXPoly nor outputYPoly can be null");
        }
        if (outputXPoly.length != outputYPoly.length) {
            throw new IllegalArgumentException(
                "outputXPoly and outputYPoly must be the same length");
        }
        if (polynomialCoeff == null) {
            throw new IllegalArgumentException(
                "polynomialCoeff cannot be null");
        }
        if (polynomialCoeff.length != 3) {
            throw new IllegalArgumentException(
                "polynomialCoeff.length has to be 3");
        }
        if (outputXPoly.length != (2*x.length + 1)) {
            throw new IllegalArgumentException("outputXPoly.length must be " +
                " (2 * x.length) + 1");
        }
        
        /*
        y = c0*1 + c1*x[i] + c2*x[i]*x[i]
        
        dy/dx = c1 + 2*c2*x[i]
        */
        
        int n = outputXPoly.length;
        
        /*
        want them in order so using count0 and count1
        n=21
        
        0,20    1    2    3    4    5    6    7    8    9
                                                  
         19   18   17   16   15   14   13   12   11   10
        */
        int count0 = 0;
        int count1 = n - 2;
        
        for (int i = 0; i < x.length; i++) {
            
            double dydx = polynomialCoeff[1] + (2. * polynomialCoeff[2] * x[i]);
            
            if (dydx == 0) {
                // same x, y's are +- dist
                outputXPoly[count0] = x[i];
                outputYPoly[count0] = y[i] + dist;
                count0++;
                outputXPoly[count1] = x[i];
                outputYPoly[count1] = y[i] - dist;
                count1--;
                continue;
            }
            
            double tangentSlope = -1./dydx;
            
            double theta = Math.atan(tangentSlope);
            
            double dy = dist * Math.sin(theta);
            
            double dx = dist * Math.cos(theta);
            
            if ((count0 == 0) && (x[i] == 0)) {
                dy = dist;
                dx = 0;
            } else if (
                ((count0 == ((x.length/2) - 1) || (count0 == (x.length/2))))
                && (x[i] == (imgWidth - 1))) {
                dy = dist;
                dx = 0;
            }
            
            //System.out.println("i=" + i + " theta=" + theta
            //    + " x[i]=" + x[i] + " dx=" + dx + " dy=" + dy);
          
            float xHigh = (float)(x[i] + dx);
            float yHigh = (float)(y[i] + dy);
            
            float xLow = (float)(x[i] - dx);
            float yLow = (float)(y[i] - dy);
              
            if (theta < 0) {
                float tx = xHigh;
                float ty = yHigh;
                xHigh = xLow;
                yHigh = yLow;
                xLow = tx;
                yLow = ty;
            }
            
            outputXPoly[count0] = xHigh;
            outputYPoly[count0] = yHigh;
            count0++;
            outputXPoly[count1] = xLow;
            outputYPoly[count1] = yLow;
            count1--;
        }
        
        outputXPoly[outputXPoly.length - 1] = outputXPoly[0];
        outputYPoly[outputXPoly.length - 1] = outputYPoly[0];
    }

    protected int nPointsInPolygon(Set<PairInt> rainbowPoints, 
        float[] xPoly, float[] yPoly) {
        
        PointInPolygon pIn = new PointInPolygon();
        
        int nInside = 0;
        
        for (PairInt p : rainbowPoints) {
            float x = p.getX();
            float y = p.getY();
            
            boolean ans = pIn.isInSimpleCurve(x, y, xPoly, yPoly, xPoly.length);
            
            if (ans) {
                nInside++;
            }
        }
        
        return nInside;
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
    
    private String debugPlot(Set<PairInt> points, float[] xPoly, float[] yPoly,
        int plotXMax, int plotYMax, String plotLabel) {
        
        float xMin = Float.MAX_VALUE;
        float xMax = Float.MIN_VALUE;
        float[] xP = new float[points.size()];
        float[] yP = new float[xP.length];
        int i = 0;
        for (PairInt p : points) {
            float x = p.getX();
            if (x < xMin) {
                xMin = x;
            }
            if (x > xMax) {
                xMax = x;
            }
            xP[i] = x;
            yP[i] = p.getY();
            i++;
        }
        
        try {
            PolygonAndPointPlotter plotter = new PolygonAndPointPlotter(0, 
                plotXMax, 0, plotYMax);
            plotter.addPlot(xP, yP, xPoly, yPoly, plotLabel);
            
            String fileName = plotter.writeFile(Integer.valueOf(9182));
            
            return fileName;
            
        } catch (IOException e) {
            Logger.getLogger(this.getClass().getName()).severe(e.getMessage());
        }
        
        return null;
    }
}
