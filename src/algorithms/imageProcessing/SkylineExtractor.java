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
import algorithms.misc.AverageUtil;
import algorithms.misc.Histogram;
import algorithms.misc.HistogramHolder;
import algorithms.misc.MiscMath;
import algorithms.util.ArrayPair;
import algorithms.util.Errors;
import algorithms.util.PairFloat;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.PolygonAndPointPlotter;
import algorithms.util.ScatterPointPlotterPNG;
import java.awt.Color;
import java.io.IOException;
import java.security.NoSuchAlgorithmException;
import java.util.*;
import java.util.logging.Logger;

/**
 * class to find the sky points in an image.  TODO: there needs
 * to be another class that is usable to combine information (sky masks) from
 * multiple images that are registered (aligned) in order to better define 
 * the skyline when there are optically thick clouds for instance.
 * 
 * @author nichole
 */
public class SkylineExtractor {

    private Logger log = Logger.getLogger(this.getClass().getName());
    
    private boolean useAlternateSkySeeds = false;
    
    private SunFinder sunFinder = null;
    private RainbowFinder rainbowFinder = null;
    
    private List<PairIntArray> skylineEdges = null;
    
    public static String debugName = "";
    public static void setDebugName(String name) {
        debugName = name;
    }
    
    /**
     * NOT YET USABLE.  This method is a placeholder for a flag that can
     * be used to better select the sky location, perhaps by user
     * interaction with UI.
     */
    public void doUseAlternateSkySelection() {
        
        useAlternateSkySeeds = true;
        
        throw new UnsupportedOperationException("Not supported yet.");
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
            
            // TODO: NOTE, could just find the border pixels of the mask's sky 
            // points which are done already within createBestSkyMask
            // instead of running the mask through the CannyEdgeFilter

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
                
        //TODO: a placeholder for an alternate method to determine seed sky
        // points has been made but not implemented and that's presumably
        // needed for complex cases, such as smooth foreground which
        // has sky colors and a sky that has alot of cloud structure.
        // This method should allow user to see first guess at sky if in
        // interactive mode then allow them to select.
        Set<PairInt> points = extractSkyStarterPoints(theta, gradientXY, 
            originalColorImage, edgeSettings, outputSkyCentroid,
            removedSets);
       
        if (points.isEmpty()) {
            
            GreyscaleImage mask = theta.createWithDimensions();
               
            // return an image of all 1's
            mask.fill(1);
            
            return mask;
        }
        
        int xOffset = theta.getXRelativeOffset();
        int yOffset = theta.getYRelativeOffset();
      
        GroupPixelColors allSkyColor = new GroupPixelColors(points,
            originalColorImage, xOffset, yOffset);
        
        boolean skyIsDarkGrey = skyIsDarkGrey(allSkyColor);
        
        log.info("skyIsDarkGrey=" + skyIsDarkGrey + " " + debugName);
        
        sunFinder = new SunFinder();
        sunFinder.findSunPhotosphere(originalColorImage, xOffset, yOffset, 
            skyIsDarkGrey);

        rainbowFinder = new RainbowFinder();
        
        // should not see sun and rainbow in same image
        if (sunFinder.getSunPoints().isEmpty()) {
                        
            rainbowFinder.findRainbowInImage(
                points, 
                removedSets.getReflectedSunRemoved(), originalColorImage, 
                xOffset, yOffset, theta.getWidth(), theta.getHeight(),
                skyIsDarkGrey, removedSets);
        }

        int nSkyPointsBeforeFindClouds = points.size();
        
        findClouds(points, rainbowFinder.getPointsToExcludeInHull(), 
            originalColorImage, theta);  

        if (!rainbowFinder.getPointsToExcludeInHull().isEmpty()) {
            // addRainbow to Hull, but only if there are sky points adjacent to hull
            rainbowFinder.addRainbowToSkyPoints(points,
                theta.getWidth() - 1, theta.getHeight() - 1);
        }
             
        addEmbeddedIfSimilarToSky(points, originalColorImage, 
            xOffset, yOffset, removedSets);

        Set<PairInt> exclude = new HashSet<PairInt>();
        exclude.addAll(removedSets.getHighContrastRemoved());
        exclude.addAll(sunFinder.getSunPoints());

        // look for patches of sky that are not yet found and are on the image
        // boundaries
        int nAdded = addImageBoundaryEmbeddedSkyIfSimilar(points, exclude, 
            originalColorImage, xOffset, yOffset, removedSets);
        
        log.info("number of sunPoints=" + sunFinder.getSunPoints().size() + " "
        + "reflectedSunRemoved.size()=" + removedSets.getReflectedSunRemoved().size());

        boolean addEmbedded = false;
        
        //if (sunPoints.isEmpty()) {
        
            growForLowContrastLimits(points, exclude, originalColorImage,
                xOffset, yOffset, 
                determineBinFactorForSkyMask(theta.getNPixels()));

            addEmbedded = true;

        //}

        if (!sunFinder.getSunPoints().isEmpty()) {
            
            sunFinder.correctSkylineForSun(points, originalColorImage, xOffset, 
                yOffset, gradientXY);
            
            addEmbedded = true;
        }
        
        if (addEmbedded) {
            addEmbeddedIfSimilarToSky(points, originalColorImage, 
                xOffset, yOffset, removedSets);
        }
                
        points.addAll(sunFinder.getSunPoints());
        
        this.skylineEdges = extractAndSmoothSkylinePoints(points, gradientXY);
        
debugPlot(points, originalColorImage, theta.getXRelativeOffset(), theta.getYRelativeOffset(), 
"final");

        GreyscaleImage mask = gradientXY.createWithDimensions();
        mask.fill(1);
        for (PairInt p : points) {
            int x = p.getX();
            int y = p.getY(); 
            mask.setValue(x, y, 0);
        }
        
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();        
        double[] xycen = curveHelper.calculateXYCentroids(points);
        outputSkyCentroid.add((int)Math.round(xycen[0]), (int)Math.round(xycen[1]));

        ImageProcessor imageProcessor = new ImageProcessor();
        imageProcessor.removeSpurs(mask);
   
        return mask;
    }
    
    /**
     * extract starter points for sky from the theta image and further filtering
     * color and contrast and the gradientXY image.
     * @param theta
     * @param gradientXY
     * @param originalColorImage
     * @param edgeSettings
     * @param outputSkyCentroid
     * @param removedSets
     * @return 
     */
    public Set<PairInt> extractSkyStarterPoints(final GreyscaleImage theta,
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
                
                int idx = originalColorImage.getInternalIndex(ox, oy);
        
                int r = originalColorImage.getR(idx);
                int g = originalColorImage.getG(idx);
                int b = originalColorImage.getB(idx);
                
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
                
                int idx = originalColorImage.getInternalIndex(ox, oy);
        
                int r = originalColorImage.getR(idx);
                int g = originalColorImage.getG(idx);
                int b = originalColorImage.getB(idx);
                
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
            
            int idx = originalColorImage.getInternalIndex(x, y);
        
            int r = originalColorImage.getR(idx);
            int g = originalColorImage.getG(idx);
            int b = originalColorImage.getB(idx);
                
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
                
                float cieX = originalColorImage.getCIEX(vIdx);
                float cieY = originalColorImage.getCIEY(vIdx);
                double diffCIEX = Math.abs(cieX - localSky.getAverageCIEX());
                double diffCIEY = Math.abs(cieY - localSky.getAverageCIEY());
                
                float saturation = originalColorImage.getSaturation(vIdx);
                
                boolean isBrown = (Math.abs(rPercentV - 0.5) < 0.4)
                    && (Math.abs(gPercentV - 0.32) < 0.1)
                    && (Math.abs(bPercentV - 0.17) < 0.1);

                if (isBrown) {
                    
                    // trying to skip over foreground such as land or sunset + water
                    
                    if (((colorDiffV/skyStDevColorDiff) > 15) && (saturation < 0.5)) {
                        
                        if ((saturation <= 0.4) && ((colorDiffV/skyStDevColorDiff) > 50) 
                            ) {

                            continue;
                        }
                         if ((saturation <= 0.4) && 
                            ((Math.abs(contrastV)/Math.abs(skyStDevContrast)) > 10.)
                            ) {
                                                    
                            continue;
                        }
                    }
                }
                
                if ( // no contrast or color change, should be sky
                    (Math.abs(contrastV) < 0.01)
                    && (colorDiffV < 10)
                    && (diffCIEX < 0.009) && (diffCIEY < 0.009)) {

                } else if (
                    (skyStDevContrast != 0.)
                    && ((Math.abs(contrastV)/skyStDevContrast) > 10.)
                    ) {

                    continue;
                } else if (
                    (skyStDevContrast != 0.)
                    && ((Math.abs(contrastV) > 0.1) 
                    && ((Math.abs(contrastV)/skyStDevContrast) > (1.5 + (Math.abs(contrastV)-0.5)*(-2.0))))
                    && ((Math.abs(colorDiffV)/skyStDevColorDiff) > 2.5)
                    ) {

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

                    continue;
                } else if (skyIsRed) {
                                        
                    if (
                        // contrast is defined by luma, so might be weak near
                        // skyline near sun for example
                        (skyStDevContrast != 0.)
                        && ((Math.abs(colorDiffV)/skyStDevColorDiff) > 15.*diffCIEX)
                        && (diffCIEX > 0.03) 
                        && ((diffCIEX/localSky.getStdDevCIEX()) > 15.*diffCIEX)
                        ) {
                        continue;
                    } else if (
                        // contrast is defined by luma, so might be weak near
                        // skyline near sun for example
                        (skyStDevContrast != 0.)
                        && ((Math.abs(colorDiffV)/skyStDevColorDiff) > 15.*diffCIEY)
                        && (diffCIEY > 0.03) 
                        && ((diffCIEY/localSky.getStdDevCIEY()) > 15.*diffCIEY)
                        ) {
                        continue;
                    } else if (skyStDevContrast == 0.) {
                        if (contrastV >= 0.) {
                            doNotAddToStack = true;
                        }
                        
                    } else {
                        //TODO:  if there are sun points, need a zone of
                        // avoidance to not erode the foreground 
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
                        continue;
                    } else if (
                        (skyStDevContrast > 0.0)
                        && (Math.abs(contrastV) > 0.05)
                        && ((Math.abs(colorDiffV)/skyStDevColorDiff) > 1.5)
                        && ((diffCIEX/localSky.getStdDevCIEX()) > 0.9) // ?
                        && ((diffCIEY/localSky.getStdDevCIEY()) > 0.9) // ?                       
                        && (skyStDevColorDiff > 0.)
                    ){
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
                        continue;
                    } else {
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

                if (isBrown) {
                    // trying to skip over foreground such as land or sunset + water                    
                    if ((colorDiffV > 15*skyStDevColorDiff) && (saturation < 0.5)) {
                        if ((saturation <= 0.4) && (colorDiffV > 50*skyStDevColorDiff) 
                            ) {

                            continue;
                        }
                         if ((saturation <= 0.4) && 
                            (Math.abs(contrastV) > 10.*Math.abs(skyStDevContrast))
                            ) {

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
                            cieC.getGreenishYellowThroughOrangePolynomial();
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

                //} else if (isSolarYellow) {
                    // this is a sky point
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
        
        /*try {
            h[0].plotHistogram("sky brightness", 235);
        } catch(Exception e) {
    
        }*/    

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

        int idx = origColorImg.getInternalIndex(col, row);
        
        int r = origColorImg.getR(idx);
        int g = origColorImg.getG(idx);
        int b = origColorImg.getB(idx);
                    
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
        Set<PairInt> excludePoints, ImageExt colorImg, int xOffset, int yOffset,
        int binFactor) {

        // it tries to avoid adding snowy mountain tops to hazy sky pixels,
        // for example.

        int cWidth = colorImg.getWidth();
        int cHeight = colorImg.getHeight();
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
                        
                        int idx = colorImg.getInternalIndex(x, y);
        
                        int r = colorImg.getR(idx);
                        int g = colorImg.getG(idx);
                        int b = colorImg.getB(idx);
                        
                        float contrast = gpc.calcContrastToOther(r, g, b);
                        //double contrastDivStDev = contrast/gpc.getStdDevContrast();
                        //float colorDiff = gpc.calcColorDiffToOther(r, g, b);
                        double colorDiffForSkyColor = skyIsRed ? 
                            (r - gpc.getAvgRed()) : (b - gpc.getAvgBlue());
                        //double colorDiffForSkyColorDivStDev = skyIsRed ?
                        //    colorDiffForSkyColor/gpc.getStdDevRed() :
                        //    colorDiffForSkyColor/gpc.getStdDevBlue();
                        double diffCIEX = colorImg.getCIEX(x, y) - gpc.getAverageCIEX();
                        double diffCIEY = colorImg.getCIEY(x, y) - gpc.getAverageCIEY();
                        //double diffCIEXDivStDev = Math.abs(diffCIEX)/gpc.getStdDevCIEX();
                        //double diffCIEYDivStDev = Math.abs(diffCIEY)/gpc.getStdDevCIEY();
                           
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
            
            int idx = colorImg.getInternalIndex(x, y);

            int rr = colorImg.getR(idx);
            int gg = colorImg.getG(idx);
            int bb = colorImg.getB(idx);
                        
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
                
                int idx = colorImg.getInternalIndex(col, row);
        
                rSum += colorImg.getR(idx);
                gSum += colorImg.getG(idx);
                bSum += colorImg.getB(idx);
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
                
                int idx = originalColorImage.getInternalIndex(ox, oy);
                
                int r = originalColorImage.getR(idx);
                int g = originalColorImage.getG(idx);
                int b = originalColorImage.getB(idx);
                
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
        
        //TODO. This method should one day allow user to see first guess at sky 
        // if in interactive mode then allow them to select sky location.
                
        if (useAlternateSkySeeds) {
            
            zeroPointLists = selectFromSkySeeds(zeroPointLists, colorImg, 
                thetaImg);
            
        } else {
            
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
        }
        
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
        ImageExt originalColorImage, int xOffset, int yOffset,
        HistogramHolder[] brightnessHistogram, 
        GroupPixelColors[] skyPartitionedByBrightness) {
        
        Set<PairInt> excludePoints = highContrastRemoved;
        
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
                skyPartitionedByBrightness,
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

    protected void addEmbeddedIfSimilarToSky(Set<PairInt> points, 
        ImageExt originalColorImage, int xOffset, int yOffset, 
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
            xOffset, yOffset, brightnessHistogram);

        // no need to update rowColRanges as this was not an external "grow"
        //add embedded pixels if they're near existing sky colors
        addIfSimilarToSky(embeddedPoints, points, 
            removedSets.getHighContrastRemoved(),
            originalColorImage, xOffset, yOffset,
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
                
                int idx = colorImg.getInternalIndex(x, y);
                
                int r = colorImg.getR(idx);
                int g = colorImg.getG(idx);
                int b = colorImg.getB(idx);
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

            log.info("fileName=" + fileName);
        
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
                                
                int idx = colorImg.getInternalIndex(x, y);
                
                int r = colorImg.getR(idx);
                int g = colorImg.getG(idx);
                int b = colorImg.getB(idx);
                
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
            return false;
        }
        
        // Looking for differences between a snow field and clouds.
        // The Warren 1982 paper suggests that snow will have a redder slope
        // over the optical bandpass than optically thick clouds.  
        // The optically thick clouds will be grey if not slightly blue.
        
        HistogramHolder bHist = Histogram.createSimpleHistogram(3,
            bDIVRGB, Errors.populateYErrorsBySqrt(bDIVRGB));
        log.info("bHist y=" + Arrays.toString(bHist.getYHist()) 
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

    protected int addImageBoundaryEmbeddedSkyIfSimilar(Set<PairInt> skyPoints, 
        Set<PairInt> exclude, ImageExt originalColorImage, int xOffset,
        int yOffset, RemovedSets removedSets) {
        
        int width = originalColorImage.getWidth();
        int height = originalColorImage.getHeight();
        
        PerimeterFinder perimeterFinder = new PerimeterFinder();

        int imageMaxColumn = width - 1;
        int imageMaxRow = height - 1;
       
        int[] skyRowMinMax = new int[2];
        
        Set<PairInt> embeddedPoints = new HashSet<PairInt>();
        Map<Integer, List<PairInt>> skyRowColRanges = perimeterFinder.find(
            skyPoints, skyRowMinMax, imageMaxColumn, embeddedPoints);
        
        int nEmbeddedBefore = embeddedPoints.size();
        
        fillInImageBoundaryGapsIfBoundedBySky(skyRowColRanges, skyRowMinMax,
            width, height);
        
        embeddedPoints = perimeterFinder.findEmbeddedGivenRowData(skyRowMinMax, 
            width - 1, skyRowColRanges);
        
        if (embeddedPoints.isEmpty()) {
            return 0;
        }
        
        //TODO: revisit this with more test images
        float maxPossibleFracNewDivSky =
            ((float)(skyPoints.size() + embeddedPoints.size())/((float)(width * height)));
        if (maxPossibleFracNewDivSky > 0.75) {
            return 0;
        }
        
        float fracNewDivSky = (float)embeddedPoints.size()/(float)skyPoints.size();
        //AZ        0.0070258044  <== small amount helpful
        //rainbow   0.048877604
        //patagonia 0.45363995    <== HELPFUL
        //rainier   0.080801226   <== adding covers the entire image.  DO NOT USE FOR THIS
        log.info("fracNewDivSky=" + fracNewDivSky 
            + " nEmbeddedBefore=" + nEmbeddedBefore 
            + " embedded=" + embeddedPoints.size() 
            + " nSky/nTotal=" + ((float)skyPoints.size()/((float)(width * height)))
            + " (nSky + nEmbeded)/nTotal=" 
            + ((float)(skyPoints.size() + embeddedPoints.size())/((float)(width * height)))
            + " " + debugName);
        
        int nSkyBefore = skyPoints.size();
        
        HistogramHolder[] brightnessHistogram = new HistogramHolder[1];

        // brightest sky is in bin [2], and dimmest is in [0]
        GroupPixelColors[] skyPartitionedByBrightness = 
            partitionInto3ByBrightness(skyPoints, originalColorImage, 
            xOffset, yOffset, brightnessHistogram);

        // no need to update rowColRanges as this was not an external "grow"
        //add embedded pixels if they're near existing sky colors
        addIfSimilarToSky(embeddedPoints, skyPoints, 
            removedSets.getHighContrastRemoved(),
            originalColorImage, xOffset, yOffset,
            brightnessHistogram, skyPartitionedByBrightness);
        
        return (skyPoints.size() - nSkyBefore);
    }

    private void fillInImageBoundaryGapsIfBoundedBySky(
        Map<Integer, List<PairInt>> skyRowColRanges, int[] skyRowMinMax, 
        int width, int height) {
        
        int imageMaxColumn = width - 1;
        int imageMaxRow = height - 1;
        
        // follow skyRowColRanges around image boundary to see if there gaps
        // on the bounds of the image.
        // Use a "fudge" to gather those image boundary touching pixels that
        //    are not in skyPoints and find the embedded among them then
        //    those that look like sky points:
        
        // ---- fill in any gaps in columns in top row if top row is 0 -----
        if (skyRowMinMax[0] == 0) {
            
            List<PairInt> rowColRanges = skyRowColRanges.get(Integer.valueOf(0));
            
            if (rowColRanges.size() > 1) {
                rowColRanges.get(0).setY(rowColRanges.get(
                    rowColRanges.size() - 1).getY());
                for (int i = rowColRanges.size() - 1; i > 0; i--) {
                    rowColRanges.remove(i);
                }
            }
            
            if (!rowColRanges.isEmpty()) {
                // fill in corners if the gap along that image boundary column
                // is bounded at a higher row.
                boolean fillInLeftCorner = false;
                if (rowColRanges.get(0).getX() > 0) {
                    // ascend rows to see if any have first colRange starting at 0
                    for (int row = skyRowMinMax[0] + 1; row <= skyRowMinMax[1]; row++) {
                        List<PairInt> rcr = skyRowColRanges.get(Integer.valueOf(row));
                        
                        if (!rcr.isEmpty() && (rcr.get(0).getX() == 0)) {
                            fillInLeftCorner = true;
                            break;
                        }
                    }
                }
                if (fillInLeftCorner) {
                    rowColRanges.get(0).setX(0);
                }
            
                boolean fillInRightCorner = false;
                if (rowColRanges.get(rowColRanges.size() - 1).getY() < imageMaxColumn) {
                    // ascend rows to see if any have last colRange == imageMaxColumn
                    for (int row = skyRowMinMax[0] + 1; row <= skyRowMinMax[1]; row++) {
                        List<PairInt> rcr = skyRowColRanges.get(Integer.valueOf(row));                        
                        if (!rcr.isEmpty() && 
                            (rcr.get(rcr.size() - 1).getY() == imageMaxColumn)) {
                            fillInRightCorner = true;
                            break;
                        }
                    }
                }
                if (fillInRightCorner) {
                    rowColRanges.get(rowColRanges.size() - 1).setY(imageMaxColumn);
                }
            }
        }
        
        // ---- fill in any gaps in columns in bottom row if top row is imageMaxRow -----
        if (skyRowMinMax[1] == imageMaxRow) {
            
            List<PairInt> rowColRanges = skyRowColRanges.get(
                Integer.valueOf(imageMaxRow));
            
            if (rowColRanges.size() > 1) {
                rowColRanges.get(0).setY(rowColRanges.get(
                    rowColRanges.size() - 1).getY());
                for (int i = rowColRanges.size() - 1; i > 0; i--) {
                    rowColRanges.remove(i);
                }
            }
            
            if (!rowColRanges.isEmpty()) {
                
                // fill in corners if the gap along that image boundary column
                // is bounded at a higher row.
                boolean fillInLeftCorner = false;
                if (rowColRanges.get(0).getX() > 0) {
                    // ascend rows to see if any have first colRange starting at 0
                    for (int row = skyRowMinMax[1] - 1; row > skyRowMinMax[0]; row--) {
                        List<PairInt> rcr = skyRowColRanges.get(Integer.valueOf(row));                        
                        if (!rcr.isEmpty() && (rcr.get(0).getX() == 0)) {
                            fillInLeftCorner = true;
                            break;
                        }
                    }
                }
                if (fillInLeftCorner) {
                    rowColRanges.get(0).setX(0);
                }
            
                boolean fillInRightCorner = false;
                if (rowColRanges.get(rowColRanges.size() - 1).getY() < imageMaxColumn) {
                    // ascend rows to see if any have last colRange == imageMaxColumn
                    for (int row = skyRowMinMax[1] - 1; row > skyRowMinMax[0]; row--) {
                        List<PairInt> rcr = skyRowColRanges.get(Integer.valueOf(row));
                        if (!rcr.isEmpty() && (
                            rcr.get(rcr.size() - 1).getY() == imageMaxColumn)) {
                            fillInRightCorner = true;
                            break;
                        }
                    }
                }
                if (fillInRightCorner) {
                    rowColRanges.get(rowColRanges.size() - 1).setY(imageMaxColumn);
                }
            }
        }
        
        // --- fill in leftmost and rightmost image column gaps if bounded by sky ---
        int minRowForFirstColumn = Integer.MAX_VALUE;
        int maxRowForFirstColumn = Integer.MIN_VALUE;
        int minRowForLastColumn = Integer.MAX_VALUE;
        int maxRowForLastColumn = Integer.MIN_VALUE;
        for (int row = skyRowMinMax[0]; row <= skyRowMinMax[1]; row++) {
            List<PairInt> rowColRanges = skyRowColRanges.get(Integer.valueOf(row));
            if (rowColRanges.isEmpty()) {
                continue;
            }
            int firstColRangeX = rowColRanges.get(0).getX();
            int lastColRangeY = rowColRanges.get(rowColRanges.size() - 1).getY();
            if (firstColRangeX == 0) {
                if (row < minRowForFirstColumn) {
                    minRowForFirstColumn = row;
                }
                if (row > maxRowForFirstColumn) {
                    maxRowForFirstColumn = row;
                }
            }
            if (lastColRangeY == imageMaxColumn) {
                if (row < minRowForLastColumn) {
                    minRowForLastColumn = row;
                }
                if (row > maxRowForLastColumn) {
                    maxRowForLastColumn = row;
                }
            }
        }
        if ((minRowForFirstColumn != Integer.MAX_VALUE) &&
            (maxRowForFirstColumn != Integer.MIN_VALUE)) {
            
            for (int row = minRowForFirstColumn; row <= maxRowForFirstColumn; row++) {
                
                List<PairInt> rowColRanges = skyRowColRanges.get(Integer.valueOf(row));
                
                if (rowColRanges.isEmpty()) {
                    rowColRanges.add(new PairInt(0, 0));
                } else {
                    rowColRanges.get(0).setX(0);
                }
            }
        }
        if ((minRowForLastColumn != Integer.MAX_VALUE) &&
            (maxRowForLastColumn != Integer.MIN_VALUE)) {
            
            for (int row = minRowForLastColumn; row <= maxRowForLastColumn; row++) {
                
                List<PairInt> rowColRanges = skyRowColRanges.get(Integer.valueOf(row));
                
                if (rowColRanges.isEmpty()) {
                    rowColRanges.add(new PairInt(imageMaxColumn, imageMaxColumn));
                } else {
                    rowColRanges.get(rowColRanges.size() - 1).setY(imageMaxColumn);
                }
            }
        }
    }

    /**
     * a placeholder for ability to select the seeds of the sky from the
     * zeroPointLists.  This might allow user interaction in the future.
     * 
     * 
     * @param zeroPointLists
     * @param colorImg
     * @param thetaImg
     * @return 
     */
    protected List<PairIntArray> selectFromSkySeeds(
        List<PairIntArray> zeroPointLists, ImageExt colorImg, 
        GreyscaleImage thetaImg) {
        
        throw new UnsupportedOperationException("Not supported yet.");
    }

    protected List<PairIntArray> extractAndSmoothSkylinePoints(
        Set<PairInt> skyPoints, GreyscaleImage gradientXY) {
        
        /*
        could use border points directly followed by a line thinner 
        and then the EdgeExtractorWithJunctions, but
        a smoother curve is produced using a gradient produced from the
        mask (that is use a CannyEdgeDetector) followed by a line thinner 
        and then the EdgeExtractorWithJunctions.
        */

        boolean useGradient = false;
        
        if (useGradient) {
            
            GreyscaleImage mask = gradientXY.createWithDimensions();
            mask.fill(250);
            for (PairInt p : skyPoints) {
                int x = p.getX();
                int y = p.getY(); 
                mask.setValue(x, y, 0);
            }

            CannyEdgeFilter cFilter = new CannyEdgeFilter();
            cFilter.applyFilter(mask);

            IEdgeExtractor edgeExtractor = new EdgeExtractorWithJunctions(mask);
            edgeExtractor.removeShorterEdges(true);
            List<PairIntArray> edges = edgeExtractor.findEdges();

            return edges;
        }
        
        Set<PairInt> outputBorderPoints = new HashSet<PairInt>();
        
        Set<PairInt> outputEmbeddedGapPoints = new HashSet<PairInt>();
                
        getEmbeddedAndBorderPoints(skyPoints, gradientXY.getWidth(),
            gradientXY.getHeight(), outputEmbeddedGapPoints,
            outputBorderPoints);
       
        PostLineThinnerCorrections pslt = new PostLineThinnerCorrections();
        pslt.correctForArtifacts(outputBorderPoints, gradientXY.getWidth(), 
            gradientXY.getHeight());
         
        GreyscaleImage output = gradientXY.createWithDimensions();
        for (PairInt p : outputBorderPoints) {
            output.setValue(p.getX(), p.getY(), 1);
        }
        
        AbstractEdgeExtractor edgeExtractor = 
            new EdgeExtractorWithJunctions(output);
        
        //edgeExtractor.removeShorterEdges(true);
        
        List<PairIntArray> edges = edgeExtractor.findEdges();
        
        // smooth by large number.  TODO: consider a resoltuion dep factor
        AverageUtil avgUtil = new AverageUtil();
        int k = 15;
        
        output = gradientXY.createWithDimensions();       
        
        for (int i = 0; i < edges.size(); ++i) {
            PairIntArray edge = edges.get(i);
            if (edge.getN() >= k) {
                edge = avgUtil.calculateBoxCarAverage(edges.get(i), k);
                for (int j = 0; j < edge.getN(); ++j) {
                    output.setValue(edge.getX(j), edge.getY(j), 1);
                }
            }
        }
        
        pslt = new PostLineThinnerCorrections();
        pslt.correctForArtifacts(output);        
        edgeExtractor = new EdgeExtractorWithJunctions(output);
        edgeExtractor.removeShorterEdges(true);
        edges = edgeExtractor.findEdges();     
        
        return edges;
    }

    public class RemovedSets {
        private Set<PairInt> removedNonCloudColors = new HashSet<PairInt>();
        private Set<PairInt> highContrastRemoved = new HashSet<PairInt>();
        private Set<PairInt> reflectedSunRemoved = new HashSet<PairInt>();
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
    
    public SunFinder getSunFinderResults() {
        return sunFinder;
    }
    
    public RainbowFinder getRainbowFinderResults() {
        return rainbowFinder;
    }
    
    public List<PairIntArray> getSkylineEdges() {
        return skylineEdges;
    }
}
