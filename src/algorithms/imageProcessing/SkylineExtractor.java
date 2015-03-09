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
        GreyscaleImage gradientXY, Image originalImage,
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
        GreyscaleImage gradientXY,
        Image originalColorImage, CannyEdgeFilterSettings edgeSettings,
        PairIntArray outputSkyCentroid) throws 
        IOException, NoSuchAlgorithmException {
        
        if (theta == null) {
            throw new IllegalArgumentException("theta cannot be null");
        }
                
        Image colorImg = originalColorImage;
        GreyscaleImage thetaImg = theta;
        GreyscaleImage gXYImg = gradientXY;
        
        int binFactor = determineBinFactorForSkyMask(theta.getNPixels());

        log.info("binFactor=" + binFactor);
        
        ImageProcessor imageProcessor = new ImageProcessor();

        if (binFactor > 1) {
            thetaImg = imageProcessor.binImage(theta, binFactor);
            colorImg = imageProcessor.binImage(originalColorImage, binFactor);
            gXYImg = imageProcessor.binImage(gradientXY, binFactor);
        }

        List<PairIntArray> zeroPointLists = getSortedContiguousZeros(thetaImg);
        
        if (zeroPointLists.isEmpty()) {
            
            GreyscaleImage mask = theta.createWithDimensions();
               
            // return an image of all 1's
            mask.fill(1);
            
            return mask;
        }
      
        List<PairIntArray> removedNonCloudColors = 
            removeSetsWithNonCloudColors(zeroPointLists, colorImg, thetaImg,
            true, 0);
        
        if (zeroPointLists.isEmpty()) {
            
            GreyscaleImage mask = theta.createWithDimensions();
               
            // return an image of all 1's
            mask.fill(1);
            
            return mask;
        }
        
        //TODO:  assumes that largest smooth component of image is sky.
        // if sky is small and a foreground object is large and featureless
        // and not found as dark, this will fail. 
        // will adjust for that one day, possibly with color validation
        reduceToLargest(zeroPointLists);
        
        if (zeroPointLists.isEmpty()) {
            
            GreyscaleImage mask = theta.createWithDimensions();
               
            // return an image of all 1's
            mask.fill(1);
            
            return mask;
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
            GreyscaleImage mask = theta.createWithDimensions();
            // return an image of all 1's
            mask.fill(1);
            return mask;
        }
        
        int nAfterHighContrastRemoval = count(zeroPointLists);
                
        Set<PairInt> reflectedSunRemoved = removeReflectedSun(zeroPointLists, 
            colorImg, thetaImg);
        
        if (zeroPointLists.isEmpty()) {
            GreyscaleImage mask = theta.createWithDimensions();
            mask.fill(1);
            return mask;
        }
        
        Set<PairInt> points = combine(zeroPointLists);
       
        int valueToSubtract = extractSkyFromGradientXY(gXYImg, points,
            highContrastRemoved);
        
        if (binFactor > 1) {
            
            thetaImg = theta;
            colorImg = originalColorImage;
            gXYImg = gradientXY;
            
            points = imageProcessor.unbinZeroPointLists(points, binFactor);
            
            removedNonCloudColors = imageProcessor.unbinZeroPointLists(removedNonCloudColors, 
                binFactor);
            
            highContrastRemoved = imageProcessor.unbinZeroPointLists(highContrastRemoved,
                binFactor);
            
            reflectedSunRemoved = imageProcessor.unbinZeroPointLists(reflectedSunRemoved,
                binFactor);
            
            if (!highContrastRemoved.isEmpty()) {
                removeBinFactorArtifacts(points, binFactor,
                    thetaImg.getWidth(), thetaImg.getHeight(),
                    thetaImg.getXRelativeOffset(), thetaImg.getYRelativeOffset());
            }
        }
        
        //now the coordinates in zeroPointLists are w.r.t. thetaImg
        
        PerimeterFinder perimeterFinder = new PerimeterFinder();
        int[] skyRowMinMax = new int[2];
        Map<Integer, PairInt> skyRowColRange = perimeterFinder.find(points, 
            skyRowMinMax);
        fillInRightAndLowerBoundarySkyPoints(points, binFactor, 
            skyRowColRange, skyRowMinMax, colorImg,
            thetaImg.getWidth(), thetaImg.getHeight(),
            thetaImg.getXRelativeOffset(), thetaImg.getYRelativeOffset());
        
        GreyscaleImage mask = gradientXY.createWithDimensions();
        mask.fill(1);
        for (PairInt p : points) {
            int x = p.getX();
            int y = p.getY(); 
            mask.setValue(x, y, 0);
        }
        
        int nSkyPointsBeforeFindClouds = points.size();
        
        Map<Integer, PixelColors> pixelColorsMap = new 
            HashMap<Integer, PixelColors>();
        
        Map<PairInt, Set<PixelColors> > skyColorsMap = new HashMap<PairInt, 
            Set<PixelColors> >();
        
        //TODO: this will be revised when have narrowed down which color spaces
        // to use.
                
        populatePixelColorMaps(points, colorImg, mask, pixelColorsMap, 
            skyColorsMap);
        
        findClouds(points, new HashSet<PairInt>(), colorImg, mask, 
            pixelColorsMap, skyColorsMap);

        mergeOverlappingPreviouslyRemovedNonCloudColors(points, 
            removedNonCloudColors);

        GroupPixelColors allSkyColor = new GroupPixelColors(points,
            colorImg, thetaImg.getXRelativeOffset(), 
            thetaImg.getYRelativeOffset());
        
        boolean skyIsDarkGrey = skyIsDarkGrey(allSkyColor);

        Set<PairInt> sunPoints = new HashSet<PairInt>();
        double[] ellipFitParams = findSunConnectedToSkyPoints(points, 
            reflectedSunRemoved, colorImg, thetaImg.getXRelativeOffset(), 
            thetaImg.getYRelativeOffset(), skyIsDarkGrey, sunPoints);
        
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
            findRainbowPoints(points, reflectedSunRemoved, colorImg, 
                thetaImg.getXRelativeOffset(), thetaImg.getYRelativeOffset(),
                skyIsDarkGrey) :
            new HashSet<PairInt>();
        
        // find embedded non-sky pixels that should be sky:
        Set<PairInt> embeddedPoints = findEmbeddedNonPointsExtendedToImageBounds(
            points, mask, 1);
        
debugPlot(points, colorImg, mask.getXRelativeOffset(), mask.getYRelativeOffset(), 
"after_first_0");

        if (!embeddedPoints.isEmpty()) {

            HistogramHolder[] brightnessHistogram = new HistogramHolder[1];

            // brightest sky is in bin [2], and dimmest is in [0]
            GroupPixelColors[] skyPartitionedByBrightness = 
                partitionInto3ByBrightness(points, originalColorImage, 
                thetaImg.getXRelativeOffset(), thetaImg.getYRelativeOffset(), 
                brightnessHistogram);

            findEmbeddedClouds(embeddedPoints, points, highContrastRemoved,
                colorImg, mask,
                pixelColorsMap, skyColorsMap,
                brightnessHistogram, skyPartitionedByBrightness);

debugPlot(points, colorImg, mask.getXRelativeOffset(), mask.getYRelativeOffset(), 
"after_0");

//plotSkyColor(points, colorImg, mask);

            if ((reflectedSunRemoved.size() < 50) &&
                ((nBeforeHighContrastRemoval - nAfterHighContrastRemoval) <= 
                (int)(((float)nBeforeHighContrastRemoval)*0.1f))) {
            
                // also, do not perform this if we can tell that the sky is mostly
                // found except the nBinFactor pixels near the border.
                // (wanting to avoid overrunning low contrast skylines such as
                // snow covered peaks under a hazy sky, etc)
                if (embeddedPoints.size() > 100) {
            
                    // histograms of hue and saturation
 
                    int nBefore = points.size();
                    //Set<PairInt> copyOfPoints = new HashSet<PairInt>(points);
                    //GreynMascaleImage copyOfMask = mask.copyImage();

                    findSeparatedClouds(sunPoints, points, highContrastRemoved,
                        colorImg, mask,
                        pixelColorsMap, skyColorsMap,
                        brightnessHistogram, skyPartitionedByBrightness);
                    
                    int nPointsAdded = points.size() - nBefore;
                    float frac = (float)nPointsAdded/(float)points.size();
                    log.info("number of points added for separated search=" +
                        nPointsAdded + " to " + points.size() + " points (" +
                        frac + "%)");
                       
debugPlot(points, colorImg, mask.getXRelativeOffset(), mask.getYRelativeOffset(), 
"after_1");          

                    if (reflectedSunRemoved.size() > 0) {
                        
                        reduceToLargestContiguousGroup(points, mask, 10, 
                            new HashSet<PairInt>());                        
                    }

debugPlot(points, colorImg, mask.getXRelativeOffset(), mask.getYRelativeOffset(), 
"after_added_separated_clouds");

                }
            }
        }

        Set<PairInt> exclude = new HashSet<PairInt>();
        exclude.addAll(highContrastRemoved);
        exclude.addAll(rainbowPoints);
        exclude.addAll(sunPoints);

        if (embeddedPoints.size() > 0) {

            HistogramHolder[] brightnessHistogram = new HistogramHolder[1];

            // brightest sky is in bin [2], and dimmest is in [0]
            GroupPixelColors[] skyPartitionedByBrightness = 
                partitionInto3ByBrightness(points, originalColorImage, 
                thetaImg.getXRelativeOffset(), 
                thetaImg.getYRelativeOffset(), brightnessHistogram);

            findEmbeddedClouds(embeddedPoints, points, exclude,
                colorImg, mask,
                pixelColorsMap, skyColorsMap,
                brightnessHistogram, skyPartitionedByBrightness);

debugPlot(points, colorImg, mask.getXRelativeOffset(), mask.getYRelativeOffset(), 
"after_2");

        }

log.info("number of sunPoints=" + sunPoints.size() + " "
    + "reflectedSunRemoved.size()=" + reflectedSunRemoved.size());

        if (sunPoints.isEmpty()) {
            
            growForLowContrastLimits(points, exclude, colorImg, mask);
        
            // TODO: might need to expand these a little
            Set<PairInt> placeholdingPoints = new HashSet<PairInt>();
            placeholdingPoints.addAll(sunPoints);
            placeholdingPoints.addAll(rainbowPoints);

            reduceToLargestContiguousGroup(points, mask, 1000, 
                placeholdingPoints);

        }
       
        if (!sunPoints.isEmpty()) {
            correctSkylineForSun(sunPoints, points, colorImg, mask, gradientXY);
        }/* else if (!rainbowPoints.isEmpty()) {
            correctSkylineForRainbow(rainbowPoints, points, colorImg, mask);
        }*/

debugPlot(points, colorImg, mask.getXRelativeOffset(), 
    mask.getYRelativeOffset(), "final");

        //TODO: refine this number comparison
        float f = (float)points.size()/(float)nSkyPointsBeforeFindClouds;
        if (f > 2.0) {
        //    findMoreClouds(points, colorImg, mask);
        }
        
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

        if (mask != null) {
            imageProcessor.removeSpurs(mask);
        }         
   
        return mask;
    }
    
    public Set<PairInt> combine(List<PairIntArray> points) {
        Set<PairInt> set = new HashSet<PairInt>();
        for (PairIntArray p : points) {
            for (int i = 0; i < p.getN(); i++) {
                int x = p.getX(i);
                int y = p.getY(i);
                PairInt pi = new PairInt(x, y);
                set.add(pi);
            }
        }
        return set;
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
        
        if (nMaxGroupN > 10000000) {
            MultiArrayMergeSort.sortByDecr(groupN, groupIndexes);
        } else {
            CountingSort.sortByDecr(groupN, groupIndexes, nMaxGroupN);
        }
        
        List<Integer> groupIds = new ArrayList<Integer>();
        groupIds.add(Integer.valueOf(groupIndexes[0]));
        
        if (nGroups > 1) {
                        
            float n0 = (float)groupN[0];
            
            for (int i = 1; i < groupN.length; i++) {
                
                if (limitToLargest) {
                    float number = groupN[i];

                    float frac = number/n0;
    System.out.println(number);    
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
    System.out.println(number + " n0=" + n0);    
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
     * any removed connected sky as long as the majority of sky points remain.
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
     * @param gradientXY
     * @param skyPoints
     * @param excludeThesePoints
     * @return 
     */
    public int extractSkyFromGradientXY(GreyscaleImage gradientXY,
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
            
            //TODO: remove this section or move it into an aspect
                
            // === count number of embedded groups of non-zeros in skyPoints ===
            
            Set<PairInt> embeddedPoints = new HashSet<PairInt>();
            
            int nCorrectedEmbeddedGroups = extractEmbeddedGroupPoints(
                skyPoints, gXY2, embeddedPoints);
            
            log.info("nIter=" + nIter + ")" 
                + " nCorrectedEmbeddedGroups=" + nCorrectedEmbeddedGroups
                + " nEmbeddedPixels=" + embeddedPoints.size()
                + " out of " + skyPoints.size()
                + " (level=" + ((float)gXYValues.getY(lastHistIdx)/originalMaxValue) 
                + " subtract=" + subtract + " out of max=" + gXYValues.getX(0)
                + ")"
            );
                        
            nIter++;
        }
/*
try {
    Image img1 = gradientXY.copyImageToGreen();
    ImageIOHelper.addToImage(skyPoints, 0, 0, img1);
    ImageDisplayer.displayImage("sky points subtract=" + subtract, img1);
} catch (IOException ex) {
    log.severe(ex.getMessage());
}
*/         
        return subtract;
    }

    private void growZeroValuePoints(Set<PairInt> points, 
        Set<PairInt> excludeThesePoints, GreyscaleImage gradientXY) {
  
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
     * water, their location in x,y must be fittable by an ellipse, else they 
     * may not be found as sun points.
     * @param skyPoints
     * @param reflectedSunRemoved previously removed points assumed to be
     * reflected light
     * @param clr
     * @param xOffset
     * @param yOffset 
     * @param skyIsDarkGrey 
     * @param outputYellowPoints 
     * @return the extracted sun points that are connected to skyPoints
     */
    protected double[] findSunConnectedToSkyPoints(Set<PairInt> skyPoints, 
        Set<PairInt> reflectedSunRemoved,
        Image clr, int xOffset, int yOffset, boolean skyIsDarkGrey,
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

                int r = clr.getR(xx, yy);
                int g = clr.getG(xx, yy);
                int b = clr.getB(xx, yy);

                if (sunColors.isSunCenterColor(r, g, b)) {
                    
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

                int r = clr.getR(vX, vY);
                int g = clr.getG(vX, vY);
                int b = clr.getB(vX, vY);

                if (sunColors.isSunCenterColor(r, g, b)) {
                    
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
     * this assumes findSunConnectedToSkyPoints has been run to populate
     * sunPoints and that those have not been filtered for a non-circular or
     * elliptical shape.  If there are no points in sunPoints (which are bright
     * yellow), then this method will not search for a rainbow.
     * @param skyPoints
     * @param reflectedSunRemoved
     * @param colorImg
     * @param xRelativeOffset
     * @param yRelativeOffset
     * @return 
     */
    private Set<PairInt> findRainbowPoints(Set<PairInt> skyPoints, 
        Set<PairInt> reflectedSunRemoved,
        Image colorImg, int xOffset, int yOffset, boolean skyIsDarkGrey) {
                
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

    private float[] calculateAvgAndStDevOfDiffs(Map<PairInt, PairFloat> 
        contrastAndColorMap) {
        
        // average difference from neighbors
        double avgContrast = 0;
        
        double avgColor = 0;
        
        int count = 0;
        
        Iterator<Map.Entry<PairInt, PairFloat> > iter = 
            contrastAndColorMap.entrySet().iterator();
        
        while (iter.hasNext()) {
            
            Map.Entry<PairInt, PairFloat> entry = iter.next();
            
            PairInt i = entry.getKey();
            int x = i.getX();
            int y = i.getY();
            
            for (int xx = (x - 1); xx <= (x + 1); xx++) {
                for (int yy = (y - 1); yy <= (y + 1); yy++) {
                    
                    PairInt j = new PairInt(xx, yy);
                    
                    if (contrastAndColorMap.containsKey(j)) {
                        
                        PairFloat iCC = entry.getValue();
                        
                        PairFloat jCC = contrastAndColorMap.get(j);
                        
                        float diffContrast = Math.abs(iCC.getX() - jCC.getX());
                        
                        float diffColor = Math.abs(iCC.getY() - jCC.getY());
                        
                        avgContrast += diffContrast;
                        
                        avgColor += diffColor;
                        
                        count++;
                    }
                }
            }            
        }
        
        avgContrast /= (double)count;
        
        avgColor /= (double)count;
        
        // standard deviation of avg difference from neighbors
        double stDevContrast = 0;
        
        double stDevColor = 0;
        
        iter = contrastAndColorMap.entrySet().iterator();
        
        while (iter.hasNext()) {
            
            Map.Entry<PairInt, PairFloat> entry = iter.next();
            
            PairInt i = entry.getKey();
            int x = i.getX();
            int y = i.getY();
            
            for (int xx = (x - 1); xx <= (x + 1); xx++) {
                for (int yy = (y - 1); yy <= (y + 1); yy++) {
                    
                    PairInt j = new PairInt(xx, yy);
                    
                    if (contrastAndColorMap.containsKey(j)) {
                        
                        PairFloat iCC = entry.getValue();
                        
                        PairFloat jCC = contrastAndColorMap.get(j);
                        
                        float diffContrast = Math.abs(iCC.getX() - jCC.getX());
                        diffContrast -= avgContrast;
                        
                        float diffColor = Math.abs(iCC.getY() - jCC.getY());
                        diffColor -= avgColor;
                        
                        stDevContrast += (diffContrast * diffContrast);
                        
                        stDevColor += (diffColor * diffColor);
                    }
                }
            } 
        }
        
        stDevContrast = Math.sqrt(stDevContrast/((double)count - 1));
        
        stDevColor = Math.sqrt(stDevColor/((double)count - 1));
        
        return new float[]{(float)avgContrast, (float)stDevContrast, 
            (float)avgColor, (float)stDevColor};
    }
    
    private void addBackMissingZeros(Set<PairInt> zeroPoints,
        GreyscaleImage gXYImg, int binFactor, int valueToSubtract) {
        
        int width = gXYImg.getWidth();
        int height = gXYImg.getHeight();
        
        GreyscaleImage img = gXYImg.copyImage();
        MatrixUtil.add(img.getValues(), -1*valueToSubtract);
        
        Set<PairInt> addPoints = new HashSet<PairInt>();
        
        for (PairInt p : zeroPoints) {

            int x = p.getX();
            int y = p.getY();

            for (int c = (x - binFactor); c <= (x + binFactor); c++) {
                if ((c < 0) || (c > (width - 1))) {
                    continue;
                }
                for (int r = (y - binFactor); r <= (y + binFactor); r++) {
                    if ((r < 0) || (r > (height - 1))) {
                        continue;
                    }
                    if ((c == x) && (r == y)) {
                        continue;
                    }

                    int neighborIdx = img.getIndex(c, r);

                    Integer index = Integer.valueOf(neighborIdx);

                    if (addPoints.contains(index)) {
                        continue;
                    }
                    if (zeroPoints.contains(index)) {
                        continue;
                    }

                    int v = img.getValue(c, r);
                    if (v == 0) {
                        addPoints.add(new PairInt(c, r));
                    }
                }
            }
        }
        
        zeroPoints.addAll(addPoints);
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

    public int extractEmbeddedGroupPoints(
        List<Set<PairInt>> embeddedGroups, Map<Integer, PairInt> rowColRange, 
        int[] rowMinMax, Set<PairInt> outputEmbeddedPoints,
        int xMinImage, int xMaxImage, int yMinImage, int yMaxImage) {
        
        int nCorrectedEmbeddedGroups = 0;
        
        PerimeterFinder finder = new PerimeterFinder();
        
        for (int gId = 0; gId < embeddedGroups.size(); gId++) {

            Set<PairInt> groupPoints = embeddedGroups.get(gId);

            int[] gRowMinMax = new int[2];
            Map<Integer, PairInt> gRowColRange = finder.find(groupPoints, 
                gRowMinMax);

            boolean unbounded = isPerimeterUnbound(gRowColRange, gRowMinMax,
                rowColRange, rowMinMax, xMinImage, xMaxImage, yMinImage, yMaxImage);

            if (!unbounded) {
                
                nCorrectedEmbeddedGroups++;
                
                outputEmbeddedPoints.addAll(groupPoints);
            }
        }
        
        return nCorrectedEmbeddedGroups;
    }
    
    private int extractEmbeddedGroupPoints(Set<PairInt> skyPoints, 
        GreyscaleImage gXY2, Set<PairInt> outputEmbeddedPoints) {
        
        PerimeterFinder finder = new PerimeterFinder();
        int[] rowMinMax = new int[2];
        Map<Integer, PairInt> rowColRange = finder.find(skyPoints, rowMinMax);
        DFSContiguousValueFinder contiguousNonZeroFinder = 
            new DFSContiguousValueFinder(gXY2);
        contiguousNonZeroFinder.findEmbeddedGroupsNotThisValue(0, rowColRange,
            rowMinMax);
        int nEmbeddedGroups = contiguousNonZeroFinder.getNumberOfGroups();

        List<Set<PairInt>> embeddedGroups = new ArrayList<Set<PairInt>>();
        
        for (int gId = 0; gId < nEmbeddedGroups; gId++) {
            
            Set<PairInt> groupPoints = new HashSet<PairInt>();

            contiguousNonZeroFinder.getXY(gId, groupPoints);
            
            embeddedGroups.add(groupPoints);
            
            outputEmbeddedPoints.addAll(groupPoints);
        }
                
        return nEmbeddedGroups;
    }

    /**
     * fill in the rightmost sky points that would be missing due to down sizing
     * the image to find sky points then up sizing the image to current scale
     * and do the same for the lowest sky points in the image (highest y).
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
    private void fillInRightAndLowerBoundarySkyPoints(Set<PairInt> skyPoints, 
        int binFactor, Map<Integer, PairInt> skyRowColRange, int[] skyRowMinMax,
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
            
            PairInt cRange = skyRowColRange.get(rowIndex);
            
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
                
                skyRowColRange.put(rowIndex, new PairInt(cRange.getX(), x));
            }
        }

        // lower boundary
        int diff = imageHeight - skyRowMinMax[1];
        if ((diff <= binFactor) && (diff > 0)) {
            
            PairInt cRange = skyRowColRange.get(Integer.valueOf(skyRowMinMax[1]));
            
            for (int row = (skyRowMinMax[1] + 1); row < imageHeight; row++) {
                PairInt colRange = new PairInt(cRange.getX(), cRange.getY());
                skyRowColRange.put(Integer.valueOf(row), colRange);
                
                for (int col = colRange.getX(); col <= colRange.getY(); col++) {
                    skyPoints.add(new PairInt(col, row));
                }
            }
        }
    }

    private void addPointsAndUpdateRowColRange(Set<PairInt> skyPoints, 
        Set<PairInt> sunPoints, Map<Integer, PairInt> skyRowColRange, 
        int[] skyRowMinMax) {
        
        if (skyPoints == null) {
            throw new IllegalArgumentException("skyPoints cannot be null");
        }
        if (sunPoints == null) {
            throw new IllegalArgumentException("sunPoints cannot be null");
        }
        if (skyRowColRange == null) {
            throw new IllegalArgumentException("skyRowColRange cannot be null");
        }
        if (skyRowMinMax == null) {
            throw new IllegalArgumentException("skyRowMinMax cannot be null");
        }
        
        if (sunPoints.isEmpty()) {
            return;
        }
        
        skyPoints.addAll(sunPoints);
        
        PerimeterFinder finder = new PerimeterFinder();
        Map<Integer, PairInt> skyRowColRange2 = finder.find(skyPoints, 
            skyRowMinMax);
        
        skyRowColRange.clear();
        
        skyRowColRange.putAll(skyRowColRange2);
        
    }

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
     * @param pixelColorsMap
     * @param skyColorsMap 
     */
    private void findClouds(Set<PairInt> skyPoints, Set<PairInt> excludePoints,
        Image originalColorImage, GreyscaleImage mask,
        Map<Integer, PixelColors> pixelColorsMap,
        Map<PairInt, Set<PixelColors> > skyColorsMap
        ) {
        
        int maskWidth = mask.getWidth();
        int maskHeight = mask.getHeight();
        
        Stack<PairInt> cloudStack = new Stack<PairInt>();
        
        Set<PairInt> candidateCloudPoints = new HashSet<PairInt>();
       
        int xOffset = mask.getXRelativeOffset();
        int yOffset = mask.getYRelativeOffset();
              
        int[] dxs = new int[]{-1, -1,  0,  1, 1, 1, 0, -1};
        int[] dys = new int[]{ 0, -1, -1, -1, 0, 1, 1,  1};
        
        // TODO: some adjustment should be made for the perimeter of the
        // sun points to not erode the skyline.  
                
        for (PairInt skyPoint : skyPoints) {
            
            int x = skyPoint.getX();
            int y = skyPoint.getY();
            Integer index = Integer.valueOf(mask.getIndex(x, y));
        
            for (int k = 0; k < dxs.length; k++) {
            
                int xx = x + dxs[k];
                int yy = y + dys[k];
                
                if ((xx < 0) || (xx > (maskWidth - 1)) || (yy < 0) || 
                    (yy > (maskHeight - 1))) {
                    continue;
                }
                
                PairInt p = new PairInt(xx, yy);
                
                if (skyPoints.contains(p) || excludePoints.contains(p)) {
                    continue;
                }
                
                //add colors for index to the local sky points entry for p
                populatePixelColorMaps(p, index, originalColorImage, 
                    xOffset, yOffset, pixelColorsMap, skyColorsMap);
        
                if (candidateCloudPoints.contains(p)) {
                    continue;
                }
                
                // include perimeter points (they're missing from  
                // skyPoints because those are derived from theta and theta is
                // the product of 2 gaussian convolutions.  the convolutions
                // widen features in gradientXY.  see caveat above sunPoints.
                // it's safe to add these without a color or contrast check:
                candidateCloudPoints.add(p);
                
                cloudStack.add(p);
            }
        }
    
        GroupPixelColors allSkyColor = new GroupPixelColors(skyPoints,
            originalColorImage, xOffset, yOffset);
        
        double rDivB = allSkyColor.getAvgRed() / allSkyColor.getAvgBlue();
        boolean skyIsRed = (rDivB > 1);
       
        log.info("==> r/b=" + rDivB
            + " redStdev=" + allSkyColor.getStdDevRed()
            + " blueStDev=" + allSkyColor.getStdDevBlue());
        
        Map<PairInt, GroupPixelColors> localSkyColors = new
            HashMap<PairInt, GroupPixelColors>();
        
        Set<PairInt> visited = new HashSet<PairInt>();
        visited.add(cloudStack.peek());
       
        //dxs = new int[]{-1,  0, 1, 0};
        //dys = new int[]{ 0, -1, 0, 1};
        
        while (!cloudStack.isEmpty()) {

            PairInt uPoint = cloudStack.pop();
            
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

                if (visited.contains(vPoint) || skyPoints.contains(vPoint)
                    || candidateCloudPoints.contains(vPoint) || 
                    excludePoints.contains(vPoint)) {
                    continue;
                }

                visited.add(vPoint);
                    
                GroupPixelColors localSky = localSkyColors.get(uPoint);
                if (localSky == null) {
                    // this should never be null:
                    Set<PixelColors> skyColors = skyColorsMap.get(uPoint);
                    localSky = new GroupPixelColors(skyColors);
                    localSkyColors.put(uPoint, localSky);
                }

                int rV = originalColorImage.getR(vX + xOffset, vY + yOffset);
                int gV = originalColorImage.getG(vX + xOffset, vY + yOffset);
                int bV = originalColorImage.getB(vX + xOffset, vY + yOffset);

                float totalRGBV = rV + gV + bV;
                
                /*is this within color and contrast limits?  if yes,
                   add to candidateCloudPoints 
                   add to cloudStack
                   add skyColors to candidateSkyColorsMap for vPoint
                   add localSky to localSkyColors for vPoint
                */

                double contrastV = localSky.calcContrastToOther(rV, gV, bV);

                //TODO: this should use CIE xy chromaticity instead of rgb for
                // similar color independent of brightness
                
                double colorDiffV = localSky.calcColorDiffToOther(rV, gV, bV);

                double skyStDevContrast = localSky.getStdDevContrast();

                double skyStDevColorDiff = localSky.getStdDevColorDiff();

                boolean doNotAddToStack = false;
                 
 //TODO: this needs adjustments...
                float rPercentV = (float)rV/totalRGBV;
                float gPercentV = (float)gV/totalRGBV;
                float bPercentV = (float)bV/totalRGBV;
                
                boolean isBrown = (Math.abs(rPercentV - 0.5) < 0.4)
                    && (Math.abs(gPercentV - 0.32) < 0.1)
                    && (Math.abs(bPercentV - 0.17) < 0.1);
              
                if (isBrown) {
                    
                    // trying to skip over foreground such as land or sunset + water
                    
                    float[] hsb = new float[3];
                    Color.RGBtoHSB(rV, gV, bV, hsb);
                    
                    if ((colorDiffV > 15*skyStDevColorDiff) && (hsb[2] < 0.5)) {
                        
                        if (hsb[2] > 0.4) {
                                                        
                            continue;
                            
                        } else if ((colorDiffV > 50*skyStDevColorDiff) || 
                            (Math.abs(contrastV) > 10.*Math.abs(skyStDevContrast))
                            ) {
                                                        
                            continue;
                        }
                    }
                }
                
                if (Double.isInfinite(skyStDevContrast)) {
                    continue;
                }
                
                if (skyIsRed) {
                    
                    // if contrast is '+' and large, may be the skyline boundary
                    if ((skyStDevContrast != 0.)
                        &&
                        (
                            ((contrastV > 0.1) && (Math.abs(contrastV) > 3.*skyStDevContrast))
                            ||
                            ((contrastV > 0) && (colorDiffV > 15*skyStDevColorDiff))
                        )
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
/* 
if (((vX + xOffset) >= 588) && ((vX + xOffset) <= 588) && ((vY + yOffset) == 266)) {
    float[] hsb = new float[3];
    Color.RGBtoHSB(rV, gV, bV, hsb);
    String str = String.format(
        "*(%d, %d) contrastV=%f skyStDevContrast=%f colorDiffV=%f skyStDevColorDiff=%f hsb=%s",
        (vX + xOffset), (vY + yOffset),
        contrastV, skyStDevContrast, colorDiffV, skyStDevColorDiff,
        Arrays.toString(hsb)
    );
    log.info(str);
int z = 1;
}*/
                    //blue filters
                    if (
                        (
                            (contrastV < 0.05)
                            || 
                            (
                                (skyStDevContrast == 0.0) && (skyStDevColorDiff == 0.0)
                                && (contrastV < 0.01)
                            )
                        )
                        &&
                        (colorDiffV < 16)
                        ) {
                       
                    } else {

                        continue;
                    }
                }
                
                candidateCloudPoints.add(vPoint);

                if (!doNotAddToStack) {
                    cloudStack.add(vPoint);
                }

                localSkyColors.put(vPoint, localSky);
                skyColorsMap.put(vPoint, skyColorsMap.get(uPoint));

            }
        }
                
        skyPoints.addAll(candidateCloudPoints);
        
        for (PairInt p : skyPoints) {
            int x = p.getX();
            int y = p.getY();            
            mask.setValue(x, y, 0);
        }
    }

    private void findMoreClouds(Set<PairInt> skyPoints, Image originalColorImage, 
        GreyscaleImage mask) {
        
        /* for cloudy skies, findClouds can find the large majority of sky,
        but sometimes there are dark cloud bands separating brighter
        clouds or sky from the majority of connected sky pixels.
        findMoreClouds attempts to connect the nearby, but unconnected sky
        with the large majority of connected sky by looking at the color
        properties of skyPoints as 3 separate groups.
        */
        
        int xOffset = mask.getXRelativeOffset();
        int yOffset = mask.getYRelativeOffset();
        
        java.util.Stack<PairInt> stack = new java.util.Stack<PairInt>();
        
        //O(N_sky)
        for (PairInt p : skyPoints) {
            stack.add(p);
        }
        Set<PairInt> visited = new HashSet<PairInt>();
        visited.add(stack.peek());
        
        int[] counts = new int[3];
        GroupPixelColors[] allSkyColors = partitionInto3ByColorDifference(skyPoints,
            originalColorImage, xOffset, yOffset, counts);
        int maxPartionIdx = MiscMath.findYMaxIndex(counts);
        
        int maskWidth = mask.getWidth();
        int maskHeight = mask.getHeight();
        
        //int[] dxs = new int[]{-1,  0, 1, 0};
        //int[] dys = new int[]{ 0, -1, 0, 1};
        int[] dxs = new int[]{-1, -1,  0,  1, 1, 1, 0, -1};
        int[] dys = new int[]{ 0, -1, -1, -1, 0, 1, 1,  1};
        
        while (!stack.isEmpty()) {
            
            PairInt skyPoint = stack.pop();
            
            int x = skyPoint.getX();
            int y = skyPoint.getY();
        
            for (int k = 0; k < dxs.length; k++) {
                    
                int dx = dxs[k];
                int dy = dys[k];
                
                int xx = x + dx;
                int yy = y + dy;
                
                if ((xx < 0) || (xx > (maskWidth - 1)) || (yy < 0) || 
                    (yy > (maskHeight - 1))) {
                    continue;
                }
                               
                PairInt p = new PairInt(xx, yy);
                
                if (skyPoints.contains(p) || visited.contains(p)) {
                    continue;
                }
               
                visited.add(p);
                
                xx += xOffset;
                yy += yOffset;
                
                int rV = originalColorImage.getR(xx, yy);
                int gV = originalColorImage.getG(xx, yy);
                int bV = originalColorImage.getB(xx, yy);
  
                int ii = maxPartionIdx;
                float clrDiff = allSkyColors[ii].calcColorDiffToOther(rV, gV, bV);
                float avgClrDiff = allSkyColors[ii].getAvgColorDiff();
                double stDev = allSkyColors[ii].getStdDevColorDiff();
                // between 10.0 and 15.0
                if (Math.abs(clrDiff - avgClrDiff) <= 15.0*stDev) {
                    float contrast = allSkyColors[ii].calcContrastToOther(rV, gV, bV);
                    float avgContrast = allSkyColors[ii].getAvgContrast();
                    stDev = allSkyColors[ii].getStdDevContrast();
                    if (Math.abs(contrast - avgContrast) <= 15.0*stDev) {

                        skyPoints.add(p);

                        mask.setValue(xx, yy, 0);

                        stack.add(p);
                        
                        break;
                    }
                }
         
                if (false)
                log.info(String.format("(%d, %d) dClr=%f,%f,%f   dContrast=%f,%f,%f   (%f,%f,%f) (%f,%f,%f)",
                    xx, yy, 
                    allSkyColors[0].calcColorDiffToOther(rV, gV, bV),
                    allSkyColors[1].calcColorDiffToOther(rV, gV, bV),
                    allSkyColors[2].calcColorDiffToOther(rV, gV, bV),
                    
                    allSkyColors[0].calcContrastToOther(rV, gV, bV),
                    allSkyColors[1].calcContrastToOther(rV, gV, bV),
                    allSkyColors[2].calcContrastToOther(rV, gV, bV),
                    
                    (allSkyColors[0].calcColorDiffToOther(rV, gV, bV) 
                        - allSkyColors[0].getAvgColorDiff() ) /
                        allSkyColors[0].getStdDevColorDiff(),
                    (allSkyColors[1].calcColorDiffToOther(rV, gV, bV) 
                        - allSkyColors[1].getAvgColorDiff() ) /
                        allSkyColors[1].getStdDevColorDiff(),
                    (allSkyColors[2].calcColorDiffToOther(rV, gV, bV) 
                        - allSkyColors[2].getAvgColorDiff() ) /
                        allSkyColors[2].getStdDevColorDiff(),
                    
                    (allSkyColors[0].calcContrastToOther(rV, gV, bV)
                        - allSkyColors[0].getAvgContrast()) /
                        allSkyColors[0].getStdDevContrast(),
                    (allSkyColors[1].calcContrastToOther(rV, gV, bV)
                        - allSkyColors[1].getAvgContrast()) /
                        allSkyColors[1].getStdDevContrast(),
                    (allSkyColors[2].calcContrastToOther(rV, gV, bV)
                        - allSkyColors[2].getAvgContrast()) /
                        allSkyColors[2].getStdDevContrast()
                    )
                );
            }
        }
        
    }

    private GroupPixelColors[] partitionInto3ByColorDifference(Set<PairInt> skyPoints, 
        Image originalColorImage, int xOffset, int yOffset, int[] outputCounts) {
         
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
    
    private GroupPixelColors[] partitionInto3ByBrightness(
        Set<PairInt> skyPoints, 
        Image originalColorImage, int xOffset, int yOffset, HistogramHolder[] h) {
         
        int i = 0;
        float[] brightness = new float[skyPoints.size()];
        PairInt[] ps = new PairInt[brightness.length];
        for (PairInt p : skyPoints) {
            int x = p.getX() + xOffset;
            int y = p.getY() + yOffset;
            int r = originalColorImage.getR(x, y);
            int g = originalColorImage.getG(x, y);
            int b = originalColorImage.getB(x, y);
            
            float[] hsb = new float[3];
            Color.RGBtoHSB(r, g, b, hsb);
        
            brightness[i] = hsb[2];
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
    
    private void findEmbeddedClouds(Set<PairInt> embeddedPoints, 
        Set<PairInt> skyPoints, Set<PairInt> excludeThesePoints,
        Image originalColorImage, GreyscaleImage mask,
        Map<Integer, PixelColors> pixelColorsMap,
        Map<PairInt, Set<PixelColors> > skyColorsMap,
        HistogramHolder[] brightnessHistogram,
        GroupPixelColors[] skyPartitionedByBrightness) {
        
        findClouds2(embeddedPoints, skyPoints, excludeThesePoints,
            originalColorImage, mask,
            pixelColorsMap, skyColorsMap,
            brightnessHistogram, skyPartitionedByBrightness, true, true);
        
debugPlot(skyPoints, originalColorImage, 
    mask.getXRelativeOffset(), mask.getYRelativeOffset(), 
    "sky_added_embedded_clouds");

    }
    
    private void findSeparatedClouds(Set<PairInt> sunPoints, 
        Set<PairInt> skyPoints, Set<PairInt> excludeThesePoints,
        Image originalColorImage, GreyscaleImage mask,
        Map<Integer, PixelColors> pixelColorsMap,
        Map<PairInt, Set<PixelColors> > skyColorsMap,
        HistogramHolder[] brightnessHistogram,
        GroupPixelColors[] skyPartitionedByBrightness) {
        
        findClouds2(sunPoints, skyPoints, excludeThesePoints,
            originalColorImage, mask,
            pixelColorsMap, skyColorsMap,
            brightnessHistogram, skyPartitionedByBrightness, false, false);
    }

    /**
     * Find clouds using given data and one of 2 methods:
     * (1) if "findOnlyConnected" is true, pointsToConnectOrExclude should hold 
     * only the points one wants to find if they pass color and contrast filters
     * and from there, any adjacent neighbors passing the filters will be
     * found and also added to the skyPoints;
     * or (2) if "findOnlyConnected" is false, pointsToConnectOrExclude should
     * hold any points that should be excluded if found as a candidate 
     * neighbor while the search is over each pixel in originalColorImage.
     * 
     * Note that in case (2), skyPoints and pointsToConnectOrExclude might
     * be the same point set.
     * @param pointsToConnectOrExclude
     * @param skyPoints
     * @param excludeThesePoints
     * @param originalColorImage
     * @param mask
     * @param brightnessHistogram
     * @param skyPartitionedByBrightness
     * @param pointsAreEmbeddedInSky 
     * @param pointIsEmbeddedInSky if true, the pixel filter will add a pass
     * for bluish or pinkish white clouds, else false.  it should be safe to
     * use this as true when searching for clouds embedded already in sky that
     * could not be confused with snowy mountain peaks or foreground.
     */
    private void findClouds2(Set<PairInt> pointsToConnectOrExclude, 
        final Set<PairInt> skyPoints, Set<PairInt> excludeThesePoints,
        Image originalColorImage, GreyscaleImage mask,
        Map<Integer, PixelColors> pixelColorsMap,
        Map<PairInt, Set<PixelColors> > skyColorsMap,
        HistogramHolder[] brightnessHistogram,
        GroupPixelColors[] skyPartitionedByBrightness,
        boolean findOnlyConnected, boolean pointsAreEmbeddedInSky) {
  
        if (pointsToConnectOrExclude == null) {
            if (!findOnlyConnected) {
                pointsToConnectOrExclude = new HashSet<PairInt>();
            } else {
                throw new IllegalArgumentException(
                    "if findOnlyConnected is true, then " + 
                    "pointsToConnectOrExclude has to contain points that " +
                    "hold the points to connect to");
            }
        }
        
        int xOffset = mask.getXRelativeOffset();
        int yOffset = mask.getYRelativeOffset();
     
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
     
        // find all image pixels that are not in skyPoints or sunPoints
        // that are within the color range of brightest and next brightest sky
        
        Set<PairInt> cloudPoints = new HashSet<PairInt>();
        
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
       
        Stack<PairInt> stack = new Stack<PairInt>();
        
        if (!findOnlyConnected) {

            findSeparatedClouds(pointsToConnectOrExclude, skyPoints,
                excludeThesePoints, originalColorImage, mask,
                brightnessHistogram, skyPartitionedByBrightness,
                skyIsRed, skyIsPurple, hasDarkGreyClouds, useKEqualsZero,
                pointsAreEmbeddedInSky, cloudPoints, stack);
            
            int n = cloudPoints.size();

            log.info("NUMBER of separated cloud points = " + n);

            DFSConnectedGroupsFinder groupsFinder = new DFSConnectedGroupsFinder();
            groupsFinder.setMinimumNumberInCluster(100);
            groupsFinder.findConnectedPointGroups(skyPoints, mask.getWidth(), 
                mask.getHeight());

            log.info("NUMBER of groups of connected skyPoints=" 
                + groupsFinder.getNumberOfGroups());

            if (groupsFinder.getNumberOfGroups() > 2) {
                                
                //only grow contiguous groups of cloudPoints.
                groupsFinder = new DFSConnectedGroupsFinder();
                groupsFinder.setMinimumNumberInCluster(30);
                groupsFinder.findConnectedPointGroups(cloudPoints, 
                    mask.getWidth(), mask.getHeight());
                
                cloudPoints.clear();
                for (int ii = 0; ii < groupsFinder.getNumberOfGroups(); ii++) {
                    Set<PairInt> g = groupsFinder.getXY(ii);
                    cloudPoints.addAll(g);
                }
                
                int nBefore = cloudPoints.size();
                                
                findClouds(cloudPoints, skyPoints, originalColorImage, mask, 
                    pixelColorsMap, skyColorsMap);

                int nAfter = cloudPoints.size();

                log.info("after 'add separated clouds', findClouds finds " 
                    + nBefore + ":" + nAfter);
                
                skyPoints.addAll(cloudPoints);

                for (PairInt p : cloudPoints) {
                    int x = p.getX();
                    int y = p.getY();            
                    mask.setValue(x, y, 0);
                }

debugPlot(skyPoints, originalColorImage, mask.getXRelativeOffset(), 
    mask.getYRelativeOffset(), "after_find_clouds_in_findClouds2");

                if ((((nAfter - nBefore)/nBefore) > 2) || 
                    ((nAfter - nBefore) > 1000)) {
 
                    // one more round of search for embedded but only in the new points
                    Set<PairInt> embeddedPoints = findEmbeddedNonPoints(
                        cloudPoints);

                    skyPoints.addAll(embeddedPoints);

                    for (PairInt p : embeddedPoints) {
                        int x = p.getX();
                        int y = p.getY();
                        mask.setValue(x, y, 0);
                    }

debugPlot(skyPoints, originalColorImage, mask.getXRelativeOffset(), 
    mask.getYRelativeOffset(), "after_findEmbeddedNonPoints");
int z = 1;

                    //TODO: cases where sky is not contiguous at this point?
                    embeddedPoints = findEmbeddedNonPointsExtendedToImageBounds(
                        skyPoints, mask, 1);

                    skyPoints.addAll(embeddedPoints);

                    for (PairInt p : embeddedPoints) {
                        int x = p.getX();
                        int y = p.getY();
                        mask.setValue(x, y, 0);
                    }
                }
                
            } else {
                
                Set<PairInt> embeddedPoints = findEmbeddedContiguousUnincludedPoints(
                    skyPoints, mask);

                skyPoints.addAll(embeddedPoints);

                for (PairInt p : embeddedPoints) {
                    int x = p.getX();
                    int y = p.getY();
                    mask.setValue(x, y, 0);
                }
            }

        } else {
        
            // == case (2)  find only points connected to 
            // pointsToConnectOrExclude and add them to skyPoints,
            // afterwards find non-sky points embedded within skyPoints
            
            // find the points in and connected to pointsToConnectOrExclude
            Set<PairInt> excludePoints = new HashSet<PairInt>();
            
            for (PairInt uPoint : pointsToConnectOrExclude) {
                
                // populate stack and cloud points with the initial points
                // that pass the color and contrast filter:
                 
                if (excludePoints.contains(uPoint)) {
                    continue;
                }
                
                filterToAddCloudPixel(excludePoints, skyPoints, 
                    originalColorImage, xOffset, yOffset, uPoint,
                    brightnessHistogram, skyPartitionedByBrightness,
                    skyIsRed, skyIsPurple, hasDarkGreyClouds, useKEqualsZero,
                    pointsAreEmbeddedInSky,
                    cloudPoints, stack
                );
            }
            
            if (stack.isEmpty()) {
                return;
            }
            
            Set<PairInt> visited = new HashSet<PairInt>();
            visited.add(stack.peek());
            
            // find the neighbors of the points in the stack
            
            while (!stack.isEmpty()) {
                
                PairInt uPoint = stack.pop();
                
                int uX = uPoint.getX();
                int uY = uPoint.getY();

                for (int vX = (uX - 1); vX <= (uX + 1); vX++) {

                    for (int vY = (uY - 1); vY <= (uY + 1); vY++) {
                    
                        PairInt vPoint = new PairInt(vX, vY);
                        
                        if (vPoint.equals(uPoint) || 
                            excludePoints.contains(vPoint) ||
                            visited.contains(vPoint)) {
                            continue;
                        }

                        filterToAddCloudPixel(excludePoints, skyPoints, 
                            originalColorImage, xOffset, yOffset, vPoint,
                            brightnessHistogram, skyPartitionedByBrightness,
                            skyIsRed, skyIsPurple, hasDarkGreyClouds, 
                            useKEqualsZero, pointsAreEmbeddedInSky,
                            cloudPoints, stack
                        );

                        visited.add(vPoint);
                    }
                }
            }
                    
            skyPoints.addAll(cloudPoints);

            for (PairInt p : cloudPoints) {
                int x = p.getX();
                int y = p.getY();            
                mask.setValue(x, y, 0);
            }
            
        } // end case (2)
        
    }

    private GroupPixelColors calculateColorOfOuterEdgeOfSun(
        Set<PairInt> sunPoints, Image originalColorImage, int xOffset, int yOffset) {
       
        EllipseHelper ellipseHelper = new EllipseHelper();
        
        double[] params = ellipseHelper.fitEllipseToPoints(sunPoints);
        
        if (params == null) {
            // not close to an ellipse
            return null;
        }
        
        float xc = (float)params[0];
        float yc = (float)params[1];
        float a = (float)params[2];
        float b = (float)params[3];
        float alpha = (float)params[4];
        
        int xMin = Integer.MAX_VALUE;
        int xMax = Integer.MIN_VALUE;
        int yMin = Integer.MAX_VALUE;
        int yMax = Integer.MIN_VALUE;
       
        for (PairInt p : sunPoints) {
            int x = p.getX();
            int y = p.getY();
            if (x < xMin) {
                xMin = x;
            }
            if (y < yMin) {
                yMin = y;
            }
            if (x > xMax) {
                xMax = x;
            }
            if (y > yMax) {
                yMax = y;
            }
        }
        
        double radius = (alpha < 0.5) ? 0.5*(xMax - xMin) : Math.sqrt(a*a - b*b);
        
        double limit = radius - 1;
        
        Set<PairInt> outer = new HashSet<PairInt>();
        for (PairInt p : sunPoints) {
            int x = p.getX() + xOffset;
            int y = p.getY() + yOffset;
            float dx = xc - x;
            float dy = yc - y;
                
            double dist = Math.sqrt(dx*dx + dy*dy);
            
            if (dist >= limit) {
                outer.add(p);
            }
        }
        
        return new GroupPixelColors(outer, originalColorImage, xOffset, yOffset);
    }

    private int count(List<PairIntArray> zeroPointLists) {
        
        int n = 0;
        for (PairIntArray p : zeroPointLists) {
            n += p.getN();
        }
        
        return n;
    }

    private void debugPlot(Set<PairInt> extSkyPoints, Image originalColorImage, 
        int xOffset, int yOffset, String outputPrefixForFileName) {
        
        //plot is made in aspects
        
    }
    
    private void debugPlot(Set<PairInt> r0, Set<PairInt> r1, Set<PairInt> r2, 
        Set<PairInt> r3, Image originalColorImage, int xOffset, int yOffset, 
        String filtered_out_of_clouds) {
        
        //plot is made in aspects
    }

    private void filterToAddCloudPixel(Set<PairInt> excludeFromThesePoints, 
        Set<PairInt> skyPoints, Image origColorImg, 
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
        
        // because of Rayleigh scattered illumination, need large tolerance
        int limit = 55;
        boolean isGrey = (Math.abs(r - g) < limit) && (Math.abs(g - b) < limit) 
            && (Math.abs(r - b) < 1.8*limit);

        float[] hsb = new float[3];
        Color.RGBtoHSB(r, g, b, hsb);

        for (int k = 2; k >= 0; k--) {                    

            if ((k == 0) && !useKEqualsZero) {
                continue;
            } 

            double contrast = skyBinsByBrightness[k].calcContrastToOther(r, g, b);

            double diffContrast = contrast - 
                skyBinsByBrightness[k].getAvgContrast();
            
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

if (((p.getX() + xOffset) >= 518) && ((p.getX() + xOffset) <= 518) && ((p.getY() + yOffset) >= 245) && ((p.getY() + yOffset) <= 245)) {
CIEChromaticity cieC = new CIEChromaticity();
float[] xychr = cieC.rgbToXYChromaticity(r, g, b);  //0.3136, 0.33;     0.3136, 0.33
float[] xychrsky = cieC.rgbToXYChromaticity(        //0.3039, 0.3136;   0.3130, 0.3218
    Math.round(skyBinsByBrightness[k].getAvgRed()),
    Math.round(skyBinsByBrightness[k].getAvgGreen()),
    Math.round(skyBinsByBrightness[k].getAvgBlue()));

boolean t1 = (pointIsEmbeddedInSky);
boolean t2 = (diffR <= skyBinsByBrightness[k].getStdDevRed());
boolean t3 = (skyIsRed && (r > skyBinsByBrightness[k].getAvgRed()));
boolean t4 = (skyIsRed && (r > 155) && (diffR <= 3.5*skyBinsByBrightness[k].getStdDevRed()));
boolean t5 = (!skyIsRed  && (diffR <= 2.0*skyBinsByBrightness[k].getStdDevRed()));
                 
boolean t6 = (diffG <= 1.5*skyBinsByBrightness[k].getStdDevGreen());
boolean t7 = (skyIsRed && (contrast < 0.) && (g > 130) && (diffG <= 2.0*skyBinsByBrightness[k].getStdDevGreen()));
boolean t8 = ((skyIsPurple || skyIsRed) && (contrast > 0.) && (g < 130) && 
                        (diffG <= 2.0*skyBinsByBrightness[k].getStdDevGreen()));
boolean t9 = (!skyIsRed && hasDarkGreyClouds && 
                        (diffG <= 2.0*skyBinsByBrightness[k].getStdDevGreen()));
boolean t10 = (diffB <= 1.2*skyBinsByBrightness[k].getStdDevBlue());
                    
boolean t11 = (skyIsPurple && (contrast > 0.) && (b < 130) && 
                        (diffB <= 2.5*skyBinsByBrightness[k].getStdDevBlue()));
boolean t12 = (!skyIsRed && 
                        (diffB <= 1.5*skyBinsByBrightness[k].getStdDevBlue()));
log.info(String.format("(%d,%d) t1=%b, t2=%b, t3=%b, t4=%b, t5=%b", col, row, t1, t2, t3, t4, t5));
log.info(String.format("t6=%b, t7=%b, t8=%b, t9=%b, t10=%b", t6, t7, t8, t9, t10));
log.info(String.format("t11=%b, t12=%b", t11, t12));

log.info(String.format("r,g,b=(%d,%d,%d) R/G=%f R/B=%f, G/B=%f  contrast=%f  diffRGB=(%f,%f,%f) stdDevRGB=(%f,%f,%f)  skyIsRed=%b hasDarkGreyClouds=%b", 
r, g, b,
(float)r/(float)g, (float)r/(float)b, (float)g/(float)b,
contrast,
diffR, diffG, diffB, 
skyBinsByBrightness[k].getStdDevRed(),
skyBinsByBrightness[k].getStdDevGreen(),
skyBinsByBrightness[k].getStdDevBlue(),
skyIsRed, hasDarkGreyClouds
));

}

            if (
                (pointIsEmbeddedInSky /*&& isGrey*/) ||
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

    private void populatePixelColorMaps(Set<PairInt> points, 
        Image originalColorImage, GreyscaleImage mask, 
        Map<Integer, PixelColors> pixelColorsMap, 
        Map<PairInt, Set<PixelColors>> skyColorsMap) {
        
        int xOffset = mask.getXRelativeOffset(); 
        int yOffset = mask.getYRelativeOffset();
        
        for (PairInt p : points) {
            
            int index = mask.getIndex(p.getX(), p.getY());
            
            populatePixelColorMaps(p, index, 
                originalColorImage, xOffset, yOffset, 
                pixelColorsMap, skyColorsMap);
        }
    }

    /**
     * populate pixelColorsMap and skyColorsMap by adding the colors from
     * pixelIndex to skyColorsMap for key p while also adding the colors
     * from pixelIndex to pixelColorsMap.
     * @param p
     * @param pixelIndex
     * @param originalColorImage
     * @param xOffset
     * @param yOffset
     * @param pixelColorsMap
     * @param skyColorsMap 
     */
    private void populatePixelColorMaps(PairInt p, int pixelIndex,
        Image originalColorImage, int xOffset, int yOffset, 
        Map<Integer, PixelColors> pixelColorsMap, 
        Map<PairInt, Set<PixelColors>> skyColorsMap) {
        
        int x = p.getX();
        int y = p.getY();
        
        Set<PixelColors> skyColors = skyColorsMap.get(p);
        if (skyColors == null) {
            skyColors = new HashSet<PixelColors>();
            skyColorsMap.put(p, skyColors);
        }
        PixelColors skyPC = pixelColorsMap.get(pixelIndex);
        if (skyPC == null) {
            int r = originalColorImage.getR(x + xOffset, y + yOffset);
            int g = originalColorImage.getG(x + xOffset, y + yOffset);
            int b = originalColorImage.getB(x + xOffset, y + yOffset);
            skyPC = new PixelColors(r, g, b);
            pixelColorsMap.put(pixelIndex, skyPC);
        }
        skyColors.add(skyPC);
    }

    private Set<PairInt> removeReflectedSun(List<PairIntArray> zeroPointLists, 
        Image colorImg, GreyscaleImage thetaImg) {
        
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

    private void mergeOverlappingPreviouslyRemovedNonCloudColors(
        Set<PairInt> points, List<PairIntArray> removedNonCloudColors) {
               
        log.info("before merge: " + points.size() + " points");
        
        for (PairIntArray pai : removedNonCloudColors) {
            int nOverlapping = 0;
            for (int i = 0; i < pai.getN(); i++) {
                int x = pai.getX(i);
                int y = pai.getY(i);
                PairInt p = new PairInt(x, y);
                if (points.contains(p)) {
                    nOverlapping++;
                }
            }
            //TODO: this needs more tests to find a limit empirically:
            float f = (float)nOverlapping/(float)pai.getN();
            
            log.fine("f=" + f + " nOverlapping=" + nOverlapping + " nPoints=" + pai.getN());

            if ((f > 0.001) && (nOverlapping < pai.getN()) && (nOverlapping > 3)) {
                for (int i = 0; i < pai.getN(); i++) {
                    PairInt p = new PairInt(pai.getX(i), pai.getY(i));
                    points.add(p);
                }
            }
        }
        
        log.info("after merge: " + points.size() + " points");
    }

    private int[] determineChangeTowardsSkyline(Set<PairInt> points, 
        int imgWidth, int imgHeight) {
        
        int[] skyMinMaxXY = MiscMath.findMinMaxXY(points);
        int dSkylineX = 0;
        int dSkylineY = 0;
        if ((skyMinMaxXY[0] == 0) && 
            (skyMinMaxXY[1] == (imgWidth - 1)) &&
            (skyMinMaxXY[2] == 0) && 
            (skyMinMaxXY[3] < (imgHeight - 1))
        ) {
            // this is most common. skyline is horizontal and below
            // sky in image
            dSkylineY = 1;
            
        } else if ((skyMinMaxXY[0] > 0) && 
            (skyMinMaxXY[1] == (imgWidth - 1)) &&
            (skyMinMaxXY[2] == 0) && 
            (skyMinMaxXY[3] == (imgHeight - 1))
        ) {
            // upside down image
            dSkylineY = -1;
          
        } else if ((skyMinMaxXY[0] == 0) && 
            (skyMinMaxXY[1] < (imgWidth - 1)) &&
            (skyMinMaxXY[2] == 0) && 
            (skyMinMaxXY[3] == (imgHeight - 1))
        ) {
            // landscape, w/ sky leftward
            dSkylineX = 1;
              
        } else if ((skyMinMaxXY[0] > 0) && 
            (skyMinMaxXY[1] == (imgWidth - 1)) &&
            (skyMinMaxXY[2] == 0) && 
            (skyMinMaxXY[3] == (imgHeight - 1))
        ) {
            // landscape, w/ sky rightward
            dSkylineX = -1;
            
        } else if ((skyMinMaxXY[0] > 0) && 
            (skyMinMaxXY[1] == (imgWidth - 1)) &&
            (skyMinMaxXY[2] == 0) && 
            (skyMinMaxXY[3] < (imgHeight - 1))
        ) {
            // skyline is lower left to main sky
            dSkylineX = -1;
            dSkylineY = 1;
            
        } else if ((skyMinMaxXY[0] > 0) && 
            (skyMinMaxXY[1] == (imgWidth - 1)) &&
            (skyMinMaxXY[2] > 0) && 
            (skyMinMaxXY[3] == (imgHeight - 1))
        ) {
            // skyline is upper left to main sky
            dSkylineX = -1;
            dSkylineY = -1;
            
        } else if ((skyMinMaxXY[0] == 0) && 
            (skyMinMaxXY[1] < (imgWidth - 1)) &&
            (skyMinMaxXY[2] > 0) && 
            (skyMinMaxXY[3] == (imgHeight - 1))
        ) {
            // skyline is upper right to main sky
            dSkylineX = 1;
            dSkylineY = -1;
            
        } else if ((skyMinMaxXY[0] == 0) && 
            (skyMinMaxXY[1] < (imgWidth - 1)) &&
            (skyMinMaxXY[2] == 0) && 
            (skyMinMaxXY[3] < (imgHeight - 1))
        ) {
            // skyline is lower right to main sky
            dSkylineX = 1;
            dSkylineY = 1;
            
        }
        
        return new int[]{dSkylineX, dSkylineY};
    }

    private void growForLowContrastLimits(Set<PairInt> points, 
        Set<PairInt> excludePoints, Image colorImg, GreyscaleImage mask) {
        
        // it tries to avoid adding snowy mountain tops to hazy sky pixels,
        // for example.
                
        int cWidth = colorImg.getWidth();
        int cHeight = colorImg.getHeight();
        int mWidth = mask.getWidth();
        int mHeight = mask.getHeight();
        int xOffset = mask.getXRelativeOffset();
        int yOffset = mask.getYRelativeOffset();
        
        GroupPixelColors allSkyColor = new GroupPixelColors(points,
            colorImg, xOffset, yOffset);
        
        PixelColors foregroundColor = getAveragedColorOfDifference(points,
            colorImg, xOffset, yOffset);
        
        double rDivB = allSkyColor.getAvgRed()/allSkyColor.getAvgBlue();
        
        boolean skyIsRed = (rDivB > 1);
        
        CIEChromaticity cieC = new CIEChromaticity();
        
        int[] dxs = new int[]{-1, -1,  0,  1, 1, 1, 0, -1};
        int[] dys = new int[]{ 0, -1, -1, -1, 0, 1, 1,  1};
       
        int nAdded = 0;
        int nIter = 0;
        int nMaxIter = 10;
        
        while ((nIter == 0) || ((nIter < nMaxIter) && (nAdded > 0))) {

            nAdded = 0;
            
            Map<PairInt, List<Double>> localSkyCIEXY = new 
                HashMap<PairInt, List<Double>>();

            Map<PairInt, List<Double>> localSkyRGB = new HashMap<PairInt,
                List<Double>>();

            Set<PairInt> borderPixels = new HashSet<PairInt>();

            for (PairInt uPoint : points) {

                int uX = uPoint.getX() + xOffset;
                int uY = uPoint.getY() + yOffset;

                for (int k = 0; k < dxs.length; k++) {

                    int vX = uX + dxs[k];
                    int vY = uY + dys[k];

                    if ((vX < 0) || (vX > (cWidth - 1)) || (vY < 0) || 
                        (vY > (cHeight - 1))) {
                        continue;
                    }

                    PairInt vPoint = new PairInt(vX - xOffset, vY - yOffset);

                    if (uPoint.equals(vPoint) || borderPixels.contains(uPoint) ||
                        points.contains(vPoint) || excludePoints.contains(vPoint)) {
                        continue;
                    }
                    
                    // these vPoint pixels are not currently in points

                    borderPixels.add(uPoint);

if ((vX >= 308) && (vX <= 310) && (vY >= 181) && (vY <= 182)) {
int z  = 1;
}

                    if (!localSkyCIEXY.containsKey(uPoint)) {

                        Set<PairInt> neighbors = getThe8NeighborPixelsWithin(
                            uPoint, points, mWidth, mHeight);
                        neighbors.add(uPoint);         
                        
                        long rSum = 0;
                        long gSum = 0;
                        long bSum = 0;

                        int[] r = new int[neighbors.size()];
                        int[] g = new int[neighbors.size()];
                        int[] b = new int[neighbors.size()];
                        int i = 0;
                        for (PairInt p : neighbors) {
                            int x = p.getX() + xOffset;
                            int y = p.getY() + yOffset;
                            r[i] = colorImg.getR(x, y);
                            g[i] = colorImg.getG(x, y);
                            b[i] = colorImg.getB(x, y);
                            rSum += r[i];
                            gSum += g[i];
                            bSum += b[i];
                            i++;
                        }
                        double avgR = (double)rSum/(double)r.length;
                        double avgG = (double)gSum/(double)r.length;
                        double avgB = (double)bSum/(double)r.length;

                        rSum = 0;
                        gSum = 0;
                        bSum = 0;
                        for (i = 0; i < r.length; i++) {
                            double rDiff = r[i] - avgR;
                            double gDiff = g[i] - avgG;
                            double bDiff = b[i] - avgB;
                            rSum += (rDiff * rDiff);
                            gSum += (gDiff * gDiff);
                            bSum += (bDiff * bDiff);
                            i++;
                        }
                        double stDevR = Math.sqrt(rSum/((double)r.length - 1.0));
                        double stDevG = Math.sqrt(gSum/((double)r.length - 1.0));
                        double stDevB = Math.sqrt(bSum/((double)r.length - 1.0));

                        List<Double> list = cieC.calcAvgAndStdDevXY(r, g, b);
                        localSkyCIEXY.put(uPoint, list);

                        list = new ArrayList<Double>();
                        list.add(Double.valueOf(avgR));
                        list.add(Double.valueOf(avgG));
                        list.add(Double.valueOf(avgB));
                        list.add(Double.valueOf(stDevR));
                        list.add(Double.valueOf(stDevG));
                        list.add(Double.valueOf(stDevB));

                        localSkyRGB.put(uPoint, list);
                    }        
                }
            }
             
            nAdded = growForLowContrastLimits(points, 
                excludePoints, colorImg, mask, skyIsRed, borderPixels, 
                localSkyRGB, foregroundColor);
        
            nIter++;
        }
      
    }
        
    private int growForLowContrastLimits(Set<PairInt> points, 
        Set<PairInt> excludePoints, Image colorImg, GreyscaleImage mask,
        boolean skyIsRed, Set<PairInt> borderPixels,
        Map<PairInt, List<Double>> localSkyRGB, PixelColors foregroundColor) {
        
        if (borderPixels.isEmpty()) {
            return 0;
        }
        
        Set<PixelColors> clrs = new HashSet<PixelColors>();
        Iterator<Map.Entry<PairInt, List<Double>>> iter = localSkyRGB.entrySet().iterator();
        while (iter.hasNext()) {
            Map.Entry<PairInt, List<Double>> entry = iter.next();
            List<Double> rgbAvgAndStDev = entry.getValue();
            PixelColors pc = new PixelColors(
                rgbAvgAndStDev.get(0).intValue(), 
                rgbAvgAndStDev.get(1).intValue(), 
                rgbAvgAndStDev.get(2).intValue());
            clrs.add(pc);
        }
        
        GroupPixelColors gpc = new GroupPixelColors(clrs);
        
        Stack<PairInt> stack = new Stack<PairInt>();
        stack.addAll(borderPixels);
        
        Set<PairInt> visited = new HashSet<PairInt>();
        visited.add(stack.peek());
        
        int xOffset = mask.getXRelativeOffset();
        int yOffset = mask.getYRelativeOffset();
        int width = mask.getWidth();
        int height = mask.getHeight();
        
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
        
        Stack<PairInt> addedSky = new Stack<PairInt>();
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

if ((vX >= 308) && (vX <= 310) && (vY >= 181) && (vY <= 182)) {
int z  = 1;
}
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
                
if ((vX >= 308) && (vX <= 310) && (vY >= 181) && (vY <= 182)) {
int z  = 1;
}                
               
            }
        }
        
        if (addedSky.isEmpty()) {
            return nAdded;
        }
        
        // search the neighbors of the just added sky points and add them if
        // similar in color.  
        // TODO: this doesn't use a contrast check so needs more testing.
        
        Set<PixelColors> tmpColors = new HashSet<PixelColors>(addedSkyColorToAdd);
        
        visited.clear();
        visited.add(addedSky.peek());
        
        while (!addedSky.isEmpty()) {

            PairInt uPoint = addedSky.pop();

            int uX = uPoint.getX() + xOffset;
            int uY = uPoint.getY() + yOffset;

            PixelColors uColor = addedSkyColorToAdd.pop();
            
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
                
                if ((Math.abs(r - uColor.getRed()) < firstSubtrTol) 
                    && (Math.abs(g - uColor.getGreen()) < firstSubtrTol) &&
                    (Math.abs(b - uColor.getBlue()) < firstSubtrTol)) {
                        
                    nAdded++;
                        
                    points.add(vPoint);
                        
                    mask.setValue(vX - xOffset, vY - yOffset, 0);
                    
                    addedSky.add(vPoint);
                    addedSkyColorToAdd.add(uColor);
                    
                } else {
                    
                    for (PixelColors c : tmpColors) {
                        
                        if ((Math.abs(r - c.getRed()) < secondSubtrTol)
                            && (Math.abs(g - c.getGreen()) < secondSubtrTol)
                            && (Math.abs(b - c.getBlue()) < secondSubtrTol)) {

                            nAdded++;

                            points.add(vPoint);

                            mask.setValue(vX - xOffset, vY - yOffset, 0);

                            addedSky.add(vPoint);
                            addedSkyColorToAdd.add(c);
                            
                            break;
                        }
                    }
                }                
            }
        }
        
        log.info("nAdded=" + nAdded);
        
        return nAdded;
    }

    private void adjustPerimeterToIncludeImageBoundaries(int[] rowMinMax, 
        Map<Integer, PairInt> rowColRange, int[] dSkylineXY, int imageWidth,
        int imageHeight) {
         
        if ((dSkylineXY[0] == 0) && (dSkylineXY[1] == 1)) {
            
            // this is most common. skyline is horizontal and below sky in image
            
            // set all colRange.x to 0 and all colRange.y to imageWidth-1
            // add rows from min y to 0 to the map   
            Iterator<Map.Entry<Integer, PairInt> > iter = rowColRange.entrySet().iterator();
            while (iter.hasNext()) {
                Map.Entry<Integer, PairInt> entry = iter.next();
                entry.getValue().setX(0);
                entry.getValue().setY(imageWidth - 1);
            }
            
            for (int row = 0; row < rowMinMax[1]; row++) {
                PairInt colRange = new PairInt(0, imageWidth - 1);
                rowColRange.put(Integer.valueOf(row), colRange);
            }
            
        } else if ((dSkylineXY[0] == 0) && (dSkylineXY[1] == -1)) {
          
            // upside down image
            
            // set all colRange.x to 0 and all colRange.y to imageWidth-1
            // add rows from max y to imageHeight - 1 to the map 
        
            // set all colRange.x to 0 and all colRange.y to imageWidth-1
            // add rows from max y to imageHeight to the map   
            Iterator<Map.Entry<Integer, PairInt> > iter = rowColRange.entrySet().iterator();
            while (iter.hasNext()) {
                Map.Entry<Integer, PairInt> entry = iter.next();
                entry.getValue().setX(0);
                entry.getValue().setY(imageWidth - 1);
            }
            
            for (int row = (rowMinMax[0] + 1); row < imageHeight; row++) {
                PairInt colRange = new PairInt(0, imageWidth - 1);
                rowColRange.put(Integer.valueOf(row), colRange);
            }
            
        } else if ((dSkylineXY[0] == 1) && (dSkylineXY[1] == -1)) {
            
            // main sky is lower left
            Iterator<Map.Entry<Integer, PairInt> > iter = rowColRange.entrySet().iterator();
            while (iter.hasNext()) {
                Map.Entry<Integer, PairInt> entry = iter.next();
                entry.getValue().setX(0);
            }
            
            PairInt maxRowColRange = rowColRange.get(Integer.valueOf(rowMinMax[1]));
            for (int row = (rowMinMax[1] + 1); row < imageHeight; row++) {
                PairInt colRange = new PairInt(maxRowColRange.getX(), 
                    maxRowColRange.getY());
                rowColRange.put(Integer.valueOf(row), colRange);
            }
            
        } else if ((dSkylineXY[0] == 1) && (dSkylineXY[1] == 0)) {
          
            // landscape, w/ sky leftward
            
            // set all colRange.x to 0
            // add rows from min y to 0 to the map
            // add rows from max y to imageHeight to the map
            Iterator<Map.Entry<Integer, PairInt>> iter = rowColRange.entrySet().iterator();
            while (iter.hasNext()) {
                Map.Entry<Integer, PairInt> entry = iter.next();
                entry.getValue().setX(0);
            }
            
            PairInt minRowColRange = rowColRange.get(Integer.valueOf(rowMinMax[0]));
            for (int row = 0; row < rowMinMax[1]; row++) {
                PairInt colRange = new PairInt(minRowColRange.getX(), 
                    minRowColRange.getY());
                rowColRange.put(Integer.valueOf(row), colRange);
            }
            for (int row = (rowMinMax[0] + 1); row < imageHeight; row++) {
                PairInt colRange = new PairInt(minRowColRange.getX(), 
                    minRowColRange.getY());
                rowColRange.put(Integer.valueOf(row), colRange);
            }
            
        } else if ((dSkylineXY[0] == 1) && (dSkylineXY[1] == 1)) {
          
            // skyline is upper right to main sky
            // (main sky is lower left)
            
            Iterator<Map.Entry<Integer, PairInt> > iter = rowColRange.entrySet().iterator();
            while (iter.hasNext()) {
                Map.Entry<Integer, PairInt> entry = iter.next();
                entry.getValue().setX(0);
            }
            PairInt maxRowColRange = rowColRange.get(Integer.valueOf(rowMinMax[1]));
            for (int row = (rowMinMax[1] + 1); row < imageHeight; row++) {
                PairInt colRange = new PairInt(maxRowColRange.getX(), 
                    maxRowColRange.getY());
                rowColRange.put(Integer.valueOf(row), colRange);
            }    
            
        } else if ((dSkylineXY[0] == -1) && (dSkylineXY[1] == -1)) {
            
            // skyline is upper left to main sky
            // (main sky is lower right)
            
            Iterator<Map.Entry<Integer, PairInt> > iter = rowColRange.entrySet().iterator();
            while (iter.hasNext()) {
                Map.Entry<Integer, PairInt> entry = iter.next();
                entry.getValue().setY(imageWidth - 1);
            }
            
            PairInt maxRowColRange = rowColRange.get(Integer.valueOf(rowMinMax[1]));
            for (int row = (rowMinMax[1] + 1); row < imageHeight; row++) {
                PairInt colRange = new PairInt(maxRowColRange.getX(), 
                    maxRowColRange.getY());
                rowColRange.put(Integer.valueOf(row), colRange);
            }
                   
        } else if ((dSkylineXY[0] == -1) && (dSkylineXY[1] == 0)) {
            
            // landscape, w/ sky rightward
            
            Iterator<Map.Entry<Integer, PairInt> > iter = rowColRange.entrySet().iterator();
            while (iter.hasNext()) {
                Map.Entry<Integer, PairInt> entry = iter.next();
                entry.getValue().setY(imageWidth - 1);
            }
            
            PairInt minRowColRange = rowColRange.get(Integer.valueOf(rowMinMax[0]));
            for (int row = 0; row < rowMinMax[1]; row++) {
                PairInt colRange = new PairInt(minRowColRange.getX(), 
                    minRowColRange.getY());
                rowColRange.put(Integer.valueOf(row), colRange);
            }
            for (int row = (rowMinMax[0] - 1); row < imageHeight; row++) {
                PairInt colRange = new PairInt(minRowColRange.getX(), 
                    minRowColRange.getY());
                rowColRange.put(Integer.valueOf(row), colRange);
            }
       
        } else if ((dSkylineXY[0] == -1) && (dSkylineXY[1] == 1)) {
            
            // main sky is upper right
            
            Iterator<Map.Entry<Integer, PairInt> > iter = rowColRange.entrySet().iterator();
            while (iter.hasNext()) {
                Map.Entry<Integer, PairInt> entry = iter.next();
                entry.getValue().setY(imageWidth - 1);
            }
            
            PairInt minRowColRange = rowColRange.get(Integer.valueOf(rowMinMax[0]));
            for (int row = 0; row < rowMinMax[0]; row++) {
                PairInt colRange = new PairInt(minRowColRange.getX(), 
                    minRowColRange.getY());
                rowColRange.put(Integer.valueOf(row), colRange);
            }
            
        } else if ((dSkylineXY[0] == 0) && (dSkylineXY[1] == 0)) {
            
            // the image may be tilted and sky might be filling
            // most of the image
            
            /*
            sky              sky
            sky              sky
            sky             ....   
            sky            .   .
            ....          .    .
            .   .        .     .
            .....sky    ........
            
            fill in bounds for sky to sky along boundaries of the
            corners that are determined to be sky.
            
            need to avoid doing that for the foreground, non-sky features.
            
            -- will average the color from corner to corner, and if it's
            cloud like rather than foreground like, will add those
            to boundaries?
            */
            
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
    public Set<PairInt> findSunPoints(Image colorImg, int xOffset, int yOffset,
        boolean skyIsDarkGrey) {
        
        Set<PairInt> set = new HashSet<PairInt>();
        
        SunColors sunColors = new SunColors();
        
        if (skyIsDarkGrey) {
            sunColors.useDarkSkiesLogic();
        }
        
        for (int col = 0; col < colorImg.getWidth(); col++) {
            for (int row = 0; row < colorImg.getHeight(); row++) {
                
                int r = colorImg.getR(col, row);
                int g = colorImg.getG(col, row);
                int b = colorImg.getB(col, row);

                if (sunColors.isSunCenterColor(r, g, b)) {
                    set.add(new PairInt(col - xOffset, row - yOffset));
                }
            }
        }
        
debugPlot(set, colorImg, xOffset, yOffset, 
"sunpoints");

        return set;
    }

    Set<PairInt> findRainbowColoredPoints(Image colorImg, 
        Set<PairInt> reflectedSunRemoved,
        int xOffset, int yOffset, boolean skyIsDarkGrey) {
        
        AbstractSkyRainbowColors colors = skyIsDarkGrey ?
            new DarkSkyRainbowColors() : new BrightSkyRainbowColors();
        
        Set<PairInt> set = new HashSet<PairInt>();
        
        for (int col = 0; col < colorImg.getWidth(); col++) {
            for (int row = 0; row < colorImg.getHeight(); row++) {
                
                if (reflectedSunRemoved.contains(new PairInt(col - xOffset,
                    row - yOffset))) {
                    continue;
                }
                
                int r = colorImg.getR(col, row);
                int g = colorImg.getG(col, row);
                int b = colorImg.getB(col, row);
                
                if (colors.isInRedThroughPurplishRed(r, g, b)) {
                    
                    set.add(new PairInt(col - xOffset, row - yOffset));

                } else if (colors.isInOrangeRed(r, g, b)) {
                    
                    set.add(new PairInt(col - xOffset, row - yOffset));
                
                } else if (colors.isInGreenishYellowOrange(r, g, b)) {
                 
                    set.add(new PairInt(col - xOffset, row - yOffset));

                }
            }
        }
        
        return set;
    }
    
    Set<PairInt> findRainbowColoredPointsConnectedToSkyPoints(
        Set<PairInt> skyPoints, Set<PairInt> reflectedSunRemoved,
        Image colorImg, int xOffset, int yOffset, boolean skyIsDarkGrey) {
        
        Set<PairInt> rainbowPoints = new HashSet<PairInt>();
        
        AbstractSkyRainbowColors rainbowColors = skyIsDarkGrey ?
            new DarkSkyRainbowColors() : new BrightSkyRainbowColors();
        
        int width = colorImg.getWidth();
        int height = colorImg.getHeight();
        
        int[] dxs = new int[]{-1, -1,  0,  1, 1, 1, 0, -1};
        int[] dys = new int[]{ 0, -1, -1, -1, 0, 1, 1,  1};
        
        java.util.Stack<PairInt> stack = new java.util.Stack<PairInt>();

        for (PairInt p : skyPoints) {
            
            int x = p.getX() - xOffset;
            int y = p.getY() - yOffset;
            
            int r = colorImg.getR(x, y);
            int g = colorImg.getG(x, y);
            int b = colorImg.getB(x, y);

            if (rainbowColors.isInGreenishYellowOrange(r, g, b)
                || rainbowColors.isInOrangeRed(r, g, b)
                || rainbowColors.isInRedThroughPurplishRed(r, g, b)) {
                
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

                int r = colorImg.getR(vX, vY);
                int g = colorImg.getG(vX, vY);
                int b = colorImg.getB(vX, vY);
                        
                if (rainbowColors.isInGreenishYellowOrange(r, g, b) ||
                    rainbowColors.isInOrangeRed(r, g, b) ||
                    rainbowColors.isInRedThroughPurplishRed(r, g, b)) {
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

    private boolean skyIsDarkGrey(GroupPixelColors allSkyColor) {
        
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
                
                if (uPoint.equals(vPoint) || !points.contains(vPoint)) {
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
    private void plotSkyColor(Set<PairInt> points, Image colorImg, 
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

    private void findSeparatedClouds(Set<PairInt> pointsToConnectOrExclude, 
        Set<PairInt> skyPoints, Set<PairInt> excludeThesePoints,
        Image originalColorImage, GreyscaleImage mask,
        HistogramHolder[] brightnessHistogram, 
        GroupPixelColors[] skyPartitionedByBrightness, 
        boolean skyIsRed, boolean skyIsPurple, boolean hasDarkGreyClouds, 
        boolean useKEqualsZero, boolean pointsAreEmbeddedInSky, 
        Set<PairInt> cloudPoints, Stack<PairInt> stack) {
        
        int xOffset = mask.getXRelativeOffset();
        int yOffset = mask.getYRelativeOffset();
        
        for (int col = 0; col < originalColorImage.getWidth(); col++) {
            for (int row = 0; row < originalColorImage.getHeight(); row++) {

                PairInt p = new PairInt(col - xOffset, row - yOffset);
    
                if (excludeThesePoints.contains(p)) {
                    continue;
                }
             
                filterToAddCloudPixel(pointsToConnectOrExclude, skyPoints, 
                    originalColorImage, xOffset, yOffset, p,
                    brightnessHistogram, skyPartitionedByBrightness,
                    skyIsRed, skyIsPurple, hasDarkGreyClouds, useKEqualsZero,
                    pointsAreEmbeddedInSky,
                    cloudPoints, stack
                );
            }
        }

        if (stack.isEmpty()) {
            return;
        }

        // this block does not use the stack further because it has already
        // searched each pixel of the image

        // find the number of groups of connected pixels within skyPoints and
        // extSkyPoints

        for (PairInt p : cloudPoints) {
            int x = p.getX();
            int y = p.getY();            
            mask.setValue(x, y, 0);
        }

        skyPoints.addAll(cloudPoints);
    }

    private void reduceToLargestContiguousGroup(Set<PairInt> points,
        GreyscaleImage mask, int minNumberPointsInGroup,
        Set<PairInt> placeholdingPoints) {
        
        int width = mask.getWidth();
        int height = mask.getHeight();
        
        Set<PairInt> temp = new HashSet<PairInt>(points);
        temp.addAll(placeholdingPoints);
        
        // only keep the contiguous sky
        DFSConnectedGroupsFinder groupsFinder = new DFSConnectedGroupsFinder();
        groupsFinder.setMinimumNumberInCluster(minNumberPointsInGroup);
        groupsFinder.findConnectedPointGroups(temp, width, height);

        int nMax = Integer.MIN_VALUE;
        int maxIdx = -1;
        for (int ii = 0; ii < groupsFinder.getNumberOfGroups(); ii++) {
            int n = groupsFinder.getNumberofGroupMembers(ii);
            if (n > nMax) {
                nMax = n;
                maxIdx = ii;
            }
        }
        if (maxIdx == -1) {
            throw new IllegalStateException("error while finding"
                + " contiguous sky points");
        }

        mask.fill(1);
        
        Set<PairInt> g = groupsFinder.getXY(maxIdx);
        points.clear();
        
        for (PairInt p : g) {
            if (!placeholdingPoints.contains(p)) {
                
                points.add(p);
                
                int x = p.getX();
                int y = p.getY();
                mask.setValue(x, y, 0);
            }
        }

    }

    private List<PairIntArray> removeSetsWithNonCloudColors(List<PairIntArray> 
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
        
        List<PairIntArray> removed = new ArrayList<PairIntArray>();
        
        if (!remove.isEmpty()) {
            
            for (int i = (remove.size() - 1); i > -1; i--) {
                
                PairIntArray r = zeroPointLists.remove(remove.get(i).intValue());
                
                removed.add(r);
            }
        }
        
        return removed;
    }

    private Set<PairInt> findEmbeddedContiguousUnincludedPoints(
        Set<PairInt> points, GreyscaleImage mask) {
        
        Set<PairInt> embeddedPoints = new HashSet<PairInt>();
        
        // adjust perimeter to enclose the edges that are touching the image
        // boundaries
        PerimeterFinder finder = new PerimeterFinder();
        int[] rowMinMax = new int[2];
        Map<Integer, PairInt> rowColRange = finder.find(points, rowMinMax);
       
        DFSContiguousValueFinder contiguousNonZeroFinder = new DFSContiguousValueFinder(mask);
        contiguousNonZeroFinder.setMinimumNumberInCluster(1);
        contiguousNonZeroFinder.findEmbeddedGroupsNotThisValue(0,
            rowColRange, rowMinMax);

        for (int i = 0; i < contiguousNonZeroFinder.getNumberOfGroups(); i++) {
            contiguousNonZeroFinder.getXY(i, embeddedPoints);
        }
        
        return embeddedPoints;
    }
    
    private Set<PairInt> findEmbeddedNonPoints(Set<PairInt> points) {
        
        Set<PairInt> embeddedPoints = new HashSet<PairInt>();
        
        PerimeterFinder finder = new PerimeterFinder();
        int[] rowMinMax = new int[2];
        Map<Integer, PairInt> rowColRange = finder.find(points, rowMinMax);

        for (int row = rowMinMax[0]; row <= rowMinMax[1]; row++) {            
            
            PairInt colRange = rowColRange.get(Integer.valueOf(row));
            
            for (int col = colRange.getX(); col <= colRange.getY(); col++) {
                
                PairInt p = new PairInt(col, row);
                
                if (!points.contains(p)) {
                    embeddedPoints.add(p);
                }
            }
        }
        
        return embeddedPoints;
    }
    
    private Set<PairInt> findEmbeddedNonPointsExtendedToImageBounds(
        Set<PairInt> points, GreyscaleImage mask, int minGroupSize) {
        
        // determine the direction towards the skyline from the main sky
        int[] dSkylineXY = determineChangeTowardsSkyline(points, 
            mask.getWidth(), mask.getHeight());

        // adjust perimeter to enclose the edges that are touching the image
        // boundaries
        PerimeterFinder finder = new PerimeterFinder();
        int[] rowMinMax = new int[2];
        Map<Integer, PairInt> rowColRange = finder.find(points, rowMinMax);
        adjustPerimeterToIncludeImageBoundaries(rowMinMax, rowColRange, 
            dSkylineXY, mask.getWidth(), mask.getHeight());

        DFSContiguousValueFinder contiguousNonZeroFinder = new DFSContiguousValueFinder(mask);
        contiguousNonZeroFinder.setMinimumNumberInCluster(minGroupSize);
        contiguousNonZeroFinder.findEmbeddedGroupsNotThisValue(0,
            rowColRange, rowMinMax);

        Set<PairInt> embeddedPoints = new HashSet<PairInt>();
        for (int i = 0; i < contiguousNonZeroFinder.getNumberOfGroups(); i++) {
            contiguousNonZeroFinder.getXY(i, embeddedPoints);
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

}
