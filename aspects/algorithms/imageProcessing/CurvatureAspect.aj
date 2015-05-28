package algorithms.imageProcessing;

import algorithms.imageProcessing.util.*;
import java.io.*;
import java.util.logging.Logger;
import algorithms.util.*;
import algorithms.misc.*;
import java.awt.Color;
import java.util.*;
import java.util.Map.Entry;
import org.ejml.ops.*;

/**
class to create debugging output from the main source code
*/
public aspect CurvatureAspect {

    private PolygonAndPointPlotter plotter = null;

    // for the ContourFinder plots:
    private PolygonAndPointPlotter plotter9 = null;
    private PolygonAndPointPlotter plotter10 = null;

    private Logger log2 = Logger.getAnonymousLogger();

    private String currentEdgeStr = "";

    private static int outImgNum = -1;

    after() :
        target(algorithms.imageProcessing.ContourFinder) 
        && execution(public ContourFinder.new(..)) {

        log2.fine("===> in aspect for ContourFinder");

        if (plotter9 == null) {
            try {
            plotter9 = new PolygonAndPointPlotter();
            plotter10 = new PolygonAndPointPlotter();
            } catch (IOException e) {
                log2.severe("ERROR: " + e.getMessage());
            }
        }
    }

    before(ScaleSpaceCurveImage scaleSpaceImage, int edgeNumber) :
        execution(List<CurvatureScaleSpaceContour> algorithms.imageProcessing.ContourFinder.findContours(ScaleSpaceCurveImage, int))
        && args(scaleSpaceImage, edgeNumber) 
	    && target(algorithms.imageProcessing.ContourFinder) {

        Object obj = thisJoinPoint.getThis();

        if (!(obj instanceof ContourFinder)) {
            return;
        }

        ContourFinder instance = (ContourFinder)obj;

        Object[] args = (Object[])thisJoinPoint.getArgs();

        ScaleSpaceCurveImage scaleSpaceImageIn = (ScaleSpaceCurveImage)args[0];
        int edgeNumberIn = ((Integer)args[1]).intValue();

        String filePath = printImage(scaleSpaceImageIn, plotter9, 9);

        log2.fine("===> in aspect for ContourFinder, filePath=" + filePath);
    }

    before() :
        execution(protected void algorithms.imageProcessing.AbstractCurvatureScaleSpaceMapper*.extractSkyline())
        && args() 
	    && target(algorithms.imageProcessing.AbstractCurvatureScaleSpaceMapper) {

        Object obj = thisJoinPoint.getThis();

        if (!(obj instanceof AbstractCurvatureScaleSpaceMapper)) {
            return;
        }

        AbstractCurvatureScaleSpaceMapper instance = 
            (AbstractCurvatureScaleSpaceMapper)obj;

        GreyscaleImage gXY2 = instance.getGradientXY().copyImage();
        PairIntArray gXYValues = Histogram.createADescendingSortByKeyArray(gXY2);

//TODO: for images where significant number of pixels were "put back" in the
// sky to correct things like repetitive structure at same scale as convolution diffs,
// the critical fraction should probably be 0.5
// it might be that very blue skies need f>= 0.5, else 0.15?
// need a test image of blue skies with as many clouds as the az or nm images
        //find where the contour value stays above 0.3 of the count for min value
        float v0 = gXYValues.getY(gXYValues.getN() - 1);
        int min = gXYValues.getX(gXYValues.getN() - 1);
        for (int i = (gXYValues.getN() - 1); i > -1; i--) {
            float f = (float)gXYValues.getY(i)/v0;
            if (f >= 0.5) {
                min = gXYValues.getX(i);
            } else {
                break;
            }
        }
        
        StringBuilder sb = new StringBuilder("==>");
        sb.append("(").append(Integer.valueOf(gXYValues.getN())).append(")");
        for (int i = (gXYValues.getN() - 1); i > -1; i--) {
            int pixelValue = gXYValues.getX(i);
            int count = gXYValues.getY(i);
            sb.append(" (").append(pixelValue).append(",")
                .append(count).append(",")
                .append((float)count/v0).append(")");
        }
        log2.info(sb.toString());

        for (int i = 0; i < gXY2.getNPixels(); i++) {
            int v = gXY2.getValue(i);
            v -= min;
            if (v < 0) {
                v = 0;
            }
            gXY2.setValue(i, v);
        }

        try {
            String dirPath = ResourceFinder.findDirectory("bin");
            ImageIOHelper.writeOutputImage(
                dirPath + "/gXY_thresh_" + outImgNum + ".png", gXY2);

        } catch (Exception e) {
            log2.severe("ERROR: " + e.getMessage());
        }
    }

    after(float[] rainbowCoeff, 
        Set<PairInt> rainbowPoints, ImageExt originalColorImage, int xOffset, 
        int yOffset) returning(RainbowFinder.Hull rainbowHull) :
        execution(private RainbowFinder.Hull RainbowFinder.createRainbowHull(float[],
        Set<PairInt>, ImageExt, int, int))
        && args(rainbowCoeff, rainbowPoints, originalColorImage, xOffset, 
        yOffset)
	    && target(algorithms.imageProcessing.RainbowFinder) {

        Object obj = thisJoinPoint.getThis();

        if (!(obj instanceof RainbowFinder)) {
            return;
        }

        Set<PairInt> hullPoints = ((RainbowFinder)obj).getPointsToExcludeInHull();

        Image clr = originalColorImage.copyImage();

        try {
            String dirPath = ResourceFinder.findDirectory("bin");

            ImageIOHelper.addToImage(hullPoints, xOffset, yOffset, clr);

            ImageIOHelper.writeOutputImage(
                dirPath + "/rainbow_hull_" + outImgNum + ".png", 
                clr);

        } catch (IOException e) {
            log2.severe("ERROR: " + e.getMessage());
        }
 
        try {
            PolygonAndPointPlotter plotter = new PolygonAndPointPlotter(0, 
                clr.getWidth(), 0, clr.getHeight());
            plotter.addPlot(rainbowHull.xHull, rainbowHull.yHull, 
                null, null, "rainbow hull " + outImgNum);
            
            String fileName = plotter.writeFile(Integer.valueOf(9182));
                        
        } catch (IOException e) {
            Logger.getLogger(this.getClass().getName()).severe(e.getMessage());
        }

    }

    after(ImageExt colorImg, int xOffset, int yOffset, boolean skyIsDarkGrey) :
        execution(public void findSunPhotosphere(ImageExt, int, int, boolean))
        && args(colorImg, xOffset, yOffset, skyIsDarkGrey)
	    && target(algorithms.imageProcessing.SunFinder) {

        Object obj = thisJoinPoint.getThis();

        if (!(obj instanceof SunFinder)) {
            return;
        }

        Set<PairInt> points = ((SunFinder)obj).getSunPoints();

        Image clr = colorImg.copyImage();

        try {
            String dirPath = ResourceFinder.findDirectory("bin");

            ImageIOHelper.addToImage(points, xOffset, yOffset, clr);

            ImageIOHelper.writeOutputImage(
                dirPath + "/sun_" + outImgNum + ".png", clr);

        } catch (IOException e) {
            log2.severe("ERROR: " + e.getMessage());
        }

    }

    after(Set<PairInt> extSkyPoints, Image originalColorImage, 
        int xOffset, int yOffset, String outputPrefixForFileName) returning() :
        execution(void debugPlot(Set<PairInt>, Image,
        int, int, String))
        && args(extSkyPoints, originalColorImage, xOffset, yOffset, outputPrefixForFileName) {

        Image clr = originalColorImage.copyImage();

        try {
            String dirPath = ResourceFinder.findDirectory("bin");

            ImageIOHelper.addToImage(extSkyPoints, xOffset, yOffset, clr);

            ImageIOHelper.writeOutputImage(
                dirPath + "/" + outputPrefixForFileName + "_" + outImgNum + ".png", clr);

        } catch (IOException e) {
            log2.severe("ERROR: " + e.getMessage());
        }
    }

    after(Set<PairInt> r0, Set<PairInt> r1, Set<PairInt> r2, Set<PairInt> r3, 
        Image originalColorImage, 
        int xOffset, int yOffset, String outputPrefixForFileName) returning() :
        execution(void SkylineExtractor.debugPlot(Set<PairInt>,
        Set<PairInt>, Set<PairInt>, Set<PairInt>, Image, int, int, String))
        && args(r0, r1, r2, r3, originalColorImage, xOffset, yOffset, 
        outputPrefixForFileName)
	    && target(algorithms.imageProcessing.SkylineExtractor) {

        Image clr = originalColorImage.copyImage();

        try {
            String dirPath = ResourceFinder.findDirectory("bin");

            List<Set<PairInt>> sets = new ArrayList<Set<PairInt>>();
            sets.add(r0);
            sets.add(r1);
            sets.add(r2);
            sets.add(r3);
            ImageIOHelper.addAlternatingColorPointSetsToImage(sets,
                xOffset, yOffset, clr);

            ImageIOHelper.writeOutputImage(
                dirPath + "/" + outputPrefixForFileName + "_" + outImgNum + ".png", clr);

        } catch (IOException e) {
            log2.severe("ERROR: " + e.getMessage());
        }
    }

    after(Set<PairInt> points, Set<PairInt> excludeThesePoints, 
        GreyscaleImage gradientXY) returning() :
        execution(private void SkylineExtractor.growZeroValuePoints(Set<PairInt>, 
        Set<PairInt>, GreyscaleImage))
        && args(points, excludeThesePoints, gradientXY)
	    && target(algorithms.imageProcessing.SkylineExtractor) {

        Object[] args = (Object[])thisJoinPoint.getArgs();
        GreyscaleImage gXY = (GreyscaleImage)args[2];

        GreyscaleImage mask = gXY.createWithDimensions();
        mask.fill(250);

        for (PairInt p : points) {
            int x = p.getX();
            int y = p.getY();
            mask.setValue(x, y, 0);
        }

        try {
            String dirPath = ResourceFinder.findDirectory("bin");
            ImageIOHelper.writeOutputImage(
                dirPath + "/sky_gxy_" + outImgNum + ".png", mask);
        } catch (IOException e) {
            log2.severe("ERROR: " + e.getMessage());
        }

    }

    before(Set<PairInt> skyPoints, Set<PairInt> reflectedSunRemoved, 
        ImageExt clrImage, int xOffset, int yOffset, boolean skyIsDarkGrey,
        Set<PairInt> sunPoints) :
        call(protected double[] SkylineExtractor.findSunConnectedToSkyPoints(
            Set<PairInt>, Set<PairInt>, ImageExt, int, int, boolean,  Set<PairInt>))
        && args(skyPoints, reflectedSunRemoved, clrImage, xOffset, yOffset,
            skyIsDarkGrey, sunPoints)
	    && target(algorithms.imageProcessing.SkylineExtractor) {

        ImageExt clr = (ImageExt)clrImage.copyImage();

        try {
            String dirPath = ResourceFinder.findDirectory("bin");

            ImageIOHelper.addToImage(skyPoints, xOffset, yOffset, clr);

            ImageIOHelper.writeOutputImage(
                dirPath + "/sky_points_before_sun_search_" + outImgNum + ".png", clr);

        } catch (IOException e) {
            log2.severe("ERROR: " + e.getMessage());
        }
    }

    after(Set<PairInt> skyPoints, Set<PairInt> reflectedSunRemoved, 
        ImageExt clrImage, int xOffset, int yOffset, boolean skyIsDarkGrey,
        Set<PairInt> sunPoints) 
        returning(double[] params) :
        execution(protected double[] SkylineExtractor.findSunConnectedToSkyPoints(
            Set<PairInt>, Set<PairInt>, ImageExt, int, int, boolean, Set<PairInt>))
        && args(skyPoints, reflectedSunRemoved, clrImage, xOffset, yOffset,
        skyIsDarkGrey, sunPoints)
	    && target(algorithms.imageProcessing.SkylineExtractor) {

        ImageExt clr = (ImageExt)clrImage.copyImage();

        log2.info("plotting " + sunPoints.size() + " sun points");

        try {
            String dirPath = ResourceFinder.findDirectory("bin");

            ImageIOHelper.addToImage(sunPoints, xOffset, yOffset, clr);

            ImageIOHelper.writeOutputImage(
                dirPath + "/sky_sun_points_" + outImgNum + ".png", clr);

        } catch (IOException e) {
            log2.severe("ERROR: " + e.getMessage());
        }
    }

    before(Set<PairInt> skyPoints, Set<PairInt> excludePoints, 
        ImageExt originalColorImage, GreyscaleImage mask
        ) 
        : call(void SkylineExtractor.findClouds(Set<PairInt>, Set<PairInt>,
        ImageExt, GreyscaleImage)) 
        && args(skyPoints, excludePoints, originalColorImage, mask) 
        && target(algorithms.imageProcessing.SkylineExtractor) {

        ImageExt clr = (ImageExt)originalColorImage.copyImage();

        int xOffset = mask.getXRelativeOffset();
        int yOffset = mask.getYRelativeOffset();

        try {
            String dirPath = ResourceFinder.findDirectory("bin");

            ImageIOHelper.addToImage(skyPoints, xOffset, yOffset, clr);

            ImageIOHelper.writeOutputImage(
                dirPath + "/sky_before_find_clouds_" + outImgNum + ".png", clr);

        } catch (IOException e) {
            log2.severe("ERROR: " + e.getMessage());
        }
    }

    after(Set<PairInt> skyPoints, Set<PairInt> excludePoints, 
        ImageExt originalColorImage, GreyscaleImage mask
        )
        returning() :
        execution(void SkylineExtractor.findClouds(
            Set<PairInt>, Set<PairInt>, ImageExt, GreyscaleImage
        ) )
        && args(skyPoints, excludePoints, originalColorImage, mask)
	    && target(algorithms.imageProcessing.SkylineExtractor) {

        ImageExt clr = (ImageExt)originalColorImage.copyImage();

        int xOffset = mask.getXRelativeOffset();
        int yOffset = mask.getYRelativeOffset();

        try {
            String dirPath = ResourceFinder.findDirectory("bin");

            ImageIOHelper.addToImage(skyPoints, xOffset, yOffset, clr);

            ImageIOHelper.writeOutputImage(
                dirPath + "/sky_after_find_clouds_" + outImgNum + 
                "_" + n3 + ".png", clr);

            ImageIOHelper.writeOutputImage(
                dirPath + "/mask_after_find_clouds_" + outImgNum + 
                "_" + n3 + ".png", mask);

            n3++;

        } catch (IOException e) {
            log2.severe("ERROR: " + e.getMessage());
        }
    }

private static int n2 = 0;
private static int n3 = 0;

    after(Set<PairInt> skyPoints, Set<PairInt> reflectedSunRemoved, 
        ImageExt colorImg, int xOffset, int yOffset,
        int imageWidth, int imageHeight,
        boolean skyIsDarkGrey, SkylineExtractor.RemovedSets removedSets)
        returning() :
        execution(public void RainbowFinder.findRainbowInImage(
        Set<PairInt>, Set<PairInt>, ImageExt, int, int,
        int, int, boolean, SkylineExtractor.RemovedSets) )
        && args(skyPoints, reflectedSunRemoved, colorImg, xOffset, yOffset, 
        imageWidth, imageHeight, skyIsDarkGrey, removedSets)
	    && target(algorithms.imageProcessing.RainbowFinder) {

        Object obj = thisJoinPoint.getThis();

        if (!(obj instanceof RainbowFinder)) {
            return;
        }

        Set<PairInt> rainbowPoints = ((RainbowFinder)obj).getRainbowPoints();

        if (rainbowPoints.isEmpty()) {
            return;
        }

        ImageExt clr = (ImageExt)colorImg.copyImage();

        try {
            String dirPath = ResourceFinder.findDirectory("bin");

            ImageIOHelper.addToImage(rainbowPoints, xOffset, yOffset, clr);

            ImageIOHelper.writeOutputImage(
                dirPath + "/rainbow_" + outImgNum + ".png", clr);

            float[] cieX = new float[rainbowPoints.size()];
            float[] cieY = new float[rainbowPoints.size()];
            CIEChromaticity cieC = new CIEChromaticity();
            int count = 0;
            for (PairInt p : rainbowPoints) {
                int x = p.getX();
                int y = p.getY();
                int r = colorImg.getR(x + xOffset, y + yOffset);
                int g = colorImg.getG(x + xOffset, y + yOffset);
                int b = colorImg.getB(x + xOffset, y + yOffset);
                float[] cie = cieC.rgbToXYChromaticity(r, g, b);
                cieX[count] = cie[0];
                cieY[count] = cie[1];
                count++;
            }
            float cieXMin = MiscMath.findMin(cieX);
            float cieXMax = MiscMath.findMax(cieX);
            float cieYMin = MiscMath.findMin(cieY);
            float cieYMax = MiscMath.findMax(cieY);
        
            PolygonAndPointPlotter plotter = new PolygonAndPointPlotter(
                0.f, 0.8f, 0.0f, 0.9f);
            plotter.addPlot(cieX, cieY, cieX, cieY, "CIE for rainbow points");

            plotter.writeFile(outImgNum*10);

        } catch (IOException e) {
            log2.severe("ERROR: " + e.getMessage());
        }
    }

    before(List<PairIntArray> zeroPointLists, Image originalColorImage, 
        GreyscaleImage theta, boolean addAlongX, int addAmount) :
        execution(private Set<PairInt> SkylineExtractor.removeSetsWithNonCloudColors(
            List<PairIntArray>, Image, GreyscaleImage, boolean, int) )
        && args(zeroPointLists, originalColorImage, theta, addAlongX, addAmount)
	    && target(algorithms.imageProcessing.SkylineExtractor) {

        Image clr = originalColorImage.copyImage();

        int xOffset = theta.getXRelativeOffset();
        int yOffset = theta.getYRelativeOffset();

        try {
            String dirPath = ResourceFinder.findDirectory("bin");

            ImageIOHelper.addAlternatingColorCurvesToImage(zeroPointLists, clr);

            ImageIOHelper.writeOutputImage(
                dirPath + "/sky_before_noncloudcolors_removed_" + outImgNum + ".png", clr);

        } catch (IOException e) {
            log2.severe("ERROR: " + e.getMessage());
        }
    }

    after(List<PairIntArray> zeroPointLists, ImageExt colorImg, 
        GreyscaleImage thetaImg)
        returning(Set<PairInt> reflectedSun) :
        execution(private Set<PairInt> SkylineExtractor.removeReflectedSun(
            List<PairIntArray>, ImageExt, GreyscaleImage) )
        && args(zeroPointLists, colorImg, thetaImg)
	    && target(algorithms.imageProcessing.SkylineExtractor) {

        ImageExt clr = (ImageExt)colorImg.copyImage();

        int xOffset = thetaImg.getXRelativeOffset();
        int yOffset = thetaImg.getYRelativeOffset();

        try {
            String dirPath = ResourceFinder.findDirectory("bin");

            ImageIOHelper.addAlternatingColorCurvesToImage(zeroPointLists, clr);

            ImageIOHelper.writeOutputImage(
                dirPath + "/sky_reflectedsun_removed_" + outImgNum + ".png", clr);

        } catch (IOException e) {
            log2.severe("ERROR: " + e.getMessage());
        }
    }

    after(List<PairIntArray> zeroPointLists, Image originalColorImage, 
        GreyscaleImage theta, boolean addAlongX, int addAmount)
        returning(Set<PairInt> removedNonCloudColors) :
        execution(private Set<PairInt> SkylineExtractor.removeSetsWithNonCloudColors(
            List<PairIntArray>, Image, GreyscaleImage, boolean, int) )
        && args(zeroPointLists, originalColorImage, theta, addAlongX, addAmount)
	    && target(algorithms.imageProcessing.SkylineExtractor) {

        Image clr = originalColorImage.copyImage();

        int xOffset = theta.getXRelativeOffset();
        int yOffset = theta.getYRelativeOffset();

        try {
            String dirPath = ResourceFinder.findDirectory("bin");

            ImageIOHelper.addAlternatingColorCurvesToImage(zeroPointLists, clr);

            ImageIOHelper.writeOutputImage(
                dirPath + "/sky_after_noncloudcolors_removed_" + outImgNum + ".png", clr);

        } catch (IOException e) {
            log2.severe("ERROR: " + e.getMessage());
        }
    }

    before(List<PairIntArray> zeroPointLists, Image originalColorImage, 
        GreyscaleImage theta, double avgY,
        Set<PairInt> outputRemovedPoints) :
        execution(private void SkylineExtractor.removeHighContrastPoints(
            List<PairIntArray>, Image, GreyscaleImage, double,
            Set<PairInt>) )
        && args(zeroPointLists, originalColorImage, theta, avgY, 
        outputRemovedPoints)
	    && target(algorithms.imageProcessing.SkylineExtractor) {

        Image clr = originalColorImage.copyImage();

        int xOffset = theta.getXRelativeOffset();
        int yOffset = theta.getYRelativeOffset();

        try {
            String dirPath = ResourceFinder.findDirectory("bin");

            ImageIOHelper.addAlternatingColorCurvesToImage(zeroPointLists, clr);

            ImageIOHelper.writeOutputImage(
                dirPath + "/sky_before_high_contrast_points_removed_" + outImgNum + ".png", clr);

        } catch (IOException e) {
            log2.severe("ERROR: " + e.getMessage());
        }
    }

    after(List<PairIntArray> zeroPointLists, Image originalColorImage, 
        GreyscaleImage theta, double avgY, Set<PairInt> outputRemovedPoints)
        returning() :
        execution(private void SkylineExtractor.removeHighContrastPoints(
            List<PairIntArray>, Image, GreyscaleImage, double,
            Set<PairInt>) )
        && args(zeroPointLists, originalColorImage, theta, avgY, 
        outputRemovedPoints)
	    && target(algorithms.imageProcessing.SkylineExtractor) {

        Image clr = originalColorImage.copyImage();

        int xOffset = theta.getXRelativeOffset();
        int yOffset = theta.getYRelativeOffset();

        try {
            String dirPath = ResourceFinder.findDirectory("bin");

            ImageIOHelper.addAlternatingColorCurvesToImage(zeroPointLists, clr);

            ImageIOHelper.writeOutputImage(
                dirPath + "/sky_high_contrast_points_removed_" + outImgNum + ".png", clr);

        } catch (IOException e) {
            log2.severe("ERROR: " + e.getMessage());
        }
    }

    after(GreyscaleImage gradientXY, Set<PairInt> skyPoints, 
        Set<PairInt> excludeThesePoints)
        returning() :
        execution(public GreyscaleImage SkylineExtractor.extractSkyFromGradientXY(
            GreyscaleImage, Set<PairInt>, Set<PairInt>) )
        && args(gradientXY, skyPoints, excludeThesePoints)
	    && target(algorithms.imageProcessing.SkylineExtractor) {

        Image clr = gradientXY.copyImageToGreen();

        try {
            String dirPath = ResourceFinder.findDirectory("bin");

            ImageIOHelper.addToImage(skyPoints, 0, 0, clr);

            ImageIOHelper.writeOutputImage(
                dirPath + "/sky_from_gradientXY_" + outImgNum + ".png", clr);

        } catch (IOException e) {
            log2.severe("ERROR: " + e.getMessage());
        }
    }

    after(ScaleSpaceCurveImage scaleSpaceImage, int sigmaIndex, int tIndex) 
        returning() :
        execution(void algorithms.imageProcessing.ContourFinder.removeContourFromImage(ScaleSpaceCurveImage, int, int))
        && args(scaleSpaceImage, sigmaIndex, tIndex) 
	    && target(algorithms.imageProcessing.ContourFinder) {

        Object obj = thisJoinPoint.getThis();

        if (!(obj instanceof ContourFinder)) {
            return;
        }

        Object[] args = (Object[])thisJoinPoint.getArgs();

        ScaleSpaceCurveImage scaleSpaceImageIn = (ScaleSpaceCurveImage)args[0];

        String filePath = printImage(scaleSpaceImageIn, plotter10, 10);
    }

    after() returning() 
        : execution(void algorithms.imageProcessing.CurvatureScaleSpaceInflectionMapper*.createMatchedPointArraysFromContourPeaks() ) 
        && args()
        && target(algorithms.imageProcessing.CurvatureScaleSpaceInflectionMapper) {
    
        Object obj = thisJoinPoint.getThis();

        if (!(obj instanceof CurvatureScaleSpaceInflectionMapper)) {
            return;
        }

        CurvatureScaleSpaceInflectionMapper instance = (CurvatureScaleSpaceInflectionMapper)obj;

        // these are in the reference frame of the original image
        PairIntArray xyc1 = instance.getMatchedXY1().copy();
        PairIntArray xyc2 = instance.getMatchedXY2().copy();

        log2.info("n matched contour points in image1 = " + xyc1.getN());
        log2.info("n matched contour points in image2 = " + xyc2.getN());

        try {

            Image img1 = instance.getOriginalImage1().copyImage();
            
            ImageIOHelper.addCurveToImage(xyc1, img1, 2, 255, 0, 0);
            String dirPath = ResourceFinder.findDirectory("bin");
            ImageIOHelper.writeOutputImage(
                dirPath + "/matched_contour_peaks1.png", img1);

            Image img2 = instance.getOriginalImage2().copyImage();
            
            ImageIOHelper.addCurveToImage(xyc2, img2, 2, 255, 0, 0);
            ImageIOHelper.writeOutputImage(
                dirPath + "/matched_contour_peaks2.png", img2);

        } catch (IOException e) {
            log2.severe("ERROR: " + e.getMessage());
        }

    }

    after() returning() 
        : execution(void algorithms.imageProcessing.CurvatureScaleSpaceInflectionMapperForOpenCurves*.createMatchedPointArraysFromContourPeaks() ) 
        && args()
        && target(algorithms.imageProcessing.CurvatureScaleSpaceInflectionMapperForOpenCurves) {
    
        Object obj = thisJoinPoint.getThis();

        if (!(obj instanceof CurvatureScaleSpaceInflectionMapperForOpenCurves)) {
            return;
        }

        CurvatureScaleSpaceInflectionMapperForOpenCurves instance = 
            (CurvatureScaleSpaceInflectionMapperForOpenCurves)obj;

        if (instance.getMatchedXY1() == null) {
             return;
        }

        // these are in the reference frame of the original image
        PairIntArray xyc1 = instance.getMatchedXY1().copy();
        PairIntArray xyc2 = instance.getMatchedXY2().copy();

        log2.info("n matched contour points in image1 = " + xyc1.getN());
        log2.info("n matched contour points in image2 = " + xyc2.getN());

        try {

            Image img1 = instance.getOriginalImage1().copyImage();
            
            ImageIOHelper.addCurveToImage(xyc1, img1, 2, 255, 0, 0);
            String dirPath = ResourceFinder.findDirectory("bin");
            ImageIOHelper.writeOutputImage(
                dirPath + "/matched_contour_peaks1.png", img1);

            Image img2 = instance.getOriginalImage2().copyImage();
            
            ImageIOHelper.addCurveToImage(xyc2, img2, 2, 255, 0, 0);
            ImageIOHelper.writeOutputImage(
                dirPath + "/matched_contour_peaks2.png", img2);

        } catch (IOException e) {
            log2.severe("ERROR: " + e.getMessage());
        }

    }

    after() returning() 
        : execution(public void algorithms.imageProcessing.AbstractCurvatureScaleSpaceInflectionMapper*.initialize() ) 
        && args()
        && target(algorithms.imageProcessing.AbstractCurvatureScaleSpaceInflectionMapper) {
    
        Object obj = thisJoinPoint.getThis();

        if (!(obj instanceof AbstractCurvatureScaleSpaceInflectionMapper)) {
            return;
        }

        AbstractCurvatureScaleSpaceInflectionMapper instance = 
            (AbstractCurvatureScaleSpaceInflectionMapper)obj;

        List<CurvatureScaleSpaceContour> c1 = instance.getContours1();
        List<CurvatureScaleSpaceContour> c2 = instance.getContours2();
        PairIntArray xy1 = new PairIntArray();
        PairIntArray xy2 = new PairIntArray();

        for (int i = 0; i < c1.size(); i++) {
            CurvatureScaleSpaceImagePoint[] cp = c1.get(i).getPeakDetails();
            for (int j = 0; j < cp.length; j++) {
                int x = cp[j].getXCoord() + instance.getOffsetImageX1();
                int y = cp[j].getYCoord() + instance.getOffsetImageY1();
                xy1.add(x, y);
            }
        }
        for (int i = 0; i < c2.size(); i++) {
            CurvatureScaleSpaceImagePoint[] cp = c2.get(i).getPeakDetails();
            for (int j = 0; j < cp.length; j++) {
                int x = cp[j].getXCoord() + instance.getOffsetImageX2();
                int y = cp[j].getYCoord() + instance.getOffsetImageY2();
                xy2.add(x, y);
            }
        }

        try {

            Image img1 = instance.getOriginalImage1().copyImage();
            
            ImageIOHelper.addCurveToImage(xy1, img1, 2, 255, 0, 0);
            String dirPath = ResourceFinder.findDirectory("bin");
            ImageIOHelper.writeOutputImage(
                dirPath + "/contour_peaks1.png", img1);

            Image img2 = instance.getOriginalImage2().copyImage();
            
            ImageIOHelper.addCurveToImage(xy2, img2, 2, 255, 0, 0);
            ImageIOHelper.writeOutputImage(
                dirPath + "/contour_peaks2.png", img2);

        } catch (IOException e) {
            log2.severe("ERROR: " + e.getMessage());
        }
    }

    before(final PairIntArray edge, 
        final Map<SIGMA, ScaleSpaceCurve> scaleSpaceCurves, int edgeNumber) 
        : call(PairIntArray CurvatureScaleSpaceCornerDetector*.findCornersInScaleSpaceMap( 
        PairIntArray, Map<SIGMA, ScaleSpaceCurve>, int) ) 
        && args(edge, scaleSpaceCurves, edgeNumber) 
        && target(algorithms.imageProcessing.CurvatureScaleSpaceCornerDetector) {

        try {
            plotter = new PolygonAndPointPlotter();
        } catch (IOException ex) {
            ex.printStackTrace();
            throw new RuntimeException("ERROR: l64" + ex.getMessage());
        }
    }

	after(PairIntArray curve, int cIndex, float[] g, final boolean calcX) 
        returning(double curvature) :
        call(double Kernel1DHelper.convolvePointWithKernel(PairIntArray, int, 
        float[], boolean))
        && args(curve, cIndex, g, calcX) 
	    && target(algorithms.imageProcessing.Kernel1DHelper) {

        /*Object[] args = (Object[])thisJoinPoint.getArgs();

        PairIntArray c = (PairIntArray)args[0];
        int cIdx = (Integer)args[1];

		String str = String.format("k=%f  x,y=(%d, %d)", curvature,
            c.getX(cIdx), c.getY(cIdx));

        log2.info(str);*/
	}

    after(GreyscaleImage img, int col, int row, float[] g, final boolean calcX) returning(double curvature) :
        call(double Kernel1DHelper*.convolvePointWithKernel(GreyscaleImage, int, int, float[], boolean))
        && args(img, col, row, g, calcX) 
	    && target(algorithms.imageProcessing.Kernel1DHelper) {		

        /*Object[] args = (Object[])thisJoinPoint.getArgs();
        Integer c = (Integer)args[1];
        Integer r = (Integer)args[2];

		String str = String.format("k=%f  x,y=(%d, %d)", curvature, c, r);
        log2.info(str);*/
	}

    after(GreyscaleImage input, HistogramHolder imgHistogram) returning(ImageStatistics stats) :
        execution(public ImageStatistics LowIntensityRemovalFilter+.removeLowIntensityPixels(GreyscaleImage, HistogramHolder))
        && args(input, imgHistogram) 
	    && target(algorithms.imageProcessing.LowIntensityRemovalFilter) {		

        Object[] args = (Object[])thisJoinPoint.getArgs();

        GreyscaleImage img = (GreyscaleImage)args[0];

        try {
            String dirPath = ResourceFinder.findDirectory("bin");
            ImageIOHelper.writeOutputImage(dirPath + "/filters_1.png", img);
        } catch(IOException e) {
            log2.severe(e.getMessage());
        }
	}

    // note, execution binds to CannyEdgeFilter which is good,
    // but call binds to the invoker CurvatureScaleSpaceCornerDetector
    after(GreyscaleImage img) returning() :
        execution(public void CannyEdgeFilter.applyFilter(GreyscaleImage))
        && args(img)
	    && target(algorithms.imageProcessing.CannyEdgeFilter) {

        Object obj = thisJoinPoint.getThis();

        log2.fine("after CannyEdgeFilter.applyFilter instance=" 
            + obj.getClass().getName());

        if (!(obj instanceof CannyEdgeFilter)) {
            return;
        }

        CannyEdgeFilter instance = (CannyEdgeFilter)obj;

        Object[] args = (Object[])thisJoinPoint.getArgs();
        GreyscaleImage image = (GreyscaleImage)args[0];

        debugDisplay(image, "after edge thinning");

        GreyscaleImage cX = instance.getGradientX();
        GreyscaleImage cY = instance.getGradientY();
        GreyscaleImage cXY = instance.getGradientXY();
        GreyscaleImage cTheta = instance.getTheta();

        outImgNum++;
        n2 = 0;
        n3 = 0;

        try {
            String dirPath = ResourceFinder.findDirectory("bin");
            ImageIOHelper.writeOutputImage(
                dirPath + "/gX_" + outImgNum + ".png", cX);
            ImageIOHelper.writeOutputImage(
                dirPath + "/gY_" + outImgNum + ".png", cY);
            ImageIOHelper.writeOutputImage(
                dirPath + "/gXY_" + outImgNum + ".png", cXY);
            ImageIOHelper.writeOutputImage(
                dirPath + "/gTheta_" + outImgNum + ".png", cTheta);

            ImageIOHelper.writeOutputImage(dirPath + "/edgethinned.png", image);

        } catch (IOException e) {
            log2.severe(e.getMessage());
        }
   
    }

    after(GreyscaleImage input, float sigma) returning() :
        call(public void ImageProcessor.blur(GreyscaleImage, float))
        && args(input, sigma)
	    && target(algorithms.imageProcessing.ImageProcessor) {

        Object[] args = (Object[])thisJoinPoint.getArgs();
        GreyscaleImage output = (GreyscaleImage)args[0];

        debugDisplay(output, "blurred by " + sigma);

        try {
            String dirPath = ResourceFinder.findDirectory("bin");
            ImageIOHelper.writeOutputImage(dirPath + "/blurred.png", output);
        } catch (IOException e) {
            log2.severe(e.getMessage());
        }
    }

    after(GreyscaleImage input) returning() :
        execution(void CannyEdgeFilter*.apply2LayerFilter(GreyscaleImage))
        && args(input)
	    && target(algorithms.imageProcessing.CannyEdgeFilter) {

        log2.fine("after CannyEdgeFilter*.apply2LayerFilter");

        Object[] args = (Object[])thisJoinPoint.getArgs();
        GreyscaleImage output = (GreyscaleImage)args[0];

        debugDisplay(output, "after 2-threshold level filtering");

        try {
            String dirPath = ResourceFinder.findDirectory("bin");
            ImageIOHelper.writeOutputImage(dirPath + "/filters2.png", output);
        } catch (IOException e) {
            log2.severe(e.getMessage());
        }
    }

    after(final GreyscaleImage theta, 
        GreyscaleImage gradientXY, ImageExt originalColorImage, 
        CannyEdgeFilterSettings edgeSettings, PairIntArray outputSkyCentroid) 
        returning(GreyscaleImage mask) :
        execution(GreyscaleImage SkylineExtractor*.createBestSkyMask(
        GreyscaleImage, GreyscaleImage, ImageExt, CannyEdgeFilterSettings, PairIntArray))
        && args(theta, gradientXY, originalColorImage, edgeSettings, outputSkyCentroid) 
        && target(algorithms.imageProcessing.SkylineExtractor) {

        GreyscaleImage mask2 = mask.copyImage();
        MatrixUtil.multiply(mask2.getValues(), 250);

        try {
            String dirPath = ResourceFinder.findDirectory("bin");
            ImageIOHelper.writeOutputImage(
                dirPath + "/mask_" + outImgNum + ".png", mask2);
        } catch (IOException e) {
            log2.severe(e.getMessage());
        }
    }
    
    after() returning : 
	    target(algorithms.imageProcessing.CurvatureScaleSpaceCornerDetector) 
        && execution(protected void extractSkyline()) {

        Object obj = thisJoinPoint.getThis();

        if (!(obj instanceof CurvatureScaleSpaceCornerDetector)) {
            return;
        }

        CurvatureScaleSpaceCornerDetector instance = 
            (CurvatureScaleSpaceCornerDetector)obj;

        Image img3 = instance.getOriginalImage().copyImage();

        List<PairIntArray> edges = instance.getSkylineEdgesInOriginalReferenceFrame();

        if (edges.isEmpty()) {
            return; 
        }

        try {
            
            for (PairIntArray edge : edges) {
                ImageIOHelper.addCurveToImage(edge, img3, 1, 255, 255, 0);
            }
            for (PairIntArray edge : edges) {
                ImageIOHelper.addCurveToImage(edge, img3, 0, 255, 255, 255);
            }
            
            debugDisplay(edges, "skyline edges extracted", true, img3.getWidth(), 
                img3.getHeight());
            
            String dirPath = ResourceFinder.findDirectory("bin");
            ImageIOHelper.writeOutputImage(
                dirPath + "/image_with_skyline_" + outImgNum + ".png", img3);
 
        } catch (IOException ex) {
            throw new RuntimeException("ERROR: l423" + ex.getMessage());
        }
    }

    after(GreyscaleImage input, int k) returning() :
        call(public void ImageProcessor.applyImageSegmentation(GreyscaleImage, int))
        && args(input, k)
	    && target(algorithms.imageProcessing.ImageProcessor) {

        Object[] args = (Object[])thisJoinPoint.getArgs();
        GreyscaleImage image = (GreyscaleImage)args[0];
        String kStr = ((Integer)args[1]).toString();

        debugDisplay(image, "segmented by k=" + kStr);

        try {
            String dirPath = ResourceFinder.findDirectory("bin");
            ImageIOHelper.writeOutputImage(dirPath + "/segmented.png", image);
        } catch (IOException e) {
            log2.severe(e.getMessage());
        }
    }

    before() 
        : call(void CurvatureScaleSpaceCornerDetector*.initialize() ) 
        && args() 
        && target(algorithms.imageProcessing.CurvatureScaleSpaceCornerDetector) {

        Object obj = thisJoinPoint.getThis();

        if (!(obj instanceof CurvatureScaleSpaceCornerDetector)) {
            return;
        }

        CurvatureScaleSpaceCornerDetector instance = 
            (CurvatureScaleSpaceCornerDetector)obj;

        if (!instance.getInitialized()) {
            debugDisplay(instance.getImage(), "original image");
        }
    }

    after(HistogramHolder h) returning(int[] peakAndMinIndexes) :
        call(protected int[] CurvatureScaleSpaceCornerDetector.findFirstPeakAndMinimum(HistogramHolder))
        && args(h) 
	    && target(algorithms.imageProcessing.CurvatureScaleSpaceCornerDetector) { 		

        Object[] args = (Object[])thisJoinPoint.getArgs();

        HistogramHolder hist = (HistogramHolder)args[0];
        
        try {

            float xMin = MiscMath.findMax(hist.getXHist());
            float xMax = MiscMath.findMax(hist.getXHist());

            float yMax = MiscMath.findMax(hist.getYHist());

            if ((xMax > 0) && (yMax > 0)) {

                plotter.addPlot(
                    xMin, xMax, 0, yMax,
                    hist.getXHist(), hist.getYHistFloat(), 
                    hist.getXErrors(), hist.getYErrors(), 
                    currentEdgeStr);

                String filePath = plotter.writeFile(Integer.valueOf(10));
            }

            log2.fine(currentEdgeStr + ") x=" + Arrays.toString(hist.getXHist()));

            log2.fine(currentEdgeStr + ") y=" + Arrays.toString(hist.getYHist()));

        } catch(IOException e) {
           log2.severe(e.getMessage());
        }
	}

    // execution works here, but call does not
    after() returning : 
	    target(algorithms.imageProcessing.CurvatureScaleSpaceCornerDetector) 
        && execution(public void findCorners())
    {
    
        Object obj = thisJoinPoint.getThis();

        log2.info("after findCorners() " + obj.getClass().getSimpleName());

        if (!(obj instanceof CurvatureScaleSpaceCornerDetector)) {
            return;
        }

        CurvatureScaleSpaceCornerDetector instance = 
            (CurvatureScaleSpaceCornerDetector)obj;

        Image img3 = instance.getOriginalImage().copyImage();

        List<PairIntArray> edges = instance.getEdgesInOriginalReferenceFrame();

        PairIntArray corners = instance.getCornersInOriginalReferenceFrame();

        try {
            // multi-colored edges alone:
            debugDisplay(edges, "edges extracted", true, img3.getWidth(), img3.getHeight());

            
            for (PairIntArray edge : edges) {
                ImageIOHelper.addCurveToImage(edge, img3, 0, 255, 255, 0);
            }
            ImageIOHelper.addCurveToImage(corners, img3, 1, 255, 0, 0);
            
            /*ImageIOHelper.addCurveToImage(
                instance.getCornersForMatchingInOriginalReferenceFrame(),
                img3, 2, 255, 0, 255);
            */

            debugDisplay(img3, "original image w/ edges and corners overplotted");
            String dirPath = ResourceFinder.findDirectory("bin");
            ImageIOHelper.writeOutputImage(
                dirPath + "/image_with_edges_corners.png", img3);
 
        } catch (IOException ex) {
            throw new RuntimeException("ERROR: l298" + ex.getMessage());
        }
    }

    before(ScaleSpaceCurve scaleSpace, int edgeNumber, boolean rmFlseCrnrs, boolean isAClosedCurve) :
        target(algorithms.imageProcessing.CurvatureScaleSpaceCornerDetector)
        && call(protected PairFloatArray CurvatureScaleSpaceCornerDetector.findCornersInScaleSpaceMap(
        ScaleSpaceCurve, int, boolean, boolean))
        && args(scaleSpace, edgeNumber, rmFlseCrnrs, isAClosedCurve) {

        log2.fine("before findCornersInScaleSpaceMap for edge " 
            + Integer.toString(edgeNumber) + ":");

        Object[] args = (Object[])thisJoinPoint.getArgs();
        Integer eNumber = (Integer)args[1];

        currentEdgeStr = eNumber.toString();

    }

    after(PairIntArray xyCurve, List<Integer> maxCandidateCornerIndexes, 
        boolean isAClosedCurve) returning(PairIntArray jaggedLines) :
        call(protected PairIntArray CurvatureScaleSpaceCornerDetector.removeFalseCorners(
        PairIntArray, List<Integer>, boolean))
        && args(xyCurve, maxCandidateCornerIndexes, isAClosedCurve) 
        && target(algorithms.imageProcessing.CurvatureScaleSpaceCornerDetector) {

        //disable for most uses
        if (true) {
            return;
        }
        Object obj = thisJoinPoint.getThis();

        if (!(obj instanceof CurvatureScaleSpaceCornerDetector)) {
            return;
        }

        CurvatureScaleSpaceCornerDetector instance = (CurvatureScaleSpaceCornerDetector)obj;

        Object[] args = (Object[])thisJoinPoint.getArgs();
        PairIntArray xy = (PairIntArray)args[0];
        List<Integer> cI = (List<Integer>)args[1];

        Image image2 = instance.getImage().copyImageToGreen();
        
        // draw a blue line between jagged line end points
        int nExtraForDot = 2;
        int rClr = 0; int bClr = 250; int gClr = 0;
        for (int i = 0; i < jaggedLines.getN(); i++) {

            int idx0 = jaggedLines.getX(i);
            int idx1 = jaggedLines.getY(i);

            int x0 = xy.getX(idx0);
            int y0 = xy.getY(idx0);

            int x1 = xy.getX(idx1);
            int y1 = xy.getY(idx1);

            for (int dx = (-1*nExtraForDot); dx < (nExtraForDot + 1); dx++) {
                float xx = x0 + dx;
                if ((xx > -1) && (xx < (image2.getWidth() - 1))) {
                    for (int dy = (-1*nExtraForDot); dy < (nExtraForDot + 1); dy++) {
                        float yy = y0 + dy;
                        if ((yy > -1) && (yy < (image2.getHeight() - 1))) {
                            image2.setRGB((int)xx, (int)yy, 0, 250, 250);
                        }
                    }
                }
            }
            for (int dx = (-1*nExtraForDot); dx < (nExtraForDot + 1); dx++) {
                float xx = x1 + dx;
                if ((xx > -1) && (xx < (image2.getWidth() - 1))) {
                    for (int dy = (-1*nExtraForDot); dy < (nExtraForDot + 1); dy++) {
                        float yy = y1 + dy;
                        if ((yy > -1) && (yy < (image2.getHeight() - 1))) {
                            image2.setRGB((int)xx, (int)yy, 0, 0, 250);
                        }
                    }
                }
            }
        }
        nExtraForDot = 1;
        rClr = 250; bClr = 0; gClr = 0;
        for (Integer idx : cI) {

            int x = xy.getX(idx.intValue());
            int y = xy.getY(idx.intValue());

            for (int dx = (-1*nExtraForDot); dx < (nExtraForDot + 1); dx++) {
                float xx = x + dx;
                if ((xx > -1) && (xx < (image2.getWidth() - 1))) {
                    for (int dy = (-1*nExtraForDot); dy < (nExtraForDot + 1); dy++) {
                        float yy = y + dy;
                        if ((yy > -1) && (yy < (image2.getHeight() - 1))) {
                            image2.setRGB((int)xx, (int)yy, rClr, gClr, bClr);
                        }
                    }
                }
            }
        }

        try {
            ImageDisplayer.displayImage("corners and jagged lines", image2);
        } catch (IOException e) {
            log2.severe(e.getMessage());
        }
      
    }

    after(ScaleSpaceCurve scaleSpace, int edgeNumber, boolean rmFlseCrnrs,
        boolean isAClosedCurve) 
        returning(PairFloatArray xy) :
        target(algorithms.imageProcessing.CurvatureScaleSpaceCornerDetector)
        && call(protected PairFloatArray CurvatureScaleSpaceCornerDetector.findCornersInScaleSpaceMap(
        ScaleSpaceCurve, int, boolean, boolean))
        && args(scaleSpace, edgeNumber, rmFlseCrnrs, isAClosedCurve) {

        log2.fine("after findCornersInScaleSpaceMap for edge " 
            + Integer.toString(edgeNumber) + ":");

        for (int i = 0; i < xy.getN(); i++) {
            log2.fine(String.format("(%f, %f)", xy.getX(i), xy.getY(i)));
        }
    }

    before() 
        : call(List<PairIntArray> EdgeExtractor*.findEdges() ) 
        && args() 
        && target(algorithms.imageProcessing.EdgeExtractor) {

        log2.fine("before EdgeExtractor*.findEdges()");

        Object obj = thisJoinPoint.getThis();

        if (!(obj instanceof EdgeExtractor)) {
            return;
        }

        EdgeExtractor instance = (EdgeExtractor)obj;

        try {
            String dirPath = ResourceFinder.findDirectory("bin");

            ImageIOHelper.writeOutputImage(dirPath + "/edges0.png", 
                instance.getImage());

            plotHistogram(instance.getImage(), "edges before dfs");

        } catch (IOException ex) {

            log2.severe("ERROR: l359" + ex.getMessage());
        }
    }

    after() returning(List<PairIntArray> edges) : 
	    target(algorithms.imageProcessing.EdgeExtractor) 
        && execution(List<PairIntArray> EdgeExtractor*.findEdges()) {

        Object obj = thisJoinPoint.getThis();

        if (!(obj instanceof EdgeExtractor)) {
            return;
        }

        EdgeExtractor instance = (EdgeExtractor)obj;

        GreyscaleImage img3 = instance.getImage();

        try {
            String dirPath = ResourceFinder.findDirectory("bin");

            // multi-colored edges alone:
            debugDisplay(edges, "edges extracted (not in original frame)", 
                true, img3.getWidth(), img3.getHeight());

        } catch (IOException ex) {

            log2.severe("ERROR: l359" + ex.getMessage());
        }
    }

    before(ImageStatistics stats, HistogramHolder originalImageHistogram) 
        : call(int LowIntensityRemovalFilter*.determineLowThreshold( 
        ImageStatistics,  HistogramHolder))
        && args(stats, originalImageHistogram) 
        && target(algorithms.imageProcessing.LowIntensityRemovalFilter) {

        log2.fine("before LowIntensityRemovalFilter*.determineLowThreshold");

        plotHistogram0(stats.getHistogram(),  "border histogram");
    }

    private int debugfindEdgeWithRange(List<PairIntArray> edges,
        int xLo, int xHi, int yLo, int yHi) {

        for (int i = 0; i < edges.size(); i++) {
            boolean c = debugContainsAPoint(edges.get(i), xLo, xHi, yLo, yHi);
            if (c) {
                return i;
            }
        }
        return -1;
    }
    private boolean debugContainsAPoint(PairIntArray edge, 
        int xLo, int xHi, int yLo, int yHi) {
        for (int i = 0; i < edge.getN(); i++) {
            int x = edge.getX(i);
            int y = edge.getY(i);
            if (x >= xLo && x <= xHi && y >= yLo && y <= yHi) {
                return true;
            }
        }
        return false;
    }

    private void debugDisplay(Image input, String label) throws IOException {
     
        ImageDisplayer.displayImage(label, input);
    }
 
    private void debugDisplay(List<PairIntArray> curves, String label,
        boolean writeImage, int width, int height) throws IOException {
        
        log2.fine("create multicolor curve debug image");

        Image img2 = new Image(width, height);

        ImageIOHelper.addAlternatingColorCurvesToImage(curves, 
            "edges_" + outImgNum + ".png", writeImage, img2);        
        
        ImageDisplayer.displayImage(label, img2);
       
    }

    public void debugDisplay(GreyscaleImage input, String label) {
     
        try {
            ImageDisplayer.displayImage(label, input);

        } catch (IOException e) {

            log2.severe(e.getMessage());
        }
    }
    
    private void debugAddCurveToImage(float[] xPoints, float[] yPoints, 
        Image input, int nExtraForDot, int rClr, int gClr, int bClr) {
        
        for (int i = 0; i < xPoints.length; i++) {
            int x = (int) xPoints[i];
            int y = (int) yPoints[i];
            
            for (int dx = (-1*nExtraForDot); dx < (nExtraForDot + 1); dx++) {
                float xx = x + dx;
                if ((xx > -1) && (xx < (input.getWidth() - 1))) {
                    for (int dy = (-1*nExtraForDot); dy < (nExtraForDot + 1); dy++) {
                        float yy = y + dy;
                        if ((yy > -1) && (yy < (input.getHeight() - 1))) {
                            input.setRGB((int)xx, (int)yy, rClr, gClr, bClr);
                        }
                    }
                }
            }
        }
    }
    
    private void debugDisplayScaleSpaceImages(PairIntArray edge,
        Map<SIGMA, ScaleSpaceCurve> map, GreyscaleImage img) {
                
        try {
            //plot the edge and all scale space images
            Image debugImg = img.copyImageToGreen();
            
            ImageIOHelper.addCurveToImage(edge, debugImg, 0, 0, 255, 0);
            
//the convolutions should be acting upon data that starts at (0,0)
// and an offset applied?
// check an edge which is near x=0 and y=0
            
            Iterator<Entry<SIGMA, ScaleSpaceCurve>> iter = map.entrySet().iterator();
            while (iter.hasNext()) {
                Entry<SIGMA, ScaleSpaceCurve> entry = iter.next();
                ScaleSpaceCurve dss = entry.getValue();
                ImageIOHelper.addCurveToImage(dss, debugImg, 0, 0, 0, 255);
            }
            
            ImageDisplayer.displayImage("scale space images", debugImg);
                            
            int z = 1;
            
        } catch (IOException ex) {
            throw new RuntimeException("ERROR: l548" + ex.getMessage());
        }
    }

    private void debugDisplayScaleSpace(ScaleSpaceCurve maxScaleSpace, 
        List<Integer> maxCandidateCornerIndexes, PolygonAndPointPlotter plotter, 
        float[] xc, float[] yc, SIGMA sigma, int edgeNumber, 
        int width, int height) {
        
        float[] tCornersForPlot = new float[xc.length];
        float[] kCornersForPlot = new float[xc.length];
        for (int ii = 0; ii < maxCandidateCornerIndexes.size(); ii++) {
            int idx = maxCandidateCornerIndexes.get(ii);
            tCornersForPlot[ii] = maxScaleSpace.getT()[idx];
            kCornersForPlot[ii] = maxScaleSpace.getK(idx);
            if (kCornersForPlot[ii] < 0) {
                kCornersForPlot[ii] *= -1;
            }
        }
        try {
            float[] tmpX = maxScaleSpace.getT();
            float[] k = maxScaleSpace.getK();
            float[] tmpk = Arrays.copyOf(k, k.length);
            for (int ii = 1; ii < tmpk.length; ii++) {
                if (tmpk[ii] < 0) {
                    tmpk[ii] *= -1;
                }
            }
            
            plotter.addPlot(0.0f, 1.f, 0, MiscMath.findMax(tmpk), 
                tmpX, tmpk, tCornersForPlot, kCornersForPlot,
                "k for sigma=" + sigma.toString() 
                + " edge=" + edgeNumber);

            String filePath = plotter.writeFile();

            Image img3 = new Image(width, height);

            ImageIOHelper.addCurveToImage(maxScaleSpace, img3, 0, 255, 255, 
                255);

            debugAddCurveToImage(xc, yc, img3, 1, 255, 0, 0);

            ImageDisplayer.displayImage(
                "corner candidates (maxSigma=" + sigma.toString() 
                + ")", img3);

        } catch (Exception e) {
            throw new RuntimeException("ERROR: l595" + e.getMessage());
        }                
    }

    private void debugDisplay(ScaleSpaceCurve scaleSpace, List<Integer> 
        candidateCornerIndexes, String label, int width, int height) 
        throws IOException {
        
        Image img2 = new Image(width, height);
        
        ImageIOHelper.addCurveToImage(scaleSpace, img2, 0, 255, 255, 255);
      
        for (Integer idx : candidateCornerIndexes) {
            float x = scaleSpace.getX(idx.intValue());
            float y = scaleSpace.getY(idx.intValue());
            
            img2.setRGB((int)x, (int)y, 255, 0, 0);
            // make a plus symbol because the dot is too small
            for (int d = -1; d < 2; d++) {
                float xx = x + d;
                if ((xx > -1) && (xx < (width - 1))) {
                    img2.setRGB((int)xx, (int)y, 255, 0, 0);
                }
            }
            for (int d = -1; d < 2; d++) {
                float yy = y + d;
                if ((yy > -1) && (yy < (height - 1))) {
                    img2.setRGB((int)x, (int)yy, 255, 0, 0);
                }
            }
        }
        
        ImageDisplayer.displayImage(label, img2);
    }
      
    public double debugAreaUnderTheCurve(float[] x, float[] y) {
        
        // using trapezoidal rule
        
        if (x.length == 0) {
            return 0;
        } else if (x.length == 1) {
            // this isn't correct, but is fine for the one case it's used for
            return (y[0]);
        }
        
        double sum = 0;
        
        for (int i = 0; i < (x.length - 1); i++) {
            float yTerm = y[i + 1] + y[i];
            float xLen = x[i + 1] - x[i];
            sum += (yTerm * xLen);
        }
        
        sum *= 0.5;
        
        return sum;
    }
    
    private void debugDisplay(PairIntArray edge, PairIntArray corners,
        String label, int width, int height) throws IOException {
        
        Image img2 = new Image(width, height);
        
        ImageIOHelper.addCurveToImage(edge, img2, 0, 255, 255, 255);
        
        ImageIOHelper.addCurveToImage(corners, img2, 1, 255, 0, 0);
        
        ImageDisplayer.displayImage(label, img2);
     
        int z = 1;
    }
    
    private void plotHistogram0(HistogramHolder h, String label) {
                
        float[] xh = h.getXHist();
        float[] yh = h.getYHistFloat();
        
        float yMin = MiscMath.findMin(yh);
        int yMaxIdx = MiscMath.findYMaxIndex(yh);
        float yMax = yh[yMaxIdx];
        
        float xMin = MiscMath.findMin(xh);
        float xMax = MiscMath.findMax(xh);        
        
        log2.fine("MODE=" + xh[yMaxIdx]);
        
        try {

            PolygonAndPointPlotter plotter = new PolygonAndPointPlotter();
            
            plotter.addPlot(
                xMin, xMax, yMin, yMax,
                xh, yh, xh, yh, label);

            String filePath = plotter.writeFile3();
            
            log2.fine("created histogram at " + filePath);
            
        } catch (IOException e) {
            e.printStackTrace(); 
            log2.severe("ERROR while plotting histogram: " + e.getMessage());
        }
    }

    private void plotHistogram(GreyscaleImage img, String label) throws 
    IOException {
        
        float[] v = new float[img.getWidth() * img.getHeight()];
        float[] ve = new float[v.length];
        
        int count = 0;
        for (int i = 0; i < img.getWidth(); i++) {
            for (int j = 0; j < img.getHeight(); j++) {
                v[count] = img.getValue(i, j);
                ve[count] = (float)Math.sqrt(v[count]);
                count++;
            }
        }
        
        HistogramHolder h = Histogram.createSimpleHistogram(v, ve);
        //HistogramHolder h = Histogram.defaultHistogramCreator(v, ve);
        
        float[] xh = h.getXHist();
        float[] xhe = h.getXErrors();
        float[] yh = h.getYHistFloat();
        float[] yhe = h.getYErrors();
        int firstZero = -1;
        for (int i = (yh.length - 1); i > -1; i--) {
            if (yh[i] > 0) {
                break;
            } else {
                firstZero = i;
            }
        }
        if (firstZero > -1) {
            int n = xh.length - firstZero;
            xh = Arrays.copyOf(xh, n);
            xhe = Arrays.copyOf(xhe, n);
            yh = Arrays.copyOf(yh, n);
            yhe = Arrays.copyOf(yhe, n);
        }
        
        float xMax = MiscMath.findMax(xh);
        float yMax = MiscMath.findMax(yh);
        
        int yMaxIdx = MiscMath.findYMaxIndex(yh);
        yMax = yh[yMaxIdx + 3];
        
        PolygonAndPointPlotter plotter2 = new PolygonAndPointPlotter();
        
        plotter2.addPlot(0, xMax, 0, yMax, xh, yh, xhe, yhe, label);
        
        String filePath = plotter2.writeFile3();
        
        log2.info("filePath=" + filePath);
    }

    public void printContents(ScaleSpaceCurveImage spaceImage) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < spaceImage.getImageSigmas().length; i++) {
            sb.append("sigma=").append(spaceImage.getImageSigmas()[i]).append(": ");
            float[] t = spaceImage.getScaleSpaceImage()[i];
            for (int j = 0; j < t.length; j++) {
                sb.append("t=").append(t[j]).append(" ");
            }
            sb.append("\n");
        }
        log2.fine(sb.toString());
    }
   
    public String printImage(ScaleSpaceCurveImage spaceImage, 
        PolygonAndPointPlotter plotterC, int fileNumber) {
        
        if (plotterC == null) {
            return null;
        }
        int nTotalPoints = 0;
        for (int i = 0; i < spaceImage.getImageSigmas().length; i++) {
            float[] t = spaceImage.getScaleSpaceImage()[i];
            nTotalPoints+= t.length;
        }
        
        try {
            int count = 0;
            float[] x = new float[nTotalPoints];
            float[] y = new float[x.length];
            for (int i = 0; i < spaceImage.getImageSigmas().length; i++) {
                float sigma = spaceImage.getImageSigmas()[i];
                float[] t = spaceImage.getScaleSpaceImage()[i];
                for (int j = 0; j < t.length; j++) {
                    y[count] = sigma;
                    x[count] = t[j];
                    count++;
                }
            }
            float xMin = 0;
            float xMax = 1.0f;
            float yMin = MiscMath.findMin(y);
            if (yMin < 0) {
                yMin *= 1.1;
            } else {
                yMin *= 0.9;
            }
            float yMax = MiscMath.findMax(y);

            plotterC.addPlot(xMin, xMax,yMin, yMax, x, y, null, null,
                "");

            String filePathC = plotterC.writeFile(fileNumber);
            
            return filePathC;
            
        } catch(Exception e) {
            log2.severe(e.getMessage());
        }
        
        return null;
    }
}
