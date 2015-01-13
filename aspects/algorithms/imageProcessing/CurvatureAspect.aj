package algorithms.imageProcessing;

import java.io.IOException;
import java.util.logging.Logger;
import algorithms.util.*;
import algorithms.misc.Histogram;
import algorithms.misc.HistogramHolder;
import algorithms.misc.MiscMath;
import algorithms.util.Errors;
import java.awt.Color;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

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

        try {
            String dirPath = ResourceFinder.findDirectory("bin");

            ImageIOHelper.writeOutputImage(dirPath + "/edgethinned.png", image);

        } catch (IOException e) {
            log2.severe(e.getMessage());
        }
    }

    after(GreyscaleImage input, float sigma) returning() :
        call(public void ImageProcesser.blur(GreyscaleImage, float))
        && args(input, sigma)
	    && target(algorithms.imageProcessing.ImageProcesser) {

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
                dirPath + "/image_with_skyline.png", img3);
 
        } catch (IOException ex) {
            throw new RuntimeException("ERROR: l405" + ex.getMessage());
        }
    }

    after(GreyscaleImage input, int k) returning() :
        call(public void ImageProcesser.applyImageSegmentation(GreyscaleImage, int))
        && args(input, k)
	    && target(algorithms.imageProcessing.ImageProcesser) {

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

    after(GreyscaleImage convolvedX, GreyscaleImage convolvedY) returning(GreyscaleImage output) :
        call(public GreyscaleImage ImageProcesser.computeTheta(GreyscaleImage, GreyscaleImage))
        && args(convolvedX, convolvedY)
	    && target(algorithms.imageProcessing.ImageProcesser) {

        log2.fine("after ImageProcesser.computeTheta");

        debugDisplay(output, "theta");

        try {
            String dirPath = ResourceFinder.findDirectory("bin");
            ImageIOHelper.writeOutputImage(dirPath + "/gTheta.png", output);
        } catch (IOException e) {
            log2.severe(e.getMessage());
        }
    }

    after(GreyscaleImage convolvedX, GreyscaleImage convolvedY) returning(GreyscaleImage output) :
        call(public GreyscaleImage ImageProcesser.combineConvolvedImages(GreyscaleImage, GreyscaleImage))
        && args(convolvedX, convolvedY)
	    && target(algorithms.imageProcessing.ImageProcesser) {

        log2.fine("after ImageProcesser*.combineConvolvedImages");

        Object[] args = (Object[])thisJoinPoint.getArgs();
        GreyscaleImage cX = (GreyscaleImage)args[0];
        GreyscaleImage cY = (GreyscaleImage)args[1];

        debugDisplay(cX, "gradient X");
        debugDisplay(cY, "gradient Y");
        debugDisplay(output, "gradient combined");

        try {
            String dirPath = ResourceFinder.findDirectory("bin");
            ImageIOHelper.writeOutputImage(dirPath + "/gX.png", cX);
            ImageIOHelper.writeOutputImage(dirPath + "/gY.png", cY);
            ImageIOHelper.writeOutputImage(dirPath + "/gXY.png", output);
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

        ImageIOHelper.addAlternatingColorCurvesToImage(curves, "edges.png",
            writeImage, img2);        
        
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
