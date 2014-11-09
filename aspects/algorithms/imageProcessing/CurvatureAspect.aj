package algorithms.imageProcessing;

import java.io.IOException;
import java.util.logging.Logger;
import algorithms.PolygonAndPointPlotter;
import algorithms.ResourceFinder;
import algorithms.misc.Histogram;
import algorithms.misc.HistogramHolder;
import algorithms.misc.MiscMath;
import algorithms.util.Errors;
import java.awt.Color;
import java.io.IOException;
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

    private Logger log2 = Logger.getAnonymousLogger();

    private int debugEdgeNumber = -1;

    private String currentEdgeStr = "";

    after() :
        target(algorithms.imageProcessing.ContourFinder) 
        && execution(public ContourFinder.new(..)) {

        log2.info("===> in aspect for ContourFinder");

        if (plotter9 == null) {
            try {
            plotter9 = new PolygonAndPointPlotter();
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

        log2.info("===> in aspect for ContourFinder, filePath=" + filePath);
    }

    after(ScaleSpaceCurveImage scaleSpaceImage, int edgeNumber) 
        returning(List<CurvatureScaleSpaceContour> contours) :
        execution(List<CurvatureScaleSpaceContour> algorithms.imageProcessing.ContourFinder.findContours(ScaleSpaceCurveImage, int))
        && args(scaleSpaceImage, edgeNumber) 
	    && target(algorithms.imageProcessing.ContourFinder) {

        Object obj = thisJoinPoint.getThis();

        if (!(obj instanceof ContourFinder)) {
            return;
        }

        ContourFinder instance = (ContourFinder)obj;

        for (int i = 0; i < contours.size(); i++) {

            CurvatureScaleSpaceContour contour = contours.get(i);

            float sigma = contour.getPeakSigma();
            float t = contour.getPeakScaleFreeLength();

            String str = String.format("edge=%d contour=%d (%f, %f)", 
                contour.getEdgeNumber(), i, sigma, t);

            log2.info(str);
        }
    }

    after() 
        returning(TransformationParameters params) :
        execution(public TransformationParameters algorithms.imageProcessing.CurvatureScaleSpaceInflectionMapper.createEuclideanTransformation())
	    && target(algorithms.imageProcessing.CurvatureScaleSpaceInflectionMapper) {

        Object obj = thisJoinPoint.getThis();

        if (!(obj instanceof CurvatureScaleSpaceInflectionMapper)) {
            return;
        }

        if (params == null) {
            return;
        }

        CurvatureScaleSpaceInflectionMapper instance = (CurvatureScaleSpaceInflectionMapper)obj;

        List<CurvatureScaleSpaceContour> contours1 = instance.getMatchedContours1();

        List<CurvatureScaleSpaceContour> contours2 = instance.getMatchedContours2();

        log2.info("number of matched contours1=" + contours1.size() +
            " number of matched contours2=" + contours2.size());

        for (int i = 0; i < contours1.size(); i++) {

            CurvatureScaleSpaceContour c1 = contours1.get(i);
            CurvatureScaleSpaceContour c2 = contours2.get(i);

            float sigma1 = c1.getPeakSigma();
            float t1 = c1.getPeakScaleFreeLength();
            int i1 = c1.getEdgeNumber();

            float sigma2 = c2.getPeakSigma();
            float t2 = c2.getPeakScaleFreeLength();
            int i2 = c2.getEdgeNumber();

            String str = String.format(
                "matched edge=%d (%f, %f) (%f, %f) edge=%d", i1, sigma1, t1, 
                sigma2, t2, i2);

            log2.info(str);
        }
        
        PairIntArray matchedXY1 = instance.getMatchedXY1();
        
        PairIntArray matchedXY2 = instance.getMatchedXY2();

        log2.info("number of matched points1=" + matchedXY1.getN() +
            " number of matched points2=" + matchedXY2.getN());

        for (int i = 0; i < matchedXY1.getN(); i++) {
            int x1 = matchedXY1.getX(i);
            int y1 = matchedXY1.getY(i);
            int x2 = matchedXY2.getX(i);
            int y2 = matchedXY2.getY(i);

            String str = String.format("matched points (%d, %d) (%d, %d)", 
                x1, y1, x2, y2);

            log2.info(str);
        }

        try {

            Image img1 = instance.getImage1().copyImageToGreen();
            
            ImageIOHelper.addCurveToImage(matchedXY1, img1, 1, 255, 0, 0);
            String dirPath = ResourceFinder.findDirectory("bin");
            ImageIOHelper.writeOutputImage(
                dirPath + "/contour_peaks1.png", img1);

            Image img2 = instance.getImage2().copyImageToGreen();
            
            ImageIOHelper.addCurveToImage(matchedXY2, img2, 1, 255, 0, 0);
            ImageIOHelper.writeOutputImage(
                dirPath + "/contour_peaks2.png", img2);

        } catch (IOException e) {
            log2.info("ERROR: " + e.getMessage());
        }

    }

    after() returning : 
	    target(algorithms.imageProcessing.CurvatureScaleSpaceCornerDetector) 
        && execution(protected void extractEdgeContours()) {
        
        Object obj = thisJoinPoint.getThis();

        if (!(obj instanceof CurvatureScaleSpaceCornerDetector)) {
            return;
        }

        CurvatureScaleSpaceCornerDetector instance = (CurvatureScaleSpaceCornerDetector)obj;

        List<PairIntArray> edges = instance.getEdges();

        debugEdgeNumber = debugfindEdgeWithRange(edges, 164, 168, 58, 63);

        log2.info("debugEdgeNumber=" + debugEdgeNumber);
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

        log2.info("after CannyEdgeFilter.applyFilter instance=" 
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

    after(GreyscaleImage input) returning() :
        call(void CannyEdgeFilter*.apply2LayerFilter(GreyscaleImage))
        && args(input)
	    && target(algorithms.imageProcessing.CannyEdgeFilter) {

        log2.info("after CannyEdgeFilter*.apply2LayerFilter");

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

    after(GreyscaleImage convolvedX, GreyscaleImage convolvedY) returning(GreyscaleImage output) :
        call(public GreyscaleImage ImageProcesser.computeTheta(GreyscaleImage, GreyscaleImage))
        && args(convolvedX, convolvedY)
	    && target(algorithms.imageProcessing.ImageProcesser) {

        log2.info("after ImageProcesser.computeTheta");

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

        log2.info("after ImageProcesser*.combineConvolvedImages");

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

            log2.info(currentEdgeStr + ") x=" + Arrays.toString(hist.getXHist()));

            log2.info(currentEdgeStr + ") y=" + Arrays.toString(hist.getYHist()));

        } catch(IOException e) {
           log2.severe(e.getMessage());
        }
	}

    after() returning : 
        execution(protected void algorithms.imageProcessing.AbstractCurvatureScaleSpaceMapper+.extractEdgeContours()) {

        Object obj = thisJoinPoint.getThis();

        if (!(obj instanceof AbstractCurvatureScaleSpaceMapper)) {
            return;
        }

        log2.info("==>in aspect of AbstractCurvatureScaleSpaceMapper");

        AbstractCurvatureScaleSpaceMapper instance = 
            (AbstractCurvatureScaleSpaceMapper)obj;

        try {
            List<PairIntArray> edges = instance.getEdges();

            Image img3 = instance.getOriginalImage().copyImageToGreen();

            debugDisplay(edges, "edges extracted", true, img3.getWidth(), 
                img3.getHeight());

        } catch (IOException e) {
            log2.severe(e.getMessage());
        }
    }

    // for some reason, execution works here, but call does not
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

        Image img3 = instance.getOriginalImage().copyImageToGreen();

        List<PairIntArray> edges = instance.getEdgesInOriginalReferenceFrame();

        PairIntArray corners = instance.getCornersInOriginalReferenceFrame();

        try {
            // multi-colored edges alone:
            debugDisplay(edges, "edges extracted", true, img3.getWidth(), img3.getHeight());

            for (PairIntArray edge : edges) {
                ImageIOHelper.addCurveToImage(edge, img3, 0, 255, 255, 0);
            }
            ImageIOHelper.addCurveToImage(corners, img3, 1, 255, 0, 0);
            debugDisplay(img3, "original image w/ edges and corners overplotted");
            String dirPath = ResourceFinder.findDirectory("bin");
            ImageIOHelper.writeOutputImage(
                dirPath + "/lab_with_edges_corners.png", img3);
 
        } catch (IOException ex) {
            throw new RuntimeException("ERROR: l298" + ex.getMessage());
        }
    }

    before(ScaleSpaceCurve scaleSpace, int edgeNumber, boolean rmFlseCrnrs, boolean isAClosedCurve) :
        target(algorithms.imageProcessing.CurvatureScaleSpaceCornerDetector)
        && call(protected PairFloatArray CurvatureScaleSpaceCornerDetector.findCornersInScaleSpaceMap(
        ScaleSpaceCurve, int, boolean, boolean))
        && args(scaleSpace, edgeNumber, rmFlseCrnrs, isAClosedCurve) {

        log2.info("before findCornersInScaleSpaceMap for edge " 
            + Integer.toString(edgeNumber) + ":");

        Object[] args = (Object[])thisJoinPoint.getArgs();
        Integer eNumber = (Integer)args[1];

        currentEdgeStr = eNumber.toString();

    }

    after(ScaleSpaceCurve scaleSpace, int edgeNumber, boolean rmFlseCrnrs,
        boolean isAClosedCurve) 
        returning(PairFloatArray xy) :
        target(algorithms.imageProcessing.CurvatureScaleSpaceCornerDetector)
        && call(protected PairFloatArray CurvatureScaleSpaceCornerDetector.findCornersInScaleSpaceMap(
        ScaleSpaceCurve, int, boolean, boolean))
        && args(scaleSpace, edgeNumber, rmFlseCrnrs, isAClosedCurve) {

        log2.info("after findCornersInScaleSpaceMap for edge " 
            + Integer.toString(edgeNumber) + ":");

        for (int i = 0; i < xy.getN(); i++) {
            log2.info(String.format("(%f, %f)", xy.getX(i), xy.getY(i)));
        }

    }

    before() 
        : call(List<PairIntArray> EdgeExtractor*.findEdges() ) 
        && args() 
        && target(algorithms.imageProcessing.EdgeExtractor) {

        log2.info("before EdgeExtractor*.findEdges()");

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

    before(ImageStatistics stats, HistogramHolder originalImageHistogram) 
        : call(int LowIntensityRemovalFilter*.determineLowThreshold( 
        ImageStatistics,  HistogramHolder))
        && args(stats, originalImageHistogram) 
        && target(algorithms.imageProcessing.LowIntensityRemovalFilter) {

        log2.info("before LowIntensityRemovalFilter*.determineLowThreshold");

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
        
        log2.info("create multicolor curve debug image");

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
        
        log2.info("MODE=" + xh[yMaxIdx]);
        
        try {

            PolygonAndPointPlotter plotter = new PolygonAndPointPlotter();
            
            plotter.addPlot(
                xMin, xMax, yMin, yMax,
                xh, yh, xh, yh, label);

            String filePath = plotter.writeFile3();
            
            log2.info("created histogram at " + filePath);
            
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
        log2.info(sb.toString());
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
