package algorithms.imageProcessing;

import algorithms.util.PolygonAndPointPlotter;
import algorithms.util.ResourceFinder;
import algorithms.misc.MiscMath;
import algorithms.util.PairIntArray;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import junit.framework.TestCase;

/**
 * @author nichole
 */
public class CurvatureScaleSpaceCornerDetectorTest extends TestCase {
    
    public CurvatureScaleSpaceCornerDetectorTest(String testName) {
        super(testName);
    }

    @Override
    protected void setUp() throws Exception {
        super.setUp();
    }

    @Override
    protected void tearDown() throws Exception {
        super.tearDown();
    }
     
    public void testConvolve_image() throws Exception {
        
        System.out.println("testConvolve_image");
        
        // this is a by-eye check of the convolution w/ the results of 
        // the paper's Figure 2
        
        String filePath = ResourceFinder.findFileInTestResources("africa2.png");
        
        GreyscaleImage img = ImageIOHelper.readImageAsGrayScaleG(filePath);
        
        Kernel1DHelper kernelDHelper = new Kernel1DHelper();
        
        CurvatureScaleSpaceCornerDetector detector = new
            CurvatureScaleSpaceCornerDetector(img);
        
        detector.useLineDrawingMode();
                
        detector.initialize();
        
        List<PairIntArray> edges = detector.getEdges();
                
        float[] g = Gaussian1D.getKernel(SIGMA.EIGHT);
        
        // if there's a 'knot' in the edge where points are out of order
        // in terms of the direction of the line, that is an endpoint is
        // actually out of order and placed in the middle of the points in the
        // array, the convolution will create spurious features due to
        // convolving points that are not spatially adjacent, their only
        // adjacent in the points array
        
        GreyscaleImage eImg = new GreyscaleImage(img.getWidth(), img.getHeight());
        GreyscaleImage cImg = new GreyscaleImage(img.getWidth(), img.getHeight());
        
        for (int idx = 0; idx < edges.size(); idx++) {
        
            PairIntArray edge = edges.get(idx);
            
            float minXValue = MiscMath.findMin(edge.getX());
        
            float minYValue = MiscMath.findMin(edge.getY());
                
            //PairIntArray edge = edges.get(0);
            PairIntArray convolvedX = new PairIntArray();

            for (int i = 0; i < edge.getN(); i++) {
                int x = edge.getX(i);
                int y = edge.getY(i);

                double lx = kernelDHelper.convolvePointWithKernel(edge, i, g, 
                    true);
                int lxRound = (int) Math.round(lx);
                convolvedX.add(lxRound, y);
                //cImg.setRGB((int)lx, y, 0, 255, 0);
                eImg.setValue(x, y, 255);
            }

            PairIntArray convolvedXY = new PairIntArray();
            for (int i = 0; i < edge.getN(); i++) {
                int x = convolvedX.getX(i);
                double ly = kernelDHelper.convolvePointWithKernel(convolvedX, 
                    i, g, false);
                convolvedXY.add(x, (int)ly);
                cImg.setValue(x, (int)ly, 255);
            }
        }
        
        ImageDisplayer.displayImage("edge", eImg);
        ImageDisplayer.displayImage("convolved edge", cImg);        
    }
    
    public void testScaleSpaceImagesFigure2() throws Exception {
        
        System.out.println("testScaleSpaceImagesFigure2");
                        
        String dirPath = ResourceFinder.findDirectory("bin");
        
        // IEEE 'TRANSACTIONS ON PATTERN ANALYSIS AND MACHINE INTELLIGENCE, 
        // VOL. PAMI-8, NO. 1. JANUARY 1986
        // Scale-Based Description and Recognition of Planar Curves and 
        // Two-Dimensional Shapes by FARZIN MOKHTARIAN AND ALAN MACKWORTH
        String filePath = ResourceFinder.findFileInTestResources("africa2.png");
        
        GreyscaleImage img = ImageIOHelper.readImageAsGrayScaleG(filePath);
        
        CurvatureScaleSpaceCornerDetector detector = new
            CurvatureScaleSpaceCornerDetector(img);
                
        detector.useLineDrawingMode();
        
        detector.initialize();
        
        Map<PairIntArray, Map<SIGMA, ScaleSpaceCurve> > map = 
            detector.findCornersInScaleSpaceMaps();
        
        
        
        Iterator<Entry<PairIntArray, Map<SIGMA, ScaleSpaceCurve> > > iter = 
            map.entrySet().iterator();
        
        PolygonAndPointPlotter plotterC = new PolygonAndPointPlotter();
        PolygonAndPointPlotter plotterX = new PolygonAndPointPlotter();
        PolygonAndPointPlotter plotterY = new PolygonAndPointPlotter();
        
        while (iter.hasNext()) {
            
            Entry<PairIntArray, Map<SIGMA, ScaleSpaceCurve> > entry = iter.next();
            
            PairIntArray edge = entry.getKey();
            
            //======== draw the sigma=4  k vs t figure, x vs t, and y vs t ======
            
            Map<SIGMA, ScaleSpaceCurve> mapOfScaleSpacesForAnEdge = entry.getValue();
            ScaleSpaceCurve scaleSpaceCurveSigma = 
                mapOfScaleSpacesForAnEdge.get(SIGMA.FOUR);
            float[] xPoints = new float[scaleSpaceCurveSigma.getSize()];
            float[] yPoints = new float[xPoints.length];
            for (int ii = 0; ii < xPoints.length; ii++) {
                xPoints[ii] = ii;
                yPoints[ii] = scaleSpaceCurveSigma.getK(ii);
            }
            float xMin = 0;
            float xMax = algorithms.misc.MiscMath.findMax(xPoints);
            xMax *= 1.1f;
            float yMin = algorithms.misc.MiscMath.findMin(yPoints);
            if (yMin < 0) {
                yMin *= 1.1;
            } else {
                yMin *= 0.9;
            }
            float yMax = algorithms.misc.MiscMath.findMax(yPoints);
            // ==== k vs t
            if (xPoints.length > 1) {
                plotterC.addPlot(
                    xMin, xMax,
                    yMin, yMax,
                    //-0.2f, 0.4f,
                    //-0.015f, 0.02f,
                    //-0.00015f, 0.0001f,
                    null, null, xPoints, yPoints,
                    "t vs. curvature");

                String filePathC = plotterC.writeFile();
            }
            // ============ draw X(t,sigma) =============
            Arrays.fill(yPoints, 0);
            for (int ii = 0; ii < xPoints.length; ii++) {
                yPoints[ii] = scaleSpaceCurveSigma.getX(ii);
            }

            xMin = 0;
            xMax = algorithms.misc.MiscMath.findMax(xPoints);

            xMax *= 1.1f;
            yMin = algorithms.misc.MiscMath.findMin(yPoints);
            if (yMin < 0) {
                yMin *= 1.1;
            } else {
                yMin *= 0.9;
            }
            yMax = algorithms.misc.MiscMath.findMax(yPoints);
            yMax *= 1.1f;

            if (xPoints.length > 1) {
                plotterX.addPlot(
                    0, xMax,
                    yMin, yMax,
                    null, null, xPoints, yPoints, 
                    "t vs. X(t, sigma)");
                String filePathX = plotterX.writeFile2();
            }

            // ============ draw Y(t,sigma) =============
            Arrays.fill(yPoints, 0);
            for (int ii = 0; ii < xPoints.length; ii++) {
                yPoints[ii] = scaleSpaceCurveSigma.getY(ii);
            }

            xMin = 0;
            xMax = algorithms.misc.MiscMath.findMax(xPoints);

            xMax *= 1.1f;
            yMin = algorithms.misc.MiscMath.findMin(yPoints);
            if (yMin < 0) {
                yMin *= 1.1;
            } else {
                yMin *= 0.9;
            }
            yMax = algorithms.misc.MiscMath.findMax(yPoints);
            yMax *= 1.1f;

            if (xPoints.length > 1) {
                plotterY.addPlot(
                    xMin, xMax,
                    yMin, yMax,
                    xPoints, yPoints, xPoints, yPoints, 
                    "t vs. Y(t, sigma)");

                String filePath3 = plotterY.writeFile3();
            }
            
            break;
        }
    }
    
    public void testScaleSpaceImagesFigure3() throws Exception {
        
        System.out.println("testScaleSpaceImagesFigure3");
                        
        String dirPath = ResourceFinder.findDirectory("bin");
        
        // IEEE 'TRANSACTIONS ON PATTERN ANALYSIS AND MACHINE INTELLIGENCE, 
        // VOL. PAMI-8, NO. 1. JANUARY 1986
        // Scale-Based Description and Recognition of Planar Curves and 
        // Two-Dimensional Shapes by FARZIN MOKHTARIAN AND ALAN MACKWORTH
        String filePath = ResourceFinder.findFileInTestResources("africa2.png");
        
        GreyscaleImage img = ImageIOHelper.readImageAsGrayScaleG(filePath);
        
        CurvatureScaleSpaceImageMaker imgMaker
            = new CurvatureScaleSpaceImageMaker(img);
                
        imgMaker.useLineDrawingMode();
        
        imgMaker.initialize();
        
        assertFalse(imgMaker.getClosedCurves().isEmpty());
        
        PairIntArray edge = imgMaker.getClosedCurves().get(0);
                
        Map<Float, ScaleSpaceCurve> scaleSpaceMap
            = imgMaker.createScaleSpaceMetricsForEdge2(edge);

        ScaleSpaceCurveImage scaleSpaceImage
            = imgMaker.convertScaleSpaceMapToSparseImage(scaleSpaceMap, 0,
                edge.getN());
        
        float[] sigmas = scaleSpaceImage.getImageSigmas();
        
        float[][] tVsSigma = scaleSpaceImage.getScaleSpaceImage();
        
        int n = 0;
        for (int col = 0; col < tVsSigma.length; col++) {
            n += tVsSigma[col].length;
        }
        
        float[] x = new float[n];
        float[] y = new float[n];
        n = 0;
        for (int col = 0; col < tVsSigma.length; col++) {
            float sigma = sigmas[col];
            for (int row = 0; row < tVsSigma[col].length; row++) {
                y[n] = sigma;
                x[n] = tVsSigma[col][row];
                n++;
            }
        }
        
        float xMin = 0;
        float xMax = 1.f;
        float yMin = 0;
        float yMax = 1.1f * algorithms.misc.MiscMath.findMax(y);
            
        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter();
        
        plotter.addPlot(xMin, xMax, yMin, yMax,
            x, y, null, null, 
            "t vs. sigma for inflection points");
        
        plotter.writeFile(5);
        
        // make the start of the t -axis t=0.58 and wrap around, then
        // reverse the t-axis and it matches Figure 3 of their paper
    }
    
    public static void main(String[] args) {
        try {
            
            CurvatureScaleSpaceCornerDetectorTest test = 
                new CurvatureScaleSpaceCornerDetectorTest("StackTest");
            
            test.testConvolve_image();
            
            test.testScaleSpaceImagesFigure2();
            
            test.testScaleSpaceImagesFigure3();
            
        } catch (Exception e) {
            e.printStackTrace();
            System.err.println("ERROR: " + e.getMessage());
        }
    }
}
