package algorithms.imageProcessing.scaleSpace;

import algorithms.imageProcessing.Gaussian1D;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.ImageDisplayer;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.Kernel1DHelper;
import algorithms.imageProcessing.SIGMA;
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
        
        ImageExt img = ImageIOHelper.readImageExt(filePath);
        
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
        
String dirPath = ResourceFinder.findDirectory("bin");
String flPath = dirPath + "/sp0.png";
ImageIOHelper.writeOutputImage(flPath, eImg);
ImageIOHelper.writeOutputImage(dirPath + "/sp1.png", cImg);
        ImageDisplayer.displayImage("edge", eImg);
        ImageDisplayer.displayImage("convolved edge", cImg);        
    
    }
    
    public void testScaleSpaceImagesFigure3() throws Exception {
        
        System.out.println("testScaleSpaceImagesFigure3");
                        
        String dirPath = ResourceFinder.findDirectory("bin");
        
        // IEEE 'TRANSACTIONS ON PATTERN ANALYSIS AND MACHINE INTELLIGENCE, 
        // VOL. PAMI-8, NO. 1. JANUARY 1986
        // Scale-Based Description and Recognition of Planar Curves and 
        // Two-Dimensional Shapes by FARZIN MOKHTARIAN AND ALAN MACKWORTH
        String filePath = ResourceFinder.findFileInTestResources("africa2.png");
        
        ImageExt img = ImageIOHelper.readImageExt(filePath);
        
        CurvatureScaleSpaceImageMaker imgMaker
            = new CurvatureScaleSpaceImageMaker(img, true);
                
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
                new CurvatureScaleSpaceCornerDetectorTest(
                    "CurvatureScaleSpaceCornerDetectorTest");
            
            test.testConvolve_image();
                        
            test.testScaleSpaceImagesFigure3();
            
        } catch (Exception e) {
            e.printStackTrace();
            System.err.println("ERROR: " + e.getMessage());
        }
    }
}
