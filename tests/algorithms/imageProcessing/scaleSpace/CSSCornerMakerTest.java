package algorithms.imageProcessing.scaleSpace;

import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.SIGMA;
import algorithms.util.PolygonAndPointPlotter;
import algorithms.util.ResourceFinder;
import algorithms.util.CornerArray;
import algorithms.util.PairIntArray;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import junit.framework.TestCase;

/**
 * @author nichole
 */
public class CSSCornerMakerTest extends TestCase {
    
    public CSSCornerMakerTest(String testName) {
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
     
    public void testScaleSpaceImagesFigure2() throws Exception {
        
        System.out.println("testScaleSpaceImagesFigure2");
                                
        // IEEE 'TRANSACTIONS ON PATTERN ANALYSIS AND MACHINE INTELLIGENCE, 
        // VOL. PAMI-8, NO. 1. JANUARY 1986
        // Scale-Based Description and Recognition of Planar Curves and 
        // Two-Dimensional Shapes by FARZIN MOKHTARIAN AND ALAN MACKWORTH
        String filePath = ResourceFinder.findFileInTestResources("africa2.png");
        
        ImageExt img = ImageIOHelper.readImageExt(filePath);
        
        int imageWidth = img.getWidth();
        int imageHeight = img.getHeight();
        
        CurvatureScaleSpaceCornerDetector detector = new
            CurvatureScaleSpaceCornerDetector(img);
                        
        detector.useLineDrawingMode();
        
        detector.initialize();
        
        CSSCornerMaker cornerMaker = new CSSCornerMaker(imageWidth, imageHeight);
        
        List<Map<SIGMA, ScaleSpaceCurve>> outEdgeScaleSpaceMaps =
            new ArrayList<Map<SIGMA, ScaleSpaceCurve>>();
        
        //NOTE: this method is not efficient because it is expected the use
        // of CornerRegions will be obsolete after the next major refactoring.
        List<CornerArray> cornersList = cornerMaker.findCornersInScaleSpaceMaps(
            detector.getEdges(), outEdgeScaleSpaceMaps);
        
        PolygonAndPointPlotter plotterC = new PolygonAndPointPlotter();
        PolygonAndPointPlotter plotterX = new PolygonAndPointPlotter();
        PolygonAndPointPlotter plotterY = new PolygonAndPointPlotter();
        
        for (Map<SIGMA, ScaleSpaceCurve> mapOfScaleSpacesForAnEdge : outEdgeScaleSpaceMaps) {
                                    
            //======== draw the sigma=4  k vs t figure, x vs t, and y vs t ======
            
            ScaleSpaceCurve scaleSpaceCurveSigma = 
                mapOfScaleSpacesForAnEdge.get(SIGMA.FOUR);
            
            if (scaleSpaceCurveSigma == null || scaleSpaceCurveSigma.getSize() < 2) {
                continue;
            }
                        
            float[] xPoints = new float[scaleSpaceCurveSigma.getSize()];
            float[] yPoints = new float[xPoints.length];
            for (int ii = 0; ii < xPoints.length; ii++) {
                xPoints[ii] = ii;
            }
            float xMin = 0;
            float xMax = 1.1f * algorithms.misc.MiscMath.findMax(xPoints);
            float yMin, yMax;
            
            // ============ draw X(t,sigma) =============
            for (int ii = 0; ii < xPoints.length; ++ii) {
                yPoints[ii] = scaleSpaceCurveSigma.getX(ii);
            }
            yMin = 0.9f * algorithms.misc.MiscMath.findMin(yPoints);
            yMax = 1.1f * algorithms.misc.MiscMath.findMax(yPoints);

            plotterX.addPlot(
                0, xMax, yMin, yMax,
                null, null, xPoints, yPoints, 
                "t vs. X(t, sigma)");
            String filePathX = plotterX.writeFile2();

                        
            for (int ii = 0; ii < xPoints.length; ii++) {
                yPoints[ii] = scaleSpaceCurveSigma.getK(ii);
            }
            yMin = algorithms.misc.MiscMath.findMin(yPoints);
            if (yMin < 0) {
                yMin *= 1.1;
            } else {
                yMin *= 0.9;
            }
            yMax = algorithms.misc.MiscMath.findMax(yPoints);
            // ==== k vs t
            plotterC.addPlot(
                xMin, xMax, yMin, yMax,
                null, null, xPoints, yPoints,
                "t vs. curvature");

            String filePathC = plotterC.writeFile();
            
            // ============ draw Y(t,sigma) =============
            Arrays.fill(yPoints, 0);
            for (int ii = 0; ii < xPoints.length; ii++) {
                yPoints[ii] = scaleSpaceCurveSigma.getY(ii);
            }

            yMin = algorithms.misc.MiscMath.findMin(yPoints);
            if (yMin < 0) {
                yMin *= 1.1;
            } else {
                yMin *= 0.9;
            }
            yMax = 1.1f * algorithms.misc.MiscMath.findMax(yPoints);

            plotterY.addPlot(
                xMin, xMax, yMin, yMax,
                xPoints, yPoints, xPoints, yPoints, 
                "t vs. Y(t, sigma)");

            String filePath3 = plotterY.writeFile3();
            
            break;
        }
    }
    
    public void testPrintScaleSpaceCurves() throws Exception {
                                    
        String filePath = ResourceFinder.findFileInTestResources("lab.gif");
        
        ImageExt img = ImageIOHelper.readImageExt(filePath);
        
        int imageWidth = img.getWidth();
        int imageHeight = img.getHeight();
        
        CurvatureScaleSpaceCornerDetector detector = new
            CurvatureScaleSpaceCornerDetector(img);
                                
        detector.initialize();
        
       /* List<PairIntArray> edges = detector.getEdgesInOriginalReferenceFrame();
        Image img2 = img.copyImage();
        ImageIOHelper.addAlternatingColorCurvesToImage(edges, "_lab_EDGES_", false, img2);
        MiscDebug.writeImage(img2, "_lab_EDGES_");
        */
        
        CSSCornerMaker cornerMaker = new CSSCornerMaker(imageWidth, imageHeight);
        
        List<Map<SIGMA, ScaleSpaceCurve>> outEdgeScaleSpaceMaps =
            new ArrayList<Map<SIGMA, ScaleSpaceCurve>>();
        
        //NOTE: this method is not efficient because it is expected the use
        // of CornerRegions will be obsolete after the next major refactoring.
        List<CornerArray> cornersList = cornerMaker.findCornersInScaleSpaceMaps(
            detector.getEdges(), outEdgeScaleSpaceMaps);
                
        int xTest = 95;
        int yTest = 299;
        
        for (int lIdx = 0; lIdx < outEdgeScaleSpaceMaps.size(); ++lIdx) {
             
            Map<SIGMA, ScaleSpaceCurve> mapOfScaleSpacesForAnEdge = 
                outEdgeScaleSpaceMaps.get(lIdx);
                
            PairIntArray edge = detector.getEdges().get(lIdx);
            
            boolean found = false;
            for (int i = 0; i < edge.getN(); ++i) {
                if ((Math.abs(xTest - edge.getX(i)) < 5) && Math.abs(yTest - edge.getY(i)) < 5) {
                    found = true;
                    break;
                }
            }
            if (!found) {
                continue;
            }
            
            CornerArray ca = cornersList.get(lIdx);
                          
            /*
            need to print i (x,y) sigma k
            */
            for (Entry<SIGMA, ScaleSpaceCurve> entry2 : mapOfScaleSpacesForAnEdge.entrySet()) {
                SIGMA sigma = entry2.getKey();
                ScaleSpaceCurve scaleSpaceCurve = entry2.getValue();
                
                for (int i = 0; i < scaleSpaceCurve.getT().length; ++i) {
                    
                    float k = scaleSpaceCurve.getK(i);
                    
                    if (Math.abs(k) < 0.002) {
                        continue;
                    }
                    
                    String str = String.format("(%d,%d) k=%.4f sigma=%.3f", 
                        Math.round(scaleSpaceCurve.getX(i)),
                        Math.round(scaleSpaceCurve.getY(i)), k, 
                        SIGMA.getValue(sigma));
                    
                    System.out.println(str);
                }
            }
            System.out.flush();
            break;
        }
    }
    
    public static void main(String[] args) {
        try {
            
            CSSCornerMakerTest test = 
                new CSSCornerMakerTest("CSSCornerMakerTest");
                                    
            //test.testScaleSpaceImagesFigure2();
            
        } catch (Exception e) {
            e.printStackTrace();
            System.err.println("ERROR: " + e.getMessage());
        }
    }
}
