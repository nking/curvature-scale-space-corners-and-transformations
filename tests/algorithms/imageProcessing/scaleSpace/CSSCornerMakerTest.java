package algorithms.imageProcessing.scaleSpace;

import algorithms.imageProcessing.scaleSpace.*;
import algorithms.imageProcessing.Gaussian1D;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.ImageDisplayer;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.Kernel1DHelper;
import algorithms.imageProcessing.SIGMA;
import algorithms.imageProcessing.util.PairIntWithIndex;
import algorithms.util.PolygonAndPointPlotter;
import algorithms.util.ResourceFinder;
import algorithms.misc.MiscMath;
import algorithms.util.PairIntArray;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
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
                        
        String dirPath = ResourceFinder.findDirectory("bin");
        
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
        
        PairIntArray outputCorners = new PairIntArray();
        Map<Integer, Set<PairIntWithIndex>> emptyJunctionsMap = 
            new HashMap<Integer, Set<PairIntWithIndex>>();
        
        CSSCornerMaker cornerMaker = new CSSCornerMaker(imageWidth, imageHeight);
        
        Map<PairIntArray, Map<SIGMA, ScaleSpaceCurve> > map = 
            cornerMaker.findCornersInScaleSpaceMaps(detector.getEdges(), 
                emptyJunctionsMap, outputCorners);
        
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
    
    public static void main(String[] args) {
        try {
            
            CSSCornerMakerTest test = 
                new CSSCornerMakerTest("CSSCornerMakerTest");
                                    
            test.testScaleSpaceImagesFigure2();
            
        } catch (Exception e) {
            e.printStackTrace();
            System.err.println("ERROR: " + e.getMessage());
        }
    }
}
