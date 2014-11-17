package algorithms.imageProcessing;

import algorithms.ResourceFinder;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class ContourFinderTest {
    
    public ContourFinderTest() {
    }
    
    @Before
    public void setUp() {
    }
    
    @After
    public void tearDown() {
    }

    @Test
    public void testFindContours() throws Exception {
        
        String filePath = ResourceFinder.findFileInTestResources(
            "closed_curve.png");
        
        GreyscaleImage img = ImageIOHelper.readImageAsGrayScaleG(filePath);
        
        CurvatureScaleSpaceImageMaker instance = 
            new CurvatureScaleSpaceImageMaker(img);
        
        instance.useLineDrawingMode();
                
        instance.initialize();
        
        assertTrue(instance.getInitialized());
                
        List<PairIntArray> curves = instance.getClosedCurves();
        
        //Collections.sort(curves, new PairIntArrayComparator());
        
        assertTrue(curves.size() == 1);
        
        Map<Float, ScaleSpaceCurve> scaleSpaceMap = 
                instance.createScaleSpaceMetricsForEdge2(curves.get(0));
         
        ScaleSpaceCurveImage scaleSpaceImage = 
             instance.convertScaleSpaceMapToSparseImage(scaleSpaceMap, 0);
        
        ContourFinder contourFinder = new ContourFinder();
         
        List<CurvatureScaleSpaceContour> result = contourFinder.findContours(
            scaleSpaceImage, 0);
        
        assertTrue(result.size() == 3);
                
        CurvatureScaleSpaceContour contour = result.get(0);
        assertTrue(Math.abs(contour.getPeakScaleFreeLength() - 0.48f) < 0.1f);
        assertTrue(Math.abs(contour.getPeakSigma() - 32.0f) < 1.0f);
        assertTrue(contour.getEdgeNumber() == 0);
        assertTrue(contour.getPeakDetails().length == 2);
        
        contour = result.get(1);
        assertTrue(Math.abs(contour.getPeakScaleFreeLength() - 0.13f) < 0.1f);
        assertTrue(Math.abs(contour.getPeakSigma() - 26.9f) < 1.0f);
        assertTrue(contour.getEdgeNumber() == 0);
        assertTrue(contour.getPeakDetails().length == 2);
        
        contour = result.get(2);
        assertTrue(Math.abs(contour.getPeakScaleFreeLength() - 0.825f) < 0.1f);
        assertTrue(Math.abs(contour.getPeakSigma() - 26.9f) < 1.0f);
        assertTrue(contour.getEdgeNumber() == 0);
        assertTrue(contour.getPeakDetails().length == 2);
    }
    
    @Test
    public void testFindContours2() throws Exception {
        
        String filePath = ResourceFinder.findFileInTestResources(
            "blox.gif");
            //"closed_curve_translate.png");
        
        GreyscaleImage img = ImageIOHelper.readImageAsGrayScaleG(filePath);
        
        CurvatureScaleSpaceImageMaker instance = 
            new CurvatureScaleSpaceImageMaker(img);
                
        instance.initialize();
        
        assertTrue(instance.getInitialized());
                
        List<PairIntArray> curves = instance.getClosedCurves();
        
        //Collections.sort(curves, new PairIntArrayComparator());
                
         Map<Float, ScaleSpaceCurve> scaleSpaceMap = 
                instance.createScaleSpaceMetricsForEdge2(curves.get(0));
         
         ScaleSpaceCurveImage scaleSpaceImage = 
             instance.convertScaleSpaceMapToSparseImage(scaleSpaceMap, 0);
         
        ContourFinder contourFinder = new ContourFinder();
        List<CurvatureScaleSpaceContour> result = contourFinder.findContours(
            scaleSpaceImage, 0);
                        
        for (int i = 0; i < result.size(); i++) {
            
            CurvatureScaleSpaceContour c = result.get(i);
            
            CurvatureScaleSpaceImagePoint[] peakPoints = c.getPeakDetails();
            
            int[] x = new int[peakPoints.length];
            int[] y = new int[x.length];
            for (int j = 0; j < x.length; j++) {
                x[j] = peakPoints[j].getXCoord();
                y[j] = peakPoints[j].getYCoord();
            }
            
            System.out.println("CONTOUR: (" + c.getPeakSigma() + ", " 
                + c.getPeakScaleFreeLength() + ") at x=" + Arrays.toString(x)
                + " y=" + Arrays.toString(y));
        }
    }
    
    public static void main(String[] args) {
        
        try {
            
            ContourFinderTest test = new ContourFinderTest();
            
            test.testFindContours2();
            
        } catch (Exception e) {
            e.printStackTrace();
            System.err.println("ERROR: " + e.getMessage());
        }
    }
}
