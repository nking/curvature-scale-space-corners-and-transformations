package algorithms.imageProcessing;

import algorithms.ResourceFinder;
import java.util.ArrayList;
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
public class CurvatureScaleSpaceContourMatcherTest {
    
    public CurvatureScaleSpaceContourMatcherTest() {
    }
    
    @Before
    public void setUp() {
    }
    
    @After
    public void tearDown() {
    }
    
    protected List<CurvatureScaleSpaceContour> getImageContours(String fileName)
        throws Exception {

        List<CurvatureScaleSpaceContour> contours = 
            new ArrayList<CurvatureScaleSpaceContour>();

        String filePath = ResourceFinder.findFileInTestResources(fileName);

        GreyscaleImage img = ImageIOHelper.readImageAsGrayScaleG(filePath);

        CurvatureScaleSpaceImageMaker instance
            = new CurvatureScaleSpaceImageMaker(img, true);

        List<PairIntArray> curves = instance.getClosedCurves();

        for (int i = 0; i < curves.size(); i++) {

            PairIntArray curve = curves.get(i);

            Map<Float, ScaleSpaceCurve> scaleSpaceMap
                = instance.createScaleSpaceMetricsForEdge2(curve);

            ScaleSpaceCurveImage scaleSpaceImage
                = instance.convertScaleSpaceMapToSparseImage(scaleSpaceMap, i);
                 
            ContourFinder contourFinder = new ContourFinder();

            List<CurvatureScaleSpaceContour> result = contourFinder.findContours(
                scaleSpaceImage, i);
        
            contours.addAll(result);
        }
        
        return contours;
    }

    @Test
    public void testMapScaleSpaceImageContours() throws Exception {
        
        List<CurvatureScaleSpaceContour> contours1 = 
            getImageContours("closed_curve.png");
        
        List<CurvatureScaleSpaceContour> contours2 = 
            getImageContours("closed_curve_translate.png");
        
        CurvatureScaleSpaceContourMatcher matcher = new
            CurvatureScaleSpaceContourMatcher(contours1, contours2);
        
        for (int i = 0; i < contours1.size(); i++) {            
            CurvatureScaleSpaceContour c = contours1.get(i);
            System.out.println("c1: " + i + ") sigma=" + c.getPeakSigma() 
                + " t=" + c.getPeakScaleFreeLength());
        }
        for (int i = 0; i < contours2.size(); i++) {
            CurvatureScaleSpaceContour c = contours2.get(i);
            System.out.println("c2: " + i + ") sigma=" + c.getPeakSigma() 
                + " t=" + c.getPeakScaleFreeLength());
        }
        
        matcher.matchContours();
       
        assertTrue(Math.abs(matcher.getSolvedShift()) < 0.1);
       
        assertTrue(Math.abs(matcher.getSolvedScale()) < 1.01);
    }

}
