package algorithms.imageProcessing.scaleSpace;

import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.util.PairInt;
import algorithms.util.ResourceFinder;
import algorithms.util.PairIntArray;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import junit.framework.TestCase;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class ContourFinderTest extends TestCase {
    
    public ContourFinderTest() {
    }
    
    public void testFindContours() throws Exception {
        
        String filePath = ResourceFinder.findFileInTestResources(
            "closed_curve.png");
        
        ImageExt img = ImageIOHelper.readImageExt(filePath);
        
        CurvatureScaleSpaceImageMaker instance = 
            new CurvatureScaleSpaceImageMaker(img);
                                                
        List<PairIntArray> curves = instance.getClosedCurves();
        
        //Collections.sort(curves, new PairIntArrayComparator());
        
        assertTrue(curves.size() == 1);
        
        Map<Float, ScaleSpaceCurve> scaleSpaceMap = 
                instance.createScaleSpaceMetricsForEdge2(curves.get(0));
         
        ScaleSpaceCurveImage scaleSpaceImage = 
             instance.convertScaleSpaceMapToSparseImage(scaleSpaceMap, 0,
             curves.get(0).getN());
        
        ContourFinder contourFinder = new ContourFinder();
         
        List<CurvatureScaleSpaceContour> result = contourFinder.findContours(
            scaleSpaceImage, 0);
        
        debugPlot(result, ImageIOHelper.readImageExt(filePath), 0, 0);
        
        assertTrue(result.size() >= 3);
        
        /*
        there are 3 main large peaks around the inner curves of the clover,
        and with a better line thinner, there are now peaks for the inflection
        points around the large lobes, but unevenly so.
        
        so, will assert the 3 large inner curve peaks
        */
        List<Float> expectedScaleFreeLength = new ArrayList<Float>();
        List<Float> expectedPeakSigma = new ArrayList<Float>();
        List<Integer> expectedEdgeNumber = new ArrayList<Integer>();
        List<PairInt> expectedPeakCoords0 = new ArrayList<PairInt>();
        List<PairInt> expectedPeakCoords1 = new ArrayList<PairInt>();
        
        expectedScaleFreeLength.add(Float.valueOf(0.51f));
        expectedPeakSigma.add(Float.valueOf(32));
        expectedEdgeNumber.add(Integer.valueOf(0));
        expectedPeakCoords0.add(new PairInt(24, 59));
        expectedPeakCoords1.add(new PairInt(23, 71));
         
        expectedScaleFreeLength.add(Float.valueOf(0.17f));
        expectedPeakSigma.add(Float.valueOf(27.5f));
        expectedEdgeNumber.add(Integer.valueOf(0));
        expectedPeakCoords0.add(new PairInt(58, 43));
        expectedPeakCoords1.add(new PairInt(51, 37));
          
        expectedScaleFreeLength.add(Float.valueOf(0.87f));
        expectedPeakSigma.add(Float.valueOf(27.5f));
        expectedEdgeNumber.add(Integer.valueOf(0));
        expectedPeakCoords0.add(new PairInt(52, 88));
        expectedPeakCoords1.add(new PairInt(60, 83));
        
        Set<Integer> found = new HashSet<Integer>();
        
        for (int i = 0; i < result.size(); i++) {
            
            CurvatureScaleSpaceContour contour = result.get(i);
            float sfl = contour.getPeakScaleFreeLength();
            float sigma = contour.getPeakSigma();
            int edgeNumber = contour.getEdgeNumber();
            CurvatureScaleSpaceImagePoint[] peakDetails = contour.getPeakDetails();
            
            for (int j = 0; j < expectedScaleFreeLength.size(); j++) {
                Integer index = Integer.valueOf(j);
                if (found.contains(index)) {
                    continue;
                }
                if (Math.abs(sfl - expectedScaleFreeLength.get(j)) < 0.15f) {
                    if (Math.abs(sigma - expectedPeakSigma.get(j)) < 5f) {
                        if (edgeNumber == expectedEdgeNumber.get(j)) {
                            PairInt exp0 = expectedPeakCoords0.get(j);
                            PairInt exp1 = expectedPeakCoords1.get(j);
                            
                            if (
                                ((Math.abs(peakDetails[0].getXCoord() - exp0.getX()) < 5)
                                && (Math.abs(peakDetails[0].getYCoord() - exp0.getY()) < 5))) {
                                if (
                                    ((Math.abs(peakDetails[1].getXCoord() - exp1.getX()) < 5)
                                    && (Math.abs(peakDetails[1].getYCoord() - exp1.getY()) < 5))) {
                                    
                                    found.add(index);
                                    break;
                                }
                            } else if (
                                ((Math.abs(peakDetails[1].getXCoord() - exp0.getX()) < 5)
                                && (Math.abs(peakDetails[1].getYCoord() - exp0.getY()) < 5))) {
                                if (
                                    ((Math.abs(peakDetails[0].getXCoord() - exp1.getX()) < 5)
                                    && (Math.abs(peakDetails[0].getYCoord() - exp1.getY()) < 5))) {
                                    
                                    found.add(index);
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }
        
        assertTrue(found.size() == 3);
    }
    
    public void estFindContours2() throws Exception {
        
        String filePath = ResourceFinder.findFileInTestResources(
            "blox.gif");
        
        ImageExt img = ImageIOHelper.readImageExt(filePath);
        
        CurvatureScaleSpaceImageMaker instance = 
            new CurvatureScaleSpaceImageMaker(img);
                        
        List<PairIntArray> curves = instance.getClosedCurves();
        
        //Collections.sort(curves, new PairIntArrayComparator());
                
         Map<Float, ScaleSpaceCurve> scaleSpaceMap = 
                instance.createScaleSpaceMetricsForEdge2(curves.get(0));
         
         ScaleSpaceCurveImage scaleSpaceImage = 
             instance.convertScaleSpaceMapToSparseImage(scaleSpaceMap, 0,
             curves.get(0).getN());
         
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
            
            test.testFindContours();
            
        } catch (Exception e) {
            e.printStackTrace();
            System.err.println("ERROR: " + e.getMessage());
        }
    }

    private void debugPlot(List<CurvatureScaleSpaceContour> result, ImageExt 
        img, int xOffset, int yOffset) throws IOException {
        
        if (result.isEmpty()) {
            return;
        }
        
        int nExtraForDot = 1;
        int rClr = 255;
        int gClr = 0;
        int bClr = 0;
        
        for (int i = 0; i < result.size(); i++) {
            
            CurvatureScaleSpaceContour cssC = result.get(i);
            
            CurvatureScaleSpaceImagePoint[] peakDetails = cssC.getPeakDetails();
            
            for (CurvatureScaleSpaceImagePoint peakDetail : peakDetails) {
                int x = peakDetail.getXCoord() + xOffset;
                int y = peakDetail.getYCoord() + yOffset;
                for (int dx = (-1*nExtraForDot); dx < (nExtraForDot + 1); dx++) {
                    float xx = x + dx;
                    if ((xx > -1) && (xx < (img.getWidth() - 1))) {
                        for (int dy = (-1*nExtraForDot); dy < (nExtraForDot + 1); 
                            dy++) {
                            float yy = y + dy;
                            if ((yy > -1) && (yy < (img.getHeight() - 1))) {
                                img.setRGB((int)xx, (int)yy, rClr, gClr, bClr);
                            }
                        }
                    }
                }
            }            
        }
        
        String dirPath = ResourceFinder.findDirectory("bin");
        
        ImageIOHelper.writeOutputImage(dirPath + "/contours.png", img);
    }
}
