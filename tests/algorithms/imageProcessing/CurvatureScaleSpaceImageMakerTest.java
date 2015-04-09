package algorithms.imageProcessing;

import algorithms.util.PolygonAndPointPlotter;
import algorithms.util.ResourceFinder;
import algorithms.misc.MiscMath;
import algorithms.util.PairIntArray;
import java.io.IOException;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class CurvatureScaleSpaceImageMakerTest {
    
    public CurvatureScaleSpaceImageMakerTest() {
    }
    
    @Before
    public void setUp() {
    }
    
    @After
    public void tearDown() {
    }
    
    private void addToPlot(PolygonAndPointPlotter plotter, float[]x, float[] y) 
        throws IOException {
        
        float xMin = 0;
        float xMax = MiscMath.findMax(x);
        xMax *= 1.1f;
        float yMin = MiscMath.findMin(y);
        if (yMin < 0) {
            yMin *= 1.1;
        } else {
            yMin *= 0.9;
        }
        float yMax = MiscMath.findMax(y);
        
        plotter.addPlot(xMin, xMax,yMin, yMax, x, y, null, null,
            "");

        String filePathC = plotter.writeFile();
        
        int z = 1;
    }
    
    @Test
    public void testCreateScaleSpaceMetricsForForEdge() throws Exception {
        
        // plot t vs sigma for a somewhat simple known curve
        
        String filePath = ResourceFinder.findFileInTestResources(
            "closed_curve.png");
        
        ImageExt img = ImageIOHelper.readImageExt(filePath);
        
        CurvatureScaleSpaceImageMaker instance = 
            new CurvatureScaleSpaceImageMaker(img);
        
        instance.useLineDrawingMode();
                
        instance.initialize();
        
        assertTrue(instance.getInitialized());
        
        PolygonAndPointPlotter plotterC = new PolygonAndPointPlotter();
        
        List<PairIntArray> curves = instance.getClosedCurves();
        
        //Collections.sort(curves, new PairIntArrayComparator());
        
        assertTrue(curves.size() == 1);
        
        boolean plotAs2 = false;
        
        if (plotAs2) {
            PairIntArray c = curves.get(0);
            int n0 = c.getN() >> 1;
            int n1 = c.getN() - n0;
            PairIntArray half0 = new PairIntArray(n0);
            PairIntArray half1 = new PairIntArray(n1);
            for (int i = 0; i < n0; i++) {
                half0.add(c.getX(i), c.getY(i));
            }
            for (int i = n0; i < c.getN(); i++) {
                half1.add(c.getX(i), c.getY(i));
            }
            curves.clear();
            curves.add(half0);
            curves.add(half1);
            
            System.out.println("curve0=" + half0);
            System.out.println("curve1=" + half1);
        }
        
        //assertTrue(curves.size() == 1);
                       
        for (int ii = 0; ii < curves.size(); ii++) {
            
            PairIntArray edge = curves.get(ii);
            
            Map<SIGMA, ScaleSpaceCurve> result = 
                instance.createScaleSpaceMetricsForEdge(edge);
            
            
            Iterator<Entry<SIGMA, ScaleSpaceCurve> > iter = 
                result.entrySet().iterator();

            int nTotalPoints = 0;
            while (iter.hasNext()) {
                Entry<SIGMA, ScaleSpaceCurve> entry = iter.next();
                ScaleSpaceCurve scaleSpaceCurve = entry.getValue();
                nTotalPoints += scaleSpaceCurve.getKIsZeroIdxSize();
            }

            if (nTotalPoints == 0) {
                throw new IllegalStateException("edge extractor + inflection mapper "
                + " returned no results");
            }

            iter = result.entrySet().iterator();

            // for each inflection point, plot idx/nPoints vs sigma
            float[] x = new float[nTotalPoints];
            float[] y = new float[x.length];
            int count = 0;
            while (iter.hasNext()) {

                Entry<SIGMA, ScaleSpaceCurve> entry = iter.next();

                SIGMA sigma = entry.getKey();

                ScaleSpaceCurve scaleSpaceCurve = entry.getValue();

                int nZeros = scaleSpaceCurve.getKIsZeroIdxSize();

                if (nZeros == 0) {
                    continue;
                }
                

                int nPoints = scaleSpaceCurve.getSize();

                for (int i = 0; i < nZeros; i++) {
                    int zIdx = scaleSpaceCurve.getKIsZeroIdx()[i];
                    x[count] = (float)zIdx/(float)nPoints;
                    y[count] = SIGMA.getValue(sigma);
                    //System.out.println("SIGM=" + sigma 
                    //    + "(" + scaleSpaceCurve.getX(zIdx) + "," 
                    //    + scaleSpaceCurve.getY(zIdx) + ") t=" + x[count] + " ii=" + ii);
                    count++;
                }
            }

            addToPlot(plotterC, x, y);
        }
    }
    
    @Test
    public void testCreateScaleSpaceMetricsForForEdge2() throws Exception {

        // plot t vs sigma for a somewhat simple known curve

        String filePath = ResourceFinder.findFileInTestResources(
            //"closed_curve.png");
            "closed_curve_translate_scale_rotate.png");
        
        ImageExt img = ImageIOHelper.readImageExt(filePath);
        
        CurvatureScaleSpaceImageMaker instance = 
            new CurvatureScaleSpaceImageMaker(img);
        
        instance.useLineDrawingMode();
        
        instance.initialize();
        
        assertTrue(instance.getInitialized());
                
        PolygonAndPointPlotter plotterC = new PolygonAndPointPlotter();
        
        List<PairIntArray> curves = instance.getClosedCurves();
        
        //Collections.sort(curves, new PairIntArrayComparator());
        
        assertTrue(curves.size() == 1);
        
        boolean plotAs2 = false;
        
        if (plotAs2) {
            PairIntArray c = curves.get(0);
            int n0 = c.getN() >> 1;
            int n1 = c.getN() - n0;
            PairIntArray half0 = new PairIntArray(n0);
            PairIntArray half1 = new PairIntArray(n1);
            for (int i = 0; i < n0; i++) {
                half0.add(c.getX(i), c.getY(i));
            }
            for (int i = n0; i < c.getN(); i++) {
                half1.add(c.getX(i), c.getY(i));
            }
            curves.clear();
            curves.add(half0);
            curves.add(half1);
            
            System.out.println("curve0=" + half0);
            System.out.println("curve1=" + half1);
        }
        
        //assertTrue(curves.size() == 1);
                       
        for (int ii = 0; ii < curves.size(); ii++) {
            
            PairIntArray edge = curves.get(ii);
            
            Map<Float, ScaleSpaceCurve> result = 
                instance.createScaleSpaceMetricsForEdge2(edge);
            
            ScaleSpaceCurveImage scaleSpaceCurveImage = 
                instance.convertScaleSpaceMapToSparseImage(result, ii, 
                edge.getN());
            
            float[] imageSigmas = scaleSpaceCurveImage.getImageSigmas();
            float[][] scaleSpaceImage = scaleSpaceCurveImage.getScaleSpaceImage();
            assertTrue(imageSigmas.length == scaleSpaceImage.length);
            for (int j = 0; j < imageSigmas.length; j++) {
                assertTrue(imageSigmas[j] != 0);
                for (int jj = 0; jj < scaleSpaceImage[j].length; jj++) {
                    assertTrue(scaleSpaceImage[j][jj] != 0);
                }
            }
            
            Iterator<Entry<Float, ScaleSpaceCurve> > iter = 
                result.entrySet().iterator();

            int nTotalPoints = 0;
            while (iter.hasNext()) {
                Entry<Float, ScaleSpaceCurve> entry = iter.next();
                ScaleSpaceCurve scaleSpaceCurve = entry.getValue();
                nTotalPoints += scaleSpaceCurve.getKIsZeroIdxSize();
            }

            if (nTotalPoints == 0) {
                throw new IllegalStateException(
                    "edge extractor + inflection mapper "
                    + " returned no results");
            }

            iter = result.entrySet().iterator();

            // for each inflection point, plot idx/nPoints vs sigma
            float[] x = new float[nTotalPoints];
            float[] y = new float[x.length];
            int count = 0;
            while (iter.hasNext()) {

                Entry<Float, ScaleSpaceCurve> entry = iter.next();

                Float sigma = entry.getKey();

                ScaleSpaceCurve scaleSpaceCurve = entry.getValue();

                int nZeros = scaleSpaceCurve.getKIsZeroIdxSize();

                if (nZeros == 0) {
                    continue;
                }

                int nPoints = scaleSpaceCurve.getSize();

                for (int i = 0; i < nZeros; i++) {
                    
                    int zIdx = scaleSpaceCurve.getKIsZeroIdx()[i];
                    
                    x[count] = (float)zIdx/(float)nPoints;
                    y[count] = sigma.floatValue();
                    
                    /*System.out.println("SIGMA=" + sigma 
                        + "(" + scaleSpaceCurve.getX(zIdx) + "," 
                        + scaleSpaceCurve.getY(zIdx) + ") t=" + x[count] 
                        + " ii=" + ii);
                    */
                    count++;
                }
            }

            addToPlot(plotterC, x, y);
        }
        int z = 1;
    }
    
    @Test
    public void testAlternateConstructor() throws Exception {

        // assert that using the image maker with already extracted edges,
        //   produces same results as when it invokes the edge extractor.

        String filePath = ResourceFinder.findFileInTestResources(
            //"closed_curve.png");
            "closed_curve_translate_scale_rotate.png");
        
        ImageExt img = ImageIOHelper.readImageExt(filePath);
        
        CurvatureScaleSpaceImageMaker instance = 
            new CurvatureScaleSpaceImageMaker(img);
        
        instance.useLineDrawingMode();
        
        instance.initialize();
        
        assertTrue(instance.getInitialized());
                        
        List<PairIntArray> curves = instance.getClosedCurves();
        
        
        CurvatureScaleSpaceImageMaker instance2 = 
            new CurvatureScaleSpaceImageMaker(img, curves);
        
        instance2.useLineDrawingMode();
        
        instance2.initialize();
        
        assertTrue(instance2.getInitialized());
                
        List<PairIntArray> curves2 = instance2.getClosedCurves();
        
        assertTrue(curves.size() == curves2.size());
        
        for (int ii = 0; ii < curves.size(); ii++) {
            
            PairIntArray edge = curves.get(ii);
            
            Map<Float, ScaleSpaceCurve> result = 
                instance.createScaleSpaceMetricsForEdge2(edge);
            
            ScaleSpaceCurveImage scaleSpaceCurveImage = 
                instance.convertScaleSpaceMapToSparseImage(result, ii, 
                edge.getN());
            
            
            PairIntArray edge2 = curves2.get(ii);
            
            Map<Float, ScaleSpaceCurve> result2 = 
                instance2.createScaleSpaceMetricsForEdge2(edge2);
            
            ScaleSpaceCurveImage scaleSpaceCurveImage2 = 
                instance2.convertScaleSpaceMapToSparseImage(result2, ii, 
                edge2.getN());
            
            
            float[] sigmas = scaleSpaceCurveImage.getImageSigmas();
            float[] sigmas2 = scaleSpaceCurveImage2.getImageSigmas();
            
            float[][] scaleImage = scaleSpaceCurveImage.getScaleSpaceImage();
            float[][] scaleImage2 = scaleSpaceCurveImage2.getScaleSpaceImage();
            
            assertTrue(Arrays.equals(sigmas, sigmas2));
            
            assertTrue(scaleImage.length == scaleImage2.length);
            
            for (int i = 0; i < scaleImage.length; i++) {
                assertTrue(Arrays.equals(scaleImage[i], scaleImage2[i]));
            }
        }
    }
    
    public static void main(String[] args) {
        
        try {
            
            CurvatureScaleSpaceImageMakerTest test = 
                new CurvatureScaleSpaceImageMakerTest();
            
            
            test.testCreateScaleSpaceMetricsForForEdge2();
            
        } catch (Exception e) {
            e.printStackTrace();
            System.err.println("ERROR: " + e.getMessage());
        }
    }
    
}
