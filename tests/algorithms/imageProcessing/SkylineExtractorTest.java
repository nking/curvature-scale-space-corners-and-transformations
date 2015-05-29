package algorithms.imageProcessing;

import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.ResourceFinder;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;
import static junit.framework.Assert.assertTrue;
import junit.framework.TestCase;
import org.junit.Test;

/**
 *
 * @author nichole
 */
public class SkylineExtractorTest extends TestCase {
    
    private Logger log = Logger.getLogger(this.getClass().getName());
    
    public SkylineExtractorTest() {
    }
    
    public void estOrderByProximity() throws Exception {
                
        int x0 = 468;
        int y0 = 98;
        List<PairInt> expected = new ArrayList<PairInt>();
        expected.add(new PairInt(x0, y0));
        expected.add(new PairInt(469, 97));
        expected.add(new PairInt(469, 96));
        expected.add(new PairInt(469, 95));
        expected.add(new PairInt(470, 95));
        expected.add(new PairInt(471, 94));
        expected.add(new PairInt(472, 93));
        expected.add(new PairInt(473, 92));
        expected.add(new PairInt(474, 92));
        expected.add(new PairInt(475, 91));
        expected.add(new PairInt(476, 92));
        expected.add(new PairInt(477, 91));
        expected.add(new PairInt(478, 92));
        expected.add(new PairInt(477, 93));        
        expected.add(new PairInt(476, 94));
        expected.add(new PairInt(475, 95));
        expected.add(new PairInt(474, 95));
        expected.add(new PairInt(473, 95));
        expected.add(new PairInt(472, 96));
        expected.add(new PairInt(471, 97));
        expected.add(new PairInt(470, 98));
        expected.add(new PairInt(469, 98));
        
        Set<PairInt> points = new HashSet<PairInt>(expected);
        
        SkylineExtractor instance = new SkylineExtractor();
        List<PairInt> result = instance.orderByProximity(points);
        
        assertTrue(result.size() == expected.size());
        
        for (int i = 0; i < result.size(); ++i) {
            PairInt r = result.get(i);
            PairInt e = expected.get(i);
            assertTrue(r.equals(e));
        }
        
    }
    
    public void testOrderByIncreasingXThenY() throws Exception {
        
        int[] x = new int[]{7,  6,   5, 4,  4, 3, 2,  1, 0};
        int[] y = new int[]{10, 10, 10, 0, 10, 4, 10, 10, 10};
        
        int[] expectedX = new int[]{0,   1,  2,  3,  4,  4,  5,  6,  7};
        int[] expectedY = new int[]{10, 10, 10,  4,  0, 10, 10, 10, 10};
        
        Set<PairInt> points = new HashSet<PairInt>();
        for (int i = 0; i < x.length; ++i) {
            points.add(new PairInt(x[i], y[i]));
        }
        
        SkylineExtractor skylineExtractor = new SkylineExtractor();
        
        DoubleLinkedCircularList outputList = new DoubleLinkedCircularList();
        Map<Integer, HeapNode> outputMap = new
            HashMap<Integer, HeapNode>();
        
        skylineExtractor.orderByIncreasingXThenY(points, outputList, outputMap);
        
        int i = 0;
        HeapNode cn = outputList.getSentinel().getLeft();
        while (cn != outputList.getSentinel()) {
            PairInt p = (PairInt)cn.getData();
            assertTrue(p.getX() == expectedX[i]);
            assertTrue(p.getY() == expectedY[i]);
            ++i;
            cn = cn.getLeft();
        }
        
        assertTrue(outputMap.size() == 8);
    }
    
    @Test
    public void testCreateBestSkyMask() throws Exception {
        
        String[] fileNames = new String[] {
            /*"brown_lowe_2003_image1.jpg",
            //"brown_lowe_2003_image1_rot.jpg",
            //"brown_lowe_2003_image2.jpg",
            "venturi_mountain_j6_0001.png",
            //"venturi_mountain_j6_0010.png",
            "seattle.jpg",
            "arches.jpg",
            "stinson_beach.jpg",
            "cloudy_san_jose.jpg",  
            "stonehenge.jpg",
            "norwegian_mtn_range.jpg",
            "halfdome.jpg",
            "costa_rica.jpg",
            "new-mexico-sunrise_w725_h490.jpg",
            "arizona-sunrise-1342919937GHz.jpg",
            "sky_with_rainbow.jpg",
            "sky_with_rainbow2.jpg",*/
            "patagonia_snowy_foreground.jpg",
            "mt_rainier_snowy_field.jpg",    
            "klein_matterhorn_snowy_foreground.jpg"
            //"30.jpg",
            //"arches_sun_01.jpg",
            //"stlouis_arch.jpg", 
            //"contrail.jpg"
        };
        
        for (String fileName : fileNames) {
            
            log.info("fileName=" + fileName);
            
            // revisit infl points.  is there a threshold removing points?
            String filePath1 = ResourceFinder.findFileInTestResources(fileName);
            ImageExt img1 = ImageIOHelper.readImageExt(filePath1);
     
            SkylineExtractor.setDebugName(fileName);

            ImageHelperForTests helper = new ImageHelperForTests(img1, true);
            
            SkylineExtractor skylineExtractor = new SkylineExtractor();
                    
            PairIntArray outputSkyCentroid = new PairIntArray();
            
            // sky are the zeros in this:
            GreyscaleImage resultMask = skylineExtractor.createBestSkyMask(
                helper.getTheta(), helper.getGradientXY(), img1,
                helper.getCannyEdgeFilterSettings(), outputSkyCentroid);
            
            int idx = fileName.lastIndexOf(".");
            String fileNameRoot = fileName.substring(0, idx);
            String filePathSkyMask = ResourceFinder.findFileInTestResources(
                 fileNameRoot + "_sky.png");

            // sky are non-zeros here (255, 255, 255)
            GreyscaleImage skyMask = ImageIOHelper.readImageAsBinary(
                filePathSkyMask);
            
            int xOffset = resultMask.getXRelativeOffset();
            int yOffset = resultMask.getYRelativeOffset();
            
            long nExpected = 0;
            long nFound = 0;
            long nOverrun = 0;
            
            for (int i = 0; i < resultMask.getWidth(); i++) {
                for (int j = 0; j < resultMask.getHeight(); j++) {
                    
                    int x0 = i + xOffset;
                    int y0 = j + yOffset;
                    
                    int vExpected = skyMask.getValue(i, j);
                    int vResult = resultMask.getValue(x0, y0);
                    
                    if (vExpected > 0) {
                        nExpected++;
                        // only count found if expected
                        if (vResult == 0) {
                            nFound++;
                        }
                    } else {
                        if (vResult == 0) {
                            nOverrun++;
                        }
                    }
                }
            }
        
            double fracFound = (double)nFound/(double)nExpected;
            double fracOverrun = (double)nOverrun/(double)nExpected;
            
            //assertTrue(Math.abs(fracFound - 1) < 0.05);
            //assertTrue(fracOverrun < 0.05);
            System.out.println(fileName + " fraction of expected sky found = "
                + fracFound + "  fracOverrun = " + fracOverrun);
        }
    }
   
    public void testPopulatePolygon2() throws Exception {
        
        RainbowFinder rFinder = new RainbowFinder();
        
        /*
        y = x^2 + 3*x
        
        */
        
        float[] polynomialCoeff = new float[] {0.0f, 3.0f, 1.0f};
        
        float[] x = new float[]{-3.5f, -3.f,  -1.5f, 0.f, 0.5f};
        float[] y = new float[]{ 2.f,   0.f, -2.25f, 0.f, 2.0f};
        
        float dist = 1;
        
        float[] outputXPoly = new float[(2 * x.length) + 1];
        float[] outputYPoly = new float[outputXPoly.length];

        rFinder.populatePolygon(x, y, dist, outputXPoly, outputYPoly, 
            polynomialCoeff, 200, 200);

        for (int i = 0; i < outputXPoly.length; i++) {
            System.out.println("i=" + i 
                + " (" + outputXPoly[i] + " , " + outputYPoly[i] + ")");
        }
        
        int z = 1;
    }
    
    public static void main(String[] args) {
        
        try {
            SkylineExtractorTest test = new SkylineExtractorTest();

            test.testCreateBestSkyMask();
            //test.estOrderByProximity();
            //test.testOrderByIncreasingXThenY();
            //test.testPopulatePolygon();
            test.testPopulatePolygon2();
        
        } catch(Exception e) {
            e.printStackTrace();
            System.out.println(e.getMessage());
            fail(e.getMessage());
        }
    }

}
