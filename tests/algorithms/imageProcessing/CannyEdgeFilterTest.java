package algorithms.imageProcessing;

import algorithms.util.PairInt;
import algorithms.util.ResourceFinder;
import java.io.IOException;
import java.util.HashSet;
import java.util.Set;
import java.util.logging.Logger;
import junit.framework.TestCase;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class CannyEdgeFilterTest extends TestCase {
    
    protected Logger log = Logger.getLogger(this.getClass().getName());
    
    public CannyEdgeFilterTest() {
    }
    
    public void testApplyFilter() throws Exception {
                
        String[] fileNames = new String[] {
            "blox.gif",
            "house.gif",
            "lab.gif",
            "africa.png",
            "susan-in.gif",
            "valve_gaussian.png",
            "lena.jpg"
        };
        
        for (String fileName : fileNames) {
                        
            String filePath1 = ResourceFinder.findFileInTestResources(fileName);
            
            GreyscaleImage img = ImageIOHelper.readImageAsGrayScaleG(filePath1);
            
            CannyEdgeFilter filter = new CannyEdgeFilter();
        
            if (fileName.equals("africa.png")) {
                filter.useLineDrawingMode();
            }
            
            int idx = fileName.lastIndexOf(".");
            String fileNameRoot = fileName.substring(0, idx);
            
            filter.applyFilter(img);
            GreyscaleImage img2 = img.copyImage();
            img.multiply(200);
            String dirPath = ResourceFinder.findDirectory("bin");
            ImageIOHelper.writeOutputImage(dirPath 
                + "/linethinned_" + fileNameRoot +".png", img);
        
            GreyscaleImage gXY = filter.getGradientXY();
            
            // ----- assert line width is 1 pix everywhere exception junctions.
            // ---- this tests results of the line thinner and post-line
            // ---- thinner corrections            
            gXY.multiply(100);
            ImageIOHelper.writeOutputImage(dirPath + "/gXY+" 
                + fileNameRoot + ".png", gXY);
            
            assertThinnedLines(img2);
        }
             
    }
    
    public void testApplyHistogramEqualization() {
        
        GreyscaleImage input = new GreyscaleImage(10, 10);
        
        int v = 0;
        for (int col = 0; col < input.getWidth(); col++) {
            for (int row = 0; row < input.getHeight(); row++) {
                if (v > 127) {
                    v = 0;
                }
                input.setValue(col, row, v);
                v++;
            }
        }
        
        CannyEdgeFilter instance = new CannyEdgeFilter();
        
        instance.applyHistogramEqualization(input);
        
        // rough check that min and max are filling 0 to 255.  the more detailed
        // tests for the histogram method are in HistogramEqualizationTest
        int min = Integer.MAX_VALUE;
        int max = Integer.MIN_VALUE;
        for (int col = 0; col < input.getWidth(); col++) {
            for (int row = 0; row < input.getHeight(); row++) {
                
                v = input.getValue(col, row);
                if (v < min) {
                    min = v;
                }
                if (v > max) {
                    max = v;
                }
            }
        }
        System.out.println("min=" + min + " max=" + max);
        assertTrue(min == 0);
        assertTrue(max == 255);
    }
    
    public void testCreateGradientProducts() throws IOException {
        
        double dTheta = 10.0;
        
        int n = (int)(360.f/dTheta);
        
        float r = (float)(((float)n) * 10.f/(2. * Math.PI));//10.0f;
        
        float xc = r + 5;
        float yc = xc;
        
        int w = (int)(xc * 2.f);
        int h = w;
        
        GreyscaleImage img = new GreyscaleImage(w, h);
        
        float expectedCurvature = (1.f/r);
                         
        int pointFactor = 10;
        
        float rDiffMax = 5.0f;
        
        for (int i = 0; i < n*pointFactor; i++) {
            /*
            (x-xc)^2 + (y-yc)^2 = r
            x = xc + r*cos(theta)
            y = yc + r*sin(theta)
            */
            double thetaRadians = (Math.PI*(i*dTheta)/180.)/pointFactor;
            
            double cos = Math.cos(thetaRadians);
            double sin = Math.sin(thetaRadians);
            
            int x = (int)(xc + (r * cos));
            int y = (int)(yc + (r * sin));
            
            img.setValue(x, y, 127);
            
            for (float dr = 0.5f; dr < rDiffMax; dr += 0.25f) {
                img.setValue((int)(xc + ((r - dr) * cos)),
                    (int)(yc + ((r - dr) * sin)), 127);
            }            
        }
        
        //String dirPath = ResourceFinder.findDirectory("bin");
        
        //ImageIOHelper.writeOutputImage(dirPath + "/tmp.png", img);
        
        CannyEdgeFilter instance = new CannyEdgeFilter();
                
        GreyscaleImage[] result = instance.createGradientProducts(img);
        
        /*
        ImageIOHelper.writeOutputImage(dirPath + "/tmpGX.png", result[0]);
        ImageIOHelper.writeOutputImage(dirPath + "/tmpGY.png", result[1]);
        ImageIOHelper.writeOutputImage(dirPath + "/tmpGXY.png", result[2]);
        ImageIOHelper.writeOutputImage(dirPath + "/tmpTheta.png", result[3]);
        */
        
        GreyscaleImage theta = result[3];
        for (int col = 0; col < theta.getWidth(); col++) {
            for (int row = 0; row < theta.getHeight(); row++) {
                
                int pix = theta.getValue(col, row);
                
                // gradX and gradY are both 0 ==> theta = 0
                if (pix != 0) {
                    
                    float radius = (float) Math.sqrt(Math.pow(col - xc, 2) 
                        + Math.pow(row - yc, 2));
                    
                    float diff = r - radius;
                    if (diff < 0) {
                        diff *= -1;
                    }
                    
                    //System.out.println(diff + " : " + (rDiffMax + 1));
                    
                    assertTrue(diff < (rDiffMax + 2));
                }
            }
        }
    }

    public static void main(String[] args) {
        
        try {
            
            CannyEdgeFilterTest test = new CannyEdgeFilterTest();
            
            /*test.testApplyHistogramEqualization();
            
            test.testLowIntensityFilter();
            
            test.testApply2LayerFilter();
            
            test.testCreateGradientProducts();
            */
            test.testApplyFilter();
            
        } catch (Exception e) {
            e.printStackTrace();
            System.err.println("ERROR: " + e.getMessage());
        }
    }

    private void assertThinnedLines(GreyscaleImage thinnedImg) {
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        final int w = thinnedImg.getWidth();
        final int h = thinnedImg.getHeight();
        
        Set<PairInt> points = imageProcessor.readNonZeroPixels(thinnedImg);
        
        assertFalse(points.isEmpty());
        
        PostLineThinnerCorrections pltc = new PostLineThinnerCorrections();
        
        /*    0         2
           0  #  #  0   1
        0  #* #  0      0
           #  0        -1
       -1  0  1  2  3
        */
        Set<PairInt> zeroes = new HashSet<PairInt>();
        Set<PairInt> ones = new HashSet<PairInt>();
        zeroes.add(new PairInt(-1, 0));
        zeroes.add(new PairInt(0, -1));
        zeroes.add(new PairInt(1, 1)); zeroes.add(new PairInt(1, -2));
        zeroes.add(new PairInt(2, 0));
        zeroes.add(new PairInt(3, -1));
        
        ones.add(new PairInt(0, 1));
        ones.add(new PairInt(1, 0));
        ones.add(new PairInt(1, -1));
        ones.add(new PairInt(2, -1));
        
        for (int nRot = 0; nRot < 4; ++nRot) {
        
            switch(nRot) {
                case 1:
                    pltc.reverseXs(zeroes, ones);
                    break;
                case 2:
                    pltc.reverseYs(zeroes, ones);
                    break;
                case 3:
                    pltc.reverseXs(zeroes, ones);
                    break;
                default:
                    break;
            }
            
            for (PairInt p : points) {

                int col = p.getX();
                int row = p.getY();           

                int patternCount = 0;

                for (PairInt p2 : ones) {
                    int x = col + p2.getX();
                    int y = row + p2.getY();
                    if ((x < 0) || (y < 0) || (x > (w - 1)) || (y > (h - 1))) {
                        patternCount = 0;
                        break;
                    }
                    PairInt p3 = new PairInt(x, y);
                    if (points.contains(p3)) {
                        patternCount++;
                    }
                }
                
                for (PairInt p2 : zeroes) {
                    int x = col + p2.getX();
                    int y = row + p2.getY();
                    if ((x < 0) || (y < 0) || (x > (w - 1)) || (y > (h - 1))) {
                        patternCount = 0;
                        break;
                    }
                    PairInt p3 = new PairInt(x, y);
                    if (!points.contains(p3)) {
                        patternCount++;
                    }
                }
                
                if (patternCount == (ones.size() + zeroes.size())) {
                    
                    pltc.debugPrint(points, 
                        new HashSet<PairInt>(), new HashSet<PairInt>(),
                        col - 2, col + 2, row - 2, row + 2);
        
                    fail("col=" + (col + thinnedImg.getXRelativeOffset())
                        + " row=" + (row + thinnedImg.getYRelativeOffset()));
                }
            }
        }
    }
}
