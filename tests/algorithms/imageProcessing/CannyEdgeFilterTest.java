package algorithms.imageProcessing;

import algorithms.util.ResourceFinder;
import java.io.IOException;
import junit.framework.TestCase;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class CannyEdgeFilterTest extends TestCase {
    
    public CannyEdgeFilterTest() {
    }
        
    public void testLowIntensityFilter() throws Exception {
        
        String fileName = "house.gif";
        //String fileName = "lab.gif";
        //String fileName = "susan-in_plus.png";
        //String fileName = "africa.png";
        //String fileName = "valve_gaussian.png";
        //String fileName = "lena.jpg";
                
        String filePath = ResourceFinder.findFileInTestResources(fileName);
        
        GreyscaleImage img = ImageIOHelper.readImageAsGrayScaleG(filePath);
        
        CannyEdgeFilter filter = new CannyEdgeFilter();
        
        if (fileName.equals("susan-in_plus.png") || fileName.equals("africa.png")) {
            filter.useLineDrawingMode();
        }
        
        filter.applyFilter(img);
        
        //Image[] gradientProducts = filter.createGradientProducts(img);                
        //Image g = gradientProducts[2];

        int z = 1;
    }

    public void estApply2LayerFilter() throws Exception {
                
        String fileName = "house.gif";
        //String fileName = "lab.gif";
        //String fileName = "susan-in.gif";
        //String fileName = "susan-in_plus.png";
        //String fileName = "africa.png";
        //String fileName = "valve_gaussian.png";
        //String fileName = "lena.jpg";
                
        String filePath = ResourceFinder.findFileInTestResources(fileName);
        
        GreyscaleImage img = ImageIOHelper.readImageAsGrayScaleG(filePath);
                
        CannyEdgeFilter filter = new CannyEdgeFilter();
                                
        GreyscaleImage[] gradientProducts = filter.createGradientProducts(img);
        
        assertTrue(gradientProducts.length == 4);
        
        GreyscaleImage g = gradientProducts[2];

        filter.apply2LayerFilter(g);
               
        //TODO: add qualitative tests...
    }
    
    public void testApplyFilter() throws Exception {
                
        //String fileName = "house.gif";
        String fileName = "lab.gif";
        //String fileName = "susan-in.gif";
        //String fileName = "susan-in_plus.png";
        //String fileName = "africa.png";
        //String fileName = "valve_gaussian.png";
        //String fileName = "lena.jpg";
                
        String filePath = ResourceFinder.findFileInTestResources(fileName);
        
        GreyscaleImage img = ImageIOHelper.readImageAsGrayScaleG(filePath);
                
        CannyEdgeFilter filter = new CannyEdgeFilter();
        
        if (fileName.equals("africa.png")) {
            filter.useLineDrawingMode();
        } else if (fileName.contains("susan")) {
            filter.useLineDrawingMode();
        }
        
        filter.applyFilter(img);
        
              
        img.multiply(255);
        //ImageDisplayer.displayImage("canny edge filtered", img);
        String dirPath = ResourceFinder.findDirectory("bin");
        ImageIOHelper.writeOutputImage(dirPath + "/linethinned.png", img);

        GreyscaleImage img2 = filter.getGradientXY();
        img2.multiply(100);
        ImageIOHelper.writeOutputImage(dirPath + "/gXY.png", img2);
        
        //TODO: add qualitative tests...
        int z = 1;        
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
}
