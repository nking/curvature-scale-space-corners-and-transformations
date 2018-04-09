package algorithms.imageProcessing.features;

import algorithms.imageProcessing.GreyscaleImage;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class GradientIntegralHistogramsTest extends TestCase {
    
    public GradientIntegralHistogramsTest() {
    }
    
    /*
          10      220   100   10
       0  1  2  3  4  5  6  7  8 
    */
    public void test0() {
        
        GreyscaleImage g = new GreyscaleImage(2, 2);
        GreyscaleImage t = new GreyscaleImage(2, 2);
        
        g.setValue(0, 0, 10);  t.setValue(0, 0, 20);  
        g.setValue(0, 1, 220); t.setValue(0, 1, 90); 
        g.setValue(1, 0, 100); t.setValue(1, 0, 130); 
        g.setValue(1, 1, 10);  t.setValue(1, 1, 170); 
        
        GradientIntegralHistograms gh = new GradientIntegralHistograms();
        
        int[][] histograms = gh.createHistograms(g, t, 9);
        
        int[] outHist = new int[9];
        int[] outN = new int[1];
        
        assertTest0(gh, outHist, outN, histograms, g.getWidth(), g.getHeight());
    }
    
    private void assertTest0(GradientIntegralHistograms gh, int[] outHist,
        int[] outN, int[][] histograms, int w, int h) {
        
        gh.extractWindow(histograms, 0, 1, 0, 1, w, h, outHist, outN);
        
        //     10      220   100   10
        //  0  1  2  3  4  5  6  7  8 
        
        assertEquals(5, outHist[0]);
        assertEquals(5, outHist[1]);
        assertEquals(220, outHist[4]);
        assertEquals(100, outHist[6]);
        assertEquals(10, outHist[8]);
        assertEquals(4, outN[0]);
        
        //     10      220   100   10
        //  0  1  2  3  4  5  6  7  8 
        gh.extractWindow(histograms, 0, 1, 0, 0, w, h, outHist, outN);
        
        assertEquals(5, outHist[0]);
        assertEquals(5, outHist[1]);
        assertEquals(0, outHist[4]);
        assertEquals(100, outHist[6]);
        assertEquals(0, outHist[8]);
        assertEquals(2, outN[0]);
                
        //     10      220   100   10
        //  0  1  2  3  4  5  6  7  8 
        
        gh.extractWindow(histograms, 0, 0, 0, 1, w, h, outHist, outN);
        /*
        g.setValue(0, 0, 10);  t.setValue(0, 0, 30);  
        g.setValue(0, 1, 220); t.setValue(0, 1, 90); 
        g.setValue(1, 0, 100); t.setValue(1, 0, 130); 
        g.setValue(1, 1, 10);  t.setValue(1, 1, 170); 
        */
        assertEquals(5, outHist[0]);
        assertEquals(5, outHist[1]);
        assertEquals(220, outHist[4]);
        assertEquals(0, outHist[6]);
        assertEquals(0, outHist[8]);
        assertEquals(2, outN[0]);
        
        gh.extractWindow(histograms, 0, 0, 0, 0, w, h, outHist, outN);
        /*
        g.setValue(0, 0, 10);  t.setValue(0, 0, 30);  
        g.setValue(0, 1, 220); t.setValue(0, 1, 90); 
        g.setValue(1, 0, 100); t.setValue(1, 0, 130); 
        g.setValue(1, 1, 10);  t.setValue(1, 1, 170); 
        */
        assertEquals(5, outHist[0]);
        assertEquals(5, outHist[1]);
        assertEquals(0, outHist[4]);
        assertEquals(0, outHist[6]);
        assertEquals(0, outHist[8]);
        assertEquals(1, outN[0]);
        
        gh.extractWindow(histograms, 1, 1, 1, 1, w, h, outHist, outN);
        /*
        g.setValue(0, 0, 10);  t.setValue(0, 0, 30);  
        g.setValue(0, 1, 220); t.setValue(0, 1, 90); 
        g.setValue(1, 0, 100); t.setValue(1, 0, 130); 
        g.setValue(1, 1, 10);  t.setValue(1, 1, 170); 
        */
        assertEquals(0, outHist[1]);
        assertEquals(0, outHist[4]);
        assertEquals(0, outHist[6]);
        assertEquals(10, outHist[8]);
        assertEquals(1, outN[0]);
    }
    
    /*
          10      225
          10      220   100    10
       0  1  2  3  4  5  6  7  8 
    */
    public void test1() {
        
        GreyscaleImage g = new GreyscaleImage(2, 3);
        GreyscaleImage t = new GreyscaleImage(2, 3);
        
        /*
        binWidth = 20
        
          10   30   50  70   90   110  130  150  170
        0   20   40   60   80  100  120  140  160  180
        */
        
        g.setValue(0, 0, 10);  t.setValue(0, 0, 20);  
        g.setValue(0, 1, 220); t.setValue(0, 1, 90); 
        g.setValue(1, 0, 100); t.setValue(1, 0, 130); 
        g.setValue(1, 1, 10);  t.setValue(1, 1, 170); 
        g.setValue(0, 2, 10);  t.setValue(0, 2, 30); 
        g.setValue(1, 2, 225); t.setValue(1, 2, 90); 
        
        int w = g.getWidth();
        int h = g.getHeight();
        
        GradientIntegralHistograms gh = new GradientIntegralHistograms();
        
        int[][] histograms = gh.createHistograms(g, t, 9);
        
        int[] outHist = new int[9];
        int[] outN = new int[1];
        
        assertTest0(gh, outHist, outN, histograms, w, h);
    
        gh.extractWindow(histograms, 0, 1, 0, 2, w, h, outHist, outN);
        /*
        g.setValue(0, 0, 10);  t.setValue(0, 0, 30);  
        g.setValue(0, 1, 220); t.setValue(0, 1, 90); 
        g.setValue(1, 0, 100); t.setValue(1, 0, 130); 
        g.setValue(1, 1, 10);  t.setValue(1, 1, 170); 
        g.setValue(0, 2, 10);  t.setValue(0, 2, 30); 
        g.setValue(1, 2, 225); t.setValue(1, 2, 90);
              10      225
              10      220   100    10
           0  1  2  3  4  5  6  7  8 
        */
        assertEquals(5, outHist[0]);
        assertEquals(15, outHist[1]);
        assertEquals(225 + 220, outHist[4]);
        assertEquals(100, outHist[6]);
        assertEquals(10, outHist[8]);
        assertEquals(6, outN[0]);
        
        gh.extractWindow(histograms, 0, 1, 2, 2, w, h, outHist, outN);
        /*
        g.setValue(0, 0, 10);  t.setValue(0, 0, 30);  
        g.setValue(0, 1, 220); t.setValue(0, 1, 90); 
        g.setValue(1, 0, 100); t.setValue(1, 0, 130); 
        g.setValue(1, 1, 10);  t.setValue(1, 1, 170); 
        g.setValue(0, 2, 10);  t.setValue(0, 2, 30); 
        g.setValue(1, 2, 225); t.setValue(1, 2, 90);
              10      225
              10      220   100    10
           0  1  2  3  4  5  6  7  8 
        */
        assertEquals(10, outHist[1]);
        assertEquals(225, outHist[4]);
        assertEquals(0, outHist[6]);
        assertEquals(0, outHist[8]);
        assertEquals(2, outN[0]);
    }
    
}
