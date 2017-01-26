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
        
        g.setValue(0, 0, 10);  t.setValue(0, 0, 30);  
        g.setValue(0, 1, 220); t.setValue(0, 1, 90); 
        g.setValue(1, 0, 100); t.setValue(1, 0, 130); 
        g.setValue(1, 1, 10);  t.setValue(1, 1, 170); 
        
        GradientIntegralHistograms gh = new GradientIntegralHistograms(g, t, 9);
        
        int[] outHist = new int[9];
        int[] outN = new int[1];
        
        assertTest0(gh, outHist, outN);
    }
    
    private void assertTest0(GradientIntegralHistograms gh, int[] outHist,
        int[] outN) {
        
        gh.extractWindow(0, 1, 0, 1, outHist, outN);
        
        //     10      220   100   10
        //  0  1  2  3  4  5  6  7  8 
        
        assertEquals(10, outHist[1]);
        assertEquals(220, outHist[4]);
        assertEquals(100, outHist[6]);
        assertEquals(10, outHist[8]);
        assertEquals(4, outN[0]);
        
        //     10      220   100   10
        //  0  1  2  3  4  5  6  7  8 
        gh.extractWindow(0, 1, 0, 0, outHist, outN);
        
        assertEquals(10, outHist[1]);
        assertEquals(0, outHist[4]);
        assertEquals(100, outHist[6]);
        assertEquals(0, outHist[8]);
        assertEquals(2, outN[0]);
                
        //     10      220   100   10
        //  0  1  2  3  4  5  6  7  8 
        
        gh.extractWindow(0, 0, 0, 1, outHist, outN);
        /*
        g.setValue(0, 0, 10);  t.setValue(0, 0, 30);  
        g.setValue(0, 1, 220); t.setValue(0, 1, 90); 
        g.setValue(1, 0, 100); t.setValue(1, 0, 130); 
        g.setValue(1, 1, 10);  t.setValue(1, 1, 170); 
        */
        assertEquals(10, outHist[1]);
        assertEquals(220, outHist[4]);
        assertEquals(0, outHist[6]);
        assertEquals(0, outHist[8]);
        assertEquals(2, outN[0]);
        
        gh.extractWindow(0, 0, 0, 0, outHist, outN);
        /*
        g.setValue(0, 0, 10);  t.setValue(0, 0, 30);  
        g.setValue(0, 1, 220); t.setValue(0, 1, 90); 
        g.setValue(1, 0, 100); t.setValue(1, 0, 130); 
        g.setValue(1, 1, 10);  t.setValue(1, 1, 170); 
        */
        assertEquals(10, outHist[1]);
        assertEquals(0, outHist[4]);
        assertEquals(0, outHist[6]);
        assertEquals(0, outHist[8]);
        assertEquals(1, outN[0]);
        
        gh.extractWindow(1, 1, 1, 1, outHist, outN);
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
        
        g.setValue(0, 0, 10);  t.setValue(0, 0, 30);  
        g.setValue(0, 1, 220); t.setValue(0, 1, 90); 
        g.setValue(1, 0, 100); t.setValue(1, 0, 130); 
        g.setValue(1, 1, 10);  t.setValue(1, 1, 170); 
        g.setValue(0, 2, 10);  t.setValue(0, 2, 30); 
        g.setValue(1, 2, 225); t.setValue(1, 2, 90); 
        
        GradientIntegralHistograms gh = new GradientIntegralHistograms(g, t, 9);
        
        int[] outHist = new int[9];
        int[] outN = new int[1];
        
        assertTest0(gh, outHist, outN);
    
        gh.extractWindow(0, 1, 0, 2, outHist, outN);
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
        assertEquals(20, outHist[1]);
        assertEquals(225 + 220, outHist[4]);
        assertEquals(100, outHist[6]);
        assertEquals(10, outHist[8]);
        assertEquals(6, outN[0]);
        
        gh.extractWindow(0, 1, 2, 2, outHist, outN);
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
