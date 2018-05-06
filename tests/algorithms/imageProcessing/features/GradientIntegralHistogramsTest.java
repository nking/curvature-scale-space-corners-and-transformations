package algorithms.imageProcessing.features;

import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.GreyscaleImage.Type;
import java.security.SecureRandom;
import java.util.Arrays;
import junit.framework.TestCase;
import static junit.framework.TestCase.assertEquals;

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
    
    private void assertTest0(GradientIntegralHistograms gh, long[] outHist,
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
    /*
          10      225
          10      220   100    10
       0  1  2  3  4  5  6  7  8 
    */
    public void test2() {
        
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
        
        long[] outHist = new long[9];
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
    
    public void testRandomExtract() throws Exception {
        
        SecureRandom random = SecureRandom.getInstance("SHA1PRNG");
        
        long seed = System.currentTimeMillis();
        //seed = 1525577808222L;
        System.out.println("SEED=" + seed);
        random.setSeed(seed);
        
        int nBins = 9;
        int w = 16;
        int h = 25;
        int nT = 100;//1<<29;
        
        // -255 to 255
        GreyscaleImage gradient = new GreyscaleImage(w, h, Type.Bits64Signed);
        // 0 to 180
        GreyscaleImage theta = new GreyscaleImage(w, h);
        
        for (int i = 0; i < theta.getNPixels(); ++i) {
            if (random.nextBoolean()) { continue;}
            
            //int r = random.nextInt(512);
            //r -= 256;
            int r = random.nextInt(256);
            gradient.setValue(i, r);
            
            r = random.nextInt(180);
            theta.setValue(i, r);
        }
        
        GradientIntegralHistograms gih = new GradientIntegralHistograms();
        int[][] hists = gih.createHistograms(gradient, theta, nBins);
        int[] outN = new int[1];
        int[] outH = new int[nBins];
        int[] sumH = new int[nBins];
        double[] bilinearParams = new double[5];
        double f0, f1;
        int t, v, b0, b1;
        
        gih._printHistograms_xy(hists, w, h);
        
        for (int i = 0; i < nT; ++i) {
            
            int x2 = random.nextInt(w - 1);
            x2++;
            int x1 = random.nextInt(x2);
            
            int y2 = random.nextInt(h - 1);
            y2++;
            int y1 = random.nextInt(y2);
            
            gih.extractWindow(hists, x1, x2, y1, y2, w, h, outH, outN);
            
            Arrays.fill(sumH, 0);
            for (int ii = x1; ii <= x2; ++ii) {
                for (int jj = y1; jj <= y2; ++jj) {
                    v = gradient.getValue(ii, jj);
                    
                    // t, b0, b1, f0, f1;
                    gih.calculateBinsAndFractions(ii, jj, theta, nBins, 
                        bilinearParams);
                    
                    t = (int)bilinearParams[0];
                    b0 = (int)bilinearParams[1];
                    b1 = (int)bilinearParams[2];
                    f0 = bilinearParams[3];
                    f1 = bilinearParams[4];
                    
                    sumH[b0] += (f0 * v);
                    sumH[b1] += (f1 * v);
                }
            }
            //System.out.println("try " + i + " outH=" + Arrays.toString(outH));
            //System.out.println("try " + i + " sumH=" + Arrays.toString(sumH));
            for (int kk = 0; kk < sumH.length; ++kk) {
                assertTrue(Math.abs(sumH[kk] - outH[kk]) <= 1);
            }
        }
    }
}
