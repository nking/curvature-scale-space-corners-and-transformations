package algorithms.imageProcessing.features;

import algorithms.imageProcessing.GreyscaleImage;
import junit.framework.TestCase;
import static junit.framework.TestCase.assertEquals;

/**
 *
 * @author nichole
 */
public class PolarThetaIntegralHistogramsTest extends TestCase {
    
    public PolarThetaIntegralHistogramsTest() {
    }
    
    /*
            1     1         1        1
      16 32 48 64 80 96 112 128 144 160 176 192 208 224 240 256
      0  1  2  3  4  5   6  7   8    9  10  11  12  13  14   15
     */
    public void test0() {

        GreyscaleImage t = new GreyscaleImage(2, 2);

        t.setValue(0, 0, 47);
        t.setValue(0, 1, 127);
        t.setValue(1, 0, 159);
        t.setValue(1, 1, 79);

        PolarThetaIntegralHistograms gh = new PolarThetaIntegralHistograms();

        int[][] histograms = gh.createHistograms(t, 16);

        int[] outHist = new int[16];
        int[] outN = new int[1];

        assertTest0(gh, outHist, outN, histograms, t.getWidth(), t.getHeight());
    }

    private void assertTest0(PolarThetaIntegralHistograms gh, int[] outHist,
        int[] outN, int[][] histograms, int w, int h) {

        /*
                1     1         1        1
          16 32 48 64 80 96 112 128 144 160 176 192 208 224 240 256
          0  1  2  3  4  5   6  7   8    9  10  11  12  13  14   15
               0,0   1,1       0,1      1,0
        */
        
        gh.extractWindow(histograms, 0, 1, 0, 1, w, h, outHist, outN);
        //System.out.println("out=" + Arrays.toString(outHist));
        assertEquals(0, outHist[0]);
        assertEquals(0, outHist[1]);
        assertEquals(1, outHist[2]);
        assertEquals(1, outHist[4]);
        assertEquals(1, outHist[7]);
        assertEquals(1, outHist[9]);
        assertEquals(4, outN[0]);

        gh.extractWindow(histograms, 0, 1, 0, 0, w, h, outHist, outN);
        assertEquals(0, outHist[0]);
        assertEquals(0, outHist[1]);
        assertEquals(1, outHist[2]);
        assertEquals(0, outHist[4]);
        assertEquals(0, outHist[7]);
        assertEquals(1, outHist[9]);
        assertEquals(2, outN[0]);
        
        /*
                1     1         1        1
          16 32 48 64 80 96 112 128 144 160 176 192 208 224 240 256
          0  1  2  3  4  5   6  7   8    9  10  11  12  13  14   15
               0,0   1,1       0,1      1,0
        */
        gh.extractWindow(histograms, 0, 0, 0, 1, w, h, outHist, outN);
        assertEquals(0, outHist[0]);
        assertEquals(0, outHist[1]);
        assertEquals(1, outHist[2]);
        assertEquals(0, outHist[4]);
        assertEquals(1, outHist[7]);
        assertEquals(0, outHist[9]);
        assertEquals(2, outN[0]);
        

        /*
                1     1         1        1
          16 32 48 64 80 96 112 128 144 160 176 192 208 224 240 256
          0  1  2  3  4  5   6  7   8    9  10  11  12  13  14   15
               0,0   1,1       0,1      1,0
        */
        gh.extractWindow(histograms, 0, 0, 0, 0, w, h, outHist, outN);
        assertEquals(0, outHist[0]);
        assertEquals(0, outHist[1]);
        assertEquals(1, outHist[2]);
        assertEquals(0, outHist[4]);
        assertEquals(0, outHist[7]);
        assertEquals(0, outHist[9]);
        assertEquals(1, outN[0]);
        

        gh.extractWindow(histograms, 1, 1, 1, 1, w, h, outHist, outN);
        assertEquals(0, outHist[0]);
        assertEquals(0, outHist[1]);
        assertEquals(0, outHist[2]);
        assertEquals(1, outHist[4]);
        assertEquals(0, outHist[7]);
        assertEquals(0, outHist[9]);
        assertEquals(1, outN[0]);
    }

    /*
                  1         1 
            1     1         1        1
      16 32 48 64 80 96 112 128 144 160 176 192 208 224 240 256
      0  1  2  3  4  5   6  7   8    9  10  11  12  13  14   15
     */
    public void test1() {

        GreyscaleImage t = new GreyscaleImage(2, 3);

        t.setValue(0, 0, 47);
        t.setValue(0, 1, 127);
        t.setValue(1, 0, 159);
        t.setValue(1, 1, 79);
        t.setValue(0, 2, 120);
        t.setValue(1, 2, 70);

        int w = t.getWidth();
        int h = t.getHeight();

        PolarThetaIntegralHistograms gh = new PolarThetaIntegralHistograms();

        int[][] histograms = gh.createHistograms(t, 16);

        int[] outHist = new int[16];
        int[] outN = new int[1];

        assertTest0(gh, outHist, outN, histograms, w, h);

        /*
                      1         1 
                1     1         1        1
          16 32 48 64 80 96 112 128 144 160 176 192 208 224 240 256
          0  1  2  3  4  5   6  7   8    9  10  11  12  13  14   15
               0,0   1,1       0,1      1,0 
                     1,2       0,2
        */ 
        gh.extractWindow(histograms, 0, 1, 0, 2, w, h, outHist, outN);
        assertEquals(0, outHist[0]);
        assertEquals(0, outHist[1]);
        assertEquals(1, outHist[2]);
        assertEquals(2, outHist[4]);
        assertEquals(2, outHist[7]);
        assertEquals(1, outHist[9]);
        assertEquals(6, outN[0]);

        gh.extractWindow(histograms, 0, 1, 2, 2, w, h, outHist, outN);
        assertEquals(0, outHist[0]);
        assertEquals(0, outHist[1]);
        assertEquals(0, outHist[2]);
        assertEquals(1, outHist[4]);
        assertEquals(1, outHist[7]);
        assertEquals(0, outHist[9]);
        assertEquals(2, outN[0]);
    }

}
