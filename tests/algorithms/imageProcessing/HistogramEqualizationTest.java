package algorithms.imageProcessing;

import algorithms.util.ResourceFinder;
import java.util.Arrays;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class HistogramEqualizationTest {
    
    public HistogramEqualizationTest() {
    }
    
    @Before
    public void setUp() {
    }
    
    @After
    public void tearDown() {
    }
    
    private GreyscaleImage getTestImage() {
        
        int[][] d = new int[8][8];
        d[0] = new int[]{52, 63, 62, 63, 67, 79, 85, 87};
        d[1] = new int[]{55, 59, 59, 58, 61, 65, 71, 79};
        d[2] = new int[]{61, 55, 68, 71, 68, 60, 64, 69};
        d[3] = new int[]{66, 90, 113, 122, 104, 70, 59, 68};
        d[4] = new int[]{70, 109, 144, 154, 126, 77, 55, 65};
        d[5] = new int[]{61, 85, 104, 106, 88, 68, 61, 76};
        d[6] = new int[]{64, 69, 66, 70, 68, 58, 65, 78};
        d[7] = new int[]{73, 72, 73, 69, 70, 75, 83, 94};
        
        GreyscaleImage img = new GreyscaleImage(8, 8);
        for (int i = 0; i < img.getWidth(); i++) {
            for (int j = 0; j < img.getWidth(); j++) {
                img.setValue(i, j, d[i][j]);
            }
        }
        
        return img;
    }

    @Test
    public void testCalculateHistogram() {
                
        GreyscaleImage img = getTestImage();
        
        HistogramEqualization instance = new HistogramEqualization(img);
        
        instance.calculateHistogram();
        
        long[] expectedCounts = new long[256];
        expectedCounts[52] = 1;
        expectedCounts[55] = 3;
        expectedCounts[58] = 2;
        expectedCounts[59] = 3;
        expectedCounts[60] = 1;
        expectedCounts[61] = 4;
        expectedCounts[62] = 1;
        expectedCounts[63] = 2;
        
        expectedCounts[64] = 2;
        expectedCounts[65] = 3;
        expectedCounts[66] = 2;
        expectedCounts[67] = 1;
        expectedCounts[68] = 5;
        expectedCounts[69] = 3;
        expectedCounts[70] = 4;
        expectedCounts[71] = 2;
        
        expectedCounts[72] = 1;
        expectedCounts[73] = 2;
        expectedCounts[75] = 1;
        expectedCounts[76] = 1;
        expectedCounts[77] = 1;
        expectedCounts[78] = 1;
        expectedCounts[79] = 2;
        expectedCounts[83] = 1;
        
        expectedCounts[85] = 2;
        expectedCounts[87] = 1;
        expectedCounts[88] = 1;
        expectedCounts[90] = 1;
        expectedCounts[94] = 1;
        expectedCounts[104] = 2;
        expectedCounts[106] = 1;
        expectedCounts[109] = 1;
        
        expectedCounts[113] = 1;
        expectedCounts[122] = 1;
        expectedCounts[126] = 1;
        expectedCounts[144] = 1;
        expectedCounts[154] = 1;
        
        assertTrue(Arrays.equals(instance.getHist(), expectedCounts));
    }

    @Test
    public void testCalculateCumulativeHistogram() {
        
        GreyscaleImage img = getTestImage();
        
        HistogramEqualization instance = new HistogramEqualization(img);
        
        instance.calculateHistogram();
        
        instance.calculateCumulativeHistogram();
        
        assertTrue(instance.getHistCMin() == 1);
        
        long[] expectedCounts = new long[256];
        expectedCounts[52] = 1;
        expectedCounts[55] = 4;
        expectedCounts[58] = 6;
        expectedCounts[59] = 9;
        expectedCounts[60] = 10;
        expectedCounts[61] = 14;
        expectedCounts[62] = 15;
        expectedCounts[63] = 17;
        
        expectedCounts[64] = 19;
        expectedCounts[65] = 22;
        expectedCounts[66] = 24;
        expectedCounts[67] = 25;
        expectedCounts[68] = 30;
        expectedCounts[69] = 33;
        expectedCounts[70] = 37;
        expectedCounts[71] = 39;
        
        expectedCounts[72] = 40;
        expectedCounts[73] = 42;
        expectedCounts[75] = 43;
        expectedCounts[76] = 44;
        expectedCounts[77] = 45;
        expectedCounts[78] = 46;
        expectedCounts[79] = 48;
        expectedCounts[83] = 49;
        
        expectedCounts[85] = 51;
        expectedCounts[87] = 52;
        expectedCounts[88] = 53;
        expectedCounts[90] = 54;
        expectedCounts[94] = 55;
        expectedCounts[104] = 57;
        expectedCounts[106] = 58;
        expectedCounts[109] = 59;
        
        expectedCounts[113] = 60;
        expectedCounts[122] = 61;
        expectedCounts[126] = 62;
        expectedCounts[144] = 63;
        expectedCounts[154] = 64;
        
        for (int i = 0; i < expectedCounts.length; i++) {
            //System.out.println(i + ") " + expectedCounts[i]
            //    + ":" + instance.getRHistC()[i]);
            // didn't include 'in between' points above, so only assert the above
            if (expectedCounts[i] > 0) {
                assertTrue(expectedCounts[i] == instance.getHistC()[i]);
            }
        }        
    }

    @Test
    public void testApplyTransformationFunction() {
        
        GreyscaleImage img = getTestImage();
        
        GreyscaleImage imgOrig = img.copyImage();
        
        HistogramEqualization instance = new HistogramEqualization(img);
        
        instance.applyFilter();
        
        long[] expectedCounts = new long[256];
        expectedCounts[52] = 1;
        expectedCounts[55] = 13;
        expectedCounts[58] = 21;
        expectedCounts[59] = 33;
        expectedCounts[60] = 37;
        expectedCounts[61] = 53;
        expectedCounts[62] = 57;
        expectedCounts[63] = 66;
        
        expectedCounts[64] = 74;
        expectedCounts[65] = 86;
        expectedCounts[66] = 94;
        expectedCounts[67] = 98;
        expectedCounts[68] = 118;
        expectedCounts[69] = 130;
        expectedCounts[70] = 146;
        expectedCounts[71] = 154;
        
        expectedCounts[72] = 158;
        expectedCounts[73] = 166;
        expectedCounts[75] = 170;
        expectedCounts[76] = 174;
        expectedCounts[77] = 178;
        expectedCounts[78] = 182;
        expectedCounts[79] = 190;
        expectedCounts[83] = 195;
        
        expectedCounts[85] = 203;
        expectedCounts[87] = 207;
        expectedCounts[88] = 211;
        expectedCounts[90] = 215;
        expectedCounts[94] = 219;
        expectedCounts[104] = 227;
        expectedCounts[106] = 231;
        expectedCounts[109] = 235;
        
        expectedCounts[113] = 239;
        expectedCounts[122] = 243;
        expectedCounts[126] = 247;
        expectedCounts[144] = 251;
        expectedCounts[154] = 255;
                      
        for (int i = 0; i < imgOrig.getWidth(); i++) {
            for (int j = 0; j < imgOrig.getHeight(); j++) {
                int v = imgOrig.getValue(i, j);
                long vExp = expectedCounts[v];
                long vT = img.getValue(i, j);
                assertTrue(vT == vExp);
            }
        }
    }
    
    @Test
    public void testAppearance() throws Exception {
        
        String fileName = "house.gif";
        //String fileName = "lab.gif";
        //String fileName = "susan-in.gif";
        //String fileName = "africa.png";
        //String fileName = "valve_gaussian.png";
        //String fileName = "lena.jpg";
        
        String filePath = ResourceFinder.findFileInTestResources(fileName);
        
        GreyscaleImage img = ImageIOHelper.readImageAsGrayScaleG(filePath);
        //Image img = ImageIOHelper.readImage(filePath);
    
        ImageDisplayer.displayImage("normal image", img);
        
        HistogramEqualization heq = new HistogramEqualization(img);
        heq.applyFilter();
        
        ImageDisplayer.displayImage("histogram equalized", img);
        
    }
    
}
