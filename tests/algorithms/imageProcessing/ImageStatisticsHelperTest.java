package algorithms.imageProcessing;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class ImageStatisticsHelperTest {
    
    public ImageStatisticsHelperTest() {
    }
    
    @Before
    public void setUp() {
    }
    
    @After
    public void tearDown() {
    }

    @Test
    public void testGetMean() {
        
        int v = 1;
        GreyscaleImage img = new GreyscaleImage(2, 2);
        for (int col = 0; col < img.getWidth(); col++) {
            for (int row = 0; row < img.getHeight(); row++) {
                img.setValue(col, row, v);
                v++;
            }
        }
        //1 2 3 4   10
        
        int result = ImageStatisticsHelper.getMean(img);
        
        assertTrue(2 == result);
        
        assertTrue(3 == ImageStatisticsHelper.getMedian(img));
        
    }

    @Test
    public void testGetQuartiles() {
        /*
                      median
             min        .         max
               .        .         .
               .   |    .    |    .
                q1   q2   q3   q4
        */
        
        int v = 3;
        GreyscaleImage img = new GreyscaleImage(4, 4);
        for (int col = 0; col < img.getWidth(); col++) {
            for (int row = 0; row < img.getHeight(); row++) {
                img.setValue(col, row, v);
            }
            v+=3;
        }
        // 3 3 3 3
        // 6 6 6 6
        //*9 9 9 9
        // 12 12 12 12
        
        /*
        int medianIdx = c.length >> 1;
        int q12Idx = (medianIdx - 1) >> 1;                      3    
        int q34Idx = (c.length + (medianIdx + 1))/2;            12
        return new float[]{c[q12Idx], c[medianIdx], c[q34Idx], c[c.length - 1]};
        */
        
        int i = 0;
        float[] imgValues = new float[img.getWidth() * img.getHeight()];
        for (int col = 0; col < img.getWidth(); col++) {
            for (int row = 0; row < img.getHeight(); row++) {
                imgValues[i] = img.getValue(col, row);
                i++;
            }
        }
        
        float[] quartiles = ImageStatisticsHelper.getQuartiles(imgValues);
        assertTrue(quartiles[0] == 3);
        assertTrue(quartiles[1] == 9);
        assertTrue(quartiles[2] == 12);
        assertTrue(quartiles[3] == 12);
   
        int borderWidth = 1;
        boolean useSturges = false;
        
        ImageStatistics result = ImageStatisticsHelper.examineImageBorders(img, 
            borderWidth, useSturges);
        
        assertTrue(result.getMean() == 7.5);
        assertTrue(result.getMin() == 3.0);
        assertTrue(result.getMax() == 12.0);
        
        result = ImageStatisticsHelper.examineImage(img, useSturges);
        assertTrue(result.getMean() == 7.5);
        assertTrue(result.getMin() == 3.0);
        assertTrue(result.getMax() == 12.0);
    }

}
