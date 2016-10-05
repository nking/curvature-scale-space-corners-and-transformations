package algorithms.imageProcessing;

import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class SummedAreaTableTest extends TestCase {
    
    public SummedAreaTableTest() {
    }
    
    public void testSummedAreaTable() {
        
        /*
        column major notation is used
        
        2  5  2  7   10  8 18   10 18  36
        1  3  1  2    5  6 11    5 11  22
        0  2  5  9    2  5  9    2  7  16
           0  1  2
         */
        GreyscaleImage img = new GreyscaleImage(3, 3);
        img.setValue(0, 0, 2); img.setValue(0, 1, 3); img.setValue(0, 2, 5);
        img.setValue(1, 0, 5); img.setValue(1, 1, 1); img.setValue(1, 2, 2);
        img.setValue(2, 0, 9); img.setValue(2, 1, 2); img.setValue(2, 2, 7);
       
        SummedAreaTable sumTable = new SummedAreaTable();
        GreyscaleImage sImg = sumTable.createAbsoluteSummedAreaTable(img);
        
        int[][] expected = new int[3][];
        expected[0] = new int[]{2, 5, 10};
        expected[1] = new int[]{7, 11, 18};
        expected[2] = new int[]{16, 22, 36};
        
        for (int i = 0; i < img.getWidth(); ++i) {
            for (int j = 0; j < img.getHeight(); ++j) {
                assertEquals(expected[i][j], sImg.getValue(i, j));
            }
        }
        
        /*
        img              sImg
        values           values
        2  5  2  7     10 18  36
        1  3  1  2      5 11  22
        0  2  5  9      2  7  16
           0  1  2
         
        extract areas from table.                
        */
        int expectedSum;
        int[] output = new int[2];
        
        
        expectedSum = 36;
        sumTable.extractWindowFromSummedAreaTable(sImg, 
            0, 0, 5, output);
        assertEquals(output[0], expectedSum);
        assertEquals(output[1], 9);
        
        expectedSum = 11;
        sumTable.extractWindowFromSummedAreaTable(sImg, 
            0, 0, 3, output);
        assertEquals(output[0], expectedSum);
        assertEquals(output[1], 4);
        
        expectedSum = 2;
        sumTable.extractWindowFromSummedAreaTable(sImg, 
            0, 0, 1, output);
        assertEquals(output[0], expectedSum);
        assertEquals(output[1], 1);
        
        
        /*
        img              sImg
        values           values
        2  5  2  7     10 18  36
        1  3  1  2      5 11  22
        0  2  5  9      2  7  16
           0  1  2
        */
        
        expectedSum = 36;
        sumTable.extractWindowFromSummedAreaTable(sImg, 
            0, 1, 5, output);
        assertEquals(expectedSum, output[0]);
        assertEquals(9, output[1]);
        
        expectedSum = 18;
        sumTable.extractWindowFromSummedAreaTable(sImg, 
            0, 1, 3, output);
        assertEquals(expectedSum, output[0]);
        assertEquals(6, output[1]);
        
        expectedSum = 3;
        sumTable.extractWindowFromSummedAreaTable(sImg, 
            0, 1, 1, output);
        assertEquals(expectedSum, output[0]);
        assertEquals(1, output[1]);
        
        
       /*
        img              sImg
        values           values
        2  5  2  7     10 18  36
        1  3  1  2      5 11  22
        0  2  5  9      2  7  16
           0  1  2
        */
        expectedSum = 36;
        sumTable.extractWindowFromSummedAreaTable(sImg, 
            0, 2, 5, output);
        assertEquals(expectedSum, output[0]);
        assertEquals(9, output[1]);
        
        expectedSum = 18 - 7;
        sumTable.extractWindowFromSummedAreaTable(sImg, 
            0, 2, 3, output);
        assertEquals(expectedSum, output[0]);
        assertEquals(4, output[1]);
        
        expectedSum = 5;
        sumTable.extractWindowFromSummedAreaTable(sImg, 
            0, 2, 1, output);
        assertEquals(expectedSum, output[0]);
        assertEquals(1, output[1]);
        
        /*
        img              sImg
        values           values
        2  5  2  7     10 18  36
        1  3  1  2      5 11  22
        0  2  5  9      2  7  16
           0  1  2
        */
        expectedSum = 36;
        sumTable.extractWindowFromSummedAreaTable(sImg, 
            1, 0, 5, output);
        assertEquals(output[0], expectedSum);
        assertEquals(output[1], 9);
        
        expectedSum = 22;
        sumTable.extractWindowFromSummedAreaTable(sImg, 
            1, 0, 3, output);
        assertEquals(output[0], expectedSum);
        assertEquals(output[1], 6);
        
        expectedSum = 5;
        sumTable.extractWindowFromSummedAreaTable(sImg, 
            1, 0, 1, output);
        assertEquals(output[0], expectedSum);
        assertEquals(output[1], 1);
        
        /*
        img              sImg
        values           values
        2  5  2  7     10 18  36
        1  3  1  2      5 11  22
        0  2  5  9      2  7  16
           0  1  2
        */
        expectedSum = 36;
        sumTable.extractWindowFromSummedAreaTable(sImg, 
            1, 1, 5, output);
        assertEquals(output[0], expectedSum);
        assertEquals(output[1], 9);
        
        expectedSum = 36;
        sumTable.extractWindowFromSummedAreaTable(sImg, 
            1, 1, 3, output);
        assertEquals(output[0], expectedSum);
        assertEquals(output[1], 9);
        
        expectedSum = 1;
        sumTable.extractWindowFromSummedAreaTable(sImg, 
            1, 1, 1, output);
        assertEquals(output[0], expectedSum);
        assertEquals(output[1], 1);
        
        //0,0  0,1  0,2  1,0  1,1  1,2
        
        /*
        img              sImg
        values           values
        2  5  2  7     10 18  36
        1  3  1  2      5 11  22
        0  2  5  9      2  7  16
           0  1  2
        */
        expectedSum = 36;
        sumTable.extractWindowFromSummedAreaTable(sImg, 
            1, 2, 5, output);
        assertEquals(output[0], expectedSum);
        assertEquals(output[1], 9);
        
        expectedSum = 20;
        sumTable.extractWindowFromSummedAreaTable(sImg, 
            1, 2, 3, output);
        assertEquals(output[0], expectedSum);
        assertEquals(output[1], 6);
        
        expectedSum = 2;
        sumTable.extractWindowFromSummedAreaTable(sImg, 
            1, 2, 1, output);
        assertEquals(output[0], expectedSum);
        assertEquals(output[1], 1);
        
        
        /*
        img              sImg
        values           values
        2  5  2  7     10 18  36
        1  3  1  2      5 11  22
        0  2  5  9      2  7  16
           0  1  2
        */
        expectedSum = 36;
        sumTable.extractWindowFromSummedAreaTable(sImg, 
            2, 0, 5, output);
        assertEquals(output[0], expectedSum);
        assertEquals(output[1], 9);
        
        expectedSum = 17;
        sumTable.extractWindowFromSummedAreaTable(sImg, 
            2, 0, 3, output);
        assertEquals(expectedSum, output[0]);
        assertEquals(4, output[1]);
        
        expectedSum = 9;
        sumTable.extractWindowFromSummedAreaTable(sImg, 
            2, 0, 1, output);
        assertEquals(expectedSum, output[0]);
        assertEquals(1, output[1]);
        
        /*
        img              sImg
        values           values
        2  5  2  7     10 18  36
        1  3  1  2      5 11  22
        0  2  5  9      2  7  16
           0  1  2
        */
        expectedSum = 36;
        sumTable.extractWindowFromSummedAreaTable(sImg, 
            2, 1, 5, output);
        assertEquals(output[0], expectedSum);
        assertEquals(output[1], 9);
        
        expectedSum = 26;
        sumTable.extractWindowFromSummedAreaTable(sImg, 
            2, 1, 3, output);
        assertEquals(output[0], expectedSum);
        assertEquals(output[1], 6);
        
        expectedSum = 2;
        sumTable.extractWindowFromSummedAreaTable(sImg, 
            2, 1, 1, output);
        assertEquals(expectedSum, output[0]);
        assertEquals(1, output[1]);
        
        /*
        img              sImg
        values           values
        2  5  2  7     10 18  36
        1  3  1  2      5 11  22
        0  2  5  9      2  7  16
           0  1  2
        */
        expectedSum = 36;
        sumTable.extractWindowFromSummedAreaTable(sImg, 
            2, 2, 5, output);
        assertEquals(output[0], expectedSum);
        assertEquals(output[1], 9);
        
        expectedSum = 12;
        sumTable.extractWindowFromSummedAreaTable(sImg, 
            2, 2, 3, output);
        assertEquals(expectedSum, output[0]);
        assertEquals(4, output[1]);
        
        expectedSum = 7;
        sumTable.extractWindowFromSummedAreaTable(sImg, 
            2, 2, 1, output);
        assertEquals(output[0], expectedSum);
        assertEquals(output[1], 1);
        
        // --- assert the mean of window size 5 ----
        GreyscaleImage mImg = sumTable.
            applyMeanOfWindowFromSummedAreaTable(sImg, 5);
        for (int i = 0; i < mImg.getNPixels(); ++i) {
            assertEquals(36/9, mImg.getValue(i));
        }
        
        // --- assert the mean of window size 3 ----
        mImg = sumTable.applyMeanOfWindowFromSummedAreaTable(sImg, 3);
        assertEquals(11/4, mImg.getValue(0, 0));
        assertEquals(18/6, mImg.getValue(0, 1));
        assertEquals(11/4, mImg.getValue(0, 2));
        assertEquals(22/6, mImg.getValue(1, 0));
        assertEquals(36/9, mImg.getValue(1, 1));
        assertEquals(20/6, mImg.getValue(1, 2));
        assertEquals(17/4, mImg.getValue(2, 0));
        assertEquals(26/6, mImg.getValue(2, 1));
        assertEquals(12/4, mImg.getValue(2, 2));
        
        // --- assert the mean of window size 1 ----
        mImg = sumTable.applyMeanOfWindowFromSummedAreaTable(sImg, 1);
        assertEquals(img.getValue(0, 0), mImg.getValue(0, 0));
        assertEquals(img.getValue(0, 1), mImg.getValue(0, 1));
        assertEquals(img.getValue(0, 2), mImg.getValue(0, 2));
        assertEquals(img.getValue(1, 0), mImg.getValue(1, 0));
        assertEquals(img.getValue(1, 1), mImg.getValue(1, 1));
        assertEquals(img.getValue(1, 2), mImg.getValue(1, 2));
        assertEquals(img.getValue(2, 0), mImg.getValue(2, 0));
        assertEquals(img.getValue(2, 1), mImg.getValue(2, 1));
        assertEquals(img.getValue(2, 2), mImg.getValue(2, 2));
    }
    
}
