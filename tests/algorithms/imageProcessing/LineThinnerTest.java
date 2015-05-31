package algorithms.imageProcessing;

import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class LineThinnerTest {
    
    public LineThinnerTest() {
    }
    
    @Before
    public void setUp() {
    }
    
    @After
    public void tearDown() {
    }

    @Test
    public void testApplyFilter() {
        
        /*                                    ideally
          9 _ _ _ _ _ _ _ _ _ _ _ _
          8   @ @ @ @ @ @ @ @ @ @ _
          7   @ @ @ @ @ @ @ @ @ @ _         @ @ @ @ @ @ @ @
          6   @ @ @ @ @ @ @ @ @ @ _         @             @
          5   @ @ @         @ @ @ _         @             @
          4   @ @ @         @ @ @ _         @             @
          3   @ @ @ @ @ @ @ @ @ @ _         @             @
          2   @ @ @ @ @ @ @ @ @ @ _         @ @ @ @ @ @ @ @
          1   @ @ @ @ @ @ @ @ @ @ _         
          0                       _
            0 1 2 3 4 5 6 7 8 9         0 1 2 3 4 5 6 7 8 9  
        
        */
               
        GreyscaleImage input = getTestRectangle();
        
        GreyscaleImage input2 = input.copyImage();
                
        ErosionFilter instance = new ErosionFilter();
        
        instance.applyFilter(input);
        
        //System.out.println("erosion:");
        //printImage(input);
            
        GreyscaleImage input3 = input2.copyImage();
        
        ZhangSuenLineThinner instance2 = new ZhangSuenLineThinner();
        
        instance2.applyFilter(input2);
        
        //System.out.println("zhang-suen:");
        //printImage(input2);
        
        Set<PairInt> expected = getExpectedThinnedTestRectangle();
        int nExpectedFound = 0;
        int nNotExpectedFound = 0;
        for (int col = 0; col < input2.getWidth(); col++) {
            for (int row = 0; row < input2.getHeight(); row++) {
                
                int v = input2.getValue(col, row);
                PairInt p = new PairInt(col, row);
                
                if (v == 1) {
                    if (expected.contains(p)) {
                        nExpectedFound++;
                    } else {
                        nNotExpectedFound++;
                    }
                }
            }
        }
        
        // there are small errors, so need a tolerance
        float eps = 0.1f * expected.size();
        assertTrue(Math.abs(nExpectedFound - expected.size()) < eps);
        assertTrue(nNotExpectedFound < eps);
        
    }
    
    @Test
    public void testApplyFilter_shell() {
        
        GreyscaleImage input = getTestShell();
        
        GreyscaleImage input2 = input.copyImage();
        
        GreyscaleImage input3 = input.copyImage();
        
        ErosionFilter instance = new ErosionFilter();
        
        instance.applyFilter(input);
        
        //System.out.println("erosion:");
        //printImage(input);
       
        //-----------------------------------------------------
        ZhangSuenLineThinner instance2 = new ZhangSuenLineThinner();
        
        instance2.applyFilter(input2);
        
        //System.out.println("zhang-suen:");
        //printImage(input2);
        
        //System.out.println("summed:");
        //printSummed(input3);
        
        Set<PairInt> expected = getExpectedThinnedTestShell();
        int nExpectedFound = 0;
        int nNotExpectedFound = 0;
        for (int col = 0; col < input2.getWidth(); col++) {
            for (int row = 0; row < input2.getHeight(); row++) {
                
                int v = input2.getValue(col, row);
                PairInt p = new PairInt(col, row);
                
                if (v == 1) {
                    if (expected.contains(p)) {
                        nExpectedFound++;
                    } else {
                        nNotExpectedFound++;
                    }
                }
            }
        }
        
        // there are small errors, so need a tolerance
        float eps = 0.1f * expected.size();
        assertTrue(Math.abs(nExpectedFound - expected.size()) < eps);
        assertTrue(nNotExpectedFound < eps);
    }
    
    private void printSummed(GreyscaleImage img) {
        
        ZhangSuenLineThinner lineThinner = new ZhangSuenLineThinner();
        
        GreyscaleImage input4 = lineThinner.sumOver8Neighborhood(img);
        
        printImage(input4);
    }
    
    private GreyscaleImage getTestRectangleLeftUpper() {
        
        GreyscaleImage input = getTestImage00();
        
        for (int row = 0; row < 5; row++) {
            for (int col = 0; col < input.getWidth(); col++) {
                input.setValue(col, row, 0);
            }
        }
        
        printImage(input);
        
        return input;
    }
    
    private GreyscaleImage getTestRectangleLeftBottom() {
        
        GreyscaleImage input = getTestRectangleLeft();
        
        for (int row = 5; row < input.getHeight(); row++) {
            for (int col = 0; col < input.getWidth(); col++) {
                input.setValue(col, row, 0);
            }
        }
        
        printImage(input);
        
        return input;
    }
    
    private GreyscaleImage getTestRectangleLeft() {
        
        GreyscaleImage input = getTestImage00();
        
        for (int col = 6; col < input.getWidth(); col++) {
            for (int row = 0; row < input.getHeight(); row++) {
                input.setValue(col, row, 0);
            }
        }
       
        //printImage(input);
        
        return input;
    }
    
    private GreyscaleImage getTestRectangleRightBottom() {
        
        GreyscaleImage input = getTestRectangleRight();
        
        for (int row = 5; row < input.getHeight(); row++) {
            for (int col = 0; col < input.getWidth(); col++) {
                input.setValue(col, row, 0);
            }
        }
        
        printImage(input);
        
        return input;
    }
    
    private GreyscaleImage getTestRectangleRightUpper() {
        
        GreyscaleImage input = getTestRectangleRight();
        
        for (int row = 0; row < 5; row++) {
            for (int col = 0; col < input.getWidth(); col++) {
                input.setValue(col, row, 0);
            }
        }
        
        printImage(input);
        
        return input;
    }
    
    private GreyscaleImage getTestRectangleRight() {
        
        GreyscaleImage input = getTestImage00();
        
        for (int col = 0; col < 6; col++) {
            for (int row = 0; row < input.getHeight(); row++) {
                input.setValue(col, row, 0);
            }
        }
       
        //printImage(input);
        
        return input;
    }
    
    private void printImage(final GreyscaleImage input) {
        for (int row = (input.getHeight() - 1); row > -1; row--) {
            StringBuilder sb = new StringBuilder(String.format("row %2d:  ", row));
            for (int col = 0; col < input.getWidth(); col++) {
                int v = input.getValue(col, row);
                String str = (v == 0) ? String.format("  ") : String.format("%d ", v);
                sb.append(str);
            }
            System.out.println(sb.toString());
        }
        StringBuilder sb = new StringBuilder(String.format("        "));
        for (int col = 0; col < input.getWidth(); col++) {
            sb.append(String.format("%2d", col));
        }
        System.out.println(sb.toString());
        System.out.println("\n");
    }
    
    private Set<PairInt> getExpectedThinnedTestRectangle() {
        
        Set<PairInt> expected = new HashSet<PairInt>();
        
        for (int row = 2; row <= 7; row++) {
            expected.add(new PairInt(2, row));
            expected.add(new PairInt(9, row));
        }
        for (int col = 2; col <= 9; col++) {
            expected.add(new PairInt(col, 2));
            expected.add(new PairInt(col, 7));
        }
        
        return expected;
    }
    
    private GreyscaleImage getTestRectangle() {
        
        GreyscaleImage input = getTestImage00();
                        
        //printImage(input);
        
        return input;
    }
    
    private Set<PairInt> getExpectedThinnedTestShell() {
    
        Set<PairInt> expected = new HashSet<PairInt>();
        
        List<String> str = new ArrayList<String>();
        str.add("               11111             ");
        str.add("           11111   11111         ");
        str.add("        1111           1111      ");
        str.add("       1                   1     ");
        str.add("      1                     1    ");
        str.add("     1                       1   ");
        str.add("    1                         1  ");
        str.add("    1                         1  ");
        str.add("    1                         1  ");
        str.add("     1                       1   ");
        str.add("      1                     1    ");
        str.add("       1                   1     ");
        str.add("        1111          1111       ");
        str.add("           11111  11111          ");
        str.add("               1111              ");
        
        int w = str.get(0).length();
        int h = str.size();
        
        //GreyscaleImage img = new GreyscaleImage(w, h);
        
        for (int row = 0; row < str.size(); row++) {
            String rowStr = str.get(row);
            for (int col = 0; col < w; col++) {
                char c = rowStr.charAt(col);
                if (c != ' ') {
                    expected.add(new PairInt(col, row));
                    //img.setValue(col, row, 1);
                }
            }
        }
        
        //printImage(img);
        
        return expected;
    }
    
    private GreyscaleImage getTestShell() {
            
        /*
        http://ascii.co.uk/art/circle
                        11111
                    1111111111111
                 1111           1111
               111                 111
              111                   111
             111                     111
            111                       111
            111                       111
            111                       111
             111                     111
              111                   111
               111                 111
                 1111          1111
                    111111111111
                        1111        
        */
     
        List<String> str = new ArrayList<String>();
        str.add("               11111             ");
        str.add("           1111111111111         ");
        str.add("        1111           1111      ");
        str.add("      111                 111    ");
        str.add("     111                   111   ");
        str.add("    111                     111  ");
        str.add("   111                       111 ");
        str.add("   111                       111 ");
        str.add("   111                       111 ");
        str.add("    111                     111  ");
        str.add("     111                   111   ");
        str.add("      111                 111    ");
        str.add("        1111          1111       ");
        str.add("           111111111111          ");
        str.add("               1111              ");
        
        int w = str.get(0).length();
        int h = str.size();
        
        GreyscaleImage img = new GreyscaleImage(w, h);
        for (int row = 0; row < str.size(); row++) {
            String rowStr = str.get(row);
            for (int col = 0; col < w; col++) {
                char c = rowStr.charAt(col);
                if (c != ' ') {
                    img.setValue(col, row, 1);
                }
            }
        }
        
        //printImage(img);
        
        return img;
    }
    
    /*http://ascii.co.uk/art/circle
                   ooo OOO OOO ooo
               oOO                 OOo
           oOO                         OOo
        oOO                               OOo
      oOO                                   OOo
    oOO                                       OOo
   oOO                                         OOo
  oOO                                           OOo
 oOO                                             OOo
 oOO                                             OOo
 oOO                                             OOo
 oOO                                             OOo
 oOO                                             OOo
  oOO                                           OOo
   oOO                                         OOo
    oOO                                       OOo
      oOO                                   OOo
        oO                                OOo
           oOO                         OOo
               oOO                 OOo
                   ooo OOO OOO ooo

    */
    private GreyscaleImage getTestImage00() {
        /*
          9 _ _ _ _ _ _ _ _ _ _ _ _
          8   @ @ @ @ @ @ @ @ @ @ _
          7   @ @ @ @ @ @ @ @ @ @ _   
          6   @ @ @ @ @ @ @ @ @ @ _   
          5   @ @ @         @ @ @ _   
          4   @ @ @         @ @ @ _   
          3   @ @ @ @ @ @ @ @ @ @ _ 
          2   @ @ @ @ @ @ @ @ @ @ _
          1   @ @ @ @ @ @ @ @ @ @ _
          0                       _
            0 1 2 3 4 5 6 7 8 9
        */
        GreyscaleImage input = new GreyscaleImage(12, 10);
        for (int col = 1; col < 11; col++) {
            input.setValue(col, 1, 1);
        }
        for (int col = 1; col < 11; col++) {
            input.setValue(col, 2, 1);
        }
        for (int col = 1; col < 11; col++) {
            input.setValue(col, 3, 1);
        }
        input.setValue(1, 4, 1);input.setValue(2, 4, 1);input.setValue(3, 4, 1);
        input.setValue(8, 4, 1);input.setValue(9, 4, 1);input.setValue(10, 4, 1);
        
        input.setValue(1, 5, 1);input.setValue(2, 5, 1);input.setValue(3, 5, 1);
        input.setValue(8, 5, 1);input.setValue(9, 5, 1);input.setValue(10, 5, 1);
        
        for (int col = 1; col < 11; col++) {
            input.setValue(col, 6, 1);
        }
        for (int col = 1; col < 11; col++) {
            input.setValue(col, 7, 1);
        }
        for (int col = 1; col < 11; col++) {
            input.setValue(col, 8, 1);
        }
      
        return input;
    }
    
    private GreyscaleImage getTestCircleImage() {
        /*
         0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0     4
         0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0     3
         0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0     2
         0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0     1
         0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0    20
         0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0     9
         0 0 0 0 1 1 1 1 0 0 0 0 0 0 0 0 1 1 1 1 0 0 0 0 0     8
         0 0 0 1 1 1 1 0 0 0 0 0 0 0 0 0 0 1 1 1 1 0 0 0 0     7
         0 0 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 0 0 0 0     6
         0 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 0 0 0     5
         0 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 0 0 0     4
         0 0 1 1 1 0 0 0 0 0 0       0 0 0 0 0 1 1 1 0 0 0     3
         0 0 1 1 1 0 0 0 0 0 0   C   0 0 0 0 0 1 1 1 1 0 0  *  2
         0 0 1 1 1 0 0 0 0 0 0       0 0 0 0 0 1 1 1 0 0 0     1
         0 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 0 0 0    10
         0 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 0 0 0     9
         0 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 0 0 0     8
         0 0 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 0 0 0 0     7
         0 0 0 1 1 1 1 0 0 0 0 0 0 0 0 0 0 1 1 1 1 0 0 0 0     6
         0 0 0 0 1 1 1 1 0 0 0 0 0 0 0 0 1 1 1 1 0 0 0 0 0     5
         0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0     4
         0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0     3
         0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0     2
         0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0     1
         0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0     0
                                 *
         0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4
                             1                   2
          
        EXPECTED when different are '+'.  spot checks only.
        So, can see that when a line thickness is an even number,
        the erosion filter alone doesn't know to make corrections
        for a larger unknown shape.  It's not ideal, but it's fast
        and the difference from where one would like the remaining
        pixel to have been is at most += 1 pixel due to the
        even number of pixels.
        
         0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0   4
         0 0 0 0 0 0 0 0 0                   0 0 0 0 0 0 0   3
         0 0 0 0 0 0 0 0                       0 0 0 0 0 0   2
         0 0 0 0 0 0 0                 1         0 0 0 0 0   1
         0 0 0 0 0 0     1 1 1 1 1 1 1 + 1 1       0 0 0 0  20 
         0 0 0 0 0     1                     1       0 0 0   9
         0 0 0 0     1                         1       0 0   8
         0 0 0       1                           1     0 0   7
         0 0       1                           + 1     0 0   6
         0 0     1                               + 1   0 0   5
         0 0   1                                 1     0 0   4   
         0 0   1                                 1     0 0   3
         0 0   1                 C               1     0 0 * 2
         0 0   1                                 1     0 0   1
         0 0   1                                 1     0 0  10
         0 0   1                                 1     0 0   9
         0 0     1                               1     0 0   8
         0 0 0     1                           1     0 0 0   7
         0 0 0 0     1                     1 1       0 0 0   6
         0 0 0 0 0     1                 1       0 0 0 0 0   5
         0 0 0 0 0 0     1             1       0 0 0 0 0 0   4
         0 0 0 0 0 0 0     1 1 1 1 1 1       0 0 0 0 0 0 0   3
         0 0 0 0 0 0 0 0                   0 0 0 0 0 0 0 0   2
         0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0   1
         0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0   0
                                 *
         0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4
                             1                   2
        
        */
        int rMax = 10;
        GreyscaleImage input = new GreyscaleImage((rMax<<1) + 5, (rMax<<1) + 5);
        
        float r1 = rMax-2;
        float r2 = rMax-1;
        float r3 = rMax;
        int n = 80;
        int xc = rMax + 2;
        int yc = xc;
                        
        for (int i = 0; i < n; i++) {
            /*
            (x-xc)^2 + (y-yc)^2 = r
            x = xc + r*cos(theta)
            y = yc + r*sin(theta)
            */
            double thetaRadians = ((double)i) * 2.*Math.PI/((double)n);
            
            double c = Math.cos(thetaRadians);
            double s = Math.sin(thetaRadians);
            
            input.setValue((int) (xc + (r1 * c)), (int) (yc + (r1 * s)), 1);
            input.setValue((int) (xc + (r2 * c)), (int) (yc + (r2 * s)), 1);
            input.setValue((int) (xc + (r3 * c)), (int) (yc + (r3 * s)), 1);
        }
        
        // touchups
        input.setValue(xc, 1, 1);
        input.setValue(1, yc, 1);
        
        System.out.println("circle:");
        //printImage(input);
        
        return input;
    }
    
    private PairIntArray getTestCircleExpectedErosion() {
       
        PairIntArray xy = new PairIntArray();
        
        int rMax = 10;
        
        float r2 = rMax-1;
        int n = 80;
        int xc = rMax + 2;
        int yc = xc;
                        
        for (int i = 0; i < n; i++) {
            /*
            (x-xc)^2 + (y-yc)^2 = r
            x = xc + r*cos(theta)
            y = yc + r*sin(theta)
            */
            double thetaRadians = ((double)i) * 2.*Math.PI/((double)n);
            
            double c = Math.cos(thetaRadians);
            double s = Math.sin(thetaRadians);
            
            int xt = (int) Math.round(xc + (r2 * c));
            int yt = (int) Math.round(yc + (r2 * s)); 
            
            xy.add(xt, yt);
        }
        
        return xy;
    }

    private GreyscaleImage getTestImage0() {
         /*                                   ideally
          7                           7    
          6   @ @ @ @ @ @ @ @ @ @     6
          5   @ @ @ @ @ @ @ @ @ @     5   @ @ @ @ @ @ @ @ 
          4   @ @ @ @ @ @ @ @ @ @     4   @             @
          3   @ @ @         @ @ @     3   @             @  
          2   @ @ @         @ @ @     2   @             @  
          1   @ @ @ @ @ @ @ @ @ @     1   @ @ @ @ @ @ @ @
          0                           0
            0 1 2 3 4 5 6 7 8 9       0 1 2 3 4 5 6 7 8 9
        */
        /*      result
         6   
         5         1 1 1 1 1 1 
         4       1             1   
         3       1             1   
         2         1         1   
         1           1 1 1 1     
         0   
        
             0 1 2 3 4 5 6 7 8 9 0 1 
        
        */
        GreyscaleImage input = new GreyscaleImage(12, 8);
        for (int col = 1; col < 11; col++) {
            input.setValue(col, 1, 1);
        }
        input.setValue(1, 2, 1);input.setValue(2, 2, 1);input.setValue(3, 2, 1);
        input.setValue(8, 2, 1);input.setValue(9, 2, 1);input.setValue(10, 2, 1);
        
        input.setValue(1, 3, 1);input.setValue(2, 3, 1);input.setValue(3, 3, 1);
        input.setValue(8, 3, 1);input.setValue(9, 3, 1);input.setValue(10, 3, 1);
        
        for (int col = 1; col < 11; col++) {
            input.setValue(col, 4, 1);
        }
        for (int col = 1; col < 11; col++) {
            input.setValue(col, 5, 1);
        }
        for (int col = 1; col < 11; col++) {
            input.setValue(col, 6, 1);
        }
        
        for (int row = (input.getHeight() - 1); row > -1; row--) {
            StringBuilder sb = new StringBuilder();
            for (int col = 0; col < input.getWidth(); col++) {
                sb.append(input.getValue(col, row)).append(" ");
            }
            System.out.println(sb.toString());
        }
        System.out.println("\n");
        
        return input;
    }

    @Test
    public void testPrefillGapsSurroundedByNeighbors() {
        /*2   @
          1 @ C @
          0   @
            0 1 2
        */
        GreyscaleImage input = new GreyscaleImage(4, 4);
        input.setValue(0, 1, 1);
        input.setValue(1, 2, 1);
        input.setValue(2, 1, 1);
        input.setValue(1, 0, 1);
        ErosionFilter instance = new ErosionFilter();
        instance.prefillGapsSurroundedByNeighbors(input);
        assertTrue(input.getValue(1, 1) == 1);
    }

    @Test
    public void testHasImmediateFourNeighbors() {
        
        /*2   @
          1 @ C @
          0   @
            0 1 2
        */
        
        GreyscaleImage input = new GreyscaleImage(4, 4);
        input.setValue(1, 0, 1);
        input.setValue(0, 1, 1);
        input.setValue(1, 1, 1);
        input.setValue(2, 1, 1);
        input.setValue(1, 2, 1);
        int col = 1;
        int row = 1;
        ErosionFilter instance = new ErosionFilter();
        boolean expResult = true;
        boolean result = instance.hasImmediateFourNeighbors(input, col, row);
        assertEquals(expResult, result);
        
        /*2   @
          1   C @
          0   @
            0 1 2
        */
        input = new GreyscaleImage(4, 4);
        input.setValue(1, 0, 1);
        input.setValue(1, 1, 1);
        input.setValue(2, 1, 1);
        input.setValue(1, 2, 1);
        col = 1;
        row = 1;
        instance = new ErosionFilter();
        expResult = false;
        result = instance.hasImmediateFourNeighbors(input, col, row);
        assertEquals(expResult, result);
    }

    @Test
    public void testIsTheEndOfAOnePixelLine() {
        /*
            #
          . @ .
          . C .
          . . .
        */
        GreyscaleImage input = new GreyscaleImage(4, 4);
        input.setValue(1, 1, 1);
        input.setValue(1, 2, 1);
        input.setValue(1, 3, 1);
        int col = 1;
        int row = 1;
        ErosionFilter instance = new ErosionFilter();
        boolean expResult = true;
        boolean result = instance.isTheEndOfAOnePixelLine(input, col, row);
        assertEquals(expResult, result);
        
        input.setValue(0, 2, 1);
        input.setValue(2, 2, 1);
        instance = new ErosionFilter();
        expResult = false;
        result = instance.isTheEndOfAOnePixelLine(input, col, row);
        assertEquals(expResult, result);
    }
    
}
