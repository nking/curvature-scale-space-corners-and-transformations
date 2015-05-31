package algorithms.imageProcessing;

import algorithms.CountingSort;
import algorithms.util.PairInt;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.logging.Logger;

/**
 *
 * @author nichole
 */
public class ZhangSuenLineThinner extends AbstractLineThinner {
    
    private static final int[][] nbrs = 
        {{0, -1}, {1, -1}, {1, 0}, {1, 1}, 
        {0, 1}, {-1, 1}, {-1, 0}, {-1, -1}, {0, -1}};
        
    private static int[][][] nbrGroups = {{{0, 2, 4}, {2, 4, 6}}, 
        {{0, 2, 6}, {0, 4, 6}}};

    protected boolean useLineDrawingMode = false;
    
    protected boolean debug = false;
    
    private Logger log = Logger.getLogger(this.getClass().getName());
    
    /**
     * for images which are already line drawings, that is images such as
     * maps with only lines, or for block images, use this to avoid a gap filling
     * stage that fills single pixel gaps surrounded by non-zero pixels.  
     * (Else, the filter applies such a gap filling algorithm to help avoid 
     * creating bubbles in thick lines).
     */
    @Override
    public void useLineDrawingMode() {
        useLineDrawingMode = true;
    }

    @Override
    public void setDebug(boolean setToDebug) {
        debug = setToDebug;
    }
    
    @Override
    public void applyFilter(final GreyscaleImage input) {
        
        GreyscaleImage summed = sumOver8Neighborhood(input);
        
        Set<PairInt> points = new HashSet<PairInt>();
        for (int col = 0; col < input.getWidth(); col++) {
            for (int row = 0; row < input.getHeight(); row++) {
                if (input.getValue(col, row) > 0) {
                    points.add(new PairInt(col, row));
                }
            }
        }
        applyLineThinner(points, 0, input.getWidth(),
            0, input.getHeight());
        input.fill(0);
        for (PairInt p : points) {
            input.setValue(p.getX(), p.getY(), 1);
        }
        
        // make corrections for artifacts created for inclined lines
        /*
        can see a 3x3 hollow square.  the correction should replace it
        with the diagonal of the exterior connecting slope.
        Is the artifact a product of the gaussian convolution size?
        If yes, would need to know the blur radius.
        */
        correctForArtifacts(input);
                
        //correctForMinorOffsetsByIntensity(input, summed);
    }
    
    public void applyLineThinner(Set<PairInt> points, int minX, int maxX,
        int minY, int maxY) {
        
        // adapted from code at http://rosettacode.org/wiki/Zhang-Suen_thinning_algorithm#Java
         
        boolean firstStep = false;
        boolean hasChanged;
        
        Set<PairInt> remove = new HashSet<PairInt>();
        
        do {
            hasChanged = false;
            firstStep = !firstStep;
             
            for (int r = minY + 1; r < maxY - 1; r++) {
                for (int c = minX + 1; c < maxX - 1; c++) {
                     
                    PairInt uPoint = new PairInt(c, r);
                    
                    if (!points.contains(uPoint)) {
                        continue;
                    }
 
                    int nn = numNeighbors(r, c, points);
                    if (nn < 2 || nn > 6) {
                        continue;
                    }
 
                    int nt = numTransitions(r, c, points);
                    if (nt != 1) {
                        continue;
                    }
 
                    if (!atLeastOneIsVacant(r, c, firstStep ? 0 : 1, points)) {
                        continue;
                    }
                     
                    remove.add(uPoint);
                }
            }
 
            if (!remove.isEmpty()) {
                
                for (PairInt p : remove) {
                    points.remove(p);
                    hasChanged = true;
                }

                remove.clear();
            }
 
        } while (hasChanged || firstStep);
    }
    
    private int numNeighbors(int r, int c, Set<PairInt> points) {
        int count = 0;
        for (int i = 0; i < nbrs.length - 1; i++) {
            int x = c + nbrs[i][0];
            int y = r + nbrs[i][1];
            PairInt p = new PairInt(x, y);
            if (points.contains(p)) {
                count++;
            } else {
                int z = 1;
            }
        }
        return count;
    }
 
    /**
     * visits neighbors in counter-clockwise direction and looks for the
     * pattern 0:1 in the current and next neighbor.  each such pattern
     * is counted as a transition.
     * @param r
     * @param c
     * @param points
     * @return 
     */
    static int numTransitions(int r, int c, Set<PairInt> points) {
        
        /* 
         5  4  3    1
         6     2    0
         7  0  1   -1
 
         -1  0  1
         */
        
        int count = 0;
        for (int i = 0; i < nbrs.length - 1; i++) {
            int x = c + nbrs[i][0];
            int y = r + nbrs[i][1];
            PairInt p = new PairInt(x, y);
            if (!points.contains(p)) {
                int x2 = c + nbrs[i + 1][0];
                int y2 = r + nbrs[i + 1][1];
                PairInt p2 = new PairInt(x2, y2);
                if (points.contains(p2)) {
                    count++;
                }
            }
        }
        return count;
    }
 
    /**
     * looking for 2 zeroes within the 4 neighborhood pattern of point (c,r).
     * @param r
     * @param c
     * @param step
     * @param points
     * @return 
     */
    static boolean atLeastOneIsVacant(int r, int c, int step, Set<PairInt> points) {
        int count = 0;
        
        int[][] group = nbrGroups[step];
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < group[i].length; j++) {
                int[] nbr = nbrs[group[i][j]];

                int x = c + nbr[0];
                int y = r + nbr[1];

                PairInt p = new PairInt(x, y);

                if (!points.contains(p)) {
                    count++;
                    break;
                }
            }
        }
        
        return (count > 1);
    }

    private void correctForArtifacts(GreyscaleImage input) {
                
        //correctForHoleArtifacts3(input);
        
        //correctForHoleArtifacts2(input);
        
        //correctForHoleArtifacts1(input);
                 
        correctForZigZag0(input);
                        
        correctForZigZag1(input);
        
        correctForZigZag2(input);
        
        correctForZigZag3(input);
        
        correctForZigZag4(input);
        
        correctForLine0(input);
        
    }
    
    private void correctForZigZag0(GreyscaleImage input) {
       
        /*
        looking for pattern
       
           0  0         2
           0  #  #      1
           #* #<=0      0
        #     0         -1
        
       -1  0  1  2
        
        and removing the topmost left #'s
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
       
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(0, -1));
        zeroes.add(new PairInt(0, -2));
        zeroes.add(new PairInt(1, 1));
        zeroes.add(new PairInt(1, -2));
        zeroes.add(new PairInt(2, 0));
        
        ones.add(new PairInt(-1, 1));
        ones.add(new PairInt(1, -1));
        ones.add(new PairInt(1, 0));
        ones.add(new PairInt(2, -1));
        
        changeToZeroes.add(new PairInt(1, 0));
        
        int startValue = 1;
        
        replacePattern(input, zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
        
        rotate90ThreeTimes(input, zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
        
    }
    
    private void correctForZigZag1(GreyscaleImage input) {
        
        /*
        looking for pattern
       
                 #      3
           0     #      2
        0  0  #  0      1
        0  #* #  0      0
        #               -1
        
       -1  0  1  2 
        
        and removing the topmost left #'s
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
        
        /*
        looking for pattern
       
                 #      3
           0     #      2
        0  0  #  0      1
        0  #* #  0      0
        #               -1
        
       -1  0  1  2 
        
        and removing the topmost left #'s
        */
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(-1, 0));
        zeroes.add(new PairInt(-1, -1));
        zeroes.add(new PairInt(0, -1));
        zeroes.add(new PairInt(0, -2));
        zeroes.add(new PairInt(2, 0));
        zeroes.add(new PairInt(2, -1));
        
        ones.add(new PairInt(-1, 1));
        ones.add(new PairInt(1, 0));
        ones.add(new PairInt(1, -1));
        ones.add(new PairInt(2, -2));
        ones.add(new PairInt(2, -3));
        
        changeToZeroes.add(new PairInt(1, 0));
        
        int startValue = 1;
        
        replacePattern(input, zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
        
        rotate90ThreeTimes(input, zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
        
    }
    
    private void correctForZigZag2(GreyscaleImage input) {
        
        /*
        looking for pattern
       
           0  0  0  #      1
           0  #* # 0       0
           0  #  0         -1
           #  0            -2
        
       -2 -1  0  1  2
        
        and removing the topmost left #'s
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
       
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(-1, 1));
        zeroes.add(new PairInt(-1, 0));
        zeroes.add(new PairInt(-1, -1));
        zeroes.add(new PairInt(0, 2));
        zeroes.add(new PairInt(0, -1));
        zeroes.add(new PairInt(1, 1));
        zeroes.add(new PairInt(1, -1));
        zeroes.add(new PairInt(2, 0));
      
        ones.add(new PairInt(-1, 2));
        ones.add(new PairInt(0, 1));
        ones.add(new PairInt(1, 0));
        ones.add(new PairInt(2, -1));
        
        changeToZeroes.add(new PairInt(0, 0));
        
        int startValue = 1;
        
        replacePattern(input, zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
        
        rotate90ThreeTimes(input, zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
        
    }
    
    private void correctForZigZag3(GreyscaleImage input) {
        
        /*
        looking for pattern
                 0  #      3
              0  #         2
           0  #  #  0      1
           0  #* 0         0
           #  0  0        -1
        
       -2 -1  0  1  2
        
        and removing the topmost left #'s
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
       
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(-1, 0));
        zeroes.add(new PairInt(-1, -1));
        zeroes.add(new PairInt(0, 1));
        zeroes.add(new PairInt(0, -2));
        zeroes.add(new PairInt(1, 1));
        zeroes.add(new PairInt(1, 0));
        zeroes.add(new PairInt(1, -3));
        zeroes.add(new PairInt(2, -1));
      
        ones.add(new PairInt(-1, 1));
        ones.add(new PairInt(0, -1));
        ones.add(new PairInt(1, -1));
        ones.add(new PairInt(1, -2));
        ones.add(new PairInt(2, -3));
        
        changeToZeroes.add(new PairInt(0, -1));
        
        int startValue = 1;
        
        replacePattern(input, zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
        
        rotate90ThreeTimes(input, zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
        
    }
    
    private void correctForZigZag4(GreyscaleImage input) {
        
        /*
        looking for pattern
                           3
        #  #  0  0         2
           0  #  0  0      1
           0  #* #  #      0
           0  0  0  0  #  -1
        
       -2 -1  0  1  2
        
        and removing the topmost left #'s
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
       
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(-1, 1));
        zeroes.add(new PairInt(-1, 0));
        zeroes.add(new PairInt(-1, -1));
        zeroes.add(new PairInt(0, 1));
        zeroes.add(new PairInt(0, -2));
        zeroes.add(new PairInt(1, 1));
        zeroes.add(new PairInt(1, -1));
        zeroes.add(new PairInt(1, -2));
        zeroes.add(new PairInt(2, 1));
        zeroes.add(new PairInt(2, -1));
      
        ones.add(new PairInt(-2, -2));
        ones.add(new PairInt(-1, -2));
        ones.add(new PairInt(0, -1));
        ones.add(new PairInt(1, 0));
        ones.add(new PairInt(2, 0));
        ones.add(new PairInt(3, 1));
        
        changeToZeroes.add(new PairInt(0, 0));
        
        int startValue = 1;
        
        replacePattern(input, zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
        
        rotate90ThreeTimes(input, zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
        
    }
    
    private void rotate90ThreeTimes(GreyscaleImage input, 
        final LinkedHashSet<PairInt> zeroes, final LinkedHashSet<PairInt> ones, 
        LinkedHashSet<PairInt> changeToZeroes, 
        final LinkedHashSet<PairInt> changeToOnes, final int startCenterValue) {
        
        // ----- change the sign of x to handle other direction -----
        for (PairInt p : zeroes) {
            p.setX(-1 * p.getX());
        }
        
        for (PairInt p : ones) {
            p.setX(-1 * p.getX());
        }
          
        for (PairInt p : changeToZeroes) {
            p.setX(-1 * p.getX());
        }
        
        for (PairInt p : changeToOnes) {
            p.setX(-1 * p.getX());
        }
                    
        replacePattern(input, zeroes, ones, changeToZeroes, changeToOnes,
            startCenterValue);
             
        // ----- change the sign of y to handle other direction -----
        for (PairInt p : zeroes) {
            p.setY(-1 * p.getY());
        }
        
        for (PairInt p : ones) {
            p.setY(-1 * p.getY());
        }
        
        for (PairInt p : changeToZeroes) {
            p.setY(-1 * p.getY());
        }
        
        for (PairInt p : changeToOnes) {
            p.setY(-1 * p.getY());
        }
        
        replacePattern(input, zeroes, ones, changeToZeroes, changeToOnes,
            startCenterValue);
        
        // ----- change the sign of x to handle another direction -----
        for (PairInt p : zeroes) {
            p.setX(-1 * p.getX());
        }
        
        for (PairInt p : ones) {
            p.setX(-1 * p.getX());
        }
        
        for (PairInt p : changeToZeroes) {
            p.setX(-1 * p.getX());
        }
        
        for (PairInt p : changeToOnes) {
            p.setX(-1 * p.getX());
        }
                    
        replacePattern(input, zeroes, ones, changeToZeroes, changeToOnes,
            startCenterValue);
    }
    
    private void replacePattern(GreyscaleImage input, 
        final LinkedHashSet<PairInt> zeroes, final LinkedHashSet<PairInt> ones, 
        final LinkedHashSet<PairInt> changeToZeroes, final LinkedHashSet<PairInt> changeToOnes, 
        final int startCenterValue) {
        
        int w = input.getWidth();
        int h = input.getHeight();
        
        for (int col = 0; col < w; col++) {
            
            for (int row = 0; row < h; row++) {
                
                int v = input.getValue(col, row);
              
                if (v != startCenterValue) {
                    continue;
                }
                
                boolean foundPattern = true;
                
                for (PairInt p : zeroes) {
                    int x = col + p.getX();
                    int y = row + p.getY();
                    if ((x < 0) || (y < 0) || (x > (w - 1)) || (y > (h - 1))) {
                        //TODO: revisit this
                        foundPattern = false;
                        break;
                    }
                    int vz = input.getValue(x, y);
                    if (vz != 0) {
                        foundPattern = false;
                        break;
                    }
                }
                
                if (!foundPattern) {
                    continue;
                }
                
                for (PairInt p : ones) {
                    int x = col + p.getX();
                    int y = row + p.getY();
                    if ((x < 0) || (y < 0) || (x > (w - 1)) || (y > (h - 1))) {
                        foundPattern = false;
                        break;
                    }
                    int vz = input.getValue(x, y);
                    if (vz != 1) {
                        foundPattern = false;
                        break;
                    }
                }
                
                if (!foundPattern) {
                    continue;
                }
                
                for (PairInt p : changeToZeroes) {
                    int x = col + p.getX();
                    int y = row + p.getY();
                    if ((x < 0) || (y < 0) || (x > (w - 1)) || (y > (h - 1))) {
                        continue;
                    }
                    input.setValue(x, y, 0);
                }
                
                for (PairInt p : changeToOnes) {
                    int x = col + p.getX();
                    int y = row + p.getY();
                    if ((x < 0) || (y < 0) || (x > (w - 1)) || (y > (h - 1))) {
                        continue;
                    }
                    input.setValue(x, y, 1);
                }
            }
        }
    }
    
    private void debugPrint(GreyscaleImage input, int xStart, int xStop,
        int yStart, int yStop) {
        
        StringBuilder sb = new StringBuilder();
                    
        for (int row = yStart; row <= yStop; row++) {
            sb.append(String.format("%3d: ", row));
            for (int col = xStart; col <= xStop; col++) {
                sb.append(String.format(" %3d ", input.getValue(col, row)));
            }
            sb.append(String.format("\n"));
        }
        
        System.out.println(sb.toString());
    }

    /**
     * removes a hole artifact in inclined lines.  note that this should
     * probably be adjusted for gaussian convolution combined radius
     * if used outside of the gradientXY image produced by the
     * CannyEdgeFilter.
     * @param input 
     */
    private void correctForHoleArtifacts3(GreyscaleImage input) {
        
        /*
        looking for pattern
        
            0    0    0    0    0    0           3
            0    0    0    0    1    1    1      2
            0    0    0    1    0    1    0      1
            0    0    1    0*   1    1    0      0
            0    1    0    1    1    0    0     -1
            0    1    1    1    0    0    0     -2
            0    1    0    0    0    0    0     -3
        
           -3   -2   -1    0    1    2    3     
        
        and removing the topmost left #'s
        
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
        
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(-1, -1));
        zeroes.add(new PairInt(-1, 1)); 
        zeroes.add(new PairInt(1, -1)); 
        zeroes.add(new PairInt(0, -2)); 
        zeroes.add(new PairInt(1, 2)); 
        zeroes.add(new PairInt(2, 2)); 
        zeroes.add(new PairInt(2, 1)); 
        zeroes.add(new PairInt(0, -3)); 
        zeroes.add(new PairInt(1, -3)); 
        zeroes.add(new PairInt(-1, -2)); 
        zeroes.add(new PairInt(-1, -3)); 
        zeroes.add(new PairInt(2, -3));
        zeroes.add(new PairInt(-2, -3));
        zeroes.add(new PairInt(-2, -2));
        zeroes.add(new PairInt(-2, -1));
        zeroes.add(new PairInt(-2, 0));
        zeroes.add(new PairInt(-1, 3));
        zeroes.add(new PairInt(0, 3));
        zeroes.add(new PairInt(1, 3));
        zeroes.add(new PairInt(2, 3));
        zeroes.add(new PairInt(3, 3));
        zeroes.add(new PairInt(3, 2));
        zeroes.add(new PairInt(3, 1));
        zeroes.add(new PairInt(3, 0));
        zeroes.add(new PairInt(3, -1));
        zeroes.add(new PairInt(-3, -3));
        zeroes.add(new PairInt(-3, 2));
        zeroes.add(new PairInt(-3, -1));
        zeroes.add(new PairInt(-3, 0));
        zeroes.add(new PairInt(-3, -1));
        zeroes.add(new PairInt(-3, -2));
        zeroes.add(new PairInt(-3, -3));
        
        ones.add(new PairInt(0, 1));
        ones.add(new PairInt(0, 2));
        ones.add(new PairInt(0, -1));
        ones.add(new PairInt(-1, 0));
        ones.add(new PairInt(-1, 2));
        ones.add(new PairInt(1, -2));
        ones.add(new PairInt(1, 0));
        ones.add(new PairInt(1, 1));
        ones.add(new PairInt(2, -2));
        ones.add(new PairInt(2, -1));
        ones.add(new PairInt(2, 0));
        ones.add(new PairInt(3, -2));
        ones.add(new PairInt(-2, 1));
        ones.add(new PairInt(-2, 2));
        ones.add(new PairInt(-2, 3));
    
        changeToZeroes.add(new PairInt(-2, 2));
        changeToZeroes.add(new PairInt(-2, 1));
        changeToZeroes.add(new PairInt(-1, 0));
        changeToZeroes.add(new PairInt(0, -1));
        changeToZeroes.add(new PairInt(1, -2));
        changeToZeroes.add(new PairInt(2, -2));
        changeToZeroes.add(new PairInt(2, 0));
        changeToZeroes.add(new PairInt(1, 1));
        changeToZeroes.add(new PairInt(0, 2));
              
        int centralValue = 0;
        
        replacePattern(input, zeroes, ones, changeToZeroes, changeToOnes, 
            centralValue);
        
        rotate90ThreeTimes(input, zeroes, ones, changeToZeroes, changeToOnes, 
            centralValue);
       
    }
    
    /**
     * removes a hole artifact in inclined lines.  note that this should
     * probably be adjusted for gaussian convolution combined radius
     * if used outside of the gradientXY image produced by the
     * CannyEdgeFilter.
     * @param input 
     */
    private void correctForHoleArtifacts2(GreyscaleImage input) {
        
        /*
        looking for pattern
        
         0    0    0    0    0    0     3
         0    0    0    1    1    1     2
         0    0    1    0    1    0     1
         0    1    0*   1    1    0     0
         0    1    1    1    0    0    -1
         0    1    0    0    0    0    -2
        
        -2   -1    0    1    2    3                
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
      
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(0, 2));
        zeroes.add(new PairInt(0, -2));
        zeroes.add(new PairInt(-1, -1));
        zeroes.add(new PairInt(-1, -2));
        zeroes.add(new PairInt(-1, -3));
        zeroes.add(new PairInt(0, -3));
        zeroes.add(new PairInt(1, -3));
        zeroes.add(new PairInt(1, -1));
        zeroes.add(new PairInt(1, 2));
        zeroes.add(new PairInt(2, -3));
        zeroes.add(new PairInt(2, 1));
        zeroes.add(new PairInt(2, 2));
        zeroes.add(new PairInt(-2, -3));
        zeroes.add(new PairInt(-2, -2));
        zeroes.add(new PairInt(-2, -1));
        zeroes.add(new PairInt(-2, 0));
        zeroes.add(new PairInt(-2, 1));
        zeroes.add(new PairInt(-2, 2));
        zeroes.add(new PairInt(3, -3));
        zeroes.add(new PairInt(3, -1));
        zeroes.add(new PairInt(3, 0));
        zeroes.add(new PairInt(3, 1));
        zeroes.add(new PairInt(3, 2));
        
        ones.add(new PairInt(0, -1));
        ones.add(new PairInt(0, 1));
        ones.add(new PairInt(-1, 0));
        ones.add(new PairInt(-1, 1));
        ones.add(new PairInt(-1, 2));
        ones.add(new PairInt(1, -2));
        ones.add(new PairInt(1, 0));
        ones.add(new PairInt(1, 1));
        ones.add(new PairInt(2, -2));
        ones.add(new PairInt(2, -1));
        ones.add(new PairInt(2, 0));
        ones.add(new PairInt(3, -2));
    
        changeToZeroes.add(new PairInt(-1, 0));
        changeToZeroes.add(new PairInt(-1, 1));
        changeToZeroes.add(new PairInt(0, -1));
        changeToZeroes.add(new PairInt(1, -2));
        changeToZeroes.add(new PairInt(1, 1));
        changeToZeroes.add(new PairInt(2, -2));
        changeToZeroes.add(new PairInt(2, 0));
        
        int centralValue = 0;
        
        replacePattern(input, zeroes, ones, changeToZeroes, changeToOnes, 
            centralValue);
        
        rotate90ThreeTimes(input, zeroes, ones, changeToZeroes, changeToOnes, 
            centralValue);
       
    }
    
    private void correctForLine0(GreyscaleImage input) {
        
        /*
        looking for pattern
       
           #            2
           #  0         1
        0  0* #  0      0
           #  0        -1
           #
        
       -1  0  1  2
        
        and removing the topmost left #'s
        */
        
        LinkedHashSet<PairInt> ones = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> zeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToZeroes = new LinkedHashSet<PairInt>();
        LinkedHashSet<PairInt> changeToOnes = new LinkedHashSet<PairInt>();
       
        /*
        looking for pattern
       
           #            2
           #  0         1
        0  0* #  0      0
           #  0        -1
           #           -2
        
       -1  0  1  2
        
        and removing the topmost left #'s
        */
        
        // y's are inverted here because sketch above is top left is (0,0)
        zeroes.add(new PairInt(-1, 0));
        zeroes.add(new PairInt(1, 1));
        zeroes.add(new PairInt(1, -1));
        zeroes.add(new PairInt(2, 0));
        
        ones.add(new PairInt(0, 2));
        ones.add(new PairInt(0, 1));
        ones.add(new PairInt(0, -1));
        ones.add(new PairInt(0, -2));
        ones.add(new PairInt(1, 0));
        
        changeToZeroes.add(new PairInt(1, 0));
        changeToOnes.add(new PairInt(0, 0));
        
        int startValue = 0;
        
        replacePattern(input, zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
        
        rotate90ThreeTimes(input, zeroes, ones, changeToZeroes, changeToOnes, 
            startValue);
        
    }
    
    private void correctForHoleArtifacts1(GreyscaleImage input) {
        
        /*     
        
        look for pattern with hole in middle,
        fill the hole,
        then use 
        boolean nullable = erosionFilter.process(input, input, 
            neighborX, neighborY);
                                          
                      1               1
                 1    0*   1          0     
                      1              -1
                                     -2
        
           -2   -1    0    1    2
        */  
        
        ErosionFilter erosionFilter = new ErosionFilter();
        
        int w = input.getWidth();
        int h = input.getHeight();
      
        int[] nbX = new int[]{-1, -1, 0, 1, 1, 1,  0, -1};
        int[] nbY = new int[]{ 0,  1, 1, 1, 0, -1,-1, -1};

        Set<PairInt> ones = new HashSet<PairInt>();
       
        // y's are inverted here because sketch above is top left is (0,0)
        
        ones.add(new PairInt(-1, 0));
        ones.add(new PairInt(0, -1));
        ones.add(new PairInt(0, 1));
        ones.add(new PairInt(1, 0));
        
        int centralValue = 0;
            
        for (int col = 0; col < w; col++) {

            for (int row = 0; row < h; row++) {

                int v = input.getValue(col, row);

                if (v != centralValue) {
                    continue;
                }

                boolean foundPattern = true;

                for (PairInt p : ones) {
                    int x = col + p.getX();
                    int y = row + p.getY();
                    if ((x < 0) || (y < 0) || (x > (w - 1)) || (y > (h - 1))) {
                        continue;
                    }
                    int vz = input.getValue(x, y);
                    if (vz != 1) {
                        foundPattern = false;
                        break;
                    }
                }

                if (!foundPattern) {
                    continue;
                }

                // set the center pixel to '1' and visit each in 8 neighbor
                // hood to determine if can null it

                input.setValue(col, row, 1);

                for (int k = 0; k < nbX.length; k++) {
                    int x = col + nbX[k];
                    int y = row + nbY[k];
                    if ((x < 0) || (y < 0) || (x > (w - 1)) || (y > (h - 1))) {
                        continue;
                    }

                    boolean nullable = erosionFilter.process(input, input, 
                        x, y);

                    if (nullable) {
                        input.setValue(x, y, 0);
                    }
                }
            }
        }
        
    }

    protected GreyscaleImage sumOver8Neighborhood(GreyscaleImage img) {
        
        GreyscaleImage summed = img.copyImage();
        
        int[] dxs = new int[]{-1, -1,  0,  1, 1, 1, 0, -1};
        int[] dys = new int[]{ 0, -1, -1, -1, 0, 1, 1,  1};
        
        int w = img.getWidth();
        int h = img.getHeight();
        
        // for each pixel, sum it's neighbors
        for (int col = 0; col < w; col++) {
            for (int row = 0; row < h; row++) {
                
                int sum = 0;
                
                for (int nIdx = 0; nIdx < dxs.length; nIdx++) {
                    
                    int x = dxs[nIdx] + col;
                    int y = dys[nIdx] + row;
                    
                    if ((x<0) || (y<0) || (x>(w-1)) || (y>(h-1))) {
                        continue;
                    }
                    int v = img.getValue(x, y);
                    
                    sum += v;                    
                }
                summed.setValue(col, row, sum);
            }
        }
        
        return summed;
    }

    private void correctForMinorOffsetsByIntensity(GreyscaleImage input, 
        GreyscaleImage summed) {
        
        /*
        for each pixel in input, put the 8 neighbors in sorted order,
        keeping track of their location.
        iterate over the list of descending sorted sums:
            while centerSum < neighborSum
                if (neighbor value in input is 0) 
                    if swap does not disconnect any lines
                        swap pixels
                        break
        */
        
        int[] dxs = new int[]{-1, -1,  0,  1, 1, 1, 0, -1};
        int[] dys = new int[]{ 0, -1, -1, -1, 0, 1, 1,  1};
       
        ErosionFilter erosionFilter = new ErosionFilter();
        
        int w = input.getWidth();
        int h = input.getHeight();
            
        for (int col = 1; col < (w - 1); col++) {
            for (int row = 1; row < (h - 1); row++) {
                
                int v = input.getValue(col, row);
                if (v == 0) {
                    continue;
                }
                
                int vSum = summed.getValue(col, row);
                
                if (vSum == 8) {
                    continue;
                }
                                
                if (erosionFilter.doesDisconnect(input, col, row)) {
                    continue;
                }
                
                int maxSum = vSum;
                int maxIdx = -1;
                
                for (int nIdx = 0; nIdx < dxs.length; nIdx++) {
                    
                    int x = dxs[nIdx] + col;
                    int y = dys[nIdx] + row;
                    
                    if ((x<0) || (y<0) || (x>(w-1)) || (y>(h-1))) {
                        continue;
                    }
                    
                    // only compare the neighbors which are swappable
                    if (input.getValue(x, y) == 0) {
                        int sum = summed.getValue(x, y);
                        if (sum > maxSum) {
                            maxIdx = nIdx;
                            maxSum = sum;
                        }
                    }
                }
                                
                if (maxIdx > -1) {
                    int x = dxs[maxIdx] + col;
                    int y = dys[maxIdx] + row;
                    input.setValue(x, y, 1);
                    input.setValue(col, row, 0);
                }
                
            }    
        }
    }
    
}
