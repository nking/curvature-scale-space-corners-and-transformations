package algorithms.imageProcessing;

import algorithms.util.PairInt;
import java.util.HashSet;
import java.util.Set;
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
    
    @Override
    public void applyFilter(GreyscaleImage input) {
                
        boolean hasABorderPixel = hasAtLeastOneBorderPixel(input);
        
        GreyscaleImage input2 = hasABorderPixel ? addOnePixelBorders(input) :
            input;
        
        int w2 = input2.getWidth();
        int h2 = input2.getHeight();
        
        Set<PairInt> points = new HashSet<PairInt>();
        for (int col = 0; col < w2; col++) {
            for (int row = 0; row < h2; row++) {
                if (input2.getValue(col, row) > 0) {
                    points.add(new PairInt(col, row));
                }
            }
        }
        applyLineThinner(points, 0, w2, 0, h2);
        input2.fill(0);
        for (PairInt p : points) {
            input2.setValue(p.getX(), p.getY(), 1);
        }

        PostLineThinnerCorrections pltc = new PostLineThinnerCorrections();
        if (edgeGuideImage != null) {
            pltc.setEdgeGuideImage(edgeGuideImage);
        }
        pltc.correctForArtifacts(input2);
                
        GreyscaleImage input3 = hasABorderPixel ? removeOnePixelBorders(input2)
            : input2;
        
        input.resetTo(input3);
                
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

}
