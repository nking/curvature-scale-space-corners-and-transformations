package algorithms.imageProcessing;

import algorithms.util.PairInt;
import java.util.HashSet;
import java.util.Set;

/**
 *
 * @author nichole
 */
public class ZhangSuenLineThinner {
    
    private static final int[][] nbrs = {{0, -1}, {1, -1}, {1, 0}, {1, 1}, 
        {0, 1}, {-1, 1}, {-1, 0}, {-1, -1}, {0, -1}};
        
    private static final int[][][] nbrGroups = {{{0, 2, 4}, {2, 4, 6}}, 
        {{0, 2, 6}, {0, 4, 6}}};
            
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
 
                    if (numTransitions(r, c, points) != 1) {
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
            }
        }
        return count;
    }
 
    static int numTransitions(int r, int c, Set<PairInt> points) {
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
        return count > 1;
    }
}
