package algorithms.imageProcessing;

import algorithms.CountingSort;
import algorithms.util.PairInt;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
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
    
    private static final int[][] nbrs = {{0, -1}, {1, -1}, {1, 0}, {1, 1}, 
        {0, 1}, {-1, 1}, {-1, 0}, {-1, -1}, {0, -1}};
        
    private static final int[][][] nbrGroups = {{{0, 2, 4}, {2, 4, 6}}, 
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

    private void correctForArtifacts(GreyscaleImage input) {
        
        correctForHoleArtifacts(input);
    }
    
    private void correctForHoleArtifacts(GreyscaleImage input) {
        /*
        can see a 3x3 hollow square.  the correction should replace it
        with the diagonal of the exterior connecting slope.
        Is the artifact a product of the gaussian convolution size?
        If yes, would need to know the blur radius.
        */
        
        // key = location number, value = offsets from center
        Map<Integer, List<PairInt>> extSlopeLocations = 
            new HashMap<Integer, List<PairInt>>();
        /*
         66700
         6###0
         5# #1
         4###2
         44322   count pixels in spaces 0 thru 7
         */
        for (int i = 0; i < 8; i++) {
             extSlopeLocations.put(Integer.valueOf(i), new ArrayList<PairInt>());
        }
        extSlopeLocations.get(Integer.valueOf(0)).add(new PairInt(1, 2));
        extSlopeLocations.get(Integer.valueOf(0)).add(new PairInt(2, 2));
        extSlopeLocations.get(Integer.valueOf(0)).add(new PairInt(2, 1));
        extSlopeLocations.get(Integer.valueOf(1)).add(new PairInt(2, 0));
        extSlopeLocations.get(Integer.valueOf(2)).add(new PairInt(2, -1));
        extSlopeLocations.get(Integer.valueOf(2)).add(new PairInt(2, -2));
        extSlopeLocations.get(Integer.valueOf(2)).add(new PairInt(1, -2));
        extSlopeLocations.get(Integer.valueOf(3)).add(new PairInt(0, -2));
        extSlopeLocations.get(Integer.valueOf(4)).add(new PairInt(-1, -2));
        extSlopeLocations.get(Integer.valueOf(4)).add(new PairInt(-2, -2));
        extSlopeLocations.get(Integer.valueOf(4)).add(new PairInt(-2, -1));
        extSlopeLocations.get(Integer.valueOf(5)).add(new PairInt(-2, 0));
        extSlopeLocations.get(Integer.valueOf(6)).add(new PairInt(-2, 1));
        extSlopeLocations.get(Integer.valueOf(6)).add(new PairInt(-2, 2));
        extSlopeLocations.get(Integer.valueOf(6)).add(new PairInt(-1, 2));
        extSlopeLocations.get(Integer.valueOf(7)).add(new PairInt(0, 2));
        
        int[] locationNumbers = new int[8];
        int[] counts = new int[8];
        
        int radius = 3;
        // looking for 3x3 hollow squares
        for (int col = (0 + (radius/2)); col < (input.getWidth() - (radius/2)); 
            col++) {
            
            for (int row = (0 + (radius/2)); row < (input.getHeight() - 
                (radius/2)); row++) {
                
                int v = input.getValue(col, row);
                
                if (v != 0) {
                    continue;
                }
                
                boolean foundPattern = true;
                
                for (int c0 = (col - 1); c0 <= (col + 1); c0++) {
                    for (int r0 = (row - 1); r0 <= (row + 1); r0++) {
                        if (c0 == col && r0 == row) {
                            continue;
                        }
                        if (input.getValue(c0, r0) == 0) {
                            foundPattern = false;
                            break;
                        }
                    }
                    if (!foundPattern) {
                        break;
                    }
                }
                if (foundPattern) {
                    // null the pixels in the box pattern:
                    for (int c0 = (col - 1); c0 <= (col + 1); c0++) {
                        for (int r0 = (row - 1); r0 <= (row + 1); r0++) {
                            if (c0 == col && r0 == row) {
                                input.setValue(c0, r0, 1);
                            } else {
                                input.setValue(c0, r0, 0);
                            }
                        }
                    }
                    
                    // determine slope of exterior lines
                    /*
                      66700
                      6###0
                      5# #1
                      4###2
                      44322   count pixels in spaces 0 thru 7
                    */
                    for (int i = 0; i < locationNumbers.length; i++) {
                        locationNumbers[i] = i;
                        counts[i] = 0;
                    }
                    
                    Iterator<Entry<Integer, List<PairInt>>> iter = 
                        extSlopeLocations.entrySet().iterator();
                    while (iter.hasNext()) {
                        Entry<Integer, List<PairInt>> entry = iter.next();
                        for (PairInt p : entry.getValue()) {
                            int x = col + p.getX();
                            int y = row + p.getY();
                            
                            if ((x < 0) || (y < 0) || (x > (input.getWidth() - 1)) ||
                                (y > (input.getHeight() - 1))) {
                                continue;
                            }
                            
                            if (input.getValue(x, y) > 0) {
                                int locationNumber = entry.getKey().intValue();
                                counts[locationNumber]++;
                            }
                        }
                    }
                    CountingSort.sortByDecr(counts, locationNumbers, 3);
                    
                    /*
                      66700
                      6###0
                      5# #1
                      4###2
                      44322   count pixels in spaces 0 thru 7
                    */
                    // depending on top 2 results in c, fill in the n
                    for (int ii = 0; ii < 2; ii++) {
                        int n0 = locationNumbers[ii];
                        int x = col;
                        int y = row;
                        switch(n0) {
                            case 0:
                                if ((x + 1) < input.getWidth()) {
                                    x += 1;
                                }
                                if ((y + 1) < input.getHeight()) {
                                    y += 1;
                                }
                                break;
                            case 1:
                                if ((x + 1) < input.getWidth()) {
                                    x += 1;
                                }
                                break;
                            case 2:
                                if ((x + 1) < input.getWidth()) {
                                    x += 1;
                                }
                                if ((y - 1) > -1) {
                                    y -= 1;
                                }
                                break;
                            case 3:
                                if ((y - 1) > -1) {
                                    y -= 1;
                                }
                                break;
                            case 4:
                                if ((x - 1) > -1) {
                                    x -= 1;
                                }
                                if ((y - 1) > -1) {
                                    y -= 1;
                                }
                                break;
                            case 5:
                                if ((x - 1) > -1) {
                                    x -= 1;
                                }
                                break;
                            case 6:
                                if ((x - 1) > -1) {
                                    x -= 1;
                                }
                                if ((y + 1) < input.getHeight()) {
                                    y += 1;
                                }
                                break;
                            case 7:
                                if ((y + 1) < input.getHeight()) {
                                    y += 1;
                                }
                                break;
                            default:
                                break;
                        }
                        input.setValue(x, y, 1);
                    }
                }
            }
        }
    }    
    
}
