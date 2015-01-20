package algorithms.compGeometry;

import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

/**
 *
 * @author nichole
 */
public class PerimeterFinder {
    
    /**
     * for the bounds of rows present in points, find the min and max of columns
     * and return the result as a map with key = row number, value = pair of
     * ints with x being the start column for that row and y being the stop
     * column (inclusive) for that row.  note that rows without points were
     * interpreted from their surrounding rows.
     * 
     * @param points
     * @param outputRowMinMax output populated as the min and max of rows are 
     * determined.
     * @return 
     */
    public Map<Integer, PairInt> find(PairIntArray points, int[] outputRowMinMax) {
        
        int minY = Integer.MAX_VALUE;
        int maxY = Integer.MIN_VALUE;
        
        Map<Integer, List<Integer> > rowCols = new 
            HashMap<Integer, List<Integer> >();
        
        // O(N)
        for (int i = 0; i < points.getN(); i++) {
            
            int x = points.getX(i);
            int y = points.getY(i);
            if (y < minY) {
                minY = y;
            }
            if (y > maxY) {
                maxY = y;
            }
            
            Integer row = Integer.valueOf(y);
            
            List<Integer> cols = rowCols.get(row);
            if (cols == null) {
                cols = new ArrayList<Integer>();
            }
            cols.add(Integer.valueOf(x));
            
            rowCols.put(row, cols);
        }
        
        outputRowMinMax[0] = minY;
        outputRowMinMax[1] = maxY;
        
        return find(outputRowMinMax, rowCols);
    }
    
    /**
     * for the bounds of rows present in points, find the min and max of columns
     * and return the result as a map with key = row number, value = pair of
     * ints with x being the start column for that row and y being the stop
     * column (inclusive) for that row.  note that rows without points were
     * interpreted from their surrounding rows.
     * 
     * @param points
     * @param outputRowMinMax output populated as the min and max of rows are 
     * determined.
     * @return 
     */
    public Map<Integer, PairInt> find(Set<PairInt> points, int[] outputRowMinMax) {
        
        int minY = Integer.MAX_VALUE;
        int maxY = Integer.MIN_VALUE;
        
        Map<Integer, List<Integer> > rowCols = new 
            HashMap<Integer, List<Integer> >();
        
        // O(N)
        for (PairInt p : points) {
            
            int x = p.getX();
            int y = p.getY();
            if (y < minY) {
                minY = y;
            }
            if (y > maxY) {
                maxY = y;
            }
            
            Integer row = Integer.valueOf(y);
            
            List<Integer> cols = rowCols.get(row);
            if (cols == null) {
                cols = new ArrayList<Integer>();
            }
            cols.add(Integer.valueOf(x));
            
            rowCols.put(row, cols);
        }
        
        outputRowMinMax[0] = minY;
        outputRowMinMax[1] = maxY;
        
        return find(outputRowMinMax, rowCols);
    }
    
    protected Map<Integer, PairInt> find(int[] minMax, 
        Map<Integer, List<Integer> > rowCols) {
        
        int minY = minMax[0];
        int maxY = minMax[1];
        
        // first find min and max cols for each row, then interpolate those
        // for missing rows
        Map<Integer, PairInt> rowColRange = new HashMap<Integer, PairInt>();
        
        Iterator<Entry<Integer, List<Integer> > > iter = rowCols.entrySet().iterator();
        
        // a little more than O(N)
        while (iter.hasNext()) {
            Entry<Integer, List<Integer> > entry = iter.next();
            Integer row = entry.getKey();
            List<Integer> cols = entry.getValue();
            Collections.sort(cols);
            
            PairInt colRange = new PairInt(cols.get(0).intValue(),
                cols.get(cols.size() - 1));
            
            rowColRange.put(row, colRange);
        }
        
        // a little more than O(N)
        
        // find missing rows and interpolate with the surrounding rows
        for (int rowIdx = (minY + 1); rowIdx < maxY; rowIdx++) {
            
            Integer row = Integer.valueOf(rowIdx);
            
            if (!rowColRange.containsKey(row)) {
                
                // previous row will already have been filled
                Integer prevRow = Integer.valueOf(rowIdx - 1);
                Integer nextRow = null;
                
                for (int nextRowIdx = (rowIdx + 1); nextRowIdx <= maxY; 
                    nextRowIdx++) {
                    
                    nextRow = Integer.valueOf(nextRowIdx);
                    
                    if (rowColRange.containsKey(nextRow)) {
                        break;
                    }
                }
                
                PairInt prevColRange = rowColRange.get(prevRow);
                PairInt nextColRange = rowColRange.get(nextRow);
                int nRows = nextRow.intValue() - prevRow.intValue() + 1;
                
                for (int r = rowIdx; r < nextRow.intValue(); r++) {
                    
                    /*
                     row 0:    @  . . @ . . @ . .
                     row 1:    .  . . . . . . . .
                     row 2:    .  . . . . . . . .
                     row 3:    .  @ . . . . . @ .                    
                    */
                    float deltaStart = (nextColRange.getX() - 
                        prevColRange.getX())/(float)(nRows - 1);
                    
                    int nR = r - prevRow.intValue();
                    
                    int cStart = (int)(prevColRange.getX() + (nR * deltaStart));
                    
                    float deltaStop = (nextColRange.getY() - 
                        prevColRange.getY())/(float)(nRows - 1);
                    
                    int cStop = (int)(prevColRange.getY() + (nR * deltaStop));
                    
                    PairInt colRange = new PairInt(cStart, cStop);
                    
                    rowColRange.put(Integer.valueOf(r), colRange);
                }
                
                // skip past the next row
                rowIdx = nextRow.intValue();
            }
        }
        
        return rowColRange;
    }
}
