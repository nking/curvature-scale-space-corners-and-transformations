package algorithms.compGeometry;

import algorithms.util.PairInt;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

/**
 * class to create a map of rows with start and stop column bounds (inclusive)
 * for each row that bounding the area occupied by given points.
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
    public Map<Integer, PairInt> find(Set<PairInt> points, int[] outputRowMinMax) {
        
        if (points == null) {
	    	throw new IllegalArgumentException("points cannot be null");
        }
        if (outputRowMinMax == null) {
	    	throw new IllegalArgumentException("outputRowMinMax cannot be null");
        }
      
        // key holds row number
        // value holds (first column number, last column number) for points in the row
        Map<Integer, PairInt> rowColRange = new HashMap<Integer, PairInt>();
        
        int minY = Integer.MAX_VALUE;
        int maxY = Integer.MIN_VALUE;
            
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
            
            Integer key = Integer.valueOf(y);
            PairInt value = rowColRange.get(key);
            if (value == null) {
                value = new PairInt(x, x);
                rowColRange.put(key, value);
            } else {
                if (x < value.getX()) {
                    value.setX(x);
                } 
                if (x > value.getY()) {
                    value.setY(x);
                }
            }
        }
        
        outputRowMinMax[0] = minY;
        outputRowMinMax[1] = maxY;
        
        for (int row = minY; row < maxY; row++) {
            
            Integer key = Integer.valueOf(row);
            PairInt value = rowColRange.get(key);
                        
            if ((value != null) && (value.getX() < value.getY())) {
                continue;
            }
            if ((value != null) && (value.getX() == value.getY())
                && (row == minY)) {
                continue;
            }
            
            // else, find values from closest previous and next rows
            Integer prevKey = Integer.valueOf(row - 1);
            PairInt prevValue = rowColRange.get(prevKey);
            
            int nextRow = -1;
            Integer nextKey = null;
            PairInt nextValue = null;
            for (nextRow = row + 1; nextRow <= maxY; nextRow++) {
                nextKey = Integer.valueOf(nextRow);
                nextValue = rowColRange.get(nextKey);
                if (nextValue != null) {
                    break;
                }
            }
           
            float dRowX = (float)(nextValue.getX() - prevValue.getX())/
                (float)(nextRow - prevKey.intValue());
            
            float x = (float)prevValue.getX() + 
                (float)(row - prevKey.intValue()) * dRowX;
                
            float dRowY = (float)(nextValue.getY() - prevValue.getY())/
                (float)(nextRow - prevKey.intValue());
            
            float y = (float)prevValue.getY() + 
                (float)(row - prevKey.intValue()) * dRowY;
            
            value = new PairInt(Math.round(x), Math.round(y));
            
            rowColRange.put(key, value);
        }
       
        return rowColRange;
    }
   
}
