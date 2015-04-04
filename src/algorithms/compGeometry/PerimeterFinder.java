package algorithms.compGeometry;

import algorithms.misc.MiscMath;
import algorithms.util.PairInt;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.List;
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
   
    /**
     * For the given points, find the ranges of columns that bound the points
     * that are contiguous and the points that are completely 
     * enclosed within points but not part of the set.
     * This returns an outline of the points attempting to correct for
     * concave portions of the hull, that is, it is roughly a concave hull
     * that includes embedded PairInts that are not in the set points.
     * 
     * @param points
     * @param outputRowMinMax output populated as the min and max of rows are 
     * determined.
     * @param outputEmbeddedGapPoints a return variable holding the found 
     * embedded points that are not in the "points" set, but are enclosed by it.
     * @return 
     */
    public Map<Integer, List<PairInt>> find2(Set<PairInt> points, 
        int[] outputRowMinMax, Set<PairInt> outputEmbeddedGapPoints) {
        
        if (points == null) {
	    	throw new IllegalArgumentException("points cannot be null");
        }
        if (outputRowMinMax == null) {
	    	throw new IllegalArgumentException("outputRowMinMax cannot be null");
        }
        
        //== O(N):
        int[] minMaxXY = MiscMath.findMinMaxXY(points);
        int minX = minMaxXY[0];
        int maxX = minMaxXY[1];
        int minY = minMaxXY[2];
        int maxY = minMaxXY[3];
        
        outputRowMinMax[0] = minY;
        outputRowMinMax[1] = maxY;       
        
        //== O((maxX-minX+1)*(maxY-minY+1)):
        // key holds row number
        // value holds (first column number, last column number) for points in the row
        Map<Integer, List<PairInt>> rowColRanges = findRowColRanges(points, 
            minX, maxX, minY, maxY);        
        
        /*
        TODO: alter algorithm to look at whether the contiguous gap pixels
        are unbounded, rather than looking at each gap along a row.
        TODO: consider where can make changes to iterate over data by
        point in "points" instead of minX, maxX, minY, maxY to reduce the
        runtime and keep it easier to make polynomial estimate.
        
        -- visit each row in rowColRanges and store the gaps if any:
           Gap: row, startGap, stopGapInclusive
           store in a stack so that the last pushed are the smallest rows.
           * runtime complexity is
             O((maxY-minY+1)*k) where k is the number of contig ranges per row
           * storage is a stack, inserts are O(1) rt and later pops are O(1) rt
        -- visit each stack member using a dfs style w/ a push for true neighbors.
           visit of gap u:
              neighbors of u are gaps immediately above in row that are "connected".
              for each neighbor v of u:
                  process the pair as a new group or add to existing group for u or v.
           * runtime complexity is > O(m) where m is the number of contig gap ranges by row
           * storage is
             reverse lookup: Gap -> group number.  needs fast search and insert
                 so Map<Gap, Integer>
             forward storage: need to be able to aggregate Gaps for a group,
                 in a random indexed container.
                 so List<List<Gap>>
        -- each contiguous gap is then present in List<List<Gap>> as first 
           level item
        
        -- for each Gap in a gap group:
              have row and gap.  use the check bounds for a row
              -- if any are not bounded, this invalidates the entire
                 group, so skip over it
                 -- else add the entire group to a List<List<Gap>>
                    that will be used for merging items in rowColRanges
           * runtime complexity is > O(m) where m is the number of contig gap ranges by row
           * storage is a List<List<Gap>>
        -- visit each item in the new gaps List<List<Gap>>
           condense the affected colRange for the row in rowColRanges
           * runtime complexity depends on the dataset...still roughly O(m) or less
        */ 
        
        // now, want to find gaps in the colRanges where there are points
        // in "points" above it and below it (i.e. the gap is completely embedded
        // in "points")
        
        //== O((maxY-minY+1)*k) where k is the number of contig ranges per row:
        for (int row = minY; row <= maxY; row++) {
            
            Integer key = Integer.valueOf(row);
            
            List<PairInt> colRanges = rowColRanges.get(key);
            
            List<Integer> mergeIndexes = new ArrayList<Integer>();
            
            for (int i = 1; i < colRanges.size(); i++) {
                
                int gapStart = colRanges.get(i - 1).getY() + 1;
                int gapStop = colRanges.get(i).getX() - 1;
                
                boolean bounded = boundedByPointsInHigherRows(row, gapStart,
                    gapStop, maxY, rowColRanges);
                
                if (!bounded) {
                    continue;
                }
                
                bounded = boundedByPointsInLowerRows(row, gapStart,
                    gapStop, minY, rowColRanges);
                
                if (!bounded) {
                    continue;
                }
                
                // the gap is bounded above and below, so
                // combine colRange before with current colRange
                mergeIndexes.add(Integer.valueOf(i));
            }
            
            if (!mergeIndexes.isEmpty()) {
                
                for (int i = (mergeIndexes.size() - 1); i > -1; i--) {
                    
                    int idx = mergeIndexes.get(i).intValue();
                    
                    PairInt edit = colRanges.get(idx - 1);
                    
                    int gapStart = edit.getY() + 1;
                    
                    PairInt current = colRanges.get(idx);
                    
                    int gapStop = current.getX() - 1;
                    
                    edit.setY(current.getY());
                                        
                    colRanges.remove(idx);
                    
                    for (int cIdx = gapStart; cIdx <= gapStop; cIdx++) {
                        outputEmbeddedGapPoints.add(new PairInt(cIdx, row));
                    }
                }
            }
        }
        
        return rowColRanges;
    }
    
    boolean boundedByPointsInHigherRows(int row, int gapStart, int gapStop,
        int maxRow, Map<Integer, List<PairInt>> rowColRange) {
        
        LinkedHashSet<Integer> gapRange = createRange(gapStart, gapStop);
        List<Integer> remove = new ArrayList<Integer>();
        
        for (int r = (row + 1); r <= maxRow; r++) {
            
            List<PairInt> colRanges = rowColRange.get(Integer.valueOf(r));
            
            boolean bounded = boundedByPoints(colRanges, gapRange, remove);
            
            if (bounded) {
                return true;
            }
        }
        
        return gapRange.isEmpty();
    }
    
    boolean boundedByPointsInLowerRows(int row, int gapStart, int gapStop,
        int minRow, Map<Integer, List<PairInt>> rowColRanges) {
        
        LinkedHashSet<Integer> gapRange = createRange(gapStart, gapStop);
        List<Integer> remove = new ArrayList<Integer>();
        
        for (int r = (row - 1); r >= minRow; r--) {
            
            List<PairInt> colRanges = rowColRanges.get(Integer.valueOf(r));

            boolean bounded = boundedByPoints(colRanges, gapRange, remove);
            
            if (bounded) {
                return true;
            }
        }
        
        return gapRange.isEmpty();
    }
    
    private boolean boundedByPoints(List<PairInt> colRanges, 
        LinkedHashSet<Integer> gapRange, List<Integer> remove) {
                            
        for (PairInt colRange : colRanges) {

            int col0 = colRange.getX();
            int col1 = colRange.getY();

            for (Integer gapColumn : gapRange) {

                int col = gapColumn.intValue();

                if ((col >= col0) && (col <= col1)) {
                    remove.add(gapColumn);
                }
            }
            for (Integer rm : remove) {
                gapRange.remove(rm);
            }
            remove.clear();

            if (gapRange.isEmpty()) {
                return true;
            }
        }

        return gapRange.isEmpty();
    }
    
    LinkedHashSet<Integer> createRange(int start, int stopInclusive) {
        
        int n = stopInclusive - start + 1;
        
        LinkedHashSet<Integer> range = new LinkedHashSet<Integer>(n);
        
        for (int i = 0; i < n; i++) {
            
            Integer value = Integer.valueOf(start + i);
            
            range.add(value);
        }
        
        return range;
    }
    
    /**
     * for the given points, find the ranges of contiguous columns and return
     * that by row. 
     * runtime complexity is O((maxX-minX+1)*(maxY-minY+1)), so this is
     * larger than O(N) for points datasets that are less dense than the min
     * max range.
     * 
     * @param points
     * @param outputRowMinMax output populated as the min and max of rows are 
     * determined.
     * @return a map with key = row, value = list of contiguous points in the
     * row.
     */
    Map<Integer, List<PairInt>> findRowColRanges(Set<PairInt> points, 
        int minX, int maxX, int minY, int maxY) {
        
        if (points == null) {
	    	throw new IllegalArgumentException("points cannot be null");
        }
        
        //TODO: consider whether this needs to have an entry for every row
        // between minY and maxY, and if not, then note that invoker
        // needs to check for null, even within expected range,
        // and change the algorithm here to iterate over points instead of 
        // row and col.
        
        // key holds row number
        // value holds (first column number, last column number) for points in the row
        Map<Integer, List<PairInt>> rowColRange = new HashMap<Integer, List<PairInt>>();
       
        //O(N):
        for (int row = minY; row <= maxY; row++) {
            
            List<PairInt> colRanges = new ArrayList<PairInt>();
            Integer key = Integer.valueOf(row);
            rowColRange.put(key, colRanges);
            
            PairInt currentColRange = null;
            
            for (int col = minX; col <= maxX; col++) {
                
                boolean contains = points.contains(new PairInt(col, row));
                
                if (currentColRange == null) {
                    if (contains) {
                        currentColRange = new PairInt(col, col);
                    }
                } else {
                    if (contains) {
                        if (col == (currentColRange.getY() + 1)) {
                            currentColRange.setY(col);
                        } else {
                            colRanges.add(currentColRange);
                            currentColRange = null;
                        }
                    } else {
                        colRanges.add(currentColRange);
                        currentColRange = null;
                    }
                }
            }
            
            //store last currentColRange
            if (currentColRange != null) {
                
                if (colRanges.isEmpty()) {
                    
                    colRanges.add(currentColRange);
                    
                } else if ( 
                    !colRanges.get(colRanges.size() - 1).equals(currentColRange)) {
                    
                    colRanges.add(currentColRange);
                }
            }
        }
        
        return rowColRange;
    }

    private class Gap {
        
        private final int row;
        
        private final int start;
        
        private final int stopInclusive;
        
        public Gap(int rowNumber, int startColumn, int stopColumnInclusive) {
            row = rowNumber;
            start = startColumn;
            stopInclusive = stopColumnInclusive;
        }
        
        /**
         * @return the row
         */
        public int getRow() {
            return row;
        }

        /**
         * @return the start
         */
        public int getStart() {
            return start;
        }

        /**
         * @return the stopInclusive
         */
        public int getStopInclusive() {
            return stopInclusive;
        }

        @Override
        public boolean equals(Object obj) {
            
            if (!(obj instanceof Gap)) {
                return false;
            }
            
            Gap other = (Gap)obj;
            
            if ((other.getRow() == row) && (other.getStart() == start) &&
                (other.getStopInclusive() == stopInclusive)) {
                return true;
            }
            
            return false;
        }

        @Override
        public int hashCode() {

            int hash = fnvHashCode(this.row, this.start, this.stopInclusive);

            return hash;
        }

        int fnv321aInit = 0x811c9dc5;
        int fnv32Prime = 0x01000193;

        protected int fnvHashCode(int i0, int i1, int i2) {

            /*
             * hash = offset_basis
             * for each octet_of_data to be hashed
             *     hash = hash xor octet_of_data
             *     hash = hash * FNV_prime
             * return hash
             *
             * Public domain:  http://www.isthe.com/chongo/src/fnv/hash_32a.c
             */
            int hash = 0;

            int sum = fnv321aInit;

            // xor the bottom with the current octet.
            sum ^= i0;

            // multiply by the 32 bit FNV magic prime mod 2^32
            sum *= fnv32Prime;

            sum ^= i1;

            sum *= fnv32Prime;
            
            sum ^= i2;

            sum *= fnv32Prime;

            hash = sum;

            return hash;
        }

    }
}
