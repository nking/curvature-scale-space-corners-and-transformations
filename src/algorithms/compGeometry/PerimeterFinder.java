package algorithms.compGeometry;

import algorithms.misc.MiscMath;
import algorithms.util.PairInt;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Stack;

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
        
        // runtime complexity is > O(m) where m is the number of contig gap ranges by row
        List<List<Gap>> contiguousGaps = findContiguousGaps(rowColRanges,
            minX, maxX, minY, maxY);
        
        //runtime complexity is > O(m) where m is the number of contig gap ranges by row
        Set<Gap> embeddedGaps = findBoundedGaps(contiguousGaps, minY, maxY, 
            rowColRanges);
        
        for (Gap gap : embeddedGaps) {
            
            int row = gap.getRow();
            
            List<PairInt> colRanges = rowColRanges.get(Integer.valueOf(row));
            
            for (int i = (colRanges.size() - 1); i > 0; i--) {
                
                int gapStart = colRanges.get(i - 1).getY() + 1;
                int gapStop = colRanges.get(i).getX() - 1;
                
                if ((gap.getStart() == gapStart) && (gap.getStopInclusive() == gapStop)) {
                    
                    PairInt edit = colRanges.get(i - 1);
                                        
                    PairInt current = colRanges.get(i);
                                        
                    edit.setY(current.getY());
                                        
                    colRanges.remove(i);
                    
                    for (int cIdx = gapStart; cIdx <= gapStop; cIdx++) {
                        outputEmbeddedGapPoints.add(new PairInt(cIdx, row));
                    }
                }
            }    
        }
        
        /*
        TODO: consider where can make changes to iterate over data by
        point in "points" instead of minX, maxX, minY, maxY to reduce the
        runtime and keep it easier to make polynomial estimate.
        */ 
        
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

    /**
     * Find the gaps in rowColRanges and put them in same group if they are 
     * connected.  Note that diagonal pixel are not considered connected
     * (though this may change if test cases show they should be).
     * 
     * @param rowColRanges
     * @param minX
     * @param maxX
     * @param minY
     * @param maxY
     * @return 
     */
    protected List<List<Gap>> findContiguousGaps(Map<Integer, List<PairInt>> 
        rowColRanges, int minX, int maxX, int minY, int maxY) {
        
        // ---------------- store the gaps in a stack ---------------------
        // runtime complexity is
        // O((maxY-minY+1)*k) where k is the number of contig ranges per row
        
        // stack is lifo, so push in reverse order of desired use
        Stack<Gap> stack = findGaps(rowColRanges, minX, maxX, minY, maxY);
        
        // ----------------- find contiguous gaps ------------------------        
        List<List<Gap>> gapGroups = new ArrayList<List<Gap>>();
        
        if (stack.isEmpty()) {
            return gapGroups;
        }
        
        Map<Gap, Integer> gapToIndexMap = new HashMap<Gap, Integer>();
        
        Map<Gap, Boolean> visited = new HashMap<Gap, Boolean>();
                
        while (!stack.isEmpty()) {
            
            Gap uNode = stack.pop();
            
            if (visited.containsKey(uNode)) {
                continue;
            }
            
            visited.put(uNode, Boolean.TRUE);
            
            int row = uNode.getRow();
            int uStart = uNode.getStart();
            int uStopIncl = uNode.getStopInclusive();
            
            // search in row above for a neighbor.
            // to save space, just looking in rowColRanges
            
            int vRow = row + 1;
            
            if (vRow > maxY) {
                continue;
            }
            
            // ------ find the connected neighbors of u, below u ---------
            
            List<Gap> vNodes = new ArrayList<Gap>();
            
            // these are ordered by increasing range:
            List<PairInt> colRanges = rowColRanges.get(Integer.valueOf(vRow));
            
            /*
                  |----|    |----|
                @@@@@@@@@            case 0
               @@@@@                 case 1
                      @@@@@          case 2
                    @@               case 3
            */
            for (int i = 1; i < colRanges.size(); i++) {
                int gapStart = colRanges.get(i - 1).getY() + 1;
                int gapStop = colRanges.get(i).getX() - 1;
                
                Gap vNode = null;
                
                // is any portion of uStart is within col0:col2
                if ((uStart <= gapStart) && (uStopIncl >= gapStop)) {
                    vNode = new Gap(vRow, gapStart, gapStop);
                } else if ((uStopIncl >= gapStart) && (uStopIncl <= gapStop)) {
                    vNode = new Gap(vRow, gapStart, gapStop);
                } else if ((uStart >= gapStart) && (uStart <= gapStop)) {
                    vNode = new Gap(vRow, gapStart, gapStop);
                }
                                
                if (vNode != null) {
                    
                    vNodes.add(vNode);                    
                }
            }
            
            // ---------------- process the neighbors -----------------
            
            for (Gap vNode : vNodes) {
                // process each node.  add to existing group or start a new one
                Integer uIdx = gapToIndexMap.get(uNode);
                Integer vIdx = gapToIndexMap.get(vNode);
                Integer groupIdx = null;
                if ((uIdx != null) && (vIdx != null)) {
                    if (uIdx.intValue() == vIdx.intValue()) {
                        // these are already in the same group... should not be visited again though
                        groupIdx = uIdx;
                    } else {
                        Integer moveTo;
                        Integer moveFrom;
                        if (uIdx.intValue() < vIdx.intValue()) {
                            moveTo = uIdx;
                            moveFrom = vIdx;
                        } else {
                            moveTo = vIdx;
                            moveFrom = uIdx;
                        }
                        List<Gap> moveFromG = gapGroups.get(moveFrom);
                        for (Gap g : moveFromG) {
                            gapToIndexMap.put(g, moveTo);
                        }
                        gapGroups.get(moveTo).addAll(moveFromG);
                        gapGroups.get(moveFrom).clear();
                        
                        int z = 1;
                    }
                } else if (uIdx != null) {
                    groupIdx = uIdx;
                    List<Gap> group = gapGroups.get(groupIdx);
                    group.add(vNode);
                    gapToIndexMap.put(vNode, groupIdx);
                } else if (vIdx != null) {
                    groupIdx = vIdx;
                    List<Gap> group = gapGroups.get(groupIdx);
                    group.add(uNode);
                    gapToIndexMap.put(uNode, groupIdx);
                } else {
                    // both are null
                    List<Gap> group = new ArrayList<Gap>();
                    group.add(uNode);
                    group.add(vNode);
                    groupIdx = Integer.valueOf(gapGroups.size());
                    gapGroups.add(group);
                    gapToIndexMap.put(uNode, groupIdx);
                    gapToIndexMap.put(vNode, groupIdx);
                }
                
                System.out.println(groupIdx + " ==> u " + uNode.toString() 
                    + " v " + vNode.toString());
                
                stack.push(vNode);
            }
            
            if (vNodes.isEmpty()) {
                // store u alone
                Integer groupIdx = gapToIndexMap.get(uNode);
                if (groupIdx == null) {
                    groupIdx = Integer.valueOf(gapGroups.size());
                    List<Gap> group = new ArrayList<Gap>();
                    group.add(uNode);
                    gapGroups.add(group);
                    gapToIndexMap.put(uNode, groupIdx);
                }
            }
        }
        
        // condense:
        boolean hasEmpty = false;
        for (List<Gap> group : gapGroups) {
            if (group.isEmpty()) {
                hasEmpty = true;
                break;
            }
        }
        if (hasEmpty) {
            List<List<Gap>> tmp = new ArrayList<List<Gap>>();
            for (List<Gap> group : gapGroups) {
                if (!group.isEmpty()) {
                    tmp.add(group);
                }
            }
            gapGroups = tmp;
        }
        /*
        for (int i = 0; i < gapGroups.size(); i++) {
            List<Gap> group = gapGroups.get(i);
            System.out.println("group: " + i);
            for (Gap gap : group) {
                System.out.println("  " + gap.toString());
            }
        }
        */
        
        return gapGroups;
    }
    
    Stack<Gap> findGaps(Map<Integer, List<PairInt>> 
        rowColRanges, int minX, int maxX, int minY, int maxY) {
        
        // ------- store the gaps in a stack --------
        // runtime complexity is
        // O((maxY-minY+1)*k) where k is the number of contig ranges per row
        
        // stack is lifo, so push in reverse order of desired use
        Stack<Gap> stack = new java.util.Stack<Gap>();
        
        //for (int row = minY; row <= maxY; row++) {
        for (int row = maxY; row >= minY; row--) {
            
            Integer key = Integer.valueOf(row);
            
            List<PairInt> colRanges = rowColRanges.get(key);
                        
            //for (int i = 1; i < colRanges.size(); i++) {
            for (int i = (colRanges.size() - 1); i > 0; i--) {
                
                int gapStart = colRanges.get(i - 1).getY() + 1;
                int gapStop = colRanges.get(i).getX() - 1;
                
                Gap gap = new Gap(row, gapStart, gapStop);
                
                stack.push(gap);
            }
        }
        
        return stack;
    }

    protected Set<Gap> findBoundedGaps(List<List<Gap>> contiguousGapGroups, 
        int minY, int maxY, Map<Integer, List<PairInt>> rowColRanges) {
        
        Set<Gap> embeddedGaps = new HashSet<Gap>();
        
        for (List<Gap> contiguousGap : contiguousGapGroups) {
            
            boolean notBounded = false;
            
            for (Gap gap : contiguousGap) {
                
                boolean bounded = boundedByPointsInHigherRows(gap.getRow(), 
                    gap.getStart(), gap.getStopInclusive(), maxY, rowColRanges);
                
                if (!bounded) {
                    notBounded = true;
                    break;
                }
                
                bounded = boundedByPointsInLowerRows(gap.getRow(), 
                    gap.getStart(), gap.getStopInclusive(), minY, rowColRanges);
                
                if (!bounded) {
                    notBounded = true;
                    break;
                }
            }
            
            if (!notBounded) {
                embeddedGaps.addAll(contiguousGap);
            }
        }
        
        return embeddedGaps;
    }

    static class Gap {
        
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

        @Override
        public String toString() {
            
            StringBuilder sb = new StringBuilder();
            sb.append("row=").append(Integer.toString(row))
                .append(" cols=").append(Integer.toString(start))
                .append(":").append(Integer.toString(stopInclusive));
            
            return sb.toString();
        }
    }
}
