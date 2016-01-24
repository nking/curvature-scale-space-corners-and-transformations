package algorithms.compGeometry;

import algorithms.imageProcessing.DFSConnectedGroupsFinder;
import algorithms.imageProcessing.MiscellaneousCurveHelper;
import algorithms.imageProcessing.ZhangSuenLineThinner;
import algorithms.util.PathStep;
import algorithms.misc.MiscMath;
import algorithms.util.BitVectorRepresentation;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.PairIntArrayWithColor;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Stack;
import java.util.logging.Logger;

/**
 * class to create a map of rows with start and stop column bounds (inclusive)
 * for each row to bound the area occupied by given points while including also
 * groups of points that are completely embedded in the point set though not
 * part of the point set.  Those embedded points are retrievable separately.
 * The data structure is somewhat like a concave hull.  It's used to help
 * quickly scan a region and not include the tops of mountains that the sky
 * concave hull is bent around, for example.
 *
 * @author nichole
 */
public class PerimeterFinder {

    private Logger log = Logger.getLogger(this.getClass().getName());

    private boolean debug = true;

    /**
     * For the given points, find the ranges of columns that bound the points
     * that are contiguous and the points that are completely
     * enclosed within points but not part of the set.
     * This returns an outline of the points attempting to correct for
     * concave portions of the hull, that is, it is roughly a concave hull
     * that includes embedded PairInts that are not in the set points.
     * runtime complexity is approx O(N).
     * @param points
     * @param outputRowMinMax output populated as the min and max of rows are
     * determined.
     * @param imageMaxColumn maximum column index of image (used to understand
     * when a pixel is in the last column)
     * @param outputEmbeddedGapPoints a return variable holding the found
     * embedded points that are not in the "points" set, but are enclosed by it.
     * @return map w/ key being row number, value being a list of column ranges
     * that inclusively bound the points in that row.
     */
    public Map<Integer, List<PairInt>> find(Set<PairInt> points,
        int[] outputRowMinMax, int imageMaxColumn,
        Set<PairInt> outputEmbeddedGapPoints) {

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

        /*
        TODO: consider where can make changes to iterate over data by
        point in "points" instead of minX, maxX, minY, maxY to reduce the
        runtime and keep it easier to make polynomial estimate.
        */

        //== O((maxX-minX+1)*(maxY-minY+1)):
        // key holds row number
        // value holds (first column number, last column number) for points in the row
        Map<Integer, List<PairInt>> rowColRanges = findRowColRanges(points,
            minX, maxX, minY, maxY);

        // runtime complexity is > O(m) where m is the number of contig gap ranges by row
        List<List<Gap>> contiguousGaps = findContiguousGaps(rowColRanges,
            minX, maxX, minY, maxY);

        //runtime complexity is > O(m) where m is the number of contig gap ranges by row
        List<List<Gap>> embeddedGapGroups = findBoundedGaps(contiguousGaps,
            minY, maxY, imageMaxColumn, rowColRanges);

        // update the rowColRanges to encapsulate the truly embedded points too
        for (List<Gap> embeddedGroup : embeddedGapGroups) {

            updateRowColRangesForVerifiedEmbedded(rowColRanges,
                embeddedGroup, outputEmbeddedGapPoints);

        }

        return rowColRanges;
    }

    /**
     * For the given points, find the ranges of columns that bound the points
     * that are contiguous and the points that are completely
     * enclosed within points but not part of the set.
     * This returns an outline of the points attempting to correct for
     * concave portions of the hull, that is, it is roughly a concave hull
     * that includes embedded PairInts that are not in the set points.
     * @param rowMinMax
     * @param imageMaxColumn
     * @param rowColRanges
     * @return 
     */
    public Set<PairInt> findEmbeddedGivenRowData(
        int[] rowMinMax, int imageMaxColumn,
        Map<Integer, List<PairInt>> rowColRanges) {

        if (rowColRanges == null) {
	    	throw new IllegalArgumentException("rowColRanges cannot be null");
        }
        if (rowMinMax == null) {
	    	throw new IllegalArgumentException("outputRowMinMax cannot be null");
        }

        int minX = Integer.MAX_VALUE;
        int maxX = Integer.MIN_VALUE;
        for (int row = rowMinMax[0]; row <= rowMinMax[1]; row++) {
            List<PairInt> rcr = rowColRanges.get(Integer.valueOf(row));
            if (rcr.isEmpty()) {
                continue;
            }
            int x0 = rcr.get(0).getX();
            int xf = rcr.get(rcr.size() - 1).getY();
            if (x0 < minX) {
                minX = x0;
            }
            if (xf > maxX) {
                maxX = xf;
            }
        }

        // runtime complexity is > O(m) where m is the number of contig gap ranges by row
        List<List<Gap>> contiguousGaps = findContiguousGaps(rowColRanges,
            minX, maxX, rowMinMax[0], rowMinMax[1]);

        //runtime complexity is > O(m) where m is the number of contig gap ranges by row
        List<List<Gap>> embeddedGapGroups = findBoundedGaps(contiguousGaps,
            rowMinMax[0], rowMinMax[1], imageMaxColumn, rowColRanges);

        Set<PairInt> outputEmbeddedGapPoints = new HashSet<PairInt>();

        // update the rowColRanges to encapsulate the truly embedded points too
        for (List<Gap> embeddedGroup : embeddedGapGroups) {

            updateRowColRangesForVerifiedEmbedded(rowColRanges,
                embeddedGroup, outputEmbeddedGapPoints);
        }

        return outputEmbeddedGapPoints;
    }

    boolean boundedByPointsInHigherRows(int row, int gapStart, int gapStop,
        int maxRow, int imageMaxColumn, Map<Integer, List<PairInt>> rowColRange) {

        LinkedHashSet<Integer> gapRange = createRange(gapStart, gapStop);

        List<Integer> remove = new ArrayList<Integer>();

        for (int r = (row + 1); r <= maxRow; r++) {

            List<PairInt> colRanges = rowColRange.get(Integer.valueOf(r));

            boolean bounded = boundedByPoints(colRanges, imageMaxColumn,
                gapRange, remove);

            if (bounded) {
                return true;
            }
        }

        return gapRange.isEmpty();
    }

    boolean boundedByPointsInLowerRows(int row, int gapStart, int gapStop,
        int minRow, int imageMaxColumn,
        Map<Integer, List<PairInt>> rowColRanges) {

        LinkedHashSet<Integer> gapRange = createRange(gapStart, gapStop);
        List<Integer> remove = new ArrayList<Integer>();

        for (int r = (row - 1); r >= minRow; r--) {

            List<PairInt> colRanges = rowColRanges.get(Integer.valueOf(r));

            boolean bounded = boundedByPoints(colRanges, imageMaxColumn,
                gapRange, remove);

            if (bounded) {
                return true;
            }
        }

        return gapRange.isEmpty();
    }

    private boolean boundedByPoints(List<PairInt> colRanges, int imageMaxColumn,
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
     * @return a map with key = row, value = list of contiguous points in the
     * row.
     */
    Map<Integer, List<PairInt>> findRowColRanges(Set<PairInt> points,
        int minX, int maxX, int minY, int maxY) {

        if (points == null) {
	    	throw new IllegalArgumentException("points cannot be null");
        }

        /*
        TODO: consider where can make changes to iterate over data by
        point in "points" instead of minX, maxX, minY, maxY to reduce the
        runtime and keep it easier to make polynomial estimate.
        */

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
     * connected.  Note that diagonal pixels are not considered connected
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
        Stack<Gap> stack = findGaps(rowColRanges, minX, maxX, minY, maxY);

        // ----------------- find contiguous gaps ------------------------
        List<List<Gap>> gapGroups = new ArrayList<List<Gap>>();

        if (stack.isEmpty()) {
            return gapGroups;
        }

        Map<Gap, Integer> gapToIndexMap = new HashMap<Gap, Integer>();

        Map<Gap, Boolean> visited = new HashMap<Gap, Boolean>();

        //TODO: add minimum and maximum runtime estimates here.

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

                //System.out.println(groupIdx + " ==> u " + uNode.toString()
                //    + " v " + vNode.toString());

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
        boolean hasAnEmpty = false;
        for (List<Gap> group : gapGroups) {
            if (group.isEmpty()) {
                hasAnEmpty = true;
                break;
            }
        }
        if (hasAnEmpty) {
            List<List<Gap>> tmp = new ArrayList<List<Gap>>(gapGroups.size() - 1);
            for (List<Gap> group : gapGroups) {
                if (!group.isEmpty()) {
                    tmp.add(group);
                }
            }
            gapGroups = tmp;
        }

        return gapGroups;
    }

    Stack<Gap> findGaps(Map<Integer, List<PairInt>>
        rowColRanges, int minX, int maxX, int minY, int maxY) {

        // ------- store the gaps in a stack --------
        // runtime complexity is
        // O((maxY-minY+1)*k) where k is the number of contig ranges per row

        // stack is lifo, so push in reverse order of desired use
        Stack<Gap> stack = new java.util.Stack<Gap>();

        for (int row = maxY; row >= minY; row--) {

            Integer key = Integer.valueOf(row);

            List<PairInt> colRanges = rowColRanges.get(key);

            for (int i = (colRanges.size() - 1); i > 0; i--) {

                int gapStart = colRanges.get(i - 1).getY() + 1;
                int gapStop = colRanges.get(i).getX() - 1;

                Gap gap = new Gap(row, gapStart, gapStop);

                stack.push(gap);
            }
        }

        return stack;
    }

    /**
     * given contiguousGapGroups, return the subset that are completely bounded
     * by points (note that the bounds are given by rowColRanges).
     * @param contiguousGapGroups
     * @param minY the minimum row of the points set used to construct rowColRanges
     * @param maxY the maximum row of the points set used to construct rowColRanges
     * @param imageMaxColumn the maximum column index of the image
     * @param rowColRanges map with key = row number and value = list of column
     * ranges for the presence of points.
     * @return subset of contiguousGapGroups that are completely bounded by
     * rowColRanges
     */
    protected List<List<Gap>> findBoundedGaps(List<List<Gap>> contiguousGapGroups,
        int minY, int maxY, int imageMaxColumn,
        Map<Integer, List<PairInt>> rowColRanges) {

        List<List<Gap>> contiguousBoundedGapGroups = new ArrayList<List<Gap>>();

        /* need a data structure to access all Gaps by row by number.
        Will use a Map with key=integer and value = set of Gaps.
        Since there are not usually very many Gaps per row, will not use
        a sorted list of Gaps as the value, but that might be something to
        consider in the future with stats of the total number of gaps in
        contiguousBoundedGapGroups compared to the number of pixels in sky.
        */
        Map<Integer, Set<Gap>> rowGapsMap = createRowMap(contiguousGapGroups);

        for (List<Gap> contiguousGap : contiguousGapGroups) {

            boolean notBounded = false;

            for (Gap gap : contiguousGap) {

                int row = gap.getRow();

                /*
                 check that the gap above it if any is not in a gap that
                 continues to the image boundary
                */
                boolean adjRowGapIsUnbounded =
                    adjacentGapIsConnectedToImageBoundary(gap.getStart(),
                        gap.getStopInclusive(),
                        rowGapsMap.get(Integer.valueOf(row + 1)),
                        rowColRanges.get(Integer.valueOf(row + 1)),
                        imageMaxColumn);

                if (adjRowGapIsUnbounded) {
                    notBounded = true;
                    break;
                }

                /*
                 same check for row below
                */
                adjRowGapIsUnbounded = adjacentGapIsConnectedToImageBoundary(
                    gap.getStart(), gap.getStopInclusive(),
                    rowGapsMap.get(Integer.valueOf(row - 1)),
                    rowColRanges.get(Integer.valueOf(row - 1)),
                    imageMaxColumn);

                if (adjRowGapIsUnbounded) {
                    notBounded = true;
                    break;
                }

                boolean bounded = boundedByPointsInHigherRows(row,
                    gap.getStart(), gap.getStopInclusive(), maxY,
                    imageMaxColumn, rowColRanges);

                if (!bounded) {
                    notBounded = true;
                    break;
                }

                bounded = boundedByPointsInLowerRows(row, gap.getStart(),
                    gap.getStopInclusive(), minY, imageMaxColumn, rowColRanges);

                if (!bounded) {
                    notBounded = true;
                    break;
                }
            }

            if (!notBounded) {
                contiguousBoundedGapGroups.add(contiguousGap);
            }
        }

        return contiguousBoundedGapGroups;
    }

    protected boolean updateForAddedPoints(
        Map<Integer, List<PairInt>> rowColRanges, int[] rowMinMax,
        Collection<PairInt> addedPoints) {

        boolean allAreWithinExistingRange = true;

        // ------------------- update rowColRanges --------------------------
        for (PairInt p : addedPoints) {

            int row = p.getY();
            Integer rowKey = Integer.valueOf(row);

            int col = p.getX();

            List<PairInt> colRanges = rowColRanges.get(rowKey);

            /*
             cases:
                no row or empty list for row in rowColRanges

                point is before first range in colRanges
                   adjacent to it or not

                point is after last range in colRanges
                   adjacent to it or not

                point is between ranges in colRanges
                   adjacent to a range or not

                point is within range in colRanges

            Note that when the insert is adjacent to a range, have to consider
            whether it is adjacent on both sides, in which case it is a merge.
            */

            if ((colRanges == null) || colRanges.isEmpty()) {

                // case: no row or empty list for row in rowColRanges

                allAreWithinExistingRange = false;

                if (colRanges == null) {
                    colRanges = new ArrayList<PairInt>();
                    rowColRanges.put(rowKey, colRanges);
                }
                colRanges.add(new PairInt(col, col));
                if (row < rowMinMax[0]) {
                    rowMinMax[0] = row;
                }
                if (row > rowMinMax[1]) {
                    rowMinMax[1] = row;
                }
                continue;
            }

            int n = colRanges.size();

            if (col < colRanges.get(0).getX()) {

                //case: point is before first range

                allAreWithinExistingRange = false;

                PairInt first = colRanges.get(0);
                if (col == (first.getX() - 1)) {
                    first.setX(col);
                } else {
                    PairInt add = new PairInt(col, col);
                    colRanges.add(0, add);
                }
            } else if (col > colRanges.get(n - 1).getY()) {

                //case: point is after last range

                allAreWithinExistingRange = false;

                PairInt last = colRanges.get(n - 1);

                if (col == (last.getY() + 1)) {
                    last.setY(col);
                } else {
                    PairInt add = new PairInt(col, col);
                    colRanges.add(add);
                }
            } else if (n == 1) {
                // was not before or after the only range, so check within
                PairInt current = colRanges.get(0);
                if (!((col >= current.getX()) && (col <= current.getY()))) {
                    // not in range.  this should not happen
                    throw new IllegalStateException("error in algorithm: point "
                        + col + ", " + row + " was not added to rowColRanges");
                }
            } else {

                //case: point is between ranges in colRanges
                // or
                //case: point is within range in colRanges

                boolean added = false;

                for (int i = (colRanges.size() - 1); i > 0; i--) {

                    PairInt current = colRanges.get(i);

                    PairInt prev = colRanges.get(i - 1);

                    //case: point is within range in colRanges
                    if ((col >= current.getX()) && (col <= current.getY())) {
                        // no need to change range
                        added = true;
                        break;
                    } else if ((i == 1) &&
                        (col >= prev.getX()) && (col <= prev.getY())) {
                        added = true;
                        break;
                    }

                    //case: point is between ranges in colRanges

                    allAreWithinExistingRange = false;

                    if ((col >= prev.getY()) && (col <= current.getX())) {

                        if (col == (prev.getY() + 1)) {

                            //if adjacent to current range too, merge them
                            if (col == (current.getX() - 1)) {
                                // extend prev to current end and remove current
                                prev.setY(current.getY());
                                colRanges.remove(current);
                            } else {
                                prev.setY(col);
                            }
                            added = true;
                            break;
                        } else if (col == (current.getX() - 1)) {
                            current.setX(col);
                            added = true;
                            break;
                        } else {
                            PairInt add = new PairInt(col, col);
                            colRanges.add(i, add);
                            added = true;
                            break;
                        }
                    }
                }
                if (!added) {

                    // find min and max of ranges:
                    int[] minMaxCols = findMinMaxColumns(rowColRanges,
                        rowMinMax);

                     List<PairInt> colRanges0 = rowColRanges.get(
                         Integer.valueOf(p.getY()));
                     StringBuffer sb = new StringBuffer();
                     for (PairInt cr : colRanges0) {
                         sb.append("colRange=" + cr.getX() + ":" + cr.getY());
                         sb.append("\n");
                     }

                    throw new IllegalStateException("point " + p.toString() +
                        " was not added to a colRange. " +
                        " colRanges minX=" + minMaxCols[0] +
                        " maxX=" + minMaxCols[1] + " minRow=" + rowMinMax[0] +
                        " maxRow=" + rowMinMax[1] + " " + sb.toString());
                }
            }
        }

        return allAreWithinExistingRange;
    }

    /**
     * runtime complexity is < O(N) where N is ~ imageMaxColumn * rowRange.
     * 
     * @param rowColRanges
     * @param rowMinMax
     * @param imageMaxColumn
     * @param addedPoints 
     */
    public void updateRowColRangesForAddedPoints(
        Map<Integer, List<PairInt>> rowColRanges, int[] rowMinMax,
        int imageMaxColumn, Collection<PairInt> addedPoints) {

        /*
        - update rowColRanges and rowMinMax for individual pixels.
          This is faster to update rowColRanges rather than create it anew
          from all points (not given in arguments).
          * runtime: O(N_addedPoints)
        - findContiguousGaps
          * runtime: > O(m) where m is the number of contig gap ranges by row
        - findBoundedGaps
          * runtime: > O(m) where m is the number of contig gap ranges by row
        - update rowColRanges to include the truly embedded gaps just verified
          * runtime:

        Note that can avoid the last 3 steps if the addedPoints all exist within
        existing ranges in rowColRanges.
        */

        boolean allAreWithinExistingRange = updateForAddedPoints(
            rowColRanges, rowMinMax, addedPoints);

        if (allAreWithinExistingRange) {
            return;
        }

        // ------- find minX and maxX -------
        // runtime is O(rowColRanges.size)
        int minY = rowMinMax[0];
        int maxY = rowMinMax[1];
        int minX = Integer.MAX_VALUE;
        int maxX = Integer.MIN_VALUE;
        for (int row = minY; row <= maxY; row++) {
            List<PairInt> colRanges = rowColRanges.get(Integer.valueOf(row));
            if ((colRanges == null) || colRanges.isEmpty()) {
                continue;
            }
            int n = colRanges.size();
            PairInt cr = colRanges.get(0);
            if (cr.getX() < minX) {
                minX = cr.getX();
            }
            if (n > 1) {
                cr = colRanges.get(n - 1);
            }
            if (cr.getY() > maxX) {
                maxX = cr.getY();
            }
        }

        // runtime complexity is > O(m) where m is the number of contig gap ranges by row
        List<List<Gap>> contiguousGaps = findContiguousGaps(rowColRanges,
            minX, maxX, minY, maxY);

        //runtime complexity is > O(m) where m is the number of contig gap ranges by row
        List<List<Gap>> embeddedGapGroups = findBoundedGaps(contiguousGaps, minY,
            maxY, imageMaxColumn, rowColRanges);

        // update the rowColRanges to encapsulate the truly embedded points too
        Set<PairInt> outputEmbeddedGapPoints = new HashSet<PairInt>();

        for (List<Gap> embeddedGroup : embeddedGapGroups) {

            updateRowColRangesForVerifiedEmbedded(rowColRanges,
                embeddedGroup, outputEmbeddedGapPoints);
        }
    }

    private void updateRowColRangesForVerifiedEmbedded(
        Map<Integer, List<PairInt>> rowColRanges,
        Collection<Gap> embeddedGaps, Set<PairInt> outputEmbeddedGapPoints) {

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
    }
    
    /**
     * find the pixels which are the borders of rowColRanges including concave
     * pixels, but excluding any pixels that are the image border pixels.
     * Note that rowColRanges has to represent a contiguous point set.
     *
     * @param rowColRanges the column bounds for each row of a contiguous
     * point set.
     * @param rowMinMax the minimum and maximum rows present in the contiguous
     * point set.
     * @param imageMaxColumn the maximum column in the image from which the
     * point set was derived.
     * @param imageMaxRow the maximum row in the image from which the point
     * set was derived.
     * @return
     */
    public Set<PairInt> getBorderPixels(Map<Integer, List<PairInt>> rowColRanges,
        int[] rowMinMax, int imageMaxColumn, int imageMaxRow) {

        Set<PairInt> borderPixels0 = getBorderPixels0(rowColRanges, rowMinMax, 
            imageMaxColumn, imageMaxRow);
        
        //remove any pixels on image border.  TODO: review why this is useful to do
        Set<PairInt> remove = new HashSet<PairInt>();
        for (PairInt p : borderPixels0) {
            if (p.getX() == 0 || p.getY() == 0 || p.getX() == imageMaxColumn ||
                p.getY() == imageMaxRow) {
                remove.add(p);
            }
        }
        
        for (PairInt p : remove) {
            borderPixels0.remove(p);
        }
        
        return borderPixels0;
    }

    /**
     * find the pixels which are the borders of rowColRanges including concave
     * pixels, and including any pixels that are the image border pixels.
     * Note that rowColRanges has to represent a contiguous point set.
     * runtime complexity is approx O(N) where N is row_range * col_range
     * @param rowColRanges the column bounds for each row of a contiguous
     * point set.
     * @param rowMinMax the minimum and maximum rows present in the contiguous
     * point set.
     * @param imageMaxColumn the maximum column in the image from which the
     * point set was derived.
     * @param imageMaxRow the maximum row in the image from which the point
     * set was derived.
     * @return
     */
    public Set<PairInt> getBorderPixels0(Map<Integer, List<PairInt>> rowColRanges,
        int[] rowMinMax, int imageMaxColumn, int imageMaxRow) {

        Set<PairInt> borderPixels = new HashSet<PairInt>();
        
        if (rowColRanges.isEmpty()) {
            return borderPixels;
        }

        if (debug) {
            algorithms.misc.MiscDebug.assertAllRowsPopulated(rowColRanges,
                rowMinMax, imageMaxColumn, imageMaxRow);
        }

        /*
        Need to handle concave bounds:
           [][][][][][][][]
           [][]    [][]  [][]
              [][][][][][]
             []     []
              [][]
           [][]  [][]  [][] <--- the 2nd to last point should be found as border too
              [][]  [][][]  <--- same for 3rd to last here and 2nd in row

          ___________________
          |[][][][][][][][]  |  Find the min and max column bounds for the
          |[][]    [][]  [][]|  region.  make a point set of all points that
          |   [][][][][][]   |  are not in a colRange.
          |  []     []       |  Find the contiguous groups among those points.
          |   [][]           |  Then iterate over the boundaries of the region
          |[][]  [][]  [][]  |  to test whether a point is found n a group,
          |   [][]  [][][]   |  and when it is, put it in a set called
          --------------------  connectedToBounds.
                                        Then iterate over each point in colRanges
        ----------------------------    and for each test if it is adjacent to
        |     ___________________       a point in connectedToBounds.
        |     |[][][][][][][][]  |      If the point is connected, put it in the
        |     |[][]    [][]  [][]|  border points set.
        |     |   [][][][][][]   |
        |     |  []     [][][]   |  Then, if the top row is > 0, add the top
        |     |   [][]           |  row colRange pixels to border points set.
        |     |[][]  []  [][][]  |
        |     |   [][]  [][][]   |  Then, if the bottow row is < max of image rows,
        |     --------------------  add the bottom row colRange pixels to
        |                           border points.
        ----------------------------
                                    Then if any colRanges are equal to the leftmost
                                    region bounds and that is > 0,
                                    those points should be added to border pixels set.

        Then if any colRanges are equal to the rightmost region bounds and that
        is < imageMaxColumn, those points should be added to border pixels set.

        */

        int[] colMinMax = getMinMaxColumnsInRanges(rowColRanges, rowMinMax);

        // ---- find nonMembers that are connected to the bounds of the -----
        // ---- region bounded by rowMinMax and colMinMax               -----

        Set<PairInt> nonMembersConnectedToBounds =
            findNonMembersConnectedToBounds(rowColRanges, rowMinMax, colMinMax,
                imageMaxColumn, imageMaxRow);

        //int[] dxs = new int[]{-1, -1,  0,  1, 1, 1, 0, -1};
        //int[] dys = new int[]{ 0, -1, -1, -1, 0, 1, 1,  1};
        int[] dxs = new int[]{-1, 0,  1, 0};
        int[] dys = new int[]{ 0, -1, 0, 1};

        // --- for each member in colRanges, if it's adjacent to a point
        // --- in nonMembersConnectedToBounds it's a border point
        for (int row = rowMinMax[0]; row <= rowMinMax[1]; ++row) {

            List<PairInt> colRanges = rowColRanges.get(Integer.valueOf(row));

            if (colRanges == null || colRanges.isEmpty()) {
                throw new IllegalStateException(
                "each row should have a point in it, else not contiguous");
            }

            for (PairInt colRange : colRanges) {
                for (int col = colRange.getX(); col <= colRange.getY(); ++col) {

                    // test if adjacent to a point in nonMembersConnectedToBounds
                    for (int idx = 0; idx < dxs.length; ++idx) {

                        int x = col + dxs[idx];
                        int y = row + dys[idx];

                        PairInt t = new PairInt(x, y);
                        if (nonMembersConnectedToBounds.contains(t)) {
                            borderPixels.add(new PairInt(col, row));
                            break;
                        }
                    }
                }
            }
        }

        // --- then add any pixels on the bounds
        int[] rows;
        if (rowMinMax[0] >= 0) {
            if (rowMinMax[1] <= imageMaxRow) {
                rows = new int[]{rowMinMax[0], rowMinMax[1]};
            } else {
                rows = new int[]{rowMinMax[0]};
            }
        } else if (rowMinMax[1] <= imageMaxRow) {
            rows = new int[]{rowMinMax[1]};
        } else {
            rows = new int[]{};
        }

        for (int row : rows) {

            List<PairInt> colRanges = rowColRanges.get(Integer.valueOf(row));

            if (colRanges == null || colRanges.isEmpty()) {
                throw new IllegalStateException(
                "each row should have a point in it, else not contiguous");
            }

            for (PairInt colRange : colRanges) {
                for (int col = colRange.getX(); col <= colRange.getY(); ++col) {
                    borderPixels.add(new PairInt(col, row));
                }
            }
        }

        for (int row = rowMinMax[0]; row <= rowMinMax[1]; ++row) {

            List<PairInt> colRanges = rowColRanges.get(Integer.valueOf(row));

            if (colRanges == null || colRanges.isEmpty()) {
                throw new IllegalStateException(
                "each row should have a point in it, else not contiguous");
            }

            int firstX = colRanges.get(0).getX();

            if ((colMinMax[0] >= 0) && (firstX >= colMinMax[0])) {
                borderPixels.add(new PairInt(firstX, row));
            }

            int n = colRanges.size();

            int lastX = colRanges.get(n - 1).getY();

            if ((colMinMax[1] <= imageMaxColumn) && (lastX <= colMinMax[1])) {
                borderPixels.add(new PairInt(lastX, row));
            }
        }

        return borderPixels;
    }

    /**
     * get a point set of the points not in column ranges for the region bounded
     * by min row, max row, and the minimum of columns and the maximum of columns.
     * runtime complexity is approx O(N) where N is n_row_range * n_col_range.
     * @param rowColRanges
     * @param rowMinMax
     * @param colMinMax
     * @param imageMaxColumn
     * @param imageMaxRow
     * @return
     */
    public Set<PairInt> getVoidsInRectangularRegion(Map<Integer, List<PairInt>> rowColRanges,
        int[] rowMinMax, int[] colMinMax, int imageMaxColumn, int imageMaxRow) {

        Set<PairInt> set = new HashSet<PairInt>();

        for (int row = rowMinMax[0]; row <= rowMinMax[1]; row++) {

            List<PairInt> colRanges = rowColRanges.get(Integer.valueOf(row));

            if ((colRanges == null) || colRanges.isEmpty()) {
                continue;
            }

            int n = colRanges.size();

            PairInt colRange = colRanges.get(0);

            for (int x = colMinMax[0]; x < colRange.getX(); ++x) {
                set.add(new PairInt(x, row));
            }

            for (int i = 1; i < n; ++i) {
                int lx = colRange.getY();
                colRange = colRanges.get(i);
                int rx = colRange.getX();
                for (int x = (lx + 1); x < rx; ++x) {
                    set.add(new PairInt(x, row));
                }
            }

            for (int x = (colRange.getY() + 1); x <= colMinMax[1]; ++x) {
                set.add(new PairInt(x, row));
            }
        }

        return set;
    }

    public int[] getMinMaxColumnsInRanges(Map<Integer, List<PairInt>> rowColRanges,
        int[] rowMinMax) {

        int minX = Integer.MAX_VALUE;
        int maxX = Integer.MIN_VALUE;

        for (int row = rowMinMax[0]; row <= rowMinMax[1]; row++) {

            List<PairInt> colRanges = rowColRanges.get(Integer.valueOf(row));

            if ((colRanges == null) || colRanges.isEmpty()) {
                continue;
            }

            int n = colRanges.size();

            int tc = colRanges.get(0).getX();

            if (tc < minX) {
                minX = tc;
            }
            tc = colRanges.get(n - 1).getY();
            if (tc > maxX) {
                maxX = tc;
            }
        }

        return new int[]{minX, maxX};
    }

    protected Gap findLastGap(Set<Gap> gaps) {
        Gap lastGap = null;
        for (Gap gap : gaps) {
            if (lastGap == null) {
                lastGap = gap;
            } else {
                if (gap.getStopInclusive() > lastGap.getStopInclusive()) {
                    lastGap = gap;
                }
            }
        }
        return lastGap;
    }

    protected Gap findFirstGap(Set<Gap> gaps) {
        Gap firstGap = null;
        for (Gap gap : gaps) {
            if (firstGap == null) {
                firstGap = gap;
            } else {
                if (gap.getStart() < firstGap.getStart()) {
                    firstGap = gap;
                }
            }
        }
        return firstGap;
    }

    protected boolean adjacentGapIsConnectedToImageBoundary(
        int startGap, int stopGapInclusive, Set<Gap> adjacentGaps,
        List<PairInt> adjacentColRanges, int imageMaxColumn) {

        if ((adjacentColRanges == null) || adjacentColRanges.isEmpty()) {
            return true;
        }

        /*
        check to see if the startGap:stopGapInclusive is adjacent to a
        row which has a leading or trailing gap which is connected to the image
        boundaries.
        */
        PairInt lastColRange = adjacentColRanges.get(adjacentColRanges.size() - 1);
        if ((adjacentGaps == null) || adjacentGaps.isEmpty()) {
            if (lastColRange.getY() < imageMaxColumn) {
                if (stopGapInclusive >= (lastColRange.getY() + 1)) {
                    return true;
                }
            }
        } else {
            Gap lastGap = findLastGap(adjacentGaps);
            if (lastColRange.getX() > lastGap.getStopInclusive()) {
                if (lastColRange.getY() < imageMaxColumn) {
                    if (stopGapInclusive >= (lastColRange.getY() + 1)) {
                        return true;
                    }
                }
            }
        }
        PairInt firstColRange = adjacentColRanges.get(0);
        if ((adjacentGaps == null) || adjacentGaps.isEmpty()) {
            if (firstColRange.getX() > 0) {
                if (startGap <= (firstColRange.getX() - 1)) {
                    return true;
                }
            }
        } else {
            Gap firstGap = findFirstGap(adjacentGaps);
            if (firstColRange.getY() < firstGap.getStart()) {
                if (firstColRange.getX() > 0) {
                    if (startGap <= (firstColRange.getX() - 1)) {
                        return true;
                    }
                }
            }
        }

        if ((adjacentGaps == null) || adjacentGaps.isEmpty()) {
            return false;
        }

        for (int col = startGap; col <= stopGapInclusive; col++) {
            for (Gap gap : adjacentGaps) {
                if ((col >= gap.getStart()) && (col <= gap.getStopInclusive())) {
                    if ((gap.getStart() == 0) || (gap.getStopInclusive() == imageMaxColumn)) {
                        return true;
                    }
                }
            }
        }

        return false;
    }

    protected Map<Integer, Set<Gap>> createRowMap(List<List<Gap>> gapLists) {

        Map<Integer, Set<Gap>> rowSetsMap = new HashMap<Integer, Set<Gap>>();

        for (List<Gap> gaps : gapLists) {

            for (Gap gap : gaps) {

                Integer row = Integer.valueOf(gap.getRow());

                Set<Gap> set = rowSetsMap.get(row);

                if (set == null) {
                    set = new HashSet<Gap>();
                    rowSetsMap.put(row, set);
                }

                set.add(gap);
            }
        }

        return rowSetsMap;
    }

    protected int findIndexOfOverlappingRange(List<PairInt> colRanges,
        PairInt findColRange) {

        if ((colRanges == null) || colRanges.isEmpty()) {
            return -1;
        }
        int fc0 = findColRange.getX();
        int fc1 = findColRange.getY();

        for (int i = 0; i < colRanges.size(); i++) {
            PairInt colRange = colRanges.get(i);
            int c0 = colRange.getX();
            int c1 = colRange.getY();

            if ((fc0 <= c0) && (fc1 >= c0)) {
                return i;
            } else if ((fc0 >= c0) && (fc0 <= c1)) {
                return i;
            }
        }
        return -1;
    }

    private int[] findMinMaxColumns(Map<Integer, List<PairInt>> rowColRanges,
        int[] rowMinMax) {

        // runtime is O(rowColRanges.size)
        int minY = rowMinMax[0];
        int maxY = rowMinMax[1];
        int minX = Integer.MAX_VALUE;
        int maxX = Integer.MIN_VALUE;
        for (int row = minY; row <= maxY; row++) {
            List<PairInt> colRanges = rowColRanges.get(Integer.valueOf(row));
            if ((colRanges == null) || colRanges.isEmpty()) {
                continue;
            }
            int n = colRanges.size();
            PairInt cr = colRanges.get(0);
            if (cr.getX() < minX) {
                minX = cr.getX();
            }
            if (n > 1) {
                cr = colRanges.get(n - 1);
            }
            if (cr.getY() > maxX) {
                maxX = cr.getY();
            }
        }

        return new int[]{minX, maxX};
    }

    private Map<Integer, List<PairInt>> copy(
        Map<Integer, List<PairInt>> rowColRanges, int[] rowMinMax) {

         Map<Integer, List<PairInt>> output = new HashMap<Integer, List<PairInt>>();
         for (int row = rowMinMax[0]; row <= rowMinMax[1]; ++row) {

             Integer key = Integer.valueOf(row);

             List<PairInt> colRanges = rowColRanges.get(key);

             List<PairInt> outputColRanges = new ArrayList<PairInt>();

             for (PairInt p : colRanges) {
                 outputColRanges.add(new PairInt(p.getX(), p.getY()));
             }

             output.put(key, outputColRanges);
         }

         return output;
    }

    /**
     * runtime complexity is approx O(N) where N is row_range * col_range
     * @param rowColRanges
     * @param rowMinMax
     * @param colMinMax
     * @param imageMaxColumn
     * @param imageMaxRow
     * @return 
     */
    public Set<PairInt> findNonMembersConnectedToBounds(
        Map<Integer, List<PairInt>> rowColRanges, int[] rowMinMax,
        int[] colMinMax, int imageMaxColumn, int imageMaxRow) {

        //approx O(N) where N is row_range * col_range
        Set<PairInt> nonMembers = getVoidsInRectangularRegion(rowColRanges,
            rowMinMax, colMinMax, imageMaxColumn, imageMaxRow);

        // approx O(N_non_members) 
        //Find the contiguous groups among nonMembers.
        DFSConnectedGroupsFinder contigFinder = new DFSConnectedGroupsFinder();
        contigFinder.setMinimumNumberInCluster(1);
        contigFinder.findConnectedPointGroups(nonMembers, imageMaxColumn,
            imageMaxRow);

        Set<Set<PairInt>> contigNonMembers = new HashSet<Set<PairInt>>();
        for (int i = 0; i < contigFinder.getNumberOfGroups(); ++i) {
            Set<PairInt> group = contigFinder.getXY(i);
            contigNonMembers.add(group);
        }

        Set<PairInt> contigNonMembersConnectedToBounds = new HashSet<PairInt>();

        // --- find the top and bottom row pixels not in colRanges and test memberships ---
        int[] rows = new int[]{rowMinMax[0], rowMinMax[1]};

        for (int row : rows) {

            List<PairInt> colRanges = rowColRanges.get(Integer.valueOf(row));

            if ((colRanges == null) || colRanges.isEmpty()) {
                continue;
            }

            int n = colRanges.size();

            PairInt colRange = colRanges.get(0);

            for (int x = colMinMax[0]; x < colRange.getX(); ++x) {
                PairInt t = new PairInt(x, row);
                Set<PairInt> keep = null;
                for (Set<PairInt> set : contigNonMembers) {
                    if (set.contains(t)) {
                        keep = set;
                        break;
                    }
                }
                if (keep != null) {
                    contigNonMembersConnectedToBounds.addAll(keep);
                    contigNonMembers.remove(keep);
                }
            }

            for (int i = 1; i < n; ++i) {
                int lx = colRange.getY();
                colRange = colRanges.get(i);
                int rx = colRange.getX();
                for (int x = (lx + 1); x < rx; ++x) {
                    PairInt t = new PairInt(x, row);
                    Set<PairInt> keep = null;
                    for (Set<PairInt> set : contigNonMembers) {
                        if (set.contains(t)) {
                            keep = set;
                            break;
                        }
                    }
                    if (keep != null) {
                        contigNonMembersConnectedToBounds.addAll(keep);
                        contigNonMembers.remove(keep);
                    }
                }
            }

            for (int x = (colRange.getY() + 1); x <= colMinMax[1]; ++x) {
                PairInt t = new PairInt(x, row);
                Set<PairInt> keep = null;
                for (Set<PairInt> set : contigNonMembers) {
                    if (set.contains(t)) {
                        keep = set;
                        break;
                    }
                }
                if (keep != null) {
                    contigNonMembersConnectedToBounds.addAll(keep);
                    contigNonMembers.remove(keep);
                }
            }
        }

        // --- scan the pixels in first and last columns and test membership ---
        // --- in contigNonMembers
        int[] cols = new int[]{colMinMax[0], colMinMax[1]};

        for (int col : cols) {

            for (int row = rowMinMax[0]; row <= rowMinMax[1]; ++row) {

                PairInt t = new PairInt(col, row);
                Set<PairInt> keep = null;
                for (Set<PairInt> set : contigNonMembers) {
                    if (set.contains(t)) {
                        keep = set;
                        break;
                    }
                }

                if (keep != null) {
                    contigNonMembersConnectedToBounds.addAll(keep);
                    contigNonMembers.remove(keep);
                }
            }
        }

        return contigNonMembersConnectedToBounds;
    }

    /* NOT READY FOR USE
     * Walk the border to find the perimeter border pixels.
     * Note that rowColRanges has to represent the contiguous point set called points.
     * Note also that this method removes spurs (single pixel width paths
     * without a turn around).
     * If wanted to keep the spurs and still have a connected sequential curve,
     * could copy this method and instead of back-tracking, add the previous
     * point in the linked list at the dead end and continue (for the spurs,
     * there will be the same points leading out and back composing the spur
     * in the final curve - special care is needed for starter points which
     * are spurs and storing "visited" for that altered algorithm).
     *
     * @param rowColRanges the column bounds for each row of a contiguous
     * point set.
     * @param rowMinMax the minimum and maximum rows present in the contiguous
     * point set.
     * @param imageMaxColumn the maximum column in the image from which the
     * point set was derived.
     * @param imageMaxRow the maximum row in the image from which the point
     * set was derived.
     * @return
    public PairIntArray getOrderedBorderPixels(Set<PairInt> points,
        Map<Integer, List<PairInt>> rowColRanges,
        int[] rowMinMax, int imageMaxColumn, int imageMaxRow) {

        PerimeterFinderMemo memo = new PerimeterFinderMemo();

        BitVectorRepresentation bitVectorHelper = new BitVectorRepresentation(points);

        //the algorithm uses rowColRanges to make a starting point and add it
        //to the path of steps.
        //
        //for each currentStep, the state at the time of instantiation is stored
        //in it's node.  the state is it's pixel location, a summary of all visited
        //nodes at that moment, and the available moves for the pixel location.
        //
        //it's available moves are adjacent points that have at least one non-point
        //neighbor (the goal is to find border points, so has to be on perimeter
        //of points) and have not been visited yet.
        //
        //at each iteration, the currentStep is first checked against memo's list
        //of known "non-solution states".  if that exists, it backtracks (without
        //storing backtrack state again).
        //If the currentStep was not found by memo, its nextSteps are examined to
        //choose the next step if any is possible, else backtrack if none are
        //possible. (Note that backtracking stores the entire path in memo to
        //capture the state of nodes that led to "no solution" and then removes
        //currentStep from previous step's nextSteps, removes currentStep from visited
        //and border,
        //then removes currentStep from path and sets currentStep to last node in path).
        //(Note, that the responsibility of updating a step's nextNodes to remove
        //a point should only occur at backtrack stage after memoization in order
        //to allow correct recognition of state later to avoid "non-solution" states).
        // 
        //If a next step is possible (did not backtrack), a PathStep is created
        //for the nextStep and it is added to the path, border, and visited.
        //
        //then iteration continues until an end state is reached.
        //end state is currentStep is adjacent to starter point after have left
        //the start region, or all points have been visited or all paths tried <--- review this

        //example for memoization and backtracking.
        //
        //In this example, can see that when the path reaches node 20, needs to
        //backtrack to node 1 and choose 1->4 instead of 1->2->3->4.
        //
        //Can see that nodes 5 to 18 are outside of an interaction distance with
        //path nodes 1 through 4, so then, the paths from 5 to 18
        //can be saved and restored if needed later instead of re-computing them 
        //when '4' is next selected as the next node in a path.  the "next move" 
        //after '4' would skip forward to '18' in a choice of one of the saved 
        //paths from 5 to 18.
        //
        //have not implemented that below, but have implemented memoization
        //for states that lead to no-solution, and avoidance of continuing on
        //those paths in the future.
        //
        //see PerimeterFinderMemo.java, in progress...
        //
        //Reuse of successful paths is needed for this algorithm if the point sets
        //are larger than a couple dozen points, UNLESS the method is given 
        //points that are already an edge thinned to single pixel widths at all
        //locations, including the diagonals.
        //
        //55
        //54                     11   12
        //53            9   10   14   13
        //52            8   15
        //51            7   16
        //50            6   17
        //49            5   20   18
        //48            4   19
        //47       1    2    3
        //46       0
        //   125  126  127  128  129  130

        PairIntArray border = new PairIntArray();

        if (rowColRanges.isEmpty()) {
            return border;
        }

        if (points.size() < 5) {
            throw new IllegalArgumentException(
            "points.size() needs to be a least 5 points");
        }

        if (debug) {
            algorithms.misc.MiscDebug.assertAllRowsPopulated(rowColRanges,
                rowMinMax, imageMaxColumn, imageMaxRow);
        }

        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();

        //xc,yc
        double[] xyCentroid = curveHelper.calculateXYCentroids(points);

        Set<PairInt> visited = new HashSet<PairInt>();
        Set<PairInt> neighbors = new HashSet<PairInt>();

        // ------- initialize with step 0 ------
        //step 0 should be the smallest row and col,
        //but because it may be a single pixel width spur (which will not end in
        //a closed curve), have to make sure that the added point has at least
        //2 neighbors, else discard it and continue forward until have found a
        //starter point with 2 neighbors.

        int y0 = rowMinMax[0];
        List<PairInt> colRanges = rowColRanges.get(Integer.valueOf(y0));
        PairInt colRange = colRanges.get(0);
        int x0 = colRange.getX();

        //storing the state of choices in a linked list.
        //   each choice holds the col,row of that pixel and the neighbor points
        //      for paths not chosen.
        //backtracking removes last node and chooses from next to last node
        
        LinkedList<PathStep> path = new LinkedList<PathStep>();
        PathStep currentStep = null;

        curveHelper.findNeighborsWithAtLeastOneNonPoint(x0, y0, neighbors,
            points, visited, imageMaxColumn + 1, imageMaxRow + 1);

        if (neighbors.size() < 2) {

            while (neighbors.size() < 2) {
                PairInt p = neighbors.iterator().next();
                x0 = p.getX();
                y0 = p.getY();
                curveHelper.findNeighborsWithAtLeastOneNonPoint(x0, y0,
                    neighbors, points, visited,
                    imageMaxColumn + 1, imageMaxRow + 1);
            }
            border.add(x0, y0);
            visited.add(new PairInt(x0, y0));

            currentStep = new PathStep(new PairInt(x0, y0), neighbors,
                bitVectorHelper.createBitstring(visited));

            path.add(currentStep);

            log.fine("start point=(" + x0 + "," + y0 + ")");

        } else {

            border.add(x0, y0);
            visited.add(new PairInt(x0, y0));

            currentStep = new PathStep(new PairInt(x0, y0), neighbors,
                bitVectorHelper.createBitstring(visited));

            path.add(currentStep);

            log.fine("start point=(" + x0 + "," + y0 + ")");

            // --- take the 2nd step if colRange.getY() is > col, to make start
            //     direction clockwise ---

            if (colRange.getY() > colRange.getX()) {

                PairInt p = new PairInt(x0 + 1, y0);
                border.add(x0 + 1, y0);
                visited.add(p);
                log.fine("adding (" + (x0 + 1) + "," + y0 + ")");

                curveHelper.findNeighborsWithAtLeastOneNonPoint(x0 + 1, y0, 
                    neighbors, points, visited,
                    imageMaxColumn + 1, imageMaxRow + 1);

                currentStep = new PathStep(p, neighbors,
                    bitVectorHelper.createBitstring(visited));

                path.add(currentStep);
            }
        }

        // consider each move in all 8 directions in which move point is
        //   !visited.constains(point) && points.contains(point)
        //
        //if more than one choice, prefer choice with fewer neighbors.
        //
        //if there are zero choices, remove the latest point from border and try again.
        //
        //repeat while col,row is not adjacent to x0,y0 and visited.size < points.size. <=== reviewing this when testing

        boolean search = true;

        int earliestBacktrack = Integer.MAX_VALUE;

        while (search) {

            boolean currentIsNoSolution = memo.isANoSolutionState(currentStep);
            
            if (currentStep.getNextSteps().isEmpty() || currentIsNoSolution) {

                while (currentStep.getNextSteps().isEmpty() || currentIsNoSolution) {

                    log.fine("   removing (" + currentStep.getCoords().getX() + ","
                        + currentStep.getCoords().getY() + ")"
                        + " border.n=" + border.getN() + " path.size=" + path.size()
                    );

                    if (border.getN() == 0) {
                        return null;
                    }

                    // no solution for this path, so store state and backtrack

                    if (!currentIsNoSolution) {
                        // if currentIsNoSolution, state has already been stored
                        memo.storeNoSolutionState(path);
                    }

                    path.removeLast();
                    
                    PathStep nextStep = path.peekLast();

                    if (nextStep != null) {
                        nextStep.getNextSteps().remove(currentStep.getCoords());
                    }

                    border.removeRange(border.getN() - 1, border.getN() - 1);

                    visited.remove(currentStep.getCoords());

                    if (path.size() < earliestBacktrack) {
                        earliestBacktrack = path.size();
                        log.info("earliestBacktrack=" + earliestBacktrack);
                    }

                    currentStep = nextStep;

                    if (currentStep == null) {
                        return border;
                    }

                    currentIsNoSolution = memo.isANoSolutionState(currentStep);
                }

            } else if (currentStep.getNextSteps().size() == 1) {

                PairInt p = currentStep.getNextSteps().iterator().next();
                border.add(p.getX(), p.getY());
                visited.add(p);
                log.fine("adding (" + p.getX() + "," + p.getY() + ")"
                    + " border.n=" + border.getN() + " path.size=" + path.size()
                );

                curveHelper.findNeighborsWithAtLeastOneNonPoint(p.getX(), 
                    p.getY(), neighbors, points,
                    visited, imageMaxColumn + 1, imageMaxRow + 1);

                currentStep = new PathStep(p, neighbors,
                    bitVectorHelper.createBitstring(visited));

                path.add(currentStep);

            } else {

                // choose the neighbor with the fewest neighbors and when there
                // are ties, choose the one furthest from the centroid

                PairInt p = null;
                int minNNeigbors = Integer.MAX_VALUE;

                for (PairInt pn : currentStep.getNextSteps()) {

                    int x2 = pn.getX();
                    int y2 = pn.getY();

                    log.fine("   compare " + " (" + x2 + "," + y2 + ")  to"
                        + " (" + currentStep.getCoords().getX() + ","
                        + currentStep.getCoords().getY() + ")");

                    int nNeighbors = curveHelper.countNeighbors(x2, y2, points,
                        imageMaxColumn + 1, imageMaxRow + 1);

                    if (p == null) {

                        minNNeigbors = nNeighbors;
                        p = pn;

                    } else if (nNeighbors < minNNeigbors) {

                        minNNeigbors = nNeighbors;
                        p = pn;

                    } else if (nNeighbors == minNNeigbors) {

                        double distC = Math.sqrt(
                            (p.getX() - xyCentroid[0])*(p.getX() - xyCentroid[0])
                            +
                            (p.getY() - xyCentroid[1])*(p.getY() - xyCentroid[1]));

                        double dist2 = Math.sqrt(
                            (x2 - xyCentroid[0])*(x2 - xyCentroid[0])
                            +
                            (y2 - xyCentroid[1])*(y2 - xyCentroid[1]));

                        if (dist2 > distC) {
                            p = pn;
                        }
                    }
                }

                border.add(p.getX(), p.getY());
                visited.add(p);
                log.fine("adding (" + p.getX() + "," + p.getY() + ")"
                    + " border.n=" + border.getN() + " path.size=" + path.size()
                );

                curveHelper.findNeighborsWithAtLeastOneNonPoint(p.getX(), 
                    p.getY(), neighbors, points,
                    visited, imageMaxColumn + 1, imageMaxRow + 1);

                currentStep = new PathStep(p, neighbors,
                    bitVectorHelper.createBitstring(visited));

                path.add(currentStep);

            }

            // only eval if have left the region of start point.
            //TODO: refine this limit for very small point sets (e.g. 9).
            //TODO: need a check on the innermost block that the points are
            //      bounded, excepting spurs which should be stored. do not want
            //      a path which turns around early, missing most of the points.
            //      can use point in polygon test...or an area test or centroid test...
            if (border.getN() > 8) {
                if ((Math.abs(currentStep.getCoords().getX() - x0) < 2)
                    && (Math.abs(currentStep.getCoords().getY() - y0) < 2)) {
                    search = false;
                }
            }
        }

        return border;
    }
    */
    
    /**
     * given the perimeter of points bounding the blob and an expected distance
     * of srchRadius between points, order the points in counter clockwise
     * manner (the CCW is consistent w/ earlier curve orderings in contour matcher).
     * approx O(N_perimeter), but has factors during searches that could be improved
     * @param perimeter
     * @param blob
     * @param srchRadius expected distance between points in the perimeter.
     * For example, a densely populated blob would have a perimeter with 
     * srchRadius of Math.sqrt(2) to include diagonals.
     * @return 
     */
    public PairIntArray orderThePerimeter(Set<PairInt> perimeter, 
        Set<PairInt> blob, float srchRadius) {
        
        ZhangSuenLineThinner lt = new ZhangSuenLineThinner();
        
        Set<PairInt> skeleton = new HashSet<PairInt>(blob);
                        
        int[] minMaxXY = MiscMath.findMinMaxXY(skeleton);
        lt.applyLineThinner(skeleton, minMaxXY[0], minMaxXY[1], minMaxXY[2],
            minMaxXY[3]);
        int[] xPoints = new int[skeleton.size()];
        int[] yPoints = new int[skeleton.size()];

        int count = 0;
        for (PairInt p : skeleton) {
            xPoints[count] = p.getX();
            yPoints[count] = p.getY();
            count++;
        }
            
        // order skeleton by x and then y
        Map<Integer, List<Integer>> xSkeletonMap = makeXMap(xPoints, yPoints);
        Map<Integer, List<Integer>> ySkeletonMap = makeYMap(xPoints, yPoints);

        // approx O(N_perimeter), but has some factors during searches that could be improved
        // dfs walk through points, adding the clockwise point
        PairIntArray orderedEdge = orderThePerimeter(perimeter, 
            xSkeletonMap, ySkeletonMap, srchRadius);
        
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        
        boolean isAdjacent = curveHelper.isAdjacent(orderedEdge, 0, 
            orderedEdge.getN() - 1);
        
        if (isAdjacent) {
            PairIntArrayWithColor p = new PairIntArrayWithColor(orderedEdge);
            p.setColor(1);
            orderedEdge = p;
        }
            
        return orderedEdge;
    }
    
    private PairIntArray orderThePerimeter(Set<PairInt> perimeter, 
        Map<Integer, List<Integer>> xSkeletonMap, 
        Map<Integer, List<Integer>> ySkeletonMap, float srchRadius) {
        
        /*
        TODO: need better data structures for finding nearest skeleton point
        and need to improve NearestPoints.
        TODO: the code below should be simplified if possible.
        */
        
        int[] xPoints = new int[perimeter.size()];
        int[] yPoints = new int[xPoints.length];
        int count = 0;
        for (PairInt p : perimeter) {
            xPoints[count] = p.getX();
            yPoints[count] = p.getY();
            count++;
        }
        NearestPoints np = new NearestPoints(xPoints, yPoints);
                
        PairInt prev = np.getSmallestXY();

        PairIntArray orderedEdge = new PairIntArray(perimeter.size());
        
        int n = perimeter.size();
        
        orderedEdge.add(prev.getX(), prev.getY());
        
        while (!perimeter.isEmpty()) {
            
            boolean rmvd = perimeter.remove(prev);
            assert(rmvd == true);
            
            if (perimeter.isEmpty()) {
                break;
            }
            
            int x = prev.getX();
            int y = prev.getY();
         
            Set<PairInt> neighborsC = np.findNeighbors(x, y, srchRadius);
            Set<PairInt> neighbors = new HashSet<PairInt>();
            for (PairInt p : neighborsC) {
                if (perimeter.contains(p) && !p.equals(prev)) {
                    neighbors.add(p);
                }
            }
            
            if (neighbors.isEmpty()) {
                // last point may have been a spur so next point is a multiple
                // of search radius away
                int nIterMax = n;
                int nIter = 0;
                while (neighbors.isEmpty() && (nIter < nIterMax)) {
                    neighborsC = np.findNeighbors(x, y, (nIter + 2)*srchRadius);
                    for (PairInt p : neighborsC) {
                        if (perimeter.contains(p) && !p.equals(prev)) {
                            neighbors.add(p);
                        }
                    }
                    nIter++;
                }
                if (neighbors.isEmpty()) {
                    // some corners may have been removed so check whether
                    // is adjacent to start point
                    if (orderedEdge.getN() > 1) {
                        int diffX = orderedEdge.getX(orderedEdge.getN() - 2) - x;
                        int diffY = orderedEdge.getY(orderedEdge.getN() - 2) - y;
                        if (Math.sqrt(diffX*diffX + diffY*diffY) < srchRadius) {
                            break;
                        }
                    }
                    if (perimeter.size() == 1) {
                        log.info("remaining point in perimeters=" + perimeter.iterator().next());
                    }
                    String err = "1) no neighbors for point (" + x + "," + y + ").  "
                        + " remaining number of points to add is "
                        + (perimeter.size() - 1);
                    err = err + "\norderedEdge=" + orderedEdge.toString();
                    throw new IllegalStateException(err);
                }
            }
             
            if (neighbors.size() > 1) {
                
                // find the neighbor closest skeleton point then find the 
                // neighbor closest to it in counter clockwise direction.
                
                List<Integer> ys = xSkeletonMap.get(Integer.valueOf(x));
                int ySkel1 = Integer.MAX_VALUE;
                if (ys != null) {
                    int minDistY2 = Integer.MAX_VALUE;
                    for (Integer y0 : ys) {
                        int distYSq = y0.intValue() - y;
                        distYSq *= distYSq;
                        if (distYSq < minDistY2) {
                            minDistY2 = distYSq;
                            ySkel1 = y0.intValue();
                        }
                    }
                }
                int xSkel1 = x;
                List<Integer> xs = ySkeletonMap.get(Integer.valueOf(x));
                int xSkel2 = Integer.MAX_VALUE;
                if (xs != null) {
                    int minDistX2 = Integer.MAX_VALUE;
                    for (Integer x0 : xs) {
                        int distXSq = x0.intValue() - x;
                        distXSq *= distXSq;
                        if (distXSq < minDistX2) {
                            minDistX2 = distXSq;
                            xSkel2 = x0.intValue();
                        }
                    }
                }
                int ySkel2 = y;
                int xSkel, ySkel;
                if ((((xSkel1 - x)*(xSkel1 - x)) + ((ySkel1 - y)*(ySkel1 - y)))
                    < 
                    (((xSkel2 - x)*(xSkel2 - x)) + ((ySkel2 - y)*(ySkel2 - y)))) {
                    xSkel = xSkel1;
                    ySkel = ySkel1;
                } else {
                    xSkel = xSkel2;
                    ySkel = ySkel2;
                }
                
                /*
                can use the skeleton to determine orientation, and hence, the
                point in neighbors which is counter clock wise to the current
                point.
                
                   (xi, yi)
                    /
                 (x,y) --- (xSkel, ySkel)
                    \
                    (xi, yi)
                     
                */
                
                if ((x == xSkel) && (y == ySkel)) {
                    // compare angles w/ previous to find the one which is most CCW
                    int minDistSq = Integer.MAX_VALUE;
                    PairInt pNext = null;
                    double minDirection = Double.MIN_VALUE;// should be a (+) number
                    for (PairInt p : neighbors) {
                        int diffX = p.getX() - xSkel;
                        int diffY = p.getY() - ySkel;
                        int distSq = diffX * diffX + diffY * diffY;
                        if (pNext == null) {
                            pNext = p;
                            minDistSq = distSq;
                            continue;
                        }
                        double direction = LinesAndAngles.direction(x, y, 
                            pNext.getX(), pNext.getY(), p.getX(), p.getY());
                        if (direction == 0) {
                            if (distSq < minDistSq) {
                                minDistSq = distSq;
                                pNext = p;
                                minDirection = direction;
                            } else if (distSq == minDistSq && (orderedEdge.getN() > 1)) {
                                // use direction w.r.t. prev prev point
                                double direction2 = LinesAndAngles.direction(
                                    orderedEdge.getX(orderedEdge.getN() - 2),
                                    orderedEdge.getY(orderedEdge.getN() - 2),
                                    x, y, p.getX(), p.getY());
                                if (direction2 <= 0) {
                                    continue;
                                }
                                minDistSq = distSq;
                                pNext = p;
                                minDirection = direction;
                            }
                        } else if (direction > 0) {
                            pNext = p;
                            minDirection = direction;
                            minDistSq = distSq;
                        }
                    }
                    if (pNext == null) {
                        // some corners may have been removed so check whether
                        // is adjacent to start point
                        if (orderedEdge.getN() > 1) {
                            int diffX = orderedEdge.getX(orderedEdge.getN() - 2) - x;
                            int diffY = orderedEdge.getY(orderedEdge.getN() - 2) - y;
                            if (Math.sqrt(diffX*diffX + diffY*diffY) < srchRadius) {
                                break;
                            }
                        }
                        String err = "2) no neighbors for point (" + x + "," + y + ").  "
                            + " remaining number of points to add is "
                            + (perimeter.size() - 1);
                        err = err + "\norderedEdge=" + orderedEdge.toString();
                        throw new IllegalStateException(err);
                    }
                    prev = pNext;
                } else {
                    // use direction as given by cross product
                    
                    // only keep those with direction >= 0
                    int minDistSq = Integer.MAX_VALUE;
                    PairInt pNext = null;
                    double minDirection = Double.MIN_VALUE;// should be a (+) number
                    for (PairInt p : neighbors) {
                        double direction = LinesAndAngles.direction(x, y, xSkel, 
                            ySkel, p.getX(), p.getY());
                        double direction2 = Double.MIN_VALUE;
                        if (pNext != null) {
                            direction2 = LinesAndAngles.direction(x, y,
                                pNext.getX(), pNext.getY(), p.getX(), p.getY());
                            if (direction2 < 0) {
                                continue;
                            }
                        }
                        if (direction < 0) {
                            continue;
                        }
                        int diffX = p.getX() - xSkel;
                        int diffY = p.getY() - ySkel;
                        int distSq = diffX * diffX + diffY * diffY;
                        if (pNext == null) {
                            pNext = p;
                            minDistSq = distSq;
                            minDirection = direction;
                            continue;
                        }       
                        if (direction > minDirection) {
                            minDistSq = distSq;
                            pNext = p;
                            minDirection = direction;
                        } else if (direction == minDirection || direction == 0) {
                            if (direction2 <= 0) {
                                continue;
                            } else {
                                minDistSq = distSq;
                                pNext = p;
                                minDirection = direction;
                            }
                        }
                    }
                    if (pNext == null) {
                        // some corners may have been removed so check whether
                        // is adjacent to start point
                        if (orderedEdge.getN() > 1) {
                            int diffX = orderedEdge.getX(orderedEdge.getN() - 2) - x;
                            int diffY = orderedEdge.getY(orderedEdge.getN() - 2) - y;
                            if (Math.sqrt(diffX*diffX + diffY*diffY) < srchRadius) {
                                break;
                            }
                        }
                        String err = "3) no neighbors for point (" + x + "," + y + ").  "
                            + " remaining number of points to add is "
                            + (perimeter.size() - 1);
                        err = err + "\norderedEdge=" + orderedEdge.toString();
                        throw new IllegalStateException(err);
                    }
                    prev = pNext;
                }
            } else {
                prev = neighbors.iterator().next();
            }
            
            orderedEdge.add(prev.getX(), prev.getY());
        }
        
        return orderedEdge;
    }
    
    /**
     * map w/ key = x and y = ascending ordered values for that x
     * @param x
     * @param y
     * @return 
     */
    private Map<Integer, List<Integer>> makeXMap(int[] x, int[] y) {
        
        Map<Integer, List<Integer>> map = new HashMap<Integer, List<Integer>>();
        for (int i = 0; i < x.length; ++i) {
            Integer key = Integer.valueOf(x[i]);
            List<Integer> list = map.get(key);
            if (list == null) {
                list = new ArrayList<Integer>();
                map.put(key, list);
            }
            list.add(Integer.valueOf(y[i]));
        }
        
        for (Map.Entry<Integer, List<Integer>> entry : map.entrySet()) {
            List<Integer> ys = entry.getValue();
            if (ys.size() > 1) {
                Collections.sort(ys);
            }
        }
        
        return map;
    }
    
    /**
     * map w/ key = y and x = ascending ordered values for that y
     * @param x
     * @param y
     * @return 
     */
    private Map<Integer, List<Integer>> makeYMap(int[] x, int[] y) {
        
        Map<Integer, List<Integer>> map = new HashMap<Integer, List<Integer>>();
        for (int i = 0; i < y.length; ++i) {
            Integer key = Integer.valueOf(y[i]);
            List<Integer> list = map.get(key);
            if (list == null) {
                list = new ArrayList<Integer>();
                map.put(key, list);
            }
            list.add(Integer.valueOf(x[i]));
        }
        
        for (Map.Entry<Integer, List<Integer>> entry : map.entrySet()) {
            List<Integer> xs = entry.getValue();
            if (xs.size() > 1) {
                Collections.sort(xs);
            }
        }
        
        return map;
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
