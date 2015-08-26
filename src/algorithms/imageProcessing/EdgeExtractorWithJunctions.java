package algorithms.imageProcessing;

import algorithms.CountingSort;
import algorithms.MultiArrayMergeSort;
import algorithms.QuickSort;
import algorithms.Rotate;
import algorithms.imageProcessing.util.MatrixUtil;
import algorithms.misc.Misc;
import algorithms.misc.MiscDebug;
import algorithms.util.PairIntArray;
import algorithms.util.PairInt;
import algorithms.util.PairIntArrayComparator;
import algorithms.util.PairIntArrayDescendingComparator;
import algorithms.util.PairIntArrayWithColor;
import algorithms.util.PointPairInt;
import java.io.IOException;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.Stack;
import java.util.logging.Level;
import thirdparty.HungarianAlgorithm;

/**
 * Edge extractor operates on an image that has already been reduced to
 * single pixel width lines and extracts edges from it, attempting to make
 * the longest edges it can - this specific edge extractor stores and uses
 * junction information, so is not as fast as EdgeExtractor.java.
 *
 Edge extraction
    Local Methods:
        (1) DFS walk through adjacent non-zero pixels to form a sequence of
            pixels called an edge.

        (2) find join points

        (3) join edges using join points

        (4) find junction points

        (5) using the junction points, splice together the edges to improve
            the edges.  Currently, improve means to make the edges longer
            (which is good for an image expected to have one contour for
             example, and doesn't harm corner finding for non-contour goals).

        (7) remove edges shorter than a minimum length (this may of may not
            be currently enabled)

  @see AbstractEdgeExtractor
  *
 * @author nichole
 */
public class EdgeExtractorWithJunctions extends AbstractEdgeExtractor {

    private boolean debug = false;

    /**
     * map with key = center of junction pixel coordinates;
     * value = set of adjacent pixels when there are more than the preceding
     * and next.
     */
    private Map<Integer, Set<Integer>> junctionMap = new HashMap<Integer, Set<Integer>>();

    /**
     * map with key = pixel coordinates of all pixels involved in junctions;
     * value = PairInt holding index of edge that pixel is located in and
     * holding the index within that edge of the pixel.
     * for example, a pixel located in edges(0) at offset=100
     * would have PairInt(0, 100) as a value.
     */
    private Map<Integer, PairInt> junctionLocationMap = new HashMap<Integer, PairInt>();

    private int defaultMaxIterJunctionJoin = 1;

    private boolean singleClosedEdge = false;

    // ordered so that closest are first, then the diagonals
    private static int[] dxs8 = new int[]{-1, 0, 1,  0, -1,   1, 1, -1};
    private static int[] dys8 = new int[]{ 0, 1, 0, -1, -1,  -1, 1,  1};

    /**
     * NOTE:  input should have a black (empty) background and edges should
     * have values > 125 counts.  Edges should also have width of 1 and no larger.
     *
     * @param input
     */
    public EdgeExtractorWithJunctions(GreyscaleImage input) {

        super(input);
    }

    /**
     * NOTE:  input should have a black (empty) background and edges should
     * have values > 125 counts.  Edges should also have width of 1 and no larger.
     * The guide image is used to alter the extracted edges back towards the
     * highest intensity pixels of the guide image.  The guide image is expected
     * to be the combined X and Y gradient image from earlier processing
     * stages.
     *
     * @param input
     * @param anEdgeGuideImage
     */
    public EdgeExtractorWithJunctions(GreyscaleImage input,
        final GreyscaleImage anEdgeGuideImage) {

        super(input, anEdgeGuideImage);
    }

    public void overrideMaxNumberIterationsJunctionSplice(int nMaxIter) {
        if (nMaxIter < 2) {
            return;
        }
        defaultMaxIterJunctionJoin = nMaxIter;
    }

    /**
     * get a hashmap of the junctions within the edges, that is the regions
     * where there are more than 2 adjacent pixels to a pixel.
     * @return hashmap with key = image pixel indexes, values = set of adjacent
     * pixels as pixel indexes.
     */
    public Map<Integer, Set<Integer> > getJunctionMap() {
        return junctionMap;
    }

    /**
     * contains information needed to find any of the pixels that are part of
     * a junction.
     * @return hashmap with key = image pixel indexes, values = a PairInt with
     * x holding the edge index and y holding the offset index of the point
     * within the edge.
     */
    public Map<Integer, PairInt> getLocatorForJunctionAssociatedPoints() {
        return junctionLocationMap;
    }

    /**
     * find the edges and return as a list of points.  The method uses a
     * DFS search through all points in the image with values > 0 to link
     * adjacent sequential points into edges.
     * As a side effect, the method also populates
     * member variables edgeJunctionMap and outputIndexLocatorForJunctionPoints.
     *
     * @return
     */
    @Override
    public List<PairIntArray> findEdgesIntermediateSteps(List<PairIntArray> edges) {

        List<PairIntArray> output = edges;

        //O(N)
        Map<PairInt, PairInt> joinPoints = findJoinPoints(output);

        if (debug) {
            algorithms.misc.MiscDebug.printJoinPoints(joinPoints, output);
            algorithms.misc.MiscDebug.writeJoinPointsImage(joinPoints, output,
                img);
        }

        log.fine("edges.size()=" + output.size() + " before join-points");

        output = joinOnJoinPoints(joinPoints, output);

        log.fine("edges.size()=" + output.size() + " after join-points");

        final int nMaxIter = defaultMaxIterJunctionJoin;
        int nIter = 0;
        int nSplices = 0;

        while ((nIter == 0) || ((nIter < nMaxIter) && (nSplices > 0))) {

            findJunctions(output);

        if (output.size() > 10000) {
            return output;
        }

            if (debug) {
                algorithms.misc.MiscDebug.printJunctions(this.junctionMap,
                    output, img);
            }

            //TODO: could improve thisfor nIter > 0 by only visiting the
            // edges which had changed (been spliced and fused on previous iter)

            nSplices = spliceEdgesAtJunctionsIfImproves(output);

            if (nSplices > 0) {
                removeEdgesShorterThan(output, 1);
            }

            ++nIter;
        }

        findJunctions(output);

        log.info("nIter=" + nIter);

        if (!this.singleClosedEdge) {
            return output;
        }

        nIter = 0;
        nSplices = 0;

        while ((nIter == 0) || ((nIter < nMaxIter) && (nSplices > 0))) {

            if (nSplices > 0) {

                removeEdgesShorterThan(output, 1);

                //TODO: this is handled in the methods invoked now so consider reducing number of times used
                findJunctions(output);

                /*
                for edges that are not closed curves, try a DFS w/ preference
                for distance, and if the result is a closed curve, replace
                */
                replaceOpenCurvesIfPossibleToClose(output);
            }

            nSplices = 0;

            int nIns = insertAdjacentForClosedCurve(output);

            //nIns += insertBetweenAdjacentForClosedCurve(output);

            nSplices += nIns;

            ++nIter;
        }

        log.info("nIter=" + nIter);

        return output;
    }

    /**
     * Find the join points of edges which are adjacent to one another and to
     * no other edges.  This method is intended to follow the first creation
     * of edges from the line thinned edge image to find the unambiguously
     * connectable edges.  (Note that the method findJunctions(...) is later
     * used to find where there are more than 2 adjacent pixels to any pixel.)
     *
     * runtime complexity is O(N).
     *
     * @param edges
     * @return hashmap with key = edge index and index within edge of pixel;
     * value = the adjacent pixel's edge index and index within its edge
     */
    protected Map<PairInt, PairInt> findJoinPoints(List<PairIntArray> edges) {

        int n = edges.size();

        int w = img.getWidth();
        int h = img.getHeight();

        //with key = x, y of point
        //with value = edge index, index within edge
        Map<PairInt, PairInt> endPointMap = new HashMap<PairInt, PairInt>(2 * n);

        // holds x,y points of first and last points of the edges
        ArrayDeque<PairInt> endPointQueue = new ArrayDeque<PairInt>(2 * n);

        //with key = edge index, index within edge of one join point
        //                 (the join point with the smaller edge number)
        //       with value = edge index, index within edge of other join point
        Set<PointPairInt> theJoinPoints = new HashSet<PointPairInt>(2 * n);

        //with key = edge number
        //with value = set of the entries in theJoinPoints for which one has a
        //first end point in this edge
        Map<Integer, Set<PointPairInt>> edgeFirstEndPointMap =
            new HashMap<Integer, Set<PointPairInt>>();

        Map<Integer, Set<PointPairInt>>  edgeLastEndPointMap =
            new HashMap<Integer, Set<PointPairInt>> ();

        // 2* O(N_edges)
        for (int edgeIdx = 0; edgeIdx < n; edgeIdx++) {

            PairIntArray edge = edges.get(edgeIdx);

            int nEdge = edge.getN();

            if (nEdge == 0) {
                continue;
            }

            PairInt xy = new PairInt(edge.getX(0), edge.getY(0));
            PairInt loc = new PairInt(edgeIdx, 0);

            endPointMap.put(xy, loc);
            endPointQueue.add(xy);

            if (nEdge == 1) {
                continue;
            }

            xy = new PairInt(edge.getX(nEdge - 1), edge.getY(nEdge - 1));
            loc = new PairInt(edgeIdx, nEdge - 1);

            endPointMap.put(xy, loc);
            endPointQueue.add(xy);

            if (nEdge == 2) {
                continue;
            }

            // add the second and second to last to endPointMap only
            xy = new PairInt(edge.getX(1), edge.getY(1));
            loc = new PairInt(edgeIdx, 1);
            endPointMap.put(xy, loc);

            if (nEdge == 3) {
                continue;
            }

            xy = new PairInt(edge.getX(nEdge - 2), edge.getY(nEdge - 2));
            loc = new PairInt(edgeIdx, nEdge - 2);
            endPointMap.put(xy, loc);
        }

        assert(endPointQueue.size() >= n);
        assert(endPointMap.size() >= 2*n);

        while (!endPointQueue.isEmpty()) {

            PairInt uXY = endPointQueue.poll();
            int uX = uXY.getX();
            int uY = uXY.getY();

            PairInt uLoc = endPointMap.get(uXY);

            assert(uLoc != null);

            int uEdgeIdx = uLoc.getX();

            Set<Integer> neighborEdgeNumbers = new HashSet<Integer>();

            PairInt vClosestXY = null;
            PairInt vClosestLoc = null;
            int closestDistSq = Integer.MAX_VALUE;
            int vClosestCanBeReordered = -99;

            for (int nIdx = 0; nIdx < dxs8.length; nIdx++) {
                int vX = uX + dxs8[nIdx];
                int vY = uY + dys8[nIdx];
                if ((vX < 0) || (vX > (w - 1)) || (vY < 0) || (vY > (h - 1))) {
                    continue;
                }
                PairInt vXY = new PairInt(vX, vY);

                PairInt vLoc = endPointMap.get(vXY);

                if ((vLoc == null) || (vLoc.getX() == uEdgeIdx)) {
                    continue;
                }

                int vEdgeIdx = vLoc.getX();

                neighborEdgeNumbers.add(Integer.valueOf(vEdgeIdx));

                if (neighborEdgeNumbers.size() > 1) {
                    // this is a junction, so break and continue at next uXY
                    break;
                }

                PairIntArray vEdge = edges.get(vEdgeIdx);

                int canBeReordered = canBeReordered(vLoc, vEdge);

                if (canBeReordered == -1) {
                    continue;
                }

                int diffX = uX - vX;
                int diffY = uY - vY;
                int distSq = (diffX * diffX) + (diffY * diffY);

                if ((distSq < closestDistSq) || ((distSq == closestDistSq)
                    && (canBeReordered == 0))) {

                    closestDistSq = distSq;
                    vClosestLoc = vLoc;
                    vClosestXY = vXY;
                    vClosestCanBeReordered = canBeReordered;
                }
            }

            if ((neighborEdgeNumbers.size() != 1) || (vClosestXY == null)) {
                continue;
            }

            // make the smaller edge number as variables '0'
            PairInt xy0 = uXY;
            PairInt xy1 = vClosestXY;
            PairInt loc0 = uLoc;
            PairInt loc1 = vClosestLoc;
            boolean smallestEdgeIdxIsV = false;
            if (loc0.getX() > loc1.getX()) {
                PairInt swap = loc0;
                loc0 = loc1;
                loc1 = swap;
                smallestEdgeIdxIsV = true;
                swap = xy0;
                xy0 = xy1;
                xy1 = swap;
            }
            int n0 = edges.get(loc0.getX()).getN();
            int n1 = edges.get(loc1.getX()).getN();

            theJoinPoints.add(new PointPairInt(loc0, loc1));

            // add the location points to the edgeFirstEndPointMap and
            // edgeLastEndPointMap maps
            for (int ii = 0; ii < 2; ++ii) {
                int locX = loc0.getX();
                int locY = loc0.getY();
                int nEdge = n0;

                if (ii == 1) {
                    locX = loc1.getX();
                    locY = loc1.getY();
                    nEdge = n1;
                }

                boolean addToFirst = (locY < 2);
                boolean addToLast = (locY > (nEdge - 3));
                if (nEdge == 2) {
                    if (locY == 1) {
                        addToFirst = false;
                        addToLast = true;
                    } else {
                        addToFirst = true;
                        addToLast = false;
                    }
                }

                if ((ii == 0) && smallestEdgeIdxIsV && (vClosestCanBeReordered == 1) && (n0 == 3)) {
                    log.fine(String.format(
                    "creating a join point for an edge with size 3: xy=(%d,%d) loc=%d:%d",
                    vClosestLoc.getX(), vClosestLoc.getY(), vClosestLoc.getX(),
                    vClosestLoc.getY()));
                    // an entry is added to both edge maps because we don't know yet
                    // which to use.
                }

                if (addToFirst) {
                    Integer key = Integer.valueOf(locX);
                    Set<PointPairInt> joinPointSet = edgeFirstEndPointMap.get(key);
                    if (joinPointSet == null) {
                        joinPointSet = new HashSet<PointPairInt>();
                    }
                    PointPairInt joinPoint = new PointPairInt(loc0, loc1);
                    joinPointSet.add(joinPoint);
                    edgeFirstEndPointMap.put(key, joinPointSet);
                }

                if (addToLast) {
                    Integer key = Integer.valueOf(locX);
                    Set<PointPairInt> joinPointSet = edgeLastEndPointMap.get(key);
                    if (joinPointSet == null) {
                        joinPointSet = new HashSet<PointPairInt>();
                    }
                    PointPairInt joinPoint = new PointPairInt(loc0, loc1);
                    joinPointSet.add(joinPoint);
                    edgeLastEndPointMap.put(key, joinPointSet);
                }
            }
        }

        reduceMultipleEndpointsForEdge(edges, edgeFirstEndPointMap,
            edgeLastEndPointMap, endPointMap, theJoinPoints);

        if (debug) {
            // 2 of the same join points for edges of size 3 in which
            // the join point is movable to the first or last point in edge.
            // the next block removes that possibility
            boolean skipForSize3 = true;
            algorithms.misc.MiscDebug.assertConsistentJoinPointStructures(edges,
                edgeFirstEndPointMap, edgeLastEndPointMap, theJoinPoints,
                skipForSize3);
        }

        // edges that are size 3 and have join points that are not in the
        // first or last position, can be moved to either, so keep
        // track of them to avoid moving to the same position twice
        Map<PairInt, Set<Integer>> edgeSize3NonEndPoints = new HashMap<PairInt,
            Set<Integer>>();
        for (PointPairInt entry : theJoinPoints) {
            PairInt loc0 = entry.getKey();
            PairInt loc1 = entry.getValue();
            int n0 = edges.get(loc0.getX()).getN();
            int n1 = edges.get(loc1.getX()).getN();
            if ((n0 == 3) && (loc0.getY() == 1)) {
                if (!edgeSize3NonEndPoints.containsKey(loc0)) {
                    edgeSize3NonEndPoints.put(loc0, new HashSet<Integer>());
                }
            }
            if ((n1 == 3) && (loc1.getY() == 1)) {
                if (!edgeSize3NonEndPoints.containsKey(loc1)) {
                    edgeSize3NonEndPoints.put(loc1, new HashSet<Integer>());
                }
            }
        }

        // 2 * O(N_join_points)
        // re-order the endpoints which are not first or last position in edge
        for (PointPairInt entry : theJoinPoints) {
            PairInt loc0 = entry.getKey();
            PairInt loc1 = entry.getValue();

            PairIntArray edge0 = edges.get(loc0.getX());
            PairIntArray edge1 = edges.get(loc1.getX());

            PairInt origLoc0 = new PairInt(loc0.getX(), loc0.getY());

            int swapIdx = -1;
            if (edgeSize3NonEndPoints.containsKey(loc0)) {
                swapIdx = 0;
                Set<Integer> takenPositions = edgeSize3NonEndPoints.get(loc0);
                if (takenPositions.size() > 1) {
                    throw new IllegalStateException(
                        "then are more than 2 join points for the same edge of size 3." +
                        " edge " + loc0.getX() + " xy=(" + edge0.getX(loc0.getY())
                        + "," + edge0.getY(loc0.getY()) +")");
                }
                if (!takenPositions.isEmpty()) {
                    swapIdx = 2;
                }
                takenPositions.add(Integer.valueOf(swapIdx));
            } else  {
                int canBeReordered = canBeReordered(loc0, edge0);
                if (canBeReordered == 1) {
                    if (loc0.getY() < 2) {
                        swapIdx = 0;
                    } else {
                        swapIdx = edge0.getN() - 1;
                    }
                }
            }
            if (swapIdx > -1) {
                int swapX = edge0.getX(swapIdx);
                int swapY = edge0.getY(swapIdx);
                edge0.set(swapIdx, edge0.getX(1), edge0.getY(1));
                edge0.set(1, swapX, swapY);
                loc0.setY(swapIdx);
            }

            swapIdx = -1;
            if (edgeSize3NonEndPoints.containsKey(loc1)) {
                swapIdx = 0;
                Set<Integer> takenPositions = edgeSize3NonEndPoints.get(loc1);
                if (takenPositions.size() > 1) {
                    throw new IllegalStateException(
                        "then are more than 2 join points for the same edge of size 3." +
                        " edge " + loc1.getX() + " xy=(" + edge1.getX(loc1.getY())
                        + "," + edge1.getY(loc1.getY()) +")");
                }
                if (!takenPositions.isEmpty()) {
                    swapIdx = 2;
                }
                takenPositions.add(Integer.valueOf(swapIdx));
            } else  {
                int canBeReordered = canBeReordered(loc1, edge1);
                if (canBeReordered == 1) {
                    if (loc1.getY() < 2) {
                        swapIdx = 0;
                    } else {
                        swapIdx = edge1.getN() - 1;
                    }
                }
            }
            if (swapIdx > -1) {
                int swapX = edge1.getX(swapIdx);
                int swapY = edge1.getY(swapIdx);
                edge1.set(swapIdx, edge1.getX(1), edge1.getY(1));
                edge1.set(1, swapX, swapY);
                loc1.setY(swapIdx);
            }
        }

        if (debug) {
            boolean skipForSize3 = false;
            algorithms.misc.MiscDebug.assertConsistentJoinPointStructures(edges,
                edgeFirstEndPointMap, edgeLastEndPointMap, theJoinPoints,
                skipForSize3);
        }

        Map<PairInt, PairInt> result = new HashMap<PairInt, PairInt>();
        for (PointPairInt entry : theJoinPoints) {
            PairInt loc0 = entry.getKey();
            PairInt loc1 = entry.getValue();
            result.put(loc0, loc1);
        }

        return result;
    }

    /**
     * Iterate over each point looking for its neighbors and noting when
     * there are more than 2, and store the found junctions as member variables.
     *
     * runtime complexity is O(N)
     *
     * @param edges
     */
    protected void findJunctions(List<PairIntArray> edges) {

        // key = center of junction pixel coordinates
        // value = set of adjacent pixels when there are more than the preceding
        //         and next.
        Map<Integer, Set<Integer> > theJunctionMap = new HashMap<Integer, Set<Integer>>();

        // key = pixel coordinates of all pixels involved in junctions
        // value = PairInt holding index of edge that pixel is located in and
        //         holding the index within that edge of the pixel.
        //         for example, a pixel located in edges(0) at offset=100
        //         would have PairInt(0, 100) as a value.
        Map<Integer, PairInt> theJunctionLocationMap = new HashMap<Integer, PairInt>();

        int w = img.getWidth();
        int h = img.getHeight();

        findJunctions(edges, theJunctionMap, theJunctionLocationMap, w, h);

        junctionMap = theJunctionMap;
        junctionLocationMap = theJunctionLocationMap;
    }

    /**
     * Iterate over each point looking for its neighbors and noting when
     * there are more than 2, and store the found junctions as member variables.
     *
     * runtime complexity is O(N)
     *
     * @param edges
     */
    protected static void findJunctions(List<PairIntArray> edges,
        Map<Integer, Set<Integer> > theJunctionMap,
        Map<Integer, PairInt> theJunctionLocationMap,
        int imageWidth, int imageHeight) {

        int n = edges.size();

        // key = image pixel index,
        // value = pairint of edge index, and index within edge
        Map<Integer, PairInt> pointLocator = new HashMap<Integer, PairInt>();

        // key = center of junction pixel coordinates
        // value = number of times this point is a value in the first theJunctionMap
        Map<Integer, Integer> theJunctionFrequencyMap = new HashMap<Integer, Integer>();

        int w = imageWidth;
        int h = imageHeight;

        // O(N)
        for (int edgeIdx = 0; edgeIdx < n; ++edgeIdx) {

            PairIntArray edge = edges.get(edgeIdx);

            for (int uIdx = 0; uIdx < edge.getN(); ++uIdx) {

                int col = edge.getX(uIdx);
                int row = edge.getY(uIdx);

                int uPixIdx = (row * w) + col;

                pointLocator.put(Integer.valueOf(uPixIdx), new PairInt(edgeIdx,
                    uIdx));
            }
        }

        // 8 * O(N)
        for (int edgeIdx = 0; edgeIdx < n; ++edgeIdx) {

            PairIntArray edge = edges.get(edgeIdx);

            for (int uIdx = 0; uIdx < edge.getN(); ++uIdx) {

                int col = edge.getX(uIdx);
                int row = edge.getY(uIdx);

                int uPixIdx = (row * w) + col;

                Set<PairInt> neighbors = new HashSet<PairInt>();

                for (int nIdx = 0; nIdx < dxs8.length; ++nIdx) {

                    int x = col + dxs8[nIdx];
                    int y = row + dys8[nIdx];

                    if ((x < 0) || (x > (w - 1)) || (y < 0) || (y > (h - 1))) {
                        continue;
                    }

                    int vPixIdx = (y * w) + x;

                    PairInt vLoc = pointLocator.get(Integer.valueOf(vPixIdx));

                    if (vLoc != null) {
                        neighbors.add(vLoc);
                    }
                }

                if (neighbors.size() > 2) {

                    Set<Integer> indexes = new HashSet<Integer>();

                    for (PairInt p : neighbors) {

                        int edge2Idx = p.getX();
                        int vIdx = p.getY();

                        PairIntArray vEdge = edges.get(edge2Idx);

                        int vPixIdx = (vEdge.getY(vIdx) * w) + vEdge.getX(vIdx);

                        Integer key = Integer.valueOf(vPixIdx);

                        theJunctionLocationMap.put(key, p);

                        indexes.add(key);

                        Integer freq = theJunctionFrequencyMap.get(key);
                        if (freq == null) {
                            theJunctionFrequencyMap.put(key, Integer.valueOf(1));
                        } else {
                            theJunctionFrequencyMap.put(key,
                                Integer.valueOf(freq.intValue() + 1));
                        }
                    }

                    theJunctionMap.put(Integer.valueOf(uPixIdx), indexes);

                    theJunctionLocationMap.put(Integer.valueOf(uPixIdx),
                        new PairInt(edgeIdx, uIdx));
                }

                // if (neighbors.size() == 2) is handled in findJoinPoints
            }
        }

        // visit each junction and make sure the real center is the only one
        // stored
        Set<Integer> remove = new HashSet<Integer>();
        Set<Integer> doNotRemove = new HashSet<Integer>();
        for (Entry<Integer, Set<Integer>> entry : theJunctionMap.entrySet()) {

            Integer pixelIndex = entry.getKey();

            Set<Integer> neighborIndexes = entry.getValue();
            int nN = neighborIndexes.size();

            int maxN = Integer.MIN_VALUE;

            for (Integer pixelIndex2 : neighborIndexes) {

                if (remove.contains(pixelIndex2)) {
                    continue;
                }

                Set<Integer> neighborIndexes2 = theJunctionMap.get(pixelIndex2);

                if (neighborIndexes2 == null) {
                    continue;
                }

                if (neighborIndexes2.size() > nN) {
                    if (neighborIndexes2.size() > maxN) {
                        maxN = neighborIndexes2.size();
                    }
                }
            }

            if (maxN > Integer.MIN_VALUE) {

                if (!doNotRemove.contains(pixelIndex)) {
                    remove.add(pixelIndex);
                }

            } else {
                doNotRemove.add(pixelIndex);
                // remove the neighbors from the junction map
                for (Integer pixelIndex2 : neighborIndexes) {
                    int mapFreq2 = theJunctionFrequencyMap.get(pixelIndex2);

                    if (!doNotRemove.contains(pixelIndex2) && (mapFreq2 == 1)) {
                        remove.add(pixelIndex2);
                    }
                }
            }
        }

        for (Integer pixelIndex : remove) {
            theJunctionMap.remove(pixelIndex);
        }

    }

    private Map<Integer, PairIntArray> createIndexedMap(List<PairIntArray> edges) {

        Map<Integer, PairIntArray> output = new HashMap<Integer, PairIntArray>();
        for (int i = 0; i < edges.size(); i++) {
            output.put(Integer.valueOf(i), edges.get(i));
        }

        return output;
    }

    /**
     * join edges using information in the member variable joinPoints
     * and update the junction and joinPoint information after the
     * changes.  Note that the join points must be first or last positions
     * in an edge in edges.
     * @param joinPoints map with key = the PairInt with x being the edge index
     * and y being the index of the first join point within the edge, value =
     * the PairInt with x being the edge index
     * and y being the index of the other join point within the edge
     * @param edges
     * @return
     */
    protected List<PairIntArray> joinOnJoinPoints(Map<PairInt, PairInt>
        joinPoints, List<PairIntArray> edges) {

        //order the join points by the edge number of the first join point

        int[] indexes = new int[joinPoints.size()];
        PairInt[][] edgeJoins = new PairInt[joinPoints.size()][2];
        int count = 0;
        for (Entry<PairInt, PairInt> entry : joinPoints.entrySet()) {
            PairInt loc0 = entry.getKey();
            PairInt loc1 = entry.getValue();
            indexes[count] = loc0.getX();
            edgeJoins[count] = new PairInt[]{loc0, loc1};
            count++;
        }
        QuickSort.sortBy1stArg(indexes, edgeJoins);

long nPointsBefore = countPixelsInEdges(edges);

        // put the edges in a map to remove and search updates faster
        Map<Integer, PairIntArray> edgesMap = createIndexedMap(edges);

        int n = edgeJoins.length;

        for (int i = (n - 1); i > -1; --i) {

//algorithms.misc.MiscDebug.printJoinPoints(edgeJoins, 0, i, edgesMap);

            PairInt[] entry = edgeJoins[i];

            PairInt loc0 = entry[0];

            PairInt loc1 = entry[1];

            // if they've already been merged
            if (loc0.getX() == loc1.getX()) {
                continue;
            }

            // the smaller edge index should be loc1 and the larger loc0,
            if (loc1.getX() > loc0.getX()) {
                PairInt tmp = loc0;
                loc0 = loc1;
                loc1 = tmp;
            }
            final int edge0Idx = loc0.getX();
            final int edge1Idx = loc1.getX();

            // edge to move
            int removedEdgeIdx = edge0Idx;
            PairIntArray edge0 = edgesMap.remove(Integer.valueOf(removedEdgeIdx));

            int n0 = edge0.getN();

            if ((loc0.getY() != 0) && (loc0.getY() != (n0 - 1))){
                /*
                this can happen if there is an error in the join point algorithm
                resulting in the same join point to 2 different edges.  an update
                in the location will eventually be an error in one of them.

                After more testing, will change this to discard the join point
                and warn of error.
                */
                throw new IllegalStateException("ERROR in the updates? " +
                    " loc0=" + edge0Idx + "," + loc0.getY() + " n=" + n0 +
                    " i=" + i + " (nedges=" + n + ") to append edge "
                    + edge0Idx + " to edge " + edge1Idx);
            }

            // join point should be at the beginning, so reverse if not
            if (loc0.getY() == (n0 - 1)) {

                edge0.reverse();
                loc0.setY(0);

                // everything with smaller index than i in edgeJoins with
                // edgeIndex==loc0.getX() needs to be updated for this reversal.
                // idx becomes n-idx-1
                for (int j = (i - 1); j > -1; --j) {
                    PairInt[] vEntry = edgeJoins[j];
                    for (int k = 0; k < 2; k++) {
                        PairInt vLoc = vEntry[k];
                        if (vLoc.getX() == edge0Idx) {
                            int idxRev = n0 - vLoc.getY() - 1;
                            vLoc.setY(idxRev);
                        }
                    }
                }
            }
            // edge to receive new edge
            PairIntArray edge1 = edgesMap.get(Integer.valueOf(loc1.getX()));
            int n1 = edge1.getN();
            // endpoint should be at the end, so reverse if not
            if (loc1.getY() == 0) {

                edge1.reverse();
                loc1.setY(n1 - 1);

                // everything with smaller index than i in edgeJoins that has
                // edgeIndex==loc1.getX() needs to be updated for this reversal.
                // idx becomes n-idx-1
                for (int j = (i - 1); j > -1; --j) {
                    PairInt[] vEntry = edgeJoins[j];
                    for (int k = 0; k < 2; k++) {
                        PairInt vLoc = vEntry[k];
                        if (vLoc.getX() == edge1Idx) {
                            int idxRev = n1 - vLoc.getY() - 1;
                            vLoc.setY(idxRev);
                        }
                    }
                }
            }

            final int indexWithinEdge0 = loc0.getY();
            final int indexWithinEdge1 = loc1.getY();

            // --- append edge0 to edge1 ----
            edge1.addAll(edge0);

            // for earlier items in array edgeJoins
            // need to update all edge indexes and indexes within edge.

            // loc0 got appended to loc1 --> [edge1][edge0]
            // points in edge0 need n1 added to them

            // first, will only update the indexes within the edge for
            // edgeIndex == loc0.getX()
            //     the new offset index = n0 + index
            for (int j = (i - 1); j > -1; --j) {
                PairInt[] vEntry = edgeJoins[j];
                for (int k = 0; k < 2; k++) {
                    PairInt vLoc = vEntry[k];
                    if (vLoc.getX() == edge0Idx) {
                        int idxEdit = vLoc.getY() + n1;
                        vLoc.setY(idxEdit);
                    }
                }
            }

            // any edge index > loc0 gets reduced by one, but if == it gets replaced by loc1.getX()
            for (int j = (i - 1); j > -1; --j) {
                PairInt[] vEntry = edgeJoins[j];
                for (int k = 0; k < 2; k++) {
                    PairInt vLoc = vEntry[k];
                    if (vLoc.getX() > edge0Idx) {
                        int editX = vLoc.getX() - 1;
                        vLoc.setX(editX);
                    } else if (vLoc.getX() == edge0Idx) {
                        int editX = edge1Idx;
                        vLoc.setX(editX);
                    }
                }
            }

            // the output map keeps 0 to loc1.getx(),
            // but loc0.getX() + 1 gets moved to loc0.getX()
            //    and on until have reached size() - 1
            for (int j = (removedEdgeIdx + 1); j <= edgesMap.size(); ++j) {
                PairIntArray v = edgesMap.remove(Integer.valueOf(j));
                assert(v != null);
                edgesMap.put(Integer.valueOf(j - 1), v);
            }

        }

        List<PairIntArray> output = new ArrayList<PairIntArray>();
        for (int i = 0; i < edgesMap.size(); i++) {
            Integer key = Integer.valueOf(i);
            PairIntArray v = edgesMap.get(key);
            assert(v != null);
            output.add(v);
        }

        return output;
    }

    /**
     * return whether the point can be swapped with the nearest endpoint.
     * <pre>
     * returns
     *     0: for no need to reorder
     *     1: for it not being an endpoint but can be re-ordered
     *     -1: for cannot be re-ordered
     * </pre>
     * @param pointEdgeLocation
     * @param edge
     * @return code
       <pre> returns
           0: for no need to reorder
           1: for it not being an endpoint but can be re-ordered
           -1: for cannot be re-ordered,
       </pre>
     */
    private int canBeReordered(PairInt pointEdgeLocation, PairIntArray edge) {

        int n = edge.getN();

        if ((pointEdgeLocation.getY() == 0) || (pointEdgeLocation.getY() == (n - 1))) {
            return 0;
        }

        if ((pointEdgeLocation.getY() > 2) && (pointEdgeLocation.getY() < (n - 3))) {
            return -1;
        }

        int idxOrig = pointEdgeLocation.getY();

        int check0Idx = 0;
        int check1Idx = 2;

        if ((n > 3) && (idxOrig == (n - 2))) {
            check0Idx = n - 3;
            check1Idx = n - 1;
        }

        int diffX = edge.getX(check0Idx) - edge.getX(check1Idx);
        int diffY = edge.getY(check0Idx) - edge.getY(check1Idx);
        if ((Math.abs(diffX) > 1) || (Math.abs(diffY) > 1)) {
            return -1;
        }

        return 1;
    }

    /**
     * If pointEdgeLocation is near the beginning or end of edge,
     * swap it with the point that is the endpoint and update the given
     * data structures with the updated information.
     * @param pointEdgeLocation
     * @param edge
     * @return 1 for did re-order, 0 for no need to re-order, else -1 for cannot
     * re-order
     */
    private int reorderIfNearEnd(PairInt pointEdgeLocation,
        PairIntArray edge,
        PairInt connectingEdgeLocation, PairIntArray connectingEdge) {

        int n = edge.getN();

        if ((pointEdgeLocation.getY() == 0) || (pointEdgeLocation.getY() == (n - 1))) {
            return 0;
        }

        /*
        looks like max offset from end is 2
         @ . [....]

         @ .
         . [....]
        */

        if ((pointEdgeLocation.getY() > 2) && (pointEdgeLocation.getY() < (n - 3))) {
            return 0;
        }

        int idxOrig = pointEdgeLocation.getY();
        int idxSwap;
        if (idxOrig == 1) {
            idxSwap = 0;
        } else if (idxOrig == (n - 2)) {
            idxSwap = n - 1;
        } else {
            /*
            can just swap it with first or last point if the swapped point
            is adjacent to the point right before the point to be swapped,
            in other words, does not break a connection in the edge
            */
            if (idxOrig == 2) {
                idxSwap = 0;
            } else {
                idxSwap = n - 1;
            }
            int prevX = edge.getX(idxOrig - 1);
            int prevY = edge.getY(idxOrig - 1);
            int endX = edge.getX(0);
            int endY = edge.getY(0);
            if ((Math.abs(prevX - endX) > 1) || (Math.abs(prevY - endY) > 1)) {
                return -1;
            }
        }

        /*
        if swap position is still adjacent to the connecting point,
        can complete the change.
        */
        int connectedX = connectingEdge.getX(connectingEdgeLocation.getY());
        int connectedY = connectingEdge.getY(connectingEdgeLocation.getY());

        int swapX = edge.getX(idxSwap);
        int swapY = edge.getY(idxSwap);

        if ((Math.abs(connectedX - swapX) > 1) || (Math.abs(connectedY - swapY) > 1)) {
            return -1;
        }

        pointEdgeLocation.setY(idxSwap);

        return 1;
    }

    private int spliceEdgesAtJunctionsIfImproves(List<PairIntArray> edges) {

        /*
        The main goal is to make better contours.

        Edges with junctions can sometimes be spliced and rejoined with another
        junction edge to make a better edge where better edge may be a longer
        edge or a closed contour useful for determining transformations
        between images containing the contour.

        could be used with shape templates...

        For now, will make a simple algorithm which tries to increase the
        length of the longest edges in a junction.
        */

        int nSplices = 0;

        // key = edge index.
        // value = pixel indexes.
        //   the pixel indexes are used to find values in junctionLocatorMap
        //   to update it as points are moved to and from edges.
        // key = edge index.
        // value = pixel indexes.
        //   the pixel indexes are used to find values in junctionLocatorMap
        //   to update it as points are moved to and from edges.
        Map<Integer, Set<Integer>> theEdgeToPixelIndexMap = createEdgeToPixelIndexMap();

        if (debug) {
            algorithms.misc.MiscDebug.assertConsistentEdgeCapacity(
                theEdgeToPixelIndexMap, junctionLocationMap, junctionMap,
                edges);
        }

        /*
        TODO: consider handling splices of edges of size 1 after this
        block.  They can be inserted in a way that does not affect line
        continuity

        For example:

            0
             1 < where '<' can be inserted after 1 and before 2 or vice versa
             2
             3
        */

        for (Entry<Integer, Set<Integer>> entry : junctionMap.entrySet()) {

            Integer centerPixelIndex = entry.getKey();
            PairInt centerLoc = junctionLocationMap.get(centerPixelIndex);
            assert(centerLoc != null);

            /*
            if (debug) {

                log.info("processing junction w/ center pixel index=" + centerPixelIndex
                    + " and loc=" + centerLoc.getX() + ":" + centerLoc.getY());

                String str = algorithms.misc.MiscDebug.printJunctionsToString(
                    junctionLocationMap, junctionMap, edges, img);
                log.info(str);

            }
            */

            Set<Integer> adjIndexes = entry.getValue();

            int[] pixIndexes = new int[adjIndexes.size() + 1];

            // lengths holds the edge up until the junction.  splices the edge
            // at the junction figuratively and counts number of pixels before
            // and after splice and keeps the longest.
            int[] lengths = new int[pixIndexes.length];

            pixIndexes[0] = centerPixelIndex.intValue();
            lengths[0] = (new Splice(edges.get(centerLoc.getX())))
                .getLengthOfLongestSide(centerLoc.getY());

            int maxN = lengths[0];
            int maxNPixIdx = pixIndexes[0];
            boolean foundAnotherEdge = false;

            int count = 1;

            // for junctions on the same edge, need to count the portions
            // of their edge which are on oppossite sides of the junction,
            // so need to track the maxN, then correct the lengths of splice
            // from same edge afterwards

            for (Integer pixIndex : adjIndexes) {

                pixIndexes[count] = pixIndex.intValue();

                PairInt loc = junctionLocationMap.get(pixIndex);

                if ((centerLoc.getX() == loc.getX())) {
                    count++;
                    continue;
                }

                foundAnotherEdge = true;

                Splice splice = new Splice(edges.get(loc.getX()));

                lengths[count] = splice.getLengthOfLongestSide(loc.getY());

                if (lengths[count] > maxN) {
                    maxN = lengths[count];
                    maxNPixIdx = pixIndex;
                }

                count++;
            }

            if (!foundAnotherEdge) {
                continue;
            }

            //-- correct the lengths for the junctions from edge maxNEdgeIdx --
            PairInt maxNLoc = junctionLocationMap.get(maxNPixIdx);
            Splice maxNSplice = new Splice(edges.get(maxNLoc.getX()));
            int maxNSpliceOtherSideN = maxNSplice.splice(new int[]{maxNLoc.getY()})[1].getN();
            // correct the splice lengths if they are from the same edge
            for (int ii = 0; ii < count; ++ii) {
                int pixIdx = pixIndexes[ii];
                PairInt loc = junctionLocationMap.get(pixIdx);
                if ((loc.getX() != maxNLoc.getX()) || (pixIdx == maxNPixIdx)) {
                    continue;
                }
                lengths[ii] = maxNSpliceOtherSideN;
            }

            if ((maxN > lengths.length) || (maxN > 10000000)) {
                MultiArrayMergeSort.sortByDecr(lengths, pixIndexes);
            } else {
                CountingSort.sortByDecr(lengths, pixIndexes, maxN);
            }

            int pixIdx0 = pixIndexes[0];

            int pixIdx1 = pixIndexes[1];

            PairInt loc0 = junctionLocationMap.get(Integer.valueOf(pixIdx0));
            final int edge0Idx = loc0.getX();
            final int indexWithinEdge0 = loc0.getY();

            PairInt loc1 = junctionLocationMap.get(Integer.valueOf(pixIdx1));
            final int edge1Idx = loc1.getX();
            final int indexWithinEdge1 = loc1.getY();

            if (edge0Idx != edge1Idx && (lengths[0] != 0) && (lengths[1] != 0)) {

                Set<Integer> edgePixelIndexes =
                    theEdgeToPixelIndexMap.get(Integer.valueOf(edge0Idx));
                assert(edgePixelIndexes != null);

                int[] splice0Y = new int[]{indexWithinEdge0};
                Splice splice0 = new Splice(edges.get(edge0Idx));

                Map<Integer, PairInt> smallerSpliceLocations0 = new
                    HashMap<Integer, PairInt>();
                Map<Integer, PairInt> largerSpliceLocations0 = new
                    HashMap<Integer, PairInt>();
                // splice splice0 into 2 edges and put updated i
                PairIntArray[] spliced0 = splice0.splice(splice0Y,
                    edgePixelIndexes, junctionLocationMap,
                    smallerSpliceLocations0, largerSpliceLocations0);

                int[] splice1Y = new int[]{indexWithinEdge1};
                Splice splice1 = new Splice(edges.get(edge1Idx));

                edgePixelIndexes =
                    theEdgeToPixelIndexMap.get(Integer.valueOf(edge1Idx));
                assert(edgePixelIndexes != null);
                Map<Integer, PairInt> smallerSpliceLocations1 = new
                    HashMap<Integer, PairInt>();
                Map<Integer, PairInt> largerSpliceLocations1 = new
                    HashMap<Integer, PairInt>();
                PairIntArray[] spliced1 = splice1.splice(splice1Y,
                    edgePixelIndexes, junctionLocationMap,
                    smallerSpliceLocations1, largerSpliceLocations1);

                // if not spliced, continue.
                if ((spliced0[0].getN() == 0) || (spliced1[0].getN() == 0)) {
                    continue;
                }

                int splice01InsertedIdx = edges.size();

                // add the smaller part of spliced0 and spliced1 to edges
                edges.add(spliced0[1]);

                int splice11InsertedIdx = edges.size();

                edges.add(spliced1[1]);

                int splice00N = spliced0[0].getN();
                int splice10N = spliced1[0].getN();

                // if splice0Y[0] is first point, reverse the edge before append
                if (splice0Y[0] == 0) {
                    spliced0[0].reverse();
                    splice0Y[0] = splice00N - 1;

                    // --- update largerSpliceLocations0 to reverse index within edge
                    for (Entry<Integer, PairInt> sEntry : largerSpliceLocations0.entrySet()) {
                        PairInt sLoc = sEntry.getValue();
                        int idxRev = splice00N - sLoc.getY() - 1;

                        sLoc.setY(idxRev);
                        assert(idxRev < splice00N);
                    }
                }

                // if splice1Y[0] is not the first point, reverse the edge before append
                if (splice1Y[0] != 0) {
                    spliced1[0].reverse();
                    splice1Y[0] = 0;

                    // --- update largerSpliceLocations1 to reverse index within edge
                    for (Entry<Integer, PairInt> sEntry : largerSpliceLocations1.entrySet()) {
                        PairInt sLoc = sEntry.getValue();
                        int idxRev = splice10N - sLoc.getY() - 1;

                        sLoc.setY(idxRev);
                        assert(idxRev < splice10N);
                    }
                }

                // append splice1 to splice0
                spliced0[0].addAll(spliced1[0]);

                ++nSplices;

                PairIntArray edge0 = edges.get(edge0Idx);
                edge0.swapContents(spliced0[0]);
                spliced0[0] = null;

                // --- update location map for information in
                // --- smallerSpliceLocations0, largerSpliceLocations0 and
                // --- smallerSpliceLocations1, largerSpliceLocations1
                for (Entry<Integer, PairInt> sEntry : largerSpliceLocations0.entrySet()) {
                    Integer sPixelIndex = sEntry.getKey();
                    PairInt sLoc = sEntry.getValue();
                    PairInt loc = junctionLocationMap.get(sPixelIndex);
                    assert (loc != null);

                    // edge index remains same, but the index within edge may have changed
                    loc.setY(sLoc.getY());
                    assert(sLoc.getY() < edge0.getN());
                }
                for (Entry<Integer, PairInt> sEntry : smallerSpliceLocations0.entrySet()) {
                    Integer sPixelIndex = sEntry.getKey();
                    PairInt sLoc = sEntry.getValue();
                    PairInt loc = junctionLocationMap.get(sPixelIndex);
                    assert (loc != null);

                    // location is new edge with edited index within edge
                    loc.setX(splice01InsertedIdx);
                    loc.setY(sLoc.getY());

                    //-- move item from theEdgeToPixelIndexMap too
                    //-- remove sPixelIndex from edge loc0.getX() and add it to splice01InsertedIdx
                    Set<Integer> sPixelIndexes = theEdgeToPixelIndexMap.get(Integer.valueOf(edge0Idx));
                    boolean removed = sPixelIndexes.remove(sPixelIndex);
                    assert(removed == true);

                    sPixelIndexes = theEdgeToPixelIndexMap.get(Integer.valueOf(splice01InsertedIdx));
                    if (sPixelIndexes == null) {
                        sPixelIndexes = new HashSet<Integer>();
                    }
                    sPixelIndexes.add(sPixelIndex);
                    theEdgeToPixelIndexMap.put(Integer.valueOf(splice01InsertedIdx), sPixelIndexes);
                }
                for (Entry<Integer, PairInt> sEntry : smallerSpliceLocations1.entrySet()) {
                    Integer sPixelIndex = sEntry.getKey();
                    PairInt sLoc = sEntry.getValue();
                    PairInt loc = junctionLocationMap.get(sPixelIndex);
                    assert (loc != null);

                    // location is new edge with edited index within edge
                    loc.setX(splice11InsertedIdx);
                    loc.setY(sLoc.getY());

                    //-- move item from theEdgeToPixelIndexMap too
                    //-- remove sPixelIndex from edge loc1.getX() and add it to splice11InsertedIdx
                    Set<Integer> sPixelIndexes =
                        theEdgeToPixelIndexMap.get(Integer.valueOf(edge1Idx));
                    boolean removed = sPixelIndexes.remove(sPixelIndex);
                    assert(removed == true);
                    theEdgeToPixelIndexMap.put(Integer.valueOf(edge1Idx), sPixelIndexes);

                    sPixelIndexes = theEdgeToPixelIndexMap.get(Integer.valueOf(splice11InsertedIdx));
                    if (sPixelIndexes == null) {
                        sPixelIndexes = new HashSet<Integer>();
                    }
                    sPixelIndexes.add(sPixelIndex);
                    theEdgeToPixelIndexMap.put(Integer.valueOf(splice11InsertedIdx), sPixelIndexes);
                }
                for (Entry<Integer, PairInt> sEntry : largerSpliceLocations1.entrySet()) {
                    Integer sPixelIndex = sEntry.getKey();
                    PairInt sLoc = sEntry.getValue();
                    PairInt loc = junctionLocationMap.get(sPixelIndex);
                    assert (loc != null);

                    // location after append is loc0.getX() with offset of splice00N
                    loc.setX(edge0Idx);
                    loc.setY(sLoc.getY() + splice00N);
                    assert(loc.getY() < edge0.getN());

                    //-- move item from theEdgeToPixelIndexMap too
                    //-- remove sPixelIndex from edge loc1.getX() and add it to loc0.getX()
                    Integer key = Integer.valueOf(edge1Idx);
                    Set<Integer> sPixelIndexes = theEdgeToPixelIndexMap.get(key);
                    boolean removed = sPixelIndexes.remove(sPixelIndex);
                    assert(removed == true);
                    theEdgeToPixelIndexMap.put(key, sPixelIndexes);

                    key = Integer.valueOf(edge0Idx);
                    sPixelIndexes = theEdgeToPixelIndexMap.get(key);
                    if (sPixelIndexes == null) {
                        sPixelIndexes = new HashSet<Integer>();
                    }
                    sPixelIndexes.add(sPixelIndex);
                    theEdgeToPixelIndexMap.put(key, sPixelIndexes);
                }

                // remove edge1Idx from edges as it's now part of edge edge0Idx
                PairIntArray removedEdge = edges.remove(edge1Idx);
                assert(removedEdge != null);

                if (debug) {
                    // assert that theEdgeToPixelIndexMap does not have any
                    // entries for the edge that was just removed.
                    // all entries in theEdgeToPixelIndexMap and the junction
                    // maps should have been updated above for edge1Idx
                    Set<Integer> pixelIndexes = theEdgeToPixelIndexMap.get(edge1Idx);
                    assert(pixelIndexes == null || pixelIndexes.isEmpty());
                }

                // ---- update all entries in junctionLocationMap and pixelIndexes for the removal
                // ---- there shouldn't be any entries in for loc1.getX() at this point though
                /*
                0                0
                ...             ...
                49               49
                50 <-- rm        prev 51, now 50.
                51
                */
                for (int eIdx = (edge1Idx + 1); eIdx <= edges.size(); ++eIdx) {

                    Integer edgeIndex = Integer.valueOf(eIdx);

                    Set<Integer> pixelIndexes = theEdgeToPixelIndexMap.get(edgeIndex);

                    if (pixelIndexes != null) {
                        for (Integer pixelIndex : pixelIndexes) {
                            PairInt loc = junctionLocationMap.get(pixelIndex);
                            loc.setX(loc.getX() - 1);
                        }

                        theEdgeToPixelIndexMap.remove(edgeIndex);

                        Integer edgeIndexEarlier = Integer.valueOf(eIdx - 1);

                        theEdgeToPixelIndexMap.put(edgeIndexEarlier, pixelIndexes);
                    }
                }

                if (debug) {
                    algorithms.misc.MiscDebug.assertConsistentEdgeCapacity(
                        theEdgeToPixelIndexMap, junctionLocationMap,
                        junctionMap, edges);
                }
            }
        }

        return nSplices;
    }

    protected void reduceMultipleEndpointsForEdge(
        List<PairIntArray> edges,
        Map<Integer, Set<PointPairInt>> edgeFirstEndPointMap,
        Map<Integer, Set<PointPairInt>> edgeLastEndPointMap,
        Map<PairInt, PairInt> endPointMap, Set<PointPairInt> theJoinPoints) {

        Map<Integer, Set<PointPairInt>> tmpFirstRemoveMap = new HashMap<Integer, Set<PointPairInt>>();

        Map<Integer, Set<PointPairInt>> tmpLastRemoveMap = new HashMap<Integer, Set<PointPairInt>>();

        for (int type = 0; type < 2; ++type) {

            Map<Integer, Set<PointPairInt>> edgeMap = edgeFirstEndPointMap;
            if (type == 1) {
                edgeMap = edgeLastEndPointMap;
            }

            for (Entry<Integer, Set<PointPairInt>> entry : edgeMap.entrySet()) {

                Set<PointPairInt> joinPointSet = entry.getValue();

                if (joinPointSet.size() < 2) {
                    continue;
                }

                // has more than one joint point for this endpoint, though some
                // might have already been removed and are in tmpFirstRemoveMap
                // if so

                if (log.isLoggable(Level.FINE)) {
                    // construct a warning
                    StringBuilder sb = new StringBuilder("edge ");
                    sb.append(entry.getKey().toString())
                        .append(" has ").append(Integer.toString(joinPointSet.size()))
                        .append(" first endpoints, so deciding between them:");

                    for (PointPairInt joinPoint : entry.getValue()) {
                        PairInt loc0 = joinPoint.getKey();
                        PairInt loc1 = joinPoint.getValue();
                        PairIntArray edge0 = edges.get(loc0.getX());
                        PairIntArray edge1 = edges.get(loc1.getX());
                        int x0 = edge0.getX(loc0.getY());
                        int y0 = edge0.getY(loc0.getY());
                        int x1 = edge1.getX(loc1.getY());
                        int y1 = edge1.getY(loc1.getY());
                        sb.append("\n")
                            .append(" ").append(Integer.toString(loc0.getX())).append(":")
                            .append(Integer.toString(loc0.getY()))
                            .append(" (").append(Integer.toString(x0)).append(",")
                            .append(Integer.toString(y0)).append(")<-->")
                            .append(" ").append(Integer.toString(loc1.getX())).append(":")
                            .append(Integer.toString(loc1.getY()))
                            .append(" (").append(Integer.toString(x1)).append(",")
                            .append(Integer.toString(y1)).append(")");
                    }
                    log.warning(sb.toString());
                }

                // choose closest join point pair, and break ties with those that
                // have higher number of members that are already endpoints.

                PointPairInt closestJoinPoint = null;
                int closestDistSq = Integer.MAX_VALUE;
                int closestCanBeReordered0 = -99;
                int closestCanBeReordered1 = -99;

                for (PointPairInt joinPoint : joinPointSet) {

                    PairInt loc0 = joinPoint.getKey();
                    PairIntArray edge0 = edges.get(loc0.getX());
                    int n0 = edge0.getN();

                    Set<PointPairInt> removedFirstSet =
                        tmpFirstRemoveMap.get(Integer.valueOf(loc0.getX()));
                    if (removedFirstSet == null) {
                        removedFirstSet = new HashSet<PointPairInt>();
                        tmpFirstRemoveMap.put(Integer.valueOf(loc0.getX()), removedFirstSet);
                    }

                    Set<PointPairInt> removedLastSet =
                        tmpLastRemoveMap.get(Integer.valueOf(loc0.getX()));
                    if (removedLastSet == null) {
                        removedLastSet = new HashSet<PointPairInt>();
                        tmpLastRemoveMap.put(Integer.valueOf(loc0.getX()), removedLastSet);
                    }

                    if ((n0 != 3) && (loc0.getY() < 2) && removedFirstSet.contains(joinPoint)) {
                        continue;
                    }
                    if ((n0 != 3) && (loc0.getY() > (n0 - 3)) && removedLastSet.contains(joinPoint)) {
                        continue;
                    }

                    PairInt loc1 = joinPoint.getValue();
                    PairIntArray edge1 = edges.get(loc1.getX());

                    assert(edge0 != null);
                    assert(edge1 != null);

                    int diffX = edge0.getX(loc0.getY()) - edge1.getX(loc1.getY());
                    int diffY = edge0.getY(loc0.getY()) - edge1.getY(loc1.getY());
                    int distSq = (diffX * diffX) + (diffY * diffY);
                    if (distSq < closestDistSq) {
                        closestDistSq = distSq;
                        closestJoinPoint = joinPoint;
                        closestCanBeReordered0 = canBeReordered(loc0, edge0);
                        closestCanBeReordered1 = canBeReordered(loc1, edge1);
                    } else if (distSq == closestDistSq) {
                        int cr0 = canBeReordered(loc0, edge0);
                        int cr1 = canBeReordered(loc1, edge1);;
                        if (closestCanBeReordered0 == 0 && closestCanBeReordered1 == 0) {
                            continue;
                        }
                        if (((cr0 == 0) && (cr1 == 0)) ||
                            (closestCanBeReordered0 != 0 && closestCanBeReordered1 != 0
                                && ((cr0 == 0) || cr1 == 0))
                            ) {
                            closestDistSq = distSq;
                            closestJoinPoint = joinPoint;
                            closestCanBeReordered0 = cr0;
                            closestCanBeReordered1 = cr1;
                        }
                    }
                }

                //TODO: revisit this
                if (closestJoinPoint == null) {
                    continue;
                }
                //assert (closestJoinPoint != null);

                for (PointPairInt joinPoint : joinPointSet) {

                    if (joinPoint.equals(closestJoinPoint)) {
                        continue;
                    }

                    PairInt loc0 = joinPoint.getKey();
                    PairIntArray edge0 = edges.get(loc0.getX());
                    int n0 = edge0.getN();

                    Set<PointPairInt> removedFirstSet0 =
                        tmpFirstRemoveMap.get(Integer.valueOf(loc0.getX()));

                    Set<PointPairInt> removedLastSet0 =
                        tmpLastRemoveMap.get(Integer.valueOf(loc0.getX()));

                    if ((n0 != 3) && (loc0.getY() < 2) && removedFirstSet0.contains(joinPoint)) {
                        continue;
                    }
                    if ((n0 != 3) && (loc0.getY() > (n0 - 3)) && removedLastSet0.contains(joinPoint)) {
                        continue;
                    }

                    PairInt loc1 = joinPoint.getValue();
                    PairIntArray edge1 = edges.get(loc1.getX());
                    int n1 = edge1.getN();

                    assert(edge0 != null);
                    assert(edge1 != null);

                    if ((n0 == 3) && (loc0.getY() == 1)) {
                        log.fine(String.format(
                            "removing point (%d,%d) but it is in the middle of an edge of size 3",
                            edge0.getX(1), edge0.getY(1)));
                    }

                    Set<PointPairInt> removedFirstSet1 =
                        tmpFirstRemoveMap.get(Integer.valueOf(loc1.getX()));
                    if (removedFirstSet1 == null) {
                        removedFirstSet1 = new HashSet<PointPairInt>();
                        tmpFirstRemoveMap.put(Integer.valueOf(loc1.getX()), removedFirstSet1);
                    }

                    Set<PointPairInt> removedLastSet1 =
                        tmpLastRemoveMap.get(Integer.valueOf(loc1.getX()));
                    if (removedLastSet1 == null) {
                        removedLastSet1 = new HashSet<PointPairInt>();
                        tmpLastRemoveMap.put(Integer.valueOf(loc1.getX()), removedLastSet1);
                    }

                    boolean isFirstEndPoint0 = false;
                    if (loc0.getY() == 0) {
                        isFirstEndPoint0 = true;
                    } else if ((loc0.getY() != (n0 - 1)) && (loc0.getY() < 2)) {
                        isFirstEndPoint0 = true;
                    }
                    boolean isFirstEndPoint1 = false;
                    if (loc1.getY() == 0) {
                        isFirstEndPoint1 = true;
                    } else if ((loc1.getY() != (n1 - 1)) && (loc1.getY() < 2)) {
                        isFirstEndPoint1 = true;
                    }

                    if (isFirstEndPoint0 && removedFirstSet0.contains(joinPoint)) {
                        continue;
                    } else if (!isFirstEndPoint0 && removedLastSet0.contains(joinPoint)) {
                        continue;
                    }

                    if (isFirstEndPoint0) {
                        removedFirstSet0.add(joinPoint);
                    } else {
                        removedLastSet0.add(joinPoint);
                    }
                    if (isFirstEndPoint1) {
                        removedFirstSet1.add(joinPoint);
                    } else {
                        removedLastSet1.add(joinPoint);
                    }

                    boolean removed = theJoinPoints.remove(joinPoint);
                    assert(removed);

                    tmpFirstRemoveMap.put(Integer.valueOf(loc0.getX()), removedFirstSet0);
                    tmpLastRemoveMap.put(Integer.valueOf(loc0.getX()), removedLastSet0);

                    tmpFirstRemoveMap.put(Integer.valueOf(loc1.getX()), removedFirstSet1);
                    tmpLastRemoveMap.put(Integer.valueOf(loc1.getX()), removedLastSet1);
                }
            }

            // --- update the edge maps for the changes ---

            Map<Integer, Set<PointPairInt>> tmpMap = new HashMap<Integer, Set<PointPairInt>>();

            for (Entry<Integer, Set<PointPairInt>> entry : edgeFirstEndPointMap.entrySet()) {

                Set<PointPairInt> rmJoinPointSet = tmpFirstRemoveMap.get(entry.getKey());

                if (rmJoinPointSet == null || rmJoinPointSet.isEmpty()) {
                    tmpMap.put(entry.getKey(), entry.getValue());
                    continue;
                }

                Set<PointPairInt> tmpJoinPointSet = new HashSet<PointPairInt>();
                for (PointPairInt joinPoint : entry.getValue()) {
                    if (!rmJoinPointSet.contains(joinPoint)) {
                        tmpJoinPointSet.add(joinPoint);
                    }
                }
                tmpMap.put(entry.getKey(), tmpJoinPointSet);
            }
            edgeFirstEndPointMap.clear();
            edgeFirstEndPointMap.putAll(tmpMap);

            tmpMap = new HashMap<Integer, Set<PointPairInt>>();

            for (Entry<Integer, Set<PointPairInt>> entry : edgeLastEndPointMap.entrySet()) {

                Set<PointPairInt> rmJoinPointSet = tmpLastRemoveMap.get(entry.getKey());

                if (rmJoinPointSet == null || rmJoinPointSet.isEmpty()) {
                    tmpMap.put(entry.getKey(), entry.getValue());
                    continue;
                }

                Set<PointPairInt> tmpJoinPointSet = new HashSet<PointPairInt>();
                for (PointPairInt joinPoint : entry.getValue()) {
                    if (!rmJoinPointSet.contains(joinPoint)) {
                        tmpJoinPointSet.add(joinPoint);
                    }
                }
                tmpMap.put(entry.getKey(), tmpJoinPointSet);
            }
            edgeLastEndPointMap.clear();
            edgeLastEndPointMap.putAll(tmpMap);
        }

    }

    public PairIntArray findAsSingleClosedEdge() {

        this.singleClosedEdge = true;

        List<PairIntArray> output = connectPixelsViaDFSForBounds();

if (debug) {
Image img0 = ImageIOHelper.convertImage(img);
for (int i = 0; i < output.size(); ++i) {
    PairIntArray pa = output.get(i);
    for (int j = 0; j < pa.getN(); ++j) {
        int x = pa.getX(j);
        int y = pa.getY(j);
        if (i == 0) {
            if (j == 0 || (j == (pa.getN() - 1))) {
                ImageIOHelper.addPointToImage(x, y, img0, 0, 200, 100, 0);
            } else {
                ImageIOHelper.addPointToImage(x, y, img0, 0, 255, 0, 0);
            }
        } else if (i == 1) {
            ImageIOHelper.addPointToImage(x, y, img0, 0, 0, 255, 0);
        } else {
            ImageIOHelper.addPointToImage(x, y, img0, 0, 0, 0, 255);
        }
    }
}
MiscDebug.writeImageCopy(img0, "output_before_merges_" + MiscDebug.getCurrentTimeFormatted() + ".png");
}
        output = findEdgesIntermediateSteps(output);

if (debug) {
Image img2 = ImageIOHelper.convertImage(img);
for (int i = 0; i < output.size(); ++i) {
    PairIntArray pa = output.get(i);
    for (int j = 0; j < pa.getN(); ++j) {
        int x = pa.getX(j);
        int y = pa.getY(j);
        if (i == 0) {
            if (j == 0 || (j == (pa.getN() - 1))) {
                ImageIOHelper.addPointToImage(x, y, img2, 0, 200, 100, 0);
            } else {
                ImageIOHelper.addPointToImage(x, y, img2, 0, 255, 0, 0);
            }
        } else if (i == 1) {
            ImageIOHelper.addPointToImage(x, y, img2, 0, 0, 255, 0);
        } else {
            ImageIOHelper.addPointToImage(x, y, img2, 0, 0, 0, 255);
        }
    }
}
MiscDebug.writeImageCopy(img2, "output_after_merges_" + MiscDebug.getCurrentTimeFormatted() + ".png");
}
        if (output.size() > 1) {
            Collections.sort(output, new PairIntArrayDescendingComparator());
        }

        PairIntArray out = output.get(0);

        reorderEndpointsIfNeeded(out);

        if (isAdjacent(out, 0, out.getN() - 1)) {
            PairIntArrayWithColor p = new PairIntArrayWithColor(out);
            p.setColor(1);
            out = p;
        }
if (debug) {
Image img2 = ImageIOHelper.convertImage(img);
for (int j = 0; j < out.getN(); ++j) {
    int x = out.getX(j);
    int y = out.getY(j);
    if (j == 0 || (j == (out.getN() - 1))) {
        ImageIOHelper.addPointToImage(x, y, img2, 0, 200, 100, 0);
    } else {
        ImageIOHelper.addPointToImage(x, y, img2, 0, 255, 0, 0);
    }
}
MiscDebug.writeImageCopy(img2, "output_after_reorder_endpoints_" + MiscDebug.getCurrentTimeFormatted() + ".png");
}
        return out;
    }

    private Map<Integer, Set<Integer>> createEdgeToPixelIndexMap() {

        // key = edge index.
        // value = pixel indexes.
        //    the pixel indexes are used to find values in junctionLocatorMap
        //    to update it as points are moved to and from edges.
        Map<Integer, Set<Integer>> theEdgeToPixelIndexMap = new HashMap<Integer, Set<Integer>>();

        for (Entry<Integer, PairInt> entry : junctionLocationMap.entrySet()) {
            Integer pixelIndex = entry.getKey();
            PairInt loc = entry.getValue();
            Integer edgeIndex = Integer.valueOf(loc.getX());
            Set<Integer> pixelIndexes = theEdgeToPixelIndexMap.get(edgeIndex);
            if (pixelIndexes == null) {
                pixelIndexes = new HashSet<Integer>();
                theEdgeToPixelIndexMap.put(edgeIndex, pixelIndexes);
            }
            pixelIndexes.add(pixelIndex);
        }

        return theEdgeToPixelIndexMap;
    }

    private int insertAdjacentForClosedCurve(List<PairIntArray> edges) {

        if (edges == null || edges.isEmpty()) {
            return 0;
        }

        Collections.sort(edges, new PairIntArrayDescendingComparator());

        int maxEdgeIdx = 0;

        findJunctions(edges);

        // key = edge index.
        // value = pixel indexes.
        //     the pixel indexes are used to find values in junctionLocatorMap
        //     to update it as points are moved to and from edges.
        Map<Integer, Set<Integer>> theEdgeToPixelIndexMap = createEdgeToPixelIndexMap();

        int nChanged = 0;

        for (int i = (edges.size() - 1); i > -1; --i) {

            if (i == maxEdgeIdx || (edges.get(i).getN() == 0)) {
                continue;
            }

            // i may be merged into maxEdgeIdx
            int nInserted = insert(theEdgeToPixelIndexMap, edges, maxEdgeIdx, i);

            if (nInserted > 0) {

                findJunctions(edges);

                theEdgeToPixelIndexMap = createEdgeToPixelIndexMap();
            }

            nChanged += nInserted;
        }

        int nMerged = 0;
        int nIter = 0;
        int nMaxIter = 10;
        if (edges.size() > nMaxIter) {
            nMaxIter = edges.size();
        }

        while ((nIter == 0) || ((nIter < nMaxIter) && (nMerged > 0))) {

            if (nMerged > 0) {
                findJunctions(edges);
                theEdgeToPixelIndexMap = createEdgeToPixelIndexMap();
            }

            nMerged = 0;

            for (int i = 0; i < edges.size(); ++i) {

                //TODO: consider skipping single pixel edges because they're
                // handled in insertAdjacentForClosedCurve?
                if (edges.get(i).getN() == 0) {
                    continue;
                }

                PairIntArray curve = edges.get(i);

                int x1 = curve.getX(0);
                int y1 = curve.getY(0);
                int xn = curve.getX(curve.getN() - 1);
                int yn = curve.getY(curve.getN() - 1);

                // uses the junctionLocationMap to search for neighbors of given point
                Set<PairInt> adjToFirstLoc = findAdjacentInOtherEdge(x1, y1, i);

                Set<PairInt> adjToLastLoc = findAdjacentInOtherEdge(xn, yn, i);

                /*
                See if there's an unambiguous closest among the 2 sets,
                and choose that if there is
                */

                double minDist = Double.MAX_VALUE;
                PairInt closest = null;
                for (PairInt loc : adjToFirstLoc) {
                    int x = edges.get(loc.getX()).getX(loc.getY());
                    int y = edges.get(loc.getX()).getY(loc.getY());
                    int diffX = x1 - x;
                    int diffY = y1 - y;
                    double dist = Math.sqrt(diffX*diffX + diffY*diffY);
                    if (dist == minDist) {
                        // discard because of ambiguity
                        closest = null;
                        break;
                    } else if (dist < minDist) {
                        minDist = dist;
                        closest = loc;
                    }
                }

                if (closest != null) {
                    //edge i is the one being inserted into by edge 2
                    int ins = insert(theEdgeToPixelIndexMap, edges, i, closest.getX());

                    /*
                    TODO: consider this instead:
                    int nInserted = insertForEdge1FirstLocation(edges, i, closest.getX(), set);
                    */
                    if (ins > 0) {
                        return ins;
                    }
                } else {
                    for (PairInt loc : adjToFirstLoc) {
                        int edge2Idx = loc.getX();
                        int ins = insert(theEdgeToPixelIndexMap, edges, i, edge2Idx);
                        if (ins > 0) {
                            return ins;
                        }
                    }
                }

                closest = null;
                minDist = Double.MAX_VALUE;
                for (PairInt loc : adjToLastLoc) {
                    int x = edges.get(loc.getX()).getX(loc.getY());
                    int y = edges.get(loc.getX()).getY(loc.getY());
                    int diffX = xn - x;
                    int diffY = yn - y;
                    double dist = Math.sqrt(diffX * diffX + diffY * diffY);
                    if (dist == minDist) {
                        // discard because of ambiguity
                        closest = null;
                        break;
                    } else if (dist < minDist) {
                        minDist = dist;
                        closest = loc;
                    }
                }

                if (closest != null) {
                    //edge i is the one being inserted into by edge 2
                    int ins = insert(theEdgeToPixelIndexMap, edges, i, closest.getX());
                    /*
                    TODO: consider instead
                    int nInserted = insertForEdge1LastLocation(edges, i, closest.getX(), set);
                    */
                    if (ins > 0) {
                        return ins;
                    }
                } else {
                    for (PairInt loc : adjToLastLoc) {
                        int edge2Idx = loc.getX();
                        int ins = insert(theEdgeToPixelIndexMap, edges, i, edge2Idx);
                        if (ins > 0) {
                            return ins;
                        }
                    }
                }
            }

            ++nIter;
        }

        return nChanged;
    }

    private int insertBetweenAdjacentForClosedCurve(List<PairIntArray> output) {

        //NOT YET IMPLEMENTED

        if (output == null || output.isEmpty()) {
            return 0;
        }

        int nChanged = 0;

        /*
        unlike insertAdjacentForClosedCurve, this method does not attempt to
        join all curves exclusively to the largest curve.
        instead it attempts to join a curve to the adjacent curves at its
        endpoints.

        previous join points and junction splice operations should have joined
        most curves, but there may be the case where two curves
        have an open curve in between them and neither of the two would ever gain
        length by adding it, so all 3 remain separate.

        Here's a case:
          1 1 1
                1
                  1
          1 1 1 1   3 3
                  2     3 3
                    2 3     3
                        3  3

        For each curve, looks at it's endpoints and searches the junction
        location map to find adjacent edge points (junctions).

        It finds the closest for both endpoints.

        Then it looks to see if those 2 edges have a junction.
        (If there is only one edge it is adjacent to, that's handled in the
        other method, insertAdjacentForClosedCurve).

        If those 2 edges are adjacent, and if that adjacency is adjacent
        to one of the endpoints of the open curve in between,
        the 3 are joined.

        */

        for (int i = 0; i < output.size(); ++i) {

            //TODO: consider skipping single pixel edges because they're handled in insertAdjacentForClosedCurve?
            if (output.get(i).getN() == 0) {
                continue;
            }

            PairIntArray curve = output.get(i);

            int x1 = curve.getX(0);
            int y1 = curve.getY(0);
            int xn = curve.getX(curve.getN() - 1);
            int yn = curve.getY(curve.getN() - 1);
            /*
            PairInt[] adjLoc = findTwoAdjacent(i, x1, y1, xn, yn);

            if (firstLastAdjLoc == null) {
                continue;
            }

            // skipping for insert into same edge because it's handled in another method
            if (firstLastAdjLoc[0].getX() == firstLastAdjLoc[1].getX()) {
                continue;
            }

            check that adjToFirstPointLoc is adjacent to adjToLastPointLoc
            */
            /*
            Here's a case:
              1 1 1
                    1
                      1
              1 1 1 1   3 3
                      2     3 3
                        2 3     3
                            3  3
            */

        }

        return nChanged;
    }

    /**
     * populate the output lists with adjacent pixel locations for
     * the edge next to the reference edge (where reference edge is usually the
     * one with largest number of points).
     *
     * @param theEdgeToPixelIndexMap
     * @param edge1Idx
     * @param edge2Idx
     * @param outputAdjacentToEdge1LocMap key is a location in edge1, value is
     * a set of locations in edge2 adjacent to the key
     */
    private void populateWithAdjacentLocations(
        Map<Integer, Set<Integer>> theEdgeToPixelIndexMap,
        final int edge1Idx, final int edge2Idx,
        Map<PairInt, Set<PairInt>> outputAdjacentToEdge1LocMap) {

        Set<Integer> edge1PixIndexes = theEdgeToPixelIndexMap.get(Integer.valueOf(edge1Idx));

        for (Integer edge1PixIndex : edge1PixIndexes) {

            PairInt edge1Loc = junctionLocationMap.get(edge1PixIndex);

            // search the 8 neighbors to see if any are in maxEdgePixelIndexes
            int x = img.getCol(edge1PixIndex);
            int y = img.getRow(edge1PixIndex);

            for (int i = 0; i < dxs8.length; ++i) {

                int x2 = x + dxs8[i];
                int y2 = y + dys8[i];
                if (x2 < 0 || (y2 < 0) || (x2 > (img.getWidth() - 1) || (y2 > (img.getHeight() - 1)))) {
                    continue;
                }

                int pixIdx2 = img.getIndex(x2, y2);

                // if this is in edge2, add it to the set
                PairInt adjacentToEdge1Loc = junctionLocationMap.get(pixIdx2);

                if ((adjacentToEdge1Loc != null) && (adjacentToEdge1Loc.getX() == edge2Idx)) {

                    // edge1JunctionLoc is next to edge1Loc
                    Set<PairInt> edge2LocSet = outputAdjacentToEdge1LocMap.get(edge1Loc);
                    if (edge2LocSet == null) {
                        edge2LocSet = new HashSet<PairInt>();
                        outputAdjacentToEdge1LocMap.put(edge1Loc, edge2LocSet);
                    }
                    edge2LocSet.add(adjacentToEdge1Loc);
                }
            }
        }
    }

    private void reorderEndpointsIfNeeded(PairIntArray out) {

        if (out == null || out.getN() < 3) {
            return;
        }

        if (isAdjacent(out, 0, out.getN() - 1)) {
            return;
        }

        int x1 = out.getX(0);
        int y1 = out.getY(0);

        int xn = out.getX(out.getN() - 1);
        int yn = out.getY(out.getN() - 1);

        int closestPivotIdx = -1;

        for (int idx = (out.getN() - 2); idx > 0; --idx) {

            int x = out.getX(idx);
            int y = out.getY(idx);

            if ((Math.abs(x - x1) < 2) && (Math.abs(y - y1) < 2)) {
                closestPivotIdx = idx;
                break;
            }
        }

        if (closestPivotIdx == -1) {
            // it's not a closed curve
            return;
        }

        if (closestPivotIdx > (out.getN() - closestPivotIdx - 1)) {
            // closest to end
            /*
            If point before closestPivotIdx is next to (x1, y1), can just
            reverse the points from closestPivotIdx to n-1.
            */
            int xPrev = out.getX(closestPivotIdx - 1);
            int yPrev = out.getY(closestPivotIdx - 1);
            if ((Math.abs(xPrev - xn) > 1) || (Math.abs(yPrev - yn) > 1)) {
                return;
            }

            int count = 0;
            int nSep = (out.getN() - closestPivotIdx) >> 1;
            for (int idx = closestPivotIdx; idx < (closestPivotIdx + nSep); ++idx) {
                int idx2 = out.getN() - count - 1;
                int swapX = out.getX(idx);
                int swapY = out.getY(idx);
                out.set(idx, out.getX(idx2), out.getY(idx2));
                out.set(idx2, swapX, swapY);
                count++;
            }

        } else {
            // closest to beginning
            /*
            If point after closestPivotIdx is next to (xn, yn), can reverse the
            points from 0 to closestPivotIdx.
            */

            closestPivotIdx = -1;

            for (int idx = 1; idx < (out.getN() - 1); ++idx) {

                int x = out.getX(idx);
                int y = out.getY(idx);

                if ((Math.abs(x - xn) < 2) && (Math.abs(y - yn) < 2)) {
                    closestPivotIdx = idx;
                    break;
                }
            }

            if (closestPivotIdx == -1) {
                // it's not a closed curve
                return;
            }

            int xNext = out.getX(closestPivotIdx + 1);
            int yNext = out.getY(closestPivotIdx + 1);
            if ((Math.abs(xNext - xn) > 1) || (Math.abs(yNext - yn) > 1)) {
                return;
            }

            out.reverse0toIdx(closestPivotIdx);
        }
    }

    private int insert(Map<Integer, Set<Integer>> theEdgeToPixelIndexMap,
        List<PairIntArray> output, int edge1Idx, int edge2Idx) {

        /*
        edge 1 is the one being inserted into by edge 2
        */

        PairIntArray edge1 = output.get(edge1Idx);

        PairIntArray edge2 = output.get(edge2Idx);

        /*
        key is a location in edge1,
        value is a set of locations in edge2 adjacent to the key
        */
        Map<PairInt, Set<PairInt>> adjToEdge1LocMap = new HashMap<PairInt, Set<PairInt>>();
        populateWithAdjacentLocations(theEdgeToPixelIndexMap, edge1Idx, edge2Idx,
            adjToEdge1LocMap);

        if (adjToEdge1LocMap.isEmpty()) {
            return 0;
        }

        // if can insert at index = 0 or last index for edge1, prefer that
        Set<PairInt> edge2LocForEdge1FirstLoc = null;
        Set<PairInt> edge2LocForEdge1LastLoc = null;

        for (Entry<PairInt, Set<PairInt>> entry : adjToEdge1LocMap.entrySet()) {
            PairInt edge1Loc = entry.getKey();
            if (edge1Loc.getY() == 0) {
                edge2LocForEdge1FirstLoc = entry.getValue();
            } else if (edge1Loc.getY() == (edge1.getN() - 1)) {
                edge2LocForEdge1LastLoc = entry.getValue();
            }
        }

        if ((edge2LocForEdge1FirstLoc != null) && !edge2LocForEdge1FirstLoc.isEmpty()) {

            int nIns = insertForEdge1FirstLocation(output, edge1Idx, edge2Idx,
                edge2LocForEdge1FirstLoc);

            if (nIns > 0) {
                return nIns;
            }
        }

        if ((edge2LocForEdge1LastLoc != null) && !edge2LocForEdge1LastLoc.isEmpty()) {

            int nIns = insertForEdge1LastLocation(output, edge1Idx, edge2Idx,
                edge2LocForEdge1LastLoc);

            if (nIns > 0) {
                return nIns;
            }
        }

        /*
        if have arrived here, the locations in edge1 are not at its endpoints.

        can insert another curve if both insertion points in edge1 are next to
        one another.

        edge1 is '@''s
                               @ @ @ @
                             @
                           @
                    0 1 2 @
                    5 4 3 @
                            @ @
        */

        if (adjToEdge1LocMap.size() < 2) {
            // there aren't 2 points to
            return 0;
        }

        /*
        key is a location in edge1,
        value is a set of locations in edge2 adjacent to the key
        Map<PairInt, Set<PairInt>> adjToEdge1LocMap
        */

        // ---- looking for 2 keys in adjToEdge1LocMap that are adjacent to
        //      one another

        List<PairInt> edge1ALocKeys = new ArrayList<PairInt>();
        List<PairInt> edge1BLocKeys = new ArrayList<PairInt>();

        for (Entry<PairInt, Set<PairInt>> entryA : adjToEdge1LocMap.entrySet()) {

            PairInt edge1ALoc = entryA.getKey();

            assert(edge1ALoc.getX() == edge1Idx);

            int xA = edge1.getX(edge1ALoc.getY());
            int yA = edge1.getY(edge1ALoc.getY());

            for (Entry<PairInt, Set<PairInt>> entryB : adjToEdge1LocMap.entrySet()) {

                PairInt edge1BLoc = entryB.getKey();

                assert(edge1BLoc.getX() == edge1Idx);

                if (edge1ALoc.equals(edge1BLoc)) {
                    continue;
                }

                int xB = edge1.getX(edge1BLoc.getY());
                int yB = edge1.getY(edge1BLoc.getY());

                int diffX = Math.abs(xA - xB);
                if (diffX > 1) {
                    continue;
                }
                int diffY = Math.abs(yA - yB);
                if (diffY > 1) {
                    continue;
                }

                int idxA = edge1ALocKeys.indexOf(edge1BLoc);
                if (idxA > -1) {
                    if (edge1BLocKeys.get(idxA).equals(edge1ALoc)) {
                        continue;
                    }
                }

                edge1ALocKeys.add(edge1ALoc);
                edge1BLocKeys.add(edge1BLoc);
            }
        }

        if (edge1ALocKeys.isEmpty()) {
            return 0;
        }

        /*
        if edge2 are endpoints can simply insert
        else have to see if can re-order edge 2 points so
        that the closest points are both placable at endpoints
        */

        List<PairInt> edge1ALocClosestEdge2Loc = new ArrayList<PairInt>();
        List<PairInt> edge1BLocClosestEdge2Loc = new ArrayList<PairInt>();

        for (int i = 0; i < edge1ALocKeys.size(); ++i) {
            PairInt edge1ALoc = edge1ALocKeys.get(i);
            int x1 = edge1.getX(edge1ALoc.getY());
            int y1 = edge1.getY(edge1ALoc.getY());

            Set<PairInt> edge2LocSet = adjToEdge1LocMap.get(edge1ALoc);
            double minDist = Double.MAX_VALUE;
            PairInt closest = null;
            for (PairInt edge2Loc : edge2LocSet) {
                int x2 = edge2.getX(edge2Loc.getY());
                int y2 = edge2.getY(edge2Loc.getY());
                int diffX = Math.abs(x1 - x2);
                int diffY = Math.abs(y1 - y2);
                double dist = Math.sqrt(diffX*diffX + diffY*diffY);
                if (dist < minDist) {
                    minDist = dist;
                    closest = edge2Loc;
                }
            }
            edge1ALocClosestEdge2Loc.add(closest);
        }

        for (int i = 0; i < edge1BLocKeys.size(); ++i) {
            PairInt edge1BLoc = edge1BLocKeys.get(i);
            int x1 = edge1.getX(edge1BLoc.getY());
            int y1 = edge1.getY(edge1BLoc.getY());

            Set<PairInt> edge2LocSet = adjToEdge1LocMap.get(edge1BLoc);
            double minDist = Double.MAX_VALUE;
            PairInt closest = null;
            for (PairInt edge2Loc : edge2LocSet) {
                int x2 = edge2.getX(edge2Loc.getY());
                int y2 = edge2.getY(edge2Loc.getY());
                int diffX = Math.abs(x1 - x2);
                int diffY = Math.abs(y1 - y2);
                double dist = Math.sqrt(diffX*diffX + diffY*diffY);
                if (dist < minDist) {
                    minDist = dist;
                    closest = edge2Loc;
                }
            }
            edge1BLocClosestEdge2Loc.add(closest);
        }

        for (int i = 0; i < edge1ALocKeys.size(); ++i) {

            PairInt edge1ALoc = edge1ALocKeys.get(i);
            PairInt edge1BLoc = edge1BLocKeys.get(i);

            PairInt edge2ALoc = edge1ALocClosestEdge2Loc.get(i);
            PairInt edge2BLoc = edge1BLocClosestEdge2Loc.get(i);

            if ((edge2ALoc.getY() == 0) && (edge2BLoc.getY() == (edge2.getN() - 1))) {

                if (edge1ALoc.getY() < edge1BLoc.getY()) {
                    edge1.insertAll(edge1ALoc.getY() + 1, edge2);
                    output.remove(edge2Idx);
                    return 1;
                } else {
                    edge2.reverse();
                    edge1.insertAll(edge1BLoc.getY() + 1, edge2);
                    output.remove(edge2Idx);
                    return 1;
                }

            } else if ((edge2BLoc.getY() == 0) && (edge2ALoc.getY() == (edge2.getN() - 1))) {
                if (edge1ALoc.getY() < edge1BLoc.getY()) {
                    edge2.reverse();
                    edge1.insertAll(edge1ALoc.getY() + 1, edge2);
                    output.remove(edge2Idx);
                    return 1;
                } else {
                    edge1.insertAll(edge1BLoc.getY() + 1, edge2);
                    output.remove(edge2Idx);
                    return 1;
                }

            } else {

                // --- if here, then edge1 and edge2 closest points weren't
                //     endpoints so we try to reorder points to make the
                //     curves joinable sequentially

                boolean edge1IsClosed = isAdjacent(edge1, 0, edge1.getN() - 1);

                boolean edge2IsClosed = isAdjacent(edge2, 0, edge2.getN() - 1);

                int nReversed = 0;

                if (output.size() == 2 && !edge1IsClosed) {

                    // attempt to reorder and then set edge1IsClosed to true

                    int edge1EndIdx = findAdjacentToTopAtBottom(edge1);

                    boolean didReverse = false;

                    if (edge1EndIdx > -1) {
   // TODO: invoke even if edge1EndIdx << (edge1.getN()/2) ?
                        didReverse = reverseBottomIfPossible(edge1, edge1EndIdx);
                        if (didReverse) {
                            if (edge1ALoc.getY() > edge1EndIdx) {
                                edge1ALoc.setY(edge1.getN() - edge1ALoc.getY() + edge1EndIdx - 1);
                            }
                            if (edge1BLoc.getY() > edge1EndIdx) {
                                edge1BLoc.setY(edge1.getN() - edge1BLoc.getY() + edge1EndIdx - 1);
                            }
                            edge1IsClosed = isAdjacent(edge1, 0, edge1.getN() - 1);
                            nReversed++;

                        } else {

                            if (edge1EndIdx < (edge1.getN()/2)) {

                                if (reverseTopIfPossible(edge1, edge1EndIdx - 1)) {
                                    if (edge1ALoc.getY() < (edge1EndIdx - 1)) {
                                        edge1ALoc.setY(edge1EndIdx - 1 - edge1ALoc.getY());
                                    }
                                    if (edge1BLoc.getY() < (edge1EndIdx - 1)) {
                                        edge1BLoc.setY(edge1EndIdx - 1 - edge1BLoc.getY());
                                    }
                                    nReversed++;
                                    edge1IsClosed = isAdjacent(edge1, 0, edge1.getN() - 1);
                                }
                            }
                        }
                    }

                    if (!didReverse) {

                        int edge1BeginIdx = findAdjacentToBottomAtTop(edge1);

                        if (edge1BeginIdx > -1) {
                            didReverse = reverseTopIfPossible(edge1, edge1BeginIdx);
                            if (didReverse) {
                                if (edge1ALoc.getY() < edge1BeginIdx) {
                                    edge1ALoc.setY(edge1BeginIdx - edge1ALoc.getY());
                                }
                                if (edge1BLoc.getY() < edge1BeginIdx) {
                                    edge1BLoc.setY(edge1BeginIdx - edge1BLoc.getY());
                                }
                                edge1IsClosed = isAdjacent(edge1, 0, edge1.getN() - 1);
                                nReversed++;
                            } else {

                                if (edge1BeginIdx > (edge1.getN()/2)) {
                                    if (reverseBottomIfPossible(edge1, edge1BeginIdx + 1)) {
                                        if (edge1ALoc.getY() > (edge1BeginIdx + 1)) {
                                            edge1ALoc.setY(edge1.getN() - edge1ALoc.getY() + edge1BeginIdx);
                                        }
                                        if (edge1BLoc.getY() > (edge1BeginIdx + 1)) {
                                            edge1BLoc.setY(edge1.getN() - edge1BLoc.getY() + edge1BeginIdx);
                                        }
                                        nReversed++;
                                        edge1IsClosed = isAdjacent(edge1, 0, edge1.getN() - 1);
                                    }
                                }
                            }
                        }
                    }
                }

                if (output.size() == 2 && edge1IsClosed && !edge2IsClosed) {

                    // attempt to reorder and then set edge2IsClosed to true

                    int edge2EndIdx = findAdjacentToTopAtBottom(edge2);

                    boolean didReverse = false;

                    if (edge2EndIdx > -1) {

                        didReverse = reverseBottomIfPossible(edge2, edge2EndIdx);
                        if (didReverse) {
                            if (edge2ALoc.getY() > edge2EndIdx) {
                                edge2ALoc.setY(edge2.getN() - edge2ALoc.getY() + edge2EndIdx - 1);
                            }
                            if (edge2BLoc.getY() > edge2EndIdx) {
                                edge2BLoc.setY(edge2.getN() - edge2BLoc.getY() + edge2EndIdx - 1);
                            }
                            edge2IsClosed = isAdjacent(edge2, 0, edge2.getN() - 1);
                            nReversed++;
                        } else {

                            if (edge2EndIdx < (edge2.getN()/2)) {

                                if (reverseTopIfPossible(edge2, edge2EndIdx - 1)) {
                                    if (edge2ALoc.getY() < (edge2EndIdx - 1)) {
                                        edge2ALoc.setY(edge2EndIdx - 1 - edge2ALoc.getY());
                                    }
                                    if (edge2BLoc.getY() < (edge2EndIdx - 1)) {
                                        edge2BLoc.setY(edge2EndIdx - 1 - edge2BLoc.getY());
                                    }
                                    nReversed++;
                                    edge2IsClosed = isAdjacent(edge2, 0, edge2.getN() - 1);
                                }
                            }
                        }
                    }

                    if (!didReverse) {

                        int edge2BeginIdx = findAdjacentToBottomAtTop(edge2);

                        if (edge2BeginIdx > -1) {
                            didReverse = reverseTopIfPossible(edge2, edge2BeginIdx);
                            if (didReverse) {
                                if (edge2ALoc.getY() < edge2BeginIdx) {
                                    edge2ALoc.setY(edge2BeginIdx - edge2ALoc.getY());
                                }
                                if (edge2BLoc.getY() < edge2BeginIdx) {
                                    edge2BLoc.setY(edge2BeginIdx - edge2BLoc.getY());
                                }
                                edge2IsClosed = isAdjacent(edge2, 0, edge2.getN() - 1);
                                nReversed++;
                            } else {
                                if (edge2BeginIdx > (edge2.getN()/2)) {
                                    if (reverseBottomIfPossible(edge2, edge2BeginIdx + 1)) {
                                        if (edge2ALoc.getY() > (edge2BeginIdx + 1)) {
                                            edge2ALoc.setY(edge2.getN() - edge2ALoc.getY() + edge2BeginIdx);
                                        }
                                        if (edge2BLoc.getY() > (edge2BeginIdx + 1)) {
                                            edge2BLoc.setY(edge2.getN() - edge2BLoc.getY() + edge2BeginIdx);
                                        }
                                        nReversed++;
                                        edge2IsClosed = isAdjacent(edge2, 0, edge2.getN() - 1);
                                    }
                                }
                            }
                        }
                    }
                }

                if (edge1IsClosed && edge2IsClosed) {

                    if (edge1ALoc.getY() > edge1BLoc.getY()) {
                        if (edge2ALoc.getY() > edge2BLoc.getY()) {
                            circularlyShift(edge1, (edge1.getN() - edge1BLoc.getY() - 1));
                            circularlyShift(edge2, (edge2.getN() - edge2BLoc.getY() - 1));
                            edge2.reverse();
                            edge1.addAll(edge2);
                            output.remove(edge2Idx);
                            return 1;
                        } else {
                            circularlyShift(edge1, (edge1.getN() - edge1BLoc.getY() - 1));
                            circularlyShift(edge2, (edge2.getN() - edge2ALoc.getY() - 1));
                            edge1.addAll(edge2);
                            output.remove(edge2Idx);
                            return 1;
                        }

                    } else {
                        if (edge2ALoc.getY() > edge2BLoc.getY()) {
                            circularlyShift(edge1, (edge1.getN() - edge1ALoc.getY() - 1));
                            circularlyShift(edge2, (edge2.getN() - edge2BLoc.getY() - 1));
                            edge1.addAll(edge2);
                            output.remove(edge2Idx);
                            return 1;
                        } else {
                            circularlyShift(edge1, (edge1.getN() - edge1ALoc.getY() - 1));
                            circularlyShift(edge2, (edge2.getN() - edge2ALoc.getY() - 1));
                            edge2.reverse();
                            edge1.addAll(edge2);
                            output.remove(edge2Idx);
                            return 1;
                        }
                    }
                } else if (nReversed > 0) {
                    return 1;
                }
            }
        }

        return 0;
    }

    private Set<PairInt> findAdjacentInOtherEdge(final int x, final int y,
        final int edgeIdx) {

        Set<PairInt> locs = new HashSet<PairInt>();

        for (int i = 0; i < dxs8.length; ++i) {

            int x2 = x + dxs8[i];
            int y2 = y + dys8[i];

            if ((x2 < 0) || (y2 < 0) || (x2 > (img.getWidth() - 1)) ||
                (y2 > (img.getHeight() - 1))) {
                continue;
            }

            Integer pixIndex2 = img.getIndex(x2, y2);

            PairInt loc2 = junctionLocationMap.get(pixIndex2);

            if (loc2 != null && (loc2.getX() != edgeIdx)) {
                locs.add(loc2);
            }
        }

        return locs;
    }

    private boolean isAdjacent(PairIntArray edge, int idx1, int idx2) {

        int x1 = edge.getX(idx1);
        int y1 = edge.getY(idx1);

        int x2 = edge.getX(idx2);
        int y2 = edge.getY(idx2);

        int diffX = Math.abs(x1 - x2);
        int diffY = Math.abs(y1 - y2);

        if ((diffX < 2) && (diffY < 2)) {
            return true;
        }

        return false;
    }

    private boolean reverseTopIfPossible(PairIntArray edge, int cIdx) {

        /*
            reversing beginning of edge:

                                  0          Can reverse segment 0 to cIdx if
                                  1          the point at index 0 is
                                 cIdx - 1    adjacent to cIdx + 1
            edge1 0      edge2   cIdx
                                 cIdx + 1
        */

        boolean isAdjacent = isAdjacent(edge, 0, cIdx + 1);

        if (isAdjacent) {

            edge.reverse0toIdx(cIdx);

            return true;
        }

        return false;
    }

    private boolean reverseBottomIfPossible(PairIntArray edge, int cIdx) {

        /*            ----  cIdx - 1 ---
            edge1 0      edge2 cIdx          Can reverse segment cIdx to end if
                               cIdx + 1      the point at n-1 is adjacent to
                               cIdx + 2      the point at cIdx - 1
                               n - 1
        */

        boolean isAdjacent = isAdjacent(edge, cIdx - 1, edge.getN() - 1);

        if (isAdjacent) {

            edge.reverseIdxtoEnd(cIdx);

            return true;
        }

        return false;
    }

    private void circularlyShift(PairIntArray edge, int positiveNumber) {

        if (positiveNumber == 0) {
            return;
        } else if (positiveNumber < 0) {
            throw new IllegalArgumentException(
            "method currently requires positiveNumber to be a positive number");
        }

        Rotate rotate = new Rotate();

        //TODO: implement the faster logic for shift rights in rotate2
        rotate.rotate2(edge.getX(), edge.getN(), -1*positiveNumber);

        rotate.rotate2(edge.getY(), edge.getN(), -1*positiveNumber);

    }

    private int insertForEdge1LastLocation(List<PairIntArray> edges,
        final int edge1Idx, final int edge2Idx, Set<PairInt> edge2LocForEdge1LastLoc) {

        PairIntArray edge1 = edges.get(edge1Idx);

        PairIntArray edge2 = edges.get(edge2Idx);

         /*
         have found locations in edge2 close to the last point in edge1

         (1) if edge2 location is an edge2 endpoint, the insert is simply
         an append or a reverse and append
         (2) ELSE try 2 tests for re-ordering edge2

         (2a) reversing end of edge2:

         ----  cIdx - 1 ---
         edge1 (n-1)    edge2 cIdx          Can reverse segment cIdx to end if
                              cIdx + 1      the point at n-1 is adjacent to
                              cIdx + 2      the point at cIdx - 1
                              n - 1

         (2a) reversing beginning of edge2:

                               0          Can reverse segment 0 to cIdx if
                               1          the point at index 0 is
                              cIdx - 1      adjacent to cIdx + 1
       edge1 (n-1)    edge2   cIdx
                              cIdx + 1
         */
        for (PairInt edge2Loc : edge2LocForEdge1LastLoc) {
            if (edge2Loc.getY() == 0) {
                edge1.addAll(edge2);
                edges.remove(edge2Idx);
                return 1;
            } else if (edge2Loc.getY() == (edge2.getN() - 1)) {
                edge2.reverse();
                edge1.addAll(edge2);
                edges.remove(edge2Idx);
                return 1;
            }
        }

        // prefer the closest to the edge1 point
        int xn = edge1.getX(edge1.getN() - 1);
        int yn = edge1.getY(edge1.getN() - 1);
        PairInt closest = null;
        double closestDist = Integer.MAX_VALUE;
        for (PairInt edge2Loc : edge2LocForEdge1LastLoc) {
            assert(edge2Loc.getX() == edge2Idx);
            int x2 = edge2.getX(edge2Loc.getY());
            int y2 = edge2.getY(edge2Loc.getY());
            int diffX = Math.abs(xn - x2);
            int diffY = Math.abs(yn - y2);
            double dist = Math.sqrt(diffX * diffX + diffY * diffY);
            if (dist < closestDist) {
                closestDist = dist;
                closest = edge2Loc;
            }
        }

        boolean didReverse = reverseTopIfPossible(edge2, closest.getY());

        if (didReverse) {
            edge1.addAll(edge2);
            edges.remove(edge2Idx);
            return 1;
        } else {
            didReverse = reverseBottomIfPossible(edge2, closest.getY());
            if (didReverse) {
                edge2.reverse();
                edge1.addAll(edge2);
                edges.remove(edge2Idx);
                return 1;
            }
        }

        // if arrived here, try the points in set that are not same as 'closest'
        for (PairInt edge2Loc : edge2LocForEdge1LastLoc) {

            if (edge2Loc.equals(closest)) {
                continue;
            }

            int idx = edge2Loc.getY();

            // same checks for whether can re-order to have an endpoint
            didReverse = reverseTopIfPossible(edge2, idx);

            if (didReverse) {
                edge1.addAll(edge2);
                edges.remove(edge2Idx);
                return 1;
            } else {
                didReverse = reverseBottomIfPossible(edge2, idx);
                if (didReverse) {
                    edge2.reverse();
                    edge1.addAll(edge2);
                    edges.remove(edge2Idx);
                    return 1;
                }
            }
        }

        return 0;
    }

    private int insertForEdge1FirstLocation(List<PairIntArray> edges,
        final int edge1Idx, final int edge2Idx, Set<PairInt> edge2LocForEdge1FirstLoc) {

        PairIntArray edge1 = edges.get(edge1Idx);

        PairIntArray edge2 = edges.get(edge2Idx);

        /*
        have found locations in edge2 close to the first point in edge1

        (1) if edge2 location is an edge2 endpoint, the insert is simply
            and append or a reverse and append
        (2) ELSE try 2 tests for re-ordering edge2

        (2a) reversing end of edge2:

                     ----  cIdx - 1 ---
        edge1 0      edge2 cIdx          Can reverse segment cIdx to end if
                           cIdx + 1      the point at n-1 is adjacent to
                           cIdx + 2      the point at cIdx - 1
                           n - 1

        (2a) reversing beginning of edge2:

                              0          Can reverse segment 0 to cIdx if
                              1          the point at index 0 is
                             cIdx - 1    adjacent to cIdx + 1
        edge1 0      edge2   cIdx
                             cIdx + 1
        */

        for (PairInt edge2Loc : edge2LocForEdge1FirstLoc) {
            if (edge2Loc.getY() == 0) {
                edge2.reverse();
                edge1.insertAll(0, edge2);
                edges.remove(edge2Idx);
                return 1;
            } else if (edge2Loc.getY() == (edge2.getN() - 1)) {
                edge2.addAll(edge1);
                edges.set(edge1Idx, edge2);
                edges.remove(edge2Idx);
                return 1;
            }
        }

        // prefer the closest to the edge1 point
        int x1 = edge1.getX(0);
        int y1 = edge1.getY(0);
        PairInt closest = null;
        double closestDist = Integer.MAX_VALUE;
        for (PairInt edge2Loc : edge2LocForEdge1FirstLoc) {
            assert(edge2Loc.getX() == edge2Idx);
            int x2 = edge2.getX(edge2Loc.getY());
            int y2 = edge2.getY(edge2Loc.getY());
            int diffX = Math.abs(x1 - x2);
            int diffY = Math.abs(y1 - y2);
            double dist = Math.sqrt(diffX*diffX + diffY*diffY);
            if (dist < closestDist) {
                closestDist = dist;
                closest = edge2Loc;
            }
        }

        boolean didReverse = reverseTopIfPossible(edge2, closest.getY());

        if (didReverse) {
            edge2.reverse();
            edge1.insertAll(0, edge2);
            edges.remove(edge2Idx);
            return 1;
        } else {
            didReverse = reverseBottomIfPossible(edge2, closest.getY());
            if (didReverse) {
                edge2.addAll(edge1);
                edges.set(edge1Idx, edge2);
                edges.remove(edge2Idx);
                return 1;
            }
        }

        // if arrived here, try the points in set that are not same as 'closest'
        for (PairInt edge2Loc : edge2LocForEdge1FirstLoc) {

            if (edge2Loc.equals(closest)) {
                continue;
            }

            int idx = edge2Loc.getY();

            // same checks for whether can re-order to have an endpoint
            didReverse = reverseTopIfPossible(edge2, idx);

            if (didReverse) {
                edge2.reverse();
                edge1.insertAll(0, edge2);
                edges.remove(edge2Idx);
                return 1;
            } else {
                didReverse = reverseBottomIfPossible(edge2, idx);
                if (didReverse) {
                    edge2.addAll(edge1);
                    edges.set(edge1Idx, edge2);
                    edges.remove(edge2Idx);
                    return 1;
                }
            }
        }
        return 0;
    }

    private boolean find(List<PairIntArray> output, int x, int y) {

        for (int i = 0; i < output.size(); ++i) {

            PairIntArray p = output.get(i);

            for (int idx = 0; idx < p.getN(); ++idx) {
                if ((p.getX(idx) == x) && (p.getY(idx) == y)) {
                    return true;
                }
            }
        }
        return false;
    }

    /**
     * find a point starting from the last index that is adjacent to the first
     * index.  the method is used to learn whether there is a reversible section.
     * @param edge
     * @return
     */
    private int findAdjacentToTopAtBottom(PairIntArray edge) {

        int x1 = edge.getX(0);
        int y1 = edge.getY(0);

        for (int idx = (edge.getN() - 1); idx > 0; --idx) {

            int x2 = edge.getX(idx);
            int y2 = edge.getY(idx);

            int diffX = Math.abs(x1 - x2);
            int diffY = Math.abs(y1 - y2);

            if ((diffX < 2) && (diffY < 2)) {
                return idx;
            }
        }

        return -1;
    }

    private int findAdjacentToBottomAtTop(PairIntArray edge) {

        int xn = edge.getX(edge.getN() - 1);
        int yn = edge.getY(edge.getN() - 1);

        for (int idx = 0; idx < (edge.getN() - 2); ++idx) {

            int x2 = edge.getX(idx);
            int y2 = edge.getY(idx);

            int diffX = Math.abs(xn - x2);
            int diffY = Math.abs(yn - y2);

            if ((diffX < 2) && (diffY < 2)) {
                return idx;
            }
        }

        return -1;
    }

    private void replaceOpenCurvesIfPossibleToClose(List<PairIntArray> edges) {

        /*
        TODO: may move the code in insert() which handles
        reordering to join edge1ALoc to edge2ALoc and edge1BLoc to edge2BLoc
        to here.

        It's a means to untie knots too, but does not try to untie every knot.

        The knots can be found as junctions in an edge where the junction has
        more than 2 points from that edge in it.

        can use this structure to make lookups faster:
        Map<Integer, Set<Integer>> theEdgeToPixelIndexMap = createEdgeToPixelIndexMap();
        */

    }
}
