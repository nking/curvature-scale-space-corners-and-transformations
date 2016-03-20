package algorithms.imageProcessing;

import algorithms.misc.Misc;
import algorithms.misc.MiscDebug;
import algorithms.misc.MiscMath;
import algorithms.shortestPath.AStar;
import algorithms.util.PairInt;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 *
 * @author nichole
 */
class UntraversableLobeRemover {

    /**
     * NOTE: very important that spurs be removed before using this because
     * a short spur found as a junction can lead to a false finding of the
     * shortest segment as a segment that isn't the spur.
     *
     * @param closedCurve
     * @return true if points were removed
     */
    public boolean applyFilter(Set<PairInt> closedCurve, Set<PairInt> exclude,
        int imageWidth, int imageHeight) {
        
        /*
        Example:
                         -3
                 #       -2
              -  #       -1    The center pattern can be found, then if each
           #  #  -  -     0    of the 3 has a neighbor that isn't adjacent
              -  #        1    to the others and those spokes are single pixel,
                    #          widths, it is a junction that cannot be traversed
                               one way and then the other.
       -2 -1  0  1  2

        To further locate the untraversable section, can find the longest path
        from each spoke back to a pixel adjacent to itself and
        the shortest of the longest paths can be trimmed.

        For example:
                        #         -3
             #          #         -2
          #     #     S2          -1
          #     #  S1              0
             #        S3           1
                         #         2
                            #      3
         -3 -2 -1  0  1  2  3

        Of the longest paths for S1, S2, and S3, that of S1 would be the
        shortest so could be trimmed.
        */

        int nChanges = 0;
        int nIter = 0;
        int nIterMax = 10;

        boolean wasChanged = false;

        while ((nIter == 0) || ((nChanges > 0) && (nIter < nIterMax))) {

            nChanges = 0;

if (false) {
Image img1 = new Image(imageWidth, imageHeight);
for (PairInt p : closedCurve) {
    img1.setRGB(p.getX(), p.getY(), 255, 0, 0);
}
MiscDebug.writeImage(img1, "before_removed_untr_lobe_"
    + MiscDebug.getCurrentTimeFormatted());
}
            for (PairInt p : closedCurve) {

                if (exclude.contains(p)) {
                    continue;
                }

                Set<PairInt> junctionPoints = findJunction(p.getX(), p.getY(),
                    closedCurve);

                Set<PairInt> diagJunctionPoints = null;

                boolean isDefaultJunction = true;
                
                if (junctionPoints == null) {

                    diagJunctionPoints = findDiagJunction(p.getX(), p.getY(),
                        closedCurve);

                    if (diagJunctionPoints == null) {
                        continue;
                    }
                    
                    isDefaultJunction = false;
                }

                assert((isDefaultJunction && junctionPoints.size() == 3) ||
                    (!isDefaultJunction && diagJunctionPoints.size() == 4));

                Junction spokes = null;

                if (isDefaultJunction) {
                    spokes = findSpokesOfJunction(junctionPoints,
                        closedCurve);
                } else {
                    spokes = findSpokesOfJunction2(diagJunctionPoints,
                        closedCurve);
                }

                if (spokes == null) {
                    continue;
                }

                List<Set<PairInt>> longestPaths = findShortestPathBetweenSpokes(
                    spokes, closedCurve, isDefaultJunction);

                //assert(longestPaths.size() == 3);

                removeShortestPath(spokes, longestPaths, closedCurve);
if (false) {
Image img3 = new Image(imageWidth, imageHeight);
for (PairInt p2 : closedCurve) {
    img3.setRGB(p2.getX(), p2.getY(), 255, 0, 0);
}
MiscDebug.writeImage(img3, "removed_untr_lobe_"
+ MiscDebug.getCurrentTimeFormatted());
}
                assert(closedCurve.size() > 0);

                SpurRemover spurRm = new SpurRemover();
                spurRm.remove(closedCurve, imageWidth, imageHeight);

                assert(closedCurve.size() > 0);

                nChanges++;

                wasChanged = true;

                break;
            }

            nIter++;
        }

        // if any have changed, use contiguous pix finder to keep only
        // the contiguous
        if (wasChanged) {

            int[] minMaxXY = MiscMath.findMinMaxXY(closedCurve);
            DFSConnectedGroupsFinder finder = new DFSConnectedGroupsFinder();
            finder.findConnectedPointGroups(closedCurve);
            int nMax = Integer.MIN_VALUE;
            int maxIdx = -1;
            for (int i = 0; i < finder.getNumberOfGroups(); ++i) {
                int n = finder.getNumberofGroupMembers(i);
                if (n > nMax) {
                    nMax = n;
                    maxIdx = i;
                }
            }
            if (maxIdx > -1) {
                closedCurve.clear();
                closedCurve.addAll(finder.getXY(maxIdx));

                SpurRemover spurRm = new SpurRemover();
                spurRm.remove(closedCurve, minMaxXY[1] + 1, minMaxXY[3] + 1);
            }
        }

        return wasChanged;
    }

    private Set<PairInt> findJunction(int x, int y, Set<PairInt> points) {

        Pattern pattern = getHorizPattern();

        Set<PairInt> patternPoints = findPattern(x, y, points, pattern);

        if (patternPoints != null) {
            return patternPoints;
        }

        swapXDirection(pattern);

        patternPoints = findPattern(x, y, points, pattern);

        if (patternPoints != null) {
            return patternPoints;
        }

        pattern = getVertPattern();

        patternPoints = findPattern(x, y, points, pattern);

        if (patternPoints != null) {
            return patternPoints;
        }

        swapYDirection(pattern);

        patternPoints = findPattern(x, y, points, pattern);

        return patternPoints;
    }

    private Set<PairInt> findDiagJunction(int x, int y, Set<PairInt> points) {

        Pattern pattern = getDiagPattern();

        Set<PairInt> patternPoints = findPattern(x, y, points, pattern);

        if (patternPoints != null) {
            return patternPoints;
        }

        swapXDirection(pattern);

        patternPoints = findPattern(x, y, points, pattern);

        if (patternPoints != null) {
            return patternPoints;
        }

        swapYDirection(pattern);

        patternPoints = findPattern(x, y, points, pattern);

        if (patternPoints != null) {
            return patternPoints;
        }

        swapXDirection(pattern);

        patternPoints = findPattern(x, y, points, pattern);

        if (patternPoints != null) {
            return patternPoints;
        }

        // ----- diag pattern 2 ---
        pattern = getDiagPattern2();

        patternPoints = findPattern(x, y, points, pattern);

        if (patternPoints != null) {
            return patternPoints;
        }

        swapXDirection(pattern);

        patternPoints = findPattern(x, y, points, pattern);

        if (patternPoints != null) {
            return patternPoints;
        }

        pattern = getDiagPattern3();

        patternPoints = findPattern(x, y, points, pattern);

        if (patternPoints != null) {
            return patternPoints;
        }

        swapYDirection(pattern);

        patternPoints = findPattern(x, y, points, pattern);

        if (patternPoints != null) {
            return patternPoints;
        }

        return patternPoints;
    }

    protected Pattern getHorizPattern() {

        /*
                         -3
                 .       -2
              -  #       -1
           .  #  -  -     0
              -  #        1
                    .

       -2 -1  0  1  2
        */

        Pattern pr = new Pattern();
        pr.ones = new HashSet<PairInt>();
        pr.zeroes = new HashSet<PairInt>();

        pr.zeroes.add(new PairInt(0, 1)); pr.zeroes.add(new PairInt(0, -1));
        pr.zeroes.add(new PairInt(1, 0));
        pr.zeroes.add(new PairInt(2, 0));

        pr.ones.add(new PairInt(1, 1)); pr.ones.add(new PairInt(1, -1));

        return pr;
    }

    protected Pattern getVertPattern() {

        /*
              .          -1
           -  #  -        0
           #  -  #        1
              -           2
       -2 -1  0  1  2
        */
        Pattern pr = new Pattern();
        pr.ones = new HashSet<PairInt>();
        pr.zeroes = new HashSet<PairInt>();

        pr.zeroes.add(new PairInt(-1, 0));
        pr.zeroes.add(new PairInt(0, 2)); pr.zeroes.add(new PairInt(0, 1));
        pr.zeroes.add(new PairInt(1, 0));

        pr.ones.add(new PairInt(-1, 1)); pr.ones.add(new PairInt(1, 1));

        return pr;
    }

    protected DiagPattern getDiagPattern() {

        /*
                         -3
              -  #  -    -2
              -  #  #    -1
              #  -  -     0
                          1
       -2 -1  0  1  2
        */

        DiagPattern pr = new DiagPattern();
        pr.ones = new HashSet<PairInt>();
        pr.zeroes = new HashSet<PairInt>();

        pr.zeroes.add(new PairInt(0, -1)); pr.zeroes.add(new PairInt(0, -2));
        pr.zeroes.add(new PairInt(1, 0));
        pr.zeroes.add(new PairInt(2, 0)); pr.zeroes.add(new PairInt(2, -2));

        pr.ones.add(new PairInt(1, -1)); pr.ones.add(new PairInt(1, -2));
        pr.ones.add(new PairInt(2, -1));

        return pr;
    }

    protected DiagPattern getDiagPattern2() {

        /*
                         -3
              -  #  -    -2
              -  #  #    -1
              -  #  -     0
                          1
          -2 -1  0  1  2
        */

        DiagPattern pr = new DiagPattern();
        pr.ones = new HashSet<PairInt>();
        pr.zeroes = new HashSet<PairInt>();

        pr.zeroes.add(new PairInt(-1, 0)); pr.zeroes.add(new PairInt(-1, -1)); pr.zeroes.add(new PairInt(-1, -2));
        pr.zeroes.add(new PairInt(1, 0)); pr.zeroes.add(new PairInt(1, -2));

        pr.ones.add(new PairInt(0, -1)); pr.ones.add(new PairInt(0, -2));
        pr.ones.add(new PairInt(1, -1));

        return pr;
    }

    protected DiagPattern getDiagPattern3() {

        /*
                         -2
              -  #  -    -1
              #  #  #     0
              -  -  -     1

          -1  0  1  2
        */

        DiagPattern pr = new DiagPattern();
        pr.ones = new HashSet<PairInt>();
        pr.zeroes = new HashSet<PairInt>();

        pr.zeroes.add(new PairInt(0, 1)); pr.zeroes.add(new PairInt(0, -1));
        pr.zeroes.add(new PairInt(1, 1));
        pr.zeroes.add(new PairInt(2, 1)); pr.zeroes.add(new PairInt(2, -1));

        pr.ones.add(new PairInt(1, 0)); pr.ones.add(new PairInt(1, -1));
        pr.ones.add(new PairInt(2, 0));

        return pr;
    }

    private void swapYDirection(Pattern pattern) {
        // ----- change the sign of y  -----
        for (PairInt p : pattern.zeroes) {
            p.setY(-1 * p.getY());
        }
        for (PairInt p : pattern.ones) {
            p.setY(-1 * p.getY());
        }
    }

    private void swapXDirection(Pattern pattern) {
        // ----- change the sign of x  -----
        for (PairInt p : pattern.zeroes) {
            p.setX(-1 * p.getX());
        }
        for (PairInt p : pattern.ones) {
            p.setX(-1 * p.getX());
        }
    }

    private Set<PairInt> findPattern(int x, int y, Set<PairInt> points,
        Pattern pattern) {

        Set<PairInt> patternPoints = new HashSet<PairInt>();

        for (PairInt p : pattern.zeroes) {
            PairInt p2 = new PairInt(x + p.getX(), y + p.getY());
            if (points.contains(p2)) {
                return null;
            }
        }
        for (PairInt p : pattern.ones) {
            PairInt p2 = new PairInt(x + p.getX(), y + p.getY());
            if (!points.contains(p2)) {
                return null;
            }
            patternPoints.add(p2);
        }

        patternPoints.add(new PairInt(x, y));

        return patternPoints;
    }

    private Junction findSpokesOfJunction(Set<PairInt> junctionPoints,
        Set<PairInt> closedCurve) {

        assert(junctionPoints.size() == 3);

        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();

        /*
        Example:
                 S       -2
              -  #       -1    The center pattern can be found, then if each
           S  #  -  -     0    of the 3 has a neighbor that isn't adjacent
              -  #        1    to the others and those spokes are single pixel,
                    S          widths, it is a junction that cannot be traversed
                               one way and then the other.
       -2 -1  0  1  2
        */

        List<Set<PairInt>> spokes = null;

        Junction junction = new Junction(3);

        int count = 0;

        for (PairInt p : junctionPoints) {

            Set<PairInt> neighbors = curveHelper.findNeighbors(p.getX(),
                p.getY(), closedCurve);

            neighbors.removeAll(junctionPoints);

            if (neighbors.size() > 2 || neighbors.isEmpty()) {
                return null;
            }

            if (count == 0) {
                spokes = new ArrayList<Set<PairInt>>();
            }
            junction.setJ(count, p);
            Set<PairInt> spoke = new HashSet<PairInt>();
            spoke.addAll(neighbors);
            spokes.add(spoke);

            count++;
        }

        assert(spokes.size() == 3);

        // assert that the points in one spoke are not adjacent to others
        // and not adjacent to the oppossing junction points
        for (int i = 0; i < spokes.size(); ++i) {
            Set<PairInt> spokeI = spokes.get(i);
            assert(!spokeI.isEmpty());
            for (PairInt pI : spokeI) {
                for (int j = (i + 1); j < spokes.size(); ++j) {
                    Set<PairInt> spokeJ = spokes.get(j);
                    assert(!spokeJ.isEmpty());
                    for (PairInt pJ : spokeJ) {
                        if (areAdjacent(pI, pJ)) {
                            return null;
                        }
                    }
                }
                if (i == 0) {
                    if (areAdjacent(pI, junction.getJ(1))) {
                        return null;
                    }
                    if (areAdjacent(pI, junction.getJ(2))) {
                        return null;
                    }
                } else if (i == 1) {
                    if (areAdjacent(pI, junction.getJ(0))) {
                        return null;
                    }
                    if (areAdjacent(pI, junction.getJ(2))) {
                        return null;
                    }
                } else {
                    if (areAdjacent(pI, junction.getJ(0))) {
                        return null;
                    }
                    if (areAdjacent(pI, junction.getJ(1))) {
                        return null;
                    }
                }
            }
        }

        junction.setS0(spokes.get(0).iterator().next());
        junction.setS1(spokes.get(1).iterator().next());
        junction.setS2(spokes.get(2).iterator().next());

        return junction;
    }

    /**
     * find spokes for the diagonal patterns
     * @param junctionPoints
     * @param closedCurve
     * @return
     */
    private Junction findSpokesOfJunction2(Set<PairInt> junctionPoints,
        Set<PairInt> closedCurve) {

        assert(junctionPoints.size() == 4);

        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();

        /*
        Example:
              S          -2
              #          -1    The center pattern can be found, then if each
              #  #  S     0    of the 3 has a neighbor that isn't adjacent
              #           1    to the others and those spokes are single pixel,
              S                widths, it is a junction that cannot be traversed
                               one way and then the other.
       -2 -1  0  1  2
        */

        List<Set<PairInt>> spokes = null;

        Junction junction = new Junction(4);

        int count = 0;

        // find the center junction by centroid
        double[] xyCen = curveHelper.calculateXYCentroids(junctionPoints);
        PairInt jC = new PairInt((int)Math.round(xyCen[0]), (int)Math.round(xyCen[1]));
        assert(junctionPoints.contains(jC));

        for (PairInt p : junctionPoints) {

            Set<PairInt> neighbors = curveHelper.findNeighbors(p.getX(),
                p.getY(), closedCurve);

            neighbors.removeAll(junctionPoints);

            if (!p.equals(jC)) {
                if (neighbors.isEmpty() || neighbors.size() > 2) {
                    return null;
                }
                if (neighbors.size() > 1) {
                    // keep the closest to jC
                    int minDistSq = Integer.MAX_VALUE;
                    PairInt minDistP = null;
                    for (PairInt pN : neighbors) {
                        int distSq = distSq(jC, pN);
                        if (distSq < minDistSq) {
                            minDistSq = distSq;
                            minDistP = pN;
                        }
                    }
                    neighbors.clear();
                    neighbors.add(minDistP);
                }
            }

            if (count == 0) {
                spokes = new ArrayList<Set<PairInt>>();
            }
            junction.setJ(count, p);
            if (!p.equals(jC)) {
                Set<PairInt> spoke = new HashSet<PairInt>();
                spoke.addAll(neighbors);
                spokes.add(spoke);
            }

            count++;
        }

        assert(spokes.size() == 3);

        // assert that the points in one spoke are not adjacent to others
        // and not adjacent to the oppossing junction points
        for (int i = 0; i < spokes.size(); ++i) {
            Set<PairInt> spokeI = spokes.get(i);
            assert(!spokeI.isEmpty());
            for (PairInt pI : spokeI) {
                for (int j = (i + 1); j < spokes.size(); ++j) {
                    Set<PairInt> spokeJ = spokes.get(j);
                    assert(!spokeJ.isEmpty());
                    for (PairInt pJ : spokeJ) {
                        if (areAdjacent(pI, pJ)) {
                            return null;
                        }
                    }
                }
                // it will be adjacent to 1 junction point
                int nAdj = 0;
                for (int k = 0; k < junction.getJLength(); ++k) {
                    if (areAdjacent(pI, junction.getJ(k))) {
                        nAdj++;
                    }
                }
                if (nAdj > 1) {
                    return null;
                }
            }
        }

        junction.setS0(spokes.get(0).iterator().next());
        junction.setS1(spokes.get(1).iterator().next());
        junction.setS2(spokes.get(2).iterator().next());

        return junction;
    }

    private void removeShortestPath(Junction spokes, List<Set<PairInt>> longestPaths,
        Set<PairInt> closedCurve) {

        int minLen = Integer.MAX_VALUE;
        int minIdx = -1;
        for (int i = 0; i < longestPaths.size(); ++i) {
            int n = longestPaths.get(i).size();
            if (n < minLen) {
                minLen = n;
                minIdx = i;
            }
        }
        assert(minIdx != -1);
        Set<PairInt> rm = longestPaths.get(minIdx);

        // make sure the junction points don't get removed
        for (int j = 0; j < spokes.getJLength(); ++j) {
            rm.remove(spokes.getJ(j));
        }

        closedCurve.removeAll(rm);
    }

    /**
     * find the shortest paths between spokes, excluding the junction points.
     * runtime complexity is O(|V|) + O(|E| + |V|lg2|V|).
     * @param spokes
     * @param closedCurve
     * @param isDefaultJunction if horiz or vert patterns were found, this
     * is the default junction pattern, else is one of the diag patterns
     * @return
     */
    private List<Set<PairInt>> findShortestPathBetweenSpokes(Junction spokes,
        Set<PairInt> closedCurve, boolean isDefaultJunction) {

        int n = closedCurve.size();

        /*
        A* needs
        -- point list as an array (excluding junctions)
        -- adjacency list as a parallel array of linked lists holding indexes
        -- heuristics (or let the algorithm determine the line of sight from
           a point to the destination.  it's only needed for the points which
           are visited).
        -- the adjacency lists should not contain the junction points

        For local bookkeeping, need to keep track of where the spokes
        variables are and need to keep track of where the points are in the
        points array in order to be able to make the adjacency list.
        */

        Map<PairInt, Integer> coordsIndexMap = new HashMap<PairInt, Integer>();

        PairInt[] points = isDefaultJunction ?
            new PairInt[n - 3] : new PairInt[n - 4];

        int s0Idx = -1;
        int s1Idx = -1;
        int s2Idx = -1;

        int count = 0;
        for (PairInt p : closedCurve) {
            boolean skip = false;
            for (int j = 0; j < spokes.getJLength(); ++j) {
                if (p.equals(spokes.getJ(j))) {
                    skip = true;
                    break;
                }
            }
            if (skip) {
                continue;
            }
            if (p.equals(spokes.getS0())) {
                s0Idx = count;
            } else if (p.equals(spokes.getS1())) {
                s1Idx = count;
            } else if (p.equals(spokes.getS2())) {
                s2Idx = count;
            }
            points[count] = p;

            coordsIndexMap.put(p, Integer.valueOf(count));

            count++;
        }

        assert(s0Idx > -1);
        assert(s1Idx > -1);
        assert(s2Idx > -1);
        assert((isDefaultJunction && (count == (n - 3))) ||
            (!isDefaultJunction && (count == (n - 4))));

        int[] dxs8 = Misc.dx8;
        int[] dys8 = Misc.dy8;

        List<LinkedList<Integer>> adjList = new ArrayList<LinkedList<Integer>>();

        for (int i = 0; i < points.length; ++i) {

            PairInt p = points[i];

            LinkedList<Integer> list = new LinkedList<Integer>();

            for (int j = 0; j < dxs8.length; ++j) {
                PairInt p2 = new PairInt(p.getX() + dxs8[j], p.getY() + dys8[j]);
                boolean skip = false;
                for (int k = 0; k < spokes.getJLength(); ++k) {
                    if (p2.equals(spokes.getJ(k))) {
                        skip = true;
                        break;
                    }
                }
                if (skip) {
                    continue;
                }
                Integer index2 = coordsIndexMap.get(p2);
                if (index2 == null) {
                    continue;
                }
                list.add(index2);
            }
            adjList.add(list);
        }

        /*
        shortest path from
           src0Idx to src1Idx
           src0Idx to src2Idx
           src1Idx to src2Idx
        */ 
        int[] path01;
        int[] path02;
        int[] path12;
        
        /*
        We want to find the smallest loop which only contains one of the
        spikes (or two in some case) to trim that leaving the longer
        traversable path.
        A combination of using AStar and reading the distance matrix
        finds this lobe.
        */
        
        AStar aStar = new AStar(points, adjList, s0Idx, s1Idx);
        path01 = aStar.search();
        int[] path01Alt = aStar.createSourceToSourcePath();
        if (path01 == null) {
            path01 = path01Alt;
        } else if ((path01Alt != null) && (path01.length > path01Alt.length)) {
            path01 = path01Alt;
        }
    
        aStar = new AStar(points, adjList, s0Idx, s2Idx);
        path02 = aStar.search();
        int[] path02Alt = aStar.createSourceToSourcePath();
        if (path02 == null) {
            path02 = path02Alt;
        } else if ((path02Alt != null) && (path02.length > path02Alt.length)) {
            path02 = path02Alt;
        }
     
        aStar = new AStar(points, adjList, s1Idx, s2Idx);
        path12 = aStar.search();
        int[] path12Alt = aStar.createSourceToSourcePath();
        if (path12 == null) {
            path12 = path12Alt;
        } else if ((path12Alt != null) && (path12.length > path12Alt.length)) {
            path12 = path12Alt;
        }
        boolean allSameLength = false;
        if ((path01 != null) && (path02 != null) && (path12 != null)
            && (path01.length == path02.length) && (path01.length == path12.length)) {
            allSameLength = true;
        } else if ((path01 != null) && (path02 != null)
            && (path01.length == path02.length)) {
            allSameLength = true;
        } else if ((path01 != null) && (path12 != null)
            && (path01.length == path12.length)) {
            allSameLength = true;
        } else if ((path02 != null) && (path12 != null)
            && (path02.length == path12.length)) {
            allSameLength = true;
        }
        if (allSameLength) {
            aStar = new AStar(points, adjList, s2Idx, s0Idx);
            path12 = aStar.search();
            path12Alt = aStar.createSourceToSourcePath();
            if (path12 == null) {
                path12 = path12Alt;
            } else if ((path12Alt != null) && (path12.length > path12Alt.length)) {
                path12 = path12Alt;
            }
        }
       
        List<Set<PairInt>> paths = new ArrayList<Set<PairInt>>();
        if (path01 != null) {
            Set<PairInt> set = new HashSet<PairInt>();
            for (int i = 0; i < path01.length; ++i) {
                PairInt p = points[path01[i]];
                if (p.equals(spokes.getS0()) || p.equals(spokes.getS1())) {
                    continue;
                }
                set.add(p);
            }
            paths.add(set);
        }

        if (path02 != null) {
            Set<PairInt> set = new HashSet<PairInt>();
            for (int i = 0; i < path02.length; ++i) {
                PairInt p = points[path02[i]];
                if (p.equals(spokes.getS0()) || p.equals(spokes.getS2())) {
                    continue;
                }
                set.add(p);
            }
            paths.add(set);
        }

        if (path12 != null) {
            Set<PairInt> set = new HashSet<PairInt>();
            for (int i = 0; i < path12.length; ++i) {
                PairInt p = points[path12[i]];
                if (p.equals(spokes.getS1()) || p.equals(spokes.getS2())) {
                    continue;
                }
                set.add(p);
            }
            paths.add(set);
        }

        return paths;
    }

    private boolean areAdjacent(PairInt p0, PairInt p1) {
        int diffX = Math.abs(p0.getX() - p1.getX());
        int diffY = Math.abs(p0.getY() - p1.getY());
        if ((diffX < 2) && (diffY < 2)) {
            return true;
        }
        return false;
    }
    
    private int distSq(PairInt p0, PairInt p1) {
        int diffX = Math.abs(p0.getX() - p1.getX());
        int diffY = Math.abs(p0.getY() - p1.getY());
        return (diffX*diffX) + (diffY*diffY);
    }

    public static class Pattern {
        Set<PairInt> ones;
        Set<PairInt> zeroes;
    }

    public static class DiagPattern extends Pattern {
    }

    public static class Junction {
        final PairInt[] s = new PairInt[3];
        final PairInt[] j;
        public Junction(int nJunctions) {
            j = new PairInt[nJunctions];
        }
        public void setS0(PairInt p) {
            s[0] = p;
        }
        public void setS1(PairInt p) {
            s[1] = p;
        }
        public void setS2(PairInt p) {
            s[2] = p;
        }
        public PairInt getS0() {
            return s[0];
        }
        public PairInt getS1() {
            return s[1];
        }
        public PairInt getS2() {
            return s[2];
        }
        public void setJ(int index, PairInt p) {
            if (index < 0 || (index > (j.length - 1))) {
                throw new IllegalArgumentException("index is out of bounds");
            }
            j[index] = p;
        }
        public PairInt getJ(int index) {
            if (index < 0 || (index > (j.length - 1))) {
                throw new IllegalArgumentException("index is out of bounds");
            }
            return j[index];
        }
        public int getJLength() {
            return j.length;
        }
    }
}
