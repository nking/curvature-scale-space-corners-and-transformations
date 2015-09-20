package algorithms.imageProcessing;

import algorithms.LongestPath;
import algorithms.misc.MiscDebug;
import algorithms.misc.MiscMath;
import algorithms.util.PairInt;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
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
    public boolean applyFilter(Set<PairInt> closedCurve, Set<PairInt> exclude) {

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

            for (PairInt p : closedCurve) {
                
                if (exclude.contains(p)) {
                    continue;
                }

                Set<PairInt> junctionPoints = findJunction(p.getX(), p.getY(),
                    closedCurve);

                if (junctionPoints == null) {
                    continue;
                }

                assert(junctionPoints.size() == 3);

                Junction spokes = findSpokesOfJunction(junctionPoints,
                    closedCurve);

                if (spokes == null) {
                    continue;
                }

                List<Set<PairInt>> longestPaths = findLongestPaths(spokes,
                    closedCurve);

                assert(longestPaths.size() == 3);

                removeShortestPath(spokes, longestPaths, closedCurve);

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
            finder.findConnectedPointGroups(closedCurve, minMaxXY[1] + 1, 
                minMaxXY[3] + 1);
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

        if (patternPoints != null) {
            return patternPoints;
        }

        return null;
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

        Junction junction = new Junction();

        int count = 0;

        for (PairInt p : junctionPoints) {

            Set<PairInt> neighbors = curveHelper.findNeighbors(p.getX(),
                p.getY(), closedCurve);

            // remove the points that are other junction points
            for (PairInt jp : junctionPoints) {
                neighbors.remove(jp);
            }

            if (neighbors.size() > 2 || neighbors.isEmpty()) {
                return null;
            }

            if (count == 0) {
                junction.j0 = p;
                spokes = new ArrayList<Set<PairInt>>();
            } else if (count == 1) {
                junction.j1 = p;
            } else {
                junction.j2 = p;
            }
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
                        int diffX = Math.abs(pI.getX() - pJ.getX());
                        int diffY = Math.abs(pI.getY() - pJ.getY());
                        if ((diffX < 2) && (diffY < 2)) {
                            return null;
                        }
                    }
                }
                if (i == 0) {
                    int diffX = Math.abs(pI.getX() - junction.j1.getX());
                    int diffY = Math.abs(pI.getY() - junction.j1.getY());
                    if ((diffX < 2) && (diffY < 2)) {
                        return null;
                    }
                    diffX = Math.abs(pI.getX() - junction.j2.getX());
                    diffY = Math.abs(pI.getY() - junction.j2.getY());
                    if ((diffX < 2) && (diffY < 2)) {
                        return null;
                    }
                } else if (i == 1) {
                    int diffX = Math.abs(pI.getX() - junction.j0.getX());
                    int diffY = Math.abs(pI.getY() - junction.j0.getY());
                    if ((diffX < 2) && (diffY < 2)) {
                        return null;
                    }
                    diffX = Math.abs(pI.getX() - junction.j2.getX());
                    diffY = Math.abs(pI.getY() - junction.j2.getY());
                    if ((diffX < 2) && (diffY < 2)) {
                        return null;
                    }
                } else {
                    int diffX = Math.abs(pI.getX() - junction.j0.getX());
                    int diffY = Math.abs(pI.getY() - junction.j0.getY());
                    if ((diffX < 2) && (diffY < 2)) {
                        return null;
                    }
                    diffX = Math.abs(pI.getX() - junction.j1.getX());
                    diffY = Math.abs(pI.getY() - junction.j1.getY());
                    if ((diffX < 2) && (diffY < 2)) {
                        return null;
                    }
                }
            }
        }

        junction.s0 = spokes.get(0).iterator().next();
        junction.s1 = spokes.get(1).iterator().next();
        junction.s2 = spokes.get(2).iterator().next();

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
        rm.remove(spokes.j0);
        rm.remove(spokes.j1);
        rm.remove(spokes.j2);

        closedCurve.removeAll(rm);
    }

    private List<Set<PairInt>> findLongestPaths(Junction spokes,
        Set<PairInt> closedCurve) {

        /*
        each path source is spokes.si and dest is one of spokes.jj or a point
        on the path for that spoke.
        */

        List<Set<PairInt>> paths = new ArrayList<Set<PairInt>>();

        LongestPath longestPath = new LongestPath();

        Set<PairInt> dest = new HashSet<PairInt>();
        dest.add(spokes.j1);
        dest.add(spokes.j2);
        Set<PairInt> path0 = longestPath.findMaxCostPath(closedCurve,
            spokes.s0, dest);
        paths.add(path0);

        dest = new HashSet<PairInt>();
        dest.add(spokes.j0);
        dest.add(spokes.j2);
        Set<PairInt> path1 = longestPath.findMaxCostPath(closedCurve,
            spokes.s1, dest);
        paths.add(path1);

        dest = new HashSet<PairInt>();
        dest.add(spokes.j0);
        dest.add(spokes.j1);
        Set<PairInt> path2 = longestPath.findMaxCostPath(closedCurve,
            spokes.s2, dest);
        paths.add(path2);

        return paths;
    }

    public static class Pattern {
        Set<PairInt> ones;
        Set<PairInt> zeroes;
    }

    public static class Junction {
        PairInt j0 = null;
        PairInt j1 = null;
        PairInt j2 = null;
        PairInt s0 = null;
        PairInt s1 = null;
        PairInt s2 = null;
    }
}
