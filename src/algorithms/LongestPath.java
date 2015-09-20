package algorithms;

import algorithms.misc.Misc;
import algorithms.util.PairInt;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.Stack;
import java.util.logging.Logger;

/**
 * A version of a path builder that prefers the longest path between src
 * and either destination node or back to itself.
 * The algorithm is not currently optimal and is meant for use to quickly
 * compare knots in closed curves.  It may be improved over time as needed.
 *
 *
 * @author nichole
 */
public class LongestPath {

    protected Map<PairInt, Double> dist = null;
    protected Map<PairInt, PairInt> prev = null;
    protected Map<PairInt, Double> heur = null;
    protected Set<PairInt> unVisited = null;

    protected Set<PairInt> dest = null;

    protected final static PairInt sentinel = new PairInt(-1, -1);

    public LongestPath() {
        dist = new HashMap<PairInt, Double>();
        prev = new HashMap<PairInt, PairInt>();
        heur = new HashMap<PairInt, Double>();
        unVisited = new HashSet<PairInt>();

        dest = null;
    }

    /**
     * find the path that maximizes the cost from the source to one of the
     * given destinations or the longest path if those destinations are not
     * in the path (one use case is for a portion of a curve which may double
     * back on itself and another is for a portion of a curve which continues
     * from one side to the other and reaches the destination).
     *
     * runtime is O(|V||E|).
     *
     * @return
     */
    public Set<PairInt> findMaxCostPath(Set<PairInt> points, PairInt src,
        Set<PairInt> potentialDest) {

        dest = potentialDest;

        /*
        goal is to visit nodes from src to dest or until nodes run out to
        find the longest path for that segment.

        a DFS style visit pattern is applied with modifications.

        the maximum distance of v among all of u's neighbors is found where
           maximum distance is the largest distance from u to v,
        else if equal is the max distance from v to the src,
        else if equal is the larger of the heuristic of potential v's if any.

        this method invokes search for the traversal.
        if search returns that the number of adjacent visited nodes
        is less than all nodes expected in path, this method
        sets the heuristic to the current dist map, initializes the
        data structures, and tries again.
        this has the effect of choosing the other path for equally distant
        neighbors and can help to correct crossed paths if there is only
        one in the path.

        for more complex curves, one could consider storing such choices
        and trying all combinations of them and keeping the first that
        spans from dest to src or adjacent to src...improving the heuristic
        with each combination...
        */

        PairInt nVisited = search(points, src, potentialDest);

        if (nVisited.getY() > 1) {

            heur.putAll(dist);
            dist = new HashMap<PairInt, Double>();
            prev = new HashMap<PairInt, PairInt>();
            unVisited = new HashSet<PairInt>();

            nVisited = search(points, src, potentialDest);
        }

        PairInt lastPoint = findMaxDist();
        boolean lastIsSrcOrDest = lastPoint.equals(src);
        if (!lastIsSrcOrDest) {
            for (PairInt p : dest) {
                if (lastPoint.equals(p)) {
                    lastIsSrcOrDest = true;
                    break;
                }
            }
            if (!lastIsSrcOrDest) {
                lastIsSrcOrDest = areAdjacent(lastPoint, src);
            }
        }
        if (!lastIsSrcOrDest) {
            //TODO: improve...
            Logger.getLogger(this.getClass().getName()).fine(
            "WARNING: neither src nor dest nodes reached");
        }

        Set<PairInt> pathPoints = new HashSet<PairInt>();
        pathPoints.addAll(dist.keySet());
        pathPoints.addAll(unVisited);

        return pathPoints;
    }

    private boolean areAdjacent(PairInt p0, PairInt p1) {
        int diffX = Math.abs(p0.getX() - p1.getX());
        int diffY = Math.abs(p0.getY() - p1.getY());
        if ((diffX < 2) && (diffY < 2)) {
            return true;
        }
        return false;
    }

    /**
     * find the path that maximizes the cost from the source to one of the
     * given destinations or the longest path if those destinations are not
     * in the path (one use case is for a portion of a curve which may double
     * back on itself and another is for a portion of a curve which continues
     * from one side to the other and reaches the destination).
     * runtime is O(|V||E|).
     *
     * @return
     */
    private PairInt search(Set<PairInt> points, PairInt src,
        Set<PairInt> potentialDest) {

        dest = potentialDest;

        Set<PairInt> visited = new HashSet<PairInt>();

        // source index is presumed to be index=0
        dist.put(src, Double.valueOf(0));

        Stack<PairInt> stack = new Stack<PairInt>();
        stack.add(src);

        //visited.add(src);

        int nIter = 0;

        int[] dxs8 = Misc.dx8;
        int[] dys8 = Misc.dy8;

        while (!stack.isEmpty()) {

            PairInt uPoint = stack.pop();

            double maxDist = Double.MIN_VALUE;
            PairInt maxDistVPoint = null;

            for (int vIdx = 0; vIdx < dxs8.length; ++vIdx) {

                PairInt vPoint = new PairInt(uPoint.getX() + dxs8[vIdx],
                    uPoint.getY() + dys8[vIdx]);

                if (visited.contains(vPoint) || !points.contains(vPoint)) {
                    continue;
                }

                // prefer the more distant of neighbors:
                int diffX = uPoint.getX() - vPoint.getX();
                int diffY = uPoint.getY() - vPoint.getY();
                double uvDist = Math.sqrt(diffX*diffX + diffY*diffY);

                unVisited.add(vPoint);

                if (uvDist > maxDist) {
                    maxDist = uvDist;
                    maxDistVPoint = vPoint;
                } else if (uvDist == maxDist) {
                    // choose which is more distant from src
                    int diffX0 = src.getX() - maxDistVPoint.getX();
                    int diffY0 = src.getY() - maxDistVPoint.getY();
                    double dist0 = Math.sqrt(diffX0*diffX0 + diffY0*diffY0);

                    int diffX1 = src.getX() - vPoint.getX();
                    int diffY1 = src.getY() - vPoint.getY();
                    double dist1 = Math.sqrt(diffX1*diffX1 + diffY1*diffY1);

                    if (dist1 > dist0) {
                        maxDist = uvDist;
                        maxDistVPoint = vPoint;
                    } else if (dist1 == dist0) {
                        Double vHeur = heur.get(vPoint);
                        Double maxDistVHeur = heur.get(maxDistVPoint);
                        if (vHeur == null || maxDistVHeur == null) {
                            continue;
                        }
                        if (vHeur.doubleValue() > maxDistVHeur.doubleValue()) {
                            maxDist = uvDist;
                            maxDistVPoint = vPoint;
                        }
                    }
                }
            }

            if (maxDistVPoint != null) {

                Double distU = dist.get(uPoint);

                // break if have reached destination or are adj to src
                //TODO: consider whether this can use distU > 1
                boolean reached = false;
                if (distU > 2) {
                    for (PairInt p : dest) {
                        if (maxDistVPoint.equals(p)) {
                            reached = true;
                            break;
                        }
                    }
                    if (reached) {
                        break;
                    }
                    if (areAdjacent(maxDistVPoint, src)) {
                        break;
                    }
                }

                assert(distU != null);
                double uPlusUV = distU.doubleValue() + maxDist;
                dist.put(maxDistVPoint, Double.valueOf(uPlusUV));
                prev.put(maxDistVPoint, uPoint);
                stack.add(maxDistVPoint);

                visited.add(uPoint);
                visited.add(maxDistVPoint);
                unVisited.remove(maxDistVPoint);

            } else {

                break;
            }

            nIter++;
        }

        PairInt nVisited = new PairInt(visited.size(), unVisited.size());

        return nVisited;
    }

    private PairInt findMaxDist() {

        PairInt p = null;

        double maxDist = Double.MIN_VALUE;

        for (Entry<PairInt, Double> entry : dist.entrySet()) {

            Double d = entry.getValue();

            if (d.doubleValue() > maxDist) {
                maxDist = d.doubleValue();
                p = entry.getKey();
            }
        }

        return p;
    }

}
