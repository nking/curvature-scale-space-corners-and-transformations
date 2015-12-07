package algorithms.imageProcessing;

import algorithms.misc.Misc;
import algorithms.util.PairInt;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

/**
 * a class to aggregate the results of the Hough transform within a given
 * tolerance in theta and radius for contiguous pixels.
 *
 * @author nichole
 */
public class DFSSimilarThetaRadiusGroupsFinder {

    public DFSSimilarThetaRadiusGroupsFinder() {

    }
    
    /**
     * aggregate the results of the Hough transform within a given tolerance in
     * theta and radius for contiguous pixels.
     */
    public Map<PairInt, PairInt> findConnectedPointGroups(List<PairInt> sortedThetaRadiusKeys,
        Map<PairInt, Set<PairInt>> thetaRadiusPixMap, int thetaTol,
        int radiusTol) {
        
        Map<PairInt, PairInt> pixToTRMap = findConnectedPointGroupsIterative(
            sortedThetaRadiusKeys, thetaRadiusPixMap, thetaTol, radiusTol);
        
        return pixToTRMap;
        
    }

    /**
     * aggregate the results of the Hough transform within a given tolerance in
     * theta and radius for contiguous pixels.
     */
    private Map<PairInt, PairInt> findConnectedPointGroupsIterative(List<PairInt> sortedThetaRadiusKeys,
        Map<PairInt, Set<PairInt>> thetaRadiusPixMap, int thetaTol,
        int radiusTol) {

        Map<PairInt, PairInt> pixToTRMap = new HashMap<PairInt, PairInt>();

        if (sortedThetaRadiusKeys.isEmpty()) {
            return pixToTRMap;
        }

        Map<PairInt, PairInt> pixToTRMapNotAggr = new HashMap<PairInt, PairInt>();
        for (Entry<PairInt, Set<PairInt>> entry : thetaRadiusPixMap.entrySet()) {
            for (PairInt p : entry.getValue()) {
                pixToTRMapNotAggr.put(p, entry.getKey());
            }
        }

        int[] dxs8 = Misc.dx8;
        int[] dys8 = Misc.dy8;

        // pairint is coords (x, y)
        Set<PairInt> visited = new HashSet<PairInt>();

        // contains pairint of coords(x, y)
        java.util.Stack<PairInt> stack = new java.util.Stack<PairInt>();

        // need the smallest indexes to be visited first, so add to the stack
        // (which is LIFO) in reverse
        //O(N)
        for (int i = (sortedThetaRadiusKeys.size() - 1); i > -1; --i) {
            stack.addAll(thetaRadiusPixMap.get(sortedThetaRadiusKeys.get(i)));
        }
        
PairInt chk1 = new PairInt(297, 115);
PairInt chk2 = new PairInt(296, 115);

        Set<PairInt> boundaryPoints = new HashSet<PairInt>();
        
        while (!stack.isEmpty()) {

            PairInt uPoint = stack.pop();
            
if (uPoint.equals(chk1)) {
int z = 1;
} else if (uPoint.equals(chk2)) {
int z = 1;
}
            if (visited.contains(uPoint)) {
                continue;
            }

            PairInt uTR = pixToTRMap.get(uPoint);
            if (uTR == null) {
                uTR = pixToTRMapNotAggr.get(uPoint);
            }

            // search neighbors having same (theta, radius) within tolerance
            // of uPoint for an existing solution in pixToMap.

            PairInt trStore = uTR;
            
            int minDiffT = Integer.MAX_VALUE;
            PairInt minDiffT_TR = null;
            PairInt minDiffT_Coord = null;
            int minDiffR = Integer.MAX_VALUE;
            PairInt minDiffR_TR = null;
            PairInt minDiffR_Coord = null;

            for (int i = 0; i < dxs8.length; ++i) {
                int vX = uPoint.getX() + dxs8[i];
                int vY = uPoint.getY() + dys8[i];
                PairInt vPoint = new PairInt(vX, vY);
                if (!pixToTRMapNotAggr.containsKey(vPoint)) {
                    continue;
                }
                PairInt vTR = pixToTRMap.get(vPoint);
                if (vTR == null) {
                    vTR = pixToTRMapNotAggr.get(vPoint);
                }
                if (vTR != null) {
                    if (vTR.equals(trStore)) {
                        continue;
                    }
                    int diffT = Math.abs(trStore.getX() - vTR.getX());
                    int diffR = Math.abs(trStore.getY() - vTR.getY());
                    if ((diffT <= thetaTol) && (diffR <= radiusTol)) {
                        if (diffT < minDiffT) {
                            minDiffT = diffT;
                            minDiffT_TR = vTR;
                            minDiffT_Coord = vPoint;
                        }
                        if (diffR < minDiffR) {
                            minDiffR = diffR;
                            minDiffR_TR = vTR;
                            minDiffR_Coord = vPoint;
                        }
                    }
                }

                int nU = thetaRadiusPixMap.get(uTR).size();

                if (minDiffT_TR != null && minDiffR_TR != null) {
                    if (minDiffR < minDiffT) {
                        int nV = thetaRadiusPixMap.get(minDiffR_TR).size();
                        if (nV >= nU) {
                            trStore = minDiffR_TR;
                        } else if (nV < nU) {
                            //reassign minDiffR_Coord
                            thetaRadiusPixMap.get(minDiffR_TR).remove(minDiffR_Coord);
                            thetaRadiusPixMap.get(trStore).add(minDiffR_Coord);
                            pixToTRMap.put(minDiffR_Coord, trStore);
                            stack.add(minDiffR_Coord);
                        } else {
                            boundaryPoints.add(uPoint);
                        }
                    } else {
                        int nV = thetaRadiusPixMap.get(minDiffT_TR).size();
                        if (nV >= nU) {
                            trStore = minDiffT_TR;
                        } else if (nV < nU) {
                            //reassign minDiffT_Coord
                            thetaRadiusPixMap.get(minDiffT_TR).remove(minDiffT_Coord);
                            thetaRadiusPixMap.get(trStore).add(minDiffT_Coord);
                            pixToTRMap.put(minDiffT_Coord, trStore);
                            stack.add(minDiffT_Coord);
                        } else {
                            boundaryPoints.add(uPoint);
                        }
                    }
                } else if (minDiffT_TR != null) {
                    int nV = thetaRadiusPixMap.get(minDiffT_TR).size();
                    if (nV >= nU) {
                        trStore = minDiffT_TR;
                    } else if (nV < nU) {
                        //reassign minDiffT_Coord
                        thetaRadiusPixMap.get(minDiffT_TR).remove(minDiffT_Coord);
                        thetaRadiusPixMap.get(trStore).add(minDiffT_Coord);
                        pixToTRMap.put(minDiffT_Coord, trStore);
                        stack.add(minDiffT_Coord);
                    } else {
                        boundaryPoints.add(uPoint);
                    }
                } else if (minDiffR_TR != null) {
                    int nV = thetaRadiusPixMap.get(minDiffR_TR).size();
                    if (nV >= nU) {
                        trStore = minDiffR_TR;
                    } else if (nV < nU) {
                        //reassign minDiffR_Coord
                        thetaRadiusPixMap.get(minDiffR_TR).remove(minDiffR_Coord);
                        thetaRadiusPixMap.get(trStore).add(minDiffR_Coord);
                        pixToTRMap.put(minDiffR_Coord, trStore);
                        stack.add(minDiffR_Coord);
                    } else {
                        boundaryPoints.add(uPoint);
                    }
                }
            }
            
            pixToTRMap.put(uPoint, trStore);

            if (!trStore.equals(uTR)) {
                thetaRadiusPixMap.get(trStore).add(uPoint);
                thetaRadiusPixMap.get(uTR).remove(uPoint);
            }

            visited.add(uPoint);
        }

        return pixToTRMap;
    }

}
