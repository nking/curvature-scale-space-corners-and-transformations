package algorithms.imageProcessing;

import algorithms.MultiArrayMergeSort;
import algorithms.compGeometry.ClosestPairBetweenSets;
import algorithms.imageProcessing.util.PairIntWithIndex;
import algorithms.misc.Misc;
import algorithms.util.PairInt;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.Stack;
import java.util.logging.Logger;

/**
 * a class to aggregate the results of the Hough transform within a given
 * tolerance in theta and radius for contiguous pixels.
 *
 * @author nichole
 */
public class DFSSimilarThetaRadiusGroupsFinder extends AbstractDFSConnectedGroupsFinder {

    private List<Set<PairInt>> sortedGroups = null;

    public DFSSimilarThetaRadiusGroupsFinder() {

        super();

        this.minimumNumberInCluster = 1;
    }

    @Override
    Logger constructLogger() {
        return Logger.getLogger(DFSSimilarThetaRadiusGroupsFinder.class.getName());
    }

    /**
     * aggregate the results of the Hough transform within a given tolerance in
     * theta and radius for contiguous pixels.
     * @param sortedThetaRadiusKeys
     * @param thetaRadiusPixMap
     * @param thetaTol
     * @param radiusTol
     * @param allowGaps if true, performs another operation to join groups
     * within a distance of radiusTol if the closest points have (theta, radius)
     * within tolerances.
     * @return 
     */
    public Map<PairInt, PairInt> findConnectedPointGroups(
        List<PairInt> sortedThetaRadiusKeys,
        Map<PairInt, Set<PairInt>> thetaRadiusPixMap, int thetaTol,
        int radiusTol, boolean allowGaps) {

        Map<PairInt, PairInt> pixToTRMap = findConnectedPointGroupsIterative(
            sortedThetaRadiusKeys, thetaRadiusPixMap, thetaTol, radiusTol);

        createGroups(pixToTRMap, thetaTol, radiusTol);

        prune();
        
        sortedGroups = sortResults();
        
        if (allowGaps) {
            mergeIfGapWithinTolerance(pixToTRMap, thetaTol, radiusTol);
            pruneSortedGroups();
        }

        return pixToTRMap;
    }

    /**
     * aggregate the results of the Hough transform within a given tolerance in
     * theta and radius for contiguous pixels.
     */
    private Map<PairInt, PairInt> findConnectedPointGroupsIterative(
        List<PairInt> sortedThetaRadiusKeys,
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

        while (!stack.isEmpty()) {

            PairInt uPoint = stack.pop();

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

    /**
     * get the groups of points, sorted by decreasing size, after invoking
     * findConnectedPointGroupsIterative.
     * @return
     */
    public List<Set<PairInt>> getSortedGroupsOfPoints() {
        return sortedGroups;
    }

    private List<Set<PairInt>> sortResults() {

        int n = this.getNumberOfGroups();

        int[] c = new int[n];
        int[] indexes = new int[n];

        for (int i = 0; i < n; ++i) {
            c[i] = getXY(i).size();
            indexes[i] = i;
        }

        MultiArrayMergeSort.sortByDecr(c, indexes);

        List<Set<PairInt>> sortedList = new ArrayList<Set<PairInt>>();

        for (int i = 0; i < n; ++i) {
            int idx = indexes[i];
            sortedList.add(getXY(idx));
        }

        return sortedList;
    }

    private void createGroups(Map<PairInt, PairInt> pixToTRMap, int thetaTol,
        int radiusTol) {

        int[] dxs8 = Misc.dx8;
        int[] dys8 = Misc.dy8;

        // pairint is coords (x, y)
        Set<PairInt> visited = new HashSet<PairInt>();

        // contains pairint of coords(x, y)
        java.util.Stack<PairInt> stack = new java.util.Stack<PairInt>();

        // need the smallest indexes to be visited first, so add to the stack
        // (which is LIFO) in reverse
        //O(N)
        for (Entry<PairInt, PairInt> entry : pixToTRMap.entrySet()) {
            stack.add(entry.getKey());
        }

        while (!stack.isEmpty()) {

            PairInt uPoint = stack.pop();

            if (visited.contains(uPoint)) {
                continue;
            }

            PairInt uTR = pixToTRMap.get(uPoint);

            boolean processed = false;

            for (int i = 0; i < dxs8.length; ++i) {

                int vX = uPoint.getX() + dxs8[i];
                int vY = uPoint.getY() + dys8[i];

                PairInt vPoint = new PairInt(vX, vY);

                PairInt vTR = pixToTRMap.get(vPoint);

                if (vTR == null) {
                    continue;
                }

                int diffT = Math.abs(uTR.getX() - vTR.getX());
                int diffR = Math.abs(uTR.getY() - vTR.getY());
                if ((diffT <= thetaTol) && (diffR <= radiusTol)) {
                    processPair(uPoint, vPoint);
                    processed = true;
                    stack.add(vPoint);
                }
            }

            if (!processed) {
                process(uPoint);
            }

            visited.add(uPoint);
        }
    }

    /**
     * iterate over groups to see if the closest pair of points are within
     * a radiusTol distance and have (theta, radius) values similar within 
     * tolerances.
     * <pre>
     * runtime complexity is in the 
     *    worse case: nGroups^2 * (O(N_i * lg2(N_i))
     * but it is expected that most groups will not have two sampling
     * points having (theta, radius) values within 2 times the tolerances,
     * so will be skipped for the closest pair search.
     * </pre>
     * This should be invoked after the member groups have been pruned.
     * 
     * Note that if an image line has high curvature artifacts in it,
     * the merge might not be able within tolerances to join the segments
     * currently.  A method afterwards can be used to look at overlapping
     * edge indexes of the bounds of each group to examine a further merge.
     * @param pixToTRMap
     * @param thetaTol
     * @param radiusTol 
     */
    private void mergeIfGapWithinTolerance(Map<PairInt, PairInt> pixToTRMap, 
        int thetaTol, int radiusTol) {
        
        if (sortedGroups == null) {
            return;
        }
        
        int n = sortedGroups.size();
        
        // search will skip a group when size is less than this
        int minSize = 1;//3;
        
        int distSqTol = radiusTol * radiusTol;
        
        Stack<Integer> stack = new Stack<Integer>();
        for (int i = (n - 1); i > -1; --i) {
            stack.add(Integer.valueOf(i));
        }
        
        Set<Integer> visited = new HashSet<Integer>();
        
        //TODO: walk through merge patterns to determine if the algorithm
        // always terminates (stack becomes empty).
        
        while (!stack.isEmpty()) {
            
            Integer index = stack.pop();
            
            if (visited.contains(index)) {
                continue;
            }
            
            int i = index.intValue();
            
            Set<PairInt> group = sortedGroups.get(i);
            
            if (group.size() < minSize) {
                continue;
            }
            
            PairInt s1 = group.iterator().next();
            PairInt tr1 = pixToTRMap.get(s1);
            
            visited.add(index);
            
            for (int j = (i + 1); j < n; ++j) {
                
                Set<PairInt> group2 = sortedGroups.get(j);
                
                if (group2.size() < minSize) {
                    continue;
                }
                
                PairInt s2 = group2.iterator().next();
                PairInt tr2 = pixToTRMap.get(s2);
                
                if (Math.abs(tr2.getX() - tr1.getX()) > (2 * thetaTol)) {
                    continue;
                }
                if (Math.abs(tr2.getY() - tr1.getY()) > (2 * radiusTol)) {
                    continue;
                }
                
                /*
                for the closest pair,
                adding a variable to represent the two different sets
                and altering the divide and conquer of the closest pair
                algorithm to only
                   compare the members if they are from different sets
                */
                ClosestPairBetweenSets cp = new ClosestPairBetweenSets();
                
                ClosestPairBetweenSets.ClosestPairInt result = 
                    cp.findClosestPair(group, group2);
                
                if ((result == null) || (result.getSeparationSquared() > distSqTol)) {
                    continue;
                }
                
                PairIntWithIndex p1, p2;
                if (result.getPoint0().getPixIndex() == 1) {
                    p1 = result.getPoint0();
                    p2 = result.getPoint1();
                } else {
                    p1 = result.getPoint1();
                    p2 = result.getPoint0();
                }
                
                PairInt pTr1 = pixToTRMap.get(new PairInt(p1.getX(), p1.getY()));
                PairInt pTr2 = pixToTRMap.get(new PairInt(p2.getX(), p2.getY()));
                
                int diffT = Math.abs(pTr1.getX() - pTr2.getX());
                int diffR = Math.abs(pTr1.getY() - pTr2.getY());
                
                if ((diffT > thetaTol) || (diffR > radiusTol)) {
                    continue;
                }
                
                // merge the smaller into the larger
                int nU = group.size();
                int nV = group2.size();
                if (nU >= nV) {
                    // merge v into u
                    group.addAll(group2);
                    group2.clear();
                    stack.add(Integer.valueOf(i));
                } else {
                    // merge u into v
                    group2.addAll(group);
                    group.clear();
                    stack.add(Integer.valueOf(j));
                }
                visited.clear();
            }
        }
        
    }

    private void pruneSortedGroups() {
        
        int n0 = this.sortedGroups.size();
        
        for (int i = (n0 - 1); i > -1; --i) {
            Set<PairInt> group = sortedGroups.get(i);
            if (group.size() == 0) {
                sortedGroups.remove(i);
            }
        }
        
        int n = this.sortedGroups.size();

        int[] c = new int[n];
        int[] indexes = new int[n];

        for (int i = 0; i < n; ++i) {
            c[i] = sortedGroups.get(i).size();
            indexes[i] = i;
        }

        MultiArrayMergeSort.sortByDecr(c, indexes);

        List<Set<PairInt>> sortedList = new ArrayList<Set<PairInt>>();

        for (int i = 0; i < n; ++i) {
            int idx = indexes[i];
            sortedList.add(sortedGroups.get(idx));
        }

        sortedGroups = sortedList;
    }
}
