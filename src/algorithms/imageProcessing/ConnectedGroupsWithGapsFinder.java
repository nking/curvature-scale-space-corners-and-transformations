package algorithms.imageProcessing;

import algorithms.compGeometry.ClosestPairBetweenSets;
import algorithms.imageProcessing.util.AngleUtil;
import algorithms.imageProcessing.util.PairIntWithIndex;
import algorithms.misc.Misc;
import algorithms.util.PairInt;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;

/**
 *  not yet tested or ready for use
 * 
 * if gapAllowed == 0, runtime complexity is 8 * O(N_all_points).
 * if gapAllowed > 0, runtime complexity depends upon the characteristics of the
 * data: the min and max bounds are used to reduce comparisons, and the theta
 * is used to filter also. the remaining list groups are compared using
 * a closest pair between sets algorithm in pairs visited using a DFS patter,
 * so the runtime complexity is then roughly
 *   nListPairs * O(N_points_in_both_lists * lg_2(N_points_in_both_lists))
 *   times a factor of about 6.
 * 
 * @author nichole
 */
public class ConnectedGroupsWithGapsFinder {
    
    /**
     * an array to hold each grouped set of indexes as an item.
     */
    protected List<Set<Integer>> groupMembership = new ArrayList<Set<Integer>>();
        
    protected int minimumNumberInCluster = 2;
    
    /*
     * map w/ key holding indexes for a point to the group it belongs to.
     * note that the point index is relative to indexer.x and indexer.y
     */
    protected Map<Integer, Integer> inputIndexToGroupIndexMap = new HashMap<Integer, Integer>();
    
    protected Logger log = Logger.getLogger(this.getClass().getName());
    
    protected boolean debug = false;
    
    public ConnectedGroupsWithGapsFinder() {
                
    }
    
    public void setMinimumNumberInCluster(int n) {
        this.minimumNumberInCluster = n;
    }
    
    public void findConnectedGroups(List<Set<PairInt>> inputList, 
        List<Integer> thetas, double allowedGap) {
                
        findClustersIterative(inputList, thetas, allowedGap);
        
        prune();        
    }
      
    protected void findClustersIterative(List<Set<PairInt>> inputList, 
        List<Integer> thetas, double allowedGap) {
        
        if (inputList.isEmpty()) {
            return;
        }
        
        int n = inputList.size();
        List<Bounds> bounds = new ArrayList<Bounds>();
        for (int i = 0; i < n; ++i) {
            Collection<PairInt> groupI = inputList.get(i);
            //xMin, xMax, yMin, yMax
            int[] minMaxXY = algorithms.misc.MiscMath.findMinMaxXY(groupI);
            Bounds b = new Bounds();
            b.minX = minMaxXY[0];
            b.maxX = minMaxXY[1];
            b.minY = minMaxXY[2];
            b.maxY = minMaxXY[3];
            bounds.add(b);
        }
        
        Set<Integer> visited = new HashSet<Integer>();
        
        java.util.Stack<Integer> stack = new java.util.Stack<Integer>();
        
        boolean foundANeighbor = false;
        
        for (int i = (n - 1); i > 0; i--) {
            stack.add(Integer.valueOf(i));
        }
               
        while (!stack.isEmpty()) {

            Integer uIndex = stack.pop();
            
            if (visited.contains(uIndex)) {
                continue;
            }
            
            int uIdx = uIndex.intValue();

            Set<PairInt> groupI = inputList.get(uIdx);
            Integer thetaI = thetas.get(uIdx);
            Bounds boundsI = bounds.get(uIdx);
            
            for (int vIdx = 0; vIdx < n ; ++vIdx) {
                if (uIdx == vIdx) {
                    continue;
                }
                Integer vIndex = Integer.valueOf(vIdx);
                
                Bounds boundsJ = bounds.get(vIdx);
                if ((boundsI.maxX <  (boundsJ.minX - allowedGap)) || 
                    (boundsI.minX >  (boundsJ.maxX + allowedGap)) ||
                    (boundsI.maxY <  (boundsJ.minY - allowedGap)) || 
                    (boundsI.minY >  (boundsJ.maxY + allowedGap))
                    ) {
                    continue;
                }
                Set<PairInt> groupJ = inputList.get(vIdx);
                Integer thetaJ = thetas.get(vIdx);
                double diff = AngleUtil.getAngleDifference(thetaI.floatValue(), 
                    thetaJ.floatValue());
                if (Math.abs(diff) > 44) {
                    continue;
                }
                ClosestPairBetweenSets cpFinder = new ClosestPairBetweenSets();
                ClosestPairBetweenSets.ClosestPairInt cp = cpFinder.findClosestPair(groupI, groupJ);
                PairIntWithIndex p1 = cp.getPoint0();
                if (p1 == null) {
                    continue;
                }
                if (Math.sqrt(cp.getSeparationSquared()) > allowedGap) {
                    continue;
                }
            
                processPair(uIdx, vIdx);

                stack.add(vIndex);

                foundANeighbor = true;
            }
            
            if (!foundANeighbor && (minimumNumberInCluster == 1)) {                
                process(uIdx);
            }
            
            visited.add(uIndex);
        }
    }
    
    protected void processPair(Integer uIndex, Integer vIndex) {
        Integer uGroupId = inputIndexToGroupIndexMap.get(uIndex);
        Integer vGroupId = inputIndexToGroupIndexMap.get(vIndex);
        if (uGroupId == null) {
            if (vGroupId == null) {
                uGroupId = Integer.valueOf(groupMembership.size());
                inputIndexToGroupIndexMap.put(uIndex, uGroupId);
                inputIndexToGroupIndexMap.put(vIndex, uGroupId);
                Set<Integer> set = new HashSet<Integer>();
                set.add(uIndex);
                set.add(vIndex);
                groupMembership.add(set);
            } else {
                groupMembership.get(vGroupId).add(uIndex);
                inputIndexToGroupIndexMap.put(uIndex, vGroupId);
            }
        } else {
            if (vGroupId == null) {
                groupMembership.get(uGroupId).add(vIndex);
                inputIndexToGroupIndexMap.put(vIndex, uGroupId);
            } else {
                // else vGroupId == uGroupId
                if (!vGroupId.equals(uGroupId)) {
                    // merge groups
                    Set<Integer> uGroup = groupMembership.get(uGroupId);
                    Set<Integer> vGroup = groupMembership.get(vGroupId);
                    int nU = uGroup.size();
                    int nV = vGroup.size();
                    if (nU >= nV) {
                        // merge v into u
                        uGroup.addAll(vGroup);
                        for (Integer p : vGroup) {
                            inputIndexToGroupIndexMap.put(p, uGroupId);
                        }
                        vGroup.clear();
                    } else {
                        // merge u into v
                        vGroup.addAll(uGroup);
                        for (Integer p : uGroup) {
                            inputIndexToGroupIndexMap.put(p, vGroupId);
                        }
                        uGroup.clear();
                    }
                }
            }
        }
    }
     
    protected void process(Integer uIndex) {
        Integer groupId = inputIndexToGroupIndexMap.get(uIndex);
        if (groupId == null) {
            groupId = Integer.valueOf(groupMembership.size());
            inputIndexToGroupIndexMap.put(uIndex, groupId);
            Set<Integer> set = new HashSet<Integer>();
            set.add(uIndex);
            groupMembership.add(set);
        }
    }
    
    /**
     * remove groups smaller than minimumNumberInCluster
     */
    protected void prune() {
        log.finest("number of groups before prune=" + groupMembership.size());
        //TODO: the data structures used could be written at the expense
        // of space complexity to reduce changes needed when group number
        // changes
        /*
         * [------] 0
         * [------] 1 <---- too few
         * [------] 2
         */
        // iterate backwards so can move items up without conflict with iterator
        for (int i = groupMembership.size() - 1; i > -1; i--) {
            Set<Integer> group = groupMembership.get(i);
            int count = group.size();
            log.finest("  group " + i + " has " + count + " members before prune (min=" + minimumNumberInCluster + ")");
            if (count < minimumNumberInCluster) {
                // remove this group and move up all groups w/ index > i by one index
                for (int j = i + 1; j < groupMembership.size(); j++) {
                    int newGroupId = j - 1;
                    // update members in pointToGroupIndex
                    Set<Integer> latest = groupMembership.get(j);
                    for (Integer p : latest) {
                        inputIndexToGroupIndexMap.put(p, Integer.valueOf(newGroupId));
                    }
                }
                Set<Integer> removed = groupMembership.remove(i);
            }
        }
        log.finest("number of groups after prune=" + groupMembership.size());
    }
    
    public int getNumberOfGroups() {
        return groupMembership.size();
    }
    
    public Set<Integer> getGroupedIndexes(int groupId) {
        if (groupMembership.isEmpty()) {
            return new HashSet<Integer>();
        }
        if (groupId > (groupMembership.size() - 1) || (groupId < 0)) {
            throw new IllegalArgumentException("groupId=" + groupId + " is outside of range of nGroups=" + groupMembership.size());
        }
        Set<Integer> set = groupMembership.get(groupId);
        return set;
    }
    
    private class Bounds {
        public int minX;
        public int maxX;
        public int minY;
        public int maxY;
    }
}
