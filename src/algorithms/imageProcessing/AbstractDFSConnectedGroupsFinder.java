package algorithms.imageProcessing;

import algorithms.util.PairInt;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;

/**
 *
 * @author nichole
 */
public abstract class AbstractDFSConnectedGroupsFinder {
    
    /**
     * an array to hold each group as an item.  each item contains a key which is an index
     * to arrays indexer.x, indexer.y and this.pointToGroupIndex
     */
    protected List<Set<PairInt>> groupMembership = new ArrayList<Set<PairInt>>();
        
    protected int minimumNumberInCluster = 3;
    
    /*
     * map w/ key holding indexes for a point to the group it belongs to.
     * note that the point index is relative to indexer.x and indexer.y
     */
    //protected int[] pointToGroupIndex = null;
    protected Map<PairInt, Integer> pointToGroupMap = new HashMap<PairInt, Integer>();
    
    protected Logger log = null;
    
    protected boolean debug = false;

    public AbstractDFSConnectedGroupsFinder() {
        log = constructLogger();
    }
    
    abstract Logger constructLogger();

    public void setDebug(boolean setDebugToTrue) {
        this.debug = setDebugToTrue;
    }

    protected void processPair(PairInt uPoint, PairInt vPoint) {
        Integer uGroupId = pointToGroupMap.get(uPoint);
        Integer vGroupId = pointToGroupMap.get(vPoint);
        if (uGroupId == null) {
            if (vGroupId == null) {
                uGroupId = Integer.valueOf(groupMembership.size());
                pointToGroupMap.put(uPoint, uGroupId);
                pointToGroupMap.put(vPoint, uGroupId);
                Set<PairInt> set = new HashSet<PairInt>();
                set.add(uPoint);
                set.add(vPoint);
                groupMembership.add(set);
            } else {
                groupMembership.get(vGroupId).add(uPoint);
                pointToGroupMap.put(uPoint, vGroupId);
            }
        } else {
            if (vGroupId == null) {
                groupMembership.get(uGroupId).add(vPoint);
                pointToGroupMap.put(vPoint, uGroupId);
            } else {
                // else vGroupId == uGroupId
                if (!vGroupId.equals(uGroupId)) {
                    // merge groups
                    Set<PairInt> uGroup = groupMembership.get(uGroupId);
                    Set<PairInt> vGroup = groupMembership.get(vGroupId);
                    int nU = uGroup.size();
                    int nV = vGroup.size();
                    if (nU >= nV) {
                        // merge v into u
                        uGroup.addAll(vGroup);
                        for (PairInt p : vGroup) {
                            pointToGroupMap.put(p, uGroupId);
                        }
                        vGroup.clear();
                        int z = 1;
                    } else {
                        // merge u into v
                        vGroup.addAll(uGroup);
                        for (PairInt p : uGroup) {
                            pointToGroupMap.put(p, vGroupId);
                        }
                        uGroup.clear();
                        int z = 1;
                    }
                }
            }
        }
    }

    protected void process(PairInt uPoint) {
        Integer groupId = pointToGroupMap.get(uPoint);
        if (groupId == null) {
            groupId = Integer.valueOf(groupMembership.size());
            pointToGroupMap.put(uPoint, groupId);
            Set<PairInt> set = new HashSet<PairInt>();
            set.add(uPoint);
            groupMembership.add(set);
        }
    }

    public void setMinimumNumberInCluster(int n) {
        this.minimumNumberInCluster = n;
    }
  
    public List<Set<PairInt>> getGroupMembershipList() {
        return new ArrayList<Set<PairInt>>(groupMembership);
    }

    public int getNumberOfGroups() {
        return groupMembership.size();
    }

    public Map<PairInt, Integer> getPointToGroupIndexes() {
        return pointToGroupMap;
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
            Set<PairInt> group = groupMembership.get(i);
            int count = group.size();
            log.finest("  group " + i + " has " + count + " members before prune (min=" + minimumNumberInCluster + ")");
            if (count < minimumNumberInCluster) {
                // remove this group and move up all groups w/ index > i by one index
                for (int j = i + 1; j < groupMembership.size(); j++) {
                    int newGroupId = j - 1;
                    // update members in pointToGroupIndex
                    Set<PairInt> latest = groupMembership.get(j);
                    for (PairInt p : latest) {
                        pointToGroupMap.put(p, Integer.valueOf(newGroupId));
                    }
                }
                Set<PairInt> removed = groupMembership.remove(i);
            }
        }
        log.finest("number of groups after prune=" + groupMembership.size());
    }

    public int getNumberofGroupMembers(int groupId) {
        if (groupMembership.isEmpty()) {
            return 0;
        }
        if (groupId > (groupMembership.size() - 1) || (groupId < 0)) {
            throw new IllegalArgumentException("groupId=" + groupId + " is outside of range of nGroups=" + groupMembership.size());
        }
        return groupMembership.get(groupId).size();
    }

    public Set<PairInt> getXY(int groupId) {
        if (groupMembership.isEmpty()) {
            return new HashSet<PairInt>();
        }
        if (groupId > (groupMembership.size() - 1) || (groupId < 0)) {
            throw new IllegalArgumentException("groupId=" + groupId + " is outside of range of nGroups=" + groupMembership.size());
        }
        Set<PairInt> set = groupMembership.get(groupId);
        return set;
    }
    
}
