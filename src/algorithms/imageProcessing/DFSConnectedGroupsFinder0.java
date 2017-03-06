package algorithms.imageProcessing;

import algorithms.misc.Misc;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

/**
 * adapted from
 * https://code.google.com/p/two-point-correlation/source/browse/src/main/java/algorithms/compGeometry/clustering/twopointcorrelation/AbstractGroupFinder.java
 * and
 * https://code.google.com/p/two-point-correlation/source/browse/src/main/java/algorithms/compGeometry/clustering/twopointcorrelation/DFSGroupFinder.java
 * under MIT License (MIT), Nichole King 2013
 * 
 * @author nichole
 */
public class DFSConnectedGroupsFinder0 {
    
    protected boolean notValue = false;
    
    /**
     * uses the 4 neighbor region if true, else the 8-neighbor region
     */
    protected boolean use4Neighbors = true;
    
    /**
     * an array to hold each group as an item.  each item contains a key which is an index
     * to arrays indexer.x, indexer.y and this.pointToGroupIndex
     */
    protected List<TIntSet> groupMembership = new ArrayList<TIntSet>();
    
    protected int minimumNumberInCluster = 3;
    
    /*
     * map w/ key holding indexes for a point to the group it belongs to.
     * note that the point index is relative to indexer.x and indexer.y
     */
    //protected int[] pointToGroupIndex = null;
    
    // key = pixIndex, value = group index
    protected TIntIntMap pointToGroupMap = new TIntIntHashMap();
    
    private final int imgWidth;
    
    protected Logger log = Logger.getLogger(this.getClass().getName());
    
    protected boolean debug = false;

    public DFSConnectedGroupsFinder0(int imageWidth) {
        imgWidth = imageWidth;
    }
    
    public void setDebug(boolean setDebugToTrue) {
        this.debug = setDebugToTrue;
    }

    public void setMinimumNumberInCluster(int n) {
        this.minimumNumberInCluster = n;
    }
    
    public void setToUse8Neighbors() {
        use4Neighbors = false;
    }
    
    public void findConnectedPointGroups(TIntSet pixIdxs) {
        
        //TODO: make an adjacency map at this point and replace
        // the loop over neighbor offsets with loop over adjacent points
        // to speed up the algorithm.        
        
        findClustersIterative(pixIdxs);
        
        prune();        
    }

    protected void findClustersIterative(TIntSet pixIdxs) {
        
        if (pixIdxs.isEmpty()) {
            return;
        }
        
        int[] dxs;
        int[] dys;
        if (use4Neighbors) {
            dxs = Misc.dx4;
            dys = Misc.dy4;
        } else {
            dxs = Misc.dx8;
            dys = Misc.dy8;
        }
        
        TIntSet visited = new TIntHashSet();
        
        java.util.Stack<Integer> stack = new java.util.Stack<Integer>();
        
        //O(N)
        if (use4Neighbors) {
            TIntIterator iter = pixIdxs.iterator();
            while (iter.hasNext()) {
                stack.add(Integer.valueOf(iter.next()));
            }
        } else {
            stack.add(Integer.valueOf(pixIdxs.iterator().next()));
        }
        
        while (!stack.isEmpty()) {

            int uPoint = stack.pop().intValue();
            
            if (visited.contains(uPoint)) {
                continue;
            }

            int uY = uPoint/imgWidth;
            int uX = uPoint - (uY * imgWidth);
            
            boolean foundANeighbor = false;
            
            //(1 + frac)*O(N) where frac is the fraction added back to stack
            
            for (int i = 0; i < dxs.length; ++i) {
                
                int vX = uX + dxs[i];
                int vY = uY + dys[i];
            
                int vPoint = (vY * imgWidth) + vX;

                if (vPoint == uPoint) {
                    continue;
                }

                if (!pixIdxs.contains(vPoint)) {
                    continue;
                }

                processPair(uPoint, vPoint);

                if (!use4Neighbors) {
                    stack.add(vPoint);
                }

                foundANeighbor = true;
            }
            
            if (!foundANeighbor && (minimumNumberInCluster == 1)) {                
                process(uPoint);
            }
            
            visited.add(uPoint);
        }
    }
  
    protected void processPair(int uPoint, int vPoint) {
        int uGroupId = pointToGroupMap.containsKey(uPoint) ?
            pointToGroupMap.get(uPoint) : -1;
        int vGroupId = pointToGroupMap.containsKey(vPoint) ?
            pointToGroupMap.get(vPoint) : -1;
        if (uGroupId == -1) {
            if (vGroupId == -1) {
                uGroupId = groupMembership.size();
                pointToGroupMap.put(uPoint, uGroupId);
                pointToGroupMap.put(vPoint, uGroupId);
                TIntSet set = new TIntHashSet();
                set.add(uPoint);
                set.add(vPoint);
                groupMembership.add(set);
            } else {
                groupMembership.get(vGroupId).add(uPoint);
                pointToGroupMap.put(uPoint, vGroupId);
            }
        } else {
            if (vGroupId == -1) {
                groupMembership.get(uGroupId).add(vPoint);
                pointToGroupMap.put(vPoint, uGroupId);
            } else {
                // else both are not null
                if (vGroupId != uGroupId) {
                    // merge groups
                    //TODO: this is where use of disjoint forest could help
                    TIntSet uGroup = groupMembership.get(uGroupId);
                    TIntSet vGroup = groupMembership.get(vGroupId);
                    int nU = uGroup.size();
                    int nV = vGroup.size();
                    if (nU >= nV) {
                        // merge v into u
                        uGroup.addAll(vGroup);
                        TIntIterator iter = vGroup.iterator();
                        while (iter.hasNext()) {
                            int pixIdx = iter.next();
                            pointToGroupMap.put(pixIdx, uGroupId);
                        }
                        vGroup.clear();
                    } else {
                        // merge u into v
                        vGroup.addAll(uGroup);
                        TIntIterator iter = uGroup.iterator();
                        while (iter.hasNext()) {
                            int pixIdx = iter.next();
                            pointToGroupMap.put(pixIdx, vGroupId);
                        }
                        uGroup.clear();
                    }
                }
            }
        }
    }

    protected void process(int uPoint) {
        int groupId = pointToGroupMap.get(uPoint);
        if (!pointToGroupMap.containsKey(uPoint)) {
            groupId = groupMembership.size();
            pointToGroupMap.put(uPoint, groupId);
            TIntSet set = new TIntHashSet();
            set.add(uPoint);
            groupMembership.add(set);
        }
    }
  
    public List<TIntSet> getGroupMembershipList() {
        return groupMembership;
    }

    public int getNumberOfGroups() {
        return groupMembership.size();
    }

    public TIntIntMap getPointToGroupIndexes() {
        return pointToGroupMap;
    }

    /**
     * remove groups smaller than minimumNumberInCluster
     */
    protected void prune() {
        log.finest("number of groups before prune=" + groupMembership.size());
        //TODO: this is where use of a disjoint forest could help
        /*
         * [------] 0
         * [------] 1 <---- too few
         * [------] 2
         */
        // iterate backwards so can move items up without conflict with iterator
        for (int i = groupMembership.size() - 1; i > -1; i--) {
            TIntSet group = groupMembership.get(i);
            int count = group.size();
            log.finest("  group " + i + " has " + count + " members before prune (min=" + minimumNumberInCluster + ")");
            if (count < minimumNumberInCluster) {
                // remove this group and move up all groups w/ index > i by one index
                for (int j = i + 1; j < groupMembership.size(); j++) {
                    int newGroupId = j - 1;
                    // update members in pointToGroupIndex
                    TIntSet latest = groupMembership.get(j);
                    TIntIterator iter2 = latest.iterator();
                    while (iter2.hasNext()) {
                        int pixIdx = iter2.next();
                        pointToGroupMap.put(pixIdx, newGroupId);
                    }
                }
                TIntSet removed = groupMembership.remove(i);
            }
        }
        log.finest("number of groups after prune=" + groupMembership.size());
    }
    
    public TIntIntMap createPointIndexMap() {
        
        TIntIntMap ptIdxMap = new TIntIntHashMap();
        
        int n = getNumberOfGroups();
        for (int i = 0; i < n; ++i) {
            
            TIntSet set = getXY(i);
            
            TIntIterator iter = set.iterator();
            while (iter.hasNext()) {
                int pixIdx = iter.next();
                ptIdxMap.put(pixIdx, i);
            }
        }
        
        return ptIdxMap;
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

    public TIntSet getXY(int groupId) {
        if (groupMembership.isEmpty()) {
            return new TIntHashSet();
        }
        if (groupId > (groupMembership.size() - 1) || (groupId < 0)) {
            throw new IllegalArgumentException("groupId=" + groupId + " is outside of range of nGroups=" + groupMembership.size());
        }
        TIntSet set = groupMembership.get(groupId);
        return set;
    }
    
}
