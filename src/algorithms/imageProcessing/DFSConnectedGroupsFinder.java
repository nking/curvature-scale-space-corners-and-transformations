package algorithms.imageProcessing;

import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.SimpleLinkedListNode;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
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
public class DFSConnectedGroupsFinder {
    
    /**
     * an array to hold each group as an item.  each item contains a key which is an index
     * to arrays indexer.x, indexer.y and this.pointToGroupIndex
     */
    protected List<Set<PairInt> > groupMembership = new ArrayList<Set<PairInt> >();
    
    protected Logger log = null;
    
    protected Set<PairInt> visited = new HashSet<PairInt>();
    
    /*
     * map w/ key holding indexes for a point to the group it belongs to.
     * note that the point index is relative to indexer.x and indexer.y
     */
    //protected int[] pointToGroupIndex = null;
    protected Map<PairInt, Integer> pointToGroupMap = new
        HashMap<PairInt, Integer >();
    
    protected int minimumNumberInCluster = 3;
    
    protected boolean notValue = false;
    
    protected boolean debug = false;
        
    public DFSConnectedGroupsFinder() {
                
        this.log = Logger.getLogger(this.getClass().getName());
    }
        
    public void setMinimumNumberInCluster(int n) {
        this.minimumNumberInCluster = n;
    }
    
    public void setDebug(boolean setDebugToTrue) {
        this.debug = setDebugToTrue;
    }
  
    public void findConnectedPointGroups(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
                
        findClustersIterative(points, imageWidth, imageHeight);
        
        prune();        
    }

    protected void findClustersIterative(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
        
        java.util.Stack<PairInt> stack = new java.util.Stack<PairInt>();
        
        //O(N)
        for (PairInt p : points) {
            stack.add(p);
        }
               
        visited.add(stack.peek());

        while (!stack.isEmpty()) {

            PairInt uPoint = stack.pop();
            
            int uX = uPoint.getX();
            int uY = uPoint.getY();

            //(1 + frac)*O(N) where frac is the fraction added back to stack
            
            for (int vX = (uX - 1); vX <= (uX + 1); vX++) {
                
                if ((vX < 0) || (vX > (imageWidth - 1))) {
                    continue;
                }
                
                for (int vY = (uY - 1); vY <= (uY + 1); vY++) {
                    if ((vY < 0) || (vY > (imageHeight - 1))) {
                        continue;
                    }
                    
                    PairInt vPoint = new PairInt(vX, vY);
                    
                    if (vPoint.equals(uPoint)) {
                        continue;
                    }
                    
                    if (!points.contains(vPoint)) {
                        continue;
                    }
                    
                    if (visited.contains(vPoint)) {
                        continue;
                    }
                    
                    visited.add(vPoint);
                    
                    processPair(uPoint, vPoint);
                
                    stack.add(vPoint);
                }
            }
        }
       
    }

    protected void processPair(PairInt uPoint, PairInt vPoint) {
              
        Integer groupId = pointToGroupMap.get(uPoint);
        
        if ((groupId != null) && (pointToGroupMap.get(vPoint) == null)) {
                    
            groupMembership.get(groupId).add(vPoint);
            
            pointToGroupMap.put(vPoint, groupId);
                        
        } else if ((groupId == null) && (pointToGroupMap.get(vPoint) != null)) {

            groupId = pointToGroupMap.get(vPoint);

            groupMembership.get(groupId).add(uPoint);
            
            pointToGroupMap.put(uPoint, groupId);
            
        } else if ((groupId == null) && (pointToGroupMap.get(vPoint) == null)) {
                        
            groupId = Integer.valueOf(groupMembership.size());
            
            pointToGroupMap.put(uPoint, groupId);
            
            pointToGroupMap.put(vPoint, groupId);
            
            Set<PairInt> set = new HashSet<PairInt>();
            set.add(uPoint);
            set.add(vPoint);
            
            groupMembership.add(set);
                      
        }
       
        int z = 1;
    }

    public List<Set<PairInt> > getGroupMembershipList() {
        return new ArrayList<Set<PairInt> >(groupMembership);
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
        for (int i = (groupMembership.size() - 1); i > -1; i--) {
            
            Set<PairInt> group = groupMembership.get(i);
            
            int count = group.size();
            
            log.finest("  group " + i + " has " + count 
                + " members before prune (min=" + minimumNumberInCluster + ")");
            
            if (count < minimumNumberInCluster) {
            
                // remove this group and move up all groups w/ index > i by one index
                for (int j = (i + 1); j < groupMembership.size(); j++) {
                    
                    int newGroupId = j - 1;
                    
                    groupMembership.set(newGroupId, groupMembership.get(j));
                    
                    // update members in pointToGroupIndex
                    Set<PairInt> latest = groupMembership.get(j);
                    
                    for (PairInt p : latest) {
                        pointToGroupMap.put(p, Integer.valueOf(newGroupId));
                    }
                }
                
                groupMembership.remove(i);
            }
        }
   
        log.finest("number of groups after prune=" + groupMembership.size());
    }
  
    public int getNumberofGroupMembers(int groupId) {
        
        if (groupMembership.isEmpty()) {
            return 0;
        }
        if (groupId > (groupMembership.size() - 1) || (groupId < 0)) {
            throw new IllegalArgumentException("groupId=" + groupId 
            + " is outside of range of nGroups=" + groupMembership.size());
        }
                
        return groupMembership.get(groupId).size();
    }
    
    public Set<PairInt> getXY(int groupId) {
        
        if (groupMembership.isEmpty()) {
            return new HashSet<PairInt>();
        }
        if (groupId > (groupMembership.size() - 1) || (groupId < 0)) {
            throw new IllegalArgumentException("groupId=" + groupId 
            + " is outside of range of nGroups=" + groupMembership.size());
        }
        
        Set<PairInt> set = groupMembership.get(groupId);
       
        return set;
    }
  
}
