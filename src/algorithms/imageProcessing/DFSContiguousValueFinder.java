package algorithms.imageProcessing;

import algorithms.misc.Misc;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import java.util.ArrayList;
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
public class DFSContiguousValueFinder {
    
   /**
     * an array to hold each group as an item.  each item contains a key which is an index
     * to arrays indexer.x, indexer.y and this.pointToGroupIndex
     */
    protected List<Set<Integer> > groupMembership = new ArrayList<Set<Integer> >();
    
    protected List<Boolean> groupIsUnbound = new ArrayList<Boolean>();
    
    protected Logger log = null;
    
    protected Set<Integer> visited = new HashSet<Integer>();
    
    /*
     * map w/ key holding indexes for a point to the group it belongs to.
     * note that the point index is relative to indexer.x and indexer.y
     */
    protected Map<Integer, Integer> pointToGroupMap = new
        HashMap<Integer, Integer >();
    
    protected int minimumNumberInCluster = 3;
    
    protected boolean notValue = false;
    
    protected boolean debug = false;
    
    protected final GreyscaleImage img;
    
    protected final Set<PairInt> exclude;
    
    /**
     * uses the 4 neighbor region if true, else the 8-neighbor region
     */
    protected boolean use4Neighbors = true;
    
    public DFSContiguousValueFinder(final GreyscaleImage input) {
        
        this.img = input;
        
        this.log = Logger.getLogger(this.getClass().getName());
        
        this.exclude = new HashSet<PairInt>();
    }
    
    public DFSContiguousValueFinder(final GreyscaleImage input, Set<PairInt>
        mask) {
        
        this.img = input;
        
        this.log = Logger.getLogger(this.getClass().getName());
        
        this.exclude = new HashSet<PairInt>(mask);
    }
        
    public void setMinimumNumberInCluster(int n) {
        this.minimumNumberInCluster = n;
    }
    
    public void setDebug(boolean setDebugToTrue) {
        this.debug = setDebugToTrue;
    }
    
    public void setToUse8Neighbors() {
        use4Neighbors = false;
    }
       
    /**
     * NOTE: to keep the performance reasonable, bin the image so that
     * the number of pixels is below 87000 pixels.
     * 
     * @param pixelValue 
     */
    public void findGroups(final int pixelValue) {
            
        findClusters(pixelValue);
        
        prune();        
    }
    
    /**
     * NOTE: to keep the performance reasonable, bin the image so that
     * the number of pixels is below 87000 pixels.
     * 
     * @param pixelValue 
     */
    public void findGroupsNotThisValue(int pixelValue) {
     
        notValue = true;
                
        findClusters(pixelValue);
        
        prune();        
    }

    protected void findClusters(final int pixelValue) {
        
        // traverse the data by ordered x values 
        //    so can break when exceed critical distance
        
        findClustersIterative(pixelValue);
        
    }
    
    /**
     * NOTE: to keep the performance reasonable, the use of jvm stack has to
     * be reduced.  For default jvm settings, the size of the local array variable
     * within this method frame should be kept to having fewer that 128k items
     * in it, therefore the GreyscaleImage should be binned before use to keep
     * the number of pixels below 128k/1.5 (roughly keep below 87000 pixels).
     * If all pixels are connected, that limit has to be lowered.
     * 
     * @param pixelValue 
     */
    protected void findClustersIterative(final int pixelValue) {
        
        int width = img.getWidth();
        int height = img.getHeight();
        
        int[] dxs;
        int[] dys;
        if (use4Neighbors) {
            dxs = Misc.dx4;
            dys = Misc.dy4;
        } else {
            dxs = Misc.dx8;
            dys = Misc.dy8;
        }
        
        java.util.Stack<Integer> stack = new java.util.Stack<Integer>();
        
        //O(N)
        for (int uIndex = (img.getNPixels() - 1); uIndex > -1; uIndex--) {
            PairInt p = new PairInt(img.getCol(uIndex), img.getRow(uIndex));
            if (!exclude.contains(p)) {
                stack.add(Integer.valueOf(uIndex));
            }
        }
                
        while (!stack.isEmpty()) {

            int uIndex = stack.pop().intValue();
            
            Integer uKey = Integer.valueOf(uIndex);
            
            if (visited.contains(uKey)) {
                continue;
            }
            
            int uPixValue = img.getValue(uIndex);
            
            if ((notValue && (uPixValue == pixelValue)) ||
                (!notValue && (uPixValue != pixelValue))) {
                
                visited.add(uKey);
                
                continue;
            }
            
            int uY = uIndex/width;
            int uX = uIndex - (uY * width);

            //(1 + frac)*O(N) where frac is the fraction added back to stack
            
            for (int i = 0; i < dxs.length; ++i) {
                
                int vX = uX + dxs[i];
                int vY = uY + dys[i];
            
                if ((vX < 0) || (vX > (width - 1)) || (vY < 0) || (vY > (height - 1))) {
                    continue;
                }
                
                if (exclude.contains(new PairInt(vX, vY))) {
                    continue;
                }

                int vIndex = (vY * width) + vX;

                Integer vKey = Integer.valueOf(vIndex);

                int vPixValue = img.getValue(vIndex);

                if ((notValue && (vPixValue == pixelValue)) ||
                    (!notValue && (vPixValue != pixelValue))) {

                    continue;
                }

                processPair(uKey, vKey, false);
                
                // inserting back at the top of the stack assures that the 
                // search continues next from an associated point
                stack.add(vKey);
            }
            
            visited.add(uKey);
        }
    }
    
    protected void processPair(Integer uIdx, Integer vIdx, boolean isBeyondBounds) {
                
        Integer groupId = pointToGroupMap.get(uIdx);
        
        if ((groupId != null) && (pointToGroupMap.get(vIdx) == null)) {
                    
            groupMembership.get(groupId).add(vIdx);
            
            pointToGroupMap.put(vIdx, groupId);
                        
        } else if ((groupId == null) && (pointToGroupMap.get(vIdx) != null)) {

            // usually, u will have been added before v is visited, so this
            // block is rarely used
            
            groupId = pointToGroupMap.get(vIdx);

            groupMembership.get(groupId).add(uIdx);
            
            pointToGroupMap.put(uIdx, groupId);
            
        } else if ((groupId == null) && (pointToGroupMap.get(vIdx) == null)) {
                        
            groupId = Integer.valueOf(groupMembership.size());
            
            pointToGroupMap.put(uIdx, groupId);
            
            pointToGroupMap.put(vIdx, groupId);
            
            Set<Integer> set = new HashSet<Integer>();
            set.add(uIdx);
            set.add(vIdx);
            
            groupMembership.add(set);
            
            groupIsUnbound.add(Boolean.FALSE);
                      
        }
        
        if (isBeyondBounds) {
            groupIsUnbound.set(groupId, Boolean.TRUE);
        }

    }
    
    protected void process(Integer uIdx, boolean isBeyondBounds) {
                
        Integer groupId = pointToGroupMap.get(uIdx);
        
        if (groupId == null) {
                        
            groupId = Integer.valueOf(groupMembership.size());
            
            pointToGroupMap.put(uIdx, groupId);
            
            Set<Integer> set = new HashSet<Integer>();
            set.add(uIdx);
            
            groupMembership.add(set);
            
            groupIsUnbound.add(Boolean.FALSE);
                      
        }
        
        if (isBeyondBounds) {
            groupIsUnbound.set(groupId, Boolean.TRUE);
        }

    }
            
    public List<Set<Integer> > getGroupMembershipList() {
        return groupMembership;
    }

    public int getNumberOfGroups() {
        return groupMembership.size();
    }

    public Map<Integer, Integer> getPointToGroupIndexes() {
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

        log.fine("PRUNE: nGroups before prune =" + groupMembership.size());

        // iterate backwards so can move items up without conflict with iterator
        for (int i = (groupMembership.size() - 1); i > -1; i--) {
            
            Set<Integer> group = groupMembership.get(i);
            
            int count = group.size();
            
            log.finest("  group " + i + " has " + count 
                + " members before prune (min=" + minimumNumberInCluster + ")" 
                + ".  doPrune=" +
                ((count < minimumNumberInCluster) || groupIsUnbound.get(i).booleanValue())
                + "   count=" + count + " minimumNumberInCluster=" + minimumNumberInCluster
                + " isUnbound=" + groupIsUnbound.get(i).booleanValue()
            );
            
            if ((count < minimumNumberInCluster) || groupIsUnbound.get(i).booleanValue()) {
                             
                // remove this group and move up all groups w/ index > i by one index
                for (int j = (i + 1); j < groupMembership.size(); j++) {
                    
                    int newGroupId = j - 1;
                                        
                    // update members in pointToGroupIndex
                    Set<Integer> latest = groupMembership.get(j);
                    //
                    for (Integer p : latest) {
                        pointToGroupMap.put(p, Integer.valueOf(newGroupId));
                    }
                } 
               
                Set<Integer> removed = groupMembership.remove(i);
                
            }
        }
   
        log.fine("number of groups after prune=" + groupMembership.size());
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
    
    public PairIntArray getXY(int groupId) {
        
        if (groupMembership.size() == 0) {
            return new PairIntArray();
        }
        if (groupId > (groupMembership.size() - 1) || (groupId < 0)) {
            throw new IllegalArgumentException("groupId=" + groupId 
            + " is outside of range of nGroups=" + groupMembership.size());
        }
        
        Set<Integer> indexes = getIndexes(groupId);
        
        PairIntArray xy = new PairIntArray(indexes.size());
                
        for (Integer index : indexes) {
                        
            img.getXY(xy, index.intValue());
        }
        
        return xy;
    }
    
    public void getXY(final int groupId, final Set<PairInt> output) {
        
        if (groupMembership.isEmpty()) {
            return;
        }
        if (groupId > (groupMembership.size() - 1) || (groupId < 0)) {
            throw new IllegalArgumentException("groupId=" + groupId 
            + " is outside of range of nGroups=" + groupMembership.size());
        }
        
        Set<Integer> indexes = groupMembership.get(groupId);

        for (Integer index : indexes) {
            
            int idx = index.intValue();
            
            int x = img.getCol(idx);
            int y = img.getRow(idx);
            
            PairInt p = new PairInt(x, y);
                       
            output.add(p);
        }        
    }
    
    public Set<Integer> getIndexes(int groupId) {
        
        if (groupMembership.isEmpty()) {
            return new HashSet<Integer>();
        }
        if (groupId > (groupMembership.size() - 1) || (groupId < 0)) {
            throw new IllegalArgumentException("groupId=" + groupId 
            + " is outside of range of nGroups=" + groupMembership.size());
        }
               
        return groupMembership.get(groupId);
    }
 
}
