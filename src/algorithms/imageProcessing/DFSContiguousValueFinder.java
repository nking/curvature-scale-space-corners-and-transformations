package algorithms.imageProcessing;

import algorithms.util.PairIntArray;
import algorithms.util.SimpleLinkedListNode;
import algorithms.util.StackInt;
import java.util.Arrays;
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
    protected SimpleLinkedListNode[] groupMembership = null;

    protected int nGroups = 0;
    
    protected Logger log = null;
    
    // 0 = unvisited, 1 = processing, 2 = visited
    protected int[] color = null;
    
    /*
     * array holding indexes for a point to the group it belongs to.
     * note that the point index is relative to indexer.x and indexer.y
     */
    protected int[] pointToGroupIndex = null;
    
    protected int minimumNumberInCluster = 3;
    
    protected boolean notValue = false;
    
    protected boolean debug = false;
    
    protected final GreyscaleImage img;
    
    public DFSContiguousValueFinder(final GreyscaleImage input) {
        
        this.img = input;
        
        this.log = Logger.getLogger(this.getClass().getName());
    }
        
    public void setMinimumNumberInCluster(int n) {
        this.minimumNumberInCluster = n;
    }
    
    public void setDebug(boolean setDebugToTrue) {
        this.debug = setDebugToTrue;
    }
    
    protected void initializeVariables() {

        pointToGroupIndex = new int[img.getNPixels()];
        
        Arrays.fill(pointToGroupIndex, -1);

        groupMembership = new SimpleLinkedListNode[10];
        
        for (int i = 0; i < groupMembership.length; i++) {
            
            groupMembership[i] = new SimpleLinkedListNode();
        }
    }
        
    public void findGroups(int pixelValue) {
                
        initializeVariables();
        
        findClusters(pixelValue);
        
        prune();        
    }
    
    public void findGroupsNotThisValue(int pixelValue) {
     
        notValue = true;
        
        initializeVariables();
        
        findClusters(pixelValue);
        
        prune();        
    }

    protected void findClusters(int pixelValue) {
        
        // traverse the data by ordered x values 
        //    so can break when exceed critical distance
                                
        // 0 = unvisited, 1 = processing, 2 = visited
        color = new int[img.getNPixels()];
        
        findClustersIterative(pixelValue);
        
    }
    
    protected void findClustersIterative(int pixelValue) {
        
        int width = img.getWidth();
        int height = img.getHeight();
        
        StackInt stack = new StackInt();
        
        /*for (int uIndex = (img.getNPixels() - 1); uIndex > -1; uIndex--) {
            stack.insert(uIndex);
        }*/
        
        for (int col = (width - 1); col > -1; col--) {
            for (int row = (img.getHeight() - 1); row > -1; row--) {
                int idx = (row * width) + col;
                stack.insert(idx);
            }
        }
               
        color[stack.peek()] = 2;
        
        while (!stack.isEmpty()) {
            
            SimpleLinkedListNode uNode = stack.pop();
            
            int uIndex = uNode.getKey();
            
            int uPixValue = img.getValue(uIndex);
            
            if ((notValue && (uPixValue == pixelValue)) ||
                (!notValue && (uPixValue != pixelValue))) {
                
                color[uIndex] = 2;
                continue;
            }
            
            /*
            note to speed this up a little, have copied out the index to
            row and col relationship from GreyscaleImage.
            */            
            int uY = uIndex/width;
            int uX = uIndex - (uY * width);
            
            for (int vX = (uX - 1); vX <= (uX + 1); vX++) {
                if ((vX < 0) || (vX > (width - 1))) {
                    continue;
                }
                for (int vY = (uY - 1); vY <= (uY + 1); vY++) {
                    if ((vY < 0) || (vY > (height - 1))) {
                        continue;
                    }
                    
                    int vIndex = (vY * width) + vX;
                    
                    if (color[vIndex] != 0 || (uIndex == vIndex)) {
                        continue;
                    }
                    
                    //TODO: could consider not allowing diagonal connections
                    
                    color[vIndex] = 2;
                
                    int vPixValue = img.getValue(vIndex);
            
                    if ((notValue && (vPixValue == pixelValue)) ||
                        (!notValue && (vPixValue != pixelValue))) {
                        
                        continue;
                    }
                    
                    processPair(uIndex, vIndex);
                
                    // inserting back at the top of the stack assures that the 
                    // search continues next from an associated point
                    stack.insert(vIndex);
                }
            }
        }        
    }

    protected void processPair(int uIdx, int vIdx) {
        
        //log.finest("processPair " + uSortedXIndex + ":" + vSortedXIndex);           
        
        int groupId;
        
        if ((pointToGroupIndex[uIdx] > -1) && (pointToGroupIndex[vIdx] == -1)) {
        
            groupId = pointToGroupIndex[uIdx];
            
            pointToGroupIndex[vIdx] = groupId;
            
            groupMembership[groupId].insert(vIdx);
            
        } else if ((pointToGroupIndex[vIdx] > -1) && (pointToGroupIndex[uIdx] == -1)) {
            
           groupId = pointToGroupIndex[vIdx];
            
            pointToGroupIndex[uIdx] = groupId;
            
            groupMembership[groupId].insert(uIdx);
            
        } else if ((pointToGroupIndex[uIdx] == -1) && (pointToGroupIndex[vIdx] == -1)) {
            
            checkAndExpandGroupMembershipArray();
            
            groupId = nGroups;
            
            pointToGroupIndex[uIdx] = groupId;
            
            pointToGroupIndex[vIdx] = groupId; 
            
            groupMembership[groupId].insert(uIdx);
            
            groupMembership[groupId].insert(vIdx);
     
            nGroups++;
            
        } else {
            
            groupId = -1;
            
            log.finest("not reading " + uIdx + ":" + vIdx );
        }
       
    }
            
    public SimpleLinkedListNode[] getGroupMembershipList() {
        return Arrays.copyOf(groupMembership, nGroups);
        //return groupMembership;
    }

    public int getNumberOfGroups() {
        return nGroups;
    }

    public int[] getPointToGroupIndexes() {
        return pointToGroupIndex;
    }
    
    /**
     * remove groups smaller than minimumNumberInCluster
     */
    protected void prune() {
        
        log.finest("number of groups before prune=" + nGroups);
        
        /*
         * [------] 0
         * [------] 1 <---- too few
         * [------] 2
         */
        // iterate backwards so can move items up without conflict with iterator
        for (int i = (nGroups - 1); i > -1; i--) {
            
            SimpleLinkedListNode group = groupMembership[i];
            
            int count = 0;
            SimpleLinkedListNode latest = group;
            while (latest != null) {
                count++;
                latest = latest.getNext();
            }
            
            log.finest("  group " + i + " has " + count 
                + " members before prune (min=" + minimumNumberInCluster + ")");
            
            if (count < minimumNumberInCluster) {
                
                // remove this group and move up all groups w/ index > i by one index
                for (int j = (i + 1); j < nGroups; j++) {
                    
                    int newGroupId = j - 1;
                    
                    groupMembership[newGroupId] = groupMembership[j];
                    
                    // update members in pointToGroupIndex
                    latest = groupMembership[j];
                    
                    while (latest != null) {
                        
                        int idx = latest.getKey();
                        
                        pointToGroupIndex[idx] = newGroupId;
                        
                        latest = latest.getNext();
                        
                    }
                }
                
                nGroups--;
            }
        }
        
        log.finest("number of groups after prune=" + nGroups);
    }

    protected void checkAndExpandGroupMembershipArray() {
        
        if (groupMembership == null) {
            throw new IllegalStateException("groupMembership cannot be null");
        }

        if (groupMembership.length < (nGroups + 1)) {
            int oldN = groupMembership.length;
            int n = (oldN == 0) ? 10 : (int) (1.5f * oldN);

            groupMembership = Arrays.copyOf(groupMembership, n);
            for (int k = oldN; k < n; k++) {
                groupMembership[k] = new SimpleLinkedListNode();
            }
        }
    }
    
    /**
     * estimate the amount of memory used by this class and its instance and class variables
     * @return the memory in bytes used by this class and its instance and class variables.
    */
    public long approximateMemoryUsed() {
        
        String arch = System.getProperty("sun.arch.data.model");

        boolean is32Bit = ((arch != null) && arch.equals("64")) ? false : true;

        int nbits = (is32Bit) ? 32 : 64;

        int overheadBytes = 16;

        int intBytes = (is32Bit) ? 4 : 8;
        int arrayBytes = 32/8;
        int refBytes = nbits/8;
        
        long sumBytes = 4*intBytes + 2*arrayBytes;
        if (groupMembership != null) {
            long nodeInBytes = SimpleLinkedListNode.approximateMemoryUsed();
            sumBytes += (groupMembership.length * nodeInBytes);
        }
        if (pointToGroupIndex != null) {
            sumBytes += (pointToGroupIndex.length * intBytes);
        }
        
        sumBytes += overheadBytes;
        
        long padding = (sumBytes % 8);
        
        sumBytes += padding;
        
        return sumBytes;
    }
    
    public int getNumberofGroupMembers(int groupId) {
        
        if (nGroups == 0) {
            return 0;
        }
        if (groupId > (nGroups - 1) || (groupId < 0)) {
            throw new IllegalArgumentException("groupId=" + groupId 
            + " is outside of range of nGroups=" + nGroups);
        }
        
        int[] indexes = getIndexes(groupId);
        
        return indexes.length;
    }
    
    public PairIntArray getXY(int groupId) {
        
        if (nGroups == 0) {
            return new PairIntArray();
        }
        if (groupId > (nGroups - 1) || (groupId < 0)) {
            throw new IllegalArgumentException("groupId=" + groupId 
            + " is outside of range of nGroups=" + nGroups);
        }
        
        int[] indexes = getIndexes(groupId);
        
        PairIntArray xy = new PairIntArray();
                
        for (int i = 0; i < indexes.length; i++) {
            
            int index = indexes[i];
            
            img.getXY(xy, index);
        }
        
        return xy;
    }
    
    public int[] getIndexes(int groupId) {
        
        if (nGroups == 0) {
            return new int[0];
        }
        if (groupId > (nGroups - 1) || (groupId < 0)) {
            throw new IllegalArgumentException("groupId=" + groupId 
            + " is outside of range of nGroups=" + nGroups);
        }
        
        int[] indexes = new int[10];
        
        int count = 0;
        
        SimpleLinkedListNode group = groupMembership[groupId];
        
        while (group != null) {
            
            int idx = group.getKey();
            
            indexes = expandIfNeeded(indexes, count + 1);
            
            indexes[count] = idx;
            
            group = group.getNext();
            
            count++;
        }
        
        return Arrays.copyOf(indexes, count);
    }
 
    protected int[] expandIfNeeded(int[] a, int nTotal) {
        if (nTotal > a.length) {
            int n = a.length + 10;
            if (nTotal > n) {
                n = nTotal;
            }
            return Arrays.copyOf(a, n);
        }
        return a;
    }
}
