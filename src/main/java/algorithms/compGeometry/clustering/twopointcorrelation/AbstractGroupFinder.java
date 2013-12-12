package algorithms.compGeometry.clustering.twopointcorrelation;

import java.util.Arrays;
import java.util.logging.Logger;

public abstract class AbstractGroupFinder implements IGroupFinder {

    /**
     * an array to hold each group as an item.  each item contains keys which hold the <xy>Point index.
     */
    protected SimpleLinkedListNode[] groupMembership = null;

    protected int nGroups = 0;
    
    protected Logger log = null;
    
    /*
     * array holding indexes for a point to the group it belongs to.
     * note that the point index is relative to indexer.getXSortedByY
     */
    protected int[] pointToGroupIndex = null;
            
    protected final float threshhold;
    
    protected final float threshholdFactor;
    
    protected int minimumNumberInCluster = 3;
    
    protected boolean debug = false;
    
    public AbstractGroupFinder(float threshhold, float threshholdFactor) {
        this.threshhold = threshhold;
        this.threshholdFactor = threshholdFactor;
    }
    
    public void setMinimumNumberInCluster(int n) {
        this.minimumNumberInCluster = n;
    }
    
    public void setDebug(boolean setDebugToTrue) {
        this.debug = setDebugToTrue;
    }
    
    protected void initializeVariables(DoubleAxisIndexer indexer) {

        pointToGroupIndex = new int[indexer.getNumberOfPoints()];
        
        Arrays.fill(pointToGroupIndex, -1);

        groupMembership = new SimpleLinkedListNode[10];
        
        for (int i = 0; i < groupMembership.length; i++) {
            
            groupMembership[i] = new SimpleLinkedListNode();
        }
    }
    
    protected abstract void findClusters(DoubleAxisIndexer indexer);
    
    public abstract void constructLogger();
    
    public void findGroups(DoubleAxisIndexer indexer) {
                
        constructLogger();
        
        initializeVariables(indexer);
        
        findClusters(indexer);
        
        prune();        
    }

    @Override
    public SimpleLinkedListNode[] getGroupMembershipList() {
        return groupMembership;
    }

    @Override
    public int getNumberOfGroups() {
        return nGroups;
    }

    @Override
    public int[] getPointToGroupIndexes() {
        return pointToGroupIndex;
    }
    
    /**
     * remove groups smaller than minimumNumberInCluster
     */
    protected void prune() {
        
        log.info("number of groups before prune=" + nGroups);
        
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
                latest = latest.next;
            }
            
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
                        
                        latest = latest.next;
                    }
                }
                
                nGroups--;
                
            }
        }
        
        log.info("number of groups after prune=" + nGroups);
    }

    protected void checkAndExpandGroupMembershipArray() {

        if (groupMembership.length < (nGroups + 1)) {
            int oldN = groupMembership.length;
            int n = (int) (1.5f * oldN);
            if (n < (oldN + 1)) {
                n = oldN + 1;
            }

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
    
    public float[] getX(int groupId, DoubleAxisIndexer indexer) {
        
        if (nGroups == 0) {
            return new float[0];
        }
        if (groupId > (nGroups - 1) || (groupId < 0)) {
            throw new IllegalArgumentException("groupId=" + groupId + " is outside of range of nGroups=" + nGroups);
        }
        
        int[] indexes = getIndexes(groupId);
        
        float[] x = new float[indexes.length];
        
        for (int i = 0; i < indexes.length; i++) {
            
            x[i] = indexer.getX()[ indexes[i] ];
        }
        
        return x;
    }
    
    public float[] getY(int groupId, DoubleAxisIndexer indexer) {
        
        if (nGroups == 0) {
            return new float[0];
        }
        if (groupId > (nGroups - 1) || (groupId < 0)) {
            throw new IllegalArgumentException("groupId=" + groupId + " is outside of range of nGroups=" + nGroups);
        }
        
        int[] indexes = getIndexes(groupId);
        
        float[] y = new float[indexes.length];
        
        for (int i = 0; i < indexes.length; i++) {
            
            y[i] = indexer.getY()[ indexes[i] ];
        }
        
        return y;
    }
    
    public int[] getIndexes(int groupId) {
        
        if (nGroups == 0) {
            return new int[0];
        }
        if (groupId > (nGroups - 1) || (groupId < 0)) {
            throw new IllegalArgumentException("groupId=" + groupId + " is outside of range of nGroups=" + nGroups);
        }
        
        int[] indexes = new int[10];
        
        int count = 0;
        
        SimpleLinkedListNode group = groupMembership[groupId];
        
        while (group != null) {
            
            int idx = group.getKey();
            
            indexes = expandIfNeeded(indexes, count + 1);
            
            indexes[count] = idx;
            
            group = group.next;
            
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
    
    public void printMembership(DoubleAxisIndexer indexer) {
        
        System.out.println(nGroups + " Groups:");
        
        for (int i = 0; i < nGroups; i++) {
            
            System.out.println("  group " + i);
            
            SimpleLinkedListNode group = groupMembership[i];
        
            while (group != null) {
                
                int idx = group.getKey();
                
                float x = indexer.getX()[idx];
                
                float y = indexer.getY()[idx];
                
                System.out.println("    (" + x + "," + y + ")");
                
                group = group.next;
            }
        }
    }

}
