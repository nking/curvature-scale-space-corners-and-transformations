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
    
    public abstract void findClusters(DoubleAxisIndexer indexer);
    
    public abstract void constructLogger();
    
    public void findGroups(DoubleAxisIndexer indexer) {
        
        constructLogger();
        
        initializeVariables(indexer);
        
        findClusters(indexer);
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
}
