package algorithms.compGeometry.clustering.twopointcorrelation;

public abstract class AbstractGroupFinder implements IGroupFinder {

    protected SimpleLinkedListNode[] groupMembership = null;
    
    protected int[] pointToGroupIndex = null;
        
    protected int nGroups = 0;
    
    protected final float threshhold;
    
    protected final float threshholdFactor;
    
    protected int minimumNumberInCluster = 3;
    
    public AbstractGroupFinder(float threshhold, float threshholdFactor) {
        this.threshhold = threshhold;
        this.threshholdFactor = threshholdFactor;
    }
    
    public void setMinimumNumberInCluster(int n) {
        this.minimumNumberInCluster = n;
    }
    
    public abstract void findGroups(DoubleAxisIndexer indexer);

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

}
