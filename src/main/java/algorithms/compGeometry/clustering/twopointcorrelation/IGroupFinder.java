package algorithms.compGeometry.clustering.twopointcorrelation;

public interface IGroupFinder {
    
    public void findGroups(DoubleAxisIndexer indexer);
    
    public SimpleLinkedListNode[] getGroupMembershipList();
    
    public int getNumberOfGroups();
    
    public int[] getPointToGroupIndexes();
    
    public void setMinimumNumberInCluster(int n);
}
