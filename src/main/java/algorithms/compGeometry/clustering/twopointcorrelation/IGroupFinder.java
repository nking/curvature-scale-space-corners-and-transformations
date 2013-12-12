package algorithms.compGeometry.clustering.twopointcorrelation;

public interface IGroupFinder {
    
    public void findGroups(DoubleAxisIndexer indexer);
    
    public SimpleLinkedListNode[] getGroupMembershipList();
    
    public int getNumberOfGroups();
    
    public int[] getPointToGroupIndexes();
    
    public void setMinimumNumberInCluster(int n);
    
    public float[] getX(int groupId, DoubleAxisIndexer indexer);
    
    public float[] getY(int groupId, DoubleAxisIndexer indexer);
    
    public int[] getIndexes(int groupId);
    
    public void printMembership(DoubleAxisIndexer indexer);
    
    public long approximateMemoryUsed();
    
    public void setDebug(boolean setDebugToTrue);
}
