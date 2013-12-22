package algorithms.compGeometry.clustering.twopointcorrelation;

public interface IGroupFinder {
    
    public void findGroups(AxisIndexer indexer);
    
    public SimpleLinkedListNode[] getGroupMembershipList();
    
    public int getNumberOfGroups();
    
    public int[] getPointToGroupIndexes();
    
    public void setMinimumNumberInCluster(int n);
    
    public float[] getX(int groupId, AxisIndexer indexer);
    
    public float[] getY(int groupId, AxisIndexer indexer);
    
    public int[] getIndexes(int groupId);
    
    public void printMembership(AxisIndexer indexer);
    
    public long approximateMemoryUsed();
    
    public void setDebug(boolean setDebugToTrue);
}
