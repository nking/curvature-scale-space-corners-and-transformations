package algorithms.compGeometry.clustering.twopointcorrelation;

/**
 * interface for classes that find groups of points based upon the background threshhold density.
 * 
 * @author nichole
 *
 */
public interface IGroupFinder {
    
    /**
     *
     * @param indexer
     */
    public void findGroups(AxisIndexer indexer);
    
    /**
     *
     * @return
     */
    public SimpleLinkedListNode[] getGroupMembershipList();
    
    /**
     *
     * @return
     */
    public int getNumberOfGroups();
    
    /**
     *
     * @return
     */
    public int[] getPointToGroupIndexes();
    
    /**
     *
     * @param n
     */
    public void setMinimumNumberInCluster(int n);
    
    /**
     *
     * @param groupId
     * @param indexer
     * @return
     */
    public float[] getX(int groupId, AxisIndexer indexer);
    
    /**
     *
     * @param groupId
     * @param indexer
     * @return
     */
    public float[] getY(int groupId, AxisIndexer indexer);
    
    /**
     *
     * @param groupId
     * @return
     */
    public int[] getIndexes(int groupId);
    
    /**
     *
     * @param indexer
     */
    public void printMembership(AxisIndexer indexer);
    
    /**
     *
     * @return
     */
    public long approximateMemoryUsed();
    
    /**
     *
     * @param setDebugToTrue
     */
    public void setDebug(boolean setDebugToTrue);
}
