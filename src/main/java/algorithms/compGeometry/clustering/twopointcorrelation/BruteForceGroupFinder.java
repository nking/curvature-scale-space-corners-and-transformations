package algorithms.compGeometry.clustering.twopointcorrelation;

/**
 * iterate over every pair to find those with distances less than threshhold*threshholdFactor
 * from one another.
 * 
 * runtime complexity O(N^2)
 * 
 * @author nichole
 */
public class BruteForceGroupFinder extends AbstractGroupFinder {

    public BruteForceGroupFinder(float threshhold, float threshholdFactor) {
        super(threshhold, threshholdFactor);
    }
    
    @Override
    public void findGroups(DoubleAxisIndexer indexer) {
        // TODO Auto-generated method stub

    }
}
