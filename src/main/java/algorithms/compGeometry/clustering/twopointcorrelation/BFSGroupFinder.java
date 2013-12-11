package algorithms.compGeometry.clustering.twopointcorrelation;

import java.util.logging.Logger;

/**
 * use a Breadth First Search algorithm to visit nodes within distance <
 * threshhold*threshholdFactor from a node to find the groups of nodes.
 * 
 * average runtime complexity is approx O(|E|), worst case runtime: O(|V| + |E|)
 * 
 * @author nichole
 */
public class BFSGroupFinder extends AbstractGroupFinder {

    public BFSGroupFinder(float threshhold, float threshholdFactor) {
        super(threshhold, threshholdFactor);
    }
    
    public void constructLogger() {
        this.log = Logger.getLogger(this.getClass().getName());
    }
    
    @Override
    public void findClusters(DoubleAxisIndexer indexer) {
        // TODO Auto-generated method stub

    }

}
