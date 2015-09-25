package algorithms.compGeometry.clustering.distanceTransform;

import algorithms.util.PairInt;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;

/**
 *
 * @author nichole
 */
public class DTGroupFinder {
    /**
     * an array to hold each group as an item.  each item contains a key which is an index
     * to arrays indexer.x, indexer.y and this.pointToGroupIndex
     */
    protected List<Set<PairInt> > groupMembership = new ArrayList<Set<PairInt> >();
    
    protected Logger log = null;
    
    protected Set<PairInt> visited = new HashSet<PairInt>();
    
     /*
     * map w/ key holding indexes for a point to the group it belongs to.
     * note that the point index is relative to indexer.x and indexer.y
     */
    //protected int[] pointToGroupIndex = null;
    protected Map<PairInt, Integer> pointToGroupMap = new
        HashMap<PairInt, Integer >();
    
    protected int minimumNumberInCluster = 3;
    
    protected boolean notValue = false;
    
    protected boolean debug = false;
    
    protected float threshholdFactor = 3.0f;
        
    public DTGroupFinder() {                
        this.log = Logger.getLogger(this.getClass().getName());
    }
        
    public void setMinimumNumberInCluster(int n) {
        this.minimumNumberInCluster = n;
    }
    
    public void setDebug(boolean setDebugToTrue) {
        this.debug = setDebugToTrue;
    }
    
    public void setThreshholdFactor(float factor) {
        this.threshholdFactor = factor;
    }

    /**
     * given the critical density and having the threshold factor, find the
     * groups of points within a critical distance of one another.
     * runtime complexity is O(N_points * lg2(N_points)).
     * @param criticalDensity
     * @param points
     * @param width
     * @param height 
     */
    void calculateGroups(float criticalDensity, Set<PairInt> points, int width,
        int height) {
        
        float thrsh = criticalDensity * threshholdFactor;
        
        findGroups(thrsh, points, width, height);
        
        prune(); 
    }
    
    private void findGroups(float thrsh, Set<PairInt> points,
        int width, int height) {
        
        if (points.isEmpty()) {
            return;
        }
        
        // association of 2 points for separation <= critSeparation
        float critSep = 2.f/thrsh;
        
        if (critSep < 1) {
            // each point is a group
            setEachPointAsAGroup();
            return;
        }
        
        //TODO: consider data structures with point location as part of their
        // structure... spatial indexing, RTrees...
        /*
        need to sort the points by x then y to be able to make small searches
        around a point as traverse the stack.
        The runtime complexity is O(N*lg2(N)).
        Note that if the dataset were sparse, could assume only need to sort
        on x and use the O(N) counting sort (if N=100 for example, and the
        dataset were sparse and there are fewer than 7 points with the same y
        for each x, there's no need to sort on y).
        */
        
       
        
    }
    
    int getNumberOfGroups() {
        throw new UnsupportedOperationException("not yet implemented");
    }
    
    Set<PairInt> getGroup(int idx) {
        throw new UnsupportedOperationException("not yet implemented");
    }

    private void prune() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    private void setEachPointAsAGroup() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

}
