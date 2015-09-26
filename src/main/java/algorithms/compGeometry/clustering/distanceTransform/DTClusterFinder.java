package algorithms.compGeometry.clustering.distanceTransform;

import algorithms.util.PairInt;
import java.util.Set;

/**
 * main class to cluster finder whose logic is based upon distance transform,
 * density threshold, and a signal to noise argument of ~ 3 as a factor.
 * 
 * @author nichole
 */
public class DTClusterFinder {
    
    private final Set<PairInt> points;
    private final int width;
    private final int height;
    
    private float critDens = Float.POSITIVE_INFINITY;
    
    private DTGroupFinder groupFinder = null;
    
    private enum STATE {
        INIT, HAVE_CLUSTER_DENSITY, HAVE_GROUPS
    }
    
    private STATE state = null;
    
    private boolean debug = false;
    
    public DTClusterFinder(Set<PairInt> points, int width, int height) {
        
        this.points = points;
        this.width = width;
        this.height = height;
        
        state = STATE.INIT;
    }
    
    public void setToDebug() {
        debug = true;
    }
    
    public void calculateCriticalDensity() {
        
        if (state.compareTo(STATE.HAVE_CLUSTER_DENSITY) > -1) {
            return;
        }
        
        DistanceTransform dtr = new DistanceTransform();
        int[][] dt = dtr.applyMeijsterEtAl(points, width, height);
        
        CriticalDensitySolver densSolver = new CriticalDensitySolver();
        
        if (debug) {
            densSolver.setToDebug();
        }
        
        this.critDens = densSolver.findCriticalDensity(dt, points.size(), width, height);               
    }
    
    public void setCriticalDensity(float dens) {
        
        if (state.compareTo(STATE.HAVE_CLUSTER_DENSITY) > -1) {
            throw new IllegalStateException("cluster density is already set");
        }
        
        this.critDens = dens;
    }
    
    public void findClusters() {
        
        if (state.compareTo(STATE.HAVE_GROUPS) < 1) {
            calculateCriticalDensity();
        } else if (state.compareTo(STATE.HAVE_CLUSTER_DENSITY) > -1) {
            return;
        }
        
        groupFinder = new DTGroupFinder();
        
        groupFinder.calculateGroups(critDens, points);
        
    }
    
    public int getNumberOfClusters() {
        
        if (groupFinder == null) {
            return 0;
        }
        
        return groupFinder.getNumberOfGroups();
    }
    
    public Set<PairInt> getCluster(int idx) {
        
        if (groupFinder == null) {
            throw new IllegalArgumentException(
                "findClusters was not successfully invoked");
        }
        
        if ((idx < 0) || (idx > (groupFinder.getNumberOfGroups() - 1))) {
            throw new IllegalArgumentException("idx is out of bounds");
        }
        
        return groupFinder.getGroup(idx);
    }
    
    public float getCriticalDensity() {
        return critDens;
    }
    
}
