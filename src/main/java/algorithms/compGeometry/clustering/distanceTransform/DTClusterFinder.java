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
    
    public DTClusterFinder(Set<PairInt> points, int width, int height) {
        
        this.points = points;
        this.width = width;
        this.height = height;
    }
    
    void calculateCriticalDensity() {
        
        DistanceTransform dtr = new DistanceTransform();
        int[][] dt = dtr.applyMeijsterEtAl(points, width, height);
        
        CriticalDensitySolver densSolver = new CriticalDensitySolver();
        
        this.critDens = densSolver.findCriticalDensity(dt, points.size(), width, height);               
    }
    
    void setCriticalDensity(float dens) {
        this.critDens = dens;
    }
    
    void findClusters() {
        
        groupFinder = new DTGroupFinder();
        
        groupFinder.calculateGroups(critDens, points, width, height);
        
    }
    
    int getNumberOfClusters() {
        
        if (groupFinder == null) {
            return 0;
        }
        
        return groupFinder.getNumberOfGroups();
    }
    
    Set<PairInt> getCluster(int idx) {
        
        if (groupFinder == null) {
            throw new IllegalArgumentException(
                "findClusters was not successfully invoked");
        }
        
        if ((idx < 0) || (idx > (groupFinder.getNumberOfGroups() - 1))) {
            throw new IllegalArgumentException("idx is out of bounds");
        }
        
        return groupFinder.getGroup(idx);
    }
    
    float getCriticalDensity() {
        return critDens;
    }
    
}
