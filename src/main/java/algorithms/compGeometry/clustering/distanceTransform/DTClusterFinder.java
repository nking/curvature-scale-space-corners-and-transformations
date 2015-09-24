package algorithms.compGeometry.clustering.distanceTransform;

import algorithms.util.PairInt;
import java.util.Set;

/**
 * main class to cluster finder based upon distance transform.
 * 
 * @author nichole
 */
public class DTClusterFinder {
    
    private final int[][] xy;
    private final Set<PairInt> points;
    private final int width;
    private final int height;
    
    private float critDens = Float.POSITIVE_INFINITY;
    
    public DTClusterFinder(Set<PairInt> points, int width, int height) {
        this.xy = new int[width][];
        for (int i = 0; i < width; ++i) {
            xy[i] = new int[height];
        }
        this.points = points;
        this.width = width;
        this.height = height;
    }
    
    void calculateCriticalDensity() {
        throw new UnsupportedOperationException("not yet implemented");
    }
    
    void setCriticalDensity(float dens) {
        this.critDens = dens;
    }
    
    void findClusters() {
        throw new UnsupportedOperationException("not yet implemented");
    }
    
    int getNumberOfClusters() {
        throw new UnsupportedOperationException("not yet implemented");
    }
    
    Set<PairInt> getCluster(int idx) {
        throw new UnsupportedOperationException("not yet implemented");
    }
    
    float getCriticalDensity() {
        return critDens;
    }
    
}
