package algorithms.imageProcessing;

import java.util.ArrayList;
import java.util.List;
import java.util.TreeSet;

/**
 * A class to serve the purpose of an adjacency matrix (via rules to choose
 * the "neighbor corner") and the state of "un-visited" for the search
 * algorithm.
 * 
 * @author nichole
 */
public class NextCorner<T extends CornerRegion> {
    
    /**
     * a list of corners extracted from a closed curve in a digital image
     */
    protected final TreeSet<T> origCorners;
    
    /**
     * A list of indexes to origCorners carrying corner regions from the curve, 
     * that is, closed curve edge, given to the corner matcher.
     * 
     * Note that the lists are ordered by descending peak sigma.   Note also 
     * that the List<Integer> indexes are referred to as corner indexes.
     */
    protected final TreeSet<T> remainingCorners;
    
    private final List<Integer> matchedCornerIndexes1;
    private final List<Integer> matchedCornerIndexes2;
    
    protected final List<FeatureComparisonStat> matchedStats;
    protected final List<Double> matchedDist;
    
    protected int matchedEdgeNumber1 = -1;
    protected int matchedEdgeNumber2 = -1;
    
    public NextCorner(final List<T> corners) {
                
        origCorners = new TreeSet<T>(new DescendingKComparator());
                        
        // populate with contours that haven't been visited
        remainingCorners = new TreeSet<T>(
            new DescendingKComparator());
        
        matchedStats = new ArrayList<FeatureComparisonStat>();
        
        matchedDist = new ArrayList<Double>();
        
        for (int i = 0; i < corners.size(); i++) {
            
            T corner = corners.get(i);
            
            origCorners.add(corner);
            
            remainingCorners.add(corner);
        }

        matchedCornerIndexes1 = new ArrayList<Integer>();
        matchedCornerIndexes2 = new ArrayList<Integer>();
        
    }
    
    public void addMatchedCorners(T corner1, T corner2, Integer cornerIndex1, 
        Integer cornerIndex2, FeatureComparisonStat stat) {
        
        double dist = Double.POSITIVE_INFINITY;
        
        addMatchedCorners(corner1, corner2, cornerIndex1, cornerIndex2, stat,
            dist);
    }
    
    public void addMatchedCorners(T corner1, T corner2, Integer cornerIndex1, 
        Integer cornerIndex2, FeatureComparisonStat stat, double dist) {
        
        markAsVisited(corner1);
        
        if (corner2 == null) {
            return;
        }
        
        matchedCornerIndexes1.add(cornerIndex1);
        matchedCornerIndexes2.add(cornerIndex2);
        matchedStats.add(stat);
        matchedDist.add(Double.valueOf(dist));
        
        if (matchedEdgeNumber1 == -1) {
            matchedEdgeNumber1 = corner1.getEdgeIdx();
        } else {
            assert(matchedEdgeNumber1 == corner1.getEdgeIdx());
        }
        
        if (matchedEdgeNumber2 == -1) {
            matchedEdgeNumber2 = corner2.getEdgeIdx();
        } else {
            assert(matchedEdgeNumber2 == corner2.getEdgeIdx());
        }    
    }
    
    /**
     * get cost as normalized sum of SSDs times normalized sums of distances
     * @param maxSSD
     * @param maxDistance
     * @return
     */
    public double getNormalizedCost(double maxSSD, double maxDistance) {

        double sumSSD = 0;
        double sumDist = 0;
        
        for (int i = 0; i < matchedStats.size(); ++i) {
            sumSSD += matchedStats.get(i).getSumIntensitySqDiff();
            sumDist += matchedDist.get(i).doubleValue();
        }
        sumSSD /= (double)matchedStats.size();
        sumDist /= (double)matchedStats.size();
        
        sumSSD /= maxSSD;
        sumDist /= maxDistance;

        double normalizedCost = sumSSD * sumDist;
        
        return normalizedCost;
    }
    
    /**
     * find the largest remaining corner within the list of unvisited.
     * side-effect of removing the returned corner index from the look-up data
     * structures.
     * 
     * @return 
     */
    public T findStrongestRemainingCorner() {
        
        if (remainingCorners.isEmpty()) {
            return null;
        }
        
        T corner = remainingCorners.first();
        
        boolean removed = remainingCorners.remove(corner);
                
        assert(removed == true);
        
        return corner;
    }
    
    public void markAsVisited(CornerRegion corner) {        
        remainingCorners.remove(corner);
    }
    
    public int getMatchedEdgeNumber2() {
        return matchedEdgeNumber2;
    }
    public int getMatchedEdgeNumber1() {
        return matchedEdgeNumber1;
    }
    
    /**
     * find the next smallest sigma peak out of all contours. Note that
     * this method has the side effect of removing the returned contour
     * from the internal look-up data structures.
     * 
     * @param target the contour or which to find the next smallest contour
     * 
     * @return 
     */
    public T findTheNextSmallestUnvisitedSibling(T target) { 

        if (target == null) {
            return null;
        }
        
        T nextLower = origCorners.higher(target);
        
        while ((nextLower != null) && !remainingCorners.contains(nextLower)) {
            nextLower = origCorners.higher(nextLower);
        }
                
        if (nextLower != null) {
            
            boolean removed = remainingCorners.remove(nextLower);
            
            assert(removed);
            
            return nextLower;
        }
        
        return null;
    }
    
    public List<Integer> getMatchedCornerIndexes1() {
        return matchedCornerIndexes1;
    }
    
    public List<Integer> getMatchedCornerIndexes2() {
        return matchedCornerIndexes2;
    }
    
    public List<FeatureComparisonStat> getMatchedFeatureComparisonStats() {
        return matchedStats;
    }
    
    public List<Double> getMatchedDistances() {
        return matchedDist;
    }
}
