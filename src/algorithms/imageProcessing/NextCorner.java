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
public class NextCorner {
    
    /**
     * a list of corners extracted from a closed curve in a digital image
     */
    protected final TreeSet<CornerRegion> origCorners;
    
    /**
     * A list of indexes to origCorners carrying corner regions from the curve, 
     * that is, closed curve edge, given to the corner matcher.
     * 
     * Note that the lists are ordered by descending peak sigma.   Note also that the 
     * List<Integer> indexes are referred to as corner indexes.
     */
    protected final TreeSet<CornerRegion> remainingCorners;
    
    protected final List<CornerRegion> matchedCorners1;
    protected final List<CornerRegion> matchedCorners2;
    
    protected int matchedEdgeNumber1 = -1;
    protected int matchedEdgeNumber2 = -1;
    
    public NextCorner(final List<CornerRegion> corners,
        List<CornerRegion> alreadyVisited) {
                
        origCorners = new TreeSet<CornerRegion>(
            new DescendingKComparator());
                        
        // populate with contours that haven't been visited
        remainingCorners = new TreeSet<CornerRegion>(
            new DescendingKComparator());
        
        for (int i = 0; i < corners.size(); i++) {
            
            CornerRegion corner = corners.get(i);
            
            origCorners.add(corner);
            
            if (alreadyVisited.contains(corner)) {
                continue;
            }
            
            remainingCorners.add(corner);
        }

        matchedCorners1 = new ArrayList<CornerRegion>();
        
        matchedCorners2 = new ArrayList<CornerRegion>();
    }
    
    public void addMatchedCorners(CornerRegion corner1, CornerRegion corner2) {
        
        markAsVisited(corner1);
        
        if (corner2 == null) {
            return;
        }
        
        matchedCorners1.add(corner1);
        matchedCorners2.add(corner2);

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
     * find the largest sigma peak within the remaining un-searched contours.  
     * Note that the method has the 
     * side-effect of removing the returned contour from the look-up data
     * structures.
     * 
     * @return 
     */
    public CornerRegion findStrongestRemainingCorner() {
        
        if (remainingCorners.isEmpty()) {
            return null;
        }
        
        CornerRegion corner = remainingCorners.first();
        
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
    public CornerRegion findTheNextSmallestUnvisitedSibling(CornerRegion target) { 

        if (target == null) {
            return null;
        }
        
        CornerRegion nextLower = origCorners.higher(target);
        
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
    
    public List<CornerRegion> getMatchedCorners1() {
        return matchedCorners1;
    }
    
    public List<CornerRegion> getMatchedCorners2() {
        return matchedCorners2;
    }
}
