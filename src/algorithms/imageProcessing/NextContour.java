package algorithms.imageProcessing;

import java.util.ArrayList;
import java.util.List;
import java.util.TreeSet;

/**
 * A class to serve the purpose of an adjacency matrix (via rules to choose
 * the "neighbor contour") and the state of "un-visited" for the A* search
 * algorithm.
 * 
 * @author nichole
 */
public class NextContour {
    
    /**
     * a list of closed contours extracted from a digital image
     */
    protected final TreeSet<CurvatureScaleSpaceContour> origContours;
    
    /**
     * A list of indexes to origContours carrying contours from the curve, that
     * is, closed curve edge, given to the contour matcher.
     * 
     * Note that the lists are ordered by descending peak sigma.   Note also that the 
     * List<Integer> indexes are referred to as contour indexes.
     */
    protected final TreeSet<CurvatureScaleSpaceContour> remainingContours;
    
    protected final List<CurvatureScaleSpaceContour> matchedContours1;
    protected final List<CurvatureScaleSpaceContour> matchedContours2;
    
    protected int matchedEdgeNumber1 = -1;
    protected int matchedEdgeNumber2 = -1;
    
    public NextContour(final List<CurvatureScaleSpaceContour> contours,
        List<CurvatureScaleSpaceContour> alreadyVisited) {
                
        origContours = new TreeSet<CurvatureScaleSpaceContour>(
            new DescendingSigmaComparator());
                        
        // populate with contours that haven't been visited
        remainingContours = new TreeSet<CurvatureScaleSpaceContour>(
            new DescendingSigmaComparator());
        
        for (int i = 0; i < contours.size(); i++) {
            
            CurvatureScaleSpaceContour contour = contours.get(i);
            
            origContours.add(contour);
            
            if (alreadyVisited.contains(contour)) {
                continue;
            }
            
            remainingContours.add(contour);
        }

        matchedContours1 = new ArrayList<CurvatureScaleSpaceContour>();
        
        matchedContours2 = new ArrayList<CurvatureScaleSpaceContour>();
    }
    
    public void addMatchedContours(CurvatureScaleSpaceContour contour1,
        CurvatureScaleSpaceContour contour2) {
        
        markAsVisited(contour1);
        
        matchedContours1.add(contour1);
        matchedContours2.add(contour2);

        if (matchedEdgeNumber1 == -1) {
            matchedEdgeNumber1 = contour1.getEdgeNumber();
        } else {
            assert(matchedEdgeNumber1 == contour1.getEdgeNumber());
        }
        
        if (matchedEdgeNumber2 == -1) {
            matchedEdgeNumber2 = contour2.getEdgeNumber();
        } else {
            assert(matchedEdgeNumber2 == contour2.getEdgeNumber());
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
    public CurvatureScaleSpaceContour findTallestContourWithinScaleSpace() {
        
        if (remainingContours.isEmpty()) {
            return null;
        }
        
        CurvatureScaleSpaceContour contour = remainingContours.first();
        
        boolean removed = remainingContours.remove(contour);
                
        assert(removed == true);
        
        return contour;
    }
    
    public void markAsVisited(CurvatureScaleSpaceContour contour) {
        
        remainingContours.remove(contour);
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
    public CurvatureScaleSpaceContour findTheNextSmallestUnvisitedSibling(
        CurvatureScaleSpaceContour target) { 

        if (target == null) {
            return null;
        }
        
        CurvatureScaleSpaceContour nextLower = origContours.higher(target);
        
        while ((nextLower != null) && !remainingContours.contains(nextLower)) {
            nextLower = origContours.higher(nextLower);
        }
                
        if (nextLower != null) {
            
            boolean removed = remainingContours.remove(nextLower);
            
            assert(removed);
            
            return nextLower;
        }
        
        return null;
    }
    
    public List<CurvatureScaleSpaceContour> getMatchedContours1() {
        return matchedContours1;
    }
    
    public List<CurvatureScaleSpaceContour> getMatchedContours2() {
        return matchedContours2;
    }
}
