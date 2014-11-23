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
    protected final TreeSet<CurvatureScaleSpaceContour> unvisitedContours;
    
    protected final List<CurvatureScaleSpaceContour> matchedContours1;
    protected final List<CurvatureScaleSpaceContour> matchedContours2;
    
    public NextContour(final List<CurvatureScaleSpaceContour> contours,
        final List<CurvatureScaleSpaceContour> alreadyVisited) {
        
        unvisitedContours = new TreeSet<CurvatureScaleSpaceContour>(
            new DescendingSigmaComparator());
        
        unvisitedContours.addAll(contours);
        
        matchedContours1 = new ArrayList<CurvatureScaleSpaceContour>();
        
        matchedContours2 = new ArrayList<CurvatureScaleSpaceContour>();
    }
    
    /**
     * find the largest sigma peak within the remaining un-searched contours.
     * Note that the method has the side-effect of marking the contour as
     * visited by removing it from the internal look-up data structures.
     * 
     * @return 
     */
    public CurvatureScaleSpaceContour findTallestContourWithinAScaleSpace() {
        
        // returns null if empty
        CurvatureScaleSpaceContour contour = unvisitedContours.pollFirst();
                                
        return contour;
    }
    
    public void unVisitMatchedContours(int idx) {
        
        CurvatureScaleSpaceContour cm1 = matchedContours1.get(idx);
        CurvatureScaleSpaceContour cm2 = matchedContours2.get(idx);
        matchedContours1.remove(cm1);
        matchedContours2.remove(cm2);
                
        unvisitedContours.add(cm1);
    }
    
    public void addMatchedContours(CurvatureScaleSpaceContour contour1,
        CurvatureScaleSpaceContour contour2) {
        
        matchedContours1.add(contour1);
        matchedContours2.add(contour2);
    }
    
    /**
     * find the next smallest sigma peak out of all contours.  the sigma is
     * found by looking up the contour that target references.  Note that
     * this method has the side effect of removing the returned contour
     * from the internal look-up data structures.
     * 
     * @param c contour for which to find the next lowest item in 
     * unvisitedContours.
     * 
     * @return 
     */
    public CurvatureScaleSpaceContour findTheNextSmallestUnvisitedSibling(
        CurvatureScaleSpaceContour c) { 

        CurvatureScaleSpaceContour contour = unvisitedContours.floor(c);

        return contour;
    }
  
    public List<CurvatureScaleSpaceContour> getMatchedContours1() {
        return matchedContours1;
    }
    
    public List<CurvatureScaleSpaceContour> getMatchedContours2() {
        return matchedContours2;
    }
}
