package algorithms.imageProcessing;

import algorithms.compGeometry.HoughTransform;
import algorithms.compGeometry.PerimeterFinder2;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.TIntList;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

/**
 * a class to hold various methods for determining vanishing points
 * and to hold the resulting vanishing points.
 * 
 * @author nichole
 */
public class VanishingPoints {
    
    public Map<PairInt, Set<PairInt>> polarMap;
    
    /**
     * points and orientations to use for calculating vanishing lines and points.
     * NOTE that the points should probably be edge points or
     * the boundaries of segmentation contiguous labels as filtered edges.
     */
    public void find(List<Set<PairInt>> listOfContigousLabels) {
        
        // -- extract the boundaries of the sets
        // -- use hough transform to find the lines
        // -- find the parallel and overlapping lines
        // -- calculate intersection of non-parallel and non-overlapping lines.
        // -- use mode and average to combine results to
        //    local vanishing points
        // -- assert that geometry is consistent when possible
        // -- assert that have no more than 3 vanishing points
       
        List<PairIntArray> listOfBounds = new ArrayList<PairIntArray>();
        
        // TODO: extract ordered bounds
        
        HoughTransform ht = new HoughTransform();
        
        // make map w/ key = pairint(theta, radius), value = keypoints
        Map<PairInt, Set<PairInt>> thetaRadiusMap = ht.findLines(
            listOfBounds);
                
        polarMap = thetaRadiusMap;
    }
    
}
