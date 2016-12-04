package algorithms.imageProcessing;

import algorithms.compGeometry.HoughTransform;
import algorithms.compGeometry.PerimeterFinder2;
import algorithms.imageProcessing.matching.LineFinder;
import algorithms.imageProcessing.matching.LinesFinder;
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
        
    /**
     * points and orientations to use for calculating vanishing lines and points.
     * NOTE that the points should probably be edge points or
     * the boundaries of segmentation contiguous labels as filtered edges.
     */
    public void find(List<Set<PairInt>> listOfContigousLabels) {
        
        LinesFinder finder = new LinesFinder();
        
        finder.find(listOfContigousLabels);
    
        throw new UnsupportedOperationException("not yet implemented");
    }
    
}
