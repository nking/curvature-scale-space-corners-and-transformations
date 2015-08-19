package algorithms.compGeometry;

import algorithms.util.PairInt;
import java.util.HashMap;
import java.util.Map;

/**
 * Class to store backtracking nodes used by PerimeterFinder's
 * getOrderedBorderPixels algorithm and return a re-usable sublist when
 * possible.  The class is specific to that method and assumes an interaction
 * distance is (diffX lte 1 and diffY lte 1).
 * 
 * @author nichole
 */
public class PerimeterFinderMemo {
    
    private int currentBaseIndex = -1;
    
    private int currentBackTrackLevel = -1;
    
    /*
    a method is given a PathStep node and it's index in the invoker's path list.  
    If currentBaseIndex here is not set, set it and currentBackTrackLevel to 
    the index.
    Else baseIndex is set, 
       and if index < currentBackTrackLevel, update currentBackTrackLevel to index
       else if index >= currentBackTrackLevel, set currentBaseIndex and 
       currentBackTrackLevel to index.
    
    the PathStep is added to Map<baseIndex, List<PathStep>> backTrackNodes
                             Map<baseIndex, List<Index>> backTrackLevels
    
    
    another method can return a sublist of nodes which start with the requested
    key and extend as far up (to smaller indexes) the list as possible where
    each node from key to the smallest possible is tested that it is beyond
    an interaction distance of all nodes it is given in a argument to test
    against.  the interaction distance is defined as being adjacent to a node.
    in other words, if the backtrack branch sublist couldn't have been affected
    by the given nodes, a copy of the sublist is returned for use by invoker.
    
    Note that finding a key in backTrackNodes should be made faster by
    creating a lookup structure:
        Map<baseIndex, Map<PairInt, Integer>> where PairInt is the node's col, row
        and the Integer values are the indexes in backTrackNodes where the node
        is.
    
    */
    
}
