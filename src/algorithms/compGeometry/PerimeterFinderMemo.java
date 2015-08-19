package algorithms.compGeometry;

import algorithms.compGeometry.PerimeterFinder.PathStep;
import algorithms.util.PairInt;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

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
    
    /**
     * key = baseIndex (which is the back track level at which highest backtracking begins
     *       for a given list of nodes); 
     * value = list of nodes added for each backtracking for this baseIndex.
     */
    private Map<Integer, List<PathStep>> backTrackNodes = new HashMap<Integer, List<PathStep>>();
                                 
    /**
     * key = baseIndex (which is the back track level at which highest backtracking begins
     *       for a given list of nodes); 
     * value = the position of the node in the invoker's path list (it's the 
     *       back track level for it).
     */
    private Map<Integer, List<Integer> > backTrackLevels = new HashMap<Integer, List<Integer>>();
           
    /**
     * key = baseIndex (which is the back track level at which highest backtracking begins
     *       for a given list of nodes); 
     * value = map with key = PairInt of col, row of the node and value being 
     *         the position of the node in the backTrackNodes's list.
     */
    private Map<Integer, Map<PairInt, Integer>> backTrackLookupMap =
        new HashMap<Integer, Map<PairInt, Integer>>();
    
    public void storeNode(int pathIndex, PathStep node) {
        
        if (pathIndex < 0) {
            throw new IllegalStateException("pathIndex must be a positive number");
        }
        if (node == null) {
            throw new IllegalStateException("node cannot be null");
        }
        
        if (currentBaseIndex == -1) {
            currentBaseIndex = pathIndex;
            currentBackTrackLevel = pathIndex;
        } else {
            if (pathIndex < currentBackTrackLevel) {
                currentBackTrackLevel = pathIndex;
            } else {
                if (backTrackNodes.get(Integer.valueOf(pathIndex)) != null) {
                    //not expecting this to happen
                    throw new IllegalStateException(
                    "Error in algorithm: not expecting backtrack at an existing level at a later time");
                }
                currentBaseIndex = pathIndex;
                currentBackTrackLevel = pathIndex;
            }
        }
        
        Integer key0 = Integer.valueOf(currentBaseIndex);
        
        List<PathStep> backTrackNodesList = backTrackNodes.get(key0);
        if (backTrackNodesList == null) {
            backTrackNodesList = new ArrayList<PathStep>();
            backTrackNodes.put(key0, backTrackNodesList);
        }
        int idx = backTrackNodesList.size();
        backTrackNodesList.add(node);
        
        List<Integer> backTrackLevelList = backTrackLevels.get(key0);
        if (backTrackLevelList == null) {
            backTrackLevelList = new ArrayList<Integer>();
            backTrackLevels.put(key0, backTrackLevelList);
        }
        backTrackLevelList.add(Integer.valueOf(currentBackTrackLevel));
        
        Map<PairInt, Integer> backTrackLookupIndexes = backTrackLookupMap.get(key0);
        if (backTrackLookupIndexes == null) {
            backTrackLookupIndexes = new HashMap<PairInt, Integer>();
            backTrackLookupMap.put(key0, backTrackLookupIndexes);
        }
        backTrackLookupIndexes.put(node.coords, Integer.valueOf(idx));
        
    }
    
    /**
     * given the current node, the node just removed, and the next node,
     * look for a stored sublist starting at next node which can be appended
     * to the current path after next node.  The removed node and next node
     * are used as a range for which to look to make sure that they are not
     * adjacent to the sublist (and hence would affect their future state).
     * @param currentNode
     * @param removedNode
     * @param nextNode
     * @return sublist that can be appended to path after nextNode, else null
     * if a non-interacting sublist was not found.
     */
    public List<PathStep> getMemoizedSublist(PathStep currentNode, 
        PathStep removedNode, PathStep nextNode) {
        
        if (currentBaseIndex == -1) {
            return null;
        }
        
        /*
        if nextNode has been memoized, then
        -- the "interaction test nodes" are then removedNode through nextNode(excl)
           in the list from backTrackNodes 
           where removedNode will be at a higher index near the end of the list
           and nextNode's prev will be at lower index than removedNode
           at (nextNode's index + 1)
        -- test for each node in the list from backTrackNodes starting
           at (nextNode's index - 1) to index -1 whether they are adjacent to
           the "interaction test nodes" adding to output list until an adjacent
           is found or until beginning of list is reached.
        */
        
        Integer key0 = Integer.valueOf(currentBaseIndex);
        
        Map<PairInt, Integer> backTrackLookupIndexes = backTrackLookupMap.get(key0);
        
        Integer nextNodeBTLIndex = backTrackLookupIndexes.get(nextNode.coords);
        
        if (nextNodeBTLIndex == null) {
            return null;
        }
        
        Integer rmNodeBTLIndex = backTrackLookupIndexes.get(removedNode.coords);
        
        if (rmNodeBTLIndex == null) {
            throw new IllegalStateException("removedNode was not found.  It's" +
                " expected that method storeNode was used just before this.");
        }
        if (nextNodeBTLIndex.intValue() >= rmNodeBTLIndex.intValue()) {
            throw new IllegalStateException("removedNode is expected to have " +
                " been a more recent backtrack than nextNode.");
        }
        
        List<PathStep> btList = backTrackNodes.get(key0);
        
        Set<PairInt> interactionTestNodes = new HashSet<PairInt>();        
        
        for (int idx = (nextNodeBTLIndex.intValue() + 1); 
            idx <= rmNodeBTLIndex.intValue(); ++idx) {
            
            PathStep node = btList.get(idx);
            interactionTestNodes.add(node.coords);
        }
        
        List<PathStep> output = new ArrayList<PathStep>();
        
        for (int idx = (nextNodeBTLIndex.intValue() - 1); idx > -1; --idx) {
            
            PathStep node = btList.get(idx);
            
            int x = node.coords.getX();
            int y = node.coords.getY();
            
            boolean adj = false;
            for (PairInt t : interactionTestNodes) {
                int diffX = Math.abs(t.getX() - x);
                int diffY = Math.abs(t.getY() - y);
                if ((diffX < 2) && (diffY < 2)) {
                    adj = true;
                    break;
                }
            }
            if (adj) {
                break;
            }
            output.add(node);
        }
        
        if (output.isEmpty()) {
            return null;
        }
        
        return output;
    }
    
}
