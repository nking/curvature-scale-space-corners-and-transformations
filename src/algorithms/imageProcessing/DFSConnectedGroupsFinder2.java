package algorithms.imageProcessing;

import algorithms.util.PairInt;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.logging.Logger;

/**
 * a class to find contiguous pixels near one another in value and including
 * logic to wrap around a set of values (0 to 360, for example, should see
 * 0 and 359 is being within a tolerance of '1' from one another). Note that
 * the minimum possible value is always 0.
 * 
 * @author nichole
 */
public class DFSConnectedGroupsFinder2 extends AbstractDFSConnectedGroupsFinder {

    private enum State {
        INITIALIZED, GROUPS_FOUND, GROUPS_PRUNED, POST_GROUP_CORRECTED
    }
    
    private State state = null;
    
    public DFSConnectedGroupsFinder2() {
        
        minimumNumberInCluster = 1;
        
        state = State.INITIALIZED;
    }
    
    Logger constructLogger() {
        return Logger.getLogger(DFSConnectedGroupsFinder2.class.getName());
    }
    
    /**
     * find contiguous pixels near one another in value and including logic to
     * wrap around a set of values (0 to 360, for example, should see 0 and 359
     * ss being within a tolerance of '1' from one another). Note that the
     * minimum possible value is always 0.
     * To correct for a group that has wandered from a total range of tolerance,
     * setTheCorrectForWandering to true
     *
     * @param pointValueMap
     * @param maxValueForWrapAround
     * @param toleranceInValue
     * @param imageWidth
     * @param imageHeight
     * @param correctForWandering if true, uses another algorithm after the
     * groups are found to make subsets within each if the range is larger 
     * than tolerance.
     */
    public void findConnectedPointGroups(Map<PairInt, Float> pointValueMap, 
        int maxValueForWrapAround, int toleranceInValue,
        int imageWidth, int imageHeight, boolean correctForWandering) {
          
        findClustersIterative(pointValueMap, maxValueForWrapAround, 
            toleranceInValue, imageWidth, imageHeight);
        
        prune();
                
        if (false && correctForWandering) {
            correctRangesIfNeeded(pointValueMap, maxValueForWrapAround,
                toleranceInValue, imageWidth, imageHeight);
        }
    }

    @Override
    protected void prune() {
        
        if (!state.equals(State.GROUPS_FOUND)) {
            throw new IllegalStateException(
                "findClustersIterative must be used before this");
        }
        
        super.prune();
        
        state = State.GROUPS_PRUNED;
    }
  
    protected void findClustersIterative(Map<PairInt, Float> pointValueMap, 
        int maxValueForWrapAround, int toleranceInValue,
        int imageWidth, int imageHeight) {
        
        if (pointValueMap.isEmpty()) {
            state = State.GROUPS_FOUND;
            return;
        }
        
        Set<PairInt> visited = new HashSet<PairInt>();
        
        java.util.Stack<PairInt> stack = new java.util.Stack<PairInt>();
        
        //O(N)
        for (Entry<PairInt, Float> entry : pointValueMap.entrySet()) {
            stack.add(entry.getKey());
        }
               
        visited.add(stack.peek());

        while (!stack.isEmpty()) {

            PairInt uPoint = stack.pop();
            
            int uX = uPoint.getX();
            int uY = uPoint.getY();

            boolean foundANeighbor = false;
            
            float uValue = pointValueMap.get(uPoint).floatValue();
            
            //(1 + frac)*O(N) where frac is the fraction added back to stack
            
            for (int vX = (uX - 1); vX <= (uX + 1); vX++) {
                if ((vX < 0) || (vX > (imageWidth - 1))) {
                    continue;
                }
                
                for (int vY = (uY - 1); vY <= (uY + 1); vY++) {
                    if ((vY < 0) || (vY > (imageHeight - 1))) {
                        continue;
                    }
                    
                    PairInt vPoint = new PairInt(vX, vY);
                    
                    if (vPoint.equals(uPoint)) {
                        continue;
                    }
                    
                    if (!pointValueMap.containsKey(vPoint)) {
                        continue;
                    }
                    
                    if (visited.contains(vPoint)) {
                        continue;
                    }
                    
                    boolean similar = false;
                    
                    float vValue = pointValueMap.get(vPoint).floatValue();
                    if (Math.abs(uValue - vValue) <= toleranceInValue) {
                        similar = true;
                    }
                    
                    // test wrap around
                    if (!similar) {
                        if (Math.abs(uValue - (vValue - maxValueForWrapAround)) <= toleranceInValue) {
                            similar = true;
                        }
                        /*
                            0        360
                            |  u    v |
                          v |  u      |
                        */
                    }
                    if (!similar) {
                        if (Math.abs(vValue - (uValue - maxValueForWrapAround)) <= toleranceInValue) {
                            similar = true;
                        }
                        /*
                            0        360
                            |  v    u |
                          u |  v      |
                        */
                    }
                    
                    if (!similar) {
                        continue;
                    }
                    
                    visited.add(vPoint);
                    
                    processPair(uPoint, vPoint);
                
                    stack.add(vPoint);
                    
                    foundANeighbor = true;
                }
            }
            
            if (!foundANeighbor && (minimumNumberInCluster == 1)) {
                
                process(uPoint);
            }
        }
        
        state = State.GROUPS_FOUND;
    }

    private void correctRangesIfNeeded(Map<PairInt, Float> pointValueMap, 
        int maxValueForWrapAround, int toleranceInValue, int imageWidth, 
        int imageHeight) {
          
        if (!state.equals(State.GROUPS_PRUNED)) {
            throw new IllegalStateException(
                "prune() must be used before this");
        }
        
        throw new UnsupportedOperationException("not yet implemented");
            
        /*
        // these are ordered by index
        List<Integer> groupsLargerThanRange = findGroupsLargerThanRange(
            maxValueForWrapAround, toleranceInValue);
        
        if (groupsLargerThanRange.size() == 0) {
            return;
        }
        
        // traverse the range from largest index to smallest so can append
        // to lists without affecting smaller indexes
        for (int i = (groupsLargerThanRange.size() - 1); i > -1; --i) {
            
            int idx = groupsLargerThanRange.get(i);
            
            Map<PairInt, Float> thetaMap = createSubset(groupMembership.get(idx),
                pointValueMap);
            
            List<Set<PairInt>> subsetsWithinTolerance = 
                findSubsetsOfRangesWithinTolerance(thetaMap, maxValueForWrapAround,
                    toleranceInValue, imageWidth, imageHeight);
            
            assert(subsetsWithinTolerance.size() > 1);
            
            int nBefore = groupMembership.get(idx).size();
            
            // replace current group:
            groupMembership.set(idx, subsetsWithinTolerance.get(0));
            
            int nAfterTot = groupMembership.get(idx).size();
            
            // append remaining and update pointToGroupMap
            for (int j = 1; j < subsetsWithinTolerance.size(); ++j) {
                
                Set<PairInt> group = subsetsWithinTolerance.get(j);
                
                int idx2 = groupMembership.size();
                
                subsetsWithinTolerance.add(group);
                
                // update Map<PairInt, Integer> pointToGroupMap
                updatePointToGroupMap(group, idx2);
            }
        }
       
        state = State.POST_GROUP_CORRECTED;
        */
    }
    
    private List<Set<PairInt>> findSubsetsOfRangesWithinTolerance( 
        Map<PairInt, Float> thetaMap,
        int maxValueForWrapAround, int toleranceInValue, 
        int imageWidth, int imageHeight) {
        
        DFSConnectedGroupsFinder3 finder = new DFSConnectedGroupsFinder3();
        
        finder.findConnectedPointGroups(thetaMap, maxValueForWrapAround, 
            toleranceInValue, imageWidth, imageHeight);

        int nGroups = finder.getNumberOfGroups();
        
        List<Set<PairInt>> sets = new ArrayList<Set<PairInt>>(nGroups);
        
        for (int i = 0; i < nGroups; ++i) {
            sets.add(finder.getXY(i));
        }
        
        return sets;
    }

}
