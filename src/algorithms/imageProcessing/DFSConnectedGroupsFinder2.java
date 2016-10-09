package algorithms.imageProcessing;

import algorithms.imageProcessing.util.MiscStats;
import algorithms.misc.Misc;
import algorithms.util.PairInt;
import java.util.Arrays;
import java.util.HashMap;
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
    
    /**
     * uses the 4 neighbor region if true, else the 8-neighbor region
     */
    protected boolean use4Neighbors = false;
    
    private State state = null;
    
    public DFSConnectedGroupsFinder2() {
    
        super();
        
        minimumNumberInCluster = 1;
        
        state = State.INITIALIZED;
    }
    
    public void setToUse4Neighbors() {
        use4Neighbors = true;
    }
    
    @Override
    Logger constructLogger() {
        return Logger.getLogger(DFSConnectedGroupsFinder2.class.getName());
    }
    
    /**
     * find contiguous pixels near one another in value and including logic to
     * wrap around a set of values (0 to 360, for example, should see 0 and 359
     * as being within a tolerance of '1' from one another). Note that the
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
                
        if (correctForWandering) {
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
        
        int[] dxs;
        int[] dys;
        if (use4Neighbors) {
            dxs = Misc.dx4;
            dys = Misc.dy4;
        } else {
            dxs = Misc.dx8;
            dys = Misc.dy8;
        }        

        while (!stack.isEmpty()) {

            PairInt uPoint = stack.pop();
            
            int uX = uPoint.getX();
            int uY = uPoint.getY();

            boolean foundANeighbor = false;
            
            float uValue = pointValueMap.get(uPoint).floatValue();
            
            //(1 + frac)*O(N) where frac is the fraction added back to stack
            
            for (int i = 0; i < dxs.length; ++i) {
                
                int vX = uX + dxs[i];
                int vY = uY + dys[i];
                
                if ((vX < 0) || (vX > (imageWidth - 1))) {
                    continue;
                }
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
        
        // traverse the range from largest index to smallest so can append
        // to lists without affecting smaller indexes
        for (int idx = (getNumberOfGroups() - 1); idx > -1; --idx) {
            
            if (groupMembership.get(idx).size() < 2) {
                continue;
            }
             
            Map<PairInt, Float> thetaMap = createSubMap(groupMembership.get(idx),
                pointValueMap);
            
            int[] startEndValues = MiscStats.determineStartEndValues(thetaMap, 
                maxValueForWrapAround, toleranceInValue);
            
            int range = (startEndValues[0] <= startEndValues[1]) ?
                (startEndValues[1] - startEndValues[0]) :
                (startEndValues[1] + (maxValueForWrapAround - startEndValues[0]));
            
            if (range <= toleranceInValue) {
                continue;
            }
            
            int nBefore = groupMembership.get(idx).size();
            
            List<Set<PairInt>> subsetsWithinTolerance = 
                findSubsetsOfRangesWithinTolerance(thetaMap, 
                maxValueForWrapAround, toleranceInValue, imageWidth, 
                imageHeight);
            
            assert(subsetsWithinTolerance.size() > 0);            
            
            // replace current group:
            groupMembership.set(idx, subsetsWithinTolerance.get(0));
            
            int nAfterTot = groupMembership.get(idx).size();
            
            // append remaining and update pointToGroupMap
            for (int j = 1; j < subsetsWithinTolerance.size(); ++j) {
                
                Set<PairInt> group = subsetsWithinTolerance.get(j);
                
                nAfterTot += group.size();
                
                int idx2 = groupMembership.size();
                
                groupMembership.add(group);
                
                // update Map<PairInt, Integer> pointToGroupMap
                updatePointToGroupMap(group, idx2);
            }
            if (nBefore != nAfterTot) {
                int z = 1;
            }
            assert(nBefore == nAfterTot);
        }
       
        state = State.POST_GROUP_CORRECTED;
    }
    
    private void updatePointToGroupMap(Set<PairInt> group, int idx) {
        
        Integer index = Integer.valueOf(idx);
        
        for (PairInt p : group) {
            this.pointToGroupMap.put(p, index);
        }
    }
    
    private Map<PairInt, Float> createSubMap(Set<PairInt> points, 
        Map<PairInt, Float> pointValueMap) {
        
        Map<PairInt, Float> subMap = new HashMap<PairInt, Float>();
        
        for (PairInt p : points) {
            subMap.put(p, pointValueMap.get(p));
        }
        
        return subMap;
    }

    private List<Set<PairInt>> findSubsetsOfRangesWithinTolerance( 
        Map<PairInt, Float> thetaMap,
        int maxValueForWrapAround, int toleranceInValue, 
        int imageWidth, int imageHeight) {
        
        DFSConnectedGroupsFinder3 finder = new DFSConnectedGroupsFinder3();
        if (use4Neighbors) {
            finder.setToUse4Neighbors();
        }
        List<Set<PairInt>> sets = finder.findConnectedPointGroups(thetaMap, 
            maxValueForWrapAround, toleranceInValue, imageWidth, imageHeight);
        
        return sets;
    }
  
}
