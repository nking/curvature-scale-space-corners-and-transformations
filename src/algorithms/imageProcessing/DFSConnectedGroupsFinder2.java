package algorithms.imageProcessing;

import algorithms.util.PairInt;
import java.util.Map;
import java.util.Map.Entry;

/**
 * a class to find contiguous pixels near one another in value and including
 * logic to wrap around a set of values (0 to 360, for example, should see
 * 0 and 359 is being within a tolerance of '1' from one another). Note that
 * the minimum possible value is always 0.
 * 
 * @author nichole
 */
public class DFSConnectedGroupsFinder2 extends DFSConnectedGroupsFinder {

    public DFSConnectedGroupsFinder2() {
        minimumNumberInCluster = 1;
    }
        
    public void findConnectedPointGroups(Map<PairInt, Float> pointValueMap, 
        int maxValueForWrapAround, int toleranceInValue,
        int imageWidth, int imageHeight) {
          
        findClustersIterative(pointValueMap, maxValueForWrapAround, 
            toleranceInValue, imageWidth, imageHeight);
        
        prune();
    }

    protected void findClustersIterative(Map<PairInt, Float> pointValueMap, 
        int maxValueForWrapAround, int toleranceInValue,
        int imageWidth, int imageHeight) {
        
        if (pointValueMap.isEmpty()) {
            return;
        }
        
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
    }
    
}
