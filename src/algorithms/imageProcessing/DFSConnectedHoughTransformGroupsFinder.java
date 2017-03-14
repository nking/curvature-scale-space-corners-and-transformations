package algorithms.imageProcessing;

import algorithms.misc.Misc;
import algorithms.util.PairInt;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;

/**
 * @author nichole
 */
public class DFSConnectedHoughTransformGroupsFinder extends AbstractDFSConnectedGroupsFinder {
        
    protected boolean notValue = false;
        
    /**
     * uses the 4 neighbor region if true, else the 8-neighbor region
     */
    protected boolean use4Neighbors = false;
    
    public DFSConnectedHoughTransformGroupsFinder() {
                
    }
    
    Logger constructLogger() {
        return Logger.getLogger(this.getClass().getName());
    }
    
    public void setMinimumNumberInCluster(int n) {
        this.minimumNumberInCluster = n;
    }
    
    public void setToUse4Neighbors() {
        use4Neighbors = true;
    }
    
    /**
     * given a map with points have theta between 0 and 359, inclusive,
     * find the groups of adjacent points having same theta and radius
     * within tolerances.
     * @param pointTRMap map w/ key = pixel coordinates (x,y) and value = 
     *        pairint of theta in degrees and radius from origin.
     * @param thetaTolerance 
     * @param radiusTolerance 
     */
    public void findConnectedPointGroups(Map<PairInt, PairInt> pointTRMap, 
        int thetaTolerance, int radiusTolerance, int wrapAroundValue) {
                
        findClustersIterative(pointTRMap, thetaTolerance, radiusTolerance,
            wrapAroundValue);
        
        prune();        
    }

    protected void findClustersIterative(Map<PairInt, PairInt> pointTRMap, 
        int thetaTolerance, int radiusTolerance, int wrapAroundValue) {
        
        if (pointTRMap.isEmpty()) {
            return;
        }
        
        int[] dxs;
        int[] dys;
        if (use4Neighbors) {
            dxs = Misc.dx4;
            dys = Misc.dy4;
        } else {
            dxs = Misc.dx8;
            dys = Misc.dy8;
        }
                
        Set<PairInt> visited = new HashSet<PairInt>();
        
        java.util.Stack<PairInt> stack = new java.util.Stack<PairInt>();
        
        //O(N)
        stack.addAll(pointTRMap.keySet());
               
        while (!stack.isEmpty()) {

            PairInt uPoint = stack.pop();
            
            if (visited.contains(uPoint)) {
                continue;
            }

            int uX = uPoint.getX();
            int uY = uPoint.getY();
            
            PairInt uTR = pointTRMap.get(uPoint);

            boolean foundANeighbor = false;
            
            //(1 + frac)*O(N) where frac is the fraction added back to stack
            
            for (int i = 0; i < dxs.length; ++i) {
                
                int vX = uX + dxs[i];
                int vY = uY + dys[i];
                
                PairInt vPoint = new PairInt(vX, vY);
                
                PairInt vTR = pointTRMap.get(vPoint);
                
                if (vTR == null) {
                    continue;
                }
                
                // if v has same theta within tolerance and same radius within
                // tolerance, add to same group.
                
                if (Math.abs(uTR.getY() - vTR.getY()) > radiusTolerance) {
                    continue;
                }
                if (Math.abs(uTR.getX() - vTR.getX()) > thetaTolerance) {
                    if (Math.abs(uTR.getX() - (vTR.getX() - wrapAroundValue)) > thetaTolerance) {
                        if (Math.abs(vTR.getX() - (uTR.getX() - wrapAroundValue)) > thetaTolerance) {
                            continue;
                        }
                    }
                }

                processPair(uPoint, vPoint);

                stack.add(vPoint);

                foundANeighbor = true;
            }
            
            if (!foundANeighbor && (minimumNumberInCluster == 1)) {                
                process(uPoint);
            }
            
            visited.add(uPoint);
        }
    }
  
}
