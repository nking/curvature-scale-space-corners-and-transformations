package algorithms.imageProcessing;

import algorithms.misc.Misc;
import algorithms.util.PairInt;
import java.util.HashSet;
import java.util.Set;
import java.util.logging.Logger;

/**
 * adapted from
 * https://code.google.com/p/two-point-correlation/source/browse/src/main/java/algorithms/compGeometry/clustering/twopointcorrelation/AbstractGroupFinder.java
 * and
 * https://code.google.com/p/two-point-correlation/source/browse/src/main/java/algorithms/compGeometry/clustering/twopointcorrelation/DFSGroupFinder.java
 * under MIT License (MIT), Nichole King 2013
 * 
 * @author nichole
 */
public class DFSConnectedGroupsFinder extends AbstractDFSConnectedGroupsFinder {
        
    protected boolean notValue = false;
        
    /**
     * uses the 4 neighbor region if true, else the 8-neighbor region
     */
    protected boolean use4Neighbors = true;
    
    public DFSConnectedGroupsFinder() {
                
    }
    
    Logger constructLogger() {
        return Logger.getLogger(DFSConnectedGroupsFinder2.class.getName());
    }
    
    public void setMinimumNumberInCluster(int n) {
        this.minimumNumberInCluster = n;
    }
    
    public void setToUse8Neighbors() {
        use4Neighbors = false;
    }
    
    public void findConnectedPointGroups(Set<PairInt> points) {
        
        //TODO: make an adjacency map at this point and replace
        // the loop over neighbor offsets with loop over adjacent points
        // to speed up the algorithm.        
        
        findClustersIterative(points);
        
        prune();        
    }

    protected void findClustersIterative(Set<PairInt> points) {
        
        if (points.isEmpty()) {
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
        stack.addAll(points);
               
        while (!stack.isEmpty()) {

            PairInt uPoint = stack.pop();
            
            if (visited.contains(uPoint)) {
                continue;
            }

            int uX = uPoint.getX();
            int uY = uPoint.getY();

            boolean foundANeighbor = false;
            
            //(1 + frac)*O(N) where frac is the fraction added back to stack
            
            for (int i = 0; i < dxs.length; ++i) {
                
                int vX = uX + dxs[i];
                int vY = uY + dys[i];
            
                PairInt vPoint = new PairInt(vX, vY);

                if (vPoint.equals(uPoint)) {
                    continue;
                }

                if (!points.contains(vPoint)) {
                    continue;
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
