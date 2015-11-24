package algorithms.imageProcessing;

import algorithms.util.PairInt;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
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
        
    public DFSConnectedGroupsFinder() {
                
        this.log = Logger.getLogger(this.getClass().getName());
    }
        
    public void setMinimumNumberInCluster(int n) {
        this.minimumNumberInCluster = n;
    }
  
    public void findConnectedPointGroups(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
                
        findClustersIterative(points, imageWidth, imageHeight);
        
        prune();        
    }

    protected void findClustersIterative(Set<PairInt> points, int imageWidth, 
        int imageHeight) {
        
        if (points.isEmpty()) {
            return;
        }
        
        java.util.Stack<PairInt> stack = new java.util.Stack<PairInt>();
        
        //O(N)
        stack.addAll(points);
               
        visited.add(stack.peek());

        while (!stack.isEmpty()) {

            PairInt uPoint = stack.pop();
            
            int uX = uPoint.getX();
            int uY = uPoint.getY();

            boolean foundANeighbor = false;
            
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
                    
                    if (!points.contains(vPoint)) {
                        continue;
                    }
                    
                    if (visited.contains(vPoint)) {
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
