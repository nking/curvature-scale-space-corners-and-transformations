package algorithms.compGeometry.clustering.twopointcorrelation;

import java.util.logging.Logger;

import algorithms.compGeometry.LinesAndAngles;
import algorithms.util.Stack;

/**
 * <pre>
 * An implementation of IGroupFinder that uses a depth first search algorithm
 * to visit all nodes in a dataset and put into groups those within 
 * distance < threshhold * threshholdFactor from each other.
 * 
 * the runtime complexity is a little larger, but on the order of O(N^2)
 * 
 * Note that the recursive solution has fewer visits, but the recursion requires a method
 * frame to be loaded for each invocation. the method frame load and unload is expensive computationally,
 * and the other caveat is that if the recursion is deep enough it may cause the program to
 * become memory bound and that further decreases performance or can halt the program.
 * </pre>
 * 
 * @author nichole
 */
public class DFSGroupFinder extends AbstractGroupFinder {

    // 0 = unvisited, 1 = processing, 2 = visited
    protected int[] color = null;
    
    protected int[] sortedXIndexes = null;
    
    protected final float thrsh;
    
    protected AxisIndexer indexer = null;
    
    public DFSGroupFinder(float threshhold, float threshholdFactor) {
        
        super(threshhold, threshholdFactor);
        
        thrsh = threshhold * threshholdFactor;        
    }
    
    public void constructLogger() {
        this.log = Logger.getLogger(this.getClass().getName());
    }   
    
    @Override
    protected void findClusters(AxisIndexer indexer) {
             
        log.info("using a critical density of " + threshhold + " * " 
            + threshholdFactor + " = " + thrsh);
        
        // traverse the data by ordered x values using 
        //     indexer.getSortedXIndexes() so can break when exceed critical distance
        
        this.indexer = indexer;
        
        sortedXIndexes = indexer.getSortedXIndexes();
                
        // 0 = unvisited, 1 = processing, 2 = visited
        color = new int[sortedXIndexes.length];

        //findClustersRecursive();
        
        findClustersIterative();
        
    }
    
    protected void findClustersIterative() {

        Stack stack = new Stack();
        
        for (int uSortedXIndex = (sortedXIndexes.length - 1); uSortedXIndex > -1; uSortedXIndex--) {
            stack.insert(uSortedXIndex);
        }
                
        color[stack.peek()] = 2;
        
        while (!stack.isEmpty()) {
            
            SimpleLinkedListNode uNode = stack.pop();
            
            int uSortedXIndex = uNode.getKey();
            
            // we process the pair when their point density is higher than thrsh:
            //  
            //   2/sep_u_v  > thrsh  where thrsh is 2.5*the background linear density
            //
            //   if want to stop the search along x axis when have surpassed an association distance,
            //      we can see  2/thrsh > sep_u_v
            //                  (2/thrsh) > u - v 
            //                  (2/thrsh) + v > u
            
            // association of 2 points for separation <= critSeparation
            
            float critSeparation = 2.f/thrsh;

            float uX = indexer.getX()[ sortedXIndexes[uSortedXIndex] ];
            float minXAssoc = uX - critSeparation;
            float maxXAssoc = uX + critSeparation;
            float uY = indexer.getY()[ sortedXIndexes[uSortedXIndex] ];
            float minYAssoc = uY - critSeparation;
            float maxYAssoc = uY + critSeparation;
            
            // for each neighbor v of u
            for (int vSortedXIndex = 0; vSortedXIndex < sortedXIndexes.length; 
                vSortedXIndex++) {

                if (color[vSortedXIndex] != 0 || (uSortedXIndex == vSortedXIndex)) {
                    continue;
                }
                
                float vX = indexer.getX()[ sortedXIndexes[vSortedXIndex] ];
                
                if (vX < minXAssoc) {
                    continue;
                }
                
                if (vX > maxXAssoc) {
                    // we've past the distance in the ordered x array of points associated with u,                 
                    break;
                }
                
                float vY = indexer.getY()[ sortedXIndexes[vSortedXIndex] ];
                
                if ((vY < minYAssoc) || (vY > maxYAssoc)) {
                    continue;
                }
                
                // one last check using the true separation
                
                double sep = Math.sqrt(LinesAndAngles.distSquared(uX, uY, vX, vY));
                  
                if (sep > critSeparation) {
                    continue;
                }
                
                log.fine("  comparing: (" + uX + "," + uY + ")  (" 
                    + vX + "," + vY + ")  sep=" + sep + " cr=" + critSeparation);
      
                color[vSortedXIndex] = 2;
                
                processPair(indexer, uSortedXIndex, vSortedXIndex);
                
                // inserting back at the top of the stack assures that the search continues next from an associated point
                stack.insert(vSortedXIndex);
            }
        }
    }
    
    protected void processPair(AxisIndexer indexer, int uSortedXIndex, int vSortedXIndex) {
        
        //log.finest("processPair " + uSortedXIndex + ":" + vSortedXIndex);           
        
        int uIdx = indexer.getSortedXIndexes()[uSortedXIndex];
        
        int vIdx = indexer.getSortedXIndexes()[vSortedXIndex];
        
        int groupId;
        
        /*
        float tmpux = indexer.getX()[ sortedXIndexes[uSortedXIndex] ];
        float tmpuy = indexer.getY()[ sortedXIndexes[uSortedXIndex] ];
        float tmpvx = indexer.getX()[ sortedXIndexes[vSortedXIndex] ];
        float tmpvy = indexer.getY()[ sortedXIndexes[vSortedXIndex] ];
        */
        
        if ((pointToGroupIndex[uIdx] > -1) && (pointToGroupIndex[vIdx] == -1)) {
        
            groupId = pointToGroupIndex[uIdx];
            
            pointToGroupIndex[vIdx] = groupId;
            
            groupMembership[groupId].insert(vIdx);
            
        } else if ((pointToGroupIndex[vIdx] > -1) && (pointToGroupIndex[uIdx] == -1)) {
            
            groupId = pointToGroupIndex[vIdx];
            
            pointToGroupIndex[uIdx] = groupId;
            
            groupMembership[groupId].insert(uIdx);
            
        } else if ((pointToGroupIndex[uIdx] == -1) && (pointToGroupIndex[vIdx] == -1)) {
            
            checkAndExpandGroupMembershipArray();
            
            groupId = nGroups;
            
            pointToGroupIndex[uIdx] = groupId;
            
            pointToGroupIndex[vIdx] = groupId; 
            
            groupMembership[groupId].insert(uIdx);
            
            groupMembership[groupId].insert(vIdx);
     
            nGroups++;
            
        } else {
            
            groupId = -1;
            
            log.finest("not reading " + uIdx + ":" + vIdx );
        }
        
        /*
        log.fine(groupId + ") processPair  (" + tmpux + "," + tmpuy + ")");
        log.fine(groupId + ") processPair  (" + tmpvx + "," + tmpvy + ")");
        */
    }
 
}
