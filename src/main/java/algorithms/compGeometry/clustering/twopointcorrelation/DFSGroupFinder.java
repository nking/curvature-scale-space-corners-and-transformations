package algorithms.compGeometry.clustering.twopointcorrelation;

import java.util.logging.Logger;

import algorithms.compGeometry.LinesAndAngles;
import algorithms.util.Stack;

/**
 * use a Depth First Search algorithm to visit all nodes in a dataset and
 * put into groups those within  distance < threshhold * threshholdFactor from each other.
 * 
 * the number of times visit is called is < O(N^2) and the execution steps within that
 * method are only partially completed as soon as separation is found to be too large.
 * 
 * one test of 55 datapoints shows runtime complexity of O(N^1.8) for the recursive DFS
 * and O(N^1.88) for iterative DFS.  The advantage to the iterative DFS is in the reduced
 * amount of memory used because the recursive loads a method frame for each invocation.
 * 
 * @author nichole
 */
public class DFSGroupFinder extends AbstractGroupFinder {

    // 0 = unvisited, 1 = processing, 2 = visited
    protected int[] color = null;
    
    protected int[] sortedXIndexes = null;
    
    protected final float thrsh;
    
    protected DoubleAxisIndexer indexer = null;
    
    public DFSGroupFinder(float threshhold, float threshholdFactor) {
        
        super(threshhold, threshholdFactor);
        
        thrsh = threshhold * threshholdFactor;
    }
    
    public void constructLogger() {
        this.log = Logger.getLogger(this.getClass().getName());
    }   
    
    @Override
    protected void findClusters(DoubleAxisIndexer indexer) {
             
        // traverse the data by ordered x values using indexer.getSortedXIndexes() so can break when exceed critical distance
        
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
            
            //  the thrsh is a density
            //   2./(uX-vX)  > thrsh ==>   2/thrsh > (uX-vX)
            //
            
            float cr = 2.f/thrsh;

            float uX = indexer.getX()[ sortedXIndexes[uSortedXIndex] ];
            float minXAssoc = uX - cr;
            float maxXAssoc = uX + cr;
            float uY = indexer.getY()[ sortedXIndexes[uSortedXIndex] ];
            float minYAssoc = uY - cr;
            float maxYAssoc = uY + cr;
            
            // for each neighbor v of u
            for (int vSortedXIndex = 0; vSortedXIndex < sortedXIndexes.length; vSortedXIndex++) {

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
                
                if (sep > cr) {
                    continue;
                }
                                          
                color[vSortedXIndex] = 2;
                
                processPair(indexer, uSortedXIndex, vSortedXIndex);
                
                // inserting back at the top of the stack assures that the search continues next from an associated point
                stack.insert(vSortedXIndex);
            }
        }
    }
    
    protected void findClustersRecursive() {

        for (int uSortedXIndex = 0; uSortedXIndex < sortedXIndexes.length; uSortedXIndex++) {
            
            if (color[uSortedXIndex] == 0) {
                visit(uSortedXIndex);
            }
        } 
    }
    
    protected void visit(int uSortedXIndex) {
        
        color[uSortedXIndex] = 1;
        
        // we process the pair when their point density is higher than thrsh:  that is  (2 points/distance) > thrsh
        //  
        //  2 points  <  thrsh * ( (ux-vx)^2 + (uy-vy)^2 )^(0.5) 
        //  to see the max extent along y, set diff in x's to zero:  
        //     2 points  <  thrsh * (uy-vy) so differences in y smaller than 2./thrsh are in a group
        //
        float cr = 2.f/thrsh;

        float uX = indexer.getX()[ sortedXIndexes[uSortedXIndex] ];
        float minXAssoc = uX - cr;
        float maxXAssoc = uX + cr;
        float uY = indexer.getY()[ sortedXIndexes[uSortedXIndex] ];
        float minYAssoc = uY - cr;
        float maxYAssoc = uY + cr;
        
        // iterate over uNode neighbors v.  can skip calc if color[v] != 0
        for (int vSortedXIndex = 1; vSortedXIndex < sortedXIndexes.length; vSortedXIndex++) {

            if (color[vSortedXIndex] != 0) {
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
            
            if (sep > cr) {
                continue;
            }
                        
            processPair(indexer, uSortedXIndex, vSortedXIndex);
            
            visit(vSortedXIndex);
        }
       
        color[uSortedXIndex] = 2;
    }
    
    protected void processPair(DoubleAxisIndexer indexer, int uSortedXIndex, int vSortedXIndex) {
        
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
            
        } else if ((pointToGroupIndex[vIdx] > -1) && (pointToGroupIndex[uIdx] == -1)){
            
            groupId = pointToGroupIndex[vIdx];
            
            pointToGroupIndex[uIdx] = groupId;
            
            groupMembership[groupId].insert(uIdx);
            
        } else if ((pointToGroupIndex[uIdx] == -1) && (pointToGroupIndex[vIdx] == -1) ){
            
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
        System.out.println(groupId + ") processPair  (" + tmpux + "," + tmpuy + ")");
        System.out.println(groupId + ") processPair  (" + tmpvx + "," + tmpvy + ")");
        */
    }
 
}
