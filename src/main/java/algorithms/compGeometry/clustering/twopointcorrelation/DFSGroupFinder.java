package algorithms.compGeometry.clustering.twopointcorrelation;

import java.util.logging.Logger;

/**
 * use a Depth First Search algorithm to visit all nodes in a dataset and
 * put into groups those within  distance < threshhold * threshholdFactor from each other.
 * 
 * the number of times visit is called is < O(N^2) and the execution steps within that
 * method are only partially completed as soon as separation is found to be too large.
 * 
 * one test of 55 datapoints shows runtime complexity of O(N^1.8)
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

        for (int uSortedXIndex = 0; uSortedXIndex < sortedXIndexes.length; uSortedXIndex++) {
            
            if (color[uSortedXIndex] == 0) {
                visit(uSortedXIndex);
            }
        }  
    }
    
    protected void visit(int uSortedXIndex) {
        
        color[uSortedXIndex] = 1;
        
        float uX = indexer.getX()[ sortedXIndexes[uSortedXIndex] ];
        float minXAssoc = uX - thrsh;
        float maxXAssoc = uX + thrsh;
        float uY = indexer.getY()[ sortedXIndexes[uSortedXIndex] ];
        float minYAssoc = uY - thrsh;
        float maxYAssoc = uY + thrsh;
        
        // iterate over uNode neighbors v.  can skip calc if color[v] != 0
        for (int vSortedXIndex = 1; vSortedXIndex < sortedXIndexes.length; vSortedXIndex++) {
            

            if (color[vSortedXIndex] != 0) {
                continue;
            }
            
            float vX = indexer.getX()[ sortedXIndexes[vSortedXIndex] ];
            
            if (vX < minXAssoc) {
                //TODO:  I think I have to add to queue for later processing
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
               
            // if we're here, vX is within assoc distance of u and so is vY so this is a neighbor
                        
            processPair(indexer, uSortedXIndex, vSortedXIndex);
            
            visit(vSortedXIndex);
        }
       
        color[uSortedXIndex] = 2;
    }
    
    protected void processPair(DoubleAxisIndexer indexer, int uSortedXIndex, int vSortedXIndex) {
        
        //log.finest("processPair " + uSortedXIndex + ":" + vSortedXIndex);           

        // u and v should be in a group together
        // v has not been assigned to a group yet, but u may have
        // if u has a group id already, add v
        // else create a group and add u and v
        
        int uIdx = indexer.getSortedXIndexes()[uSortedXIndex];
        
        int vIdx = indexer.getSortedXIndexes()[vSortedXIndex];
        
        int groupId;
        
        /*float tmpux = indexer.getX()[ sortedXIndexes[uSortedXIndex] ];
        float tmpuy = indexer.getY()[ sortedXIndexes[uSortedXIndex] ];
        float tmpvx = indexer.getX()[ sortedXIndexes[vSortedXIndex] ];
        float tmpvy = indexer.getY()[ sortedXIndexes[vSortedXIndex] ];*/
        
        if (pointToGroupIndex[uIdx] > -1) {
        
            groupId = pointToGroupIndex[uIdx];
            
            pointToGroupIndex[vIdx] = groupId;
            
            groupMembership[groupId].insert(vIdx);
            
        } else if (pointToGroupIndex[vIdx] > -1) {
            
            groupId = pointToGroupIndex[vIdx];
            
            pointToGroupIndex[uIdx] = groupId;
            
            groupMembership[groupId].insert(uIdx);
            
        } else {
            
            checkAndExpandGroupMembershipArray();
            
            groupId = nGroups;
            
            pointToGroupIndex[uIdx] = groupId;
            
            pointToGroupIndex[vIdx] = groupId; 
            
            groupMembership[groupId].insert(uIdx);
            
            groupMembership[groupId].insert(vIdx);
            
            nGroups++;
        }
        
        //System.out.println(groupId + ") processPair  (" + tmpux + "," + tmpuy + ")");
        //System.out.println(groupId + ") processPair  (" + tmpvx + "," + tmpvy + ")");
    }

}
