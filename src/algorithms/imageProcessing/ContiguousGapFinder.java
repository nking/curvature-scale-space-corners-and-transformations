package algorithms.imageProcessing;

import algorithms.misc.Misc;
import algorithms.misc.MiscMath;
import algorithms.util.PairInt;
import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.logging.Logger;

/**
 * 
 * @author nichole
 */
public class ContiguousGapFinder {
    
   /**
     * an array to hold each group as an item.  each item contains a key which is an index
     * to arrays indexer.x, indexer.y and this.pointToGroupIndex
     */
    protected List<Set<PairInt> > groupMembership = new ArrayList<Set<PairInt> >();
        
    protected Logger log = null;
        
    /*
     * map w/ key =pairint coord,
       value = group it belongs to.
     */
    protected TObjectIntMap<PairInt> pointToGroupMap = new
        TObjectIntHashMap<PairInt>();
    
    protected int minimumNumberInCluster = 3;
    
    protected boolean debug = false;
    
    protected final Set<PairInt> points;
    
    private final int minX;
    private final int maxX;
    private final int minY;
    private final int maxY;
    
    /**
     * uses the 4 neighbor region if true, else the 8-neighbor region
     */
    protected boolean use4Neighbors = true;
    
    public ContiguousGapFinder(final Set<PairInt> thePoints) {
        
        this.points = thePoints;
        
        this.log = Logger.getLogger(this.getClass().getName());
        
        int[] minMaxXY = MiscMath.findMinMaxXY(points);
        
        this.minX = minMaxXY[0];
        this.maxX = minMaxXY[1];
        this.minY = minMaxXY[2];
        this.maxY = minMaxXY[3];
    }
   
    public void setMinimumNumberInCluster(int n) {
        this.minimumNumberInCluster = n;
    }
    
    public void setDebug(boolean setDebugToTrue) {
        this.debug = setDebugToTrue;
    }
    
    public void setToUse8Neighbors() {
        use4Neighbors = false;
    }
    
    public void findGaps() {
            
        findGapClusters();
        
        prune();        
    }
    
    protected void findGapClusters() {
        
        int[] dxs;
        int[] dys;
        if (use4Neighbors) {
            dxs = Misc.dx4;
            dys = Misc.dy4;
        } else {
            dxs = Misc.dx8;
            dys = Misc.dy8;
        }
            
        for (int i = minX; i <= maxX; ++i) {
            for (int j = minY; j <= maxY; ++j) {
                PairInt p = new PairInt(i, j);
                if (points.contains(p)) {
                    continue;
                }
                for (int k = 0; k < dxs.length; ++k) {
                    int x2 = i + dxs[k];
                    int y2 = j + dys[k];
                    if (x2 < minX || y2 < minY || x2 > maxX ||
                        y2 > maxY) {
                        continue;
                    }
                    PairInt p2 = new PairInt(x2, y2);
                    if (!points.contains(p2)) {
                        processPair(p, p2);
                    }
                }
            }
        }
    }
    
    protected void processPair(PairInt u, PairInt v) {
                
        int groupId = -1;
        if (pointToGroupMap.containsKey(u)) {
            groupId = pointToGroupMap.get(u);
        }
        
        boolean containsV = pointToGroupMap.containsKey(v);
        
        if ((groupId > -1) && !containsV) {
                    
            groupMembership.get(groupId).add(v);
            
            pointToGroupMap.put(v, groupId);
                        
        } else if ((groupId == -1) && containsV) {

            // usually, u will have been added before v is visited, so this
            // block is rarely used
            
            groupId = pointToGroupMap.get(v);

            groupMembership.get(groupId).add(u);
            
            pointToGroupMap.put(u, groupId);
            
        } else if ((groupId == -1) && !containsV) {
                        
            groupId = groupMembership.size();
            
            pointToGroupMap.put(u, groupId);
            
            pointToGroupMap.put(v, groupId);
            
            Set<PairInt> set = new HashSet<PairInt>();
            set.add(u);
            set.add(v);
            
            groupMembership.add(set);
                                  
        }
        
    }
       
    public List<Set<PairInt> > getGapGroupMembershipList() {
        return groupMembership;
    }

    public int getNumberOfGapGroups() {
        return groupMembership.size();
    }

    public TObjectIntMap<PairInt> getPointToGroupIndexes() {
        return pointToGroupMap;
    }
    
    /**
     * remove groups smaller than minimumNumberInCluster
     */
    protected void prune() {
                
        log.finest("number of groups before prune=" + groupMembership.size());
        
        /*
         * [------] 0
         * [------] 1 <---- too few
         * [------] 2
         */

        log.fine("PRUNE: nGroups before prune =" + groupMembership.size());

        // iterate backwards so can move items up without conflict with iterator
        for (int i = (groupMembership.size() - 1); i > -1; i--) {
            
            Set<PairInt> group = groupMembership.get(i);
            
            int count = group.size();
            
            if (count < minimumNumberInCluster) {
                             
                // remove this group and move up all groups w/ index > i by one index
                for (int j = (i + 1); j < groupMembership.size(); j++) {
                    
                    int newGroupId = j - 1;
                                        
                    // update members in pointToGroupIndex
                    Set<PairInt> latest = groupMembership.get(j);
                    //
                    for (PairInt p : latest) {
                        pointToGroupMap.put(p, newGroupId);
                    }
                } 
               
                Set<PairInt> removed = groupMembership.remove(i);
            }
        }
   
        log.fine("number of groups after prune=" + groupMembership.size());
    }

    public int getNumberofGroupMembers(int groupId) {
        
        if (groupMembership.isEmpty()) {
            return 0;
        }
        if (groupId > (groupMembership.size() - 1) || (groupId < 0)) {
            throw new IllegalArgumentException("groupId=" + groupId 
            + " is outside of range of nGroups=" + groupMembership.size());
        }
                
        return groupMembership.get(groupId).size();
    }
    
    public void getXY(final int groupId, final Set<PairInt> output) {
        
        if (groupMembership.isEmpty()) {
            return;
        }
        if (groupId > (groupMembership.size() - 1) || (groupId < 0)) {
            throw new IllegalArgumentException("groupId=" + groupId 
            + " is outside of range of nGroups=" + groupMembership.size());
        }
        
        Set<PairInt> indexes = groupMembership.get(groupId);

        output.addAll(indexes);
    }
    
    public Set<PairInt> getXY(final int groupId) {
        
        if (groupMembership.isEmpty()) {
            return null;
        }
        if (groupId > (groupMembership.size() - 1) || (groupId < 0)) {
            throw new IllegalArgumentException("groupId=" + groupId 
            + " is outside of range of nGroups=" + groupMembership.size());
        }
        
        return groupMembership.get(groupId);
    }

    /**
     * @return the minX
     */
    public int getMinX() {
        return minX;
    }

    /**
     * @return the maxX
     */
    public int getMaxX() {
        return maxX;
    }

    /**
     * @return the minY
     */
    public int getMinY() {
        return minY;
    }

    /**
     * @return the maxY
     */
    public int getMaxY() {
        return maxY;
    }
}
