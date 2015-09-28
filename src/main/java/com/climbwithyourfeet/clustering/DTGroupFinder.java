package com.climbwithyourfeet.clustering;

import com.climbwithyourfeet.clustering.util.PairInt;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.logging.Logger;

/**
 *
 * @author nichole
 */
public class DTGroupFinder {
    /**
     * an array to hold each group as an item. Note that the original point
     * given in points is preserved so any specialization information available
     * to the copy() method is present in the grouped points too.
     * One use case is the point being scaled CIE XY colors and a specialization 
     * of PairInt that has a field holding the pixel index
     */
    protected List<Set<PairInt> > groupMembership = new ArrayList<Set<PairInt> >();
    
    protected Logger log = null;
        
     /*
     * map w/ key holding indexes for a point to the group it belongs to.
     * note that the point index is relative to indexer.x and indexer.y
     */
    protected Map<PairInt, Integer> pointToGroupMap = new
        HashMap<PairInt, Integer >();
    
    /**
      key = group receiving merges, value = set of integers to merge into this
      key.  note that values are all greater than key
     */
    protected TreeMap<Integer, Set<Integer>> mergeGroupsBeforePrune = new
        TreeMap<Integer, Set<Integer>>();
    
    protected int minimumNumberInCluster = 3;
    
    protected boolean notValue = false;
    
    protected boolean debug = false;
    
    protected float threshholdFactor = 2.5f;
    
    protected float critDensity = 0;
            
    public DTGroupFinder() {                
        this.log = Logger.getLogger(this.getClass().getName());
    }
        
    public void setMinimumNumberInCluster(int n) {
        this.minimumNumberInCluster = n;
    }
    
    public void setDebug(boolean setDebugToTrue) {
        this.debug = setDebugToTrue;
    }
    
    public void setThreshholdFactor(float factor) {
        this.threshholdFactor = factor;
    }

    /**
     * given the critical density and having the threshold factor, find the
     * groups of points within a critical distance of one another.
     * runtime complexity is O(N_points * lg2(N_points)).
     * @param criticalDensity
     * @param points
     * @param width
     * @param height 
     */
    void calculateGroups(float criticalDensity, Set<PairInt> points) {
        
        this.critDensity = criticalDensity;
        
        /*
        TODO: revisit the scale
        if critDensity == 1, clustered points are adjcent to one another,
        and the final crit sep needs to be '1',
        so a correction to thrsh is made here for that.
        */
        
        float thrsh = criticalDensity * threshholdFactor;
        
        if (critDensity == 1) {
            // fudge to result in critical separation of 1
            thrsh = 2.f;
        }
        
        findGroups(thrsh, points);
        
        merge();
        
        prune(); 
    }
    
    private void findGroups(float thrsh, Set<PairInt> points) {
        
        if (points.isEmpty()) {
            return;
        }
        
        // association of 2 points for separation <= critSeparation
        float critSep = 2.f/thrsh;
        
        if (critSep < 1) {
            // each point is a group
            if (minimumNumberInCluster < 2) {
                setEachPointAsAGroup(points);
            }
            return;
        }
        
        log.info("critSep=" + critSep);
        
        //TODO: consider data structures with point location as part of their
        // structure... spatial indexing, RTrees...
        /*
        Sort the points by x then y to be able to make small searches around a
        point as traverse the sorted array.
        The runtime complexity is O(N*lg2(N)).
        Note that if the dataset were sparse, could assume only need to sort
        on x and could hence use the O(N) counting sort (if N=100 for example, 
        and the dataset were sparse and there are fewer than 7 points with the 
        same y for each x, O(N) and search all points with a specific x is 
        faster than O(N*lg2(N))).
        */
        
        PairInt[] sorted = new PairInt[points.size()];
        int count = 0;
        for (PairInt p : points) {
            sorted[count] = p;
            count++;
        }
        PISort.mergeSortByXThenY(sorted);
                
        for (int i = 0; i < sorted.length; ++i) {
            
            PairInt uPoint = sorted[i];
            
            // process the pair when their point density is higher than thrsh:
            //  
            //   2/sep_u_v  > thrsh  where thrsh is 2.5*the background linear density
            
            // association of 2 points for separation <= critSeparation
            
            float uX = uPoint.getX();
            float minXAssoc = uX - critSep;
            float maxXAssoc = uX + critSep;
            float uY = uPoint.getY();
            float minYAssoc = uY - critSep;
            float maxYAssoc = uY + critSep;
            
            boolean assocFound = false;
            
            /*
            uX = sorted[0].x
                search backward while x >= minXAssoc
                    while x == minXAssoc, only proceed for y >= minYAssoc
                search forward while x <= maxXAssoc
                    while x == maxXAssoc, only proceed for y <= maxYAssoc
            */
            
            // search backward within radius critSep
            for (int j = (i - 1); j > -1; --j) {
                
                PairInt vPoint = sorted[j];
                
                float vX = vPoint.getX();
                
                if (vX < minXAssoc) {
                    break;
                }
                
                if (vX == minXAssoc) {
                    if (vPoint.getY() < minYAssoc) {
                        break;
                    }
                }
                // for given x, if y < minYAssoc can skip, but not break
                if (vPoint.getY() < minYAssoc) {
                    continue;
                }
                
                // check the diagonal distances:
                double sep = Math.sqrt((vX - uX)*(vX - uX) + 
                    (vPoint.getY() - uY)*(vPoint.getY() - uY));
                  
                if (sep > critSep) {
                    continue;
                }
                
                // if arrive here, vX is within an assoc radius and so is vY
                processPair(uPoint, vPoint);
                                
                assocFound = true;
            }
            
            // search forward within radius critSep
            for (int j = (i + 1); j < sorted.length; ++j) {
                
                PairInt vPoint = sorted[j];
                
                float vX = vPoint.getX();
                
                if (vX > maxXAssoc) {
                    break;
                }
                
                if (vX == maxXAssoc) {
                    if (vPoint.getY() > maxYAssoc) {
                        break;
                    }
                }
                // for given x, if y > maxYAssoc can skip, but not break
                if (vPoint.getY() > maxYAssoc) {
                    //TODO: if knew where next x was in sorted, could increment j to that
                    continue;
                }
                
                // check the diagonal distances:
                double sep = Math.sqrt((vX - uX)*(vX - uX) + 
                    (vPoint.getY() - uY)*(vPoint.getY() - uY));
                  
                if (sep > critSep) {
                    continue;
                }
                
                // if arrive here, vX is within an assoc radius and so is vY
                processPair(uPoint, vPoint);
                                
                assocFound = true;
            }
            
            // if no points were within assoc distance, this becomes its
            // own group
            if ((minimumNumberInCluster < 2) && !assocFound) {
                process(uPoint);
            }            
        }
        
    }
    
    int getNumberOfGroups() {
        return groupMembership.size();
    }
    
    Set<PairInt> getGroup(int groupIdx) {
        
        if (groupMembership.isEmpty()) {
            return new HashSet<PairInt>();
        }
        if (groupIdx > (groupMembership.size() - 1) || (groupIdx < 0)) {
            throw new IllegalArgumentException("groupIdx=" + groupIdx 
            + " is outside of range of nGroups=" + groupMembership.size());
        }
        
        Set<PairInt> set = groupMembership.get(groupIdx);
       
        return set;
    }

    private void setEachPointAsAGroup(Set<PairInt> points) {
        
        for (PairInt p : points) {
            
            int sz = groupMembership.size();
            
            pointToGroupMap.put(p, Integer.valueOf(sz));
            
            Set<PairInt> set = new HashSet<PairInt>();
            
            set.add(p);
            
            groupMembership.add(set);
        }
    }

    private void processPair(PairInt uPoint, PairInt vPoint) {
        
        Integer groupId = pointToGroupMap.get(uPoint);
        
        Integer vGroupId = pointToGroupMap.get(vPoint);
        
        if ((groupId != null) && (vGroupId == null)) {
                    
            groupMembership.get(groupId).add(vPoint);
            
            pointToGroupMap.put(vPoint, groupId);
                        
        } else if ((groupId == null) && (vGroupId != null)) {

            groupId = vGroupId;

            groupMembership.get(groupId).add(uPoint);
            
            pointToGroupMap.put(uPoint, groupId);
            
        } else if ((groupId == null) && (vGroupId == null)) {
            
            groupId = Integer.valueOf(groupMembership.size());
            
            pointToGroupMap.put(uPoint, groupId);
            
            pointToGroupMap.put(vPoint, groupId);
            
            Set<PairInt> set = new HashSet<PairInt>();
            set.add(uPoint);
            set.add(vPoint);
            
            groupMembership.add(set);
                      
        } else if ((groupId != null) && (vGroupId != null)) {
            
            Integer uIdx = groupId;
            
            Integer vIdx = vGroupId;
            
            if (uIdx != vIdx) {
                if (uIdx > vIdx) {
                    Integer swap = uIdx;
                    uIdx = vIdx;
                    vIdx = swap;
                }
                
                // now uIdx is < vIdx
                
                /*
                key = group receiving merges, value = set of integers to merge into this
                key.  note that values are all greater than key
                */
                
                Set<Integer> group = mergeGroupsBeforePrune.get(uIdx);
                if (group == null) {
                    group = new HashSet<Integer>();
                    mergeGroupsBeforePrune.put(uIdx, group);
                }
                group.add(vIdx);
            }
        }
    }
    
    protected void process(PairInt uPoint) {
                
        Integer groupId = pointToGroupMap.get(uPoint);
        
        if (groupId == null) {
                        
            groupId = Integer.valueOf(groupMembership.size());
            
            pointToGroupMap.put(uPoint, groupId);
            
            Set<PairInt> set = new HashSet<PairInt>();
            set.add(uPoint);
            
            groupMembership.add(set);
        }
    }
    
    private void prune() {
        
        log.finest("number of groups before prune=" + groupMembership.size());
        
        //TODO: the data structures used could be written at the expense
        // of space complexity to reduce changes needed when group number
        // changes
        
        /*
         * [------] 0
         * [------] 1 <---- too few
         * [------] 2
         */
        // iterate backwards so can move items up without conflict with iterator
        for (int i = (groupMembership.size() - 1); i > -1; i--) {
            
            Set<PairInt> group = groupMembership.get(i);
            
            int count = group.size();
            
            log.finest("  group " + i + " has " + count 
                + " members before prune (min=" + minimumNumberInCluster + ")");
            if (count < minimumNumberInCluster) {
            
                // remove this group and move up all groups w/ index > i by one index
                for (int j = (i + 1); j < groupMembership.size(); j++) {
                    
                    int newGroupId = j - 1;
                                        
                    // update members in pointToGroupIndex
                    Set<PairInt> latest = groupMembership.get(j);
                    
                    for (PairInt p : latest) {
                        pointToGroupMap.put(p, Integer.valueOf(newGroupId));
                    }
                }
                
                Set<PairInt> removed = groupMembership.remove(i);
            }
        }
        
        log.finest("number of groups after prune=" + groupMembership.size());
    }

    private void merge() {
        
        log.finest("merge " + mergeGroupsBeforePrune.size() + " groups");
    
        /*
        merge values in tree first:
          9->10           is move 10 into group 9
          7->{8, 10}      is move 10 which is now 9) into group 7
        */
        Map<Integer, Set<Integer>> merged = new HashMap<Integer, Set<Integer>>();
        
        for (Integer moveToGroupIndex : mergeGroupsBeforePrune.descendingKeySet()) {
            
            Set<Integer> moveIndexes = mergeGroupsBeforePrune.get(moveToGroupIndex);
            
            // if any value in moveIndexes is a key in merged, add its values
            // to moveIndexes
            Set<Integer> add = new HashSet<Integer>();
            for (Integer moveGroupIndex : moveIndexes) {
                Set<Integer> set2 = merged.get(moveGroupIndex);
                if (set2 != null) {
                    add.addAll(set2);
                }
            }
            moveIndexes.addAll(add);
            
            for (Integer moveGroupIndex : moveIndexes) {
                        
                Set<PairInt> moveGroup = groupMembership.get(moveGroupIndex);
            
                for (PairInt p : moveGroup) {
                    pointToGroupMap.put(p, moveToGroupIndex);
                }

                groupMembership.get(moveToGroupIndex).addAll(moveGroup);

                groupMembership.get(moveGroupIndex).clear();
                
                Set<Integer> set2 = merged.get(moveGroupIndex);
                if (set2 == null) {
                    set2 = new HashSet<Integer>();
                    merged.put(moveGroupIndex, set2);
                }
                set2.add(moveToGroupIndex);
            }
        }
        
        mergeGroupsBeforePrune.clear();
    }

}
