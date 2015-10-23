package algorithms.compGeometry.clustering;

import algorithms.CountingSort;
import algorithms.MultiArrayMergeSort;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 *
 * @author nichole
 */
public class FixedDistanceGroupFinder {

    protected final boolean[][] visited;

    protected final int[] groupNumber;

    protected final List<Set<Integer>> groupList = new ArrayList<Set<Integer>>();

    protected List<Set<Integer>> sortedGroupList = null;

    protected final float[] x;

    protected final float[] y;

    public FixedDistanceGroupFinder(final float[] xPoints, final float[] yPoints) {

        if (xPoints == null) {
            throw new IllegalArgumentException("xPoints cannot be null");
        }
        if (yPoints == null) {
            throw new IllegalArgumentException("yPoints cannot be null");
        }

        int n = xPoints.length;

        visited = new boolean[n][];
        for (int i = 0; i < n; ++i) {
            visited[i] = new boolean[n];
        }

        groupNumber = new int[n];
        Arrays.fill(groupNumber, -1);

        x = xPoints;

        y = yPoints;

    }

    public int[] getGroupNumbers() {
        return groupNumber;
    }
    
    public int getNumberOfGroups() {
        return groupList.size();
    }
    
    /**
     * get a list of the x,y indexes in the given group number
     * @param groupNumber
     * @return
     */
    public Set<Integer> getGroupIndexes(final int groupNumber) {
        
        if (groupNumber < 0 || groupNumber > (groupList.size() - 1)) {
            throw new IllegalArgumentException("groupNumber is out of bounds");
        }
        
        Set<Integer> set = groupList.get(groupNumber);
        
        return set;
    }

    /**
     * get a list of the x,y indexes in each group in a list sorted by 
     * descending numbers of members (i.e., the largest group membership is 
     * index 0 of the returned list).
     * @return
     */
    public List<Set<Integer>> getDescendingSortGroupList() {
        
        if (sortedGroupList == null) {

            sortedGroupList = sortGroupList();
        }

        return sortedGroupList;
    }
    
    public void findGroupsOfPoints(final double maxPointSeparation) {
        
        double maxDistSq = maxPointSeparation * maxPointSeparation;
        
        java.util.Stack<Integer> stack = new java.util.Stack<Integer>();
        
        /*
        using a stack for the outer loop will increase the runtime complexity
        for this stage, but leads to un-fractured groups.
        */
        
        //O(N)
        for (int i = 0; i < x.length; ++i) {
            Integer index = Integer.valueOf(i);
            stack.add(index);
        }
        
        while (!stack.isEmpty()) {

            int uIdx = stack.pop().intValue();
            
            float uX = x[uIdx];
            float uY = y[uIdx];
                        
            for (int vIdx = 0; vIdx < x.length; ++vIdx) {
                
                if (visited[uIdx][vIdx]) {
                    continue;
                }
                
                visited[uIdx][vIdx] = true;
                visited[vIdx][uIdx] = true;
                    
                float vX = x[vIdx];
                float vY = y[vIdx];
                
                float diffX = uX - vX;
                float diffY = uY - vY;
                
                double distSq = (diffX * diffX) + (diffY * diffY);
                
                if (distSq <= maxDistSq) {
                                        
                    processPair(uIdx, vIdx);
                    
                    stack.add(Integer.valueOf(vIdx));
                }
            }
        }
    }

    protected void processPair(final int uIdx, final int vIdx) {
              
        int groupIdx = groupNumber[uIdx];
        
        if ((groupIdx != -1) && (groupNumber[vIdx] == -1)) {
                    
            groupList.get(Integer.valueOf(groupIdx)).add(Integer.valueOf(vIdx));
            
            groupNumber[vIdx] = groupIdx;
                        
        } else if ((groupIdx == -1) && (groupNumber[vIdx] != -1)) {

            groupIdx = groupNumber[vIdx];

            groupList.get(Integer.valueOf(groupIdx)).add(Integer.valueOf(uIdx));
            
            groupNumber[uIdx] = groupIdx;
            
        } else if ((groupIdx == -1) && (groupNumber[vIdx] == -1)) {
                        
            groupIdx = groupList.size();
            
            groupNumber[uIdx] = groupIdx;
            
            groupNumber[vIdx] = groupIdx;
            
            Set<Integer> set = new HashSet<Integer>();
            set.add(Integer.valueOf(uIdx));
            set.add(Integer.valueOf(vIdx));
            
            groupList.add(set);
                      
        } 
    }
    
    private List<Set<Integer>> sortGroupList() {

        int[] nMembers = new int[groupList.size()];
        int[] groupListIndexes = new int[groupList.size()];

        int nMax = Integer.MIN_VALUE;
        for (int i = 0; i < groupList.size(); ++i) {

            groupListIndexes[i] = i;

            nMembers[i] = groupList.get(i).size();

            if (nMembers[i] > nMax) {
                nMax = nMembers[i];
            }
        }

        if ((nMax > groupListIndexes.length) || (nMax > 10000000)) {
            MultiArrayMergeSort.sortByDecr(nMembers, groupListIndexes);
        } else {
            CountingSort.sortByDecr(nMembers, groupListIndexes, nMax);
        }

        List<Set<Integer>> sorted = new ArrayList<Set<Integer>>();

        for (int idx : groupListIndexes) {

            Set<Integer> set = groupList.get(idx);

            sorted.add(set);
        }

        return sorted;
    }

}
