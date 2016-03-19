package algorithms.imageProcessing;

import algorithms.imageProcessing.util.MiscStats;
import algorithms.util.PairInt;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.logging.Logger;

/**
 * a class to find contiguous pixels having a value within
 * tolerance of the given value and using
 * logic to wrap around a set of values (0 to 360, for example, should see
 * 0 and 359 is being within a tolerance of '1' from one another). Note that
 * the minimum possible value is always 0.
 * 
 * @author nichole
 */
class DFSConnectedGroupsFinder3 {
    
    /**
     * find contiguous pixels having a value within tolerance of the given value and 
     * using logic to
     * wrap around a set of values (0 to 360, for example, should see 0 and 359
     * as being within a tolerance of '1' from one another). Note that the
     * minimum possible value is always 0.
     * To correct for a group that has wandered from a total range of tolerance,
     * setTheCorrectForWandering to true
     *
     * @param pointValueMap
     * @param maxValueForWrapAround
     * @param toleranceInValue
     * @param imageWidth
     * @param imageHeight
     */
    public List<Set<PairInt>> findConnectedPointGroups(Map<PairInt, Float> pointValueMap, 
        int maxValueForWrapAround, int toleranceInValue,
        int imageWidth, int imageHeight) {
        
        return findContiguousGroupsWithRangesWithinTolerance(
            pointValueMap, maxValueForWrapAround, toleranceInValue, 
            imageWidth, imageHeight);
    }
  
    private List<Set<PairInt>> findContiguousGroupsWithRangesWithinTolerance(
        Map<PairInt, Float> thetaMap,
        int maxValueForWrapAround, int toleranceInValue, 
        int imageWidth, int imageHeight) {
        
        int n0 = thetaMap.size();
        
        if (n0 <= 1) {
            List<Set<PairInt>> list = new ArrayList<Set<PairInt>>();
            list.add(thetaMap.keySet());
            return list;
        }
        
        float[] values = new float[n0];
        int count = 0;
        for (Entry<PairInt, Float> entry : thetaMap.entrySet()) {
            values[count] = entry.getValue().floatValue();
            count++;
        }
        Arrays.sort(values);
        
        int[] startEndIndexes = MiscStats.determineStartEndIndexes(values, 
            maxValueForWrapAround, toleranceInValue);
        
        boolean isWrapAround = (startEndIndexes[0] > startEndIndexes[1]);
                
        float range = isWrapAround ?
            (values[startEndIndexes[1]] 
                + (maxValueForWrapAround - values[startEndIndexes[0]])) :
            (values[startEndIndexes[1]] - values[startEndIndexes[0]]);
        
        if (range == toleranceInValue) {
            List<Set<PairInt>> list = new ArrayList<Set<PairInt>>();
            list.add(thetaMap.keySet());
            return list;
        }
                      
        int nDividers = (int)(range/toleranceInValue);
        
        int nSeeds = nDividers + 1;
        
        float binWidth = (range/(float)nSeeds)/2.f;
                
        float[] seeds = new float[nSeeds];
        for (int i = 0; i < nSeeds; ++i) {
            float v = values[startEndIndexes[0]] + (((2 * i) + 1) * binWidth);
            if (v > maxValueForWrapAround) {
                v -= maxValueForWrapAround;
            }
            seeds[i] = v;
        }
        
        assert(seeds[seeds.length -1] != 0);
    
        int n2 = 0;
        
        List<Set<PairInt>> output = new ArrayList<Set<PairInt>>();
        
        Map<PairInt, List<Integer>> pointGroupMap = new HashMap<PairInt, List<Integer>>();
        
        for (int i = 0; i < seeds.length; ++i) {
            
            float seed = seeds[i];
            
            ContiguousFinder finder = new ContiguousFinder();
            
            finder.findConnectedPointGroups(thetaMap, seed,
                maxValueForWrapAround, toleranceInValue, imageWidth, 
                imageHeight);
            
            int n = finder.getNumberOfGroups();
            
            for (int j = 0; j < n; ++j) {
                
                Set<PairInt> group = finder.getXY(j);
              
                n2 += group.size();
                
                Integer groupIndex = Integer.valueOf(output.size());
                
                output.add(group);
                
                for (PairInt p : group) {
                    List<Integer> list = pointGroupMap.get(p);
                    if (list == null) {
                        list = new ArrayList<Integer>();
                        pointGroupMap.put(p, list);
                    }
                    list.add(groupIndex);
                }
            }
        }
        
        assert(n2 >= n0);
        
        resolveMultipleMappings(output, pointGroupMap, thetaMap);
        
        return output;
        /*
        Because groups that were found by DFSConnectedGroupsFinder2
        were found by an offset from each pixel's value,
        those groups may have a range of min value to max value that is larger
        than tolerance.
        
        This method finds the "seeds" that result in contiguous subsets 
        whose min to max values are less than or equal to tolerance.
        
        Ideally, the algorithm finds the fewest subsets that each have a 
        range <= tolerance.  the total subsets equal the total membership
        of the super set.
        
        a solution:
           create seeds that are within tolerance of one another and span
        the range
        
        The optimal solution would be able to determine the "seed" values such
        that the fewest contiguous point groups meet the goal.
        
        /*
        The Voronoi algorithm is useful for 1D data, but this is 3D data
        and the other 2Ds require adjacency for all members within a
        seed's cell, so will not use Voronoi.
        
        The optimal value of median seeds is not obvious because of the need
        for adjacency for the membership of pixels.  In other words, it isn't
        a partition problem either because of the need for contiguous 
        membership.
        
        A solution which may not be optimal, would be to make connected groups
        for each seed ignoring whether the point is already present in another
        seed group (but noting dual membership).
        Then to resolve the best membership for the ambiguous points that are
        present in more than one group.
           -- the ambiguous points can choose the seed they are closer to in
              value, and ties can be broken by avg rgb difference of the group
              and the ambiguous point, and further ties can be broken by
              proximity to group centroid.
              -- a caveat is to avoid breaking the connection for an adjacent
                 pixel which does not have another group it could be a member of.
        This solution thus far would be nSeeds times 
        semi-linear in npoints of group, then less than semi-linear for resolving
        ambigous point memberships.
        
        It's possible to repeat that solution for all possible combinations of
        seeds within each side of the median divider(s) and find the best
        solution among those. 
        For 4 seeds on one side and 3 on another, that would be 4*3 combinations.
        For the worse case scenario the range of the entire group would be 360
        and the tolerance would be very small so that the number of seeds
        possible on each side of a divider would factor to a large number of
        combinations n0*n1*n2*n3...
       
        leaning towards using the non ideal solution as it does separate 
        different value sets as defined by tolerance,
        but the solution might not provide the smallest number of contiguous
        groups.
        
        looking at graph cuts algorithms next...
        */

    }

    private void resolveMultipleMappings(List<Set<PairInt>> output, 
        Map<PairInt, List<Integer>> pointGroupMap, 
        Map<PairInt, Float> thetaMap) {
        
        //TODO: consider using rgb to resolve best group for ambiguous matchings
        // and a tie in thetaAvg difference
        
        float[] thetaAvg = new float[output.size()];
        PairInt[] centroids = new PairInt[output.size()];
        Arrays.fill(thetaAvg, -1);
        
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        
        // for first implementation, do not perform check for whether moving
        // a point out of a group disconnects another point.
        for (Entry<PairInt, List<Integer>> entry : pointGroupMap.entrySet()) {
            List<Integer> indexes = entry.getValue();
            if (indexes.size() > 1) {
                
                PairInt p = entry.getKey();
                
                // choose between groups
                float minDiffTheta = Float.MAX_VALUE;
                int minDiffThetaIdx = -1;
                int minDiffThetaDistSq = Integer.MAX_VALUE;
                for (int i = 0; i < indexes.size(); ++i) {
                    int idx = indexes.get(i).intValue();
                    if (thetaAvg[idx] == -1) {
                        thetaAvg[idx] = calculateAvg(output.get(idx), thetaMap);
                    }
                    if (centroids[idx] == null) {
                        double[] xyCen = curveHelper.calculateXYCentroids(output.get(idx));
                        centroids[idx] = new PairInt((int)Math.round(xyCen[0]),
                            (int)Math.round(xyCen[1]));
                    }
                    float theta = thetaMap.get(p).floatValue();
                    float diffT = Math.abs(thetaAvg[idx] - theta);
                    int diffX = p.getX() - centroids[idx].getX();
                    int diffY = p.getY() - centroids[idx].getY();
                    int distSq = (diffX * diffX) + (diffY * diffY);
                    if ((diffT < minDiffTheta) ||  ((diffT == minDiffTheta) &&
                        (distSq < minDiffThetaDistSq))) {
                        minDiffTheta = diffT;
                        minDiffThetaIdx = idx;
                        minDiffThetaDistSq = distSq;
                    }
                }
                for (int i = (indexes.size() - 1); i > -1 ; --i) {
                    int idx = indexes.get(i).intValue();
                    if (idx == minDiffThetaIdx) {
                        continue;
                    }
                    output.get(idx).remove(p);
                    indexes.remove(indexes.get(i));
                }
int z = 1;                
            }
        }
int z = 1;
    }

    private float calculateAvg(Set<PairInt> points, Map<PairInt, Float> thetaMap) {
        
        float sum = 0;
        for (PairInt p : points) {
            sum += thetaMap.get(p).floatValue();
        }
        sum /= (float)points.size();
        
        return sum;
    }

    protected static class ContiguousFinder extends AbstractDFSConnectedGroupsFinder {

        public ContiguousFinder() {
            minimumNumberInCluster = 1;        
        }
    
        public void findConnectedPointGroups(Map<PairInt, Float> pointValueMap,
            float value, int maxValueForWrapAround, int toleranceInValue,
            int imageWidth, int imageHeight) {

            findClustersIterative(pointValueMap, value, maxValueForWrapAround,
                toleranceInValue, imageWidth, imageHeight);

            prune();
        }
        
        private void findClustersIterative(Map<PairInt, Float> pointValueMap,
            float value, int maxValueForWrapAround, int toleranceInValue,
            int imageWidth, int imageHeight) {

            Set<PairInt> visited = new HashSet<PairInt>();

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

                boolean similar = false;
                if (Math.abs(uValue - value) <= toleranceInValue) {
                    similar = true;
                }
                if (!similar) {
                    if (Math.abs(uValue - (value - maxValueForWrapAround)) <= toleranceInValue) {
                        similar = true;
                    }
                }
                if (!similar) {
                    if (Math.abs(value - (uValue - maxValueForWrapAround)) <= toleranceInValue) {
                        similar = true;
                    }
                }
                if (!similar) {
                    continue;
                } 

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

                        float vValue = pointValueMap.get(vPoint).floatValue();
                        similar = false;
                        if (Math.abs(vValue - value) <= toleranceInValue) {
                            similar = true;
                        }
                        if (!similar) {
                            if (Math.abs(vValue - (value - maxValueForWrapAround)) <= toleranceInValue) {
                                similar = true;
                            }
                        }
                        if (!similar) {
                            if (Math.abs(value - (vValue - maxValueForWrapAround)) <= toleranceInValue) {
                                similar = true;
                            }
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

        @Override
        Logger constructLogger() {
            return Logger.getLogger(DFSConnectedGroupsFinder2.class.getName());
        }
    
    }
}
