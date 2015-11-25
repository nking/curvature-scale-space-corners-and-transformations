package algorithms.imageProcessing;

import algorithms.util.PairInt;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;

/**
 * a class to find contiguous pixels near one another having a value within
 * tolerance of the given value and using
 * logic to wrap around a set of values (0 to 360, for example, should see
 * 0 and 359 is being within a tolerance of '1' from one another). Note that
 * the minimum possible value is always 0.
 * 
 * @author nichole
 */
class DFSConnectedGroupsFinder3 extends AbstractDFSConnectedGroupsFinder {

    private enum State {
        INITIALIZED, GROUPS_FOUND, GROUPS_PRUNED, POST_GROUP_CORRECTED
    }
    
    private State state = null;
    
    public DFSConnectedGroupsFinder3() {
        
        minimumNumberInCluster = 1;
        
        state = State.INITIALIZED;
    }
    
    Logger constructLogger() {
        return Logger.getLogger(DFSConnectedGroupsFinder2.class.getName());
    }
    
    /**
     * near one another having a value within tolerance of the given value and 
     * using logic to
     * wrap around a set of values (0 to 360, for example, should see 0 and 359
     * ss being within a tolerance of '1' from one another). Note that the
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
    public void findConnectedPointGroups(Map<PairInt, Float> pointValueMap, 
        int maxValueForWrapAround, int toleranceInValue,
        int imageWidth, int imageHeight) {
        
        throw new UnsupportedOperationException("not yet implemented");
    }

    @Override
    protected void prune() {
        
        if (!state.equals(State.GROUPS_FOUND)) {
            throw new IllegalStateException(
                "findClustersIterative must be used before this");
        }
        
        super.prune();
        
        state = State.GROUPS_PRUNED;
    }
  
    private List<Set<PairInt>> findContiguousGroupsWithRangesWithinTolerance(
        Map<PairInt, Float> thetaMap,
        int maxValueForWrapAround, int toleranceInValue, 
        int imageWidth, int imageHeight) {
        
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
        use the median or several medians (e.g. quartiles) as dividers of values.
        
        The number of median dividers is 1 if the range < 2*tolerance 
        and > tolerance, and is 2 if the range is < 3*tolerance, etc.
        
        Then the "seed" values for subset membership are the mid points of
        the region above and below a divider and resulting groups will meet the
        goal.
        
        The optimal solution would be able to determine the "seed" values such
        that the fewest contiguous point groups meet the goal.
        
        /*
        range is 9
        tolerance=8, nDividers=int(range/tolerance) = 1
        1
        1
        3  
        4  median of median = 3 or 4 (avg=3.5)
        5
        6  
        6 median 
        6
        7
        8  median of median = 8
        9
        10
        10
        
        That provides seeds of values, but the result still needs to have
        the points contiguous and the fewest number of contiguous groups.
        
        The Voronoi algorithm is useful for 1D data, but this is 3D data
        and the other 2Ds require adjacency for all members within a
        seed's cell, so will not use Voronoi.
        
        The optimal value of median seeds is not obvious because of the need
        for adjacency for the membership of pixels.  In other words, it isn't
        a partition problem either because of the need for contiguous 
        membership.
        
        A solution which may not be optimal, would be to make connected groups
        for each seed ignoring whether the seed is already present in another
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
        
        It's possible to repeat that solution for all possible combination of
        seeds within each side of the median divider(s) and find the best
        solution among those. 
        For 4 seeds on one side and 3 on another, that would be 4*3 combinations.
        For the worse case scenario the range of the entire group would be 360
        and the tolerance would be very small so that the number of seeds
        possible on each side of a divider would factor to a large number of
        combinations n0*n1*n2*n3...
       
        Will use the non ideal solution as it does better separate truly 
        different value sets as defined by tolerance,
        but the solution might not provide the smallest number of contiguous
        groups.
        
        */

        throw new UnsupportedOperationException("not yet implemented");
    }

}
