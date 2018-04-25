package packing;

/**
 * Class to return a set of points for the intersection of 2 point sets where 
 * the set of points are separated by a given interval in x and y, that is,
 * the bins are squares.
 * The two point sets can be irregular in shape.  
 * The solution is not guaranteed to be optimal.
 * 
 * The runtime complexity is expected to be:
 *     O(min(N_points) * lg_2(N_points))
 * 
 * Considering implementing these:
 * 
 * Option 1:
 *     resembles 2D strips:
 *     -- find the intersection of the point sets by iterating over the smallest
 *        point set.
 *     -- sort the intersection data structure by x, then y
 *     -- fill the intersection space with rectangles of given x and y size.
 *        -- use IntervaleRangeSearch for collision checks
 *        
 * Option 2:
 *     find the intersection of the point sets by iterating over the smallest
 *        point set.
 *     use a medial axis built from a voronoi diagram
 *     (see algorithms.compGeometry.MedialAxis), 
 *     then fill the space using order based upon the most connected medial 
 *     axis points and then the space connected to those assigned bins, 
 *     iteratively.
 *     It also uses IntervaleRangeSearch for collision checks.
 * 
 * Option 3:
 *    find the intersection of the point sets by iterating over the smallest
 *        point set.
 *     create ordered, connected boundary points of the
 *     intersection.
 *     then walk along the border, filling in bins
 *        at cell x,y spacings and adding the connected
 *        intervals to a stack iteratively to continue filling
 *        the intersection with bins.
 
 * 
 * @author nichole
 */
public class Intersection2DPacking {
    
}
