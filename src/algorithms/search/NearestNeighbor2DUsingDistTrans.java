package algorithms.search;

/**
 * a nearest neighbor based on the distance
 * transform.
 * 
 * @author nichole
 */
public class NearestNeighbor2DUsingDistTrans {
    
    /*
    
    NOTE: can only be used to search for a point
    that does not exist in the original set of points.
    
    - accepts integer arrays of x and y
    - transforms if needed to make sure there is at least
      a gap of 1 between each point.
      - stores the scale factpr if any
      - *NOTE that the runtime complexity for calculating the
         distance transform is then max x times max y,
         so this method should not be used on very large axes data.
    - calculate distance transform.
    
    - nn queries:
      - transform query coords to scaled reference frame.
      - scan the 9 neighors to find minimum.  the neearest point
       should be off in that direction offset by that amount.
      - transform result back to original reference frame.
    */
}
