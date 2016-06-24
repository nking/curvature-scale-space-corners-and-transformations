package algorithms.search;

/**
 * a nearest neighbor based on the distance
 * transform.
 * 
 * @author nichole
 */
public class NearestNeighbor2DUsingDistTrans {
    
    /*
    - accepts integer arrays of x and y
    - transforms if needed to make sure there is at least
      a gap of 1 between each point.
      - stores the scale factpr if any
    - calculate distance transform.
    
    - nn queries:
      - transform query coords to scaled reference frame.
      - scan the 9 neighors to find minimum.  the neearest point
       should be off in that direction offset by that amount.
      - transform result back to original reference frame.
    */
}
