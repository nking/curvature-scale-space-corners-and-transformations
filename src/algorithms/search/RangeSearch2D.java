package algorithms.search;

/**
 * adapted from 
 * https://www.cise.ufl.edu/class/cot5520fa09/CG_RangeTrees.pdf
 * which is the same as
 * http://www.cs.uu.nl/docs/vakken/ga/slides5b.pdf
 * 
 * @author nichole
 */
public class RangeSearch2D {
   
    /*
    Algorithm Build2DRangeTree(P)
        1. Construct the associated structure: 
           Build a binary search tree Tassoc on the 
           set Py of y-coordinates in P
           (can this be a threaded search tree?
            should add a split to the class)
         2. if P contains only one point
         3. then Create a leaf ν storing this point, and make
                Tassoc the associated structure of ν.
         4. else Split P into Pleft and Pright, the subsets ≤ and >
                the median x-coordinate xmid
         5. νleft ← Build2DRangeTree(Pleft)
         6. νright ← Build2DRangeTree(Pright)
         7. Create a node ν storing xmid, make νleft the left
                child of ν, make νright the right child of ν, and
                make Tassoc the associated structure of ν
         8. return ν

    Improve:
        Suppose we pre-sort P on y-coordinate, and whenever we split
        P into Pleft and Pright, we keep the y-order
        For a sorted set, the associated structure can be built in linear
        time.
    Then this construction is O(N * lg_2(N))
    
    */
   
    
    /*
    Algorithm 2DRangeQuery(T,[x : x0]×[y : y0])
       1. νsplit ←FindSplitNode(T,x,x0)
       2. if νsplit is a leaf
       3. then report the point stored at νsplit, if an answer
       4. else ν ← lc(νsplit)
       5. while ν is not a leaf
       6. do if x ≤ xν
       7. then 1DRangeQ(Tassoc(rc(ν)),[y : y0])
       8. ν ← lc(ν)
       9. else ν ← rc(ν)
      10. Check if the point stored at ν must be reported.
      11. Similarly, follow the path from rc(νsplit) to x0

    query time is O(log_2(n)^d)
    
    this can be improved to O(log_2(n)^(d-1)) with 
        fractional cascading
   
    Fractional Cascading:
    
    */
}
