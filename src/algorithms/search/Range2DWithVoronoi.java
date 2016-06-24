package algorithms.search;

/**
 *
 * @author nichole
 */
public class Range2DWithVoronoi {

    /*
     VoronoiFortunesSweep voronoi = 
         new VoronoiFortunesSweep();
        
     // O(N * log_2(N)) to build
     voronoi.generateVoronoi(x, y, xmin, xmax, ymin, ymax, 
        minDist);
    
     then have a list of edges containing site pairs
     and can retrieve the list of sites
    
     make a map w/ key site indexes and value = pther site index.
    
     The sites are already ordered by x and y so use
     those for input to anything else.
        
     -> to find nearest site, could make an alternate constructor
     for kd tree to accept sorted x,y.
     (see divide and conquer for closest points...there,
     the sorts for x and y partitions are minimized)
    
     when find closest site,
       then examine all of its edges
           the edges have site1 and site2
           -- for a range query with a distance limit,
              can use that here to exclude those
              with dist larger than limit.
           -- put the distances in a minheap and
              extract k of them.
              (note, that if a dist limit is used,
              and the limit is small, can use
              counting sort instead)
    
    */
}
