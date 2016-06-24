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
     and have a sorted list of sites
    
     make a map w/ key site indexes and value = pther site index.
    
     The sites are already ordered by x and y so use
     those for input to anything else.
     
     to find nearest site, can use the KDTreeFloat.
     (build tree is O(N*lg_2(N))
     (NOTE: could consider transfoming data to integers with
      at least a gap of 1 between each point and then
      create a distance transform...to interpret results one
      needs to then scan the immediate 8 neigbhors to find
      minima for the answer...so the distance transform 
      method should only be used when N is larger than 1 << 9.)
     
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
