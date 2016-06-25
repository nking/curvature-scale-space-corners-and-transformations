package algorithms.search;

import algorithms.compGeometry.voronoi.VoronoiFortunesSweep;
import algorithms.compGeometry.voronoi.VoronoiFortunesSweep.GraphEdge;
import algorithms.compGeometry.voronoi.VoronoiFortunesSweep.Site;
import algorithms.imageProcessing.FixedSizeSortedVector;
import algorithms.util.PairFloat;
import algorithms.util.PairInt;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * A class to find the k nearest neighbors of a given
 * query point.  Internally, it uses a voronoi diagram
 * to find the neighbors of the nearest point to the
 * query point and returns the k closest to the query point.
 * 
  <pre>
  constructor, one time cost, runtime complexity:
      O(N * log_2(N))
  kNN query best runtime complexity:
      O(log_2(N)) + O(n_nearest_edges * log_2(k))
  kNN query worse runtime complexity:
      O(N) + O(n_nearest_edges * log_2(k))
  </pre>
  
 * @author nichole
 */
public class KNearestNeighbors {

    private VoronoiFortunesSweep voronoi = null;
    
    private Map<PairFloat, Set<Integer>> siteIndexesMap = null;
    
    // TODO: may want to swap this out for a faster means of
    //    finding the nearest (x[i], y[i]) for a query point
    private KDTreeFloat kdTree = null;
    
    public KNearestNeighbors(int[] x, int[] y) {
        init(x, y);
    }
    
    public KNearestNeighbors(float[] x, float[] y) {
        init(x, y);
    }
    
    private void init(int[] x, int[] y) {
        
        int n = x.length;
        
        float[] x2 = new float[n];
        float[] y2 = new float[n];
    
        for (int i = 0; i < n; ++i) {
            x2[i] = x[i];
            y2[i] = y[i];
        }
        
        init(x2, y2);
    }
    
    private void init(float[] x, float[] y) {
        
        int n = x.length;
        
        float xmin = Float.MAX_VALUE;
        float xmax = Float.MIN_VALUE;
        float ymin = Float.MAX_VALUE;
        float ymax = Float.MIN_VALUE;
        for (int i = 0; i < n; ++i) {
            float xp = x[i];
            float yp = y[i];
            if (xp < xmin) {
                xmin = xp;
            }
            if (xp > xmax) {
                xmax = xp;
            }
            if (yp < ymin) {
                ymin = yp;
            }
            if (yp > ymax) {
                ymax = yp;
            }
        }
        
        int minDist = 1;

        voronoi = new VoronoiFortunesSweep();

        // O(N * log_2(N)) to build
        voronoi.generateVoronoi(x, y, xmin, xmax, ymin, ymax,
            minDist);
        
        LinkedList<GraphEdge> edges = voronoi.getAllEdges();
        
        Site[] sites = voronoi.getSites();
        
        siteIndexesMap = new HashMap<PairFloat, Set<Integer>>();
    
        for (GraphEdge edge : edges) {
            int s1 = edge.site1;
            int s2 = edge.site2;
            
            PairFloat p1 = sites[s1].getCoord();
            PairFloat p2 = sites[s2].getCoord();
            
            Set<Integer> indexes = siteIndexesMap.get(p1);
            if (indexes == null) {
                indexes = new HashSet<Integer>();
                siteIndexesMap.put(p1, indexes);
            }
            indexes.add(s2);
            
            indexes = siteIndexesMap.get(p2);
            if (indexes == null) {
                indexes = new HashSet<Integer>();
                siteIndexesMap.put(p2, indexes);
            }
            indexes.add(s1);
        }
        assert(siteIndexesMap.size() == x.length);
        
        // retrieve the points from voronoi sites because they are sorted
        n = sites.length;
        float[] xp = new float[n];
        float[] yp = new float[n];
        for (int i = 0; i < n; ++i) {
            PairFloat p = sites[i].getCoord();
            xp[i] = p.getX();
            yp[i] = p.getY();
        }
        
        //O(N*lg_2(N)
        kdTree = new KDTreeFloat(xp, yp, true);
    }

    private float dist(float x, float y, PairFloat p) {
        float diffX = x - p.getX();
        float diffY = y - p.getY();
        return (float)Math.sqrt(diffX * diffX) + (diffY * diffY);
    }
    
    private class PairDist implements Comparable<PairDist>{
        PairFloat s1;
        float dist;
        @Override
        public int compareTo(PairDist other) {
            if (dist < other.dist) {
                return -1;
            } else if (dist > other.dist) {
                return 1;
            }
            return 0;
        }
    }
    
    public List<PairFloat> findNearest(int k, float x, float y) {
        return findNearest(k, x, y, Float.MAX_VALUE);
    }
    
    public List<PairFloat> findNearest(int k, float x, float y,
        float maxDistance) {
             
        // O(log_2(N) at best, but some extreme queries are O(N)
        // this returns one or equidistant multiple answers
        Set<PairFloat> nearest = kdTree.findNearestNeighbor(x, y);
      
        if (nearest == null) {
            return null;
        }
        
        // less than O(n_nearest)edges*log_2(k))
        FixedSizeSortedVector<PairDist> vec = 
            new FixedSizeSortedVector<PairDist>(k, PairDist.class);
        
        Site[] sites = voronoi.getSites();
        
        Set<PairFloat> added = new HashSet<PairFloat>();
        
        // sites are in siteIndexesMap
        for (PairFloat site : nearest) {
            
            float dist = dist(x, y, site);
            
            if (dist > maxDistance) {
                continue;
            }
            
            if (!added.contains(site)) {
                PairDist pd = new PairDist();
                pd.s1 = site;
                pd.dist = dist;
                vec.add(pd);
                
                added.add(site);
            }
            
            Set<Integer> siteIndexes = siteIndexesMap.get(site);
            
            if (siteIndexes == null) {
                throw new IllegalStateException("error in algorithm:"
                    + " voronoi diagram has no neighbors for "
                    + " (" + site.getX() + "," + site.getY() + ")");
            }
            
            for (Integer index2 : siteIndexes) {
                PairFloat site2 = sites[index2.intValue()].getCoord();
                dist = dist(x, y, site2);
                if (dist > maxDistance) {
                    continue;
                }
                if (!added.contains(site2)) {
                    PairDist pd = new PairDist();
                    pd.s1 = site2;
                    pd.dist = dist;
                    vec.add(pd);
                
                    added.add(site2);
                }
            }
        }
        
        List<PairFloat> output = new ArrayList<PairFloat>();
        PairDist[] a = vec.getArray();
        for (int i = 0; i < vec.getNumberOfItems(); ++i) {
            output.add(a[i].s1);
        }
        
        return output;
    }
    
}
