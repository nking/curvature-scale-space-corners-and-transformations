package algorithms.search;

import algorithms.compGeometry.voronoi.VoronoiFortunesSweep;
import algorithms.compGeometry.voronoi.VoronoiFortunesSweep.GraphEdge;
import algorithms.compGeometry.voronoi.VoronoiFortunesSweep.Site;
import algorithms.imageProcessing.FixedSizeSortedVector;
import algorithms.util.PairFloat;
import algorithms.util.PolygonAndPointPlotter;
import java.io.IOException;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;

/**
 * A class to find the k nearest neighbors of a given
 * query point.  Internally, it uses a voronoi diagram
 * to find the neighbors of the nearest point to the
 * query point and returns the k closest to the query 
 * point by adjacent voronoi sites.
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
    
    private float xmin = Float.MAX_VALUE;
    private float xmax = Float.MIN_VALUE;
    private float ymin = Float.MAX_VALUE;
    private float ymax = Float.MIN_VALUE;
        
    public KNearestNeighbors(int[] x, int[] y) {
        init(x, y);
    }
    
    public KNearestNeighbors(float[] x, float[] y) {
        init(x, y);
    }
    
    private void init(int[] x, int[] y) {
        
        if (x == null || x.length == 0 || y == null
            || y.length == 0 || (x.length != y.length)) {
            throw new IllegalArgumentException(
                "x and y cannot be null or empty "
                    + "and must be same lengths");
        }
        
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
        
        if (x == null || x.length == 0 || y == null
            || y.length == 0 || (x.length != y.length)) {
            throw new IllegalArgumentException(
                "x and y cannot be null or empty "
                    + "and must be same lengths");
        }
        if (x.length < 3) {
            throw new IllegalArgumentException("x and y "
            + " must be at least length 3 in size");
        }
        
        int n = x.length;
        
        xmin = Float.MAX_VALUE;
        xmax = Float.MIN_VALUE;
        ymin = Float.MAX_VALUE;
        ymax = Float.MIN_VALUE;
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
        
        int minDist = 0;

        voronoi = new VoronoiFortunesSweep();
    
        // O(N * log_2(N)) to build
        voronoi.generateVoronoi(x, y, xmin, xmax, ymin, ymax,
            minDist);
        
        LinkedList<GraphEdge> edges = voronoi.getAllEdges();
        
        Site[] sites = voronoi.getSites();
        
        siteIndexesMap = new HashMap<PairFloat, Set<Integer>>();
    
        assert (!edges.isEmpty());
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
            indexes.add(Integer.valueOf(s2));

            indexes = siteIndexesMap.get(p2);
            if (indexes == null) {
                indexes = new HashSet<Integer>();
                siteIndexesMap.put(p2, indexes);
            }
            indexes.add(Integer.valueOf(s1));
        }

        assert(sites.length == x.length);

        // points closer than minDist are not present,
        // so the map is possibly smaller than all points.
        // therefore, for this use of voronoi, need minDist=0.
        {//DEBUG
            if (siteIndexesMap.size() != x.length) {
                Logger.getLogger(this.getClass().getName())
                    .warning("siteMap.size=" + siteIndexesMap.size()
                    + " x.length=" + x.length);
            }
        }
        //assert(siteIndexesMap.size() == x.length);
        
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
       
        // O(log_2(N) at best, but some extreme queries are O(N).
        // nearest site(s). (if same distances, returns more than one).
        Set<PairFloat> nearest = kdTree.findNearestNeighbor(x, y);
      
        if (nearest == null) {
            return null;
        }
        
        /*
        a fixed vector of size k tracks the nearest and nearest
        adjacent.
        
        the search for k nearest continues in the adjacent sites
        as long as the adjacent site (whose neighbors should be
        searched) is nearer than the last item in the fixed vector.
        
        */
        
        // each fixed size vector comparison on insert is O(log_2(k))
        FixedSizeSortedVector<PairDist> vec = 
            new FixedSizeSortedVector<PairDist>(k, PairDist.class);
        
        Site[] sites = voronoi.getSites();
        
        Set<PairFloat> visited = new HashSet<PairFloat>();
        
        ArrayDeque<PairFloat> queue = new ArrayDeque<PairFloat>();
        queue.addAll(nearest);
        
        while (!queue.isEmpty()) {
            
            PairFloat site = queue.pop();
                 
            if (visited.contains(site)) {
                continue;
            }
            visited.add(site);
            
            float dist = dist(x, y, site);

            if (dist > maxDistance) {
                continue;
            }
            
            // if vec is not full or if site is closer than
            //  last full vec member, add site and add it's neighbors
            //  to queue
            
            int nV = vec.getNumberOfItems();
            if ((nV < k) || ((nV > 0) && 
                (dist < vec.getArray()[nV-1].dist))) {
            
                PairDist pd = new PairDist();
                pd.s1 = site;
                pd.dist = dist;
                vec.add(pd);                
            
                // add neighbors to queue
                
                Set<Integer> siteIndexes = siteIndexesMap.get(site);
            
                if (siteIndexes == null) {
                    throw new IllegalStateException("error in algorithm:"
                        + " voronoi diagram has no neighbors for "
                        + " (" + site.getX() + "," + site.getY() + ")");
                }
            
                for (Integer index2 : siteIndexes) {
                    PairFloat site2 = sites[index2.intValue()].getCoord();
                    queue.add(site2);
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
    
    public void debug(int fileNumber) throws IOException {
        
        Site[] sites = voronoi.getSites();
        int n = sites.length;
        float[] x = new float[n];
        float[] y = new float[n];
        for (int i = 0; i < n; ++i) {
            x[i] = sites[i].getCoord().getX();
            y[i] = sites[i].getCoord().getY();
        }
        
        LinkedList<GraphEdge> edges = voronoi.getAllEdges();
        
        PolygonAndPointPlotter plotter = 
            new PolygonAndPointPlotter(xmin, xmax, ymin, ymax);
        
        float[] xPolygon = null;
        float[] yPolygon = null;
        
        plotter.addPlot(x, y, xPolygon, yPolygon, "points");
        
        int n2 = edges.size();
        xPolygon = new float[2*n2];
        yPolygon = new float[2*n2];
        int count = 0;
        for (GraphEdge edge : edges) {
            xPolygon[count] = edge.x1;
            yPolygon[count] = edge.y1;
            xPolygon[count + 1] = edge.x2;
            yPolygon[count + 1] = edge.y2;
            count += 2;
        }
        plotter.addPlotWithLines(x, y, xPolygon, yPolygon, 
            "edges");
        String filePath = plotter.writeFile("debug_voron_" + fileNumber);
        System.out.println("wrote file=" + filePath);
    }
    
}
