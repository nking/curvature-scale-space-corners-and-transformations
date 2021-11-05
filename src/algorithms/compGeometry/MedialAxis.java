package algorithms.compGeometry;

import thirdparty.voronoi.VoronoiFortunesSweep;
import thirdparty.voronoi.VoronoiFortunesSweep.*;
import algorithms.imageProcessing.BresenhamsLine;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.PostLineThinnerCorrections;
import algorithms.misc.MiscMath;
import algorithms.msts.PrimsMST;
import algorithms.search.NearestNeighbor2D;
import algorithms.util.PairInt;
import algorithms.util.PolygonAndPointPlotter;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.list.TIntList;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.logging.Logger;

/**
   A class to create a 2D medial axis given shape points
   and the boundary.
       from https://en.wikipedia.org/wiki/Medial_axis
       "the medial axis of a subset S which is bounded 
        by planar curve C is the locus of the centers 
        of circles that are tangent to curve C in two 
        or more points."

   The class uses a Voronoi diagram, removes the 
   boundary connecting edges, and then tries to
   remove points that are not medial axis points
   that are artifacts of small bumps in the
   shape boundary.

   The class is a work in progress, so there may still
   be non-medial axis points that have not been removed
   yet.

 * @author nichole
 */
public class MedialAxis {
    
    private Logger log = Logger.getLogger(this.getClass().getName());
    
    private final Set<PairInt> points;
    private Set<PairInt> boundary = null;
    
    /**
     * if a line thinner was applied, the points
     * removed from the boundary are kept here.
     * they are needed when excluding "exterior"
     * points from the medial axis edges.
     */
    private final Set<PairInt> removedPoints = 
        new HashSet<PairInt>();
    
    //xMin, xMax, yMin, yMax
    private final int[] minMaxXY;
    
    private List<GraphEdge> edges = null;
    
    private Set<PairInt> edgePoints = null;
    
    private boolean applyLT = false;
    
    /**
     * constructor containing all points in the area
     * and the bounding points.
     * 
     * @param shapePoints
     * @param boundaryPoints 
     */
    public MedialAxis(final Set<PairInt> shapePoints,
        final Set<PairInt> boundaryPoints) {

        this.points = new HashSet<PairInt>(shapePoints);
        
        this.boundary = new HashSet<PairInt>(boundaryPoints);
            
        points.removeAll(boundary);
        
        minMaxXY = MiscMath.findMinMaxXY(boundary);         
    }
    
    public void setToApplyLineThinner() {
        applyLT = true;
    }
    
    /**
     * find the medial axis, but with the warning that 
     * it may contain extra spikes out to false medial
     * axis points that are artifacts of bumps in the
     * boundary of the shape.
     * The artifacts are not harmful for the main reason
     * the class was built, so this faster method is
     * offered as an option.
     * The runtime complexity is purely that of Voronoi
     * Fortune Sweep, O(N * log_2(N)) where N is the
     * number of boundary points.
     * 
     */
    public void fastFindMedialAxis() {
    
        if (edges != null) {
            throw new IllegalStateException(
                "find... has already been inboked");
        }
        
        if (applyLT) {
            applyLineThinner();
        }
        
        edges = findVoronoiInteriorEdges();
        
        //plotVoronoi();        
    }
    
    public void findMedialAxis() {
        
        if (edges != null) {
            throw new IllegalStateException(
                "find... has already been inboked");
        }
        
        applyLineThinner();
        
        edges = findVoronoiInteriorEdges2();
                
        //plotVoronoi();

    }
    
    private List<GraphEdge> findVoronoiInteriorEdges() {
        
        float xmin = minMaxXY[0];
        float xmax = minMaxXY[1];
        float ymin = minMaxXY[2];
        float ymax = minMaxXY[3];
        
        int n = boundary.size();
        float[] x = new float[n];
        float[] y = new float[n];
        
        int count = 0;
        for (PairInt p : boundary) {
            float xp = p.getX();
            float yp = p.getY();
            x[count] = xp;
            y[count] = yp;
            count++;
        }
                
        int minDist = 0;
        int offset = 2;
        
        VoronoiFortunesSweep voronoi = new VoronoiFortunesSweep();
        
        voronoi.generateVoronoi(x, y, 
            xmin - offset, xmax + offset, 
            ymin - offset, ymax + offset, 
            minDist);
        
        //voronoi.plot(1234);
                
        LinkedList<GraphEdge> edges = voronoi.getAllEdges();
        
        List<GraphEdge> output = new ArrayList<GraphEdge>();
       
        count = 0;
        for (GraphEdge edge : edges) {
            int x1 = Math.round(edge.x1);
            int y1 = Math.round(edge.y1);
            int x2 = Math.round(edge.x2);
            int y2 = Math.round(edge.y2);

            PairInt p1 = new PairInt(x1, y1);
            PairInt p2 = new PairInt(x2, y2);
            if (p1.equals(p2)) {
                continue;
            }
        
            if (removedPoints.contains(p1) ||
                removedPoints.contains(p2) ||
                boundary.contains(p1) ||
                boundary.contains(p2)) {
                continue;
            }
            
            if (points.contains(p1) && points.contains(p2)) {
                output.add(edge);
            }
        }
        
        if (output.isEmpty() && boundary.size() > 3 &&
            boundary.size() < 12) {
            
            // small space, and the thinning stage before
            // using this class may have removed interior
            // points so the points.contains(p) fails.
            // this does the more expensive point in polygon
            // test too.
            // TODO: find a fast correction for the line thinning
            // subsequent removal of points from shape points.
            
            x = Arrays.copyOf(x, x.length + 1);
            y = Arrays.copyOf(y, y.length + 1);
            x[x.length - 1] = x[0];
            y[y.length - 1] = y[0];
            
            PointInPolygon pip = new PointInPolygon();
            for (GraphEdge edge : edges) {
                int x1 = Math.round(edge.x1);
                int y1 = Math.round(edge.y1);
                int x2 = Math.round(edge.x2);
                int y2 = Math.round(edge.y2);

                PairInt p1 = new PairInt(x1, y1);
                PairInt p2 = new PairInt(x2, y2);
                if (p1.equals(p2)) {
                    continue;
                }
                
                // if p1 is in points or is interior to boundary
                // and same for p2, can keep it
                if ((points.contains(p1) ||
                    pip.isInSimpleCurve(x1, y1, x, y, x.length)) && 
                    (points.contains(p2) ||
                    pip.isInSimpleCurve(x2, y2, x, y, x.length))) {
                    
                    output.add(edge);
                }
            }
        }
        
        return output;
    }
    
    private List<GraphEdge> findVoronoiInteriorEdges2() {
               
        //TODO: improve this:
        Set<PairInt> c = findBoundaryProblems();
        
        NearestNeighbor2D nn = new NearestNeighbor2D(
            c, minMaxXY[1], minMaxXY[3]);
        
        float xmin = minMaxXY[0];
        float xmax = minMaxXY[1];
        float ymin = minMaxXY[2];
        float ymax = minMaxXY[3];
        
        int n = boundary.size();
        float[] x = new float[n];
        float[] y = new float[n];
        
        int count = 0;
        for (PairInt p : boundary) {
            float xp = p.getX();
            float yp = p.getY();
            x[count] = xp;
            y[count] = yp;
            count++;
        }
                
        int minDist = 0;
        
        VoronoiFortunesSweep voronoi = 
            new VoronoiFortunesSweep();
        
        voronoi.generateVoronoi(x, y, 
            xmin - 1, xmax + 1, ymin - 1, ymax + 1, 
            minDist);
        
        LinkedList<GraphEdge> edges = voronoi.getAllEdges();
        
        List<GraphEdge> output = new ArrayList<GraphEdge>();
       
        Site[] sites = voronoi.getSites();
        
        /*
        need to store edge points and make an adjacency
        map for them.
        while storing the edge points, need to find the
        edge points which are very near the set rm.
        
        then will make an mst for the adjacency map.
        
        then will remove the rm points and the branch
        they are on up until they reach a parent with
        another child.
        */
        
        TIntObjectMap<TIntIntMap>
            adjCostMap = new TIntObjectHashMap<TIntIntMap>();
       
        /*
        edge: v1, v2
           each vertex in edge gets a vertex index
              that may already exist
        each edge is stored in map w/ key=pairint(v1,v2)
           where v1<v2.
        
        rm are the verex indexes to remove
        */
        TIntSet rm = new TIntHashSet();
        TObjectIntMap<PairInt> vertexIndexes 
            = new TObjectIntHashMap<PairInt>();
        Map<PairInt, GraphEdge> vertexEdgeMap = 
            new HashMap<PairInt, GraphEdge>();
        
        count = 0;
        for (GraphEdge edge : edges) {
            int x1 = Math.round(edge.x1);
            int y1 = Math.round(edge.y1);
            int x2 = Math.round(edge.x2);
            int y2 = Math.round(edge.y2);

            PairInt p1 = new PairInt(x1, y1);
            PairInt p2 = new PairInt(x2, y2);
           
            if (p1.equals(p2)) {
                continue;
            }
            if (removedPoints.contains(p1) ||
                removedPoints.contains(p1) ||
                boundary.contains(p1) ||
                boundary.contains(p2)) {
                continue;
            }
            
            if (points.contains(p1) && points.contains(p2)) {

                int idx1;
                if (vertexIndexes.containsKey(p1)) {
                    idx1 = vertexIndexes.get(p1);
                } else {
                    idx1 = vertexIndexes.size();
                    vertexIndexes.put(p1, idx1);
                }
                        
                Set<PairInt> nearest1 = 
                    nn.findClosest(x1, y1, 3);
                if (!nearest1.isEmpty()) {
                    rm.add(idx1);
                    //System.out.println("rm: " + p1);
                }
                
                int idx2;
                if (vertexIndexes.containsKey(p2)) {
                    idx2 = vertexIndexes.get(p2);
                } else {
                    idx2 = vertexIndexes.size();
                    vertexIndexes.put(p2, idx2);
                }
                Set<PairInt> nearest2 = 
                    nn.findClosest(x2, y2, 3);
                if (!nearest2.isEmpty()) {
                    rm.add(idx2);
                    //System.out.println("rm: " + p2);
                }
                
                PairInt eKey;
                if (idx1 < idx2) {
                    eKey = new PairInt(idx1, idx2);
                } else {
                    eKey = new PairInt(idx2, idx1);
                }
                
                if (vertexEdgeMap.containsKey(eKey)) {
                    continue;
                }
                vertexEdgeMap.put(eKey, edge);
                
                TIntIntMap map = adjCostMap.get(idx1);
                if (map == null) {
                    map = new TIntIntHashMap();
                    adjCostMap.put(idx1, map);
                }
                map.put(idx2, 1);
                
                map = adjCostMap.get(idx2);
                if (map == null) {
                    map = new TIntIntHashMap();
                    adjCostMap.put(idx2, map);
                }
                map.put(idx1, 1);
            }
        }
        
       // System.out.println("nVertexes=" + vertexIndexes.size());
       // System.out.println("nEdges=" + vertexEdgeMap.size());
        int maxCost = PrimsMST.maxEdgeCost(adjCostMap);
        PrimsMST mst = new PrimsMST();
        mst.calculateMinimumSpanningTree(adjCostMap, maxCost);
        
        int[] prev = mst.getPredeccessorArray();
        Map<Integer, LinkedList<Integer>> treeMap = mst.makeTreeFromPrev();
                
        TIntSet rm2 = new TIntHashSet();
        
        Set<GraphEdge> rmEdges = new HashSet<GraphEdge>();
        TIntIterator iter = rm.iterator();
        while (iter.hasNext()) {
            
            int idx = iter.next();
                                    
            // keep walking up tree until a parent has children.n>2
            while (true) {
                
                if (rm2.contains(idx)) {
                    break;
                }
                
                int prevIdx = prev[idx];
                
                PairInt rmKey;
                if (idx < prevIdx) {
                    rmKey = new PairInt(idx, prevIdx);
                } else {
                    rmKey = new PairInt(prevIdx, idx);
                }

                rmEdges.add(vertexEdgeMap.get(rmKey));
                
                LinkedList<Integer> pC = treeMap.get(prevIdx);
                if (pC == null) {
                    break;
                }
                //boolean removed = pC.remove(idx);
                //rm2.add(idx);
                //assert(removed);
                if (pC.size() > 1) {
                    break;
                }
                idx = prevIdx;
            }
        }
        
        //System.out.println("nEdges to remove=" + rmEdges.size());
        
        Iterator<Entry<PairInt, GraphEdge>> iter2 = 
            vertexEdgeMap.entrySet().iterator();
        while (iter2.hasNext()) {
            Entry<PairInt, GraphEdge> entry = iter2.next();
            GraphEdge edge = entry.getValue();
            if (!rmEdges.contains(edge)) {
                output.add(edge);
            }
        }
        
        System.out.println("nEdges=" + output.size());
        
        return output;
    }
    
    private void plotVoronoi() {
        
        float xmin = minMaxXY[0];
        float xmax = minMaxXY[1];
        float ymin = minMaxXY[2];
        float ymax = minMaxXY[3];
        
        int n = boundary.size();
        float[] x = new float[n];
        float[] y = new float[n];
        
        int count = 0;
        for (PairInt p : boundary) {
            float xp = p.getX();
            float yp = p.getY();
            x[count] = xp;
            y[count] = yp;
            count++;
        }
               
        if (edges == null) {
            fastFindMedialAxis();
        }
                
        try {
        PolygonAndPointPlotter plotter = 
            new PolygonAndPointPlotter(xmin - 1, xmax + 1, 
                ymin - 1, ymax + 1);
        
        float[] xPolygon = null;
        float[] yPolygon = null;
                
        n = edges.size();
        xPolygon = new float[2*n];
        yPolygon = new float[2*n];
        count = 0;
        for (GraphEdge edge : edges) {
            float x1 = edge.x1;
            float y1 = edge.y1;
            float x2 = edge.x2;
            float y2 = edge.y2;
            
            xPolygon[count] = x1;
            yPolygon[count] = y1;
            xPolygon[count + 1] = x2;
            yPolygon[count + 1] = y2;
            count += 2;
        }
        plotter.addPlotWithLines(x, y, xPolygon, yPolygon, 
            "edges");
        
        count = 0;
        for (GraphEdge edge : edges) {
            int x1 = Math.round(edge.x1);
            int y1 = Math.round(edge.y1);
            int x2 = Math.round(edge.x2);
            int y2 = Math.round(edge.y2);

            PairInt p1 = new PairInt(x1, y1);
            PairInt p2 = new PairInt(x2, y2);

            if (points.contains(p1) && points.contains(p2)) {
                xPolygon[count] = x1;
                yPolygon[count] = y1;
                xPolygon[count + 1] = x2;
                yPolygon[count + 1] = y2;
                count += 2;
            }
        }
        xPolygon = Arrays.copyOf(xPolygon, count);
        yPolygon = Arrays.copyOf(yPolygon, count);
        
        plotter.addPlotWithLines(x, y, xPolygon, yPolygon, 
            "edited for medial axes");
        
        //-----
        Set<PairInt> pts = getMedialAxisPoints();
        n = pts.size();
        x = new float[n];
        y = new float[n];
        count = 0;
        for (PairInt p : pts) {
            float xp = p.getX();
            float yp = p.getY();
            x[count] = xp;
            y[count] = yp;
            count++;
        }
        plotter.addPlotWithLines(x, y, xPolygon, yPolygon, 
            "med axis pts");
        
        String filePath = plotter.writeFile(1000);
        System.out.println("wrote file=" + filePath);
        } catch (Throwable t) {
            
        }
    }

    /**
     * return the results as the graph edges.
     * The output is copied so modifications
     * won't affect the instance variables.
     * @return 
     */
    protected List<GraphEdge> getMedAxisAsEdges() {
        return new ArrayList<GraphEdge>(edges);
    }
    
    /**
     * return the results as undirected point sets.
     * Note that the method runtime complexity is
     * O(N_edges) because it creates points between
     * the edge endpoints.
     * @return 
     */
    public Set<PairInt> getMedialAxisPoints() {
        
        if (edgePoints != null) {
            return edgePoints;
        }
       
        Set<PairInt> output = new HashSet<PairInt>();
        
        for (GraphEdge edge : edges) {
        
            int x1 = Math.round(edge.x1);
            int y1 = Math.round(edge.y1);
            int x2 = Math.round(edge.x2);
            int y2 = Math.round(edge.y2);

            BresenhamsLine.createLinePoints(x1, y1, x2, y2,
                output);
        }
        
        this.edgePoints = output;
        
        return output;
    }
    
    private Set<PairInt> findBoundaryProblems() {
                
        PostLineThinnerCorrections pltc = new 
            PostLineThinnerCorrections();
        
        return pltc.findBoundaryPattern(boundary, minMaxXY[1], minMaxXY[3]);
        
    }

    private void applyLineThinner() {
        
        Set<PairInt> b = new HashSet<PairInt>(boundary);
        
        ImageProcessor imp = new ImageProcessor();
        imp.applyThinning(b, minMaxXY[1] + 1, 
            minMaxXY[3] + 1);
        
        //SpurRemover spurRm = new SpurRemover();
        //spurRm.remove(b, minMaxXY[1] + 3, 
        //    minMaxXY[3] + 3);
        
        /*
        ZhangSuenLineThinner lt = new ZhangSuenLineThinner();
        lt.applyLineThinner(b, 
            minMaxXY[0] - 1, minMaxXY[1] + 1, 
            minMaxXY[2] - 1, minMaxXY[3] + 1);
        */
            
        // store the removed points
        removedPoints.addAll(boundary);
        removedPoints.removeAll(b);
        
        System.out.println("line thinning removed " +
            removedPoints.size() + " from the boundary");
        
        boundary = b;
    }

    public Set<PairInt> getBoundary() {
        return boundary;
    }
}
