package algorithms.compGeometry;

import algorithms.compGeometry.voronoi.VoronoiFortunesSweep;
import algorithms.compGeometry.voronoi.VoronoiFortunesSweep.GraphEdge;
import algorithms.imageProcessing.MiscellaneousCurveHelper;
import algorithms.misc.MiscMath;
import algorithms.search.NearestNeighbor2D;
import algorithms.util.PairInt;
import algorithms.util.PolygonAndPointPlotter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
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
    private final Set<PairInt> boundary;
    
    //xMin, xMax, yMin, yMax
    private final int[] minMaxXY;
    
    private List<GraphEdge> edges = null;
    
    private Set<PairInt> edgePoints = null;
    
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
                
        edges = findVoronoiInteriorEdges();
        
        //plotVoronoi();        
    }
    
    public void findMedialAxis() {
        
        plotVoronoi();

        NearestNeighbor2D np = new NearestNeighbor2D(boundary, 
            minMaxXY[1], minMaxXY[3]);
        
        throw new UnsupportedOperationException("not yet impl");
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
        
        VoronoiFortunesSweep voronoi = 
            new VoronoiFortunesSweep();
        
        voronoi.generateVoronoi(x, y, 
            xmin - 1, xmax + 1, ymin - 1, ymax + 1, 
            minDist);
        
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

            if (points.contains(p1) && points.contains(p2)) {
                output.add(edge);
            }
        }
        
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
        
        plotter.addPlot(x, y, xPolygon, yPolygon, "points");
        
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
        
        MiscellaneousCurveHelper curveHelper =
            new MiscellaneousCurveHelper();
        
        Set<PairInt> output = new HashSet<PairInt>();
        
        for (GraphEdge edge : edges) {
        
            int x1 = Math.round(edge.x1);
            int y1 = Math.round(edge.y1);
            int x2 = Math.round(edge.x2);
            int y2 = Math.round(edge.y2);

            curveHelper.createLinePoints(x1, y1, x2, y2,
                output);
        }
        
        this.edgePoints = output;
        
        return output;
    }
}
