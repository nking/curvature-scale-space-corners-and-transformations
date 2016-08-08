package algorithms.compGeometry;

import algorithms.compGeometry.voronoi.VoronoiFortunesSweep;
import algorithms.compGeometry.voronoi.VoronoiFortunesSweep.GraphEdge;
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
    private final NearestNeighbor2D np;
    
    /** as circles within points are searched, they are placed
    in processed
    */
    private final Set<PairInt> processed;
        
    //xMin, xMax, yMin, yMax
    private final int[] minMaxXY;
    
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
         
        this.np = new NearestNeighbor2D(boundary, 
            minMaxXY[1], minMaxXY[3]);
        
        processed = new HashSet<PairInt>(points.size());
        
    }
    
    public void findMedialAxis() {
        plotVoronoi();
        throw new UnsupportedOperationException("not yet impl");
    }
    
    private void plotVoronoi() {
        
        float xmin = Float.MAX_VALUE;
        float xmax = Float.MIN_VALUE;
        float ymin = Float.MAX_VALUE;
        float ymax = Float.MIN_VALUE;
        
        for (PairInt p : boundary) {
            float xp = p.getX();
            float yp = p.getY();
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

            if ((points.contains(p1) || processed.contains(p1))
                && (points.contains(p2) || processed.contains(p2))) {
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
        
        String filePath = plotter.writeFile(1000);
        System.out.println("wrote file=" + filePath);
        } catch (Throwable t) {
            
        }
    }

    // TODO: consider offering the edges
    // or points, whichever user prefers
    protected List<MedialAxisPoint> getMedAxisList() {
        return medAxisList;
    }
    
    public Set<PairInt> getMedialAxisPoints() {
        Set<PairInt> mAPs = new HashSet<PairInt>();
        for (MedialAxis1.MedialAxisPoint mp : medAxisList) {
            mAPs.add(mp.getCenter());
        }
        return mAPs;
    }

}
