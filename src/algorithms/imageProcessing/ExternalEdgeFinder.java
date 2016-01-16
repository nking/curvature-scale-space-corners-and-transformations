package algorithms.imageProcessing;

import algorithms.compGeometry.PointInPolygon;
import algorithms.compGeometry.convexHull.GrahamScan;
import algorithms.compGeometry.convexHull.GrahamScanTooFewPointsException;
import algorithms.util.PairIntArray;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * class to find where edges in a set of edges from one list are external to the 
 * entire bounds of the reference edges.
 * 
 * @author nichole
 */
public class ExternalEdgeFinder {
    
    /**
     * find the edges in edgesToSearch that are outside the region occupied by
     * edgesAsReferenceFrame + buffer and return the result as a set of indexes
     * with respect to edgesToSearch.
     * @param edgesToSearch
     * @param edgesAsReferenceFrame
     * @param buffer
     * @return 
     * @throws algorithms.compGeometry.convexHull.GrahamScanTooFewPointsException 
     */
    public Set<Integer> findExteriorEdges(List<PairIntArray> edgesToSearch,
        List<PairIntArray> edgesAsReferenceFrame, int buffer) throws 
        GrahamScanTooFewPointsException {
        
        // find the convex hull of edgesAsReferenceFrame
        GrahamScan grahamScan = calculateConvexHull(edgesAsReferenceFrame);
        
        Set<Integer> skipSet = new HashSet<Integer>();
        
        if (grahamScan == null) {
            return skipSet;
        }
        
        // hull has start and end points repeated:
        float[] xHull = grahamScan.getXHull();
        float[] yHull = grahamScan.getYHull();
        expandHull(xHull, yHull, buffer);        
        
        PointInPolygon pointInPolygonChecker = new PointInPolygon();
        
        for (int edgeNumber = 0; edgeNumber < edgesToSearch.size(); edgeNumber++) {
            
            PairIntArray edge = edgesToSearch.get(edgeNumber);
                        
            for (int i = 0; i < edge.getN(); i++) {
                
                int x = edge.getX(i);
                int y = edge.getY(i);
                
                boolean isWithinCurve = pointInPolygonChecker.isInSimpleCurve(
                    x, y, xHull, yHull, xHull.length);
                
                if (!isWithinCurve) {
                    
                    // a point in edge is outside the hull
                    skipSet.add(Integer.valueOf(edgeNumber));
                    
                    break;
                }
            }
        }
        
        return skipSet;
    }
    
    /**
     * find the edges in edgesToSearch that are outside the region of the 
     * image frame.
     * 
     * @param edgesToSearch
     * @param imageReferenceWidth
     * @param imageReferenceHeight
     * @return 
     */
    public Set<Integer> findExteriorEdges(List<PairIntArray> edgesToSearch,
        int imageReferenceWidth, int imageReferenceHeight) {
        
        Set<Integer> skipSet = new HashSet<Integer>();
                
        for (int edgeNumber = 0; edgeNumber < edgesToSearch.size(); edgeNumber++) {
            
            PairIntArray edge = edgesToSearch.get(edgeNumber);
                        
            boolean isWithinFrame = true;
            
            for (int i = 0; i < edge.getN(); i++) {
                
                int x = edge.getX(i);
                int y = edge.getY(i);
                
                if ((x < 0) || (x > (imageReferenceWidth - 1))) {
                    isWithinFrame = false;
                    break;
                }
                
                if ((y < 0) || (y > (imageReferenceHeight - 1))) {
                    isWithinFrame = false;
                    break;
                }
            }
            
            if (!isWithinFrame) {
                skipSet.add(Integer.valueOf(edgeNumber));
            }
        }
        
        return skipSet;
    }
    
    private GrahamScan calculateConvexHull(List<PairIntArray> edges) throws 
        GrahamScanTooFewPointsException {
        
        int n2 = 0;
        for (PairIntArray edge : edges) {
            n2 += edge.getN();
        }
        if (n2 < 3) {
            return null;
        }
        int count = 0;
        float[] x2 = new float[n2];
        float[] y2 = new float[x2.length];
        for (PairIntArray edge : edges) {
            for (int i = 0; i < edge.getN(); i++) {
                x2[count] = edge.getX(i);
                y2[count] = edge.getY(i);
            }
        }
        
        GrahamScan scan = new GrahamScan();
        
        scan.computeHull(x2, y2);
        
        return scan;
    }

    private void expandHull(float[] xHull, float[] yHull, int buffer) {
        
        if (buffer == 0) {
            return;
        }
        
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        
        double[] xyCentroids = curveHelper.calculateXYCentroids(xHull, yHull);
        
        float xcen = (float)xyCentroids[0];
        float ycen = (float)xyCentroids[1];
        for (int i = 0; i < xHull.length; i++) {
            
            if (xHull[i] < xcen) {
                xHull[i] -= buffer;
            } else if (xHull[i] > xcen) {
                xHull[i] += buffer;
            }
            
            if (yHull[i] < ycen) {
                yHull[i] -= buffer;
            } else if (yHull[i] > ycen) {
                yHull[i] += buffer;
            }
        }    
    }
}
