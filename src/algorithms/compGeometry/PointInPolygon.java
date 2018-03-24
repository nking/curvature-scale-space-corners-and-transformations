package algorithms.compGeometry;

import algorithms.misc.MiscMath;

/**
 * test whether a point is within a polygon.
 *
 * adapted from http://rosettacode.org/wiki/Ray-casting_algorithm impl for Go that
 * counts the number of times a ray intersects the edges of the polygon starting
 * from the point and extending to higher x.  This class also does a quick
 * check for whether the point is in the line segment first.
 *
 * @author nichole
 */
public class PointInPolygon {

    public PointInPolygon() {
    }

    float eps = 0.00001f;

    /**
     * given a polygon (xPolygon, yPolygon) that has the same start and end
     * points, determine whether the point (xpt, ypt) is in the polygon.
     * 
     * @param xPt
     * @param yPt
     * @param xPolygon x values of a polygon in which each x should be smaller
     * than Integer.MAX_VALUE so that cross products do not overflow int type.
     * @param yPolygon y values of a polygon in which each y should be smaller
     * than Integer.MAX_VALUE so that cross products do not overflow int type.
     * @param nPolygonPoints
     * @return whether (xPt, yPt) is on or within {xPolygon, yPolygon}
     */
    public boolean isInSimpleCurve(int xPt, int yPt, int[] xPolygon,
        int[] yPolygon, int nPolygonPoints) {
        
        int xMax = MiscMath.findMax(xPolygon, nPolygonPoints) + 1;
        
        int n2 = 0;
        int n3 = 0;
        int n0 = 0;
        
        int sumIntersectingRays = 0;
        
        for (int i = 0; i < (nPolygonPoints - 1); i++) {
            
            int x2 = xPolygon[i];
            int y2 = yPolygon[i];
            int x3 = xPolygon[i + 1];
            int y3 = yPolygon[i + 1];
            
            if (LinesAndAngles.pointIsInLine(xPt, yPt, x2, y2, x3, y3)) {
                return true;
            } 
            
            boolean intersects =
                LinesAndAngles.linesIntersect(xPt, yPt, xMax, yPt,
                x2, y2, x3, y3);
            
            if (intersects) {
                // avoid counting intersection with a vertex twice
                if (yPt != y3) {
                    sumIntersectingRays++;
                }
                
                if (yPt == y2) {
                    n2++;
                } 
                if (yPt == y3) {
                    n3++;
                }
                if ((yPt != y2) && (yPt != y3)) {
                    n0++;
                }
            }
        }
        
        boolean odd = ((sumIntersectingRays & 1) == 1);
        
        // for complex concave and convex curves:
        if (!odd && ((n0 & 1) == 0) && ((n2 & 1) == 1) && ((n3 & 1) == 1) ) {
            odd = true;
        }
        
        return odd;
    }
    
    /**
     * given a polygon (xPolygon, yPolygon) that has the same start and end
     * points, determine whether the point (xpt, ypt) is in the polygon.
     * 
     * @param xPt
     * @param yPt
     * @param xPolygon x values of a polygon in which each x should be smaller
     * than Integer.MAX_VALUE so that cross products do not overflow int type.
     * @param yPolygon y values of a polygon in which each y should be smaller
     * than Integer.MAX_VALUE so that cross products do not overflow int type.
     * @param nPolygonPoints
     * @return whether (xPt, yPt) is on or within {xPolygon, yPolygon}
     */
    public boolean isInSimpleCurve(float xPt, float yPt, float[] xPolygon,
        float[] yPolygon, int nPolygonPoints) {
        
        float xMax = MiscMath.findMax(xPolygon, nPolygonPoints) + 1;
        
        int sumIntersectingRays = 0;
        
        int n2 = 0;
        int n3 = 0;
        int n0 = 0;
                
        for (int i = 0; i < (nPolygonPoints - 1); i++) {
            
            float x2 = xPolygon[i];
            float y2 = yPolygon[i];
            float x3 = xPolygon[i + 1];
            float y3 = yPolygon[i + 1];
            
            if (LinesAndAngles.pointIsInLine(xPt, yPt, x2, y2, x3, y3)) {
                return true;
            }
            
            boolean intersects =
                LinesAndAngles.linesIntersect(xPt, yPt, xMax, yPt,
                x2, y2, x3, y3);
            
            if (intersects) {
                // avoid counting intersection with a vertex twice
                if (yPt != y3) {
                    sumIntersectingRays++;
                }
                
                if (yPt == y2) {
                    n2++;
                } 
                if (yPt == y3) {
                    n3++;
                }
                if ((yPt != y2) && (yPt != y3)) {
                    n0++;
                }
            }
        }
                
        boolean odd = ((sumIntersectingRays & 1) == 1);
        
        if (!odd && ((n0 & 1) == 0) && ((n2 & 1) == 1) && ((n3 & 1) == 1) ) {
            odd = true;
        }
        
        return odd;
    }
    
}
