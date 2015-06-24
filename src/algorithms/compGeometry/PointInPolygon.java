package algorithms.compGeometry;

/**
 * test whether a point is within a polygon.
 *
 * follows http://rosettacode.org/wiki/Ray-casting_algorithm impl for Go that
 * counts the number of times a ray intersects the edges of the polygon starting
 * from the point and going ANY fixed direction. This algorithm is sometimes
 * also known as the crossing number algorithm or the even-odd rule algorithm.
 * The implementation here is different, and includes logic for testing whether
 * the point is in the polygon lines too.
 *
 * @author nichole
 */
public class PointInPolygon {

    public PointInPolygon() {
    }

    float eps = 0.00001f;

    public boolean isInSimpleCurve(float xPt, float yPt, float[] xPolygon,
        float[] yPolygon, int nPolygonPoints) {
        
        int sumIntersectingRays = 0;
        
        boolean possibleDoubleCount = false;
        
        for (int i = 0; i < nPolygonPoints; i++) {

            if ((i + 2) > nPolygonPoints) {
                continue;
            }

            /*
            returns 0 if does not intersect, returns 1 if it does intersect,
             returns 2 if the point lies on a line directly (and in that
             case, the invoker should not test further),
             returns 3 for yPt being equal to one of the segments and xPt being
             to the right of both segments (and in that case, the invoker should
             not count the result twice).
            */
            int result = rayIntersects(xPt, yPt, xPolygon, yPolygon, i, i + 1);
        
            if (result == 0) {
                continue;
            } else if (result == 2) {
                return true;
            }
            
            if (result == 3) {
                if (possibleDoubleCount) {
                    // do not count this, but do reset flag
                    possibleDoubleCount = false;
                    continue;
                }
                possibleDoubleCount = true;                
            }
          
            sumIntersectingRays++;
        }
        
        return ((sumIntersectingRays & 1) == 1);
    }

    /**
     * pt intersects edge by projecting that the point will travel horizontally.
     * returns 0 if does not intersect, returns 1 if it does intersect,
     * returns 2 if the point lies on a line directly (and in that
     * case, the invoker should not test further),
     * returns 3 for yPt being equal to one of the segments and xPt being
     * to the right of both segments (and in that case, the invoker should
     * not count the result twice).
     *
     * @return
     */
    int rayIntersects(float xPt, float yPt, float[] xPolygon,
        float[] yPolygon, int index1, int index2) {

        float ax = xPolygon[index1];
        float ay = yPolygon[index1];
        float bx = xPolygon[index2];
        float by = yPolygon[index2];
        
        if (ay >= by) {
            float xtmp = ax;
            float ytmp = ay;
            ax = bx;
            ay = by;
            bx = xtmp;
            by = ytmp;
        }
        
        if ((xPt == ax) && (yPt == ay)) {
            return 2;
        }
        if ((xPt == bx) && (yPt == by)) {
            return 2;
        }
        
        if ((yPt < ay) || (yPt > by)) {
            return 0;
        }
        
        // test for vertical a->b,  "is on line" and "right of line"
        if (ax == bx) {
            if (xPt == ax) {
                if ((yPt > ay) && (yPt < by)) {
                    return 2;
                }
                return 0;
            } else if (xPt < ax) {
                if ((yPt > ay) && (yPt < by)) {
                    return 1;
                }
                return 0;
            } else {
                return 0;
            }
        }
        
        // since (xPt,yPt) is within y range ay:by
        if ((xPt <= ax) && (xPt <= bx)) {
            if ((yPt == by) || (yPt == ay)) {
                // we don't want to count the point twice as "intersects"
                return 3;
            }
            return 1;
        }
        if ((xPt > ax) && (xPt > bx)) {
            return 0;
        }
        
        // (xPt,yPt) is within y range ay:by  and bounded by ax and bx (but != ax nor bx)
        
        // test "is on line"
        float slopeAB = (by - ay)/(bx - ax);
        
        float slopeAPt = (yPt - ay)/(xPt - ax);
        
        float slopePtB = (by - yPt)/(bx - xPt);
        
        // TODO: may need to widen eps to include 1 pixel rounding
        if ((Math.abs(slopeAB - slopeAPt) < eps) && 
            (Math.abs(slopeAB - slopePtB) < eps)) {
            
            return 2;
        }
        
        if ((xPt > ax) && (bx > xPt)) {
            float xLine = ax + ((yPt - ay)/slopeAB);
            if (xPt < xLine) {
                return 1;
            }
        }
        
        return 0;
    }

}
