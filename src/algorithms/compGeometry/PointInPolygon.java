package algorithms.compGeometry;

/**
 * test whether a point is within a polygon.
 *
 * follows http://rosettacode.org/wiki/Ray-casting_algorithm impl for Go that
 * counts the number of times a ray intersects the edges of the polygon starting
 * from the point and going ANY fixed direction. This algorithm is sometimes
 * also known as the crossing number algorithm or the even-odd rule algorithm
 *
 * @author nichole
 */
public class PointInPolygon {

    public PointInPolygon() {
    }

    float eps = 0.00001f;

    /**
     * check whether pt is in simple curve. note that the calculation returns
     * false if a pt is on the boundary of the polygon.
     *
     * adapted from the Go version:
     * http://rosettacode.org/wiki/Ray-casting_algorithm
     */
    public boolean isInSimpleCurve(float xPt, float yPt, float[] xPolygon,
        float[] yPolygon, int nPolygonPoints) {
        
        if (isInPolyPoints(xPt, yPt, xPolygon, yPolygon, 0, nPolygonPoints)) {
            return true;
        }

        int sumIntersectingRays = 0;
        for (int i = 0; i < nPolygonPoints; i++) {

            if ((i + 2) > nPolygonPoints) {
                continue;
            }

            boolean does = rayIntersects(xPt, yPt, xPolygon, yPolygon, i, i + 1);

            if (does) {
                sumIntersectingRays++;
            }
        }
        
        return ((sumIntersectingRays & 1) == 1);
    }

    /**
     * check whether pt is in simple curve. note that the calculation. returns
     * false if a pt is on the boundary of the polygon.
     *
     * Note that this method can be used when (xPolygon, yPolygon) do not have
     * the same starting and ending point.
     *
     * adapted from the Go version:
     * http://rosettacode.org/wiki/Ray-casting_algorithm
     */
    public boolean isInSimpleCurve_calcForPolygonWithoutDoubleStartEndPoint(
        float xPt, float yPt, float[] xPolygon, float[] yPolygon,
        int nPolygonPoints) {

        if (isInPolyPoints(xPt, yPt, xPolygon, yPolygon, 0, nPolygonPoints)) {
            return true;
        }
        
        int sumIntersectingRays = 0;
        for (int i = 0; i < nPolygonPoints; i++) {

            int i2 = i + 1;
            if (i2 > (nPolygonPoints - 1)) {
                i2 = 0;
            }

            boolean does = rayIntersects(xPt, yPt, xPolygon, yPolygon, i, i2);

            if (does) {
                sumIntersectingRays++;
            }
        }

        return ((sumIntersectingRays & 1) == 1);
    }

    /**
     * check whether pt is in simple curve. note that the calculation. returns
     * false if a pt is on the boundary of the polygon.
     *
     * adapted from the Go version:
     * http://rosettacode.org/wiki/Ray-casting_algorithm
     */
    public boolean isInSimpleCurve(float xPt, float yPt, float[] xPolygon,
        float[] yPolygon, int startPolygon, int endPolygon) {

        if (isInPolyPoints(xPt, yPt, xPolygon, yPolygon, startPolygon, endPolygon)) {
            return true;
        }
        
        int sumIntersectingRays = 0;
        for (int i = startPolygon; i < endPolygon; i++) {

            if ((i + 2) > endPolygon) {
                continue;
            }

            boolean does = rayIntersects(xPt, yPt, xPolygon, yPolygon, i, i + 1);

            if (does) {
                sumIntersectingRays++;
            }
        }

        return ((sumIntersectingRays % 2) == 1);
    }

    /**
     * pt intersects edge by projecting that the point will travel horizontally.
     *
     * @return
     */
    boolean rayIntersects(float xPt, float yPt, float[] xPolygon,
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
        if ((yPt == ay) || (yPt == by)) {
            yPt += eps;
        }

        if ((yPt < ay) || (yPt > by)) {
            return false;
        }
        
        if (ax > bx) {
            if (xPt > ax) {
                return false;
            }
            if (xPt < bx) {
                return true;
            }
        } else {
            if (xPt > bx) {
                return false;
            }
            if (xPt < ax) {
                return true;
            }
            //TODO: Revisit this.  it will miss for right side
            if ((xPt == ax) || (xPt == bx)) {
                xPt += eps;
            }
        }

        boolean intersects = (yPt - ay)/(xPt - ax) >= (by - ay)/(bx - ax);
        
        return intersects;
    }

    protected boolean isInPolyPoints(float xPt, float yPt, float[] xPolygon, 
        float[] yPolygon, int startPolygonIdx, int stopPolygonIdxExcl) {
        
        for (int i = startPolygonIdx; i < stopPolygonIdxExcl; i++) {
            
            if (Math.abs(xPt - xPolygon[i]) < eps) {
                if (Math.abs(yPt - yPolygon[i]) < eps) {
                    return true;
                }
            }
        }
        
        return false;
    }

}
