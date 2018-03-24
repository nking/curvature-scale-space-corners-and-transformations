package algorithms.compGeometry;

import algorithms.imageProcessing.MiscellaneousCurveHelper;

/**
 *
 * @author nichole
 */
public class PointInPolygon2 {
    
    private final int[] xPolygon;
    private final int[] yPolygon;
    private final float[] xPolygonF;
    private final float[] yPolygonF;
    
    private final double[] xyCentroid;
    private final double minDistFromCenter;
    private final double maxDistFromCenter;
    // xmin, xmax, ymin, ymax
    private final double[] xyMinMax;
    
    /**
     * 
     * @param xPoly x coordinates of a polygon where first and last point are
     * the same
     * @param yPoly y coordinates of a polygon where first and last point are
     * the same
     */
    public PointInPolygon2(int[] xPoly, int[] yPoly) {
        
        this.xPolygon = xPoly;
        this.yPolygon = yPoly;
        
        this.xPolygonF = null;
        this.yPolygonF = null;
        
        MiscellaneousCurveHelper mch = new MiscellaneousCurveHelper();
        xyCentroid = mch.calculateXYCentroids(xPolygon, yPolygon);
        
        xyMinMax = new double[]{
            Double.MAX_VALUE, Double.MIN_VALUE, Double.MAX_VALUE, Double.MIN_VALUE
        };
        
        double diffX, diffY, diff;
        double min = Double.MAX_VALUE;
        double max = Double.MIN_VALUE;
        for (int i = 0; i < xPoly.length; ++i) {
            int xp = xPoly[i];
            int yp = yPoly[i];
            diffX = xp - xyCentroid[0];
            diffY = yp - xyCentroid[1];
            diff = Math.sqrt(diffX*diffX + diffY*diffY);
            if (diff < min) {
                min = diff;
            }
            if (diff > max) {
                max = diff;
            }
            if (xp < xyMinMax[0]) {
                xyMinMax[0] = xp;
            }
            if (yp < xyMinMax[2]) {
                xyMinMax[2] = yp;
            }
            if (xp > xyMinMax[1]) {
                xyMinMax[1] = xp;
            }
            if (yp > xyMinMax[3]) {
                xyMinMax[3] = yp;
            }
        }
        this.maxDistFromCenter = max;
        this.minDistFromCenter = min;
    }
    
    /**
     * 
     * @param xPoly x coordinates of a polygon where first and last point are
     * the same
     * @param yPoly y coordinates of a polygon where first and last point are
     * the same
     */
    public PointInPolygon2(float[] xPoly, float[] yPoly) {
        
        this.xPolygonF = xPoly;
        this.yPolygonF = yPoly;
        
        this.xPolygon = null;
        this.yPolygon = null;
        
        MiscellaneousCurveHelper mch = new MiscellaneousCurveHelper();
        xyCentroid = mch.calculateXYCentroids(xPolygon, yPolygon);
        
        xyMinMax = new double[]{
            Double.MAX_VALUE, Double.MIN_VALUE, Double.MAX_VALUE, Double.MIN_VALUE
        };
        
        double diffX, diffY, diff;
        double min = Double.MAX_VALUE;
        double max = Double.MIN_VALUE;
        for (int i = 0; i < xPoly.length; ++i) {
            float xp = xPoly[i];
            float yp = yPoly[i];
            diffX = xp - xyCentroid[0];
            diffY = yp - xyCentroid[1];
            diff = Math.sqrt(diffX*diffX + diffY*diffY);
            if (diff < min) {
                min = diff;
            }
            if (diff > max) {
                max = diff;
            }
            if (xp < xyMinMax[0]) {
                xyMinMax[0] = xp;
            }
            if (yp < xyMinMax[2]) {
                xyMinMax[2] = yp;
            }
            if (xp > xyMinMax[1]) {
                xyMinMax[1] = xp;
            }
            if (yp > xyMinMax[3]) {
                xyMinMax[3] = yp;
            }
        }
        this.maxDistFromCenter = max;
        this.minDistFromCenter = min;
    }
    
    public boolean isOutsideMaxRadius(float xp, float yp) {
        double diffX = xp - xyCentroid[0];
        double diffY = yp - xyCentroid[1];
        double diff = Math.sqrt(diffX*diffX + diffY*diffY);
    
        return (diff > this.maxDistFromCenter);
    }
    
    public boolean isOutsideBoundingBox(float xp, float yp) {
        if (xp < xyMinMax[0]) {
            return true;
        }
        if (xp > xyMinMax[1]) {
            return true;
        }
        if (yp < xyMinMax[2]) {
            return true;
        }
        if (yp > xyMinMax[3]) {
            return true;
        }
        return false;
    }
}
