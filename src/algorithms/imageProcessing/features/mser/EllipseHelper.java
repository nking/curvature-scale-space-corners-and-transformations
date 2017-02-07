package algorithms.imageProcessing.features.mser;

import algorithms.QuickSort;
import algorithms.misc.MiscMath;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import java.util.HashSet;
import java.util.Set;

/**
 *
 * @author nichole
 */
public class EllipseHelper {
    
    // note that this may include values that extend beyond the original image
    // dimensions
    private final PairIntArray ellipse;
    
    private final int x;
    
    private final int y;
    
    private final Bounds bounds;
    
    private final double minor;
    private final double major;
    private final double eccentricity;
    
    /**
     * constructor w/ center x, y and coefficients 
     * v0x, v1x, v0y, v1y, e0sq, e1sq that are
     * the eigenvectors and values of the ellipse 
     * calculated from the spatial moments of the ellipse.
     * @param x
     * @param y
     * @param coeffs 
     */
    public EllipseHelper(int x, int y, double[] coeffs) {
        
        double v0x = coeffs[0];
        double v1x = coeffs[1];
        double v0y = coeffs[2];
        double v1y = coeffs[3];
        
        ellipse = new PairIntArray();
        
        this.major = 2. * coeffs[4];
        this.minor = 2. * coeffs[5];
        
        this.x = x;
        this.y = y;
        
        for (double t = 0.0; t < 2.0 * Math.PI; t += 0.001) {
            
            double mc = Math.cos(t);
            double ms = Math.sin(t);
            
            int x2 = (int)Math.round(x + (mc * v0x + ms * v1x) * 2.0 + 0.5);
            int y2 = (int)Math.round(y + (mc * v0y + ms * v1y) * 2.0 + 0.5);

            ellipse.add(x2, y2);
        }
        
        eccentricity = Math.sqrt(major*major - minor*minor)/major;
        
        bounds = new Bounds(ellipse);
    }
    
    /**
     * constructor w/ center x, y and points of the ellipse
     * @param xy points of the ellipse
     */
    public EllipseHelper(int x, int y, PairIntArray xy) {
        
        ellipse = xy;
        this.x = x;
        this.y = y;
        bounds = new Bounds(ellipse);
        
        // use a Region to calculate the major and minor axes
        Region r = new Region();
        for (int i = 0; i < xy.getN(); ++i) {
            int x2 = xy.getX(i);
            int y2 = xy.getY(i);
            r.accumulate(x2, y2);
        }
        double[] coeffs = r.calcParamTransCoeff();
        this.major = 2. * coeffs[4];
        this.minor = 2. * coeffs[5];
    
        eccentricity = Math.sqrt(major*major - minor*minor)/major;
    }
    
    public double getEccentricity() {
        return eccentricity;
    }
    
    public PairIntArray getEllipse() {
        return ellipse;
    }
    
    /**
     * test whether the point (xPoint, yPoint) is within or on the ellipse
     * boundary.
     * 
     * @param xPoint
     * @param yPoint
     * @return 
     */
    public boolean isWithin(int xPoint, int yPoint) {
                
        return bounds.isWithin(xPoint, yPoint);
    }
    
    public double getMajorTimesMinor() {
        return major * minor;
    }
    
    public boolean intersects(EllipseHelper other) {
        
        if (other.bounds.maxRow < bounds.minRow || 
            other.bounds.minRow > bounds.maxRow) {
            return false;
        }
        if (other.bounds.maxCol < bounds.minCol || 
            other.bounds.minCol > bounds.maxCol) {
            return false;
        }
        PairIntArray otherXY = other.getEllipse();
        for (int i = 0; i < otherXY.getN(); ++i) {
            if (bounds.isWithin(otherXY.getX(i), otherXY.getY(i))) {
                return true;
            }
        }
        
        // check for whether one is embedded within the other
        // and the coordinates of the ellipses do not intersect
        // so wasn't detected yet
        if (surrounds(other)) {
            return true;
        }
        if (other.surrounds(this)) {
            return true;
        }
        
        return false;
    }
    
    public boolean surrounds(EllipseHelper other) {
        if (bounds.minRow < other.bounds.minRow &&
            bounds.maxRow > other.bounds.maxRow &&
            bounds.minCol < other.bounds.minCol &&
            bounds.maxCol > other.bounds.maxCol) {
            return true;
        }
        return false;
    }
    
    private class Bounds {
        
        final int minRow;
        final int maxRow;
        final TIntObjectMap<PairInt> rowBounds;
    
        final int minCol;
        final int maxCol;
        
        public Bounds(PairIntArray xy) {
            
            int[] minMaxXY = MiscMath.findMinMaxXY(xy);
            
            this.minRow = minMaxXY[2];
            this.maxRow = minMaxXY[3];
            this.rowBounds = new TIntObjectHashMap<PairInt>();
            
            int minX = Integer.MAX_VALUE;
            int maxX = Integer.MIN_VALUE;
            
            for (int i = 0; i < xy.getN(); ++i) {
                
                int col = xy.getX(i);
                
                int row = xy.getY(i);
                
                PairInt xMinMax = rowBounds.get(row);
                if (xMinMax == null) {
                    xMinMax = new PairInt(col, col);
                    rowBounds.put(row, xMinMax);
                } else {
                    if (col < xMinMax.getX()) {
                        xMinMax.setX(col);
                    } else if (col > xMinMax.getY()) {
                        xMinMax.setY(col);
                    }
                }
                
                if (col < minX) {
                    minX = col;
                }
                if (col > maxX) {
                    maxX = col;
                }
            }
            this.minCol = minX;
            this.maxCol = maxX;
        }
        
        boolean isWithin(int xPoint, int yPoint) {
            
            if (yPoint < minRow || yPoint > maxRow) {
                return false;
            }
            if (xPoint < minCol || xPoint > maxCol) {
                return false;
            }
            
            PairInt xMinMax = rowBounds.get(yPoint);
            assert(xMinMax != null);
            
            if (xPoint < xMinMax.getX() || xPoint > xMinMax.getY()) {
                return false;
            }
            
            return true;
        }
    
        Set<PairInt> createEllipseFillingPoints() {
            
            Set<PairInt> points = new HashSet<PairInt>();
        
            for (int r = minRow; r <= maxRow; ++r) {
                PairInt xMinMax = rowBounds.get(r);
                for (int c = xMinMax.getX(); c <= xMinMax.getY(); ++c) {
                    points.add(new PairInt(c, r));
                }
            }
        
            return points;
        }
    }
    
    public Set<PairInt> createEllipseFillingPoints() {
        return bounds.createEllipseFillingPoints();
    }
    
    public int getSemiMinor() {
        return (int)Math.round(minor);
    }
    
    public int getSemiMajor() {
        return (int)Math.round(major);
    }
    
    public int[] getXYCenter() {
        return new int[]{x, y};
    }
}
