package algorithms.compGeometry;

import algorithms.compGeometry.convexHull.GrahamScan;
import algorithms.compGeometry.convexHull.GrahamScanTooFewPointsException;
import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageDisplayer;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.MiscellaneousCurveHelper;
import algorithms.misc.MiscMath;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

/**
 *
 * @author nichole
 */
public class PerimeterFinder {
    
    /**
     * for the bounds of rows present in points, find the min and max of columns
     * and return the result as a map with key = row number, value = pair of
     * ints with x being the start column for that row and y being the stop
     * column (inclusive) for that row.  note that rows without points were
     * interpreted from their surrounding rows.
     * 
     * @param points
     * @param outputRowMinMax output populated as the min and max of rows are 
     * determined.
     * @return 
     */
    public Map<Integer, PairInt> find(PairIntArray points, int[] outputRowMinMax) {
        
        int minY = Integer.MAX_VALUE;
        int maxY = Integer.MIN_VALUE;
        
        float[] xP = new float[points.getN()];
        float[] yP = new float[points.getN()];
        
        // O(N)
        for (int i = 0; i < points.getN(); i++) {
            
            int x = points.getX(i);
            int y = points.getY(i);
            if (y < minY) {
                minY = y;
            }
            if (y > maxY) {
                maxY = y;
            }
            
            xP[i] = x;
            yP[i] = y;
        }
        
        outputRowMinMax[0] = minY;
        outputRowMinMax[1] = maxY;
        
        return find(outputRowMinMax, xP, yP);
    }
    
    /**
     * for the bounds of rows present in points, find the min and max of columns
     * and return the result as a map with key = row number, value = pair of
     * ints with x being the start column for that row and y being the stop
     * column (inclusive) for that row.  note that rows without points were
     * interpreted from their surrounding rows.
     * 
     * @param points
     * @param outputRowMinMax output populated as the min and max of rows are 
     * determined.
     * @return 
     */
    public Map<Integer, PairInt> find(Set<PairInt> points, int[] outputRowMinMax) {
        
        int minY = Integer.MAX_VALUE;
        int maxY = Integer.MIN_VALUE;
        
        float[] xP = new float[points.size()];
        float[] yP = new float[points.size()];
            
        int i = 0;
        // O(N)
        for (PairInt p : points) {
            
            int x = p.getX();
            int y = p.getY();
            if (y < minY) {
                minY = y;
            }
            if (y > maxY) {
                maxY = y;
            }
            
            xP[i] = x;
            yP[i] = y;
            i++;
        }
        
        outputRowMinMax[0] = minY;
        outputRowMinMax[1] = maxY;
        
        return find(outputRowMinMax, xP, yP);
    }
    
    protected Map<Integer, PairInt> find(int[] minMax, float[] xP, float[] yP) {
        
        if (xP == null) {
	    	throw new IllegalArgumentException("xP cannot be null");
        }
	    if (yP == null) {
	    	throw new IllegalArgumentException("yP cannot be null");
        }
        if (minMax == null) {
	    	throw new IllegalArgumentException("minMax cannot be null");
        }
                
        Map<Integer, PairInt> rowColRange = new HashMap<Integer, PairInt>();
        
        int minY = minMax[0];
        int maxY = minMax[1];
        
        try {
            
            GrahamScan scan = new GrahamScan();
            scan.computeHull(xP, yP);
        
            float[] xHull = scan.getXHull();
            float[] yHull = scan.getYHull();
            
            int len = (maxY - minY) + 1;
            
            // the indexes correspond to idx + yMin so that the first is for yMin
            float[] startCols = new float[len];
            float[] stopCols = new float[len];
            
            populateStartAndStopColsForAllRows(xHull, yHull, minY, maxY, 
                startCols, stopCols);
            
            for (int i = 0; i < startCols.length; i++) {
                int startCol = Math.round(startCols[i]);
                int stopCol = Math.round(stopCols[i]);
                PairInt p = new PairInt(startCol, stopCol);
                
                Integer row = Integer.valueOf(i + minY);
                
                rowColRange.put(row, p);
            }
            
        } catch(GrahamScanTooFewPointsException e) {
            
            int minX = Math.round(MiscMath.findMin(xP));
            int maxX = Math.round(MiscMath.findMax(xP));
            
            for (int row = minY; row <= maxY; row++) {
                Integer r = Integer.valueOf(row);
                PairInt p = new PairInt(minX, maxX);
                rowColRange.put(r, p);
            }
            
            return rowColRange;
        }
        
        return rowColRange;
    }

    /**
     * given the convex hull, populate start and stop columns for each row
     * between the minimum and maximum rows.  Note that the method expects
     * a hull ordered clockwise and having the same first and last point
     * (these characteristics are present in the hull created by GrahamScan).
     * @param xHull
     * @param yHull
     * @param minY
     * @param maxY
     * @param startCols
     * @param stopCols 
     */
    private void populateStartAndStopColsForAllRows(float[] xHull, float[] yHull, 
        int minY, int maxY, float[] startCols, float[] stopCols) {
        
        if (xHull == null) {
            throw new IllegalArgumentException("xHull cannot be null");
        }
        if (yHull == null) {
            throw new IllegalArgumentException("yHull cannot be null");
        }
        if (startCols == null) {
            throw new IllegalArgumentException("startCols cannot be null");
        }
        if (stopCols == null) {
            throw new IllegalArgumentException("stopCols cannot be null");
        }
        if (xHull.length != yHull.length) {
            throw new IllegalArgumentException(
            "xHull and yHull must have same length");
        }
        if (stopCols.length != startCols.length) {
            throw new IllegalArgumentException(
            "stopCols and startCols must have same length");
        }
        
        int len = (maxY - minY) + 1;
        
        if (startCols.length != len) {
            throw new IllegalArgumentException(
            "startCols length is expected to be " + len + " from minY and maxY");
        }
        
        Arrays.fill(startCols, -1);
        Arrays.fill(stopCols, -1);
       
        /*
         walk the hulls and fill in each row's start or stopCol between
         it and next point.

         NOTE: the hull created by GrahamScan is ordered clockwise and includes
         the same first and last point.
         */

        boolean onRightHull = true;

        for (int hIdx = 0; hIdx < (xHull.length - 1); hIdx++) {

            int compIdx = (hIdx + 1);
            
            float xh = xHull[hIdx];
            float yh = yHull[hIdx];

            int yi = Math.round(yh);
            int yj = Math.round(yHull[compIdx]);
            
            if (yi == maxY) {
                onRightHull = false;
            }
            
            int nRows = Math.abs(yj - yi) + 1;
            
            float dx = (nRows == 1) ? 0 :
                (xHull[compIdx] - xh) / (float) (nRows - 1);

            float xStart = xh;
            
            if (yj < yi) {
                int swap = yi;
                yi = yj;
                yj = swap;
                dx *= -1;
                xStart = xHull[compIdx];
            }
            
            for (int row = yi; row <= yj; row++) {
                int nr = row - yi;
                float x = xStart + (nr * dx);
                int idx = row - minY;
                if (onRightHull) {
                    stopCols[idx] = x;
                } else {
                    startCols[idx] = x;
                }
            }
        }
        
        // assert that all are filled
        for (int i = 0; i < startCols.length; i++) {
            if (startCols[i] == -1 || stopCols[i] == -1) {
                throw new IllegalStateException(
                "walk the hull to fill start and stop is not yet correct");
            }
        }
    }
}
