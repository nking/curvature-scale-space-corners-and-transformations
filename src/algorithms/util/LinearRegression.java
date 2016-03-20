package algorithms.util;

import algorithms.QuickSort;
import algorithms.misc.MiscMath;
import java.io.IOException;
import java.util.Arrays;
import java.util.logging.Logger;

/**
 *
 * @author nichole
 */
public class LinearRegression {
    
    private Logger log = Logger.getLogger(this.getClass().getName());
    
    void calculateXYDifferences(PairIntArray xy1, PairIntArray xy2,
        int[] dxOutput, int[] dyOutput) {
        
        int n = xy1.getN();
        
        for (int i = 0; i < n; i++) {
            int diffX = xy1.getX(i) - xy2.getX(i);
            int diffY = xy1.getY(i) - xy2.getY(i);            
            dxOutput[i] = diffX;
            dyOutput[i] = diffY;         
        }
    }
    
    /**
     * calculate the theil sen estimator for the set of points and return
     * the yIntercept and slope that can be used to plot a line that is the
     * linear regression of the x and y points.
     * NOTE: a side effect of the method is that x and y become partially
     * sorted.
     * @param x
     * @param y
     * @return 
     */
    public float[] calculateTheilSenEstimatorParams(int[] x, int[] y) {
        
        int n = x.length;
        
        /*      
        for 1000 points, for each possible pair w/ image 2 points,
        the real solution would be looking for a match within 
        2.5*stdev or 3 * stdev      
        */
        
        /* linear regression w/ theil sen estimator:
        http://en.wikipedia.org/wiki/Theil%E2%80%93Sen_estimator
        
        median m of the slopes (yj − yi)/(xj − xi) determined by all pairs of 
        sample points. 
        */
        int count = 0;
        float[] s = new float[n*n];
        for (int i = 0; i < n; i++) {
            for (int j = (i + 1); j < n; j++) {
                if ((i == j) || (x[j] - x[i]) == 0) {
                    continue;
                }
                s[count] = (float)(y[j] - y[i])/((float)x[j] - x[i]);
                count++;
            }
        }
        
        s = Arrays.copyOf(s, count);
        Arrays.sort(s);
        int idx = s.length/2;
        float median;
        if ((idx & 1) == 0) {
            median = (s[idx] + s[idx - 1])/2.f;
        } else {
            median = s[idx];
        }
        
        log.info("thiel sen beta=" + median);
       
        // find the y-intercept as the median of the values y[i] − median * x[i]
        float[] s2 = new float[x.length];
        for (int i = 0; i < x.length; i++) {
            s2[i] = y[i] - median * x[i];
        }
        QuickSort.sort(s2, x, y, 0, s2.length - 1);
        int medianIdx = s2.length/2;
        
        /*
           (y1 - y0)/(x1 - x0) = slope
            y1 - y0 = slope*(x1 - x0);
            y1 = y0 + slope*(x1 - x0);
            y1 = (y0 - slope*x0) + slope*x1
            y1 =  yIntercept     + slope*x1
        */
        
        float yIntercept = y[medianIdx] - median * x[medianIdx];
        
        //the estimation of yIntercept needs to be improved:
        int np = 10;
        while (((medianIdx - np) < 0) || ((medianIdx + np) > (x.length - 1))) {
            np--;
            if (np < 0 || np == 0) {
                break;
            }
        }
        if (np > 0) {
            float sum = 0;
            for (int j = (medianIdx - np); j <= (medianIdx + np); j++) {
                sum += (y[j] - median * x[j]);
            }
            yIntercept = sum/((float)(2*np + 1));
        }
        
        return new float[]{yIntercept, median};
    }
    
    public float[] calculateParamsForLinearRegression(PairIntArray xy1, 
        PairIntArray xy2) {
        
        int n = xy1.getN();
        
        int[] dx = new int[n];
        int[] dy = new int[n];
        
        calculateXYDifferences(xy1, xy2, dx, dy);
        
        return calculateTheilSenEstimatorParams(dx, dy);
    }
    
    public void plotTheLinearRegression(PairIntArray xy1, PairIntArray xy2) {
        
        int n = xy1.getN();
        
        int[] dx = new int[n];
        int[] dy = new int[n];
        
        calculateXYDifferences(xy1, xy2, dx, dy);
        
        plotTheLinearRegression(dx, dy);
    }
    
    /**
     * make a plot of the linear regression of arrays x and y.
     * NOTE: a side effect of the method is that x and y become partially
     * sorted.
     * @param x
     * @param y 
     */
    public void plotTheLinearRegression(int[] x, int[] y) {
                        
        float[] tsbParams = calculateTheilSenEstimatorParams(x, y);
        
        float yIntercept = tsbParams[0];
        
        float slope = tsbParams[1];
        
        /*
        plot dx, dy
        and plot a line generated from the yIntercept and median: yIntercept − median*x_i
        */        
        int xMin = MiscMath.findMin(x);
        int xMax = MiscMath.findMax(x);
        int len = xMax - xMin + 1;
        int[] tsbX = new int[len];
        int[] tsbY = new int[len];
        int count = 0;
        for (int xCoord = xMin; xCoord <= xMax; xCoord++) {
            float yCoord = yIntercept + slope * (float)xCoord;
            tsbX[count] = xCoord;
            tsbY[count] = Math.round(yCoord);
            count++;
        }
        
        int yMin = MiscMath.findMin(y);
        int yMax = MiscMath.findMax(y);
       
        try {
            PolygonAndPointPlotter plotter = new PolygonAndPointPlotter();
            plotter.addPlot(
                xMin, xMax, yMin, yMax,
                x, y, 
                tsbX, tsbY,
                "X vs Y and thiel sen beta linear regression line");

            plotter.writeFile();
            
        } catch(IOException e) {
            
            log.severe("ERROR while trying to write plot: " + e.getMessage());
        }
    }
    
    // ======================================================================
    void calculateXYDifferences(PairFloatArray xy1, PairFloatArray xy2,
        float[] dxOutput, float[] dyOutput) {
        
        int n = xy1.getN();
        
        for (int i = 0; i < n; i++) {
            float diffX = xy1.getX(i) - xy2.getX(i);
            float diffY = xy1.getY(i) - xy2.getY(i);            
            dxOutput[i] = diffX;
            dyOutput[i] = diffY;         
        }
    }
    
    /**
     * calculate the theil sen estimator for the set of points and return
     * the yIntercept and slope that can be used to plot a line that is the
     * linear regression of the x and y points.
     * NOTE: a side effect of the method is that x and y become partially
     * sorted.
     * @param x
     * @param y
     * @return 
     */
    public float[] calculateTheilSenEstimatorParams(float[] x, float[] y) {
        
        int n = x.length;
        
        /*      
        for 1000 points, for each possible pair w/ image 2 points,
        the real solution would be looking for a match within 
        2.5*stdev or 3 * stdev      
        */
        
        /* linear regression w/ theil sen estimator:
        http://en.wikipedia.org/wiki/Theil%E2%80%93Sen_estimator
        
        median m of the slopes (yj − yi)/(xj − xi) determined by all pairs of 
        sample points. 
        */
        int count = 0;
        float[] s = new float[n*n];
        for (int i = 0; i < n; i++) {
            for (int j = (i + 1); j < n; j++) {
                if ((i == j) || (x[j] - x[i]) == 0) {
                    continue;
                }
                s[count] = (y[j] - y[i])/(x[j] - x[i]);
                count++;
            }
        }
        
        if (count == 0) {
            // this can happen for vertical lines
            return new float[]{Float.NaN, Float.MAX_VALUE};
        }
        
        float median;
        s = Arrays.copyOf(s, count);
        Arrays.sort(s);
        int idx = s.length/2;
        if ((idx & 1) == 0) {
            median = (s[idx] + s[idx - 1])/2.f;
        } else {
            median = s[idx];
        }
        
        log.info("thiel sen beta=" + median);
       
        // find the y-intercept as the median of the values y[i] − median * x[i]
        float[] s2 = new float[x.length];
        for (int i = 0; i < x.length; i++) {
            s2[i] = y[i] - median * x[i];
        }
        QuickSort.sort(s2, x, y, 0, s2.length - 1);
        int medianIdx = s2.length/2;
        
        /*
           (y1 - y0)/(x1 - x0) = slope
            y1 - y0 = slope*(x1 - x0);
            y1 = y0 + slope*(x1 - x0);
            y1 = (y0 - slope*x0) + slope*x1
            y1 =  yIntercept     + slope*x1
        */
        
        float yIntercept = y[medianIdx] - median * x[medianIdx];
        
        //the estimation of yIntercept needs to be improved:
        int np = 10;
        while (((medianIdx - np) < 0) || ((medianIdx + np) > (x.length - 1))) {
            np--;
            if (np < 0 || np == 0) {
                break;
            }
        }
        if (np > 0) {
            float sum = 0;
            for (int j = (medianIdx - np); j <= (medianIdx + np); j++) {
                sum += (y[j] - median * x[j]);
            }
            yIntercept = sum/((float)(2*np + 1));
        }
        
        return new float[]{yIntercept, median};
    }
    
    /**
     * calculate the theil sen estimator for the set of points and return
     * the yIntercept and slope that can be used to plot a line that is the
     * linear regression of the x and y points.
     * NOTE: a side effect of the method is that x and y become partially
     * sorted.
     * @param x
     * @param y
     * @return 
     */
    public float[] calculateTheilSenEstimatorMedian(float[] x, float[] y) {
        
        int n = x.length;
        
        /*      
        for 1000 points, for each possible pair w/ image 2 points,
        the real solution would be looking for a match within 
        2.5*stdev or 3 * stdev      
        */
        
        /* linear regression w/ theil sen estimator:
        http://en.wikipedia.org/wiki/Theil%E2%80%93Sen_estimator
        
        median m of the slopes (yj − yi)/(xj − xi) determined by all pairs of 
        sample points. 
        */
        int count = 0;
        float[] s = new float[n*n];
        for (int i = 0; i < n; i++) {
            for (int j = (i + 1); j < n; j++) {
                if ((i == j) || (x[j] - x[i]) == 0) {
                    continue;
                }
                s[count] = (y[j] - y[i])/(x[j] - x[i]);
                count++;
            }
        }
        
        s = Arrays.copyOf(s, count);
        Arrays.sort(s);
        int idx = s.length/2;
        float median;
        if ((idx & 1) == 0) {
            median = (s[idx] + s[idx - 1])/2.f;
        } else {
            median = s[idx];
        }
        
        log.info("thiel sen beta=" + median);
       
        // find the y-intercept as the median of the values y[i] − median * x[i]
        float[] s2 = new float[x.length];
        for (int i = 0; i < x.length; i++) {
            s2[i] = y[i] - median * x[i];
        }
        QuickSort.sort(s2, x, y, 0, s2.length - 1);
        int medianIdx = s2.length/2;
       
        float xMedian = x[medianIdx];
        float yMedian = y[medianIdx];
        // improve the vlue over several points
        int np = 10;
        while (((medianIdx - np) < 0) || ((medianIdx + np) > (x.length - 1))) {
            np--;
            if (np < 0 || np == 0) {
                break;
            }
        }
        if (np > 0) {
            float sumX = 0;
            float sumY = 0;
            for (int j = (medianIdx - np); j <= (medianIdx + np); j++) {
                sumX += x[j];
                sumY += y[j];
            }
            xMedian = sumX/((float)(2*np + 1));
            yMedian = sumY/((float)(2*np + 1));
        }
        
        return new float[]{xMedian, yMedian};
    }
    
    public float[] calculateParamsForLinearRegression(PairFloatArray xy1, 
        PairFloatArray xy2) {
        
        int n = xy1.getN();
        
        float[] dx = new float[n];
        float[] dy = new float[n];
        
        calculateXYDifferences(xy1, xy2, dx, dy);
        
        return calculateTheilSenEstimatorParams(dx, dy);
    }
    
    public void plotTheLinearRegression(PairFloatArray xy1, PairFloatArray xy2) {
        
        int n = xy1.getN();
        
        float[] dx = new float[n];
        float[] dy = new float[n];
        
        calculateXYDifferences(xy1, xy2, dx, dy);
        
        plotTheLinearRegression(dx, dy);
    }
    
    /**
     * make a plot of the linear regression of arrays x and y.
     * NOTE: a side effect of the method is that x and y become partially
     * sorted.
     * @param x
     * @param y 
     */
    public void plotTheLinearRegression(float[] x, float[] y) {
                        
        float[] tsbParams = calculateTheilSenEstimatorParams(x, y);
        
        float yIntercept = tsbParams[0];
        
        float slope = tsbParams[1];
        
        /*
        plot dx, dy
        and plot a line generated from the yIntercept and median: yIntercept − median*x_i
        */        
        int xMin = (int)Math.floor(MiscMath.findMin(x)) - 1;
        int xMax = (int)Math.ceil(MiscMath.findMax(x)) + 1;
        int len = xMax - xMin + 1;
        float[] tsbX = new float[len];
        float[] tsbY = new float[len];
        int count = 0;
        for (int xCoord = xMin; xCoord <= xMax; xCoord++) {
            float yCoord = yIntercept + slope * (float)xCoord;
            tsbX[count] = xCoord;
            tsbY[count] = yCoord;
            count++;
        }
        
        int yMin = (int)Math.floor(MiscMath.findMin(y)) - 1;
        int yMax = (int)Math.ceil(MiscMath.findMax(y)) + 1;
       
        try {
            PolygonAndPointPlotter plotter = new PolygonAndPointPlotter();
            plotter.addPlot(
                xMin, xMax, yMin, yMax,
                x, y, 
                tsbX, tsbY,
                "X vs Y and thiel sen beta linear regression line");

            plotter.writeFile();
            
        } catch(IOException e) {
            
            log.severe("ERROR while trying to write plot: " + e.getMessage());
        }
    }
}
