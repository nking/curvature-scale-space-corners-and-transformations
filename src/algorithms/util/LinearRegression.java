package algorithms.util;

import algorithms.KSelect;
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
                float diffX = x[j] - x[i];
                if (diffX == 0) {
                    continue;
                }
                float diffY = y[j] - y[i];
                s[count] = diffY/diffX;
                count++;
            }            
        }
        
        KSelect kSelect = new KSelect();
        
        float median = kSelect.findMedianOfMedians(s, 0, count - 1);
        
        //y_i − median*x_i
        log.info("thiel sen beta=" + median);
       
        // find the y-intercept as the median of the values y[i] − median * x[i]
        float[] s2 = new float[x.length];
        for (int i = 0; i < x.length; i++) {
            s2[i] = y[i] - median * x[i];
        }
                    
        int s2MedianIdx = kSelect.findMedianOfMediansIdx(s2, 0, s2.length - 1);
        /*
           (y1 - y0)/(x1 - x0) = slope
            y1 - y0 = slope*(x1 - x0);
            y1 = y0 + slope*(x1 - x0);
            y1 = (y0 - slope*x0) + slope*x1
            y1 =  yIntercept     + slope*x1
        */
        
        float yIntercept = y[s2MedianIdx] - median * x[s2MedianIdx];
        
        //the estimation of yIntercept needed to be improved
        int np = 10;
        while (((s2MedianIdx - np) < 0) || ((s2MedianIdx + np) > (x.length - 1))) {
            np--;
            if (np < 0 || np == 0) {
                break;
            }
        }
        if (np > 0) {
            float sum = 0;
            for (int j = (s2MedianIdx - np); j <= (s2MedianIdx + np); j++) {
                sum += (y[j] - median * x[j]);
            }
            yIntercept = sum/((float)(2*np + 1));
        }
        
        return new float[]{yIntercept, median};
    }
    
    public float[] calculateParamsForLinearRegression(PairIntArray xy1, PairIntArray xy2) {
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
    
    public void plotTheLinearRegression(int[] dx, int[] dy) {
                        
        float[] tsbParams = calculateTheilSenEstimatorParams(dx, dy);
        
        float yIntercept = tsbParams[0];
        
        float slope = tsbParams[1];
        
        /*
        plot dx, dy
        and plot a line generated from the yIntercept and median: yIntercept − median*x_i
        */        
        int xMin = MiscMath.findMin(dx);
        int xMax = MiscMath.findMax(dx);
        int len = xMax - xMin + 1;
        int[] tsbX = new int[len];
        int[] tsbY = new int[len];
        int count = 0;
        for (int x = xMin; x <= xMax; x++) {
            float y = yIntercept + slope * x;
            tsbX[count] = x;
            tsbY[count] = (int)Math.round(y);
            count++;
        }
        
        int yMin = MiscMath.findMin(dy);
        int yMax = MiscMath.findMax(dy);
       
        try {
            PolygonAndPointPlotter plotter = new PolygonAndPointPlotter();
            plotter.addPlot(
                xMin, xMax, yMin, yMax,
                dx, dy, 
                tsbX, tsbY,
                "diff X vs diff Y and thiel sen beta linear regression line");

            plotter.writeFile();
            
        } catch(IOException e) {
            
            log.severe("ERROR while trying to write plot: " + e.getMessage());
        }
    }
}
