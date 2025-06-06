package algorithms.compGeometry;

import algorithms.sort.MultiArrayMergeSort;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.util.AngleUtil;
import algorithms.misc.Misc;
import algorithms.statistics.LinearRegression;
import algorithms.util.PairFloatArray;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.PairIntArrayWithColor;
import gnu.trove.list.TIntList;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Random;
import java.util.Set;

/**
 * a class for Hough transforms of simple geometric shapes.  currently a line
 * is implemented.
 * 
 * Helpful in starting this was a look at the code available from
 * http://vase.essex.ac.uk/software/HoughTransform/
 * and 
 * http://homepages.inf.ed.ac.uk/rbf/HIPR2/flatjavasrc/Hough.java
 * though this implementation is different.
 * 
 * @author nichole
 */
public class HoughTransform {

    public HoughTransform() {
    }

    /**
     * given an edge of points, computes the Hough
     * transform of lines and returns results as an associate array with 
     * key = pair with x = polar theta in degrees and y = distance from
     * the origin in pixels; value = number of transformation points having
     * the key.   Note that the angle is calculated for expectations of a
     * counter clockwise ordered curve and the vector of the angle is 
     * perpendicular to p1 (direction given by right hand rule).
     * The angles are 0 to 360.
     * 
     * Note that if the edge has less than 3 points, an empty map is returned.
     * 
     * runtime complexity is O(N_edge_pts), but includes transcendental operations.
     * 
     * @param edge a curve defined by the points within
     * @param imageWidth
     * @param imageHeight
     * @return thetaRadiusPixCoords mappings
     */
    public Map<PairInt, Set<PairInt>> calculateLineGivenEdge(PairIntArray edge,
        int imageWidth, int imageHeight) {
        
        Map<PairInt, Set<PairInt>> outputPolarCoordsPixMap = new HashMap<PairInt, Set<PairInt>>();
        
        if (edge.getN() < 3) {
            return outputPolarCoordsPixMap;
        }
        
        // theta is 0 to 360
        Map<Integer, Double> cosineMap = Misc.getCosineThetaMapForTwoPI();
        Map<Integer, Double> sineMap = Misc.getSineThetaMapForTwoPI();
        
        boolean curveIsClosed = (edge instanceof PairIntArrayWithColor) &&
            (((PairIntArrayWithColor)edge).isClosedCurve()) &&
            (edge.getN() > 2);

        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
                    
        int n = edge.getN();
        
        for (int i = 0; i < n; ++i) {
            
            int x = edge.getX(i);
            int y = edge.getY(i);
            
            int xp, yp, xn, yn;
            
            if (i == 0) {
                if (curveIsClosed) {
                    xp = edge.getX(n - 1);
                    yp = edge.getY(n - 1);
                } else {
                    // use replication for boundary
                    xp = x;
                    yp = y;
                }
                xn = edge.getX(i + 1);
                yn = edge.getY(i + 1);
            } else if (i == (n - 1)) {
                xp = edge.getX(i - 1);
                yp = edge.getY(i - 1);
                if (curveIsClosed) {
                    xn = edge.getX(0);
                    yn = edge.getY(0);
                } else {
                    xn = x;
                    yn = y;
                }
            } else {
                xp = edge.getX(i - 1);
                yp = edge.getY(i - 1);
                
                xn = edge.getX(i + 1);
                yn = edge.getY(i + 1);
            }
            
            // note, this is not the angle along the edge, it's perpendicular
            // to it, but the calculation is consistent
            double t = curveHelper.calculateAngleTangentToMidpoint(xp, yp, x, y, 
                xn, yn); 
            
            double tDegrees = t * 180./Math.PI;
            
            int tInt = (int)Math.round(tDegrees);
            
            if (tInt > 359) {
                tInt = tInt - 360;
            }
            
            Integer theta = Integer.valueOf(tInt);

            double ct = cosineMap.get(theta).doubleValue();
            double st = sineMap.get(theta).doubleValue();

            double r = Math.sqrt(Math.pow((x * ct), 2) + Math.pow((y * st), 2));

            PairInt p = new PairInt(tInt, (int)Math.round(r));

            Set<PairInt> set = outputPolarCoordsPixMap.get(p);
            if (set == null) {
                set = new HashSet<PairInt>();
                outputPolarCoordsPixMap.put(p, set);
            }
            set.add(new PairInt(x, y));
        }
        
        return outputPolarCoordsPixMap;
    }
   
    /**
     * given the theta values and a set of point coordinates, returns a map
     * of the angle and distance from the image origin.
     * 
     * @param points
     * @param theta360
     * @return 
     */
    public Map<PairInt, PairInt> calculateHoughTransforms(
        Set<PairInt> points, GreyscaleImage theta360) {
        
        // theta is 0 to 360
        Map<Integer, Double> cosineMap = Misc.getCosineThetaMapForTwoPI();
        Map<Integer, Double> sineMap = Misc.getSineThetaMapForTwoPI();
                        
        Map<PairInt, PairInt> pointTRMap = new HashMap<PairInt, PairInt>();
                
        for (PairInt p : points) {
        
            int x = p.getX();
            int y = p.getY();
            int t = theta360.getValue(x, y);
            
            /*if ((t - 90) > 0) {
                t -= 90;
            } else if ((t + 90) < 360) {
                t += 90;
            }*/
           
            Integer theta = Integer.valueOf(t);
            
            double ct = cosineMap.get(theta).doubleValue();
            double st = sineMap.get(theta).doubleValue();

            double r = Math.sqrt(Math.pow((x * ct), 2) + Math.pow((y * st), 2));

            PairInt pTR = new PairInt(t, (int) Math.round(r));

            pointTRMap.put(p, pTR);
        }
        
        return pointTRMap;
    }
    
    /**
     * 
     * @param xCoords
     * @param yCoords
     * @param theta orientations of xCoords, yCoords in range 0 to 359,
     * inclusive.
     * @return 
     */
    public Map<PairInt, PairInt> calculateHoughTransforms(
        TIntList xCoords, TIntList yCoords, TIntList theta) {
        
        if (xCoords.size() == yCoords.size() || xCoords.size() == 
            theta.size()) {
            throw new IllegalArgumentException("xCoords, yCoords, and theta "
                + "must be same lengths");
        }
        
        // theta is 0 to 360
        Map<Integer, Double> cosineMap = Misc.getCosineThetaMapForTwoPI();
        Map<Integer, Double> sineMap = Misc.getSineThetaMapForTwoPI();
                        
        Map<PairInt, PairInt> pointTRMap = new HashMap<PairInt, PairInt>();
                
        for (int i = 0; i < xCoords.size(); ++i) {
            int x = xCoords.get(i);
            int y = yCoords.get(i);
            int t = theta.get(i);
            
            Integer th = Integer.valueOf(t);
            
            double ct = cosineMap.get(th).doubleValue();
            double st = sineMap.get(th).doubleValue();

            double r = Math.sqrt(Math.pow((x * ct), 2) + Math.pow((y * st), 2));

            PairInt pTR = new PairInt(t, (int) Math.round(r));

            pointTRMap.put(new PairInt(x, y), pTR);
        }
        
        return pointTRMap;
    }
           
    protected PairIntArray createLine(int len, int xOff, int yOff) {
        PairIntArray a = new PairIntArray(len);
        for (int i = 0; i < len; ++i) {
            a.add(xOff + i, yOff + i);
        }
        return a;
    }
    
    /**
     * runtime complexity is O(N * lg_2(N)).
     * @param thetaRadiusPixMap
     * @return 
     */
    public List<PairInt> sortByVotes(Map<PairInt, Set<PairInt>> thetaRadiusPixMap) {
        
        int[] votes = new int[thetaRadiusPixMap.size()];
        int[] indexes = new int[votes.length];
        PairInt[] keys = new PairInt[votes.length];
        
        int count = 0;
        for (Entry<PairInt, Set<PairInt>> entry : thetaRadiusPixMap.entrySet()) {
            votes[count] = entry.getValue().size();
            keys[count] = entry.getKey();
            indexes[count] = count;
            count++;
        }
        
        MultiArrayMergeSort.sortByDecr(votes, indexes);

        List<PairInt> outSortedKeys = new ArrayList<PairInt>();
        
        for (int i = 0; i < indexes.length; ++i) {
            int idx = indexes[i];            
            outSortedKeys.add(keys[idx]);
        }
        
        return outSortedKeys;
    }

    /**
     * sets true and false in present array for presence in points set.
       <pre>
          19 20 21 22 23
          14 15 16 17 18
          10 11  @ 12 13    
           5  6  7  8  9       
           0  1  2  3  4
       </pre>
     * @param x
     * @param y
     * @param points
     * @param present 
    */
    private void populateNeighborhoodPresence(int x, int y, Set<PairInt> points, 
        boolean[] present) {
        
        int count = 0;
        for (int dy = -2; dy <= 2; ++dy) {
            int y2 = y + dy;
            for (int dx = -2; dx <= 2; ++dx) {
                if (dx == 0 && dy == 0) {
                    continue;
                }
                int x2 = x + dx;
                PairInt p = new PairInt(x2, y2);
                if (points.contains(p)) {
                    present[count] = true;
                } else {
                    present[count] = false;
                }
                ++count;
            }
        }
    }

    /**
     *  uses thiel sen estimator to calculate the slope and y intercept,
     * then calculates the average distance of the points from the line
     * and then the standard deviation of that average.  then removes
     * outliers that are further from the line than meanTolerance
     * and then recalculates.  Note that the argument set points is modified
     * to remove outliers.
     * @param linePoints input and output modified set
     * @param roughTheta
     * @return
     */
    private float[] calculateReducedLinePolarCoordsAndStats(
        Set<PairInt> linePoints, Integer roughTheta) {
        
        PairFloatArray xy = new PairFloatArray();
        
        for (PairInt p : linePoints) {
            xy.add(p.getX(), p.getY());
        }
            
        LinearRegression lReg = new LinearRegression();
        //lReg.plotTheLinearRegression(x, y);
        
        float[] yInterceptAndSlope;
        
        if (xy.getN() > 3162) {
        //if (xy.getN() > 46340) {
            // need to reduce the number of points
            //TODO: replace with faster sampling
            Random random = new Random(System.nanoTime());
            int n = 3162;
            float[] xs = new float[n];
            float[] ys = new float[n];
            Set<Integer> chosen = new HashSet<Integer>();
            while (chosen.size() < n) {
                int idx = random.nextInt(n);
                Integer index = Integer.valueOf(idx);
                if (!chosen.contains(index)) {
                    xs[chosen.size()] = xy.getX(idx);
                    ys[chosen.size()] = xy.getY(idx);
                    chosen.add(index);
                }
            }
            yInterceptAndSlope = lReg.calculateTheilSenEstimatorParams(xs, ys);
        } else {            
            float[] x = Arrays.copyOf(xy.getX(), xy.getN());
            float[] y = Arrays.copyOf(xy.getY(), xy.getN());
            yInterceptAndSlope = lReg.calculateTheilSenEstimatorParams(x, y);
        }
        
        double thetaRadians;
        double rSum = 0;
        
        boolean isVertical = (yInterceptAndSlope[1] == Float.MAX_VALUE);
        boolean isHorizontal = (yInterceptAndSlope[1] == 0); 
        
        if (isVertical) {
            thetaRadians = Math.PI/2.;
            int count = 0;
            for (PairInt p : linePoints) {
                if ((count & 1) == 1) {
                    continue;
                }
                rSum += p.getX();
                ++count;
            }
            rSum /= count;
        } else if (isHorizontal) {
            thetaRadians = 0;
            int count = 0;
            for (PairInt p : linePoints) {
                if ((count & 1) == 1) {
                    continue;
                }
                rSum += p.getY();
                ++count;
            }
            rSum /= count;
        } else {
            thetaRadians = Math.atan(yInterceptAndSlope[1]);
            double ct = Math.cos(thetaRadians);
            double st = Math.sin(thetaRadians);
            int count = 0;
            for (PairInt p : linePoints) {
                if ((count & 1) == 1) {
                    continue;
                }
                rSum += (p.getX() * ct) + (p.getY() * st);
                ++count;
            }
            rSum /= count;
        }
        
        int t = (int)Math.round(thetaRadians*180./Math.PI);
        if (t < 0) {
            t += 360;
        }
        
        // ----- evaluate ------
        
        /*
        have y = m * x + k  where m is slope and k is yIntercept

        point (x0, y0)

        d = sqrt( 
            (((x0 + m*y0 - m*k)/(m*m + 1)) - x0)^2 +
            (((m*(x0 + m*y0 - m*k))/(m*m + 1)) + k - y0)^2
        )
        */
        
        float m = yInterceptAndSlope[1];
        float k = yInterceptAndSlope[0];
        float mn, stdv;
            
        double sum = 0;
        int count = 0;
        double[] d = new double[linePoints.size()];
        List<PairInt> points = new ArrayList<PairInt>(linePoints);
        for (PairInt p : points) {
            int x0 = p.getX();
            int y0 = p.getY();
            if (isVertical) {
                d[count] = x0 - rSum;
            } else if (isHorizontal) {
                d[count] = y0 - rSum;
            } else {
                /*
                m*x - y + yIntercept = 0
                */
                d[count] = Math.sqrt(Math.abs(m*x0 - y0 + k)/(m*m + 1));  
            }
            sum += d[count];
            count++;
        }
        sum /= (double)linePoints.size();
        mn = (float)sum;
        sum = 0;
        for (int i = 0; i < d.length; ++i) {
            double diff = d[i] - mn;
            sum += (diff * diff);
            count++;
        }
        stdv = (float)Math.sqrt(sum/((double)linePoints.size() - 1.));        
        
        if (stdv == 0) {
            return new float[]{t, (float)Math.abs(rSum), Math.abs(mn), stdv};
        }
        
        // TODO: might want to allow a tolerance to be passed in
        double limit = stdv;
        // if it's highly inclined line, tolerance has to allow for step width
        if ((Math.abs(AngleUtil.getAngleDifference(t, 22.5f)) < 12) || 
            (Math.abs(AngleUtil.getAngleDifference(t, 67.5f)) < 12) ||
            (Math.abs(AngleUtil.getAngleDifference(t, 112.5f)) < 12) ||
            (Math.abs(AngleUtil.getAngleDifference(t, 157.5f)) < 12) ||
            (Math.abs(AngleUtil.getAngleDifference(t, 202.5f)) < 12) ||
            (Math.abs(AngleUtil.getAngleDifference(t, 247.5f)) < 12) ||
            (Math.abs(AngleUtil.getAngleDifference(t, 292.5f)) < 12) ||
            (Math.abs(AngleUtil.getAngleDifference(t, 337.5f)) < 12)
            ) {
            //TODO:
            // this picks up the meeting of lines at a corner, so
            // may need revision
            limit = 2*stdv;
        }
        if (limit < 1) {
            limit = 1;
        }
        
        Set<Integer> remove = new HashSet<Integer>();
        for (int i = 0; i < d.length; ++i) {
            double diff = Math.abs(d[i] - mn);
            if (diff > limit) {
                remove.add(Integer.valueOf(i));
            }
        }
        
        // --- recalc mean and stdev w/o outliers ----
  //TODO: consider recalculating the slope and intercept w/ thiel sen
        sum = 0;
        count = 0;
        for (int i = 0; i < d.length; ++i) {
            if (remove.contains(Integer.valueOf(i))) {
                linePoints.remove(points.get(i));
                continue;
            }
            sum += d[i];
            count++;
        }
        sum /= (double)count;
        mn = (float)sum;
        count = 0;
        sum = 0;
        for (int i = 0; i < d.length; ++i) {
            if (remove.contains(Integer.valueOf(i))) {
                continue;
            }
            double diff = Math.abs(d[i] - mn);
            sum += (diff * diff);
            count++;
        }
        stdv = (float)Math.sqrt(sum/((double)count - 1.));        
        return new float[]{t, (float)Math.abs(rSum), Math.abs(mn), stdv};
    }
    
    public class HoughTransformLines {
        
        private final Map<PairInt, PairInt> pixelToPolarCoordMap;
        
        private final List<Set<PairInt>> sortedLineGroups;
        
        public HoughTransformLines(Map<PairInt, PairInt> pixToTRMap,
            List<Set<PairInt>> sortedGroups) {
            this.pixelToPolarCoordMap = pixToTRMap;
            this.sortedLineGroups = sortedGroups;
        }
        
        public Map<PairInt, PairInt> getPixelToPolarCoordMap() {
            return pixelToPolarCoordMap;
        }
        
        public List<Set<PairInt>> getSortedLineGroups() {
            return sortedLineGroups;
        }
    }
   
}
