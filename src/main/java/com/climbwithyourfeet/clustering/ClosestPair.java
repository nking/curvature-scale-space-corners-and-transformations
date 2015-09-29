package com.climbwithyourfeet.clustering;

import com.climbwithyourfeet.clustering.util.PairFloat;
import com.climbwithyourfeet.clustering.util.Sorter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Find the closest pair of points in a set of points by using divide and 
 * conquer to reduce the number of comparisons
 * 
 * From pseudocode in Intro to Algorithms by Cormen et al.
 * 
 * runtime complexity is ~ O(N lg N).
 *
 * @author nichole
 */
public class ClosestPair {

    private boolean debug = false;
    
    /**
     *
     */
    public ClosestPair() {
    }
    
    /**
     *
     */
    public void setDebug() {
        this.debug = true;
    }

    /**
     * find the closest pair within the x,y set of points.
     * note that the given arrays will have their item orders altered.
     *
     * @param xPoints
     * @param yPoints
     * @return
     *     returns 5 items in array:  shortest distance between pairs, pt1 x, pt1 y, pt2 x, pt2 y
     */
    public ClosestPairFloat findClosestPair(List<Float> xPoints, List<Float> yPoints) {
        
        if (xPoints == null) {
	    	throw new IllegalArgumentException("xpoints cannot be null");
        }
	    if (yPoints == null) {
	    	throw new IllegalArgumentException("ypoints cannot be null");
        }
	    if (xPoints.size() != yPoints.size()) {
	    	throw new IllegalArgumentException(
                "xpoints must have the same number of items as y");
        }
        
        List<PairFloat> p = new ArrayList<PairFloat>();
        List<PairFloat> x = new ArrayList<PairFloat>();
        List<PairFloat> y = new ArrayList<PairFloat>();
        
        for (int i = 0; i < xPoints.size(); i++) {
            float xPoint = xPoints.get(i);
            float yPoint = yPoints.get(i);
            // since the points xy are not going to be modified, can reuse xy
            PairFloat xy = new PairFloat(xPoint, yPoint);
            p.add(xy);
            y.add(xy);
        }
        Sorter.mergeSortByXThenY(p);
        x.addAll(p);
        
        Sorter.mergeSortByYThenX(y);
        
        return divide(p, x, y);
    }
                
    /**
     * use divide and conquer to find the pair of points with the smallest 
     * separation.  note that x and y have to contain objects that are the
     * same instance as the objects in p.
     * 
     * @param p should be ordered left to right already
     * @param x
     * @param y
     * @return 
     */
    protected ClosestPairFloat divide(List<PairFloat> p, List<PairFloat> x, 
        List<PairFloat> y) {
        
        if (p.size() <= 3) {
            return bruteForceMinDistance(p);
        }
                            
        int q = p.size() >> 1; // 0 1 *2*  3 4 5

        List<PairFloat> pL = new ArrayList<PairFloat>(q);
        List<PairFloat> pR = new ArrayList<PairFloat>(q);
        
        // need a datastructure to search for pL members.  for a hash table, search is O(1).
        //  will assume that the java HashSet implementation has O(1) insert and search
        Set<PairFloat> isInL = new HashSet<PairFloat>();
        
        List<PairFloat> xL = new ArrayList<PairFloat>(q);
        List<PairFloat> xR = new ArrayList<PairFloat>(q);
        List<PairFloat> yL = new ArrayList<PairFloat>(q);
        List<PairFloat> yR = new ArrayList<PairFloat>(q);
        
        // p is already sorted by x, so fill in the subsets for p and x
        for (int i = 0; i < q; ++i) {
            PairFloat pi = p.get(i);
            pL.add(pi);
            xL.add(pi);
            isInL.add(pi);
        }
        for (int i = q; i < p.size(); ++i) {
            PairFloat pi = p.get(i);
            pR.add(pi);
            xR.add(pi);
        }
        
        for (int i = 0; i < y.size(); ++i) {
            PairFloat yi = y.get(i);
            if (isInL.contains(yi)) {
                yL.add(yi);
            } else {
                yR.add(yi);
            }
        }
        
        // use 2 recursive calls for conquer

        // find closest pair in pL w/ pL, xL, and yL
        ClosestPairFloat cpL = divide(pL, xL, yL);

        // find closest pair in pP w/ pR, xR, and yR
        ClosestPairFloat cpR = divide(pR, xR, yR);

        ClosestPairFloat d = (cpL.separation <= cpR.separation) ? cpL : cpR;
        
        return combine(yL, yR, d, xL.get(xL.size() - 1), xR.get(0));            
    }
    
    /**
     * runtime complexity:
     *    O(N) plus a small fraction of points contributing O(m lg m) for a sort
     * @param yL
     * @param yR
     * @param d
     * @param leftMostL
     * @param rightMostR
     * @return 
     */
    protected ClosestPairFloat combine(List<PairFloat> yL, List<PairFloat> yR, 
        ClosestPairFloat d,  PairFloat leftMostL, PairFloat rightMostR) {
        
        List<PairFloat> yPrime = new ArrayList<PairFloat>();
        
        float delta = d.separation;
        
        // find the points in yL that have x within delta from leftMostL
        // traversing yL (rather than xL) to keep yPrime ordered
        for (int i = 0; i < yL.size(); ++i) {
            float xPoint = yL.get(i).getX();
            if ((leftMostL.getX() - xPoint) <= delta) {
                yPrime.add(yL.get(i));
            }
        }
        
        // find the points in yR that have x within delta from rightMostR
        // traversing yR (rather than xR) to keep yPrime ordered
        for (int i = 0; i < yR.size(); ++i) {
            float xPoint = yR.get(i).getX();
            if ((xPoint - rightMostR.getX()) <= delta) {
                yPrime.add(yR.get(i));
            }
        }
        
        int idx0 = -1;
        int idx1 = -1;
        double minDistSq = Double.MAX_VALUE;
        
        for (int i = 0; i < yPrime.size(); i++) {
            
            PairFloat yi = yPrime.get(i);
            
            for (int j = (i + 1); j < (i + 8); j++) {
                if (j > (yPrime.size() - 1)) {
                    break;
                }
                
                PairFloat yj = yPrime.get(j);
                
                float diffX = yi.getX() - yj.getX();
                float diffY = yi.getY() - yj.getY();
                
                double distSq = (diffX * diffX) + (diffY * diffY);
                if (distSq < minDistSq) {
                    minDistSq = distSq;
                    idx0 = i;
                    idx1 = j;
                }
            }
        }
       
        if (minDistSq == Double.MAX_VALUE) {
            return d;
        }
        
        float dist = (float)Math.sqrt(minDistSq);
        
        if (dist < delta) {
            return new ClosestPairFloat(yPrime.get(idx0), yPrime.get(idx1), dist);
        } else {
            return d;
        }        
    }
    
    /**
     *
     * @param p
     * @return
     */
    protected ClosestPairFloat bruteForceMinDistance(List<PairFloat> p) {
        
        if (p.size() < 1) {
            return new ClosestPairFloat(null, null, Float.MAX_VALUE);
        }
        
        double minDistSq = Double.MAX_VALUE;
        int idx0 = -1;
        int idx1 = -1;

        for (int i = 0; i < p.size(); i++) {
            
            PairFloat pi = p.get(i);
            
            for (int j = i; j < p.size(); j++) {
                if (i == j) {
                    continue;
                }
                
                PairFloat pj = p.get(j);
                
                float diffX = pi.getX() - pj.getX();
                float diffY = pi.getY() - pj.getY();
                double distSq = (diffX * diffX) + (diffY * diffY);

                if (distSq < minDistSq) {
                    minDistSq = distSq;
                    idx0 = i;
                    idx1 = j;
                }
            }
        }
        
        float minDist = (float)Math.sqrt(minDistSq);
        
        return new ClosestPairFloat(p.get(idx0), p.get(idx1), minDist);
    }

    /**
     *
     */
    public class ClosestPairFloat {
        PairFloat point0;
        PairFloat point1;
        float separation;

        /**
         *
         * @param p0
         * @param p1
         * @param sep
         */
        public ClosestPairFloat(PairFloat p0, PairFloat p1, float sep) {
            this.point0 = p0;
            this.point1 = p1;
            this.separation = sep;
        }
    }
    
}
