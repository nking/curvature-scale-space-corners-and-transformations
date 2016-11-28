package algorithms.compGeometry;

import algorithms.compGeometry.convexHull.GrahamScan;
import algorithms.compGeometry.convexHull.GrahamScanTooFewPointsException;
import algorithms.imageProcessing.util.PairIntWithIndex;
import algorithms.util.PairInt;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * find the furthest pair between points by finding the upper and lower
 * convex hull and then compare the distances between upper and lower hull
 * points.
 *
 * @author nichole
 */
public class FurthestPair {
    
    public PairInt[] find(Set<PairInt> points) {
        
        int n = points.size();
        
        if (n < 2) {
            return null;
        } else if (n == 2) {
            Iterator<PairInt> iter = points.iterator();
            return new PairInt[]{iter.next(), iter.next()};
        }
        
        float[] x = new float[n];
        float[] y = new float[n];
        
        int count = 0;
        for (PairInt p : points) {
            x[count] = p.getX();
            y[count] = p.getY();
            count++;
        }
        
        GrahamScan scan = new GrahamScan();
        
        try {
            
            scan.computeHull(x, y);
            
        } catch (GrahamScanTooFewPointsException ex) {
            
            // this can happen for colinear points reducting sorted size < 3
            //Logger.getLogger(FurthestPair.class.getName()).log(Level.SEVERE, null, ex);
            
            return findWithBruteForce(points);
        }
        
        float[] xHull = scan.getXHull();
        float[] yHull = scan.getYHull();
                
        float maxDistSq = Float.MIN_VALUE;
        int idx1 = -1;
        int idx2 = -1;
    
        for (int i = 0; i < xHull.length; ++i) {
            for (int j = (i + 1); j < xHull.length; ++j) {
            
                float distSq = distanceSq(xHull[i], yHull[i], xHull[j], yHull[j]);

                if (distSq > maxDistSq) {
                    maxDistSq = distSq;
                    idx1 = i;
                    idx2 = j;
                }
            }
        }
        
        if (idx1 != -1) {
            PairInt p0 = new PairInt(Math.round(xHull[idx1]), Math.round(yHull[idx1]));
            PairInt p1 = new PairInt(Math.round(xHull[idx2]), Math.round(yHull[idx2]));
            return new PairInt[]{p0, p1};
        }
        
        return null;
    }
    
    int distanceSq(PairIntWithIndex p1, PairIntWithIndex p2) {
        
        int dx = p1.getX() - p2.getX();
        int dy = p1.getY() - p2.getY();
        
        return (dx*dx + dy*dy);
    }
    int distanceSq(PairInt p1, PairInt p2) {
        
        int dx = p1.getX() - p2.getX();
        int dy = p1.getY() - p2.getY();
        
        return (dx*dx + dy*dy);
    }
    float distanceSq(float x1, float y1, float x2, float y2) {
        
        float dx = x1 - x2;
        float dy = y1 - y2;
        
        return (dx*dx + dy*dy);
    }
    
    PairInt[] findWithBruteForce(Set<PairInt> points) {
        
        if (points.size() < 2) {
            return null;
        }
        
        float max = Float.MIN_VALUE;
        PairInt maxP1 = null;
        PairInt maxP2 = null;
        
        List<PairInt> list = new ArrayList<PairInt>(points);
        for (int i = 0; i < list.size(); ++i) {
            PairInt p1 = list.get(i);
            for (int j = (i + 1); j < list.size(); ++j) {
                PairInt p2 = list.get(j);
                float d = distanceSq(p1, p2);
                if (d > max) {
                    max = d;
                    maxP1 = p1;
                    maxP2 = p2;
                }
            }
        }
        
        return new PairInt[]{maxP1, maxP2};
    }
}
