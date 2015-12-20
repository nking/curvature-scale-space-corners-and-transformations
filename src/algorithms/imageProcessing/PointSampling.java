package algorithms.imageProcessing;

import algorithms.compGeometry.PointPartitioner;
import algorithms.compGeometry.PointPartitioner.Bounds;
import algorithms.misc.MiscMath;
import algorithms.util.PairInt;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

/**
 *
 * @author nichole
 */
public class PointSampling {
    
    /**
     * Given a set of points, determine cells across the data (roughly k cells 
     * for each dimension or roughly equal area) to create seeds, then create a 
     * list of distance based values (distance of points from the nearest seeds) 
     * that is intended to be used for random selection of points from the total 
     * set in a spatially stratified manner, biased towards points closer to the 
     * seed centers.
     * 
     * @param points 
     * @param numCellsPerDimensions 
     * @return  
     */
    public PointValueDistr createSpatialInvDistBasedValues(Set<PairInt> points, 
        int numCellsPerDimensions) {
        
        return createSpatialDistBasedValues(points, numCellsPerDimensions, true);
    }
    
    /**
     * Given a set of points, determine cells across the data (roughly k cells 
     * for each dimension or roughly equal area) to create seeds, then create a 
     * list of distance based values (distance of points from the nearest seeds) 
     * that is intended to be used for random selection of points from the total 
     * set in a spatially stratified manner, biased towards points closer to the 
     * seed centers.
     * 
     * @param points 
     * @param numCellsPerDimensions 
     * @return  
     */
    public PointValueDistr createSpatialDistBasedValues(Set<PairInt> points, 
        int numCellsPerDimensions) {
        
        return createSpatialDistBasedValues(points, numCellsPerDimensions, false);
    }
    
    
    /**
     * Given a set of points, determine cells across the data (roughly k cells 
     * for each dimension or roughly equal area) to create seeds, then create a 
     * list of distance based values (distance of points from the nearest seeds) 
     * that is intended to be used for random selection of points from the total 
     * set in a spatially stratified manner, biased towards points closer to the 
     * seed centers.
     * 
     * @param points 
     * @param numCellsPerDimensions 
     * #param inverse if true, uses maxDist - distance for distances to weight
     * for points closer to the seeds.
     * @return  
     */
    protected PointValueDistr createSpatialDistBasedValues(Set<PairInt> points, 
        int numCellsPerDimensions, boolean inverse) {
        
        if (numCellsPerDimensions < 1) {
            throw new IllegalArgumentException(
                "numCellsPerDimensions must be > 0");
        }
        
        PointPartitioner pp = new PointPartitioner();
        
        List<Bounds> bounds = pp.findCells(numCellsPerDimensions, points);
        
        int[] maxXY = findMaxXAndY(bounds);
        
        List<PairInt> seeds = centers(bounds);
        
        // make an image with points being "1" and seeds being "0" and add
        // the seeds to points
        int w = maxXY[0] + 1;
        int h = maxXY[1] + 1;
        
        Map<PairInt, Integer> pointDistMap = useDistanceTransform(points, w, h, 
            seeds);
        
        int maxDist = Integer.MIN_VALUE;
        if (inverse) {
            for (Entry<PairInt, Integer> entry : pointDistMap.entrySet()) {
                int dist = entry.getValue().intValue();
                if (dist > maxDist) {
                    maxDist = dist;
                }
            }
        }
        
        int n = pointDistMap.size();
        
        BigInteger[] cumulativeCount = new BigInteger[pointDistMap.size()];
        PairInt[] pointsArray = new PairInt[pointDistMap.size()];
        
        int count = 0;
        BigInteger prev = null;
        for (Entry<PairInt, Integer> entry : pointDistMap.entrySet()) {
            
            int value = entry.getValue().intValue();
            if (inverse) {
                value = maxDist + 1 - value;
            }
            
            pointsArray[count] = entry.getKey();
            if (prev == null) {
                prev = BigInteger.ZERO;
            }
            
            BigInteger current = new BigInteger(MiscMath.writeToBigEndianBytes(value));
            
            cumulativeCount[count] = prev.add(current);
            
            prev = cumulativeCount[count];
            
            count++;
        }
        
        PointValueDistr pv = new PointValueDistr(prev, pointsArray, cumulativeCount);
        
        return pv;
    }
    
    protected List<PairInt> centers(List<Bounds> bounds) {
        
        List<PairInt> centers = new ArrayList<PairInt>();
        
        for (Bounds b : bounds) {
            
            float avgX = (b.lowerLeft.getX() + b.lowerRight.getX() 
                + b.upperLeft.getX() + b.upperRight.getX())/4.f;
            
            float avgY = (b.lowerLeft.getY() + b.lowerRight.getY() 
                + b.upperLeft.getY() + b.upperRight.getY())/4.f;
            
            PairInt p = new PairInt(Math.round(avgX), Math.round(avgY));
            
            centers.add(p);
        }
        
        return centers;
    }

    private int[] findMaxXAndY(List<Bounds> bounds) {
        
        int maxX = Integer.MIN_VALUE;
        
        int maxY = Integer.MIN_VALUE;
        
        for (Bounds b : bounds) {
            if (b.lowerRight.getX() > maxX) {
                maxX = b.lowerRight.getX();
            }
            if (b.upperRight.getX() > maxX) {
                maxX = b.upperRight.getX();
            }
            if (b.upperLeft.getY() > maxY) {
                maxY = b.upperLeft.getY();
            }
            if (b.upperRight.getY() > maxY) {
                maxY = b.upperRight.getY();
            }
        }
        
        return new int[]{maxX, maxY};
    }

    private Map<PairInt, Integer> useDistanceTransform(Set<PairInt> points, 
        int w, int h, List<PairInt> seeds) {
             
        com.climbwithyourfeet.clustering.DistanceTransform dt = 
            new com.climbwithyourfeet.clustering.DistanceTransform();
        
        Set<com.climbwithyourfeet.clustering.util.PairInt> seedPoints
            = new HashSet<com.climbwithyourfeet.clustering.util.PairInt>();
        
        for (PairInt p : seeds) {
            com.climbwithyourfeet.clustering.util.PairInt p2 = new
                com.climbwithyourfeet.clustering.util.PairInt(p.getX(), p.getY());
            seedPoints.add(p2);
        }
        
        int[][] distances = dt.applyMeijsterEtAl(seedPoints, w, h);
        
        Map<PairInt, Integer> distMap = new HashMap<PairInt, Integer>();
        for (PairInt p : points) {
            int distSq = distances[p.getX()][p.getY()];
            Integer dist = Integer.valueOf((int)Math.round(Math.sqrt(distSq)));
            distMap.put(p, dist);
        }
        
        return distMap;
    }    
    
}
