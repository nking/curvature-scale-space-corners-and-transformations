package algorithms.imageProcessing.features;

import algorithms.compGeometry.PerimeterFinder;
import algorithms.imageProcessing.MiscellaneousCurveHelper;
import algorithms.imageProcessing.ZhangSuenLineThinner;
import algorithms.misc.MiscMath;
import algorithms.util.PairInt;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * class to encapsulate methods to make a rough skeletonization of a set of
 * points (using line thinning) and methods to find the closest skeleton
 * point to a query point.  The class currently expects to handle a list of
 * such points to make an internal resulting list of skeletons 
 * (a.k.a. medial axes for each point set in a list).
 * 
 * (for example, this instance is constructed from points which are a larger
 * set of blobs in a hierarchy of segmentation.  this instance and the bounds
 * are used to merge the groups of points from a finer segmentation into 
 * the groups defined by larger bounds.  the resulting merged list is then
 * ordered to have same order as the internal data here and then both are
 * filtered to remove empty sets.  thereafter, the medial axes can be used
 * complementarily with the grouped point lists for functions that need to know
 * the direction of "inward" for the group of points.)
 * 
 * @author nichole
 */
public class BlobMedialAxes {
    
    private List<Map<Integer, List<Integer>>> skeletonXMapList = null;
    private List<Map<Integer, List<Integer>>> skeletonYMapList = null;
    private double[][] xyCentroids = null;
    private final List<Double> lColorList;
    private final List<Double> aColorList;
    private final List<Double> bColorList;
    
    public BlobMedialAxes(final List<Set<PairInt>> blobs, 
        final List<Double> lClrList, final List<Double> aClrList,
        final List<Double> bClrList) {
        
        int n = blobs.size();    
        xyCentroids = new double[n][2];
        skeletonXMapList = new ArrayList<Map<Integer, List<Integer>>>(n);
        skeletonYMapList = new ArrayList<Map<Integer, List<Integer>>>(n);
        
        this.lColorList = new ArrayList<Double>(lClrList);
        this.aColorList = new ArrayList<Double>(aClrList);
        this.bColorList = new ArrayList<Double>(bClrList);
        
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        
        ZhangSuenLineThinner lt = new ZhangSuenLineThinner();
        
        // make skeleton for each blob for detailed perimeter "inward" directions
        for (int i = 0; i < blobs.size(); ++i) {
            
            Set<PairInt> blob = blobs.get(i);
            Set<PairInt> skeleton = new HashSet<PairInt>(blob);
            
            int[] minMaxXY = MiscMath.findMinMaxXY(skeleton);
            lt.applyLineThinner(skeleton, minMaxXY[0], minMaxXY[1], minMaxXY[2],
                minMaxXY[3]);
            int[] xPoints = new int[skeleton.size()];
            int[] yPoints = new int[skeleton.size()];

            int count = 0;
            for (PairInt p : skeleton) {
                xPoints[count] = p.getX();
                yPoints[count] = p.getY();
                count++;
            }
            
            double[] xyCen = curveHelper.calculateXYCentroids(blob);
            
            // order skeleton by x and then y
            Map<Integer, List<Integer>> xSkeletonMap = makeXMap(xPoints, yPoints);
            Map<Integer, List<Integer>> ySkeletonMap = makeYMap(xPoints, yPoints);
            
            xyCentroids[i] = xyCen;
            skeletonXMapList.add(xSkeletonMap);
            skeletonYMapList.add(ySkeletonMap);
        }
    }
    
    public int getNumberOfItems() {
        return xyCentroids.length;
    }
    
    public PairInt findClosestPoint(int index, int x, int y) {
        
        if (index < 0 || index > (skeletonXMapList.size() - 1)) {
            throw new IllegalArgumentException("index is out of bounds");
        }
        
        return findClosestPoint(x, y, skeletonXMapList.get(index), skeletonYMapList.get(index));
    }
    
    protected PairInt findClosestPoint(int x, int y, 
        Map<Integer, List<Integer>> xSkeletonMap,
        Map<Integer, List<Integer>> ySkeletonMap) {
        
        List<Integer> ys = xSkeletonMap.get(Integer.valueOf(x));
        int ySkel1 = Integer.MAX_VALUE;
        if (ys != null) {
            int minDistY2 = Integer.MAX_VALUE;
            for (Integer y0 : ys) {
                int distYSq = y0.intValue() - y;
                distYSq *= distYSq;
                if (distYSq < minDistY2) {
                    minDistY2 = distYSq;
                    ySkel1 = y0.intValue();
                }
            }
        }
        int xSkel1 = x;
        List<Integer> xs = ySkeletonMap.get(Integer.valueOf(x));
        int xSkel2 = Integer.MAX_VALUE;
        if (xs != null) {
            int minDistX2 = Integer.MAX_VALUE;
            for (Integer x0 : xs) {
                int distXSq = x0.intValue() - x;
                distXSq *= distXSq;
                if (distXSq < minDistX2) {
                    minDistX2 = distXSq;
                    xSkel2 = x0.intValue();
                }
            }
        }
        int ySkel2 = y;
        int xSkel, ySkel;
        if ((((xSkel1 - x)*(xSkel1 - x)) + ((ySkel1 - y)*(ySkel1 - y)))
            <
            (((xSkel2 - x)*(xSkel2 - x)) + ((ySkel2 - y)*(ySkel2 - y)))) {
            xSkel = xSkel1;
            ySkel = ySkel1;
        } else {
            xSkel = xSkel2;
            ySkel = ySkel2;
        }
        
        return new PairInt(xSkel, ySkel);
    }
    
    public double[] getOriginalBlobXYCentroid(int index) {
        
        if (index < 0 || index > (skeletonXMapList.size() - 1)) {
            throw new IllegalArgumentException("index is out of bounds");
        }
        
        double[] c = Arrays.copyOf(xyCentroids[index], 2);
        
        return c;
    }
    
    /**
     * get the LAB color space averaged colors for the points within the 
     * bounds at index in lists.
     * @param index
     * @return 
     */
    public float[] getLABColors(int index) {
        
        if (index < 0 || index > (skeletonXMapList.size() - 1)) {
            throw new IllegalArgumentException("index is out of bounds");
        }
        
        float[] c = new float[3];
        c[0] = lColorList.get(index).floatValue();
        c[1] = aColorList.get(index).floatValue();
        c[2] = bColorList.get(index).floatValue();
        
        return c;
    }
    
    /**
     * map w/ key = x and y = ascending ordered values for that x
     * @param x
     * @param y
     * @return 
     */
    protected Map<Integer, List<Integer>> makeXMap(int[] x, int[] y) {
        
        Map<Integer, List<Integer>> map = new HashMap<Integer, List<Integer>>();
        for (int i = 0; i < x.length; ++i) {
            Integer key = Integer.valueOf(x[i]);
            List<Integer> list = map.get(key);
            if (list == null) {
                list = new ArrayList<Integer>();
                map.put(key, list);
            }
            list.add(Integer.valueOf(y[i]));
        }
        
        for (Map.Entry<Integer, List<Integer>> entry : map.entrySet()) {
            List<Integer> ys = entry.getValue();
            if (ys.size() > 1) {
                Collections.sort(ys);
            }
        }
        
        return map;
    }
    
    /**
     * map w/ key = y and x = ascending ordered values for that y
     * @param x
     * @param y
     * @return 
     */
    protected Map<Integer, List<Integer>> makeYMap(int[] x, int[] y) {
        
        Map<Integer, List<Integer>> map = new HashMap<Integer, List<Integer>>();
        for (int i = 0; i < y.length; ++i) {
            Integer key = Integer.valueOf(y[i]);
            List<Integer> list = map.get(key);
            if (list == null) {
                list = new ArrayList<Integer>();
                map.put(key, list);
            }
            list.add(Integer.valueOf(x[i]));
        }
        
        for (Map.Entry<Integer, List<Integer>> entry : map.entrySet()) {
            List<Integer> xs = entry.getValue();
            if (xs.size() > 1) {
                Collections.sort(xs);
            }
        }
        
        return map;
    }

    /**
     * rewrite the internal data structures to remove the subset of indexes
     * given.  Note that the internal skeleton
     * maps are for now, not trimmed to the smaller subset of points as the
     * larger set - that is not harmful to the internal functions.
     * @param removeIndexes an ascending list of unique indexes to remove
     */
    public void removeIndexes(List<Integer> removeIndexes) {
        
        // NOTE: should consider retaining blobs in constructor so can
        // rebuild the skeleton maps here from the reduced set of points.
        // for now, deciding that the extra information in the skeleton maps 
        // doesn't harm the use of the structures.
        
        Set<Integer> exclude = new HashSet<Integer>(removeIndexes);
        
        int n2 = xyCentroids.length - removeIndexes.size();
        
        List<Double> ell = new ArrayList<Double>();
        List<Double> a = new ArrayList<Double>();
        List<Double> b = new ArrayList<Double>();
        
        double[][] xyCentroid2 = new double[n2][2];
        int count = 0;
        for (int i = 0; i < xyCentroids.length; ++i) {
            
            Integer index = Integer.valueOf(i);
            if (!exclude.contains(index)) {
                xyCentroid2[count] = Arrays.copyOf(this.xyCentroids[i], 2);
                ell.add(lColorList.get(i));
                a.add(aColorList.get(i));
                b.add(bColorList.get(i));
                count++;
            }
        }   
        
        this.xyCentroids = xyCentroid2;
        
        this.lColorList.clear();
        this.lColorList.addAll(ell);
        
        this.aColorList.clear();
        this.aColorList.addAll(a);
        
        this.bColorList.clear();
        this.bColorList.addAll(b);
    }
    
}
