package algorithms.compGeometry.clustering;

import algorithms.imageProcessing.GroupPixelRGB0;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.MiscellaneousCurveHelper;
import algorithms.util.PairInt;
import algorithms.util.TrioInt;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import java.awt.Color;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.logging.Logger;

/**
 * k-means clustering is a method of cluster analysis to partition n
 * observations into k clusters in which each observation belongs to the cluster
 * with the nearest mean.
 * This results in a partitioning of the data space into Voronoi cells, which
 * for this single parameter analysis, is 1-D.
 * 
 * This version of k-means needs to be given seeds and hence k from those.
 * 
 * Useful reading:
 * http://en.wikipedia.org/wiki/K-means_clustering
 * 
 * @author nichole
 */
public class KMeansHSV {
    
    protected Logger log = Logger.getLogger(this.getClass().getName());
           
    /**
     * this is k and is chosen by the user
     */
    protected final int nSeeds;
    
    // After final iteration, standard deviations are stored in seedVariances 
    // instead of variances
    protected final static int nMaxIter = 100;
    protected int nIter = 0;
    protected final ImageExt img;
    
    // [nSeeds][h, s, v, x, y]
    protected final float[][] seeds;
    
    /**
     * note that an internal method binPoints makes an assumption that the
     * minimum values in an img pixel is 0 and a maximum is 255.
     * 
     * This class is experimental and specific to a use case of merging
     * segmented cells.
     * 
     * @param starterSeeds
     * @param starterSeedColors
     * @param img
     * 
     */
    public KMeansHSV(List<PairInt> starterSeeds,
        List<GroupPixelRGB0> starterSeedColors, ImageExt img) {
        
        if (starterSeeds.size() != starterSeedColors.size()) {
            throw new IllegalArgumentException("starterSeeds and"
                + " starterSeedColors must be same size");
        }
        
        this.nSeeds = starterSeeds.size();
        this.nIter = 0;
        
        // starter seeds, sorted by increasing value
        // [index][r,g,b, x, y]
        this.seeds = createStartSeeds(img, starterSeeds, starterSeedColors);
        
        this.img = img;
        
        System.out.println("k=" + nSeeds);
    }
        
    /**
     *
     * @param unitSets sets to be merged using k-means.  note that
     * a set is never split into a smaller unit.
     */
    public TIntList computeMeans(List<Set<PairInt>> unitSets) {
    
        // x, y, n
        List<TrioInt> unitCentroids = new ArrayList<TrioInt>();
        
        // r, g, b, n
        List<GroupPixelRGB0> unitClrs = new ArrayList<GroupPixelRGB0>();
        
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        
        for (int i = 0; i < unitSets.size(); ++i) {
            
            Set<PairInt> set = unitSets.get(i);
            
            GroupPixelRGB0 clr = new GroupPixelRGB0();
            clr.calculateColors(set, img, 0, 0);
            unitClrs.add(clr);
        
            double[] xyCen = curveHelper.calculateXYCentroids(set);
            
            unitCentroids.add(new TrioInt((int)Math.round(xyCen[0]),
                (int)Math.round(xyCen[1]), set.size()));
        }
        
        boolean hasConverged = false;
                
        float normXY = img.getNPixels()/nSeeds;
        double maxError = 0.01 * Math.sqrt(this.nSeeds * (3. + normXY));
        
        TIntList unitLabels = null;
        
        while (!hasConverged && (nIter < nMaxIter)) {
            
            unitLabels = binPoints(unitSets, unitClrs, unitCentroids);
            
            double l2Norm = calculateMeanOfSeedPoints(unitSets, unitClrs, 
                unitCentroids, unitLabels);
           
            System.out.println("nIter=" + nIter + " l2Norm=" 
                + l2Norm + " maxError=" + maxError);
            
            if (l2Norm == Double.POSITIVE_INFINITY) {
                break;
            } else if (l2Norm < maxError) {
                break;
            }
    
            nIter++;
        }
        
        return unitLabels;
    }
    
    /**
     * calculate the mean value of all points within a seed bin and return them
     *   as new seed bin centers.  note that if there is a bin without points
     *   in it, null is returned.
     *
     */
    private double calculateMeanOfSeedPoints(List<Set<PairInt>> unitSets, 
        List<GroupPixelRGB0> unitClrs, List<TrioInt> unitCentroids,
        TIntList unitLabels) {

        int[][] sum = new int[nSeeds][5];
        for (int i = 0; i < 3; ++i) {
            sum[i] = new int[5];
        }
        int[] count = new int[nSeeds];
        
        float[] hsv = new float[3];

        for (int i = 0; i < unitSets.size(); i++) {

            /*
            GroupPixelRGB0 clr = unitClrs.get(i);
            Color.RGBtoHSB(Math.round(clr.getAvgRed()), 
                Math.round(clr.getAvgGreen()), 
                Math.round(clr.getAvgBlue()), hsv);
            
            TrioInt xy = unitCentroids.get(i);
            int seedIdx = unitLabels.get(i);

            sum[seedIdx][0] += (xy.getZ() * hsv[0]);
            sum[seedIdx][1] += (xy.getZ() * hsv[1]);
            sum[seedIdx][2] += (xy.getZ() * hsv[2]);
            sum[seedIdx][3] += (xy.getZ() * xy.getX());
            sum[seedIdx][4] += (xy.getZ() * xy.getY());
        
            count[seedIdx] += xy.getZ();
            */
            
            int seedIdx = unitLabels.get(i);
            
            for (PairInt p : unitSets.get(i)) {
                int r = img.getR(p);
                int g = img.getG(p);
                int b = img.getB(p);
                Color.RGBtoHSB(r, g, b, hsv);
                sum[seedIdx][0] += hsv[0];
                sum[seedIdx][1] += hsv[1];
                sum[seedIdx][2] += hsv[2];
                sum[seedIdx][3] += p.getX();
                sum[seedIdx][4] += p.getY();
            }
            
            count[seedIdx] += unitSets.get(i).size();
        }
        
        double l2Norm = 0;
        
        // TODO: adding experimental logic here to stop the merging when
        // a cluster count == 0.
        // the invoker must recognize return value as specific to end of merging.
        
        for (int i = 0; i < nSeeds; i++) {
            for (int j = 0; j < 5; ++j) {
                if (count[i] == 0) {
                    assert(sum[i][j] == 0);
                    return Double.POSITIVE_INFINITY;
                } else {
                    sum[i][j] /= count[i];
                }
                double diff = seeds[i][j] - sum[i][j];
                l2Norm += (diff * diff);
                seeds[i][j] = sum[i][j];
            }
        }
        
        return Math.sqrt(l2Norm);
    }
    
    protected TIntList binPoints(List<Set<PairInt>> unitSets,
        List<GroupPixelRGB0> unitColors, List<TrioInt> unitCentroids) {
        
        TIntList unitLabels = new TIntArrayList();
        
        for (int i = 0; i < unitSets.size(); ++i) {
            
            Set<PairInt> set = unitSets.get(i);
            
            GroupPixelRGB0 clrs = unitColors.get(i);
            
            int seedIdx = findNearestNeighbor(clrs, unitCentroids.get(i));
            
            unitLabels.add(seedIdx);
        }

        return unitLabels;
    }
 
    private int findNearestNeighbor(GroupPixelRGB0 unitColors, 
        TrioInt unitCentroid) {
        
        /*
        [nSeeds][h, s, v, x, y]
        float[][] seeds;
        */
        
        float[] hsv = new float[3];
        Color.RGBtoHSB(Math.round(unitColors.getAvgRed()), 
            Math.round(unitColors.getAvgGreen()), 
            Math.round(unitColors.getAvgBlue()), hsv);
        
        int x = unitCentroid.getX();
        int y = unitCentroid.getY();
        
        double minDistSq = Double.MAX_VALUE;
        int minDistIdx = -1;
        
        // TODO: need normalization for 
        // diffX * diffX + diffY * diffY.
        // slic super pixels uses the super pixel area sXs.
        // will use average as area/k
        float norm = img.getNPixels()/nSeeds;
        
        for (int i = 0; i < seeds.length; ++i) {            
            float diffH = seeds[i][0] - hsv[0];
            float diffS = seeds[i][1] - hsv[1];
            float diffV = seeds[i][2] - hsv[2];
            float diffX = seeds[i][3] - x;
            float diffY = seeds[i][4] - y;
            
            double distSq = 
                diffH * diffH + diffS * diffS + diffV * diffV +
                ((diffX * diffX + diffY * diffY)/(norm*norm));
            
            if (distSq < minDistSq) {
                minDistSq = distSq;
                minDistIdx = i;
            }
        }
        
        return minDistIdx;
    }

    /**
     [nSeeds][h, s, v, x, y]
     int[][] seeds;
    
     * @param img
     * @param starterSeeds
     * @param starterSeedColors
     * @return 
     */
    private float[][] createStartSeeds(ImageExt img, List<PairInt> starterSeeds, 
        List<GroupPixelRGB0> starterSeedColors) {
        
        int k = starterSeeds.size();
        
        float[] hsv = new float[3];
        
        float[][] s = new float[k][];
        for (int i = 0; i < k; ++i) {
            GroupPixelRGB0 rgb = starterSeedColors.get(i);
            PairInt xy = starterSeeds.get(i);
            
            Color.RGBtoHSB((int)Math.round(rgb.getAvgRed()), 
                (int)Math.round(rgb.getAvgGreen()), 
                (int)Math.round(rgb.getAvgBlue()), hsv);
            
            s[i] = new float[5];
            s[i][0] = hsv[0];
            s[i][1] = hsv[1];
            s[i][2] = hsv[2];
            s[i][3] = xy.getX();
            s[i][4] = xy.getY();
        }
        
        return s;
    }

}
