package algorithms.imageProcessing.optimization.segmentation;

import algorithms.imageProcessing.MiscellaneousCurveHelper;
import algorithms.imageProcessing.util.MatrixUtil;
import algorithms.util.PairInt;
import java.util.Arrays;
import java.util.List;
import java.util.Set;
import thirdparty.HungarianAlgorithm;

/**
 *
 * @author nichole
 */
public class SegmentationResults {
    
    private final int[] xCentroids;
    private final int[] yCentroids;
    private final int[] nPoints;
    
    public SegmentationResults(List<Set<PairInt>> segmentedSets) {
        
        int n = segmentedSets.size();
        
        xCentroids = new int[n];
        yCentroids = new int[n];
        nPoints = new int[n];
        
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        
        for (int i = 0; i < segmentedSets.size(); ++i) {
            
            Set<PairInt> set = segmentedSets.get(i);
            
            double[] xyCen = curveHelper.calculateXYCentroids(set);
            
            xCentroids[i] = (int)Math.round(xyCen[0]);
            
            yCentroids[i] = (int)Math.round(xyCen[1]);
            
            nPoints[i] = set.size();
        }
    }
    
    public SegmentationResults(int[][] xyCenN) {
        
        if (xyCenN.length != 3) {
            throw new IllegalArgumentException("expecting xyCenN.length == 3");
        }
        
        int n = xyCenN[0].length;
        
        xCentroids = new int[n];
        yCentroids = new int[n];
        nPoints = new int[n];
                
        for (int i = 0; i < n; ++i) {
            
            xCentroids[i] = xyCenN[0][i];
            
            yCentroids[i] = xyCenN[1][i];
            
            nPoints[i] = xyCenN[2][i];
        }
    }
    
    public int getNumberOfItems() {
        return xCentroids.length;
    }
    
    /**
     * calculate the difference between other and this instance
     * using an optimal O(N^3) bipartite matching and a penalty for haivng
     * different sizes of solutions.  Note that if the number of items
     * in the solution is large, this method needs to be altered to 
     * determine when to use a greedy O(N^2) method.
     * @param expected
     * @return 
     */
    public double calculateDifference(SegmentationResults expected) {
        
        /*
        TODO: the cost function needs to use the dimensions of the cells
           or define a difference between perimeters for optimization.
        */
        
        boolean doNormalize = false;
        
        /*
        find the closest match in other for each item and sum the
        differences.
        then add a term that is a penalty sum for not having the
        right number of items.
        
        The segmentation is expected to result in not very many items
        so am chosing to use bipartite matching, which has at best a 
        runtime complexity of O(N^3) though there are methods that 
        use O(m*n*lg2(n)).
        
        NOTE: if this method needs to be adjusted for large n in the future,
        one could make greedy O(N^2) approach and start with the 
        maximum of nPoints.
        */
        
        //TODO: consider adding a term for the number of points in each
        // item
        
        int n1 = this.xCentroids.length;
        int n2 = expected.xCentroids.length;
        
        if (n1 == 0 || n2 == 0) {
            return Double.MAX_VALUE;
        }
        double nSum = 0;
        double[] weights = new double[n2];
        for (int i = 0; i < n2; ++i) {
            weights[i] += expected.nPoints[i];
            nSum += weights[i];
        }
        /*
        // weights that normalize by number of points inversely
        double tot = 0;
        for (int i = 0; i < n2; ++i) {
            double div = (nSum - weights[i]) / ((n2 - 1.) * nSum);
            weights[i] = (float) div;
            tot += div;
        }*/
        // weights that iven items with larger number of points a higher proportion of difference
        double tot = 0;
        for (int i = 0; i < n2; ++i) {
            weights[i] /= nSum;
            tot += weights[i];
        }
        
        /*
        2
        3
        5    
        wanting the portion of "5"'s results to be 5/10 of the total differences
            and if "5"'s results are missing in this instance, then 
            it should be counted as it's maximum differnce * npoints["5"]
            where it's maximum difference will be approximated using it's x an y values.
        
        2 d=1   
        3 d=1
        5 d=1
            nSum=10  (2./10.)*1 + (3./10.)*1 + (5./10)*1.
        */
        
        float[][] cost = new float[n1][n2];
        for (int i1 = 0; i1 < n1; ++i1) {
            
            cost[i1] = new float[n2];
            Arrays.fill(cost[i1], Float.MAX_VALUE);
            
            int x1 = xCentroids[i1];
            int y1 = yCentroids[i1];
            
            for (int i2 = 0; i2 < n2; ++i2) {
                
                int x2 = expected.xCentroids[i2];
                int y2 = expected.yCentroids[i2];
                                
                int diffX = x1 - x2;
                int diffY = y1 - y2;
                
                cost[i1][i2] = diffX * diffX + diffY * diffY;
            }
        }
        
        boolean transposed = false;
        if (cost.length > cost[0].length) {
            cost = MatrixUtil.transpose(cost);
            transposed = true;
        }
        
        HungarianAlgorithm b = new HungarianAlgorithm();
        int[][] match = b.computeAssignments(cost);
        
        double diffSum = 0;
        double maxDiff = Double.MIN_VALUE;
        int maxDiffIdx2 = -1;
        
        boolean[] matched = new boolean[n2];
        
        for (int i = 0; i < match.length; i++) {
            int idx1 = match[i][0];
            int idx2 = match[i][1];
            if (idx1 == -1 || idx2 == -1) {
                continue;
            }
            if (transposed) {
                int swap = idx1;
                idx1 = idx2;
                idx2 = swap;
            }            
            int x1 = xCentroids[idx1];
            int y1 = yCentroids[idx1];
            int x2 = expected.xCentroids[idx2];
            int y2 = expected.yCentroids[idx2];
            
            matched[idx2] = true;
                                
            int diffX = x1 - x2;
            int diffY = y1 - y2;
                
            double diff = Math.sqrt(diffX * diffX + diffY * diffY);
            
            diff *= weights[idx2];
            
            diffSum += diff;
            
            if (diff > maxDiff) {
                maxDiff = diff;
            }
        }
        
        // for missing expected items, estimating a cost as the centroid times
        // the number of point
        int nc = 0;
        for (int i = 0; i < n2; ++i) {
            if (!matched[i]) {
                int x2 = expected.xCentroids[i];
                int y2 = expected.yCentroids[i];
                double diff = Math.sqrt(x2 * x2 + y2 * y2) * expected.nPoints[i]/2.;
                diff *= weights[i];
                diffSum += diff;
                nc++;
            }
        }
        
        //TODO: may revise this.  penalty is the difference in number of items,
        // times max difference in matched items.
        double penalty = Math.abs(n1 - (n2 - nc)) * maxDiff;
        
        diffSum += penalty;
        
        return diffSum;
    }
    
}
