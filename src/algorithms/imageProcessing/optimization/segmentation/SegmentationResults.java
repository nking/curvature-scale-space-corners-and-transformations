package algorithms.imageProcessing.optimization.segmentation;

import algorithms.compGeometry.PerimeterFinder;
import algorithms.imageProcessing.MiscellaneousCurveHelper;
import algorithms.imageProcessing.ZhangSuenLineThinner;
import algorithms.imageProcessing.util.MatrixUtil;
import algorithms.misc.Misc;
import algorithms.misc.MiscMath;
import algorithms.search.KNearestNeighbors;
import algorithms.util.PairInt;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
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
    private List<Set<PairInt>> perimeters;
    
    private int maxY = Integer.MIN_VALUE;
    private int maxX = Integer.MIN_VALUE;
    
    public SegmentationResults(List<Set<PairInt>> segmentedSets) {
        
        int n = segmentedSets.size();
        
        int yMax = Integer.MIN_VALUE;
        int xMax = Integer.MIN_VALUE;
        
        {//debug
            for (int i = 0; i < n; ++i) {
                Set<PairInt> set = segmentedSets.get(i);
                int[] xMinMaxYMinMax = MiscMath.findMinMaxXY(set);            
                if (xMinMaxYMinMax[1] > xMax) {
                    xMax = xMinMaxYMinMax[1];
                }
                if (xMinMaxYMinMax[3] > yMax) {
                    yMax = xMinMaxYMinMax[3];
                }
            }
            /*
            if (n > 0) {
                Image img = new Image(xMax + 1, yMax + 1);
                long ts = MiscDebug.getCurrentTimeFormatted();
                MiscDebug.writeAlternatingColor(
                    img, segmentedSets, "seg_" + ts);
            }*/
        }
                
        xCentroids = new int[n];
        yCentroids = new int[n];
        nPoints = new int[n];
        perimeters = new ArrayList<Set<PairInt>>();
        
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        
        for (int i = 0; i < n; ++i) {
            
            Set<PairInt> set = segmentedSets.get(i);
            
            double[] xyCen = curveHelper.calculateXYCentroids(set);
            
            xCentroids[i] = (int)Math.round(xyCen[0]);
            
            yCentroids[i] = (int)Math.round(xyCen[1]);
            
            nPoints[i] = set.size();
        
            Set<PairInt> border = extractBorder(set);
            
            perimeters.add(border);
            
            //int[]{xMin, xMax, yMin, yMax}
            int[] xMinMaxYMinMax = MiscMath.findMinMaxXY(border);
            
            if (xMinMaxYMinMax[1] > xMax) {
                xMax = xMinMaxYMinMax[1];
            }
            if (xMinMaxYMinMax[3] > yMax) {
                yMax = xMinMaxYMinMax[3];
            }
        }
    }
    
    public int getNumberOfItems() {
        return xCentroids.length;
    }
    
    public List<Set<PairInt>> getPointsList() {
        return this.perimeters;
    }
    
    /**
     * see notes in BenchmarkMeasurer.java
     * 
     * @param expected
     * @return 
     */
    public double evaluate(SegmentationResults expected,
        int dMax) {
            
        if (perimeters.size() == 0) {
            return 0;
        }
        
        BenchmarkMeasurer measurer = new BenchmarkMeasurer();
        
        float fMeasure = measurer.evaluate(this, expected, dMax);
        
        return fMeasure;
    }
    
    /**
     * calculate the difference between other and this instance.
     * edits are in progress.
     * 
     * @param expected
     * @return 
     * @deprecated 
     */
    public double calculateDifference(SegmentationResults expected) {
        
        if (true) {
            throw new UnsupportedOperationException(
                "chaging the cost function to the precisiona nd recall");
        }
        
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
        // weights that have items with larger number 
        // of points a higher proportion of difference
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
        /*int nc = 0;
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
        */
        
        /*
        now calculate the difference between points.
        creating a penalty for missing points and additional
        points as the sum of the normalized distance from
        expected centroid.
        */
        double penalty2 = 0;
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
            int x2 = expected.xCentroids[idx2];
            int y2 = expected.yCentroids[idx2];
            
            Set<PairInt> set2 = expected.getPointsList().get(idx2);
            Set<PairInt> set1 = expected.getPointsList().get(idx1);
                
            // xMin, xMax, yMin, yMax
            int[] minMaxXY = MiscMath.findMinMaxXY(set2);
    
            int wX = minMaxXY[1] - minMaxXY[0];
            int wY = minMaxXY[3] - minMaxXY[2];
            double normalization = Math.sqrt(wX * wX + wY * wY);
            
            double sumDist = 0;
            for (PairInt p : set2) {
                if (set1.contains(p)) {
                    continue;
                }
                int diffX = p.getX() - x2;
                int diffY = p.getY() - y2;
                sumDist += Math.sqrt(diffX * diffX + diffY * diffY);
            }
            for (PairInt p : set1) {
                if (set2.contains(p)) {
                    continue;
                }
                int diffX = p.getX() - x2;
                int diffY = p.getY() - y2;
                sumDist += Math.sqrt(diffX * diffX + diffY * diffY);
            }
            sumDist /= normalization;
            penalty2 += sumDist;
        }
        
        diffSum += penalty2;
                
        return diffSum;
    }

    private Set<PairInt> extractBorder(Set<PairInt> set) {

        // until PerimeterFinder is fixed,
        // will use a simple search for points that have
        // no neighbors, and then apply a line thinner
        // to that.
        //    may want to consider a use of the blob medial
        //    axes to determine "inward" and hence fill in
        //    embedded holes and gaps.
        
        Set<PairInt> border = new HashSet<PairInt>();
        
        int xMin = Integer.MAX_VALUE;
        int xMax = Integer.MIN_VALUE;
        int yMin = Integer.MAX_VALUE;
        int yMax = Integer.MIN_VALUE;
        
        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;
        
        for (PairInt p : set) {
            int x = p.getX();
            int y = p.getY();
            for (int i = 0; i < dxs.length; ++i) {
                int x2 = x + dxs[i];
                int y2 = y + dys[i];
                PairInt p2 = new PairInt(x2, y2);
                if (!set.contains(p2)) {
                    border.add(p);
                    break;
                }
            }
            if (x < xMin) {
                xMin = x;
            }
            if (y < yMin) {
                yMin = y;
            }
            if (x > xMax) {
                xMax = x;
            }
            if (y > yMax) {
                yMax = y;
            }
        }
        
        xMax++;
        yMax++;
        yMin--;
        xMin--;
        if (xMin < 0) {
            xMin = 0;
        }
        if (yMin < 0) {
            yMin = 0;
        }
        
        ZhangSuenLineThinner lt = new ZhangSuenLineThinner();
        lt.applyLineThinner(border, xMin, xMax, yMin, yMax);
      
        return border;
    }

    public int sumNPerimeters() {

        int n = 0;
        
        for (Set<PairInt> perimeter : perimeters) {
            n += perimeter.size();
        }
        
        return n;
    }    

    public Set<PairInt> getAllPoints() {

        Set<PairInt> output = new HashSet<PairInt>();
        
        for (Set<PairInt> set : perimeters) {
            output.addAll(set);
        }
        
        return output;
    }
    
    public List<Set<PairInt>> getPerimeters() {
        return perimeters;
    }
    
    public KNearestNeighbors createKNN() {
        
        int n = sumNPerimeters();
        
        if (n == 0) {
            throw new IllegalStateException("perimeters "
                + "cannot be empty");
        }
        
        float[] x = new float[n];
        float[] y = new float[n];
        int count = 0;
        for (Set<PairInt> perimeter : perimeters) {
            for (PairInt p : perimeter) {
                x[count] = p.getX();
                y[count] = p.getY();
                count++;
            }
        }
        KNearestNeighbors kNN = new KNearestNeighbors(x, y);
        return kNN;
    }

}
