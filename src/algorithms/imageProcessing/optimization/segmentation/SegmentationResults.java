package algorithms.imageProcessing.optimization.segmentation;

import algorithms.compGeometry.PerimeterFinder;
import algorithms.imageProcessing.FixedSizeSortedIntVector;
import algorithms.imageProcessing.MiscellaneousCurveHelper;
import algorithms.imageProcessing.util.MatrixUtil;
import algorithms.misc.Misc;
import algorithms.misc.MiscMath;
import algorithms.util.PairInt;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import thirdparty.HungarianAlgorithm;
import thirdparty.ods.Integerizer;
import thirdparty.ods.XFastTrie;
import thirdparty.ods.XFastTrieNode;

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
        
        xCentroids = new int[n];
        yCentroids = new int[n];
        nPoints = new int[n];
        perimeters = new ArrayList<Set<PairInt>>();
        
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        
        int yMax = Integer.MIN_VALUE;
        int xMax = Integer.MIN_VALUE;
        
        for (int i = 0; i < segmentedSets.size(); ++i) {
            
            Set<PairInt> set = segmentedSets.get(i);
            
            double[] xyCen = curveHelper.calculateXYCentroids(set);
            
            xCentroids[i] = (int)Math.round(xyCen[0]);
            
            yCentroids[i] = (int)Math.round(xyCen[1]);
            
            nPoints[i] = set.size();
        
            Set<PairInt> border = extractBorder(set);
            
            perimeters.add(border);
            
            //int[]{xMin, xMax, yMin, yMax}
            int[] xMinMaxYMinMax = MiscMath.findMinMaxXY(border);
            
            if (xMinMaxYMinMax[1] < xMax) {
                xMax = xMinMaxYMinMax[1];
            }
            if (xMinMaxYMinMax[3] < yMax) {
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
    public double evaluate(SegmentationResults expected) {
        
        if (true) {
            throw new UnsupportedOperationException(
                "chaging the cost function to the precisiona nd recall");
        }
        
        /*
        for each point in perimeters, find the 6 nearest
        neighbors in expected.perimeters.

        discard those larger than dMax=2
        
        after have the set of candidates,
            discard those with no candidate neighbors
        
        int[] dist  x0 -> y0
        */
        
        int nPerimeterPoints = sumNPerimeters();
        
        int dMax = 2;
        int dMaxSq = dMax * dMax;
           
        /*
        need a k nearest neighbors search to make the cost
        matrix input for the min cost bipartite matching.
        
        choices are:
        (1) for each point in a boundary in perimeters,
           iteratate over a region x += dMax and y +- dMax
           and test membership in expected.perimeters set.
           can use a fixedsortedintvector to keep the top k
           (smallest distances).
           For dMax = 2, the scan over neighbors is O(24)
           and each insert into sorted vector is O(lg_2(k)).
           so the maximum runtime complexity is O(24 * lg_2(k)).
           for k=3 --> O(38)
        (2) use XFastTries for predecessor and successor
            searches over x and y separately, then
            look for which results within dMax are a valid
            (x,y) pair in expected.perimeters.
            The runtime complexity is dependent upon the
            size of the image and upon k.
            For 2048, w, the max number of bits needed is
            12 (signed numbers...) for example.
            maximum runtime complexity is
              O(2 * k * lg_2(w)) + O(k*k)
            so for k=3 and large image of 2048 --> O(22)
            and for an image closer to 400 x 600 --> O(20) so
            does not change much with increasing max dimension.
        
        Note that both of the above have a factor of
        N_perimeter_points not included in notes above.
        */
        
        Set<PairInt> allExpectedPoints = getAllPoints(expected);
          
        int n1 = nPerimeterPoints;
        int n2 = allExpectedPoints.size();
        float[][] cost = new float[n1][n2];        
        
        XFastTrie<XFastTrieNode<Integer>, Integer> xbt
            = loadWithXPoints();
        
        XFastTrie<XFastTrieNode<Integer>, Integer> ybt
            = loadWithYPoints();
        
        // r.t. complexity is m*m + O(2*m*lg_2(w)), 
        // so approx O(21) for maxDimension = 2048
        // or 
        int m = 3;
        int[] xIdxs = new int[m];
        int[] yIdxs = new int[m];
        
        for (Set<PairInt> perimeter : expected.perimeters) {
        
            for (PairInt p : perimeter) {
                int x = p.getX();
                int y = p.getY();
                
                int nX = findClosest(x, xbt, dMax, xIdxs);
                int nY = findClosest(y, ybt, dMax, yIdxs);
                
                for (int i = 0; i < nX; ++i) {
                    int x2 = xIdxs[i];
                    for (int j = 0; j < nY; ++j) {
                        int y2 = yIdxs[i];
                        PairInt p2 = new PairInt(x2, y2);
                        if (allExpectedPoints.contains(p2)) {
        
                            int diffX = x2 - x;
                            int diffY = y2 - y;
                            int distSq = (diffX * diffX) + (diffY * diffY);
                            if (distSq > dMaxSq) {
                                continue;
                            }
                            
                            /*
                            store in cost matrix... need indexed points
                            */
                        }
                    }
                }
            }
        }
        
        //TODO: return fMeasure
        return 1;
    }
    
    /**
     * calculate the difference between other and this instance.
     * edits are in progress.
     * 
     * @param expected
     * @return 
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
        
        PerimeterFinder perimeterFinder = new PerimeterFinder();
        
        Set<PairInt> outputEmbeddedGapPoints = new HashSet<PairInt>();

        //int[]{xMin, xMax, yMin, yMax}
        int[] xMinMaxYMinMax = MiscMath.findMinMaxXY(set);
        int imageMaxColumn = xMinMaxYMinMax[1] - 1;
        int imageMaxRow = xMinMaxYMinMax[3] - 1;
       
        int[] rowMinMax = new int[2];
                
        Map<Integer, List<PairInt>> rowColRanges = 
            perimeterFinder.find(set, rowMinMax, imageMaxColumn, 
            outputEmbeddedGapPoints);

        if (!outputEmbeddedGapPoints.isEmpty()) {
            // update the perimeter for "filling in" embedded points
            perimeterFinder.updateRowColRangesForAddedPoints(rowColRanges, 
                rowMinMax, imageMaxColumn, outputEmbeddedGapPoints);
        }
       
        Set<PairInt> borderPixels = perimeterFinder.getBorderPixels(
            rowColRanges, rowMinMax, imageMaxColumn, imageMaxRow);
        
        /*
        ZhangSuenLineThinner lt = new ZhangSuenLineThinner();
        lt.applyLineThinner(borderPixels, 0, imageMaxColumn, 0, imageMaxRow);
        PostLineThinnerCorrections pltc = new PostLineThinnerCorrections();
        pltc.correctForExtCorner(borderPixels, imageMaxColumn + 1, imageMaxRow + 1);
        */ 
        
        return borderPixels;
    }

    private int sumNPerimeters() {

        int n = 0;
        
        for (Set<PairInt> perimeter : perimeters) {
            n += perimeter.size();
        }
        
        return n;
    }

    private XFastTrie<XFastTrieNode<Integer>, Integer> loadWithXPoints() {
        
        int xW = 1 + (int)(Math.floor(Math.log(maxX)/Math.log(2)));
        
        Integerizer<Integer> it = new Integerizer<Integer>() {
            @Override
            public int intValue(Integer x) {
                return x;
            }
        };
        
        XFastTrieNode<Integer> node = new XFastTrieNode<Integer>();
        
        XFastTrie<XFastTrieNode<Integer>, Integer> xbt = 
            new XFastTrie<XFastTrieNode<Integer>, Integer>(
            node, it, xW);
        
        for (Set<PairInt> set : perimeters) {
            for (PairInt p : set) {
                int x = p.getX();
                xbt.add(Integer.valueOf(x));
            }
        }
        
        return xbt;
    }
    
    private XFastTrie<XFastTrieNode<Integer>, Integer> 
        loadWithYPoints() {
        
        int yW = 1 + (int)(Math.floor(Math.log(maxY)/Math.log(2)));
        
        Integerizer<Integer> it = new Integerizer<Integer>() {
            @Override
            public int intValue(Integer x) {
                return x;
            }
        };
        
        XFastTrieNode<Integer> node = new XFastTrieNode<Integer>();
        
        XFastTrie<XFastTrieNode<Integer>, Integer> ybt = 
            new XFastTrie<XFastTrieNode<Integer>, Integer>(
            node, it, yW);
        
        for (Set<PairInt> set : perimeters) {
            for (PairInt p : set) {
                int y = p.getY();
                ybt.add(Integer.valueOf(y));
            }
        }
        
        return ybt;
    }

    private int findClosest(int value, 
        XFastTrie<XFastTrieNode<Integer>, Integer> bt, 
        int dMax, int[] output) {
        
        Integer vIndex = Integer.valueOf(value);
       
        int k = output.length;
        
        int n = 0;
        Integer v0 = bt.find(vIndex);
        if (v0 != null) {
            output[n] = v0.intValue();
            n++;
        }
        
        v0 = vIndex;
        for (int i = 0; i < k/2; ++k) {
            v0 = bt.predecessor(v0);
            if (v0 == null || v0.intValue() > dMax) {
                break;
            } 
            output[n] = v0.intValue();
            n++;
        }
        v0 = vIndex;
        for (int i = 0; i < (k - n); ++k) {
            v0 = bt.successor(v0);
            if (v0 == null || v0.intValue() > dMax) {
                break;
            } 
            output[n] = v0.intValue();
            n++;
        }
        return n;
    }

    private Set<PairInt> getAllPoints(SegmentationResults expected) {

        Set<PairInt> output = new HashSet<PairInt>();
        
        for (Set<PairInt> set : expected.perimeters) {
            output.addAll(set);
        }
        
        return output;
    }
    
}
