package algorithms.imageProcessing.matching;

import algorithms.QuickSort;
import algorithms.compGeometry.PerimeterFinder2;
import algorithms.imageProcessing.ColorHistogram;
import algorithms.imageProcessing.MiscellaneousCurveHelper;
import algorithms.util.PairIntArray;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.set.TIntSet;
import java.util.List;
import algorithms.imageProcessing.matching.PartialShapeMatcher.Result;
import algorithms.imageProcessing.transform.EuclideanEvaluator;
import algorithms.imageProcessing.transform.EuclideanTransformationFit;
import algorithms.imageProcessing.transform.TransformationParameters;
import algorithms.util.PairInt;
import algorithms.util.TwoDIntArray;
import algorithms.util.VeryLongBitString;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TObjectDoubleMap;
import gnu.trove.map.hash.TObjectDoubleHashMap;
import gnu.trove.set.hash.TIntHashSet;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import thirdparty.edu.princeton.cs.algs4.Interval;
import thirdparty.edu.princeton.cs.algs4.Interval2D;
import thirdparty.edu.princeton.cs.algs4.QuadTree;

/**
 * uses PartialShapeMatcher and search patterns to
 * find the best fitting match between a
 * group of adjacent segmented cells
 * to a template shape.  Any additional information
 * that helps limit the search should be done before
 * this stage and included or excluded in the
 * adjacency map.  Note that this class could change
 * in the future to use weights in the adjacency map,
 * but for now it's a binary relationship map.
 *
 * Also note that the class uses internal cache which requires single
 * threaded use of this class.
 * 
 * @author nichole
 */
public class ShapeFinder {
    
    /*
    TODO: 
       will add a method to include the orb keypoints and oriented half 
          descriptors.
       the global 2D bins for starting searches will be present for this new
         method.
       the descriptors will be used before partial shape matching
          to evaluate whether adding the cell reduces the cost of the
          path to the cell.
       will also change the floyd warshall aggregation/search to a dijkstra's
          search.
          -- both the descriptors and partial shape matching 
          might be necessary to determine the best aggregation at each step.
    */
    
    private double maxDiffChordSum = Float.MAX_VALUE;
    
    private double maxDistTransformSum = Float.MAX_VALUE;
    
    private boolean useEuclideanInCost = false;
    private float[] euclideanScaleRange = new float[]{0.9f, 1.1f};

    private final float areaFactor = 2.f;
    private final TIntSet outOfSizeRange = new TIntHashSet();

    // partial shape matcher sampling distance
    private final int dp = 1;
    private final boolean useSameNSampl = true;
    private final int innerTolerance = 3;
    
    private final Map<VeryLongBitString, PairIntArray> aggregatedBoundaries =
        new HashMap<VeryLongBitString, PairIntArray>();
    private final Map<VeryLongBitString, Result> aggregatedResultMap =
        new HashMap<VeryLongBitString, Result>();
    private final Map<VeryLongBitString, Double> aggregatedCostMap =
        new HashMap<VeryLongBitString, Double>();
    
    private int[] minMaxXY = null;
    
    private final int[] dimensionsT;
    private int areaT;
    private final int nB;
    
    private final List<PairIntArray> orderedBoundaries;
    private final List<Set<PairInt>> pointsList;
    private final TIntObjectMap<TIntSet> adjacencyMap; 
    private final PairIntArray template;
    private final int[][] templateHSVHist;
    private final List<TwoDIntArray> imgHSVHist;
    
    private QuadTree<Integer, Integer> centroidQT = null;
    
    private final ColorHistogram ch = new ColorHistogram();
    private final int[][] cache1;
    private final float chHSVLimit1 = 0.1f;
    private final float chHSVLimit2 = 0.25f;
    private final float chHSVLimit3 = 0.35f;
    private final Set<VeryLongBitString> skipSet = new HashSet<VeryLongBitString>();
        
    private final int nTop = 100;

    //DEBUG
    TIntSet expectedIndexes = new TIntHashSet();
    VeryLongBitString tBS = null;
    int tBSIdx = -1;
    
    /**
     * NOT READY FOR USE.
     */
    public ShapeFinder(List<PairIntArray> orderedBoundaries,
        List<Set<PairInt>> pointsList, TIntObjectMap<TIntSet> adjacencyMap, 
        PairIntArray template, int[][] templateHSVHist, 
        List<TwoDIntArray> imgHSVHist) {
        
        this.orderedBoundaries = orderedBoundaries;
        this.pointsList = pointsList;
        this.adjacencyMap = adjacencyMap;
        this.template = template;
        this.templateHSVHist = templateHSVHist;
        this.imgHSVHist = imgHSVHist;
        
        cache1 = new int[templateHSVHist.length][templateHSVHist[0].length];
        for (int i = 0; i < cache1.length; ++i) {
            cache1[i] = new int[templateHSVHist[0].length];
        }
        
        assert(pointsList.size() == orderedBoundaries.size());
        
        this.dimensionsT = calcDimensions(template);
        this.areaT = (int)Math.round(Math.sqrt(dimensionsT[0] * dimensionsT[0] +
            dimensionsT[1] * dimensionsT[1]));
        this.nB = orderedBoundaries.size();
    }
    
    /**
     * NOT READY FOR USE.
     * uses PartialShapeMatcher and search patterns to
     * find the best fitting match between a
     * group of adjacent segmented cells
     * to a template shape.  Any additional information
     * that helps limit the search should be done before
     * this stage and included or excluded in the
     * adjacency map.  Note that this class could change
     * in the future to use weights in the adjacency map,
     * but for now it's a binary relationship map.
     * NOTE also that methods to return the top k results
     * or the equivalent best within a tolerance might be
     * made in the future.
     *
     * The runtime is dependent upon which of the 3 search
     * patterns is kept in the end, so that will be
     * filled in here after implementation and testing.
     *
     * @return
     */
    public Result[] findMatchingCells() {

        /*
        some notes that will be updated when return to the search
        task.  higher priority before that is improvement of input so that true match
        is top result consistently.
        
        step (1) O(N_cells) + O(N_cells_in_bin * log_2(N_cells_in_bin)):
           the bin dimensions will be max of template dimensions.
           actually, should probably be the diagonal for max dist
           between any two points but that could be done more
           accurately.
           - 2D bins, that is x and y.
           - sort the cells by x centroid and ties by y centroid
        step (2) search pattern, 3 are presented.
            The bigger picture is a sliding window of the 2D bin
            where the core of the search and aggregation is done
            within the 2D bin, but the aggregation can continue
            into the next bin to the right, next bin above,
            or above and to the right.
            The search of the cells could proceed in 3 different
            ways.  will implement all 3 and keep the fastest of
            the robust of them.
            (2a)
               try all combinations in a bin, then all combinations
               that can be aggregated on to them that are in
               the adjacent "next" bins.
               this method is thorough brute force, but is
               combinatorial within bounds of area 3*binsize (or smaller
               if adj map is sparse) in runtime,
               so a faster but still accurate method would be better.
            (2b and c)
               over the whole image, use the partial shape matcher
               with cell boundaries that have area less than
               two or so times the dimension of the template.
               store the top k of those matches rather than just
               the best matching presumably.  this is the first
               step of 2 different searches, 2b and 2c.
               (2b)
                   visit each bin determining whether a match
                   to the whole template is findable within the
                   top k search results of each, aggregating
                   adjacent where it is matching or where it is in between
                   such matches (that is, embedded) ...
                   this is solving the problem on a cell level and combining
                   results to make a better solution.
                   difficult to see this as robust, but
                   if significant fraction of matching cells were
                   within a small top k of 5 or less of individul cell shape matching,
                   this might be feasible and it would be polynomial
                   search, bounded by the segmentation fragmentation
                   template size, and the choice of k.
               (2c)
                   for each bin, choose as the start of a search
                   the best fitting cell, then proceed with a path
                   like aggregation search (dijkstra's, for example).
                   more robust would be to let each cell with
                   a reasonable matching cost, be a start seed in
                   a separate path search.
                   if the later is done, need to keep track of
                   already searched groups across path solutions
                   and reuse those result.
                   Dijkstra's runtime for a single path search is O(E lg V).
                   **NOTE that this pattern of a dijkstra's search for each
                   feasible start could be done over the whole image's cells
                   instead of a bin's cells and could presumably be done
                   over a smaller number of starts than all cells by sorting
                   the start list by the partial match cost (score)
                   and reuse of the aggregated cell partial results.

               there are then 4 search patterns to implement here.
        */
        
        // ----- init state ----
        List<PairInt> orderedBoundaryCentroids = calculateCentroids();
        assert(orderedBoundaryCentroids.size() == nB);
        
        centroidQT = new QuadTree<Integer, Integer>();
        for (int i = 0; i < nB; ++i) {
            PairInt pCen = orderedBoundaryCentroids.get(i);
            centroidQT.insert(pCen.getX(), pCen.getY(), Integer.valueOf(i));
        }
        
        this.minMaxXY = calculateMinMaxXY();
        
        assert(minMaxXY[0] >= 0);
        assert(minMaxXY[1] >= 0);
        assert(minMaxXY[2] >= 0);
        assert(minMaxXY[3] >= 0);

        matchIndividually();
        
        // ----------
        
        { // DEBUG
            //DEBUG specific to a test
            Set<PairInt> allPoints = new HashSet<PairInt>();
            expectedIndexes.clear();
            tBS = new VeryLongBitString(orderedBoundaries.size());
            PairInt[] cens = new PairInt[8];
            cens[0] = new PairInt(184,86);
            cens[1] = new PairInt(187,79);
            cens[2] = new PairInt(184,60);
            cens[3] = new PairInt(177,41);
            cens[4] = new PairInt(174,76);
            cens[5] = new PairInt(173,52);
            cens[6] = new PairInt(189,45);
            cens[7] = new PairInt(184,35);

            for (int i = 0; i < orderedBoundaries.size(); ++i) {
                PairInt pCen = orderedBoundaryCentroids.get(i);
                System.out.println("pCen=" + pCen + " idx=" + i);
                for (int j = 0; j < cens.length; ++j) {
                   int dx = Math.abs(pCen.getX() - cens[j].getX());
                   int dy = Math.abs(pCen.getY() - cens[j].getY());
                   if ((dx < 3 && dy < 3) /*|| 
                       (pCen.getX() > 170 && 190 < pCen.getX() &&
                       pCen.getY() > 30 && 90 < pCen.getY())*/) {
                       tBS.setBit(i);
                       expectedIndexes.add(i);
                       if (expectedIndexes.size() == 4) {
                           tBSIdx = i;
                       }
                       allPoints.addAll(pointsList.get(i));
                       System.out.println("*pCen=" + pCen + " idx=" + i + " "
                          + " adj=" + adjacencyMap.get(i));
                       break;
                   }
                }
            }
            System.out.println("found " + expectedIndexes.size() + " of " 
                + cens.length + " expected.  tBS=" + Arrays.toString(tBS.getSetBits()));
               /* 
                try {
                   PairIntArray b =  mergeAdjacentOrderedBorders(tBS);
                   int[] xPolygon = null;
                   int[] yPolygon = null;
                   PolygonAndPointPlotter plotter = new PolygonAndPointPlotter();
                   int[] xminmxyminmiac = MiscMath.findMinMaxXY(allPoints);
                   int[] xp, yp;
                   int n2 = allPoints.size();
                   xp = new int[n2];
                   yp = new int[n2];
                   int i = 0;
                   for (PairInt p : allPoints) {
                       xp[i] = p.getX();
                       yp[i] = p.getY();
                       i++;
                   }
                   plotter.addPlot(xminmxyminmiac[0], xminmxyminmiac[1],
                       xminmxyminmiac[2], xminmxyminmiac[3],
                       xp, yp, xPolygon, yPolygon, "shape tbs");
                   plotter.writeFile2();
               } catch (Throwable t) { }        
               int z = 1;
               */
        }// end DEBUG
        
        
        if (true) {
            return wideFWSearch(orderedBoundaryCentroids);
        }

        throw new UnsupportedOperationException("not yet implemented");
    }

    private Result[] wideFWSearch(List<PairInt> orderedBoundaryCentroids) {    
        
        //TODO: looks like need to retain the best of each bin
        //      and further compare those to one another.
        //      -- might need to consider a fast inner pattern comparison.
        //      -- might need to consider a quick way to trim the aggregated
        //         shapes by the found correspondence because that correspondence
        //         step of partial shape matching excludes occluded sections and
        //         external non-object inclusions.
        //         currently, that information of the subset of the aggregated
        //         shape that does match is not being used after the partial 
        //         shape match,
        //         instead the entire aggregated shape is being used in further
        //         evalutations such as the color histograms.
        //         this change then implies that using half descriptors 
        //         (meaning, only the half of the descriptor facing inward of
        //         shape) of keypoints would be up for consideration again.
        //         - such descriptors probably need to be color in the global search
        //           stage (as a replacement or supplement to color histograms)
        //           and would only include the keypoints matched by correspondence
        //           (w/ additional logic to keep the comparison to intersection
        //           of shapes...excluding the occlusion and additions of non-object
        //           to the aggregated shape)
        //         - such descriptors, if used at the stage of euclidean evaluation
        //           in partial shape matching, would need to be few in number 
        //           to keep the partial shape matching fast...)
        //TODO: consider if another pattern of color evaluation would
        //      be quicker than the color histograms.  they're fast because
        //      the sizes are only 32 bins and comparisons are the same. 
        //      but might be able to
        //      use average colors in quadrants...oriented averages,
        //      but with small number of steps to combine color averages.
        //      (might require standard deviations).  a step towards descriptors
        //      in terms of using spatial location with respect to same quadrant
        //      in template object, but fewer comparisons...
        //TODO: once the search over the entire image is robust, will make 
        //      another method attempting a faster search pattern than the
        //      nBins * O(N^3) Floyd Warshall (excluding partial shape matching
        //      and color histograms in the complexity summary just now).
        //TODO: consider making an alternative method for the outer search
        //      pattern, that is, the global search pattern of bins across the 
        //      image.  could instead use a kmeans alteration of the bin centers
        //      and sizes similarly as super-pixels does, and then only search 
        //      the top m of those ... 
        //      the means is the intersection of color histograms...
        
        int maxDim = Math.max(dimensionsT[0], dimensionsT[1]);
        
        int w = minMaxXY[1] + 1;
        int h = minMaxXY[3] + 1;
        
        int nPPB = 1;
        
        int nX = nPPB * ((int)Math.floor(w/maxDim) + 1);
        int nY = nPPB * ((int)Math.floor(h/maxDim) + 1);
        
        int binCount = 0;
        
        double minCost = Double.MAX_VALUE;
        VeryLongBitString minBitString = null;
        
        long tShapeSum = 0;
         
        for (int j = 0; j < nY; ++j) {
            int startY = (maxDim/nPPB) * j;
            if (startY > (h - 1)) {
                continue;
            }
            int stopY = startY + maxDim;
            if (stopY >= h) {
                stopY = h - 1;
            }
            Interval<Integer> intY = new Interval<Integer>(startY, stopY);
            for (int i = 0; i < nX; ++i) {
                
                long t0 = System.currentTimeMillis();
                
                int startX = (maxDim/nPPB) * i;
                if (startX > (w - 1)) {
                    continue;
                }
                int stopX = startX + maxDim;
                if (stopX >= w) {
                    stopX = w - 1;
                }
                Interval<Integer> intX = new Interval<Integer>(startX, stopX);
                Interval2D<Integer> rect = new Interval2D<Integer>(intX, intY);
               
                List<Integer> indexes = centroidQT.query2D(rect);
                
                System.out.println(String.format(
                    "shape finder bin (%d:%d, %d:%d)  nSegments=", startX, stopX, 
                    startY, stopY, indexes.size()));
                
                if (indexes.isEmpty()) {
                    // no segmented cells in this bin
                    continue;
                }
                
                // chose the smallest x, smallest y in bin
                int index = findSmallestXYCentroid(indexes, orderedBoundaryCentroids);
           
                if (index < 0) {
                    continue;
                }
                
                binCount++;
                
                PairInt pCen = orderedBoundaryCentroids.get(index);
                System.out.println("pCen=" + pCen + " i=" + index);

                VeryLongBitString bitString = minCostAggregationFW(index,
                    orderedBoundaryCentroids);

                if (bitString == null) {
                    continue;
                }

                double cost = aggregatedCostMap.get(bitString).doubleValue();
                if (cost < minCost) {
                    minCost = cost;
                    minBitString = bitString;
                }
               
                long t1 = System.currentTimeMillis();
                long t1Sec = (t1 - t0) / 1000;
                System.out.println(t1Sec + " sec to match " + template.getN()
                    + " points to " + aggregatedBoundaries.get(bitString).getN());
                tShapeSum += t1Sec;
            }
        }
        
        System.out.println("time spent in shape matching=" + tShapeSum);
        
        System.out.println("binCount=" + binCount);

        if (minBitString == null) {
            return null;
        }

        long t0 = System.currentTimeMillis();
        
        // ---- choosing nTop results from below to return as results -------

        Result r = aggregatedResultMap.get(minBitString);
        Object[] data = new Object[2];
        data[0] = minBitString;
        data[1] = aggregatedBoundaries.get(minBitString);
        r.setData(data);
       
        // because of the blur performed on the boundaries, sometimes, there
        // are slightly different indexes that can result in the same boundaries.
        // so, find and remove redundant entries.
        // TODO: may be able to remove this
        {   int rmvd = aggregatedCostMap.size();
            Set<VeryLongBitString> key2Set = aggregatedCostMap.keySet();
            VeryLongBitString[] key2 = key2Set.toArray(new 
                VeryLongBitString[aggregatedCostMap.size()]);
            for (int i = 0; i < key2.length; ++i) {
                VeryLongBitString bs1 = key2[i];
                if (!key2Set.contains(bs1)) {
                    continue;
                }
                PairIntArray p1 = aggregatedBoundaries.get(bs1);
                Set<PairInt> set1 = new HashSet<PairInt>();
                for (int j = 0; j < p1.getN(); ++j) {
                    set1.add(new PairInt(p1.getX(j), p1.getY(j)));
                }
                Set<VeryLongBitString> rm = new HashSet<VeryLongBitString>();
                for (VeryLongBitString bs2 : key2Set) {
                    if (bs1.equals(bs2)) {
                        continue;
                    }
                    PairIntArray p2 = aggregatedBoundaries.get(bs2);
                    Set<PairInt> set2 = new HashSet<PairInt>();
                    for (int j = 0; j < p2.getN(); ++j) {
                        set2.add(new PairInt(p2.getX(j), p2.getY(j)));
                    }
                    Set<PairInt> s1Minus2 = new HashSet<PairInt>(set1);
                    s1Minus2.removeAll(set2);
                    Set<PairInt> s2Minus1 = new HashSet<PairInt>(set2);
                    s2Minus1.removeAll(set1);
                    if (s1Minus2.size() == 0 && s2Minus1.size() == 0) {
                        rm.add(bs2);
                    }
                }
                for (VeryLongBitString bs : rm) {
                    key2Set.remove(bs);
                    aggregatedCostMap.remove(bs);
                    aggregatedResultMap.remove(bs);
                    aggregatedBoundaries.remove(bs);
                }
            }
            rmvd -= aggregatedCostMap.size();
            System.out.println("removed " + rmvd + " redundant results"); 
        }
        
        // could use FixedSizeSortedVector to make this faster for only nTop objects
        int nc = aggregatedCostMap.size();
        float[] costs = new float[nc];
        VeryLongBitString[] keys = new VeryLongBitString[nc];
        int count = 0;
        for (Entry<VeryLongBitString, Double> entry : aggregatedCostMap.entrySet()) {
            VeryLongBitString key = entry.getKey();
            if (containsASkipItem(key)) {
                continue;
            }
            int[] bits = key.getSetBits();
            if (bits.length == 1 && outOfSizeRange.contains(bits[0])) {
                continue;
            }
            
            fillCache1(bits);
            float intersection = ch.intersection(cache1, templateHSVHist);
            //System.out.println("-> intersection=" + intersection + " bits=" + 
            //    Arrays.toString(bits));
            if (intersection < chHSVLimit3) {
                skipSet.add(key);
                continue;
            }
            assert(aggregatedResultMap.containsKey(key));
            keys[count] = key;
            costs[count] = entry.getValue().floatValue();
            count++;
        }
        costs = Arrays.copyOf(costs, count);
        keys = Arrays.copyOf(keys, count);
        
        QuickSort.sortBy1stArg(costs, keys);
        
        int nt = (count >= nTop) ? nTop : count;
        
        System.out.println("returning " + nt + " results");
        Result[] results = new Result[nt + 1];
        for (int i = 0; i < nt; ++i) {
            VeryLongBitString key = keys[i];
            r = aggregatedResultMap.get(key);
            PairIntArray b = aggregatedBoundaries.get(key);
            data = new Object[2];
            data[0] = key;
            data[1] = b;
            r.setData(data);
            assert(b != null);
            results[i + 1] = r;
            int[] idxs = key.getSetBits();
            
            fillCache1(idxs);
            float intersection = ch.intersection(cache1, templateHSVHist);
            
            StringBuilder sb = new StringBuilder(", sizes=");
            for (int idx : idxs) {
                int sz = pointsList.get(idx).size();
                sb.append(sz).append(" ");
            }
            String costStr = String.format("cost=%.4f",
                aggregatedCostMap.get(key).floatValue());
            System.out.println(i + " " + costStr + " " + Arrays.toString(idxs) 
                + sb.toString() + " intersection=" + intersection);
        }
        Double cost = aggregatedCostMap.get(tBS);
        if (cost == null) {
            PairIntArray b =  mergeAdjacentOrderedBorders(tBS);
            cost = calcAndStoreMatchCost(b, tBS);
        }
        System.out.println("expected=" + Arrays.toString(tBS.getSetBits()));
        if (tBS.getSetBits().length > 0) {
            results[0] = aggregatedResultMap.get(tBS);
            data = new Object[2];
            data[0] = tBS;
            data[1] = aggregatedBoundaries.get(tBS);
            results[0].setData(data);
        } else {
            results[0] = results[1];
        }
       
        long t1 = System.currentTimeMillis();
        long t1Sec = (t1 - t0)/1000;
        System.out.println("post shape mtching filter time = " + t1Sec);
        
        return results;
    }    

    private void matchIndividually() {

        double maxChord = Double.MIN_VALUE;
        double maxDist = Double.MIN_VALUE;
        
        for (int i = 0; i < nB; ++i) {

            PairIntArray p = orderedBoundaries.get(i);

            if (p.getN() < 6) {
                continue;
            }

            int[] dimensions = calcDimensions(p);

            int area = (int)Math.round(Math.sqrt(dimensions[0] * dimensions[0] +
                dimensions[1] * dimensions[1]));

            if (area > (areaFactor * areaT)) {
                // NOTE: this assumes the segmentation of the object from
                // a large background or foreground is complete separation.
                // If that isn't true, then a search pattern which is
                // more like a skyline partial matching pattern must be
                // used instead.
                // TODO: need a way to recognize that or at least mention
                // it in the documentation. the skyline search pattern will
                // be implemented soon.
                outOfSizeRange.add(i);
                continue;
            }

            if (area < (areaT/areaFactor)) {
                //TODO: may need to revise this limit
                outOfSizeRange.add(i);
                continue;
            }

            PartialShapeMatcher matcher = new PartialShapeMatcher();
            matcher.overrideSamplingDistance(dp);
            if (useSameNSampl) {
                matcher.setToUseSameNumberOfPoints();
            }
            Result r = matcher.match(template, p);
            if (r == null) {
                continue;
            }
            
            double d = r.getChordDiffSum();
            if (d > maxChord) {
                maxChord = d;
            }
            if (useEuclideanInCost) {
                d = calcTransformationDistanceSum(r, template, 
                    orderedBoundaries.get(i), false);
                if (d > maxDist) {
                    maxDist = d;
                }
            }
        }
        
        maxDiffChordSum = maxChord;
        
        if (useEuclideanInCost) {
            maxDistTransformSum = maxDist;
        }
    }

    private int[] calcDimensions(PairIntArray a) {

        int minX = Integer.MAX_VALUE;
        int maxX = Integer.MIN_VALUE;
        int minY = Integer.MAX_VALUE;
        int maxY = Integer.MIN_VALUE;

        for (int i = 0; i < a.getN(); ++i) {
            int x = a.getX(i);
            int y = a.getY(i);
            if (x < minX) {
                minX = x;
            }
            if (x > maxX) {
                maxX = x;
            }
            if (y < minY) {
                minY = y;
            }
            if (y > maxY) {
                maxY = y;
            }
        }

        return new int[]{maxX - minX + 1, maxY - minY + 1};
    }

    private List<PairInt> calculateCentroids() {

        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();

        List<PairInt> centroids = new ArrayList<PairInt>();

        for (int i = 0; i < nB; ++i) {

            double[] xyCen = curveHelper.calculateXYCentroids(
                orderedBoundaries.get(i));

            PairInt xy = new PairInt((int)Math.round(xyCen[0]),
                (int)Math.round(xyCen[1]));

            centroids.add(xy);
        }

        return centroids;
    }
    
    private int[] calculateMinMaxXY() {

        int minX = Integer.MAX_VALUE;
        int maxX = Integer.MIN_VALUE;
        int minY = Integer.MAX_VALUE;
        int maxY = Integer.MIN_VALUE;
                    
        for (int i = 0; i < nB; ++i) {
            
            PairIntArray a = orderedBoundaries.get(i);
            
            for (int j = 0; j < a.getN(); ++j) {
                int x = a.getX(j);
                int y = a.getY(j);
                if (x < minX) {
                    minX = x;
                }
                if (y < minY) {
                    minY = y;
                }
                if (x > maxX) {
                    maxX = x;
                }
                if (y > maxY) {
                    maxY = y;
                }
            }
        }

        return new int[]{minX, maxX, minY, maxY};
    }

    /**
     * a Floyd Warshall search to find min-cost aggregation of segmented
     * cells within a limited distance of adjacency.
     *
     * @param index
     * @param orderedBoundaryCentroids
     * @return
     */
    private VeryLongBitString minCostAggregationFW(int index, List<PairInt> centroids) {

        assert(centroids.size() == orderedBoundaries.size());
        
        TObjectDoubleMap<PairInt> distMap = new TObjectDoubleHashMap<PairInt>();
        Map<PairInt, VeryLongBitString> indexesMap = new HashMap<PairInt,
            VeryLongBitString>();

        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();

        TIntList indexes = new TIntArrayList();

        PairInt indexXY = centroids.get(index);

        int nT1 = template.getN();

        float factor = 2.f;
        float maxDist = factor * (float)Math.max(dimensionsT[0], dimensionsT[1]);
        
        // ------- initialize the local maps and store partial results in
        //            output vars too

        // insert items within distance of index centroid and areaFactor times
        // the dimensionsT
        for (int i = 0; i < nB; ++i) {

            if (orderedBoundaries.get(i).getN() < 6) {
                continue;
            }

            PairInt xy = centroids.get(i);
            double dist = distance(indexXY, xy);
            if (dist > maxDist) {
                continue;
            }

            int[] dimensions = calcDimensions(orderedBoundaries.get(i));
            int area = (int)Math.round(Math.sqrt(dimensions[0] * dimensions[0] +
                dimensions[1] * dimensions[1]));
            if (area > areaFactor*areaT || (dimensions[0] > dimensionsT[0]) ||
                (dimensions[1] > dimensionsT[1])) {
                continue;
            }

            indexes.add(i);
        }

        System.out.println("indexes=" + indexes);

        for (int i = 0; i < indexes.size(); ++i) {

            int idx1 = indexes.get(i);

            if (!adjacencyMap.containsKey(idx1)) {
                continue;
            }

            VeryLongBitString bs1 = new VeryLongBitString(nB);
            bs1.setBit(idx1);

            Double cost1 = aggregatedCostMap.get(bs1);

            for (int j = 0; j < indexes.size(); ++j) {

                int idx2 = indexes.get(j);

                if (!adjacencyMap.containsKey(idx1) ||
                    !adjacencyMap.get(idx1).contains(idx2)) {
                    continue;
                }

                if (idx1 == idx2) {
                    PairInt key = new PairInt(idx1, idx1);
                    distMap.put(key, cost1);
                    indexesMap.put(key, bs1);
                    continue;
                }

                VeryLongBitString bs2 = bs1.copy();
                bs2.setBit(idx2);
                
                if (containsASkipItem(bs2)) {
                    continue;
                }
                
                int[] bs2Bits = bs2.getSetBits();
                
                // ---- filter out very different color histogram groups ---
                fillCache1(bs2Bits);
                float intersection = ch.intersection(cache1, templateHSVHist);
                if (intersection < chHSVLimit1) {
                    // skip this in future
                    skipSet.add(bs2);
                    continue;
                }
                // --------
                
                PairInt key2;
                key2 = new PairInt(idx1, idx2);
                indexesMap.put(key2, bs2);

                PairIntArray boundary12 = aggregatedBoundaries.get(bs2);
                Double cost12 = aggregatedCostMap.get(bs2);

                if (cost12 == null) {

                    if (boundary12 == null) {
                        boundary12 = mergeAdjacentOrderedBorders(bs2);
                    }

                    cost12 = calcAndStoreMatchCost(boundary12, bs2);
                }

                if (cost12 != null) {

                    distMap.put(key2, cost12);

                    //DEBUG
                    if (expectedIndexes.contains(idx1) ||
                        expectedIndexes.contains(idx2)) {
                        System.out.println("contains an expected: " +
                            " key2=" + key2 + " cost=" + cost12);
                    }
                }
            }
        }
        // these skipSet items that could be filtered before input is
        // given to this class instance
        System.out.println("0 skipSet.size=" + skipSet.size());
        // ---- end initializing local maps

        PairInt minCostIdx = null;
        double minCost = Double.MAX_VALUE;

        /*
        s0 = costIJ;
        s1 = costIK + costKJ;
        */

        for (int i0 = 0; i0 < indexes.size(); ++i0) {
            int k = indexes.get(i0);
            for (int i1 = 0; i1 < indexes.size(); ++i1) {
                int i = indexes.get(i1);
                for (int i2 = 0; i2 < indexes.size(); ++i2) {
                    int j = indexes.get(i2);

                    // set these by adjacency.  false if not adjacent
                    boolean tIJ = true;
                    boolean tIK = true;
                    boolean tKJ = true;

                    if (!adjacencyMap.containsKey(i) ||
                        !adjacencyMap.get(i).contains(j)) {
                        tIJ = false;
                    }
                    if (!adjacencyMap.containsKey(i) ||
                        !adjacencyMap.get(i).contains(k)) {
                        tIK = false;
                    }
                    if (!adjacencyMap.containsKey(k) ||
                        !adjacencyMap.get(k).contains(j)) {
                        tKJ = false;
                    }

                    if (!tIJ && (!tIK || !tKJ)) {
                        continue;
                    }

                    boolean setPrev = true;
                    if (i == j) {
                        setPrev = false;
                    }

                    PairInt keyIJ;
                    keyIJ = new PairInt(i, j);

                    if (i == j) {
                        // this was set in initialization block above
                        if (tIJ) {
                            assert(distMap.containsKey(keyIJ));
                            assert(indexesMap.containsKey(keyIJ));
                            if (distMap.get(keyIJ) < minCost) {
                                minCost = distMap.get(keyIJ);
                                minCostIdx = keyIJ;
                                assert(indexesMap.get(minCostIdx) != null);
                            }
                        }
                        continue;
                    }

                    PairInt keyIK, keyKJ;
                    keyIK = new PairInt(i, k);
                    keyKJ = new PairInt(j, k);

                    Double s0 = null;
                    Double s1 = null;
                    //s0 = costIJ;
                    //s1 = costIK + costKJ;

                    VeryLongBitString bs0 = null;
                    VeryLongBitString bs1 = null;

                    if (tIJ) {
                        if (distMap.containsKey(keyIJ)) {
                            s0 = distMap.get(keyIJ);
                            bs0 = indexesMap.get(keyIJ);
                        } else {
                            PairInt p = new PairInt(keyIJ.getX(), keyIJ.getX());
                            if (indexesMap.get(p) == null) {
                                bs0 = new VeryLongBitString(nB);
                                bs0.setBit(keyIJ.getX());
                                //indexesMap.put(p, bs0.copy());
                            } else {
                                bs0 = indexesMap.get(p).copy();
                            }
                            bs0.setBit(keyIJ.getY());
                            if (!aggregatedCostMap.containsKey(bs0)) {
                                PairIntArray boundary = orderedBoundaries.get(keyIJ.getX());
                                if (keyIJ.getX() != keyIJ.getY()) {
                                    boundary = mergeAdjacentOrderedBorders(bs0);
                                }
                                Double s = calcAndStoreMatchCost(boundary, bs0);
                                if (s == null) {
                                    tIJ = false;
                                }
                            }
                            if (tIJ) {
                                s0 = aggregatedCostMap.get(bs0);
                            }
                        }
                        if (tIJ) {
                            if (containsASkipItem(bs0)) {
                                tIJ = false;
                            }
                        }
                        if (tIJ) {
                            // ---- filter out very different color histogram groups ---
                            int[] bs0Bits = bs0.getSetBits();
                            fillCache1(bs0Bits);
                            float intersection = ch.intersection(cache1, templateHSVHist);
                            if (intersection < chHSVLimit2) {
                                // skip this in future
                                skipSet.add(bs0);
                                tIJ = false;
                            }
                            // --------
                        }
                    }
                    
                    VeryLongBitString bsIK = indexesMap.get(keyIK);
                    VeryLongBitString bsKJ = indexesMap.get(keyKJ);
                    if (tIK && tKJ) {
                        //s1 = costIK + costKJ;
                        if (bsIK == null) {
                            bsIK = new VeryLongBitString(nB);
                            bsIK.setBit(keyIK.getX());
                            bsIK.setBit(keyIK.getY());
                        }
                        if (bsKJ == null) {
                            bsKJ = new VeryLongBitString(nB);
                            bsKJ.setBit(keyKJ.getX());
                            bsKJ.setBit(keyKJ.getY());
                        }
                        if (containsASkipItem(bsIK)) {
                            tIK = false;
                        }
                    }
                    if (tIK && tKJ) {
                        if (containsASkipItem(bsKJ)) {
                            tKJ = false;
                        }
                    }
                    if (tIK && tKJ) {
                        // ---- filter out very different color histogram groups ---
                        int[] bsIKBits = bsIK.getSetBits();
                        fillCache1(bsIKBits);
                        float intersection = ch.intersection(cache1, templateHSVHist);
                        if (intersection < chHSVLimit2) {
                            // skip this in future
                            skipSet.add(bsIK);
                            tIK = false;
                        }
                    }
                    if (tIK && tKJ) {
                        int[] bsKJBits = bsKJ.getSetBits();
                        fillCache1(bsKJBits);
                        float intersection = ch.intersection(cache1, templateHSVHist);
                        if (intersection < chHSVLimit2) {
                            // skip this in future
                            skipSet.add(bsKJ);
                            tKJ = false;
                        }
                        // --------
                    }
                    
                    if (tIK && tKJ) {
                        bs1 = bsIK.or(bsKJ);
                        if (containsASkipItem(bs1)) {
                            tIK = false;
                            tKJ = false;
                        }
                    }
                    if (tIK && tKJ) {   
                        PairIntArray bIKKJ = aggregatedBoundaries.get(bs1);
                        s1 = aggregatedCostMap.get(bs1);
                        if (s1 == null) {
                            bIKKJ = mergeAdjacentOrderedBorders(bs1);
                            Result r = null;
                            if (bIKKJ.getN() > 6) {
                                PartialShapeMatcher matcher = new PartialShapeMatcher();
                                matcher.overrideSamplingDistance(dp);
                                if (useSameNSampl) {
                                    matcher.setToUseSameNumberOfPoints();
                                }
                                r = matcher.match(template, bIKKJ);
                            }

                            if (r == null) {
                                tIK = false;
                                tKJ = false;
                            } else {
                                 
                                int nI = r.getNumberOfMatches();
                                float f = 1.f - ((float)nI/(float)nT1);
                                double d = r.getChordDiffSum()/maxDiffChordSum;
                                double c = f * f + d * d;

                                if (useEuclideanInCost) {
                                    double dist = calcTransformationDistanceSum(r,
                                        template, bIKKJ, true)/maxDistTransformSum;
                                    c += (dist * dist);
                                }
                                
                                s1 = Double.valueOf(c);
                                
                                aggregatedBoundaries.put(bs1, bIKKJ);
                                aggregatedCostMap.put(bs1, s1);
                                aggregatedResultMap.put(bs1, r);
                            }
                        }
                    }

                    if (tIJ && s0 != null && s1 != null &&
                        ((s0.doubleValue() <= s1.doubleValue()) || (!tIK || !tKJ))) {
                        assert(tIJ);

                        if ((distMap.containsKey(keyIJ) && (distMap.get(keyIJ) > s0))
                            || !distMap.containsKey(keyIJ)) {
                            distMap.put(keyIJ, s0);
                            indexesMap.put(keyIJ, bs0);
                        }

                        //DEBUG
                        {
                            if (expectedIndexes.contains(i) ||
                                expectedIndexes.contains(j)) {
                                if (keyIJ.getX() <= keyIJ.getY()) {
                                    System.out.println("contains an expected: " +
                                        " keyIJ=" + keyIJ.getX() + "," 
                                        + keyIJ.getY() + " cost=" + s0);
                                } else {
                                    System.out.println("contains an expected: " +
                                        " keyIJ=" + keyIJ.getY() + "," 
                                        + keyIJ.getX() + " *cost=" + s0);
                                }
                            }
                        }
                        
                    } else if (tIK && tKJ && (s1 != null)) {
                        assert(tIK);
                        assert(tKJ);
                        assert(bs1 != null);
                        
                        if ((distMap.containsKey(keyIJ) && (distMap.get(keyIJ) > s1))
                            || !distMap.containsKey(keyIJ)) {
                            distMap.put(keyIJ, s1);
                            indexesMap.put(keyIJ, bs1);
                        }
                        
                        
                        {//DEBUG
                            if (expectedIndexes.contains(i) ||
                                expectedIndexes.contains(j)
                                || expectedIndexes.contains(j)) {
                                if (keyIJ.getX() <= keyIJ.getY()) {
                                    System.out.println("contains an expected: " +
                                        " keyIJ=" + keyIJ.getX() + "," 
                                        + keyIJ.getY() + " ik+kh, keyIJ cost=" + s1);
                                } else {
                                    System.out.println("contains an expected: " +
                                        " keyIJ=" + keyIJ.getY() + "," 
                                        + keyIJ.getX() + " *ik+kh, keyIJcost=" + s1);
                                }
                            }
                        }

                    } else {
                        continue;
                    }

                    if (indexesMap.containsKey(keyIJ)) {
                        assert (indexesMap.get(keyIJ) != null);
                        if (distMap.get(keyIJ) < minCost) {
                            minCost = distMap.get(keyIJ);
                            minCostIdx = keyIJ;
                        }
                    }
                }
            }
        }

        if (minCostIdx == null) {
            return null;
        }

        VeryLongBitString minBitString = indexesMap.get(minCostIdx);
        
        {//DEBUG
            System.out.println("skipSet.size=" + skipSet.size());
            // expected:
            Double cost = aggregatedCostMap.get(tBS);
            if (cost == null) {
                PairIntArray b =  mergeAdjacentOrderedBorders(tBS);
                cost = calcAndStoreMatchCost(b, tBS);
            }
            fillCache1(tBS.getSetBits());
            float intersection = ch.intersection(cache1, templateHSVHist);
            
            System.out.println("expected true answer cost=" + cost +
                " " + Arrays.toString(tBS.getSetBits()) + " intersection=" +
                intersection);
            for (Entry<PairInt, VeryLongBitString> entry : indexesMap.entrySet()) {
                if (entry.getValue().equals(tBS)) {
                    System.out.println("Expected found in distMap is " +
                        distMap.get(entry.getKey()));
                }
            }

            //found
            int[] minidxs = indexesMap.get(minCostIdx).getSetBits();
            for (int idx : minidxs) {
                PairIntArray pts = orderedBoundaries.get(idx);
                System.out.println("CEN=" +
                    Arrays.toString
                    (curveHelper.calculateXYCentroids(pts))
                    + " n=" + pts.getN()
                );
            }
            
            fillCache1(minBitString.getSetBits());
            intersection = ch.intersection(cache1, templateHSVHist);

            System.out.println("minCostIdx=" + minCostIdx
                + " minCost=" + minCost + " bs="
                + Arrays.toString(minBitString.getSetBits()) + " intersection="
                + intersection);
        }
        
        return minBitString;
    }

    private double distance(PairInt xy1, PairInt xy2) {

        int diffX = xy1.getX() - xy2.getX();
        int diffY = xy1.getY() - xy2.getY();
        double dist = Math.sqrt(diffX * diffX + diffY * diffY);

        return dist;
    }

    private Double calcAndStoreMatchCost(PairIntArray boundary, 
        VeryLongBitString bs) {

        aggregatedBoundaries.put(bs, boundary);

        if (boundary.getN() < 6) {
            return null;
        }

        PartialShapeMatcher matcher = new PartialShapeMatcher();
        matcher.overrideSamplingDistance(dp);
        if (useSameNSampl) {
            matcher.setToUseSameNumberOfPoints();
        }
        Result result12 = matcher.match(template, boundary);
        if (result12 == null) {
            return null;
        }
        aggregatedResultMap.put(bs, result12);

        int nT1 = template.getN();
        // NOTE: this may need to be revised.  using the max diff chord
        // sum from only the single segmented cell matches as
        // the normalization for the Salukwdze distance^2.
        // May need to re-examine the bounds values and adjust the
        // normalizations.
        int nI = result12.getNumberOfMatches();
        float f = 1.f - ((float) nI / (float) nT1);
        double d = result12.getChordDiffSum() / maxDiffChordSum;
        double s = (f * f + d * d);
        
        if (useEuclideanInCost) {
            double dist = calcTransformationDistanceSum(result12,
                template, boundary, true) / maxDistTransformSum;
            s += (dist * dist);
        }
        
        Double cost12 = Double.valueOf(s);

        aggregatedCostMap.put(bs, cost12);

        return cost12;
    }

    private PairIntArray mergeAdjacentOrderedBorders(VeryLongBitString bs) {

        int[] setBits = bs.getSetBits();
        Set<PairInt> allPoints = new HashSet<PairInt>();
        for (int idx : setBits) {
            Set<PairInt> set = pointsList.get(idx);
            allPoints.addAll(set);
        }

        System.out.println("allPoints.size=" + allPoints.size()
        + " setBits=" + Arrays.toString(setBits));

        PerimeterFinder2 f2 = new PerimeterFinder2();

        return f2.extractOrderedBorder(allPoints);
    }

    /**
     * uses the euclidean transformation on the correspondence list
     * in r.  if useLimits is set, and if the calculated transformation
     * parameters' scale is out of range, then maxDistTransformSum is returned,
     * else, the summed differences are returned.
     * @param r
     * @param p
     * @param q
     * @param useLimits
     * @return 
     */
    private double calcTransformationDistanceSum(Result r, PairIntArray p,
        PairIntArray q, boolean useLimits) {
        
        if (r.getTransformationParameters() == null) {
            return Double.MAX_VALUE;
        }
        
        PairIntArray left = new PairIntArray(r.getNumberOfMatches());
        PairIntArray right = new PairIntArray(r.getNumberOfMatches());
        for (int i = 0; i < r.getNumberOfMatches(); ++i) {
            int idx = r.getIdx1(i);
            left.add(p.getX(idx), p.getY(idx));
            idx = r.getIdx2(i);
            right.add(q.getX(idx), q.getY(idx));
        }
        
        int tolerance = 5;
        
        EuclideanEvaluator evaluator = new EuclideanEvaluator();
        EuclideanTransformationFit fit = evaluator.evaluate(left,
            right, r.getTransformationParameters(), tolerance);
        
        if (useLimits) {
            TransformationParameters params = fit.getTransformationParameters();
            float scale = params.getScale();
            if (scale < euclideanScaleRange[0] || scale > euclideanScaleRange[1]) {
                return maxDistTransformSum;
            }
        }
        
        List<Double> distances = fit.getErrors();
        
        double sum = 0;
        for (Double d : distances) {
            sum += d;
        }
        
        return sum;
    }
    
    private void fillCache1(int[] indexes) {
        clear(cache1);
        for (int idx : indexes) {
            ch.add2To1(cache1, imgHSVHist.get(idx).a);
        }
    }
    
    private void clear(int[][] cache) {
        for (int i = 0; i < cache.length; ++i) {
            Arrays.fill(cache[i], 0);
        }
    }
    
    /**
     * checks whether skipSet contains bs or whether bs is composed
     * of any member in skipSet.
     * 
     * @param bs
     * @return 
     */
    private boolean containsASkipItem(VeryLongBitString bs) {
        
        if (skipSet.contains(bs)) {
            return true;
        }
        
        // TODO: this could probably be improved, but for now, expecting
        // that the number of items in skipSet is small.
        // if that's not the case, then more filtering could be done to the
        // input to this instance.
        
        for (VeryLongBitString skip : skipSet) {
            
            VeryLongBitString intersection = bs.and(skip);
            
            if (intersection.getNSetBits() == skip.getNSetBits()) {
                return true;
            }
        }
        
        return false;
    }

    private int findSmallestXYCentroid(List<Integer> indexes,
        List<PairInt> centroids) {
        
        int minX = Integer.MAX_VALUE;
        int minY = Integer.MAX_VALUE;
        int minIdx = -1;
        
        for (Integer index : indexes) {
            PairInt pCen = centroids.get(index.intValue());
            int x = pCen.getX();
            int y = pCen.getY();
            if (x < minX) {
                minX = x;
                minY = y;
                minIdx = index.intValue();
            } else if ((x == minX) && (y < minY)) {
                minX = x;
                minY = y;
                minIdx = index.intValue();
            }
        }
        
        return minIdx;
    }
}
