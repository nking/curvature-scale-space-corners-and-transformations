package algorithms.imageProcessing.matching;

import algorithms.imageProcessing.MiscellaneousCurveHelper;
import algorithms.util.PairIntArray;
import gnu.trove.map.TIntObjectMap;
import java.util.List;
import algorithms.imageProcessing.matching.PartialShapeMatcher.Result;
import algorithms.misc.Misc;
import algorithms.misc.MiscMath;
import algorithms.util.OneDIntArray;
import algorithms.util.PairInt;
import algorithms.util.QuadInt;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Set;
import javax.naming.directory.SearchResult;
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
 * adjacency map.  
 *
 * Also note that the class uses internal cache which requires single
 * threaded use of this class.
 * 
 * @author nichole
 */
public class ShapeFinder {
    
    // partial shape matcher sampling distance
    private final int dp = 1;
    private final boolean useSameNSampl = false;

    private final PairIntArray bounds1;
    private final int[] ch1;
    private final float scale1;
    private final float sz1;
    private final List<Set<PairInt>> listOfSets2;
    private final List<OneDIntArray> listOfCH2s;
    private final float scale2;
    private final TObjectIntMap<int[]> keyIndexMap;
    private final TIntObjectMap<PairIntArray> indexBoundsMap;
    private final float intersectionLimit;
    
    private final TObjectIntMap<PairInt> pointIndexes2Map;
    private final TIntObjectMap<TIntSet> adj2Map;
    private final List<PairInt> xyCen2List;
    private final List<QuadInt> xyMinMax2List;
    
    private final int xMax1;
    private final int yMax1;
    private final int xMax2;
    private final int yMax2;
    
    /**
     * NOT READY FOR USE.
     * 
     * note that this is the data of one image in pyramid1 compared to the data
     * in one image of pyramid2 and it is assumed that one such pair will have
     * template and true match of nearly the same scale.  
     * On other words, this instance of ShapeFinder might not hold the true 
     * match at same scale, but another instance (holding other scles of the
     * pyramids) should if the true match is
     * present.
     * 
     * @param bounds1 clockwise ordered boundaries of template object.
     * @param ch1 1-D color histogram of template object
     * @param scale1 a bookkeeping number for the octave downsampling scale factor
     * with respect to the complete pyramid1.
     * @param listOfSets2 labeled point sets of dataset2 to be searched for template
     * object.listOfCH2s the color histograms for listOfSets2
     * @param scale2 a bookkeeping number for the octave downsampling scale factor
     * with respect to the complete pyramid2.
     * @param keyIndexMap an in/out variable for use in combination wtth
     * indexBoundsMap to cache the clockwise boundaries of aggregated adjacent
     * labeled cells.
     * The key = sorted indexes of listOfSets2, value = the lookup key for
     * indexBoundsMap,
     * @param indexBoundsMap an in/out variable for use in combination wtth
     * keyIndexMap to cache the clockwise boundaries of aggregated adjacent
     * labeled cells.
     * The key = the lookup key from keyIndexMap, value = clockwise ordered
     * bounding points of the aggregated sets of labels from the key in
     * indexBoundsMap,
     */
    public ShapeFinder(PairIntArray bounds1, int[] ch1, float scale1, float sz1,
        int xMax1, int yMax1,
        List<Set<PairInt>> listOfSets2, List<OneDIntArray> listOfCH2s,
        float scale2, TObjectIntMap<int[]> keyIndexMap,
        TIntObjectMap<PairIntArray>  indexBoundsMap,
        int xMax2, int yMax2,
        float intersectionLimit) {
        
        this.bounds1 = bounds1;
        this.ch1 = ch1;
        this.scale1 = scale1;
        this.sz1 = sz1;
        this.listOfSets2 = listOfSets2;
        this.listOfCH2s = listOfCH2s;
        this.scale2 = scale2;
        this.keyIndexMap = keyIndexMap;
        this.indexBoundsMap = indexBoundsMap;
        this.intersectionLimit = intersectionLimit;
        this.xMax1 = xMax1;
        this.yMax1 = yMax1;
        this.xMax2 = xMax2;
        this.yMax2 = yMax2;
        
        MiscellaneousCurveHelper curveHelper = new MiscellaneousCurveHelper();
        
        this.xyMinMax2List = new ArrayList<QuadInt>();
        this.xyCen2List = new ArrayList<PairInt>();
        this.pointIndexes2Map = new TObjectIntHashMap<PairInt>();
        for (int i = 0; i < listOfSets2.size(); ++i) {
            Set<PairInt> set = listOfSets2.get(i);
            for (PairInt p : set) {
                pointIndexes2Map.put(p, i);
            }
            double[] xyCen = curveHelper.calculateXYCentroids(set);
            xyCen2List.add(new PairInt(
                (int)Math.round(xyCen[0]), (int)Math.round(xyCen[1])));
            int[] minMaxXY = MiscMath.findMinMaxXY(set);
            xyMinMax2List.add(new QuadInt(minMaxXY[0], minMaxXY[1], minMaxXY[2],
                minMaxXY[3]));
        }
        
        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;
        this.adj2Map = new TIntObjectHashMap<TIntSet>();
        for (int idx1 = 0; idx1 < listOfSets2.size(); ++idx1) {
            Set<PairInt> set = listOfSets2.get(idx1);
            for (PairInt p : set) {
                int x = p.getX();
                int y = p.getY();
                // find adjacent labels that are not idx1
                for (int k = 0; k < dxs.length; ++k) {
                    int x2 = x + dxs[k];
                    int y2 = y + dys[k];
                    PairInt p2 = new PairInt(x2, y2);
                    int idx2 = pointIndexes2Map.get(p2);
                    if (idx2 == idx1) {
                        continue;
                    }
                    TIntSet set2 = adj2Map.get(idx1);
                    if (set2 == null) {
                        set2 = new TIntHashSet();
                        adj2Map.put(idx1, set2);
                    }
                    set2.add(idx2);
                }
            }
        }
    }

    /**
     * NOT READY FOR USE.
     * uses PartialShapeMatcher and search patterns to
     * find the best fitting match between a
     * group of adjacent segmented cells
     * to a template shape. 
     *
     * The runtime is dependent upon which of the 3 search
     * patterns is kept in the end, so that will be
     * filled in here after implementation and testing.
     *
     * @return
     */
    public ShapeFinderResult findAggregated() {

        /*
        some notes that will be updated when return to the search
        task.
        
        basically, have a global grid search,
        wrapping local FloydWarshall searches where search
            is the best fitting bounds of aggregated labels sets and
            their color histograms
        
        step (1) O(N_cells) + O(N_cells_in_bin * log_2(N_cells_in_bin)):
           - 2D bins, that is x and y of size sz1.
           - sort the cells within each bin by x centroid and ties by y centroid
        step (2) local search pattern
            The bigger picture is a sliding window of the 2D bin
            where the core of the search and aggregation is done
            within the 2D bin, but the aggregation can continue
            into the next bin to the right, next bin above,
            or above and to the right.
            
           this is the part that is changing.
           would like to use Dijkstra's within each bin to find the best 
             matching aggregation of cells, but unfortunately, the 
             aggregation of cell boundaries is not inductive, that is
             two wrong steps of aggregation might have better results than
             two correct steps.
           therefore, using the more time consuming but more complete
           Floyd-Warshal search within each 2D bin.
           Note that this might be improvable in the future.               
        */
        
        globalGridSearch();
           
        throw new UnsupportedOperationException("not yet implemented");
    }

    private void globalGridSearch() {
        
        int w = xMax2 + 1;
        int h = yMax2 + 1;
        
        int nPPB = 1;
        
        int nX = nPPB * ((int)Math.floor(w/sz1) + 1);
        int nY = nPPB * ((int)Math.floor(h/sz1) + 1);
        int binSz = (int)Math.ceil(sz1);
        
        QuadTree<Integer, Integer> centroidQT = new QuadTree<Integer, Integer>();
            
        // populated bin numbers, x, y
        Set<PairInt> binNumbers = new HashSet<PairInt>();
        
        // store the segmented cell centroids into quadtree
        for (int i = 0; i < xyCen2List.size(); ++i) {
            PairInt xy = xyCen2List.get(i);
            
            centroidQT.insert(xy.getX(), xy.getY(), Integer.valueOf(i));
            
            int xBin = xy.getX()/binSz;
            int yBin = xy.getY()/binSz;
            binNumbers.add(new PairInt(xBin, yBin));
        }
        
        // visit the bins by smallest x,y first,
        //   and start each search with the labeled set which has the
        //   smallest x,y in that bin.
        //   any search can agregate from an adjacent bin, but not beyond
        //   that so a restricted list of set indexes is created
        //   for each search start.
        
        // sort so that all x for same y are listed, then next x for larger y
        List<PairInt> sortedBinNumbers = new ArrayList<PairInt>(binNumbers);
        Collections.sort(sortedBinNumbers, new XYSort());
        
        for (int i = 0; i < sortedBinNumbers.size(); ++i) {
            PairInt xyBin = sortedBinNumbers.get(i);
            int startX = xyBin.getX() * binSz;
            int startY = xyBin.getY() * binSz;
            
            int stopX = startX + binSz;
            if (stopX > (w - 1)) {
                stopX = w - 1;
            }
            
            int stopY = startY + binSz;
            if (stopY > (h - 1)) {
                stopY = h - 1;
            }
            
            Interval<Integer> intX = new Interval<Integer>(startX, stopX);
            Interval<Integer> intY = new Interval<Integer>(startY, stopY);
            Interval2D<Integer> rect = new Interval2D<Integer>(intX, intY);
             
            // centroid list indexes (these are same indexes are for label2 sets)
            List<Integer> indexes = centroidQT.query2D(rect);
            assert(!indexes.isEmpty());
                
            // chose the smallest x, smallest y in bin
            int srchIdx = findSmallestXYCentroid(indexes, xyCen2List);
            
            QuadInt srchXYMinMax = xyMinMax2List.get(srchIdx);
            
            // aggregate adjacent indexes of larger x or larger y
            TIntSet adjIdxs = new TIntHashSet();
            adjIdxs.addAll(indexes);
            
            PairInt nextBin = new PairInt(xyBin.getX() + 1, xyBin.getY());
            if (binNumbers.contains(nextBin)) {
                Interval<Integer> intX2 = new Interval<Integer>(startX + binSz, 
                    stopX + binSz);
                Interval<Integer> intY2 = new Interval<Integer>(startY, stopY);
                Interval2D<Integer> rect2 = new Interval2D<Integer>(intX2, intY2);
                // centroid list indexes
                List<Integer> indexes2 = centroidQT.query2D(rect2);
                if (indexes2 != null) {
                    for (int idx2 : indexes2) {
                        // if extent + srch extent is > 1.15*sz1, exclude
                        int maxDim = minMaxXY(srchXYMinMax, 
                            xyMinMax2List.get(idx2));
                        if (maxDim <= 1.15 * sz1) {
                            adjIdxs.add(idx2);
                        }
                    }
                }
            }
            nextBin = new PairInt(xyBin.getX() + 1, xyBin.getY() + 1);
            if (binNumbers.contains(nextBin)) {
                Interval<Integer> intX2 = new Interval<Integer>(startX + binSz, 
                    stopX + binSz);
                Interval<Integer> intY2 = new Interval<Integer>(startY + binSz, 
                    stopY + binSz);
                Interval2D<Integer> rect2 = new Interval2D<Integer>(intX2, intY2);
                // centroid list indexes
                List<Integer> indexes2 = centroidQT.query2D(rect2);
                if (indexes2 != null) {
                    for (int idx2 : indexes2) {
                        int maxDim = minMaxXY(srchXYMinMax, 
                            xyMinMax2List.get(idx2));
                        if (maxDim <= 1.15 * sz1) {
                            adjIdxs.add(idx2);
                        }
                    }
                }
            }
            nextBin = new PairInt(xyBin.getX(), xyBin.getY() + 1);
            if (binNumbers.contains(nextBin)) {
                Interval<Integer> intX2 = new Interval<Integer>(startX, stopX);
                Interval<Integer> intY2 = new Interval<Integer>(startY + binSz, 
                    stopY + binSz);
                Interval2D<Integer> rect2 = new Interval2D<Integer>(intX2, intY2);
                // centroid list indexes
                List<Integer> indexes2 = centroidQT.query2D(rect2);
                if (indexes2 != null) {
                    for (int idx2 : indexes2) {
                        int maxDim = minMaxXY(srchXYMinMax, 
                            xyMinMax2List.get(idx2));
                        if (maxDim <= 1.15 * sz1) {
                            adjIdxs.add(idx2);
                        }
                    }
                }
            }
            
            // TODO: invoke the local search with restricted label range
            SearchResult r = searchUsingFloydWarshall(srchIdx, adjIdxs);
            
            //TODO: when replace aspectj w/ another AOP library, assert that 
            //   all list 2 set labels are covered in adjIdxs and srchIdx,
            //   after all FW invocations.
        }
        
        throw new UnsupportedOperationException(
            "Not supported yet."); 
    }

    private SearchResult searchUsingFloydWarshall(int srchIdx, TIntSet adjIdxs) {
        
        
        throw new UnsupportedOperationException(
            "Not supported yet.");
    }

    private int findSmallestXYCentroid(List<Integer> indexes, 
        List<PairInt> xyCenList) {
        
        int minX = Integer.MAX_VALUE;
        int minX_Y = Integer.MAX_VALUE;
        int minIdx = -1;
        
        for (int i = 0; i < indexes.size(); ++i) {
            int idx = indexes.get(i);
            int x = xyCenList.get(idx).getX();
            int y = xyCenList.get(idx).getY();
            if (x < minX) {
                minX = x;
                minX_Y = y;
                minIdx = idx;
            } else if (x == minX) {
                if (y < minX_Y) {
                    minX_Y = y;
                    minIdx = idx;
                }
            }
        }
        
        assert(minIdx > -1);
        
        return minIdx;
    }

    private int minMaxXY(QuadInt srchXYMinMax, QuadInt otherXYMinMax) {
        int diffX = Math.max(
            Math.abs(srchXYMinMax.getA() - otherXYMinMax.getB()),
            Math.abs(srchXYMinMax.getB() - otherXYMinMax.getA()));
        
        int diffY = Math.max(
            Math.abs(srchXYMinMax.getC() - otherXYMinMax.getD()),
            Math.abs(srchXYMinMax.getD() - otherXYMinMax.getC()));
        
        // consider taking the diagonal here
        return Math.max(diffX, diffY);
    }
    
    public static class ShapeFinderResult extends Result {
        
        protected float intersection;
        protected final PairIntArray bounds1;
        protected final PairIntArray bounds2;
        
        public ShapeFinderResult(int n1, int n2, int offset, PairIntArray bounds1,
            PairIntArray bounds2) {
            
            super(n1, n2, offset);
            
            this.bounds1 = bounds1;
            this.bounds2 = bounds2;
        }
        
        public ShapeFinderResult(int n1, int n2, int nOriginal, 
            int offset, PairIntArray bounds1,
            PairIntArray bounds2) {
            
            super(n1, n2, nOriginal, offset);
            
            this.bounds1 = bounds1;
            this.bounds2 = bounds2;
        }
    }
    
    // sorted so that all X of same Y are first, then all X of larger Y
    private class XYSort implements Comparator<PairInt> {

        @Override
        public int compare(PairInt o1, PairInt o2) {

            int comp = Integer.compare(o1.getY(), o2.getY());
            if (comp != 0) {
                return comp;
            }
            comp = Integer.compare(o1.getX(), o2.getX());
            return comp;
        }

    }
}
