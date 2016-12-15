package algorithms.imageProcessing.matching;

import algorithms.MultiArrayMergeSort;
import algorithms.QuickSort;
import algorithms.compGeometry.LinesAndAngles;
import algorithms.imageProcessing.SummedColumnTable;
import algorithms.imageProcessing.features.RANSACEuclideanSolver;
import algorithms.imageProcessing.transform.EuclideanTransformationFit;
import algorithms.imageProcessing.transform.MatchedPointsTransformationCalculator;
import algorithms.imageProcessing.transform.TransformationParameters;
import algorithms.imageProcessing.transform.Transformer;
import algorithms.misc.Misc;
import algorithms.misc.MiscDebug;
import algorithms.search.KNearestNeighbors;
import algorithms.util.CorrespondencePlotter;
import algorithms.util.PairFloat;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import gnu.trove.iterator.TObjectFloatIterator;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TObjectFloatMap;
import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.map.hash.TObjectFloatHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.TreeMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import thirdparty.edu.princeton.cs.algs4.Interval;
import thirdparty.edu.princeton.cs.algs4.IntervalRangeSearch;

/**
 NOTE: NOT READY FOR USE YET... still testing.

<pre>
based upon algorithm in paper
 "Efficient Partial Shape Matching
    of Outer Contours: by Donoser
     - called IS-Match, integral shape match
     - a silhouette of ordered points are sampled
         making it an "order preserved assignment problem".
       - a chord angle descriptor is local and global and
         is invariant to similarity transformations.
       - the method returns partial sub matches
         so works with articulated data and occluded shapes
       - uses an efficient integral image based matching algorithm
       - the multi-objective optimization uses principles of
         Paretto efficiency, defined with the fraction of the
         total matched and the summed differences of angles.
       - the final result returned is the sequences and
         the total fraction matched and summed absolute differences,
         instead of the Salukwadze distance of a Paretto frontier.

       * point sampling:
         (a) same number of points over each contour
             - can handle similarity transforms.
               disadvantage is that it is less able to
               handle occlusion or extraneous shapes in
               the shape, though the euclidean projection
               used within, works around this.
             - to use in this mode:
                   setToUseSameNumberOfPoints()
                      and the number of points can be set 
                      using dp, with nSample = min(p.n, q.n)/dp
                      and
                   overrideSamplingDistance(dp)
         (b) OR, equidistant points
             - can handle occlusion.
               disadvantage in matching when there are
               scale differences between the shapes is
               handled in part by an internal euclidean
               projection.
             - the default sampling of pixels is 3.
               this can be changed using
                   overrideSamplingDistance(dp)

       The runtime complexity for building the integral
       image is O(m*n) where n and m are the number of sampled
       points on the input shapes.

       The runtime complexity for the search of the
       integral image of summed differences and analysis
       will be added here:
 </pre>
 <em>NOTE: You may need to pre-process the shape points
     for example, smooth the boundary.</em>
 <pre>
     This method:  
        PairIntArray p = imageProcessor
            .extractSmoothedOrderedBoundary()
        uses a Gaussian smoothing of 2 sigma,
        but a smaller sigma can be specified.
  </pre>
  @author nichole
 */
public class PartialShapeMatcher2 {

    //TODO: edit to read only the block width of 
    //      coluns along columns of matrix...reading
    //      over the block size square 
    //      includes a less strict range of indexes.
    //      and in that case, can replace the area summed table
    //      with a column summed table.  that results in
    //      2 steps for extrcting the sum instead of four.
    //      - could also edit to use the range search intervals
    //        while reading from the matrix and use the
    //        salukwzde cost for decisions about intersections.
    //        it's the same logic as throughout the code but
    //        moved up front and condensed
    //TODO: articulated match needs improvements
    
    /**
     * in sampling the boundaries of the shapes, one can
     * choose to use the same number for each (which can result
     * in very different spacings for different sized curves)
     * or one can choose a set distance between sampling
     * points.
     * dp is the set distance between sampling points.
       The authors of the paper use 3 as an example.
     */
    protected int dp = 3;

    private boolean useSameNumberOfPoints = false;

    private boolean srchForArticulatedParts = false;

    // this helps to remove points far from
    // euclidean transformations using RANSAC.
    // it should probably always be true.
    private boolean performEuclidTrans = true;

    private float pixTolerance = 20;

    // 10 degrees is 0.1745 radians
    // for a fit to a line, consider 1E-9
    private float thresh = 0.01f;//(float)(Math.PI/180.) * 10.f;

    private int minLength = 3;
    
    protected Logger log = Logger.getLogger(this.getClass().getName());

    private boolean debug = false;

    public void setToArticulatedMatch() {
        srchForArticulatedParts = true;
    }
    
    /**
     * turn off the euclidean transformation filter and addition of points.
     * NOTE that this should probably not normally be used.
     */
    public void _overrideToDisableEuclideanMatch() {
        performEuclidTrans = false;
    }
    
    /**
     * override the threshold for using a chord differernce value
     * to this.   By default it is set to a generous 10 degree, but
     * to fit a perfect line or similar, one may want to reduce this
     * threshold to 1E-9 or so.
     * @param t 
     */
    public void _overrideToThreshhold(float t) {
        this.thresh = t;
    }

    /**
    if this is set, the same number of points
    are used to sample both shapes.
    The number of points is min(p.n, q.n)/dp.
    You can change dp from the default of
    3 by using the method overrideSamplingDistance(dp).
    */
    public void setToUseSameNumberOfPoints() {
        useSameNumberOfPoints = true;
    }

    public void overrideSamplingDistance(int d) {
        this.dp = d;
    }

    public void setToDebug() {
        debug = true;
        log.setLevel(Level.FINE);
    }

    /**
      NOT READY FOR USE... still testing...

      A shape is defined as the clockwise ordered sequence
      of points P_1...P_N
      and the shape to match has points Q_1...Q_N.
      The spacings used within this method are equidistant
      unless changed using method setToUseSameNumberOfPoints().
      The default spacing is 3, 
      so override that if a different number
      is needed.
      
     <em>NOTE: You may need to pre-process the shape points
     for example, smooth the boundary.</em>
     <pre>
     This method:  
        PairIntArray p = imageProcessor
            .extractSmoothedOrderedBoundary()
        uses a Gaussian smoothing of 2 sigma,
        but a smaller sigma can be specified.
      </pre>
     @param p
     @param q
    */
    public Result match(PairIntArray p, PairIntArray q) {

        log.info("p.n=" + p.getN() + " q.n=" + q.getN()
            + " useSameNumberOfPoints=" + useSameNumberOfPoints
            + " dp=" + dp);

        if (p.getN() < 2 || q.getN() < 2) {
            throw new IllegalArgumentException("p and q must "
            + " have at least dp*2 points = " + (dp * 2));
        }

        if (useSameNumberOfPoints) {
            return matchSameNumber(p, q);
        }
        
        if (dp == 1) {
            return match0(p, q);
        }
     
        PairIntArray pSub = new PairIntArray(p.getN()/dp);
        PairIntArray qSub = new PairIntArray(q.getN()/dp);
    
        for (int i = 0; i < p.getN(); i += dp) {
            pSub.add(p.getX(i), p.getY(i));
        }
        
        for (int i = 0; i < q.getN(); i += dp) {
            qSub.add(q.getX(i), q.getY(i));
        }
        
        log.info("pSub.n=" + pSub.getN() + " qSub.n=" + qSub.getN());
        
        Result rSub = match0(pSub, qSub);
        
        if (rSub == null) {
            return null;
        } 
        
        // -- put results back into frame of p and q --
          
        // TODO: consider an option to further match the
        // points between correspondence if dp > 1
      
        Result r = new Result(p.getN(), q.getN(), rSub.getOriginalN1());

        for (int i = 0; i < rSub.idx1s.size(); ++i) {
            int idx1 = rSub.idx1s.get(i);
            int idx2 = rSub.idx2s.get(i);
            idx1 *= dp;
            idx2 *= dp;
            if (idx1 >= p.getN()) {
                int z = 0;
            }
            if (idx2 >= q.getN()) {
                int z = 0;
            }
            r.idx1s.add(idx1);
            r.idx2s.add(idx2);
        }
        assert(assertIndexesWithinBounds(r, p.getN(), q.getN()));
        r.chordDiffSum = rSub.chordDiffSum;
        r.distSum = rSub.distSum;
        
        if (rSub.articulatedSegment != null) {
            r.articulatedSegment = new TIntArrayList();
            for (int i = 0; i < rSub.idx1s.size(); ++i) {
                int segIdx = rSub.articulatedSegment.get(i);
                r.articulatedSegment.add(segIdx);
            }
        }
        
        if (rSub.getTransformationParameters() != null) {
       
            /*
            transX = xt0 -
                (xc*scale + (((x0-xc)*scale*math.cos(theta))
                + ((y0-yc)*scale*math.sin(theta)))

            transY = yt0 -
                (yc*scale + ((-(x0-xc)*scale*math.sin(theta))
                + ((y0-yc)*scale*math.cos(theta)))            
            */
            
            // translation increases by scale factor change, dp
            // rotation shouldn't change
            // scale changes by factor dp
            
            TransformationParameters params = rSub.getTransformationParameters().copy();
            params.setScale(params.getScale());
            params.setTranslationX(params.getTranslationX());
            params.setTranslationY(params.getTranslationY());
            r.setTransformationParameters(params);
            
        }
        
        return r;
    }

    private Result matchSameNumber(PairIntArray p, PairIntArray q) {

        log.fine("p.n=" + p.getN() + " q.n=" + q.getN());

        if (p.getN() < 2 || q.getN() < 2) {
            throw new IllegalArgumentException("p and q must "
            + " have at least dp*2 points = " + (dp * 2));
        }

        int nSampl = Math.min(p.getN(), q.getN())/dp;

        PairIntArray pSub = new PairIntArray(nSampl);
        PairIntArray qSub = new PairIntArray(nSampl);

        int pDp = p.getN()/nSampl;
        int qDp = q.getN()/nSampl;
    
        for (int i = 0; i < p.getN(); i += pDp) {
            pSub.add(p.getX(i), p.getY(i));
        }
        
        for (int i = 0; i < q.getN(); i += qDp) {
            qSub.add(q.getX(i), q.getY(i));
        }
        
        log.fine("pSub.n=" + pSub.getN() + " qSub.n=" + qSub.getN());
        
        Result rSub = match0(pSub, qSub);
        
        if (rSub == null) {
            return null;
        }
        
        // -- put results back into frame of p and q --
          
        // unrealistic value to ensure it's not used.
        int offset = Integer.MAX_VALUE;

        Result r = new Result(p.getN(), q.getN(), rSub.getOriginalN1());
            
        for (int i = 0; i < rSub.idx1s.size(); ++i) {
            int idx1 = rSub.idx1s.get(i);
            int idx2 = rSub.idx2s.get(i);
            idx1 *= pDp;
            idx2 *= qDp;
            if (idx1 >= p.getN()) {
                int z = 0;
            }
            if (idx2 >= q.getN()) {
                int z = 0;
            }
            r.idx1s.add(idx1);
            r.idx2s.add(idx2);
        }
        r.chordDiffSum = rSub.chordDiffSum;
        
        if (rSub.articulatedSegment != null) {
            r.articulatedSegment = new TIntArrayList();
            for (int i = 0; i < rSub.idx1s.size(); ++i) {
                int segIdx = rSub.articulatedSegment.get(i);
                r.articulatedSegment.add(segIdx);
            }
        }
        
        if (rSub.getTransformationParameters() != null) {
            
            /*
            two different scale factors need to be applied to the
                parameters, pDp and qDp.
            
            transX = xt0 -
                (xc*scale + (((x0-xc)*scale*math.cos(theta))
                + ((y0-yc)*scale*math.sin(theta)))

            transY = yt0 -
                (yc*scale + ((-(x0-xc)*scale*math.sin(theta))
                + ((y0-yc)*scale*math.cos(theta)))            
            
            if use 0,0 for origin, the equations simplify to:
                transX = xt0 - (x0*scale)
                transY = yt0 - (y0*scale)
            
            transforming x0,y0 by qDp and xt0,yt0 by pDp:                 
                new transX = (xt0 * pDp) - (x0 * qDp * scale)
                new transY = (yt0 * pDp) - (y0 * qDp * scale)
            
            looks like need to recalc transformation
            */
            
            PairIntArray left = new PairIntArray(r.idx1s.size());
            PairIntArray right = new PairIntArray(r.idx2s.size());
            assert(r.idx1s.size() == r.idx2s.size());
            for (int i = 0; i < r.idx1s.size(); ++i) {
                int idx = r.idx1s.get(i);
                left.add(p.getX(idx), p.getY(idx));
                idx = r.idx2s.get(i);
                right.add(q.getX(idx), q.getY(idx));
            }
            
            MatchedPointsTransformationCalculator tc = 
                new MatchedPointsTransformationCalculator();
            
            TransformationParameters params = tc.calulateEuclideanWithoutFilter(
                left, right, 0, 0);
            
            r.setTransformationParameters(params);            
        }

        return r;
    }
    
    private Result match0(PairIntArray p, PairIntArray q) {

        if (p == null || p.getN() < 2) {
            throw new IllegalArgumentException("p must have at "
                + "least 2 points");
        }
        
        if (q == null || q.getN() < 2) {
            throw new IllegalArgumentException("q must have at "
                + "least 2 points");
        }
        
        // --- make difference matrices ---

        //md[0:n2-1][0:n1-1][0:n1-1]
        float[][][] md;
        Result r;
        if (p.getN() <= q.getN()) {
            md = createDifferenceMatrices(p, q);
            r = match0(md, p, q);
            if (r != null) {
                System.out.println("not transposed");
                assert(assertIndexesWithinBounds(r, p.getN(), q.getN()));
            }
        } else {
            md = createDifferenceMatrices(q, p);
            r = match0(md, q, p);
            if (r != null) {
                System.out.println("transpose");
                r = r.transpose();
                assert(assertIndexesWithinBounds(r, p.getN(), q.getN()));
            }
        }

        return r;
    }
        
    private Result match0(float[][][] md, PairIntArray p, PairIntArray q) {

        if (p == null || p.getN() < 2) {
            throw new IllegalArgumentException("p must have at "
                + "least 2 points");
        }
        
        if (q == null || q.getN() < 2) {
            throw new IllegalArgumentException("q must have at "
                + "least 2 points");
        }
        
        if (p.getN() > q.getN()) {
            throw new IllegalArgumentException(
            "q.n must be >= p.n");
        }
        
        int n1 = p.getN();
        int n2 = q.getN();

        /*
        TODO: edit these comments for the latest code
        
        the matrices in md can be analyzed for best
        global solution and/or separately for best local
        solution.

        This method will return results for a local
        solution to create the point correspondence list.

        Note that the local best could be two different
        kinds of models, so might write two
        different methods for the results.
        (1) the assumption of same object but with some
            amount of occlusion and maybe additional
            shapes present due to segmentation not being
            able to isolate the object completely.
        (2) the assumption of same object but with
           some parts being differently oriented, for
           an example, the scissors opened versus closed.
           The occlusion should be handled for this one too.

        For the multi-objective optimization cost,
        need the sum of differences of chords and the fraction
        of the whole to calculate the Salukwdze distance 
        of the Paretto frontier.
        */
        
        List<SR> minima = findMinima(md, n1, n2);
        
        int topK = 20;
        if (minima.size() < topK) {
            topK = minima.size();
        }
        
        minima = minima.subList(0, topK);
                
        /*
        TODO: need to restore the options
            espec euclid transformation
        need to also add for the best result a check whether
            projection effects may be causing missing points
            due to foreshortening, etc, and an attempt to
            correct those.
        
        NOTE: tests show that most of the objects have the correct main
        offset (and hence correspoendence) for the top result, but not
        always (sometimes due to order...for articulated objects w/ 
        patterns that are similar enough to look like they repeat, 
        an articulated section that is the wrong match may have equal
        cost as the similar pattern which is the right match, but the
        adjacent points which are articulated differently do not
        lead to a longer match for the true pattern segment...
        so for articulated patterns like a clover leaf for example,
        2 leafs which are not adjacent might be unchanged 
        in q w.r.t. shape p, and they are
        separated by 2 articulated leafs, so matching the overall pattern
        requires an analysis that prefers a geometry (rigid or rigid + tolerance
        for projection) and then finds
        the matching "interleaving" patterns that remain unmatched after the
        geometric match.
        
        NOTE also that the results depend upon the threshold too.  For the scissors
        tests, a smaller threshold near 0.2 works better that twice that size,
        but for the noise gingerbread man profiles which have even been smoothed,
        a larger threshold has better results.
        
        Therefore, can see that an evaulation of the topK is how to find the best
        solution among the top intervals and then search for remaining results
        consistent with the user set options and remaining unmatched points.
        
        TODO: for the andr 02 statue test, the larger projection
        requires some edits in this....
        */
        
        if (srchForArticulatedParts) {
        
            // NOT YET IMPL
            
            OrderedClosedCurveCorrespondence occ = 
                new OrderedClosedCurveCorrespondence();
        
            occ.setMinimumLength(minLength);

            //TODO: for articulated, may want to consider an option here
            //   where a partial interval is rejected if more than the endpoints
            //   are trimmed.  currently, the interval, if not entirely within
            //   an unmatched and consistent range, is filtered to the subset of 
            //   the interval which does fit.

            occ.addIntervals(minima.subList(0, 10), n1, n2);

            List<SR> results = occ.getResultsAsList();

            Result best = new Result(n1, n2, n1);
            for (int i = 0; i < results.size(); ++i) {

                SR sr = results.get(i);

                for (int idx1 = sr.startIdx1; idx1 <= sr.stopIdx1; ++idx1) {
                    int idx2 = idx1 + sr.offsetIdx2;
                    if (idx2 > (n2 - 1)) {
                        idx2 -= n2;
                    }
                    best.idx1s.add(idx1);
                    best.idx2s.add(idx2);
                }
            }

            //throw new UnsupportedOperationException("not yet implmented");

            return best;

        } else if (performEuclidTrans) {

            // solve for transformation, add points near projection,
            // return sorted solutions, best is at top.
            // note that RANSAC has been used to remove outliers
            // from the already matched points too

            // the added points from the transformation
            List<PairIntArray> addedPoints = new ArrayList<PairIntArray>(topK);
                
            List<Result> results = transformAndEvaluate(minima, p, q,
                md, pixTolerance, topK, addedPoints);

            if (results != null) {
                for (int i = 0; i < topK; ++i) {
                //for (int i = 0; i < 1; ++i) {
                    Result result = results.get(i);
                    populateWithChordDiffs(result, md, p.getN(), q.getN());
                    log.fine("calc'ed chords: " + result.toStringAbbrev());
                    assert(assertIndexesWithinBounds(result, p.getN(), q.getN()));
                }
                
                double maxChordDiff = findMaxChordDiff(results);
                
                Collections.sort(results, new ResultsComparator(maxChordDiff, n1));
            }
            
            if (results == null) {
                return null;
            }
            
            return results.get(0);

        }
        
        throw new UnsupportedOperationException("not yet implemented");        
    }

    private List<SR> findMinima(float[][][] md, int n1, int n2) {

        //q.n is >= p.n, that is n2 >= n1
            
        // reading over a range of window sizes to keep the 
        // sum/nPix below thresh and keeping the mincost solutions.

        // find the intervals of contiguous minima and assign 
        // curve indexes to the largest segments.
        // (note that the objective formula for the cost
        // is the Salukwzde distance).

        List<SR> allResults = new ArrayList<SR>();
        
        Interval<Integer> interval = null;
        
        // this is learned from the first search at offset=0
        double maxChordSum = Double.MIN_VALUE;
        float[] outC = new float[2];
        double ds1, ds2;
        
        SummedColumnTable sct = new SummedColumnTable();
        
        /*
            MXM              NXN
                         30 31 32 33
         20 21 22        20 21 22 23
         10 11 12        10 11 12 13
         00 01 02        00 01 02 03   p_i_j - q_i_j

                         01 02 03 00
         20 21 22        31 32 33 30
         10 11 12        21 22 23 20
         00 01 02        11 12 13 10  p_i_j - q_(i+1)_(j+1)

                         12 13 10 11
         20 21 22        02 03 00 01
         10 11 12        32 33 30 31
         00 01 02        22 23 20 21  p_i_j - q_(i+2)_(j+2)

                         23 20 21 22
         20 21 22        13 10 11 12
         10 11 12        03 00 01 02
         00 01 02        33 30 31 32  p_i_j - q_(i+3)_(j+3)
        */
           
        for (int offset = 0; offset < md.length; offset++) {
            
            float[][] a = md[offset];

            List<Interval<Integer>> outputIntervals = new ArrayList<Interval<Integer>>();
            List<SR> outputValues = new ArrayList<SR>();
            
            search(a, outputIntervals, outputValues, offset);
            
            if (offset == 0) {
                maxChordSum = findMaxDiffChordSum(outputValues);
            }
            
            // the intervals are in the reference frame of p and of q shifted by
            // offset.
                        
            allResults.addAll(outputValues);
        }
        
        for (SR sr: allResults) {
            sr.maxChordSum = maxChordSum;
        }
        
        // sort by salukwzde distance
        Collections.sort(allResults, new SRComparator());
        
        System.out.println("nIntervals to proces=" + allResults.size());
       
        return allResults;
        
    }
    
    private double calcSalukDist(double compChord, double maxChord,
        int length, int maxLength) {
        double d = compChord/maxChord;
        double f = 1. - ((double)length/(double)maxLength);
        return f*f + d*d;
    }

    private Result smallNumberEuclidean(SR sr, PairIntArray p, PairIntArray q, 
        float pixTol, PairIntArray outputAddedPoints) {

        int len = sr.mLen;
        if (len < 2) {
            return null;
        }
        
        PairIntArray leftXY = new PairIntArray(len);
        PairIntArray rightXY = new PairIntArray(len);
        PairIntArray leftUnmatchedXY = new PairIntArray(p.getN() - len);
        PairIntArray rightUnmatchedXY = new PairIntArray(q.getN() - len);

        populatePointArrays(sr, p, q, leftXY, rightXY,
            leftUnmatchedXY, rightUnmatchedXY);

        
        MatchedPointsTransformationCalculator tc = 
            new MatchedPointsTransformationCalculator();        
                
        // NOTE: could use the method which removes outliers also, but the
        //    the contiguous points in sr are a single unit essentially
        
        TransformationParameters params = tc.calulateEuclideanWithoutFilter(
            leftXY, rightXY, 0, 0);
        
        if (params == null) {
            return null;
        }
        
        // since this class is using equidistant
        // points on shape boundary, need to compare
        // for same scale.
        if (params.getScale() < 0.9 || params.getScale() > 1.1) {
            log.fine("WARNING: " +
                " euclidean transformation scale: "  + params);
        }
        
        // ---- apply transformation -----
        Result result = applyTransformation(params, p, q,
            leftXY, rightXY, leftUnmatchedXY, rightUnmatchedXY, 
            pixTol, outputAddedPoints);
    
        if (debug) {
            String debugTag = Integer.toString(MiscDebug.getCurrentTimeFormatted());
            try {
                CorrespondencePlotter plotter = new CorrespondencePlotter(p, q);
                for (int i = 0; i < result.getNumberOfMatches(); ++i) {
                    int idx1 = result.getIdx1(i);
                    int idx2 = result.getIdx2(i);
                    int x1 = p.getX(idx1);
                    int y1 = p.getY(idx1);
                    int x2 = q.getX(idx2);
                    int y2 = q.getY(idx2);
                    if ((i % 4) == 0) {
                        plotter.drawLineInAlternatingColors(
                            x1, y1, x2, y2, 0);
                    }
                }
                String filePath = plotter.writeImage("_"
                    + "_debug4_" + debugTag);
            } catch (Throwable t) {
            }
        }
        
        return result;
    }
    
    private Result applyTransformation(TransformationParameters params,
        PairIntArray p, PairIntArray q,
        PairIntArray leftXY, PairIntArray rightXY,
        PairIntArray leftUnmatchedXY, PairIntArray rightUnmatchedXY, 
        float pixTol, PairIntArray outputAddedPoints) {
        
        Transformer transformer = new Transformer();
        PairIntArray leftTr = transformer.applyTransformation(
            params, leftUnmatchedXY);
        PairIntArray leftTr0 = transformer.applyTransformation(
            params, leftXY);

        String debugTag = Integer.toString(MiscDebug.getCurrentTimeFormatted());
        if (false && debug) {
            try {
                CorrespondencePlotter plotter = new CorrespondencePlotter(p, q);
                for (int i = 0; i < leftXY.getN(); ++i) {
                    int x1 = leftXY.getX(i);
                    int y1 = leftXY.getY(i);
                    int x2 = rightXY.getX(i);
                    int y2 = rightXY.getY(i);
                    if ((i % 5) == 0) {
                        plotter.drawLineInAlternatingColors(x1, y1,
                            x2, y2, 0);
                    }
                }
                String filePath = plotter.writeImage("_"
                    + "_debug1_" + debugTag);
                plotter = new CorrespondencePlotter(p, q);
                for (int i = 0; i < leftXY.getN(); ++i) {
                    int x1 = leftXY.getX(i);
                    int y1 = leftXY.getY(i);
                    int x2 = leftTr0.getX(i);
                    int y2 = leftTr0.getY(i);
                    if ((i % 5) == 0) {
                        plotter.drawLine(x1, y1, x2, y2,
                            255, 0, 0, 0);
                    }
                }
                filePath = plotter.writeImage("_"
                    + "_debug2_" + debugTag);

                plotter = new CorrespondencePlotter(p, q);
                for (int i = 0; i < leftUnmatchedXY.getN(); ++i) {
                    int x1 = leftUnmatchedXY.getX(i);
                    int y1 = leftUnmatchedXY.getY(i);
                    int x2 = leftTr.getX(i);
                    int y2 = leftTr.getY(i);
                    if ((i % 5) == 0) {
                        plotter.drawLine(x1, y1, x2, y2,
                            255, 0, 0, 0);
                    }
                }

                filePath = plotter.writeImage("_"
                    + "_debug3_" + debugTag);
            } catch (Throwable t) {
            }
        }

        log.fine(debugTag + " params=" + params);

        // find the best matches to the unmatched in
        // q

        // optimal is currently hungarian, but
        //     may replace w/ another bipartite matcher
        //     in future

        //boolean useOptimal = false;

        TObjectFloatMap<PairInt> idxMap;
        //if (useOptimal) {
        //    idxMap = optimalMatch(leftTr, rightUnmatched,
        //        pixTol);
        //} else {
            idxMap = nearestMatch(leftTr, rightUnmatchedXY, pixTol);
        //}

        log.fine(debugTag + "transformation nearest matches=" +
            idxMap.size());

        /*
        combine results from leftXY-rightXY
        with idxMap where idxMap is
            idx1, idx2 of leftTr and rightUnmatched, resp.
                          left
        */

        TObjectIntMap<PairInt> pPoints = Misc.createPointIndexMap(p);
        TObjectIntMap<PairInt> qPoints = Misc.createPointIndexMap(q);

        Result result = new Result(p.getN(), q.getN(), p.getN());
            result.setTransformationParameters(params);
        
        TObjectFloatIterator<PairInt> iter = idxMap.iterator();
        for (int i = 0; i < idxMap.size(); ++i) {
            iter.advance();
            PairInt idxIdx = iter.key();
            float dist = iter.value();
            PairInt ell = new PairInt(
                leftUnmatchedXY.getX(idxIdx.getX()),
                leftUnmatchedXY.getY(idxIdx.getX())
            );
            int pIdx = pPoints.get(ell);
            assert(pIdx > -1);
            PairInt ar = new PairInt(
                rightUnmatchedXY.getX(idxIdx.getY()),
                rightUnmatchedXY.getY(idxIdx.getY())
            );
            int qIdx = qPoints.get(ar);
            assert(qIdx > -1);

            result.insert(pIdx, qIdx, dist);

            outputAddedPoints.add(pIdx, qIdx);
        }

        for (int i = 0; i < rightXY.getN(); ++i) {
            int diffX = leftTr0.getX(i) - rightXY.getX(i);
            int diffY = leftTr0.getY(i) - rightXY.getY(i);
            float dist = (float)Math.sqrt(diffX * diffX +
                diffY * diffY);

            PairInt ell = new PairInt(leftXY.getX(i), leftXY.getY(i));
            int pIdx = pPoints.get(ell);
            assert(pIdx > -1);

            PairInt ar = new PairInt(rightXY.getX(i), rightXY.getY(i));
            int qIdx = qPoints.get(ar);
            assert(qIdx > -1);

            result.insert(pIdx, qIdx, dist);
        }

        if (debug) {
            try {
                CorrespondencePlotter plotter = new CorrespondencePlotter(p, q);
                for (int i = 0; i < result.getNumberOfMatches(); ++i) {
                    int idx1 = result.getIdx1(i);
                    int idx2 = result.getIdx2(i);
                    int x1 = p.getX(idx1);
                    int y1 = p.getY(idx1);
                    int x2 = q.getX(idx2);
                    int y2 = q.getY(idx2);
                    if ((i % 4) == 0) {
                        plotter.drawLineInAlternatingColors(
                            x1, y1, x2, y2, 0);
                    }
                }
                String filePath = plotter.writeImage("_"
                    + "_debug4_" + debugTag);
            } catch (Throwable t) {
            }
        }

        return result;
    }
    
    private TObjectFloatMap<PairInt> nearestMatch(
        PairIntArray xy1, PairIntArray xy2,
        double tolerance) {

        KNearestNeighbors knn = null;

        if (xy1.getN() > 3) {
            knn = new KNearestNeighbors(
                xy2.getX(), xy2.getY());
        }

        int k = 3;

        TObjectIntMap<PairInt> indexMap2 =
            Misc.createPointIndexMap(xy2);

        float[] dist = new float[3*xy1.getN()];
        int[] idx2s = new int[dist.length];
        int[] idx1s = new int[dist.length];
        int count = 0;
        for (int i = 0; i < xy1.getN(); ++i) {
            int x = xy1.getX(i);
            int y = xy1.getY(i);
            List<PairFloat> nearest = null;
            if (knn != null) {
                nearest = knn.findNearest(k, x, y, (float)tolerance);
            } else {
                nearest = kNearestBruteForce(k, x, y,
                    (float)tolerance, xy2);
            }
            if (nearest == null) {
                continue;
            }
            // note nearest has already been filtered by tolerance
            for (PairFloat p : nearest) {
                int x2 = Math.round(p.getX());
                int y2 = Math.round(p.getY());
                double d = Math.sqrt(
                    distanceSqEucl(x, y, x2, y2));
                dist[count] = (float)d;
                idx1s[count] = i;
                idx2s[count] = indexMap2.get(new PairInt(x2, y2));
                count++;
            }
        }

        idx1s = Arrays.copyOf(idx1s, count);
        idx2s = Arrays.copyOf(idx2s, count);
        dist = Arrays.copyOf(dist, count);
        int[] ilu = new int[count];
        for (int i = 0; i < count; ++i) {
            ilu[i] = i;
        }
        QuickSort.sortBy1stArg(dist, ilu);

        TIntSet a1 = new TIntHashSet();
        TIntSet a2 = new TIntHashSet();

        TObjectFloatMap<PairInt> costMap =
            new TObjectFloatHashMap<PairInt>();

        for (int i = 0; i < dist.length; ++i) {
            int idx = ilu[i];
            int idx1 = idx1s[idx];
            int idx2 = idx2s[idx];
            if (a1.contains(idx1) || a2.contains(idx2)) {
                continue;
            }
            PairInt p = new PairInt(idx1, idx2);
            costMap.put(p, dist[i]);
            a1.add(idx1);
            a2.add(idx2);
        }

        return costMap;
    }

    private boolean assertIndexesWithinBounds(Result result, 
        int n1, int n2) {
        
        for (int i = 0; i < result.idx1s.size(); ++i) {
            int idx1 = result.idx1s.get(i);
            if (idx1 >= n1) {
                return false;
            }
            int idx2 = result.idx2s.get(i);
            if (idx2 >= n2) {
                return false;
            }
        }
        
        return true;
    }

    private double findMaxChordDiff(List<Result> results) {

        double max = Double.MIN_VALUE;
        for (Result r : results) {
            double d = r.chordDiffSum;
            if (d > max) {
                max = d;
            }
        }
        return max;
    }
    
    public static class ResultsComparator implements Comparator<Result> {

        final double maxChordDiffSum;
        final int maxNMatchable;
        
        public ResultsComparator(double maxChordDiff, int maxNumberMatchable) {
            this.maxChordDiffSum = maxChordDiff;
            this.maxNMatchable = maxNumberMatchable;
        }
        
        @Override
        public int compare(Result o1, Result o2) {
            double d1 = o1.calculateSalukwdzeDistanceSquared(maxChordDiffSum, 
                maxNMatchable);
            double d2 = o2.calculateSalukwdzeDistanceSquared(maxChordDiffSum, 
                maxNMatchable);
            if (d1 < d2) {
                return -1;
            } else if (d1 > d2) {
                return 1;
            }
            return 0;
        }
    
    }

    public static class SRComparator implements Comparator<SR> {

        @Override
        public int compare(SR o1, SR o2) {
            double d1 = o1.calcSalukDist();
            double d2 = o2.calcSalukDist();
            if (d1 < d2) {
                return -1;
            } else if (d1 > d2) {
                return 1;
            }
            return 0;
        }
    
    }

    private double findMaxDiffChordSum(List<SR> list) {
        double max = Double.MIN_VALUE;
        for (SR sr : list) {
            double v = sr.maxChordSum;
            if (v > max) {
                max = v;
            }
        }
        return max;
    }

    public static class SR implements Comparable<SR> {
        // range of match is idx1 thru idx2, inclusive and if idx2 < idx1,
        //   the interval has wrapped around the closed curve
        int startIdx1;
        int stopIdx1;
        int offsetIdx2;
        int row;
        int mLen;
        int nMax;
        double diffChordSum;
        double maxChordSum;
        boolean chordSumNeedsUpdate = true;

        @Override
        public boolean equals(Object obj) {
            if (!(obj instanceof SR)) {
                return false;
            }
            SR other = (SR)obj;
            return (other.startIdx1 == startIdx1 && other.stopIdx1 == stopIdx1 &&
                other.mLen == mLen && other.diffChordSum == diffChordSum &&
                other.row == row && other.offsetIdx2 == offsetIdx2);
        }
        
        @Override
        public int hashCode() {
            int hash = fnvHashCode();
            return hash;
        }
        //Public domain:  http://www.isthe.com/chongo/src/fnv/hash_32a.c
        protected static int fnv321aInit = 0x811c9dc5;
        protected static int fnv32Prime = 0x01000193;
        protected int fnvHashCode() {
            int sum = fnv321aInit;
            sum ^= startIdx1;
            sum *= fnv32Prime;
            sum ^= stopIdx1;
            sum *= fnv32Prime;
            sum ^= offsetIdx2;
            sum *= fnv32Prime;
            sum ^= row;
            sum *= fnv32Prime;
            sum ^= mLen;
            sum *= fnv32Prime;
            sum ^= Double.hashCode(maxChordSum);
            sum *= fnv32Prime;
            return sum;
        }
        
        @Override
        public int compareTo(SR other) {
            double d1 = calcSalukDist(diffChordSum, maxChordSum, mLen, nMax);
            double d2 = calcSalukDist(other.diffChordSum, other.maxChordSum, 
                other.mLen, other.nMax);
            if (d1 < d2) {
                return -1;
            } else if (d1 > d2) {
                return 1;
            }
            return 0;
        }
        
        public void setChordSumNeedsUpdate(boolean needsUpdate) {
            chordSumNeedsUpdate = needsUpdate;
        }
        
        double calcSalukDist() {
            return calcSalukDist(diffChordSum, maxChordSum, 
                (stopIdx1 - startIdx1 + 1), mLen);
        }
        
        double calcSalukDist(double compChord, double maxChord,
            int length, int maxLength) {
            double d = compChord/maxChord;
            double f = 1. - ((double)length/(double)maxLength);
            return f*f + d*d;
        }
    }
    
    /**
     * search the difference chord sum matrix a for minimum Salukwzde cost
     * and put results in outputIntervals and the sum of the chord differences
     * in outputValues.
     * @param a
     * @param outputIntervals
     * @param outputValues 
     */
    private void search(float[][] a, List<Interval<Integer>> outputIntervals,
        List<SR> outputValues, int offset) {
    
        //q.n is >= p.n, that is, n2 >= n1
        int n1 = a.length;
        int n2 = a[0].length;

        // storing the interval of consecutive indexes, each below threshold,
        //   and storing as the value, the key to entry in intervalMap
        IntervalRangeSearch<Integer, SR> rangeSearch =
            new IntervalRangeSearch<Integer, SR>();

        Interval<Integer> interval = null;

        Set<PairInt> added = new HashSet<PairInt>();
 
        double maxChordSum = Double.MIN_VALUE;

        float[] outC = new float[2];
        int stop, start;
        
        SummedColumnTable sct = new SummedColumnTable();
        
        // using row major notation of a[row][col]
        for (int r = (n1 - 1); r >= minLength; --r) {
            for (int row = 0; row < n1; ++row) {
                for (int col = 0; col < n2; col += r) {
                    stop = col + r;
                    if (stop > (n2 - 1)) {
                        stop = n2 - 1;
                    }
                    sct.extractWindowInColumn(a, col, stop, row, outC);
                    if (outC[1] < 1) {
                        continue;
                    }
                    float d = outC[0]/outC[1];
                    if (debug) {
                        //System.out.println(String.format(
                        //"len=%d i=%d d=%.2f", (stop - j + 1), i, d));
                    }
                    if (
                        d > thresh 
                        //|| outC[0] > thresh3
                        //outC[0] > thresh
                        //|| d > (5.f * Math.PI/180.)
                        ) {
                        continue;
                    }
                    if (d > maxChordSum) {
                        maxChordSum = d;
                    }
                    start = col;
                    if (stop == start) {
                        // do not store single index matches
                        continue;
                    }

                    int ni = stop - start + 1;
                    if (ni < minLength) {
                        continue;
                    }

                    interval =  new Interval<Integer>(start, stop);
                    PairInt s = new PairInt(start, stop);
                    if (added.contains(s)) {
                        continue;
                    }

                    //a, col, stop, row, outC
                    SR sr = new SR();
                    sr.startIdx1 = col;
                    sr.stopIdx1 = stop;
                    sr.offsetIdx2 = offset;
                    sr.row = row;
                    sr.diffChordSum = d;
                    sr.setChordSumNeedsUpdate(false);
                    sr.maxChordSum = maxChordSum;
                    sr.mLen = ni;
                    sr.nMax = n2;
                
                    boolean didIns = rangeSearch.putIfLessThan(interval, sr, sr);

                    added.add(s);
                }
            } // end loop j
        } // end r
                
        rangeSearch.getAllIntervals(outputIntervals, outputValues);
    
        assert(outputIntervals.size() == outputValues.size());
    }

    /**
     * create the matrices of differences between p
     * and q.  Note that the matrix differences are
     * absolute differences.
     * index0 is rotations of q,  index1 is p.n, index2 is q.n
      returns a[0:q.n-1][0:p.n-1][0:p.n-1]
    */
    protected float[][][] createDifferenceMatrices(
        PairIntArray p, PairIntArray q) {

        if (p.getN() > q.getN()) {
            throw new IllegalArgumentException(
            "q.n must be >= p.n");
        }

        /*
        | a_1_1...a_1_N |
        | a_2_1...a_2_N |
               ...
        | a_N_1...a_N_N |
           elements on the diagonal are zero

           to shift to different first point as reference,
           can shift down k-1 rows and left k-1 columns.
        */

        //log.fine("a1:");
        float[][] a1 = createDescriptorMatrix(p, p.getN());

        //log.fine("a2:");
        float[][] a2 = createDescriptorMatrix(q, q.getN());

        /*
        - find rxr sized blocks similar to one another
          by starting at main diagonal element
          A_1(s,s) and A_2(m,m)
          which have a small angular difference value

          //s range 0 to M-1
          //m range 0 to M-1, M<=N
          //r range 2 to sqrt(min(N, M))

          - to calculate all D_a(s,m,r)
               uses concept of "integral image"
                   by Viola and Jones
               (for their data, I(x,y)=i(x,y)+I(x-1,y)+I(x,y-1)-I(x-1,y-1))
             - N integral images int_1...int_N of size MXM
               are built for N descriptor
               difference matrices M_D^n
               where the number of sampled points on the two
               shapes is N and M, respectively

               where M_D^n
                   = A_1(1:M,1:M) - A_2(n:n+M-1,n:n+M-1)

               then, all matching triplets {s,m,r} which
               provide a difference value D_a(s,m,r) below
               a fixed threshold are calculated.
        ---------------------------

        (1) make difference matrices.
            there will be N A_2 matrices in which each
            is shifted left and up by 1 (or some other value).

            M_D^n = A_1(1:M,1:M) - A_2(n:n+M-1,n:n+M-1)
                shifting A_2 by 0 through N covering all
                orientation angles.
        (2) make Summary Area Tables of the N M_D^m matrices.
        (3) search: starting on the diagonals of the integral images
            made from the N M_D^n matrices,
            D_α(s, m, r) can be calculated for every block of any
            size starting at any point on the diagonal in
            constant time.
        */

        /*
            MXM              NXN
                         30 31 32 33
         20 21 22        20 21 22 23
         10 11 12        10 11 12 13
         00 01 02        00 01 02 03   p_i_j - q_i_j

                         01 02 03 00
         20 21 22        31 32 33 30
         10 11 12        21 22 23 20
         00 01 02        11 12 13 10  p_i_j - q_(i+1)_(j+1)

                         12 13 10 11
         20 21 22        02 03 00 01
         10 11 12        32 33 30 31
         00 01 02        22 23 20 21  p_i_j - q_(i+2)_(j+2)

                         23 20 21 22
         20 21 22        13 10 11 12
         10 11 12        03 00 01 02
         00 01 02        33 30 31 32  p_i_j - q_(i+3)_(j+3)
        */

        // --- make difference matrices ---
        int n1 = p.getN();
        int n2 = q.getN();
        float[][][] md = new float[n2][][];
        float[][] prevA2Shifted = null;
        for (int i = 0; i < n2; ++i) {
            float[][] shifted2;
            if (prevA2Shifted == null) {
                shifted2 = copy(a2);
            } else {
                // shifts by 1 to left and up by 1
                rotate(prevA2Shifted);
                shifted2 = prevA2Shifted;
            }
            // NOTE: absolute values are stored.
            //M_D^n = A_1(1:M,1:M) - A_2(n:n+M-1,n:n+M-1)
            md[i] = subtract(a1, shifted2);
            assert(md[i].length == n1);
            assert(md[i][0].length == n1);
            prevA2Shifted = shifted2;
        }

        // ---- make summary area table for md-----
        for (int i = 0; i < md.length; ++i) {
            applySummedColumnTableConversion(md[i]);
        }

        return md;
    }

    /**
     given the shape points for p and q,
     create a matrix of descriptors, describing the difference
     in chord angles.

     The chord descriptor is invariant to translation, rotation,
     and scale:
       - a chord is a line joining 2 region points
       - uses the relative orientation between 2 chords
         angle a_i_j is from chord P_i_P_j to reference
         point P_i
         to another sampled point and chord P_j_P_(j-d) and P_j

         d is the number of points before j in the sequence of points P.

         a_i_j is the angle between the 2 chords P_i_P_j and P_j_P_(j-d)
    */
    protected float[][] createDescriptorMatrix(PairIntArray p, int n) {
        
        int dp1 = 1;

        float[][] a = new float[n][];
        for (int i = 0; i < n; ++i) {
            a[i] = new float[n];
        }

        /*
             P1      Pmid

                  P2
        */

        log.fine("n=" + n);

        for (int i1 = 0; i1 < n; ++i1) {
            int start = i1 + 1 + dp1;
            for (int ii = start; ii < (start + n - 1 - dp1); ++ii) {
                int i2 = ii;

                int imid = i2 - dp1;
                // wrap around
                if (imid > (n - 1)) {
                    imid -= n;
                }

                // wrap around
                if (i2 > (n - 1)) {
                    i2 -= n;
                }

                //log.fine("i1=" + i1 + " imid=" + imid + " i2=" + i2);

                double angleA = LinesAndAngles.calcClockwiseAngle(
                    p.getX(i1), p.getY(i1),
                    p.getX(i2), p.getY(i2),
                    p.getX(imid), p.getY(imid)
                );
                
                if (Double.isNaN(angleA)) {
                    if (i2 < i1 && i1 < imid) {
                        angleA = LinesAndAngles
                            .calcClockwiseAngle(
                            p.getX(i2), p.getY(i2),
                            p.getX(i1), p.getY(i1),
                            p.getX(imid), p.getY(imid));
                        if (Double.isNaN(angleA)) {
                            angleA = LinesAndAngles
                            .calcClockwiseAngle(
                            p.getX(i2), p.getY(i2),
                            p.getX(i1), p.getY(i1),
                            p.getX(imid), p.getY(imid));
                        }
                    } else {
                        System.out.println(
                        "SKIP i1=" + i1 + " imid=" + imid + " i2=" + i2);
                        continue;
                    }
                }

                /*
                String str = String.format(
                    "[%d](%d,%d) [%d](%d,%d) [%d](%d,%d) a=%.4f",
                    i1, p.getX(i1), p.getY(i1),
                    i2, p.getX(i2), p.getY(i2),
                    imid, p.getX(imid), p.getY(imid),
                    (float) angleA * 180. / Math.PI);
                log.fine(str);
                */

                a[i1][i2] = (float)angleA;
                
                if (i2 == (i1 + 2)) {
                    // fill in missing point, assume same value
                    if (a[i1][i2] == a[i1][ii]) {
                        a[i1][i1 + 1] = a[i1][i2];
                    }
                }
            }
        }

        return a;
    }

    protected int distanceSqEucl(int x1, int y1, int x2, int y2) {
        int diffX = x1 - x2;
        int diffY = y1 - y2;
        return (diffX * diffX + diffY * diffY);
    }

    private float[][] copy(float[][] a) {
        float[][] a2 = new float[a.length][];
        for (int i = 0; i < a2.length; ++i) {
            a2[i] = Arrays.copyOf(a[i], a[i].length);
        }
        return a2;
    }

    private void rotate(float[][] prevShifted) {

         // shift x left by 1 first
         for (int y = 0; y < prevShifted[0].length; ++y) {
             float tmp0 = prevShifted[0][y];
             for (int x = 0; x < (prevShifted.length- 1); ++x){
                 prevShifted[x][y] = prevShifted[x + 1][y];
             }
             prevShifted[prevShifted.length - 1][y] = tmp0;
         }

         // shift y down by 1
         for (int x = 0; x < prevShifted.length; ++x) {
             float tmp0 = prevShifted[x][0];
             for (int y = 0; y < (prevShifted[x].length - 1); ++y){
                 prevShifted[x][y] = prevShifted[x][y + 1];
             }
             prevShifted[x][prevShifted[x].length - 1] = tmp0;
         }
    }

    /**
     * subtract the portion of a2 that is same size as
     * a1 from a1.
     * @param a1
     * @param a2
     * @return
     */
    private float[][] subtract(float[][] a1, float[][] a2) {

        /*
         MXM     NXN
                 20 21 22
         10 11   10 11 12
         00 01   00 01 02

                 01 02 00
         10 11   21 22 20
         00 01   11 12 10

                 12 10 11
         10 11   02 00 01
         00 01   22 20 21

        subtracting only the MXM portion
        */

        assert(a1.length == a1[0].length);
        assert(a2.length == a2[0].length);

        int n1 = a1.length;
        int n2 = a2.length;

        assert(n1 <= n2);

        float[][] output = new float[n1][];
        for (int i = 0; i < n1; ++i) {
            output[i] = new float[n1];
            for (int j = 0; j < n1; ++j) {
                float v = a1[i][j] - a2[i][j];
                if (v < 0) {
                    v *= -1;
                }
                output[i][j] = v;
            }
        }

        return output;
    }

    private void print(String label, float[][] a) {

        StringBuilder sb = new StringBuilder(label);
        sb.append("\n");

        for (int j = 0; j < a[0].length; ++j) {
            sb.append(String.format("row: %3d", j));
            for (int i = 0; i < a.length; ++i) {
                sb.append(String.format(" %.4f,", a[i][j]));
            }
            log.fine(sb.toString());
            sb.delete(0, sb.length());
        }
    }

    protected void applySummedColumnTableConversion(float[][] mdI) {

        int w = mdI.length;
        int h = mdI[0].length;
        
        // sum along columns, that is a[i][*]
        for (int i = 0; i < w; ++i) {
            for (int j = 0; j < h; ++j) {
                if (j > 0) {
                    mdI[i][j] += mdI[i][j - 1];
                }
            }
        }
    }

    private List<PairFloat> kNearestBruteForce(
        int k, int x, int y, float tolerance,
        PairIntArray xy2) {

        int nn = 0;
        float[] dist = new float[xy2.getN()];
        int[] indexes = new int[xy2.getN()];
        for (int i = 0; i < xy2.getN(); ++i) {
            int diffX = xy2.getX(i) - x;
            int diffY = xy2.getY(i) - y;
            dist[i] = (float) Math.sqrt(
                diffX * diffX + diffY * diffY);
            if (dist[i] <= tolerance) {
                nn++;
            }
            indexes[i] = i;
        }

        if (nn == 0) {
            return null;
        }

        QuickSort.sortBy1stArg(dist, indexes);

        if (k < nn) {
            nn = k;
        }

        List<PairFloat> list = new ArrayList<PairFloat>(nn);
        for (int i = 0; i < nn; ++i) {
            int idx = indexes[i];
            int x2 = xy2.getX(idx);
            int y2 = xy2.getY(idx);
            PairFloat p = new PairFloat(x2, y2);
            list.add(p);
        }

        return list;
    }

    public static class Result {
        protected double distSum = 0;
        protected double chordDiffSum = 0;
        protected boolean chordsNeedUpdates = true;
        protected TransformationParameters params = null;
        // indexes for the correspondence from shape 1
        protected TIntList idx1s = new TIntArrayList();
        // indexes for the correspondence from shape 2
        protected TIntList idx2s = new TIntArrayList();
        // if articulated search, this contains a segment number
        //   for each gap filled with contiguous correspondence.
        protected TIntList articulatedSegment = null;
        protected final int n1;
        protected final int n2;
        protected final int origN1;
        protected Object[] data = null;
        public Result(int n1, int n2, int origN1) {
            this.n1 = n1;
            this.n2 = n2;
            this.origN1 = origN1;
        }
        
        public int getArticulatedSegment(int idx) {
            return articulatedSegment.get(idx);
        }
        public void insert(int idx1, int idx2, float dist) {
            idx1s.add(idx1);
            idx2s.add(idx2);
            distSum += dist;
        }

        public int getIdx1(int index) {
            return idx1s.get(index);
        }
        public int getIdx2(int index) {
            return idx2s.get(index);
        }

        public int getNumberOfMatches() {
            return idx1s.size();
        }

        public float getFractionOfWhole() {
            return (float)idx1s.size()/(float)origN1;
        }

        protected double getNormalizedChordDiff(double maxChordSum) {
            double d = chordDiffSum/maxChordSum;
            return d;
        }
        
        public double getChordDiffSum() {
            return chordDiffSum;
        }

        /**
         * The Salukwdze distance is the metric used as a
         * cost in comparisons and this returns the square of
         * that (sqrt operation not performed).
         * (note that
         * the maximum chord sum that was used to determine
         * a best solution is not always stored, because the
         * Result correspondence lists grow afterwards depending
         * upon options).
         * @param maxChordSum
         * @return
         */
        public float calculateSalukwdzeDistanceSquared(double maxChordSum) {
            float f = 1.f - getFractionOfWhole();
            double d = getNormalizedChordDiff(maxChordSum);
            float s = (float)(f * f + d * d);
            return s;
        }
        
        /**
         * The Salukwdze distance is the metric used as a
         * cost in comparisons and this returns the square of
         * that (sqrt operation not performed).
         * (note that
         * the maximum chord sum that was used to determine
         * a best solution is not always stored, because the
         * Result correspondence lists grow afterwards depending
         * upon options)
         * @param maxChordSum
         * @return
         */
        public float calculateSalukwdzeDistanceSquared(double maxChordSum,
            int maxNumberOfMatchable) {
            float f = 1.f - ((float)getNumberOfMatches()/
                (float)maxNumberOfMatchable);
            double d = getNormalizedChordDiff(maxChordSum);
            float s = (float)(f * f + d * d);
            return s;
        }

        void addToChordDifferenceSum(float diff) {
            chordDiffSum += diff;
        }

        public int getOriginalN1() {
            return origN1;
        }
        
        public Object[] getData() {
            return data;
        }
        
        public void setData(Object[] theData) {
            this.data = theData;
        }
        
        public void setTransformationParameters(TransformationParameters
            euclidParams) {
            this.params = euclidParams;
        }
        
        public TransformationParameters getTransformationParameters() {
            return params;
        }
        
        // NOTE: need to complete this for some cases such as
        // same number of points sampling.
        public double getDistSum() {
            return distSum;
        }
        
        protected boolean assertIndexesWithinBounds() {
            for (int i = 0; i < idx1s.size(); ++i) {
                int idx1 = idx1s.get(i);
                if (idx1 >= n1) {
                    return false;
                }
                int idx2 = idx2s.get(i);
                if (idx2 >= n2) {
                    return false;
                }
            }
            return true;
        }

        /**
         * reverse the mappings from list 1 to list 2
         * to the reference frame of list 2 to list 1.
         * @return
         */
        public Result transpose() {

            assert(assertIndexesWithinBounds());
            
            Result t = new Result(n2, n1, n1);
            t.idx1s.addAll(idx2s);
            t.idx2s.addAll(idx1s);
            t.distSum = distSum;
            t.chordDiffSum = chordDiffSum;
            t.data = data;
            
            assert(t.assertIndexesWithinBounds());
            
            if (articulatedSegment != null) {
                t.articulatedSegment = new TIntArrayList(articulatedSegment);
            }
            
            if (params != null) {
                MatchedPointsTransformationCalculator tc =
                    new MatchedPointsTransformationCalculator();
                TransformationParameters params2 = 
                    tc.swapReferenceFrames(params);
                t.setTransformationParameters(params2);
            }

            return t;
        }

        @Override
        public String toString() {

            StringBuilder sb = new StringBuilder();
            sb.append(String.format(
                "nMatched=%d frac=%.4f distSum=%.4f dChordSum=%.4f",
                idx1s.size(), getFractionOfWhole(),
                (float)distSum, (float)chordDiffSum));
            if (params != null) {
                sb.append("\nparams=").append(params.toString());
            }
            sb.append("\n");

            for (int i = 0; i < idx1s.size(); ++i) {
                sb.append(String.format("%d (%d, %d)\n",
                    i, idx1s.get(i), idx2s.get(i)));
            }

            return sb.toString();
        }

        public String toStringAbbrev() {

            StringBuilder sb = new StringBuilder();
            sb.append(String.format(
                "nMatched=%d frac=%.4f dChordSum=%.4f distSum=%.4f ",
                idx1s.size(), getFractionOfWhole(),
                (float)chordDiffSum, (float)distSum));

            return sb.toString();
        }

        public void sortByIdx1() {
            QuickSort.sortBy1stArg(idx1s, idx2s);
        }

        /**
         * find within array idxs1 gaps of values
         * and store those in the output while also
         * storing the value and idx1 of the items
         * before and after the gap into
         * outputValueIdx1Map.  (outputValueIdx1Map
         * is for gap related predecessor and successor
         * lookups).
         * @param minGapSize
         * @param outputValueIdx1Map
         * @return
         */
        PairIntArray findGapRanges(int minGapSize,
            TreeMap<Integer, Integer> outputValueIdx1Map) {

            assert(minGapSize > 0);

            int n = idx1s.size();

            PairIntArray ranges = new PairIntArray();

            if (n < 2) {
                return ranges;
            }

            sortByIdx1();

            //  cases:
            //  0 1 - - 4
            //  - - 2 3
            //  0 1 - - -
            int prevI = idx1s.get(0);
            if (prevI > (minGapSize - 1)) {
                ranges.add(0, prevI - 1);
                outputValueIdx1Map.put(Integer.valueOf(prevI),
                    Integer.valueOf(0));
                outputValueIdx1Map.put(Integer.valueOf(idx1s.get(1)),
                    Integer.valueOf(1));
            }
            for (int i = 1; i < idx1s.size(); ++i) {
                int cI = idx1s.get(i);
                if (cI > (prevI + minGapSize)) {
                   // startGap = prevI + 1 and end is cI - 1
                   ranges.add(prevI + 1, cI - 1);
                   outputValueIdx1Map.put(Integer.valueOf(prevI),
                      Integer.valueOf(i - 1));
                   outputValueIdx1Map.put(Integer.valueOf(cI),
                      Integer.valueOf(i));
                }
                prevI = cI;
            }
            if (prevI < (n1 - 1 - minGapSize)) {
                ranges.add(prevI + 1, n1 - 1);
                outputValueIdx1Map.put(Integer.valueOf(prevI),
                    Integer.valueOf(idx1s.size() - 1));
            }

            return ranges;
        }

        /**
         * find the gaps in idx1s values and return
         * the start and end of the gap values in
         * a range as pairs in PairIntArray.  Note that
         * the returned PairIntArray has been sorted
         * so that the longest gap ranges are at the
         * smallest indexes.
         * @param minGapSize
         * @param outputValueIdx1Map
         * @return 
         */
        private PairIntArray findGapRangesDescBySize(
            int minGapSize, TreeMap<Integer, Integer>
                outputValueIdx1Map) {

            PairIntArray gaps = findGapRanges(minGapSize,
                outputValueIdx1Map);

            int[] len = new int[gaps.getN()];
            int[] indexes = new int[len.length];
            for (int i = 0; i < len.length; ++i) {
                len[i] = gaps.getY(i) - gaps.getX(i) + 1;
                indexes[i] = i;
            }
            MultiArrayMergeSort.sortByDecr(len, indexes);

            PairIntArray out = new PairIntArray(gaps.getN());
            for (int i = 0; i < len.length; ++i) {
                int idx = indexes[i];
                out.add(gaps.getX(idx), gaps.getY(idx));
            }

            return out;
        }

        public TIntIntMap createValueIndexMap1() {
            TIntIntMap map = new TIntIntHashMap();
            for (int i = 0; i < idx1s.size(); ++i) {
                map.put(idx1s.get(i), i);
            }
            return map;
        }
    }

    private List<Result> transformAndEvaluate(
        List<SR> intervals, PairIntArray p, PairIntArray q,
        float[][][] md, float pixTol, int topK,
        List<PairIntArray> addedPointLists) {

        if (intervals.size() == 0) {
            return null;
        }

        List<Result> results = new ArrayList<Result>(topK);
        for (int i = 0; i < topK; ++i) {
            
            SR sr = intervals.get(i);
            
            System.out.println("srLen=" + sr.mLen);
            
            PairIntArray added = new PairIntArray();
            Result result = null;
            
            if (sr.mLen < 7) {
                //too few for the RANSAC algorithm
                result = smallNumberEuclidean(sr, p, q, pixTol, added);
            } else { 
                result = addByTransformation(sr, p, q, pixTol, added);
            }
            
            if (result != null) {
                assert(assertIndexesWithinBounds(result, p.getN(), q.getN()));
                results.add(result);
                addedPointLists.add(added);
            }
        }

        if (results.isEmpty()) {
            return null;
        }
        
        { //DEBUG, recalc the saluk score to see if top has changed
            //TODO:   
        }
     
        return results;
    }
    
    /**
     * 
     * @param result
     * @param md
     * @param n1
     * @param n2 
     */
    private void populateWithChordDiffs(Result result,
        float[][][] md, int n1, int n2) {

        assert(md[0][0].length == n1);
        
        //md[0:n2-1][0:n1-1][0:n1-1]
        
        // NOTE: this update method doesn't effectively use the summed column
        //    structure of md, so will consider changing
        //    the intermediate results to remain as SR intervals,
        //    even if the interval is a single pixel.
        //    and this argument receives a list of SR instead of a Result

        result.chordDiffSum = 0;
        
        SummedColumnTable sct = new SummedColumnTable();
        float[] output = new float[2];
        int row = 0;
        int offset = 0;
        int idx1, idx2;
        
        for (int i = 0; i < result.getNumberOfMatches(); ++i) {
            
            idx1 = result.getIdx1(i);
            idx2 = result.getIdx2(i);
            
            offset = 0;
            
            //NOTE: this is not the same as reading a block w.r.t. the original
            //   offset that was used
            
            if (idx2 > (n1 - 1)) {
                // idx2 is outside of the default row 0 of md[0] array, so
                // need to calculate the offset.
                
                // if within bounds, will add 4 to the offset to get the pixel
                // away from bounds
                
                offset = (idx2 - (n1 - 1)) + 4;
                idx2 -= offset;
                if (idx2 < 0) {
                    offset = idx2 - (n1 - 1);
                    idx2 -= offset;
                }
            }

            if (idx2 > (md[offset][0].length - 1)) {
                System.out.println("idx1=" + idx1 +
                    " idx2=" + idx2 + " n1=" + n1 + " n2=" + n2
                    + " md.length=" + md.length
                    + " md[0].length=" + md[0].length
                    + " md[0][0].length=" + md[0][0].length);
            } else if (idx1 > (md[offset][0].length - 1)) {
                System.out.println("**idx1=" + idx1 +
                    " idx2=" + idx2 + " n1=" + n1 + " n2=" + n2
                    + " md.length=" + md.length
                    + " md[0].length=" + md[0].length
                    + " md[0][0].length=" + md[0][0].length);
            }
            
            sct.extractWindowInColumn(md[offset], idx1, idx2, row, output);
            
            float d = output[0]/output[1];
            
            result.chordDiffSum += d;
        }
        
        result.chordsNeedUpdates = false;
    }

    protected void populatePointArrays(SR sr, PairIntArray p, PairIntArray q,
        PairIntArray pOut, PairIntArray qOut,
        PairIntArray pUnmatchedOut, PairIntArray qUnmatchedOut) {

        int offset = sr.offsetIdx2;
        int len = sr.mLen;

        int n1 = p.getN();
        int n2 = q.getN();

        TIntSet pIdxs = new TIntHashSet();
        TIntSet qIdxs = new TIntHashSet();

        for (int i = 0; i < len; ++i) {

            int pIdx = i;
            assert(pIdx < n1);

            pIdxs.add(pIdx);
            pOut.add(p.getX(pIdx), p.getY(pIdx));

            int qIdx = pIdx + offset;
            if (qIdx >= n2) {
                qIdx -= n2;
            }
            
            qIdxs.add(qIdx);
            qOut.add(q.getX(qIdx), q.getY(qIdx));
        }

        assert(pOut.getN() == len);
        assert(qOut.getN() == len);

        // -- populate unmatched p ----
        for (int i = 0; i < p.getN(); ++i) {
            if (!pIdxs.contains(i)) {
                pUnmatchedOut.add(p.getX(i), p.getY(i));
            }
        }
        
        for (int i = 0; i < q.getN(); ++i) {
            if (!qIdxs.contains(i)) {
                qUnmatchedOut.add(q.getX(i), q.getY(i));
            }
        }
        
        assert(pUnmatchedOut.getN() + pOut.getN() == p.getN());
        assert(qUnmatchedOut.getN() + qOut.getN() == q.getN());
    }
    
    private Result addByTransformation(SR sr, PairIntArray p, PairIntArray q,
        float pixTol, PairIntArray outputAddedPoints) {

        int len = sr.mLen;
        if (len < 2) {
            return null;
        }
        
        PairIntArray leftXY = new PairIntArray(len);
        PairIntArray rightXY = new PairIntArray(len);
        PairIntArray leftUnmatchedXY = new PairIntArray(p.getN() - len);
        PairIntArray rightUnmatchedXY = new PairIntArray(q.getN() - len);

        populatePointArrays(sr, p, q, leftXY, rightXY,
            leftUnmatchedXY, rightUnmatchedXY);

        log.fine(" lft.n=" + leftXY.getN()
            + " rgt.n=" + rightXY.getN());

        return addByTransformation(p, q,
            leftXY, rightXY, leftUnmatchedXY, rightUnmatchedXY,
            pixTol, outputAddedPoints);
    }
    
    private Result addByTransformation(
        PairIntArray p, PairIntArray q,
        PairIntArray left, PairIntArray right,
        PairIntArray leftUnmatched, PairIntArray rightUnmatched,
        float pixTol, PairIntArray outAddedPoints) {

        String debugTag = "offset=" + Integer.toString(
            MiscDebug.getCurrentTimeFormatted());

        PairIntArray outLeft = new PairIntArray(left.getN());
        PairIntArray outRight = new PairIntArray(left.getN());
        
        RANSACEuclideanSolver euclid =
            new RANSACEuclideanSolver();
        EuclideanTransformationFit fit = euclid.calculateEuclideanTransformation(
            left, right, outLeft, outRight);

        TransformationParameters params = (fit != null) ?
            fit.getTransformationParameters() : null;

        if (params == null) {
            //TODO: reconsider whether to package up
            // the given sequence s and return it here
            log.fine(debugTag + " no euclidean fit");
            return null;
        }

        // since this class is using equidistant
        // points on shape boundary, need to compare
        // for same scale.
        if (params.getScale() < 0.9 || params.getScale() > 1.1) {
            log.fine("WARNING: " + debugTag +
                " euclidean transformation scale: "  + params);
        }
        left = outLeft;
        right = outRight;

        log.fine(debugTag + " partial fit=" + fit.toString()
            + " params=" + params
            + " reset left.n=" + left.getN()
            + " right.n=" + right.getN());

        log.fine("dp=" + dp + " pixTol=" + pixTol);

        // ---- apply transformation -----
        Result result = applyTransformation(params, p, q,
            left, right, leftUnmatched, rightUnmatched, 
            pixTol, outAddedPoints);
        
        return result;        
    }

}
