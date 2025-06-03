package algorithms.imageProcessing.matching;

import algorithms.compGeometry.LinesAndAngles;
import algorithms.imageProcessing.SummedAreaTable;
import algorithms.misc.MiscMath;
import algorithms.signalProcessing.CurveResampler;
import algorithms.util.CorrespondencePlotter;
import algorithms.util.PairFloatArray;
import algorithms.util.PairIntArray;
import algorithms.util.PolygonAndPointPlotter;

import java.io.IOException;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
<pre>
based upon algorithm in paper
 "Efficient Partial Shape Matching
    of Outer Contours: by Donoser
  
     - called IS-Match, integral shape match
     - finds the best matching segments between 2 closed curves.
     - NOTE: to find the for a query and multiple targets, use the
       ...editing here to consider vector embedding approaches...
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
         using the Salukwadze distance of a Paretto frontier.

       * point sampling:
         For best results, the curves should have at least 30 points.
         (a) same number of points over each contour
             - can handle similarity transforms.
               disadvantage is that it is less able to
               handle occlusion or extraneous shapes in
               the shape.
               using the euclidean option works around this.
             - 
               to use in this mode:
                   setToUseSameNumberOfPoints()
                      and the number of points can be set 
                      using dp, with nSample = min(p.n, q.n)/dp
                      and
                   overrideSamplingDistance(dp)
               and optionally:
                   setToUseEuclidean()
         (b) OR, equidistant points (default)
             - can handle occlusion.
               disadvantage in matching when there are
               scale differences between the shapes is
               handled in part by an internal euclidean
               projection.
             - the default sampling of pixels is 3.
               this can be changed using
                   overrideSamplingDistance(dp)
               unit test using dp of 1 and 2 work well

       The runtime complexity for building the integral
       image is O(n*m*n) where m and n are the number of sampled
       points on the reduced input shapes.  The reduced size is the size after any transformations like
       resampling.

       The runtime complexity for the search depends upon the search choice:
          FAST: r.t.c. O(M*N*log(topK)).  this is the default.
          ALL_BLOCK_SIZES: r.t.c. O(M*N*DR * log(topK)).
          EXHAUSTIVE: the largest r.t.c. is exponential from recursively assembling all full curve matches
              composed of different block sizes and different number of gaps.  for 1 GB RAM max heap,
              one should use this only if number of points in the curve is around 20 or less.  If you have
              more RAM, use jvm args to set heap size min and max higher and then you'll have a higher point
              limit.  This is interesting for identifying the best complete partial matches of the two curves.

       <emph>The runtime complexity for the search of the
       integral image of summed differences and analysis,
       is n2 * (O(n1 * lg2(n1))</emph> where n1 and n2 are the number of points
       in the 2 shapes after the spacing has been considered.
       
       The algorithm runtime complexity could be reduced more, but with
       more loss in accuracy by selecting a discrete number of diagonal read 
       block sizes.
       For instance, if wanted to read only a single block size of
       n1/10, the total runtime complexity would be approx n2 * n1/10.
       (and in this case, could build the summed difference chord matrix
       smaller and more quickly than O(n2 * n1 * n1) before that.
       Essentially, would use summed column table and the rows would
       only be offsets of n1/10 so the total runtime complexity for building
       would be O(n1*(n1/10)) and reads would be smaller...it's effectively
       the result of comparing the 2 curves at offets of n1/10 intervals which
       would give fast but very rough results that would not handle articultion 
       well if at all.)
 </pre>
 
 <em>NOTE: You may need to pre-process the shape points
     for example, smooth the boundary.</em>
 <pre>
     This method:  
        PairIntArray p = imageProcessor
            .extractSmoothedOrderedBoundary()
        uses a Gaussian smoothing of 1 sigma,
  </pre>
  @author nichole
 */
public class PartialShapeMatcher {

    public static final float TWO_PI = (float)(2.*Math.PI);

    // the range of difference from TWO_PI for which an angle will be considered 0.
    public static final float ZERO_TOL = 0.00001f;

    /**
     * in sampling the boundaries of the shapes, one can
     * choose to use the same number for each (which can result
     * in very different spacings for different sized curves)
     * or one can choose a set distance between sampling
     * points.
     * dp is the set distance between sampling points.
       The authors of the paper use 3 as an example.

     The use of dp in this code is not well tested and could be improved.
    */
    protected int dp = 1;

    private boolean useSameNumberOfPoints = false;

    // 10 degrees is 0.1745 radians
    // for a fit to a line, consider 1E-9
    private float thresh = 1.f;//(float)(Math.PI/180.) * 10.f;

    private boolean overrideMinLength = false;
    private final int DEFAULT_MIN_LENGTH = 7;
    private int minLength = 7;

    private int topK = 1;

    /**
     * FAST: r.t.c. O(M*N*log(topK)).  this is the default.
     * ALL_BLOCK_SIZES: O(M*N*DR*log(topK))
     * EXHAUSTIVE: the largest r.t.c. from recursively assembling all full curve matches
     * composed of different block sizes and different number of gaps.  for 1 GB RAM max heap,
     * one should use this only if number of points in the curve is around 20 or less.  If you have
     * more RAM, use jvm args to set heap size min and max higher and then you'll have a higher point
     * limit.
     */
    private static enum SEARCH {
        FAST, ALL_BLOCK_SIZES,
        EXHAUSTIVE
    }
    private SEARCH search = SEARCH.FAST;

    protected Logger log = Logger.getLogger(this.getClass().getName());

    private boolean debug = false;

    /**
     * override the threshold for using a chord difference value
     * for the average value.   
     * By default the threshhold is set to 1.
     * @param t  threshhold to use
     */
    public void _overrideToThreshhold(float t) {
        this.thresh = t;
    }
    
    /**
     * override the default minimum length of 7.
     * @param length minimum length to use in matching blocks
     */
    public void overrideMinimumLength(int length) {
        this.minLength = length;
        this.overrideMinLength = true;
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

    /**
     * the default sampling distance is 1.  use this method to override it.
     * @param d set the sampling distance between points to use for curves.
     */
    public void overrideSamplingDistance(int d) {
        this.dp = d;
    }

    /**
     * change from the default fast search of r.t.c. O(M*N) to an ALL_BLOCK_SIZES search
     * with r.t.c. roughly  O(M*N*R*DR) where the block size R spans from the minimumLength (default 7) to
     * the maximum length N and uses a delta R that is log(N).
     */
    public void overrideToSearchAllBlockSizes(){
        this.search = SEARCH.ALL_BLOCK_SIZES;
    }

    /**
     * change from the default fast search of r.t.c. O(M*N) to an exhaustive search which has
     * the largest r.t.c. from recursively assembling all full curve matches
     * composed of different block sizes and different number of gaps.  for 1 GB RAM max heap,
     * one should use this only if number of points in the curve is around 20 or less.  If you have
     * more RAM, use jvm args to increase heap size min and max higher and then you'll have a higher point
     * limit than 20.
     */
    public void overrideToSearchExhaustive(){
        this.search = SEARCH.EXHAUSTIVE;
    }

    /**
     * return topK results.  By default topK is 1.
     * @param k the number of results to return
     */
    public void overrideTopK(int k) {
        if (k < 0) {
            throw new IllegalArgumentException("k must be >= 1");
        }
        topK = k;
    }

    /**
     * use this to enable the debug log comments and plots
     */
    public void setToDebug() {
        debug = true;
        log.setLevel(Level.FINE);
    }

    /**
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
     temporarily using exhaustive search.
     @param p a closed curve to match to q.  note that the format does not expect that start == stop point.
     @param q a closed curve to match to p.  Note that the format does not expect that start == stop point.
     @return the matched intervals.
    */
    public List<Match.Points> match(PairFloatArray p, PairFloatArray q) throws Exception {

        log.info("p.n=" + p.getN() + " q.n=" + q.getN()
            + " useSameNumberOfPoints=" + useSameNumberOfPoints
            + " dp=" + dp);

        if (p.getN() < 2 || q.getN() < 2) {
            throw new IllegalArgumentException("p and q must "
            + " have at least dp*2 points = " + (dp * 2));
        }
        
        if ((p.getN()/dp) < 30 || (q.getN()/dp) < 30) {
            log.warning("WARNING: consider overriding dp to set it to 1 because these are small curves");
        }

        final int n1 = p.getN();
        final int n2 = q.getN();

        boolean interchange = n2 < n1;
        if (interchange) {
            PairFloatArray tmp = p;
            p = q;
            q = p;
        }

        // resample if needed:
        PairFloatArray p2 = null;

        if (useSameNumberOfPoints && n1 != n2) {
            float[][] pxy = new float[2][p.getN()];
            pxy[0] = Arrays.copyOf(p.getX(), p.getN());
            pxy[1] = Arrays.copyOf(p.getY(), p.getN());
            float[][] xyOut = CurveResampler.resample(pxy, n2);
            p2 = new PairFloatArray();
            for (int i = 0; i < xyOut[0].length; ++i) {
                p2.add(xyOut[0][i], xyOut[1][i]);
            }
        }

        PairFloatArray pSub = null;
        PairFloatArray qSub = null;
        if (dp > 1) {
            PairFloatArray tmpP = (p2 != null) ?  p2 : p;
            pSub = new PairFloatArray(tmpP.getN()/dp);
            for (int i = 0; i < tmpP.getN(); i += dp) {
                pSub.add(tmpP.getX(i), tmpP.getY(i));
            }
            qSub = new PairFloatArray(q.getN()/dp);
            for (int i = 0; i < q.getN(); i += dp) {
                qSub.add(q.getX(i), q.getY(i));
            }
        } else {
            pSub = (p2 != null) ?  p2 : p;
            qSub = q;
        }

        if (!overrideMinLength) {
            this.minLength = Math.max(DEFAULT_MIN_LENGTH, (int)Math.round(0.1*Math.max(pSub.getN(), qSub.getN())));
        }

        if (debug) {
            long ts = System.nanoTime();
            log.info(String.format("ts=%d, pSub.n=%d, qSub.n=%d, minLength=%d\n", ts, pSub.getN(), qSub.getN(), minLength));
            plot(pSub, "p_"+ts);
            plot(qSub, "q_"+ts);
        }

        List<Match.Points> pointsList = match0(pSub, qSub);

        if (debug && !pointsList.isEmpty()) {
            long ts = System.nanoTime();
            plotResults(pointsList, pSub, qSub, 4, ts + "_corres_", false);
        }

        // transform results back into original reference frames p and q

        if (dp > 1) {
            for (Match.Points points : pointsList) {
                for (int i = 0; i < points.pIdxs.length; ++i) {
                    points.pIdxs[i] = Math.round(points.pIdxs[i]/(float)dp);
                    points.qIdxs[i] = Math.round(points.qIdxs[i]/(float)dp);
                }
            }
        }

        if (useSameNumberOfPoints && n1 != n2) {
            // multiply by n1 and divide by n2
            float factor = (float)n1/(float)n2;
            if (interchange) {
                factor = 1.f/factor;
            }
            for (int i = 0; i < pointsList.size(); ++i) {
                Match.Points points = pointsList.get(i);
                //pIdxs are w.r.t. p2 reference frame
                for (int j = 0; j < points.pIdxs.length; ++j) {
                    int idx = Math.round(points.pIdxs[j]*factor);
                    //TODO: this could be improved
                    if (!interchange && idx >= n1){
                        idx = n1 - 1;
                    } else if (interchange && idx >= n2){
                        idx = n2 - 1;
                    }
                    points.pIdxs[j] = idx;
                }
            }
        }

        // if there is a sequential mapping of same point in p, drop all but the first.
        for (int i = 0; i < pointsList.size(); ++i) {
            Match.Points points = pointsList.get(i);
            int prev = points.pIdxs[0];
            int idx = 1;
            for (int j = 1; j < points.pIdxs.length; ++j) {
                if (points.pIdxs[j] != prev) {
                    points.pIdxs[idx] = points.pIdxs[j];
                    points.qIdxs[idx] = points.qIdxs[j];
                    ++idx;
                }
                prev = points.pIdxs[j];
            }
            if (idx < points.pIdxs.length) {
                points.pIdxs = Arrays.copyOf(points.pIdxs, idx);
                points.qIdxs = Arrays.copyOf(points.qIdxs, idx);
            }
        }

        if (interchange) {
            for (int i = 0; i < pointsList.size(); ++i) {
                Match.Points points = pointsList.get(i);
                //pIdxs are w.r.t. p2 reference frame
                int[] tmp = points.pIdxs;
                points.pIdxs = points.qIdxs;
                points.qIdxs = tmp;
            }
        }

        return pointsList;
    }

    private List<Match.Points> match0(PairFloatArray p, PairFloatArray q) throws Exception {

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
        assert(p.getN() <= q.getN());
        float[][][] md = createDifferenceMatrices(p, q);
        applySummedAreaTableConversion(md);
        List<Match.Points> points = match0(md, p, q);

        return points;
    }

    private List<Match.Points> match0(float[][][] md, PairFloatArray p, PairFloatArray q) throws Exception {

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
           an example, the scissors opened versus closed in unit tests.
           The occlusion should be handled for this one too.

        For the multi-objective optimization cost,
        need the sum of differences of chords and the fraction
        of the whole to calculate the Salukwdze distance 
        of the Paretto frontier.
        */
        List<Match> matches = findMinima(md, n1, n2);

        List<Match.Points> points = new ArrayList<>();
        for (Match m : matches) {
            points.add(new Match.Points(m));
        }

        return points;
    }

    /**
     *
     * @param md the M X N X N integral sum of chord differences following the paper.
     * @param n1 the number of points in closed curve 1
     * @param n2 the number of points in closed curve 2
     * @return the best matching of curve 1 to curve 2 for each offset chord diff image.
     * The offset diff images are md[offset].  Note that the variable maxChordSum has been updated in each.
     * @throws Exception
     */
    private List<Match> findMinima(float[][][] md, int n1, int n2) throws Exception {

        //q.n is >= p.n, that is n2 >= n1
        
        //md[0:n2-1][0:n1-1][0:n1-1]
        if (md.length != n2 || md[0].length != n1) {
            throw new IllegalArgumentException("error in algorithm arguments:"
                + " revise to set n1 and n2 from md");
        }
        if (n2 < n1) {
            throw new IllegalArgumentException("n2 must be >= n1");
        }        
        
        // reading over a range of window sizes to keep the 
        // sum/nPix below thresh and keeping the mincost solutions.

        // find the intervals of contiguous minima and assign 
        // curve indexes to the largest segments.
        // (note that the objective formula for the cost
        // is the Salukwzde distance).

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
        
        md[0:n2-1][0:n1-1][0:n1-1]
        */

        if (search.equals(SEARCH.FAST)) {
            return fastMinimaSearch(md);
        } else if (search.equals(SEARCH.ALL_BLOCK_SIZES)) {
            return allBlockSizesMinimaSearch(md, n1, n2);
        } else if (search.equals(SEARCH.EXHAUSTIVE)) {
            return exhaustiveMinimaSearch(md, n1, n2);
        } else {
            throw new IllegalStateException(String.format("%s not implemented\n", search));
        }
    }

    /**
     * find the closest matches using a linear block search.
     * the r.t.c. for the search is O(N*M) and the sort is O(M*N*log(topK))
     * so the total r.t.c. is O(M*N*log(topK)).
     * @param md
     * @return the topK best matches
     * @throws Exception
     */
    private List<Match> fastMinimaSearch(float[][][] md) throws Exception {

        List<Match> candidates = new ArrayList<>();

        int M = md.length;
        int N = md[0].length;

        SummedAreaTable st = new SummedAreaTable();
        float[] outC = new float[2];
        float weight;

        // block sizes from minLength to N.
        // dR is the number of block sizes.  setting that to log(n).  consider fewer.
        int nIntervals = Math.max(1, (int)(Math.log(N)/Math.log(2)));
        int dR = Math.max(1, (N - minLength)/nIntervals);
        if (debug) {
            log.info(String.format("dR = %d\n", dR));
        }

        Match curM;
        for (int offset = 0; offset < M; ++offset) {
            for (int iDiag = 0; iDiag < N; ++iDiag) {
                if (N - iDiag < minLength) {
                    break;
                }
                st.extractWindowFromSummedAreaTable(md[offset], iDiag, N - 1, iDiag, N - 1, outC);
                if (outC[1] < 1) {
                    continue;
                }
                weight = outC[0] / outC[1];
                if (weight > thresh) {
                    continue;
                }
                curM = new Match(N);
                curM.add(iDiag, offset, N - iDiag, weight, 0);
                candidates.add(curM);
            }
        }

        if (candidates.isEmpty()) {
            return candidates;
        }

        // set maxChordSum in all
        double maxChordSum = Double.NEGATIVE_INFINITY;
        for (Match m : candidates) {
            maxChordSum = Math.max(maxChordSum, m.diffChordSum);
        }
        for (Match m : candidates) {
            m.maxChordSum = maxChordSum;
        }

        // heap w/ restricted size to reduce comparisons to M*N*log(topK)
        TreeSet<Match> tree = new TreeSet<>();
        for (Match m : candidates){
            tree.add(m);
            if (tree.size() > topK) {
                tree.pollLast();
            }
        }
        candidates = new ArrayList<>();
        for (Match m : tree) {
            candidates.add(m);
        }

        return candidates;
    }

    private List<Match> allBlockSizesMinimaSearch(float[][][] md, int n1, int n2) throws Exception {
        List<Match> candidates = new ArrayList<>();

        int M = md.length;
        int N = md[0].length;

        SummedAreaTable st = new SummedAreaTable();
        float[] outC = new float[2];
        float weight;

        // block sizes from minLength to N.
        // dR is the number of block sizes.  setting that to log(n).  consider fewer.
        int nIntervals = Math.max(1, (int)(Math.log(N)/Math.log(2)));
        int dR = Math.max(1, (N - minLength)/nIntervals);
        System.out.printf("df = %d\n", dR);

        Match curM;
        for (int offset = 0; offset < M; ++offset) {
            for (int iDiag = 0; iDiag < N; ++iDiag) {
                for (int r = minLength; r < N; r += dR) {
                    if (iDiag + r - 1 >= N) {
                        break;
                    }
                    st.extractWindowFromSummedAreaTable(md[offset], iDiag, iDiag + r - 1, iDiag, iDiag + r - 1, outC);
                    if (outC[1] < 1) {
                        continue;
                    }
                    weight = outC[0] / outC[1];
                    if (weight > thresh) {
                        continue;
                    }
                    curM = new Match(N);
                    curM.add(iDiag, offset, r, weight, 0);
                    candidates.add(curM);
                }
            }
        }

        if (candidates.isEmpty()) {
            return candidates;
        }

        // set maxChordSum in all
        double maxChordSum = Double.NEGATIVE_INFINITY;
        for (Match m : candidates) {
            maxChordSum = Math.max(maxChordSum, m.diffChordSum);
        }
        for (Match m : candidates) {
            m.maxChordSum = maxChordSum;
        }

        TreeSet<Match> tree = new TreeSet<>();
        for (Match m : candidates){
            tree.add(m);
            if (tree.size() > topK) {
                tree.pollLast();
            }
        }
        candidates = new ArrayList<>();
        for (Match m : tree) {
            candidates.add(m);
        }

        return candidates;
    }

    private List<Match> exhaustiveMinimaSearch(float[][][] md, int n1, int n2) throws Exception {

        // starting another approach.
        // for each position on A1 diagonal i,
        //    make combinations of blocks with same position on next offset images
        //    etc, sequentially matching or skipping and moving onto next offset image.
        // a recursion can fill all of these combinations.
        // then can use salukwdze comparator to find best among them.
        int N = md[0].length;
        double nIntervals = Math.log(N)/Math.log(2);
        nIntervals = 1;
        int dR = (int)Math.ceil((float)(N - minLength)/nIntervals);
        if (dR == 0) {
            dR = 1;
        }

        int rUpper = N;
        int nR = ((N - minLength)/dR) + 1;
        List<Match> candidates = new ArrayList<>();

        System.out.printf("\nnR=%d, dR=%d\n", nR, dR);

        // n*m* log(n)
        recursion(md, 0, rUpper, minLength, dR, candidates, new Match(N),
                new SummedAreaTable(), new float[2], 0);

        if (debug) {
            Set<Match> unique = new HashSet<>(candidates);
            System.out.printf("n=%d, nR=%d, nCandidates=%d, nUnique=%d\n", N, nR, candidates.size(), unique.size());
        }

        if (candidates.isEmpty()) {
            return candidates;
        }

        // set maxChordSum in all
        double maxChordSum = Double.NEGATIVE_INFINITY;
        for (Match m : candidates) {
            maxChordSum = Math.max(maxChordSum, m.diffChordSum);
        }
        for (Match m : candidates) {
            m.maxChordSum = maxChordSum;
        }

        TreeSet<Match> tree = new TreeSet<>();
        for (Match m : candidates){
            tree.add(m);
            if (tree.size() > topK) {
                tree.pollLast();
            }
        }
        candidates = new ArrayList<>();
        for (Match m : tree) {
            candidates.add(m);
        }

        return candidates;
    }

    private void recursion(float[][][] md, int iDiag, int r, final int rMin, final int dR, List<Match> out,
        Match curM, final SummedAreaTable st, final float[] outC, int rangeNum) {

        int n = md[0].length;
        // base case:  end of the diagonal.
        // store current match
        // reduce block size and start from top of diagonal again
        if (iDiag + r > n) {
            // end state.  store copy, and reset state if possible
            if (curM.mLen > 0) {
                out.add(curM.copy());
            }
            r -= dR;
            if (r < rMin) {
                return;
            }
            iDiag = 0;
            curM = new Match(curM.N);
            rangeNum = 0;
        }

        // uses implicit backtracking to save memory

        // skip i.
        // there will be more than one curM that ends as all skipped
        recursion(md, iDiag + 1, r, rMin, dR, out, curM, st, outC, rangeNum);

        float weight;
        for (int offset = 0; offset < md.length; ++offset) {
            st.extractWindowFromSummedAreaTable(md[offset], iDiag, iDiag + r-1, iDiag, iDiag + r-1, outC);
            if (outC[1] < 1) {
                continue;
            }
            weight = outC[0]/outC[1];
            if (weight > thresh) {
                continue;
            }
            //consider that we should not go to the next point on the diagonal which is iDiag+r because
            //    it being adjacent means that a larger block starting at higher index should have covered
            //        this case.  in other words, by skipping 1, we avoid adjacent blocks on the diagonal
            //        that were already covered by previous larger block match.
            curM.add(iDiag, offset, r, weight, rangeNum);
            recursion(md, iDiag + r + 1, r, rMin, dR, out, curM, st, outC, rangeNum + 1);
        }
    }

    /**
     * create the matrices of differences between p
     * and q.  Note that the matrix differences are
     * absolute differences.
     * index0 is rotations of q,  index1 is p.n, index2 is q.n
      returns a[0:q.n-1][0:p.n-1][0:p.n-1]
    */
    protected float[][][] createDifferenceMatrices(PairFloatArray p, PairFloatArray q) {

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
         make difference matrices.
            there will be N A_2 matrices in which each
            is shifted left and up by 1 (or some other value).

            M_D^n = A_1(1:M,1:M) - A_2(n:n+M-1,n:n+M-1)
                shifting A_2 by 0 through N covering all
                orientation angles.
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
        
        //TODO: look into Toeplitz matrix and cyclic matrix

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
       
        //print("differences:", md);

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

     r.t.c. is O(n^2)
    */
    public static float[][] createDescriptorMatrix(PairFloatArray p, int n) {
        
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

                //double angleA = LinesAndAngles.calcAngle(
                double angleA = LinesAndAngles.calcClockwiseAngle(
                    p.getX(i1), p.getY(i1),
                    p.getX(i2), p.getY(i2),
                    p.getX(imid), p.getY(imid)
                );
                
                if (Double.isNaN(angleA)) {
                    if (i2 < i1 && i1 < imid) {
                        angleA = LinesAndAngles.calcAngle(
                        //angleA = LinesAndAngles.calcClockwiseAngle(
                            p.getX(i2), p.getY(i2),
                            p.getX(i1), p.getY(i1),
                            p.getX(imid), p.getY(imid));
                        if (Double.isNaN(angleA)) {
                            angleA = LinesAndAngles.calcAngle(
                            //angleA = LinesAndAngles.calcClockwiseAngle(
                                p.getX(i2), p.getY(i2),
                                p.getX(i1), p.getY(i1),
                                p.getX(imid), p.getY(imid));
                        }
                    } else {
                       // System.out.println(
                       // "SKIP i1=" + i1 + " imid=" + imid + " i2=" + i2);
                        continue;
                    }
                }

                float angleAF = (float)angleA;

                /*
                String str = String.format(
                    "[%d](%d,%d) [%d](%d,%d) [%d](%d,%d) a=%.4f",
                    i1, p.getX(i1), p.getY(i1),
                    i2, p.getX(i2), p.getY(i2),
                    imid, p.getX(imid), p.getY(imid),
                    (float) angleA * 180. / Math.PI);
                log.fine(str);
                */

                if (angleAF >= TWO_PI) {// rounding errors from int to float and double
                    angleAF -= TWO_PI;
                }
                if (Math.abs(angleAF - TWO_PI) < ZERO_TOL) {
                    angleAF = 0.f;
                }

                a[i1][i2] = angleAF;
                
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

    private float[][] copy(float[][] a) {
        float[][] a2 = new float[a.length][];
        for (int i = 0; i < a2.length; ++i) {
            a2[i] = Arrays.copyOf(a[i], a[i].length);
        }
        return a2;
    }

    protected void rotate(float[][] prevShifted) {

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
                float v = Math.abs(a1[i][j] - a2[i][j]);
                if (Math.abs(TWO_PI - v) < ZERO_TOL) {
                    v = 0.f;
                }
                output[i][j] = v;
            }
        }

        return output;
    }
    
    protected void applySummedAreaTableConversion(float[][][] md) {
        for (float[][] floats : md) {
            applySummedAreaTableConversion(floats);
        }
    }

    protected void applySummedAreaTableConversion(float[][] mdI) {

        int w = mdI.length;
        int h = mdI[0].length;
        
        // sum along columns, that is a[i][*]
        for (int i = 0; i < w; ++i) {
            for (int j = 0; j < h; ++j) {
                if (i > 0 && j > 0) {
                    mdI[i][j] += (mdI[i - 1][j] + mdI[i][j - 1]
                        - mdI[i - 1][j - 1]);
                } else if (i > 0) {
                    mdI[i][j] += mdI[i - 1][j];
                } else if (j > 0) {
                    mdI[i][j] += mdI[i][j - 1];
                }
            }
        }
    }

    private void print(String label, float[][] a) {

        StringBuilder sb = new StringBuilder(label);
        sb.append("\n");

        for (int i = 0; i < a.length; ++i) {
            sb.append(String.format("row %3d: ", i));
            for (int j = 0; j < a[i].length; ++j) {
                sb.append(String.format(" %.2f,", a[i][j]));
            }
            log.fine(sb.toString());
            System.out.println(sb.toString());
            sb.delete(0, sb.length());
        }
    }

    private String plot(PairFloatArray p, String fileLabel) throws Exception {

        float[] x = Arrays.copyOf(p.getX(), p.getN());
        float[] y = Arrays.copyOf(p.getY(), p.getN());
        float xMax = MiscMath.findMax(x) + 1;
        float yMax = MiscMath.findMax(y) + 1;

        PolygonAndPointPlotter plot = new PolygonAndPointPlotter();

        plot.addPlot(0, xMax, 0, yMax,
                x, y, x, y, "");

        return plot.writeFile(fileLabel);
    }

    private List<String> plotResults(List<Match.Points> results, PairFloatArray p, PairFloatArray q,
                                     int spacing, String fileSuffix, boolean printIndexes) {
        List<String> writtenFiles = new ArrayList<>();

        PairIntArray pInt = new PairIntArray();
        PairIntArray qInt = new PairIntArray();
        for (int i = 0; i < p.getN(); ++i) {
            pInt.add(Math.round(p.getX(i)), Math.round(p.getY(i)));
        }
        for (int i = 0; i < q.getN(); ++i) {
            qInt.add(Math.round(q.getX(i)), Math.round(q.getY(i)));
        }

        try {
            for (int i = 0; i < results.size(); ++i) {
                CorrespondencePlotter plotter = new CorrespondencePlotter(pInt, qInt);
                Match.Points result = results.get(i);
                for (int ii = 0; ii < result.pIdxs.length; ii += spacing) {
                    int idx1 = result.pIdxs[ii];
                    int idx2 = result.qIdxs[ii];
                    int x1 = pInt.getX(idx1);
                    int y1 = pInt.getY(idx1);
                    int x2 = qInt.getX(idx2);
                    int y2 = qInt.getY(idx2);

                    if (printIndexes) {
                        System.out.println(String.format(
                                "(%d, %d) <=> (%d, %d)", x1, y1, x2, y2));
                    }

                    plotter.drawLineInAlternatingColors(x1, y1, x2, y2, 0);
                }
                String filePath = plotter.writeImage(String.format("_%s_%d", fileSuffix, i));
                writtenFiles.add(filePath);
            }
        } catch (IOException ex) {
            log.severe(ex.getMessage());
        }
        return writtenFiles;
    }
}
