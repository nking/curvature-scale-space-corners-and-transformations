package algorithms.imageProcessing.matching;

import algorithms.compGeometry.LinesAndAngles;
import algorithms.misc.Histogram;
import algorithms.misc.HistogramHolder;
import algorithms.misc.MiscMath;
import algorithms.util.Errors;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.TreeMap;
import java.util.logging.Logger;

/**
 NOTE: NOT READY FOR USE YET.
 TODO: need to change read pattern of
 difference sat, and optimization.

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
             - can handle similarity transforms, but not occlusion
         (b) OR, equidistant points
             - can handle occlusion, but not scale
       ** equidistant is used here.

       NOTE: changes will be made soon to accomodate
             search of remaining points when there are
             unequal number of points.

       The runtime complexity for building the integral
       image is O(m*n) where n and m are the number of sampled
       points on the input shapes.

       The runtime complexity for the search of the
       integral image of summed differences and analysis
       will be added here:
       *
 * @author nichole
 */
public class PartialShapeMatcher {

    /**
     * in sampling the boundaries of the shapes, one can
     * choose to use the same number for each (which can result
     * in very different spacings for different sized curves)
     * or one can choose a set distance between sampling
     * points.
     * dp is the set distance between sampling points.
       The authors use 3 as an example.
     */
    protected int dp = 5;

    protected Logger log = Logger.getLogger(this.getClass().getName());

    public void overrideSamplingDistance(int d) {
        this.dp = d;
    }

    /**
     * NOT READY FOR USE.

      A shape is defined as the clockwise ordered sequence
      of points P_1...P_N
      and the shape to match has points Q_1...Q_N.
      The spacings used within this method are equidistant
      and the default is 5, so override that if a different number
      is needed.
      The fixed equidistant spacing is invariant to rotation
      and translation, but not to scale, so if the user needs to solve
      for scale, need to do so outside of this method, that is, apply
      scale changes to the datasets before use of this method..

     NOTE: You may want to pre-process the shape points by using
     PairIntArray p = ImageProcessor.extractSmoothedOrderedBoundary

     @param p
     @param q
    */
    public Sequences match(PairIntArray p, PairIntArray q) {

        log.info("p.n=" + p.getN() + " q.n=" + q.getN());

        int diffN = p.getN() - q.getN();

        // --- make difference matrices ---

        //md[0:n2-1][0:n1-1][0:n1-1]
        float[][][] md;
        int n1, n2;
        if (diffN <= 0) {
            n1 = p.getN();
            n2 = q.getN();
            md = createDifferenceMatrices(p, q);
        } else {
            n1 = q.getN();
            n2 = p.getN();
            md = createDifferenceMatrices(q, p);
        }

        log.info("p=" + p.toString());
        log.info("q=" + q.toString());

        /*
        the matrices in md can be analyzed for best
        global solution and/or separately for best local
        solution.

        This method will return results for a local
        solution to create the point correspondence list.

        Note that the local best could be two different
        kinds of models, so might write two
        different methods for the results.
        (1) the assumption of same object but with some
            amount of occlusion, hence gaps in correspondence.
        (2) the assumption of same object but with
           some parts being differently oriented, for
           an example, the scissors opened versus closed.
        */

        /*
        need sum of differences in sequence and the fraction
        of the whole.
        paretto efficiency is that all are at a compromise of best state,
        such that increasing the state of one did not worsen the
        state of another.
        prefer:
           -- smaller total difference for largest fraction of whole
           -- 2ndly, largest total coverage

        Note that for the rigid model (excepting scale transformation)
        one would want to maximize the 3nd point, coverage, first
        with a consistent transformation.

        The articulated model chooses the 2nd point, second to get
        best fits of components first.
        */

        /*
        NOTE: this may change with more testing.

        goal is to find the best chains of sequential
        matches below a threshold and then use
        multi-objective optimization to choose the
        best consistent aggregation of chains as the
        final correspondence list.

        the rigid model allowing occlusion
        (not yet implemented) will
        likely be a better solution and is similar to
        matching patterns elsewhere in this project.

        focusing first on this articulated match to
        look at the range of ability to identify a
        whole object which may be occluded and which
        may have parts which have separate rigid
        rotation (such as the scissors opened bersus
        closed) and which may include extended objects due
        to the segmentation including a part of
        the background.

        caveats to the articulated match are that
        greedy solutions built by fraction of whole
        may have many different kinds of errors,
        but composing sequences with the top k
        fraction (at every aggregation of sequential
        segments) quickly leads to an unfeasibly
        large number of sequences to evaluate.

        */

        /*
        pattern:
            a. set rMax = n1/3 and use extractSimilat.
            n. create Sequences from best of those results
               without further use of matchArticulated
            c. repeat a and b for rMax = n1/4
            d. compare results and keep best
            e. if optional (unimplemented) performEval is
               set, create a projection to match the remaining
               points which belong to the shape.
        */

        boolean performEval = false;

        int[] rs = new int[]{n1/3, n1/4};

        List<Sequences> results = new ArrayList<Sequences>();

        for (int r : rs) {

            if (r < 2) {
                // r=1 reads diagonal of 0's only
                r = 2;
            }

            List<Sequences> seqsList = new ArrayList<Sequences>();
            // build the matching sequential sequences by
            // searching md in given block size r
            extractSimilar(md, seqsList, r);

            if (!seqsList.isEmpty()) {
                Sequences result = findBest(seqsList);
                results.add(result);
            }
        }

        //TODO: an optional, but default evaluation
        // stage to match the remaining points with
        // a transformation model (and ransac)
        // that would include those separated by
        // occlusion and exclude those not part
        // of the p shape.
        Sequences best;
        if (performEval) {
            // solve for transformation, add points, and
            // return best solution.
            if (diffN <= 0) {
                best = transformAndEvaluate(results, p, q, md);
            } else {
                best = transformAndEvaluate(results, q, p, md);
            }
        } else {
            // compare fraction and sum of differences
            best = findBest(results);
        }

        if (diffN <= 0) {
            return best;
        }

        transpose(best, n1, n2);

        return best;
    }

    protected void extractSimilar(float[][][] md,
        List<Sequences> sequences, int r) {

        //md[0:n2-1][0:n1-1][0:n1-1]

        int n2 = md.length;
        int n1 = md[0].length;

        // 23 degrees is 0.4014
        double thresh = (Math.PI/180.) * 10.;//*23.;

        MinDiffs mins = new MinDiffs(n1);
        findMinDifferenceMatrix(md, r, thresh, mins);

        // 10 degrees is 0.175
        double tolerance = 0.25;//0.1;//0.25;

        DiffMatrixResults equivBest = new DiffMatrixResults(n1);
        findEquivalentBest(md, r, mins, thresh, tolerance,
            n1, n2, equivBest);

        printEquivBest(equivBest);

        /*
        DiffMatrixResults
            map0 key = offset
                 value = map with key = start i,j
                        value = Sequence i:j, w/ offset and length
        */

        ResultSequences rs = new ResultSequences(n1);

        // ---- merge adjacent and overlapping Sequence w/ same offset ---
        TIntObjectIterator<Map<PairInt, Sequence>> iter0 =
            equivBest.map0.iterator();
        for (int i0 = 0; i0 < equivBest.map0.size(); ++i0) {
            iter0.advance();

            int offset = iter0.key();

            // all sequences in this map have the same offset
            Map<PairInt, Sequence> map10 = iter0.value();
            List<Sequence> seqs = new ArrayList<Sequence>(map10.values());
            Sequence.mergeSequences(seqs);

            Collections.sort(seqs, new SequenceComparator4());

            //changed to form adjacent segments where wrap
            // around is present, so assert format here
            assert(assertNoWrapAround(seqs));

            Sequence s0 = seqs.get(0);
            Sequence s1 = seqs.get(seqs.size() - 1);

            assert(s0.getOffset() == offset);
            assert(s1.getOffset() == offset);

            // --- look for whether all points are matched ----
            if (s0.isStartSentinel() && s1.isStopSentinel()) {
                Sequences sWrap = new Sequences();
                sWrap.getSequences().addAll(seqs);
                rs.add(offset, sWrap);
                continue;
            }

            // note that the sequences have been sorted by startIdx1

            // --- create a new Sequences when there's a gap
            Sequences current = new Sequences();
            current.getSequences().add(s0);
            for (int i = 1; i < seqs.size(); ++i) {

                s0 = current.getSequences().get(
                    current.getSequences().size() - 1);
                int s0StopIdx1 = s0.getStopIdx1();

                Sequence s = seqs.get(i);
                assert(s.getOffset() == offset);

                if ((s.startIdx1 - s0StopIdx1) < 2) {
                    current.getSequences().add(s0);
                } else {
                    rs.add(offset, current);
                    current = new Sequences();
                    current.getSequences().add(s);
                }
            }
            rs.add(offset, current);
        }

        equivBest = null;

        int nIter = 0;
        boolean didAppend = false;
        do {

            didAppend = false;

            // ----- append adjacent sequences with nearly same offsets ----
            TIntObjectIterator<Map<PairInt, Sequences>> iter1 =
                rs.map0.iterator();
            for (int i = 0; i < rs.map0.size(); ++i) {
                iter1.advance();

                int offset = iter1.key();
                Map<PairInt, Sequences> map0 = iter1.value();
                Iterator<Entry<PairInt, Sequences>> iter2 =
                    map0.entrySet().iterator();

                //TODO: could make the offset before and after 
                // depend upon dp, in other words, a large dp
                // would need a large gap search.
                // if that change is ever made, one would need to
                // also need to make a sequence.precedes() and 
                // sequence.proceeds() that accepts a gap argument.
                
                int offset1 = offset + 1;
                int offsetNMinus1 = offset - 1;
                if (offset1 >= n2) {
                    offset1 -= n2;
                }
                if (offsetNMinus1 >= n2) {
                    offsetNMinus1 -=  n2;
                }
                
                int[] sOffsets = new int[]{offset1, offsetNMinus1};

                while (iter2.hasNext()) {

                    Entry<PairInt, Sequences> entry = iter2.next();

                    Sequences current = entry.getValue();

                    if (current.getSequences().isEmpty()) {
                        continue;
                    }

                    Sequence s0 = current.getSequences().get(0);
                    Sequence s1 = current.getSequences().get(
                        current.getSequences().size() - 1);

                    for (int sOffset : sOffsets) {

                        // -- looking for offset=1 adjacent sequences ---
                        Map<PairInt, Sequences> mapOff =
                            rs.map0.get(sOffset);
                        if (mapOff != null) {
                            Iterator<Entry<PairInt, Sequences>> iter10 =
                                mapOff.entrySet().iterator();
                            while (iter10.hasNext()) {
                                Entry<PairInt, Sequences> entry2 =
                                    iter10.next();
                                Sequences s2 = entry2.getValue();
                                if (s2.getSequences().isEmpty()) {
                                    continue;
                                }
                                Sequence st0 = s2.getSequences().get(0);
                                Sequence st1 = s2.getSequences().get(
                                    s2.getSequences().size() - 1);

                                // look for a sequences that precedes s0
                                boolean precedes = s0.precedes(st1);
                                if (precedes) {
                                    current.getSequences().addAll(0,
                                        s2.getSequences());
                                    //clear s2 in map0 and map1
                                    s2.getSequences().clear();
                                    s0 = st0;
                                    didAppend = true;
                                } else {
                                    // look for a sequences that follows s1
                                    boolean proceeds = s1.proceeds(st0);
                                    if (proceeds) {
                                        current.getSequences().addAll(
                                            s2.getSequences());
                                        //clear s2 in map0 and map1
                                        s2.getSequences().clear();
                                        s1 = st1;
                                        didAppend = true;
                                    }
                                }
                            }
                        }
                    }
                }
            }

        } while (didAppend);

        // ----- calc Sequences stats ----
        TIntObjectIterator<Map<PairInt, Sequences>> iter1 =
            rs.map0.iterator();
        for (int i = 0; i < rs.map0.size(); ++i) {
            iter1.advance();

            int offset = iter1.key();
            Map<PairInt, Sequences> map0 = iter1.value();
            sequences.addAll(map0.values());
        }

        populateStats(sequences);

        log.info(sequences.size() + " sequences");
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

        // TODO: look at histograms of angles.
        /*
        float binSz = (float)(Math.PI/8.f);
        HistogramHolder hist1 = createHistogram(a1, binSz);

        HistogramHolder hist2 = createHistogram(a2, binSz);

        try {
            hist1.plotHistogram("a1 hist", "a1_hist");
            hist2.plotHistogram("a2 hist", "a2_hist");
        } catch(Throwable t) {
        }
        */

        /*
        - find rxr sized blocks similar to one another
          by starting at main diagonal element
          A_1(s,s) and A_2(m,m)
          which have a small angular difference value

                         1
          D_a(s,m,r) = ---- * summation_i_0_to_(r-1)
                        r^2
                            * summation_j_0_(r-1)
                               of [A_1(s+i,s+j) - A_2(m+i,m+j)]^2

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
                shifting A_2 by 0 through n where n is (N-M+1?),
                but shifting it by N instead would cover all
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
            applySummedAreaTableConversion(md[i]);
        }

        //printDiagonal(md[0], mdCoords[0]);

        log.fine("md.length=" + md.length);

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
    protected float[][] createDescriptorMatrix(PairIntArray p,
        int n) {

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
            int start = i1 + 1 + dp;
            for (int ii = start; ii < (start + n - 1 - dp); ++ii) {
                int i2 = ii;

                int imid = i2 - dp;
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

    protected void applySummedAreaTableConversion(float[][] mdI) {

        for (int x = 0; x < mdI.length; ++x) {
            for (int y = 0; y < mdI[x].length; ++y) {

                if (x > 0 && y > 0) {
                    mdI[x][y] += (mdI[x - 1][y] + mdI[x][y - 1]
                        - mdI[x - 1][y - 1]);
                } else if (x > 0) {
                    mdI[x][y] += mdI[x - 1][y];
                } else if (y > 0) {
                    mdI[x][y] += mdI[x][y - 1];
                }
            }
        }
    }

    private void filterForConsistentClockwise(
        List<Sequences> tracks) {

        TIntList rmList = new TIntArrayList();

        for (int i = 0; i < tracks.size(); ++i) {

            Sequences sequences = tracks.get(i);

            if (!sequences.isConsistentClockwise()) {
                rmList.add(i);
            }
        }

        log.info("removing " + rmList.size()
            + " tracks from " + tracks.size());

        for (int i = (rmList.size() - 1); i > -1; --i) {
            int rmIdx = rmList.get(i);
            tracks.remove(rmIdx);
        }

    }

    private void transpose(Sequences sequences,
        int n1, int n2) {

        if (sequences == null) {
            return;
        }

        float f = sequences.getFractionOfWhole();
        f *= ((float)n1/(float)n2);
        sequences.setFractionOfWhole(f);

        List<Sequence> sqs = sequences.getSequences();
        List<Sequence> tr = new ArrayList<Sequence>();
        for (Sequence s : sqs) {
            Sequence str = s.transpose();
            tr.add(str);
        }
        sqs.clear();
        sqs.addAll(tr);
    }

    protected int[] findOffsetHistPeaks(float[] values,
        int n1, int n2) {

        if (values == null || values.length == 0) {
            return null;
        }

        float binWidth = 6.f;
        //int nBins = 10;
        HistogramHolder hist =
            Histogram.createSimpleHistogram(
                binWidth,
                //nBins,
                values,
                Errors.populateYErrorsBySqrt(values));

        // to help wrap around, moving items in last bin
        // into first bin.  then offset of 0 +- binWidth/2
        // should be present there
        int limit = n2 - (int)Math.ceil((float)binWidth/2.f);
        for (int i = (hist.getYHist().length - 1); i > -1; --i) {
            float x = hist.getXHist()[i];
            if (x < limit) {
                break;
            }
            hist.getYHist()[0] += hist.getYHist()[i];
            hist.getYHist()[i] = 0;
        }

        List<Integer> indexes =
            MiscMath.findStrongPeakIndexesDescSort(
            hist, 0.1f);

        // also looking at last non zero
        int idx = MiscMath.findLastNonZeroIndex(hist);
        if (idx > -1) {
            log.info("last non zero offset=" +
                hist.getXHist()[idx]);
        }

        int[] offsets;
        if (indexes != null && indexes.size() > 0) {
            int end = (indexes.size() < 3) ? indexes.size() : 3;
            offsets = new int[end];
            for (int j = 0; j < end; ++j) {
                offsets[j] = Math.round(hist.getXHist()
                    [indexes.get(j).intValue()]);
            }
        } else {
            int yPeakIdx = MiscMath.findYMaxIndex(
                hist.getYHistFloat());
            if (yPeakIdx == -1) {
                return null;
            }
            offsets = new int[]{
                Math.round(hist.getXHist()[yPeakIdx])
            };
        }

        try {
            hist.plotHistogram("offsets", "offsets");
        } catch (Throwable t) {

        }

        return offsets;
    }

    protected int[] findOffsets(List<Sequence> sequences,
        int n1, int n2) {

        int n = 0;
        for (int i = 0; i < sequences.size(); ++i) {
            n += sequences.get(i).length();
        }

        // histogram with bins 10 degrees in size
        float[] values = new float[n];
        n = 0;
        for (int i = 0; i < sequences.size(); ++i) {
            int offset = sequences.get(i).getOffset();
            int len = sequences.get(i).length();
            for (int j = 0; j < len; ++j) {
                values[n] = offset;
                n++;
            }
        }

        int[] offsets = findOffsetHistPeaks(values, n1, n2);

        return offsets;
    }

    private void print(String prefix, List<Sequence> sequences) {

        for (int i = 0; i < sequences.size(); ++i) {
            log.info(String.format("%d %s %s", i, prefix,
                sequences.get(i)));
        }
    }

    private void print0(List<Sequence> sequences,
        String label) {

        // --- sort a copy of sequences by fraction, diff, then
        // startIdx1 and print

        List<Sequence> copy = new ArrayList<Sequence>(sequences);

        // desc sort by fraction
        Collections.sort(copy, new SequenceComparator2());

        float maxDiff = Float.MIN_VALUE;

        for (int i = 0; i < copy.size(); ++i) {

            log.info(String.format("%s FSORT %d  %s", label, i, copy.get(i).toString()));

            if (copy.get(i).absAvgSumDiffs > maxDiff) {
                maxDiff = copy.get(i).absAvgSumDiffs;
            }
        }

        Collections.sort(copy, new SequenceComparator3(maxDiff));
        for (int i = 0; i < copy.size(); ++i) {
            log.info(String.format("%s FDSORT %d  %s", label, i,
                copy.get(i).toString()));
        }
    }

    private HistogramHolder createHistogram(float[][] a,
        float binSz) {

        float[] values = new float[a.length * a.length -
            a.length];

        int count = 0;
        for (int i = 0; i < a.length; ++i) {
            for (int j = 0; j < a[0].length; ++j) {
                if (i == j) {
                    continue;
                }
                values[count] = Math.abs(a[i][j]);
                count++;
            }
        }

        HistogramHolder hist =
            Histogram.createSimpleHistogram(binSz, values,
                Errors.populateYErrorsBySqrt(values));

        return hist;
    }

    private void printEquivBest(DiffMatrixResults equivBest) {

        StringBuilder sb = new StringBuilder();

        TIntObjectIterator<Map<PairInt, Sequence>> iter0 =
            equivBest.map0.iterator();

        for (int i0 = 0; i0 < equivBest.map0.size(); ++i0) {

            iter0.advance();

            int offset = iter0.key();
            Map<PairInt, Sequence> map10 = iter0.value();

            Iterator<Entry<PairInt, Sequence>> iter1 = map10.entrySet().iterator();

            while (iter1.hasNext()) {

                Entry<PairInt, Sequence> entry = iter1.next();

                Sequence s = entry.getValue();

                sb.append(String.format(
                    "equivBest: %s\n", s.toString()));
            }

            if (sb.length() > 0) {
                log.info(sb.toString());
                sb.delete(0, sb.length());
            }
        }
    }

    private Sequences transformAndEvaluate(List<Sequences> results,
        PairIntArray p, PairIntArray q, float[][][] md) {

        //NOTE: md is passed in to be able to add stats of
        // points found after transformation

        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    private Sequences findBest(List<Sequences> results) {

        // keep the result with higher fraction.

        // NOTE: may need to revise to include abs difference

        Sequences best = null;
        float bestFraction = Float.MIN_VALUE;
        for (Sequences result : results) {
            float frac = result.getFractionOfWhole();
            log.info("result: " + result.toString());
            if (frac > bestFraction) {
                best = result;
                bestFraction = frac;
            }
        }

        return best;
    }

    private void populateStats(List<Sequences> seedTracks) {
        // calculate the stats for each track (== Sequences)
        for (Sequences track : seedTracks) {
            int sumLen = 0;
            float sumFrac = 0;
            double sumDiffs = 0;
            for (Sequence s : track.getSequences()) {
                int len = s.length();
                float diff = s.absAvgSumDiffs * len;
                sumLen += len;
                sumDiffs += diff;
                sumFrac += s.fractionOfWhole;
            }
            track.setAbsSumDiffs(sumDiffs);
            track.setAvgSumDiffs((float)(sumDiffs/(float)sumLen));
            track.setFractionOfWhole(sumFrac);
        }
    }

    private class DiffMatrixResults {

        /**
        key = offset
        value = map with key = pairint i, j (start of sequence),
                         value = Sequence i:j, w/ offset and length
        Note that when more than one value for keys are present,
        log it to debug and keep the longest fraction
        though, may want to consider smallest avg diff
        */
        private final TIntObjectMap<Map<PairInt, Sequence>> map0;

        public DiffMatrixResults(int n) {
            map0 = new TIntObjectHashMap<Map<PairInt, Sequence>>(n);
        }

        public void add(int i, int j, int jOffset,
            float diff, int blockSize, int n1, int n2) {

            Map<PairInt, Sequence> map2 = map0.get(jOffset);

            if (map2 == null) {
                map2 = new HashMap<PairInt, Sequence>();
                map0.put(jOffset, map2);
            }

            /*
            sequence is i-blockSize through i
            */
            int startIdx2 = j - blockSize;
            if (startIdx2 < 0) {
                startIdx2 += n2;
            }

            Sequence s = new Sequence(n1, n2, jOffset);
            s.startIdx1 = i - blockSize;
            s.startIdx2 = startIdx2;
            s.stopIdx2 = j;
            s.absAvgSumDiffs = diff;
            s.fractionOfWhole = (float)s.length()/(float)n1;

            PairInt key2 = new PairInt(s.startIdx1, s.startIdx2);
            PairInt key3 = new PairInt(i, s.stopIdx2);

            if (map2.containsKey(key2)) {
                Sequence existing = map2.get(key2);
                if (existing.absAvgSumDiffs > diff) {
                    log.info("discarding s=" + existing);
                    map2.put(key2, s);
                }
            } else {
                map2.put(key2, s);
            }
        }
    }

    /**
     * similar to DiffMatrixResults, except holds aggregated
     * contents of DiffMatrixResults as Sequences, stored
     * by key = offset.
     */
    private class ResultSequences {

        /**
        key = offset
        value = map with key = pairint i, j (start of Sequences),
                         value = Sequences
        Note that when more than one value for keys are present,
        log it to debug and keep the longest fraction
        though, may want to consider smallest avg diff
        */
        private final TIntObjectMap<Map<PairInt, Sequences>> map0;

        /**
         * same as map0 but with pairint i, j being stop of sequence
        */
        private final TIntObjectMap<Map<PairInt, Sequences>> map1;

        public ResultSequences(int n) {
            map0 = new TIntObjectHashMap<Map<PairInt, Sequences>>(n);
            map1 = new TIntObjectHashMap<Map<PairInt, Sequences>>(n);
        }

        public void add(int offset, Sequences s) {

            Sequence s2 = s.getSequences().get(0);

            assert(s2.getOffset() == offset);

            Map<PairInt, Sequences> map2 = map0.get(offset);
            Map<PairInt, Sequences> map3 = map1.get(offset);

            if (map2 == null) {
                map2 = new HashMap<PairInt, Sequences>();
                map0.put(offset, map2);
                map3 = new HashMap<PairInt, Sequences>();
                map1.put(offset, map3);
            }

            PairInt key2 = new PairInt(s2.startIdx1, s2.startIdx2);

            Sequence s3 = s.getSequences().get(s.getSequences().size() - 1);
            assert(s2.getN1() == s3.getN1());
            assert(s2.getN2() == s3.getN2());
            assert(s2.getOffset() == s3.getOffset());
            int n1 = s3.getN1();
            int stopIdx1 = s3.getStopIdx1();
            if (stopIdx1 >= n1) {
                stopIdx1 -= n1;
            }
            PairInt key3 = new PairInt(stopIdx1, s3.stopIdx2);

            if (map2.containsKey(key2)) {
                Sequences existing = map2.get(key2);
                log.warning("discarding SEQUENCES=" + existing);
            }
            map2.put(key2, s);

            if (map3.containsKey(key3)) {
                Sequences existing = map3.get(key3);
                log.warning("discarding SEQUENCES=" + existing);
            }
            map3.put(key3, s);
        }
    }

    private class MinDiffs {
        // first dimension index: md[*][][]
        int[] idxs0;
        // i is the index of this array and represents the index
        //       of point in p array
        // j is i + idxs[0] and represents the index
        //       of point in q array

        float[] mins;
        public MinDiffs(int n) {
            idxs0 = new int[n];
            mins = new float[n];
            Arrays.fill(idxs0, -1);
            Arrays.fill(mins, Float.MAX_VALUE);
        }
    }

    private class TrackComparator implements
        Comparator<Sequences> {

        final int maxNPoints;

        public TrackComparator(int n1) {
            this.maxNPoints = n1;
        }

        @Override
        public int compare(Sequences o1, Sequences o2) {

            // prefer high fraction of whole, then diff

            if (o1.getFractionOfWhole() > o2.getFractionOfWhole()) {
                return -1;
            } else if (o1.getFractionOfWhole() < o2.getFractionOfWhole()) {
                return 1;
            }

            // hard wiring a minimum size of 5 for segments
            float ns = (float)(maxNPoints/5);

            float ns1 = 1.f - ((float)o1.getSequences().size()/ns);

            float ns2 = 1.f - ((float)o2.getSequences().size()/ns);

            if (ns1 > ns2) {
                return -1;
            } else if (ns1 < ns2) {
                return 1;
            }

            return 0;
        }

        public int compare2(Sequences o1, Sequences o2) {

            // adding a term to prefer the larger
            // fraction, but in a smaller number of
            // larger segments.

            // hard wiring a minimum size of 5 for segments
            float ns = (float)(maxNPoints/5);

            float ns1 = 1.f - ((float)o1.getSequences().size()/ns);

            float ns2 = 1.f - ((float)o2.getSequences().size()/ns);

            //NOTE: this may need to change for cases where,
            // for example, have one very large segment that
            // is the right answer and several smaller matches
            // that are false due to occlusion... presumably
            // other sequences have as many false matches, but
            // this needs alot more testing.

            float s1 = o1.getFractionOfWhole() * ns1;
            float s2 = o2.getFractionOfWhole() * ns2;

            if (s1 > s2) {
                return -1;
            } else if (s1 < s2) {
                return 1;
            }

            if (o1.getFractionOfWhole() > o2.getFractionOfWhole()) {
                return -1;
            } else if (o1.getFractionOfWhole() < o2.getFractionOfWhole()) {
                return 1;
            }

            if (o1.getAbsSumDiffs() < o2.getAbsSumDiffs()) {
                return -1;
            } else if (o1.getAbsSumDiffs() > o2.getAbsSumDiffs()) {
                return 1;
            }

            return 0;
        }
    }

    /**
     * comparator for a preferring high fraction and low differences,
       then descending sort of fraction,
     * then diff, then startIdx
     */
    private class SequenceComparator3 implements
        Comparator<Sequence> {

        private final float maxDiff;

        public SequenceComparator3(float maxDiff) {
            this.maxDiff = maxDiff;
        }

        @Override
        public int compare(Sequence o1, Sequence o2) {

            float d1 = 1.f - (o1.absAvgSumDiffs/maxDiff);
            float d2 = 1.f - (o2.absAvgSumDiffs/maxDiff);

            float s1 = o1.fractionOfWhole + d1;

            float s2 = o2.fractionOfWhole + d2;

            if (s1 > s2) {
                return -1;
            } else if (s1 < s2) {
                return 1;
            }

            if (o1.fractionOfWhole > o2.fractionOfWhole) {
                return -1;
            } else if (o1.fractionOfWhole < o2.fractionOfWhole) {
                return 1;
            }

            if (o1.absAvgSumDiffs < o2.absAvgSumDiffs) {
                return -1;
            } else if (o1.absAvgSumDiffs > o2.absAvgSumDiffs) {
                return 1;
            }

            if (o1.startIdx1 < o2.startIdx1) {
                return -1;
            } else if (o1.startIdx1 > o2.startIdx1) {
                return 1;
            }

            return 0;
        }

    }

    /**
     * comparator for descending sort of fraction,
     * then diff, then startIdx
     */
    private class SequenceComparator2 implements
        Comparator<Sequence> {

        @Override
        public int compare(Sequence o1, Sequence o2) {

            if (o1.fractionOfWhole > o2.fractionOfWhole) {
                return -1;
            } else if (o1.fractionOfWhole < o2.fractionOfWhole) {
                return 1;
            }

            if (o1.absAvgSumDiffs < o2.absAvgSumDiffs) {
                return -1;
            } else if (o1.absAvgSumDiffs > o2.absAvgSumDiffs) {
                return 1;
            }

            if (o1.startIdx1 < o2.startIdx1) {
                return -1;
            } else if (o1.startIdx1 > o2.startIdx1) {
                return 1;
            }
            return 0;
        }

    }

    /**
     * comparator to sort by ascending startIdx, then
     * descending fraction of whole
     */
    private class SequenceComparator implements
        Comparator<Sequence> {

        @Override
        public int compare(Sequence o1, Sequence o2) {

            if (o1.startIdx1 < o2.startIdx1) {
                return -1;
            } else if (o1.startIdx1 > o2.startIdx1) {
                return 1;
            }
            if (o1.fractionOfWhole > o2.fractionOfWhole) {
                return -1;
            } else if (o1.fractionOfWhole < o2.fractionOfWhole) {
                return 1;
            }

            if (o1.startIdx2 < o2.startIdx2) {
                return -1;
            } else if (o1.startIdx2 > o2.startIdx2) {
                return 1;
            }

            if (o1.absAvgSumDiffs < o2.absAvgSumDiffs) {
                return -1;
            } else if (o1.absAvgSumDiffs > o2.absAvgSumDiffs) {
                return 1;
            }

            // should not arrive here
            return 0;
        }
    }

    /**
     *
     * @param md 3 dimensional array of difference matrices
     * @param r block size
     * @return
     */
    private void findMinDifferenceMatrix(
        float[][][] md, int r, double threshold,
        MinDiffs output) {

        if (r < 1) {
            throw new IllegalArgumentException("r cannot be < 1");
        }

        double c = 1./(double)(r*r);

        //md[0:n2-1][0:n1-1][0:n1-1]

        int n1 = md[0].length;
        int n2 = md.length;

        int[] idxs0 = output.idxs0;
        float[] mins = output.mins;

        int count = 0;

        double lThresh = Math.sqrt(r) * threshold;

        for (int jOffset = 0; jOffset < md.length; jOffset++) {
            log.fine(String.format("block=%d md[%d]", r, jOffset));
            float[][] a = md[jOffset];
            float sum = 0;
            for (int i = r; i < a.length; i++) {
                float s1 = a[i][i] - a[i-r][i] - a[i][i-r] + a[i-r][i-r];
                log.finest(
                    String.format(
                    " [%d,%d] %.4f, %.4f, %.4f, %.4f => %.4f",
                    i, i, a[i][i], a[i-r][i], a[i][i-r],
                    a[i-r][i-r], s1*c));

                s1 *= c;

                if (s1 < 0) {
                    if (s1 < -0.016) {
                        // warn if resolution errors are 1 degree or more
                        log.warning("s1=" + s1);
                    }
                    s1 += -1;
                }

                log.fine(String.format(" [%2d,%2d<-%2d] => %.4f",
                    i,
                    ((i + jOffset) < n2) ?
                    i + jOffset : (i + jOffset) - n2,
                    ((i + jOffset - r + 1) < n2) ?
                    i + jOffset - r + 1 : (i + jOffset - r + 1) - n2,
                    s1*c));

log.info("*CHECK: i=" + i + " j=" + (i + jOffset)
+ " jOffset=" + jOffset
+ " d=" + s1 + " r=" + r);

                if (s1 > lThresh) {
                   continue;
                }

                // note, idx from q is i + jOffset
                count++;
                sum += s1;

                if (s1 < mins[i]) {
                    int idx2 = i + jOffset;
                    if (idx2 >= n2) {
                        idx2 -= n2;
                    }
                    mins[i] = s1;
                    idxs0[i] = jOffset;
                }
                /*
                // store each index i-r through i if best
                for (int ii = (i - r); ii <= i; ++ii) {
                    if (s1 < mins[ii]) {
                        int idx2 = ii + jOffset;
                        if (idx2 >= n2) {
                            idx2 -= n2;
                        }
                        mins[ii] = s1;
                        idxs0[ii] = jOffset;
                    }
                }
                */
            }
            if (count == 0) {
                sum = Integer.MAX_VALUE;
            }
            log.fine(String.format(
                "SUM=%.4f block=%d md[%d]", sum, r, jOffset));
        }


        {
            for (int i = 0; i < idxs0.length; ++i) {
                if (mins[i] == Float.MAX_VALUE) {
                    continue;
                }
                int j = i + idxs0[i];
                if (j > n2) {
                    j -= n1;
                }
                log.info("MIN i=" + i + " j="
                + j + " offset=" + idxs0[i] + "  mind=" + mins[i]);
            }
        }


        log.info("OFFSETS=" + Arrays.toString(idxs0));
        log.info("mins=" + Arrays.toString(mins));

    }

    /**
     *
     * @param md
     * @param r
     * @param mins
     * @param threshold
     * @param tolerance
     * @param n1
     * @param n2
     * @param output contains pairs of i and jOffset, where
     * j is i + jOffset
     */
    private void findEquivalentBest(float[][][] md, int r,
        MinDiffs mins, double threshold, double tolerance,
        int n1, int n2, DiffMatrixResults output) {

        //md[0:n2-1][0:n1-1][0:n1-1]

        assert(md.length == n2);

        double c = 1./(double)(r*r);

        double lThresh = Math.sqrt(r) * threshold;

        // capture all "best" within mins[i] += tolerance

        for (int jOffset = 0; jOffset < n2; jOffset++) {
            float[][] a = md[jOffset];
            for (int i = r; i < n1; ++i) {

                if (mins.idxs0[i] == -1) {
                    // there is no best for this p index
                    // this should not happen with current
                    // read pattern
                   // log.severe("ERROR: min not set for i=" + i);
                    continue;
                }

                // mins.mins[i] is the best for index i (== P_i)
                // mins.idxs0[i] is jOffset of best
                // j is index i + jOffset

                float s1 = a[i][i] - a[i-r][i] - a[i][i-r] + a[i-r][i-r];

                s1 *= c;

                if (s1 > lThresh) {
                   continue;
                }

                if (s1 < 0) {
                    s1 *= -1;
                }

                double best = mins.mins[i];

                if (Math.abs(s1 - best) > tolerance) {
                    continue;
                }

                int idx2 = jOffset + i;
                if (idx2 >= n2) {
                    idx2 -= n2;
                }

                output.add(i, idx2, jOffset, (float)s1,
                    r, n1, n2);
            }
        }
    }

    private float findMaxAvgDiff(List<Sequence> sequences) {

        float max = Float.MIN_VALUE;
        for (Sequence s : sequences) {
            float d = s.absAvgSumDiffs;
            if (d > max) {
                max = d;
            }
       }
        return max;
    }

    private boolean assertNoWrapAround(List<Sequence> sequences) {
        for (Sequence s : sequences) {
            boolean valid = s.assertNoWrapAround();
            if (!valid) {
                return false;
            }
        }
        return true;
    }

    private boolean assertNoWrapAround2(List<Sequences> sequences) {
        for (Sequences seqs : sequences) {
            boolean valid = assertNoWrapAround(seqs.getSequences());
            if (!valid) {
                return false;
            }
        }
        return true;
    }
}
