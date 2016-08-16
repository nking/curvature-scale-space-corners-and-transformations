package algorithms.imageProcessing.matching;

import algorithms.QuickSort;
import algorithms.compGeometry.LinesAndAngles;
import algorithms.misc.Histogram;
import algorithms.misc.HistogramHolder;
import algorithms.util.Errors;
import algorithms.util.IntIntDouble;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TIntDoubleMap;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.TreeMap;
import java.util.logging.Logger;
import javax.swing.JScrollPane;

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

     * @param p
     * @param q
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

        List<Sequence> sequences = new ArrayList<Sequence>();
        List<Sequence> discarded = new ArrayList<Sequence>();

        int rMax = (int)Math.sqrt(n1);
        if (rMax < 2) {
            rMax = 2;
        }
        int rMin = 2;

        // build the matching sequential sequences by
        // searching from block size 2 to size sqrt(n1)
        extractSimilar(md, sequences, discarded, rMin, rMax);

        //changed to form adjacent segments where wrap
        // around is present, so assert format here
        assert(assertNoWrapAround(sequences));

        Sequences sequences0 = matchArticulated(
            sequences, discarded, n1, n2);

        assert(assertNoWrapAround(sequences));

        //addFeasibleDiscarded(sequences0, discarded);

        if (diffN <= 0) {
            return sequences0;
        }

        transpose(sequences0, n1, n2);

        return sequences0;
    }

    protected void extractSimilar(float[][][] md,
        List<Sequence> sequences,
        List<Sequence> discarded, int rMin, int rMax) {

        //md[0:n2-1][0:n1-1][0:n1-1]

        int n2 = md.length;
        int n1 = md[0].length;

        // 23 degrees is 0.4014
        double thresh = 23. * Math.PI/180.;

        MinDiffs mins = new MinDiffs(n1);
        for (int r = rMin; r <= rMax; ++r) {
            findMinDifferenceMatrix(md, r, thresh, mins);
        }

        // 10 degrees is 0.175
        double tolerance = 0.2;//0.1;//0.25;

        DiffMatrixResults equivBest = new DiffMatrixResults(n1);
        for (int r = rMin; r <= rMax; ++r) {
            findEquivalentBest(md, r, mins, thresh, tolerance,
                n1, n2, equivBest);
        }
        
        /*
 equivBest.sortListIndexes();
        
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < n1; ++i) {
            IndexesAndDiff iad = equivBest.indexesAndDiff[i];
            if (iad != null) {
                LinkedList<IntIntDouble> list =
                    iad.list;
                for (IntIntDouble iid : list) {
                    sb.append(String.format(
                        "i=%d j=%d offset=%d d=%.4f\n", 
                        i, iid.getA(), iid.getB(), 
                        iid.getC()));
                }
            }
            log.info(sb.toString());
            sb.delete(0, sb.length());
        }*/
        

        // ----- find sequential correspondences ----
        equivBest.sortListIndexes();

        for (int idx1 = 0; idx1 < n1; ++idx1) {

            IndexesAndDiff indexesAndDiff =
                equivBest.indexesAndDiff[idx1];
            if (indexesAndDiff == null) {
                continue;
            }

            LinkedList<IntIntDouble> list =
                indexesAndDiff.list;

            double sumAbsDiff = 0;

            for (IntIntDouble node : list) {

                int idx2 = node.getA();
                double diff = node.getC();
                sumAbsDiff += Math.abs(diff);
                int offset = node.getB();

                Sequence s = new Sequence(n1, n2, offset);
                s.startIdx1 = idx1;
                s.startIdx2 = idx2;
                s.stopIdx2 = idx2;
                //search through higher index lists to aggregate
                int nextLIdx = idx1 + 1;
                while (nextLIdx < n1) {

                    IndexesAndDiff indexesAndDiff2 =
                        equivBest.indexesAndDiff[nextLIdx];
                    if (indexesAndDiff2 == null) {
                        break;
                    }

                    Map<PairInt, IntIntDouble> jLookupMap2 =
                        indexesAndDiff2.jLookupMap;

                    int idx3 = s.stopIdx2 + 1;
                    
                    PairInt key2 = new PairInt(idx3, offset);
                    IntIntDouble node2 = jLookupMap2.get(key2);
                    if (node2 != null) {
                        s.stopIdx2 = idx3;

                        diff = node2.getC();
                        sumAbsDiff += Math.abs(diff);

                        LinkedList<IntIntDouble> list2
                            = indexesAndDiff2.list;
                        list2.remove(node2);
                        jLookupMap2.remove(key2);
                    } else {
                        break;
                    }

                    nextLIdx++;
                }

                int n = s.length();
                s.fractionOfWhole = (float)n/(float)n1;
                s.absAvgSumDiffs = (float)(sumAbsDiff/(float)n);

                if (s.stopIdx2 - s.startIdx2 > 1) {
                    if (s.absAvgSumDiffs <= tolerance) {
                        sequences.add(s);
                        log.info(String.format(
                            "%d seq %d:%d to %d  frac=%.4f  avg diff=%.4f",
                            (sequences.size() - 1),
                            s.startIdx1, s.startIdx2, s.stopIdx2,
                            s.fractionOfWhole, s.absAvgSumDiffs));
                    } else if (s.absAvgSumDiffs <= 3*tolerance) {
                        discarded.add(s);
                    }
                }
            }
        }

        log.info(sequences.size() + " sequences");
    }

    protected Sequences matchArticulated(List<Sequence> sequences,
        List<Sequence> higherErrorSequences, int n1, int n2) {

        // (1) choose the topK from sequences sorted by fraction
        // and then add to those
        int topK = 10 * (1 + (Math.max(n1, n2))/250);
        int end = (topK > sequences.size()) ? sequences.size() : topK;
        
        Collections.sort(sequences, new SequenceComparator2());

        List<Sequences> tracks = new ArrayList<Sequences>();
        for (int i = 0; i < end; ++i) {
            Sequence s = sequences.get(i);
            Sequences track = new Sequences();
            tracks.add(track);
            track.getSequences().add(s.copy());
            log.info("seed " + i + " : " + s);
        }

        return matchArticulated(sequences, higherErrorSequences,
            tracks, n1, n2);
    }

    protected Sequences matchArticulated(List<Sequence> sequences,
        List<Sequence> higherErrorSequences,
        List<Sequences> seedTracks, int n1, int n2) {

        print0(sequences, "S");

        print0(higherErrorSequences, "DS");

        print("sort0:", sequences);

        // (1) combine seedTracks that have same offset
        combineIfSameOffset(seedTracks);

        Collections.sort(sequences, new SequenceComparator2());

        //TODO: revise the datastructures here to make the
        // runtime complexity smaller after logic changes are finished.

        // (2) put the sorted sequences into sets with keys
        //     being offsets
        TreeMap<Integer, List<Sequence>> seqeuncesMap =
            placeInMapByOffsets(sequences);

        // (3) add to seedTracks, the best of same offset sequences.
        for (int i = 0; i < seedTracks.size(); ++i) {
            Sequences track = seedTracks.get(i);
            int offset = track.getSequences().get(0).getOffset();
            List<Sequence> sList = seqeuncesMap.get(Integer.valueOf(offset));
            track.getSequences().addAll(sList);
            Sequence.mergeSequences(track.getSequences());
        }

        // (4) add to seedTracks, the best of nearest offsets
        //     within a tolerance, if they don't intersect
        //     a current range
        int maxDiffOffset = 5;
        // -- if does not intersect existing range
        // -- if is consistent clockwise

        assert(assertNoWrapAround2(seedTracks));

        for (int i = 0; i < seedTracks.size(); ++i) {
            Sequences track = seedTracks.get(i);
            log.info("pre-sorted track " + i + ": " + track.toString());
        }

        //filterForConsistentClockwise(tracks);

        // calculate the stats for each track (== Sequences)
        for (Sequences track : seedTracks) {
            int sumLen = 0;
            float sumFrac = 0;
            double sumDiffs = 0;
            for (Sequence s : track.getSequences()) {
                int len = s.stopIdx2 - s.startIdx2 + 1;
                float diff = s.absAvgSumDiffs * len;
                sumLen += len;
                sumDiffs += diff;
                sumFrac += s.fractionOfWhole;
            }
            track.setAbsSumDiffs(sumDiffs);
            track.setAvgSumDiffs((float)(sumDiffs/(float)sumLen));
            track.setFractionOfWhole(sumFrac);
        }

        Collections.sort(seedTracks, new TrackComparator(n1));
        for (int i = 0; i < seedTracks.size(); ++i) {
            Sequences track = seedTracks.get(i);
            log.info("track " + i + ": " + track.toString());
        }
        
        if (seedTracks.isEmpty()) {
            return null;
        }

        return seedTracks.get(0);
    }

    protected double matchRigidWithOcclusion(List<Sequence> srquences,
        int n1) {

        throw new UnsupportedOperationException("not yet implemented");
    }

    /**
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
                output[i][j] = a1[i][j] - a2[i][j];
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

    protected boolean merge(List<Sequence> sequences) {

        if (sequences.size() < 2) {
            return false;
        }

        /*
        0
        1
        2
        3  --0
        4  --
        */

        boolean didMerge = false;

        //TODO: the sequences data structure will
        // be revised to keep the sequences in order.

        LinkedHashSet<Sequence> results
            = new LinkedHashSet<Sequence>();

        for (int i = 1; i < sequences.size(); ++i) {
            Sequence s0 = sequences.get(i - 1);
            Sequence s = sequences.get(i);
            Sequence[] merged = s0.merge(s);
            if (merged == null) {
                results.add(s0);
                results.add(s);
            } else {
                for (Sequence st : merged) {
                    results.add(st);
                    didMerge = true;
                }
            }
        }

        if (didMerge) {
            sequences.clear();
            sequences.addAll(results);
        }

        return didMerge;
    }

    private boolean canAppend(List<Sequence> track,
        Sequence testS, int n1) {

        if (intersectsExistingRange(track, testS)) {
            return false;
        }

        //TODO: could improve the datastructure
        // to make this delete faster.  logic in
        // the code is still changing currently.

        track.add(testS);
        boolean isC = Sequences.isConsistentClockwise(track);
        boolean rmvd = track.remove(testS);
        assert(rmvd);

        return isC;
    }

    /**
     * assuming that seedTracks each have only one offset
     * (expecting only one sequence, more specifically),
     * combine the items in seedTracks that have the same
     * offset.  "Combine" means merge if adjacent or overlapping,
     * else append.
     * @param seedTracks
     */
    private void combineIfSameOffset(List<Sequences> seedTracks) {

        if (seedTracks.size() < 2) {
            return;
        }

        // --- make the offset lookup map ----
        TIntObjectMap<TIntSet> offsetIndexMap
            = new TIntObjectHashMap<TIntSet>();

        for (int i = 0; i < seedTracks.size(); ++i) {

            Sequences track = seedTracks.get(i);

            int offset = track.getSequences().get(0).getOffset();

            TIntSet indexes = offsetIndexMap.get(offset);
            if (indexes == null) {
                indexes = new TIntHashSet();
                offsetIndexMap.put(offset, indexes);
            }
            indexes.add(i);
        }

        // ---- combine tracks that have same offset ----

        TIntSet rmSet = new TIntHashSet();

        for (int i = 0; i < seedTracks.size(); ++i) {

            if (rmSet.contains(i)) {
                continue;
            }

            Sequences track = seedTracks.get(i);
            int offset = track.getSequences().get(0).getOffset();

            TIntSet indexes = offsetIndexMap.get(offset);
            TIntIterator iter = indexes.iterator();
            while (iter.hasNext()) {
                int oIdx = iter.next();
                if (oIdx == i || rmSet.contains(oIdx)) {
                    continue;
                }
                rmSet.add(oIdx);

                // merge all of track2 with track
                Sequences track2 = seedTracks.get(oIdx);
                track.getSequences().addAll(track2.getSequences());

                Sequence.mergeSequences(track.getSequences());
            }

            if (!indexes.isEmpty()) {

                // NOTE: stats not yet popualated, so not updating here

                // clear indexes so no others are merged from it
                indexes.clear();
            }
        }

        if (!rmSet.isEmpty()) {
            TIntList rmList = new TIntArrayList(rmSet);
            rmList.sort();
            for (int i = (rmList.size() - 1); i > -1; --i) {
                int rmIdx = rmList.get(i);
                seedTracks.remove(rmIdx);
            }
        }
    }

    /**
     * create a map with key being offset of the Sequence,
     * and value being a list of all Sequences in sequences
     * which have that offset.
     * each Sequence in a sequences item must have same
     * offset already.
     * @param sequences
     * @return
     */
    private TreeMap<Integer, List<Sequence>>
        placeInMapByOffsets(List<Sequence> sequences) {

        TreeMap<Integer, List<Sequence>> map =
            new TreeMap<Integer, List<Sequence>>();

        for (Sequence s : sequences) {

            int offset = s.getOffset();
            Integer key = Integer.valueOf(offset);

            List<Sequence> list = map.get(key);
            if (list == null) {
                list = new ArrayList<Sequence>();
                map.put(key, list);
            }
            list.add(s);
        }

        return map;
    }

    private class IndexesAndDiff {

        private final LinkedList<IntIntDouble> list;

        /**
         * key = pairint of j, jOffset, value=list node
         */
        private final Map<PairInt, IntIntDouble> jLookupMap;

        public IndexesAndDiff() {
            list = new LinkedList<IntIntDouble>();
            jLookupMap = new HashMap<PairInt, IntIntDouble>();
        }

        public void add(int j, int jOffset, double diff) {

            IntIntDouble node = new IntIntDouble(j, jOffset, diff);

            PairInt key = new PairInt(j, jOffset);

            IntIntDouble existing = jLookupMap.get(key);

            if (existing != null) {
                double diff0 = Math.abs(existing.getC());
                if (diff0 > Math.abs(node.getC())) {
                    list.add(node);
                    jLookupMap.put(key, node);
                }
            } else {            
                list.add(node);
                jLookupMap.put(key, node);
            }
        }

        void sortByJ() {
            
            if (list.size() > 1) {

                IntIntDouble[] a
                    = list.toArray(new IntIntDouble[list.size()]);
                QuickSort.sortByA(a);
                list.clear();
                for (IntIntDouble abc : a) {
                    list.add(abc);
                }
            }

            rewriteJLookupMap();
        }

        private void rewriteJLookupMap() {

            jLookupMap.clear();

            for (IntIntDouble node : list) {
                int j = node.getA();
                int jOffset = node.getB();
                PairInt key = new PairInt(j, jOffset);
                jLookupMap.put(key, node);
            }
        }
    }

    private class DiffMatrixResults {
        /**
        index = i
        value = linked list of j,jOffset, and diff
        */
        private final IndexesAndDiff[] indexesAndDiff;

        public DiffMatrixResults(int n) {
            indexesAndDiff = new IndexesAndDiff[n];
        }

        public void add(int i, int j, int jOffset,
            double diff) {
            if (indexesAndDiff[i] == null) {
                indexesAndDiff[i] = new IndexesAndDiff();
            }
            indexesAndDiff[i].add(j, jOffset, diff);
        }

        void sortListIndexes() {
            for (IndexesAndDiff id : indexesAndDiff) {
                if (id != null) {
                    id.sortByJ();
                }
            }
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

        for (int jOffset = 0; jOffset < md.length; jOffset++) {
            log.fine(String.format("block=%d md[%d]", r, jOffset));
            float[][] a = md[jOffset];
            float sum = 0;
            //for (int i = 0; i < a.length; i+=r) {
            for (int i = (r - 1); i < a.length; i++) {
                float s1;
                if ((i - r) > -1) {
                    s1 = a[i][i] - a[i-r][i] - a[i][i-r] + a[i-r][i-r];
                    log.finest(
                        String.format(
                        " [%d,%d] %.4f, %.4f, %.4f, %.4f => %.4f",
                        i, i, a[i][i], a[i-r][i], a[i][i-r],
                        a[i-r][i-r], s1*c));
                } else {
                    s1 = a[i][i];
                    log.finest(
                        String.format(" [%d,%d] %.4f => %.4f",
                        i, i, a[i][i], s1*c));
                }
                s1 *= c;

                log.fine(String.format(" [%2d,%2d<-%2d] => %.4f",
                    i,
                    ((i + jOffset) < n2) ?
                    i + jOffset : (i + jOffset) - n2,
                    ((i + jOffset - r + 1) < n2) ?
                    i + jOffset - r + 1 : (i + jOffset - r + 1) - n2,
                    s1*c));

                float absS1 = s1;
                if (absS1 < 0) {
                    absS1 *= -1;
                }
                if (absS1 > threshold) {
                   continue;
                }
               
                // note, idx from q is i + jOffset
                count++;
                sum += absS1;
                if (absS1 < Math.abs(mins[i])) {
                    int idx2 = i + jOffset;
                    if (idx2 >= n1) {
                        // idx2 - (n1-i) = offset
                        idx2 -= n1;
                    }
                    mins[i] = s1;
                    idxs0[i] = jOffset;

                    if (false)
                    // fill in the rest of the diagonal in this block
                    for (int k = (i-1); k > (i-r); k--) {
                        if (k < 0) {
                            break;
                        }
//NOTE: that the kOffset would need
// to wrap around, suggests the summed
// area table might need to be created
// in opposite direction in y and read
// in opposite direction here, (then kOffset is jOffset - 1, for example)
// ...haven't thought this through yet... 
                        int kOffset = n1 - (i - k);
                        if (kOffset < 0) {
                            continue;
                        }
                        if (absS1 < Math.abs(mins[k])) {
                            idx2 = k + jOffset;
                            if (idx2 >= n1) {
                                idx2 -= n1;
                            }
                            mins[k] = s1;
                            // for consistency between i and j, need an edited
                            // offset instead of jOffset:
                            idxs0[k] = kOffset;
log.info("CHECK: i=" + i + " j=" + (i + jOffset)
+ " jOffset=" + jOffset + " k=" + k
+ " kOffset=" + kOffset + " r=" + r);
                        }
                    }
                }
            }
            if (count == 0) {
                sum = Integer.MAX_VALUE;
            }
            log.fine(String.format(
                "SUM=%.4f block=%d md[%d]", sum, r, jOffset));
        }

        log.fine("OFFSETS=" + Arrays.toString(idxs0));
        log.fine("mins=" + Arrays.toString(mins));
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

        // capture all "best" within mins[i] += tolerance

        for (int jOffset = 0; jOffset < n2; jOffset++) {
            float[][] a = md[jOffset];
            for (int i = 0; i < n1; i+=r) {
                if (mins.idxs0[i] == -1) {
                    // there is no best for this p index
                    continue;
                }

                // mins.mins[i] is the best for index i (== P_i)
                // mins.idxs0[i] is jOffset of best
                // j is index i + jOffset

                float s1;
                if ((i - r) > -1) {
                    s1 = a[i][i] - a[i-r][i] - a[i][i-r] + a[i-r][i-r];
                } else {
                    s1 = a[i][i];
                }
                s1 *= c;

                float absS1 = s1;
                if (absS1 < 0) {
                    absS1 *= -1;
                }
                if (absS1 > threshold) {
                   continue;
                }

                double best = mins.mins[i];

                if (Math.abs(s1 - best) > tolerance) {
                    continue;
                }

                int idx2 = jOffset + i;
                if (idx2 >= n1) {
                    idx2 -= n1;
                }

                output.add(i, idx2, jOffset, s1);

                if (false)
                // fill in the rest of the diagonal in this block
                for (int k = (i-1); k > (i-r); k--) {
                    if (k < 0) {
                        break;
                    }
                    int kOffset = n1 - (i - k);
                    if (kOffset < 0) {
                        continue;
                    }
                    if ((k - r) > -1) {
                        s1 = a[k][k] - a[k-r][k] - a[k][k-r] +
                            a[k-r][k-r];
                    } else {
                        s1 = a[k][k];
                    }
                    s1 *= c;

                    absS1 = s1;
                    if (absS1 < 0) {
                        absS1 *= -1;
                    }
                    if (absS1 > threshold) {
                        continue;
                    }

                    if (Math.abs(s1 - best) > tolerance) {
                        continue;
                    }
                    idx2 = jOffset + k;
                    if (idx2 >= n1) {
                        idx2 -= n1;
                    }
log.info("CHECK: i=" + i + " j=" + (i + jOffset)
+ " jOffset=" + jOffset + " k=" + k
+ " kOffset=" + kOffset + " r=" + r);
                    // for consistency between i and j, need an edited
                    // offset instead of jOffset:
                    output.add(k, idx2, kOffset, s1);
                }
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

    private boolean intersectsExistingRange(
        List<Sequence> existingList, Sequence sTest) {

        for (Sequence s : existingList) {
            if (s.intersects(sTest)) {
                return true;
            }
        }

        return false;
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
