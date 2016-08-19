package algorithms.imageProcessing.matching;

import algorithms.MultiArrayMergeSort;
import algorithms.QuickSort;
import algorithms.compGeometry.LinesAndAngles;
import algorithms.imageProcessing.MiscellaneousCurveHelper;
import algorithms.imageProcessing.features.RANSACEuclideanSolver;
import algorithms.imageProcessing.transform.EuclideanTransformationFit;
import algorithms.imageProcessing.transform.MatchedPointsTransformationCalculator;
import algorithms.imageProcessing.transform.TransformationParameters;
import algorithms.imageProcessing.transform.Transformer;
import algorithms.imageProcessing.util.MatrixUtil;
import algorithms.misc.Histogram;
import algorithms.misc.HistogramHolder;
import algorithms.misc.Misc;
import algorithms.misc.MiscMath;
import algorithms.search.KNearestNeighbors;
import algorithms.search.NearestNeighbor2D;
import algorithms.util.CorrespondencePlotter;
import algorithms.util.Errors;
import algorithms.util.IntIntDouble;
import algorithms.util.PairFloat;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import gnu.trove.iterator.TIntIntIterator;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.iterator.TObjectFloatIterator;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.TObjectFloatMap;
import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.map.hash.TObjectFloatHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.logging.Logger;
import thirdparty.HungarianAlgorithm;

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
    public Result match(PairIntArray p, PairIntArray q) throws NoSuchAlgorithmException {

        log.info("p.n=" + p.getN() + " q.n=" + q.getN());

        if (p.getN() < 2 || q.getN() < 2) {
            throw new IllegalArgumentException("p and q must "
            + " have at least 2 points");
        }
        
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
        pattern (1):
            a. set rMax = n1/3 and use extractSimilat.
            n. create Sequences from best of those results
               without further use of matchArticulated
            c. repeat a and b for rMax = n1/4
            d. compare results and keep best
            e. if optional (unimplemented) performEval is
               set, create a projection to match the remaining
               points which belong to the shape.
        pattern (2):
            same as patthern (1) except that 1.b. uses
            matchArticulated
        */

        boolean usePattern1 = true;
        boolean performEval = false;

        int[] rs = new int[]{n1/3
            //, n1/4
        };

        List<Sequence> sequences = new ArrayList<Sequence>();

        for (int r : rs) {

            if (r < 2) {
                // r=1 reads diagonal of 0's only
                r = 2;
            }

            // build the matching sequential sequences by
            // by reading block size matches
            extractSequences(md, sequences, r);

            //changed to form adjacent segments where wrap
            // around is present, so assert format here
            assert(assertNoWrapAround(sequences));
        }

        Collections.sort(sequences, new SequenceComparator2());
        print("SORT", sequences);

        if (rs.length > 1) {
            Sequence.mergeSequences(sequences);
        }

        // NOTE: the android statues
        // and scissor tests show that correct
        // main offset is within sequences now.

        // evaluate topK or all sequences

        Result best;
        // solve for transformation, add points, and
        // return best solution.
        if (diffN <= 0) {
            best = transformAndEvaluate(sequences, p, q, md);
        } else {
            best = transformAndEvaluate(sequences, q, p, md);
        }

        if (diffN <= 0) {
            return best;
        }

        best = best.transpose();

        return best;
    }

    protected void extractSequences(float[][][] md,
        List<Sequence> sequences, int r) {

        //md[0:n2-1][0:n1-1][0:n1-1]

        int n2 = md.length;
        int n1 = md[0].length;

        // TODO:  need to revise block reading
        // and storage of results.
        // reading block at i through i-r in both
        // dimensions is not currently stored
        // as a sequence of i-r through i

        // 23 degrees is 0.4014
        double thresh = (Math.PI/180.) * 10.;//*23.;

        MinDiffs mins = new MinDiffs(n1);
        findMinDifferenceMatrix(md, r, thresh, mins);

        Sequence st = null;
        int currentOffset = -1;
        //for (int i = r; i < mins.idxs0.length; ++i) {
        for (int i = 1; i < mins.idxs0.length; ++i) {
            int offset = mins.idxs0[i];
            if (offset == -1 && st == null) {
                continue;
            }
            int blockSize = r;
            if ((i - r + 1) < 0) {
                if (i < 2) {
                    blockSize = 2;
                } else {
                    blockSize = i + 1;
                }
            }
            if (st == null) {
                assert(offset != -1);
                currentOffset = offset;
                st = new Sequence(n1, n2, offset);
                st.startIdx1 = i - blockSize + 1;
if (st.startIdx1 < 0) {
log.info("i=" + i + " r=" + r + " blockSize=" + blockSize
+ " st=" + st);
}
                assert(st.startIdx1 >= 0);
            } else {
                if (currentOffset != offset) {
                    int stopIdx1 = i - 1;
                    int len = stopIdx1 - st.startIdx1;
                    st.startIdx2 = st.startIdx1 + offset;
                    st.stopIdx2 = st.startIdx2 + len - 1;
                    assert(st.length() == len);

                    sequences.add(st);

                    currentOffset = offset;
                    if (offset != -1) {
                        st = new Sequence(n1, n2, offset);
                        st.startIdx1 = i - blockSize + 1;
                        assert(st.startIdx1 >= 0);
                    } else {
                        st = null;
                    }
                }
            }
        }
        if (st != null) {
            int stopIdx1 = n1 - 1;
            int len = stopIdx1 - st.startIdx1;
            st.startIdx2 = st.startIdx1 + st.getOffset();
            st.stopIdx2 = st.startIdx2 + len - 1;
            assert(st.length() == len);
            sequences.add(st);
        }

        log.info(sequences.size() + " sequences");
    }

    private Sequences createSequencesWithBest(
        List<Sequence> sequences, int n1, int n2) {

        if (sequences.isEmpty()) {
            return null;
        }

        Collections.sort(sequences, new SequenceComparator2());

        print("FSORT", sequences);

        Sequence best = sequences.get(0);

        Sequences track = new Sequences();
        track.getSequences().add(best);

        // calculate the stats for track
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

        return track;
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

        // NOTE: if symmetry is present, might need
        // a period analysis

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

    private void populateWithChordDiffs(Result result,
        float[][][] md, int n1, int n2) {

        //md[0:n2-1][0:n1-1][0:n1-1]

        for (int i = 0; i < result.getNumberOfMatches(); ++i) {
            
            int idx1 = result.getIdx1(i);
            int idx2 = result.getIdx2(i);
            
            float d = read(md, idx1, idx2);
                            
            result.addToChordDifferenceSum(d);
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
        private double distSum = 0;
        private double chordDiffSum = 0;
        private TIntList idx1s = new TIntArrayList();
        private TIntList idx2s = new TIntArrayList();
        private final int n1;
        private final int n2;
        private final int origOffset;
        public Result(int n1, int n2, int offset) {
            this.n1 = n1;
            this.n2 = n2;
            this.origOffset = offset;
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

        void addToChordDifferenceSum(float diff) {
            chordDiffSum += diff;
        }
        
        /**
         * reverse the mappings from list 1 to list 2
         * to the reference frame of list 2 to list 1.
         * @return 
         */
        public Result transpose() {
            
            Result t = new Result(n2, n1, n1 - origOffset);
            t.idx1s.addAll(idx2s);
            t.idx2s.addAll(idx1s);
            t.distSum = distSum;
            t.chordDiffSum = chordDiffSum;
            
            return t;
        }
        
        @Override
        public String toString() {
            
            StringBuilder sb = new StringBuilder();
            sb.append(String.format(
                "offset=%d nMatched=%d distSum=%.4f dChordSum=%.4f",
                origOffset, idx1s.size(), (float)distSum,
                (float)chordDiffSum));
            sb.append("\n");
            
            for (int i = 0; i < idx1s.size(); ++i) {
                sb.append(String.format("%d (%d, %d)\n",
                    i, idx1s.get(i), idx2s.get(i)));
            }
            
            return sb.toString();
        }
    }

    private Result transformAndEvaluate(
        List<Sequence> sequences,
        PairIntArray p, PairIntArray q,
        float[][][] md) throws NoSuchAlgorithmException {

        int topK = 20;
        if (topK > sequences.size()) {
            topK = sequences.size();
        }

        List<Result> results = new ArrayList<Result>();
        for (int i = 0; i < topK; ++i) {
            Sequence s = sequences.get(i);
            if (s.length() < 7) {
                continue;
            }
            Result result = addByTransformation(s, p, q);
            populateWithChordDiffs(result, md, p.getN(), q.getN());
            results.add(result);
        }
        
        if (results.isEmpty()) {
            return null;
        }

        Collections.sort(results,
            new ResultComparator(results, p.getN()));

        {//DEBUG
            for (int i = 0; i < results.size(); ++i) {
                Result r = results.get(i);
                log.info(i + " transform result=" + r.toString());
            }
        }

        return results.get(0);
    }

    protected void populatePointArrays(Sequence s,
        PairIntArray p, PairIntArray q,
        PairIntArray pOut, PairIntArray qOut,
        PairIntArray pUnmatchedOut,
        PairIntArray qUnmatchedOut) {

        int n1 = p.getN();
        int n2 = q.getN();
        int offset = s.getOffset();

        TIntSet pIdxs = new TIntHashSet();
        TIntSet qIdxs = new TIntHashSet();

        for (int i = s.startIdx1; i <= s.getStopIdx1();
            ++i) {

            int pIdx = i;
            if (pIdx >= n1) {
                pIdx -= n1;
            }
            pIdxs.add(pIdx);
            pOut.add(p.getX(pIdx), p.getY(pIdx));

            int qIdx = pIdx + offset;
            if (qIdx >= n2) {
                qIdx -= n2;
            } else if (qIdx < 0) {
                qIdx += n2;
            }
            qIdxs.add(qIdx);
            qOut.add(q.getX(qIdx), q.getY(qIdx));

            log.fine("pIdx=" + pIdx + " qIdx=" + qIdx +
                " offset=" + offset);
        }
        assert(pOut.getN() == s.length());
        assert(qOut.getN() == s.length());

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
        assert(pUnmatchedOut.getN() + pOut.getN() ==
            p.getN());
        assert(qUnmatchedOut.getN() + qOut.getN() ==
            q.getN());
    }

    private Result addByTransformation(Sequence s, PairIntArray p,
        PairIntArray q) throws NoSuchAlgorithmException {

        PairIntArray leftXY = new PairIntArray(s.length());
        PairIntArray rightXY = new PairIntArray(s.length());
        PairIntArray leftUnmatchedXY = new PairIntArray(p.getN() -
            s.length());
        PairIntArray rightUnmatchedXY = new PairIntArray(q.getN() -
            s.length());
        populatePointArrays(s, p, q, leftXY, rightXY,
            leftUnmatchedXY, rightUnmatchedXY);

        log.info("lft.n=" + leftXY.getN() + 
            " rgt.n=" + rightXY.getN());

        int offset = s.getOffset();

        MatchedPointsTransformationCalculator
            tc = new MatchedPointsTransformationCalculator();

        double scale = 1.;
        //TransformationParameters params =
        //    tc.calulateEuclideanGivenScale(
        //    scale, leftXY, rightXY, 0, 0);

        /*
        TransformationParameters params =
            tc.calulateEuclideanWithoutFilter(
            leftXY, rightXY, 0, 0);
        */

        /*
        float[] weights = new float[leftXY.getN()];
        Arrays.fill(weights, 1.f/(float)leftXY.getN());
        float[] outputScaleRotTransXYStDev = new float[4];
        TransformationParameters params =
            tc.calulateEuclidean(
            leftXY, rightXY, weights, 0, 0,
            outputScaleRotTransXYStDev);
        */

        PairIntArray outLeft = new PairIntArray();
        PairIntArray outRight = new PairIntArray();
        RANSACEuclideanSolver euclid =
            new RANSACEuclideanSolver();
        EuclideanTransformationFit fit = euclid.calculateEuclideanTransformation(
            leftXY, rightXY, outLeft, outRight);

        TransformationParameters params = (fit != null) ?
            fit.getTransformationParameters() : null;
        leftXY = outLeft;
        rightXY = outRight;
        if (params == null) {
            return null;
        }
        log.info("offset=" + offset + " partial fit=" + fit.toString());

        // 4 pixels or some factor of dp
        double pixTol = 20;//30;
        log.info("dp=" + dp + " pixTol=" + pixTol);

        Transformer transformer = new Transformer();
        PairIntArray leftTr = transformer.applyTransformation(
            params, leftUnmatchedXY);
        PairIntArray leftTr0 = transformer.applyTransformation(
            params, leftXY);

//if (s.getOffset() == 27) {
    try {
       CorrespondencePlotter plotter = new
           CorrespondencePlotter(p, q);
       for (int i = 0; i < leftXY.getN(); ++i) {
           int x1 = leftXY.getX(i);
           int y1 = leftXY.getY(i);
           int x2 = rightXY.getX(i);
           int y2 = rightXY.getY(i);
           if ((i % 5) == 0) {
               plotter.drawLineInAlternatingColors(x1, y1, x2, y2,
               0);
           }
       }
       String filePath = plotter.writeImage("_" +
           "_debug1_" + offset);
       plotter = new
           CorrespondencePlotter(p, q);
       for (int i = 0; i < leftXY.getN(); ++i) {
           int x1 = leftXY.getX(i);
           int y1 = leftXY.getY(i);
           int x2 = leftTr0.getX(i);
           int y2 = leftTr0.getY(i);
           if ((i % 5) == 0) {
               plotter.drawLine(x1, y1, x2, y2,
               255,0, 0, 0);
           }
       }
       filePath = plotter.writeImage("_" +
           "_debug2_" + offset);

       plotter = new
           CorrespondencePlotter(p, q);
       for (int i = 0; i < leftUnmatchedXY.getN(); ++i) {
           int x1 = leftUnmatchedXY.getX(i);
           int y1 = leftUnmatchedXY.getY(i);
           int x2 = leftTr.getX(i);
           int y2 = leftTr.getY(i);
           if ((i % 5) == 0) {
               plotter.drawLine(x1, y1, x2, y2,
               255,0, 0, 0);
           }
       }

       filePath = plotter.writeImage("_" +
           "_debug3_" + offset);
    } catch (Throwable t) {}
    log.info("offset=" + offset + " params=" + params);
     int z = 1;
//}

        // find the best matches to the unmatched in
        // q

        // optimal is currently hungarian, but
        //     may replace w/ another bipartite matcher
        //     in future

        boolean useOptimal = false;

        TObjectFloatMap<PairInt> idxMap;
        if (useOptimal) {
            idxMap = optimalMatch(leftTr, rightUnmatchedXY,
                pixTol);
        } else {
            idxMap = nearestMatch(leftTr, rightUnmatchedXY,
                pixTol);
        }

        /*
        combine results from leftXY-rightXY
        with idxMap where idxMap is
            idx1, idx2 of leftTr and rightUnmatched, resp.
                          left

        into list of IntIntDouble:
        i, j, diff where i and j are w.r.t. p and q
        */

        TObjectIntMap<PairInt> pPoints = Misc.convertToPointIndex(p);
        TObjectIntMap<PairInt> qPoints = Misc.convertToPointIndex(q);

        Result result = new Result(p.getN(), q.getN(), offset);

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
            PairInt ar = new PairInt(
                rightUnmatchedXY.getX(idxIdx.getY()),
                rightUnmatchedXY.getY(idxIdx.getY())
            );
            int qIdx = qPoints.get(ar);

            result.insert(pIdx, qIdx, dist);
        }

        for (int i = 0; i < rightXY.getN(); ++i) {
            int diffX = leftTr0.getX(i) - rightXY.getX(i);
            int diffY = leftTr0.getY(i) - rightXY.getY(i);
            float dist = (float)Math.sqrt(diffX * diffX +
                diffY * diffY);

            PairInt ell = new PairInt(leftXY.getX(i),
                leftXY.getY(i));
            int pIdx = pPoints.get(ell);

            PairInt ar = new PairInt(rightXY.getX(i),
                rightXY.getY(i));
            int qIdx = qPoints.get(ar);

            result.insert(pIdx, qIdx, dist);
        }

//if (s.getOffset() == 27) {
    try {
       CorrespondencePlotter plotter = new
           CorrespondencePlotter(p, q);
       for (int i = 0; i < result.getNumberOfMatches(); ++i) {
           int idx1 = result.getIdx1(i);
           int idx2 = result.getIdx2(i);
           int x1 = p.getX(idx1);
           int y1 = p.getY(idx1);
           int x2 = q.getX(idx2);
           int y2 = q.getY(idx2);
           if ((i % 5) == 0) {
               plotter.drawLineInAlternatingColors(
                   x1, y1, x2, y2, 0);
           }
       }
       String filePath = plotter.writeImage("_" +
           "_debug4_" + offset);
    } catch (Throwable t) {}
    log.info("offset=27 RESULTS=" + result.toString());
    int z0 = 1;
//}

        return result;
    }

    private TObjectFloatMap<PairInt>
        optimalMatch(PairIntArray xy1,
        PairIntArray xy2, double tolerance) {

        TObjectFloatMap<PairInt> costMap =
            new TObjectFloatHashMap<PairInt>();

        float[][] cost = new float[xy1.getN()][xy2.getN()];
        for (int i1 = 0; i1 < xy1.getN(); ++i1) {
            cost[i1] = new float[xy2.getN()];
            Arrays.fill(cost[i1], Float.MAX_VALUE);
            int x1Tr = xy1.getX(i1);
            int y1Tr = xy1.getY(i1);
            for (int i2 = 0; i2 < xy2.getN(); ++i2) {
                int x2 = xy2.getX(i2);
                int y2 = xy2.getY(i2);
                int diffX = x1Tr - x2;
                int diffY = y1Tr - y2;
                cost[i1][i2] = (float)Math.abs(diffX * diffX +
                    diffY * diffY);
                costMap.put(new PairInt(i1, i2),
                    cost[i1][i2]);
            }
        }

        boolean transposed = false;
        if (cost.length > cost[0].length) {
            cost = MatrixUtil.transpose(cost);
            transposed = true;
        }

        HungarianAlgorithm b = new HungarianAlgorithm();
        int[][] match = b.computeAssignments(cost);

        TObjectFloatMap<PairInt> resultMap =
            new TObjectFloatHashMap<PairInt>();

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
            PairInt p = new PairInt(idx1, idx2);
            float costIJ = costMap.get(p);
            if (costIJ <= tolerance) {
                resultMap.put(p, costIJ);
            }
        }

        return resultMap;
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
            Misc.convertToPointIndex(xy2);

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

    private class ResultComparator implements
        Comparator<Result> {

        final float maxDiff;
        final float maxChord;

        final float n;

        public ResultComparator(List<Result> list, int n1) {
            double max = Float.MIN_VALUE;
            double max2 = Float.MIN_VALUE;
            for (Result r : list) {
                if (r.distSum > max) {
                    max = r.distSum;
                }
                if (r.chordDiffSum > max2) {
                    max2 = r.chordDiffSum;
                }
            }
            maxDiff = (float)max;
            maxChord = (float)max2;
            n = n1;
        }

        @Override
        public int compare(Result o1, Result o2) {
            return compare1(o1, o2);
        }

        public int compare1(Result o1, Result o2) {

            //float d1 = 1.f - (float)(o1.distSum/maxDiff);
            //float d2 = 1.f - (float)(o2.distSum/maxDiff);
            float d1 = 1.f - (float)(o1.chordDiffSum/maxChord);
            float d2 = 1.f - (float)(o2.chordDiffSum/maxChord);

            float f1 = (float)o1.getNumberOfMatches()/n;
            float f2 = (float)o2.getNumberOfMatches()/n;

            float s1 = f1 + d1;

            float s2 = f2 + d2;

            if (s1 > s2) {
                return -1;
            } else if (s1 < s2) {
                return 1;
            }

            if (f1 > f2) {
                return -1;
            } else if (f1 < f2) {
                return 1;
            }

            if (d1 > d2) {
                return -1;
            } else if (d1 < d2) {
                return 1;
            }

            return 0;
        }

        public int compare2(Result o1, Result o2) {

            // adding a term to prefer the larger
            // fraction, but in a smaller number of
            // larger segments.

            // hard wiring a minimum size of 5 for segments
            float ns = (float)(n/5);

            float ns1 = 1.f -
                ((float)o1.getNumberOfMatches()/ns);

            float ns2 = 1.f
                - ((float)o2.getNumberOfMatches()/ns);

            //NOTE: this may need to change for cases where,
            // for example, have one very large segment that
            // is the right answer and several smaller matches
            // that are false due to occlusion... presumably
            // other sequences have as many false matches, but
            // this needs alot more testing.

            float f1 = (float)o1.getNumberOfMatches()/n;
            float f2 = (float)o2.getNumberOfMatches()/n;

            float s1 = f1 * ns1;
            float s2 = f2 * ns2;

            if (s1 > s2) {
                return -1;
            } else if (s1 < s2) {
                return 1;
            }

            if (f1 > f2) {
                return -1;
            } else if (f1 < f2) {
                return 1;
            }

            if (o1.distSum < o2.distSum) {
                return -1;
            } else if (o1.distSum > o2.distSum) {
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
            //for (int i = r; i < a.length; i++) {
            for (int i = 0; i < a.length; i++) {
                float s1;
                if ((i - r + 1) < 0) {
                    int rTmp;
                    if (i < 2) {
                        rTmp = 2;
                        s1 = a[1][1] - a[0][1] - a[1][0]
                            + a[0][0];
                    } else {
                        rTmp = i + 1;
                        s1 = a[i][i] - a[i-rTmp+1][i] 
                            - a[i][i-rTmp+1]
                            + a[i-rTmp+1][i-rTmp+1];
                    }
                    float cTmp = 1.f/(float)(rTmp*rTmp);
                    s1 *= cTmp;
                } else {
                    s1 = a[i][i] - a[i-r+1][i] - a[i][i-r+1]
                        + a[i-r+1][i-r+1];
                    s1 *= c;
                }

                if (s1 < 0) {
                    if (s1 < -0.016) {
                        // warn if resolution errors are 1 degree or more
                        log.warning("s1=" + s1);
                    }
                    s1 += -1;
                }

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
                    + j + " offset=" + idxs0[i] + "  mind="
                    + mins[i] + " r=" + r);
            }
        }

        log.info("OFFSETS=" + Arrays.toString(idxs0));
        log.info("mins=" + Arrays.toString(mins));

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
    
    protected float read(float[][][] md, int i, int j) {
        
        int n1 = md.length;
        int n2 = md[0].length;

        int r = (n1 > 3) ? 3 : n1 - 1;
        
        if ((i - r + 1) < 0) {
            r = 2;
            if (i < 2) {
                // read at i=1, block of size 2
                i = 1;
            }
        }
        
        return read(md, i, j, r);
    }
    
    /**
     * read the summed difference matrix to extract
     * the difference between descriptors for point
     * P_i and Q_j for block size r.
     * @param md
     * @param i
     * @param j
     * @param r
     * @return 
     */
    protected float read(float[][][] md, int i, int j,
        int r) {
        
        if ((i - r + 1) < 0) {
            throw new IllegalArgumentException(
            "i and block size are not readable."
                + " (i-r+1) must be >= 0");
        }
        
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
        
        the method is meant to return the values
        for a single i, j pair.
        since the diagonals are zero, one would
        want to integrate over at least a few
        points on the row, giving the chord
        differences of i and its proceeding few 
        points from j and its proceeding few points.
        */
       
        int n1 = md.length;
        int n2 = md[0].length;
        
        int offset = j - i;
        if (offset < 0) {
            offset += n2;
        } else if (offset >= n2) {
            offset -= n2;
        }
        
        float c = 1.f/((float)r*(float)r);
        
        float[][] a = md[offset];
       
        float s1 = a[i][i] - a[i-r+1][i] - a[i][i-r+1] 
            + a[i-r+1][i-r+1];
        
        s1 *= c;
        if (s1 < 0) {
            s1 += -1;
        }
        
        return s1;
    }
}
