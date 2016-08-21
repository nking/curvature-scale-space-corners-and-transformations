package algorithms.imageProcessing.matching;

import algorithms.QuickSort;
import algorithms.compGeometry.LinesAndAngles;
import algorithms.imageProcessing.features.RANSACEuclideanSolver;
import algorithms.imageProcessing.transform.EuclideanTransformationFit;
import algorithms.imageProcessing.transform.TransformationParameters;
import algorithms.imageProcessing.transform.Transformer;
import algorithms.imageProcessing.util.MatrixUtil;
import algorithms.misc.Misc;
import algorithms.search.KNearestNeighbors;
import algorithms.util.CorrespondencePlotter;
import algorithms.util.PairFloat;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.iterator.TObjectFloatIterator;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.TObjectFloatMap;
import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.map.hash.TObjectFloatHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import thirdparty.HungarianAlgorithm;

/**
 NOTE: NOT READY FOR USE YET.

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

       The runtime complexity for building the integral
       image is O(m*n) where n and m are the number of sampled
       points on the input shapes.

       The runtime complexity for the search of the
       integral image of summed differences and analysis
       will be added here:

 * @author nichole
 */
public class PartialShapeMatcher {

    /*
    TODO:
    The scissors case shows that the articulated solution
    still needs improvement, specifically the method
    combineBestDisjoint(...).
    */

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

    private boolean srchForArticulatedParts = true;

    // this helps to remove points far from
    // euclidean transformations using RANSAC.
    private boolean performEuclidTrans = true;

    private float pixTolerance = 20;

    // 10 degrees is 0.1745
    private float thresh = (float)(Math.PI/180.) * 10.f;;

    protected Logger log = Logger.getLogger(this.getClass().getName());

    private boolean debug = true;

    public void setToArticulatedMatch() {
        srchForArticulatedParts = true;
    }

    public void overrideSamplingDistance(int d) {
        this.dp = d;
    }

    public void setToDebug() {
        debug = true;
        log.setLevel(Level.FINE);
    }

    /**
      NOT READY FOR USE.

      A shape is defined as the clockwise ordered sequence
      of points P_1...P_N
      and the shape to match has points Q_1...Q_N.
      The spacings used within this method are equidistant
      and the default is 5, so override that if a different number
      is needed.  For shapes with less than 200 points
      in the perimeters, tests with override to set
      spacing dp to 1 have worked well.

      The fixed equidistant spacing is invariant to rotation
      and translation, but not to scale, so if the user needs to solve
      for scale, need to do so outside of this method, that is, apply
      scale changes to the datasets before use of this method..

     NOTE: You may want to pre-process the shape points by using
     PairIntArray p = ImageProcessor.extractSmoothedOrderedBoundary

     @param p
     @param q
    */
    public Result match(PairIntArray p, PairIntArray q) throws
        NoSuchAlgorithmException {

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
            amount of occlusion and maybe additional
            shapes present due to segmentation not being
            able to isolate the object completely.
        (2) the assumption of same object but with
           some parts being differently oriented, for
           an example, the scissors opened versus closed.
           The occlusion should be handled for this one too.

        for the multi-objective optimization score,
        need sum of differences in sequence and the fraction
        of the whole as the variables in the
        Salukwdze distance of the Paretto frontier.

        Note, there is an optional ability to use euclidean
        transformation to add points to the solutions
        (still in progress...)
        */

        //TODO: revisit this.  might need to depend
        // upon p.n
        //NOTE: there's currently an error in merge
        // so cannot use more than one r until fix that
        int[] rs = new int[]{
            //n1/2, 
            n1/5
        };
        /*
        if (srchForArticulatedParts) {
            // the scissors tests gave better results
            // with larger block size
            rs = new int[]{n1/2};
        } else {
            rs = new int[]{n1/3};
        };
        */

        int topK = 1;

        List<Sequence> sequences = new ArrayList<Sequence>();

        for (int r : rs) {

            if (r < 2) {
                // r=1 reads diagonal of 0's only
                r = 2;
            }

            // build the matching sequential sequences by
            // by reading chord difference over a block size
            // and aggregating sequential same offsets
            extractSequences(md, sequences, r, thresh);

            assert(assertNoWrapAround(sequences));
        }

        if (rs.length > 1) {
            mergeSequences(sequences);
            for (Sequence s : sequences) {
                if (s.sumDiffsIsNotFilled()) {
                    populateWithChordDiffs(s, md);
                }
            }
        }

        // NOTE: the android statues
        // and scissor tests show that correct
        // main offset is now the top item
        // when sequences are sorted by
        // salukwzde distance.

        Result best;

        if (performEuclidTrans) {

            List<Result> results;

            // solve for transformation, add points near projection,
            // return sorted solutions, best is at top.
            if (diffN <= 0) {
                results = transformAndEvaluate(sequences, p, q,
                    md, pixTolerance, topK);
            } else {
                results = transformAndEvaluate(sequences, q, p,
                    md, pixTolerance, topK);
            }

            if (srchForArticulatedParts) {
                best = combineBestDisjoint(results, md,
                    n1, n2);
            } else {
                if (results == null || results.isEmpty()) {
                    best = null;
                } else {
                    best = results.get(0);
                }
            }

        } else {

            if (srchForArticulatedParts) {
                List<Result> results =
                    createSortedResults(sequences, n1, n2, topK);
                best = combineBestDisjoint(results, md,
                    n1, n2);
                populateWithChordDiffs(best, md, n1, n2);
            } else {
                best = createResult(sequences.get(0),
                    n1, n2);
                populateWithChordDiffs(best, md, n1, n2);
            }
        }

        if (diffN <= 0) {
            return best;
        }

        best = best.transpose();

        return best;
    }

    protected void extractSequences(float[][][] md,
        List<Sequence> sequences, int r,
        float thresh) {

        //md[0:n2-1][0:n1-1][0:n1-1]

        int n2 = md.length;
        int n1 = md[0].length;

        MinDiffs mins = new MinDiffs(n1);
        findMinDifferenceMatrix(md, r, thresh, mins);

        // merge ranges of sequential offsets.
        // key=offset, data=start i, start r, stop i

        int[] offsets = new int[n1];
        int[] startIs = new int[n1];
        int[] startRs = new int[n1];
        int[] stopIs = new int[n1];

        int cIdx = -1;
        int currentOffset = -1;
        int startBlock = -1;
        int startI = -1;
        for (int i = 1; i < mins.idxs0.length; ++i) {
            int offset = mins.idxs0[i];
            if (offset == -1) {
                if (currentOffset != -1) {
                    cIdx++;
                    offsets[cIdx] = currentOffset;
                    startIs[cIdx] = startI;
                    startRs[cIdx] = startBlock;
                    stopIs[cIdx] = i - 1;
                }
                currentOffset = -1;
                startI = -1;
                startBlock = -1;
                continue;
            }
            if (currentOffset == -1) {
                startI = i;
                startBlock = mins.blockSizes[i];
                currentOffset = offset;
            } else if (currentOffset != offset) {
                // create and store sequence
                cIdx++;
                offsets[cIdx] = currentOffset;
                startIs[cIdx] = startI;
                startRs[cIdx] = startBlock;
                stopIs[cIdx] = i - 1;

                currentOffset = offset;
                startI = i;
                startBlock = mins.blockSizes[i];
            }
        }
        if (currentOffset != -1) {
            cIdx++;
            offsets[cIdx] = currentOffset;
            startIs[cIdx] = startI;
            startRs[cIdx] = startBlock;
            stopIs[cIdx] = n1 - 1;
        }

        if (cIdx == -1) {
            return;
        }

        // once more to equal length
        ++cIdx;

        offsets = Arrays.copyOf(offsets, cIdx);
        startIs = Arrays.copyOf(startIs, cIdx);
        startRs = Arrays.copyOf(startRs, cIdx);
        stopIs = Arrays.copyOf(stopIs, cIdx);

        int stopI, readI;
        int[] rUsed = new int[1];
        for (int i = 0; i < cIdx; ++i) {

            currentOffset = offsets[i];

            // read backwards from current start block
            readI = startIs[i] - startRs[i];
            while (readI > 0) {

                startBlock = startRs[i];

                if (mins.idxs0[readI] == -1) {
                    break;
                }

                float min = mins.mins[readI];

                double lThresh = Math.sqrt(startBlock) * thresh;

                float s1 = read(md, readI, currentOffset,
                    startBlock, rUsed);

                log.fine("s1=" + s1 + " lThresh="
                    + lThresh + " min=" + min +
                    " i=" + readI + " r=" + rUsed[0]);

                if (Math.abs(s1 - min) > lThresh) {
                    break;
                }

                startIs[i] = readI;
                startRs[i] = rUsed[0];
                readI = startIs[i] - startRs[i];
            }

            //TODO: add a read forward from stopI + block size
        }

        for (int i = 0; i < cIdx; ++i) {
            Sequence[] seqs = createSequences(
                n1, n2, offsets[i],
                startIs[i] - startRs[i] + 1, stopIs[i]);
            for (Sequence s : seqs) {
                assert(s.length() <= n1);
                assert(s.length() > 0);
                populateWithChordDiffs(s, md);
                sequences.add(s);
            }
        }

        // NOTE: at this stage, a salukdwze sort
        // produces the right answer as item 0

        log.fine(sequences.size() + " sequences");
    }

    protected void mergeSequences(List<Sequence> sequences) {

        // separate the sequences into lists by offset

        TIntObjectMap<List<Sequence>> offsetMap =
            new TIntObjectHashMap<List<Sequence>>();

        for (Sequence s : sequences) {
            int offset = s.getOffset();
            List<Sequence> list = offsetMap.get(offset);
            if (list == null) {
                list = new ArrayList<Sequence>();
                offsetMap.put(offset, list);
            }
            list.add(s);
        }

        // merge each offset list and place results
        // back into sequences.

        TIntObjectIterator<List<Sequence>> iter =
            offsetMap.iterator();
        sequences.clear();
        for (int i = 0; i < offsetMap.size(); ++i) {
            iter.advance();
            int offset = iter.key();

            List<Sequence> list = iter.value();
            Sequence.mergeSequences(list, log);
            sequences.addAll(list);
        }
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

    private void print(String prefix, List<Sequence> sequences) {

        for (int i = 0; i < sequences.size(); ++i) {
            log.info(String.format("%d %s %s", i, prefix,
                sequences.get(i)));
        }
    }

    private void populateWithChordDiffs(Result result,
        float[][][] md, int n1, int n2) {

        //md[0:n2-1][0:n1-1][0:n1-1]

        result.chordDiffSum = 0;

        for (int i = 0; i < result.getNumberOfMatches(); ++i) {

            int idx1 = result.getIdx1(i);
            int idx2 = result.getIdx2(i);

            float d = read(md, idx1, idx2);

            result.addToChordDifferenceSum(d);
        }
    }

    private void populateWithChordDiffs(Sequence s,
        float[][][] md) {

        //md[0:n2-1][0:n1-1][0:n1-1]

        int n1 = s.getN1();
        int offset = s.getOffset();

        // r is block size
        int r = s.length();

        assert(r <= n1);

        int i = s.getStopIdx1();

        float s1 = read(md, i, offset, r, null);

        s.sumDiffs = s1;
    }

    /**
     *
     * @param md
     * @param i
     * @param offset
     * @param blockSize
     * @param blockUsed the output variable to return the
     * actual block size used if i is small.  this variable
     * can be null.
     * @return
     */
    private float read(float[][][] md, int i, int offset,
        int blockSize, int[] blockUsed) {

        //md[0:n2-1][0:n1-1][0:n1-1]

        float[][] a = md[offset];

        // r is block size
        int r = blockSize;

        if ((i - r + 1) < 0) {
            r = i - 1;
            if (i < 2) {
                // read at i=1, block of size 2
                i = 1;
                r = 2;
            }
        }

        float c = 1.f/((float)r*r);
        float s1 = a[i][i] - a[i-r+1][i] - a[i][i-r+1]
            + a[i-r+1][i-r+1];
        s1 *= c;
        if (s1 < 0) {
            s1 *= -1;
        }

        if (blockUsed != null) {
            blockUsed[0] = r;
        }

        return s1;
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

    private Result combineBestDisjoint(List<Result> results,
        float[][][] md, int n1, int n2) {

        //TODO: this one could be improved in many ways.
        // need to look into details of scissors offset=16 test
        // and follow when the other half of the scissors
        // doesn't get merged here.

        if (results == null || results.isEmpty()) {
            return null;
        }

        /*
        if disjoint and clockwise consistent,
        combine into one result to handle
        articulated components.
        */
        Result best = results.get(0);
        best.sortByIdx1();
        for (int i = 1; i < results.size(); ++i) {

            TIntSet p1Best = new TIntHashSet(best.idx1s);
            TIntSet p2Best = new TIntHashSet(best.idx2s);

            Result r = results.get(i);

            // if both non-intersecting sets are added
            // together, assert that the results
            // are clockwise consistent,
            // that is all idx1 are > prev idx1 and same
            // for idx2, both with respect to wrapping around
            // axis.
            // at expense of space complexity will make
            // 2 parallel arrays, first sorted by idx1
            // and assert the idx1 property,
            // then find the location of the idx2 = 0
            // and read idx2 before and after for same
            // validation

            int nConflicts = 0;
            TIntList c1 = new TIntArrayList(best.idx1s);
            TIntList c2 = new TIntArrayList(best.idx2s);

            for (int j = 0; j < r.idx1s.size(); ++j) {
                int idx1 = r.idx1s.get(j);
                if (p1Best.contains(idx1)) {
                    nConflicts++;
                    continue;
                }
                int idx2 = r.idx2s.get(j);
                if (p2Best.contains(idx2)) {
                    nConflicts++;
                    continue;
                }
                c1.add(idx1);
                c2.add(idx2);
            }

            log.fine("best.n=" + best.getNumberOfMatches()
                + " r.n=" + r.getNumberOfMatches() +
                " nConflicting=" + nConflicts +
                " combined disjoint.n=" + c1.size());

            QuickSort.sortBy1stArg(c1, c2);

            // c1 and c2 have CW ordered indexes of shapes
            int jIdx1 = c1.get(0);
            int jIdx2 = c2.get(0);
            int sentinel2 = -1;
            boolean consistent = true;
            for (int j = 1; j < c1.size(); ++j) {
                int idx1 = c1.get(j);
                if (idx1 <= jIdx1) {
                    consistent = false;
                    break;
                }
                jIdx1 = idx1;

                if (sentinel2 == -1) {
                    int idx2 = c2.get(j);
                    if (idx2 <= jIdx2) {
                        sentinel2 = j;
                    }
                    jIdx2 = idx2;
                }
            }
            if (!consistent) {
                continue;
            }
            // -- read up to sentinel2 then sentinel2 and after
            // to validate idx2
            jIdx2 = c2.get(0);
            for (int j = 1; j < sentinel2; ++j) {
                int idx2 = c2.get(j);
                if (idx2 <= jIdx2) {
                    consistent = false;
                    sentinel2 = j;
                    break;
                }
                jIdx2 = idx2;
            }
            if (consistent && (sentinel2 != -1)) {
                jIdx2 = c2.get(sentinel2);
                for (int j = sentinel2 + 1; j < c2.size(); ++j) {
                    int idx2 = c2.get(j);
                    if (idx2 <= jIdx2) {
                        consistent = false;
                        break;
                    }
                    jIdx2 = idx2;
                }
            }
            if (!consistent) {

                int nAdded = sentinel2 - best.getNumberOfMatches();

                //TODO: this should be revised
                float f = (float)nAdded/(float)best.getNumberOfMatches();
                log.fine("fraction of pts to possibly "
                    + "add = " + f + " nAdded=" + nAdded);
                if (f < 0.333 || (nAdded < 1)) {
                    continue;
                }
                // TODO: assert does not include sentinel2
                c1 = c1.subList(0, sentinel2);
                c2 = c2.subList(0, sentinel2);
            }

            // -- these points are non intersecting and
            // clockwise consistent so rewrite the
            // contents of best with them.

            int nAdded = c1.size() - best.getNumberOfMatches();
            if (nAdded < 1) {
                continue;
            }

            Result tmpBest = new Result(n1, n2, best.origOffset);
            tmpBest.idx1s.addAll(c1);
            tmpBest.idx2s.addAll(c2);

            populateWithChordDiffs(tmpBest, md, n1, n2);

            log.fine("*calc'ed chords: " + tmpBest.toStringAbbrev());

            // NOTE: distSum update is not accurate, but is
            // not used after transformations,
            // so if that ever changes,
            // need to store the individual distances
            // in Result in transformAndEvaluate
            // to add these properly when
            // combined with other Result.
            double dBest = best.distSum;
            double d = r.distSum/(double)r.getNumberOfMatches();

            tmpBest.distSum = dBest + (nAdded * d);

            // --- compare the combined to the previous best ---
            List<Result> tmp = new ArrayList<Result>(results);
            tmp.add(tmpBest);
            ResultComparator rc = new ResultComparator(tmp);
            int comp = rc.compare(best, tmpBest);
            //if (comp <= 0) {
                best = tmpBest;
                //TODO: revise how many segments to
                // safely add back in.
                //break;
            //}
        }

        return best;
    }

    /**
     * create a sequence out of the given parameters,
     * and parse it into ranges such that both idx1
     * and idx2 are increasing in both (splits on the
     * wrap around indexes for the shapes).
     * Note that the method does not fill in the
     * sum of chords.
     * @param n1
     * @param n2
     * @param offset
     * @param startIdx1
     * @param stopIdx1
     * @return
     */
    private Sequence[] createSequences(int n1, int n2,
        int offset, int startIdx1, int stopIdx1) {

        assert(startIdx1 >= 0);
        assert(stopIdx1 < n1);

        Sequence s = new Sequence(n1, n2, offset);
        s.startIdx1 = startIdx1;
        int len = stopIdx1 - startIdx1 + 1;
        s.startIdx2 = s.startIdx1 + offset;
        s.stopIdx2 = s.startIdx2 + len - 1;

        Sequence[] seqs = Sequence.parse(s);

        return seqs;
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

        public float getFractionOfWhole() {
            return (float)idx1s.size()/(float)n1;
        }

        void addToChordDifferenceSum(float diff) {
            chordDiffSum += diff;
        }

        public int getOriginalOffset() {
            return origOffset;
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
                "offset=%d nMatched=%d frac=%.4f distSum=%.4f dChordSum=%.4f",
                origOffset, idx1s.size(), getFractionOfWhole(),
                (float)distSum,
                (float)chordDiffSum));
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
                "offset=%d nMatched=%d distSum=%.4f dChordSum=%.4f",
                origOffset, idx1s.size(), (float)distSum,
                (float)chordDiffSum));

            return sb.toString();
        }

        public void sortByIdx1() {
            QuickSort.sortBy1stArg(idx1s, idx2s);
        }
    }

    private List<Result> createSortedResults(
        List<Sequence> sequences, int n1, int n2, int topK) {

        if (sequences.isEmpty()) {
            return null;
        }

        // sorts by Salukwdze distance
        Collections.sort(sequences,
            new SequenceComparator3(sequences));

        if (debug) {
            print("PSORT", sequences);
        }

        if (sequences.size() < topK) {
            topK = sequences.size();
        }

        List<Result> results = new ArrayList<Result>();
        for (int i = 0; i < topK; ++i) {
            Result result = createResult(sequences.get(i),
                n1, n2);
            results.add(result);
        }

        return results;
    }

    private List<Result> transformAndEvaluate(
        List<Sequence> sequences,
        PairIntArray p, PairIntArray q,
        float[][][] md, float pixTol, int topK)
        throws NoSuchAlgorithmException {

        if (sequences.isEmpty()) {
            return null;
        }

        // debug: print sequences sorted by fraction of whole
        Collections.sort(sequences, new SequenceComparator2());

        if (debug) {
            print("FSORT", sequences);
        }

        // sorts by Salukwdze distance
        Collections.sort(sequences,
            new SequenceComparator3(sequences));

        if (debug) {
            print("PSORT", sequences);
        }

        //NOTE: at this stage, the top item in sequences
        // is the correct answer in tests so far

        if (topK > sequences.size()) {
            topK = sequences.size();
        }

        List<Result> results = new ArrayList<Result>();
        for (int i = 0; i < topK; ++i) {
            Sequence s = sequences.get(i);
            if (s.length() < 7) {
                // 7 points are needed for the RANSAC algorithm
                continue;
            }
            Result result = addByTransformation(s, p, q, pixTol);
            if (result != null) {
                populateWithChordDiffs(result, md, p.getN(), q.getN());
                log.fine("calc'ed chords: " + result.toStringAbbrev());
                results.add(result);
            }
        }

        if (results.isEmpty()) {
            return null;
        }

        Collections.sort(results,
            new ResultComparator(results));

        return results;
    }

    /**
     * create a Result instance that lacks the
     * distance and chord difference fields.
     * @param s
     * @param n1
     * @param n2
     * @return
     */
    protected Result createResult(Sequence s, int n1, int n2) {

        int offset = s.getOffset();

        Result result = new Result(n1, n2, offset);
        for (int i = s.startIdx1; i <= s.getStopIdx1();
            ++i) {

            int pIdx = i;
            assert(pIdx < n1);

            int qIdx = pIdx + offset;
            if (qIdx >= n2) {
                qIdx -= n2;
            } else if (qIdx < 0) {
                qIdx += n2;
            }
            result.idx1s.add(pIdx);
            result.idx2s.add(qIdx);
        }

        return result;
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

        log.fine("s=" + s);

        for (int i = s.startIdx1; i <= s.getStopIdx1();
            ++i) {

            int pIdx = i;
            assert(pIdx < n1);

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
        PairIntArray q, double pixTol) throws
        NoSuchAlgorithmException {

        int offset = s.getOffset();

        PairIntArray leftXY = new PairIntArray(s.length());
        PairIntArray rightXY = new PairIntArray(s.length());
        PairIntArray leftUnmatchedXY
            = new PairIntArray(p.getN() - s.length());
        PairIntArray rightUnmatchedXY = new PairIntArray(q.getN() -
            s.length());

        populatePointArrays(s, p, q, leftXY, rightXY,
            leftUnmatchedXY, rightUnmatchedXY);

        log.fine("offset=" + offset
            + " lft.n=" + leftXY.getN()
            + " rgt.n=" + rightXY.getN());

        PairIntArray outLeft = new PairIntArray();
        PairIntArray outRight = new PairIntArray();

        RANSACEuclideanSolver euclid =
            new RANSACEuclideanSolver();
        EuclideanTransformationFit fit = euclid.calculateEuclideanTransformation(
            leftXY, rightXY, outLeft, outRight);

        TransformationParameters params = (fit != null) ?
            fit.getTransformationParameters() : null;

        if (params == null) {
            //TODO: reconsider whether to package up
            // the given sequence s and return it here
            log.fine("offset=" + offset + " no euclidean fit");
            return null;
        }

        // since this class is using equidistant
        // points on shape boundary, need to compare
        // for same scale.
        if (params.getScale() < 0.9 || params.getScale() > 1.1) {
            log.warning("offset=" + offset +
                " euclidean transformation scale: "  + params);
        }
        leftXY = outLeft;
        rightXY = outRight;
        log.fine("offset=" + offset
            + " partial fit=" + fit.toString()
            + " params=" + params
            + " reset left.n=" + leftXY.getN()
            + " right.n=" + rightXY.getN());

        log.fine("dp=" + dp + " pixTol=" + pixTol);

        Transformer transformer = new Transformer();
        PairIntArray leftTr = transformer.applyTransformation(
            params, leftUnmatchedXY);
        PairIntArray leftTr0 = transformer.applyTransformation(
            params, leftXY);

        if (debug) {
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
                    + "_debug1_" + offset);
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
                    + "_debug2_" + offset);

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
                    + "_debug3_" + offset);
            } catch (Throwable t) {
            }
        }

        log.fine("offset=" + offset + " params=" + params);

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

        log.info("for offset=" + offset + "transformation nearest matches=" +
            idxMap.size());
        
        /*
        combine results from leftXY-rightXY
        with idxMap where idxMap is
            idx1, idx2 of leftTr and rightUnmatched, resp.
                          left
        */

        TObjectIntMap<PairInt> pPoints = Misc.createPointIndexMap(p);
        TObjectIntMap<PairInt> qPoints = Misc.createPointIndexMap(q);

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
            assert(pIdx > -1);
            PairInt ar = new PairInt(
                rightUnmatchedXY.getX(idxIdx.getY()),
                rightUnmatchedXY.getY(idxIdx.getY())
            );
            int qIdx = qPoints.get(ar);
            assert(qIdx > -1);

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
            assert(pIdx > -1);

            PairInt ar = new PairInt(rightXY.getX(i),
                rightXY.getY(i));
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
                    if ((i % 5) == 0) {
                        plotter.drawLineInAlternatingColors(
                            x1, y1, x2, y2, 0);
                    }
                }
                String filePath = plotter.writeImage("_"
                    + "_debug4_" + offset);
            } catch (Throwable t) {
            }
        }
        //log.info("offset=27 RESULTS=" + result.toString());

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

    private class MinDiffs {
        // first dimension index: md[*][][]
        int[] idxs0;
        // i is the index of this array and represents the index
        //       of point in p array
        // j is i + idxs[0] and represents the index
        //       of point in q array

        float[] mins;
        int[] blockSizes;
        public MinDiffs(int n) {
            idxs0 = new int[n];
            mins = new float[n];
            blockSizes = new int[n];
            Arrays.fill(idxs0, -1);
            Arrays.fill(mins, Float.MAX_VALUE);
        }
    }

    /**
     * use Salukwdze distance to compare Result
     * instances.
     */
    private class ResultComparator implements
        Comparator<Result> {

        final float maxDist;
        final float maxChord;

        /**
         *
         * @param list all results that the comparator will
         * be used on.  The maximum values of distance and
         * sum of chords are found from this and used to
         * normalize the values to be compared with the
         * compare method.
         */
        public ResultComparator(List<Result> list) {
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
            maxDist = (float)max;
            maxChord = (float)max2;
        }

        @Override
        public int compare(Result o1, Result o2) {

            float d1 = (float)(o1.chordDiffSum/maxChord);
            float d2 = (float)(o2.chordDiffSum/maxChord);

            float f1 = 1.f - o1.getFractionOfWhole();
            float f2 = 1.f - o2.getFractionOfWhole();

            // salukwdze distance squared
            float s1 = f1*f1 + d1*d1;
            float s2 = f2*f2 + d2*d2;

            // needs nMatches to be weighted more highly
            //float s1 = f1*f1 + 2*d1*d1;
            //float s2 = f2*f2 + 2*d2*d2;

            if (s1 < s2) {
                return -1;
            } else if (s1 > s2) {
                return 1;
            }

            if (f1 < f2) {
                return -1;
            } else if (f1 > f2) {
                return 1;
            }

            if (d1 < d2) {
                return -1;
            } else if (d1 > d2) {
                return 1;
            }

            return 0;
        }
    }

    /**
     * paretto frontier salukwdze comparator using
     * difference in chords and fraction of whole
     */
    private class SequenceComparator3 implements
        Comparator<Sequence> {

        final float maxChord;

        public SequenceComparator3(List<Sequence> list) {
            double max = Float.MIN_VALUE;
            for (Sequence s : list) {
                assert(!s.sumDiffsIsNotFilled());
                if (s.sumDiffs > max) {
                    max = s.sumDiffs;
                }
            }
            maxChord = (float)max;
        }

        @Override
        public int compare(Sequence o1, Sequence o2) {

            assert(!o1.sumDiffsIsNotFilled());
            assert(!o2.sumDiffsIsNotFilled());

            float d1 = (float)(o1.sumDiffs/maxChord);
            float d2 = (float)(o2.sumDiffs/maxChord);

            float f1 = 1.f - o1.getFractionOfWhole();
            float f2 = 1.f - o2.getFractionOfWhole();

            // salukwdze distance squared
            float s1 = f1*f1 + d1*d1;
            float s2 = f2*f2 + d2*d2;

            if (s1 < s2) {
                return -1;
            } else if (s1 > s2) {
                return 1;
            }

            if (f1 < f2) {
                return -1;
            } else if (f1 > f2) {
                return 1;
            }

            if (d1 < d2) {
                return -1;
            } else if (d1 > d2) {
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

            float f1 = o1.getFractionOfWhole();
            float f2 = o2.getFractionOfWhole();

            if (f1 > f2) {
                return -1;
            } else if (f1 < f2) {
                return 1;
            }

            assert(!o1.sumDiffsIsNotFilled());
            assert(!o2.sumDiffsIsNotFilled());

            if (o1.sumDiffs < o2.sumDiffs) {
                return -1;
            } else if (o1.sumDiffs > o2.sumDiffs) {
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
     *
     * @param md 3 dimensional array of difference matrices
     * @param r block size
     * @return
     */
    private void findMinDifferenceMatrix(
        float[][][] md, int r, float threshold,
        MinDiffs output) {

        if (r < 1) {
            throw new IllegalArgumentException("r cannot be < 1");
        }

        float c = 1.f/(float)(r*r);

        //md[0:n2-1][0:n1-1][0:n1-1]

        int n1 = md[0].length;
        int n2 = md.length;

        int[] idxs0 = output.idxs0;
        float[] mins = output.mins;
        int[] rs = output.blockSizes;

        int count = 0;

        double lThresh = Math.sqrt(r) * threshold;

        log.finest("lThresh=" + lThresh);

        for (int jOffset = 0; jOffset < md.length; jOffset++) {
            log.finest(String.format("block=%d md[%d]", r, jOffset));
            float[][] a = md[jOffset];
            float sum = 0;
            for (int i = 0; i < a.length; i++) {
                int r1 = r;
                float s1;
                if ((i - r + 1) < 0) {
                    if (i < 2) {
                        r1 = 2;
                        s1 = a[1][1] - a[0][1] - a[1][0]
                            + a[0][0];
                    } else {
                        r1 = i + 1;
                        s1 = a[i][i] - a[i-r1+1][i]
                            - a[i][i-r1+1]
                            + a[i-r1+1][i-r1+1];
                    }
                    float cTmp = 1.f/(float)(r1*r1);
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
                    s1 *= -1;
                }

                log.finest("*CHECK: i=" + i + " j=" + (i + jOffset)
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
                    rs[i] = r1;
                }
            }
            if (count == 0) {
                sum = Integer.MAX_VALUE;
            }
            log.fine(String.format(
                "SUM=%.4f block=%d md[%d]", sum, r, jOffset));
        }

        if (debug) {
            for (int i = 0; i < idxs0.length; ++i) {
                if (mins[i] == Float.MAX_VALUE) {
                    continue;
                }
                int j = i + idxs0[i];
                if (j >= n2) {
                    j -= n2;
                }
                log.info("MIN i=" + i + " j="
                    + j + " offset=" + idxs0[i] + "  mind="
                    + mins[i] + " r=" + rs[i]);
            }
            log.info("OFFSETS=" + Arrays.toString(idxs0));
            log.info("mins=" + Arrays.toString(mins));
        }
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
            if (s1 < -0.016) {
                // warn if resolution errors are 1 degree or more
                log.warning("s1=" + s1 + " deg="  + (s1*180./Math.PI));
            }
            s1 *= -1;
        }

        return s1;
    }
}
