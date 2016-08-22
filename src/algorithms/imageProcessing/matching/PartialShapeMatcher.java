package algorithms.imageProcessing.matching;

import algorithms.QuickSort;
import algorithms.compGeometry.LinesAndAngles;
import algorithms.imageProcessing.features.RANSACEuclideanSolver;
import algorithms.imageProcessing.transform.EuclideanTransformationFit;
import algorithms.imageProcessing.transform.TransformationParameters;
import algorithms.imageProcessing.transform.Transformer;
import algorithms.imageProcessing.util.MatrixUtil;
import algorithms.misc.Misc;
import algorithms.misc.MiscMath;
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
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import thirdparty.HungarianAlgorithm;
import thirdparty.edu.princeton.cs.algs4.Interval;
import thirdparty.edu.princeton.cs.algs4.IntervalRangeSearch;

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
    Also, cases that have projections that are not 
    euclidean need improvements in the addByTransformation
    method.
    
    when articulated works well, need to offer an
    option to perform and save results before and 
    after the articulated additions so that user
    can examine both and decide between them.
    
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

    private boolean srchForArticulatedParts = false;

    // this helps to remove points far from
    // euclidean transformations using RANSAC.
    // it should probably always be true.
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
        int[] rs = new int[]{
            n1/2, 
            n1/5
        };
        
        int topK = 1;

        MergedMinDiffs mergedMinDiffs = null;
        
        for (int r : rs) {

            if (r < 2) {
                // r=1 reads diagonal of 0's only
                r = 2;
            }

            // build the matching sequential sequences by
            // by reading chord difference over a block size
            // and aggregating sequential same offsets
            MergedMinDiffs mmd = extractSequences(md, r, thresh);
            if (mergedMinDiffs == null) {
                mergedMinDiffs = mmd;
            } else {
                mergedMinDiffs = merge(md, mergedMinDiffs, mmd);
            }
        }
      
        MergedMinDiffs2 mergedMinDiffs2 = 
            new MergedMinDiffs2(mergedMinDiffs, n1, n2);
        populateChordDifferences(md, mergedMinDiffs2);
   
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
            // note that RANSAC has been used remove outliers
            // from the already matched points too and those
            // are not tracked, just removed.
            
            // the added points are additionally added to this
            // separate list
            List<PairIntArray> addedPoints = 
                new ArrayList<PairIntArray>(topK);
            
            if (diffN <= 0) {
                results = transformAndEvaluate(mergedMinDiffs2, p, q,
                    md, pixTolerance, topK, addedPoints);
            } else {
                results = transformAndEvaluate(mergedMinDiffs2, q, p,
                    md, pixTolerance, topK, addedPoints);
            }
            
            if (results != null && !results.isEmpty()) {
                
                // TODO: calc equiv offsets of addedPoints and
                // if the mmd2 entry w/ that offset and
                // point(s) has significant number of 
                // consecutive points not already in results[i],
                // add those to results[i].
                // TODO: ideally, might want to evaluate
                // this candidate section of consecutive points 
                // with euclidean transformation before adding.
                // -- any added points from an mmd2 item should
                //    be noted or zeroed out so can't be added
                //    in next code block.
                // In summary, the methods in "performEuclidTrans"
                // start with the best main offset of matching
                // points for an mmd2 item, 
                // calculate a euclidean transformation 
                // to roughly add new unmatched points and then
                // from those points, bootstrap to their implied
                // offsets found in another item in mmd2 
                // and use euclidean transformation on those
                // new items (the disjoint consecutive portion)
                // to remove outliers from those and then add them
                // to the best results for the mmd2 item being
                // analyzed.
            }

            if (srchForArticulatedParts) {
                
                /*
                TODO:
                combinedBestDisjoint needs to receive results 
                and mergedMinDiffs[i>=topK]
                */
                
                best = combineBestDisjoint(results, md,
                    n1, n2);
            } else {
                if (results == null || results.isEmpty()) {
                    best = null;
                } else {
                    best = results.get(0);
                }
            }
            if (best != null) {
                populateWithChordDiffs(best, md, n1, n2);
            }

        } else {

            if (srchForArticulatedParts) {
                //NOTE: to combine best disjoint, need topK > 1            
                List<Result> results =
                    createSortedResults(mergedMinDiffs2, n1, n2, topK);
                best = combineBestDisjoint(results, md,
                    n1, n2);
                populateWithChordDiffs(best, md, n1, n2);
            } else {
                mergedMinDiffs2.sortBySalukwdzeDistance();
                best = createResult(mergedMinDiffs2, 0);
            }
        }

        if (diffN <= 0) {
            return best;
        }

        best = best.transpose();

        return best;
    }

    private MergedMinDiffs condense(MergedMinDiffs mmd, 
        int n1, int n2) {
    
        Map<Integer, IntervalRangeSearch<Integer, Integer>> offsetMap = 
            new HashMap<Integer, IntervalRangeSearch<Integer, Integer>>();
        
        TIntSet rmSet = new TIntHashSet();
        TIntIntMap collisionMap = new TIntIntHashMap();
        
        for (int i = 0; i < mmd.offsets.length; ++i) {
            int offset = mmd.offsets[i];
            int startI = mmd.getStartIMinusBlock(i);
            int stopI = mmd.stopIs[i];
       
            Integer key = Integer.valueOf(offset);
            
            IntervalRangeSearch<Integer, Integer> rt =
                offsetMap.get(key);
            
            if (rt == null) {
                rt = new IntervalRangeSearch<Integer, Integer>();
                offsetMap.put(key, rt);
            }
       
            // NOTE: making the interval min and max
            // 1 pixel wider in order to collide with
            // adjacent intervals too, so can merge them
            // also.
            
            Interval<Integer> interval = new Interval<Integer>(
                startI - 1, stopI + 1);
   
            Integer replaced = rt.put(interval, Integer.valueOf(i));

            if (replaced != null) {
                
                // existing interval, so merge
                
                int mergeIntoIdx = replaced.intValue();
                if (collisionMap.containsKey(mergeIntoIdx)) {
                    mergeIntoIdx = collisionMap.get(mergeIntoIdx);
                }
                collisionMap.put(i, mergeIntoIdx);
                rmSet.add(i);
                                
                /*                
                collision with "replaced" and "i"
                   - merge i with mmd[mergeIntoIdx]
                   - add i to rmSet
                   - create a pointer from i to mergeIntoIdx
                     so that the next collision can find
                     the original merge index
                */
                int mstartI = mmd.getStartIMinusBlock(mergeIntoIdx);
                int mstopI = mmd.stopIs[mergeIntoIdx];
            
                if (interval.min().intValue() <= mstartI) {
                    if (interval.max() > mstopI) {
                        // interval is larger than current start to stop
                        mmd.stopIs[mergeIntoIdx] = mmd.stopIs[i];
                    }
                    // put back into format of startI and separate r
                    mmd.startRs[mergeIntoIdx] = mmd.startRs[i];
                    mmd.startIs[mergeIntoIdx] = mmd.startIs[i];
                } else if (interval.min() <= mstopI &&
                    (interval.max() > mstopI)) {
                    mmd.stopIs[mergeIntoIdx] = mmd.stopIs[i];
                } else if (interval.min() >= mstartI 
                    && interval.max() <= mstopI) {
                    // interval is within existing mmd1 interval
                    // make sure it gets deleted
                }
            }
        }
        
        if (rmSet.isEmpty()) {
            return mmd;
        }
        
        int nTot = mmd.offsets.length - rmSet.size();
        
        MergedMinDiffs mmdm = new MergedMinDiffs(nTot, n1);
        
        int count = 0;
        for (int i = 0; i < mmd.offsets.length; ++i) {
            if (rmSet.contains(i)) {
                continue;
            }
            mmdm.offsets[count] = mmd.offsets[i];
            mmdm.startIs[count] = mmd.startIs[i];
            mmdm.startRs[count] = mmd.startRs[i];
            mmdm.stopIs[count] = mmd.stopIs[i];
            count++;
        }
        assert(count == nTot);
    
        log.fine("condensed " + rmSet.size() + " items");
        
        return mmdm;
    }

    private void populateChordDifferences(float[][][] md, 
        MergedMinDiffs2 mmd2) {

        if (!mmd2.chordsNeedUpdates) {
            return;
        }
        
        int n1 = mmd2.n1;
        int n2 = mmd2.n2;
        
        int n = mmd2.sumChordDiffs.length;
        assert(n <= n1);
        
        int[] rUsed = new int[1];
        
        for (int i = 0; i < n; ++i) {
           
            int offset = mmd2.mmd.offsets[i];

            int stopI = mmd2.mmd.stopIs[i];
            assert(stopI <= n1);
            
            int r = mmd2.getLength(i);
            assert(r <= n1);

            float s1 = read(md, stopI, offset, r, 
                rUsed);
            
            mmd2.sumChordDiffs[i] = s1;        
        }
         
        mmd2.chordsNeedUpdates = false;
    }
    
    private class MergedMinDiffs2 {
        private final MergedMinDiffs mmd;
        private final float[] sumChordDiffs;
        private final int n1;
        private final int n2;
        private boolean chordsNeedUpdates = true;
        public MergedMinDiffs2(MergedMinDiffs mmd, int n1, int n2) {
            this.mmd = mmd;
            sumChordDiffs = new float[mmd.offsets.length];
            this.n1 = n1;
            this.n2 = n2;
            assert(mmd.n1 == n1);
        }
        public int length() {
            return sumChordDiffs.length;
        }
        public int getLength(int index) {
            return mmd.getLength(index);
        }
        public float getFraction(int index) {
            return mmd.getFraction(index);
        }
        public int getStartIMinusBlock(int index) {
            return mmd.getStartIMinusBlock(index);
        }
        
        public void sortBySalukwdzeDistance() {
        
            if (chordsNeedUpdates) {
                throw new IllegalStateException(
                "the chord difference sums must be filled");
            }
            
            assert(sumChordDiffs.length == mmd.offsets.length);
            
            float maxChord = MiscMath.findMax(sumChordDiffs);
            
            // make an array of S distance and an
            // array of indexes and sort by S distance
            
            float[] sd = new float[sumChordDiffs.length];
            int[] indexes = new int[sd.length];
            
            for (int i = 0; i < sd.length; ++i) {
                float d = sumChordDiffs[i]/maxChord;
                float f = 1.f - getFraction(i);
                sd[i] = f*f + d*d;
                indexes[i] = i;
            }
            
            QuickSort.sortBy1stArg(sd, indexes);
            
            int n = sd.length;
            int[] off = Arrays.copyOf(mmd.offsets, n);
            int[] str1 = Arrays.copyOf(mmd.startIs, n);
            int[] strR1 = Arrays.copyOf(mmd.startRs, n);
            int[] stp1 = Arrays.copyOf(mmd.stopIs, n);
            float[] sumDiffs = Arrays.copyOf(sumChordDiffs, n);
            for (int i = 0; i < n; ++i) {
                int idx = indexes[i];
                mmd.offsets[i] = off[idx];
                mmd.startIs[i] = str1[idx];
                mmd.startRs[i] = strR1[idx];
                mmd.stopIs[i] = stp1[idx];
                sumChordDiffs[i] = sumDiffs[idx];
            }
            
            if (debug) {
                print("PSORT");
            }
        }
        
        public void print(String label) {
            
            int n = sumChordDiffs.length;
                        
            for (int i = 0; i < n; ++i) {
                
                log.info(String.format(
                    "%d %s offset=%d  f=%.4f  d=%.4f  startI=%d  stopI=%d len=%d", 
                     i, label, mmd.offsets[i], 
                     getFraction(i), sumChordDiffs[i],
                     getStartIMinusBlock(i), mmd.stopIs[i], getLength(i)));
            }
        }
    }
    
    private class MergedMinDiffs {
        final int n1;
        int[] offsets;
        int[] startIs;
        int[] startRs;
        int[] stopIs;
        public MergedMinDiffs(int size, int n1) {
            offsets = new int[size];
            startIs = new int[size];
            startRs = new int[size];
            stopIs = new int[size];
            this.n1 = n1;
        }
        
        public int getLength(int index) {
            int len = stopIs[index] - 
                getStartIMinusBlock(index) + 1;
            return len;
        }
           
        public float getFraction(int index) {
            float len = getLength(index);
            assert(len <= n1);
            float f = len/(float)n1;
            return f;
        }
        
        public int getStartIMinusBlock(int index) {
            return startIs[index] - startRs[index] + 1;
        }
    
        public void print(String label) {
            
            int n = offsets.length;
                        
            for (int i = 0; i < n; ++i) {
                
                log.info(String.format(
                    "%d %s offset=%d  f=%.4f  startI=%d  stopI=%d len=%d", 
                     i, label, offsets[i], 
                     getFraction(i), getStartIMinusBlock(i),
                     stopIs[i], getLength(i)));
            }
        }
    }
    
    private MergedMinDiffs merge(float[][][] md,
        MergedMinDiffs mmd1, MergedMinDiffs mmd2) {
    
        int n2 = md.length;
        int n1 = md[0].length;
        
        int nTot = mmd1.offsets.length + mmd2.offsets.length;
        
        MergedMinDiffs comb = new MergedMinDiffs(nTot, n1);
        
        int nmmd1 = mmd1.offsets.length;
        System.arraycopy(mmd1.offsets, 0, comb.offsets, 0, nmmd1);
        System.arraycopy(mmd1.startIs, 0, comb.startIs, 0, nmmd1);
        System.arraycopy(mmd1.startRs, 0, comb.startRs, 0, nmmd1);
        System.arraycopy(mmd1.stopIs, 0, comb.stopIs, 0, nmmd1);
        int count = nmmd1;
        for (int i = 0; i < mmd2.offsets.length; ++i) {
            comb.offsets[count] = mmd2.offsets[i];
            comb.startIs[count] = mmd2.startIs[i];
            comb.startRs[count] = mmd2.startRs[i];
            comb.stopIs[count] = mmd2.stopIs[i];
            count++;
        }
        assert(count == nTot);      
            
        if (debug) {        
            mmd1.print("mmd1 before condense");
        }
        
        mmd1 = condense(comb, n1, n2);
            
        if (debug) {
            mmd1.print("mmd1 after condense");
        }
        
        return mmd1;
    }

    protected MergedMinDiffs extractSequences(float[][][] md,
        int r, float thresh) {

        //md[0:n2-1][0:n1-1][0:n1-1]

        int n2 = md.length;
        int n1 = md[0].length;

        MinDiffs mins = new MinDiffs(n1);
        findMinDifferenceMatrix(md, r, thresh, mins);

        // merge ranges of sequential offsets.
        // key=offset, data=start i, start r, stop i

        MergedMinDiffs mmd = new MergedMinDiffs(n1, n1);
        int[] offsets = mmd.offsets;
        int[] startIs = mmd.startIs;
        int[] startRs = mmd.startRs;
        int[] stopIs = mmd.stopIs;

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
            return null;
        }

        // once more to equal length
        ++cIdx;

        mmd.offsets = Arrays.copyOf(offsets, cIdx);
        mmd.startIs = Arrays.copyOf(startIs, cIdx);
        mmd.startRs = Arrays.copyOf(startRs, cIdx);
        mmd.stopIs = Arrays.copyOf(stopIs, cIdx);

        // re-read matrix to find eqivalent best
        // for same block size and existing ranges
        // at given offset.
        
        int stopI, readI;
        int[] rUsed = new int[1];
        for (int i = 0; i < cIdx; ++i) {

            currentOffset = mmd.offsets[i];

            // read backwards from current start block
            readI = mmd.startIs[i] - mmd.startRs[i];
            while (readI > 0) {

                if (mins.idxs0[readI] == -1) {
                    break;
                }

                startBlock = mmd.startRs[i];
                
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
                
                mmd.startIs[i] = readI;
                mmd.startRs[i] = rUsed[0];
                readI = mmd.startIs[i] - mmd.startRs[i];
            }

            //TODO: add a read forward from stopI + block size
        
        }
        
        return mmd;
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
        int r1 = blockSize;
        
        float s1;
        if ((i - r1 + 1) < 0) {
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
            s1 = a[i][i] - a[i-r1+1][i] - a[i][i-r1+1]
                + a[i-r1+1][i-r1+1];
            float c = 1.f/((float)r1*r1);
            s1 *= c;
        }

        if (s1 < 0) {
            s1 *= -1.f;
        }  

        if (blockUsed != null) {
            blockUsed[0] = r1;
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
            float n = (n1 < n2) ? n1 : n2;
            return (float)idx1s.size()/(float)n;
        }
        
        protected double getNormalizedChordDiff(double maxChordSum) {
            double d = chordDiffSum/maxChordSum;
            return d;
        }
        
        /**
         * The Salukwdze distance is the metric used as a 
         * cost in comparisons.
         * (note that other portions of the code use
         * the square of the distance.  also note that 
         * the maximum chord sum that was used to determine
         * a best solution is not stored, because the
         * Result membership grows afterwards depending
         * upon options, though this could change in 
         * the future).
         * @param maxChordSum
         * @return 
         */
        public float calculateSalukwdzeDistance(double maxChordSum) {
            float f = 1.f - getFractionOfWhole();
            double d = getNormalizedChordDiff(maxChordSum);
            float s = (float)Math.sqrt(f * f + d * d);
            return s;
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
        MergedMinDiffs2 mmd2, int n1, int n2, int topK) {

        if (mmd2.length() == 0) {
            return null;
        }

        mmd2.sortBySalukwdzeDistance();

        //NOTE: at this stage, the top item in sequences
        // is the correct answer in tests so far

        if (topK > mmd2.length()) {
            topK = mmd2.length();
        }

        List<Result> results = new ArrayList<Result>();
        for (int i = 0; i < topK; ++i) {
            Result result = createResult(mmd2, i);
            results.add(result);
        }

        return results;
    }

    /**
     * 
     * @param mmd2 list of merged sequential minima 
     * as ranges of correspondence of p and q with
     * different offsets.  the list contains many
     * different solutions.
     * @param p array of shape p boundary coordinates
     * @param q array of shape q boundary coordinates
     * @param md the summed permuted chord difference matrix.
     * @param pixTol a tolerance for use in matching
     * projected unmatched points using derived transformation.
     * @param topK the number of top items in mmd2
     * to analyze, where top is w.r.t. the order after
     * mmd2 has been sorted by Salukwdze distance 
     * (a side effect of this method).
     * @param addedPointLists an <em>output list</em> 
     * of same size as results (which is size topK or null). 
     * each item contains a list of
     * points added to the result item because it was
     * found via projection.  the points have different
     * "offsets" than the mmd2 best item they were derived
     * from by projection.
     * @return list of topK results from the sorted
     * mmd2 list and it's projection analysis.
       Note that mmd2 and return results have same 
       ordering and indexes, that is results[0] is
       derived from mmd2 item 0.
     * @throws NoSuchAlgorithmException thrown when the
     * chosen algorithm used in ransac is not available.
     * TODO: change that so that it finds an algorithm
       and doesn't throw the exception.
     */
    private List<Result> transformAndEvaluate(
        MergedMinDiffs2 mmd2, PairIntArray p, PairIntArray q,
        float[][][] md, float pixTol, int topK,
        List<PairIntArray> addedPointLists)
        throws NoSuchAlgorithmException {

        if (mmd2.length() == 0) {
            return null;
        }

        mmd2.sortBySalukwdzeDistance();
        
        //NOTE: at this stage, the top item in sequences
        // is the correct answer for the main offset
        // of matches in tests so far

        if (topK > mmd2.length()) {
            topK = mmd2.length();
        }

        List<Result> results = new ArrayList<Result>(topK);
        for (int i = 0; i < topK; ++i) {
            int len = mmd2.getLength(i);
            if (len < 7) {
                // 7 points are needed for the RANSAC algorithm
                continue;
            }
            PairIntArray added = new PairIntArray();
            Result result = addByTransformation(mmd2, i, p, q, 
                pixTol, added);
            if (result != null) {
                populateWithChordDiffs(result, md, p.getN(), q.getN());
                log.fine("calc'ed chords: " + result.toStringAbbrev());
                results.add(result);
                addedPointLists.add(added);
            }
        }

        if (results.isEmpty()) {
            return null;
        }

        return results;
    }

    /**
     * create a Result instance that lacks the
     * distance and chord difference fields.
     * @param mmd2     
     * @param index     
     * @return
     */
    protected Result createResult(MergedMinDiffs2 mmd2, int index) {

        int offset = mmd2.mmd.offsets[index];
        int n1 = mmd2.n1;
        int n2 = mmd2.n2;
        
        int s0 = mmd2.getStartIMinusBlock(index);
        int s1 = mmd2.mmd.stopIs[index];

        Result result = new Result(n1, n2, offset);
        for (int i = s0; i <= s1; ++i) {

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
        
        if (!mmd2.chordsNeedUpdates) {
            result.chordDiffSum = mmd2.sumChordDiffs[index];
        }

        return result;
    }

    protected void populatePointArrays(MergedMinDiffs2 mmd2,
        int index, PairIntArray p, PairIntArray q,
        PairIntArray pOut, PairIntArray qOut,
        PairIntArray pUnmatchedOut,
        PairIntArray qUnmatchedOut) {

        int offset = mmd2.mmd.offsets[index];
        int len = mmd2.getLength(index);

        int n1 = p.getN();
        int n2 = q.getN();
        assert(n1 == mmd2.n1);
        assert(n2 == mmd2.n2);

        TIntSet pIdxs = new TIntHashSet();
        TIntSet qIdxs = new TIntHashSet();

        int s0 = mmd2.getStartIMinusBlock(index);
        int s1 = mmd2.mmd.stopIs[index];
        for (int i = s0; i <= s1; ++i) {

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
        assert(pUnmatchedOut.getN() + pOut.getN() ==
            p.getN());
        assert(qUnmatchedOut.getN() + qOut.getN() ==
            q.getN());
    }

    private Result addByTransformation(MergedMinDiffs2 mmd2, 
        int index, PairIntArray p, PairIntArray q, 
        double pixTol, PairIntArray outputAddedPoints) throws
        NoSuchAlgorithmException {

        int offset = mmd2.mmd.offsets[index];
        int len = mmd2.getLength(index);

        PairIntArray leftXY = new PairIntArray(len);
        PairIntArray rightXY = new PairIntArray(len);
        PairIntArray leftUnmatchedXY
            = new PairIntArray(p.getN() - len);
        PairIntArray rightUnmatchedXY = new PairIntArray(q.getN() -
            len);

        populatePointArrays(mmd2, index, p, q, leftXY, rightXY,
            leftUnmatchedXY, rightUnmatchedXY);

        log.fine("offset=" + offset
            + " lft.n=" + leftXY.getN()
            + " rgt.n=" + rightXY.getN());

        PairIntArray outLeft = new PairIntArray();
        PairIntArray outRight = new PairIntArray();

        return addByTransformation(p, q,
            leftXY, rightXY,
            leftUnmatchedXY, rightUnmatchedXY,
            pixTol, offset,
            outLeft, outRight, outputAddedPoints);
     }
    
    private Result addByTransformation(
        PairIntArray p, PairIntArray q,
        PairIntArray left, PairIntArray right,
        PairIntArray leftUnmatched, PairIntArray rightUnmatched,
        double pixTol, int origOffset,
        PairIntArray outLeft, PairIntArray outRight,
        PairIntArray outAddedPoints) throws
        NoSuchAlgorithmException {
         
        String debugTag = "offset=" + Integer.toString(origOffset);
        
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
            log.warning(debugTag +
                " euclidean transformation scale: "  + params);
        }
        left = outLeft;
        right = outRight;
        log.fine(debugTag + " partial fit=" + fit.toString()
            + " params=" + params
            + " reset left.n=" + left.getN()
            + " right.n=" + right.getN());

        log.fine("dp=" + dp + " pixTol=" + pixTol);

        Transformer transformer = new Transformer();
        PairIntArray leftTr = transformer.applyTransformation(
            params, leftUnmatched);
        PairIntArray leftTr0 = transformer.applyTransformation(
            params, left);

        if (debug) {
            try {
                CorrespondencePlotter plotter = new CorrespondencePlotter(p, q);
                for (int i = 0; i < left.getN(); ++i) {
                    int x1 = left.getX(i);
                    int y1 = left.getY(i);
                    int x2 = right.getX(i);
                    int y2 = right.getY(i);
                    if ((i % 5) == 0) {
                        plotter.drawLineInAlternatingColors(x1, y1,
                            x2, y2, 0);
                    }
                }
                String filePath = plotter.writeImage("_"
                    + "_debug1_" + debugTag);
                plotter = new CorrespondencePlotter(p, q);
                for (int i = 0; i < left.getN(); ++i) {
                    int x1 = left.getX(i);
                    int y1 = left.getY(i);
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
                for (int i = 0; i < leftUnmatched.getN(); ++i) {
                    int x1 = leftUnmatched.getX(i);
                    int y1 = leftUnmatched.getY(i);
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

        boolean useOptimal = false;

        TObjectFloatMap<PairInt> idxMap;
        if (useOptimal) {
            idxMap = optimalMatch(leftTr, rightUnmatched,
                pixTol);
        } else {
            idxMap = nearestMatch(leftTr, rightUnmatched,
                pixTol);
        }

        log.info(debugTag + "transformation nearest matches=" +
            idxMap.size());
       
        /*
        combine results from leftXY-rightXY
        with idxMap where idxMap is
            idx1, idx2 of leftTr and rightUnmatched, resp.
                          left
        */

        TObjectIntMap<PairInt> pPoints = Misc.createPointIndexMap(p);
        TObjectIntMap<PairInt> qPoints = Misc.createPointIndexMap(q);

        Result result = new Result(p.getN(), q.getN(), origOffset);

        TObjectFloatIterator<PairInt> iter = idxMap.iterator();
        for (int i = 0; i < idxMap.size(); ++i) {
            iter.advance();
            PairInt idxIdx = iter.key();
            float dist = iter.value();
            PairInt ell = new PairInt(
                leftUnmatched.getX(idxIdx.getX()),
                leftUnmatched.getY(idxIdx.getX())
            );
            int pIdx = pPoints.get(ell);
            assert(pIdx > -1);
            PairInt ar = new PairInt(
                rightUnmatched.getX(idxIdx.getY()),
                rightUnmatched.getY(idxIdx.getY())
            );
            int qIdx = qPoints.get(ar);
            assert(qIdx > -1);

            result.insert(pIdx, qIdx, dist);
            
            outAddedPoints.add(pIdx, qIdx);
            
            // quick look at properties to find these
            // added points in the min diff merged lists
            {
                if (qIdx < pIdx) {
                    qIdx += q.getN();
                }
                int offsetA = qIdx - pIdx;
                log.info(debugTag 
                    + ": added pair with implied offset=" + 
                    offsetA + " i=" + pIdx);
            }
        }

        for (int i = 0; i < right.getN(); ++i) {
            int diffX = leftTr0.getX(i) - right.getX(i);
            int diffY = leftTr0.getY(i) - right.getY(i);
            float dist = (float)Math.sqrt(diffX * diffX +
                diffY * diffY);

            PairInt ell = new PairInt(left.getX(i),
                left.getY(i));
            int pIdx = pPoints.get(ell);
            assert(pIdx > -1);

            PairInt ar = new PairInt(right.getX(i),
                right.getY(i));
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
                    + "_debug4_" + debugTag);
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
