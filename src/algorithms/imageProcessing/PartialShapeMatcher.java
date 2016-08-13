package algorithms.imageProcessing;

import algorithms.QuickSort;
import algorithms.compGeometry.LinesAndAngles;
import algorithms.misc.Histogram;
import algorithms.misc.HistogramHolder;
import algorithms.util.Errors;
import algorithms.util.PairIntArray;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TIntDoubleMap;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Map.Entry;
import java.util.Set;
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
       
      A shape is defined as the CW ordered sequence of points P_1...P_N
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
        
        Note that for the rigid model (exvepting scale transformation)
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
        
        the rigid model allowing occlusiong 
        (not yet implemented) will
        likely be a better solution and is similar to
        matching patterns elsewhere in this project.
        
        focusing first on this articulated match to
        look at the range of ability to identify a
        whole object which may be occluded and which
        may have parts which have separate rigid 
        rotation (such as the scissors opened bersus
        closed).
           
        caveats to the articulated match are that
        greedy solutions built by fraction of whole
        may have many different kinds of errors,
        but composing sequences with the top k 
        fraction quickly leads to an unfeasibly
        large number of sequences to evaluate.
        
        */
        
        List<Sequence> sequences = new ArrayList<Sequence>();
        List<Sequence> discarded = new ArrayList<Sequence>();
                
        extractSimilar(md, sequences, discarded, 3, 3);
        
        // sort by descending fraction of whole
        Collections.sort(sequences, new SequenceComparator2());
        
        Sequence s0 = sequences.get(0);
                
        int rMax = (int)Math.sqrt(n1);
        if (rMax < 2) {
            rMax = 2;
        }
        
        sequences = new ArrayList<Sequence>();
        discarded = new ArrayList<Sequence>();
                
        extractSimilar(md, sequences, discarded, 2, rMax);
        
        Sequences sequences0 = matchArticulated(
            sequences, discarded, s0, n1, n2);
        
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
        
        /*
        TODO: will apply a different pattern of reading
        the blocks and merging results next.
        */
        
        //NOTE: presumably, at least 3 consecutive
        // points are needed for a good match,
        // so will start r block size at 3
        MinDiffs mins = new MinDiffs(n1);
        for (int r = rMin; r <= rMax; ++r) {
            findMinDifferenceMatrix(md, r, thresh, mins);
        }
        
        // 10 degrees is 0.175
        double tolerance = 0.25;
          
        DiffMatrixResults equivBest = new DiffMatrixResults(n1);
        for (int r = rMin; r <= rMax; ++r) {
            findEquivalentBest(md, r, mins, thresh, tolerance,
                n1, n2, equivBest);
        }
                
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < n1; ++i) {
            sb.append(String.format("[%4d]: ", i));
            TIntList list = equivBest.indexes[i];
            if (list == null) {
                sb.append("  NA");
            } else {
                list.sort();
                for (int j = 0; j < list.size(); ++j) {
                    sb.append(Integer.toString(list.get(j)));
                    sb.append(",");
                }
            }
            sb.append(" | ");
            log.info(sb.toString());
            sb.delete(0, sb.length());
        }
               
        // ----- find sequential correspondences ----
        
        for (int idx1 = 0; idx1 < n1; ++idx1) {
            TIntList list = equivBest.indexes[idx1];
            if (list == null) {
                continue;
            }
            TDoubleList diffList = equivBest.diffs[idx1];
            QuickSort.sortBy1stArg(list, diffList);
            
            double sumAbsDiff = 0;
            
            for (int j = 0; j < list.size(); ++j) {
                
                int idx2 = list.get(j);
                double diff = diffList.get(j);
                sumAbsDiff += Math.abs(diff);
                
                Sequence s = new Sequence();
                s.startIdx1 = idx1;
                s.startIdx2 = idx2;
                s.stopIdx2 = idx2;
                //search through higher index lists to aggregate
                int nextLIdx = idx1 + 1;
                while (nextLIdx < n1) {
                    TIntList list2 = equivBest.indexes[nextLIdx];
                    if (list2 == null) {
                        break;
                    }
                    TIntSet set2 = equivBest.indexSets[nextLIdx];
                    int idx3 = s.stopIdx2 + 1;
                    if (set2.contains(idx3)) {                            
                        
                        s.stopIdx2 = idx3;
                        
                        // TODO: replace this expensive operation 
                        // with a better data type after refactor logic
                        int rmIdx = list2.indexOf(idx3);
                        
                        diff = equivBest.diffs[nextLIdx].get(rmIdx);
                        
                        sumAbsDiff += Math.abs(diff);
                        
                        list2.removeAt(rmIdx);
                        equivBest.diffs[nextLIdx].removeAt(rmIdx);
                        set2.remove(idx3);
                    } else {
                        break;
                    }
                    
                    nextLIdx++;
                }
            
                int n = s.stopIdx2 - s.startIdx2 + 1;
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
        List<Sequence> higherErrorSequences, Sequence s0, int n1, int n2) {

        List<Sequences> tracks = new ArrayList<Sequences>();
        Sequences track = new Sequences();
        tracks.add(track);
        track.sequences.add(s0.copy());
        
        return matchArticulated(sequences, higherErrorSequences, 
            tracks, n1, n2);
    }
    
    protected Sequences matchArticulated(List<Sequence> sequences,
        List<Sequence> higherErrorSequences, int n1, int n2) {
    
        // (1) choose the topK from sequences sorted by fraction
        // and then add to those
        int topK = 15;
        
        Collections.sort(sequences, new SequenceComparator2());
                
        List<Sequences> tracks = new ArrayList<Sequences>();
        for (int i = 0; i < topK; ++i) {
            Sequence s = sequences.get(i);
            Sequences track = new Sequences();
            tracks.add(track);
            track.sequences.add(s.copy());
        }
        
        return matchArticulated(sequences, higherErrorSequences, 
            tracks, n1, n2);
    }
    
    protected Sequences matchArticulated(List<Sequence> sequences,
        List<Sequence> higherErrorSequences, 
        List<Sequences> seedTracks, int n1, int n2) {
    
        print0(sequences, "S");
        
        print0(higherErrorSequences, "DS");
    
        //(1) merge overlapping consistent
        mergeOverlappingConsistent(seedTracks,
            sequences, n1, n2);
       
        print("sort0:", sequences);
        
        //TODO: revise the datastructures here to make the
        // runtime complexity smaller after logic changes are finished.
        
        // (6) add to "tracks", sequences w/ similar offset then best
        //      of those that don't conflict w/ existing
        for (int i = 0; i < seedTracks.size(); ++i) {
            
            Sequences currentTrack = seedTracks.get(i);
            Sequence s0 = currentTrack.sequences.get(0);
            
            List<Sequence> copy = new ArrayList<Sequence>(sequences);
            Collections.sort(copy, new SequenceComparator2());
            for (int j = (copy.size() - 1); j > -1; --j) {
                Sequence s = copy.get(j);
                if (s.equals(s0) || 
                    intersectsExistingRange(
                        currentTrack.sequences, s, n1)) {
                    copy.remove(j);
                }
            }
            List<Sequence> copy2 = new ArrayList<Sequence>(copy);
            
            boolean didAdd = true;
            while (didAdd) {
                didAdd = false;
                
                if (copy.isEmpty()) {
                    break;
                }
                
                //TODO: may revise this:
                int maxDiffOffset = 5;
                // prefer same offset if it exists,
                // else, choose top of sort
                Sequence ls = currentTrack.sequences
                    .get(currentTrack.sequences.size() - 1);
                    
                int offset = calcOffset12(ls, n1); 
                int minDiffOffset = Integer.MAX_VALUE;
                Sequence minOffsetS = null;
                for (int j = (copy.size() - 1); j > -1; --j) {
                    Sequence st = copy.get(j);
                    int off = calcOffset12(st, n1);
                    if (off == offset) {                        
                        st = copy.remove(j);
                        if (!canAppend(currentTrack.sequences, st, n1)) {                                           
                            // merge
                            currentTrack.sequences.add(st.copy());
                            merge(currentTrack.sequences, n1, n2);
                            continue;
                        }
                        currentTrack.sequences.add(st.copy());
                        didAdd = true;
                        break;
                    } else {
                        int diffOffset = Math.abs(Math.abs(off) - 
                            Math.abs(offset));
                        if (diffOffset < minDiffOffset) {   
                            if (!canAppend(currentTrack.sequences,
                            st, n1)) {
                                copy2.remove(j);
                                continue;
                            }
                            minDiffOffset = diffOffset;
                            minOffsetS = st;
                        } else if (diffOffset > maxDiffOffset) {
                            copy.remove(j);
                            continue;
                        }
                    }
                }
                if (!didAdd) {
                    if (minDiffOffset <= maxDiffOffset) {                    
                        // TODO: may want to revise this
                        currentTrack.sequences.add(minOffsetS.copy());
                        copy.remove(minOffsetS);
                        didAdd = true;
                    }
                }
            }
            
            Set<Sequence> exists = new HashSet<Sequence>(currentTrack.sequences);
            
            // add the top of sorted, if it fits into remaining space
            // TOOO: consider a higher quality filter here?
            // tolerance = 0.1?  that should probably be done
            // at an earlier stage
            for (int j = (copy2.size() - 1); j > -1; --j) {
                Sequence s = copy2.get(j);
                if (exists.contains(s) || 
                    intersectsExistingRange(
                        currentTrack.sequences, s, n1)) {
                    copy2.remove(j);
                }
            }
            didAdd = true;
            while (didAdd) {
                didAdd = false;
                if (copy2.isEmpty()) {
                    break;
                }
                Sequence st = copy2.get(0);
                while (!canAppend(currentTrack.sequences,
                    st, n1)) {
                    copy2.remove(0);
                    //TODO: replace w/ a linked list or other
                    // to improve runtime complexity here when
                    // logic changes are finished
                    if (copy2.isEmpty()) {
                        st = null;
                        break;
                    }
                    st = copy2.get(0);
                }
                if (st == null) {
                    break;
                }
                currentTrack.sequences.add(st.copy());
                copy2.remove(st);
                didAdd = true;
            }
          
        }
        
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
            for (Sequence s : track.sequences) {
                int len = s.stopIdx2 - s.startIdx2 + 1;
                float diff = s.absAvgSumDiffs * len;
                sumLen += len;
                sumDiffs += diff;
                sumFrac += s.fractionOfWhole;
            }
            track.absSumDiffs = sumDiffs;
            track.avgSumDiffs = (float)(sumDiffs/(float)sumLen);
            track.fractionOfWhole = sumFrac;
        }
        
        // sorting here needs to prefer higher fraction and
        // longer segments too
        Collections.sort(seedTracks, new TrackComparator(n1));
        for (int i = 0; i < seedTracks.size(); ++i) {
            Sequences track = seedTracks.get(i);
            log.info("track " + i + ": " + track.toString());
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
        // later, the matching of merged segments needs
        // a way to prefer high quality matches that
        // include curves over high quality matches
        // of straight lines for examples when aggregating
        // segments.
        // large angles for small diff between i, j
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
          //r range 2 to sqrt(min(n, n))

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

            if (!isConsistentClockwise(sequences.sequences)) {
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
    
    private boolean isConsistentClockwise(
        List<Sequence> sequences) {

        if (sequences.isEmpty()) {
            return true;
        }

        Collections.sort(sequences, new SequenceComparator4());

        /*
         all startIdx1 should be increasing,
         and wrap around should be considered.
         then, all startIdx2 should be increasing
         and wrap around whould be considered.
         */
        Sequence s0 = sequences.get(0);

        int ns = sequences.size();

        // check startIdx1 then startIdx2
        for (int check = 0; check < 2; ++check) {
            boolean wrapped = false;
            int prev = (check == 0) ? s0.startIdx1
                : s0.startIdx2;
            for (int j = 1; j <= ns; ++j) {                   
                Sequence s;
                if (j == ns) {
                    if (check == 0) {
                        break;
                    }
                    s = sequences.get(0);
                } else {
                    s = sequences.get(j);
                }
                int idx = (check == 0) ? s.startIdx1
                    : s.startIdx2;
                if (idx == prev) {
                    return false;
                } else if (idx < prev) {
                    if (wrapped) {
                        return false;
                    }
                    wrapped = true;
                    prev = idx;
                }
                prev = idx;
            } // end loop over j sequences in a track
        } // end loop over check
        
        return true;        
    }

    private void transpose(Sequences sequences, 
        int n1, int n2) {

        sequences.fractionOfWhole *= 
            ((float)n1/(float)n2);
        
        List<Sequence> sqs = sequences.getList();
        
        for (Sequence s : sqs) {
            int n = s.stopIdx2 - s.startIdx2;
            int startIdx2 = s.startIdx1;
            s.startIdx1 = s.startIdx2;
            s.startIdx2 = startIdx2;
            s.stopIdx2 = s.startIdx2 + n;
        }
    }

    private void print(String prefix, List<Sequence> sequences) {
        for (int i = 0; i < sequences.size(); ++i) {
            System.out.println(String.format(
                "%d %s %s", i, prefix, sequences.get(i)));
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

    private void mergeOverlappingConsistent(List<Sequences> 
        seedTracks, List<Sequence> sequences, int n1, int n2) {

        /*
        NOTE: if seedTracks becomes a large number of
        tracks, should use a data structure here to make
        finding an intersection of ranges faster
        (e.g. range tree).
        */
        
        for (Sequences track : seedTracks) {
            
            assert(track.sequences.size() == 1);
            
            Sequence s0 = track.sequences.get(0);
            int s0startIdx1 = s0.startIdx1;
            int s0stopIdx1 = s0startIdx1 +
                (s0.stopIdx2 - s0.startIdx2);
            int s0Offset = s0.startIdx2 - s0.startIdx1;
            
            TIntSet skip = new TIntHashSet();
            
            boolean merged = true;
            while (merged) {
                merged = false;
                for (int i = 0; i < sequences.size(); ++i) {
                    if (skip.contains(i)) {
                        continue;
                    }
                    Sequence s = sequences.get(i);
                    int sOffset = s.startIdx2 - s.startIdx1;
                    if (sOffset != s0Offset) {
                        skip.add(i);
                        continue;
                    }
                    int startIdx1 = s.startIdx1;
                    int stopIdx1 = startIdx1 +
                        (s.stopIdx2 - s.startIdx2);
                    if (startIdx1 == s0startIdx1 && 
                        stopIdx1 == s0stopIdx1) {
                        skip.add(i);
                        continue;
                    }
                    boolean doMerge = false;
                    if (s0startIdx1 >= startIdx1 &&  
                        s0startIdx1 <= stopIdx1) {
                        doMerge = true;
                    } else if (s0stopIdx1 >= startIdx1 &&
                        s0stopIdx1 <= stopIdx1) {
                        doMerge = true;
                    } else if (s0startIdx1 <= startIdx1 &&
                        s0stopIdx1 >= stopIdx1) {
                        doMerge = true;
                    }
                    if (doMerge) {                      
                        merged = merge(s0, s, n1, n2);
                        skip.add(i);
                    }
                }
            }
        }   
    }
    
    protected void merge(List<Sequence> sequences, 
        int n1, int n2) {

        if (sequences.size() < 2) {
            return;
        }
        
        /*
        0
        1
        2 
        3  --0
        4  -- 
        */
               
        for (int i = (sequences.size() - 1); i > 0; --i) {
            boolean merged = merge(sequences.get(i - 1),
                sequences.get(i), n1, n2);
            if (merged) {
                sequences.remove(i);
            }
        }
    }
    
    /**
     * TODO:  this method needs tests
     * 
     * @param mergeInto
     * @param mergeFrom
     * @param n1
     * @param n2
     * @return 
     */
    protected boolean merge(Sequence mergeInto, Sequence mergeFrom,
        int n1, int n2) {
 
        int s0stopIdx1 = mergeInto.startIdx1 +
            (mergeInto.stopIdx2 - mergeInto.startIdx2);
        int s0Offset = calcOffset12(mergeInto, n1);

        int sOffset = calcOffset12(mergeFrom, n1);
        if (sOffset != s0Offset) {
            return false;
        }
        int stopIdx1 = mergeFrom.startIdx1 +
            (mergeFrom.stopIdx2 - mergeFrom.startIdx2);
        if ((mergeFrom.startIdx1 == mergeInto.startIdx1) && 
            (stopIdx1 == s0stopIdx1)) {
            // these are same ranges, so let invoker remove mergeFrom
            return true;
        }

        if (!intersects(mergeInto, mergeFrom, n1)) {
            return false;
        }

        int len0 = mergeInto.stopIdx2 - mergeInto.startIdx2 + 1;
        float f0 = mergeInto.fractionOfWhole;
        float d0 = mergeInto.absAvgSumDiffs;
        float d0Tot = d0 * len0;

        int len = mergeFrom.stopIdx2 - mergeFrom.startIdx2 + 1;
        float f = mergeFrom.fractionOfWhole;
        float d = mergeFrom.absAvgSumDiffs;
        float dTot = d * len;

        // handle wrap around for idx1 axis, 
        // but never wrap the idx2 axis.

        if (stopIdx1 < n1 && s0stopIdx1 < n1) {

            if (mergeFrom.startIdx1 < mergeInto.startIdx1) {
                int nAdded = mergeInto.startIdx1 - mergeFrom.startIdx1;
                d0Tot += (dTot/(float)nAdded);
                mergeInto.startIdx1 = mergeFrom.startIdx1;
                mergeInto.startIdx2 = mergeFrom.startIdx2;
                len0 += nAdded;
            }
            if (stopIdx1 > s0stopIdx1) {
                int nAdded = stopIdx1 - s0stopIdx1;
                d0Tot += (dTot/(float)nAdded);
                s0stopIdx1 = stopIdx1;
                mergeInto.stopIdx2 = mergeFrom.stopIdx2;
                len0 += nAdded;
            }

            mergeInto.fractionOfWhole = (float)len0/(float)n1;
            mergeInto.absAvgSumDiffs = d0Tot/(float)len0;

            return true;
        }

        if (s0stopIdx1 > (n1 - 1)) {

            s0stopIdx1 -= n1;

            if (stopIdx1 > (n1 - 1)) {

                stopIdx1 -= n1;

                boolean didMerge = false;
                
                /*
                s : 0-->lstPt     frstPt   n-1   (lstPt overrun)

                s0: 0-->lstPt     frstPt   n-1   (lstPt overrun)
                */
                if (stopIdx1 > s0stopIdx1) {
                    if (stopIdx1 >= mergeInto.startIdx1) {
                        // the merge would make a full circle
                        // which would not be interpretable w/
                        // current logic and data structure, so
                        // keep the sequences separate.
                        // TODO: might need to add a flag indicating
                        // that the sequences represent full possible
                        // correspondence
                        return false;
                    }
                    /*
                    s : 0----->lstPt     frstPt   n-1   (lstPt overrun)

                    s0: 0-->lstPt     frstPt   n-1   (lstPt overrun)
                    */
                    assert(mergeFrom.stopIdx2 > mergeInto.stopIdx2);
                    int nAdded = stopIdx1 - s0stopIdx1;
                    d0Tot += (dTot/(float)nAdded);
                    s0stopIdx1 = stopIdx1;
                    mergeInto.stopIdx2 = mergeFrom.stopIdx2;
                    len0 += nAdded;
                    didMerge = true;
                } 
                if (mergeFrom.startIdx1 < mergeInto.startIdx1) {
                    /*
                    s : 0-->lstPt  frstPt   n-1   (lstPt overrun)

                    s0: 0-->lstPt     frstPt   n-1   (lstPt overrun)
                    */
                    if (mergeFrom.startIdx1 <= s0stopIdx1) {
                        // the merge would make a full circle
                        // which would not be interpretable w/
                        // current logic and data structure, so
                        // keep the sequences separate.
                        // TODO: might need to add a flag indicating
                        // that the sequences represent full possible
                        // correspondence
                        return false;
                    }
                    assert(mergeFrom.stopIdx2 < mergeInto.stopIdx2);
                    int nAdded = mergeInto.startIdx1 - mergeFrom.startIdx1;
                    d0Tot += (dTot/(float)nAdded);
                    mergeInto.startIdx1 = mergeFrom.startIdx1;
                    mergeInto.startIdx2 = mergeFrom.startIdx2;
                    len0 += nAdded;
                    didMerge = true;
                }
                if (!didMerge) {
                    return false;
                }
             } else {
                boolean didMerge = false;
                if (mergeFrom.startIdx1 < mergeInto.startIdx1) {
                    /*
                    s :           frstPt  lstPt   n-1

                    s0: 0-->lstPt      frstPt   n-1   (lstPt overrun)
                    */
                    if (mergeFrom.startIdx1 <= s0stopIdx1) {
                        // see merge notes in block above
                        return false;
                    }
                    assert(mergeFrom.stopIdx2 < mergeInto.stopIdx2);
                    int nAdded = mergeInto.startIdx1 - mergeFrom.startIdx1;
                    d0Tot += (dTot/(float)nAdded);
                    mergeInto.startIdx1 = mergeFrom.startIdx1;
                    mergeInto.startIdx2 = mergeFrom.startIdx2;
                    len0 += nAdded;
                    didMerge = true;
                }
                if (!didMerge) {
                    return false;
                }
            }
        } else if (stopIdx1 > (n1 - 1)) {

            stopIdx1 -= n1;

            boolean didMerge = false;
            
            if (mergeInto.startIdx1 < mergeFrom.startIdx1) {
                /*
                s : 0-->lstPt               frstPt  n-1  (lstPt overrun)

                s0:                 frstPt  lstPt
                */
                if (stopIdx1 < mergeInto.startIdx1) {
                    int nAdded = n1 - s0stopIdx1 + stopIdx1; 
                    d0Tot += (dTot/(float)nAdded);
                    s0stopIdx1 = stopIdx1;
                    mergeInto.stopIdx2 = mergeFrom.stopIdx2;
                    len0 += nAdded;
                    didMerge = true;
                }
            }
            if (!didMerge) {
                return false;
            }
        }

        mergeInto.fractionOfWhole = (float)len0/(float)n1;
        mergeInto.absAvgSumDiffs = d0Tot/(float)len0;

        return true;
    }

    private int calcOffset12(Sequence s, int n1) {
        int offset = s.startIdx2 - s.startIdx1;
        if (offset < 0) {
            // assuming the smaller is correct answer
            int offset2 = (n1 - s.startIdx1) + s.startIdx2;
            if (offset2 < -offset) {
                return offset2;
            }
        }
        return offset;
    }

    private boolean canAppend(List<Sequence> track, 
        Sequence testS, int n1) {
        
        //24:10
        if (testS.startIdx1==24 && testS.startIdx2==10){
            int z = 1;
        }
        
        if (intersectsExistingRange(track, testS, n1)) {
            return false;
        }
        
        //TODO: could improve the datastructure
        // to make this delete faster.  logic in 
        // the code is still changing currently.
        
        track.add(testS);
        boolean isC = isConsistentClockwise(track);
        boolean rmvd = track.remove(testS);
        assert(rmvd);
        
        return isC;
    }

    private class DiffMatrixResults {
        private TIntSet[] indexSets = null;
        // indexes' indexes are i and values are j
        private TIntList[] indexes = null;
        private TDoubleList[] diffs = null;
        public DiffMatrixResults(int n) {
            indexes = new TIntList[n];
            diffs = new TDoubleList[n];
            indexSets = new TIntSet[n];
        }
        public void add(int i, int j, double diff) {
            if (indexes[i] == null) {
                indexes[i] = new TIntArrayList();
                diffs[i] = new TDoubleArrayList();
                indexSets[i] = new TIntHashSet();
            }
            if (!indexSets[i].contains(j)) {
                indexSets[i].add(j);
                indexes[i].add(j);
                diffs[i].add(diff);
            }
        }
    }
    
    public static class Sequence {
        int startIdx1;
        int startIdx2 = -1;
        int stopIdx2 = -1;
        float absAvgSumDiffs;
        float fractionOfWhole;
        public int getStartIdx1() {
            return startIdx1;
        }
        public int getStartIdx2() {
            return startIdx2;
        }
        public int getStopIdx2() {
            return stopIdx2;
        }
        public Sequence copy() {
            Sequence cp = new Sequence();
            cp.startIdx1 = startIdx1;
            cp.startIdx2 = startIdx2;
            cp.stopIdx2 = stopIdx2;
            cp.absAvgSumDiffs = absAvgSumDiffs;
            cp.fractionOfWhole = fractionOfWhole;
            return cp;
        }
        @Override
        public String toString() {
            StringBuilder sb = new StringBuilder();
            sb.append(String.format(
            "(%d:%d to %d, f=%.4f d=%.4f)", 
            startIdx1, startIdx2, stopIdx2,
            fractionOfWhole, absAvgSumDiffs));
            return sb.toString();
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
    
    public static class Sequences {
        List<Sequence> sequences = new ArrayList<Sequence>();
        float fractionOfWhole;
        double absSumDiffs;
        float avgSumDiffs;
        public List<Sequence> getList() {
            return sequences;
        }
        @Override
        public String toString() {
            StringBuilder sb = new StringBuilder();
            sb.append(String.format("frac=%.4f",
                fractionOfWhole));
            sb.append(String.format(", avgDiff=%.4f,  sumDiff=%.4f",
                avgSumDiffs, absSumDiffs));
            for (Sequence s : sequences) {
                sb.append("\n").append(s.toString());
            }
            return sb.toString();
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
            
            // adding a term to prefer the larger
            // fraction, but in a smaller number of 
            // larger segments.
            
            // hard wiring a minimum size of 5 for segments
            float ns = (float)(maxNPoints/5);
           
            float ns1 = 1.f - ((float)o1.sequences.size()/ns);
            
            float ns2 = 1.f - ((float)o2.sequences.size()/ns);
            
            //NOTE: this may need to change for cases where,
            // for example, have one very large segment that
            // is the right answer and several smaller matches
            // that are false due to occlusion... presumably
            // other sequences have as many false matches, but
            // this needs alot more testing.
            
            float s1 = o1.fractionOfWhole * ns1;
            float s2 = o2.fractionOfWhole * ns2;
            
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
            
            if (o1.absSumDiffs < o2.absSumDiffs) {
                return -1;
            } else if (o1.absSumDiffs > o2.absSumDiffs) {
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
     * comparator for descending sort startIdx1,
     * then startIdx2
     */    
    private class SequenceComparator4 implements
        Comparator<Sequence> {

        @Override
        public int compare(Sequence o1, Sequence o2) {
            
            if (o1.startIdx1 < o2.startIdx1) {
                return -1;
            } else if (o1.startIdx1 > o2.startIdx1) {
                return 1;
            }
            
            if (o1.startIdx2 < o2.startIdx2) {
                return -1;
            } else if (o1.startIdx2 > o2.startIdx2) {
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
            log.info(String.format("block=%d md[%d]", r, jOffset));
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
                
                log.info(String.format(" [%2d,%2d<-%2d] => %.4f", 
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
                
                // note, idx from q is i + iOffset
                count++;        
                sum += absS1;
                if (absS1 < Math.abs(mins[i])) {
                    int idx2 = i + jOffset;
                    if (idx2 >= n1) {
                        idx2 -= n1;
                    }
                    mins[i] = s1;
                    idxs0[i] = jOffset;
                    
                    // fill in the rest of the diagonal in this block
                    for (int k = (i-1); k > (i-r); k--) {
                        if (k < 0) {
                            break;
                        }
                        if (absS1 < Math.abs(mins[k])) {
                            idx2 = k + jOffset;
                            if (idx2 >= n1) {
                                idx2 -= n1;
                            }
                            mins[k] = s1;
                            idxs0[k] = jOffset;
                        }
                    }
                }
            }
            if (count == 0) {
                sum = Integer.MAX_VALUE;
            }
            log.info(String.format(
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
                if (idx2 >= n2) {
                    idx2 -= n2;
                }
                
                output.add(i, idx2, s1);
                
                // fill in the rest of the diagonal in this block
                for (int k = (i-1); k > (i-r); k--) {
                    if (k < 0) {
                        break;
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
                    output.add(k, idx2, s1);
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
        List<Sequence> existingList, Sequence sTest,
        int n1) {

        for (Sequence s : existingList) {
            if (intersects(s, sTest, n1)) {
                return true;
            }
        }
        
        return false;
    }
    
    private boolean intersects(Sequence s, Sequence sTest,
        int n1) {
        
        int sTestStopIdx1 = sTest.startIdx1 +
            (sTest.stopIdx2 - sTest.startIdx2);
         
        int stopIdx1 = s.startIdx1
            + (s.stopIdx2 - s.startIdx2);
        
        // aggregation does not extend idx2 past n2-1
        if (sTestStopIdx1 < n1 && stopIdx1 < n1) {
            return intersectsNoWrapAroundCheck(s, sTest);
        }

        if (sTestStopIdx1 > (n1 - 1)) {
            
            sTestStopIdx1 -= n1;
            
            if (stopIdx1 > (n1 - 1)) {
            
                stopIdx1 -= n1;
                
                /*
                s    : 0-->lastPoint     firstPoint   n-1   (lastPoint overrun)
                
                sTest: 0-->lastPoint     firstPoint   n-1   (lastPoint overrun)
                */

                // test that they don't intersect
                if (sTestStopIdx1 > stopIdx1 && sTest.startIdx1 < s.startIdx1) {
                    if ((sTest.stopIdx2 < s.startIdx2)
                        || (sTest.startIdx2 > s.stopIdx2)) {
                        return false;
                    }
                }
                
            } else {
                /*
                s    :              firstPoint  lastPoint   n-1
                
                sTest: 0-->lastPoint                     firstPoint   n-1   (lastPoint overrun)
                */ 
                // test that they do not intersect
                if (sTestStopIdx1 < s.startIdx1 && sTest.startIdx1 > 
                    stopIdx1) {
                    if ((sTest.stopIdx2 < s.startIdx2)
                        || (sTest.startIdx2 > s.stopIdx2)) {
                        return false;
                    }
                }
            }
            
        } else if (stopIdx1 > (n1 - 1)) {
            
            stopIdx1 -= n1;
            
            /*
            s    : 0-->lastPoint                     firstPoint   n-1   (lastPoint overrun)
                
            sTest:              firstPoint  lastPoint 
            */
            // test that they do not intersect
            if (stopIdx1 < sTest.startIdx1 && s.startIdx1 > 
                sTestStopIdx1) {
                if ((sTest.stopIdx2 < s.startIdx2)
                    || (sTest.startIdx2 > s.stopIdx2)) {
                    return false;
                }
            }
        }
        
        return true;
    }
    
    private boolean intersectsNoWrapAroundCheck(
        Sequence s, Sequence sTest) {
        
        int sTestStopIdx1 = sTest.startIdx1 +
            (sTest.stopIdx2 - sTest.startIdx2);
         
        int stopIdx1 = s.startIdx1
            + (s.stopIdx2 - s.startIdx2);

        // test that first indexes don't intersect
        if ((sTestStopIdx1 < s.startIdx1)
            || (sTest.startIdx1 > stopIdx1)) {
            // test that second indexes don't intersect
            if ((sTest.stopIdx2 < s.startIdx2)
                || (sTest.startIdx2 > s.stopIdx2)) {
                return false;
            }
        }
        
        return true;
    }
}
