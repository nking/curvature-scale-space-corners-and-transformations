package algorithms.imageProcessing;

import algorithms.QuickSort;
import algorithms.compGeometry.LinesAndAngles;
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
        global solution and separately for best local
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
        
        List<Sequence> sequences = new ArrayList<Sequence>();
        List<Sequence> discarded = new ArrayList<Sequence>();
                
        extractSimilar(md, sequences, discarded);
       
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
                
        Sequences sequences0 = matchArticulated(
            sequences, n1, n2);
        
        //addFeasibleDiscarded(sequences0, discarded);
        
        if (diffN <= 0) {
            return sequences0;
        }
        
        transpose(sequences0, n1, n2);
        
        return sequences0;
    }
    
    protected void extractSimilar(float[][][] md,
        List<Sequence> sequences, List<Sequence> discarded) {
        
        //md[0:n2-1][0:n1-1][0:n1-1]
        
        int n2 = md.length;
        int n1 = md[0].length;
       
        int rMax = (int)Math.sqrt(n1);
        if (rMax < 1) {
            rMax = 1;
        }
     
        // 23 degrees is 0.4014
        double thresh = 23. * Math.PI/180.;
        
        /*
        TODO: will apply a different pattern of reading
        the blocks and merging results next.
        */
    
        MinDiffs mins = new MinDiffs(n1);
        for (int r = 2; r <= rMax; ++r) {
            findMinDifferenceMatrix(md, r, thresh, mins);
        }
       
        // 10 degrees is 0.175
        double tolerance = 0.25;
  
        DiffMatrixResults equivBest = new DiffMatrixResults(n1);
        for (int r = 2; r <= rMax; ++r) {
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
                        
                        // NOTE: an expensive operation which will be
                        // replaced with new block pattern reading:
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
                            "seq %d:%d to %d  frac=%.4f  avg diff=%.4f",
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
        int n1, int n2) {
        
        //(1) ascending sort ordered by startIdx1 and then
        // descending fraction of whole 
        Collections.sort(sequences, new SequenceComparator());
        
        // (1.5) descending sort of fraction, then diff, then startIdx
        List<Sequence> list2 = new ArrayList<Sequence>(sequences);
        float maxAvgDiff = findMaxAvgDiff(list2);
        System.out.println("list2.sz=" + list2.size());
        Collections.sort(list2, new SequenceComparator3(maxAvgDiff));
        //Collections.sort(list2, new SequenceComparator2());

        //(2) a lookup for items in list2
        //    belonging to >= startIdx1.
        TreeMap<Integer, TIntList> startLookup =
            new TreeMap<Integer, TIntList>();
        for (int i = 0; i < list2.size(); ++i) {
            Sequence s1 = list2.get(i);
            Integer key = Integer.valueOf(s1.startIdx1);
            for (int j = 0; j < list2.size(); ++j) {
                Sequence s2 = list2.get(j);
                if (s2.startIdx1 < key.intValue()) {
                    continue;
                }
                if (!startLookup.containsKey(key)) {
                    startLookup.put(key, new TIntArrayList());
                }                
                // TODO: revisit this...avoiding adding entire list
                if (startLookup.get(key).size() > n1/4) {
                    break;
                }
                startLookup.get(key).add(j);
            }
            log.fine("FSORT: " + i + " " + s1.toString());
        }

        //(2.5) a lookup for items in list2
        //    belonging to >= startIdx2.
        TreeMap<Integer, TIntList> startLookup2 =
            new TreeMap<Integer, TIntList>();
        for (int i = 0; i < list2.size(); ++i) {
            Sequence s1 = list2.get(i);
            Integer key = Integer.valueOf(s1.startIdx2);
            for (int j = 0; j < list2.size(); ++j) {
                Sequence s2 = list2.get(j);
                if (s2.startIdx2 < key.intValue()) {
                    continue;
                }
                if (!startLookup2.containsKey(key)) {
                    startLookup2.put(key, new TIntArrayList());
                }                
                // TODO: revisit this...avoiding adding entire list
                if (startLookup2.get(key).size() > n2/4) {
                    break;
                }
                startLookup2.get(key).add(j);
            }
        }
        
        // (3) create "tracks" of sequences
        Set<Sequence> added = new HashSet<Sequence>();
        
        List<Sequences> tracks = new ArrayList<Sequences>();
        
        for (int i = 0; i < sequences.size(); ++i) {
            
            Sequence s = sequences.get(i);
            if (added.contains(s)) {
                continue;
            }
            added.add(s);
            
            Sequences currentTrack = new Sequences();
            tracks.add(currentTrack);
            currentTrack.sequences.add(s);
        
            /*
            (next.startIdx1 >
                   startIdx1 + (stopIdx2 - stopIdx1))
            */
            Sequence lastSequence = s;
            while (true) {
                
                int nextStartIdx1 = lastSequence.startIdx1
                    + (lastSequence.stopIdx2 -
                    lastSequence.startIdx2) + 1;

                // if nextStartIdx1 is larger than n1, it will
                // be missing from startLookup and we've
                // ended the progression through idx1
                
                Entry<Integer, TIntList> entry = 
                    startLookup.ceilingEntry(
                        Integer.valueOf(nextStartIdx1));

                if (entry == null) {
                    break;
                }
               
                TIntList list2Indexes = entry.getValue();
                boolean appended = false;
                for (int j = 0; j < list2Indexes.size(); ++j) {
                    int list2Idx = list2Indexes.get(j);
                    Sequence s2 = list2.get(list2Idx);
                    if (s2.startIdx2 <= lastSequence.stopIdx2) {
                        continue;
                    }
                    /*
                    if (!verifyConsistenCW(currentTrack, s2)) {
                        continue;
                    }*/
                    currentTrack.sequences.add(s2);
                    added.add(s2);
                    appended = true;
                    lastSequence = s2;
                    break;
                }
                if (!appended) {
                    break;
                }
            }

            // if first sequence startIdx2 is > 0, need to
            // search the region before startIdx2 also
            if (currentTrack.sequences.get(0).startIdx2 > 1) {
                int nextStartIdx2 = 0;
                while (true) {
                    Entry<Integer, TIntList> entry = 
                        startLookup2.ceilingEntry(
                            Integer.valueOf(nextStartIdx2));

                    if (entry == null) {
                        break;
                    }

                    TIntList list2Indexes = entry.getValue();
                    boolean appended = false;
                    for (int j = 0; j < list2Indexes.size(); ++j) {
                        int list2Idx = list2Indexes.get(j);
                        Sequence s2 = list2.get(list2Idx);
                        if (s2.stopIdx2 >= 
                            currentTrack.sequences.get(0).startIdx2) {
                            continue;
                        }
                        
                        // check for range clash...
                        if (intersectsExistingRange1(currentTrack.sequences,
                            s2)) {
                            continue;
                        } 
                        
                        currentTrack.sequences.add(s2);
                        added.add(s2);
                        appended = true;
                        lastSequence = s2;
                        break;
                    }
                    if (!appended) {
                        break;
                    }

                    nextStartIdx2 = lastSequence.stopIdx2 + 1;
                    if (nextStartIdx2 > currentTrack.sequences.get(0).startIdx2) {
                        break;
                    }
                }
            }
        }
     
        filterForConsistentClockwise(tracks);
     
        // calculate the stats for each track (== Sequences)
        for (Sequences track : tracks) {
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
        Collections.sort(tracks, new TrackComparator(n1));
        for (int i = 0; i < tracks.size(); ++i) {
            Sequences track = tracks.get(i);
            log.info("track " + i + ": " + track.toString());
        }

        return tracks.get(0);
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

            if (sequences.sequences.isEmpty()) {
                rmList.add(i);
                continue;
            }

            Collections.sort(sequences.sequences,
                new SequenceComparator4());

            /*
             all startIdx1 should be increasing,
             and wrap around should be considered.
             then, all startIdx2 should be increasing
             and wrap around whould be considered.
             */
            Sequence s0 = sequences.sequences.get(0);

            boolean notValid = false;
            
            int ns = sequences.sequences.size();

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
                        s = sequences.sequences.get(0);
                    } else {
                        s = sequences.sequences.get(j);
                    }
                    int idx = (check == 0) ? s.startIdx1
                        : s.startIdx2;
                    if (idx == prev) {
                        rmList.add(i);
                        notValid = true;
                        break;
                    } else if (idx < prev) {
                        if (wrapped) {
                            rmList.add(i);
                            notValid = true;
                            break;
                        }
                        wrapped = true;
                        prev = idx;
                    }
                    prev = idx;
                } // end loop over j sequences in a track

                if (notValid) {
                    break;
                }
                
            } // end loop over check
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

    private boolean intersectsExistingRange1(
        List<Sequence> existingList, Sequence s) {

        int stopIdx1 = s.startIdx1 +
            (s.stopIdx2 - s.startIdx2);
        
        for (Sequence s0 : existingList) {
            
            int s0stopIdx1 = s0.startIdx1 +
                (s0.stopIdx2 - s0.startIdx2);
            if (s.startIdx1 >= s0.startIdx1 &&
                s.startIdx1 <= s0stopIdx1) {
                return true;
            }
            if (stopIdx1 >= s0.startIdx1 &&
                stopIdx1 <= s0stopIdx1) {
                return true;
            }
            if (s.startIdx1 <= s0.startIdx1 &&
                stopIdx1 >= s0stopIdx1) {
                return true;
            }
        }

        return false;
    }
}
