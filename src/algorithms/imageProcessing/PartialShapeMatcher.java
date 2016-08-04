package algorithms.imageProcessing;

import algorithms.compGeometry.LinesAndAngles;
import algorithms.util.PairIntArray;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

/**
 *
 * @author nichole
 */
public class PartialShapeMatcher {
 
    /*
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
       - inclused calculation og a global PAreto frontier for
         measuring partial simulatiry
    
       * point sampling:
         (a) same number of points over each contour
             - can handle similarity transforms, but not occlusion
         (b) OR, equidistant points
             - can handle occlusion, but not scale
       ** equidistant is used

       - a shape is defined as the CW ordered sequence of points P_1...P_N
         and the shape to match has points Q_1...Q_N
       - descriptor invariant to translation, rotation, and scale:
           - a chord is a line joining 2 region points
           - uses the relative orientation between 2 chords
             angle a_i_j is from chord P_i_P_j to reference
             point P_i
             to another sampled point and chord P_j_P_(j-d) and P_j

             d is the number of points before j in the sequence of points P.

             a_i_j is the angle between the 2 chords P_i_P_j and P_j_P_(j-d)
       - regarding the descriptor matrix:
          the angles for on reference point P_1 w.r.t N sample points
          (and a fixed d such as d=3)
          would be a_i_1...a_i_N
          then, can do the same for each reference point to
          create a matrix w/ each reference point's angles as a row
             | a_1_1...a_1_N |
             | a_2_1...a_2_N |
               ...
             | a_N_1...a_N_N |
                 and elements on the diagonal are zero

                 note that the local information is close to
                 the diagonal

             the choice of the first point in the sequence
             affects the descriptors and hence comparisons.
    
             to change to a different first point, one can transform
             the matrix by circular shifts.
             for example, if k is chosen as first in sequence,
             would need to shift matrix up by k-1 rows and then
             left by k-1 rows.

       - 2 shape contours R_1 and R_2 w/ matrixes A_1, A_2
         of sizes NXN and MXM and M <= N

         to compare the matrices to find a partial match:
            - find rxr sized blocks similar to one another
              by starting at main diagonal element
              A_1(s,s) and A_2(m,m)
              which have a small angular difference value
 
                             1
              D_a(s,m,r) = ---- * summation_i_0_to_(r-1)
                            r^2
                                * summation_j_0_(r-1)
                                   of [A_1(s+i,s+j) - A_2(m+i,m+j)]^2
              //s range 0 to N-1
              //m range 0 to M-1
              //r range 0 or 1 to min(N, M)

              finding all of the possible matches and lengths
                is infeasible with brute force for large number of points,
                so authors created a Summary and Table approach
            Summary And Table:
                 - to calculate all D_a(s,m,r)
                   uses concept of "integral image"
                       by Viola and Jones
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

                   the final set is made by merging
                   different matches to make a set of connected
                   point correspondences
                   (see scissors in figure 4).

              Similarity of Shape
                  - adopting method of Bronstein etal. [17]. 
                    They proposed to use the Pareto-framework for 
                    quantitative interpretation of partial similarity. 
                    They define two quantities: 
                       partiality λ(X′, Y′), which describes 
                       the length of the parts
                       (the higher the value the smaller the part).
                       it's the size of the parts compared to the
                       entire object, so is 0 for partiality of the "whole".
                         X'^c = X\X'
                    and dissimilarity ε(X′, Y′), 
                       which measures the dissimilarity between the parts, 
                       
                       where X′ and Y′ are two contour parts of the shape. 

                    A pair Φ(X∗, Y ∗) = (λ(X∗, Y ∗), ε(X∗, Y ∗)) of 
                    partiality and dissimilarity values, 
                    where dissim is lowest for given partiality,
                    is the Pareto optimum.
                         the Pareto optimums are curves called the 
                         Pareto frontier.
                         they calculate a global optimal Pareto frontier 
                         by returning the minimum descriptor difference 
                         for all partialities.
                         Then the single partial similarity value,
                         a.k.a.  Salukwadze distance is calculated 
                         based on the Pareto frontier by
                           d_S(X,Y) =  inf  | phi(X*, Y*)|_1
                                      (X*,Y*)

                    where |..|_1 is the L1 norm of the vector
                           and L1 norm is the "manhattan dist",
                           that is abs(diffX) + abs(diffY), and
                           is not the euclidean distance

                    so d_X(X,Y) measures dist from (0,0) to
                        the closest point on the Paretto Frontier

            runtime complexity is based on image integral analysis
               and is O(m*n) where n and m are the number of sampled
               points on the input shapes.

               for their shape database tests, they found n=30 was
               enough to find matches.

               Integral images: 
                  "Rapid object detection using a boosted 
                   cascade of simple features" by
                    Viola and Jones 2001
               Pareto Frontier:
                  "Partial similarity of objects, or how to 
                  compare a centaur to a horse"
                   by Bronstein et al. 2008

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
    
    public void overrideSamplingDistance(int d) {
        this.dp = d;
    }
    
    /**
     * NOT READY FOR USE.
     * 
      the spacings used are equidistant, so note that
        any scale factor between p and q has to be 
        applied to the data before using this method.
     this method will return a score when completed.
      
     * @param p
     * @param q 
    */
    public double match(PairIntArray p, PairIntArray q) {
       
        int minN = Math.min(p.getN(), q.getN());
       
        //System.out.println("a1:");
        float[][] a1 = createDescriptorMatrix(p, minN);
        
        //System.out.println("a2:");
        float[][] a2 = createDescriptorMatrix(q, minN);       
        
        // --- make difference matrices ---
        
        // the rotated matrix for index 0 rotations is q.  
        // TODO: the same needs to be done for p separately
        // and combine results.
        // the scissors test shows this for points 17-32.
        // the articulated solution needs one of the matrices
        // in md[index0] to start at point 17 so that the summed
        // differences don't include the large difference in
        // transition from points 15 to 17.
        // In summary, need to combine results for these
        // for operations with p to q with results for q to p
        float[][][] md = createDifferenceMatrices(a1, a2, minN);
        
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
        
        List<Sequence> sequencesPQ = extractSimilar(md);
        
        return matchArticulated(sequencesPQ, minN);
        
        //printing out results for md[0] and md[-3] and +3
        // to look at known test data while checking logic
        //print("md[0]", md[0]);
        //print("md[3]", md[3]);  // <----- can see this is all zeros as expected
        //print("md[-3]", md[md.length - 1]);
        
    }
    
    protected List<Sequence> extractSimilar(float[][][] md) {

        /*
        - choose reasonable upper limit to r, such as sqrt(n)
            - for that rMax,
              - find the best blocks and find the
                avg diff of those and eiher the min and max of reange
                or st.dev.
       - visit all blocks as currently do and find the best,
         which may be more than one index of diff matrices,
         within the avg and stdev.
         (note, can build a std dev summed area table for this image alone)
       - inspect the best results:
         - form the data as nPart / nWhole for chains of similar indexes?
         - (still reading paretto grontier analysis, but i think that's what it
           is composed of)
        */
        
        int n1 = md[0].length;
       
        int rMax = (int)Math.sqrt(n1);
        if (rMax < 1) {
            rMax = 1;
        }

        double thresh = 30;
                
        MinDiffs mins = new MinDiffs(n1);
        for (int r = 1; r <= rMax; ++r) {
            findMinDifferenceMatrix(md, r, thresh, mins);
        }
        
        double tolerance = 3.;
        
        DiffMatrixResults equivBest = new DiffMatrixResults(n1);
        for (int r = 1; r <= rMax; ++r) {
            findEquivalentBest(md, r, mins, tolerance,
                equivBest);
        }
        
        /*
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
            System.out.println(sb.toString());
            sb.delete(0, sb.length());
        }
        */
        
        List<Sequence> sequences = new ArrayList<Sequence>();
        for (int idx1 = 0; idx1 < n1; ++idx1) {
            TIntList list = equivBest.indexes[idx1];
            if (list == null) {
                continue;
            }
            list.sort();
            for (int j = 0; j < list.size(); ++j) {
                int idx2 = list.get(j);
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
                        list2.remove(idx3);
                        set2.remove(idx3);
                    } else {
                        break;
                    }
                    
                    nextLIdx++;
                }
                if (s.stopIdx2 - s.startIdx2 > 1) {                
                    sequences.add(s);
                }
            }
        }
                
        System.out.println(sequences.size() + " sequences");
        for (int i = 0; i < sequences.size(); ++i) {
            Sequence s = sequences.get(i);
            int len = s.stopIdx2 - s.startIdx2 + 2;
            float frac = (float)len/(float)n1;
            System.out.println(String.format(
            "seq %d:%d to %d  %.4f", s.startIdx1, s.startIdx2,
                s.stopIdx2, frac));
        }
        
        return sequences;     
    }
    
    protected double matchArticulated(List<Sequence> srquences,
        int n1) {
        
        throw new UnsupportedOperationException("not yet implemented");
    }
    
    protected double matchRigidWithOcclusion(List<Sequence> srquences,
        int n1) {
        
        throw new UnsupportedOperationException("not yet implemented");
    }
    
    // index0 is rotations of p.n, index1 is p.n, index2 is q.n
    protected float[][][] createDifferenceMatrices(
        float[][] a1, float[][] a2, int minN) {
     
        /*
        | a_1_1...a_1_N |
        | a_2_1...a_2_N |
               ...
        | a_N_1...a_N_N |
           elements on the diagonal are zero
        
           to shift to different first point as reference,
           can shift down k-1 rows and left k-1 columns.
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
          //m range 0 to N-1, M<=N
          //r range 1? to min(N, M)

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
        (3) starting on the diagonals of the integral images
            made from the N M_D^n matrices,
            D_α(s, m, r) can be calculated for every block of any 
            size starting at any point on the diagonal in 
            constant time.
      
        reading the pareto frontier papers...
        lower threshold...
        building correspondence list from M_D^n...
        
        */

        // --- make difference matrices ---
        float[][][] md = new float[minN][][];
        float[][] prevA2Shifted = null;
        for (int i = 0; i < md.length; ++i) {
            float[][] shifted2;
            if (prevA2Shifted == null) {
                shifted2 = copy(a2);
            } else {
                // shifts by 1 to left and down by 1
                rotate(prevA2Shifted);
                shifted2 = prevA2Shifted;
            }
            //M_D^n = A_1(1:M,1:M) - A_2(n:n+M-1,n:n+M-1)
            md[i] = subtract(a1, shifted2);
            prevA2Shifted = shifted2;
        }
        
        // ---- make summary area table for md-----
        for (int i = 0; i < md.length; ++i) {
            applySummedAreaTableConversion(md[i]);
        }
        
        System.out.println("md.length=" + md.length);
        
        return md;
    }
    
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
        
        System.out.println("n=" + n);
        
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
                
                //System.out.println("i1=" + i1 + " imid=" + imid + " i2=" + i2);
   
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
                System.out.println(str);
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
        
        assert(a1.length <= a2.length);
        assert(a1[0].length <= a2[0].length);
        
        float[][] output = new float[a1.length][];
        for (int i = 0; i < a1.length; ++i) {
            output[i] = new float[a1[i].length];
            for (int j = 0; j < a1[i].length; ++j) {
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
            System.out.println(sb.toString());
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

    private class DiffMatrixResults {
        private TIntSet[] indexSets = null;
        private TIntList[] indexes = null;
        public DiffMatrixResults(int n) {
            indexes = new TIntList[n];
            indexSets = new TIntSet[n];
        }
        public void add(int index, int value) {
            if (indexes[index] == null) {
                indexes[index] = new TIntArrayList();
                indexSets[index] = new TIntHashSet();
            }
            if (!indexSets[index].contains(value)) {
                indexSets[index].add(value);
                indexes[index].add(value);
            }
        }
    }
    
    private class Sequence {
        int startIdx1;
        int startIdx2 = -1;
        int stopIdx2 = -1;
    }
    
    private class MinDiffs {
        int[] idxs1;
        int[] idxs2;
        float[] mins;
        public MinDiffs(int n) {
            idxs1 = new int[n];
            idxs2 = new int[n];
            mins = new float[n];
            Arrays.fill(idxs1, -1);
            Arrays.fill(idxs2, -1);
            Arrays.fill(mins, Float.MAX_VALUE);
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
        
        double c = 1./(double)(r*r);
       
        int n1 = md[0].length;

        /*
        md[n2offset][n1][n2]
        
        md[0   ][0   ][0   ]
          [n1-1][0-n1][n2-1]
        */
        
        int[] idxs0 = output.idxs1;
        int[] idxs2 = output.idxs2;
        float[] mins = output.mins;
     
        for (int iOffset = 0; iOffset < md.length; iOffset++) {
            System.out.println("md[" + iOffset + "]:");
            float[][] a = md[iOffset];
            float sum = 0;
            for (int i = 0; i < a.length; i+=r) {
                float s1;
                if ((i - r) > -1) {
                    s1 = a[i][i] - a[i - r][i] - a[i][i - r] + a[r][r];
                    System.out.println(
                        String.format(
                        " [%d,%d] %.4f, %.4f, %.4f, %.4f => %.4f", 
                        i, i, a[i][i], a[i - r][i], a[i][i - r],
                        a[r][r], s1*c));
                } else {
                    s1 = a[i][i];
                    System.out.println(
                        String.format(
                        " [%d,%d] %.4f => %.4f", 
                        i, i, a[i][i], s1*c));
                }
                s1 *= c;
                float absS1 = s1;
                if (absS1 < 0) {
                    absS1 *= -1;
                }
                if (absS1 > threshold) {
                   continue;
                }
                
                // note, idx from q is i + iOffset
                
                sum += absS1;
                if (absS1 < Math.abs(mins[i])) {
                    int idx2 = i + iOffset;
                    if (idx2 > n1) {
                        idx2 -= n1;
                    }
                    mins[i] = s1;
                    idxs0[i] = iOffset;
                    idxs2[i] = idx2;
                    
                    // fill in the rest of the diagonal in this block
                    for (int k = (i-1); k > (i - r); k--) {
                        if (k < 0) {
                            break;
                        }
                        if (mins[i] < mins[k]) {
                            idx2 = k + iOffset;
                            if (idx2 > n1) {
                                idx2 -= n1;
                            }
                            mins[k] = s1;
                            idxs0[k] = iOffset;
                            idxs2[k] = idx2;
                        }
                    }
                }
            }
            System.out.println("SUM=" + sum);
        }
        
        System.out.println("OFFSETS=" + Arrays.toString(idxs0));
        System.out.println("idx2=" + Arrays.toString(idxs2));
      
    }

    private void findEquivalentBest(float[][][] md, int r, 
        MinDiffs mins, double tolerance,
        DiffMatrixResults output) {
    
        int n1 = mins.idxs1.length;
         
        double c = 1./(double)(r*r);
               
        // capture all "best" within minSigns[i] += 2*variances[i]
        for (int iOffset = 0; iOffset < md.length; iOffset++) {
            float[][] a = md[iOffset];
            for (int i = 0; i < a.length; i+=r) {
                if (mins.idxs1[i] == -1) {
                    continue;
                }
                double s1;
                if ((i - r) > -1) {
                    s1 = a[i][i] - a[i - r][i] - a[i][i - r] + a[r][r];
                } else {
                    s1 = a[i][i];
                }
                s1 *= c;
                
                double avg = mins.mins[i];
    
                if (Math.abs(s1 - avg) > tolerance) {
                    continue;
                }
                
                int idx2 = iOffset + i;
                if (idx2 > n1) {
                    idx2 -= n1;
                }
                
                output.add(i, idx2);
                
                // fill in the rest of the diagonal in this block
                /*for (int k = (i-1); k > (i - r); k--) {
                    if (k < 0) {
                        break;
                    }
                    if ((k - r) > -1) {
                        s1 = a[k][k] - a[k - r][k] - a[k][k - r] 
                            + a[r][r];
                    } else {
                        s1 = a[k][k];
                    }
                    s1 *= c;
                    if (Math.abs(s1 - avg) > tolerance) {
                        continue;
                    }
                    idx2 = iOffset + k;
                    if (idx2 > n1) {
                        idx2 -= n1;
                    }
                    output.add(k, idx2);
                }*/
            }
        }
    }
    
}
