package algorithms.imageProcessing;

import algorithms.compGeometry.LinesAndAngles;
import algorithms.util.PairIntArray;
import java.util.Arrays;

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
                
        if (p.getN() > q.getN()) {
            throw new IllegalArgumentException(
            "q must be <= p in length.");
        }
       
        // --- make difference matrices ---
        float[][][] md = createDifferenceMatrices(p, q);
        
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
        
        return matchArticulated(md);
                
        //printing out results for md[0] and md[-3] and +3
        // to look at known test data while checking logic
        //print("md[0]", md[0]);
        //print("md[3]", md[3]);  // <----- can see this is all zeros as expected
        //print("md[-3]", md[md.length - 1]);
        
    }
    
    protected double matchArticulated(float[][][] md) {
       
        int rMax = (int)Math.sqrt(md[0].length);
        if (rMax < 1) {
            rMax = 1;
        }
        
        int[][] mins = new int[rMax][];
        for (int r = 1; r <= rMax; ++r) {
            mins[r - 1] = findMinDifferenceMatrix(md, r);
        }
        
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < md[0].length; ++i) {
            sb.append(String.format("[%4d]: ", i));
            for (int r = 0; r < mins.length; ++r) {
                String str;
                if (mins[r][i] == -1) {
                    str = " NA ";
                } else {
                    str = String.format("%3d", mins[r][i]);
                }
                sb.append(str).append(" ");
            }
            System.out.println(sb.toString());
            sb.delete(0, sb.length());
        }
        
        //printBlocks("md[0]", md[0]);
        //printBlocks("md[3]", md[3]);
        //printBlocks("md[-3]", md[md.length - 4]);
        
        throw new UnsupportedOperationException("not yet implemented");
    }
    
    protected double matchRigidWithOcclusion(float[][][] md) {
        
        printBlocks("md[0]", md[0]);
        printBlocks("md[3]", md[3]);
        printBlocks("md[-3]", md[md.length - 4]);
        
        throw new UnsupportedOperationException("not yet implemented");
    }
    
    protected float[][][] createDifferenceMatrices(
        PairIntArray p, PairIntArray q) {
                
        if (p.getN() > q.getN()) {
            throw new IllegalArgumentException(
            "q must be <= p in length.");
        }
        
        /*
        | a_1_1...a_1_N |
        | a_2_1...a_2_N |
               ...
        | a_N_1...a_N_N |
           elements on the diagonal are zero
        
           to shift to different first point as reference,
           can shift up k-1 rows and left k-1 columns.
        */
        
        float[][] a1 = createDescriptorMatrix(p);
        
        float[][] a2 = createDescriptorMatrix(q);
        
        int n1 = a1.length;
        
        int n2 = a2.length;
        
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
                         1
          D_a(s,m,r) = ---- * summation_i_0_to_(r-1)
                        r^2
                            * summation_j_0_(r-1)
                               of [A_1(s+i,s+j) - A_2(m+i,m+j)]^2
          //s range 0 to M-1
          //m range 0 to N-1, M<=N
          //r range 1? to min(N, M)

        Looking at D_a in detail for sets of s,m,r to 
        better understand the integral images of the
        difference matrices:

        s=0,m=0,r=2:
        (1/4) * ((A_1(0,0) - A_2(0,0)) 
               + (A_1(0,1) - A_2(0,1))
               + (A_1(0,2) - A_2(0,2)) +
                 (A_1(1,0) - A_2(1,0)) 
               + (A_1(1,1) - A_2(1,1))
               + (A_1(1,2) - A_2(1,2)) +
                 (A_1(2,0) - A_2(2,0)) 
               + (A_1(2,1) - A_2(2,1))
               + (A_1(2,2) - A_2(2,2)))
        (1/4) * ((A_1(0,0) + A_1(0,1) + A_1(0,2))
               + (A_1(1,0) + A_1(1,1) + A_1(1,2))
               + (A_1(2,0) + A_1(2,1) + A_1(2,2))
               - (A_2(0,0) + A_2(0,1) + A_2(0,2))
               - (A_2(1,0) + A_2(1,1) + A_2(1,2))
               - (A_2(2,0) + A_2(2,1) + A_2(2,2))
        (1/4) * (- A_1INT(1,1) + A_1INT(1,2)
                 + A_1INT(2,1) + A_1INT(2,2))
               - (- A_2INT(1,1) + A_2INT(1,2)
                  + A_2INT(2,1) + A_2INT(2,2))
        can use the summary area table method to make
        array summed along x and y.
        from wikipedia: "the value at any point (x, y) 
        in the summed area table is just the sum of all 
        the pixels above and to the left of (x, y), inclusive"
  
        I(x,y)=i(x,y)+I(x-1,y)+I(x,y-1)-I(x-1,y-1))
                        sum x      sum y
        2 7  6  4      7 13 17   14 18 30 
        1 3  0  3      3  0  6    7  5 13
        0 4  1  2      4  5  7    4  5  7
          0  1  2             
        
        s=0,m=0,r=2:
            D_a(s,m,r) = (1/4) * (- A_1INT(1,1) + A_1INT(1,2)
                                  + A_1INT(2,1) + A_1INT(2,2))
                               - (- A_2INT(1,1) + A_2INT(1,2)
                                  + A_2INT(2,1) + A_2INT(2,2))
        s=0,m=1,r=2:
            D_a(s,m,r) = (1/4) * (- A_1INT(1,1) + A_1INT(1,2)
                                  + A_1INT(2,1) + A_1INT(2,2))
                               - (- A_2INT(2,1) + A_2INT(2,2)
                                  + A_2INT(3,1) + A_2INT(3,2))
        
        rewritten using the integral images:
        D_a(s,m,r) = 
            (1/r^2) * (- A_1INT(s+r-1,s+r-1) + A_1INT(s+r-1,s+r)
                       + A_1INT(s+r,  s+r-1) + A_1INT(s+r,  s+r))
                    - (- A_2INT(m+r-1,m+r-1) + A_2INT(m+r-1,m+r)
                       + A_2INT(m+r,  m+r-1) + A_2INT(m+r,  m+r))
        
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
        float[][][] md = new float[n1][][];
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
        
        // ---- make summary area tables from md -----
        for (int i = 0; i < md.length; ++i) {
            float[][] mdI = md[i];
            // I(x,y)=i(x,y)+I(x-1,y)+I(x,y-1)-I(x-1,y-1))
            for (int x = 0; x < mdI.length; ++x) {
                for (int y = 0; y < mdI[x].length; ++y) {
                    if (x > 0 && y > 0) {
                        mdI[x][y] += (mdI[x-1][y] + mdI[x][y-1]
                            - mdI[x-1][y-1]);
                    } else if (x > 0) {
                        mdI[x][y] += mdI[x-1][y];
                    } else if (y > 0) {
                        mdI[x][y] += mdI[x][y-1];
                    }
                }
            }
        }
       
        System.out.println("md.length=" + md.length);
        
        return md;
    }
    
    protected float[][] createDescriptorMatrix(PairIntArray p) {
        
        int n = (int)Math.ceil((double)p.getN()/dp);
        
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
    
    /**
     * 
     * @param md 3 dimensional array of difference matrices
     * @param r block size
     * @return 
     */
    private int[] findMinDifferenceMatrix(float[][][] md,
        int r) {
        
        double c = (1./(r*r));
        
        int[] idxs = new int[md[0].length];
        double[] mins = new double[md[0].length];
        Arrays.fill(idxs, -1);
        Arrays.fill(mins, Double.MAX_VALUE);
        
        for (int i0 = 0; i0 < md.length; i0++) {
            System.out.println("md[" + i0 + "]:");
            float[][] a = md[i0];
            for (int i = 0; i < a.length; i+=r) {
                double d;
                if ((i - r) > -1) {
                    d = a[i][i] - a[i - r][i] - a[i][i - r] + a[r][r];
                    System.out.println(
                        String.format(
                        " [%d,%d] %.4f, %.4f, %.4f, %.4f => %.4f", 
                        i, i, a[i][i], a[i - r][i], a[i][i - r],
                        a[r][r], d*c));
                } else {
                    d = a[i][i];
                }
                if (d < 0) {
                    d *= -1;
                }
                if (d > 50) {
                    continue;
                }
                d *= c;
                if (d < mins[i]) {
                    mins[i] = d;
                    idxs[i] = i0;
                }
            }
        }
        
        return idxs;
    }

    private void printBlocks(String label, float[][] a) {

        System.out.println(label);
        
        /*
        Find quadratic sub-areas within reference and query 
        descriptor matrices starting at main diagonal which 
        are most similar.
        
        */
        // try a set block size
        int r = 2;
        
        double c = (1./(r*r));
        
        double prev = Double.MAX_VALUE;
        
        for (int i = r; i < a.length; i+=r) {
            
            double d;
            
            if ((i - r) > -1) {
                d = c * (a[i][i] - a[i - r][i - r] +
                    a[i - r][i] + a[i][i - r]);
                System.out.println(
                    String.format(
                    " [%d,%d] %.4f, %.4f, %.4f, %.4f => %.4f", 
                    i, i, -a[i - r][i - r],
                    a[i - r][i], a[i][i - r],
                    a[i][i], d));
            } else {
                d = c * a[i][i];
            }
            
            /*
            d = a[i][i];
            if (prev < Double.MAX_VALUE) {
                d -= prev;
            }
            prev = d;
                    
            System.out.println(
                String.format(" [%d,%d]%.4f,", i, i, d));
            */
        }
        
    }

}
