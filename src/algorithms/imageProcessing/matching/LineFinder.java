package algorithms.imageProcessing.matching;

import algorithms.compGeometry.LinesAndAngles;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import thirdparty.edu.princeton.cs.algs4.Interval;
import thirdparty.edu.princeton.cs.algs4.IntervalRangeSearch;

/**
NOTE: NOT READY FOR USE YET... still testing.

Adapted the chord descriptor difference matrix from PartialShapeMatcher.java 
* which was made following the paper
"Efficient Partial Shape Matching of Outer Contours" by Donoser

The differences in chords of the implied second shape, a line is always 0,
so edits are present specific to a one dimensional line instead of a closed
curve.

@author nichole
*/
public class LineFinder {

    /**
     * in sampling the boundaries of the shapes, one can
     * choose to use the same number for each (which can result
     * in very different spacings for different sized curves)
     * or one can choose a set distance between sampling
     * points.
     * dp is the set distance between sampling points.
       The authors of the paper use 3 as an example.
     */
    protected int dp = 1;

    // 10 degrees is 0.1745 radians
    // for a fit to a line, consider 1E-9
    private float thresh = (float)(1.e-7);

    protected Logger log = Logger.getLogger(this.getClass().getName());

    private boolean debug = false;

    /**
     * override the threshhold for using a chord differernce value
     * to this.   By default it is set to 1.e-7 radians.
     * @param t 
     */
    public void _overrideToThreshhold(float t) {
        this.thresh = t;
    }

    /**
     * override the sampling distance along the boundary which is by default
     * set to 1.  If a larger value is used, the curve is sampled at
     * the given spacing, and the results are interpolated as filled in
     * between a sampled range.
     * @param d 
     */
    public void overrideSamplingDistance(int d) {
        this.dp = d;
    }

    public void setToDebug() {
        debug = true;
        log.setLevel(Level.FINE);
    }

    /**
      NOT READY FOR USE... still testing...

      A shape is defined as the clockwise ordered sequence
      of points P_1...P_N.
      The spacings used within this method are equidistant
      The default spacing is 1, 
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
     @param p
    */
    public LineResult match(PairIntArray p) {

        log.fine("p.n=" + p.getN());

        if (p.getN() < 2) {
            throw new IllegalArgumentException("p must "
            + " have at least dp*2 points = " + (dp * 2));
        }

        if (dp == 1) {
            return match0(p);
        }
     
        PairIntArray pSub = new PairIntArray(p.getN()/dp);
    
        for (int i = 0; i < p.getN(); i += dp) {
            pSub.add(p.getX(i), p.getY(i));
        }
        
        log.fine("pSub.n=" + pSub.getN());
        
        LineResult rSub = match0(pSub);
        
        if (rSub == null) {
            return null;
        } 
        
        // -- put results back into frame of p --
        
        LineResult r = new LineResult();
        
        List<PairInt> lr = rSub.getLineIndexRanges();
        
        for (int i = 0; i < lr.size(); ++i) {
            PairInt startStop = lr.get(i);
            int x = startStop.getX() * dp;
            int y = startStop.getY() * dp;
            r.addLineRange(new PairInt(x, y));
        }
        
        return r;
    }

    private LineResult match0(PairIntArray p) {

        if (p == null || p.getN() < 2) {
            throw new IllegalArgumentException("p must have at "
                + "least 2 points");
        }
        
        // --- make difference matrices ---

        //md[0:n2-1][0:n1-1][0:n1-1]
        int n1 = p.getN();
        float[][] md = createDifferenceMatrices(p);
        
        // need to find the intervals of contiguous 0s
        // and assign member points to the largest segments.
        
        /*
        reading each row
           read start and stop cols of values < threshold
           and store each as an interval.
           if intersects with existing interval, compare size and
              keep the one that is larger in the interval tree
        
        then need a correction for wrap around if there are
           zeros at first and last index
        */
        
        IntervalRangeSearch<Integer, Integer> rangeSearch =
            new IntervalRangeSearch<Integer, Integer>();
        
        // key = map size at put, value = interval
        TIntObjectMap<Interval<Integer>> intervalMap = 
            new TIntObjectHashMap<Interval<Integer>>();
       
        Interval<Integer> interval = null;
        
        for (int i = 0; i < md.length; ++i) {
            
            int start = -1;
            
            for (int j = 0; j < md[i].length; ++j) {
                
                float d = md[i][j];
                
                boolean store = false;
                
                if (d < thresh) {
                    if (start == -1) {
                        // do not start if on diagonal
                        // TODO: revisit this soon
                        if (j == i) {
                            continue;
                        }
                        start = j;
                    }
                    if ((j == (md[i].length - 1)) && (start > -1)) {
                        store = true;
                    }
                } else if (start > -1) {
                    store = true;
                }
                
                if (store) {
                    
                    // create an interval,
                    int stop = j - 1;
                                        
                    // to prevent two intersecting lines from being merged
                    // into one, will use a start interval one
                    // index higher, and correct for it later.
                    // also, not storing single index matches
                    
                    start++;
                    if (start > (md[i].length - 1)) {
                        start = (md[i].length - 1);
                    }
                    
                    if (stop < start) {
                        // do not store single index matches
                        start = -1;
                        continue;
                    }       
                  
                    interval =  new Interval<Integer>(start, stop);
        
                    int sz = intervalMap.size();
                    
                    // store it in range search
                    Integer existing = rangeSearch.put(interval, 
                        Integer.valueOf(sz));
              
                    if (existing != null) {
                        // clashes with existing, so make sure the largest
                        // remains in range tree
                        Interval<Integer> comp = intervalMap.get(existing);
                        int nc = comp.max().intValue() - comp.min().intValue() + 2;
                        int ni = stop - start + 2;
                        if (nc > ni) {
                            Integer rmvd = rangeSearch.put(comp, existing);
                            assert(rmvd != null);
                            assert(rmvd.intValue() == sz);
                        } else {
                            intervalMap.put(Integer.valueOf(sz), interval);
                        }
                    } else {
                        intervalMap.put(Integer.valueOf(sz), interval);
                    }
                    
                    // reset startI 
                    start = -1;
                }
            }
        } // end loop j
                
        List<Interval<Integer>> list = rangeSearch.getAllIntrvals();
        
        TIntSet existing = new TIntHashSet();
        
        LineResult result = new LineResult();
        for (Interval<Integer> interval2 : list) {
            // correct for the interval start being +1
            int start = interval2.min() - 1;
            if (existing.contains(start)) {
                start++;
            }
            int stop = interval2.max();
            if (existing.contains(stop)) {
                stop--;
            }
            PairInt s = new PairInt(start, stop);
            result.addLineRange(s);
            for (int i = start; i <= stop; ++i) {
                existing.add(i);
            }
        }
        
        // TODO: add here a search to extend into gaps where possible.
        //    some of the steps above in trying to assure separated
        //    segments, have left some indexes of found lines out
        //    at the intersection of lines
                
        return result;
    }

    /**
     * create the matrices of differences between p
     * and q.  Note that the matrix differences are
     * absolute differences.
      returns a[0:p.n-1][0:p.n-1]
    */
    protected float[][] createDifferenceMatrices(
        PairIntArray p) {

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
        float[][] a2 = createLineDescriptorMatrix(p.getN());

        /*
            MXM              <XM
         20 21 22        20 21 22
         10 11 12        10 11 12 
         00 01 02        00 01 02   p_i_j - q_i_j
        */
 
        // --- make difference matrices ---
        int n1 = p.getN();
        int n2 = n1;
        float[][] md = copy(a2);
        // NOTE: absolute values are stored.
        //M_D^n = A_1(1:M,1:M) - A_2(n:n+M-1,n:n+M-1)
        md = subtract(a1, md);
        
        if (debug) {
            print("a1", a1); 
            print("a2", a2); 
            print("diff matrix", md);
        }
        
        return md;
    }
    
    protected float[][] createLineDescriptorMatrix(int n) {
        
        float[][] a = new float[n][];
        for (int i = 0; i < n; ++i) {
            a[i] = new float[n];
            Arrays.fill(a[i], (float)Math.PI);
            // diagonal is zero
            a[i][i] = 0;
        }
        
        return a;
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

                double angleA = LinesAndAngles
                    .calcClockwiseAngle(
                    p.getX(i1), p.getY(i1),
                    p.getX(i2), p.getY(i2),
                    p.getX(imid), p.getY(imid)
                );

                //System.out.println("i1=" + i1 + " imid=" + imid + " i2=" + i2 +
                //    "  angleA=" + angleA);
                
                a[i1][i2] = (float)angleA;
                
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

    private void print(String label, float[][][] a) {
        for (int i = 0; i < a.length; ++i) {
            print(label + " off " + i, a[i]);
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

    public static class LineResult {
        
        List<PairInt> lineIndexRanges = new ArrayList<PairInt>();
        
        public LineResult() {
        }
        
        /**
         * add start and stop index ranges for a found line segment
         * @param startStop 
         */
        public void addLineRange(PairInt startStop) {
            lineIndexRanges.add(startStop);
        }
        
        public List<PairInt> getLineIndexRanges() {
            return lineIndexRanges;
        }
    }
}
