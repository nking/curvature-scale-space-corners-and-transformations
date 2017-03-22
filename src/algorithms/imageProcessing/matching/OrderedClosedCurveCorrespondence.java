package algorithms.imageProcessing.matching;

import algorithms.util.PairIntArray;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import java.util.ArrayList;
import java.util.List;
import java.util.Map.Entry;
import java.util.TreeMap;

/**
 * a class to handle additions, that are sometimes merges 
 * of overlapping or
 * adjacent correspondence intervals, that are 
 * clockwise consistent.
 * The class performs checks
 * to the correspondency list as intervals are added (or attempted).
 * The class is specialized for use with PartialShapeMatcher.java
 * which has intervals it attempts to add to this structure in
 * order of increasing cost.
 * The intervals are given as ranges of indexes called idx1
 * and the shape they match to is specified as an offset
 * from the idx1 indexes (that is, the 2nd correspondence list is
 * implied by the offset from the first for a given interval range).
 * Each interval has a single offset which may be different
 * from the offset in other intervals
 * (because the PartialShapeMatchers allow occlusion and articulation).
 *
 * NOTE: this class is not "thread safe", that is, only a single thread should
 * access it because it uses an internal cache that is not guarded..
 * 
 * @author nichole
 */
class OrderedClosedCurveCorrespondence {

    // visit each interval in order of smallest cost,
    // and only add the clockwise consistent intervals to a combined output
    protected TreeMap<Integer, SR> t1 = new TreeMap<Integer, SR>();

    private int minLength = 3;

    // NOTE: this makes the code "not thread safe"
    private int[] cachedIdx2 = new int[2];
    
    private int nMatched = 0;
    
    private boolean doStopAt90Percent = true;
    
    // begin purely debug variables
    private boolean debug = false;
    public PairIntArray dbg1 = null;
    public PairIntArray dbg2 = null;
    public int dp = 1;
    // end purely debug variables
    
    public void setToDebug() {
        debug = true;
    }
    
    public void setMinimumLength(int length) {
        minLength = length;
    }
    
    public void overrideStopAt90PercentMatched() {
        this.doStopAt90Percent = false;
    }

    public void addIntervals(List<SR> intervals, int n1, int n2) {

        // rule from PartialShapeMatcher is n1 <= n2.
        float nMaxMatchable = n1;
        
        for (SR sr: intervals) {
            if (debug) {
                System.out.println("cost=" + sr.calcSalukDist() + 
                " sr=" + sr.startIdx1 + " : " + sr.stopIdx1);
            }
            
            addInterval(sr, n1, n2);
            
            if (doStopAt90Percent && (nMatched > 0.9 * nMaxMatchable)) {
                return;
            }
        }

    }

    public List<SR> getResultsAsList() {

        List<SR> list = new ArrayList<SR>();

        for (Entry<Integer, SR> entry : t1.entrySet()) {
            list.add(entry.getValue());
        }

        return list;
    }

    
   
    private void print(SR sr, String label, int n2) {
        if (debug && dbg1 != null) {
            calculateIds2s(sr, n2);
            System.out.println(label + String.format(
            "\n    -->add p: %d %d : (%d, %d) : (%d, %d) off=%d\n", 
            sr.startIdx1, sr.stopIdx1,
            dp*dbg1.getX(sr.startIdx1), dp*dbg1.getY(sr.startIdx1),
            dp*dbg1.getX(sr.stopIdx1), dp*dbg1.getY(sr.stopIdx1),
            sr.offsetIdx2)
            + String.format(
        "    idx2s: %d %d : (%d, %d) : (%d, %d) \n", 
            cachedIdx2[0], cachedIdx2[1],
            dp*dbg2.getX(cachedIdx2[0]), dp*dbg2.getY(cachedIdx2[0]),
            dp*dbg2.getX(cachedIdx2[1]), dp*dbg2.getY(cachedIdx2[1]))
            );
        }
    }

    private void addFirstInterval(SR sr) {

        if (!t1.isEmpty()) {
            throw new IllegalStateException("addFirstInterval "
                + " is meant for use with an empty tree");
        }

        Integer k1 = Integer.valueOf(sr.startIdx1);
        t1.put(k1, sr);
        
        nMatched += sr.mLen;
    }

    /**
     * add intervals to the clockwise ordered unique correspondence
     * list internal to this instance.
     * Note that each interval is expected to be clockwise consistent
     * (stopIdx1 > startIdx1) and the list of intervals is expected
     * to be sorted so that the highest priority (== lowest cost)
     * intervals are at the smallest list indexes, that is the
     * list is increasing in cost with index.
     * @param sr
     * @param n1
     * @param n2
     */
    public void addInterval(SR sr, int n1, int n2) {

        assert(sr.startIdx1 != sr.stopIdx1);

        if (t1.isEmpty()) {
            addFirstInterval(sr);
            if (debug) {
                print(sr, "first : ", n2);
            }
            return;
        }

        //assert clockwise consistent

        //NOTE: if an interval is trimmed rather than discarded here
        // because of clockwise consistency,
        // then might need to consider re-doing the interval sort...
        // (an adaptive optimal: if an interval is trimmed, might
        // edit it in allResults list, re-sort and start again...
        // will not do that here, but might consider a better way to have
        // same result in the future).

        /*
        ---------------------
        goal: to check that a new interval to insert into t1 is consistent
              with t1 existing indexes in idx1 and in idx2
              where consistency is clockwise ordering of both lists
              and unique matchings.

        first structures:
            t1 is an ordered tree map w/ key = sr.startIdx1 of interval sr
                and value = interval sr.

        NOTE: some defintions w.r.t. TreeMap are
            ceiling method returns a key-value mapping associated
                with the least key greater than or equal to the given key,
                or null if there is no such key.
            floor method returns a key-value mapping associated
                with the greatest key less than or equal to the given key,
                or null if there is no such key

        NOTE: to simplify the order checks of idx2, will add a phase
            to idx2 when idx2 < idx1.
            For example, let n1=n2=10, and one pair has idx1=2 w/ idx2=9
            then the next pair w/ idx1=3 maps to idx2=0,
            but to keep idx2 increasing, will add n2 to make it 10.

        NOTE: have added an exclusion clause that may need to be edited.
            If a candidate interval will be adjacent to an existing interval
            in t1 in terms of idx1, then idx2 must be adjacent also
             within a pixel or so.
            This is to prevent a large discontinuity.
        
        --------------------------------------------------------
        case 0: sr.startIdx1 ceiling is null, that is, there are
                no intervals in t1 at same or larger index position
                than st.startIdx1
                and there is at least 1 existing interval in t1.

            content ordered by idx1
                      t1  |  interval
               ------------------------
                          |
                      -#- | [-#-]
               startIdx1  | [sr]

        --------------------------------------------------------
        case 1: sr.startIdx1 ceiling is not null, that is, there are
                intervals in t1 at same or larger index position
                than st.startIdx1
                and there are not intervals in t1 at a
                smaller index position that st.startIdx1.

            content ordered by idx1
                      t1  |  interval
               ------------------------
                          |
               startIdx1  | [sr]
                      -#- | [-#-]

        --------------------------------------------------------
        case 2: sr.startIdx1 ceiling is not null, that is, there are
                intervals in t1 at same or larger index position
                than st.startIdx1
                and there are intervals in t1 at a
                smaller index position that st.startIdx1.

            content ordered by idx1
                      t1  |  interval
               ------------------------
                      -#- | [-#-]
               startIdx1  | [sr]
                      -#- | [-#-]

        */
        
        // possibly intersecting, so remove complete intersection,
        // or inconsistent intersection
        Entry<Integer, SR> midE = t1.ceilingEntry(
            Integer.valueOf(sr.startIdx1));
        if (isEmbeddedOrInconsistentWith(sr, midE, n2)) {
            return;
        }
       
        
        Entry<Integer, SR> above = t1.floorEntry(
            Integer.valueOf(sr.startIdx1 - 1));
        if (isEmbeddedOrInconsistentWith(sr, above, n2)) {
            return;
        }
        
        Entry<Integer, SR> below = t1.ceilingEntry(
            Integer.valueOf(sr.stopIdx1 + 1));
        
        if (isEmbeddedOrInconsistentWith(sr, below, n2)) {
            return;
        }
        
        if (sr.startIdx1 == 0) {

            if (below == null) {
                // this can happen if there's only one item in t1 and
                // sr has the same or smaller range than it.
                return;
            }
           
            // case 1

            addForCase1(sr, below, n2);

        } else if (below == null) {

            // case 0

            // no entries below sr are in t1

            assert(sr.startIdx1 > 0);
            
            if (above == null) {
                // this can happen if there's only one item in t1 and
                // sr has the same or smaller range than it.
                return;
            }

            addForCase0(sr, above, n2);

        } else {

            if (above == null) {

                // case 1

                // there are no entries above sr in t1

                addForCase1(sr, below, n2);

            } else {

                // case 2

                addForCase2(sr, above, below, n2);
            }
        }
    }

    private void addForCase0(SR sr, Entry<Integer, SR> above,
        int n2) {

        /*
        (1) find t1 floor for sr.startIdx1 - 1.
        (2) test that entire range is consistent
            -- if floor stopIdx1 is < sr.startIdx1
               -- if floor stopIdx2 is < sr.startIdx2
                  can add interval
        (3) iterate over each idx1,idx2 in sr interval
            test for each ifx1,idx2
            -- if idx1 > floor stopIdx1
               -- if idx2 > floor stopIdx2
                  can add interval
        */

        if (excludeCase0(above.getValue(), sr, n2)) {
            return;
        }
        
        if (case0AllConsistent(sr, above, n2)) {
            Integer k1 = Integer.valueOf(sr.startIdx1);
            t1.put(k1, sr);
            nMatched += sr.mLen;
            if (debug) {
                print(sr, "case 0", n2);
            }
            return;
        }

        SR floor = above.getValue();
        calculateIds2s(floor, n2);
        int floorStopIdx2 = cachedIdx2[1];

        TIntList subsetIdx1s = new TIntArrayList();

        populateCase0Idx1s(sr, floor.startIdx1, floorStopIdx2, 
            subsetIdx1s, n2);

        int ns = subsetIdx1s.size();
        if (ns < minLength) {
            return;
        }

        assert(assertContiguous(subsetIdx1s));

        sr.startIdx1 = subsetIdx1s.get(0);
        sr.stopIdx1 = subsetIdx1s.get(ns - 1);
        sr.mLen = sr.stopIdx1 - sr.startIdx1 + 1;
        sr.setChordSumNeedsUpdate(true);
        if (debug) {
            print(sr, "case 0 indiv", n2);
        }
        Integer k1 = Integer.valueOf(sr.startIdx1);
        t1.put(k1, sr);
        nMatched += sr.mLen;
    }

    private void addForCase1(SR sr, Entry<Integer, SR> below,
        int n2) {

        /*
               content ordered by idx1
                      t1  |  interval
               ------------------------
                          |
               startIdx1  | [sr]
                      -#- | [-#-]

        (1) find t1 ceiling for sr.stopIdx1 + 1.
        (2) test that entire range is consistent
            -- if ceiling startIdx2 is larger than sr.stopIdx2
               return is consistent
        (3) iterate over each idx1,idx2 in sr interval
            test for each ifx1,idx2
            -- if ceiling start idx2 is larger than idx2
               return is consistent
        */

        if (excludeCase1(below.getValue(), sr, n2)) {
            return;
        }
        
        if (case1AllConsistent(sr, below, n2)) {
            Integer k1 = Integer.valueOf(sr.startIdx1);
            t1.put(k1, sr);
            nMatched += sr.mLen;
            if (debug) {
                print(sr, "case 1", n2);
            }
            return;
        }

        SR ceil = below.getValue();
        calculateIds2s(ceil, n2);
        int ceilStrtIdx2 = cachedIdx2[0];

        TIntList subsetIdx1s = new TIntArrayList();

        populateCase1Idx1s(sr, ceilStrtIdx2, subsetIdx1s, n2);

        int ns = subsetIdx1s.size();
        if (ns < minLength) {
            return;
        }

        assert(assertContiguous(subsetIdx1s));

        sr.startIdx1 = subsetIdx1s.get(0);
        sr.stopIdx1 = subsetIdx1s.get(ns - 1);
        sr.mLen = sr.stopIdx1 - sr.startIdx1 + 1;
        sr.setChordSumNeedsUpdate(true);
        if (debug) {
            print(sr, "case 1 indev", n2);
        }
        Integer k1 = Integer.valueOf(sr.startIdx1);
        t1.put(k1, sr);
        nMatched += sr.mLen;
    }

    private boolean case0AllConsistent(SR sr,
        Entry<Integer, SR> strt1Floor, int n2) {

        SR floor = strt1Floor.getValue();
        calculateIds2s(floor, n2);
        int floorStopIdx2 = cachedIdx2[1];

        calculateIds2s(sr, n2);
        int startIdx2 = cachedIdx2[0];

        /*
               content ordered by idx1
                      t1  |  interval
               ------------------------
                          |
                      -#- | [-#-]
               startIdx1  | [sr]
        */
        if (floor.stopIdx1 < sr.startIdx1) {
            if (floorStopIdx2 < startIdx2) {
                return true;
            }
        }

        return false;
    }

    private void populateCase0Idx1s(int offset,
        TIntList inputIdx1s, int floorStopIdx1,
        int floorStopIdx2, TIntList outSubsetIdx1s, int n2) {

        for (int i = 0; i < inputIdx1s.size(); ++i) {
            int idx1 = inputIdx1s.get(i);
            int idx2 = idx1 + offset;
           
            if (idx1 > floorStopIdx1) {
                if (idx2 > floorStopIdx2) {
                    outSubsetIdx1s.add(idx1);
                }
            }
        }
    }

    private void populateCase0Idx1s(SR sr, int floorStopIdx1,
        int floorStopIdx2, TIntList outSubsetIdx1s, int n2) {

        TIntList input = new TIntArrayList();
        for (int idx1 = sr.startIdx1; idx1 <= sr.stopIdx1; ++idx1) {
            input.add(idx1);
        }

        populateCase0Idx1s(sr.offsetIdx2, input, floorStopIdx1,
            floorStopIdx2, outSubsetIdx1s, n2);
    }

    private void populateCase1Idx1s(int offset,
        TIntList inputIdx1s, int ceilStrtIdx2,
        TIntList subsetIdx1s, int n2) {

        for (int i = 0; i < inputIdx1s.size(); ++i) {
            int idx1 = inputIdx1s.get(i);
            int idx2 = idx1 + offset;
            
            /*
            content ordered by idx1
                      t1  |  interval
               ------------------------
                          |
               startIdx1  | [sr]
                      -#- | [-#-]
            */

            if (ceilStrtIdx2 > idx2) {
                subsetIdx1s.add(idx1);
            }
        }
    }

    private void populateCase1Idx1s(SR sr, int ceilStrtIdx2,
        TIntList outSubsetIdx1s, int n2) {

        TIntList input = new TIntArrayList();
        for (int idx1 = sr.startIdx1; idx1 <= sr.stopIdx1; ++idx1) {
            input.add(idx1);
        }

        populateCase1Idx1s(sr.offsetIdx2, input, ceilStrtIdx2,
            outSubsetIdx1s, n2);
    }

    private boolean assertContiguous(TIntList list) {

        if (list.size() < 2) {
            return true;
        }

        int prev = list.get(0);
        for (int i = 1; i < list.size(); ++i) {
            int v = list.get(i);
            if (v == (prev + 1)) {
                prev = v;
                continue;
            }
            return false;
        }

        return true;
    }

    private boolean case1AllConsistent(SR sr, Entry<Integer, SR> stp1Ceil,
        int n2) {

        SR ceil = stp1Ceil.getValue();
        calculateIds2s(ceil, n2);
        int ceilStrtIdx2 = cachedIdx2[0];

        calculateIds2s(sr, n2);
        int stpIdx2 = cachedIdx2[1];

        return (ceilStrtIdx2 > stpIdx2);

    }
    
    /**
     * calculate the idx2s from idx1 and offset.
     * can retrieve the result from cachedIdx2
     * @param sr 
     */
    private void calculateIds2s(SR sr, int n2) {
        cachedIdx2[0] = sr.startIdx1 + sr.offsetIdx2;
        cachedIdx2[1] = sr.stopIdx1 + sr.offsetIdx2;
    }

    private void addForCase2(SR sr, Entry<Integer, SR> strt1Floor,
        Entry<Integer, SR> stp1Ceil, int n2) {

        SR floor = strt1Floor.getValue();
        calculateIds2s(floor, n2);
        int floorStopIdx2 = cachedIdx2[1];
        
        TIntList subsetIdx1s = new TIntArrayList();

        populateCase0Idx1s(sr, floor.startIdx1, floorStopIdx2, subsetIdx1s, n2);

        int ns = subsetIdx1s.size();
        if (ns < minLength) {
            return;
        }

        assert(assertContiguous(subsetIdx1s));

        SR ceil = stp1Ceil.getValue();
        calculateIds2s(ceil, n2);
        int ceilStrtIdx2 = cachedIdx2[0];

        TIntList subsetIdx1s2 = new TIntArrayList();

        populateCase1Idx1s(sr.offsetIdx2, subsetIdx1s,
            ceilStrtIdx2, subsetIdx1s2, n2);

        ns = subsetIdx1s2.size();
        if (ns < minLength) {
            return;
        }

        assert(assertContiguous(subsetIdx1s2));

        sr.startIdx1 = subsetIdx1s2.get(0);
        sr.stopIdx1 = subsetIdx1s2.get(ns - 1);
        sr.mLen = sr.stopIdx1 - sr.startIdx1 + 1;
        sr.setChordSumNeedsUpdate(true);
        if (debug) {
            print(sr, "case 2 indev", n2);
        }
        Integer k1 = Integer.valueOf(sr.startIdx1);
        t1.put(k1, sr);
        nMatched += sr.mLen;
    }

    private boolean excludeCase0(SR above, SR sr, int n2) {
                
        // exclude if aboveStopIdx1 is near strt1
        //  and strt2 is far from abovestp2
        calculateIds2s(above, n2);
        int aboveStp2 = cachedIdx2[1];
        if ((sr.startIdx1 - above.startIdx1) < 2) {
            calculateIds2s(sr, n2);
            int strt2 = cachedIdx2[0];
            int stp2 = cachedIdx2[1];
            if ((strt2 < aboveStp2) || ((strt2 - aboveStp2) > 2)) {
                print(sr, "excluding by c0: ", n2);         
                return true;
            }
        }
        return false;
    }
    
    private boolean excludeCase1(SR below, SR sr, int n2) {
                
        // exclude if below.startIdx is nearly adj ro sr.stp1
        //  and belowstrtidx2 is far from stopIdx2
        calculateIds2s(below, n2);
        int belowStrt2 = cachedIdx2[0];
        if ((below.startIdx1 - sr.stopIdx1) < 2) {
            calculateIds2s(sr, n2);
            int strt2 = cachedIdx2[0];
            int stp2 = cachedIdx2[1];
            if ((belowStrt2 < stp2) || ((belowStrt2 - stp2) > 2)) {
                print(sr, "excluding by c1: ", n2);         
                return true;
            }
        }
        return false;
    }

    private boolean isEmbeddedOrInconsistentWith(SR sr, 
        Entry<Integer, SR> compE, int n2) {
        
        if (compE == null) {
            return false;
        }
        
        calculateIds2s(sr, n2);
        int strtIdx2 = cachedIdx2[0];
        int stpIdx2 = cachedIdx2[1];

        SR comp = compE.getValue();
        calculateIds2s(comp, n2);
        int compStrtIdx2 = cachedIdx2[0];
        int compStpIdx2 = cachedIdx2[1];

        if (comp.startIdx1 < sr.startIdx1) {
            if (!(compStrtIdx2 < stpIdx2)) {
                return true;
            }
        } else if (comp.startIdx1 > sr.startIdx1) {
            if (!(compStrtIdx2 > stpIdx2)) {
                return true;
            }
        }
        if (sr.startIdx1 >= comp.startIdx1 && sr.stopIdx1 <= comp.stopIdx1) {
            return true;
        }
        if (strtIdx2 >= compStrtIdx2 && stpIdx2 <= compStpIdx2) {
            return true;
        }

        return false;        
    }

}
