package algorithms.imageProcessing.matching;

import algorithms.imageProcessing.matching.PartialShapeMatcher2.SR;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.TreeMap;

/**
 * a class to handle additions and clockwise consistency checks
 * to a growing correspondency list.
 * The class is specialized for use with PartialShapeMatcher2.java
 * which has intervals it attempts to add to this structure in
 * order of increasing cost.
 * The intervals are given as ranges of indexes called idx1
 * and the shape they match to is specified as an offset
 * from the idx1 indexes.
 * Each interval has a single offset which may be different
 * from the offset in other intervals.
 *
 * @author nichole
 */
class OrderedClosedCurveCorrespondence {

    // visit each interval in order of smallest cost,
    // and only add the clockwise consistent intervals to a combined output
    protected TreeMap<Integer, PartialShapeMatcher2.SR> t1 = new TreeMap<Integer, PartialShapeMatcher2.SR>();

    private int minLength = 3;

    public void setMinimumLength(int length) {
        minLength = length;
    }

    public void addIntervals(List<SR> intervals, int n1, int n2) {

        for (PartialShapeMatcher2.SR sr: intervals) {
            boolean didIns = addInterval(sr, n1, n2);
        }

    }

    private void addFirstInterval(SR sr) {

        if (!t1.isEmpty()) {
            throw new IllegalStateException("addFirstInterval "
                + " is meant for use with an empty tree");
        }

        Integer k1 = Integer.valueOf(sr.startIdx1);
        t1.put(k1, sr);
    }

    /**
     * add intervals to the clockwise ordered unique correspondence
     * list intrnal to this instance.
     * Note that each interval is expected to be clockwise consistent
     * (stopIdx1 > startIdx1) and the list of intervals is expected
     * to be sorted so that the highest priority (== lowest cost)
     * intervals are at the smallest list indexes, that is the
     * list is increasing in cost with index.
     * @param sr
     * @param n1
     * @param n2
     * @return
     */
    public boolean addInterval(SR sr, int n1, int n2) {

        assert(sr.startIdx1 != sr.stopIdx1);

        if (t1.isEmpty()) {
            addFirstInterval(sr);
            return true;
        }

        //assert clockwise consistent

        int startIdx2 = sr.startIdx1 - sr.offsetIdx2;
        int stopIdx2 = sr.stopIdx1 - sr.offsetIdx2;
        if (startIdx2 < 0) {
            startIdx2 += n2;
            if (stopIdx2 < 0) {
                stopIdx2 += n2;
            }
        }

        //NOTE: if an interval is trimmed rather than discarded here
        // because of clockwise consistency,
        // then might need to consider re-doing the interval sort...
        // (an adaptive optimal: if an interval is trimmed, might
        // edit it in allResults list, re-sort and start again...
        // will not do that here, but might consider a better way to have
        // same result in the future).

        /*
        outlined here are decision trees for different use cases in order
        to find all needed datastructures.

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


            find t1 floor for sr.startIdx1 - 1.
            -- if floor stopIdx1 < sr.startIdx1
               -- if floor stopIdx2 < sr.startIdx2
                  can add interval
               -- else
                  can test each idx1, idx2 and add those
                  with idx2 > floor stopIdx2.
            -- else
               can test each idx1, idx2 in region and add those
                  with idx1 > floot stopIdx1 &&
                  with idx2 > floor stopIdx2

           then dividing the logic of the decision tree into 3 tests to make the
           logic usable within part of case 2 also.

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

        Then the methods for these invoked by the case finding
        logic are:
        
        - addForCase0(SR sr, Entry<Integer, SR> strt1Floor,
            int n2) : void
        - case0AllConsistent(SR sr,
            Entry<Integer, SR> strt1Floor, int n2) : boolean
        - populateCase0Idx1s(SR sr, int floorStopIdx1,
            int floorStopIdx2, TIntList subsetIdx1s, int n2) : void
        
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

        find t1 ceiling for sr.stopIdx1 + 1.
           -- if ceiling startIdx2 is larger than sr.stopIdx2
                 can add interval

        then dividing the logic of the decision tree into 3 tests to make
        the logic usable within part of case 2 also.

            (1) find t1 ceiling for sr.stopIdx1 + 1.

            (2) test that entire range is consistent
                -- if ceiling startIdx2 is larger than sr.stopIdx2
                   return is consistent
            (3) iterate over each idx1,idx2 in sr interval
                test for each ifx1,idx2
                -- if ceiling start idx2 is larger than idx2
                   return is consistent
        
        Then the methods for these invoked by the case finding
        logic are:
        
        - addForCase1(SR sr, Entry<Integer, SR> stp1Ceil, int n2) : void
        - case1AllConsistent(SR sr,
            Entry<Integer, SR> stp1Ceil, int n2) : boolean
        - populateCase1Idx1s(SR sr, int ceilStartIdx1,
            TIntList subsetIdx1s, int n2) : void
        
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

            for each idx1, idx2 in sr,
               test for case 1 and if returns true,
               test for case 0 and if returns true, add it

         Then the methods for these invoked by the case finding
         logic are:
        
         - addForCase2(SR sr, Entry<Integer, SR> strt1Floor, 
             Entry<Integer, SR> stp1Ceil, int n2) : void
         - populateCase0Idx1s(int offset, TIntList inputIdx1s, 
             int floorStopIdx1, int floorStopIdx2, 
             TIntList outSubsetIdx1s, int n2) : void
         - populateCase1Idx1s(int offset, TIntList inputIdx1s,
            int ceilStartIdx1, TIntList subsetIdx1s, int n2) : void
        
        --------------------------------------------------------

        */

        Integer strt1 = Integer.valueOf(sr.startIdx1);
        Integer stp1 = Integer.valueOf(sr.stopIdx1);
        int strt2Int = sr.startIdx1 - sr.offsetIdx2;
        int stp2Int = sr.stopIdx1 - sr.offsetIdx2;
        if (strt2Int < 0) {
            strt2Int += n2;
        }
        if (stp2Int < 0) {
            stp2Int += n2;
        }
        Integer strt2 = Integer.valueOf(strt2Int);
        Integer stp2 = Integer.valueOf(stp2Int);

        Entry<Integer, SR> strt1Ceil = t1.ceilingEntry(strt1);

        if (strt1Ceil == null) {

            // case 0

            assert(sr.startIdx1 > 0);

            Entry<Integer, SR> strt1Floor = t1.floorEntry(
                Integer.valueOf(sr.startIdx1 - 1));

            if (strt1Floor == null) {
                throw new IllegalStateException("Error in algorithm: "
                    + " there should be at least one entry in t1 at this point");
            }

            addForCase0(sr, strt1Floor, n2);

        } else {

            if (sr.startIdx1 == 0) {
                
                // case 1

                Entry<Integer, SR> stp1Ceil = t1.floorEntry(
                    Integer.valueOf(sr.stopIdx1 + 1));
                
                addForCase1(sr, stp1Ceil, n2);
                
            } else {

                Entry<Integer, SR> strt1Floor = t1.floorEntry(
                    Integer.valueOf(sr.startIdx1 - 1));
                
                Entry<Integer, SR> stp1Ceil = t1.floorEntry(
                    Integer.valueOf(sr.stopIdx1 + 1));
                
                if (strt1Floor == null) {
                    
                    // case 1

                    addForCase1(sr, stp1Ceil, n2);
                    
                } else {
                    
                    // case 2

                    addForCase2(sr, strt1Floor, stp1Ceil, n2);
                }
            }
        }

        throw new UnsupportedOperationException("not yet mplemented");
    }

    private void addForCase0(SR sr, Entry<Integer, SR> strt1Floor,
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

        if (case0AllConsistent(sr, strt1Floor, n2)) {
            Integer k1 = Integer.valueOf(sr.startIdx1);
            t1.put(k1, sr);
            return;
        }

        SR floor = strt1Floor.getValue();
        int floorStopIdx1 = floor.stopIdx1;
        int floorStopIdx2 = floor.stopIdx1 - floor.offsetIdx2;
        if (floorStopIdx2 < 0) {
            floorStopIdx2 += n2;
        }
        if (floorStopIdx2 < floorStopIdx1) {
            floorStopIdx2 += n2;
        }

        TIntList subsetIdx1s = new TIntArrayList();

        populateCase0Idx1s(sr, floorStopIdx1, floorStopIdx2, subsetIdx1s, n2);

        int ns = subsetIdx1s.size();
        if (ns < minLength) {
            return;
        }

        assert(assertContiguous(subsetIdx1s));

        sr.startIdx1 = subsetIdx1s.get(0);
        sr.stopIdx1 = subsetIdx1s.get(ns - 1);
        sr.setChordSumNeedsUpdate(true);

        Integer k1 = Integer.valueOf(sr.startIdx1);
        t1.put(k1, sr);
    }

    private void addForCase1(SR sr, Entry<Integer, SR> stp1Ceil,
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

        if (case1AllConsistent(sr, stp1Ceil, n2)) {
            Integer k1 = Integer.valueOf(sr.startIdx1);
            t1.put(k1, sr);
            return;
        }

        SR ceil = stp1Ceil.getValue();
        int ceilStrtIdx2 = ceil.startIdx1 - ceil.offsetIdx2;
        if (ceilStrtIdx2 < 0) {
            ceilStrtIdx2 += n2;
        }
        if (ceilStrtIdx2 < ceil.startIdx1) {
            ceilStrtIdx2 += n2;
        }

        TIntList subsetIdx1s = new TIntArrayList();

        populateCase1Idx1s(sr, ceilStrtIdx2, subsetIdx1s, n2);

        int ns = subsetIdx1s.size();
        if (ns < minLength) {
            return;
        }

        assert(assertContiguous(subsetIdx1s));

        sr.startIdx1 = subsetIdx1s.get(0);
        sr.stopIdx1 = subsetIdx1s.get(ns - 1);
        sr.setChordSumNeedsUpdate(true);

        Integer k1 = Integer.valueOf(sr.startIdx1);
        t1.put(k1, sr);
    }
    
    private boolean case0AllConsistent(SR sr,
        Entry<Integer, SR> strt1Floor, int n2) {

        SR floor = strt1Floor.getValue();
        int floorStopIdx1 = floor.stopIdx1;
        int floorStopIdx2 = floor.stopIdx1 - floor.offsetIdx2;
        if (floorStopIdx2 < 0) {
            floorStopIdx2 += n2;
        }
        if (floorStopIdx2 < floorStopIdx1) {
            floorStopIdx2 += n2;
        }

        int startIdx1 = sr.startIdx1;
        int startIdx2 = sr.startIdx1 - sr.offsetIdx2;
        if (startIdx2 < 0) {
            startIdx2 += n2;
        }
        if (startIdx2 < startIdx1) {
            startIdx2 += n2;
        }

        /*
               content ordered by idx1
                      t1  |  interval
               ------------------------
                          |
                      -#- | [-#-]
               startIdx1  | [sr]
        */
        if (floorStopIdx1 < startIdx1) {
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
            int idx2 = idx1 - offset;
            if (idx2 < 0) {
                idx2 += n2;
            }
            if (idx2 < idx1) {
                idx2 += n2;
            }

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
            int idx2 = idx1 - offset;
            if (idx2 < 0) {
                idx2 += n2;
            }
            if (idx2 < idx1) {
                idx2 += n2;
            }
            
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
            if (v == (prev - 1)) {
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
        int ceilStrtIdx2 = ceil.startIdx1 - ceil.offsetIdx2;
        if (ceilStrtIdx2 < 0) {
            ceilStrtIdx2 += n2;
        }
        if (ceilStrtIdx2 < ceil.startIdx1) {
            ceilStrtIdx2 += n2;
        }

        int stpIdx2 = sr.stopIdx1 - sr.offsetIdx2;
        if (stpIdx2 < 0) {
            stpIdx2 += n2;
        }
        if (stpIdx2 < sr.startIdx1) {
            stpIdx2 += n2;
        }
        
        return (ceilStrtIdx2 > stpIdx2);
        
    }

    private void addForCase2(SR sr, Entry<Integer, SR> strt1Floor, 
        Entry<Integer, SR> stp1Ceil, int n2) {
        
        SR floor = strt1Floor.getValue();
        int floorStopIdx1 = floor.stopIdx1;
        int floorStopIdx2 = floor.stopIdx1 - floor.offsetIdx2;
        if (floorStopIdx2 < 0) {
            floorStopIdx2 += n2;
        }
        if (floorStopIdx2 < floorStopIdx1) {
            floorStopIdx2 += n2;
        }

        TIntList subsetIdx1s = new TIntArrayList();

        populateCase0Idx1s(sr, floorStopIdx1, floorStopIdx2, subsetIdx1s, n2);

        int ns = subsetIdx1s.size();
        if (ns < minLength) {
            return;
        }

        SR ceil = stp1Ceil.getValue();
        int ceilStrtIdx2 = ceil.startIdx1 - ceil.offsetIdx2;
        if (ceilStrtIdx2 < 0) {
            ceilStrtIdx2 += n2;
        }
        if (ceilStrtIdx2 < ceil.startIdx1) {
            ceilStrtIdx2 += n2;
        }

        TIntList subsetIdx1s2 = new TIntArrayList();

        populateCase1Idx1s(sr.offsetIdx2, subsetIdx1s,
            ceilStrtIdx2, subsetIdx1s2, n2);

        ns = subsetIdx1s2.size();
        if (ns < minLength) {
            return;
        }

        assert(assertContiguous(subsetIdx1s));

        sr.startIdx1 = subsetIdx1s.get(0);
        sr.stopIdx1 = subsetIdx1s.get(ns - 1);
        sr.setChordSumNeedsUpdate(true);

        Integer k1 = Integer.valueOf(sr.startIdx1);
        t1.put(k1, sr);       
    }

}
