package algorithms.imageProcessing.matching;

import algorithms.imageProcessing.matching.PartialShapeMatcher2.SR;
import java.util.List;
import java.util.Map;
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

    public void addIntervals(List<SR> intervals, int n1, int n2) {
        
        for (PartialShapeMatcher2.SR sr: intervals) {
            boolean didIns = addInterval(sr, n1, n2);
        }
        
    }
    
    private void addFirstInterval(SR sr, int n2) {
        
        if (!t1.isEmpty()) {
            throw new IllegalStateException("addFirstInterval "
                + " is meant for use with an empty tree");
        }
        
        for (int i1 = sr.startIdx1; i1 <= sr.stopIdx1; ++i1) {
            int i2 = i1 - sr.offsetIdx2;
            if (i2 < 0) {
                i2 += n2;
            }
            Integer k1 = Integer.valueOf(i1);
            t1.put(k1, sr);
        }
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
            addFirstInterval(sr, n2);
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

        NOTE: ideally would like simpler logic such as
            t1 w/ keys of idx1 and t2 w/ keys of idx2
            being used in simplest manner to assert consistency

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

        
            find t1 floor for sr.startIdx1.
            -- if floor stopIdx1 >= sr.startIdx1
               -- if floor stopIdx1 > sr.stopIdx1
                  sr is not consistent, so do not add
               -- else 
                  a subset of sr might be added where 
                  idx1 is > floor stopIdx1 if passes idx2 checks:
                  -- if floor stopIdx2 > sr.stopIdx2
                     sr is not consistent, so do not add
                  -- else
                     subset of sr where 
                     idx1 is > floor stopIdx1 and
                     idx2 is > floor stopIdx2
                     can be added
            -- else floor stopIdx1 is < sr.startIdx1
               -- if floor stopIdx2 >= sr.startIdx2
                  -- if floor stopIdx2 > sr.stopIdx2
                     sr is not consistent, so do not add
                  -- else
                     subset of sr where 
                     idx2 is > floor stopIdx2
                     can be added
               -- else floor stopIdx2 is < sr.startIdx2
                  can add interval
        
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
       
        find t1 ceiling for sr.stopIdx1.
           -- if ceiling startIdx1 is larger than sr.stopIdx1
              -- if ceiling startIdx2 is larger than sr.stopIdx2
                 can add interval
              -- else 
                 can add subset of interval sr where 
                     idx2 < sr.startIdx2
           -- else ceiling startIdx1 is smaller than or equal to 
              sr.stopIdx1
                  can add subset of interval sr where 
                     idx1 < sr.startIdx1 and idx2 < sr.startIdx2
           
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
        
        --------------------------------------------------------

        */

        throw new UnsupportedOperationException("not yet mplemented");
    }        
}
