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

    // note: these variables are evolving while outlining the decision tree:
    /**
     * wIdx2 is an idx1 index. if an idx2==0 is in t1 and if that occurs after
     * an idx2 > 0, mapping of idx1 to idx2, wIdx2 will be set and the wrap
     * around is known precisely as this value which is idx1 (and it maps to
     * idx2 == 0).
     */
    protected int wIdx2 = -1;

    /**
     * maxIdx2BeforeW is an idx1 index. if wIdx2UpperLimit is > -1 that means a
     * wrap around has been found but the true start of it may occur at an
     * earlier index idx1. maxIdx2BeforeW is the largest idx1 index which maps
     * to an idx2 before the wrap around. maxIdx2BeforeW reaches its maximum
     * when idx1 maps to n2-1, but that mapping might not exist for this p,q
     * correspondence. This variable is used in consistency checks if wIdx2==-1
     * and maxIdx2BeforeW is greater than -1.
     */
    protected int maxIdx2BeforeW = -1;

    /**
     * wIdx2UpperLimit is an idx1 index. if an idx2 is in t1 and if that occurs
     * after an idx2 of smaller value wIdx2UpperLimit will be set and updated
     * for matchings of smaller idx1 and idx2 that approach idx2 being 0. this
     * indicates that a wrap around has occurred, but the final first location
     * of the idx2 wrap around is not necessarily determined. This means that
     * between the interval above wIdx2UpperLimit and the interval containing
     * wIdx2UpperLimit, an inserted interval might have idx1 indexes that map to
     * a region between idx1 after maxIdx2BeforeW and before n2-1 or it might
     * map to region between where idx1 maps to idx2=0 and idx1=wIdx2UpperLimit.
     */
    protected int wIdx2UpperLimit = -1;

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
            if (wIdx2 > -1) {
                // the exact location of idx2 wrap around is found
                continue;
            }
            if (i2 == (n2 - 1)) {
                maxIdx2BeforeW = i1;
                continue;
            }
            if (i2 < i1) {
                // wrap around has or is occurring
                if (i2 == 0) {
                    wIdx2 = i1;
                } else if (wIdx2UpperLimit == -1) {
                    wIdx2UpperLimit = i1;
                }
            } else {
                maxIdx2BeforeW = i1;
            }
        }
    }
    
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
                and value = interval sr

            the wrap around variables above are used to check consistency of
                idx2 indexes.

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
               
        paused here...
        
        --------------------------------------------------------
        case 1: sr.startIdx1 ceiling is not null, that is, there are
                intervals in t1 at same or larger index position
                than st.startIdx1 
                and there may or may not be intervals in t1 at a
                smaller index position that st.startIdx1.

            content ordered by idx1
                      t1  |  interval
               ------------------------
                          |      
               startIdx1  | [sr]
                      -#- | [-#-]
       
        paused here
           -- ceiling key is larger than sr.stopIdx1
           -- ceiling key is not larger than sr.stopIdx1
        
        --------------------------------------------------------

        */

        throw new UnsupportedOperationException("not yet mplemented");
    }        
}
