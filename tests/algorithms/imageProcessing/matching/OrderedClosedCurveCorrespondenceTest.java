package algorithms.imageProcessing.matching;

import algorithms.imageProcessing.matching.PartialShapeMatcher.SR;
import algorithms.util.PairInt;
import java.util.ArrayList;
import java.util.List;
import junit.framework.TestCase;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class OrderedClosedCurveCorrespondenceTest extends TestCase {
    
    public OrderedClosedCurveCorrespondenceTest() {
    }
    
    public void testAddIntervals() {
        
        List<SR> intervals = new ArrayList<SR>();
        
        List<PairInt> ranges = new ArrayList<PairInt>();
        ranges.add(new PairInt(0, 3));
        ranges.add(new PairInt(4, 6));
        ranges.add(new PairInt(7, 9));
        
        for (int i = 0; i < ranges.size(); ++i) {
            PairInt range = ranges.get(i);
            SR sr = new SR();
            sr.startIdx1 = range.getX(); sr.stopIdx1 = range.getY(); 
            sr.row = 0; 
            sr.mLen = sr.stopIdx1 - sr.startIdx1 + 1; 
            sr.nMax = 10; 
            sr.diffChordSum = 1; sr.diffChordSum = 10;
            intervals.add(sr);
        }
        
        OrderedClosedCurveCorrespondence occ = new 
            OrderedClosedCurveCorrespondence();
        occ.addIntervals(intervals, 10, 10);
        
        List<SR> results = occ.getResultsAsList();
        
        assertEquals(3, results.size());
        
        for (SR result : results) {
            PairInt found = new PairInt(result.startIdx1, result.stopIdx1);
            assertTrue(ranges.remove(found));
        }
        assertEquals(0, ranges.size());
        
        // ------ same ranges and different order ----
        ranges = new ArrayList<PairInt>();
        ranges.add(new PairInt(4, 6));
        ranges.add(new PairInt(7, 9));
        ranges.add(new PairInt(0, 3));
        
        for (int i = 0; i < ranges.size(); ++i) {
            PairInt range = ranges.get(i);
            SR sr = new SR();
            sr.startIdx1 = range.getX(); sr.stopIdx1 = range.getY(); 
            sr.row = 0; 
            sr.mLen = sr.stopIdx1 - sr.startIdx1 + 1;
            sr.nMax = 10; 
            sr.diffChordSum = 1; sr.diffChordSum = 10;
            intervals.add(sr);
        }
        
        occ = new OrderedClosedCurveCorrespondence();
        occ.addIntervals(intervals, 10, 10);
        
        results = occ.getResultsAsList();
        
        assertEquals(3, results.size());
        
        for (SR result : results) {
            PairInt found = new PairInt(result.startIdx1, result.stopIdx1);
            assertTrue(ranges.remove(found));
        }
        assertEquals(0, ranges.size());
    }
    
    public void testAddIntervals2() {
        
        List<SR> intervals = new ArrayList<SR>();
        
        /*
        0 : 3   6 : 9
        4 : 6   0 : 2
        7 : 9   3 : 5
        */
        List<PairInt> ranges = new ArrayList<PairInt>();
        ranges.add(new PairInt(4, 6));
        ranges.add(new PairInt(7, 9));
        ranges.add(new PairInt(0, 3));
        
        for (int i = 0; i < ranges.size(); ++i) {
            PairInt range = ranges.get(i);
            SR sr = new SR();
            sr.startIdx1 = range.getX(); sr.stopIdx1 = range.getY(); 
            sr.offsetIdx2 = 6;
            sr.mLen = sr.stopIdx1 - sr.startIdx1 + 1;
            sr.row = 0; sr.nMax = 10; 
            sr.diffChordSum = 1; sr.diffChordSum = 10;
            intervals.add(sr);
        }
        
        OrderedClosedCurveCorrespondence occ = new 
            OrderedClosedCurveCorrespondence();
        occ.addIntervals(intervals, 10, 10);
        
        List<SR> results = occ.getResultsAsList();
        
        assertEquals(3, results.size());
        
        
    }
    
    public void testAddIntervals3() {
        
        List<SR> intervals = new ArrayList<SR>();
        
         /*
              0  1  2  3  4  5  6  7  8  9
            9          -                    
            8       -                    
            7    -                    
            6 -         
            5                                  )
            4                               +
            3                            +
            2                   0  
            1                0   
            0             0
              0  1  2  3  4  5  6  7  8  9 10 11
         
        */
        
        /*
        0 : 3    6 : 9  offset=6
        4 : 6    0 : 2  offset=6
        9 : 11   3 : 5  offset= 3 + (n2-1-9) + 1 = 4     
        */
        int n1 = 12; int n2 = 10;
        List<PairInt> ranges = new ArrayList<PairInt>();
        ranges.add(new PairInt(4, 6));
        ranges.add(new PairInt(9, 11));
        ranges.add(new PairInt(0, 3));
        
        for (int i = 0; i < ranges.size(); ++i) {
            PairInt range = ranges.get(i);
            SR sr = new SR();
            sr.startIdx1 = range.getX(); sr.stopIdx1 = range.getY(); 
            sr.offsetIdx2 = 6;
            sr.row = 0; 
            sr.mLen = sr.stopIdx1 - sr.startIdx1 + 1;
            sr.nMax = n1; 
            sr.diffChordSum = 1; sr.diffChordSum = 10;
            if (range.getX() == 4 && range.getY() == 6) {
                sr.offsetIdx2 = 6;
            } else if (range.getX() == 9 && range.getY() == 11) {
                sr.offsetIdx2 = 4;
            }
            intervals.add(sr);
        }
        
        OrderedClosedCurveCorrespondence occ = new 
            OrderedClosedCurveCorrespondence();
        occ.addIntervals(intervals, n1, n2);
        
        List<SR> results = occ.getResultsAsList();
        
        assertEquals(3, results.size());
        
    }
    
    public void testAddIntervals4() {
        
        List<SR> intervals = new ArrayList<SR>();
        
         /*
              0  1  2  3  4  5  6  7  8  9
            9          -                    
            8       -                    
            7    -                    
            6 -         
            5                                  )
            4                               +
            3                            +
            2                   0  
            1                0   
            0             0
              0  1  2  3  4  5  6  7  8  9 10 11
         
        */
        
        /*
        0 : 3    6 : 9  offset=6
        4 : 6    0 : 2  offset=6
        9 : 11   3 : 5  offset= 3 + (n2-1-9) + 1 = 4     
        */
        int n1 = 12; int n2 = 10;
        List<PairInt> ranges = new ArrayList<PairInt>();
        ranges.add(new PairInt(4, 6));
        ranges.add(new PairInt(9, 11));
        ranges.add(new PairInt(0, 3));
        
        for (int i = 0; i < ranges.size(); ++i) {
            PairInt range = ranges.get(i);
            SR sr = new SR();
            sr.startIdx1 = range.getX(); sr.stopIdx1 = range.getY(); 
            sr.offsetIdx2 = 6;
            sr.row = 0; 
            sr.mLen = sr.stopIdx1 - sr.startIdx1 + 1;
            sr.nMax = n1; 
            sr.diffChordSum = 1; sr.diffChordSum = 10;
            if (range.getX() == 4 && range.getY() == 6) {
                sr.offsetIdx2 = 6;
            } else if (range.getX() == 9 && range.getY() == 11) {
                sr.offsetIdx2 = 4;
            }
            intervals.add(sr);
        }
        
        // add a bad interval last that overlaps the others and should
        //  not be successfullt added to the tree
        SR sr = new SR();
        sr.startIdx1 = 6; sr.stopIdx1 = 9; 
        sr.offsetIdx2 = 0;
        sr.row = 0; 
        sr.mLen = sr.stopIdx1 - sr.startIdx1 + 1;
        sr.nMax = n1; 
        sr.diffChordSum = 1; sr.diffChordSum = 10;
        intervals.add(sr);
        
        OrderedClosedCurveCorrespondence occ = new 
            OrderedClosedCurveCorrespondence();
        occ.addIntervals(intervals, n1, n2);
        
        List<SR> results = occ.getResultsAsList();
        
        assertEquals(3, results.size());
        
    }
}
