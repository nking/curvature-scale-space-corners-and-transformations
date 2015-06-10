package algorithms.imageProcessing;

import algorithms.util.PairInt;
import java.util.HashSet;
import java.util.Set;
import junit.framework.TestCase;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class PostLineThinnerCorrectionsTest extends TestCase {
    
    public PostLineThinnerCorrectionsTest() {
    }
    
    public void testCorrectForArtifacts() {
        
    }
    
    public void testCorrectForHoleArtifacts1() throws Exception {
        
        /* 
                      1                   7
                      1*              2   6
                 1    0    1          1   5 
                      1               0   4
                      1              -1   3
        
           -2   -1    0    1    2
            2   3     4    5    6
        */ 
        
        int w = 10;
        int h = 10;
        
        Set<PairInt> points = new HashSet<PairInt>();
        
        for (int col = 3; col <= 5; col++) {
            if (col != 4) {
                points.add(new PairInt(col, 5));
            }
        }
        for (int row = 3; row <= 7; row++) {
            if (row != 5) {
                points.add(new PairInt(4, row));
            }
        }
        
        assertTrue(points.size() == 6);
                
        PostLineThinnerCorrections pc = new PostLineThinnerCorrections();
        
        pc.correctForHoleArtifacts00_10(points, w, h);
        
        /* 
                      1                   7
                      1*              2   6
                 1    0    1          1   5 
                      1               0   4
                      1              -1   3
        
           -2   -1    0    1    2
            2   3     4    5    6
        
        
                      1                   7
                      1               2   6
                      1               1   5 
                      1*              0   4
                      1              -1   3
        
           -2   -1    0    1    2
            2    3    4    5    6
        */
        for (int row = 3; row <= 7; row++) {
            assertTrue(points.contains(new PairInt(4, row)));
        }
        
        assertFalse(points.contains(new PairInt(3, 5)));
        assertFalse(points.contains(new PairInt(5, 5)));
    }

}
