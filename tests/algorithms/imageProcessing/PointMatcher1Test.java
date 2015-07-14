package algorithms.imageProcessing;

import java.util.logging.Logger;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class PointMatcher1Test extends TestCase {

    private Logger log = Logger.getLogger(this.getClass().getName());

    public PointMatcher1Test() {
    }

    public void testSortByDescendingMatches3() throws Exception {
        
        PointMatcher matcher = new PointMatcher();
        
        TransformationPointFit[] fits = new TransformationPointFit[4];
        
        fits[0] = new TransformationPointFit(
            new TransformationParameters(),
            10, 10.0, 5.0, 
            10.0f, 10.0f);
        
        fits[1] = new TransformationPointFit(
            new TransformationParameters(),
            20, 5.0, 2.0, 
            10.0f, 10.0f);
        
        fits[2] = new TransformationPointFit(
            new TransformationParameters(),
            2, 25.0, 20.0, 
            10.0f, 10.0f);
        
        fits[3] = new TransformationPointFit(
            new TransformationParameters(),
            20, 7.0, 6.0, 
            10.0f, 10.0f);
      
        matcher.sortByDescendingMatches(fits, 0, fits.length - 1);
        
        assertTrue(fits[0].getNumberOfMatchedPoints() == 20);
        assertTrue(Math.abs(fits[0].getMeanDistFromModel() - 5) < 0.1);
        
        assertTrue(fits[1].getNumberOfMatchedPoints() == 20);
        assertTrue(Math.abs(fits[1].getMeanDistFromModel() - 7) < 0.1);
        
        assertTrue(fits[2].getNumberOfMatchedPoints() == 10);
        
        assertTrue(fits[3].getNumberOfMatchedPoints() == 2);
    }
    
    public static void main(String[] args) {

        try {
            PointMatcher1Test test = new PointMatcher1Test();

            test.testSortByDescendingMatches3();

        } catch(Exception e) {
            e.printStackTrace();
            System.out.println(e.getMessage());
            fail(e.getMessage());
        }
    }
}
