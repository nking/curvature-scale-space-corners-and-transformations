package algorithms.compGeometry;

import java.util.logging.Logger;
import junit.framework.TestCase;

public class LinesAndAnglesTest extends TestCase {

    protected Logger log = Logger.getLogger(this.getClass().getSimpleName());

    public LinesAndAnglesTest(String testName) {
        super(testName);
    }

    @Override
    protected void setUp() throws Exception {
        super.setUp();
    }

    @Override
    protected void tearDown() throws Exception {
        super.tearDown();
    }

    public void testCrossProduct() {

        log.info("crossProduct");

        /*         |
         *       4 -
         *         |
         *       3 -  o2
         *         |
         * 1 o   2 -
         *         |
         *       1 -
         *         |
         *       0 o--|--|--|--|--|--|
         *         0  1  2  3  4
         *
         */
        double x1 = -3;
        double y1 = 2;
        double x2 = 1;
        double y2 = 3;
        double result = LinesAndAngles.crossProduct(x1, y1, x2, y2);
        assertTrue(result < 0); // clockwise
        assertTrue(result == -11);

        x1 = 3;
        result = LinesAndAngles.crossProduct(x1, y1, x2, y2);
        assertTrue(result >= 0); // counter clockwise
        assertTrue(result == 7);
    }

    
    public void testCalculatePolarSineTheta() throws Exception {
        
    }
}
