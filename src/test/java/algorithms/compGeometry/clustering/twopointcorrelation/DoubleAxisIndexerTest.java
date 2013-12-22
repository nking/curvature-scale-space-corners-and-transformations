package algorithms.compGeometry.clustering.twopointcorrelation;

import java.security.SecureRandom;
import java.util.logging.Logger;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class DoubleAxisIndexerTest extends TestCase {

    protected Logger log = Logger.getLogger(this.getClass().getSimpleName());

    public DoubleAxisIndexerTest(String testName) {
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



     public void testFindXYMinMax() throws Exception {
        log.info("testFindXYMinMax()");

        float[] x = new float[] {
            0,  1,  2,  3,  4,  5,  6,  7,  8,
            0,              4,  5,  6,  7,  8,
            0,              4,  5,  6,  7,  8,
            0,              4,  5,  6,  7,  8,
            0,  1,  2,  3,  4,  5,  6,  7,  8,
            0,  1,  2,  3,  4,              8,
            0,  1,  2,  3,  4,              8,
            0,  1,  2,  3,  4,              8,
            0,  1,          4,  5,  6,  7,  8,
            0,  1,          4,  5,  6,  7,  8,
            0,  1,  2,  3,  4,  5,  6,  7,  8
        };
        float[] y = new float[] {
            0,  0,  0,  0,  0,  0,  0,  0,  0,
            1,              1,  1,  1,  1,  1,
            2,              2,  2,  2,  2,  2,
            3,              3,  3,  3,  3,  3,
            4,  4,  4,  4,  4,  4,  4,  4,  4,
            5,  5,  5,  5,  5,              5,
            6,  6,  6,  6,  6,              6,
            7,  7,  7,  7,  7,              7,
            8,  8,          8,  8,  8,  8,  8,
            9,  9,          9,  9,  9,  9,  9,
           10, 10, 10, 10, 10, 10, 10, 10, 10
        };


        DoubleAxisIndexer indexer = new DoubleAxisIndexer();
        indexer.sortAndIndexX(x, y, x.length);

        float[] minMaxes = indexer.findXYMinMax();

        assertTrue(minMaxes[0] == 0.0f);
        assertTrue(minMaxes[1] == 8.0f);
        assertTrue(minMaxes[2] == 0.0f);
        assertTrue(minMaxes[3] == 10.0f);
    }

    public void testSortAndIndexXThenY() throws Exception {
        log.info("testSortAndIndexXThenY()");

        int npoints = 100;

        float[] x = new float[npoints];
        float[] y = new float[npoints];

        /* make a test data set that is already ordered as expected:
         *           [3]
         *        [2]
         *     [1]
         *  [0]
         */

        for (int i = 0; i < npoints; i++) {
            x[i] = i;
            y[i] = i;
        }

        DoubleAxisIndexer indexer = new DoubleAxisIndexer();
        indexer.sortAndIndexX(x, y, npoints);

        int[] xIndexes = indexer.getSortedXIndexes();

        float[] xOriginal = indexer.getX();
        float[] yOriginal = indexer.getY();

        for (int i = 0; i < npoints; i++) {
            assertTrue(xIndexes[i] == i);

            int xIndex = xIndexes[i];
            assertTrue(xOriginal[xIndex] == i);
            assertTrue(yOriginal[xIndex] == i);
        }
    }

    public void testSortAndIndexXThenY_0() throws Exception {

        log.info("testSortAndIndexXThenY_0()");

        int npoints = 100;

        float[] x = new float[npoints];
        float[] y = new float[npoints];

        /* make a test data set that is ordered exactly oppossite for both x and y:
         *           [0]
         *        [1]
         *     [2]
         *  [3]
         */

        int count = npoints-1;
        for (int i = 0; i < npoints; i++) {
            x[i] = count;
            y[i] = count;
            count--;
        }

        DoubleAxisIndexer indexer = new DoubleAxisIndexer();
        indexer.sortAndIndexX(x, y, npoints);

        int[] xIndexes = indexer.getSortedXIndexes();

        float[] xOriginal = indexer.getX();
        float[] yOriginal = indexer.getY();

        for (int i = 0; i < npoints; i++) {
            int xIndex = xIndexes[i];
            assertTrue(xOriginal[xIndex] == i);
            assertTrue(yOriginal[xIndex] == i);
        }
    }

    public void testSortAndIndexXThenY_1() throws Exception {

        log.info("testSortAndIndexXThenY_1()");

        int npoints = 100;

        float[] x = new float[npoints];
        float[] y = new float[npoints];

        /* make a test data set that is ordered exactly opposite for x but not y:
         *  [3]
         *     [2]
         *        [1]
         *           [0]
         */

        int count = npoints-1;
        for (int i = 0; i < npoints; i++) {
            x[i] = count;
            y[i] = i;
            count--;
        }

        DoubleAxisIndexer indexer = new DoubleAxisIndexer();
        indexer.sortAndIndexX(x, y, npoints);

        int[] xIndexes = indexer.getSortedXIndexes();

        float[] xOriginal = indexer.getX();
        float[] yOriginal = indexer.getY();

        count = npoints-1;
        for (int i = 0; i < npoints; i++) {

            int xIndex = xIndexes[i];
            
            assertTrue(xIndex == count);

            assertTrue(xOriginal[xIndex] == i);
            assertTrue(yOriginal[xIndex] == count);
            
            count--;
        }
    }

    public void testSortAndIndexXThenY_2() throws Exception {
        log.info("testSortAndIndexXThenY_2()");

        int npoints = 100;

        float[] x = new float[npoints];
        float[] y = new float[npoints];

        /*
         * randomly create data then assert that value  previous is less than last for indexes
         */
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        //sr.setSeed(System.currentTimeMillis());
        sr.setSeed(123456789);

        for (int i = 0; i < npoints; i++) {
            x[i] = sr.nextFloat()*(float)(npoints - 1);
            y[i] = sr.nextFloat()*(float)(npoints - 1);
        }

        DoubleAxisIndexer indexer = new DoubleAxisIndexer();
        indexer.sortAndIndexX(x, y, npoints);

        // compare results of sorted indexes
        int[] sortedIndexes = null;
        float[] sortedValues = null;
        int caseCount = 0;
        for (int i = 0; i < 1; i++) {
            switch(caseCount) {
                case 0:
                    sortedIndexes = indexer.getSortedXIndexes();
                    sortedValues = indexer.getX();
                    break;
                default:
                    break;
            }
            float previous = sortedValues[sortedIndexes[0]];
            for (int j = 1; j < sortedValues.length; j++) {
                int index = sortedIndexes[j];
                float current = sortedValues[index];
                //String str = String.format("(%d %d) %.2f %.2f", j, index, previous, current );
                //log.info(str);
                assertTrue(previous <= current);
                previous = current;
            }
        }
        int z = 1;
    }
}
