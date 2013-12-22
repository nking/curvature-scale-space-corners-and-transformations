package algorithms.compGeometry.clustering.twopointcorrelation;

import java.io.File;
import static junit.framework.Assert.assertNotNull;
import static junit.framework.Assert.assertTrue;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class SerializerUtilTest extends TestCase {

    public SerializerUtilTest(String testName) {
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

    /**
     * Test of findResourcesDirectory method, of class SerializerUtil.
     */
    public void testFindTmpDataDirectory() throws Exception {

        String result = SerializerUtil.findTmpDataDirectory();

        assertNotNull(result);
        assertTrue(result.length() > 0);

        File dir = new File(result);
        assertTrue(dir.isDirectory());
    }

    /**
     * Test of serializeIndexer method, of class SerializerUtil.
     */
    public void testSerializeIndexer() throws Exception {
        int npoints = 100;

        float[] x = new float[npoints];
        float[] y = new float[npoints];
        for (int i = 0; i < npoints; i++) {
            x[i] = i;
            y[i] = i;
        }

        DoubleAxisIndexer indexer = new DoubleAxisIndexer();
        indexer.sortAndIndexX(x, y, npoints);

        String result = SerializerUtil.serializeIndexer(indexer);

        assertNotNull(result);
        assertTrue(result.length() > 0);

        File fl = new File(result);
        assertTrue(fl.isFile());
        assertTrue(fl.exists());

        DoubleAxisIndexer indexer2 = SerializerUtil.readPersistedPoints(result);

        assertTrue(indexer2.getNumberOfPoints() == indexer.getNumberOfPoints());

        for (int ii = 0; ii < 3; ii++) {
            float[] array1;
            float[] array2;
            switch(ii) {
                case 0:
                    array1 = indexer.getX();
                    array2 = indexer2.getX();
                    break;
                case 1:
                    array1 = indexer.getY();
                    array2 = indexer2.getY();
                    break;
                case 2:
                    array1 = indexer.getXErrors();
                    array2 = indexer2.getXErrors();
                    break;
                case 3:
                    array1 = indexer.getYErrors();
                    array2 = indexer2.getYErrors();
                    break;
                default:
                    array1 = new float[0];
                    array2 = new float[0];
            }
            if (array1 == null) {
                assertNull(array2);
            } else if (array2 == null) {
                assertNull(array1);
            } else {
                for (int i = 0; i < indexer.getNumberOfPoints(); i++) {
                    float a1 = array1[i];
                    float a2 = array2[i];
                    assertTrue(a1 == a2);
                }
            }
        }
    }
}
