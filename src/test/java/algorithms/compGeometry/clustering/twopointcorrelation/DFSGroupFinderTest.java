package algorithms.compGeometry.clustering.twopointcorrelation;

import java.util.logging.Logger;
import algorithms.util.PolygonAndPointPlotter;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class DFSGroupFinderTest extends TestCase {

    protected Logger log = Logger.getLogger(this.getClass().getSimpleName());

    public DFSGroupFinderTest(String testName) {
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
  
    public void testFindClusters() throws Exception {
        
        float[] x = new float[] {
                                5,  6,  7,  8,
                                5,  6,  7,  8,
                                5,  6,  7,  8,
                                5,  6,  7,  8,
            0,  1,  2,  3,      5,  6,  7,  8,
            0,  1,  2,  3,                 
            0,  1,  2,  3,                 
            0,  1,  2,  3,                 
            0,  1,              5,  6,  7,  8,
            0,  1,              5,  6,  7,  8,
            0,  1,          4,  5,  6,  7,  8
        };
        float[] y = new float[] {
                                0,  0,  0,  0,
                                1,  1,  1,  1,
                                2,  2,  2,  2,
                                3,  3,  3,  3,
            4,  4,  4,  4,      4,  4,  4,  4,
            5,  5,  5,  5,                 
            6,  6,  6,  6,                 
            7,  7,  7,  7,                
            8,  8,              8,  8,  8,  8,
            9,  9,              9,  9,  9,  9,
           10, 10,          10, 10, 10, 10, 10
        };

        AxisIndexer indexer = new AxisIndexer();
        indexer.sortAndIndexX(x, y, x.length);

        
        float[] xmmm = indexer.findXYMinMax();
        
        float xmin = xmmm[0];
        float xmax = xmmm[1];
        float ymin = xmmm[2];
        float ymax = xmmm[3];
        
        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter(
            xmin, xmax, ymin, ymax
        );
        plotter.addPlot(indexer.getX(), indexer.getY(), null, null, "");
        String filePath = plotter.writeFile();
        log.info("filePath=" + filePath);
        
        DFSGroupFinder groupFinder = new DFSGroupFinder(2.0f, 1.0f);
        groupFinder.findGroups(indexer);
        
        SimpleLinkedListNode[] list = groupFinder.getGroupMembershipList();
        
        int nGroups = groupFinder.getNumberOfGroups();
        
        assertTrue(nGroups == 3);
        
        int[] expected = new int[]{22, 20, 13};
        boolean[] found = new boolean[3];
        
        log.info("nGroups=" + nGroups);
        
        for (int i = 0; i < nGroups; i++) {
            int count = 0;
            SimpleLinkedListNode latest = list[i];
            while (latest != null) {
                count++;
                latest = latest.next;
            }
            log.info("group " + i + " count=" + count);
            
            for (int j = 0; j < expected.length; j++) {
                if (expected[j] == count) {
                    found[j] = true;
                    break;
                }
            }
        }
        
        for (boolean f : found) {
            assertTrue(f);
        }
        
        assertNotNull(groupFinder.getPointToGroupIndexes());
        
        groupFinder.printMembership(indexer);
        
        // test some exceptions
        boolean threwException = false;
        try {
            float[] xa = groupFinder.getX(1, null);
        } catch (IllegalArgumentException e) {
            threwException = true;
        }
        assertTrue(threwException);
        
        
        threwException = false;
        try {
            float[] xa = groupFinder.getX(100, null);
        } catch (IllegalArgumentException e) {
            threwException = true;
        }
        assertTrue(threwException);
        
        
        threwException = false;
        try {
            float[] ya = groupFinder.getY(1, null);
        } catch (IllegalArgumentException e) {
            threwException = true;
        }
        assertTrue(threwException);
        
        
        threwException = false;
        try {
            float[] ya = groupFinder.getY(100, null);
        } catch (IllegalArgumentException e) {
            threwException = true;
        }
        assertTrue(threwException);
        
        
        threwException = false;
        try {
            int[] ia = groupFinder.getIndexes(100);
        } catch (IllegalArgumentException e) {
            threwException = true;
        }
        assertTrue(threwException);
    }
    
    public void testExceptions() throws Exception {
        
        boolean threwException;
        
        DFSGroupFinder groupFinder = new DFSGroupFinder(2.0f, 1.0f);
        groupFinder.setDebug(true);
        
        threwException = false;
        try {
            groupFinder.findGroups(null);
        } catch (IllegalArgumentException e) {
            threwException = true;
        }
        assertTrue(threwException);
        
        
        threwException = false;
        try {
            groupFinder.initializeVariables(null);
        } catch (IllegalArgumentException e) {
            threwException = true;
        }
        assertTrue(threwException);
        
        
        threwException = false;
        try {
            groupFinder.printMembership(null);
        } catch (IllegalArgumentException e) {
            threwException = true;
        }
        assertTrue(threwException);
        
        float[] xa = groupFinder.getX(1, null);
        assertTrue(xa.length == 0);
        
        assertTrue(groupFinder.getY(1, null).length == 0);
        
        assertTrue(groupFinder.getIndexes(1).length == 0);
        
    }

}
