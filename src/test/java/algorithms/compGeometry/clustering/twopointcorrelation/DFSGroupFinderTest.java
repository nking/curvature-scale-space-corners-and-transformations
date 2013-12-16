package algorithms.compGeometry.clustering.twopointcorrelation;

import java.security.SecureRandom;
import java.util.logging.Logger;

import algorithms.compGeometry.clustering.twopointcorrelation.RandomClusterAndBackgroundGenerator.CLUSTER_SEPARATION;
import algorithms.misc.MiscMath;
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


        DoubleAxisIndexer indexer = new DoubleAxisIndexer();
        indexer.sortAndIndexXThenY(x, y, x.length);
        
        float xmin = indexer.getX()[ indexer.getSortedXIndexes()[0] ];
        float xmax = indexer.getX()[ indexer.getSortedXIndexes()[indexer.getNumberOfPoints() - 1] ];
        float ymin = indexer.getY()[ indexer.getSortedYIndexes()[0] ];
        float ymax = indexer.getY()[ indexer.getSortedYIndexes()[indexer.getNumberOfPoints() - 1] ];
        
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
        
    }

}
