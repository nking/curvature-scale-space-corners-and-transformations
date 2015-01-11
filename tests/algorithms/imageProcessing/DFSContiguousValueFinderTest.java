package algorithms.imageProcessing;

import algorithms.util.PairIntArray;
import algorithms.util.SimpleLinkedListNode;
import java.util.ArrayList;
import java.util.List;
import junit.framework.TestCase;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class DFSContiguousValueFinderTest extends TestCase {

    public DFSContiguousValueFinderTest(String testName) {
        super(testName);
    }
  
    public void testFindGroups() throws Exception {
        
        GreyscaleImage img = new GreyscaleImage(3, 3);
        
        /*
         0 0 0       0 1 2
         1 0 1       3 4 5
         1 1 1       6 7 8
        */
        img.setValue(0, 1, 1);
        img.setValue(0, 2, 1);
        img.setValue(1, 2, 1);
        img.setValue(2, 1, 1);
        img.setValue(2, 2, 1);
        
        for (int s = 0; s < 2; s++) {
            
            DFSContiguousValueFinder finder = new DFSContiguousValueFinder(img);

            if (s == 0) {
                finder.findGroups(0);
            } else {
                finder.findGroupsNotThisValue(1);
            }

            SimpleLinkedListNode[] groups = finder.getGroupMembershipList();

            assertNotNull(groups);

            assertTrue(finder.getNumberOfGroups() == 1);

            PairIntArray xy = finder.getXY(0);

            assertNotNull(xy);

            assertTrue(xy.getN() == 4);

            List<String> expected = new ArrayList<String>();
            expected.add(String.format("%d:%d", 0, 0));
            expected.add(String.format("%d:%d", 1, 0));
            expected.add(String.format("%d:%d", 2, 0));
            expected.add(String.format("%d:%d", 1, 1));

            for (int i = 0; i < xy.getN(); i++) {
                String v = String.format("%d:%d", xy.getX(i), xy.getY(i));
                int eIdx = expected.indexOf(v);
                assertTrue(eIdx > -1);
                expected.remove(eIdx);
            }

            assertTrue(expected.isEmpty());

            expected = new ArrayList<String>();
            expected.add(String.format("%d:%d", 0, 0));
            expected.add(String.format("%d:%d", 1, 0));
            expected.add(String.format("%d:%d", 2, 0));
            expected.add(String.format("%d:%d", 1, 1));

            SimpleLinkedListNode indexList = groups[0];
            while (indexList != null) {
                int idx = indexList.getKey();
                int x = img.getCol(idx);
                int y = img.getRow(idx);

                String v = String.format("%d:%d", x, y);
                int eIdx = expected.indexOf(v);
                assertTrue(eIdx > -1);
                expected.remove(eIdx);

                indexList = indexList.getNext();
            }

            assertTrue(expected.isEmpty());
        }
    }
}
