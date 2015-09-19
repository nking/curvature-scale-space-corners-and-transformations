package algorithms.imageProcessing;

import algorithms.compGeometry.PerimeterFinder;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.ResourceFinder;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
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

            List<Set<Integer> > groups = finder.getGroupMembershipList();

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

            Set<Integer> indexList = groups.get(0);
            for (Integer index : indexList) {
                int idx = index.intValue();
                int x = img.getCol(idx);
                int y = img.getRow(idx);

                String v = String.format("%d:%d", x, y);
                int eIdx = expected.indexOf(v);
                assertTrue(eIdx > -1);
                expected.remove(eIdx);
            }

            assertTrue(expected.isEmpty());
        }
    }
    
    public void testEmbeddedPointsInBL2003Sky() throws Exception {
                
        //TODO update this
        
        String filePath = ResourceFinder.findFileInTestResources("test_mask_0.png");
        
        GreyscaleImage img = ImageIOHelper.readImageAsBinary(filePath);
        
        DFSContiguousValueFinder finder = new DFSContiguousValueFinder(img);

        finder.findGroups(0);
                
        assertTrue(finder.getNumberOfGroups() == 1);
        
    }

    private void plotPerimeters(int xMin, int xMax, int yMin, int yMax, 
        Map<Integer, PairInt> skyRowColRange, int[] skyRowMinMax, 
        Map<Integer, PairInt> gRowColRange, int[] gRowMinMax, int gId) throws IOException {
        
        Image img = new Image(xMax + 1, yMax + 1);
        
        PairIntArray group = new PairIntArray();
        
        for (int r = gRowMinMax[0]; r <= gRowMinMax[1]; r++) {
                    
            PairInt cRange = gRowColRange.get(Integer.valueOf(r));
            
            int y = r;
            
            group.add(cRange.getX(), y);
            group.add(cRange.getY(), y);
        }
        
        PairIntArray sky = new PairIntArray();
        for (int r = skyRowMinMax[0]; r <= skyRowMinMax[1]; r++) {
                    
            PairInt cRange = skyRowColRange.get(Integer.valueOf(r));
            
            int y = r;
            
            sky.add(cRange.getX(), y);
            sky.add(cRange.getY(), y);
        }
        
        ImageIOHelper.addCurveToImage(sky, img, 
            2, 0, 0, 255);
        ImageIOHelper.addCurveToImage(group, img, 
            1, 255, 0, 0);
        
        ImageDisplayer.displayImage(Integer.toString(gId), img);
    }
    
}
