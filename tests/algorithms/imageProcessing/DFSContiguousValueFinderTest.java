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
                
        String filePath = ResourceFinder.findFileInTestResources("test_mask_0.png");
        
        GreyscaleImage img = ImageIOHelper.readImageAsBinary(filePath);
        
        DFSContiguousValueFinder finder = new DFSContiguousValueFinder(img);

        finder.findGroups(0);
                
        assertTrue(finder.getNumberOfGroups() == 1);
        
        Set<PairInt> skyPoints = new HashSet<PairInt>();
        finder.getXY(0, skyPoints);
        
        finder = new DFSContiguousValueFinder(img);

        finder.findGroups(1);
                
        assertTrue(finder.getNumberOfGroups() == 5);
        
        // ===== get perimeter of sky ====
        PerimeterFinder perimeterFinder = new PerimeterFinder();
        int[] skyRowMinMax = new int[2];
        Map<Integer, PairInt> skyRowColRange = perimeterFinder.find(skyPoints, 
            skyRowMinMax);
        
        // ===== get contiguous groups of points embedded in sky perimeter
        //       that are not in skyPoints
        finder = new DFSContiguousValueFinder(img);
        
        finder.findEmbeddedGroupsNotThisValue(0, skyRowColRange, skyRowMinMax);
        
        int nEmbeddedGroups = finder.getNumberOfGroups();
        
        List<Set<PairInt> > embeddedGroups = new ArrayList<Set<PairInt> >();
        for (int i = 0; i < nEmbeddedGroups; i++) {
            Set<PairInt> set = new HashSet<PairInt>();
            finder.getXY(i, set);
            embeddedGroups.add(set);
        }
                
        Collections.sort(embeddedGroups, new SizeComparator<Set<PairInt>>());
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        // the first 4 in embeddedGroups should be found as bound by the sky
        //    perimeter
        // while the rest should not.  the rest are concave foreground
        //    features extending into the bottom perimeter of a convex hull 
        //    around sky points.
        
        // === test the bounds method =====
        for (int gIdx = 0; gIdx < embeddedGroups.size(); gIdx++) {
            
            int[] gRowMinMax = new int[2];
            Map<Integer, PairInt> gRowColRange = perimeterFinder.find(
                embeddedGroups.get(gIdx), gRowMinMax);
        
            // plot these
            // ALSO, isPerimeterUnbound needs to compare to image
            // boundaries with == not the sky perimeter.
            // ALSO, isPerimeterUnbound needs to know when the points are
            //    aligned with sky perimeter when that's not image boundary.
            //    (currently it thinks equals with sky perimeter is embedded
            //    but that might not be true because the embedded points aren't counted beyond?
            
            plotPerimeters(0, img.getWidth() - 1, 0, img.getHeight() - 1,
                skyRowColRange, skyRowMinMax, gRowColRange, gRowMinMax, gIdx);
            
            boolean unbound = imageProcessor.isPerimeterUnbound(
                gRowColRange, gRowMinMax, 
                 skyRowColRange, skyRowMinMax,
                 0, img.getWidth() - 1, 0, img.getHeight() - 1);
            
            assertFalse(unbound);
        }

        // points that should be on or within the sky perimeter
        
        Set<PairInt> embeddedPoints = new HashSet<PairInt>();
        imageProcessor.extractEmbeddedGroupPoints(
            embeddedGroups, skyRowColRange, skyRowMinMax, embeddedPoints,
            0, img.getWidth() - 1, 0, img.getHeight() - 1);
        
        int z = 1;
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
    
    private class SizeComparator<T extends  Collection> implements Comparator<T> {

        @Override
        public int compare(T o1, T o2) {
            
            if (o1 == null && o2 != null) {
                return 1;
            } 
            if (o1 != null && o2 == null) {
                return -1;
            }
            
            return Integer.compare(o1.size(), o2.size());
        }

    }
}
