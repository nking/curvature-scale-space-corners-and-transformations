package algorithms.connected;

import algorithms.util.ResourceFinder;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.ImageProcessor;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class ConnectedValuesFinderTest extends TestCase {
    
    public ConnectedValuesFinderTest() {
    }
    
    public void testFindGroups() throws Exception {

        int w = 3;
        int h = 3;
              
        GreyscaleImage img = new GreyscaleImage(w, h);
        
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
            
            ConnectedValuesFinder finder = new ConnectedValuesFinder(img);

            if (s == 0) {
                finder.findGroups(0);
            } else {
                finder.findGroupsNotThisValue(1);
            }

            int nGroups = finder.getNumberOfGroups();
            
            assertTrue(nGroups == 1);

            TIntSet xy = finder.getXY(0);

            assertNotNull(xy);

            assertTrue(xy.size() == 4);

            TIntSet expected = new TIntHashSet(4);
            expected.add((0 * w) + 0);
            expected.add((0 * w) + 1);
            expected.add((0 * w) + 2);
            expected.add((1 * w) + 1);

            TIntIterator iter = xy.iterator();
            while (iter.hasNext()) {
                int pixIdx = iter.next();
                assertTrue(expected.remove(pixIdx));
            }
            assertTrue(expected.isEmpty());
        }
    }
    
    public void testInBL2003Sky() throws Exception {
                
        //TODO update this
        
        String filePath = ResourceFinder.findFileInTestResources("test_mask_0.png");
        
        GreyscaleImage img = ImageIOHelper.readImageAsBinary(filePath);
        
        ConnectedValuesFinder finder = new ConnectedValuesFinder(img);

        finder.findGroups(0);
                
        assertTrue(finder.getNumberOfGroups() == 1);
    }
    
    public void testInBL2003Sky_2() throws Exception {
                
        //TODO update this
        
        String filePath = ResourceFinder.findFileInTestResources("test_mask_0.png");
        
        GreyscaleImage img = ImageIOHelper.readImageAsBinary(filePath);
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        int[] imgValues = imageProcessor.convertToInt(img);
        
        ConnectedValuesFinder finder = new ConnectedValuesFinder(imgValues,
            img.getWidth(), img.getHeight());

        finder.findGroups(0);
                
        assertTrue(finder.getNumberOfGroups() == 1);
    }

}
