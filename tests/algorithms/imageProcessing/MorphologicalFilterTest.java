package algorithms.imageProcessing;

import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class MorphologicalFilterTest extends TestCase {
    
    public MorphologicalFilterTest() {
    }

    public void testBwMorphThin() {
        
        System.out.println("bwMorphThin");
        
        int[][] img = new int[4][];
        img[0] = new int[]{0, 1, 1, 1, 0};
        img[1] = new int[]{0, 1, 1, 1, 0};
        img[2] = new int[]{0, 1, 1, 1, 0};
        img[3] = new int[]{0, 1, 1, 1, 0};
        
        MorphologicalFilter mFilter = new MorphologicalFilter();
        int[][] skel = mFilter.bwMorphThin(img, Integer.MAX_VALUE);
     
        int[][] expected = new int[4][];
        expected[0] = new int[]{0, 0, 0, 0, 0};
        expected[1] = new int[]{0, 0, 1, 0, 0};
        expected[2] = new int[]{0, 0, 1, 0, 0};
        expected[3] = new int[]{0, 0, 0, 0, 0};
        
        for (int i = 0; i < expected.length; ++i) {
            for (int j = 0; j < expected[i].length; ++j) {
                assertEquals(expected[i][j], skel[i][j]);
            }
        }
    }
    
}
