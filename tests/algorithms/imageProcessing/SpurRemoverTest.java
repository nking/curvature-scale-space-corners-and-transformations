package algorithms.imageProcessing;

import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class SpurRemoverTest extends TestCase {
    
    public SpurRemoverTest() {
    }

    public void testRemove_GreyscaleImage() {

        /*  \       /    8
             \     /     7
              - - -      6
              |   |      5
              |   |      4
              - - - - -  3
              |          2
              |          1
                         0
        0 1 2 3 4 5 6 7
        */
        GreyscaleImage img = new GreyscaleImage(10, 10);
        img.setValue(3, 1, 100);  img.setValue(3, 2, 100);
        img.setValue(3, 3, 100); img.setValue(4, 3, 100);
        img.setValue(5, 3, 100); img.setValue(6, 3, 100);
        img.setValue(7, 3, 100);
        img.setValue(3, 4, 100); img.setValue(5, 4, 100);
        img.setValue(3, 5, 100); img.setValue(5, 5, 100);
        img.setValue(3, 6, 100); img.setValue(4, 6, 100);
        img.setValue(5, 6, 100);
        img.setValue(2, 7, 100); 
        img.setValue(6, 7, 100);
        img.setValue(1, 8, 100); 
        img.setValue(7, 8, 100);
        
        SpurRemover spurRemover = new SpurRemover();
        spurRemover.remove(img);
        /*  \       /    8
             \     /     7
              - - -      6
              |   |      5
              |   |      4
              - - - - -  3
              |          2
              |          1
                         0
        0 1 2 3 4 5 6 7
        */
        assertEquals(0,img.getValue(3, 1));  
        assertEquals(0,img.getValue(3, 2));
        assertEquals(0,img.getValue(3, 3)); 
        assertEquals(100,img.getValue(4, 3));
        assertEquals(0,img.getValue(5, 3)); 
        assertEquals(0,img.getValue(6, 3));
        assertEquals(0,img.getValue(7, 3));
        assertEquals(100,img.getValue(3, 4)); 
        assertEquals(100,img.getValue(5, 4));
        assertEquals(100,img.getValue(3, 5)); 
        assertEquals(100,img.getValue(5, 5));
        assertEquals(0,img.getValue(3, 6)); 
        assertEquals(100,img.getValue(4, 6));
        assertEquals(0,img.getValue(5, 6));
        assertEquals(0,img.getValue(2, 7)); 
        assertEquals(0,img.getValue(6, 7));
        assertEquals(0,img.getValue(1, 8)); 
        assertEquals(0,img.getValue(7, 8));
    }

}
