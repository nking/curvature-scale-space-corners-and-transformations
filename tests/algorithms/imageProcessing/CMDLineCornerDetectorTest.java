package algorithms.imageProcessing;

import java.io.File;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class CMDLineCornerDetectorTest {
    
    public CMDLineCornerDetectorTest() {
    }
    
    @Before
    public void setUp() {
    }
    
    @After
    public void tearDown() {
    }

    @Test
    public void test() {
        
        String[] args = new String[]{
            "-image1=testresources/closed_curve.png",
            "-mark_image",
            "-include_edges",
            "-text_output",
            "-input_is_line_drawing",
            
        };
        
        CMDLineCornerDetector.main(args);
        
        String cwd = System.getProperty("user.dir") 
            + System.getProperty("file.separator");
        
        String[] expectedFile = new String[] {
            "corners_closed_curve.tsv",
            "corners_closed_curve.png"
        };

        for (String efl : expectedFile) {
            File fl = new File(cwd + efl);
            assertTrue(fl.exists());
            assertTrue(fl.delete());
        }
    }

}
