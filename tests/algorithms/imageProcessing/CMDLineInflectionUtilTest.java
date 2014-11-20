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
public class CMDLineInflectionUtilTest {
    
    public CMDLineInflectionUtilTest() {
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
            "-image2=testresources/closed_curve_translate_scale_rotate135.png",
            "-mark_image",
            "-include_edges",
            "-text_output"
            //,"-input_is_line_drawing",
            
        };
        
        CMDLineInflectionUtil.main(args);
        
        String cwd = System.getProperty("user.dir") 
            + System.getProperty("file.separator");
        
        String[] expectedFile = new String[] {
            "transformed_closed_curve.png",
            "inflection_points_closed_curve_translate_scale_rotate135.png",
            "inflection_points_closed_curve.png",
            "transformed_edges_closed_curve_closed_curve_translate_scale_rotate135.tsv",
            "transformation_closed_curve_closed_curve_translate_scale_rotate135.txt",
            "matching_points_closed_curve_closed_curve_translate_scale_rotate135.tsv"
        };

        for (String efl : expectedFile) {
            File fl = new File(cwd + efl);
            assertTrue(fl.exists());
            assertTrue(fl.delete());
        }
    }

    public static void main(String[] args) {
        
        try {
            
            CMDLineInflectionUtilTest test = 
                new CMDLineInflectionUtilTest();
            
            test.test();
            
        } catch (Exception e) {
            e.printStackTrace();
            System.err.println("ERROR: " + e.getMessage());
        }
    }
}
