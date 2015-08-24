package algorithms.imageProcessing;

import java.io.File;
import junit.framework.TestCase;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class CMDLineInflectionUtilTest extends TestCase {
    
    public CMDLineInflectionUtilTest() {
    }
    
    public void test() {
        
        // temporarily disabled until changes are in place
        if (true) {
            return;
        }
        
        String[] args = new String[]{
            "-image1=testresources/closed_curve.png",
            "-image2=testresources/closed_curve_translate_scale_rotate135.png",
            "-mark_image",
            "-input_is_line_drawing",
            "-include_edges",
            "-text_output",
            "-refine_transformations"
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
