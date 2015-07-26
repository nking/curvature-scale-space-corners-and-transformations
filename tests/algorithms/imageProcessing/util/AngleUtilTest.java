package algorithms.imageProcessing.util;

import algorithms.util.ResourceFinder;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import junit.framework.TestCase;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class AngleUtilTest extends TestCase {
    
    public AngleUtilTest() {
    }
    
    public void testSubtract() throws Exception {
        
        AngleUtil instance = new AngleUtil();
                        
        BufferedReader br = null;
        FileReader reader = null;
        
        //1 1 1.394857 0.343024 5.231353 (8.000000, -45.000000, 56.000000, -20.000000)
        
        String pattern = "(\\d)\\s(\\d)\\s([\\-0-9\\.]+)\\s([\\-0-9\\.]+)\\s" +
            "([\\-0-9\\.]+)\\s" 
            + "\\(" 
            + "([\\-\\.0-9]+),\\s([\\-\\.0-9]+),\\s([\\-\\.0-9]+),\\s([\\-\\.0-9]+)" 
            + "\\)";
        Pattern p = Pattern.compile(pattern);
        Matcher m = null;
        
        String filePath1 = ResourceFinder.findFileInTestResources("angles.tsv");
        try {
            reader = new FileReader(new File(filePath1));
            br = new BufferedReader(reader);
                        
            String line = br.readLine();
            while (line != null) {                
                m = p.matcher(line);
                if (m.matches()) {
                    double diffX1 = Double.valueOf(m.group(6)).doubleValue();
                    double diffY1 = Double.valueOf(m.group(7)).doubleValue();
                    double diffX2 = Double.valueOf(m.group(8)).doubleValue();
                    double diffY2 = Double.valueOf(m.group(9)).doubleValue();
                    
                    double expected = Double.valueOf(m.group(5)).doubleValue();
                    
                    double result = instance.subtract(diffX1, diffY1, diffX2, 
                        diffY2);
                    
                    System.out.println("line=" + line);
                    System.out.println("expected=" + expected + " result=" 
                        + result);
                    
                    assertTrue(Math.abs(expected - result) < 0.1);
                }
                
                line = br.readLine();
            }
            
        } finally {
            if (reader != null) {
                reader.close();
            }
            if (br != null) {
                br.close();
            }
        }
        
    }
    
    public void testSubtract2() throws Exception {
        
        /*
                  +Y
                 90
        QII       |       QI
                  |     
                  |
     180-------------------- +X  0, 360
                  |   
                  |      
         QIII     |       QIV 
                 270
        */
        double x, y, expected, result;
        
        x = 10;
        y = 10;
        expected = 45.*Math.PI/180.;
        result = AngleUtil.polarAngleCCW(x, y);
        assertTrue(Math.abs(result - expected) < 0.1);
        
        x = -10;
        y = 10;
        expected = (45. + 90.)*Math.PI/180.;
        result = AngleUtil.polarAngleCCW(x, y);
        assertTrue(Math.abs(result - expected) < 0.1);
        
        x = -10;
        y = -10;
        expected = (45. + 180.)*Math.PI/180.;
        result = AngleUtil.polarAngleCCW(x, y);
        assertTrue(Math.abs(result - expected) < 0.1);
        
        x = 10;
        y = -10;
        expected = (45. + 270.)*Math.PI/180.;
        result = AngleUtil.polarAngleCCW(x, y);
        assertTrue(Math.abs(result - expected) < 0.1);
    }
}
