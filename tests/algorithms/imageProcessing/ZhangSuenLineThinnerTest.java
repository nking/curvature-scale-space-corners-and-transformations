package algorithms.imageProcessing;

import algorithms.util.PairInt;
import java.util.HashSet;
import java.util.Set;
import java.util.logging.Logger;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class ZhangSuenLineThinnerTest extends TestCase {

    private Logger log = Logger.getLogger(this.getClass().getName());
    
    public ZhangSuenLineThinnerTest(String testName) {
        super(testName);
    }

    @Override
    protected void setUp() throws Exception {
        super.setUp();
    }

    @Override
    protected void tearDown() throws Exception {
        super.tearDown();
    }
   
    public void testApplyLineThinner() throws Exception {

        // adapted from code at http://rosettacode.org/wiki/Zhang-Suen_thinning_algorithm#Java
        
        String[] image = new String[18];
        image[0]  = "                                                          "; //0
        image[1]  = " #################                   #############        "; //1
        image[2]  = " ##################               ################        "; //2
        image[3]  = " ###################            ##################        "; //3
        image[4]  = " ########     #######          ###################        "; //4
        image[5]  = "   ######     #######         #######       ######        "; //5
        image[6]  = "   ######     #######        #######                      "; //6
        image[7]  = "   #################         #######                      "; //7
        image[8]  = "   ################          #######                      "; //8
        image[9]  = "   #################         #######                      "; //9
        image[10] = "   ######     #######        #######                      "; //10
        image[11] = "   ######     #######        #######                      "; //11
        image[12] = "   ######     #######         #######       ######        "; //12
        image[13] = " ########     #######          ###################        "; //13
        image[14] = " ########     ####### ######    ################## ###### "; //14
        image[15] = " ########     ####### ######      ################ ###### "; //15
        image[16] = " ########     ####### ######         ############# ###### "; //16
        image[17] = "                                                          ";//17
        //012345678901234567890123456789012345678901234567890123456789
        //          1         2         3         4         5
        
        Set<PairInt> points = new HashSet<PairInt>(image.length);
 
        for (int r = 0; r < image.length; r++) {
            char[] cs = image[r].toCharArray();
            for (int c = 0; c < cs.length; c++) {
                if (cs[c] != ' ') {
                    PairInt p = new PairInt(c, r);
                    points.add(p);
                }
            }
        }
        
        ZhangSuenLineThinner zsLT = new ZhangSuenLineThinner();
        
        zsLT.applyLineThinner(points, 0, 80, 0, image.length - 1);
        
        for (int r = 0; r < image.length; r++) {
            
            StringBuffer sb = new StringBuffer();
            
            for (int c = 0; c < image[r].toCharArray().length; c++) {
                
                PairInt p = new PairInt(c, r);
                
                if (points.contains(p)) {
                    sb.append("#");
                } else {
                    sb.append(" ");
                }
            }
            
            System.out.println(sb.toString());
        }
    }
   
}
