package algorithms.imageProcessing;

import algorithms.imageProcessing.transform.TransformationParameters;
import algorithms.imageProcessing.transform.Transformer;
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
            
            StringBuilder sb = new StringBuilder();
            
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
    
    public void testApplyLineThinner2() throws Exception {

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
        
        GreyscaleImage img = new GreyscaleImage(image.length, image[0].length());
        for (int i = 0; i < image.length; ++i) {
            String line = image[i];
            for (int j = 0; j < line.length(); ++j) {
                if (line.charAt(j) == '#') {
                    img.setValue(i, j, 255);
                }
            }
        }
        
        printImage(img);
        
        TransformationParameters params = new TransformationParameters();
        params.setScale(1);
        params.setTranslationX(0);
        params.setTranslationY(0);
        params.setOriginX(0);
        params.setOriginY(0);
        
        params.setRotationInDegrees(270);
        params.setTranslationX(img.getHeight());
        
        //params.setRotationInDegrees(180);
        //params.setTranslationX(img.getWidth());
        //params.setTranslationY(img.getHeight());
        
        Transformer transformer = new Transformer();
        img = transformer.applyTransformation(img, params, 
                img.getHeight(), img.getWidth());
        //        img.getWidth(), img.getHeight());
 
        ImageProcessor imageProcessor = new ImageProcessor();
        
        imageProcessor.applyThinning(img);
        
        for (int r = (img.getHeight() - 1); r > -1; r--) {
            StringBuilder sb = new StringBuilder();            
            for (int c = 0; c < img.getWidth(); c++) {
                int v = img.getValue(c, r);                
                if (v != 0) {
                    sb.append("#");
                } else {
                    sb.append(" ");
                }
            }
            System.out.println(sb.toString());
        }
    }
   
    private void printImage(GreyscaleImage img) {
        for (int r = (img.getHeight() - 1); r > -1; r--) {
            StringBuilder sb = new StringBuilder();
            for (int c = 0; c < img.getWidth(); c++) {
                int v = img.getValue(c, r);                
                if (v != 0) {
                    sb.append("#");
                } else {
                    sb.append(" ");
                }
            }
            System.out.println(sb.toString());
        }
    }
}
