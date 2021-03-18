package algorithms.imageProcessing.transform;

import algorithms.matrix.MatrixUtil;
import algorithms.util.FormatArray;
import java.util.Arrays;
import junit.framework.TestCase;
import org.junit.Test;

/**
 *
 * @author nichole
 */
public class TriangulationTest extends TestCase {
    
    public TriangulationTest() {
    }

    public void testCreateIntrinsicCameraMatrix() {
    }

    /**
     * Test of calculateWCSPoint method, of class Triangulation.
     */
    public void _testCalculateWCSPoint() {
        /*
        test data from:
        https://www.microsoft.com/en-us/research/project/a-flexible-new-technique-for-camera-calibration-2/
        
        The publication is Zhengyou Zhang, 1998
        "A Flexible New Technique for Camera Calibration" 
        https://www.microsoft.com/en-us/research/wp-content/uploads/2016/02/tr98-71.pdf
        
        The model plane contains a pattern of 8×8 squares, so there are 256 corners. 
        The size of the pattern is 17cm x 17cm. The 2D coordinates (in inches) 
        of these points are available here. (We assume the plane is at Z=0.)
        
        We have taken five an off-the-shelf PULNiX CCD camera with 6 mm lens. 
        The image resolution is 640×480.
        
        There is significant lens distortion in the images.
        
        And here is what the calibration tells us about the camera: 
        The pixel is square (aspect ratio = 1); the focal length = 832.5 pixels; 
        the image center is at (303.959, 206.585); there is a significant radial 
        distortion: k1 = -0.228601, k2 = 0.190353. 
        
        The format of the calibration file is: 
            a, c, b, u0, v0, k1, k2, then 
            the rotation matrix and translation vector for the first image, 
            the rotation matrix and translation vector for the second image, etc.
            where 
              a = focal length in x, or rather a scale factor
              c = skew
              b = focal length in y, or rather a scale factor
              u0 = image center x
              v0 = image center y
              k1 = radial distortion
              k2 = radial distortion
              (Tables 1 and 2 of the publication have more details)
        
            832.5 0.204494 832.53 303.959 206.585

            -0.228601 0.190353

            0.992759 -0.026319 0.117201
            0.0139247 0.994339 0.105341
            -0.11931 -0.102947 0.987505
            -3.84019 3.65164 12.791

            0.997397 -0.00482564 0.0719419
            0.0175608 0.983971 -0.17746
            -0.0699324 0.178262 0.981495
            -3.71693 3.76928 13.1974

            0.915213 -0.0356648 0.401389
            -0.00807547 0.994252 0.106756
            -0.402889 -0.100946 0.909665
            -2.94409 3.77653 14.2456

            0.986617 -0.0175461 -0.16211
            0.0337573 0.994634 0.0977953
            0.159524 -0.101959 0.981915
            -3.40697 3.6362 12.4551

            0.967585 -0.196899 -0.158144
            0.191542 0.980281 -0.0485827
            0.164592 0.0167167 0.98622
            -4.07238 3.21033 14.3441
        
        The point data seem to in format:
            each line is the 4 corners of a square as x, y pairs.
        
        Here are the first 3 lines of each image set:
        
        Image 1:
        63.43921044061905 405.57679766845445    92.46270141677354 407.4556539075571    91.80636571669007 438.65765085408424    62.58724663945761 436.28844212118605    
116.28035530429925 409.17858333240645    146.4502700963233 410.92321756063495    145.61153657396474 442.5588539245722    115.4621289243657 440.2901448310849    
170.832696837639 412.46038130831136    201.84874167693556 414.36697339157973    201.07282342562434 446.1786800590262    169.96948547603384 444.042108267029
        
        Image 2:
        74.95173767129486 409.09267757581756    104.31259784448483 410.8263322391833    106.06048939567938 438.7099113674728    77.00428944064123 436.94253844835623    
127.82252063924655 412.2041930215117    157.94718654000116 413.6707476599427    159.23134690487825 442.0125923370657    129.4639745161323 439.99858857260546    
181.84295543141715 414.98756134254074    212.70112359460154 416.543317839545    213.47114738670456 444.7447412257028    183.11454542957057 443.2278534824685
        
        Image 3:
        137.22826265754128 394.36338179898917    160.187943549149 397.0974538410207    159.34730413574266 425.12927320814157    136.29854255272394 421.95770763177325    
180.16912538559387 399.3808729240103    204.86292900573017 402.06267557016076    203.72675697709863 431.01683243234135    179.18487675331212 427.73591150424966    
225.75104300713795 404.3630066605049    251.92825010658618 407.17957914453706    250.96538074145894 436.9423114376981    224.6283650625028 433.54478632758776
        
        
        */
        double[][] k1Intr = Triangulation.createIntrinsicCameraMatrix(830.3, 307, 206.6);
        double[][] k2Intr = Arrays.copyOf(k1Intr, k1Intr.length);
        
        double[][] k1ExtrRot = new double[3][3];
        k1ExtrRot[0] = new double[]{0.992759, -0.026319, 0.117201};
        k1ExtrRot[1] = new double[]{0.0139247, 0.994339, 0.105341};
        k1ExtrRot[2] = new double[]{-0.11931, -0.102947, 0.987505};
        double[] k1ExtrTrans = new double[]{-3.84019, 3.65164, 12.791};
        
        double[][] k2ExtrRot = new double[3][3];
        k2ExtrRot[0] = new double[]{0.997397, -0.00482564, 0.0719419};
        k2ExtrRot[1] = new double[]{0.0175608, 0.983971, -0.17746};
        k2ExtrRot[2] = new double[]{-0.0699324, 0.178262, 0.981495};
        double[] k2ExtrTrans = new double[]{-3.71693, 3.76928, 13.1974};
        
        // line 3 last point in image 1
        double[][] x1 = new double[3][1];
        x1[0] = new double[]{169.96948547603384};
        x1[1] = new double[]{444.042108267029};
        x1[2] = new double[]{1};
        
        // line 3 last point in image 2
        double[][] x2 = new double[3][1];
        x2[0] = new double[]{183.11454542957057};
        x2[1] = new double[]{443.2278534824685};
        x2[2] = new double[]{1};
        
        double[] xw = Triangulation.calculateWCSPoint(
            k1Intr, k1ExtrRot, k1ExtrTrans,
            k2Intr, k2ExtrRot, k2ExtrTrans,
            x1, x2);
        
        System.out.printf("xw=%s", FormatArray.toString(xw, "%.3e"));
            
            
    }
    
    /**
     * Test of calculateWCSPoint method, of class Triangulation.
     */
    public void testCalculateWCSPoint2() {
        /*
        test data from:
        http://www.vision.caltech.edu/bouguetj/calib_doc/htmls/example5.html
        
        choosing a point that is in left01.jpg and right01.jpg
        (307, 159)  (184, 172)
        */
        double[][] k1Intr = Triangulation.createIntrinsicCameraMatrix(-533.07, 341.6, 234.3);
        double[][] k2Intr = Triangulation.createIntrinsicCameraMatrix(-536.7, 326.5, 249.3);
        
        double[][] k1ExtrRot = MatrixUtil.createIdentityMatrix(3);
        double[] k1ExtrTrans = new double[]{0, 0, 1}; // check this
        
        double[][] k2ExtrRot = Triangulation.createRodriguesFormulaRotationMatrix(
            new double[]{0.00611, 0.00409, -0.00359});
        double[] k2ExtrTrans = new double[]{-99.85, 0.82, 0.44};
        
        System.out.printf("k2ExtrRot\n=%s\n", FormatArray.toString(k2ExtrRot, "%.3e"));
        
        // line 3 last point in image 1
        double[][] x1 = new double[3][1];
        x1[0] = new double[]{169.96948547603384};
        x1[1] = new double[]{444.042108267029};
        x1[2] = new double[]{1};
        
        // line 3 last point in image 2
        double[][] x2 = new double[3][1];
        x2[0] = new double[]{183.11454542957057};
        x2[1] = new double[]{443.2278534824685};
        x2[2] = new double[]{1};
        
        double[] xw = Triangulation.calculateWCSPoint(
            k1Intr, k1ExtrRot, k1ExtrTrans,
            k2Intr, k2ExtrRot, k2ExtrTrans,
            x1, x2);
        
        System.out.printf("xw=%s\n", FormatArray.toString(xw, "%.3e"));
            
            
    }
}
