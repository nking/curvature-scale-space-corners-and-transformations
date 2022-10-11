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
        double[][] k1Intr = Camera.createIntrinsicCameraMatrix(830.3, 307, 206.6);
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
        
        Triangulation.WCSPt wcsPt = Triangulation.calculateWCSPoint(
            k1Intr, k1ExtrRot, k1ExtrTrans,
            k2Intr, k2ExtrRot, k2ExtrTrans,
            x1, x2);

        double[] xw = wcsPt.X;

        // expecting Z values between about 310 and 500
        //  can see that in their 3D diagram on page http://robots.stanford.edu/cs223b04/JeanYvesCalib/htmls/example5.html
        // NOTE that the cameras are located at Z of about 150.
        System.out.printf("xw=%s", FormatArray.toString(xw, "%.3e"));
        System.out.println("expecting z values between 310 and 500. " +
                "  The cameras are at Z ~ 150.");
            
    }
    
    /**
     * Test of calculateWCSPoint method, of class Triangulation.
     */
    public void testCalculateWCSPoint2() {
        
        //test data from:
        //http://www.vision.caltech.edu/bouguetj/calib_doc/htmls/example5.html
        // now at http://robots.stanford.edu/cs223b04/JeanYvesCalib/
        // "Fifth calibration example - Calibrating a stereo system, stereo image rectification and 3D stereo triangulation"
        // by Jean-Yves Bouguet
        // Camera Calibration Toolbox for Matlab
        //
        // left camera: 
        //    focal length = 533.5, 533.5
        //    cc = 341.6, 235.2
        //    skew = 0, 0
        //    radial distortion k = -0.288, 0.097, 0.001, -0.0003, 0
        //
        // right camera: 
        //    focal length = 536.8, 536.5
        //    cc = 326.3, 250.1
        //    skew = 0, 0
        //    radial distortion k = -0.289, 0.107, 0.001, -0.0001, 0
        //
        // rotation vector om=0.00669, 0.00452, -0.0035
        // translation vector t = -99.80198, 1.12443, 0.05041
        //
        // note: the checkerboard squares are 30mm in WCS metric
        //
        // note: radial distortion should be corrected:  use on the original coordinates:
        //     x_corrected = x*(1 + k1*r^2 + k2r^4) where r is distance of point from cc.
        
        double[] k1IntrTrans = new double[]{-341.6, -234.3, 0};
        double[] k2IntrTrans = new double[]{-326.5, -249.3, 0};
        double[] k1IntrTransInv = Arrays.copyOf(k1IntrTrans, k1IntrTrans.length);
        MatrixUtil.multiply(k1IntrTransInv, -1);
        double[] k2IntrTransInv = Arrays.copyOf(k2IntrTrans, k2IntrTrans.length);
        MatrixUtil.multiply(k2IntrTransInv, -1);
        
        double[][] k1Intr = Camera.createIntrinsicCameraMatrix(533.07, k1IntrTrans[0], k1IntrTrans[1]);
        double[][] k2Intr = Camera.createIntrinsicCameraMatrix(536.7, k2IntrTrans[0], k2IntrTrans[1]);
        double[][] k2IntrInv = Camera.createIntrinsicCameraMatrixInverse(536.7, k2IntrTrans[0], k2IntrTrans[1]);
        
        double[][] k1ExtrRot = MatrixUtil.createIdentityMatrix(3);
        double[] k1ExtrTrans = new double[]{0, 0, 0};
        
        double[][] k2ExtrRot = Rotation.createRodriguesFormulaRotationMatrix(
            new double[]{0.00611, 0.00409, -0.00359});
        double[] k2ExtrTrans = new double[]{-99.85, 0.82, 0.44};
        
        System.out.printf("k1ExtrRot\n=%s\n", FormatArray.toString(k1ExtrRot, "%.3e"));
        System.out.printf("k1ExtrTrans\n=%s\n\n", FormatArray.toString(k1ExtrTrans, "%.3e"));
        System.out.printf("k2ExtrRot\n=%s\n", FormatArray.toString(k2ExtrRot, "%.3e"));
        System.out.printf("k2ExtrTrans\n=%s\n\n", FormatArray.toString(k2ExtrTrans, "%.3e"));
        
        //choosing a point that is in left01.jpg and right01.jpg
        //(307, 159)  (184, 172)
        double[][] x1 = new double[3][1];
        x1[0] = new double[]{307};
        x1[1] = new double[]{159};
        x1[2] = new double[]{1};
        
        double[][] x2 = new double[3][1];
        x2[0] = new double[]{184};
        x2[1] = new double[]{172};
        x2[2] = new double[]{1};
          
        double[] x1c, x2c;
        
        x1c = MatrixUtil.extractColumn(x1, 0);
        x1c = MatrixUtil.add(x1c, k1IntrTrans);
        x2c = MatrixUtil.extractColumn(x2, 0);
        x2c = MatrixUtil.add(x2c, k2IntrTrans);
        System.out.printf("x1c=%s\n", FormatArray.toString(x1c, "%.3e"));
        System.out.printf("x2c=%s\n", FormatArray.toString(x2c, "%.3e"));
        
        double[] tmp;
        //check: Xc_1_right = R * Xc_1_left + T
        // can see need to subtract center of image 1 from coords then add center of image 2.
        // can also see that the extrinsic translation is in units of image pixels,
        tmp = MatrixUtil.multiplyMatrixByColumnVector(k2ExtrRot, x1c);
        tmp = MatrixUtil.add(tmp, k2ExtrTrans);
        System.out.printf("?? %s\n", FormatArray.toString(tmp, "%.4e"));
        tmp = MatrixUtil.add(tmp, k2IntrTransInv);
        System.out.printf("==> %s\n", FormatArray.toString(tmp, "%.4e")); //1.9179e+02, 1.7495e+02, 2.0417e+00
        MatrixUtil.multiply(tmp, 1./tmp[2]);
        System.out.printf("*?? %s\n", FormatArray.toString(tmp, "%.4e"));
        
        double[] trans = MatrixUtil.subtract(x2c, x1c);            
        System.out.printf("?? trans = %s\n", FormatArray.toString(trans, "%.4e"));
        //  result = -1.0790e+02, -2.0000e+00, 0.0000e+00
        
        // correct for radial distortions:
        /*double r2, distortion;
        for (int i = 0; i < x1[0].length; ++i) {
            //x_corrected = x*(1 + k1*r^2 + k2r^4)
            r2 = x1[0][i] * x1[0][i] + x1[1][i] * x1[1][i];
            distortion = (1 + -0.288*r2 + 0.097*r2*r2);
            x1[0][i] *= distortion;
            x1[1][i] *= distortion;
            r2 = x2[0][i] * x2[0][i] + x2[1][i] * x2[1][i];
            distortion = (1 + -0.289*r2 + 0.107*r2*r2);
            x2[0][i] *= distortion;
            x2[1][i] *= distortion;
        }*/

        Triangulation.WCSPt wcsPt = Triangulation.calculateWCSPoint(
            k1Intr, k1ExtrRot, k1ExtrTrans,
            k2Intr, k2ExtrRot, k2ExtrTrans,
            x1, x2);

        double[] xw = wcsPt.X;
        
        double[][] camera1 = Camera.createCamera(k1Intr, k1ExtrRot, k1ExtrTrans);
        double[][] camera2 = Camera.createCamera(k2Intr, k2ExtrRot, k2ExtrTrans);
        double[][] camera1Inv = Camera.createCameraInverse(k1Intr, k1ExtrRot, k1ExtrTrans);
        double[][] camera2Inv = Camera.createCameraInverse(k2Intr, k2ExtrRot, k2ExtrTrans);
        
        System.out.printf("xw=%s\n\n", FormatArray.toString(xw, "%.3e"));
        
        // very rough depth from disparity without rectification
        //   z = (baseline * focalLength) / (differences of x displacements from centers of their images)
        //     = (99*534)/((341-169.97)-(326.5-183.1))
        //     = 1913
        //
        // Also, xw=-6.414e-02, -1.390e-01, -9.882e-01, 1.933e-03
        //   xw[0]/xw[3], xw[1]/xw[3], xw[2]/xw[3] =
        //     = -33.18158303155717, -71.9089498189343, -511.2260734609415
        // The value 511 fits with their figure
    
        MatrixUtil.multiply(xw, 1./xw[xw.length - 1]);
        System.out.printf("   =%s\n\n", FormatArray.toString(xw, "%.3e"));
        
        double[] expectedx1 = MatrixUtil.multiplyMatrixByColumnVector(camera1, xw);
        double[] expectedx2 = MatrixUtil.multiplyMatrixByColumnVector(camera2, xw);
        MatrixUtil.multiply(expectedx1, 1./expectedx1[expectedx1.length - 1]);
        MatrixUtil.multiply(expectedx2, 1./expectedx2[expectedx2.length - 1]);
        
        System.out.printf("camera1*xw=%s\n", FormatArray.toString(expectedx1, "%.3e"));
        System.out.printf("camera2*xw=%s\n", FormatArray.toString(expectedx2, "%.3e"));
         
        double[][] expectedX1 = MatrixUtil.multiply(camera1Inv, x1);
        double[][] expectedX2 = MatrixUtil.multiply(camera2Inv, x2);
        
        System.out.printf("camera1^-1 * x1=\n%s\n", FormatArray.toString(expectedX1, "%.3e"));
        System.out.printf("camera2^-1 * x2=\n%s\n", FormatArray.toString(expectedX2, "%.3e"));
        
        assertTrue(Math.abs(Math.abs(xw[2]) - 425) < 150);
        
        //a correspondence closer to middle of image 1: (345, 188)  (215,238)
        x1 = new double[3][1];
        x1[0] = new double[]{345};
        x1[1] = new double[]{188};
        x1[2] = new double[]{1};
        
        // line 3 last point in image 2
        x2 = new double[3][1];
        x2[0] = new double[]{215};
        x2[1] = new double[]{238};
        x2[2] = new double[]{1};
                  
        x1c = MatrixUtil.extractColumn(x1, 0);
        x1c = MatrixUtil.add(x1c, k1IntrTrans);
        x2c = MatrixUtil.extractColumn(x2, 0);
        x2c = MatrixUtil.add(x2c, k2IntrTrans);
        System.out.printf("x1c=%s\n", FormatArray.toString(x1c, "%.3e"));
        System.out.printf("x2c=%s\n", FormatArray.toString(x2c, "%.3e"));
        
        //check: Xc_1_right = R * Xc_1_left + T
        tmp = MatrixUtil.multiplyMatrixByColumnVector(k2ExtrRot, x1c);
        tmp = MatrixUtil.add(tmp, k2ExtrTrans);
        System.out.printf("?? %s\n", FormatArray.toString(tmp, "%.4e"));
        tmp = MatrixUtil.add(tmp, k2IntrTransInv);
        System.out.printf("==> %s\n", FormatArray.toString(tmp, "%.4e"));
        MatrixUtil.multiply(tmp, 1./tmp[2]);
        System.out.printf("*?? %s\n", FormatArray.toString(tmp, "%.4e"));

        wcsPt = Triangulation.calculateWCSPoint(
            k1Intr, k1ExtrRot, k1ExtrTrans,
            k2Intr, k2ExtrRot, k2ExtrTrans,
            x1, x2);

        xw = wcsPt.X;
        
        System.out.printf("xw=%s\n\n", FormatArray.toString(xw, "%.3e"));
        
        //xw=6.670e-03, -5.187e-02, -9.986e-01, 2.102e-03
        // xw[0]/xw[3], xw[1]/xw[3], xw[2]/xw[3] = (3.173168411037107, -24.676498572787818, -475.07136060894385)
    
        MatrixUtil.multiply(xw, 1./xw[xw.length - 1]);
        System.out.printf("   =%s\n\n", FormatArray.toString(xw, "%.3e"));
        
        expectedX1 = MatrixUtil.multiply(camera1Inv, x1);
        expectedX2 = MatrixUtil.multiply(camera2Inv, x2);
        
        System.out.printf("camera1^-1 * x1=\n%s\n", FormatArray.toString(expectedX1, "%.3e"));
        System.out.printf("camera2^-1 * x2=\n%s\n", FormatArray.toString(expectedX2, "%.3e"));
        
        assertTrue(Math.abs(Math.abs(xw[2]) - 425) < 150);
    }
}
