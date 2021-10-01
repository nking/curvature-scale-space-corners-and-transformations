package algorithms.imageProcessing.transform;

import static algorithms.imageProcessing.transform.Rotation.extractThetaFromZYX;
import algorithms.matrix.MatrixUtil;
import algorithms.util.FormatArray;
import java.util.Arrays;
import junit.framework.TestCase;
import static junit.framework.TestCase.assertNotNull;
import no.uib.cipr.matrix.NotConvergedException;
import org.junit.Test;

/**
 *
 * @author nichole
 */
public class ReconstructionTest extends TestCase {
    
    public ReconstructionTest() {
    }
    
      /**
     * Test of calculateUsingEssentialMatrix method, of class CameraPose.
     */
    public void testCalculateUsingEssentialMatrix() throws Exception {
        System.out.println("calculateUsingEssentialMatrix");
       
        double[][] k1 = Zhang98Data.getIntrinsicCameraMatrix();
        double[][] k2 = MatrixUtil.copy(k1);
        //x1, x2 size is 3 X 256
        double[][] x1 = Zhang98Data.getObservedFeaturesInImage(1);
        double[][] x2 = Zhang98Data.getObservedFeaturesInImage(5);
                
        Reconstruction.ReconstructionResults result = 
            Reconstruction.calculateUsingEssentialMatrix(k1, k2, x1, x2);
        
        assertNotNull(result);
        
        System.out.printf("\nresult:\nrot1=%strans1=%s\n", 
            FormatArray.toString(result.k1ExtrRot, "%.4e"),
            FormatArray.toString(result.k1ExtrTrans, "%.4e"));
        System.out.printf("\nresult:\nrot2=%strans2=%s\n", 
            FormatArray.toString(result.k2ExtrRot, "%.4e"),
            FormatArray.toString(result.k2ExtrTrans, "%.4e"));
        
        System.out.printf("\nimg1:\nrot=%strans=%s\n", 
                FormatArray.toString(Zhang98Data.getRotation(1), "%.4e"),
                FormatArray.toString(Zhang98Data.getTranslation(1), "%.4e"));
        System.out.printf("\nimg5:\nrot=%strans=%s\n", 
                FormatArray.toString(Zhang98Data.getRotation(5), "%.4e"),
                FormatArray.toString(Zhang98Data.getTranslation(5), "%.4e"));
        
        double[][] diffRSameCenter = Rotation.procrustesAlgorithmForRotation(
            Zhang98Data.getRotation(1), Zhang98Data.getRotation(5));
        
        System.out.printf("\ndifference in rot between img1 and img5=\n%s\n", 
           FormatArray.toString(diffRSameCenter, "%.4e"));
        
        double[] out = new double[3];
        extractThetaFromZYX(diffRSameCenter, out);
        
        System.out.printf("\ndifference in rot between img1 and img5 in euler angles=\n");
        for (double a : out) {
            System.out.printf("%.2f ", 180.*a/Math.PI);
        }
        System.out.println();
        double[] diffTrans = MatrixUtil.subtract(Zhang98Data.getTranslation(1), Zhang98Data.getTranslation(5));
        System.out.printf("\ndifference in trans between img1 and img5=\n%s\n",
            FormatArray.toString(diffTrans, "%.3e"));
        System.out.flush();
    }

    /**
     * Test of calculateReconstruction method, of class Reconstruction.
     */
    public void testCalculateReconstructionWithIntrinsicCamera() throws NotConvergedException {
        
        //test data from:
        //http://www.vision.caltech.edu/bouguetj/calib_doc/htmls/example5.html        
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
        // note: the checkerboard sqaures are 30mm in WCS metric
        //
        // note: radial distortion should be corrected:  use on the original coordinates:
        //     x_corrected = x*(1 + k1*r^2 + k2r^4) where r is distance of point from cc.
        
        double[][] k1Intr;
        double[][] k2Intr;
            k1Intr = Camera.createIntrinsicCameraMatrix(533.07, 341.6, 234.3);
            k2Intr = Camera.createIntrinsicCameraMatrix(536.7, 326.5, 249.3);
        
        double[][] k1ExtrRot = MatrixUtil.createIdentityMatrix(3);
        double[] k1ExtrTrans = new double[]{0, 0, 0};
        
        double[][] k2ExtrRot = Rotation.createRodriguesFormulaRotationMatrix(
            new double[]{0.00611, 0.00409, -0.00359});
        double[] k2ExtrTrans = new double[]{-99.85, 0.82, 0.44};
        //double[] k2ExtrTransRev = Arrays.copyOf(k2ExtrTrans, k2ExtrTrans.length);
        //MatrixUtil.multiply(k2ExtrTransRev, -1);
        
        System.out.printf("expected k2ExtrRot from Rodrigues formula\n=%s\n", FormatArray.toString(k2ExtrRot, "%.3e"));
        /*
        [junit] =1.000e+00, 3.577e-03, 4.101e-03 
        [junit] -3.602e-03, 1.000e+00, 6.103e-03 
        [junit] -4.079e-03, -6.117e-03, 1.000e+00 
        */
        
        double[][] x1 = new double[3][10];
        x1[0] = new double[]{129, 160, 140, 226, 225, 232, 341, 407, 532, 527};
        x1[1] = new double[]{145, 319, 361, 391, 61, 284, 289, 122, 48, 302};
        x1[2] = new double[]{1, 1, 1, 1, 1, 1, 1, 1,  1, 1};
        
        double[][] x2 = new double[3][10];
        x2[0] = new double[]{76, 110, 84, 164, 110, 124, 218, 275, 401, 402};
        x2[1] = new double[]{168, 331, 372, 401, 78, 295, 302, 134, 52, 318};
        x2[2] = new double[]{1, 1, 1, 1, 1, 1, 1, 1,  1, 1};
        
        /*
        System.out.printf("results=%s\n\n", rr.toString());
        
  //multiply by focal length?  or K?   see if triangulation in szeliski has advice
        
        System.out.println("XW points normalized by last coordinates:");
        double[] pt = new double[4];
        for (int j = 0; j < rr.XW[0].length; ++j) {    
            for (int i = 0; i < rr.XW.length; ++i) {
                pt[i] = rr.XW[i][j]/rr.XW[3][j];
            }
            System.out.printf("%s\n", FormatArray.toString(pt, "%.3e"));
        }
        */
    }
}
