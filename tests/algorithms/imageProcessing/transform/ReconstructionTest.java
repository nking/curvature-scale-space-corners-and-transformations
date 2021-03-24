package algorithms.imageProcessing.transform;

import algorithms.matrix.MatrixUtil;
import algorithms.util.FormatArray;
import java.util.Arrays;
import junit.framework.TestCase;
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
     * Test of calculateReconstruction method, of class Reconstruction.
     */
    public void testCalculateReconstructionWithIntrinsicCamera() throws NotConvergedException {
        
        //test data from:
        //http://www.vision.caltech.edu/bouguetj/calib_doc/htmls/example5.html        
        // "Fifth calibration example - Calibrating a stereo system, stereo image rectification and 3D stereo triangulation"
        // by Jean-Yves Bouguet
        // Camera Calibration Toolbox for Matlab
        
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
        x1[0] = new double[]{307, 345, 478, 476, 278, 247, 407, 408,  450, 233};
        x1[1] = new double[]{159, 188, 87, 265, 256, 191, 123, 228,    90, 283};
        x1[2] = new double[]{1, 1, 1, 1, 1, 1, 1, 1,  1, 1};
        
        double[][] x2 = new double[3][8];
        x2[0] = new double[]{184, 215, 344, 347, 161, 131, 275, 278,  402, 123};
        x2[1] = new double[]{172, 238, 95, 278, 268, 205, 134, 241,  317, 296};
        x2[2] = new double[]{1, 1, 1, 1, 1, 1, 1, 1,  1, 1};
        
        Reconstruction.ReconstructionResults rr = Reconstruction
            .calculateReconstructionWithIntrinsic(k1Intr, k2Intr, x1, x2);
        
        System.out.printf("results=%s\n\n", rr.toString());
        
        System.out.println("XW points normalized by last coordinates:");
        double[] pt = new double[4];
        for (int j = 0; j < rr.XW[0].length; ++j) {    
            for (int i = 0; i < rr.XW.length; ++i) {
                pt[i] = rr.XW[i][j]/rr.XW[3][j];
            }
            System.out.printf("%s\n", FormatArray.toString(pt, "%.3e"));
        }
        
    }
}
