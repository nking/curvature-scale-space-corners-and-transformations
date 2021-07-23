package algorithms.imageProcessing.transform;

import algorithms.imageProcessing.transform.Camera.CameraExtrinsicParameters;
import algorithms.imageProcessing.transform.Camera.CameraIntrinsicParameters;
import algorithms.imageProcessing.transform.Camera.CameraMatrices;
import algorithms.util.FormatArray;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import junit.framework.TestCase;
import no.uib.cipr.matrix.NotConvergedException;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class PNPTest extends TestCase {
    
    // TODO: test with Zhang dataset
    // TODO: test with Cambridge Landmark dataset
    
    static double eps = 1.e-5;
    
    private static final Level LEVEL = Level.INFO;
    private static final Logger log;
    static {
        log = Logger.getLogger(CameraCalibration.class.getSimpleName());
        log.setLevel(LEVEL);
    }
    
    public PNPTest() {
    }

    /**
     * Test of solveForPose method, of class PNP.
     */
    public void testSolveForPose_ZhangData() throws Exception {
        
        // see testresources/zhang1998/README.txt
        
        // they use f(r) = 1 + k1*r + k2*r^2:
        boolean useR2R4 = false;
        
        int nFeatures = 256;
        int nImages = 5;
        
        // 3 X 256
        double[][] coordsW = Zhang98Data.getFeatureWCS();
        assertEquals(3, coordsW.length);
        assertEquals(nFeatures, coordsW[0].length);
        
        //3 X (256*5)
        double[][] coordsI = Zhang98Data.getFeaturesInAllImages();
        assertEquals(3, coordsI.length);
        assertEquals(nFeatures*nImages, coordsI[0].length);
        
        System.out.printf("coordsW dimensions = [%d X %d]\ncoordsI dimensions = [%d X %d]\n",
            coordsW.length, coordsW[0].length, coordsI.length, coordsI[0].length);
        
        CameraMatrices cameraMatrices = CameraCalibration.estimateCamera(
            nFeatures, coordsI, coordsW, useR2R4);
        
        CameraIntrinsicParameters kIntr = cameraMatrices.getIntrinsics();
        List<Camera.CameraExtrinsicParameters> extrinsics = cameraMatrices.getExtrinsics();
        
        double alpha = kIntr.getIntrinsic()[0][0];
        double gamma = kIntr.getIntrinsic()[0][1];
        double u0 = kIntr.getIntrinsic()[0][2];
        double beta = kIntr.getIntrinsic()[1][1];
        double v0 = kIntr.getIntrinsic()[1][2];
        double[] kRadial = cameraMatrices.getRadialDistortCoeff();
        
        double fX = alpha;
        double fY = beta;
        double oX = u0;
        double oY = v0;
        double skew = gamma;
        
        double alphaE = 871.445;
        double gammaE = 0.2419;
        double u0E = 300.7676;
        double betaE = 871.1251;
        double v0E = 220.8684;
        double k1E = 0.1371;
        double k2E = -0.20101;        
       
        log.log(LEVEL, String.format("\n(fX, fY)=(%.3e, %.3e).  expected=(%.3e, %.3e)\n", fX, fY, alphaE, betaE));
        log.log(LEVEL, String.format("(oX, oY)=(%.3e, %.3e).  expected=(%.3e, %.3e)\n", oX, oY, u0E, v0E));
        log.log(LEVEL, String.format("skew=%.3e.  expected=%.3e\n", skew, gammaE));
        log.log(LEVEL, String.format("[kRadial]=[%.3e, %.3e].  expected=[%.3e, %.3e]\n", 
            kRadial[0], kRadial[1], k1E, k2E));
                
        Camera.CameraExtrinsicParameters ex1;
        for (int i = 0; i < nImages; ++i) {
            ex1 = extrinsics.get(i);
            log.log(LEVEL, String.format("\n"));
            log.log(LEVEL, String.format("   r%d=\n%s\n", i, FormatArray.toString(ex1.getRotation(), "%.3e")));
            log.log(LEVEL, String.format("ansR%d=\n%s\n", i, FormatArray.toString(Zhang98Data.getRotation(i), "%.3e")));
            log.log(LEVEL, String.format("   t%d=\n%s\n", i,FormatArray.toString(ex1.getTranslation(), "%.3e")));
            log.log(LEVEL, String.format("ansT%d=\n%s\n", i,FormatArray.toString(Zhang98Data.getTranslation(i), "%.3e")));
        }
        
        // now have initial parameters to refine using BundleAdjustment.java in other tests
        alphaE = 832.5010;
        gammaE = 0.2046;
        u0E = 303.9584;
        betaE = 832.5309;
        v0E = 206.5879;
        k1E = -0.228601; 
        k2E = 0.190353;
        
        final int nMaxIter = 100;
        
        List<CameraExtrinsicParameters> refinedExtr = PNP.solveForPose(
            coordsI, coordsW, 
            kIntr, cameraMatrices.getExtrinsics(),
            cameraMatrices.getRadialDistortCoeff(), nMaxIter, useR2R4); 
        
        log.log(LEVEL, String.format("\nAfter PNP\n"));
        
        assertEquals(nImages, refinedExtr.size());
                
        for (int i = 0; i < nImages; ++i) {
            log.log(LEVEL, String.format("\n"));
            log.log(LEVEL, String.format("   r%d=\n%s\n", i, 
                    FormatArray.toString(refinedExtr.get(i).getRotation(), "%.3e")));
            log.log(LEVEL, String.format("ansR%d=\n%s\n", i, 
                    FormatArray.toString(Zhang98Data.getRotation(i), "%.3e")));
            log.log(LEVEL, String.format("   t%d=\n%s\n", i,
                    FormatArray.toString(refinedExtr.get(i).getTranslation(), "%.3e")));
            log.log(LEVEL, String.format("ansT%d=\n%s\n", i,
                    FormatArray.toString(Zhang98Data.getTranslation(i), "%.3e")));
        }
    }
    
}
