package algorithms.imageProcessing.transform;

import algorithms.matrix.BlockMatrixIsometric;
import algorithms.matrix.MatrixUtil;
import algorithms.util.FormatArray;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import junit.framework.TestCase;
import static junit.framework.TestCase.assertEquals;

/**
 *
 * @author nichole
 */
public class BundleAdjustmentTest extends TestCase {
    
    static double eps = 1.e-5;
    
    private static final Level LEVEL = Level.INFO;
    private static final Logger log;
    static {
        log = Logger.getLogger(CameraCalibration.class.getSimpleName());
        log.setLevel(LEVEL);
    }
    
    public BundleAdjustmentTest() {
    }

    /**
     * Test of solveUsingSparse method, of class BundleAdjustment.
     */
    public void testSolveUsingSparse_ZhangData() throws Exception {

        // this shows that the refinement methods need improvements
        
        // see testresources/zhang1998/README.txt
        /*
        each image is 640×480  [pixels^2]
        The pixel is square (aspect ratio = 1); 
        the focal length = 832.5 pixels; 
        the image center is at (303.959, 206.585); 
        there is a significant radial distortion: k1 = -0.228601, k2 = 0.190353.
        */
        
        // they use f(r) = 1 + k1*r + k2*r^2:
        boolean useR2R4 = true;
        
        int nFeatures = 256;
        int nImages = 5;
        
        // 3 X 256
        double[][] coordsW = Zhang98Data.getFeatureWCS();
        assertEquals(3, coordsW.length);
        assertEquals(nFeatures, coordsW[0].length);
        // mean of x coords is 3.361.  stdev=2.056
        // mean of y coords is -3.361.  stdev=2.056
        
        //double[] meanT = MatrixUtil.mean(MatrixUtil.transpose(coordsW));
        //double[] stdevT = MatrixUtil.standardDeviation(MatrixUtil.transpose(coordsW));
        
        //3 X (256*5)
        double[][] coordsI = Zhang98Data.getObservedFeaturesInAllImages();
        assertEquals(3, coordsI.length);
        assertEquals(nFeatures*nImages, coordsI[0].length);
        
        log.fine(String.format("coordsW dimensions = [%d X %d]\ncoordsI dimensions = [%d X %d]\n",
            coordsW.length, coordsW[0].length, coordsI.length, coordsI[0].length));

        boolean useBouguetForRodrigues = false;

        Camera.CameraMatrices cameraMatrices = CameraCalibration.estimateCameraPlanar(nFeatures, coordsI, coordsW, useR2R4);
        
        Camera.CameraIntrinsicParameters kIntr = cameraMatrices.getIntrinsics();
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
        double gammaE = 0.2419; //<== not significantly different from 0 according to Zhang stats
        double u0E = 300.7676;
        double betaE = 871.1251;
        double v0E = 220.8684;
        double k1E = 0.1371;
        double k2E = -0.20101;

        /*log.fine("From CameraCalibration.estimateCamera:\n");
        log.fine(String.format("\n(fX, fY)=(%.3e, %.3e).  expected=(%.3e, %.3e)\n", fX, fY, alphaE, betaE));
        log.fine(String.format("(oX, oY)=(%.3e, %.3e).  expected=(%.3e, %.3e)\n", oX, oY, u0E, v0E));
        log.fine(String.format("skew=%.3e.  expected=%.3e\n", skew, gammaE));
        log.fine(String.format("[kRadial]=[%.3e, %.3e].  expected=[%.3e, %.3e]\n", 
            kRadial[0], kRadial[1], k1E, k2E));
         */

        List<String> distances = new ArrayList<String>();

        // now have initial parameters to refine using BundleAdjustment.java in other tests
        alphaE = 832.5010;
        gammaE = 0.2046;
        u0E = 303.9584;
        betaE = 832.5309;
        v0E = 206.5879;
        k1E = -0.228601; 
        k2E = 0.190353;
        
        final int nMaxIter = 100;
        int i, j;
        double[][] r;
        
        // each image has nFeatures
        TIntObjectMap<TIntSet> imageFeaturesMap = new TIntObjectHashMap<TIntSet>();
        for (i = 0; i < nImages; ++i) {
            TIntSet set = new TIntHashSet();
            for (j = 0; j < nFeatures; ++j) {
                set.add(j);
            }
            imageFeaturesMap.put(i, set);
        }

        // for each camera, add an intrinsic camera matrix.  [(3*nImages) X 3]
        // for this dataset, there is only one camera
        BlockMatrixIsometric intr = new BlockMatrixIsometric(MatrixUtil.zeros(nImages*3, 3), 3, 3);
        // for each camera, a row of radial distortion coefficients [nImages X 2]
        double[][] kRadials = new double[nImages][];
        for (i = 0; i < nImages; ++i) {
            intr.setBlock(kIntr.getIntrinsic(), i, 0);
            kRadials[i] = new double[]{kRadial[0], kRadial[1]};
            //kRadials[i] = new double[2];
        }

        //the extrinsic camera parameter rotation vectors
        //stacked along the 3 columns, that is the size is [nImages X 3] where
        //nImages is coordsI[0].length/coordsW[0].length.  each array is size [1X3]
        double[][] extrRotVecs = new double[nImages][];
        double[][] extrTrans = new double[nImages][];
        for (i = 0; i < nImages; ++i) {
            Camera.CameraExtrinsicParameters init = extrinsics.get(i);
            extrRotVecs[i] = Arrays.copyOf(init.getRodriguesVector(), init.getRodriguesVector().length);
            extrTrans[i] = Arrays.copyOf(init.getTranslation(), 3);
        }
        
        //Zhang98Data.printObservedMinusProjected_Camera_Frame();
        //Zhang98Data.printObservedMinusProjected_Camera_Frame_Eqn2();
        //Zhang98Data.printObservedMinusProjected_Image_Frame();        
        //Zhang98Data.printObservedMinusProjected_Image_Frame_Eqn2();

        BundleAdjustment ba = new BundleAdjustment();
        ba.setUseHomography();
        ba.solveSparsely(coordsI, coordsW, imageFeaturesMap,
            intr, extrRotVecs, extrTrans,
            kRadials, nMaxIter, useR2R4, useBouguetForRodrigues);
        
        log.fine(String.format("\nAfter BundleAdjustment\n"));
        
        // only f is refined:
        double fX2 = intr.getBlock(0, 0)[0][0];
        double fY2 = intr.getBlock(0, 0)[1][1];
        double oX2 = intr.getBlock(0, 0)[0][2];
        double oY2 = intr.getBlock(0, 0)[1][1];
        double skew2 = intr.getBlock(0, 0)[0][1];
        double[] kRadial2 = kRadials[0];

        System.out.printf("Compare method results to Zhang 1998:\n");

        log.fine(String.format("\n(fx0, fx1, fx expected)=(%.3e, %.3e, %.3e)\n",
            fX, fX2, alphaE));
        log.fine(String.format("\n(fy0, fy1, fy expected)=(%.3e, %.3e, %.3e)\n",
                fY, fY2, betaE));
        log.fine(String.format("\n(ox0, ox1, ox expected)=(%.3e, %.3e, %.3e)\n",
                oX, oX2, u0E));
        log.fine(String.format("\n(oy0, oy1, oy expected)=(%.3e, %.3e, %.3e)\n",
                oY, oY2, v0E));
        log.fine(String.format("\n(skew0, skew1, skew expected)=(%.3e, %.3e, %.3e)\n",
                skew, skew2, gammaE));
        log.fine(String.format("kRadial0=[%.3e, %.3e], kRadial1=[%.3e, %.3e],  kRadial expected=[%.3e, %.3e]\n",
            kRadial[0], kRadial[1], kRadial2[0], kRadial2[1], k1E, k2E));

        boolean passive = true;

        Camera.CameraExtrinsicParameters[] refined = new Camera.CameraExtrinsicParameters[nImages];
        for (i = 0; i < nImages; ++i) {
            Camera.CameraExtrinsicParameters init = extrinsics.get(i);
            double[][] xi = Zhang98Data.getObservedFeaturesInImage(i + 1);
            //this needs image coordinates because internally it is projecting X to camera then image and comparing that to x

            refined[i] = CameraPose.bouguetPoseRefine(init, kIntr, xi, coordsW, useBouguetForRodrigues);

            // print camera calib, bouguet refined rot(rVec), then zhang98
            //double[] cCalibRVec = init.getRodriguesVector();
            //double[] bRefinedRVec = refined[i].getRodriguesVector();
            double[][] cCalibRot = init.getRotation();
            double[][] bouguetRefinedRot = refined[i].getRotation();
            double[][] expectedRot = Zhang98Data.getRotation(i + 1);
            double[][] bARefinedRot = Rotation.createRotationRodriguesFormula(extrRotVecs[i], passive);

            System.out.printf("%d) rotation:\n  rot0=\n%s\n  rot refined (BA)=\n%s\n   rot refined (bouguet)=\n%s\n  rot expected=\n%s\n",
                    i,
                    FormatArray.toString(cCalibRot, "%.3e"),
                    FormatArray.toString(bARefinedRot, "%.3e"),
                    FormatArray.toString(bouguetRefinedRot, "%.3e"),
                    FormatArray.toString(expectedRot, "%.3e")
                    );

            // translations
            double[] cCalibT = init.getTranslation();
            double[] bouguetRefinedT = refined[i].getTranslation();
            double[] bARefinedT = extrTrans[i];
            double[] expectedT = Zhang98Data.getTranslation(i + 1);

            System.out.printf("\n%d) translation:\n  T0=%s\n  T refined (BA)=%s\n   T refined (Bouguet)=%s\n   T expected=%s\n",
                    i,
                    FormatArray.toString(cCalibT, "%.3e"),
                    FormatArray.toString(bARefinedT, "%.3e"),
                    FormatArray.toString(bouguetRefinedT, "%.3e"),
                    FormatArray.toString(expectedT, "%.3e"));

        }

    }

}
