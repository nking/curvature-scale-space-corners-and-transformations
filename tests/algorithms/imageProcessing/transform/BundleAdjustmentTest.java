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
    
    public void test0(){}
    
    /**
     * Test of solveUsingSparse method, of class BundleAdjustment.
     */
    public void testSolveUsingSparse_ZhangData() throws Exception {
        
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
        
        log.log(LEVEL, String.format("coordsW dimensions = [%d X %d]\ncoordsI dimensions = [%d X %d]\n",
            coordsW.length, coordsW[0].length, coordsI.length, coordsI[0].length));
        
        Camera.CameraMatrices cameraMatrices = CameraCalibration.estimateCamera(nFeatures, coordsI, coordsW, useR2R4);
        
        Camera.CameraIntrinsicParameters kIntr = cameraMatrices.getIntrinsics();
        List<Camera.CameraExtrinsicParameters> extrinsics = cameraMatrices.getExtrinsics();

        // current + delta = expected (zhang)
        double[][] rDeltas = new double[nImages][];
        double[][] tDeltas = new double[nImages][];
        for (int i = 0; i < nImages; ++i) {
            // zhang data doesn't have the Rodrigues vectors, only has the rotation matrices,
            // so this is not a unique delta for rotation:
            double[] rz = Rotation.extractRodriguesRotationVectorBouguet(Zhang98Data.getRotation(i + 1)).om;
            rDeltas[i] = MatrixUtil.subtract(rz, extrinsics.get(i).getRodriguesVector());
            tDeltas[i] = MatrixUtil.subtract(Zhang98Data.getTranslation(i + 1), extrinsics.get(i).getTranslation());
        }
        System.out.printf("delta rVecs=\n%s\n", FormatArray.toString(rDeltas, "%.3e"));
        System.out.printf("delta trans=\n%s\n", FormatArray.toString(tDeltas, "%.3e"));
        
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

        log.log(LEVEL, "From CameraCalibration.estimateCamera:\n");
        log.log(LEVEL, String.format("\n(fX, fY)=(%.3e, %.3e).  expected=(%.3e, %.3e)\n", fX, fY, alphaE, betaE));
        log.log(LEVEL, String.format("(oX, oY)=(%.3e, %.3e).  expected=(%.3e, %.3e)\n", oX, oY, u0E, v0E));
        log.log(LEVEL, String.format("skew=%.3e.  expected=%.3e\n", skew, gammaE));
        log.log(LEVEL, String.format("[kRadial]=[%.3e, %.3e].  expected=[%.3e, %.3e]\n", 
            kRadial[0], kRadial[1], k1E, k2E));

        List<String> distances = new ArrayList<String>();

        double distR, distT;
        Camera.CameraExtrinsicParameters ex1;
        for (int i = 0; i < nImages; ++i) {
            ex1 = extrinsics.get(i);
            distR = MatrixUtil.spectralNorm(MatrixUtil.pointwiseSubtract(Zhang98Data.getRotation(i+1), ex1.getRotation()));
            distT = MatrixUtil.lPSum(
                    MatrixUtil.subtract(Zhang98Data.getTranslation(i+1), ex1.getTranslation()), 2
            );
            log.log(LEVEL, String.format("\n"));
            distances.add(String.format("%d) dist R = %.3e, dist T = %.3e", i, distR, distT));
            log.log(LEVEL, String.format("%d) dist R = %.3e, dist T = %.3e", i, distR, distT));
            log.log(LEVEL, String.format("   rVec%d=\n%s\n", i, FormatArray.toString(ex1.getRodriguesVector(), "%.3e")));
            log.log(LEVEL, String.format("   r%d=\n%s\n", i, FormatArray.toString(ex1.getRotation(), "%.3e")));
            log.log(LEVEL, String.format("ansR%d=\n%s\n", i, FormatArray.toString(Zhang98Data.getRotation(i+1), "%.3e")));
            log.log(LEVEL, String.format("   t%d=\n%s\n", i,FormatArray.toString(ex1.getTranslation(), "%.3e")));
            log.log(LEVEL, String.format("ansT%d=\n%s\n", i,FormatArray.toString(Zhang98Data.getTranslation(i+1), "%.3e")));
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
            ex1 = extrinsics.get(i);
            extrRotVecs[i] = Arrays.copyOf(ex1.getRodriguesVector(), ex1.getRodriguesVector().length);
            extrTrans[i] = Arrays.copyOf(ex1.getTranslation(), 3);
        }
        
        //Zhang98Data.printObservedMinusProjected_Camera_Frame();
        //Zhang98Data.printObservedMinusProjected_Camera_Frame_Eqn2();
        //Zhang98Data.printObservedMinusProjected_Image_Frame();        
        //Zhang98Data.printObservedMinusProjected_Image_Frame_Eqn2();
        
        BundleAdjustment ba = new BundleAdjustment();
        //ba.setUseHomography();
        ba.solveSparsely(coordsI, coordsW, imageFeaturesMap,
            intr, extrRotVecs, extrTrans,
            kRadials, nMaxIter, useR2R4);
        
        log.log(LEVEL, String.format("\nAfter BundleAdjustment\n"));
        
        // only f is refined:
        fX = intr.getBlock(0, 0)[0][0];
        fY = intr.getBlock(0, 0)[1][1];
        oX = intr.getBlock(0, 0)[0][2];
        oY = intr.getBlock(0, 0)[1][1];
        skew = intr.getBlock(0, 0)[0][1];
        kRadial = kRadials[0];
        
        log.log(LEVEL, String.format("\n(fX, fY)=(%.3e, %.3e).  expected=(%.3e, %.3e)\n", 
            fX, fY, alphaE, betaE));
        log.log(LEVEL, String.format("(oX, oY)=(%.3e, %.3e).  expected=(%.3e, %.3e)\n", 
            oX, oY, u0E, v0E));
        log.log(LEVEL, String.format("skew=%.3e.  expected=%.3e\n", 
            skew, gammaE));
        log.log(LEVEL, String.format("[kRadial]=[%.3e, %.3e].  expected=[%.3e, %.3e]\n", 
            kRadial[0], kRadial[1], k1E, k2E));

        for (i = 0; i < nImages; ++i) {
            r = Rotation.createRodriguesRotationMatrixBouguet(extrRotVecs[i]).r;

            distR = MatrixUtil.spectralNorm(MatrixUtil.pointwiseSubtract(Zhang98Data.getRotation(i+1), r));
            distT = MatrixUtil.lPSum(
                    MatrixUtil.subtract(Zhang98Data.getTranslation(i+1), extrTrans[i]), 2
            );

            log.log(LEVEL, String.format("\n"));
            distances.add(String.format("%d) dist R = %.3e, dist T = %.3e", i, distR, distT));
            log.log(LEVEL, String.format("%d) dist R = %.3e, dist T = %.3e", i, distR, distT));
            log.log(LEVEL, String.format("   r%d=\n%s\n", i, 
                FormatArray.toString(r, "%.3e")));
            log.log(LEVEL, String.format("ansR%d=\n%s\n", i, 
                FormatArray.toString(Zhang98Data.getRotation(i+1), "%.3e")));
            log.log(LEVEL, String.format("   t%d=\n%s\n", i, 
                FormatArray.toString(extrTrans[i], "%.3e")));
            log.log(LEVEL, String.format("ansT%d=\n%s\n", i,
                FormatArray.toString(Zhang98Data.getTranslation(i+1), "%.3e")));
        }

        for (String str : distances) {
            log.log(LEVEL, str);
        }
        
    }

}
