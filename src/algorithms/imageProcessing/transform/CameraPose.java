package algorithms.imageProcessing.transform;

import algorithms.imageProcessing.features.RANSACSolver;
import algorithms.imageProcessing.matching.ErrorType;
import algorithms.imageProcessing.transform.Camera.CameraExtrinsicParameters;
import algorithms.imageProcessing.transform.Camera.CameraIntrinsicParameters;
import algorithms.imageProcessing.transform.Camera.CameraParameters;
import algorithms.matrix.MatrixUtil;
import algorithms.util.FormatArray;
import java.util.Arrays;
import java.util.List;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.NotConvergedException;
import no.uib.cipr.matrix.RQ;

/**
 given a set of features in image coordinates and world coordinate space with
  known camera intrinsic parameters, estimate the camera pose, that is
  extract the camera extrinsic parameters.
 
 * TODO: consider solving with M-estimators.
 * see http://research.microsoft.com/en- us/um/people/zhang/INRIA/Publis/Tutorial-Estim/node24.html
 * 
 * @author nichole
 */
public class CameraPose {
    
    public static double eps = 1e-7;
    /**
     * given correspondence between two images in image coordinates calculate 
     * the extrinsic camera parameters.
     * 
     * <pre>
     * following CMU lectures of Kris Kitani at 
     * http://www.cs.cmu.edu/~16385/s17/Slides/12.5_Reconstruction.pdf
     * 
     * add other references:
     * </pre>
     * @param k1 intrinsic camera matrix for image 1 in units of pixels.
     * @param k2 intrinsic camera matrix for image 2 in units of pixels.
     * @param x1 the image 1 set of correspondence points.  format is 3 x N where
     * N is the number of points.
     * @param x2 the image 2 set of correspondence points.  format is 3 x N where
     * N is the number of points.
     * @return 
     * @throws no.uib.cipr.matrix.NotConvergedException 
     */
    public static Camera.CameraExtrinsicParameters[] calculateUsingEssentialMatrix(
        double[][] k1, double[][] k2,
        double[][] x1, double[][] x2) throws NotConvergedException {
        
        double[][] outputXW = MatrixUtil.zeros(4, x1[0].length);
        return calculateUsingEssentialMatrix(k1, k2, x1, x2, outputXW);
    }
    
     /**
     * given correspondence between two images in image coordinates calculate 
     * the extrinsic camera parameters.
     * 
     * <pre>
     * following CMU lectures of Kris Kitani at 
     * http://www.cs.cmu.edu/~16385/s17/Slides/12.5_Reconstruction.pdf
     * 
     * add:
     * Sect 7.2 of Szeliski 2010
     * </pre>
     * @param k1 intrinsic camera matrix for image 1 in units of pixels.
     * @param k2 intrinsic camera matrix for image 2 in units of pixels.
     * @param x1 the image 1 set of correspondence points in image reference frame.  
     * format is 3 x N where N is the number of points.
     * @param x2 the image 2 set of correspondence points in image reference frame.  
     * format is 3 x N where N is the number of points.
     * @param outputXW
     * @return 
     * @throws no.uib.cipr.matrix.NotConvergedException 
     */
    public static CameraExtrinsicParameters[] calculateUsingEssentialMatrix(
        double[][] k1, double[][] k2,
        double[][] x1, double[][] x2, double[][] outputXW) throws NotConvergedException {
                
        if (x1.length != 3 || x2.length != 3) {
            throw new IllegalArgumentException("x1.length must be 3 and so must x2.length");
        }
        int n = x1[0].length;
        if (x2[0].length != n) {
            throw new IllegalArgumentException("x1 and x2 must be same dimensions");
        }
        
        /*
        http://www.cs.cmu.edu/~16385/s17/Slides/12.5_Reconstruction.pdf
        
        (1) compute fundamental mat5rix FM from the correspondence x1, x2
        (2) compute the camera matrices P1, P2 from FM.
        (3) For each point correspondence, compute the point X in 3D space (triangularization)
        */
        
        double[][] k1IntrInv = Camera.createIntrinsicCameraMatrixInverse(k1);
        double[][] k2IntrInv = Camera.createIntrinsicCameraMatrixInverse(k2);
        
        // the direction of the points is calculated by K^-1 * x
        double[][] x1Direction = MatrixUtil.multiply(k1IntrInv, x1);
        double[][] x2Direction = MatrixUtil.multiply(k2IntrInv, x2);
                
        DenseMatrix x1M = new DenseMatrix(x1Direction);
        DenseMatrix x2M = new DenseMatrix(x2Direction);
        
        // normalizing by unit standard to improve results of epipolar solution:
        EpipolarTransformer.NormalizedXY normXY1 
            = EpipolarTransformer.normalizeUsingUnitStandard(x1M);
        EpipolarTransformer.NormalizedXY normXY2 
            = EpipolarTransformer.normalizeUsingUnitStandard(x2M);
        
        DenseMatrix leftM = normXY1.getXy();
        DenseMatrix rightM = normXY2.getXy();
        
        double tolerance = 3.84; //3.84 5.99 7.82        
        boolean useToleranceAsStatFactor = true;
        ErrorType errorType = ErrorType.SAMPSONS;
        EpipolarTransformationFit fitR = null;
        boolean reCalcIterations = false;
        
        DenseMatrix normalizedE;
        
        /*EpipolarTransformer tr = new EpipolarTransformer();        
        normalizedE = tr.calculateEpipolarProjection(leftM, rightM);
        DenseMatrix vNFM = tr.validateSolution(normalizedFM, leftM, rightM);
        
        Distances distances = new Distances();
        if (useToleranceAsStatFactor) {
            fitR = distances.calculateError2(vNFM, leftM, rightM,
                    errorType, tolerance);
        } else {
            fitR = distances.calculateError(vNFM, leftM, rightM,
                    errorType, tolerance);
        }*/
        
        //if true, solves for the Essential Matrix, else solves
        //for the Fundamental Matrix.  The difference is in the diagonal used for
        //dimension reduction.
        boolean coordsAreInCameraRefFrame = true;
        
        RANSACSolver solver = new RANSACSolver();
        fitR = solver.calculateEpipolarProjection(
            leftM, rightM, errorType, useToleranceAsStatFactor, tolerance,
                reCalcIterations, coordsAreInCameraRefFrame);
        
        System.out.println("RANSAC fit=" + fitR.toString());
        
        normalizedE = fitR.getFundamentalMatrix();
        
        // this is now back in the reference frame of the x1Direction and x2Direction
        DenseMatrix essentialM = EpipolarTransformer.denormalizeTheFundamentalMatrix(
            normalizedE, normXY1.getNormalizationMatrices(),
            normXY2.getNormalizationMatrices());
        
        //=========
        double[][] _essentialMatrix = MatrixUtil.convertToRowMajor(essentialM);
        
        MatrixUtil.SVDProducts svdE = MatrixUtil.performSVD(_essentialMatrix);
        
        assert(svdE.u[0].length == 3 && svdE.u.length == 3);
        
        double detU = MatrixUtil.determinant(svdE.u);
        double detV = MatrixUtil.determinant(svdE.vT);
        
        System.out.printf("SVD.u=\n%s\n", FormatArray.toString(svdE.u, "%.3e"));
        System.out.printf("SVD.s=\n%s\n", FormatArray.toString(svdE.s, "%.3e"));
        System.out.printf("SVD.vT=\n%s\n", FormatArray.toString(svdE.vT, "%.3e"));
        System.out.printf("det(SVD.u)=%.2f\n", detU);
        System.out.printf("det(SVD.vT)=%.2f\n", detV);
        
        /*
        Szeliski 2010 chap 7:
        Once an estimate for the essential matrix E has been recovered, 
        the direction of the translation vector t can be estimated. 
        
        Note that the absolute distance between the two cameras can never 
        be recovered from pure image measurements alone without knowledge 
        about absolute camera and point positions or distances, often called ground 
        control points in photogrammetry.
        */
        
        // det(R)=1 is a proper rotation matrix.  rotation angles are counterclockwise.
        //           it's a special orthogonal matrix and provides the
        //           defining matrix representation of the group of proper n-dimensional rotations, denoted
        //           by SO(n). http://scipp.ucsc.edu/~haber/ph251/rotreflect_17.pdf
        // det(R)=-1 is an improper rotation matrix representing rotations that
        //           require mirrors.
        //           The most general improper rotation matrix is a product of a proper rotation by an
        //           angle θ about some axis nˆ and a mirror reflection through a plane that passes through
        //           the origin and is perpendicular to nˆ.  NOTE: nˆ is determined by
        //           the right hand rule.
        
        /*
        Sect 7.2 of Szeliski 2010 eqn (7.25) introduces
        R3 and R4 constructed from -U as 2 more rotation possibilities to be tested
        and that is necessary in some cases where det(R) would otherwise be -1
        (reflection).
        ...we only know both E and tˆup to a sign. Furthermore, the matrices U and V
        are not guaranteed to be rotations (you can flip both their signs and 
        still get a valid SVD).   
        For this reason, we have to generate all four possible rotation matrices
        
        R = +-U * R_Z(+-90)^T * V^T
           and keep the 2 whose determinant = 1
        */
        double[][] r1 = MatrixUtil.zeros(3, 3);
        double[][] r2 = MatrixUtil.zeros(3, 3);
        double[][] u = MatrixUtil.zeros(3, 3);
        populateWithDet1Rs(svdE, r1, r2, u);
        
        // last column in u is the second epipole and is the direction of vector t
        double[] t1 = MatrixUtil.extractColumn(u, 2);
        double[] t2 = Arrays.copyOf(t1, t1.length);
        MatrixUtil.multiply(t2, -1); 
        
        System.out.printf("R1=\n%s\n", FormatArray.toString(r1, "%.3e"));
        System.out.printf("R2=\n%s\n", FormatArray.toString(r2, "%.3e"));
        System.out.printf("t1=\n%s\n", FormatArray.toString(t1, "%.3e"));
        System.out.printf("t2=\n%s\n", FormatArray.toString(t2, "%.3e"));
        
        //then of the 4 possible choices find the one with largest number of positive Z.
          
        //NOTE: the last column vector in u is the smallest
        //    eigenvector.  it is epipole2, that is, the right image position 
        //    of the epipolar projection of the left camera center.
        //    it's int the left null space of E.
        
        //or use directions x1Direction, x2Direction
                
        x1M = extractIndices(new DenseMatrix(x1), fitR.inlierIndexes);
        x2M = extractIndices(new DenseMatrix(x2), fitR.inlierIndexes);
        x1 = MatrixUtil.convertToRowMajor(x1M);
        x2 = MatrixUtil.convertToRowMajor(x2M);
        
        // solution 1:  Rot1 and T1
        // solution 2:  Rot1 and T2
        // solution 3:  Rot2 and T2
        // solution 4:  Rot2 and T1
        double[][] rSelected = MatrixUtil.zeros(3, 3);
        double[] tSelected = new double[3];
        double[][] XW = MatrixUtil.zeros(4, x1[0].length);
        chooseRAndT(x1, x2, k1, k2, r1, r2, t1, t2, rSelected, tSelected, XW);  
        
        Camera.CameraExtrinsicParameters c1 = new Camera.CameraExtrinsicParameters();
        c1.setRotation(MatrixUtil.createIdentityMatrix(3));
        c1.setTranslation(new double[]{0, 0, 0});
          
        Camera.CameraExtrinsicParameters c2 = new Camera.CameraExtrinsicParameters();
        c2.setRotation(rSelected);
        c2.setTranslation(tSelected);
        
        return new Camera.CameraExtrinsicParameters[]{c1, c2};
    }

    /**
     * given a set of features in image space and world coordinate space with
     * known camera intrinsic parameters, estimate the camera pose, that is
     * extract the camera extrinsic parameters.
     * calibrating the camera extrinsic parameters is a.k.a. 
     * perspective-n-point-problem where n is the number of features (a.k.a. points).
     * This method uses DLT and could be followed by non-linear optimization
     * to improve the parameter estimates.
     * <pre>
     * references:
     *  Kitani lecture notes http://www.cs.cmu.edu/~16385/s17/Slides/11.3_Pose_Estimation.pdf
     * </pre>
     * @param x the image coordinates of the features in format 3 X N where
     * 3 is for x, y, 1 rows, and N columns is the number of features.  At least 3 features are needed to 
     * calculate the extrinsic parameters.
     * NOTE x and X should both be distortion-free or both should be distorted.
     * @param X the world coordinates of the features in format 3 X N where
     * 3 is for x, y, 1 rows, and N columns is the number of features.  At least 3features are needed to 
     * calculate the extrinsic parameters.
     * NOTE x and X should both be distortion-free or both should be distorted.
     * @return 
     */
    public static CameraParameters calculatePoseUsingDLT(double[][] x, double[][] X) 
        throws NotConvergedException {
                
        if (x.length != 3) {
            throw new IllegalArgumentException("x.length must be 3");
        }
        if (X.length != 3) {
            throw new IllegalArgumentException("X.length must be 3");
        }
        int n = x[0].length;
        if (n < 6) {
            throw new IllegalArgumentException("x must have at least 6 correspondences");
        }
        if (X[0].length != n) {
            throw new IllegalArgumentException("the number of columns in X must be the same as in x");
        }
        
        // normalize by last coordinate:
        /*for (int i = 0; i < x[0].length; ++i) {
            x[0][i] /= x[2][i];
            x[1][i] /= x[2][i];
        }*/
                
        // 2*n X 12       
        double xi, yi, Xi, Yi, Zi;
        double[][] ell = new double[2*n][12];
        for (int i = 0; i < n; ++i) {
            xi = x[0][i];
            yi = x[1][i];
            Xi = X[0][i];
            Yi = X[1][i];
            Zi = X[2][i];
            // http://www.cs.cmu.edu/~16385/s17/Slides/11.3_Pose_Estimation.pdf
            ell[2*i]     = new double[]{Xi, Yi, Zi, 1, 0, 0, 0, 0, -xi*Xi, -xi*Yi, -xi*Zi, -xi};
            ell[2*i + 1] = new double[]{0, 0, 0, 0, Xi, Yi, Zi, 1, -yi*Xi, -yi*Yi, -yi*Zi, -yi};
        }
        
        MatrixUtil.SVDProducts svd = MatrixUtil.performSVD(ell);
        
        // vT is 12X12.  last row in vT is the eigenvector for the smallest eigenvalue
        double[] xOrth = svd.vT[svd.vT.length - 1];
        
        // reshape into 3 X 4
        double[][] P2 = MatrixUtil.zeros(3, 4);
        System.arraycopy(xOrth, 0, P2[0], 0, 4);
        System.arraycopy(xOrth, 4, P2[1], 0, 4);
        System.arraycopy(xOrth, 8, P2[2], 0, 4);
        
        MatrixUtil.SVDProducts svdP2 = MatrixUtil.performSVD(P2);
        double[] c = MatrixUtil.extractColumn(svdP2.u, 2);
        
        // assert P2*c = 0
        double[] check0 = MatrixUtil.multiplyMatrixByColumnVector(P2, c);
        System.out.printf("check that P2*c=0:%s\n", FormatArray.toString(check0, "%.3e"));
        
        double[][] M = MatrixUtil.copySubMatrix(P2, 0, 2, 0, 2);
        RQ rq = RQ.factorize(new DenseMatrix(M));
        
        System.out.printf("RQ.R=\n%s\n", FormatArray.toString(
            MatrixUtil.convertToRowMajor(rq.getR()), "%.3e"));
        System.out.printf("RQ.Q=\n%s\n", FormatArray.toString(
            MatrixUtil.convertToRowMajor(rq.getQ()), "%.3e"));
        
        double[][] kIntr = MatrixUtil.convertToRowMajor(rq.getR());
        MatrixUtil.multiply(kIntr, 1./kIntr[2][2]);
            
        double[][] kExtrRot = MatrixUtil.convertToRowMajor(rq.getQ());
        System.out.printf("  decomposed into intrinsic=\n   %s\n", FormatArray.toString(kIntr, "%.3e"));
        System.out.printf("  decomposed into extrinsic rotation=\n   %s\n", FormatArray.toString(kExtrRot, "%.3e"));
            
        System.out.printf("  decomposed into extrinsic translation=\n   %s\n", FormatArray.toString(c, "%.3e"));
        
        CameraExtrinsicParameters extrinsics = new CameraExtrinsicParameters();
        extrinsics.setRotation(kExtrRot);
        extrinsics.setTranslation(c);
        
        CameraIntrinsicParameters intrinsics = new CameraIntrinsicParameters();
        intrinsics.setIntrinsic(kIntr);
        CameraParameters camera = new CameraParameters(intrinsics, extrinsics);
        
        return camera;
    }
    
    /**
     * given a set of features in image space and world coordinate space with
     * known camera intrinsic parameters, estimate the camera pose, that is
     * extract the camera extrinsic parameters.
     * calibrating the camera extrinsic parameters is a.k.a. 
     * perspective-n-point-problem where n is the number of features (a.k.a. points).
     * This method uses DLT and could be followed by non-linear optimization
     * to improve the parameter estimates.
     * <pre>
     * references:
     * Ma, Chen, & Moore 2003 "Camera Calibration: a USU Implementation"
     * http://www.cs.cmu.edu/~16385/s17/Slides/11.3_Pose_Estimation.pdf
     * Zhang 1999, "Flexible Camera Calibration By Viewing a Plane From Unknown Orientations"
     * Szeliski 2010 draft of "Computer Vision: Algorithms and Applications"
     * </pre>
     * @param intrinsics
     * @param x the image coordinates of the features in format 3 X N where
     * 3 is for x, y, 1 rows, and N columns is the number of features.  At least 3 features are needed to 
     * calculate the extrinsic parameters.
     * NOTE x and X should both be distortion-free or both should be distorted.
     * @param X the world coordinates of the features in format 3 X N where
     * 3 is for x, y, 1 rows, and N columns is the number of features.  At least 3features are needed to 
     * calculate the extrinsic parameters.
     * NOTE x and X should both be distortion-free or both should be distorted.
     * @return 
     */
    public static CameraExtrinsicParameters calculatePoseUsingCameraCalibration(
        Camera.CameraIntrinsicParameters intrinsics, double[][] x,
        double[][] X) throws NotConvergedException {
                
        if (x.length != 3) {
            throw new IllegalArgumentException("x.length must be 3");
        }
        if (X.length != 3) {
            throw new IllegalArgumentException("X.length must be 3");
        }
        int n = x[0].length;
        if (n < 3) {
            throw new IllegalArgumentException("x must have at least 3 correspondences");
        }
        if (X[0].length != n) {
            throw new IllegalArgumentException("the number of columns in X must be the same as in x");
        }
        
        int nImages = 1;
        
        // following Ma et al. 2003
        double[][] h = CameraCalibration.solveForHomography(x, X);
        
        CameraExtrinsicParameters kExtr = CameraCalibration.solveForExtrinsic(
            intrinsics, h);
        
        return kExtr;
    }
    
    /**
     * NOT YET IMPLEMENTED.
     * given n 3D-to-2D point correspondences, estimates the pose 
     * of a calibrated camera (a.k.a. P-n-P) with computational complexity O(n)
     * using the Moreno-Noguer et al. 2007 non-iterative algorithm.
     * This could be followed by non-linear optimization
     * to improve the parameter estimates.
     * <pre>
     * references:
     * Moreno-Noguer, Lepetite, & Fua 2007, "Accurate Non-Iterative O(n) Solution to the PnP Problem"
     * Szeliski 2010 draft of "Computer Vision: Algorithms and Applications"
     * </pre>
     * @param intrinsics
     * @param x the image coordinates of the features in format 3 X N where
     * 3 is for x, y, 1 rows, and N columns is the number of features.  At least 3 features are needed to 
     * calculate the extrinsic parameters.
     * NOTE x and X should both be distortion-free or both should be distorted.
     * @param X the world coordinates of the features in format 3 X N where
     * 3 is for x, y, 1 rows, and N columns is the number of features.  At least 3features are needed to 
     * calculate the extrinsic parameters.
     * NOTE x and X should both be distortion-free or both should be distorted.
     * @return 
     */
    public static CameraExtrinsicParameters calculatePoseUsingPNP(
        Camera.CameraIntrinsicParameters intrinsics, double[][] x,
        double[][] X) throws NotConvergedException {
                
        if (x.length != 3) {
            throw new IllegalArgumentException("x.length must be 3");
        }
        if (X.length != 3) {
            throw new IllegalArgumentException("X.length must be 3");
        }
        int n = x[0].length;
        if (n < 4) {
            throw new IllegalArgumentException("x must have at least 4 correspondences");
        }
        if (X[0].length != n) {
            throw new IllegalArgumentException("the number of columns in X must be the same as in x");
        }
        
        // Szeliski 2010 refers to perspective-n-point-problem (PnP) references  
        //   (Haralick, Lee, Ottenberg et al. 1994; Quan and Lan 1999; Moreno-Noguer, Lepetit, and Fua 2007)
        
        //port the c++ impl of  Moreno-Noguer, Lepetit, and Fua (2007)  here?
        //https://github.com/cvlab-epfl/EPnP/tree/master/cpp   
        
        throw new UnsupportedOperationException("not yet implemented");
    }
    
     /**
     * among the 4 rotation and translation combinations from R1, R1, T1, and T2, 
     * select the one with the largest number of projected Z coordinates which are
     * positive, that is, in front of both cameras.
     * NOTE that inaccuracies in this chirality are larger for points further 
     * away from the cameras and closer to the plane at infinity.
     * NOTE that the determinants of R1 and R2 should have already been checked to be +1.
     * @param x1 image 1 portion of the correspondence pairs.
     * @param x2 image 2 portion of the correspondence pairs.
     * @param k1 intrinsic camera matrix for camera 1
     * @param k2 intrinsic camera matrix for camera 2
     * @param R1 rotation matrix whose determinant is +1
     * @param R2 rotation matrix whose determinant is +1
     * @param t1 translation vector (the direction between camera centers)
     * @param t2 translation vector (the direction between camera centers)
     * @param rSelected output variable holding the R1 or R2, whichever was the 
     * first found as a valid solution.
     * @param tSelected output variable holding the t1 or t2, whichever was the 
     * first found as a valid solution.
     * @param outputX the real world coordinates of the projection of x1 and x2 using
     * triangulation. else null if no valid solution was found
     */
    private static void chooseRAndT(double[][] x1, double[][] x2, 
        double[][] k1, double[][] k2,
        double[][] R1, double[][] R2, double[] t1, double[] t2, 
        double[][] rSelected, double[] tSelected, double[][] outputX) {
    
        int n = x1[0].length;
        
        if (outputX.length != 4) {
            throw new IllegalArgumentException("outputX.length must be 4");
        }
        if (outputX[0].length != n) {
            throw new IllegalArgumentException("outputX[0].length must be the same as x1[0].length");
        }
        
        // for this model, for the first image, the camera extrinsics are
        //    R = I and t = [0], which leaves all rotation and translation in
        //    the 2nd camera extrinsics w.r.t. the first.
        double[][] k1ExtrRot = MatrixUtil.createIdentityMatrix(3);
        double[] k1ExtrTrans = new double[3];
        
        // save the first that pass the tests for Z>=0.
        double[][] bestR = null;
        double[] bestT = null;
        double[][] bestXW = null;
        String bestLabel = null;
        int bestNPosZ = Integer.MIN_VALUE;
        
        double[][] XW;
        double[] XWPt;
        String label = null;
        
        XWPt = new double[4];
        XW = new double[4][n];
        for (int i = 0; i < 4; ++i) {
            XW[i] = new double[n];
        }
            
        double[][] rTst = null;
        double[] tTst = null;
        double[][] x1Pt = new double[3][1];
        double[][] x2Pt = new double[3][1];
        int i, j, ii;
        for (i = 0; i < 3; ++i) {
            x1Pt[i] = new double[1];
            x2Pt[i] = new double[1];
        }
        
        int nPosZ; 
        
        for (j = 0; j < 4; ++j) {
            switch(j) {
                case 0: {
                    label = "R1, T1";
                    rTst = R1;
                    tTst = t1;
                    break;
                }
                case 1: {
                    label = "R1, T2";
                    rTst = R1;
                    tTst = t2;
                    break;
                }
                case 2: {
                    label = "R2, T1";
                    rTst = R2;
                    tTst = t1;
                    break;
                }
                default: {                    
                    label = "R2, T2";
                    rTst = R2;
                    tTst = t2;
                    break;
                }
            }
            nPosZ = 0;
            for (i = 0; i < n; ++i) {
                for (ii = 0; ii < 3; ++ii) {
                    x1Pt[ii][0] = x1[ii][i];
                    x2Pt[ii][0] = x2[ii][i];
                }
                //
                XWPt = Triangulation.calculateWCSPoint(
                    k1, k1ExtrRot, k1ExtrTrans, 
                    k2, rTst, tTst, 
                    x1Pt, x2Pt);
                if (XWPt[2] >= 0) {
                    nPosZ++;
                }
                for (ii = 0; ii < 4; ++ii) {
                    XW[ii][i] = XWPt[ii];
                } 
            }
            if (nPosZ > bestNPosZ) {
                bestNPosZ = nPosZ;
                bestR = rTst;
                bestT = tTst;
                bestLabel = label;
                bestXW = MatrixUtil.copy(XW);
            }
        }
        
        if (bestR == null) {
            return;
        }
        
        // copy into output variables:
        for (i = 0; i < bestR.length; ++i) {
            System.arraycopy(bestR[i], 0, rSelected[i], 0, bestR[i].length);
        }
        System.arraycopy(bestT, 0, tSelected, 0, bestT.length);
        
        System.out.println("choosing solution: " + bestLabel);
        //double estimatedRotY = Math.atan(R[0][2]/R[0][0]) * (180./Math.PI);
        double estimatedRotZ = Math.atan(-bestR[1][0]/bestR[1][1]) * (180./Math.PI);
        System.out.printf("estimated rotation in degrees about z axis from R=%.2f\n", estimatedRotZ);
        //System.out.printf("X_WCS=\n%s\n", FormatArray.toString(bestXW, "%.3e"));
        System.out.flush();
        
        for (i = 0; i < XW.length; ++i) {
            System.arraycopy(XW[i], 0, outputX[i], 0, XW[i].length);
        }
        
    }
   
    private static DenseMatrix extractIndices(DenseMatrix m, List<Integer> inlierIndexes) {
        DenseMatrix out = new DenseMatrix(m.numRows(), inlierIndexes.size());
        int r = 0;
        for (int i = 0; i < inlierIndexes.size(); ++i) {
            int idx = inlierIndexes.get(i);
            for (int j = 0; j < m.numRows(); ++j) {
                out.add(j, r, m.get(j, idx));
            }
            r++;
        }
        return out;
    }

    static void populateWithDet1Rs(MatrixUtil.SVDProducts svdE, 
        double[][] r1Out, double[][] r2Out, double[][] uOut) {
        
        //Szeliski 2010, eqn (7.25)
        
        // R_Z+90 and R_Z_-90 from 
        // Ma, Soatto, Kosecká, and Sastry, "An Invitation to 3-D Vision"
        
        //R_z_90^T  = [ [0, 1, 0], [0, -1, 0], [0, 0, 1] ]
        //R_z_-90^T = [ [0, -1, 0], [0, 1, 0], [0, 0, 1] ]
         
        //R_z_90
        double[][] r90T = new double[3][3];
        r90T[0] = new double[]{0, 1, 0};
        r90T[1] = new double[]{-1, 0, 0};
        r90T[2] = new double[]{0, 0, 1};
        double[][] r90NegT = MatrixUtil.transpose(r90T);
        
        double[][] u = svdE.u;
        double[][] uNeg = MatrixUtil.copy(u);
        MatrixUtil.multiply(uNeg, -1);
        
        double[][] rr1 = MatrixUtil.multiply(MatrixUtil.multiply(u, r90T), svdE.vT);
        double[][] rr2 = MatrixUtil.multiply(MatrixUtil.multiply(uNeg, r90T), svdE.vT);
        double[][] rr3 = MatrixUtil.multiply(MatrixUtil.multiply(u, r90NegT), svdE.vT);
        double[][] rr4 = MatrixUtil.multiply(MatrixUtil.multiply(uNeg, r90NegT), svdE.vT);
        
        double det1 = MatrixUtil.determinant(rr1);
        double det2 = MatrixUtil.determinant(rr2);
        double det3 = MatrixUtil.determinant(rr3);
        double det4 = MatrixUtil.determinant(rr4);
        
        System.out.printf("det(r1,r2,r3,r4)=%.3e,%.3e,%.3e,%.3e\n", det1, det2, det3, det4);
        System.out.printf("r1:\n%s\n", FormatArray.toString(rr1, "%.4e"));
        System.out.printf("r2:\n%s\n", FormatArray.toString(rr2, "%.4e"));
        System.out.printf("r3:\n%s\n", FormatArray.toString(rr3, "%.4e"));
        System.out.printf("r4:\n%s\n", FormatArray.toString(rr4, "%.4e"));
        
        boolean useUPos = true;
        
        int i;
        if (Math.abs(det1 - 1.) < eps) {
            System.out.printf("using +U\n");
            for (i = 0; i < 3; ++i) {
                System.arraycopy(rr1[i], 0, r1Out[i], 0, rr1[i].length);
            }
        } else if (Math.abs(det2 - 1.) < eps) {
            System.out.printf("using -U\n");
            useUPos = false;
            for (i = 0; i < 3; ++i) {
                System.arraycopy(rr2[i], 0, r1Out[i], 0, rr2[i].length);
            }
        } else {
            throw new IllegalStateException("neither rotation matrix is SO(3)");
        }
        if (Math.abs(det3 - 1.) < eps) {
            if (!useUPos) {
                throw new IllegalStateException("expecting to need +U");
            }
            System.out.printf("using +U\n");
            for (i = 0; i < 3; ++i) {
                System.arraycopy(rr3[i], 0, r2Out[i], 0, rr3[i].length);
            }
        } else if (Math.abs(det4 - 1.) < eps) {
            if (useUPos) {
                throw new IllegalStateException("expecting to need -U");
            }
            System.out.printf("using -U\n");
            for (i = 0; i < 3; ++i) {
                System.arraycopy(rr4[i], 0, r2Out[i], 0, rr4[i].length);
            }
        } else {
            throw new IllegalStateException("neither rotation matrix is SO(3)");
        }
        
        if (useUPos) {
            for (i = 0; i < 3; ++i) {
                System.arraycopy(u[i], 0, uOut[i], 0, u[i].length);
            }
        } else {
            for (i = 0; i < 3; ++i) {
                System.arraycopy(uNeg[i], 0, uOut[i], 0, uNeg[i].length);
            }
        }
    }
}
