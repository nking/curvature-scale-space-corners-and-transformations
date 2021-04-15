package algorithms.imageProcessing.transform;

import algorithms.imageProcessing.features.RANSACSolver;
import algorithms.imageProcessing.matching.ErrorType;
import algorithms.matrix.MatrixUtil;
import algorithms.util.FormatArray;
import java.util.Arrays;
import java.util.List;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.NotConvergedException;
import no.uib.cipr.matrix.RQ;

/**
 * given correspondence between two images calculate the camera
 * parameters as intrinsic and extrinsic parameters,
 * and the real world position.
 * 
 * Euler rotations:
        
        about z-axis (yaw):           about x-axis (roll):       about the y-axis (pitch):
            | cos φ   -sin φ    0 |    |    1       0       0 |  |  cos ψ    0  sin ψ |
            | sin φ    cos φ    0 |    |    0   cos θ   sin θ |  |      0    1      0 |
            |     0        0    1 |    |    0  -sin θ   cos θ |  | -sin ψ    0  cos ψ |        
        
 * useful reading:
 * <pre>
 * http://www.cs.cmu.edu/~16385/s17/Slides/12.5_Reconstruction.pdf
 * add other references here
 * </pre>
 * @author nichole
 */
public class Reconstruction {

    public static class ReconstructionResults {
        double[][] XW;
        double[][] k1Intr;
        double[][] k2Intr;
        double[][] k1ExtrRot;
        double[] k1ExtrTrans;
        double[][] k2ExtrRot;
        double[] k2ExtrTrans;
        
        @Override
        public String toString() {
            StringBuilder sb = new StringBuilder();
            sb.append("XW=\n");
            if (XW != null) {
                sb.append(FormatArray.toString(XW, "%.4e"));
            }
            sb.append("k1 intrinsic=\n");
            if (k1Intr != null) {
                sb.append(FormatArray.toString(k1Intr, "%.4e"));
            }
            sb.append("k1 extrinsic rotation=\n");
            if (k1ExtrRot != null) {
                sb.append(FormatArray.toString(k1ExtrRot, "%.4e"));
            }
            sb.append("k1 extrinsic translation=\n");
            if (k1ExtrTrans != null) {
                sb.append(FormatArray.toString(k1ExtrTrans, "%.4e"));
                sb.append("\n");
            }
            sb.append("k2 intrinsic=\n");
            if (k2Intr != null) {
                sb.append(FormatArray.toString(k2Intr, "%.4e"));
            }
            sb.append("k2 extrinsic rotation=\n");
            if (k2ExtrRot != null) {
                sb.append(FormatArray.toString(k2ExtrRot, "%.4e"));
            }
            sb.append("k2 extrinsic translation=\n");
            if (k2ExtrTrans != null) {
                sb.append(FormatArray.toString(k2ExtrTrans, "%.4e"));
                sb.append("\n");
            }
            return sb.toString();
        }
    }
    
    public static class CameraExtrinsics {
        double[][] rot;
        double[] trans;
       
        @Override
        public String toString() {
            StringBuilder sb = new StringBuilder();
            sb.append("rot=\n");
            if (rot != null) {
                sb.append(FormatArray.toString(rot, "%.4e"));
            }
            sb.append("trans=\n");
            if (trans != null) {
                sb.append(FormatArray.toString(trans, "%.4e"));
            }
            return sb.toString();
        }
    }
    
     /**
     * given correspondence between two images calculate the camera
     * parameters, and the real world position.
     * 
     * <pre>
     * following CMU lectures of Kris Kinai at 
     * http://www.cs.cmu.edu/~16385/s17/Slides/12.5_Reconstruction.pdf
     * 
     * add other references here
     * </pre>
     * @param camera1 image 1 camera matrix of intrinsic and extrinsic parameters.
     * the size is 3 x 4.
     * @param camera2 image 2 camera matrix of intrinsic and extrinsic parameters.
     * the size is 3 x 4.
     * @param x1 the image 1 set of correspondence points.  format is 3 x N where
     * N is the number of points.
     * @param x2 the image 2 set of correspondence points.  format is 3 x N where
     * N is the number of points.
     * @return 
     * @throws no.uib.cipr.matrix.NotConvergedException 
     */
    public static ReconstructionResults calculateReconstruction(
        double[][] camera1, double[][] camera2,
        double[][] x1, double[][] x2) throws NotConvergedException {
        
        if (x1.length != 3 || x2.length != 3) {
            throw new IllegalArgumentException("x1.length must be 3 and so must x2.length");
        }
        int n = x1[0].length;
        if (x2[0].length != n) {
            throw new IllegalArgumentException("x1 and x2 must be same dimensions");
        }
        
        /*
        http://www.cs.cmu.edu/~16385/s17/Slides/12.5_Reconstruction.pdf
        
        (3) For each point correspondence, compute the point X in 3D space (triangularization)
        */
        
        double[][] XW = new double[4][n];
        for (int i = 0; i < 4; ++i) {
            XW[i] = new double[n];
        }
        double[] XWPt = new double[4];
        
        double[][] x1Pt = new double[3][1];
        double[][] x2Pt = new double[3][1];
        int i, j, ii;
        for (i = 0; i < 3; ++i) {
            x1Pt[i] = new double[1];
            x2Pt[i] = new double[1];
        }
                    
        for (i = 0; i < n; ++i) {
            for (ii = 0; ii < 3; ++ii) {
                x1Pt[ii][0] = x1[ii][i];
                x2Pt[ii][0] = x2[ii][i];
            }
            //
            XWPt = Triangulation.calculateWCSPoint(
                camera1, camera2, x1Pt, x2Pt);
            for (ii = 0; ii < 4; ++ii) {
                XW[ii][i] = XWPt[ii];
            } 
        }
                
        ReconstructionResults rr = new ReconstructionResults();
        rr.XW = XW;

        return rr;
    }
    
    /**
     * given correspondence between two images calculate the camera
     * parameters, and the real world position.
     * 
     * <pre>
     * following CMU lectures of Kris Kinai at 
     * http://www.cs.cmu.edu/~16385/s17/Slides/12.5_Reconstruction.pdf
     * 
     * add other references here
     * </pre>
     * @param x1 the image 1 set of correspondence points.  format is 3 x N where
     * N is the number of points.
     * @param x2 the image 2 set of correspondence points.  format is 3 x N where
     * N is the number of points.
     * @return 
     */
    public static ReconstructionResults calculateReconstruction(
        double[][] x1, double[][] x2) throws NotConvergedException {
        
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
        
        DenseMatrix x1M = new DenseMatrix(x1);
        DenseMatrix x2M = new DenseMatrix(x2);
        
        EpipolarTransformer.NormalizedXY normXY1 = EpipolarTransformer.normalize(x1M);
        EpipolarTransformer.NormalizedXY normXY2 = EpipolarTransformer.normalize(x2M);
        DenseMatrix leftM = normXY1.getXy();
        DenseMatrix rightM = normXY2.getXy();
        
        double tolerance = 3.84; //3.84 5.99 7.82        
        boolean useToleranceAsStatFactor = true;
        ErrorType errorType = ErrorType.SAMPSONS;
        EpipolarTransformationFit fitR = null;
        boolean reCalcIterations = false;
        
        EpipolarTransformer tr = new EpipolarTransformer();
        
        /*
        DenseMatrix normalizedFM = tr.calculateEpipolarProjection(leftM, rightM);
        DenseMatrix vNFM = tr.validateSolution(normalizedFM, leftM, rightM);
        
        Distances distances = new Distances();
        if (useToleranceAsStatFactor) {
            fitR = distances.calculateError2(normalizedFM, leftM, rightM,
                    errorType, tolerance);
        } else {
            fitR = distances.calculateError(normalizedFM, leftM, rightM,
                    errorType, tolerance);
        }
        */
  
        // see sect 7.21 of szeliski 
                
        RANSACSolver solver = new RANSACSolver();
        fitR = solver.calculateEpipolarProjection(
            leftM, rightM, errorType, useToleranceAsStatFactor, tolerance,
                reCalcIterations);
        
        DenseMatrix fm = EpipolarTransformer.denormalizeTheFundamentalMatrix(
            fitR.getFundamentalMatrix(), 
            normXY1.getNormalizationMatrices(),
            normXY2.getNormalizationMatrices());
        
        double[][] _fm = MatrixUtil.convertToRowMajor(fm);
        
        x1M = extractIndices(x1M, fitR.inlierIndexes);
        x2M = extractIndices(x2M, fitR.inlierIndexes);
        x1 = MatrixUtil.convertToRowMajor(x1M);
        x2 = MatrixUtil.convertToRowMajor(x2M);
        
        System.out.println("RANSAC fit=" + fitR.toString());
        
        //(2) compute the camera matrices P1, P2 from FM.
        
        MatrixUtil.SVDProducts svd = MatrixUtil.performSVD(_fm);
        
        assert(svd.u[0].length == 3 && svd.u.length == 3);

        double detV = MatrixUtil.determinant(svd.vT);
        
        System.out.printf("SVD.u=\n%s\n", FormatArray.toString(svd.u, "%.3e"));
        System.out.printf("SVD.s=\n%s\n", FormatArray.toString(svd.s, "%.3e"));
        System.out.printf("SVD.vT=\n%s\n", FormatArray.toString(svd.vT, "%.3e"));
        System.out.printf("det(SVD.u)=%.2f\n", MatrixUtil.determinant(svd.u));
        System.out.printf("det(SVD.vT)=%.2f\n", detV);
        
        //NOTE: the last column vector in u is the smallest
        //    eigenvector.  it is epipole2, that is, the right image position 
        //    of the epipolar projection of the left camera center.
        //    it's int the left null space of E.
     
        // camera matrix P for left image
        double[][] camera1 = MatrixUtil.zeros(3, 4);
        for (int i = 0; i < 3; ++i) {
            camera1[i][i] = 1;
        }
        
        double[][] r90 = new double[3][3];
        r90[0] = new double[]{0, -1, 0};
        r90[1] = new double[]{1, 0, 0};
        r90[2] = new double[]{0, 0, 1};        
        double[][] r90T = MatrixUtil.transpose(r90);
        double[] s = Arrays.copyOf(svd.s, svd.s.length);
        s[2] = s[1];
        
        double[][] homog = MatrixUtil.multiply(svd.u, r90T);
        homog = MatrixUtil.multiplyByDiagonal(homog, s);
        homog = MatrixUtil.multiply(homog, svd.vT);
        
        // e1 is the left image view of the projected right camera nadir (= last col in V)
        // e2 is the right image view of the projected left camera nadir (= last col in U)
        double[][] epipoles = tr.calculateEpipoles(fm);
        
        double[] e2NotNormalized = MatrixUtil.transpose(svd.u)[2];
        
        System.out.printf("e1=%s\n", FormatArray.toString(epipoles[0], "%.3e"));
        System.out.printf("e2=%s\n", FormatArray.toString(epipoles[1], "%.3e"));
        
        // Szeleski: P2 =[homog ̃| e1]
        double[][] camera2 = MatrixUtil.zeros(3, 4);
        for (int i = 0; i < 3; ++i) {
            camera2[i] = new double[4];
            System.arraycopy(homog[i], 0, camera2[i], 0, 3);
            camera2[i][3] = e2NotNormalized[i];//epipoles[1][i];
        }
        
        System.out.printf("Szeliski camera2:\n%s\n", FormatArray.toString(camera2, "%.3e"));
                
        // Hartley & Zisserman: camera matrix P2 for right image = [ [e2_x] * F | e2 ]
        
        /*
        double[][] e2SkewSym = MatrixUtil.skewSymmetric(epipoles[1]);
        double[][] k2R = MatrixUtil.multiply(e2SkewSym, _fm);
        double[][] camera2 = new double[3][4];
        for (int i = 0; i < 3; ++i) {
            camera2[i] = new double[4];
            System.arraycopy(k2R[i], 0, camera2[i], 0, 3);
            camera2[i][3] = epipoles[1][i];
        }*/
        
        // Hartley and Zisserman vgg code uses:  [ -[e2_x] * F | e2 ]
        double[][] e2SkewSym 
            //= MatrixUtil.skewSymmetric(e2NotNormalized);
            = MatrixUtil.skewSymmetric(epipoles[1]);
        MatrixUtil.multiply(e2SkewSym, -1);
        double[][] k2R = MatrixUtil.multiply(e2SkewSym, _fm);
        double[][] camera2HZ = new double[3][4];
        for (int i = 0; i < 3; ++i) {
            camera2HZ[i] = new double[4];
            System.arraycopy(k2R[i], 0, camera2HZ[i], 0, 3);
            //camera2[i][3] = e2NotNormalized[i];
            camera2HZ[i][3] = epipoles[1][i];
        }
        System.out.printf("Hartley & Zisserman camera2:\n%s\n", FormatArray.toString(camera2HZ, "%.3e"));
        
            // can extract intrinsic and rotation from left 3 columns of camera2
            //    matrix using RQ decomposition
            // see http://www.cs.cmu.edu/~16385/s17/Slides/11.3_Pose_Estimation.pdf
            double[][] M = MatrixUtil.copySubMatrix(camera2, 0, 2, 0, 2);
            RQ rq = RQ.factorize(new DenseMatrix(M));
            double[][] k2Intr = MatrixUtil.convertToRowMajor(rq.getR());
            MatrixUtil.multiply(k2Intr, 1./k2Intr[2][2]);
            double[][] k2ExtrRot = MatrixUtil.convertToRowMajor(rq.getQ());
            System.out.printf("  decomposed into intrinsic=\n   %s\n", FormatArray.toString(k2Intr, "%.3e"));
            System.out.printf("  decomposed into extrinsic rotation=\n   %s\n", FormatArray.toString(k2ExtrRot, "%.3e"));
            
        return calculateReconstruction(camera1, camera2, x1, x2);
    }

     /**
     * given correspondence between two images calculate the extrinsic camera
     * parameters, and the real world position.
     * 
     * <pre>
     * following CMU lectures of Kris Kinai at 
     * http://www.cs.cmu.edu/~16385/s17/Slides/12.5_Reconstruction.pdf
     * 
     * add other references here
     * </pre>
     * @param k1 intrinsic camera matrix for image 1 in units of pixels.
     * @param k2 intrinsic camera matrix for image 2 in units of pixels.
     * @param x1 the image 1 set of correspondence points.  format is 3 x N where
     * N is the number of points.
     * @param x2 the image 2 set of correspondence points.  format is 3 x N where
     * N is the number of points.
     * @return 
     */
    public static ReconstructionResults calculateReconstructionWithIntrinsic(
        double[][] k1, double[][] k2,
        double[][] x1, double[][] x2) throws NotConvergedException {
        
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
        
        EpipolarTransformer tr = new EpipolarTransformer();
        
        /*
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
        
        RANSACSolver solver = new RANSACSolver();
        fitR = solver.calculateEpipolarProjection(
            leftM, rightM, errorType, useToleranceAsStatFactor, tolerance,
                reCalcIterations);
        
        System.out.println("RANSAC fit=" + fitR.toString());
        
        normalizedE = fitR.getFundamentalMatrix();
        
        DenseMatrix essentialM = EpipolarTransformer.denormalizeTheFundamentalMatrix(
            normalizedE, normXY1.getNormalizationMatrices(),
            normXY2.getNormalizationMatrices());
                
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
        szeliski:
        Once an estimate for the essential matrix E has been recovered, 
        the direction of the translation vector t can be estimated. 
        Note that the absolute distance between the two cameras can never 
        be recovered from pure image measurements alone, regardless of 
        how many cameras or points are used. Knowledge about absolute 
        camera and point positions or distances, of- ten called ground 
        control points in photogrammetry, is always required to 
        establish the final scale, position, and orientation.
        */
        double[] t1 = MatrixUtil.extractColumn(svdE.u, 2);
        double[] t2 = Arrays.copyOf(t1, t1.length);
        MatrixUtil.multiply(t2, -1); 
        
        //negative 90 transposed
        double[][] r90 = new double[3][3];
        r90[0] = new double[]{0, -1, 0};
        r90[1] = new double[]{1, 0, 0};
        r90[2] = new double[]{0, 0, 1};
        
  //TODO: recovering expected rotation angle... looks incorrect
          
        //positive 90 transposed
        double[][] r90T = MatrixUtil.transpose(r90);
        
        //R = ±U R^T±90^T V^T
        double[][] uNegative = MatrixUtil.copy(svdE.u);
        MatrixUtil.multiply(uNegative, -1.);
        
        double[][] R1 = MatrixUtil.multiply(svdE.u, r90T);
        R1 = MatrixUtil.multiply(R1, svdE.vT);
        
        double[][] R2 = MatrixUtil.multiply(svdE.u, r90);
        R2 = MatrixUtil.multiply(R2, svdE.vT);
        
        double[][] R3 = MatrixUtil.multiply(uNegative, r90T);
        R3 = MatrixUtil.multiply(R3, svdE.vT);
        
        double[][] R4 = MatrixUtil.multiply(uNegative, r90);
        R4 = MatrixUtil.multiply(R4, svdE.vT);
        
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
        double detR1 = MatrixUtil.determinant(R1);
        double detR2 = MatrixUtil.determinant(R2);
        double detR3 = MatrixUtil.determinant(R3);
        double detR4 = MatrixUtil.determinant(R4);
        
        double[][] R3N = MatrixUtil.copy(R3);
        MatrixUtil.multiply(R3N, 1./R3N[2][2]);
        double[][] R4N = MatrixUtil.copy(R4);
        MatrixUtil.multiply(R4N, 1./R4N[2][2]);
        
        System.out.printf("t1=\n%s\n", FormatArray.toString(t1, "%.3e"));
        System.out.printf("t2=\n%s\n", FormatArray.toString(t2, "%.3e"));
        System.out.printf("R1=\n%s\n", FormatArray.toString(R1, "%.3e"));
        System.out.printf("R2=\n%s\n", FormatArray.toString(R2, "%.3e"));
        System.out.printf("R3=\n%s\n", FormatArray.toString(R3, "%.3e"));
        System.out.printf("R3N=\n%s\n", FormatArray.toString(R3N, "%.3e"));
        System.out.printf("R4=\n%s\n", FormatArray.toString(R4, "%.3e"));
        System.out.printf("R4N=\n%s\n", FormatArray.toString(R4N, "%.3e"));
        System.out.printf("det(R1)=%.3e\n", detR1);
        System.out.printf("det(R2)=%.3e\n\n", detR2);
        System.out.printf("det(R3)=%.3e\n", detR3);
        System.out.printf("det(R4)=%.3e\n\n", detR4);
        
        double tol = 1.e-5;
        
        // keep the 2 that have det(R) == 1 
        double[][] rot1 = null; 
        double[][] rot2 = null;
        if (Math.abs(detR1 - 1.) < tol) {
            rot1 = R1;
        }
        if (Math.abs(detR2 - 1.) < tol) {
            if (rot1 == null) {
                rot1 = R2;
            } else {
                rot2 = R2;
            }
        }
        if (Math.abs(detR3 - 1.) < tol) {
            if (rot1 == null) {
                rot1 = R3;
            } else {
                rot2 = R3;
            }
        }
        if (Math.abs(detR4 - 1.) < tol) {
            if (rot1 == null) {
                rot1 = R4;
            } else {
                rot2 = R4;
            }
        }
        
        if (rot1 == null && rot2 == null) {
            return null;
        }
        
        //then of the 4 possible choices find the one with largest number of positive Z.
          
        //NOTE: the last column vector in u is the smallest
        //    eigenvector.  it is epipole2, that is, the right image position 
        //    of the epipolar projection of the left camera center.
        //    it's int the left null space of E.
           
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
        double[][] XW = chooseRAndT(x1, x2, k1, k2,
            rot1, rot2, t1, t2, rSelected, tSelected);
        
        if (XW == null) {
            return null;
        }
                
        ReconstructionResults rr = new ReconstructionResults();
        rr.XW = XW;
        rr.k2ExtrRot = rSelected;
        rr.k2ExtrTrans = tSelected;
        rr.k2Intr = k2;
        rr.k1ExtrRot = MatrixUtil.createIdentityMatrix(3);
        rr.k1ExtrTrans = new double[tSelected.length];
        rr.k1Intr = k1;

        return rr;
    }
    
     /**
     * given correspondence between two images calculate the extrinsic camera
     * parameters.
     * 
     * Can construct the essential matrix with the 2nd camera extrinsics: [R | R*t ]
     * 
     * <pre>
     * following CMU lectures of Kris Kinai at 
     * http://www.cs.cmu.edu/~16385/s17/Slides/12.5_Reconstruction.pdf
     * 
     * add other references here
     * </pre>
     * @param k1 intrinsic camera matrix for image 1 in units of pixels.
     * @param k2 intrinsic camera matrix for image 2 in units of pixels.
     * @param x1 the image 1 set of correspondence points.  format is 3 x N where
     * N is the number of points.
     * @param x2 the image 2 set of correspondence points.  format is 3 x N where
     * N is the number of points.
     * @return 
     */
    public static CameraExtrinsics[] calculateCameraExtrinsics(
        double[][] k1, double[][] k2,
        double[][] x1, double[][] x2) throws NotConvergedException {
        
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
        
        EpipolarTransformer tr = new EpipolarTransformer();
        
        /*
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
        
        RANSACSolver solver = new RANSACSolver();
        fitR = solver.calculateEpipolarProjection(
            leftM, rightM, errorType, useToleranceAsStatFactor, tolerance,
                reCalcIterations);
        
        System.out.println("RANSAC fit=" + fitR.toString());
        
        normalizedE = fitR.getFundamentalMatrix();
        
        DenseMatrix essentialM = EpipolarTransformer.denormalizeTheFundamentalMatrix(
            normalizedE, normXY1.getNormalizationMatrices(),
            normXY2.getNormalizationMatrices());
                
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
        szeliski:
        Once an estimate for the essential matrix E has been recovered, 
        the direction of the translation vector t can be estimated. 
        Note that the absolute distance between the two cameras can never 
        be recovered from pure image measurements alone, regardless of 
        how many cameras or points are used. Knowledge about absolute 
        camera and point positions or distances, of- ten called ground 
        control points in photogrammetry, is always required to 
        establish the final scale, position, and orientation.
        */
        double[] t1 = MatrixUtil.extractColumn(svdE.u, 2);
        double[] t2 = Arrays.copyOf(t1, t1.length);
        MatrixUtil.multiply(t2, -1); 
        
        //negative 90 transposed
        double[][] r90 = new double[3][3];
        r90[0] = new double[]{0, -1, 0};
        r90[1] = new double[]{1, 0, 0};
        r90[2] = new double[]{0, 0, 1};
        
  //TODO: recovering expected rotation angle... looks incorrect
          
        //positive 90 transposed
        double[][] r90T = MatrixUtil.transpose(r90);
        
        //R = ±U R^T±90^T V^T
        double[][] uNegative = MatrixUtil.copy(svdE.u);
        MatrixUtil.multiply(uNegative, -1.);
        
        double[][] R1 = MatrixUtil.multiply(svdE.u, r90T);
        R1 = MatrixUtil.multiply(R1, svdE.vT);
        
        double[][] R2 = MatrixUtil.multiply(svdE.u, r90);
        R2 = MatrixUtil.multiply(R2, svdE.vT);
        
        double[][] R3 = MatrixUtil.multiply(uNegative, r90T);
        R3 = MatrixUtil.multiply(R3, svdE.vT);
        
        double[][] R4 = MatrixUtil.multiply(uNegative, r90);
        R4 = MatrixUtil.multiply(R4, svdE.vT);
        
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
        double detR1 = MatrixUtil.determinant(R1);
        double detR2 = MatrixUtil.determinant(R2);
        double detR3 = MatrixUtil.determinant(R3);
        double detR4 = MatrixUtil.determinant(R4);
        
        double[][] R3N = MatrixUtil.copy(R3);
        MatrixUtil.multiply(R3N, 1./R3N[2][2]);
        double[][] R4N = MatrixUtil.copy(R4);
        MatrixUtil.multiply(R4N, 1./R4N[2][2]);
        
        System.out.printf("t1=\n%s\n", FormatArray.toString(t1, "%.3e"));
        System.out.printf("t2=\n%s\n", FormatArray.toString(t2, "%.3e"));
        System.out.printf("R1=\n%s\n", FormatArray.toString(R1, "%.3e"));
        System.out.printf("R2=\n%s\n", FormatArray.toString(R2, "%.3e"));
        System.out.printf("R3=\n%s\n", FormatArray.toString(R3, "%.3e"));
        System.out.printf("R3N=\n%s\n", FormatArray.toString(R3N, "%.3e"));
        System.out.printf("R4=\n%s\n", FormatArray.toString(R4, "%.3e"));
        System.out.printf("R4N=\n%s\n", FormatArray.toString(R4N, "%.3e"));
        System.out.printf("det(R1)=%.3e\n", detR1);
        System.out.printf("det(R2)=%.3e\n\n", detR2);
        System.out.printf("det(R3)=%.3e\n", detR3);
        System.out.printf("det(R4)=%.3e\n\n", detR4);
        
        double tol = 1.e-5;
        
        // keep the 2 that have det(R) == 1 
        double[][] rot1 = null; 
        double[][] rot2 = null;
        if (Math.abs(detR1 - 1.) < tol) {
            rot1 = R1;
        }
        if (Math.abs(detR2 - 1.) < tol) {
            if (rot1 == null) {
                rot1 = R2;
            } else {
                rot2 = R2;
            }
        }
        if (Math.abs(detR3 - 1.) < tol) {
            if (rot1 == null) {
                rot1 = R3;
            } else {
                rot2 = R3;
            }
        }
        if (Math.abs(detR4 - 1.) < tol) {
            if (rot1 == null) {
                rot1 = R4;
            } else {
                rot2 = R4;
            }
        }
        
        if (rot1 == null && rot2 == null) {
            return null;
        }
        
        //then of the 4 possible choices find the one with largest number of positive Z.
          
        //NOTE: the last column vector in u is the smallest
        //    eigenvector.  it is epipole2, that is, the right image position 
        //    of the epipolar projection of the left camera center.
        //    it's int the left null space of E.
           
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
        double[][] XW = chooseRAndT(x1, x2, k1, k2,
            rot1, rot2, t1, t2, rSelected, tSelected);
        
        CameraExtrinsics c1 = new CameraExtrinsics();
        c1.rot = MatrixUtil.createIdentityMatrix(3);
        c1.trans = new double[]{0, 0, 0};
          
        CameraExtrinsics c2 = new CameraExtrinsics();
        c2.rot = rSelected;
        c2.trans = tSelected;
        
        return new CameraExtrinsics[]{c1, c2};
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
     * @return the real world coordinates of the projection of x1 and x2 using
     * triangulation. else null if no valid solution was found
     */
    private static double[][] chooseRAndT(double[][] x1, double[][] x2, 
        double[][] k1, double[][] k2,
        double[][] R1, double[][] R2, double[] t1, double[] t2, 
        double[][] rSelected, double[] tSelected) {
    
        int n = x1[0].length;
        
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
            return null;
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
        System.out.printf("X_WCS=\n%s\n", FormatArray.toString(bestXW, "%.3e"));
        System.out.flush();
        
        return bestXW;
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
 
}
