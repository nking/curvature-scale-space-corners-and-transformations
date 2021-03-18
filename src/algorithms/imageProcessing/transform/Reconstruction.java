package algorithms.imageProcessing.transform;

import algorithms.imageProcessing.features.RANSACSolver;
import algorithms.imageProcessing.matching.ErrorType;
import algorithms.matrix.MatrixUtil;
import algorithms.util.FormatArray;
import java.util.List;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.Matrices;
import no.uib.cipr.matrix.NotConvergedException;

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
    }
    
    /**
     * given correspondence between two images calculate the extrinsic camera
     * parameters, and the real world position.
     * 
     * <pre>
     * following lCMU lectures of Kris Kinai at 
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
    public static ReconstructionResults calculateReconstruction(
        double[][] k1, double[][] k2,
        double[][] x1, double[][] x2) throws NotConvergedException {
        
        if (x1.length != 3 || x2.length != 3) {
            throw new IllegalArgumentException("x1.length must be 3 and so must x2.length");
        }
        int n = x1.length;
        if (x2.length != n) {
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
        
        //Essential matrix: E = K2^T * FM * K1
        double[][] k2T = MatrixUtil.transpose(k2);
        double[][] _essentialM = MatrixUtil.multiply(k2T, _fm);
        _essentialM = MatrixUtil.multiply(_essentialM, k1);
        
        MatrixUtil.SVDProducts svdE = MatrixUtil.performSVD(_essentialM);
        
        assert(svdE.u[0].length == 3 && svdE.u.length == 3);

        System.out.printf("SVD.u=\n%s", FormatArray.toString(svdE.u, "%.3e"));
        System.out.printf("det(SVD.u)=%.2f\n", MatrixUtil.determinant(svdE.u));
        System.out.printf("SVD.vT=\n%s", FormatArray.toString(svdE.vT, "%.3e"));
        System.out.printf("det(SVD.vT)=%.2f\n", MatrixUtil.determinant(svdE.vT));
        
        // same as Hartley1992's E^T
        double[][] w = new double[3][3];
        w[0] = new double[]{0, -1, 0};
        w[1] = new double[]{1, 0, 0};
        w[2] = new double[]{0, 0, 1};
        double[][] eH92 = MatrixUtil.transpose(w);
        double[][] zH92 = new double[3][3];
        zH92[0] = new double[]{0, -1, 0};
        zH92[1] = new double[]{1, 0, 0};
        zH92[2] = new double[]{0, 0, 0};

        double[][] R1 = MatrixUtil.multiply(svdE.u, w);
        R1 = MatrixUtil.multiply(R1, svdE.vT);

        double[][] R2 = MatrixUtil.multiply(svdE.u, MatrixUtil.transpose(w));
        R2 = MatrixUtil.multiply(R2, svdE.vT);

        DenseMatrix uM = new DenseMatrix(svdE.u);
        double[] t1 = Matrices.getColumn(uM, 2).getData();

        double[] t2 = Matrices.getColumn(uM, 2).getData();
        MatrixUtil.multiply(t2, -1);
        
        /*
        // R1Hartley92 is the same as R2 of Kitani's lecture notes
        // R2Hartley92 is the same as R1 of Kitani's lecture notes
        double[][] R1Hartley92 = MatrixUtil.multiply(svdE.u, eH92);
        R1Hartley92 = MatrixUtil.multiply(R1Hartley92, svdE.vT);
        double[][] R2Hartley92 = MatrixUtil.multiply(svdE.u, w);
        R2Hartley92 = MatrixUtil.multiply(R2Hartley92, svdE.vT);
        double[][] SHartley92 = MatrixUtil.multiply(MatrixUtil.transpose(svdE.vT), zH92);
        SHartley92 = MatrixUtil.multiply(SHartley92, svdE.vT);
        double[][] Q1Hartley92 = MatrixUtil.multiply(R1Hartley92, SHartley92);
        double[][] Q2Hartley92 = MatrixUtil.multiply(R2Hartley92, SHartley92);
        */
        double detR1 = MatrixUtil.determinant(R1);
        double detR2 = MatrixUtil.determinant(R2);
        System.out.printf("R1=\n%s\n", FormatArray.toString(R1, "%.3e"));
        System.out.printf("R2=\n%s\n", FormatArray.toString(R2, "%.3e"));
        System.out.printf("t1=\n%s\n", FormatArray.toString(t1, "%.3e"));
        System.out.printf("t2=\n%s\n", FormatArray.toString(t2, "%.3e"));
        System.out.printf("det(R1)=%.3e\n", detR1);
        System.out.printf("det(R2)=%.3e\n\n", detR2);
        /*System.out.printf("R1_Hartly92=\n%s\n", FormatArray.toString(R1Hartley92, "%.3e"));
        System.out.printf("R2_Hartely92=\n%s\n", FormatArray.toString(R2Hartley92, "%.3e"));
        System.out.printf("S_Hartly92=\n%s\n", FormatArray.toString(SHartley92, "%.3e"));
        System.out.printf("Q1_Hartly92=\n%s\n", FormatArray.toString(Q1Hartley92, "%.3e"));
        System.out.printf("Q2_Hartly92=\n%s\n", FormatArray.toString(Q2Hartley92, "%.3e"));
        System.out.flush();
        */
        
        // solution 1:  R1 and T1
        // solution 2:  R1 and T2
        // solution 3:  R2 and T2
        // solution 4:  R2 and T1
        double[][] rSelected = MatrixUtil.zeros(3, 3);
        double[] tSelected = new double[3];
        double[][] XW = chooseRAndT(x1, x2, k1, k2,
            R1, R2, t1, t2, rSelected, tSelected);
        
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
    
    
    /**
     * 
     * @param x1
     * @param x2
     * @param k1
     * @param k2
     * @param R1
     * @param R2
     * @param t1
     * @param t2
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
        
        // valid equation has det(R) = 1
        //   if =-1, its a reflection
        double detR1 = MatrixUtil.determinant(R1);
        double detR2 = MatrixUtil.determinant(R2);
        
        //Compute 3D point using triangulation, valid solution has positive Z value
        // (Note: negative Z means point is behind the camera )
        boolean detR1Pos1 = Math.abs(detR1 - 1.) < 1e-4;
        boolean detR2Pos1 = Math.abs(detR2 - 1.) < 1e-4;
        
        // save the first that pass the tests for Z>=0.
        double[][] R = null;
        double[] T = null;
        double[][] XW;
        double[] XWPt;
        String goodSolnLabel = null;
        String label =null;
        
        XWPt = new double[4];
        XW = new double[3][n];
        for (int i = 0; i < 3; ++i) {
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
        
        boolean goodSoln;
        for (j = 0; j < 4; ++j) {
            goodSoln = true;
            switch(j) {
                case 0: {
                    if (!detR1Pos1) {
                        goodSoln = false;
                        break;
                    }
                    label = "R1, T1";
                    rTst = R1;
                    tTst = t1;
                    break;
                }
                case 1: {
                    if (!detR1Pos1) {
                        goodSoln = false;
                        break;
                    }
                    label = "R1, T2";
                    rTst = R1;
                    tTst = t2;
                    break;
                }
                case 2: {
                    if (!detR2Pos1) {
                        goodSoln = false;
                        break;
                    }
                    label = "R2, T1";
                    rTst = R2;
                    tTst = t1;
                    break;
                }
                default: {
                    if (!detR2Pos1) {
                        goodSoln = false;
                        break;
                    }
                    label = "R2, T2";
                    rTst = R2;
                    tTst = t2;
                    break;
                }
            }
            if (!goodSoln) {
                continue;
            }
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
                if (XWPt[3] < 0) {
                    goodSoln = false;
                    // break out of i loop
                    break;
                }
                for (ii = 0; ii < 4; ++ii) {
                    XW[ii][i] = XWPt[ii];
                } 
            }
            if (goodSoln) {
                goodSolnLabel = label;
                R = rTst;
                T = tTst;
                break;
            }
        }
        
        if (R == null) {
            //TODO: wrap the solution finiding in a retry.  may need a retry loop for adjusting f
            // once a good solution is found too.
            return null;
        }
        
        // copy into output variables:
        for (i = 0; i < R.length; ++i) {
            System.arraycopy(R[i], 0, rSelected[i], 0, R[i].length);
        }
        System.arraycopy(T, 0, tSelected, 0, T.length);
        
        System.out.println("choosing solution: " + goodSolnLabel);
        //double estimatedRotY = Math.atan(R[0][2]/R[0][0]) * (180./Math.PI);
        double estimatedRotY = Math.atan(-R[2][0]/R[2][2]) * (180./Math.PI);
        System.out.printf("estimated rotation about y axis from R=%.2f\n", estimatedRotY);
        System.out.printf("X_WCS=\n%s\n", FormatArray.toString(XW, "%.3e"));
        System.out.flush();
        
        return XW;
    }
 
}
