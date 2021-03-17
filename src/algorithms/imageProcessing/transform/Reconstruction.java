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
        
        // solution 1:  R1 and T1
        // solution 2:  R1 and T2
        // solution 3:  R2 and T2
        // solution 4:  R2 and T1
        double[][] rSelected = MatrixUtil.zeros(3, 3);
        double[] tSelected = new double[3];
        double[] XW = chooseRAndT(x1, x2, k1, k2,
            R1, R2, t1, t2, rSelected, tSelected);
        editing
        if (XW == null) {
            return null;
        }
        
        double[][] i3 = new double[3][3];
        for (int i = 0; i < 3; ++i) {
            i3[i] = new double[3];
            i3[i][i] = 1;
        }
        
        ReconstructionResults rr = new ReconstructionResults();
        rr.XW = XW;
        rr.k2ExtrRot = rSelected;
        rr.k2ExtrTrans = tSelected;
        rr.k2Intr = k2;
        rr.k1ExtrRot = i3;
        rr.k1ExtrTrans = new double[tSelected.length];
        rr.k1Intr = k1;

        //Compute 3D point using triangulation, valid solution has positive Z value
        // (Note: negative Z means point is behind the camera )
        /*
        System.out.printf("R1=\n%s\n", FormatArray.toString(R1, "%.3e"));
        System.out.printf("R2=\n%s\n", FormatArray.toString(R2, "%.3e"));
        System.out.printf("t1=\n%s\n", FormatArray.toString(t1, "%.3e"));
        System.out.printf("t2=\n%s\n", FormatArray.toString(t2, "%.3e"));
        System.out.printf("det(R1)=%.3e\n", detR1);
        System.out.printf("det(R2)=%.3e\n\n", detR2);
        System.out.printf("R1_Hartly92=\n%s\n", FormatArray.toString(R1Hartley92, "%.3e"));
        System.out.printf("R2_Hartely92=\n%s\n", FormatArray.toString(R2Hartley92, "%.3e"));
        System.out.printf("S_Hartly92=\n%s\n", FormatArray.toString(SHartley92, "%.3e"));
        System.out.printf("Q1_Hartly92=\n%s\n", FormatArray.toString(Q1Hartley92, "%.3e"));
        System.out.printf("Q2_Hartly92=\n%s\n", FormatArray.toString(Q2Hartley92, "%.3e"));
        System.out.flush();
        */
        
        return null;
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
    private static double[][] chooseRAndT(double[][] x1, double[][] x2, double[][] k1, double[][] k2,
        double[][] R1, double[][] R2, double[] t1, double[] t2, 
        double[][] rSelected, double[] tSelected) {
    
        // valid equation has det(R) = 1 (rotation and reflection)
        //   if =-1, its a reflection
        double detR1 = MatrixUtil.determinant(R1);
        double detR2 = MatrixUtil.determinant(R2);
        
        // save the first that pass the tests for Z>=0.
        double[][] R = null;
        double[] T = null;
        double[][] XW;
        double[] XWPt;
        String goodSolnLabel = null;
        
        
        boolean goodSoln;
        
        editing here to iterate over each correspondence pair
        
        // check R1 and t1
        XWPt = Triangulation.calculateWCSPoints(k1, R1, t1, k2, R2, t2, x1Pt, x2Pt);
        goodSoln = ((Math.abs(detR1) - 1.) < 1e-5) && (XW[3] > 0);
        System.out.printf("Good Solution=%b\nX from inv(|[R1|t1])*x2=\n%s\n", 
            goodSoln, FormatArray.toString(XW, " %.0f"));
        if (goodSoln && R == null) {
            goodSolnLabel = "R1, T1";
            R = R1;
            T = t1;
        } 
        
        // check R1 and t2
        XWPt = Triangulation.calculateWCSPoints(k1, R1, t2, k2, R2, t2, x1Pt, x2Pt);
        goodSoln = ((Math.abs(detR1) - 1.) < 1e-5) && (XW[3] > 0);
        System.out.printf("Good Solution=%b\nX from inv(|[R1|t2])*x2=\n%s\n", 
            goodSoln, FormatArray.toString(XW, " %.0f"));
        if (goodSoln && R == null) {
            goodSolnLabel = "R1, T2";
            R = R1;
            T = t2;
        } 
        
        // check R2 and t1
        XWPt = Triangulation.calculateWCSPoints(k1, R2, t1, k2, R2, t2, x1Pt, x2Pt);
        goodSoln = ((Math.abs(detR1) - 1.) < 1e-5) && (XW[3] > 0);
        System.out.printf("Good Solution=%b\nX from inv(|[R2|t1])*x2=\n%s\n", 
            goodSoln, FormatArray.toString(XW, " %.0f"));
        if (goodSoln && R == null) {
            goodSolnLabel = "R2, T1";
            R = R2;
            T = t1;
        } 
        
        // check R2 and t2
        XW = Triangulation.calculateWCSPoints(k1, R2, t2, k2, R2, t2, x1Pt, x2Pt);
        goodSoln = ((Math.abs(detR1) - 1.) < 1e-5) && (XW[3] > 0);
        System.out.printf("Good Solution=%b\nX from inv(|[R2|t2])*x2=\n%s\n", goodSoln, FormatArray.toString(XW, " %.0f"));
        if (goodSoln && R == null) {
            goodSolnLabel = "R2, T2";
            R = R2;
            T = t2;
        } 
        
        if (R == null) {
            //TODO: wrap the solution finiding in a retry.  may need a retry loop for adjusting f
            // once a good solution is found too.
            return null;
        }
        
        for (int i = 0; i < R.length; ++i) {
            System.arraycopy(R[i], 0, rSelected[i], 0, R[i].length);
        }
        System.arraycopy(T, 0, tSelected, 0, T.length);
        
        System.out.println("choosing solution: " + goodSolnLabel);
        //double estimatedRotY = Math.atan(R[0][2]/R[0][0]) * (180./Math.PI);
        double estimatedRotY = Math.atan(-R[2][0]/R[2][2]) * (180./Math.PI);
        System.out.printf("estimated rotation about y axis from R=%.2f\n", estimatedRotY);
        System.out.flush();
        
        return XW;
    }
 
    public static boolean allZsArePositive(double[][] XW) {
        return allOfDimensionArePositive(XW, 2);
    }
    public static boolean allOfDimensionArePositive(double[][] XW, final int dimension) {
        for (int i = 0; i < XW[0].length; ++i) {
            if (XW[dimension][i] < 0) {
                return false;
            }
        }
        return true;
    }
}
