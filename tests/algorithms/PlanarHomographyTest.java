package algorithms;

import algorithms.imageProcessing.transform.Rotation;
import algorithms.matrix.MatrixUtil;
import algorithms.misc.MiscMath0;
import algorithms.statistics.UnivariateNormalDistribution;
import algorithms.imageProcessing.transform.CameraCalibration;
import algorithms.util.FormatArray;
import junit.framework.TestCase;
import no.uib.cipr.matrix.NotConvergedException;

import java.util.Arrays;
import java.util.Random;

/**
 *
 * @author nichole
 */
public class PlanarHomographyTest extends TestCase {

    /*
    4 - point pairs homography with constraint that no 3 of the 4 are collinear.
    (1) normalize the homography Lemma 5.18 of Masks, end of Section 5.3.2.
    (2) decompose the planar homography into motion and structure parameters
    (R, (1/d)*T, N) where R is the rotation matrix transforming frame 1 to frame 2,
    T is the translation vector transforming frame 1 into frame 2,
    d is the distance from the world coordinate plane P to the optical center of the first camera,
    d > 0, and N is the unit normal vector of the plane P with respect to
    the first camera frame.
    see Figure 5.10.
     */
    public void test0() throws NotConvergedException {

        long seed = System.nanoTime();
        //seed = 188903032980675L;
        System.out.println("SEED=" + seed);
        Random rand = new Random(seed);

        double angle = Math.PI/10.;
        double d = 5;
        double lambda = 4;
        // rotation is cc about the y-axis (pitch)
        double[][] r = new double[3][];
        r[0] = new double[]{Math.cos(angle), 0, Math.sin(angle)};
        r[1] = new double[]{0, 1, 0};
        r[2] = new double[]{-Math.sin(angle), 0, Math.cos(angle)};

        int i;
        int j;

        double[] t = new double[]{2, 0, 0};
        double[] n = new double[]{1, 0, 2}; // chosen to not have sqrt(sqsum))=1

        // then the plane P is tangent to n
        // which is (-2, 0, 1) or any multiple of that.

        System.out.printf("angle of rotation about y axis =%.3f, depth d=%.1f, lambda scale=%.1f\n",
                angle, d, lambda);
        System.out.printf("from camera1 to camera2, Rotation=\n%s\n", FormatArray.toString(r, "%.3e"));
        System.out.printf("from camera1 to camera2, T=\n%s\n", FormatArray.toString(t, "%.3e"));
        System.out.printf("perpendicular to plane P: N=\n%s\n", FormatArray.toString(n, "%.3e"));

        double[][] hl = MatrixUtil.outerProduct(t, n);
        MatrixUtil.multiply(hl, 1./d);
        hl = MatrixUtil.pointwiseAdd(r, hl);
        MatrixUtil.multiply(hl, lambda);

        System.out.printf("hl=\n%s\n", FormatArray.toString(hl, "%.3f"));

        MatrixUtil.SVDProducts svd = MatrixUtil.performSVD(hl);

        System.out.printf("SVD(hl).s=\n%s\n", FormatArray.toString(svd.s, "%.3f"));
        System.out.printf("SVD(hl).v=\n%s\n", FormatArray.toString(
                MatrixUtil.transpose(svd.vT), "%.3f"));
        System.out.printf("SVD(hl).u=\n%s\n", FormatArray.toString(
                svd.u, "%.3f"));

        // from lemma 5.18: H = HL/svd(HL).sigma[1] <-- makes the 2nd singular value == 1
        double[][] h = MatrixUtil.copy(hl);
        MatrixUtil.multiply(h, 1./svd.s[1]);

        System.out.printf("h=\n%s\n", FormatArray.toString(h, "%.3f"));

        // or create it using hth = V^T * sigma_normalized_by_snd_singular * V
        double[][] hth = MatrixUtil.createATransposedTimesA(h);
        System.out.printf("\nH^T*H=\n%s\n", FormatArray.toString(hth, "%.3f"));

        // check that the 2nd singular value of H^T*H is "1"
        MatrixUtil.SVDProducts svd2 = MatrixUtil.performSVD(hth);
        System.out.printf("SVD(H^T*H).s=\n%s\n", FormatArray.toString(svd2.s, "%.3f"));
        //System.out.printf("SVD(H^T*H).vT=\n%s\n", FormatArray.toString(svd2.vT, "%.3f"));
        //System.out.printf("SVD(H^T*H).u=\n%s\n", FormatArray.toString(svd2.u, "%.3f"));

        //multiply u and v by -1 if det(U) < 0
        double detU = MatrixUtil.determinant(svd2.u);
        System.out.printf("det U=%.3f\n", detU);
        if (detU < 0) {
            MatrixUtil.multiply(svd2.u, -1);
            MatrixUtil.multiply(svd2.vT, -1);
            double[][] chk = MatrixUtil.multiplyByDiagonal(svd2.u, svd2.s);
            chk = MatrixUtil.multiply(chk,svd2.vT);
            //System.out.printf("H^T*H=\n%s\n", FormatArray.toString(chk, "%.3f"));
        }
        System.out.printf("SVD(H^T*H).vT=\n%s\n", FormatArray.toString(svd2.vT, "%.3f"));
        System.out.printf("SVD(H^T*H).u=\n%s\n", FormatArray.toString(svd2.u, "%.3f"));

        double[] s2 = new double[3];
        for (i = 0; i < 3; ++i) {
            s2[i] = svd2.s[i];
        }
        double[][] vsvt = MatrixUtil.multiply(
                MatrixUtil.multiplyByDiagonal(MatrixUtil.transpose(svd2.vT), s2),
                svd2.vT);
        System.out.printf("V*S*V^T=\n%s\n", FormatArray.toString(vsvt, "%.3f"));

        //H = +- (R+ (1/d)T*N^T).

        // v2 is orthogonal to N and T
        //double[] v2 = Arrays.copyOf(svd2.vT[1], 3);

        //lambda2_j * x2_j = H * lambda1_j * x1_j

        /*
        (1) decompose H into R, T/d, and N
        (2) simulate 4 points from the homography
         */

        /* decompose: see the 4 solutions in Table 5.1.
          the 2 which have N^T*e3 = n3 > 0 pass the positive depth constraint

        U1 = [v2, u1, [v2]_x * u1]
        W1 = [H*v2, H*u1, [H*v2]_x * H * u1]]

        U2 = [v2, u2, [v2]_x * u2)
        W2 = [H*v2, H*u2, [H*v2]_x * H * u2]]

        Soln 1:
        R1 = W1 * U1^T
        N1 = [v2]_x * u1
        (1/d)*T1 = (H - R1)*N1

        Soln 2:
        R2 = W2 * U2^T
        N2 = [v2]_x * u2
        (1/d)*T2 = (H - R2)*N2

        Soln 3:
        R3 = R1
        N3 = -N1
        (1/d)*T3 = -(1/d)*T1

        Soln 4:
        R4 = R2
        N4 = -N2
        (1/d)*T4 = -(1/d)*T2

        the 2 which have N^T*e3 = n3 > 0 pass the positive depth constraint
        */

/*
        these are not the same solutions used in the book.  the resulting
        rotation matrices in example 5.20 have column vectors whose
        square sums = 1
                check the math below
*/
        // ==== MASKS eqn (5.44) ====

        double s1sq = svd2.s[0];
        double s2sq = svd2.s[1];
        double s3sq = svd2.s[2];
        // rows of vT are columns of v
        double[] v1 = Arrays.copyOf(svd2.vT[0], 3);
        double[] v2 = Arrays.copyOf(svd2.vT[1], 3);
        double[] v3 = Arrays.copyOf(svd2.vT[2], 3);

        // creating u1 and u2 from svd2.u columns 0 and 2.
        //   avoiding column 1 (=eigenvector 1) which has singular value=1?
        double f1 = Math.sqrt(1. - s3sq);
        double f2 = Math.sqrt(s1sq - 1.);
        double f3 = Math.sqrt(s1sq - s3sq);
        double[] u1 = new double[3];
        double[] u2 = new double[3];
        for (i = 0; i < 3; ++i) {
            u1[i] = (v1[i] * f1 + v3[i] * f2)/f3;
            u2[i] = (v1[i] * f1 - v3[i] * f2)/f3;
        }
        System.out.printf("u1=\n%s\n", FormatArray.toString(u1, "%.3f"));
        System.out.printf("u2=\n%s\n", FormatArray.toString(u2, "%.3f"));

        double[][] skewV2 = MatrixUtil.skewSymmetric(v2);
        double[][] skewHV2 = MatrixUtil.skewSymmetric(
                MatrixUtil.multiplyMatrixByColumnVector(h, v2)
        );

        //U1 = [v2, u1, skew(v2)*u1];
        double[][] U1 = new double[3][];
        U1[0] = Arrays.copyOf(v2, v2.length);
        U1[1] = Arrays.copyOf(u1, u1.length);
        U1[2] = MatrixUtil.multiplyMatrixByColumnVector(skewV2, u1);
        U1 = MatrixUtil.transpose(U1);

        //U2 = [v2, u2, skew(v2)*u2];
        double[][] U2 = new double[3][];
        U2[0] = Arrays.copyOf(v2, v2.length);
        U2[1] = Arrays.copyOf(u2, u2.length);
        U2[2] = MatrixUtil.multiplyMatrixByColumnVector(skewV2, u2);
        U2 = MatrixUtil.transpose(U2);

        //W1 = [H*v2, H*u1, skew(H*v2)*H*u1];
        double[][] W1 = new double[3][];
        W1[0] = MatrixUtil.multiplyMatrixByColumnVector(h, v2);
        W1[1] = MatrixUtil.multiplyMatrixByColumnVector(h, u1);
        W1[2] = Arrays.copyOf(W1[1], W1[1].length);
        W1[2] = MatrixUtil.multiplyMatrixByColumnVector(skewHV2, W1[2]);
        W1 = MatrixUtil.transpose(W1);

        //W2 = [H*v2, H*u2, skew(H*v2)*H*u2];
        double[][] W2 = new double[3][];
        W2[0] = MatrixUtil.multiplyMatrixByColumnVector(h, v2);
        W2[1] = MatrixUtil.multiplyMatrixByColumnVector(h, u2);
        W2[2] = Arrays.copyOf(W2[1], W2[1].length);
        W2[2] = MatrixUtil.multiplyMatrixByColumnVector(skewHV2, W2[2]);
        W2 = MatrixUtil.transpose(W2);

        //N1 = skew(v2)*u1;
        double[] n1 = MatrixUtil.multiplyMatrixByColumnVector(skewV2, u1);
        //N2 = skew(v2)*u2;
        double[] n2 = MatrixUtil.multiplyMatrixByColumnVector(skewV2, u2);

        // solutions with n1[2] > 0 or n2[2] > 0 are in front of the camera

        double[][] W1U1T = MatrixUtil.multiply(W1, MatrixUtil.transpose(U1));
        double[][] W2U2T = MatrixUtil.multiply(W2, MatrixUtil.transpose(U2));
        double[] HW1U1TN1 = MatrixUtil.multiplyMatrixByColumnVector(
                MatrixUtil.pointwiseSubtract(h, W1U1T), n1);
        double[] HW2U2TN2 = MatrixUtil.multiplyMatrixByColumnVector(
                MatrixUtil.pointwiseSubtract(h, W2U2T), n2);

        //Sol(:,:,1) = [W1*U1', (H - W1*U1')*N1, N1];
        //Sol(:,:,2) = [W2*U2', (H - W2*U2')*N2, N2];
        //Sol(:,:,3) = [W1*U1', -(H - W1*U1')*N1, -N1];
        //Sol(:,:,4) = [W2*U2', -(H - W2*U2')*N2, -N2];
        System.out.printf("\nSolution 1:\n");
        System.out.printf("R:\n%s\n", FormatArray.toString(W1U1T, "%.3f"));
        System.out.printf("N:\n%s\n", FormatArray.toString(n1, "%.3f"));
        System.out.printf("(1/d)*T:\n%s\n", FormatArray.toString(HW1U1TN1, "%.3f"));

        System.out.printf("\nSolution 2:\n");
        System.out.printf("R:\n%s\n", FormatArray.toString(W2U2T, "%.3f"));
        System.out.printf("N:\n%s\n", FormatArray.toString(n2, "%.3f"));
        System.out.printf("(1/d)*T:\n%s\n", FormatArray.toString(HW2U2TN2, "%.3f"));

        System.out.printf("\nSolution 3:\n");
        MatrixUtil.multiply(HW1U1TN1, -1);
        MatrixUtil.multiply(n1, -1);
        System.out.printf("R:\n%s\n", FormatArray.toString(W1U1T, "%.3f"));
        System.out.printf("N:\n%s\n", FormatArray.toString(
                n1, "%.3f"));
        System.out.printf("(1/d)*T:\n%s\n", FormatArray.toString(HW1U1TN1, "%.3f"));

        System.out.printf("\nSolution 4:\n");
        MatrixUtil.multiply(HW2U2TN2, -1);
        MatrixUtil.multiply(n2, -1);
        System.out.printf("R:\n%s\n", FormatArray.toString(W2U2T, "%.3f"));
        System.out.printf("N:\n%s\n", FormatArray.toString(
                n2, "%.3f"));
        System.out.printf("\n(1/d)*T:\n%s\n", FormatArray.toString(HW2U2TN2, "%.3f"));

        double[][] solnR = W2U2T;
        double[] solnN = n2;
        double[] soln1dT = HW2U2TN2;
        // solution 1 and solution 4 have n[2] > 0
        // assert the rotation of solution 1 and solution 4
        //     have R^T = R^−1 and det R = 1
        r = W2U2T;
        double tol = 1E-5;
        double detR = MatrixUtil.determinant(r);
        double[][] invR = MatrixUtil.inverse(r);
        double[][] rT = MatrixUtil.transpose(r);
        assertTrue(Math.abs(detR - 1) < tol);
        for (i = 0; i < 3; ++i) {
            for (j = 0; j < 3; ++j) {
                assertTrue(Math.abs(invR[i][j] - rT[i][j]) < tol);
            }
        }
        // assert orthogonal columns
        double c0c1Sum = 0;
        double c0c2Sum = 0;
        double c1c2Sum = 0;
        for (i = 0; i < 3; ++i) {
            c0c1Sum += r[i][0]*r[i][1];
            c0c2Sum += r[i][0]*r[i][2];
            c1c2Sum += r[i][1]*r[i][2];
        }
        assertTrue(Math.abs(c0c1Sum) < tol);
        assertTrue(Math.abs(c1c2Sum) < tol);
        assertTrue(Math.abs(c1c2Sum) < tol);

        // ====== create some points in image 1 and image 2 consistent with the solution
        // x2 ~ H * x1
        // and [x2]_x * H * x1 = 0

        // x1 = lambda * P * X where P=[I|0] and
        // x2 = lambda * P * X where P = [R | t]

        // draw X from a range of 512 x 512 in WCS
        int nCorres = 4;
        int size = 512;

        /*
        N = (0.447, -0.000, 0.894) for theta = 10 degrees
        or N = (1, 0, 2)  unnormalized
        d = 5
        (1) create a plane P and generate points on plane P
        (2) create image plane 1 in direction N from plane P and at a distance d
        (3) project plane P points to image plane 1
        (4) rotate and translate image plane 1 to create image plane 2
        (5) rotate and translate points from image 1 to the reference frame of image 2

        (1)
        N = n[0]i + n[1]j + n[2]k
        d = n[0]*xw_1[0] + n[1]*xw_1[1] + n[2]*xw_1[2]
        Plane passing thru xw_1 perpendicular to vector N is
           n[0]*(x - xw_1[0]) + n[1]*(y - xw_1[1]) + n[2]*(z - xw_1[2]) = 0
           OR
           n[0]*x + n[1]*y + n[2]*z = d

        generate 1 point
            xw_1[0] = rand.next(double)
            xw_1[1] = rand.next(double)
            xw_1[2] = rand.next(double)

        then using n, xw_1 and the equation for the plane:
            n[0]*(x - xw_1[0]) + n[1]*(y - xw_1[1]) + n[2]*(z - xw_1[2]) = 0
        repeat this for 3 more points, while avoiding colinear sets of any 3 points.
        x = rand.next(double)
        y = rand.next(double)
        n[2]*(z - xw_1[2]) = -n[0]*(x - xw_1[0]) - n[1]*(y - xw_1[1])
        (z - xw_1[2]) = (-n[0]*(x - xw_1[0]) - n[1]*(y - xw_1[1]))/n[2]
        z = xw_1[2] + ( (-n[0]*(x - xw_1[0]) - n[1]*(y - xw_1[1]))/n[2])

        assert:
           n[0]*(x - xw_1[0]) + n[1]*(y - xw_1[1]) + n[2]*(z - xw_1[2]) = 0

        */

        //TODO:
        //put these into a method for extracting motion and structure from 2 views from 2 calibrated images

        // the plane P is tangent to n
        // which is (-2, 0, 1) or any multiple of that.

        /*
        N = n[0]i + n[1]j + n[2]k
        d = n[0]*xw_1[0] + n[1]*xw_1[1] + n[2]*xw_1[2]
        Plane passing thru xw_1 perpendicular to vector N is
           n[0]*(x - xw_1[0]) + n[1]*(y - xw_1[1]) + n[2]*(z - xw_1[2]) = 0
           OR
           n[0]*x + n[1]*y + n[2]*z = d
           generate the first point in the plane as (0,0,0)
           then for n[0]*(x) + 0 + n[2]*(z) = 0
              generate x, and set y = 0
              calc z = -(n[0]*x)/n[2]
       */

        // generate points in plane P in WCS
        double[][] XW = new double[4][nCorres];
        for (i = 0; i < 4; ++i) {
            XW[i] = new double[nCorres];
        }
        Arrays.fill(XW[3], 1);
        XW[0][0] = 0;
        XW[1][0] = 0;
        XW[2][0] = 0;
        boolean isCollinear = false;
        do {
            // y is 0 for this example's planar homography
            System.arraycopy(
                    UnivariateNormalDistribution.randomSampleOfUnitStandard(rand, nCorres-1),
                    0, XW[0], 1, nCorres - 1);
            for (i = 1;  i < nCorres; ++i) {
                XW[2][i] = -(solnN[0]*XW[0][i])/solnN[2];
            }
            isCollinear = areCollinear(nCorres, XW, tol);
        } while (isCollinear);

        System.out.printf("XW=\n%s\n", FormatArray.toString(XW, "%.3f"));

        // adapted from MASKS book test_5_2.m

        // x1 = [3X4 projection for plane P to camera 1] * XW

        // n is the direction perpendicular to the plane P.
        // rotate the plane by n
        double[] solnNNeg = Arrays.copyOf(solnN, solnN.length);
        MatrixUtil.multiply(solnNNeg, -1);

        /*
        We need the rotation that turns direction vector (1, 0, 2) to the direction vector (0, 0, 1) so that
        it is pointing at camera 1.
        We can calculate the axis of rotation as the cross product of the
        vectors and the angle of rotation as the dot product of the vectors,
        then form the Barfoot quaternion from the angle and axis,
        then the rotation matrix from the quaternion.
        */
        // normalize the directions to make the math simpler.
        double[] nPNorm = MatrixUtil.normalizeL2(n);
        double[] nC1Norm = new double[]{0, 0, 1};
        double thetaWCSC1 = -Math.acos(nPNorm[0]*nC1Norm[0] + nPNorm[1]*nC1Norm[1] + nPNorm[2]*nC1Norm[2]);
        double[] axisWCSC1 = MatrixUtil.crossProduct(nPNorm, nC1Norm);
        double[] quatH = Rotation.createHamiltonQuaternion(thetaWCSC1, axisWCSC1);
        double[] quatB = Rotation.convertHamiltonToBarfootQuaternion(quatH);
        double[][] rWCSC1 = Rotation.createRotationMatrixFromQuaternion4(quatB);
        rWCSC1 = MatrixUtil.copySubMatrix(rWCSC1, 0, 2, 0, 2);
        System.out.printf("r*v=\n%s\n", FormatArray.toString(
                MatrixUtil.multiplyMatrixByColumnVector(rWCSC1, nPNorm),
                "%.4e"
        ));
            // expecting (0, 0, 1)
            double[] checkNC1 = MatrixUtil.multiplyMatrixByColumnVector(rWCSC1, nPNorm);
            assertTrue(Math.abs(checkNC1[0] - nC1Norm[0]) < 1E-7);
            assertTrue(Math.abs(checkNC1[1] - nC1Norm[1]) < 1E-7);
            assertTrue(Math.abs(checkNC1[2] - nC1Norm[2]) < 1E-7);

            boolean passive = false;

            // using the rodrigues formula was faster but a little less accurate:
        double[][] r2 = Rotation.createRodriguesFormulaRotationMatrix(axisWCSC1, passive);
        System.out.printf("r*v=\n%s\n", FormatArray.toString(
                MatrixUtil.multiplyMatrixByColumnVector(r2, nPNorm),
                "%.4e"
        ));

        // d is 5
        // |x1_i - xw_i| = 5  = sqrt of sq sums of translation components.
        // will translate along z only:
        System.out.println("rotation from plane P to camera 1 plane is n.  translating by z only");
        double[][] pWCSC1 = new double[3][];
        pWCSC1[0] = new double[]{rWCSC1[0][0], rWCSC1[0][1], rWCSC1[0][2], 0};
        pWCSC1[1] = new double[]{rWCSC1[1][0], rWCSC1[1][1], rWCSC1[1][2], 0};
        pWCSC1[2] = new double[]{rWCSC1[2][0], rWCSC1[2][1], rWCSC1[2][2], d};

        // [3X4] * [4XnCorres] = [3XnCorres]
        double[][] xc1 = MatrixUtil.multiply(pWCSC1, XW);

        // check the distances
        double chk;
        double dX, dY, dZ;
        for (i = 0; i < nCorres; ++i) {
            dX = xc1[0][i] - XW[0][i];
            dY = xc1[1][i] - XW[1][i];
            dZ = xc1[2][i] - XW[2][i];
            chk = Math.sqrt(dX*dX + dY*dY + dZ*dZ);
            System.out.printf("plane P to camera 1 d=%.4e\n", chk);
        }

        // to get the normal to xc1, can solve for the right null space of A*x=0 where
        // A is the stack of rows of xc points.  the right null space is the
        // last column of SVD(A).V where the last column is for eigenvalue=0.
        {
            double[][] stack1 = new double [nCorres][];
            for (i = 0; i < nCorres; ++i) {
                stack1[i] = new double[]{xc1[0][i], xc1[1][i], xc1[2][i], 1};
                //stack1[i] = new double[]{XW[0][i], XW[1][i], XW[2][i], XW[3][i]};
            }
            int idx = stack1[0].length;
            MatrixUtil.SVDProducts svd4 = MatrixUtil.performSVD(stack1);
            System.out.printf("the normal to camera 1 plane =%s\n",
                FormatArray.toString(svd4.vT[idx-1],"%.4e"));

            System.out.printf("the distance to camera 1 plane =%.4e\n", svd4.vT[idx-1][idx-1]);
            System.out.printf("(%.4e, %.4e, %.4e)\n", n[0]/svd4.vT[idx-1][0], n[1]/svd4.vT[idx-1][1],
                    n[2]/svd4.vT[idx-1][2]);
        }

        // xc2 = R*XC1 + t
        double[][] xc2 = MatrixUtil.multiply(r, xc1);
        for (i = 0; i < nCorres; ++i) {
            for (j = 0; j < 3; ++j) {
                xc2[j][i] += t[j];
            }
        }

      //  assert x0 = (h00x+h01y+h02)/(h20x + h21y + h22)
      //       and y0 = (h10x+h11y+h12)/(h20x + h21y + h22)

        // given R and t, can construct essential matrix E = [t]_x*R
        //    and then

        //lambda2 * x2 = H * lambda1 * x1

        //solnN is the unit normal vector of the plane P with respect to the first camera frame.

        // section 2.4.1 is homogeneous coordinate system
        // homogeneous coordinates append a "1" to (x,y,z) coordinates,
        // resulting in a 4-dimensional vactor (x,y,z,1).
        // This embeds Euclidean E^3 into a hyperplane in R^4 instead of R^3 (MASKS, "An Invitation to 3-D Vision")

        System.out.printf("XC1=\n%s\n", FormatArray.toString(xc1, "%.3f"));
        System.out.printf("XC2=\n%s\n", FormatArray.toString(xc2, "%.3f"));

        // with camera intrinsic matrices, we can create the image coordinates.
        double[][] intr = MatrixUtil.zeros(3, 3);
        intr[0][0] = 1200; //1200 <== these unused are from test_5_2.m
        intr[1][1] = 900; //900
        intr[0][2] = 600; //600
        intr[1][2] = 450; //450
        intr[2][2] = 1;

        // no radial distortions in these cameras
        double[][] x1 = new double[3][nCorres];
        double[][] x2 = new double[3][nCorres];
        for (i = 0; i < nCorres; ++i) {
            // normalized pinhole projection X_c/Z_c
            x1[0][i] = xc1[0][i]/xc1[2][i];
            x1[1][i] = xc1[1][i]/xc1[2][i];
            x1[2][i] = 1;
            x2[0][i] = xc2[0][i]/xc2[2][i];
            x2[1][i] = xc2[1][i]/xc2[2][i];
            x2[2][i] = 1;
        }
        // convert to image frame:
        x1 = MatrixUtil.multiply(intr, x1); // 3 x nCorres
        x2 = MatrixUtil.multiply(intr, x2);

        boolean useNormConditioning = true;

        double[][] hEst3 = CameraCalibration.solveForHomography(x1, x2, useNormConditioning);

        MatrixUtil.SVDProducts svd5 = MatrixUtil.performSVD(hEst3);

        double[][] h3 = MatrixUtil.copy(hEst3);
        MatrixUtil.multiply(h3, 1./svd5.s[1]);

        // see algorithm 5.2

        System.out.printf("h3=\n%s\n", FormatArray.toString(h3, "%.3e"));

        // ===========================
        // another approach:
        //     randomly generate xc1 in camera frame,
        //     use xc2 ~ H * xc1 <--> also x2 ~ H * x1
        //     and then triangulation to derive XW.
        //     the cameras' intrinsics are known.

    }

    private void createSkewMatrix(double[][] x2, int col, double[][] x2j) {
        x2j[0][0] = 0;
        x2j[1][1] = 0;
        x2j[2][2] = 0;
        x2j[0][1] = -x2[2][col];
        x2j[0][2] = x2[1][col];
        x2j[1][0] = x2[2][col];
        x2j[1][2] = -x2[0][col];
        x2j[2][0] = -x2[1][col];
        x2j[2][1] = x2[0][col];
    }

    private void extractColumn(double[][] x, int col, double[][] xj) {
        xj[0][0] = x[0][col];
        xj[1][0] = x[1][col];
        xj[2][0] = x[2][col];
    }

    private boolean areCollinear(int ns, double[][] X, double tol) {
        int k = 3;
        int[] selectedIndexes = new int[k];
        long nComb = MiscMath0.computeNDivKTimesNMinusK(ns, k);
        SubsetChooser chooser = new SubsetChooser(ns, k);
        int c = 0;
        while (chooser.getNextSubset(selectedIndexes) != -1) {
            double x1 = X[0][selectedIndexes[0]];
            double y1 = X[1][selectedIndexes[0]];
            double x2 = X[0][selectedIndexes[1]];
            double y2 = X[1][selectedIndexes[1]];
            double x3 = X[0][selectedIndexes[2]];
            double y3 = X[1][selectedIndexes[2]];
            if (areCollinear(x1, y1, x2, y2, x3, y3, tol)) {
                return true;
            }
            c++;
        }
        return false;
    }

    private boolean areCollinear(double tol, double x1, double y1, double x2, double y2,
                                double x3, double y3) {

        double direction = ((x2 - x1)*(y3 - y1)) - ((y2 - y1)*(x3 - x1));

        return Math.abs(direction) < tol;
    }
}
