package algorithms.imageProcessing.transform;

import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.features.RANSACSolver;
import algorithms.imageProcessing.matching.ErrorType;
import algorithms.imageProcessing.transform.Reconstruction.ReconstructionResults;
import algorithms.matrix.MatrixUtil;
import algorithms.matrix.MatrixUtil.SVDProducts;
import algorithms.misc.MiscMath0;
import algorithms.util.FormatArray;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;
import java.text.Normalizer;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.NotConvergedException;
import no.uib.cipr.matrix.SVD;

/**
 * given correspondences, pairs of points for the same objected projected into 2 images,
 * methods in this class calculate transformations necessary to make two images
 * parallel to the baseline between the camera optical centers
 * and at the focal distance.  Rectification makes the epipolar lines parallel
 * and the epipoles at infinity.
 * 
   NOTE: if the epipoles are within the images, consider using Pollefeys, 2000
   (not implemented here at this time).

   NOTE: alternatively, if have at least 7 points to create an epipolar mapping, can skip
   rectification and use stereo matching:
   "Stereo Processing by Semi-Global Matching and Mutual Information" by Hirschmuller 2008
   
 * @author nichole
 */
public class Rectification {

    // notes from Serge Belongie lectures from Computer Vision II, CSE 252B, USSD.
    // other references from Kris Kitani's lectures in 16-385 Computer Vision,
    // Carnegie Mellon University,
    // and Ma, Soatto, Kosecka,& Sastry 2012 "Invitation to Computer Vision, From Images to Geometric Models",
    //
    // homogeneous repr of a point is x_vec = (x, y, 1)^T
    // equation f a line is ell = a*x + b*y + c = 0;
    // line rewritten in homogeneous coordinatrs is x_vec^T * ell.
    // general conic in 3 dimensions is a*x^2 + b*x*y + c*y^2 + d*x*z + e*y*z + f*z^2 = 0.
    //     rewritten using 2D homogenouse coords, quadratice form: x_vec^T * C * x_vec = 0
    //                  [a   b/2   d/2 ]
    //        where C = [b/2   c   c/2 ]
    //                  [d/2 c/2     f ]
    //        C has 6 variable, 5 DOF, so need 5 points
    //     can then reformat x_vec^T * C * x_vec = 0 into
    //        the "design matrix" * "the carrier vector" = 0
    //               A * c = 0
    //        c = SVD(A).V^T[n-1], the eigenvector assoc w/ smallest eigenvalue
    //
    // lecture 5.5 Epipolar Rectification:
    //    given 2 views of a scene, the goal is to apply projective transformations
    //    to the images so that all epipolar lines correspond to the horizontal 
    //    scan lines.
    //    so need to find 2 linear transformations, H1 and H2 that map the
    //    epipoles to infinity along x-axis (coord [1, 0,0]^T).
    //
    // (1) compute E (or F) and e2.
    // (2) map e2 to infinity to make the epipolar lines parallel using H2
    //     (Hartley 1997).
    //     There is a family of H2's that will do this, parameterized by
    //     v (where v is real matrix of dimension 3x3)
    //        H = ([T]_x)^T * E + T * e^T
    //      where T is the translation vector between the 2 cameras.
    //
    //      find the H2 such that 
    //          H2*e2 ~ [1, 0, 0]^T
    //      where H2 is as close as possible to a rigid body transformattion
    //      
    //      define the translation of the image center to the origin:
    //                 [ 1  0  -o_x ]
    //           G_T = [ 0  1  -o_y ]
    //                 [ 0  0    1  ]
    //
    //      declare the rotation about the Z-axis to put the translated
    //      epipole onto the x-axis:
    //          G_R * G_T * e2 = [e_x_coord  0  1]^T
    //      
    //      define matrix G to "send the epipole to infinity":
    //              [ 1            0   0 ]
    //          G = [ 0            1   0 ]
    //              [ 1/e_x_coord  0   1 ]
    //
    //      therefore the rectification for the 2nd view is
    //          H2 = G * G_R * G_T and is a real 3X3 matrix
    //
    // (3)  To find an H compatible with E (or F), use the method of section
    //      5.4 for finding H from E.
    //      use the least squares version of 
    //          H = ([T]_x)^T * E + T * e^T
    //      for multiple points and choose the H that minimizes the
    //      distortion induced by the rectification transformation.
    //  add details here
    //
    // (4) compute the matching homography  
    //        H2 = H1 * H
    // 
    // (5) apply H1 and H2 to the left and right image respectively
    //
    // NOTE that if the camera is moving toward the image, the epipole
    // is inside the image and one must use another method.
    // (see Pollefeys)
    //
    //
    //
    // Hartley 1999, "Theory and Practice of Projective Rectification"
    // http://www.cs.ait.ac.th/~mdailey/cvreadings/Hartley-Rectify.pdf
    //
    // Mallon & Whelan 2005, "Projective Rectification from the Fundamental Matrix"
    // http://doras.dcu.ie/4662/1/JM_IVC_2005.pdf
    // http://www.cipa.dcu.ie/papers/ivc_2005_jm.pdf
    //
    // Monasse, Morel, and Tang 2011
    // "Three-step image rectification"
    // https://core.ac.uk/download/pdf/48342838.pdf
    //
    
    /**
     * NOT READY FOR USE
     * Rectify (i.e, warp) the left image correspondence points x1 and
     * right image correspondence points x2
     * so that corresponding horizontal scanlines are epipolar lines.
     
     * NOTE that the method rectify() is preferred.
     * 
     * <pre>
     * references:
     * following the algorithm of Kitani lecture 13.1, class 16-385 Computer Vision,
         Carnegie Mellon University.
       
       there are many alternative approaches depending on datasets mentioned in
       Seliski 2010, "Computer Vision: Algorithms and Applications"
    
       Ma, Soatto, Kosecka,& Sastry 2012 "Invitation to Computer Vision, 
           From Images to Geometric Models",
       Chapter 11, Section 11.5.1
     * </pre>
     * @param x1 the image 1 set of correspondence points.  format is 3 x N where
     * N is the number of points.
     * @param x2 the image 2 set of correspondence points.  format is 3 x N where
     * N is the number of points.
     * @param k1Intr intrinsic parameters for camera 1
     * @param k2Intr intrinsic parameters for camera 2
     * @return 
     */
    private static RectifiedPoints epipolar(double[][] k1Intr,
        double[][] k2Intr, double[][] x1, double[][] x2) throws NotConvergedException {
        
        /*
        from Kitani lecture:
        
       when rectified, the images are parallel, i.e. epipolar lines are horizontal:
      R = 1 0 0   t = T,0,0  [t]_x = 0 0 0
          0 1 0                      0 0 -T
          0 0 1                      0 T 0

      E = [t]_x*R = 0  0  0
                    0  0 -T
                    0  T  0

      and x^T*E*x' = 0 has to remain true

      [u v 1] [0  0  0] [u'] = 0
              [0  0 -T] [v']
              [0  T  0] [1 ]

          [u v 1] [0   ] = 0
                  [-T  ]
                  [v'*T]

   0 + -v*T + v'*T = 0
   v*T = v'*T;  y coord is always the same

   to reproject image planes onto a common plane parallel to the line
   between camera centers, need a homography (3X3 tranform) for each
   image.

   1. Rotate the right camera by R  
      (aligns camera coordinate system orientation only)
      1a. Compute E to get R
          (use reconstruction w/ intrinsic if have it, or without if not)
          Reconstruction.calculateUsingEssentialMatrix
          Reconstruction.calculateProjectiveReconstruction
       and the epipoles e1, e2
         e1 is the last row of svd.vt
            e1 is the right nullspace (in right singular vector) of F
               (F*e1 = 0 or E*e1 = 0)
        e2 is the last column of svd.u
            e2 is the left nullspace (in left singular vector) of F
                (e2^T*F = 0  or e2^T*E = 0)
            e2 when normalized by 3rd coord is in coord space of left image and
               it is the location of the right camera center.
            NOTE: translation is also the last col of svd.u
       epipolar lines l1_i and l2_i
        l2 = E*x1
        l1 = E^T*x2
        */
        
       ReconstructionResults re = Reconstruction.calculateUsingEssentialMatrix(k1Intr, k2Intr, x1, x2);
       
       double[] t = Arrays.copyOf(re.k2ExtrTrans, re.k2ExtrTrans.length);
       double[][] r = MatrixUtil.copy(re.k2ExtrRot);

       // right nullspace of F:
       double[] e1 = Arrays.copyOf(re.svdVt[2], re.svdVt[2].length);
       // left nullspace of F:
       double[] e2 = MatrixUtil.transpose(re.svdU)[2];
       //MatrixUtil.multiply(e1, 1./e1[2]);
       //MatrixUtil.multiply(e2, 1./e2[2]);
       e1 =  MatrixUtil.normalizeL2(e1);
       e2 =  MatrixUtil.normalizeL2(e2);
       
       /*
       MASKS Proposition 5.3: 
            e2^T*E = 0, E*e1 = 0.
            e2 ~ T and e1 ~ R^T * T where ~ is up to a scale factor
       
       let r_1 = e2 = T/||T|| so that epipole coincides w/ translation vector
       */
       double[] r1 = Arrays.copyOf(t, t.length);
       r1 = MatrixUtil.normalizeL2(r1);
       
       double[][] tSkewSym = MatrixUtil.skewSymmetric(t);
       double[][] rtSkewSym = MatrixUtil.multiply(r, tSkewSym);
       
       System.out.printf("t=%s\ne1=%s\ne2=%s\nr1=%s\nessentialMatrix=\n%s\n(R*(t_skewsym))=\n%s\n", 
           FormatArray.toString(t, "%.4e"),
           FormatArray.toString(e1, "%.4e"),
           FormatArray.toString(e2, "%.4e"),
           FormatArray.toString(r1, "%.4e"),
           FormatArray.toString(re.essentialMatrix, "%.4e"),
           FormatArray.toString(rtSkewSym, "%.4e"));
       
        /*
        forming rRect to transform the epipole e2 to [1,0,0]^T.
          rRect * e2 = [1,0,0]
       
          direction vector of optical axis = [0, 0, -1], pointing towards positive z
       
          r1 = rRect[0] = e2
          r2 = rRect[1] = cross product of e and the direction vector of the optical axis
          r3 = rRect[2] = r1 cross r2.
       
         cross proeduct of e2 and optical axis direction [0,0,1]
           (T_y*-1 - T_z*0)/||T|| = -T_y/||T||
           (T_z*0 - T_x*-1)/||T|| = T_x/||T||
           (T_x*0 - T_y*0)/||T||; = 0
         except for dropping the T_z^2 term from ||T|| in the normalization.
       
         r_2 = [-T_y  T_x  0]/sqrt(T_x^2 + T_y^2)
        */
        
        double td = 1./Math.sqrt(t[0]*t[0] + t[1]*t[1]);
        double[] r2 = new double[]{-t[1] * td, t[0] * td, 0};

        //let r_3 = r1Xr2 orthogonal vector
        double[] r3 = MatrixUtil.crossProduct(r1, r2);

        double[][] rRect = MatrixUtil.zeros(3, 3);
        System.arraycopy(r1, 0, rRect[0], 0, r1.length);
        System.arraycopy(r2, 0, rRect[1], 0, r2.length);
        System.arraycopy(r3, 0, rRect[2], 0, r3.length);

        System.out.printf("rRect=%s\n", FormatArray.toString(rRect, "%.4e"));

        double[][] r2Rot = MatrixUtil.copy(rRect);
        double[][] r1Rot = MatrixUtil.multiply(r, rRect);

        // MASKS eqn (11.28) where H2*e2 is r2Rot*r1 here.  assert = [1,0,0]^T.
        double[] tst = MatrixUtil.multiplyMatrixByColumnVector(r1Rot, e1);
        System.out.printf("r2Rot*e1=%s\nexpecting=[1, 0, 0]\n",
                FormatArray.toString(tst, "%.4e"));
        
        tst = MatrixUtil.multiplyMatrixByColumnVector(r2Rot, e2);
        System.out.printf("r2Rot*e2=%s\nexpecting=[1, 0, 0]\n\n",
                FormatArray.toString(tst, "%.4e"));

        /*
       2. Rotate (rectify) the left camera so that the epipole is at infinity
          [x2 y2 z2] = R1 * [x1 y1 z1] = warped left which should 
           equal [x2 y2 z2] with caveat
                             due to occlusion, etc.

       points p = (f/z2)*[x2 y2 z2]
       if have intrinsic parameters matrix K then
           points p ~ K*R1*[x1 y1 z1]
             *Kitani notes that you may need to alter f inside K to keep
              points within the original image size
              (details are in Ma et al "An Invitiation to #-D..."
               pg 400, Chapt 11)
        
       f=(W/2)*((tan(fov/2))^-1)
         */
        k1Intr = MatrixUtil.copy(k1Intr);
        k2Intr = MatrixUtil.copy(k2Intr);
        k1Intr[0][0] *= -1; 
        k1Intr[1][1] *= -1;
        k2Intr[0][0] *= -1; 
        k2Intr[1][1] *= -1;
        
        System.out.printf("R1=rRect=%s\n",
                FormatArray.toString(r1Rot, "%.4e"));
        System.out.printf("R2=R*rRect=%s\n",
                FormatArray.toString(r2Rot, "%.4e"));
        
        double[][] _h1 = MatrixUtil.multiply(k1Intr, r1Rot);
        double[][] _h2 = MatrixUtil.multiply(k2Intr, r2Rot);
        //_h1 = r1Rot;
        //_h2 = r2Rot;

        // x1 is left image points
        double[][] x1R = MatrixUtil.multiply(_h1, x1);
        // x2 is right image points
        double[][] x2R = MatrixUtil.multiply(_h2, x2);

        int i, j;
        int n = x1R[0].length;
        
        // normalize z-coords to be 1        
        for (i = 0; i < n; ++i) {
            for (j = 0; j < 3; ++j) {
                x1R[j][i] /= x1R[2][i];
                x2R[j][i] /= x2R[2][i];
            }
        }
        
        RectifiedPoints rPts = new RectifiedPoints();
        rPts.setX1(x1R);
        rPts.setX2(x2R);
        rPts.setH1(_h1);
        rPts.setH2(_h2);
        
        return rPts;
    }
    
    /**
     *
     * This is for un-calibrated cameras.  If the images have a large range of
     * depth in them or if the epipoles are inside the images, 
     * this algorithm can result in distortions.
     * If one has camera intrinsic and extrinsic parameters, Ma et al. suggest
     * use method of Fusiello et al. 1997 for Euclidean projection. 
     <pre>
     following algorithm 11.9 of Ma, Soatto, Kosecka,& Sastry (MASKS)
     * "An Invitation to Computer Vision, From Images to Geometric Models".
     also present in their code projRectify.m
     </pre>
     Note that the method is meant for use when the epipoles are outside of the image.
     In this case, one could use a nonlinear polar rectification
     as suggested by Pollefeys 2000
     ("3D model from Images" or "A simple and efficient rectification method for general motion", Pollefys, Koch, and Van Gool)
     * @param fm the fundamental matrix between image 1 and image 2.
     * @param x1 the image 1 set of correspondence points. format is 3 x N
     * where N is the number of points. NOTE: since intrinsic parameters are not
     * known, users of this method should presumably center the coordinates in
     * some manner (e.g. subtract the image center or centroid of points) since
     * internally an identity matrix is used for K.
     * @param x2 the image 2 set of correspondence points. format is 3 x N where
     * N is the number of points. NOTE: since intrinsic parameters are not
     * known, users of this method should presumably center the coordinates in
     * some manner (e.g. subtract the image center or centroid of points).
     * @param oX image center along x-axis in pixels (usually width/2).
     * @param oY image center along y-xaxis in pixels (usually height/2).
     * @return rectified points and the homography matrices used to transform them.
     */
    public static RectifiedPoints rectify(double[][] fm, double[][] x1, double[][] x2,
        double oX, double oY) throws NoSuchAlgorithmException, NotConvergedException {

        if (x1.length != 3 || x2.length != 3) {
            throw new IllegalArgumentException("x1.length must be 3 and so must x2.length");
        }
        int n0 = x1[0].length;
        if (x2[0].length != n0) {
            throw new IllegalArgumentException("x1 and x2 must be same dimensions");
        }
        if (n0 < 4) {
            throw new IllegalArgumentException("need at least 4 points for the planar homography");
        }
        if (oX < 1 || oY < 1) {
            throw new IllegalArgumentException("oX and oY must be positive integers");
        }

        // first find a transfor motion H2 that maps the second epipole e2 to infinity and aligns the epipolar lines
        // with the scanlines as H2.
        // H2*e2 ~ [1,0,0)^T, but that leaves 6 degrees of freedom in H2, so choices are made below
        // to keep H2 as close as possible to a rigid body transformation.

        int n = x1[0].length;

        DenseMatrix fmM = new DenseMatrix(fm);

        SVD svd = SVD.factorize(fmM);

        double[][] uF = MatrixUtil.convertToRowMajor(svd.getU());
        //double[][] vTF = MatrixUtil.convertToRowMajor(svd.getVt());
        double[] s = Arrays.copyOf(svd.getS(), svd.getS().length);

        double[] ep2 = MatrixUtil.extractColumn(uF, 2);
        MatrixUtil.multiply(ep2, 1./MatrixUtil.lPSum(ep2, 2));
        //System.out.printf("ep2=%s\n", FormatArray.toString(ep2, "%.3e"));
        //double[] ep1 = Arrays.copyOf(vTF[2], vTF[2].length);

        double[] ep2im = Arrays.copyOf(ep2, ep2.length);
        MatrixUtil.multiply(ep2im, 1./ep2im[2]); // redoing the normalization that MTJ toolkit already applied to u, undone in ep2

        double[] vRand = Arrays.copyOf(ep2, ep2.length);
        Random rand = new Random(System.currentTimeMillis());
        for (int i = 0; i < vRand.length; ++i) {
            vRand[i] *= rand.nextDouble();
        }

        double[][] M = MatrixUtil.multiply(MatrixUtil.transpose(MatrixUtil.skewSymmetric(ep2)), fm);
        M = MatrixUtil.pointwiseAdd(M, MatrixUtil.outerProduct(ep2, vRand));

        // % take the epipole in the second view.
        // this translates the image center (oX, oY, 1) to the origin (0, 0, 1).
        // this is G_T in Sect 11.5 of MASKS.
        double[][] Tr = new double[3][];
        Tr[0] = new double[]{1, 0, -oX};
        Tr[1] = new double[]{0, 1, -oY};
        Tr[2] = new double[]{0, 0, 1};

        // translates the normalized epipole2 to new ref frame with center (0, 0, 1)
        double[] p2T = MatrixUtil.multiplyMatrixByColumnVector(Tr, ep2im);
        //p2T = [a,  b,  1]

        //  % rotate the epipole to lie on the x-axis
        // Rr is a rotation around the z-axis that rotates the translated epipole onto the x-axis;
        // i.e. G_R*G_T*e2 = [x_e, 0, 1]^T where x_e is the x coordinate of ep2.
        //This is G_r in Sect 11.5 of MASKS.
        double theta = Math.atan(-p2T[1]/p2T[0]);
        double[][] Rr = new double[3][];
        Rr[0] = new double[]{Math.cos(theta), -Math.sin(theta), 0};
        Rr[1] = new double[]{Math.sin(theta), Math.cos(theta), 0};
        Rr[2] = new double[]{0, 0, 1};

        // rotates the normalized epipole2 so that its vector w.r.t. origin [0,0,1] is parallel to the x-axis
        double[] p2R = MatrixUtil.multiplyMatrixByColumnVector(Rr, p2T);
        // p2R = [c, 0,  1]

        // G matrix transforms the epipole2 from the x-axis in the image plane to infinity [1,0,0]^T,
        // i.e. G*[x_e, 0, 1]^T ~ [x_e, 0, -1+1] ~ [1, 0, 0]^T
        double[][] G = new double[3][];
        G[0] = new double[]{1, 0, 0};
        G[1] = new double[]{0, 1, 0};
        G[2] = new double[]{-1/p2R[0], 0, 1};

        double[] pim2r = MatrixUtil.multiplyMatrixByColumnVector(G, p2R);
        // pim2R = [-d, -0,  0]

        //% rectifying transformation for the second image
        double[][] H2 = MatrixUtil.multiply(MatrixUtil.multiply(G, Rr), Tr);

        // solve for a corresponding transformation H1 for the first view,
        // called the matching homography, H1, obtained via the fundamental matrix F.

        //% one method - compute matching homography - solve for unknown plane v so as
        //  % to minimize the disparity

        //H1 = H2*H, where H can be any homography which is compatible with the fundamental matrix F,
        // i.e. [e2]_x * H ~ F.
        // Given the two conditions H2*e2 ~ [1,0,0]^T and H1 = H2*H,

        // NOTE that the choice in H is not unique

        // the three-parameter family of homographies H = (([e2]_x)^T)*F + e2*vT
        // compatible with the fundamental matrix F,
        // since v as a member of Real^3 can be arbitrary.
        // H has to be chosen with care.
        // A common choice is to set the free parameters v € R” in such a way that the distance between
        // x2, and H*x1 for previously matched feature points is minimized.
        // ~ the algebraic error associated with the homography transfer equation (11.30)
        // [x2]_x * H * x1 = [x2]_x * ( (([e2]_x)^T)*F + e2*vT ) * x1 ~ 0
        // then can solve for v using least squares with the objective
        // function (11.31): min_v( summation_{j=1,n} ( || [x2_j]_x * ( (([e2]_x)^T)*F + e2*vT ) * x1_j ||^2  )

        /* can solve by x = A^T * b where A^T is the pseudoinverse of A for full columnrank if have exact data (no errors),
         else A*x - b = 0 least squares fit using the right nullspace of SVD(A2)
          where A2 is A with a concatenated column of -1*b.

        a00 a01  * x0 - b0 =  a00*x0 + a01*x1  - b0
        a10 a11    x1   b1    a10*x0 + a11*x1  - b1

        a00 a01 -b0  * x0
        a10 a11 -b1    x1
                       1
         */

        int i;
        double t1, t2, t3;
        t1 = ep2[0];
        t2 = ep2[1];
        t3 = ep2[2];
        double[][] A = MatrixUtil.zeros(2*n, 3);
        double[] b = new double[2*n];
        double[] row1 = new double[3];
        double[] row2 = new double[3];
        for (i = 0; i < n; ++i) {
            row1[0] = -t2 * x1[0][i] + t3 * x1[0][i] * x2[1][i];
            row1[1] = -t2*x1[1][i] + t3*x1[1][i]*x2[1][i];
            row1[2] = -t2 + t3*x2[1][i];

            row2[0] = t1*x1[0][i] - t3*x1[0][i]*x2[0][i];
            row2[1] = t1*x1[1][i] - t3*x1[1][i]*x2[0][i];
            row2[2] = t1 - t3*x2[0][i];
            System.arraycopy(row1, 0, A[i*2], 0, row1.length);
            System.arraycopy(row2, 0, A[i*2 + 1], 0, row2.length);

            b[i*2] = M[1][0]*x1[0][i]+ M[1][1]*x1[1][i]+M[1][2]-M[2][0]*x1[0][i]*x2[1][i]
                -M[2][1]*x1[1][i]*x2[1][i]- M[2][2]*x2[1][i];

            b[i*2 + 1] = -M[0][0]*x1[0][i]-M[0][1]*x1[1][i]-M[0][2]+M[2][0]*x1[0][i]*x2[0][i]
                +M[2][1]*x1[1][i]*x2[0][i]+ M[2][2]*x2[0][i];
        }

        // [3 X 2*n]
        double[][] pAInv = MatrixUtil.pseudoinverseFullColumnRank(A);
        // 3 X 2*n] * [2*n X 1] = [3 X 1]
        double[] aa = MatrixUtil.multiplyMatrixByColumnVector(pAInv, b);

        System.out.printf("aa=\n%s\n", FormatArray.toString(aa, "%.3e"));

        double[][] ep2AAt = MatrixUtil.outerProduct(ep2, aa);

        double[][] H = MatrixUtil.pointwiseAdd(M, ep2AAt);
        double[][] H1 = MatrixUtil.multiply(H2, H);

        double[][] _h1 = MatrixUtil.multiply(MatrixUtil.inverse(Tr), H1);
        double[][] _h2 = MatrixUtil.multiply(MatrixUtil.inverse(Tr), H2);

        double[][] xim1r = MatrixUtil.multiply(_h1, x1);
        double[][] xim2r = MatrixUtil.multiply(_h2, x2);
        int j;
        for (i = 0; i < n; ++i) {
            for (j = 0; j < 3; ++j) {
                xim1r[j][i] /= xim1r[2][i];
                xim2r[j][i] /= xim2r[2][i];
            }
        }

        RectifiedPoints out = new RectifiedPoints();
        out.setX1(xim1r);
        out.setX2(xim2r);
        out.setH1(_h1);
        out.setH2(_h2);

        return out;

        /*
        keeping notes here from another approach to rectification.
        not sure what the reference was.  might have been from the CMU Kitani lectures or Fusiello tutorial
        or Hartley et al.

        having neither intrinsic nor extrinsic camera parameters:
        
        P1 = K_0*[I|0]
        P2 = K_2*[R|t]
        and scalar_2*x2 = epipole_2 + scalar_1*P1[subset 3x3]*(P0[subset 3x3])^-1 * x1
                        = epipole_2 + scalar_1*K2*R*(K1^-1)*x1
        and epipole_2=K2*t
        
        which can be written in homogenous coordinates:
           x2^T * [epipole_2]_x*K2*R*(K1^-1) * x1 = 0
        and 
           F = [epipole_2]_x*K2*R*(K1^-1)  ==> x2^T * F * x1 = 0
        
        http://www.diegm.uniud.it/fusiello/teaching/mvg/elementsCV.pdf
        can determine a homography matrix 
            x2_i is approx H * x1_i for each point i.
        then using cross product:
            x2_i cross H * x1_i = 0
        
        exploit the properties of the Kronecker product and the vec operator to
        transform this into a null-space problem and then derive a linear solution:
        
            x2_i cross H * x1_i = 0
            [x2_i]_x * H * x1_i = 0
            vec( [x2_i]_x * H * x1_i ) = 0
            ( (x1_i)^T kronecker_delta [x2_i]_x) * vec(H) = 0
        
        The rank of ( (x1_i)^T kronecker_delta_product [x2_i]_x) is 2
        The number of equations in ( (x1_i)^T kronecker_delta [x2_i]_x) * vec(H) = 0 
           is 3 
        and it has 9 unknowns.
        
        Because of the rank of the kronecker product, there are only 2 independent equations
        out of the 3.
        
        let A be a factorization matrix for the 2*nPoints of equations
        A is 2*nPoints X 9
        In general A will have rank 8 and the solution is the 1-dimensional 
        right null-space of A.
        So H can be solved for nPoints .geq. 4.
        
        If the data are not exact and more than 4 points are used, the rank of A is 9
        and a least squares solution is sought.
        
        The least-squares solution for vec(H^T) is the singular vector 
        corresponding to the smallest singular value of A.
        
        NOTE: the scalars below are depths. They're the distance from the
            object in world reference to the focal place of the camera.
            Only for a special choice of the world reference frame 
            (the plane at infinity as the refence plane) does this depth
            coincide with the object's third coordinate (Z).
            (when lambda=1 the scalar is the depth of the object;
             where P=lambda*K*[R|t] and K[2][2]=1)
        
        Estimating the parallax:
            parallax * epipole_2 = (scalar_2/scala_1)*x2 - H*x1

            epipole_2, x2 and H*m1 are collinear
        
            (1/parallax) = (epipole_2 cross x2) dot (x2 cross H*x1) / || x2 cross H*x1 ||^2

            because the epipole and homography can only be determined up to a scale,
            the magnitude of the parallax can also only be estimated up to scale.

        Disparity
            Consider two identical cameras separated by a translation along a 
            direction perpendicular to principal axis (w.l.o.g. assume X axis). 
            This is the so called “normal case” for stereo (see also Sec.8.2).
        
            Since the focal planes coincide then scalar_i = scalar_2 and the 
            right epipole is at infinity: epipole_2 = [b*f, 0, 0]^T
                where f is the focal length (in pixels), 
                b is the magnitude of the traslation (in X).
            Moreover, since Ki = K2 then x1 = x1'
                 where x1'is P1[subset 3x3]*(P0[subset 3x3])^-1*x1
                 
            Eq. (29):
                 epipole_2 = scalar_2*x2 − scalar_1*x1'
`           simplifies to:
                [b*f/scala_2, 0, 0}^T = x2 - x1' (72)
        
            The difference of the coordinates of conjugate points have only one 
            non-zero component (horizontal, w.l.o.g.), 
            and this scalar value is called 
            binocular disparity. It is proportional to the reciprocal of the depth.

            if x1 = H*x1' then 1/parallax = scalar_2.
            That occurs when the reference plane is the plane at infinity,
            H_infinity = P1[subset 3x3]*(P0[subset 3x3])^-1.
            And in that case, parallax*[1, 0,0]^T = x2 - x1.
        */

    }
   
    /**
     use the homography from rectify(...) to warp the image img such that
     epipolar lines correspond to scan lines.
     <pre>
     following Chapter 11 of "Invitation to Computer Vision, From Images to Geometric Models"
     by Ma, Soatto, Kosecka,& Sastry 2012 (MASKS).
     the code is adapted from their examples_code/Hwarp.m which is freely available for non-commercial purposes.
     </pre>
     
     * @param img image to be rectified
     * @param h homography transformation matrix
     * @return 
    */
    public static double[][] hBackwardWarp(double[][] img, double[][] h) throws NotConvergedException {

        int ydim = img.length;
        int xdim = img[0].length;

        int i;
        int j;

        Grid grid = meshgridForH(h, xdim, ydim);

        assert(grid.xi.length == ydim);
        assert(grid.xi[0].length == xdim);

        ImageProcessor iP = new ImageProcessor();
        double[][] im1 = MatrixUtil.zeros(ydim, xdim);
        for (j = 0; j < xdim; ++j) {
            for (i = 0; i < ydim; ++i) {
                //im1 = interp2(double(im0),xi,yi,'bilinear');
                if (grid.xi[i][j] > 0 && grid.xi[i][j] < xdim && grid.yi[i][j] > 0 && grid.yi[i][j] < ydim) {
                    im1[i][j] = iP.biLinearInterpolation(img, grid.xi[i][j], grid.yi[i][j]);
                }
            }
        }

        return im1;
    }

    /**
     use the homography from rectify(...) to warp the image img such that
     epipolar lines correspond to scan lines.
     <pre>
     following Chapter 11 of "Invitation to Computer Vision, From Images to Geometric Models"
     by Ma, Soatto, Kosecka,& Sastry 2012 (MASKS).
     the code is adapted from their examples_code/Hwarp.m which is freely available for non-commercial purposes.
     </pre>

     * @param img image to be rectified
     * @param h homography transformation matrix
     * @return
     */
    public static RectifiedImage hBackwardWarp(GreyscaleImage img, double[][] h) throws NotConvergedException {

        int ydim = img.getHeight();
        int xdim = img.getWidth();

        int i;
        int j;

        Grid grid = meshgridForH(h, xdim, ydim);

        assert(grid.xi.length == ydim);
        assert(grid.xi[0].length == xdim);

        ImageProcessor iP = new ImageProcessor();

        RectifiedImage rImg = new RectifiedImage(xdim, ydim);

        double interp;
        int v;
        float x, y;
        for (j = 0; j < xdim; ++j) {
            for (i = 0; i < ydim; ++i) {
                //im1 = interp2(double(im0),xi,yi,'bilinear');
                x = (float)Math.round(grid.xi[i][j]); //[480 X 720]  ydim=480
                y = (float)Math.round(grid.yi[i][j]);
                if (x > 0 && x < xdim && y > 0 && y < ydim) {
                    interp = iP.biLinearInterpolation(img, x, y);
                    v = (int) Math.round(interp);
                    rImg.setRGB(j, i, v, v, v);
                }
            }
        }

        return rImg;
    }

    /**
     use the homography from rectify(...) to create the mesh grid indexes needed
     to warp the image img such that
     epipolar lines correspond to scan lines.
     <pre>
     following Chapter 11 of "Invitation to Computer Vision, From Images to Geometric Models"
     by Ma, Soatto, Kosecka,& Sastry 2012 (MASKS).
     the code is adapted from their examples_code/Hwarp.m which is freely available for non-commercial purposes.
     </pre>

     * @param h homography transformation matrix
     * @param xdim the length of the x-axis of the image to be transformed (i.e. image width)
       @param ydim the length of the y-axis of the image to be transformed (i.e. the image height)
     * @return
     */
    static Grid meshgridForH(final double[][] h, final int xdim, final int ydim) throws NotConvergedException {

        double[] ulc = MatrixUtil.multiplyMatrixByColumnVector(h, new double[]{1,1,1});
        MatrixUtil.multiply(ulc, 1./ulc[2]);

        double[] urc = MatrixUtil.multiplyMatrixByColumnVector(h, new double[]{xdim,1,1});
        MatrixUtil.multiply(urc, 1./urc[2]);

        double[] llc = MatrixUtil.multiplyMatrixByColumnVector(h, new double[]{1,ydim,1});
        MatrixUtil.multiply(llc, 1./llc[2]);

        double[] lrc = MatrixUtil.multiplyMatrixByColumnVector(h, new double[]{xdim,ydim,1});
        MatrixUtil.multiply(lrc, 1./lrc[2]);

        double xmin = MiscMath0.findMin(new double[]{ulc[0],llc[0],urc[0],lrc[0]});
        double xmax = MiscMath0.findMax(new double[]{ulc[0],llc[0],urc[0],lrc[0]});
        double ymin = MiscMath0.findMin(new double[]{ulc[1],urc[1],llc[1],lrc[1]});
        double ymax = MiscMath0.findMax(new double[]{ulc[1],urc[1],llc[1],lrc[1]});

        TDoubleList rangeX = new TDoubleArrayList();
        TDoubleList rangeY = new TDoubleArrayList();
        //range of xmin to xmax inclusive, and of size xdim
        double d = (xmax - xmin)/(xdim - 1);
        int c = 0;
        double ii = xmin;
        while (c < xdim) {
            rangeX.add(ii);
            ii += d;
            ++c;
        }
        //range of ymin to ymax inclusive, and of size ydim
        d = (ymax - ymin)/(ydim - 1);
        c = 0;
        ii = ymin;
        while (c < ydim) {
            rangeY.add(ii);
            ii += d;
            ++c;
        }
        assert(rangeX.size() == xdim);
        assert(rangeY.size() == ydim);

        int j;
        /*
        matlab:
        x = 1, 2, 3;
        y = 1, 2, 3, 4, 5;
        [X2,Y2] = meshgrid(x,y)
        X2 = 5×3
             1     2     3
             1     2     3
             1     2     3
             1     2     3
             1     2     3
        Y2 = 5×3
             1     1     1
             2     2     2
             3     3     3
             4     4     4
             5     5     5
         c=0
         meshgrid then reshape
         for j=0, j<xdim;++j
            for (i = 0; i < ydim; ++i)
                xx[c]=x[j];
                yy[c]=y[i];
                ++c
        but we are using the opposite indexing order, so use fill pattern Y2 for X2 and vice versa
        for (i = 0; i < ydim; ++i)
          for j=0, j<xdim;++j
                xx[c]=x[j];
                yy[c]=y[i];
                ++c
         */
        double[] xx = new double[ydim*xdim];
        double[] yy = new double[ydim*xdim];
        c = 0;
        int i;
        for (i = 0; i < ydim; ++i) {
            for (j = 0; j < xdim; ++j) {
                xx[c] = rangeX.get(j);
                yy[c] = rangeY.get(i);
                ++c;
            }
        }

        double[][] gg = new double[3][];
        gg[0] = Arrays.copyOf(xx, xx.length);
        gg[1] = Arrays.copyOf(yy, yy.length);
        gg[2] = new double[xx.length];
        Arrays.fill(gg[2], 1);

        // gg = [xx; yy; ones(1,ydim*xdim)];// [3 X nr*nc]
        //ww = H \ gg  ==> ww = pInv(H)*gg // [3 X 3] * [3 X nr*nc] = [3 X nr*nc]
        double[][] ww = MatrixUtil.multiply(MatrixUtil.pseudoinverseRankDeficient(h), gg);
        //System.out.printf("ww=pinv(h)*gg=\n%s\n", FormatArray.toString(ww, "%.4e"));
        for (i = 0; i < ww[0].length; ++i) {
            for (j = 0; j < 3; ++j) {
                ww[j][i] /= ww[2][i];
            }
        }
        //double[] wx = ww[0];
        //double[] wy = ww[1];

        //xi = reshape(wx,[ydim,xdim]); // [nr X nc]
        //yi = reshape(wy,[ydim,xdim]); // [nr X nc]
        double[][] xi = MatrixUtil.zeros(ydim, xdim);
        double[][] yi = MatrixUtil.zeros(ydim, xdim);
        // fill by columns
        c = 0;
        for (i = 0; i < ydim; ++i) {
            for (j = 0; j < xdim; ++j) {
                xi[i][j] = ww[0][c];
                yi[i][j] = ww[1][c];
                ++c;
            }
        }

        Grid grid = new Grid();
        grid.xi = xi;
        grid.yi = yi;

        return grid;
    }

    public static class Grid {
        public double[][] xi;
        public double[][] yi;
    }
    
    public static class RectifiedImage extends Image {
        public RectifiedImage(int theWidth, int theHeight) {
            super(theWidth, theHeight);
        }
    }

    public static class RectifiedPoints {
        private double[][] x1;
        private double[][] x2;
        private double[][] h2;
        private double[][] h1;

        /**
         * @return the x1
         */
        public double[][] getX1() {
            return x1;
        }

        /**
         * @param x1 the x1 to set
         */
        public void setX1(double[][] x1) {
            this.x1 = x1;
        }

        /**
         * @return the x2
         */
        public double[][] getX2() {
            return x2;
        }

        /**
         * @param x2 the x2 to set
         */
        public void setX2(double[][] x2) {
            this.x2 = x2;
        }

        /**
         * @return the h2
         */
        public double[][] getH2() {
            return h2;
        }

        /**
         * @param h2 the h2 to set
         */
        public void setH2(double[][] h2) {
            this.h2 = h2;
        }

        /**
         * @return the h1
         */
        public double[][] getH1() {
            return h1;
        }

        /**
         * @param h1 the h1 to set
         */
        public void setH1(double[][] h1) {
            this.h1 = h1;
        }
    }

}
