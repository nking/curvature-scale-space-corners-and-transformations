package algorithms.imageProcessing.transform;

import algorithms.imageProcessing.features.RANSACSolver;
import algorithms.imageProcessing.matching.ErrorType;
import algorithms.matrix.MatrixUtil;
import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;
import java.util.Arrays;
import java.util.List;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.NotConvergedException;

/**
 *
 * @author nichole
 */
public class Rectification {

    // notes from Serge Belongie lectures from Computer Vision II, CSE 252B, USSD
    // and Ma, Soatto, Kosecka,& Sastry "Invitation to Computer Vision, From Images to Geometric Models",
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
     * implementation of epipolar rectification following algorithm 11.9 of
     * Ma, Soatto, Kosecka,& Sastry "Invitation to Computer Vision, From Images to Geometric Models".
     * This is for un-calibrated cameras.  If the images have a large range of
     * depth in them or if the epipoles are inside the images, 
     * this algorithm can result in distortions.
     * If one has camera intrinsic and extrinsic parameters, Ma et al. suggest
     * use method of Fusiello et al. 1997 for Euclidean projection. 
     * 
     * *@param x1 the image 1 set of correspondence points. format is 3 x N
     * where N is the number of points. NOTE: since intrinsic parameters are not
     * known, users of this method should presumably center the coordinates in
     * some manner (e.g. subtract the image center or centroid of points) since
     * internally an identity matrix is used for K.
     * @param x1
     * @param x2 the image 2 set of correspondence points. format is 3 x N where
     * N is the number of points. NOTE: since intrinsic parameters are not
     * known, users of this method should presumably center the coordinates in
     * some manner (e.g. subtract the image center or centroid of points).
     * @param oX camera optical center along x-axis
     * @param oY camera optical center along y-axis
     * @return
     */
    public static RectifiedImage rectify(double[][] x1, double[][] x2, double oX, double oY) throws NoSuchAlgorithmException, NotConvergedException {

        if (x1.length != 3 || x2.length != 3) {
            throw new IllegalArgumentException("x1.length must be 3 and so must x2.length");
        }
        int n0 = x1[0].length;
        if (x2[0].length != n0) {
            throw new IllegalArgumentException("x1 and x2 must be same dimensions");
        }
        if (n0 < 7) {
            throw new IllegalArgumentException("need at least 7 points for the uncalirated camera (s) 2 view solution");
        }

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

        RANSACSolver solver = new RANSACSolver();
        fitR = solver.calculateEpipolarProjection(
                leftM, rightM, errorType, useToleranceAsStatFactor, tolerance,
                reCalcIterations, false);

        DenseMatrix fm = EpipolarTransformer.denormalizeTheFundamentalMatrix(
                fitR.getFundamentalMatrix(),
                normXY1.getNormalizationMatrices(),
                normXY2.getNormalizationMatrices());

        double[][] _fm = MatrixUtil.convertToRowMajor(fm);

        //x1M = extractIndices(x1M, fitR.inlierIndexes);
        //x2M = extractIndices(x2M, fitR.inlierIndexes);
        //x1 = MatrixUtil.convertToRowMajor(x1M);
        //x2 = MatrixUtil.convertToRowMajor(x2M);

        int n = x1[0].length;

        System.out.println("RANSAC fit=" + fitR.toString());

        EpipolarTransformer tr = new EpipolarTransformer();
        double[][] e1e2 = tr.calculateEpipoles(fm);
        double[] e2 = e1e2[1];
        
        /*
        
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
        
        If the data are not exact and more than 4 points are used, therank of A is 9
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
        
        // a translation:
        double[][] gT = MatrixUtil.createIdentityMatrix(3);
        gT[0][2] = -oX;
        gT[1][2] = -oY;
        
        // see projRectify.m in examples-code from Ma et al. supplementary book material
        // https://cs.gmu.edu/~kosecka/MASKS_book.html
        // which states:
        // "THE CODE ON THIS PAGE IS DISTRIBUTED FREE FOR NON-COMMERCIAL USE.
        // Copyright (c) MASKS, 2003."
        
        // gR is a rotation:
        // gR * gT * e2 = [e2_x, 0, 1]^T
        double[] p2T = MatrixUtil.multiplyMatrixByColumnVector(gT, e2);
        // rotate the epipole to lie on the x-axis
        double theta = Math.atan(-p2T[1]/p2T[0]);
        // rotation about z-axis:
        double[][] gR = Rotation.createEulerYawRotationMatrix(theta);
        double[] p2R = MatrixUtil.multiplyMatrixByColumnVector(gR, p2T);
       
        double[][] g = MatrixUtil.createIdentityMatrix(3);
        g[2][0] = -1./p2R[0];
        
        double[] pim2R = MatrixUtil.multiplyMatrixByColumnVector(g, p2R);
        
        double[][] h2 = MatrixUtil.multiply(g, gR);
        h2 = MatrixUtil.multiply(h2, gT);
        
        double[][] M = MatrixUtil.multiply(
            MatrixUtil.transpose(MatrixUtil.skewSymmetric(e2)), _fm);
            //  + ep2*rand(1,3);  3X1 * 1X3 = 3X3 
        double[] randV = new double[3];
        SecureRandom rand = SecureRandom.getInstanceStrong();
        randV[0] = rand.nextDouble();
        randV[1] = rand.nextDouble();
        randV[2] = rand.nextDouble();
        double[][] randEp2 = MatrixUtil.outerProduct(e2, randV);
        M = MatrixUtil.elementwiseAdd(M, randEp2);
        
        // determine H then H1 by solving for unknown plane v to minimize the 
        //   disparity in matching homography
        
        // from Chap 11 near eqn (11.30) of Ma, Soatto, et al.)
        //  algebraic error assoc w/ homography transfer:
        //     [x2]_x * H * x1 = [x2]_x * ( ([t]_x)^T * F + t*v^T ) * x1 is approx 0
        //     where [*]_x is the skew-symmetric matrix used in cross product operations
        //     and t is the 2nd epipole of the fundamental matrix.
        
        double t1 = e2[0]; 
        double t2 = e2[1]; 
        double t3 = e2[3];
        double[][] a = new double[2*n][3];
        double[] b = new double[2*n];
        for (int i = 0; i < n;++i) {
            a[2*i] = new double[]{
                -t2*x1[0][i] + t3*x1[0][i]*x2[1][i],
                -t2*x1[1][i] + t3*x1[1][i]*x2[1][i], -t2+t3*x2[1][i]};
            a[2*i + 1] = new double[]{
                t1*x1[0][i] - t3*x1[0][i]*x2[0][i],
                t1*x1[1][i] - t3*x1[1][i]*x2[0][i], t1-t3*x2[0][i]};
            b[2*i] = M[1][0]*x1[0][i] + M[1][1]*x1[1][i] + M[1][2] - M[2][0]*x1[0][i]*x2[1][i]
                  - M[2][1]*x1[1][i]*x2[1][i] - M[2][2]*x2[1][i];
            b[2*i + 1] = -M[0][0]*x1[0][i] - M[0][1]*x1[1][i] - M[0][2] + M[2][0]*x1[0][i]*x2[0][i]
                  + M[2][1]*x1[1][i]*x2[0][i] + M[2][2]*x2[0][i];
        }
        
        // 3 X n
        double[][] aInv = MatrixUtil.pseudoinverseFullRank(a);
        
        // 3 x 1
        //aa = A\b;
        double[] aa = MatrixUtil.multiplyMatrixByColumnVector(aInv, b);
        
        // 3 X 3
        //H = M + ep2*aa';
        double[][] h = MatrixUtil.elementwiseAdd(M, MatrixUtil.outerProduct(e2, aa));
        
        // 3 X 3
        //H1 = H2*H;
        double[][] h1 = MatrixUtil.multiply(h2, h);
        
        double[][] x1R = MatrixUtil.multiply(h1, x1);
        double[][] x2R = MatrixUtil.multiply(h2, x2);
        
        // normalize z-coords to be 1
        for (int i = 0; i < n; ++i) {
            x1R[0][i] /= x1R[2][i];
            x2R[0][i] /= x2R[2][i];
        }
        
        RectifiedImage out = new RectifiedImage();
        out.x1 = x1R;
        out.x2 = x2R;
        out.h1 = h1;
        out.h2 = h2;

        return out;        
    }
    
    /**
     use the homography from rectify(...) to warp the image img such that
     epipolar lines correspond to scan lines.
     following Hwarp.m in examples-code from Ma et al. supplementary book material
     https://cs.gmu.edu/~kosecka/MASKS_book.html
     which states:
     "THE CODE ON THIS PAGE IS DISTRIBUTED FREE FOR NON-COMMERCIAL USE.
     Copyright (c) MASKS, 2003."
        
     * @param img two dimension array holding pixel intensities in format
     * where x axis is along columns and y axis is along rows.
     * @param h
     * @return 
     */
    public static RectifiedImage hWarp(double[][] img, double[][] h) throws NotConvergedException {
        
        int ydim = img.length/2;
        int xdim = img[0].length/2;
        
        //NOTE: handling the zero-base coordinates offset at end
        
        // upper, lower, left and right corners
        double[] ulc = MatrixUtil.multiplyMatrixByColumnVector(h, new double[]{1, 1, 1});
        MatrixUtil.multiply(ulc, 1./ulc[2]);
        double[] urc = MatrixUtil.multiplyMatrixByColumnVector(h, new double[]{xdim, 1, 1});
        MatrixUtil.multiply(urc, 1./urc[2]);
        double[] llc = MatrixUtil.multiplyMatrixByColumnVector(h, new double[]{1, ydim, 1});
        MatrixUtil.multiply(llc, 1./llc[2]);
        double[] lrc = MatrixUtil.multiplyMatrixByColumnVector(h, new double[]{xdim, ydim, 1});
        MatrixUtil.multiply(lrc, 1./lrc[2]);
  
        // compute the new meshgrid 
        int xmin = (int)Math.min(Math.min(Math.min(ulc[0], llc[0]), urc[0]), lrc[0]);
        int xmax = (int)Math.max(Math.max(Math.max(ulc[0], llc[0]), urc[0]), lrc[0]);
        int ymin = (int)Math.min(Math.min(Math.min(ulc[1], llc[1]), urc[1]), lrc[1]);
        int ymax = (int)Math.max(Math.max(Math.max(ulc[1], llc[1]), urc[1]), lrc[1]);
                
        // generate coordinates in the new image  
        //range_x = xmin:xmax;
        //range_y = ymin:ymax;
        //[x2, y2] = meshgrid(range_x,range_y); % original image 
        // ==> x2 is a double array of size (ymax-ymin+1) X (xmax-xmin+1)
        //        where each row is [xmin, xmin+1, ... xmax]
        // ==> y2 is a double array of size ydim X xdim
        //        where each column is [ymin, ymin+1, ... ymax]
        
        //[ydim, xdim] = size(x2);
        ydim = ymax-ymin+1;
        xdim = xmax-xmin+1;
        
        int len = ydim*xdim;
        
        //xx = reshape(x2,[1,ydim*xdim]);
        //yy = reshape(y2,[1,ydim*xdim]);
        // matlab's reshape fills column by column
        double[] xx = new double[len];
        double[] yy = new double[len];
        int row, col, c;
        c = 0;
        // x2 is xmin, xmin+1, ... xmax
        //       xmin, xmin+1, ... xmax
        //       ...
        for (col = xmin; col <= xmax; ++col) {
            for (row = 0; row < ydim; ++row) {
                xx[c] = col;
                c++;
            }
        }
        // y2 is ymin,   ymin,   ...
        //       ymin+1, ymin+1, ... 
        //       ...
        //       ymax    ymax
        c = 0;
        for (row = 0; row < xdim; ++row) {
            for (col = ymin; col <= ymax; ++col) {
                yy[c] = col;
                c++;
            }
        }
        
        //gg = [xx; yy; ones(1,ydim*xdim)];
        //ww = H \ gg;
        double[][] gg = new double[3][len];
        gg[0] = xx;
        gg[1] = yy;
        gg[2] = new double[len];
        Arrays.fill(gg[2], 1.);
        
        // 3x3
        double[][] hInv = MatrixUtil.pseudoinverseRankDeficient(h);
        // 3Xlen
        double[][] ww = MatrixUtil.multiply(hInv, gg);
        
        double[] wx = ww[0];//Arrays.copyOf(ww[0], len);
        double[] wy = ww[1];//Arrays.copyOf(ww[1], len);
        for (col = 0; col < ww[0].length; ++col) {
            wx[col] /= ww[2][col];
            wy[col] /= ww[2][col];
        }
        
        //xi = reshape(wx,[ydim,xdim]);
        //yi = reshape(wy,[ydim,xdim]);        
        double[][] xi = new double[ydim][xdim];
        double[][] yi = new double[ydim][xdim];
        //  00  01   =>  00  20  01
        //  10  11       10  01  11
        //  20  21
        c = 0;
        for (col = 0; col < xdim; ++col) {
            for (row = 0; row < ydim; ++row) {
                xi[row][col] = wx[c];
                yi[row][col] = wy[c];
                c++;
            }
        }
        
        //im1 = interp2(double(im0),xi,yi,'bilinear');
        //
        // NOTE: need to handle the zero-based coordinate system
        // im0[yi[0][0][xi[0][0]]
        
        throw new UnsupportedOperationException("unfinished");
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

    public static class RectifiedImage {
        public double[][] x1;
        public double[][] x2;
        public double[][] h2;
        public double[][] h1;
    }

}
