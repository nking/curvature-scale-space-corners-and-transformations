package algorithms.imageProcessing.transform;

import algorithms.matrix.BlockMatrixIsometric;
import algorithms.matrix.MatrixUtil;
import algorithms.util.FormatArray;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.set.TIntSet;
import java.io.IOException;
import java.util.Arrays;
import java.util.logging.Logger;

import no.uib.cipr.matrix.*;

/**
 given intrinsic and extrinsic camera parameters, coordinates for points
 in a world reference frame, and observations of those points in one or more
 camera, return data structures needed by Levenberg-Marquardt algorithm
 in refinement of the intrinsic and extrinsic camera parameters.
 BundleAdjustment calculates partial derivatives of the parameters
 and calculates the re-projection error to form the parameter update steps,
 the gradient covector, and the evaluation of the objective (minimize sum of squares of
 the re-projection error).

 The main method is solveSparsely(...).
 
 * From Triggs et al. 2000, "Bundle Adjustment - A Modern Synthesis"
 * "Bundle adjustment is the problem of refining a visual reconstruction to 
 * produce jointly optimal 3D structure and viewing parameter (camera pose and/or 
 * calibration) estimates. Optimal means that the parameter estimates are found 
 * by minimizing some cost function that quantifies the model fitting error, and 
 * jointly that the solution is simultaneously optimal with respect to both 
 * structure and camera variations. The name refers to the ‘bundles’ of light 
 * rays leaving each 3D feature and converging on each camera centre, which are 
 * ‘adjusted’ optimally with respect to both feature and camera positions. 
 * Equivalently — unlike independent model methods, which merge partial 
 * reconstructions without up-dating their internal structure — all of the
 * structure and camera parameters are adjusted together ‘in one bundle’."
 * 
 * This version of Bundle-Adjustment uses Levenberg-Marquardt Algorithm (LMA) 
 * step in each iteration solved by dense Cholesky decomposition (LMA-CD).
 * 
 * TODO: consider implementing or finding an implementation of 
 * Jakob Engel, Vladlen Koltun, and Daniel Cremers. 
 * "Direct sparse odometry". 
 * IEEE transactions on pattern analysis and machine intelligence, 
 * 40(3):611–625, 2018.
 *
 * TODO:  add gauge fix.  And related to that, consider adding constraints 
 * suggested in Seliski 2010: u_0 and v_0 are close to half the image lengths 
 * and widths, respectively.  the angle between 2 image axes is close to 90.
  the focal lengths along both axes are greater than 0.     
 * 
 * TODO: consider implementing the "reduced structure system" for the cases
 * where (9^3)*mImages > (3^3)*nFeatures,  The "reduced camera system" is
 * currently implemented only, but "reduced structure system"
 * for some cases would be more efficient.
 * 
 * <pre>
 * for a bigger picture summary, see Section 
 * "Objective functions for estimating epipolar geometry" on page 169 of
 *  MASKS (Ma, Soatto, Kosecká, and Sastry 2012, 
 * "An Invitation to 3-D Vision").
 * Also see Chap 5.2, pp 127-128 of MASKS for constrained optimization.
 * </pre>
 <pre>
 Note, if need to estimate the intrinsic camera for initial conditions, one can generally
 start with:
 MASKS Algorithm 11.6, step 1:
 Guess a calibration matrix K by choosing the optical center at the center of the image,
 assuming the pixels to be square, and guessing the focal length f. For example, for
 an image plane of size (Dx X Dy) pixels, a typical guess is
     | f O Dx/2
 K = | 0 f Dy/2
     | 0 0 1
 with f = k X Dx, where k is typically chosen in the interval [0.5, 2].
 </pre>
 * @author nichole
 */
public class BundleAdjustment {

    private int useHomography = 0;
        
    //private static final Level LEVEL = Level.INFO;
    private static final Logger log = Logger.getLogger(CameraCalibration.class.getSimpleName());
    
    public BundleAdjustment() {
        useHomography = 0;
    }

    // use negative because the gradients use convention -J^T*f where f.
    private static final double updateSign = -1;

    private static final boolean useBouguetDerivs = true;
    
    /**
     * setting to use the homography matrix in the WCS projection to camera
     * 2-D reference frame.   The homography matrix is the first 2 columns of
     * the rotation matrix and the translation vector as the 3rd column
     * in the homography.  
     * The default uses, instead, the full projection
     * of Rotation * x_WCS + translation.
     * Note that the partial derivatives calculated in making
     * the Jacobian are derived from the full projection, not the homography,
     * so there will be inconsistencies until that is fixed for the homography
     * derivatives.
     */
    public void setUseHomography() {
        useHomography =1;
    }
    
    /**
     * NOTE: this method could use improvements and more testing.
     *
     * given world scene features, the features observed in images,
     * initial camera calibration and extrinsic parameters, use the 
     * iterative non-linear Levenberg-Marquardt (L-M)
     * algorithm to minimize the re-projection error by refining the values of
     * coordsW, intrinsic, and extrinsic camera parameters.
     * 
     * The L-M is an iterative non-linear optimization to minimize the 
     * objective.  For bundle adjustment, the objective is the 
     * re-projection error.
     * L-M is guaranteed to converge to eventually find an improvement, 
     * because an update with a sufficiently small magnitude and a negative scalar 
     * product with the gradient is guaranteed to do so.
     * 
     * This method exploits some of the properties of sparse matrices in the 
     * block structure of the Jacobian by using the Schur complement to form
     * reduced camera matrices which can be solved by Cholesky factoring and
     * forward and backward substitution (halving the computation time for the
     * parameter vector).
     * Regarding the bundle adjustment as refinement for extrinsic parameters:
        <pre>
            However, several researchers have noted (Fitzgibbon and 
            Zisserman, 1998, Nistér, 2001, Pollefeys, 1999) that in the 
            application of camera tracking, performing bundle adjustment each 
            time a new frame has been added to the estimation can prevent the 
            tracking process from failing over time. 
            Thus, bundle adjustment can over time in a sequential estimation 
            process have a much more dramatic impact than mere accuracy 
            improvement, since it improves the initialization for future 
            estimates, which can ultimately enable success in cases that would 
            otherwise miserably fail.
        </pre>
     <pre>
     References:
     
     http://users.ics.forth.gr/~lourakis/sba/PRCV_colloq.pdf
     lecture by Lourakis  “Bundle adjustment gone public”
     
     Zhongnan Qu's master thesis, "Efficient Optimization for Robust Bundle 
     Adjustment", 2018 Technical University of Munich
     
     Chen, Chen, & Wang 2019, "Bundle Adjustment Revisited"
     
     Engels, Stewenius, Nister 2006, “Bundle Adjustment Rules”
      
     Bill Triggs, Philip Mclauchlan, Richard Hartley, Andrew Fitzgibbon. 
     Bundle Adjustment – A Modern Synthesis. 
     International Workshop on Vision Algorithms, 
     Sep 2000, Corfu, Greece. pp.298–372, 
     10.1007/3-540-44480-7_21 . inria-00548290 

     T. Barfoot, et al., "Pose estimation using linearized rotations and 
     quaternion algebra", Acta Astronautica (2010), doi:10.1016/j.actaastro.2010.06.049
         -- using the rotation and translation update details.
         -- one of the 2 examples is interesting for the problem of pose for
            a pair of stereo-images.  it also uses cholesky factoring of block
            sparse matrix structure.
           
     Tomasi 2007,CPS 296.1 Supplementary Lecture Notes, Duke University
     
     graph partitioning:
        https://cseweb.ucsd.edu/classes/fa04/cse252c/manmohan1.pdf
        recursive partitioning w/ elimination graph and vertex cut.
        
        Triggs et al. 2000, "Bundle Adjustment – A Modern Synthesis", Section 6
        
        Skeletal graphs for efficient structure from motion
        N Snavely, SM Seitz, R Szeliski - 2008
 
     Graph partitioning in this project:
        NormalizedCuts.java which uses the Fiedler vector of the Laplacian.
        UnweightedGraphCommunityFinder.java
        
     </pre>
     TODO: review nad improve the derivatives.  need good test data.
     * @param coordsI the features observed in different images (in coordinates
     * of the image reference frame).
     * The format of coordsI is 3 X (nFeatures*nImages). Each row should
     * have nFeatures of one image, followed by nFeatures of the next image,
       etc.  The first dimension is for the x,y, and z axes.
       Note that if a feature is not present in the image, that should be
       an entry in imageMissingFeatureMap.
     * @param coordsW the features in a world coordinate system.  The format is
     * 3 X nFeatures.  The first dimension is for the x,y, and z axes.
     * @param imageFeaturesMap an associative array holding the features
     * in each image.  They key is the image number in coordsI
     * which is j/nFeatures where j is the index of the 2nd dimension,
     * that is coordsI[j].  The value is a set of feature numbers which are
     * missing from the image.  The feature numbers correspond to the
     * 2nd dimension indexes in coordsW.
     * @param intr the intrinsic camera parameter matrices stacked along rows
     * to make a tall double array of size [(mImages*3) X 3] where each image block is
     * size 3X3.  Note that only the focus parameter is refined in this method.
     * @param extrRVecs the extrinsic camera parameter (Rodrigues) rotation vectors
     * stacked along the 3 columns, that is the size is [mImages X 3].
     * each image block is size 1X3.
     * @param extrTrans the extrinsic camera parameter translation vectors
     * stacked along the 3 columns, that is the size is [nImages X 3] where
     * nImages is coordsI[0].length/coordsW[0].length.  each array is size
     * 1X3.
     * @param kRadials a double array wherein each row holds the
     * radial distortion coefficients k1 and k2 for an image.
     * NOTE: current implementation accepts values of 0 for k1 and k2.
     * @param nMaxIter
     * @param useBouguetForRodrigues if true, uses only the Bouguet algoirthms for Rodrigues rotation matrices and vectors
     * @param useR2R4 useR2R4 use radial distortion function from Ma et al. 2004 for model #4 in Table 2,
        f(r) = 1 +k1*r^2 + k2*r^4 if true,
        else use model #3 f(r) = 1 +k1*r + k2*r^2.
     * @throws no.uib.cipr.matrix.NotConvergedException
     * @throws java.io.IOException 
     */
    public void solveSparsely(
            double[][] coordsI, double[][] coordsW, TIntObjectMap<TIntSet> imageFeaturesMap,
            BlockMatrixIsometric intr, double[][] extrRVecs, double[][] extrTrans,
            double[][] kRadials, final int nMaxIter, boolean useR2R4, boolean useBouguetForRodrigues)
        throws NotConvergedException, IOException {

        int nFeatures = coordsW[0].length;
        int mImages = coordsI[0].length / nFeatures;

        if (nFeatures < 6) {
            throw new IllegalArgumentException("coordsW[0].length must be at least 6");
        }
        if ((coordsI[0].length % nFeatures) > 0) {
            throw new IllegalArgumentException("the number of images present in coordsI should"
                    + " be evenly divided by the number of features in coordsW (that is, coordsI should"
                    + " have that same number of features for each image)");
        }
        if (coordsI.length != 3) {
            throw new IllegalArgumentException("coordsI must have 3 rows.");
        }
        if (coordsW.length != 3) {
            throw new IllegalArgumentException("coordsW must have 3 rows.");
        }
        if (coordsI[0].length != nFeatures * mImages) {
            throw new IllegalArgumentException("coordsI[0].length must be evenly "
                    + "divisible by nFeatures which is coordsW[0].length");
        }
        if (intr.getA().length != 3 * mImages) {
            throw new IllegalArgumentException("intr.length must be 3*mImages");
        }
        if (intr.getA()[0].length != 3) {
            throw new IllegalArgumentException("intr[0].length must be 3");
        }
        if (kRadials.length != mImages) {
            throw new IllegalArgumentException("kRadials.length must be equal "
                    + "to the number of cameras.");
        }
        if (kRadials[0].length != 2) {
            throw new IllegalArgumentException("kRadials[0].length must be 2.");
        }
        if (extrRVecs[0].length != 3) {
            throw new IllegalArgumentException("extrRVecs[0].length must be 3");
        }
        if (extrRVecs.length != mImages) {
            throw new IllegalArgumentException("extrRVecs.length must be mImages "
                    + "where mImages = coordsI[0].length/coordsW[0].length");
        }
        if (extrTrans[0].length != 3) {
            throw new IllegalArgumentException("extrTrans[0].length must be 3");
        }
        if (extrTrans.length != mImages) {
            throw new IllegalArgumentException("extrTrans.length must be mImages "
                    + "where mImages = coordsI[0].length/coordsW[0].length");
        }
        if (imageFeaturesMap == null) {
            throw new IllegalArgumentException("imageFeaturesMap cannot be null");
        }
        if (imageFeaturesMap.size() != mImages) {
            throw new IllegalArgumentException("imageFeaturesMap size must equal "
                    + "the number of images which = coordsI[0].length/coordsW[0].length");
        }

        //TODO: as part of "fix the gauge", need to consider the first
        //    camera to have axes aligned with world axes, so that the
        //    origin of the world scene is [0,0,0] or [0,0,0,1]
        //    NLK: see MASKS Table 6.5.
                
        /*
        RCM is the reduced camera matrix in the augmented normal equation.
        
        The RCM can be solved without inverting A since it is a sparse matrix:
              https://www.cvg.ethz.ch/teaching/3dvision/2014/slides/class06eth14.pdf
              Sparse Matrix factorization:
                  (1) LU factorization: A = L*U where original equation is A*x = b
                      (a) use LU decomposition to solve L,U = luDecomp(A)
                      (b) let c = U*x and rewrite A*x=b as L*c = b.
                            then solve c from forward elimination
                      (c) solve x from back substitution in c = U*x
                  (2) QR factorization: A = Q*R
                  (3) Cholesky factorization:A = L*L^T
                  Some details for the above 3 are in Triggs et al. 2000 Appendix B and near page 23
                  ("Bundle Adjustment – A Modern Synthesis",
                  Springer-Verlag, pp.298–372, 2000, 
                  Lecture Notes in Computer Science (LNCS), 
                  10.1007/3-540-44480-7_21. inria-00590128)
                  https://hal.inria.fr/file/index/docid/590128/filename/Triggs-va99.pdf
              Or Iterative methods:
                  (1) Conjugate Gradient
                  (2) Gauss-Seidel
        
              consider Google Ceres  https://code.google.com/p/ceres-solver/
              or Lourakis SBA
        
             more suggestions in http://users.ics.forth.gr/~lourakis/sba/PRCV_colloq.pdf
               (1) Store as dense, decompose with ordinary linear algebra ◦
                    M. Lourakis, A. Argyros: SBA: A Software Package For Generic 
                       Sparse Bundle Adjustment. ACM Trans. Math. Softw. 36(1): (2009) ◦ 
                    C. Engels, H. Stewenius, D. Nister: Bundle Adjustment Rules. 
                       Photogrammetric Computer Vision (PCV), 2006. 
               (2) Store as sparse, factorize with sparse direct solvers ◦ 
                    K. Konolige: Sparse Sparse Bundle Adjustment. BMVC 2010: 1-11
               (3) Store as sparse, use conjugate gradient methods memory efficient, 
                    iterative, precoditioners necessary! ◦ 
                    S. Agarwal, N. Snavely, S.M. Seitz, R. Szeliski: 
                       Bundle Adjustment in the Large. ECCV (2) 2010: 29-42 ◦ 
                    M. Byrod, K. Astrom: Conjugate Gradient Bundle Adjustment. 
                       ECCV (2) 2010: 114-127 
               (4) Avoid storing altogether ◦ 
                    C. Wu, S. Agarwal, B. Curless, S.M. Seitz: Multicore Bundle 
                       Adjustment. CVPR 2011: 30 57-3064 ◦ 
                    M. Lourakis: Sparse Non-linear Least Squares Optimization 
                       for Geometric Vision. ECCV (2) 2010: 43-56
        
              and computer vision library in java http://boofcv.org/index.php?title=Example_Sparse_Bundle_Adjustment
              there is also a java binding for SBA http://seinturier.fr/jorigin/jsba.html
        
              http://users.ics.forth.gr/~lourakis/sparseLM/
        
            And more suggestsions in Qu:
               (1) LMA-Cholesky-Decomposition
               (2) Inexact increment step solved by Preconditioned Conjugate Gradient 
               (3) Local parameterization
               (4) IRLS
               (5) (Adaptive) DINM algorithm and DINM with local parameterization
               (6)
        */
            
        /* Triggs et al. 2000: "state updates should be evaluated using a stable 
        local parametrization based on increments from the current estimate".
        and this on 3F points:
        "3D points: Even for calibrated cameras, 
        ** ==> vision geometry and visual reconstructions are intrinsically projective. 
        If a 3D (X Y Z)^⊤ parametrization (or equivalently a homogeneous affine (X Y Z 1)^⊤ one)
        is used for very distant 3D points, large X, Y, Z displacements are 
        needed to change the image significantly. I.e., in (X Y Z) space the 
        cost function becomes very flat and steps needed for cost adjustment 
        become very large for distant points. 
            In comparison, with a homogeneous projective parametrization (X Y Z W)^⊤,
        the behaviour near infinity is natural, finite and well-conditioned so 
        long as the normalization keeps the homogeneous 4-vector finite at 
        infinity (by sending W → 0 there). In fact, there is no immediate 
        visual distinction between the images of real points near infinity and 
        virtual ones ‘beyond’ it (all camera geometries admit such virtual 
        points as bona fide projective constructs).  The optimal reconstruction 
        of a real 3D point may even be virtual in this sense, if image noise 
        happens to push it ‘across infinity’.  Also, there is nothing to stop a 
        reconstructed point wandering beyond infinity and back during the 
        optimization. This sounds bizarre at first, but it is an inescapable 
        consequence of the fact that the natural geometry and error model for 
        visual reconstruction is projective rather than affine. 
        Projectively, infinity is just like any other place. 
        ** ==> Affine parametrization (X Y Z 1)^⊤ is acceptable for points near the
        origin with close-range convergent camera geometries, but it is 
        disastrous for distant ones because it artificially cuts away half of 
        the natural parameter space, and hides the fact by sending the resulting 
        edge to infinite parameter values.
        **==> Instead, you should use a homogeneous parametrization (X Y Z W )^⊤
        for distant points, e.g. with spherical normalization summation X_i^2 = 1."
        
        on Rotations:
        "Rotations should be parametrized using either quaternions subject to 
        ∥q∥2 = 1, or local perturbations R*δR or δR*R of an existing rotation R, 
        where δR can be any well- behaved 3 parameter small rotation 
        approximation, e.g. δR = (I + [ δr ]_×), the Rodriguez formula, 
        local Euler angles, etc."

        on State Updates:
        "State updates: Just as state vectors x represent points in some 
        nonlinear space, state updates x → x+ δx represent displacements in 
        this nonlinear space that often can not be represented exactly by vector
        addition. Nevertheless, we assume that we can locally linearize the 
        state manifold, locally resolving any internal constraints and freedoms 
        that it may be subject to, to produce an unconstrained vector δx 
        parametrizing the possible local state displacements. 
        We can then, e.g., use Taylor expansion in δx to form a local cost 
        model f(x + δx)."
        
        on Error Modeling:
        "A typical ML cost function would be the summed negative log likelihoods 
        of the prediction errors of all the observed image features. For 
        Gaussian error distributions, this reduces to the sum of squared 
        covariance-weighted prediction errors (§3.2). A MAP estimator would 
        typically add cost terms giving certain structure or camera calibration 
        parameters a bias towards their expected values."
        ...
        "One of the great strengths of adjustment computations is their ability 
        to combine information from disparate sources. Assuming that the sources 
        are statistically independent of one another given the model, the total 
        probability for the model given the combined data is the product of the 
        probabilities from the individual sources. To get an additive cost 
        function we take logs, so the total log likelihood for the model given 
        the combined data is the sum of the individual source log likelihoods."
        ...
        "Information usually comes from many independent sources. In bundle 
        adjustment these include: covariance-weighted reprojection errors of 
        individual image features; other measurements such as 3D positions of 
        control points, GPS or inertial sensor readings; predictions from 
        uncertain dynamical models (for ‘Kalman filtering’ of dynamic cameras 
        or scenes); prior knowledge expressed as soft constraints (e.g. on 
        camera calibration or pose values); and supplementary sources such as 
        overfitting, regularization or description length penalties."
        
        see section 3.1 footnote 2.
        
        Regarding step control:
            recommends professional software or
            see 
               R. Fletcher. Practical Methods of Optimization. John Wiley, 1987.
               J. Nocedal and S. J. Wright. Numerical Optimization. Springer-Verlag, 1999.
               P. Gill, W. Murray, and M. Wright. Practical Optimization. Academic Press, 1981
        
        In bundle adjustment, certain well-known ambiguities 
        (poorly-controlled parameter combinations) often dominate the uncertainty. 
        Camera distance and focal length estimates, 
        and structure depth and camera baseline ones (bas-relief), 
        are both strongly correlated whenever the perspective is weak 
        (note: from wikipedia: weak perspective is used when when the depth of
        the object along the line of sight is small compared to the distance 
        from the camera, and the field of view is small.  hence, all points on 
        a 3D object are at the same distance Z_avg from the camera without 
        significant errors in the projection )
        and become strict ambiguities in the affine limit. The well-conditioned 
        diagonal blocks of the Hessian give no hint of these ambiguities: when 
        both features and cameras are free, the overall network is much less 
        rigid than it appears to be when each treats the other as fixed.
        
        ** ==> For updates involving a previously unseen 3D feature or image, 
        new variables must also be added to the system.
        (see page 33 in Section 8.1 and Section 8.2 on page 35)
        ...If these parameters are eliminated using reduction (19), the 
        observation update can be applied directly to the reduced Hessian and 
        gradient. The eliminated parameters can then be updated by simple 
        back-substitution (19) and their covariances by (17). In particular, 
        if we cease to receive new information relating to a block of parameters 
        (an image that has been fully treated, a 3D feature that has become 
        invisible), they and all the observations relating to them can be 
        subsumed once-and-for-all in a reduced Hessian and gradient on the 
        remaining parameters. If required, we can later re-estimate the 
        eliminated parameters by back-substitution. Otherwise, we do not 
        need to consider them further.
        */

        //TODO: consider adding constraints suggested in Szeliski 2010:
        // u_0 and v_0 are close to half the image lengths and widths, respectively.
        // the angle between 2 image axes is close to 90.
        // the focal lengths along both axes are greater than 0.

        //factor to raise or lower lambda.  
        //   consider using the eigenvalue spacing of J^T*J (Transtrum & Sethna, "Improvements to the Levenberg-Marquardt algorithm for nonlinear least-squares minimization")
        double lambdaF = 2;
        double eps = 1E-12;

        // (J^T*J + lambda*I) * deltas = -J^T*f where f is the reprojection residuals

        //In a single reprojection error formula, there are altogether 12 arguments 
        //   (9 camera parameters and 3 feature point positions).

        // deltas the point parameters (== world coordinate features)
        double[] outDP = new double[3 * nFeatures];
        // deltas for the camera parameters
        double[] outDC = new double[9 * mImages];

        // Qu array u for parameters is ordered: rot_0, trans_0, intr_0, ...rot_m-1, trans_m-1, intr_m-1, then x_0, ... x_n
        // but the delta parameter array for all params is ordered:
        //     dRot_0, ... dRot_m-1,  dTrans_0, ...dTrans_m-1, dIntr_0,...dIntr_m-1, dX_0, ... dX_n
        // gradient g is same length   

        //the gradient covector for point parameters.  used in calc gain ration and stopping
        double[] outGradP = new double[3 * nFeatures];
        // the gradient covector for camera parameters.  used in calc gain ration and stopping
        double[] outGradC = new double[9 * mImages];

        final double tolP = 1.e-2;
        final double tolG = 1.e-3;

        // not using these as they are estimated in calculateLMVectorsSparsely
        //initDeltaPWithQu(outDP);
        //initDeltaCWithQu(outDC);

        // evaluation of the objective re-projection error. 
        //the sum of squares of the observed feature - projected feature in camera reference frame
        final double[] outFSqSum = new double[1];

        double[] outInitLambda = new double[1];

        double lambda = 0;
        double f;
        double fBest = Double.POSITIVE_INFINITY;

        int nIter = 0;
        while (nIter < nMaxIter) {
            // for nIter = 0:
            //     use lambda=0, evaluate objective, and get the max of diagonal of (J^T*J) as the output initLambda

            try {
                // solve for the out variables, given the initial solution and lambda=0:
                calculateLMVectorsSparsely(coordsI, coordsW,
                        imageFeaturesMap, intr, extrRVecs, extrTrans, kRadials, useR2R4,
                        outDP, outDC, outGradP, outGradC, outFSqSum, lambda, outInitLambda, useBouguetForRodrigues);

                log.info(String.format("FSqSum=%.3e", outFSqSum[0]));

            } catch (NaNException e) {
                System.err.println(e.getMessage());
                return;
            }

            // sum of the squares of the re-projection errors:
            f = outFSqSum[0];

            if (nIter == 0) {
                lambda = outInitLambda[0];
                log.info(String.format("max diag of Hessian lambda=%.7e\n", lambda));
                // set to null to prevent re-calculating the max of diagonal of (J^T*J) again in calculateLMVectorsSparsely
                outInitLambda = null;
                fBest = f;
            }

            log.fine(String.format(
                    "(nIter=0) lambda=%.3e F=%.3e\n  dC=%s\n  gradC=%s\n\n",
                    lambda, outFSqSum[0],
                    FormatArray.toString(outDC, "%.3e"),
                    FormatArray.toString(outGradC, "%.3e")));
            log.fine(String.format(
                    "dP=%s\n  gradP=%s\n", FormatArray.toString(outDP, "%.3e"),
                    FormatArray.toString(outGradP, "%.3e")
            ));

            // make a test update of the parameters and calculate fsqsum for those test params
            BlockMatrixIsometric intrTest = intr.copy();
            double[][] extrRVecsTest = MatrixUtil.copy(extrRVecs);
            double[][] extrTransTest = MatrixUtil.copy(extrTrans);
            double[][] kRadialsTest = MatrixUtil.copy(kRadials);
            double[][] coordsWTest = MatrixUtil.copy(coordsW);
            /*rotation elements are indexes 0, 1, 2 of dC
            translation elements are indexes 3,4,5 of dC
            focus, radial1, radial2 are indexes 6,7,8 of dC
            3D WCS of features points are indexes 0, 1, 2 of dP
            */
            updateTranslation(extrTransTest, outDC);
            updateRotationVectors(extrRVecsTest, outDC);
            updateIntrinsic(intrTest, outDC);
            updateRadialDistortion(kRadialsTest, outDC);
            updateWorldC(coordsWTest, outDP);

            double fTest = calcReprojectionErrors(coordsI, coordsW, intrTest, extrRVecsTest, extrTransTest,
                    kRadialsTest, useR2R4, useBouguetForRodrigues);

            // calc gain ratio before potentially changing f;
            double gainRatio = calculateGainRatio(fTest, fBest, outDC, outDP, lambda, outGradC, outGradP, eps);

            // calc variables needed for stopping conditions before potentially updating the parameters
            boolean deltaCStop = isNegligible(outDC, tolP);
            boolean gradCStop = isNegligible(outGradC, tolG);
            boolean deltaPStop = isNegligible(outDP, tolP);
            boolean gradPStop = isNegligible(outGradP, tolG);

            log.info(String.format("%d) FBest=%.3e, F=%.3e, FTest=%.3e accept=%b  gain ratio=%.3e",
                    nIter, fBest, f, fTest, (fTest < fBest), gainRatio));

            boolean accept = false;
            if (fTest < fBest) {
                accept = true;
                fBest = fTest;
                //set params = test params
                intr.set(intrTest);
                MatrixUtil.copy(extrRVecsTest, extrRVecs);
                MatrixUtil.copy(extrTransTest, extrTrans);
                MatrixUtil.copy(kRadialsTest, kRadials);
                MatrixUtil.copy(coordsWTest, coordsW);
                // increase the step size by decreasing lambda
                lambda /= lambdaF;
           // } else if (nIter > 1 && fTest > f) {
                //TODO: check that this is a valid stopping condition
         //       break;
            } else {
                //log.info(String.format("%d) FSqSum=%.3e, FTestSqSum=%.3e accept=%b gain ratio=%.3e", nIter, f, fTest, accept, gainRatio));
                // decrease the step size by increasing lambda
                lambda *= lambdaF;
                lambdaF *= 2;
            }

            log.info(String.format("new lambda=%.11e\n", lambda));

            // check for convergence
            if (deltaCStop || gradCStop
                    //|| (gainRatio > 0)
            ) {
                break;
            }
            if (lambda < eps) {
                break;
            }
            ++nIter;
        }

        log.info("nIter=" + nIter);
    }

    /**
     *
     * @param coordsI
     * @param coordsW
     * @param intr
     * @param extrRVecs
     * @param extrTrans
     * @param kRadials
     * @param useR2R4
     * @param useBouguetForRodrigues if true, uses only the Bouguet algoirthms for Rodrigues rotation matrices and vectors
     * @return
     * @throws IOException
     * @throws NotConvergedException
     */
    protected double calcReprojectionErrors(double[][] coordsI, double[][] coordsW,
        BlockMatrixIsometric intr, double[][] extrRVecs, double[][] extrTrans,
        double[][] kRadials, final boolean useR2R4, boolean useBouguetForRodrigues
                                            ) throws IOException, NotConvergedException {

        int nFeatures = coordsW[0].length;
        int mImages = coordsI[0].length / nFeatures;

        if (coordsI.length != 3) {
            throw new IllegalArgumentException("coordsI.length must be 3");
        }
        if (coordsW.length != 3) {
            throw new IllegalArgumentException("coordsW.length must be 3");
        }
        if (coordsI[0].length != nFeatures * mImages) {
            throw new IllegalArgumentException("coordsI[0].length must be evenly "
                    + "divisible by nFeatures which is coordsW[0].length");
        }
        if (intr.getA().length != 3 * mImages) {
            throw new IllegalArgumentException("intr.length must be 3*mImages");
        }
        if (intr.getA()[0].length != 3) {
            throw new IllegalArgumentException("intr[0].length must be 3");
        }
        if (kRadials.length != mImages) {
            throw new IllegalArgumentException("kRadials.length must be equal "
                    + "to the number of cameras.");
        }
        if (kRadials[0].length != 2) {
            throw new IllegalArgumentException("kRadials[0].length must be 2.");
        }
        if (extrRVecs[0].length != 3) {
            throw new IllegalArgumentException("extrRVecs[0].length must be 3");
        }
        if (extrRVecs.length != mImages) {
            throw new IllegalArgumentException("extrRVecs.length must be mImages "
                    + "where mImages = coordsI[0].length/coordsW[0].length");
        }
        if (extrTrans[0].length != 3) {
            throw new IllegalArgumentException("extrTrans[0].length must be 3");
        }
        if (extrTrans.length != mImages) {
            throw new IllegalArgumentException("extrTrans.length must be mImages "
                    + "where mImages = coordsI[0].length/coordsW[0].length");
        }
        if (nFeatures < 6) {
            throw new IllegalArgumentException("need at least 6 features in an image");
        }

        double fSqSum = 0;

        double[][] auxIntr = MatrixUtil.zeros(3, 3);

        //size is [3 X 3*mImages] with each block being [3X3]
        BlockMatrixIsometric rotMatrices = createRotationMatricesFromVectors(extrRVecs, useBouguetForRodrigues);
        double[][] rotM = MatrixUtil.zeros(3, 3);

        // [2X1]
        double[] fIJ2 = new double[2];

        // i for n features, j for m images
        int i, j, k;
        for (j = 0; j < mImages; ++j) { // this is camera c in Engels pseudocode

            //intr is 3 X 3*nCameras where each block is size 3X3.
            intr.getBlock(auxIntr, j, 0);

            // get the rotation matrix rotM [3X3]
            rotMatrices.getBlock(rotM, 0, j);

            double[] omckk = Arrays.copyOf(extrRVecs[j], extrRVecs[j].length);
            double[] Tckk = Arrays.copyOf(extrTrans[j], extrTrans[j].length);

            Camera.CameraIntrinsicParameters intrinsics
                    = new Camera.CameraIntrinsicParameters(auxIntr, kRadials[j], useR2R4);

            // since the derivatives aren't needed for projection, can just use rigid motion and then transform to image frame
            CameraPose.ProjectedPoints pp = CameraPose.bouguetProjectPoints2(coordsW, omckk, Tckk, intrinsics,
                    useBouguetForRodrigues);
            double[][] x = pp.xEst;  //[2 X n] // these are in image reference frame
            //CameraPose.ProjectedPoints pp  = CameraPose.bouguetRigidMotion(coordsW, extrRVecs[j], extrTrans[j]);
            //double[][] x = pp.xEst;
            //if (useCameraFrame == 0) {
            //    x = Camera.cameraToPixelCoordinates(x, intrinsics);
            //}

            double[][] xkk = MatrixUtil.copySubMatrix(coordsI, 0, 1, nFeatures * j, nFeatures * (j + 1) - 1);

            //[2 X n]
            double[][] fj = MatrixUtil.pointwiseSubtract(xkk, x);

            for (i = 0; i < nFeatures; ++i) {
                //aIJ [2X9]  camera partial derivs (rot, trans, f, k0, k1)
                //bIJ [2X3] point partial derivs (dXW)
                for (k = 0; k < 2; ++k) {
                    fIJ2[k] = fj[k][i];
                }
                fSqSum += MatrixUtil.innerProduct(fIJ2, fIJ2);
            }
        }
        return fSqSum;
    }

    /**
     * solve for bundle adjustment data structures needed by the Levenberg-Marquardt
     * algorithm to refine the intrinsic and extrinsic camera parameters.
     * The algorithms uses the sparse structure of the jacobian to reduce
     * the computation time and memory needed.
     * The code needs initial parameter estimates of intrinsic and extrinsic
     * camera parameters (in intr, extrRot, kRadial, and extrTrans).
     The code returns results in outGradP, outGradC,
     outFSqSum, outDP, outDC for refined parameters
     and the gradient, objective and parameter update steps needed by
     code such as Levenberg-Marquardt.
     NOTE: the code does not update the intrinsic and extrinsic camera parameters, allowing
     the L-M algorithm to handle that.
     The runtime complexity is ~ O(nFeatures * mImages^2).
     Assumptions used in forming the partial derivatives of the intrinsic camera parameters
     are no skew, focal length along x is the same as focal length along y, square pixels.
     Cholesky decomposition is used with forward and back substitution
     to avoid inverting the reduced camera matrix
     and to half the runtime compared to L-U decomposition.
     Note that there can be more than one camera and should be 6 of more features per image
     (3 for rot, 3 for trans) and among those, need 3 per camera for the intrinsic parameters
     and 2 or more vantage points for the point parameters (no reference for these
     numbers, just a rough estimate from counting the number of unknowns).

     Also note that the code uses the Jacobian J = [J_P J_C] following
     Engels et al., which is reversed from the Lourakis Jacobian J = [J_c J_P].
     Qu Jacobian (and hence reduced camera matrix are consistent with Lourakis.
     See details in doc/bundle_adjustment.pdf
     <pre>
     References:

     additional information is present in directory doc as "bundle_adjustment.pdf"
     and "Levenberg-Marquardt_notes.pdf"

     * http://users.ics.forth.gr/~lourakis/sba/PRCV_colloq.pdf
     lecture by Lourakis  “Bundle adjustment gone public”

     Engels, Stewenius, Nister 2006, “Bundle Adjustment Rules”

     Bill Triggs, Philip Mclauchlan, Richard Hartley, Andrew Fitzgibbon.
     Bundle Adjustment – A Modern Synthesis.
     International Workshop on Vision Algorithms,
     Sep 2000, Corfu, Greece. pp.298–372,
     10.1007/3-540-44480-7_21 . inria-00548290

     Zhongnan Qu's master thesis, "Efficient Optimization for Robust Bundle
     Adjustment", 2018 Technical University of Munich

     Chen, Chen, & Wang 2019, "Bundle Adjustment Revisited"

     T. Barfoot, et al., "Pose estimation using linearized rotations and
     quaternion algebra", Acta Astronautica (2010), doi:10.1016/j.actaastro.2010.06.049
     -- using the rotation and translation update details.
     -- one of the 2 examples is interesting for the problem of pose for
     a pair of stereo-images.  it also uses cholesky factoring of block
     sparse matrix structure.

     The partial derivatives are from Bouguet's Camera Calibration Toolbox.
     The Bouguet toolbox webpage is currently at http://robots.stanford.edu/cs223b04/JeanYvesCalib/
     and states that the source code is freely available.
     The github repositories with forked Bouguet Matlab code do not have license
     information.  Those references are
     https://github.com/fragofer/TOOLBOX_calib
     and
     https://github.com/hunt0r/Bouguet_cam_cal_toolbox
     and the code adapted from is compute_extrinsic_refine.m which depends upon other code from the toolbox.
     </pre>
     TODO: review and improve the derivatives here.  e.g. re-do porting of the Qu derivs...
     * @param coordsI the features observed in different images (in coordinates
     * of the image reference frame).  The different images may or may not be from the same camera.
     * The format of coordsI is 3 X (nFeatures*nImages). Each row should
     * have nFeatures of one image, followed by nFeatures of the next image,
    etc.  The first dimension is for the x,y, and z axes.
    Note that if a feature is not present in the image, that should be
    an entry in imageMissingFeatureMap.
     * @param coordsW the features in a world coordinate system.  The format is
     * 3 X nFeatures.  The first dimension is for the x,y, and z axes.
     * @param imageFeaturesMap an associative array holding the features
     * present in each image.  They key is the image number in coordsI
     * which is j/nFeatures where j is the index of the 2nd dimension,
     * that is coordsI[j].  The value is a set of feature numbers which are
     * missing from the image.  The feature numbers correspond to the
     * 2nd dimension indexes in coordsW.
     * @param intr the intrinsic camera parameter matrices stacked along rows
     * to make a tall double array of size [(mImages*3) X 3] where each block is
     * size 3X3.   Note that only the focus parameter is refined in this class.
     * @param extrRVecs the extrinsic camera parameter rotation euler angles
     * stacked along the 3 columns, that is the size is nImages X 3 where
     * nImages is coordsI[0].length/coordsW[0].length.  each array is size
     * 1X3.
     * @param extrTrans the extrinsic camera parameter translation vectors
     * stacked along the 3 columns, that is the size is nImages X 3 where
     * nImages is coordsI[0].length/coordsW[0].length.  each array is size
     * 1X3.
     * @param kRadials a double array wherein each row holds the
     * radial distortion coefficients k1 and k2 for an image, so the total size is [nCameras X 2].
     * NOTE: current implementation accepts values of 0 for k1 and k2.
     * @param useR2R4 useR2R4 use radial distortion function from Ma et al. 2004 for model #4 in Table 2,
    f(r) = 1 +k1*r^2 + k2*r^4 if true,
    else use model #3 f(r) = 1 +k1*r + k2*r^2.
     * @param outDP an output array holding the update values for the point parameters.
     * The length should be 3*nFeatures.
     * @param outDC an output array holding the update values for the camera parameters.
     * The length should be 9*mImages.
     * @param outGradP an output array holding the gradient covector for point parameters
     *  (-J_P^T*(x-x_hat) as the summation of bij^T*fij).  The length should be 3 * number of features.
     * This is used by the L-M algorithm to calculate the gain ratio and evaluate stopping criteria.
     * @param outGradC an output array holding the gradient covector for camera parameters
     * (-J_C^T*(x-x_hat) as the summation of aij^T*fij).
     * The length should be 9 times the number of images.
     * This is used by the L-M algorithm to calculate the gain ratio and evaluate stopping criteria.
     * @param outFSqSum an output array holding the evaluation of the objective,
     * that is the sum of squares of the observed feature - projected feature.
     * It's the re-projection error.  NOTE that this evaluation is for the
     * given parameters, not the given parameters plus the returned update steps
     * (outDC and outDP).
     * The length should be 1.
     * @param lambda the damping parameter.  upon first use, this is 0 and
     * outInitLambda is not null so that the sparse Hessian is calculated without
     * the damping term and is used to find the initial value of
     * lambda which it places in outInitLambda.  upon all subsequent uses of
     * this method, it's expected that lambda > 0 and outInitLambda is null.
     * @param outInitLambda when not null this is the output parameter holding
     * the maximum of the diagonal of j^T*J.  the array has length 1.
     * @param useBouguetForRodrigues if true, uses only the Bouguet algoirthms for Rodrigues rotation matrices and vectors
     * @throws no.uib.cipr.matrix.NotConvergedException
     * @throws java.io.IOException
     */
    protected void calculateLMVectorsSparsely(double[][] coordsI, double[][] coordsW,
      TIntObjectMap<TIntSet> imageFeaturesMap,
      BlockMatrixIsometric intr, double[][] extrRVecs, double[][] extrTrans,
      double[][] kRadials, final boolean useR2R4,
      double[] outDP, double[] outDC, double[] outGradP, double[] outGradC,
      double[] outFSqSum, final double lambda,
      double[] outInitLambda, boolean useBouguetForRodrigues) throws NotConvergedException, IOException, NaNException {

        int nFeatures = coordsW[0].length;
        int mImages = coordsI[0].length/nFeatures;

        if (coordsI.length != 3) {
            throw new IllegalArgumentException("coordsI.length must be 3");
        }
        if (coordsW.length != 3) {
            throw new IllegalArgumentException("coordsW.length must be 3");
        }
        if (coordsI[0].length != nFeatures*mImages) {
            throw new IllegalArgumentException("coordsI[0].length must be evenly "
                    + "divisible by nFeatures which is coordsW[0].length");
        }
        if (intr.getA().length != 3*mImages) {
            throw new IllegalArgumentException("intr.length must be 3*mImages");
        }
        if (intr.getA()[0].length != 3) {
            throw new IllegalArgumentException("intr[0].length must be 3");
        }
        if (kRadials.length != mImages) {
            throw new IllegalArgumentException("kRadials.length must be equal "
                    + "to the number of cameras.");
        }
        if (kRadials[0].length != 2) {
            throw new IllegalArgumentException("kRadials[0].length must be 2.");
        }
        if (extrRVecs[0].length != 3) {
            throw new IllegalArgumentException("extrRVecs[0].length must be 3");
        }
        if (extrRVecs.length != mImages) {
            throw new IllegalArgumentException("extrRVecs.length must be mImages "
                    + "where mImages = coordsI[0].length/coordsW[0].length");
        }
        if (extrTrans[0].length != 3) {
            throw new IllegalArgumentException("extrTrans[0].length must be 3");
        }
        if (extrTrans.length != mImages) {
            throw new IllegalArgumentException("extrTrans.length must be mImages "
                    + "where mImages = coordsI[0].length/coordsW[0].length");
        }
        if (outDP.length != 3*nFeatures) {
            throw new IllegalArgumentException("outDP.length must be 3*nFeatures "
                    + "where nFeatures=coordsW[0].length");
        }
        if (outDC.length != 9*mImages) {
            throw new IllegalArgumentException("outDC.length must be 9*mImages "
                    + "where mImages=coordsI[0].length/coordsW[0].length");
        }
        if (outGradP.length != 3*nFeatures) {
            throw new IllegalArgumentException("outGradP.length must be 3 * number of features");
        }
        if (outGradC.length != 9*mImages) {
            throw new IllegalArgumentException("outGradC.length must be 9 * number of images");
        }
        if (outFSqSum.length != 1) {
            throw new IllegalArgumentException("outFSqSum.length must be 1");
        }
        if (!(lambda  >= 0.)) {
            throw new IllegalArgumentException("lambda must be a positive number");
        }
        if (imageFeaturesMap == null) {
            throw new IllegalArgumentException("imageFeaturesMap cannot be null");
        }
        if (imageFeaturesMap.size() != mImages) {
            throw new IllegalArgumentException("imageFeaturesMap size must equal "
                    + "the number of images which = coordsI[0].length/coordsW[0].length");
        }
        if (nFeatures < 6) {
            throw new IllegalArgumentException("need at least 6 features in an image");
        }

        if (outInitLambda != null) {
            // this will hold the maximum of the diagonals of hPPI and hCCJ
            outInitLambda[0] = Double.NEGATIVE_INFINITY;
        }

        Arrays.fill(outGradC, 0); // [9*mImages]
        Arrays.fill(outGradP, 0); // [3*nFeatures]
        Arrays.fill(outDC, 0);    // [9*mImages]
        Arrays.fill(outDP, 0);    //[3*nFeatures]
        Arrays.fill(outFSqSum, 0);

        double[] bC = outGradC;// [9*mImages]
        double[] bP = outGradP;// [3*nFeatures]

        double[][] auxIntr = MatrixUtil.zeros(3, 3);

        //size is [3 X 3*mImages] with each block being [3X3]
        BlockMatrixIsometric rotMatrices = createRotationMatricesFromVectors(extrRVecs, useBouguetForRodrigues);
        double[] rotAux = new double[3];
        double[][] rotM = MatrixUtil.zeros(3, 3);

        //[2X9]  camera partial derivs (rot, trans, f, k0, k1)
        double[][] aIJ = MatrixUtil.zeros(2, 9);
        //[2X3] point partial derivs (dXW)
        double[][] bIJ = MatrixUtil.zeros(2, 3);
        //aka jP_I_J^T [3X2]
        double[][] bIJT = MatrixUtil.zeros(3, 2);
        //aka jC_I_J^T  [9X2]
        double[][] aIJT = MatrixUtil.zeros(9, 2);
        // aka jP^T*JP; [3X3]
        double[][] bIJsq = MatrixUtil.zeros(3, 3);
        // [2X1]
        double[] fIJ2 = new double[2];
        //aka bP [3X1]
        double[] bIJTF = new double[3];
        //aka bC [9X1]  which is aka jC_I_J^T * fIJ
        double[] aIJTF = new double[9];

        //HPC is a.k.a. J_C^T*J_P a.k.a. (W^*)^T
        //W_i_j^T are stored here as [3X9] blocks. each W_i_j^T is b_i_j^T*a_i_j
        BlockMatrixIsometric hPCBlocks /*W^T*/= new BlockMatrixIsometric(MatrixUtil.zeros(3*nFeatures, 9*mImages), 3, 9);
        double[][] auxHPC = MatrixUtil.zeros(3, 9);
        double[][] auxHPCT = MatrixUtil.zeros(9, 3);

        //HPP is a.k.a. (J_P^T*J_P) a.k.a. V
        //The diagonals, V_i, are stored here as [3X3] blocks. each V_i is a summation
        // of b_i_j^T*b_i_j over all j images.
        BlockMatrixIsometric hPPIBlocks /*VI*/= new BlockMatrixIsometric(MatrixUtil.zeros(3*nFeatures, 3), 3, 3);
        // aka V_i; a [3X3] block
        double[][] hPPI = MatrixUtil.zeros(3, 3);

        //HCC is a.k.a. J_C^T*J_C a.k.a. U^* (and in Qu thesis is variable B).
        //The diagonals, U_j, are stored here as [9X9] blocks. each U_j is a summation of a_i_j^T*a_i_j over all i features.
        // HCC is set into matrix A as it is calculated.
        // storing hCCJBlocks in mA, while populating mA with only HCC.  later will add the negative rightsize of mA to mA
        BlockMatrixIsometric hCCJBlocks /*UJ*/= new BlockMatrixIsometric(MatrixUtil.zeros(9*mImages, 9), 9, 9);
        double[][] auxHCCJ = MatrixUtil.zeros(9, 9);

        // i for n features, j for m images
        int i, j, k;
        for (j = 0; j < mImages; ++j) { // this is camera c in Engels pseudocode

            //intr is 3 X 3*nCameras where each block is size 3X3.
            intr.getBlock(auxIntr, j, 0);

            // get the rotation matrix rotM [3X3]
            rotMatrices.getBlock(rotM, 0, j);

            double[] omckk = Arrays.copyOf(extrRVecs[j], extrRVecs[j].length);
            double[] Tckk = Arrays.copyOf(extrTrans[j], extrTrans[j].length);

            Camera.CameraIntrinsicParameters intrinsics
                    = new Camera.CameraIntrinsicParameters(auxIntr, kRadials[j], useR2R4);
            CameraPose.ProjectedPoints pp;
            if (useBouguetDerivs) {
                pp = CameraPose.bouguetProjectPoints2(coordsW, omckk, Tckk, intrinsics, useBouguetForRodrigues);
            } else {
                pp = quProjectPoints(coordsW, omckk, Tckk, intrinsics, useBouguetForRodrigues);
            }
            double[][] x = pp.xEst;  //[2 X n] // these are in image reference frame
            double[][] dxdom = pp.dxdom; // [2*n X 3]
            double[][] dxdT = pp.dxdT; // [2*n X 3]
            double[][] dxdF = pp.dxdF;//[2*n X 2]
            double[][] dxdK = pp.dxdK;//[2*n X 4]
            double[][] dP;
            if (useBouguetDerivs) {
                dP = MatrixUtil.multiply(dxdT, rotM); //[2*n X 3]
            } else {
                dP = pp.dxdX;
            }

            double[][] xkk = MatrixUtil.copySubMatrix(coordsI, 0, 1, nFeatures*j, nFeatures*(j + 1) - 1);

            //[2 X n]
            double[][] fj = MatrixUtil.pointwiseSubtract(xkk, x);

            for (i = 0; i < nFeatures; ++i) {
                //aIJ [2X9]  camera partial derivs (rot, trans, f, k0, k1)
                //bIJ [2X3] point partial derivs (dXW)
                for (k = 0; k < 2; ++k) {
                    System.arraycopy(dxdom[i*2 + k], 0, aIJ[k], 0, 3 );
                    System.arraycopy(dxdT[i*2 + k], 0, aIJ[k], 3, 3 );
                    aIJ[k][6] = dxdF[i*2 + k][0];
                    aIJ[k][7] = dxdK[i*2 + k][0];
                    aIJ[k][8] = dxdK[i*2 + k][1];
                    System.arraycopy(dP[i*2 + k], 0, bIJ[k], 0, 3 );
                    fIJ2[k] = fj[k][i];
                }
                //System.out.printf("AIJ=%s\n", FormatArray.toString(aIJ,"%.3e"));

                //bIJ^T; [3X2]  aka jP^T
                MatrixUtil.transpose(bIJ, bIJT);
                //aIJ^T; [9X2] aka jC^T
                MatrixUtil.transpose(aIJ, aIJT);

                outFSqSum[0] += MatrixUtil.innerProduct(fIJ2, fIJ2);

                /*
                gradP = bP = -JP^T * F  [3n X 1]
                ———————————————-----------------
                -B11T*F11-B12T*F12-B13T*F13
                -B21T*F21-B22T*F22-B23T*F23
                -B31T*F31-B32T*F32-B33T*F33
                -B41T*F41-B42T*F42-B43T*F43

                where each row is [3X2]*[2X1] = [3X1]
                */
                //bIJTF =  bIJT * fIJ;// [3X2]*[2X1] = [3X1]
                MatrixUtil.multiplyMatrixByColumnVector(bIJT, fIJ2, bIJTF);
                for (k = 0; k < 3; ++k) {
                    // i is the current feature
                    bP[i*3 + k] -= bIJTF[k];
                }

                /*
                gradC = bC = -JC^T * F =  [9m X 1]
                ———————————————-
                -A11T*F11-A21T*F21-A31T*F31-A41T*F41
                -A12T*F12-A22T*F22-A32T*F32-A42T*F42
                -A13T*F13-A23T*F23-A33T*F33-A43T*F43

                where each row is [9X2][2X1] = [9X1]
                */
                //aIJTF = aIJT * fIJ.  [9X2]*[2X1]=[1X9]
                MatrixUtil.multiplyMatrixByColumnVector(aIJT, fIJ2, aIJTF);
                for (k = 0; k < 9; ++k) {
                    // j is the current image
                    bC[j*9 + k] -= aIJTF[k]; // bc is [9*mImages]
                }

                // compute block (i,j) of hPC as hPC=jPTJC [3X9]
                //hPC[i][j] = bIJT * aIJ;
                MatrixUtil.multiply(bIJT, aIJ, auxHPC);
                hPCBlocks.setBlock(auxHPC, i, j);

                // bIJsq = bij^T * bij = [3X2]*[2X3] = [3X3]
                MatrixUtil.multiply(bIJT, bIJ, bIJsq);

                // HPP_i, aka V_i: for feature i, sum over all images. [3X3]
                // fetch existing block for feature i, add BIJsq to it, update the block
                hPPIBlocks.addToBlock(bIJsq, i, 0);

                //[9X9]
                //each HCC_j a.k.a. U_j, is a summation of a_i_j^T*a_i_j over all i features
                createATransposedTimesA(aIJ, auxHCCJ);
                hCCJBlocks.addToBlock(auxHCCJ, j, 0);
            } // end loop i over features
        } // end loop j over images

        log.info(String.format("outGradC=\n%s\n", FormatArray.toString(outGradC, "%.3e")));
        log.info(String.format("outGradP=%.3e, %.3e, %.3e, %.3e, %.3e, %.3e, ...",
                outGradP[0], outGradP[1], outGradP[2], outGradP[3], outGradP[4], outGradP[5]));

        // augment the diagonals of HPP and HCC by the dampening term.

        // Section 2 of Engels at al., between  eqns (4) and (5)
        // augment H_PP (and H_CC) by damping term
        //  (J_P)^T*J_P + lambda*diag((J_P)^T*J_P).
        // Each H_PP_i and each H_CC_j are the diagonal blocks of J^T*J.
        // Solomon's "Numerical Algorithms" Section 12.1.2:
        //  J^T*J is positive semi-definite and so J^T*J + λ*I_n×n must be positive definite.
        for (j = 0; j < mImages; ++j) {
            hCCJBlocks.getBlock(auxHCCJ, j, 0);
            if (outInitLambda != null) {
                if (maxDiag(auxHCCJ, outInitLambda)) {
                    log.info(String.format("max of diagonal blocks of HCC (aka U): new lambda=%.7e\n", outInitLambda[0]));
                }
            }
            for (k = 0; k < auxHCCJ.length; ++k) {
                auxHCCJ[k][k] += lambda;
            }
            hCCJBlocks.setBlock(auxHCCJ, j, 0);
        }

        for (i = 0; i < nFeatures; ++i) {
            hPPIBlocks.getBlock(hPPI, i, 0);
            if (outInitLambda != null) {
                // find maximum of the diagonal of hPP.  the diagonal of HPP is each [3X3] block hPPI.
                if (maxDiag(hPPI, outInitLambda)) {
                    log.info(String.format("max diag of hPPI (aka V_i): new lambda=%.7e\n", outInitLambda[0]));
                }
            }
            for (k = 0; k < 3; ++k) {
                hPPI[k][k] += lambda;
            }
            hPPIBlocks.setBlock(hPPI, i, 0);
        }

        // solve for dC and dP (== gradC and gradP, respectively)

        /*solve each line separately in the augmented blocks of Schur decomposition:
           |dP + (HPP^-1*HPC)*dC     | = |HPP^-1*bP           |
           |(-HPCT*HPP^-1*HPC+HCC)*dC|   |-HPCT*HPP^-1*bP + bC|

         Solving the 2nd line 1st for dC:
         mA = HCC - HPCT*HPP^-1*HPC
         vB = bC - HPCT*HPP^-1*bP
         */

        //calc HPCT * HPP^-1 which is in mA and vB
        // = W * V^-1
        //[9m X 3n]           [3n X 3n]
        //W11 W21 W31 W41  *  V1^-1 0     0     0
        //W12 W22 W32 W42     0     V2^-1 0     0
        //W13 W23 W33 W43     0     0     V3^-1 0
        //                    0     0     0     V4^-1
        //W*V^-1=       [9mX3n]
        //W11*V1^-1   W21*V2^-1   W31*V3^-1   W41*V4^-1
        //W12*V1^-1   W22*V2^-1   W32*V3^-1   W42*V4^-1
        //W13*V1^-1   W23*V2^-1   W33*V3^-1   W43*V4^-1
        //each block is [9X3]
        /*
        i=1:nFeatures
           invVI = hPPIInv
           j=1:mImages
               [row J, col I] = (hPC[i][j])^T * invVI
         */
        //tPC = HPC^T*HPP^-1 = W*V^-1 [9m X 3n]
        BlockMatrixIsometric tPCBlocks = new BlockMatrixIsometric(MatrixUtil.zeros(9*mImages, 3*nFeatures),
                9, 3);
        double[][] tPC = MatrixUtil.zeros(9, 3);

        //tP = HPP^-1*bP = V^-1*bP [3n X 1] or [1n X 3] row format
        // blocks: n rows and 1 column
        BlockMatrixIsometric tPRowBlocks = new BlockMatrixIsometric(
                MatrixUtil.zeros(nFeatures, 3), 1, 3);
        double[] tPI = new double[3];

        // matrix A is the reduced camera matrix a.k.a. Schur complement. [9m X 9m]; [mXm] block matrix with  blocks [9x9]
        // U_J is stored it in alone until i and j loops complete the first time, then
        //     the rest of mA is subtracted in.
        BlockMatrixIsometric mA = new BlockMatrixIsometric(MatrixUtil.zeros(9*mImages, 9*mImages), 9, 9);

        //aka (V_i)^-1; a [3X3] block
        double[][] invHPPI = null;
        // for each feature i
        double[] bPI = new double[3];

        // calc tPC and tP entirely
        for (i = 0; i < nFeatures; ++i) {
            //calc tP = HPP^-1*bP = V^-1*bP and set into tPBlocks(i,0) += invVI*(Σ_j(-BIJT*F1J))
            //invHPPI aka V^-1 is [3X3]
            hPPIBlocks.getBlock(hPPI, i, 0);
            invHPPI = MatrixUtil.pseudoinverseRankDeficient(hPPI);

            System.arraycopy(bP, i*3, bPI, 0, 3);
            // [3X3][3X1]=[3X1]
            MatrixUtil.multiplyMatrixByColumnVector(invHPPI, bPI, tPI);
            tPRowBlocks.addToRowBlock(tPI, i, 0);

            for (j = 0; j < mImages; ++j) { // this is camera c in Engels pseudocode
                //TODO consider how to handle feature not present in image here
                if (!imageFeaturesMap.get(j).contains(i)) {
                    continue;
                }
                //calc tPC = HPC^T*HPP^-1 = W*V^-1 and set into tPCBlocks(j,i)=WIJ*invVI//[9X3]
                // auxHPC is [3X9]
                // tPC block is [9X3][3X3]=[9X3]
                hPCBlocks.getBlock(auxHPC, i, j);
                MatrixUtil.transpose(auxHPC, auxHPCT);
                MatrixUtil.multiply(auxHPCT, invHPPI, tPC);
                tPCBlocks.setBlock(tPC, j, i);
            } // end loop over j images
        } // end loop over i images

        // mA = HCC - HPCT*HPP^-1*HPC = U - W*(V^-1)*W^T
        //    = HCC - tPC * HPC
        // first term of mA is the diagonal blocks of HCC
        for (j = 0; j < mImages; ++j) {
            hCCJBlocks.getBlock(auxHCCJ, j, 0);
            mA.setBlock(auxHCCJ, j, j);
        }
        // the rest of mA subtracts from itself tPC times HPC
        double[][] tmp = MatrixUtil.multiply(tPCBlocks.getA(), hPCBlocks.getA());
        MatrixUtil.pointwiseSubtract(mA.getA(), tmp, mA.getA());

        tmp = null;

        //vB = bC - HPCT*HPP^-1*bP   = bC - W*(V^-1)*bP
        //   = bC - HPCT * tP
        // [9m X 1] - [9*mImages X 3*nFeatures] [3*nFeatures X 1]

        //stack tPRowBlocks.getA() along rows
        double[] tmp3 = MatrixUtil.stack(MatrixUtil.transpose(tPRowBlocks.getA()));// nFeatures X 3 => 3 X nFeatures => 3*nFeatures X 1
        tmp3 = MatrixUtil.multiplyMatrixByColumnVector(MatrixUtil.transpose(hPCBlocks.getA()), tmp3);

        //mImages*9
        double[] vB = MatrixUtil.subtract(bC, tmp3);

        tmp3 = null;

        // at this point, we have calculated mA and vB

        /* TODO: (optional) Fix gauge by freezing coordinates and thereby reducing
            the linear system with a few dimensions.

           ** Section 9 of Triggs et al. 2000,
           "Bundle Adjustment – A Modern Synthesis"
               "Section 9 returns to the theoretical issue of gauge freedom
                (datum deficiency), including the theory of inner constraints."
                NLK: see MASKS Table 6.5.

            Section 9.2.1, Up vector selection, of Szeliski 2010

            Triggs 1998, "Optimal estimation of matching constraints.
              3D Structure from Multiple Images of Large-scale Environments SMILE’98,
              Lecture Notes in Computer Science
              (see Section 3.1 page 8
                   "the gauge freedom is the 3 d.o.f. choice of plane."

            N Snavely, SM Seitz, R Szeliski - 2008
            "Skeletal graphs for efficient structure from motion"

            Forstner & Wrobel refer to it as "Free Block Adjustment"

            ** Daniel D. Morris, Kenichi Kanatani and Takeo Kanade,
            "Gauge Fixing for Accurate 3D Estimation"

           Also, in this project, can see it as fixing the exrinsic parameters
              of the first camera to rotation = I and translation=0.
           Also in this project, Reconstruction.java:
              see implementation of metric constraints, after the comments
              Fig 3.1 of Tomasi & Kanade 1991 or Fig 2. of Belongie lecture notes
              Belongie Section 16.4.4 (c)
              See Step 3 - Metric Constraints

           reasons to fix the gauge:
              -- decrease drift in location accuracy
              -- smaller covariance

        gauge fix not yet included here.
        */

        // cholesky decompostion to solve for dC in mA*dC=vB
        // (using the sparsity of upper and lower triangular matrices results in
        //    half the computation time of LU decomposition in comparison)

        // mA is square [mImages*9, mImages*9]
        //    but not necessarily symmetric positive definite needed by the
        //    Cholesky decomposition, so need to find the nearest or a nearest
        //    symmetric positive definite.

        log.fine(String.format("mA=%s\n", FormatArray.toString(mA.getA(), "%.3e")));
        log.fine(String.format("vB=%s\n", FormatArray.toString(vB, "%.3e")));

        boolean useInv = false;

        double eps = 1.e-11;

        try {
            // this method attempts to find the nearest symmetric positive *definite* matrix to A:
            double[][] aPSD = MatrixUtil.nearestPositiveSemidefiniteToA(mA.getA(), eps);
            /*    double[][] b = MatrixUtil.nearestSymmetricToA(mA.getA());
                EVD evdB = EVD.factorize(new DenseMatrix(b));
                SVD svdB = SVD.factorize(new DenseMatrix(b));
            */
            DenseCholesky chol = new DenseCholesky(aPSD.length, false);
            chol = chol.factor(new LowerSPDDenseMatrix(new DenseMatrix(aPSD)));
            LowerTriangDenseMatrix _cholL = chol.getL();

            //[mImages*9, mImages*9]
            double[][] cholL = Matrices.getArray(_cholL);
            double[][] cholLT = MatrixUtil.transpose(cholL);

            log.fine(String.format("cholL=\n%s\n", FormatArray.toString(cholL, "%.3e")));
            log.fine(String.format("cholL*LT=\n%s\n", FormatArray.toString(
                    MatrixUtil.multiply(cholL, cholLT), "%.3e")));

            /* avoid inverting A by using Cholesky decomposition w/ forward and backward substitution to find dC.
            mA * dC = vB
            mA = L ﹡ L^* as Cholesky decomposition of A

            L ﹡ L^* * dC = vB

            let y = (L^* * dC)
            L ﹡ (y) = vB ==> solve for y via forward subst
            returning to y = (L^* * dC), can solve for dC via backward subst
            */

            double[] yM = MatrixUtil.forwardSubstitution(cholL,  vB);
            // temporary exit until find reasons for very large numbers in some
            //   of the arrays
            if (hasNaN(yM)) {
                throw new NaNException("Errors due to unusually large numbers");
            }
            // [[mImages*9 X mImages*9] * x = [mImages*9] ==> x is length mImages*9
            // x is dC
            MatrixUtil.backwardSubstitution(cholLT, yM, outDC);
            log.fine(String.format("yM=%s\n", FormatArray.toString(yM, "%.3e")));
            log.info(String.format("outDC from forward, backward substitution=\n%s\n", FormatArray.toString(outDC, "%.3e")));

        } catch (Throwable t) {
            // cholesky decomp of nearest psd to a failed.
            // use inverse
            useInv = true;
        }
        if (useInv) {
            // [mImages*9 X mImages*9] * y = [mImages X 9]
            // length is vB.length is [mImages*9 X 1]
            double[][] _mInv = MatrixUtil.pseudoinverseFullColumnRank(mA.getA());
            MatrixUtil.multiplyMatrixByColumnVector(_mInv, vB, outDC);
            log.info(String.format("outDC from pInv(M) * B =\n%s\n", FormatArray.toString(outDC, "%.3e")));
        }

        // tPC = HPC^T*(HPP^-1)
        // tPC^T = (HPP^-1) * HPC
        // calc outDP:  dP = invHPPI * bPI - invHPPI * HPC * dC
        //              dP = tP            - tPC^T * dC // [3nX1] - [3nX9m]*[9*mImagesX1]
        //                                              // = [3*n_features X 1]

        MatrixUtil.multiplyMatrixByColumnVector(MatrixUtil.transpose(tPCBlocks.getA()), outDC, outDP);
        for (i = 0; i < nFeatures; ++i) {
            tPRowBlocks.getRowBlock(tPI, i, 0);  // tPI is length 3
            // tPI - next 3 items of outDP
            for (k = 0; k < tPI.length; ++k) {
                outDP[i*3 + k] = tPI[k] - outDP[i*3 + k];
            }

            // Engels: compute updated point (i.e. world coord features)
            // NOTE: for this class, the updates are handled by the invoker of
            //       this method.  The updated parameters are given to the code so
            //       that when outFSqSum is calculated above, it is using
            //       the updated parameters and world coords.
        } // end loop over feature i

        // outDP can be large, so print just a few
        log.info(String.format("outDP=%.3e, %.3e, %.3e, %.3e, %.3e, %.3e, ...",
                outDP[0], outDP[1], outDP[2], outDP[3], outDP[4], outDP[5]));
        log.fine(String.format("outDP=%s\n", FormatArray.toString(outDP, "%.3e")));

        //outGradC
        //outGradP
        //outDC
        //outDP
        //outFSqSum

    }

    /**
    // for one image
         * @param useBouguetForRodrigues if true, uses only the Bouguet algoirthms for Rodrigues rotation matrices and vectors
*/
    private CameraPose.ProjectedPoints quProjectPoints(double[][] coordsW, double[] om, double[] t,
        Camera.CameraIntrinsicParameters intrinsics, boolean useBouguetForRodrigues) {

        int nFeatures = coordsW[0].length;

        double[] xWI = new double[3];
        double[] xWCI = new double[3];
        double[] xWII;

        double[][] rotM;
        if (useBouguetForRodrigues) {
            rotM = Rotation.createRodriguesFormulaRotationMatrix(om);
        } else {
            rotM = Rotation.createRodriguesFormulaRotationMatrix(om);
        }

        double[] rotAux = new double[3];

        AuxiliaryArrays aa = new AuxiliaryArrays();
        //[2X9]  camera partial derivs (rot, trans, f, k0, k1)
        double[][] aIJ = MatrixUtil.zeros(2, 9);
        //[2X3] point partial derivs (dXW)
        double[][] bIJ = MatrixUtil.zeros(2, 3);

        double[][] intr = intrinsics.getIntrinsic();
        double k1;
        double k2;
        if (intrinsics.getRadialDistortionCoeffs() != null) {
            k1 = intrinsics.getRadialDistortionCoeffs()[0];
            k2 = intrinsics.getRadialDistortionCoeffs()[1];
        } else {
            k1 = 0;
            k2 = 0;
        }

        CameraPose.ProjectedPoints p = new CameraPose.ProjectedPoints();

        //[2 X n] projected points xEst = R*X+T, where R = rodrigues(om), X is world coordinates of object, and T is translation
        p.xEst = MatrixUtil.zeros(2, nFeatures);
        //[2*n X 3] derivatives of XP w.r.t. rotation vector om
        p.dxdom = MatrixUtil.zeros(2*nFeatures, 3);
        //[2*n X 3] derivatives of XP w.r.t. translation vector
        p.dxdT = MatrixUtil.zeros(2*nFeatures, 3);
        //[2*n X 2] derivatives of XP w.r.t. camera focal length
        p.dxdF = MatrixUtil.zeros(2*nFeatures, 2);
        //[2*n X 2] derivatives of XP w.r.t. camera principal point
        p.dxdC = MatrixUtil.zeros(2*nFeatures, 2);
        //[2*n X 4] derivatives of XP w.r.t. camera distortion coefficients
        p.dxdK = MatrixUtil.zeros(2*nFeatures, 4);
        //[2*n X 1] derivatives of XP w.r.t. camera skew coefficient between x and y pixel
        p.dxdAlpha = new double[2*nFeatures];
        //[2*n X 3] derivatives of XP w.r.t. the real world point
        p.dxdX = MatrixUtil.zeros(2*nFeatures, 3);

        int i, k;
        for (i = 0; i < nFeatures; ++i) {
            //populate xWI; extract the world feature.  size [1X3]
            MatrixUtil.extractColumn(coordsW, i, xWI);

            //R * X + t
            Camera.worldToCameraCoordinates(xWI, rotM, t, rotAux, xWCI);

            aIJBIJ(xWI, xWCI, intr, k1, k2, om, rotM, t, aa, aIJ, bIJ);

            //transpose to image frame, includes radial distortion corrections
            xWII = Camera.cameraToPixelCoordinates(xWCI, intrinsics);
            p.xEst[0][i] = xWII[0];
            p.xEst[1][i] = xWII[1];

            //set aIJ and bIJ into output arrays. [2X9], [2X3] into [2*n X 3]  ...
            for (k = 0; k < 3; ++k) {
                p.dxdom[i*2][k]     = aIJ[0][k];
                p.dxdom[i*2 + 1][k] = aIJ[1][k];
                p.dxdT[i*2][k]     = aIJ[0][k + 3];
                p.dxdT[i*2 + 1][k] = aIJ[1][k + 3];

                p.dxdX[i*2][k]     = bIJ[0][k];
                p.dxdX[i*2 + 1][k]     = bIJ[1][k];
            }
            p.dxdF[i*2][0] = aIJ[0][6];
            p.dxdF[i*2 + 1][1] = aIJ[1][6];

            //[2*n X 2] derivatives of XP w.r.t. camera principal point
            //p.dxdC
            //[2*n X 1] derivatives of XP w.r.t. camera skew coefficient between x and y pixel
            //p.dxdAlpha

            //[2*n X 4] derivatives of XP w.r.t. camera distortion coefficients
            p.dxdK[i*2][0] = aIJ[0][7];
            p.dxdK[i*2 + 1][0] = aIJ[1][7];
            p.dxdK[i*2][1] = aIJ[0][8];
            p.dxdK[i*2 + 1][1] = aIJ[1][8];
        }

        return p;
    }

    /**
     * the partial derivative of the
     * final 2D re-projected coordinates of the i-th feature point
     * w.r.t. 2D coordinates of perspective projection of the i-th feature point.
     * Defined in Qu 2018 eqn (3.12).
     * 
     * @param xWCNI a world point projected to the camera reference frame and
     * normalized by it's last coordinate.
     * xWCI = column i of coordsW transformed to camera coordinates; 
     * xWCNI = xWCI/xWCI[2];
     * @param intr
     * @param k1 radial distortion coefficient 1
     * @param k2 radial distortion coefficient 2
     * @param out output array of size [2X2]
     */
    void pdCpIJCIJ(double[] xWCNI, double[][] intr,
        double k1, double k2, double[][] out) {
        
        if (out.length != 2 || out[0].length != 2) {
            throw new IllegalArgumentException("out size must be 2X2");
        }
        
        double x = -xWCNI[0];
        double y = -xWCNI[1];
        
        double x2 = x*x;
        double x4 = x2*x2;
        double y2 = y*y;
        double y4 = y2*y2;
        
        double f1 = intr[0][0];
        
        double pdxx = f1 * (1 + k1*(3*x2 + y2) + k2*(5*x4 + y4 + 6*x2*y2));
        double pdxy = 2*f1*x*y*(k1 + 2*k2*(x2 + y2));
        double pdyx = pdxy;
        double pdyy = f1*(k1*(3*y2 + x2) + k2*(5*y4 + x4 + 6*x2*y2));
        
        out[0][0] = pdxx;
        out[0][1] = pdxy;
        out[1][0] = pdyx;
        out[1][1] = pdyy;
    }
    
    /**
     * the partial derivative of the 2D coordinates of perspective projection 
     * of the i-th feature point normalized, w.r.t. to the same not normalized.
     * Defined in Qu 2018 eqn (3.16).
     * 
     * @param xWCI a world point projected to the camera reference frame.
     * xWCI = column i of coordsW transformed to camera coordinates; 
     * @param out output array of size [2X3]
     */
    void pdCIJXWIJ(double[] xWCI, double[][] out) {
        
        if (out.length != 2 || out[0].length != 3) {
            throw new IllegalArgumentException("out size must be 2X3");
        }
        
        double x = xWCI[0];
        double y = xWCI[1];
        double z = xWCI[2];
        double z2 = z*z;
        
        out[0][0] = -1./z;
        out[0][1] = 0;
        out[0][2] = x/z2;
        
        out[1][0] = 0;
        out[1][1] = -1./z;
        out[1][2] = y/z2;
        
    }
    
    /**
     * the partial derivative of the
     * final 2D re-projected coordinates of the i-th feature point
     * w.r.t. the intrinsic camera parameters.
     * Defined in Qu 2018 eqn (3.10).
     * 
     * @param xWCNI a world point projected to the camera reference frame and
     * normalized by it's last coordinate.
     * xWCI = column i of coordsW transformed to camera coordinates; 
     * xWCNI = xWCI/xWCI[2];
     * @param intr
     * @param k1 radial distortion coefficient 1
     * @param k2 radial distortion coefficient 2
     * @param out output array of size [2X3]
     */
    void pdCpIJYJ(double[] xWCNI, double[][] intr,
        double k1, double k2, double[][] out) {
        
        if (out.length != 2 || out[0].length != 3) {
            throw new IllegalArgumentException("out size must be 2X3");
        }
        
        double x = xWCNI[0];
        double y = xWCNI[1];
        
        double x2 = x*x;
        double y2 = y*y;
        double r2 = x2 + y2;
        double r4 = r2*r2;
        
        double f1 = intr[0][0];
        
        double dis = 1 + k1*r2 + k2*r4;
        
        out[0][0] = dis*x;
        out[0][1] = f1*r2*x;
        out[0][2] = f1*r4 * x;
        out[1][0] = dis*y;
        out[1][1] = f1*r2*y;
        out[1][2] = f1*r4*y;
        
    }
    
    /**
     * the partial derivative of the
     * final 2D re-projected coordinates of the i-th feature point
     * w.r.t. 2D coordinates of perspective projection of the i-th feature point.
     * Defined in Qu 2018 eqns (3.28 - 3.33).
     * @param xWI the 3-D coordinates of a world scene feature.
     * @param phi rotation angle vector of length 3 in units of radians
     * @param out output array of size [3X3]
     */
    void pdXWIJPhiJ(double[] xWI, double[] phi, double[][] out) {
        
        if (out.length != 3 || out[0].length != 3) {
            throw new IllegalArgumentException("out size must be 3X3");
        }
        
        double x = xWI[0];
        double y = xWI[1];
        double z = xWI[2];
        double pX = phi[0]; // in radians
        double pY = phi[1];
        double pZ = phi[2];
        double p = Math.sqrt(pX*pX + pY*pY + pZ*pZ);
        double p2 = p*p;
        double p3 = p2*p;
        double p4 = p2*p2;
        double c = Math.cos(p);
        double s = Math.sin(p);
        
        double dXdPxx = - (pX*x*s)/p +
            ((2*pX*x + pY*y + pZ*z)*(1.-c) + (pX*pY*z - pX*pZ*y)*c)/p2 +
            ((-pX*pY*z + pX*pZ*y + pX*pX*(pY*y + pZ*z + pX*x))*s)/p3 +
            (2.*pX*pX*(pX*x + pY*y + pZ*z)*(c - 1.))/p4;
        
        double dXdPxy = (z*s - pY*x*s)/p +
            (pX*y*(1.-c) + (pY*pY*z - pY*pZ*y)*c)/p2 +
            ((pX*pX*x + pX*pY*y + pX*pZ*z + pZ*y - pY*z)*pY*s)/p3 +
            (2.*pX*pY*(pX*x + pY*y + pZ*z)*(c - 1.))/p4;
        
        double dXdPxz = (-y*s - pZ*x*s)/p +
            (pX*z*(1.-c) + (pZ*pY*z - pZ*pZ*y)*c)/p2 +
            ((pX*pX*x + pX*pY*y + pX*pZ*z + pZ*y - pY*z)*pZ*s)/p3 +
            (2.*pX*pZ*(pX*x + pZ*z + pY*y)*(c-1.))/p4;
        
        double dXdPyx = (-z*s - pX*y*s)/p +
            (pY*x*(1.-c) + (pX*pZ*x - pX*pX*z)*c)/p2 +
            ((pY*pY*y + pY*pZ*z + pY*pX*x + pX*z - pZ*x)*pX*s)/p3 +
            (2.*pY*pX*(pY*y + pX*x + pZ*z)*(c-1.))/p4;
        
        double dXdPyy = - (pY*y*s)/p +
            ((2*pY*y + pZ*z + pX*x)*(1.-c) + (pY*pZ*x - pY*pX*z)*c)/p2 +
            ((-pY*pZ*x + pY*pX*z + pY*pY*(pZ*z + pX*x + pY*y))*s)/p3 +
            (2.*pY*pY*(pY*y + pZ*z + pX*x)*(c - 1.))/p4;
        
        double dXdPyz = (x*s - pZ*y*s)/p +
            (pY*z*(1.-c) + (pZ*pZ*x - pZ*pX*z)*c)/p2 +
            ((pY*pY*y + pY*pZ*z + pY*pX*x + pX*z - pZ*x)*pZ*s)/p3 +
            (2.*pY*pZ*(pY*y + pZ*z + pX*x)*(c - 1.))/p4;
        
        double dXdPzx = (y*s - pX*z*s)/p +
            (pZ*x*(1.-c) + (pX*pX*y - pX*pY*x)*c)/p2 +
            ((pZ*pZ*z + pZ*pX*x + pZ*pY*y + pY*x - pX*y)*pX*s)/p3 +
            (2.*pZ*pX*(pZ*z + pX*x + pY*y)*(c - 1.))/p4;
     
        double dXdPzy = (-x*s - pY*z*s)/p +
            (pZ*y*(1.-c) + (pY*pX*y - pY*pY*x)*c)/p2 +
            ((pZ*pZ*z + pZ*pX*x + pZ*pY*y + pY*x - pX*y)*pY*s)/p3 +
            (2.*pZ*pY*(pZ*z + pY*y + pX*x)*(c-1.))/p4;
        
        double dXdPzz = - (pZ*z*s)/p +
            ((2*pZ*z + pX*x + pY*y)*(1.-c) + (pZ*pX*y - pZ*pY*x)*c)/p2 +
            ((-pZ*pX*y + pZ*pY*x + pZ*pZ*(pX*x + pY*y + pZ*z))*s)/p3 +
            (2.*pZ*pZ*(pZ*z + pX*x + pY*y)*(c - 1.))/p4;
        
        out[0][0] = dXdPxx;
        out[0][1] = dXdPxy;
        out[0][2] = dXdPxz;
        out[1][0] = dXdPyx;
        out[1][1] = dXdPyy;
        out[1][2] = dXdPyz;
        out[2][0] = dXdPzx;
        out[2][1] = dXdPzy;
        out[2][2] = dXdPzz;
    }
    
    /**
     * NOTE: dFdPhi and dFdT include a partial derivative that is in the image reference frame.
     *
     * NOTE: there may be a problem if using the homography matrix [r0 r1 t] to
     * transform world scene to camera coordinates as the partial derivatives here
     * are assuming the use of the 3rd column of rotation too, that is R * T.
     * 
     for aIJ creates dF/dCameraParams which are the 9 parameters of 
     extrinsic and intrinsic,
     where the 9 parameters are the Qu notation for the variables phi_j, t_j, y_j.
     for each image = 9*nImages elements (j index is used for images).
     for bIJ creates dF/dPointParams which are the 3 parameters of the world point position.
     for each feature = 3 * mFeatures elements (i index is used for features)
     * Defined in Lourakis lecture slide 10.
     * 
     * @param xWI a world scene feature.
     * xWI = column i of coordsW
     * @param xWCI xWI projected to the camera reference frame.
     * xWCI = column i of coordsW transformed to camera coordinates, but not normalize;
     * @param intr
     * @param k1 radial distortion coefficient 1
     * @param k2 radial distortion coefficient 2
     * //@param rot extrinsic camera parameter rotation matrix.
     * @param rotAngles [1 X 3] array holding euler rotation angles.
     * @param rot [3 X 3] rotation matrix which was created with
     * Rotation.createRotationZYX(rotAngles...);
     * @param trans extrinsic camera parameter translation vector.
     * @param aa a group of arrays passed in by invoking code, re-used to avoid
     * constructing more objects.  AuxiliaryArrays aa = AuxiliaryArrays().
     * @param outAIJ output array of size [2X9]
     * @param outBIJ output array of size [2X3]
     */
    void aIJBIJ(double[] xWI, double[] xWCI, double[][] intr, double k1, double k2, 
        double[] rotAngles, double[][] rot, double[] trans, AuxiliaryArrays aa,
        double[][] outAIJ, double[][] outBIJ) {
                
        if (outAIJ.length != 2 || outAIJ[0].length != 9) {
            throw new IllegalArgumentException("outAIJ size must be 2X9");
        }
        if (outBIJ.length != 2 || outBIJ[0].length != 3) {
            throw new IllegalArgumentException("outBIJ size must be 2X3");
        }
        if (aa == null) {
            throw new IllegalArgumentException("aa cannot be null");
        }
        
        double[] xWCNI = Arrays.copyOf(xWCI, xWCI.length);
        int i;
        for (i = 0; i < xWCI.length; ++i) {
            xWCNI[i] /= xWCNI[2];
        }

        // 2X2.  Qu 2018 eqn (3.12)
        double[][] dCPdC = aa.a2X2; // this is in the image reference frame
        pdCpIJCIJ(xWCNI, intr, k1, k2, dCPdC);
        
        // 2X3.  Qu 2018 eqn (3.16)
        double[][] dCdX = aa.b2X3;
        pdCIJXWIJ(xWCI, dCdX);
        
        // 2X3.  Qu 2018 eqn (3.10)
        double[][] dCPdY = aa.c2X3;  // this is in the image reference frame
        pdCpIJYJ(xWCNI, intr, k1, k2, dCPdY);
        
        // 3X3.  Qu 2018 eqns (3.28 - 3.33)
        double[][] dXdP = aa.d3X3;
        pdXWIJPhiJ(xWI, rotAngles, dXdP);
       
        //========================================
        
        // [2X3].  Qu 2018 eqn (3.35)
        double[][] dFdT = aa.e2X3;
        MatrixUtil.multiply(dCPdC, dCdX, dFdT);
        //dFdT = dCdX;// excluding the image reference frame term

        // [2X3].  Qu 2018 eqn (3.34)
        double[][] dFdPhi = aa.f2X3;
        MatrixUtil.multiply(dFdT, dXdP, dFdPhi);

        // [2X3]. Qu 2018 eqn (3.36)
        double[][] dFdY = dCPdY;
    
        // [2X3].  Qu 2018 eqn (3.37)
        double[][] dFdX = aa.h2X3;
        MatrixUtil.multiply(dFdT, rot, dFdX);
        
        if (useHomography == 1) {
            
            // replace dFdPhi and dFdT
            
            boolean useLeftHanded = true;
            double[] h = new double[9];
            populateCameraProjectionHomography(rot, trans, h, useLeftHanded);
            
            double[][] jF = MatrixUtil.zeros(2, 9);
            PNP.calculateJF(xWI, h, jF);
            
            //J_g = dh/dp where h has 9 elements, and the number of parameters
            // jG is [9X6]
            double[][] jG = PNP.calculateJG(rotAngles);
            
            // [2X9]*[9X6] = [2X6]  this is dF/dTrans and dF/dRot combined
            double[][] j = MatrixUtil.multiply(jF, jG);  
            
            /*
            J_f: each block is [2X1] for x and y:
            ∂f1/∂h1  ∂f1/∂h2  ∂f1/∂h3  ∂f1/∂h4  ∂f1/∂h5  ∂f1/∂h6  ∂f1/∂h7  ∂f1/∂h8  ∂f1/∂h9

            J_g: each block is [1X1]
            ∂h1/∂p1  ∂h1/∂p2  ∂h1/∂p3  ∂h1/∂p4  ∂h1/∂p5  ∂h1/∂p6
            ∂h2/∂p1  ∂h2/∂p2  ∂h2/∂p3  ∂h2/∂p4  ∂h2/∂p5  ∂h2/∂p6
            ∂h3/∂p1  ∂h3/∂p2  ∂h3/∂p3  ∂h3/∂p4  ∂h3/∂p5  ∂h3/∂p6
            ∂h4/∂p1  ∂h4/∂p2  ∂h4/∂p3  ∂h4/∂p4  ∂h4/∂p5  ∂h4/∂p6
            ∂h5/∂p1  ∂h5/∂p2  ∂h5/∂p3  ∂h5/∂p4  ∂h5/∂p5  ∂h5/∂p6
            ∂h6/∂p1  ∂h6/∂p2  ∂h6/∂p3  ∂h6/∂p4  ∂h6/∂p5  ∂h6/∂p6
            ∂h7/∂p1  ∂h7/∂p2  ∂h7/∂p3  ∂h7/∂p4  ∂h7/∂p5  ∂h7/∂p6
            ∂h8/∂p1  ∂h8/∂p2  ∂h8/∂p3  ∂h8/∂p4  ∂h8/∂p5  ∂h8/∂p6
            ∂h9/∂p1  ∂h9/∂p2  ∂h9/∂p3  ∂h9/∂p4  ∂h9/∂p5  ∂h9/∂p6

            j = j_F*J_g:  [2x6] and each block is [2X1]
            (∂f1/∂h1)*(∂h1/∂p1) + (∂f1/∂h2)*(∂h2/∂p1 + (∂f1/∂h3)*(∂h3/∂p1) + ...+ (∂f1/∂h9)*(∂h9/∂p1)
            (∂f1/∂h1)*(∂h1/∂p2) + (∂f1/∂h2)*(∂h2/∂p2 + (∂f1/∂h3)*(∂h3/∂p2) + ...+ (∂f1/∂h9)*(∂h9/∂p2)
            ...
            (∂f1/∂h1)*(∂h1/∂p6) + (∂f1/∂h2)*(∂h2/∂p6 + (∂f1/∂h3)*(∂h3/∂p6) + ...+ (∂f1/∂h9)*(∂h9/∂p6)
            
            where p=(thetax,thetay,thetaz,transx,transy,transz)
            */
            // transpose to right-handed system
            j = MatrixUtil.transpose(j);
            for (i = 0; i < 2; ++i) {
                System.arraycopy(dFdPhi[i], 0, outAIJ[i], 0, j[i].length);
                System.arraycopy(dFdY[i], 0, outAIJ[i], 6, dFdY[i].length);
                System.arraycopy(dFdX[i], 0, outBIJ[i], 0, dFdX[i].length);
            }
        } else {
        
            //------
            // a holds camera parameters.  it's 2X9.  phi, t, y(=f, k1, f2)
            // b hold point parameters.    it's 2X3.  x
            for (i = 0; i < 2; ++i) {
                System.arraycopy(dFdPhi[i], 0, outAIJ[i], 0, dFdPhi[i].length);
                System.arraycopy(dFdT[i], 0, outAIJ[i], 3, dFdT[i].length);
                System.arraycopy(dFdY[i], 0, outAIJ[i], 6, dFdY[i].length);
                System.arraycopy(dFdX[i], 0, outBIJ[i], 0, dFdX[i].length);
            }
        }
    }

    /**
     *
     * @param extrRVecs
     * @param useBouguetForRodrigues if true, uses only the Bouguet algoirthms for Rodrigues rotation matrices and vectors
     * @return
     */
    private BlockMatrixIsometric createRotationMatricesFromVectors(double[][] extrRVecs, boolean useBouguetForRodrigues) {
        
        //extrRVecs is mImages*[1X3]
        
        int mImages = extrRVecs.length;
        
        double[][] rot = MatrixUtil.zeros(3, 3);
        
        BlockMatrixIsometric m = new BlockMatrixIsometric(MatrixUtil.zeros(3, 3*mImages), 3, 3);

        Rotation.RodriguesRotation rr;
        int i;
        for (i = 0; i < mImages; ++i) {
            if (useBouguetForRodrigues) {
                rr = Rotation.createRodriguesRotationMatrixBouguet(extrRVecs[i]);
                m.setBlock(rr.r, 0, i);
            } else {
                m.setBlock(Rotation.createRodriguesFormulaRotationMatrix(extrRVecs[i]), 0, i);
            }

            // SO3 still?
            //double detR = MatrixUtil.determinant(rr.r);
            //double[][] orthChk = MatrixUtil.createATransposedTimesA(rr.r);
            //log.info(String.format("orthChk=\n%s\ndet=%.3e", FormatArray.toString(orthChk,"%.3e"), detR));
        }
        return m;
    }

    /**
     * gain = (f(p) - f(p + delta p)) / ell(delta p)
             where ell(delta p) is (delta p)^T * (lambda * (delta p)) + J^T * ( b - f))
       gain = (f - fTest) / ( (delta p)^T * (lambda * (delta p) + J^T * ( b - f)) )
             
     * @param fTest
     * @param f
     * @param dC steps of change for the camera parameters in array of length
     * 9*mImages.  dC contains 3 rotation, 3 translation, 3 intrinsic parameters for one
     * image, followed by the same 9 for the next image, etc.
     * @param dP steps of change for the point parameters in array of length
     * 3*nFeatures with elements ordered as follows: dX_0, dY_0, dZ_0, ... dX_n-1, dY_n-1, dZ_n-1.
     * @param lambda
     * @param gradC
     * @param gradP
     * @param eps
     * @return 
     */
    private double calculateGainRatio(double fTest, double f, double[] dC, double[] dP, double lambda,
        double[] gradC, double[] gradP, double eps) {

        // (M. Lourakis, A. Argyros: SBA: A Software Package For Generic
        // Sparse Bundle Adjustment. ACM Trans. Math. Softw. 36(1): (2009))
        //  gain ratio = ( fPrev - fNew) / ( deltaParams^T * (lambda * deltaParams + J^T*fPrev) )
        //let s = 9*mImages + 9*mImages
        //   [1Xs]         *    ([1X1]*[sX1]             + [sX1])     = [1X1]
        //(delta params)^T *  (lambda * (delta params) + gradient)
        
        double[] dParams = new double[dC.length + dP.length];
        System.arraycopy(dC, 0, dParams, 0, dC.length);
        System.arraycopy(dP, 0, dParams, dC.length, dP.length);

        double[] gradient = new double[gradC.length + gradP.length];
        System.arraycopy(gradC, 0, gradient, 0, gradC.length);
        System.arraycopy(gradP, 0, gradient, gradC.length, gradP.length);
        MatrixUtil.multiply(gradient, -1); // our gradC and gradP are g = (-J^T * f) while Lourakis uses g = (J^T * f)

        double[] denom = Arrays.copyOf(dParams, dParams.length);
        MatrixUtil.multiply(denom, lambda);
        denom = MatrixUtil.add(denom, gradient);
      
        double d = MatrixUtil.innerProduct(dParams, denom);
            
        if (Math.abs(d) < eps) {
            return Double.NEGATIVE_INFINITY;
        }

        double gain = (f - fTest)/d;

        return gain;
    }
    
    private boolean isNegligible(double[] c, double eps) {
        for (int i = 0; i < c.length; ++i) {
            if (Math.abs(c[i]) > eps) {
                return false;
            }
        }
        return true;
    }

    /**
     * update t by deltaT
     * @param extrTrans
     * @param dC steps of camera parameters in array of length 9*mImages.
     * dC contains 3 rotation, 3 translation, 3 intrinsic parameters for one
     * image, followed by the same 9 for the next image, etc.  only the
     * translation elements are used in this method.
     */
    private void updateTranslation(double[][] extrTrans, double[] dC) {
        
        // from Danping Zou lecture notes, Shanghai Jiao Tong University,
        // EE382-Visual localization & Perception, “Lecture 08- Nonlinear least square & RANSAC”
        // http://drone.sjtu.edu.cn/dpzou/teaching/course/lecture07-08-nonlinear_least_square_ransac.pdf

        // parameter perturbations for a vector are:
        //     x + delta x
        int i, j;
        
        int mImages = extrTrans.length;
                
        for (j = 0; j < mImages; ++j) {
            // vector perturbation for translation:
            for (i = 0; i < 3; ++i) {
                //translation elements are indexes 3,4,5
                extrTrans[j][i] += updateSign*dC[j*9 + (i+3)];
            }
        }
    }

    /**
     * update the focus parameters of the intrinsic matrices.
     * 
     * @param intr the intrinsic camera parameter matrices stacked along
     * rows in a double array of size [3*nCameras X 3] where each block is
     * size 3X3.
     * @param dC steps of camera parameters in array of length 9*mImages.
     * dC contains 3 rotation, 3 translation, 3 intrinsic parameters for one
     * image, followed by the same 9 for the next image, etc.  only the
     * intrinsic camera parameters are used in this method.
     */
    private void updateIntrinsic(BlockMatrixIsometric intr, 
        double[] dC) {
        
        int nImages = intr.getA().length/3;
        
        double[][] kIntr = MatrixUtil.zeros(3, 3);
        
        //NOTE: follow up on Szeliski text stating that updating the intrinsic parameters is more involved
        
        int j;
        
        // using addition for updates for now
        for (j = 0; j < nImages; ++j) {
            // focus is parameter index 6 within the 9 for each image in dC
            intr.getBlock(kIntr, j, 0);
            kIntr[0][0] += updateSign*dC[j*9 + 6];
            kIntr[1][1] += updateSign*dC[j*9 + 6];
            intr.setBlock(kIntr, j, 0);
        }
    }
        
    private void updateRadialDistortion(double[][] kRadials, double[] dC) {
        
        int nImages = kRadials.length;
                
        // using addition for updates for now
        for (int i = 0; i < nImages; ++i) {
            kRadials[i][0] += updateSign*dC[i*9 + 7];
            kRadials[i][1] += updateSign*dC[i*9 + 8];
        }
    }

    /**
     * update rotation matrix theta vectors with steps in dC.
     * @param extrRVecs the extrinsic camera parameter rotation euler angles
     * stacked along the 3 columns, that is the size is [nImages X 3].
     * @param dC steps of camera parameters in array of length 9*mImages.
     * dC contains 3 rotation, 3 translation, 3 intrinsic parameters for one
     * image, followed by the same 9 for the next image, etc.  only the
     */
    private void updateRotationVectors(double[][] extrRVecs, double[] dC) {

        // from Danping Zou lecture notes, Shanghai Jiao Tong University,
        // EE382-Visual localization & Perception, “Lecture 08- Nonlinear least square & RANSAC”
        // http://drone.sjtu.edu.cn/dpzou/teaching/course/lecture07-08-nonlinear_least_square_ransac.pdf
        // parameter perturbations for a Lie group such as rotation are:
        //     R * (I - [delta x]_x) where [delta x]_x is the skew-symmetric matrix of delta_x
        // which is the same as eqn (29c) of Barfoot, et al. 2010,

        //double[] dRVec = new double[3];

        int i, j;
                        
        for (j = 0; j < extrRVecs.length; ++j) {
            // apply perturbation to the rotation vector
            for (i = 0; i < 3; ++i) {
                extrRVecs[j][i] += updateSign*dC[j*9 + i];
            }
        }
    }
    
    private void initDeltaPWithQu(double[] deltaP) {
        
        int nFeatures = deltaP.length/3;
                
        /*Qu thesis eqn (3.38)
        
        3 delta x ~ 1e-8
        */
        int i /*features*/, j /*parameters*/;
        for (i = 0; i < nFeatures; ++i) {
            // i*3 + 0,1,2
            for (j = 0; j < 3; ++j) {
                deltaP[i*3 + j] = 1e-8;
            }
        }
    }
    
    private void initDeltaCWithQu(double[] deltaC) {
        
        int mImages = deltaC.length/9;
                
        /*Qu thesis eqn (3.38)
        
        3 delta thetas ~ 1e-8
        3 delta translation ~1e-5
        1 delta focus ~ 1
        2 delta kRadial ~ 1e-3
        */
        int i /*parameter*/, j /*image*/;
        for (j = 0; j < mImages; ++j) {
            // j*9 + 0,1,2
            for (i = 0; i < 3; ++i) {
                // delta theta
                deltaC[j*9 + i] = 1e-8;
            }
            // j*9 + 3,4,5
            for (i = 3; i < 6; ++i) {
                // delta translation
                deltaC[j*9 + i] = 1e-5;
            }
            // delta focus
            deltaC[j*9 + 6] = 1;
            // delta radial coefficients
            deltaC[j*9 + 7] = 1e-8;
            deltaC[j*9 + 8] = 1e-8;
        }
    }
    
    /**
     * update the coordinates for the features in the world scene.
     * @param coordsW the features in a world coordinate system.  The format is
     * 3 X nFeatures.  The first dimension is for the x,y, and z axes.
     * @param deltaP an array of length nFeatures * 3 holding the steps in
     * world coordinates for the features in the world scene.
     * The elements are ordered as follows: dX_0, dY_0, dZ_0, ... dX_n-1, dY_n-1, dZ_n-1. 
     */
    private void updateWorldC(double[][] coordsW, double[] deltaP) {
        int nFeatures = coordsW[0].length;
        
        int i;
        for (i = 0; i < nFeatures; ++i) {
            coordsW[0][i] += updateSign*deltaP[i*3 + 0];
            coordsW[1][i] += updateSign*deltaP[i*3 + 1];
            coordsW[2][i] += updateSign*deltaP[i*3 + 2];
        }
    }
    
    private boolean maxDiag(double[][] m, double[] outInitLambda) {
        int i;
        boolean updated = false;
        for (i = 0; i < m.length; ++i) {
            if (outInitLambda[0] < m[i][i]) {
                outInitLambda[0] = m[i][i];
                updated = true;
            }
        }
        return updated;
    }
    private boolean max(double[][] m, double[] outInitLambda) {
        int i, j;
        boolean updated = false;
        for (i = 0; i < m.length; ++i) {
            for (j = 0; j < m[i].length; ++j) {
                if (outInitLambda[0] < m[i][j]) {
                    outInitLambda[0] = m[i][j];
                    updated = true;
                }
            }
        }
        return updated;
    }

    /**
     * @param nFeatures
     * @param coordsI the features observed in different images (in coordinates 
     * of the image reference frame).  The
     * different images may or may not be from the same camera.  
     * The format of coordsI is 3 X (nFeatures*nImages). Each row should
     * have nFeatures of one image, followed by nFeatures of the next image,
       etc.  The first dimension is for the x,y, and z axes.
       Note that if a feature is not present in the image, that should be
       an entry in imageMissingFeatureMap.
     * @param intr the intrinsic camera parameter matrices stacked along rows
     * to make a tall double array of size [(mImages*3) X 3] where each block is
     * size 3X3.   Note that only the focus parameter is refined in this class.
     * @param kRadial
     * @param useR2R4
     * @return 
     * @throws no.uib.cipr.matrix.NotConvergedException 
     * @throws java.io.IOException 
     */
    private double[][] transformPixelToCamera(int nFeatures, double[][] coordsI, 
        BlockMatrixIsometric intr, double[][] kRadial, boolean useR2R4) throws NotConvergedException, IOException {
        
        int mImages = coordsI[0].length/nFeatures;
        
        double[][] c = MatrixUtil.zeros(coordsI.length, coordsI[0].length);
        double[][] x = MatrixUtil.zeros(3, nFeatures);
        double[][] xc;
        double[][] kIntr = MatrixUtil.zeros(3, 3);
        
        int i, j, k;
        for (j = 0; j < mImages; ++j) {
            MatrixUtil.copySubMatrix(coordsI, 0, 2, j*nFeatures, ((j+1)*nFeatures)-1, x);
            
            intr.getBlock(kIntr, j, 0);
            xc = Camera.pixelToCameraCoordinates(x, kIntr, kRadial[j], useR2R4);

            for (i = 0; i < 3; ++i) {
                for (k = 0; k < nFeatures; ++k) {
                    xc[i][k] /= xc[2][k];
                }
                System.arraycopy(xc[i], 0, c[i], j*nFeatures, nFeatures);
            }
        }
        return c;
    }
    
    /**
     * populate the homography matrix to transform world screen coordinates
     * to projected 2D coordinates in the camera reference frame.
     * @param rot
     * @param translation
     * @param outH 
     */
    void populateCameraProjectionHomography(double[][] rot,
        double[] translation, double[][] outH) {
        outH[0][0] = rot[0][0];
        outH[1][0] = rot[1][0];
        outH[2][0] = rot[2][0];
        outH[0][1] = rot[0][1];
        outH[1][1] = rot[1][1];
        outH[2][1] = rot[2][1];
        outH[0][2] = translation[0];
        outH[1][2] = translation[1];
        outH[2][2] = translation[2];
    }
    /**
     * populate the homography matrix to transform world screen coordinates
     * to projected 2D coordinates in the camera reference frame.
     * method follows Wetzstein, eqn 19.
     * @param rot
     * @param translation
     * @param outH 
     */
    void populateCameraProjectionHomography(double[][] rot,
        double[] translation, double[] outH, boolean useLeftHanded) {
        if (useLeftHanded) {
            // Wetzstein lecture uses convention: "looking down the negative z-axis"
            outH[0] = rot[0][0];
            outH[1] = rot[0][1];
            outH[2] = translation[0];
            outH[3] = rot[1][0];
            outH[4] = rot[1][1];
            outH[5] = translation[1];
            outH[6] = rot[2][0]; //-
            outH[7] = rot[2][1]; //-
            outH[8] = translation[2]; //-
        } else {
            outH[0] = rot[0][0];
            outH[1] = rot[1][0];
            outH[3] = rot[2][0];
            outH[4] = rot[0][1];
            outH[6] = rot[1][1];
            outH[7] = rot[2][1];
            outH[2] = translation[0];
            outH[5] = translation[1];
            outH[8] = translation[2];
        }
    }

    private boolean hasNaN(double[] b) {
        for (double a : b) {
            if (Double.isNaN(a)) {
                return true;
            }
        }
        return false;
    }

    private boolean hasNaN(double[][] a) {
        int i, j;
        for (i = 0; i < a.length; ++i) {
            for(j = 0; j < a[0].length; ++j) {
                if (Double.isNaN(a[i][j])) {
                    return true;
                }
            }
        }
        return false;
    }

    static class AuxiliaryArrays {
        final double[][] a2X2;
        final double[][] b2X3;
        final double[][] c2X3;
        final double[][] d3X3;
        final double[][] e2X3;
        final double[][] f2X3;
        final double[][] g3X3;
        final double[][] h2X3;
        final Rotation.AuxiliaryArrays aa;
        public AuxiliaryArrays() {
            a2X2 = MatrixUtil.zeros(2, 2);
            b2X3 = MatrixUtil.zeros(2, 3);
            c2X3 = MatrixUtil.zeros(2, 3);
            d3X3 = MatrixUtil.zeros(3, 3);
            e2X3 = MatrixUtil.zeros(2, 3);
            f2X3 = MatrixUtil.zeros(2, 3);
            g3X3 = MatrixUtil.zeros(3, 3);
            h2X3 = MatrixUtil.zeros(2, 3);
            aa = new Rotation.AuxiliaryArrays();
        }
    }

    private static void createATransposedTimesA(double[][] a, double[][] out) {
        if (a == null || a.length == 0) {
            throw new IllegalArgumentException("m cannot be null or empty");
        }
        int m = a.length;
        int n = a[0].length;
        if (out.length != n || out[0].length != n) {
            throw new IllegalArgumentException("out must be size a[0].length X a[0].length");
        }
        int outCol, i, j;
        double sum;
        for (i = 0; i < n; i++) {
            for (outCol = 0; outCol < n; outCol++) {
                sum = 0;
                for (j = 0; j < m; j++) {
                    sum += (a[j][outCol] * a[j][i]);
                }
                out[outCol][i] = sum;
            }
        }
    }
}
