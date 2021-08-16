package algorithms.imageProcessing.transform;

import algorithms.matrix.BlockMatrixIsometric;
import algorithms.matrix.LinearEquations;
import algorithms.matrix.MatrixUtil;
import algorithms.util.FormatArray;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.set.TIntSet;
import java.io.IOException;
import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;
import no.uib.cipr.matrix.NotConvergedException;

/**
 given intrinsic and extrinsic camera parameters, coordinates for points
 in a world reference frame, and observations of those points in one or more
 camera, return data structures needed by Levenberg-Marquardt algorithm
 in refinement of the intrinsic and extrinsic camera parameters.
 BundleAdjustment calculates partial derivatives of the parameters
 and calculates the re-projection error to form the parameter update steps,
 the gradient covector, and the evaluation of the objective (sum of squares of
 the re-projection error).
 
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
 * reconstructions without up- dating their internal structure — all of the 
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
 * TODO: consider implementing or finding an implementation of:
 * Agarwal, Snavel, Seitz, and Szeliski
 * "Bundle Adjustment in the Large"
 * European conference on computer vision, pages 29–42. Springer, 2010
* 
 * @author nichole
 */
public class BundleAdjustment {
    
    private static final Level LEVEL = Level.INFO;
    private static final Logger log;
    static {
        log = Logger.getLogger(CameraCalibration.class.getSimpleName());
        log.setLevel(LEVEL);
    }
    
    /**
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
           
     Tomasi 2007,CPS 296.1 Supplementary Lectur Notes, Duke University
     
     graph partioning:
        https://cseweb.ucsd.edu/classes/fa04/cse252c/manmohan1.pdf
        recursive partitioning w/ elimination graph and vertex cut.
        
        Triggs et al. 2000, "Bundle Adjustment – A Modern Synthesis", Section 6
        
        Skeletal graphs for efficient structure from motion
        N Snavely, SM Seitz, R Szeliski - 2008
 
     Graph partitioning in this project:
        NormalizedCuts.java which uses the Fiedler vector of the Laplacian.
        UnweightedGraphCommunityFinder.java
        
     </pre>
     TODO: add gauge fix in coordination with invoker which provides initial
     * parameter estimates.
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
     * @param extrRotThetas the extrinsic camera parameter rotation euler angles
     * stacked along the 3 columns, that is the size is [mImages X 3].  
     * each image block is size 1X3.
     * @param extrTrans the extrinsic camera parameter translation vectors
     * stacked along the 3 columns, that is the size is [nImages X 3] where
     * nImages is coordsI[0].length/coordsW[0].length.  each array is size
     * 1X3.
     * @param kRadials a double array wherein each row holds the 
     * radial distortion coefficients k1 and k2 for an image.
     * NOTE: current implementation accepts values of 0 for k1 and k2.
     * TODO: consider whether to allow option of leaving out radial distortion
     * by allowing kRadials to be null.
     * @param useR2R4 useR2R4 use radial distortion function from Ma et al. 2004 for model #4 in Table 2,
        f(r) = 1 +k1*r^2 + k2*r^4 if true,
        else use model #3 f(r) = 1 +k1*r + k2*r^2.
     * @param nMaxIter
     * 
     * @throws no.uib.cipr.matrix.NotConvergedException 
     */
    public static void solveSparsely(
        double[][] coordsI, double[][] coordsW, TIntObjectMap<TIntSet> imageFeaturesMap,
        BlockMatrixIsometric intr, double[][] extrRotThetas, double[][] extrTrans,
        double[][] kRadials, final int nMaxIter, boolean useR2R4) 
        throws NotConvergedException, IOException {
        
        int nFeatures = coordsW[0].length;
        int mImages = coordsI[0].length/nFeatures;
                                 
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
        if (extrRotThetas[0].length != 3) {
            throw new IllegalArgumentException("extrRotThetas[0].length must be 3");
        }
        if (extrRotThetas.length != mImages) {
            throw new IllegalArgumentException("extrRotThetas.length must be mImages "
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
        //    origin of the world scence is [0,0,0] or [0,0,0,1]
                
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
                  ("Bundle Ajustment – A Modern Synthesis",
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
        If a 3D (X Y Z)⊤ parametrization (or equivalently a homogeneous affine (X Y Z 1)⊤ one) 
        is used for very distant 3D points, large X, Y, Z displacements are 
        needed to change the image significantly. I.e., in (X Y Z) space the 
        cost function becomes very flat and steps needed for cost adjustment 
        become very large for distant points. 
            In comparison, with a homogeneous projective parametrization (X Y Z W)⊤, 
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
        ** ==> Affine parametrization (X Y Z 1)⊤ is acceptable for points near the 
        origin with close-range convergent camera geometries, but it is 
        disastrous for distant ones because it artificially cuts away half of 
        the natural parameter space, and hides the fact by sending the resulting 
        edge to infinite parameter values.
        **==> Instead, you should use a homogeneous parametrization (X Y Z W )⊤ 
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
        
        //TODO: consider adding constraints suggested in Seliski 2010:
        // u_0 and v_0 are close to half the image lengths and widths, respectively.
        // the angle between 2 image axes is close to 90.
        // the focal lengths along both axes are greater than 0.
       
        //factor to raise or lower lambda.  
        //   consider using the eigenvalue spacing of J^T*J (Transtrum & Sethna, "Improvements to the Levenberg-Marquardt algorithm for nonlinear least-squares minimization")
        final double lambdaF = 2;
        double eps = 1E-12;
       
         // copy the parameter data structures into test (tentative) data structures
        BlockMatrixIsometric intrTest = intr.copy();
        double[][] extrRotThetasTest = MatrixUtil.copy(extrRotThetas);
        double[][] extrTransTest = MatrixUtil.copy(extrTrans);
        double[][] kRadialsTest = MatrixUtil.copy(kRadials);
        double[][] coordsWTest = MatrixUtil.copy(coordsW);
        
        //In a single reprojection error formula, there are altogether 12 arguments 
        //   (9 camera parameters and 3 feature point positions).
        
        // update values for the point parameters (== world coordinate features)
        double[] outDP = new double[3*nFeatures];
        // updatevalues for the camera parameters
        double[] outDC = new double[9*mImages];
                
        // Qu array u for parameters is ordered: rot_0, trans_0, intr_0, ...rot_m-1, trans_m-1, intr_m-1, then x_0, ... x_n
        // but the delta parameter array for all params is ordered:
        //     dRot_0, ... dRot_m-1,  dTrans_0, ...dTrans_m-1, dIntr_0,...dIntr_m-1, dX_0, ... dX_n
        // gradient g is same length   
        
        //the gradient covector for point parameters.  used in calc gain ration and stopping
        double[] outGradP = new double[3*nFeatures];
        // the gradient covector for camera parameters.  used in calc gain ration and stopping
        double[] outGradC = new double[9*mImages];
     
        final double tolP = 1.e-3;
        final double tolG = 1.e-3;
        
        // not using these as they are estimated in calculateLMVectorsSparsely
        //initDeltaPWithQu(outDP);
        //initDeltaCWithQu(outDC);
        
        // evaluation of the objective re-projection error. 
        //the sum of squares of the observed feature - projected feature in camera reference frame
        final double[] outFSqSum = new double[1];
        
        double[] outInitLambda = new double[1];
        
        // use lambda=0, evaluate objective, and get the max of diagnonal of (J^T*J) as the output initLambda
        calculateLMVectorsSparsely(coordsI, coordsWTest,  
            imageFeaturesMap, intrTest, extrRotThetasTest, extrTransTest, kRadialsTest, useR2R4,
            outDP, outDC, outGradP, outGradC, outFSqSum, 0., outInitLambda);
        
        double lambda = outInitLambda[0];
        
        // set to null to prevent calculating the max of diagnonal of (J^T*J) 
        outInitLambda = null;
        
        // sum the squares to evalute the re-projection error:
        double fPrev = outFSqSum[0];
        
        double fTest = Double.POSITIVE_INFINITY;
        
        log.log(LEVEL,
            String.format(
            "(nIter=0) lambda=%.3e F=%.3e\n  dC=%s\n  dP=%s\n  gradC=%s\n gradP=%s\n", 
            lambda, outFSqSum[0],
            FormatArray.toString(outDC, "%.3e"), FormatArray.toString(outDP, "%.3e"), 
            FormatArray.toString(outGradC, "%.3e"), FormatArray.toString(outGradP, "%.3e")));                 
                       
        double gainRatio;
                             
        // update the features with outDP?  see step 7 of Engels
        
        int nIter = 0;
        while ((nIter < nMaxIter) && (Math.abs(fPrev - fTest) >= eps)) {
                
            nIter++;
            
            // update the test data structures by deltaP and deltaC
            // use the 2nd three elements in outDC:
            updateTranslation(extrTransTest, outDC);
            updateRotThetas(extrRotThetasTest, outDC);
            updateIntrinsic(intrTest, outDC);
            updateRadialDistortion(kRadialsTest, outDC);
            //      updateWorldC(coordsWTest, outDP);

            calculateLMVectorsSparsely(coordsI, coordsWTest,  
                imageFeaturesMap, intrTest, extrRotThetasTest, extrTransTest, kRadialsTest, useR2R4,
                outDP, outDC, outGradP, outGradC, outFSqSum, 0, outInitLambda);
        
            fTest = outFSqSum[0];
                
            gainRatio = calculateGainRatio(fTest, fPrev, outDC, outDP, lambda, 
                outGradC, outGradP, eps);
            
            log.log(LEVEL,
                String.format(
                "(nIter=%d) lambda=%.3e fPrev=%.11e fTest=%.11e diff=%.11e\n "
                + " g.r.=%.3e\ndC=%s\ndP=%s\ngradC=%s\ngradP=%s\n", 
                nIter, lambda, fPrev, fTest, (fPrev-fTest),
                gainRatio,
                FormatArray.toString(outDC, "%.11e"), FormatArray.toString(outDP, "%.11e"), 
                FormatArray.toString(outGradC, "%.11e"), FormatArray.toString(outGradP, "%.11e")));
            
            /*
            for large values of lambda, the update is a very steep descent and
            deltaP is very small.
            If the damping term is small the approach is a nearly linear problem.
            
            NOTE: the damping term is used like a factor in the perturbation
            of a symmetric matrix.  see:
                https://nhigham.com/2021/02/16/diagonally-perturbing-a-symmetric-matrix-to-make-it-positive-definite/
            */
            if (gainRatio > 0) {
                // near the minimimum, which is good.
                // decrease lambda
                lambda /= lambdaF;
                assert(fTest < fPrev);                
                fPrev = fTest;
                
                double[][] _dt = MatrixUtil.elementwiseSubtract(extrTrans, extrTransTest);
                double[][] _rt = MatrixUtil.elementwiseSubtract(extrRotThetas, extrRotThetasTest);
                double _dts = MatrixUtil.frobeniusNorm(_dt);
                double _rts = MatrixUtil.frobeniusNorm(_rt);
                System.out.printf("delta trans=%.3e, %s\n", _dts, FormatArray.toString(_dt, "%.11e"));
                System.out.printf("delta rot=%.3e\n%s\n", _rts, FormatArray.toString(_rt, "%.11e"));
                
                //copy test data structures into the original in-out datastructures
                intr.set(intrTest);
                MatrixUtil.copy(extrRotThetasTest, extrRotThetas);
                MatrixUtil.copy(extrTransTest, extrTrans);
                MatrixUtil.copy(kRadialsTest, kRadials);
                MatrixUtil.copy(coordsWTest, coordsW);
                
                // ======= stopping conditions ============
                //   step length vanishes:  deltaParams --> 0
                //   gradient of f(x) vanishes: -J^T * (b - fgp) --> 0
                //MatrixUtil.multiply(gradientCheck, -1.);            
                if (isNegligible(outDC, tolP) || isNegligible(outDP, tolP) ||
                    isNegligible(outGradC, tolG) || isNegligible(outGradP, tolG)) {
                    break;
                }
            } else {
                // increase lambda to reduce step length and get closer to 
                // steepest descent direction
                lambda *= lambdaF;
            }
            System.out.printf("new lambda=%.11e\n", lambda);           
        }        
    }
    
    /**
     * NOT YET TESTED.
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
     </pre>
     * @param coordsI the features observed in different images (in coordinates 
     * of the image reference frame).  The
     * different images may or may not be from the same camera.
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
     * @param extrRotThetas the extrinsic camera parameter rotation euler angles
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
     * TODO: consider whether to allow option of leaving out radial distortion
     * by allowing kRadials to be null.
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
     */
    protected static void calculateLMVectorsSparsely(double[][] coordsI, double[][] coordsW,
        TIntObjectMap<TIntSet> imageFeaturesMap,
        BlockMatrixIsometric intr, double[][] extrRotThetas, double[][] extrTrans,
        double[][] kRadials, final boolean useR2R4,
        double[] outDP, double[] outDC, double[] outGradP, double[] outGradC, 
        double[] outFSqSum, final double lambda,
        double[] outInitLambda) throws NotConvergedException, IOException {
        
        int nFeatures = coordsW[0].length;
        int mImages = coordsI[0].length/nFeatures;
        double k1, k2;
        
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
        if (extrRotThetas[0].length != 3) {
            throw new IllegalArgumentException("extrRot[0].length must be 3");
        }
        if (extrRotThetas.length != mImages) {
            throw new IllegalArgumentException("extrRot.length must be mImages "
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
                
        /*
        The partial derivatives are implemented following Qu 2018.
        
        The assumptions are no skew, square pixels, f=f_x=f_y.
        
        The factorization is of the reduced camera matrix of the Schur decomposition
        as it is often assumed that (9^3)*mImages < (3^3)*nFeatures, so one
        inverts the matrix HPP (aka V*).
        
        TODO: consider branching if case (9^3)*mImages > (3^3)*nFeatures
        in order to invert matrix HCC (aka U*) instead of matrix HPP. This is
        called the "reduced structure system" in contrast to the "reduced
        camera system" solution pattern below.
        */
        
        if (outInitLambda != null) {
            // this will hold the maximum of the diagonals of hPPI and hCCJ
            outInitLambda[0] = Double.NEGATIVE_INFINITY;
        }
        
        Arrays.fill(outGradC, 0);
        Arrays.fill(outGradP, 0);
        Arrays.fill(outDC, 0);
        Arrays.fill(outDP, 0);
        Arrays.fill(outFSqSum, 0);
        
        double[] bC = outGradC;
        double[] bP = outGradP;
        
        // for each feature i
        double[] bPI = new double[3];
        // for each image j
        double[] bCJ = new double[9];
               
        //HPC is a.k.a. J_C^T*J_P a.k.a. (W^*)^T
        //W_i_j^T are stored here as [3X9] blocks. each W_i_j^T is b_i_j^T*a_i_j 
        BlockMatrixIsometric hPCBlocks /*W^T*/= new BlockMatrixIsometric(MatrixUtil.zeros(3*nFeatures, 9*mImages), 3, 9);
        
        //HPP^-1 is a.k.a. (J_P^T*J_P)^-1 a.k.a. V^-1
        //The diagonals, V_i, are stored here as [3X3] blocks. each V_i is a summation of b_i_j^T*b_i_j over all j images.
        //The inverse is each diagonal block inverted.
        BlockMatrixIsometric hPPIInvBlocks /*VI^-1*/= new BlockMatrixIsometric(MatrixUtil.zeros(3*nFeatures, 3), 3, 3);

        //HCC is a.k.a. J_C^T*J_C a.k.a. U^* (and in Qu thesis is variable B).
        //The diagonals, U_j, are stored here as [9X9] blocks. each U_j is a summation of a_i_j^T*a_i_j over all i features.
        // HCC is set into matrix A as it is calculated.
        // storing hCCJBlocks in mA, while populating mA with only HCC.  later will add the negative rightsize of mA to mA
        //BlockMatrixIsometric hCCJBlocks /*UJ*/= new BlockMatrixIsometric(MatrixUtil.zeros(9*mImages, 9), 9, 9);

        // matrix A is the reduced camera matrix a.k.a. Schur complement. [9m X 9m]; [mXm] block matrix with  blocks [9x9]
        BlockMatrixIsometric mA = new BlockMatrixIsometric(MatrixUtil.zeros(9*mImages, 9*mImages), 9, 9);
        
        // aka V_i; a [3X3] block
        double[][] hPPI = MatrixUtil.zeros(3, 3); 
        //aka (V_i)^-1; a [3X3] block
        double[][] invHPPI = null;
                
        // matrix A, aka reduced camera matrix; [9m X 9m]; [mXm] block matrix with  blocks [9x9]
        //BlockMatrixIsometric mA = new BlockMatrixIsometric(
        //    MatrixUtil.zeros(mImages*9, mImages*9), 9, 9);
        double[][] auxMA = MatrixUtil.zeros(9, 9);
        double[][] auxMA2 = MatrixUtil.zeros(9, 9);
        double[][] auxMA3 = MatrixUtil.zeros(9, 9);
        
        // vector B, on the rhs of eqn; a matrix acting as a vector with m blocks of size [9X1]
 //       double[][] vB = MatrixUtil.zeros(mImages, 9);
        
        double[][] auxHPC = MatrixUtil.zeros(3, 9);
        double[][] auxHPCT = MatrixUtil.zeros(9, 3);

        //[2X9]
        double[][] aIJ = MatrixUtil.zeros(2, 9);
        //[2X3]
        double[][] bIJ = MatrixUtil.zeros(2, 3);

        //aka jP_I_J^T [3X2]
        double[][] bIJT = MatrixUtil.zeros(3, 2);
        //aka jC_I_J^T  [9X2]
        double[][] aIJT = MatrixUtil.zeros(9, 2);
        // aka jP^T*JP; [3X3]
        double[][] bIJsq = MatrixUtil.zeros(3, 3);
        
        // [3X1]
        double[] fIJ = new double[3];
        double[] fIJ2 = new double[2];
        
        //aka bP [3X1]
        double[] bIJTF = new double[3];
        //aka bC [9X1]
        double[] aIJTF = new double[9];
 
        double[] xWI = new double[3];
        double[] xIJ = new double[3];
        double[] xIJC = new double[3];
        double[] xIJHat = new double[3];
        double[] xWCI = new double[3];
        double[] xWCNI = new double[3];
        //double[] xWCNI_2_19 = new double[3];
        double[] xWCNDI = new double[3];
        //double[] xWCNDI_2_20 = new double[3];
        double[][] xIJCs = null;  
        
        // for 0, use image (pixel) reference frame, else for 1 use camera frame.
        final int useCameraFrame = 0;
        
        if (useCameraFrame == 1) {
            xIJCs = transformPixelToCamera(nFeatures, coordsI, intr, kRadials, useR2R4);        
        }
        
        //size is [3 X 3*mImages] with each block being [3X3]
        BlockMatrixIsometric rotMatrices = createRotationMatricesFromEulerAngles(extrRotThetas);
        
        double[] rotAux = new double[3];
        double[][] rotM = MatrixUtil.zeros(3, 3);
                
        AuxiliaryArrays aa = new AuxiliaryArrays();
        double[][] auxIntr = MatrixUtil.zeros(3, 3);
        
        // i for n features, j for m images
        int i, j, k;
        
        //largest runtime complexity in the loop clauses is O(nFeatures * mImages)
        //largest runtime complexity for block multiplication is O(
                     
        for (i = 0; i < nFeatures; ++i) { // this is track p in Engels pseudocode
                        
            //reset hPPI to 0's; [3x3]// a.k.a. V*_i.
            MatrixUtil.fill(hPPI, 0);
            //reset bPI to 0's; //[3X1]
            Arrays.fill(bPI, 0);
            
            //populate xWI; extract the world feature.  size [1X3]
            MatrixUtil.extractColumn(coordsW, i, xWI);
                        
            for (j = 0; j < mImages; ++j) { // this is camera c in Engels pseudocode
                
                //TODO consider how to handle feature not present in image here.
                //   if not in image, then aIJ and bIJ are both 0, so make sure that's handled here
                if (!imageFeaturesMap.get(j).contains(i)) {
                    continue;
                }
                
                // get the rotation matrix rotM [3X3]
                rotMatrices.getBlock(rotM, 0, j);
                
                //transform to camera reference frame. size [1X3]
                Camera.worldToCameraCoordinates(xWI, rotM, extrTrans[j],
                    rotAux, xWCI);
                
                //intr is 3 X 3*nCameras where each block is size 3X3.
                intr.getBlock(auxIntr, j, 0);
          
                // normalize
                for (k = 0; k < 3; ++k) {
                    xWCNI[k] = xWCI[k] / xWCI[2];
                }
                // Qu eqn 2.19, used instead of xWCNI
                /*xWCNI_2_19[0] = -auxIntr[0][0]*xWCNI[0];
                xWCNI_2_19[1] = -auxIntr[1][1]*xWCNI[1];
                xWCNI_2_19[2] = 1;*/
                
                k1 = kRadials[j][0];
                k2 = kRadials[j][1];

                // distort results are in xWCNDI
                if (useCameraFrame == 0) {                            
                    CameraCalibration.applyRadialDistortion(xWCNI, k1, k2, useR2R4, xWCNDI);
                    // Qu eqn 2.20, used instead of xWCNDI.  world feature reprojected to camera and  distorted
                    //CameraCalibration.applyRadialDistortion(xWCNI_2_19, k1, k2, useR2R4, xWCNDI_2_20);
                }
                
                //aIJ is [2X9].  a is partial derivative of measurement vector X w.r.t. camera portion of parameter vector P
                //bIJ is [2X3]
                // populate aIJ and bIJ as output of method:
                aIJBIJ(xWI, xWCI, auxIntr, k1, k2, extrRotThetas[j], rotM, 
                    extrTrans[j], aa, aIJ, bIJ);
                         
                //populate  bIJT; [3X2]  aka jP^T
                MatrixUtil.transpose(bIJ, bIJT);

                //populate aIJT; [9X2] aka jC^T
                MatrixUtil.transpose(aIJ, aIJT);
                                
                // populate xIJ; the observed feature i in image j.  [1X3]
                MatrixUtil.extractColumn(coordsI, nFeatures*j + i, xIJ);
                               
                if (useCameraFrame == 0) {
                    // these are DISTORTED and in IMAGE reference frame
                    
                    // xWCNDI are the world scene points projected into camera frame and distorted
                   
                    // camera to pixel transformation (distortion already applied):
                    //populate xIJHat;the projected distorted feature i into image j reference frame.  [1X3]
                    MatrixUtil.multiplyMatrixByColumnVector(auxIntr, 
                        xWCNDI, xIJHat);
                    
                    //populate xIJHat;the projected feature i into image j reference frame.  [1X3]
                    //MatrixUtil.multiplyMatrixByColumnVector(auxIntr, xWCNI, xIJHat);

if (j==0 && ((i%50)==0)) {
    double[] diff = new double[xIJ.length];
    MatrixUtil.elementwiseSubtract(xIJ, xIJHat, diff);
    System.out.printf("xWCNI=%s\n", FormatArray.toString(xWCNI, "%.3e"));
    System.out.printf("xWCNDI=%s\n", FormatArray.toString(xWCNDI, "%.3e"));
    System.out.printf("*xIJHat=%s\n", FormatArray.toString(xIJHat, "%.3e"));
    System.out.printf("*xIJ=%s\n", FormatArray.toString(xIJ, "%.3e"));
    System.out.printf("*diff=%s\n", FormatArray.toString(diff, "%.7e"));
    System.out.flush();
}
                    // [1X3] - [1X3] = [1X3]
                    MatrixUtil.elementwiseSubtract(xIJ, xIJHat, fIJ);
                    
                } else {  
                    // these are DISTORTION-FREE and in CAMERA reference frame
                                        // populate xIJC; the observed feature i in image j transformed to camera reference frame.  [1X3]
                    MatrixUtil.extractColumn(xIJCs, nFeatures*j + i, xIJC);
                  
                    // xWCNI are the world scene points projected into camera frame
                    
                    // xIJC are the observed points in the camera reference frame with distortion removed
                    
//System.out.printf("(%d,%d)\nxIJC=%s\n", i,j, FormatArray.toString(xIJC, "%.3e"));
//System.out.printf("xWCNI=%s\n\n", FormatArray.toString(xWCNI, "%.3e"));
                
                    // [1X3] - [1X3] = [1X3]
                    MatrixUtil.elementwiseSubtract(xIJC, xWCNI, fIJ);
                }
                
                //sum of squares
                outFSqSum[0] += MatrixUtil.innerProduct(fIJ, fIJ);
                             
 System.out.printf("F(%d,%d)=%.3e\n", i,j,MatrixUtil.innerProduct(fIJ, fIJ));
 System.out.flush();

                //== subtract jP^T*f (aka bP) from bP ==
                 
                //aIJ is [2X9]
                //bIJ is [2X3]
                //dropping the z-axis which is value=0 (from normalized-normalized=1-1)
                System.arraycopy(fIJ, 0, fIJ2, 0, 2);
           
                //bIJTF =  bIJT * fIJ;// [3X2]*[2X1] = [3X1]
                MatrixUtil.multiplyMatrixByColumnVector(bIJT, fIJ2, bIJTF);
                
                /*
                gradP = bP = -JP^T * F  [3n X 1]
                ———————————————-----------------
                -B11T*F11-B12T*F12-B13T*F13
                -B21T*F21-B22T*F22-B23T*F23
                -B31T*F31-B32T*F32-B33T*F33
                -B41T*F41-B42T*F42-B43T*F43

                where each row is [3X2]*[2X1] = [3X1]
                */

                //[3 X 1]
                MatrixUtil.elementwiseSubtract(bPI, bIJTF, bPI);
                
                //System.out.printf("bIJTF=%s\n", FormatArray.toString(bIJTF, "%.3e"));
                //System.out.printf("gradP_%d=%s\n", i, FormatArray.toString(bPI, "%.3e"));
                //System.out.flush();

                System.arraycopy(bPI, 0, bP, i*3, 3);
                
                // populate bIJsq bij^T * bij = [3X2]*[2X3] = [3X3]
                MatrixUtil.multiply(bIJT, bIJ, bIJsq);

                // add jP^T*JP for this image to the hPPi block for feature i.
                //    HPP_i = V_i = Σ_j(BIJT*B1J) 
                // Engels: add jP^T*JP to upper triangular part of hPP aka V.
                // sum bijsq over all images and set into diagonal of hPP_i which is V*_i
                //     elementwise addition of 3X3 blocks:
                MatrixUtil.elementwiseAdd(hPPI, bIJsq, hPPI);
            
                // if camera c is free (presumably, context is fixed or free variables).
                if (true) {
                    
                    //U_j = Σ_i(AIJT*A1J) where U is jCT*JC which is HCC=U
                    //add jCT*JC aka U to upper triangular part of block (j,j) of lhs mA; 
                    //mA[j][j] += aIJT * aIJ; // [9X9]
                    MatrixUtil.multiply(aIJT, aIJ, auxMA);

/*if ((i%25)==0) {                    
System.out.printf("(ft%d,img%d)\naIJ=\n%s\nbIJ=\n%s\n", 
    i, j, FormatArray.toString(aIJ, "%.3e"), FormatArray.toString(bIJ, "%.3e"));
System.out.flush();
}*/
                    mA.getBlock(auxMA2, j, j);
                    MatrixUtil.elementwiseAdd(auxMA, auxMA2, auxMA3);
                    mA.setBlock(auxMA3, j, j);
        
                    //dont use U_j (= mA) yet, nor augment it with damping term
                    //until finished summing over all i
                    
                    // compute block (i,j) of hPC as hPC=jPTJC [3X9]
                    //hPC[i][j] = bIJT * aIJ;
                    MatrixUtil.multiply(bIJT, aIJ, auxHPC);
                    hPCBlocks.setBlock(auxHPC, i, j);
                    
                    //aIJ is [2X9]
                    //bIJ is [2X3]
                
                    // subtract aIJT*f (where bc = -aIJT*f aka -jCT*f) from block row j in vB. [9X1]
                    //aIJTF = aIJT * fIJ.  [9X2]*[2X1]=[1X9]
                    MatrixUtil.multiplyMatrixByColumnVector(aIJT, fIJ2, aIJTF);

                    /* 
                    bC = -JC^T * F =  [9m X 1]
                    ———————————————-
                    -A11T*F11-A21T*F21-A31T*F31-A41T*F41
                    -A12T*F12-A22T*F22-A32T*F32-A42T*F42
                    -A13T*F13-A23T*F23-A33T*F33-A43T*F43

                    where each row is [9X2][2X1] = [9X1]
                    */

                    // fill bCJ with last entry for image J
                    System.arraycopy(bC, j*9, bCJ, 0, 9);
                
                    MatrixUtil.elementwiseSubtract(bCJ, aIJTF, bCJ);

System.out.printf("(ft%d,img%d): bCJ=%s\n     aIJTF=%s\nj*9=%d\n", 
    i, j,  FormatArray.toString(bCJ, "%.3e"), FormatArray.toString(aIJTF, "%.3e"), j*9);
System.out.flush();
        
                    //bC length is [9*mImages]
                    System.arraycopy(bCJ, 0, bC, j*9, 9);                    
                }
            } // end image j loop
            
System.out.printf("bC=%s\n", FormatArray.toString(bC, "%.3e"));

            // now HPP_i (aka V_i) is summed over all j: V_i = Σ_j(BIJT*B1J)

            if (outInitLambda != null) {
                if (
                    maxDiag(hPPI, outInitLambda)
                ) System.out.printf("hPPI) new lambda=%.3e\n", outInitLambda[0])
                ;
            }
            
            // augment hPPI by damping term
            //  (J_P)^T*J_P + lambda*diag((J_P)^T*J_P)
            for (k = 0; k < 3; ++k) {
                hPPI[k][k] *= (1. + lambda);
            }
            //invert hPPI  which is a diagonal block of hPP aka V* // [3X3]
            invHPPI = MatrixUtil.pseudoinverseRankDeficient(hPPI);
            
            hPPIInvBlocks.setBlock(invHPPI, i, 0);
            
        }// end loop i over features
        
        if (outInitLambda != null) {
            for (j = 0; j < mImages; ++j) {
                mA.getBlock(auxMA, j, j);
                if (maxDiag(auxMA, outInitLambda)) {
                    System.out.printf("auxMa %d,%d) new lambda=%.3e\n", i, j, outInitLambda[0]);
                }
            }
        }

        // aIJ^T*aIJ is now in auxMA.  each U_j is a summation of a_i_j^T*a_i_j over all i features.
        // augment U_J with the damping term:  (J_C)^T*J_C + lambda*diag((J_C)^T*J_C)
        for (j = 0; j < mImages; ++j) {
            mA.getBlock(auxMA, j, j);
            for (k = 0; k < auxMA.length; ++k) {
                auxMA[k][k] *= (1. + lambda);
            }
            mA.setBlock(auxMA, j, j);
        }
        
        // HPC is finished and is a dense matrix
        // hPPInv is finished and the diagonals are stored
        // HCC is finished and the diagonals are stored in mA
        // bC is aIJTF and this is stored in outGradC
        // bP is bIJTF and thisis stored in outGradP
        
        System.out.printf("HPC=W^T=\n%s\n\n", FormatArray.toString(hPCBlocks.getA(), "%.3e"));
        System.out.printf("HPPInv=V^-1=\n%s\n\n", FormatArray.toString(hPPIInvBlocks.getA(), "%.3e"));
        System.out.printf("HCC=U=\n%s\n\n", FormatArray.toString(mA.getA(), "%.3e"));
        System.out.printf("gradC=bC=\n%s\n\n", FormatArray.toString(bC, "%.3e"));
        System.out.printf("gradP=bP=\n%s\n\n", FormatArray.toString(bP, "%.3e"));
        System.out.flush();
                
        // blocks: n rows and 1 column
        BlockMatrixIsometric tPRowBlocks = new BlockMatrixIsometric(
            MatrixUtil.zeros(nFeatures, 3), 1, 3);
        double[] tPI = new double[3];
        
        BlockMatrixIsometric tPCBlocks = new BlockMatrixIsometric(MatrixUtil.zeros(9*mImages, 3*nFeatures), 9, 3);
        double[][] tPC = MatrixUtil.zeros(9, 3);
        
        for (i = 0; i < nFeatures; ++i) {
            //calc tP = HPP^-1*bP = V^-1*bP and set into tPBlocks(i,0) += invVI*(Σ_j(-BIJT*F1J))
            //   invHPPI is [3X3]
            hPPIInvBlocks.getBlock(invHPPI, i, 0);
            System.arraycopy(bP, i*3, bPI, 0, 3);
            // [3X3][3X1]=[3X1]
            MatrixUtil.multiplyMatrixByColumnVector(invHPPI, bPI, tPI);
            tPRowBlocks.setRowBlock(tPI, i, 0);
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
            }
        }
        
        double[][] vB = MatrixUtil.zeros(mImages, 9); 
        double[] vBJ = new double[9];
        for (i = 0; i < nFeatures; ++i) {
            tPRowBlocks.getRowBlock(tPI, i, 0);            
            for (j = 0; j < mImages; ++j) { // this is camera c in Engels pseudocode
                //TODO consider how to handle feature not present in image here
                if (!imageFeaturesMap.get(j).contains(i)) {
                    continue;
                }
                hPCBlocks.getBlock(auxHPC, i, j);
                MatrixUtil.transpose(auxHPC, auxHPCT);
                MatrixUtil.multiplyMatrixByColumnVector(auxHPCT, tPI, vBJ);
                for (k = 0; k < vBJ.length; ++k) {
                    vB[j][k] -= vBJ[k];
                }
            }
        }
        // finish vB: vB += bC
        for (j = 0; j < mImages; ++j) { // this is camera c in Engels pseudocode
            System.arraycopy(bC, j*9, bCJ, 0, 9);
            for (k = 0; k < bCJ.length; ++k) {
                vB[j][k] += bCJ[k];
            }
        }
        
        // mA -= tPC*HPC  [9*mImages X 9*mImages]
        double[][] mARight = MatrixUtil.multiply(tPCBlocks.getA(), hPCBlocks.getA());
        mARight = MatrixUtil.elementwiseSubtract(mA.getA(), mARight);
        mA = new BlockMatrixIsometric(mARight, 9, 9);
        
        /* TODO: (optional) Fix gauge by freezing coordinates and thereby reducing 
            the linear system with a few dimensions.
        
           ** Section 9 of Triggs et al. 2000, 
           "Bundle Adjustment – A Modern Synthesis"
               "Section 9 returns to the theoretical issue of gauge freedom
                (datum deficiency), including the theory of inner constraints."
        
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
        
        //Engels step (6) (Linear Solving) 
        //    Cholesky factor the left hand side matrix B and solve for dC. 
        //    Add frozen coordinates back in.
        
        // cholesky decompostion to solve for dC in mA*dC=vB
        // (using the sparsity of upper and lower triangular matrices results in
        //    half the computation time of LU decomposition in comparison)
        
        // mA is square [mImages*9, mImages*9]
        //    but not necessarily symmetric positive definite needed by the
        //    Cholesky decomposition, so need to find the nearest or a nearest
     
System.out.printf("mA=%s\n", FormatArray.toString(mA.getA(), "%.3e"));
System.out.printf("vB=%s\n", FormatArray.toString(vB, "%.3e"));
System.out.flush();

        double eps = 1.e-7;
        // this method attempts to find the nearest symmetric positive *definite* matrix to A:
        double[][] aPSD = MatrixUtil.nearestPositiveSemidefiniteToA(mA.getA(), eps);
       
        /*
System.out.printf("aPSD=%s\n", FormatArray.toString(aPSD, "%.3e"));
        EVD evd = EVD.factorize(new DenseMatrix(aPSD));
        Arrays.sort(evd.getRealEigenvalues());
System.out.printf("EVD(aPSD).eigenvalues=%s\n", FormatArray.toString(evd.getRealEigenvalues(), "%.3e"));
SVD svd = SVD.factorize(new DenseMatrix(aPSD));
System.out.printf("SVD(aPSD).eigenvalues=%s\n", FormatArray.toString(svd.getS(), "%.3e"));
System.out.flush();
        */
        
        //[mImages*9, mImages*9]
        double[][] cholL = LinearEquations.choleskyDecompositionViaLDL(aPSD, eps/10.);
        
 System.out.printf("cholL=%s\n", FormatArray.toString(cholL, "%.3e"));

        /* avoid inverting A by using Cholesky decomposition w/ forward and 
        backward substitution to find x as dC.
            A ﹡ x = b: 
            A = L ﹡ L^*
            L ﹡ y = b ==> y via forward subst
            L^* ﹡ x = y ==> x via backward subst
        */

        // length is vB.length vectorized which is [mImages*9]
        double[] yM = MatrixUtil.forwardSubstitution(cholL, 
            MatrixUtil.reshapeToVector(vB)
        );
        
        // [[mImages*9 X mImages*9] * x = [mImages*9] ==> x is length mImages*9
        // x is dC
        MatrixUtil.backwardSubstitution(MatrixUtil.transpose(cholL), yM, outDC);
          
System.out.printf("yM=%s\n", FormatArray.toString(yM, "%.3e"));
System.out.printf("outDC=%s\n", FormatArray.toString(outDC, "%.3e"));
System.out.flush();

        // calc outDP:  dP = invHPPI * bPI - invHPPI * HPC * dC
        //              dP = tP - tPC^T * dC // [3nX1] - [3nX9m]*[9*mImagesX1] 
        //outDP : [3nX1]
        //tP : [3nX1]
        //tPCBlocks : (9*mImages, 3*nFeatures), 9, 3);
        //outDC: [9*mImagesX1]
        
        double[] dCJ = new double[9];
        double[] tmpI = new double[3]; // holds tPC^T * dC
        double[] tmp = new double[3];
        for (i = 0; i < nFeatures; ++i) {
            //[row i, col 0] = tPI - (Σ_j(tPC[J,I]^T*dCJ))
            Arrays.fill(tmpI, 0);
            for (j = 0; j < mImages; ++j) {
                //[9X3]
                tPCBlocks.getBlock(tPC, j, i);
                // tPC^T [3X9]  in auxHPC
                MatrixUtil.transpose(tPC, auxHPC);
                //dCJ [9X1].
                System.arraycopy(outDC, j*9, dCJ, 0, 9);
                // tmp [3X1]
                MatrixUtil.multiplyMatrixByColumnVector(auxHPC, dCJ, tmp);
                for (k = 0; k < tmp.length; ++k) {
                    tmpI[k] += tmp[k];
                }
            }
            tPRowBlocks.getRowBlock(tPI, i, 0);
            for (k = 0; k < tPI.length; ++k) {
                tPI[k] -= tmpI[k];
                outDP[i*3 + k] = tPI[k];
            }
           
            // Engels: compute updated point (i.e. world coord features)
            // NOTE: for this class, the updates are handled by the invoker of 
            //       this method.  The updated parameters are given to the code so
            //       that when outFSqSum is calculated above, it is using 
            //       the updated parameters and world coords.
            
        } // end loop over feature i
            
System.out.printf("outDP=%s\n", FormatArray.toString(outDP, "%.3e"));
System.out.flush();

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
    static void pdCpIJCIJ(double[] xWCNI, double[][] intr,
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
    static void pdCIJXWIJ(double[] xWCI, double[][] out) {
        
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
    static void pdCpIJYJ(double[] xWCNI, double[][] intr,
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
     * 
     * @param xWI the 3-D coordinates of a world scene feature.
     * @param phi rotation angle vector of length 3 in units of radians
     * @param out output array of size [3X3]
     */
    static void pdXWIJPhiJ(double[] xWI, double[] phi, double[][] out) {
        
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
    static void aIJBIJ(double[] xWI, double[] xWCI, double[][] intr, double k1, double k2, 
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
        
//TODO: may need sign corrections:      
        // 2X2.  Qu 2018 eqn (3.12)
        double[][] dCPdC = aa.a2X2;
        pdCpIJCIJ(xWCNI, intr, k1, k2, dCPdC);
        
        // 2X3.  Qu 2018 eqn (3.16)
        double[][] dCdX = aa.b2X3;
        pdCIJXWIJ(xWCI, dCdX);
        
        // 2X3.  Qu 2018 eqn (3.10)
        double[][] dCPdY = aa.c2X3;
        pdCpIJYJ(xWCNI, intr, k1, k2, dCPdY);
        
        // 3X3.  Qu 2018 eqns (3.28 - 3.33)
        double[][] dXdP = aa.d3X3;
        pdXWIJPhiJ(xWI, rotAngles, dXdP);
       
        //========================================
        
        // [2X3].  Qu 2018 eqn (3.35)
        double[][] dFdT = aa.e2X3;
        MatrixUtil.multiply(dCPdC, dCdX, dFdT);
        
        // [2X3].  Qu 2018 eqn (3.34)
        double[][] dFdPhi = aa.f2X3;
        MatrixUtil.multiply(dFdT, dXdP, dFdPhi);

        // [2X3]. Qu 2018 eqn (3.36)
        double[][] dFdY = dCPdY;
    
        // [2X3].  Qu 2018 eqn (3.37)
        double[][] dFdX = aa.h2X3;
        MatrixUtil.multiply(dFdT, rot, dFdX);
        
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

    private static BlockMatrixIsometric createRotationMatricesFromEulerAngles(
        double[][] extrRotThetas) {
        
        //extrRot is mImages*[1X3]
        
        int mImages = extrRotThetas.length;
        
        double[][] rot = MatrixUtil.zeros(3, 3);
        
        BlockMatrixIsometric m = new BlockMatrixIsometric(MatrixUtil.zeros(3, 3*mImages), 3, 3);
       
        Rotation.AuxiliaryArrays aa = new Rotation.AuxiliaryArrays();
        
        int i;
        for (i = 0; i < mImages; ++i) {
            Rotation.createRotationZYX(extrRotThetas[i], aa, rot);
            m.setBlock(rot, 0, i);
        }
        return m;
    }

    /**
     * gain = (f(p) - f(p + delta p)) / ell(delta p)
             where ell(delta p) is (delta p)^T * (lambda * (delta p)) + J^T * ( b - f))
       gain = (f - fPrev) / ( (delta p)^T * (lambda * (delta p) + J^T * ( b - f)) )
             
     * @param fNew
     * @param fPrev
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
    private static double calculateGainRatio(double fNew, double fPrev, 
        double[] dC, double[] dP, double lambda, 
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

        double[] denom = Arrays.copyOf(dParams, dParams.length);
        for (int i = 0; i < denom.length; ++i) {
            denom[i] *= lambda;
        }

        denom = MatrixUtil.add(denom, gradient);
      
        double d = MatrixUtil.innerProduct(dParams, denom);
            
        if (Math.abs(d) < eps) {
            return Double.NEGATIVE_INFINITY;
        }

        double gain = (fNew - fPrev)/d;

        return gain;
    }
    
    private static boolean isNegligible(double[] c, double eps) {
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
    private static void updateTranslation(double[][] extrTrans, double[] dC) {
        
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
                extrTrans[j][i] += dC[j*9 + (i+3)];
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
    private static void updateIntrinsic(BlockMatrixIsometric intr, 
        double[] dC) {
        
        int nImages = intr.getA().length/3;
        
        double[][] kIntr = MatrixUtil.zeros(3, 3);
        
        //NOTE: follow up on Szeliski text stating that updating the intrinsic parameters is more involved
        
        int j;
        
        // using addition for updates for now
        for (j = 0; j < nImages; ++j) {
            // focus is parameterindex 6 within the 9 for each image in dC
            intr.getBlock(kIntr, j, 0);
            kIntr[0][0] += dC[j*9 + 6];
            kIntr[1][1] += dC[j*9 + 6];
            intr.setBlock(kIntr, j, 0);
        }
    }
        
    private static void updateRadialDistortion(double[][] kRadials, double[] dC) {
        
        int nImages = kRadials.length;
                
        // using addition for updates for now
        for (int i = 0; i < nImages; ++i) {
            kRadials[i][0] += dC[i*9 + 7];
            kRadials[i][1] += dC[i*9 + 8];
        }
    }

    /**
     * update rotation matrix theta vectors with steps in dC.
     * @param extrRotThetas the extrinsic camera parameter rotation euler angles
     * stacked along the 3 columns, that is the size is [nImages X 3].
     * @param dC steps of camera parameters in array of length 9*mImages.
     * dC contains 3 rotation, 3 translation, 3 intrinsic parameters for one
     * image, followed by the same 9 for the next image, etc.  only the
     * rotation elements are used in this method.
     */
    private static void updateRotThetas(double[][] extrRotThetas,  
        double[] dC) {
                                 
        // from Danping Zou lecture notes, Shanghai Jiao Tong University,
        // EE382-Visual localization & Perception, “Lecture 08- Nonlinear least square & RANSAC”
        // http://drone.sjtu.edu.cn/dpzou/teaching/course/lecture07-08-nonlinear_least_square_ransac.pdf
        // parameter perturbations for a Lie group such as rotation are:
        //     R * (I - [delta x]_x) where [delta x]_x is the skew-symetric matrix of delta_x 
        
        // T. Barfoot, et al. 2010, 
        // Pose estimation using linearized rotations and quaternion algebra, 
        // Acta Astronautica (2010), doi:10.1016/j.actaastro.2010.06.049
        // eqn (31) for updating rotation:
        // C(theta) = C(deltaPhi) * previous C
        //     where C is rotation matrix r
        //           calculated as C(theta) = 
        //     and deltaPhi = sTheta * deltaTheta
        
        //    rotation matrix formed from rZ * rY * rX (yaw, pitch, and roll)
        //    which is the same convention used by Wetzstein
                
        double[] deltaTheta = new double[3];
        
        double[][] r = MatrixUtil.zeros(3, 3);
        Rotation.AuxiliaryArrays aa = new Rotation.AuxiliaryArrays();
        double[] thetaExtracted = new double[3];
        
        // ---- update theta ----        
        //extracting theta from the updated rotation would keep the theta
        //    vector consistent with the rotation matrix,
        //    but the value is different that updating theta with delta theta
        //    by addition.
        //    The difference may create a problem with convergence for delta theta.
        
        int j;
                        
        for (j = 0; j < extrRotThetas.length; ++j) {
            
            // copy delta thetas for this image into deltaTheta array
            System.arraycopy(dC, j*9, deltaTheta, 0, 3);
            
            double[] qUpdated = 
                Rotation.applySingularitySafeRotationPerturbationQuaternion(
                extrRotThetas[j], deltaTheta);
        
            // [4X4]
            double[][] qR = Rotation.createRotationMatrixFromQuaternion4(qUpdated);
            qR = MatrixUtil.transpose(qR);

            // rotation is [0:2, 0:2] of qR  

            // update in-out variable r
            for (int i = 0; i < 3; ++i) {
                System.arraycopy(qR[i], 0, r[i], 0, 3);
            }

            // ---- update theta ----        
            //extracting theta from the updated rotation would keep the theta
            //    vector consistent with the rotation matrix,
            //    but the value is a little different that updating theta with delta theta
            //    by addition.
            Rotation.extractThetaFromZYX(r, thetaExtracted);
            
            System.arraycopy(thetaExtracted, 0, extrRotThetas[j], 0, extrRotThetas[j].length);        
        }    
    }
    
    private static void initDeltaPWithQu(double[] deltaP) {
        
        int nFeatures = deltaP.length/3;
                
        /*Qu thesis eqn (3.38)
        
        delta thetas ~ 1e-8
        delta translation ~1e-5
        delta focus ~ 1
        delta kRadial ~ 1e-3
        delta x ~ 1e-8
        */
        int i /*features*/, j /*parameters*/;
        for (i = 0; i < nFeatures; ++i) {
            // i*3 + 0,1,2
            for (j = 0; j < 3; ++j) {
                deltaP[i*3 + j] = 1e-8;
            }
        }
    }
    
    private static void initDeltaCWithQu(double[] deltaC) {
        
        int mImages = deltaC.length/9;
                
        /*Qu thesis eqn (3.38)
        
        delta thetas ~ 1e-8
        delta translation ~1e-5
        delta focus ~ 1
        delta kRadial ~ 1e-3
        delta x ~ 1e-8
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
    private static void updateWorldC(double[][] coordsW, double[] deltaP) {
        int nFeatures = coordsW[0].length;
        
        int i;
        for (i = 0; i < nFeatures; ++i) {
            coordsW[0][i] += deltaP[i*3 + 0];
            coordsW[1][i] += deltaP[i*3 + 1];
            coordsW[2][i] += deltaP[i*3 + 2];
        }
    }
    
    private static boolean maxDiag(double[][] m, double[] outInitLambda) {
        int i;
        boolean updated = false;
        for (i = 0; i < m.length; ++i) {
            if (outInitLambda[0] < m[i][i]) {
                outInitLambda[0] = m[i][i];
                updated =true;
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
    private static double[][] transformPixelToCamera(int nFeatures, double[][] coordsI, 
        BlockMatrixIsometric intr, double[][] kRadial, boolean useR2R4) throws NotConvergedException, IOException {
        
        int nImages = coordsI[0].length/nFeatures;
        
        double[][] c = MatrixUtil.zeros(coordsI.length, coordsI[0].length);
        double[][] x = MatrixUtil.zeros(3, nFeatures);
        double[][] xc;
        double[][] kIntr = MatrixUtil.zeros(3, 3);
        
        int i, j, k;
        for (j = 0; j < nImages; ++j) {
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
    
}
