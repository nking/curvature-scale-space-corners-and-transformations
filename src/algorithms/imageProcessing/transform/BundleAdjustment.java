package algorithms.imageProcessing.transform;

import algorithms.matrix.BlockMatrixIsometric;
import algorithms.matrix.MatrixUtil;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.set.TIntSet;
import java.util.Arrays;
import no.uib.cipr.matrix.DenseCholesky;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.LowerTriangDenseMatrix;
import no.uib.cipr.matrix.NotConvergedException;
import no.uib.cipr.matrix.UpperTriangDenseMatrix;

/**
 given intrinsic and extrinsic camera parameters, coordinates for points
 in a world reference frame, and observations of those points in one or more
 camera, return data structures needed by Levenberg-Marquardt algorithm
 in refinement of the intrinsic and extrinsic camera parameters.
 BundleAdjustment calculates partial derivatives of the parameters
 and calculates the re-projection error to form the parameter update steps,
 the gradient vector, and the evaluation of the objective (sum of squares of
 the re-projection error).
 
 * @author nichole
 */
public class BundleAdjustment {
    
    /**
     * given initial camera calibration and extrinsic parameters, use the 
     * iterative non-linear Levenberg-Marquardt (L-M)
     * algorithm to minimize the re-projection error by finding the best values of
     * intrinsic and extrinsic camera parameters.
     * 
     * The L-M is an iterative non-linear optimization to minimize the 
     * objective.  For bundle adjustment, the objective is the 
     * re-projection error.
 * L-M is guaranteed to converge to eventually find an improvement, 
 * because an update with a sufficiently small magnitude and a negative scalar 
 * product with the gradient is guaranteed to do so.
     * 
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
           
     Tomasi 2007,CPS 296.1 Supplementary Lectur Notes, Duke University
     
     Useful mention of graph partioning in:
        https://cseweb.ucsd.edu/classes/fa04/cse252c/manmohan1.pdf
        recursive partitioning w/ elimination graph and vertex cut.
     Graph partitioning in this project:
        NormalizedCuts.java which uses the Fiedler vector of the Laplacian.
        UnweightedGraphCommunityFinder.java
        
     </pre>
     @param coordsI the features observed in different images (in coordinates 
     * of the image reference frame).  The
     * different images may or may not be from the same camera.  The image
     * to camera relationship is defined in the associative array imageToCamera.
     * The format of coordsI is 3 X (nFeatures*nImages). Each row should
     * have nFeatures of one image, followed by nFeatures of the next image,
       etc.  The first dimension is for the x,y, and z axes.
       Note that if a feature is not present in the image, that should be
       an entry in imageMissingFeatureMap.
     * @param coordsW the features in a world coordinate system.  The format is
     * 3 X nFeatures.  The first dimension is for the x,y, and z axes.
     * @param imageToCamera an associative array relating the image  of
     * coordsI to the camera in intr.  the key is the image
     * and the value is the camera.
     * Note that the number of features (nFeatures) is coordsW[0].length, so the image
     * number in coordsI is j/nFeatures where j is the index of the 2nd dimension,
     * that is coordsI[j], and the camera
     * number in intr is k/3 where k is intr[k];
     * @param imageMissingFeaturesMap an associative array holding the features
     * that are missing from an image.  They key is the image number in coordsI 
     * which is j/nFeatures where j is the index of the 2nd dimension,
     * that is coordsI[j].  The value is a set of feature numbers which are
     * missing from the image.  The feature numbers correspond to the 
     * 2nd dimension indexes in coordsW.
     * @param intr the intrinsic camera parameter matrices stacked along
     * rows in a double array of size 3 X nCameras.
     * @param extrRot the extrinsic camera parameter rotation euler angles
     * stacked along the 3 columns, that is the size is nImages X 3 where
     * nImages is coordsI[0].length/coordsW[0].length.
     * Each row is size [1X3] and that is a set of euler angles for the image number.
     * @param extrTrans the extrinsic camera parameter translation vectors
     * stacked along the 3 columns, that is the size is nImages X 3 where
     * nImages is coordsI[0].length/coordsW[0].length.
     * @param kRadial an array holding radial distortion coefficients k1 and k2.
     * NOTE: current implementation accepts values of 0 for k1 and k2.
     * TODO: consider whether to allow option of leaving out radial distortion
     * by allowing kRadial to be null.
     * @param useR2R4 useR2R4 use radial distortion function from Ma et al. 2004 for model #4 in Table 2,
        f(r) = 1 +k1*r^2 + k2*r^4 if true,
        else use model #3 f(r) = 1 +k1*r + k2*r^2.
     * @param nMaxIter
     * 
     * @throws Exception if there is an error in use of MPSolver during the
     * removal of radial distortion, a generic exception is thrown with the
     * error message from the MPSolve documentation.
     * @throws no.uib.cipr.matrix.NotConvergedException 
     */
    public static void solveSparsely(
        double[][] coordsI, double[][] coordsW,
        TIntIntMap imageToCamera,  TIntObjectMap<TIntSet> imageMissingFeaturesMap,
        double[][] intr, double[][] extrRot, double[][] extrTrans,
        double[] kRadial, final int nMaxIter, boolean useR2R4) 
        throws NotConvergedException, Exception {
        
        // number of features:
        int m = coordsW[0].length;
        
        if (m < 4) {
            throw new IllegalArgumentException("imageC[0].length must be at least 4");
        }
        
        int nImages = coordsI[0].length/m;
       
        if (nImages < 2) {
            throw new IllegalArgumentException("need at least 2 images");
        }
        
        if ((coordsI[0].length % m) > 0) {
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
        
        // TODO: look into methods in MTJ
        
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
                  Some details for the above 3 are in Triggs et al. 1999/2000 Appendix B and near page 23
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
        
        //TODO: consider adding constraints suggestded in Seliski 2010:
        // u_0 and v_0 are close to half the image lengths and widths, respectively.
        // the angle between 2 image axes is close to 90.
        // the focal lengths along both axes are greater than 0.
        
       
                
        throw new UnsupportedOperationException("not yet finished");
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
       The runtime complexity is ~ O(nFeatures * mImages^2)
       assumptions used in forming the partial derivatives of the intrinsic camera parameters
       are no skew, focal length along x is the same as focal length along y, square pixels.
       Cholesky decomposition is used with forward and back substitution
       to avoid inverting the reduced camera matrix
       and to half the runtime compared to LU decomposition
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
     * different images may or may not be from the same camera.  The image
     * to camera relationship is defined in the associative array imageToCamera.
     * The format of coordsI is 3 X (nFeatures*nImages). Each row should
     * have nFeatures of one image, followed by nFeatures of the next image,
       etc.  The first dimension is for the x,y, and z axes.
       Note that if a feature is not present in the image, that should be
       an entry in imageMissingFeatureMap.
     * @param coordsW the features in a world coordinate system.  The format is
     * 3 X nFeatures.  The first dimension is for the x,y, and z axes.
     * @param imageToCamera an associative array relating the image  of
     * coordsI to the camera in intr.  the key is the image
     * and the value is the camera.
     * Note that the number of features (nFeatures) is coordsW[0].length, so the image
     * number in coordsI is j/nFeatures where j is the index of the 2nd dimension,
     * that is coordsI[j], and the camera
     * number in intr is k/3 where k is intr[k];
     * @param imageMissingFeaturesMap an associative array holding the features
     * that are missing from an image.  They key is the image number in coordsI 
     * which is j/nFeatures where j is the index of the 2nd dimension,
     * that is coordsI[j].  The value is a set of feature numbers which are
     * missing from the image.  The feature numbers correspond to the 
     * 2nd dimension indexes in coordsW.
     * @param intr the intrinsic camera parameter matrices stacked along
     * rows in a double array of size 3 X 3*nCameras where each block is
     * size 3X3.
     * @param extrRot the extrinsic camera parameter rotation euler angles
     * stacked along the 3 columns, that is the size is nImages X 3 where
     * nImages is coordsI[0].length/coordsW[0].length.  each array is size
     * 1X3.
     * @param extrTrans the extrinsic camera parameter translation vectors
     * stacked along the 3 columns, that is the size is nImages X 3 where
     * nImages is coordsI[0].length/coordsW[0].length.  each array is size
     * 1X3.
     * @param kRadial an array holding radial distortion coefficients k1 and k2.
     * NOTE: current implementation accepts values of 0 for k1 and k2.
     * TODO: consider whether to allow option of leaving out radial distortion
     * by allowing kRadial to be null.
     * @param useR2R4 useR2R4 use radial distortion function from Ma et al. 2004 for model #4 in Table 2,
        f(r) = 1 +k1*r^2 + k2*r^4 if true,
        else use model #3 f(r) = 1 +k1*r + k2*r^2.
     * @param outDP an output array holding the update values for the point parameters.
     * The length should be 3*nFeatures.
     * @param outDC an output array holding the update values for the camera parameters.
     * The length should be 9*mImages.
     * @param outGradP an output array holding the gradient vector for point parameters
     *  (-J^T*(x-x_hat) as the summation of bij^T*fij).  The length should be 3.
     * This is used by the L-M algorithm to calculate the gain ratio and evaluate stopping criteria.
     * @param outGradC an output array holding the gradient vector for camera parameters
     * (-J^T*(x-x_hat) as the summation of aij^T*fij).
     * The length should be 9.
     * This is used by the L-M algorithm to calculate the gain ratio and evaluate stopping criteria.
     * @param outFSqSum and output array holding the evaluation of the objective,
     * that is the sum of squares of the observed feature - projected feature.
     * It's the re-projection error.
     * The length should be 1.
     */
    static void calculateLMVectorsSparsely(double[][] coordsI, double[][] coordsW,
        TIntIntMap imageToCamera,  TIntObjectMap<TIntSet> imageMissingFeaturesMap,
        BlockMatrixIsometric intr, double[][] extrRot, double[][] extrTrans,
        double[] kRadial, boolean useR2R4,
        double[] outDP, double[] outDC, double[] outGradP, double[] outGradC, double[] outFSqSum) throws NotConvergedException {
            
        int nFeatures = coordsW[0].length;
        int mImages = coordsI[0].length/nFeatures;
        int nCameras = intr.getA()[0].length/3;
        double k1 = kRadial[0];
        double k2 = kRadial[1];
        
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
        if (intr.getA()[0].length != 3) {
            throw new IllegalArgumentException("intr[0].length must be 3");
        }
        if (extrRot[0].length != 3) {
            throw new IllegalArgumentException("extrRot[0].length must be 3");
        }
        if (extrRot.length != mImages) {
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
        if (kRadial.length != 2) {
            throw new IllegalArgumentException("kRadial.length must be 2.");
        }
        
        if (outDP.length != 3*nFeatures) {
            throw new IllegalArgumentException("outDP.length must be 3*nFeatures "
                    + "where nFeatures=coordsW[0].length");
        }
        if (outDC.length != 9*mImages) {
            throw new IllegalArgumentException("outDC.length must be 9*mImages "
                    + "where mImages=coordsI[0].length/coordsW[0].length");
        }
        if (outGradP.length != 3) {
            throw new IllegalArgumentException("outGradP.length must be 3");
        }
        if (outGradC.length != 9) {
            throw new IllegalArgumentException("outGradC.length must be 9");
        }
        if (outFSqSum.length != 1) {
            throw new IllegalArgumentException("outFSqSum.length must be 1");
        } 
        if (imageToCamera == null) {
            throw new IllegalArgumentException("imageToCamera cannot be null");
        }
        if (imageMissingFeaturesMap == null) {
            throw new IllegalArgumentException("imageMissingFeaturesMap cannot be null");
        }
                
        /*
        following the Engels pseudocode for solving for the parameter vector via the
                reduced camera matrix and cholesky factoring with forward and back substitution
                to avoid inverting the reduced camera matrix.
        
        The partial derivatives are implemented following Qu 2018.
        The assumptions are no skew, square pixels, f=f_x=f_y.
        
        The factorization is of the reduced camera matrix of the Schur decomposition
        as it is often assumed that (9^3)*mImages < (3^3)*nFeatures, so one
        inverts the matrix HPP (aka V*).
        TODO: consider branching if case (9^3)*mImages > (3^3)*nFeatures
        in order to invert matrix HCC (aka U*) instead of matrix HPP.
        
        */
        
        // matrix A, aka reduced camera matrix; [9m X 9m]; [mXm] block matrix with  blocks [9x9]
        BlockMatrixIsometric mA = new BlockMatrixIsometric(
            MatrixUtil.zeros(mImages*9, mImages*9), 9, 9);
        double[][] auxMA = MatrixUtil.zeros(9, 9);
        double[][] auxMA2 = MatrixUtil.zeros(9, 9);
        double[][] auxMA3 = MatrixUtil.zeros(9, 9);
        
        // vector B, on the rhs of eqn; a matrix acting as a vector with m blocks of size [9X1]
        double[][] vB = MatrixUtil.zeros(mImages, 9);
        
        // aka V_i; a [3X3] block
        double[][] hPPI = MatrixUtil.zeros(3, 3); 
        //aka (V_i)^-1; a [3X3] block
        double[][] invHPPI;
        
        // aka jPTf; row j of bP; [3X1]
        double[] bPI = new double[3];
        
        //aka jP_I_J^T [3X2]
        double[][] bIJT = MatrixUtil.zeros(3, 2);
        //aka jC_I_J^T  [9X2]
        double[][] aIJT = MatrixUtil.zeros(9, 2);
        // aka jP^T*JP; [3X3]
        double[][] bIJsq = MatrixUtil.zeros(3, 3);
        
        //[2X9]
        double[][] aIJ = MatrixUtil.zeros(2, 9);
        //[2X3]
        double[][] bIJ = MatrixUtil.zeros(2, 3);

        //aka bP [3X1]
        double[] bIJTF = new double[3];
        //aka bC [9X1]
        double[] aIJTF = new double[9];
         
        // [2X1]
        double[] fIJ = new double[2];
        
        //m rows of blocks of size [3X9]
        BlockMatrixIsometric hPCJ = new BlockMatrixIsometric(
            MatrixUtil.zeros(mImages*3, 9), 3, 9);
        double[][] auxHPCJ = MatrixUtil.zeros(3, 9);
        double[][] auxHPCJ2 = MatrixUtil.zeros(3, 9);
        
        // [3X1]
        double[] tP = new double[3]; 
        // [9X3]
        double[][] tPC = MatrixUtil.zeros(9, 3);
        // [9X1]
        double[] tmp = new double[9];
        // [9X3]
        double[][] hPCIJT = MatrixUtil.zeros(9, 3);
        
        // n blocks of [1X3] // i is nFeatures
        double[][] tPs = MatrixUtil.zeros(nFeatures, 3);
        // nXm blocks of [3X9] // i is nFeatures, j is mImages
        BlockMatrixIsometric tPCTs = new BlockMatrixIsometric(
            MatrixUtil.zeros(nFeatures*3, mImages*9), 3, 9);
        double[][] auxTPCTs = MatrixUtil.zeros(3, 9);

        //size is [3 X 3*mImages)with each block being [3X3]
        BlockMatrixIsometric rotMatrices = createRotationMatricesFromEulerAngles(extrRot);
        
        double[] xWI = new double[3];
        double[] xIJ = new double[3];
        double[] xIJHat = new double[3];
        double[] xWCI = new double[3];
        double[] xWCNI = new double[3];
        double[] xWCNDI = new double[3];
        double[] rotAux = new double[3];
        double[][] rotM = MatrixUtil.zeros(3, 3);
                
        AuxiliaryArrays aa = new AuxiliaryArrays();
        double[][] auxIntr = MatrixUtil.zeros(3, 3);
        
        // i for n features, for m images
        int i, j, j2, k, cameraNumber;
                     
        //runtime complexity for this loop is roughly O(nFeatures * mImages^2)
        for (i = 0; i < nFeatures; ++i) { // this is variable p in Engels pseudocode
                        
            //reset hPPI to 0's; [3x3]// a.k.a. V*_i.
            MatrixUtil.fill(hPPI, 0);
            //reset bPI to 0's; //[3X1]
            Arrays.fill(bPI, 0);
            
            //populate xWI; extract the world feature.  size [1X3]
            MatrixUtil.extractColumn(coordsW, i, xWI);
            
            for (j = 0; j < mImages; ++j) {
                
 //TODO consider how to handle feature not present in image here
                
                // get the rotation matrix rotM [3X3]
                rotMatrices.getBlock(rotM, 0, j);
                
                //transform to camera reference frame. size [1X3]
                Camera.worldToCameraCoordinates(xWI, rotM, extrTrans[j],
                    rotAux, xWCI);
                
                // normalize
                for (k = 0; k < 3; ++k) {
                    xWCNI[k] = xWCI[k] / xWCI[2];
                }

                // distort results are in xWCNDI
                CameraCalibration.applyRadialDistortion(xWCNI,
                    k1, k2, useR2R4, xWCNDI);
                
                cameraNumber = imageToCamera.get(j);
                
                //intr is 3 X 3*nCameras where each block is size 3X3.
                intr.getBlock(auxIntr, 0, cameraNumber);
                
                // populate aIJ and bIJ as output of method:
                aIJBIJ(xWCI, auxIntr, k1, k2, extrRot[j], extrTrans[j], aa, aIJ, bIJ);
                
                //populate  bIJT; [3X2]  aka jP^T
                MatrixUtil.transpose(bIJ, bIJT);

                // populate bIJsq bij^T * bij = [3X2]*[2X3] = [3X3]
                MatrixUtil.multiply(bIJT, bIJ, bIJsq);

                // add jP^T*JP to upper triangular part of hPP aka V
                // sum bijsq over all images and set into diagonal of hPP_i which is V*_i
                //     elementwise addition of 3X3 blocks:
                MatrixUtil.elementwiseAdd(hPPI, bIJsq, hPPI);
                
                // populate xIJ; the observed feature i in image j
                MatrixUtil.extractColumn(coordsI, nFeatures*j + i, xIJ);
                
                //populate xIJHat;the projected feature i into image j reference frame.  [1X3]
                MatrixUtil.multiplyMatrixByColumnVector(auxIntr, xWCNDI, xIJHat);

                 // f_i_j [2X1]
                MatrixUtil.elementwiseSubtract(xIJ, xIJHat, fIJ);
                outFSqSum[0] += MatrixUtil.lPSum(fIJ, 2);

                // subtract jP^T*f (aka bP) from bP
                 
                //bIJTF =  bIJT * fIJ;// [3X1]
                MatrixUtil.multiplyMatrixByColumnVector(bIJT, fIJ, bIJTF);
                MatrixUtil.elementwiseSubtract(bPI, bIJTF, bPI);
                
                //populate aIJT; [9X2] aka jC^T
                MatrixUtil.transpose(aIJ, aIJT);

                // if camera c is free means this?
                // if image j has feature i in it ?
                if (!imageMissingFeaturesMap.get(j).contains(i)) {
                    
                    // add jCT*JC aka U to upper triangular part of block (j,j) of lhs mA; // [9X9]
                    //mA[j][j] = aIJT * aIJ;
                    MatrixUtil.multiply(aIJT, aIJ, auxMA);
                    mA.setBlock(auxMA, j, j);

                    //paused here
                    // compute block (i,j) of hPC as hPC=jPTJC [3X9]
                    //    and store until end of image j loop.
                    //hPCJ[j] = bIJT * aIJ;
                    MatrixUtil.multiply(bIJT, aIJ, auxHPCJ);
                    hPCJ.setBlock(auxHPCJ, j, 0);

                    // subtract aIJT*f (where bc = -aIJT*f aka -jCT*f) from block row j in vB. [9X1]
                    //aIJTF = aIJT * fIJ.  [1X9]
                    MatrixUtil.multiplyMatrixByColumnVector(aIJT, fIJ, aIJTF);

                    MatrixUtil.elementwiseSubtract(outGradC, aIJTF, outGradC);
               
                    MatrixUtil.elementwiseSubtract(vB[j], aIJTF, vB[j]);
                             
                }
            } // end image j loop

            //invert hPPI  which is a diagonal block of hPP aka V* // [3X3]
            invHPPI = MatrixUtil.pseudoinverseRankDeficient(hPPI);
            
            // hPPI^-1 * bPI is V^-1*bPI // [3X3][3X1] = [3X1]
            //tP = invHPPI * bPI;
            MatrixUtil.multiplyMatrixByColumnVector(invHPPI, bPI, tP);
            
            //tPs[i] = tP;
            System.arraycopy(tP, 0, tPs[i], 0, tP.length);

            // outer product of feature i
            // for each free camera means?
            for (j = 0; j < mImages; ++j) {
                if (!imageMissingFeaturesMap.get(j).contains(i)) {
                    
                    // subtract hPC^T * tP = hPC^T * (hPP^-1) * bP from part j for image j of rhs vB;
                    //                     = W * V^-1 * bP
                    // [9X1] block

                    // calc hPC^T aka W for feature i, all j images;  [9X3]
                    hPCJ.getBlock(auxHPCJ, j, 0);
                    MatrixUtil.transpose(auxHPCJ, hPCIJT);

                    //tmp = hPC^T * tP = hPC^T * (hPP^-1) * bP = W * V^-1 * bP
                    // [9X3][3X1] = [9X1]
                    MatrixUtil.multiplyMatrixByColumnVector(hPCIJT, tP, tmp);
                     
                    // (at this point, for all j's of current feature i))
                    //   subtract hPC^T * tP from element j of rhs vB.
                    /* e.g.  (W*V^-1)*(bP) =
                                                                  vector element i=1
                            element i=1,j=1                       (all j for this i already calculated)
                              \/                                    \/
                        -> | W11*V1   W21*V2  W31*V3  W41*V4 | * | B11T*F11+B12T*F12+B13T*F13 | <-
                           | W12*V1   W22*V2  W32*V3  W42*V4 |   | B21T*F21+B22T*F22+B23T*F23 |
                           | W13*V1   W23*V2  W33*V3  W43*V4 |   | B31T*F31+B32T*F32+B33T*F33 |
                                                                 | B41T*F41+B42T*F42+B43T*F43 |
                          *The V's are inverses in this matrix
                    */

                    MatrixUtil.elementwiseSubtract(vB[j], tmp, vB[j]);

                    // calc tPC = hPC^T * invHPPI (aka W*V^-1) 
                    // [9X3][3X3] = [9X3]
                    MatrixUtil.multiply(hPCIJT, invHPPI, tPC);

                     // tPC^T = invHPPI^T * hPCIJ = invHPPI * hPCIJ; // [3X9]
                     //set block i,j of tPCTs to MatrixUtil.multiply(invHPPI, hPCJ[j]);
                    MatrixUtil.multiply(invHPPI, auxHPCJ, auxTPCTs);
                    tPCTs.setBlock(auxTPCTs, i, j);
                     
                    for (j2 = j+1; j2 < mImages; ++j2) {
                        // calc tPC * hPCJ[j2] = hPC^T * invHPPI * hPCJ[j2] // [9X3][3X9]=[9X9]
                        //   subtract from block (j, j2) of lhs mA
                        hPCJ.getBlock(auxHPCJ2, j2, 0);
                        MatrixUtil.multiply(tPC, auxHPCJ2, auxMA2);
                        mA.getBlock(auxMA3, j, j2);
                        MatrixUtil.elementwiseSubtract(auxMA3, auxMA2, auxMA3);
                        mA.setBlock(auxMA3, j, j2);
                    }
                } // end if free camera
            } // end image j loop
        } // end i features loop
        
        // (optional) Fix gauge by freezing coordinates and thereby reducing the linear system with a few dimensions.

        // cholesky decompostion to solve for dC in mA*dC=vB
         // (using the sparsity of upper and lower triangular matrices results in
        //    half the computation time of LU decomposition in comparison)
        DenseCholesky chol = no.uib.cipr.matrix.DenseCholesky.factorize(new DenseMatrix(mA.getA()));
        LowerTriangDenseMatrix cholL = chol.getL();
        UpperTriangDenseMatrix cholU = chol.getU();
        double[] yM = MatrixUtil.forwardSubstitution(cholL, 
            MatrixUtil.reshapeToVector(vB)
        );
        
        // dC is [9m X 1]
        double[] dC = MatrixUtil.backwardSubstitution(cholU, yM);

        // [3nX1]
        double[] dP = new double[nFeatures];
        // [9X1]
        double[] dCJ = new double[9];
        
        double[] tmp2 = new double[3];
        for (i = 0; i < nFeatures; ++i) {
            // start with point update for feature i, dP = tP
            //dP[i] = tPs[i]; where tPs is [nFeaturesX3]
            //[1X3]
            System.arraycopy(tPs[i], 0, dP, i*3, 3);
            for (j = 0; j < mImages; ++j) {
                // subtract tPC^T*dCJ where dCJ is for image j (that is dCJ = subvector: dC[j*9:(j+1)*9)
                //[9X1]
                System.arraycopy(dC, j*9, dCJ, 0, 9);
                //tmp2 = tPCTs(i,j)*dCJ;
                //[3X9][9X1]=[3X1]
                tPCTs.getBlock(auxTPCTs, i, j);
                MatrixUtil.multiplyMatrixByColumnVector(auxTPCTs, dCJ, tmp2);
                
                //dP[i] = element wise subtract dP[i] - tmp2;
                for (k = 0; k < 3; ++k) {
                    dP[i*3 + k] -= tmp2[k];
                }
            }
             // compute updated point
        }
             
        throw new UnsupportedOperationException("not yet finished");
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
        
        double x = -xWCNI[0];
        double y = -xWCNI[1];
        
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
     * @param xWCI a world point projected to the camera reference frame.
     * xWCI = column i of coordsW transformed to camera coordinates, but not normalize; 
     * @param phi rotation angle vector of length 3 in units of radians
     * @param out output array of size [3X3]
     */
    static void pdXWIJPhiJ(double[] xWCI, double[] phi, double[][] out) {
        
        if (out.length != 3 || out[0].length != 3) {
            throw new IllegalArgumentException("out size must be 3X3");
        }
        
        double x = xWCI[0];
        double y = xWCI[1];
        double z = xWCI[2];
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
     * @param xWCI a world point projected to the camera reference frame.
     * xWCI = column i of coordsW transformed to camera coordinates, but not normalize; 
     * @param intr
     * @param k1 radial distortion coefficient 1
     * @param k2 radial distortion coefficient 2
     * //@param rot extrinsic camera parameter rotation matrix.
     * @param phi the rotation angles.  formed from 
     * phi = Rotation.extractRotation(rot);
     * @param trans extrinsic camera parameter translation vector.
     * @param aa a group of arrays passed in by invoking code, re-used to avoid
     * constructing more objects.  AuxiliaryArrays aa = AuxiliaryArrays().
     * @param outAIJ output array of size [2X9]
     * @param outBIJ output array of size [2X3]
     */
    static void aIJBIJ(double[] xWCI, double[][] intr, double k1, double k2, 
        double[] phi, double[] trans, AuxiliaryArrays aa,
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
        
        // 2X2
        double[][] dCPdC = aa.a2X2;
        pdCpIJCIJ(xWCNI, intr, k1, k2, dCPdC);
        
        // 2X3
        double[][] dCdX = aa.b2X3;
        pdCIJXWIJ(xWCI, dCdX);
        
        // 2X3
        double[][] dCPdY = aa.c2X3;
        pdCpIJYJ(xWCNI, intr, k1, k2, dCPdY);
        
        // 3X3
        double[][] dXdP = aa.d3X3;
        pdXWIJPhiJ(xWCI, phi, dXdP);
       
        //========================================
        
        // [2X3]
        double[][] dFdT = aa.e2X3;
        MatrixUtil.multiply(dCPdC, dCdX, dFdT);
        
        // [2X3]
        double[][] dFdPhi = aa.f2X3;
        MatrixUtil.multiply(dFdT, dXdP, dFdPhi);
        
        // [2X3]
        double[][] dFdY = dCPdY;
        
        // 3X3
        double[][] rot = aa.g3X3;
        Rotation.calculateRotationZYX(phi, aa.aa, rot);
        
        // [2X3]
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
        double[][] extrRot) {
        
        //extrRot is mImages*[1X3]
        
        int mImages = extrRot.length;
        
        double[][] rot = MatrixUtil.zeros(3, 3);
        
        BlockMatrixIsometric m = new BlockMatrixIsometric(MatrixUtil.zeros(3, 3*mImages), 3, 3);
       
        Rotation.AuxiliaryArrays aa = new Rotation.AuxiliaryArrays();
        
        int i;
        for (i = 0; i < mImages; ++i) {
            Rotation.calculateRotationZYX(extrRot[i], aa, rot);
            m.setBlock(rot, 0, i);
        }
        return m;
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
