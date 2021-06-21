package algorithms.imageProcessing.transform;

import algorithms.matrix.MatrixUtil;
import java.util.Arrays;

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
       assumptions used in forming the partial deivatives of the intrinsic camera parameters
       are no skew, focal length along x is the same as focal length along y, square pixels.
       Cholesky decomposition is used to half the runtime compared to LU decomposition.
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
     * @param coordsI
     * @param coordsW
     * @param intr
     * @param extrRot
     * @param extrTrans
     * @param kRadial
     * @param nMaxIter
     * @param useR2R4
     * @param outDP
     * @param outDC
     * @param outGradP
     * @param outGradC
     * @param outFSqSum 
     */
    public static void solveSparse(double[][] coordsI, double[][] coordsW,
        double[][] intr, double[][] extrRot, double[][] extrTrans,
        double[] kRadial, final int nMaxIter, boolean useR2R4,
        double[] outDP, double[] outDC, double[] outGradP, double[] outGradC, double[] outFSqSum) {
            
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
    static void pdCpIJCIJ(double[] xWCNI, Camera.CameraIntrinsicParameters intr,
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
        
        double f1 = intr.getIntrinsic()[0][0];
        
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
    static void pdCpIJYJ(double[] xWCNI, Camera.CameraIntrinsicParameters intr,
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
        
        double f1 = intr.getIntrinsic()[0][0];
        
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
    static void aIJBIJ(double[] xWCI, Camera.CameraIntrinsicParameters intr,
        double k1, double k2, 
        double[] phi,
        double[] trans, 
        AuxiliaryArrays aa,
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
