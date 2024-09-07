package algorithms.imageProcessing.transform;

import algorithms.matrix.MatrixUtil;

import java.util.Arrays;

import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.NotConvergedException;
import no.uib.cipr.matrix.SVD;

/**
 * a utility class holding rotations associated with euler angles,
 * quaternions, and angle-axis representations (using rodrigues formula for the later).
 * 
 * a rotation matrix describes a transformation in euclidean space.  it is a
 * member of  special orthogonal group SO(n) where n is usually 3, but can
 * be 4 for quaternion rotation matrix.  It has the properties R^T = R^−1 and det(R) = +-1.
 <pre>
   det(R)=1 is a proper rotation matrix.  rotation angles are counterclockwise.
       it's a special orthogonal matrix and provides the
       defining matrix representation of the group of proper n-dimensional rotations, denoted
       by SO(n). http://scipp.ucsc.edu/~haber/ph251/rotreflect_17.pdf
   det(R)=-1 is an improper rotation matrix representing rotations that
       require mirrors.
       The most general improper rotation matrix is a product of a proper rotation by an
       angle θ about some axis nˆ and a mirror reflection through a plane that passes through
       the origin and is perpendicular to nˆ.  NOTE: nˆ is determined by
       the right hand rule.
 </pre>
 * This Rotation.java class uses "passive" Euler rotations of vectors counterclockwise in a right-handed
 * coordinate system (y counterclockwise from x) by pre-multiplication 
 * (R on the left) unless otherwise specified.
 * If any one of these is changed (such as rotating axes
 * instead of vectors, a "active" transformation), then the inverse of the
 * example matrix should be used, which coincides with its transpose.
 * e.g. (A*B*C)^-1 = (C^-1) * (B^-1) * (A^-1).

 see wikipedia article "Active and Passive Transformations"

 active rotations (a.k.a. alibi transformations):
     the reference frame stays fixed and the object of interest moves. The rotation is of the object being described.
     These are clockwise rotations (left hand system) of the object of interest about the fixed reference system origin .
     active rotations are the historical system used.

 passive rotations (a.k.a. alias transformations):
     the reference frame (coordinate system) is rotated and the object of interest stays fixed.
     These are counterclockwise rotations (right hand system) of the reference frame about its own
     origin while the object of interest says fixed.
     An example use is that from the perspective of being inside a plane - the inertial reference frame appears to move
     with the opposite rotation.  This system is used when we are describing the motion of the object we are in and
     controlling.

 an active rotation R(theta) is equivalent to the passive rotation R(-theta) which is R^T.

 active transformations are often used for multiple maneuvers of a body.

 * <b><ul>RIGHT HAND SYSTEM</ul></b>
 * The equations in this section use a right hand system (== passive transformations, CCW rotations of the object while
 * reference frame is fixed).
 * The right hand system is consistent with methods in physics, engineering,
 * and computer science in general.  The <b>Hamilton quaternion</b> is consistent
 * with the right hand system and has format [scalar vector].
 * The NASA 1977 publication and Szeliski 2010
 *    define a quaternion as (qw, qx, qy, qz) where qw is a scalar and
 * [qx, qy, qz] is a vector.   That is the Hamilton Quaternion format.
 * 
 * In contrast, Barfoot, Forbes, & Furgale 2010, "Pose estimation using linearized
 * rotations and quaternion algebra", Acta Astronautica (2010), doi:10.1016/j.actaastro.2010.06.049.)
 * define a quaternion as (qx, qy, qz, qw), scalar last.
 * One might need to transform some properties to quaternions and then modify 
 * the placement of the scalar term in order
 * to compare results with other "Hamilton" "Right hand" system results with scalar first format.
 * The Barfoot equations use active left-hand rule system and a "scalar last" format.
 *         
 <pre>
 R(q) = | (w^2 + x^2 - y^2 - z^2)  2(xy - wz)               2(xz + wy)               |
        | 2(xy + wz)               w^2 - x^2 + y^2 - z^2)   2(yz - wx)               |
        | 2(xz - wy)               2(yz + wx)               (w^2 - x^2 - y^2  + z^2) |
 rewritten 
 R(q) = | 1-2*(y^2 + z^2)   2(xy - zw)       2(xz + yw)      |
        | 2(xy + zw)        1-2(x^2 + z^2)   2(yz - xw)      |
        | 2(xz - yw)        2(yz + xw)       1-2(x^2 + y^2)  |
 </pre>
 Quaternions can be derived from the axis/angle representation through the 
 formula q = (v, w) = (sin(theta)*n_hat,cos(theta)),
 where n_hat and theta are the rotation axis and angle.

 passive, CCW rotations (right hand system) of the reference frame about its own
 origin while the object of interest says fixed.
     rotation about z-axis (yaw):    about the y-axis (pitch):     about x-axis (roll):
                | cos φ   -sin φ    0 |          |  cos ψ    0  sin ψ |         |    1       0       0 |
                | sin φ    cos φ    0 |          |      0    1      0 |         |    0   cos θ  -sin θ |
                |     0        0    1 |          | -sin ψ    0  cos ψ |         |    0   sin θ   cos θ |
     Note that in the java methods below, the documentation should specify whether active or passive are used.

 applications of quaternions are intrinsic or extrinsic.

 intrinsic:
    apply transformation to the axis of the rotated coordinate reference frame (a.k.a. body frame).
    expresses F_w relative to F_b (F_b->F_w).
    vec_c = R_{b->c} * (R_{a->b} * vec_a) = R_{a->c} * vec_a
    The Barfoot et al. paper and book uses intrinsic.

 extrinsic:
    apply transformation to the axis in the World Coordinate reference frame.
    expresses F_b relative to F_w (F_w->F_b).
    Direct Cosine Matrices (DCM) used in extrinsic operations...

 The skew symmetric matrix (implemented in MatrixUtil.skewSymmetric(), is:
 [v]_x = |  0  -z   y |
         |  z   0  -x |
         | -y   x   0 |

 Its opposite is the anti-symmetric matrix and is equal to -1 * skewSymmetric:
 [[v]] = |  0   z  -y |
         | -z   0   x |
         |  y  -x   0 |

 From Shuster we have properties of the skew symmetric matrix (as -1 * anti-sym matrix [[v]]).
 <pre>
 let [u]_x be the skew symmetric of u and [v]_x be the skew-symmetric of v.
     where u and v are column vectors.

     -[u]_x * v = [v]_x * u

     -[u]_x * u = 0

     [u]_x = [ -[u]_x ]^T

     -[u]_x * -[v]_x = (u dot v) * I + v * u^T where '*' is matrix multiplication as usual.
     [u]_x * [v]_x = (u dot v) * I + v * u^T

     [u]_x * [v]_x - [v]_x * [u]_x = v * u^T - u*v^T = [u cross v]]_x

     u * v^T * -[w]_x + -[w]_x * v * u^T = -[u cross (v cross w)]_x
           where the later, triple cross product is used in making the Fundamental Matrix
           of photogrammetry, for example.
            implemented in MatrixUtil.tripleProduct() using a shorter form from Boas.
 </pre>
 
 * Eigen, ROS, and Google Ceres use Hamilton convention (active, LH, CCW rotations of object in fixed reference frame).
*  Also  Wolfram Mathematica, Matlab’s aerospace(!) and robotics toolbox, 
*  Boost, GNU Octave, NASA’s SPICE.

Note that Shuster 1993 use the active, LH, CW rotations of the object while reference frame is fixed.
        about z-axis (yaw):           about the y-axis (pitch):    about x-axis (roll): 
            | cos φ    sin φ    0 |    |  cos ψ    0 -sin ψ |      |    1       0       0 |  
            |-sin φ    cos φ    0 |    |      0    1      0 |      |    0   cos θ   sin θ |  
            |     0        0    1 |    |  sin ψ    0  cos ψ |      |    0  -sin θ   cos θ |
 Note that in the java methods below, the documentation should specify whether active or passive are used.


 <pre>
 Some references in the right hand systems:
 
 "Euler Angles and Quaternions and Transformations", NASA Mission Planning
  and Analysis Division, 1977, Shuttle Program.

  Szeleski 2010
  
  Shuster 1993, "A Survey of Attitude Representations"
    Journal of Astronautical Sciences, Vol 41, No. 4, Oct-Dec 1993, pp 439-517
    http://www.ladispe.polito.it/corsi/Meccatronica/02JHCOR/2011-12/Slides/Shuster_Pub_1993h_J_Repsurv_scan.pdf
  
  "Representing Attitude: Euler Angles, Unit Quaternions, and Rotation Vectors"
   James Diebel, 2006
 </pre>
* 
 * <pre><b><ul>LEFT HAND SYSTEM</ul></b></pre>
 * A left hand system, in contrast, is used by aerospace engineering and
 * aerospace medicine to describe the z-axis as pointing downwards (in the 
 * direction of gravity w.r.t. a geo-centric system.)
 * 
 * Chapter 4, "Human Response to Acceleration" by Banks, Brinkly, Allnut, and Harding
    in "Fundamentals of Aerospace Medicine"
    Excerpt:
    "For example, the Advisory Group for Aerospace
    Research and Development (AGARD) standard for human
    acceleration differs from the AGARD standard for aircraft
    design (in which the z-acceleration axis is reversed and
    positive downward).
        ...consistent with the AGARD standard (1), the
    Table of Equivalents for Acceleration Terminology (2), the
    Aviation Space and Environmental Medicine Standard (3),
    and the majority of the Aerospace Medicine literature, the
    positive direction of each of these axes is here described by
    ‘‘the left-hand rule.’’   That is, the x-axis dimension is an arrow
     with the positive direction forward, the y-axis dimension has
    the positive direction rightward, and"
  
  <pre>

 Rotations about the same axis are additive.
 
* Regarding rotation axes and perspectives:
*     https://en.wikipedia.org/wiki/Rotation_matrix#Ambiguities
*     The coordinates of a point P may change due to either a rotation of the 
*     coordinate system CS (alias), or a rotation of the point P (alibi). 
*     In the latter case, the rotation of P also produces a rotation of the 
*     vector v representing P. In other words, either P and v are fixed while 
*     CS rotates (alias), or CS is fixed while P and v rotate (alibi). Any 
*     given rotation can be legitimately described both ways, as vectors and 
*     coordinate systems actually rotate with respect to each other, about the 
*     same axis but in opposite directions. Throughout this article, we chose the 
*     alibi approach to describe rotations. For instance, 
*        R(theta)= | cos(theta)  -sin(theta) |
*                  | sin(theta)   cos(theta) |
*     represents a counterclockwise rotation of a vector v by an angle θ, or a 
*     rotation of the coordinate system (CS) by the same angle but in the opposite direction 
*     (i.e. clockwise). Alibi and alias transformations are also known as active 
*     and passive transformations, respectively.
* 
*     Pre-multiplication or post-multiplication
         The same point P can be represented either by a column vector v or a 
         row vector w. Rotation matrices can either pre-multiply column vectors 
         (Rv), or post-multiply row vectors (wR). However, Rv produces a 
         rotation in the opposite direction with respect to wR. 
         Throughout this article, rotations produced on column vectors are 
         described by means of a pre-multiplication. To obtain exactly the same 
         rotation (i.e. the same final coordinates of point P), the equivalent 
         row vector must be post-multiplied by the transpose of R (i.e. wR^T).
 
 skew symmetric of vector v is [v]_x:
       |    0   -v[2]  v[1]  |
       |  v[2]    0    -v[0] |
       | -v[1]  v[0]    0    |
    
    e^(i*theta) = cos(theta) + i*sin(theta)
    
    exponential of a matrix is expanded via power series:  e^A = I + A + (higher order terms AA/2! + AAA/3!...)
    
    for A = | 0  -1 |
            | 1   0 |
    
    e^(A*t) = |  cos(t)  -sin(t)  |  through diagonalizing A
              |  sin(t)   cos(t)  |
    
    for A being skew symmetric: |  0  -z  y |
                                |  z   0 -x |
                                | -y   x  0 |
    e^(A*theta) = I + A*sin(theta) +A^2*(1 - cos(theta))
    
    
    R_axis = exp(-[e_axis]_x*theta_axis)
    
    R_x = exp(-[1,0,0]_x * theta_x) = exp( 0  0  0 )*theta_x )
                                           0  0  1
                                           0 -1  0
      
   terms used describing axes of rotation, attitude or orientation:
        heading, attitude, bank,
        pitch, yaw, roll,
        pan, tilt, roll
        
   A rotation can be represented by a rotation axis nˆ and an angle ω,
   or equivalently by a 3D vector ω = θ*nˆ
   to project vector v onto the axis nˆ:
      v_parallel = nˆ*(nˆ· v) = (nˆ*nˆ^T)*v and this is not affected by rotation
      v_perpendicular = v - v_parallel
                      = (I - nˆ*nˆ^T)*v
   can rotate v by 90 degrees using the cross product:
      nˆ cross v = [nˆ]_x * v
      where [nˆ]_x is the skew-symmetric matrix of nˆ

   can rotate v by 180 degrees:
      nˆ cross (nˆ cross v) = ([nˆ]_x)^2 * v = -v_perpendicular
   which shows that can also write
      v_parallel = (I + ([nˆ]_x)^2 * v)

 examples transforming sequences:
     passive A(a)*B(b)*C(c)
     active A(-a)*B(-b)*C(-c)
                   = A(a)^T * B(b)^T * C(c)^T
                   = (C(c) * B(b) * C(c))^T

     passive intrinsic XYZ = extrinsic transposed(active)
                           = extrinsic transposed(XYZ(-1*angles))
     extrinsic passive XYZ = intrinsic active XYZ
                           = intrinsic transpose(XYZ(-1*angles))

     example:
        from scipy.spatial.transform import Rotation as R
        R.from_euler('XYZ', [-0.1, -0.2, -0.3]).as_matrix()
        array([[ 0.93629336,  0.28962948, -0.19866933],
               [-0.27509585,  0.95642509,  0.0978434 ],
               [ 0.21835066, -0.03695701,  0.97517033]])
        R.from_euler('xyz', [0.1, 0.2, 0.3]).inv().as_matrix()
        array([[ 0.93629336,  0.28962948, -0.19866933],
               [-0.27509585,  0.95642509,  0.0978434 ],
               [ 0.21835066, -0.03695701,  0.97517033]])
 </pre>

 * TODO: add a method to extract the quaternion rotation from a 4X4 rotation matrix.
 * see Barfoot et al.
 * see http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/
 *
 * TODO: add slerp
 *
 * <pre>
  also note that for optimization of rotation :
      prefer to update the rotation vectors and/or the rotation matrices.
      do not update the euler angles because eqn (28) of Barfoot et al. shows that an inverse term
      is not defined at singularities.
      do not extract euler angles or vector from a rotation matrix unless unavoidable because the
      extraction is ambiguous.
      the create rotation matrix from rotation vector or from euler angles is fine.
 * </pre>
 *
 * @author nichole
 */
public class Rotation {

    /**
     * types of Euler rotation sequences implemented in this class for the Barfoot quaternion
     * equations.
     * Note:
     * passive A(a)*B(b)*C(c)
     * active A(-a)*B(-b)*C(-c)
     *         = A(a)^T * B(b)^T * C(c)^T
     *         = (C(c) * B(b) * C(c))^T
     */
    public static enum EulerSequence {
        XYZ_ACTIVE /*1-2-3*/, ZYX_ACTIVE /*3-2-1*/
    };

    public abstract static class RotationPerturbation {
        /**
         * dPhi is the rotation vector defined near eqn (26) of Barfoot et al.
         * Barfoot, Forbes, & Furgale 2010, "Pose estimation using linearized
         * rotations and quaternion algebra", Acta Astronautica (2010), doi:10.1016/j.actaastro.2010.06.049.)
         */
        double[] dPhi;

        final EulerSequence seq;

        public RotationPerturbation(double[] dPhi, EulerSequence seq) {
            this.dPhi = dPhi;
            this.seq = seq;
        }
        public abstract RotationPerturbation copy();
    }
    public static class RotationPerturbationMatrix extends RotationPerturbation{
        /**
         * the current updated rotation
         */
        double[][] rotation;
        public RotationPerturbationMatrix(double[] dPhi, EulerSequence seq) {
            super(dPhi, seq);
        }
        public RotationPerturbationMatrix(double[] dPhi, EulerSequence seq, double[][] rotation) {
            super(dPhi, seq);
            this.rotation = rotation;
        }
        public RotationPerturbationMatrix copy() {
            return new RotationPerturbationMatrix(Arrays.copyOf(this.dPhi, dPhi.length), this.seq,
                    MatrixUtil.copy(rotation));
        }
    }
    public static class RotationPerturbationQuaternion extends RotationPerturbation{
        /**
         * the current updated quaternion
         */
        double[] quaternion;
        public RotationPerturbationQuaternion(double[] dPhi, EulerSequence seq) {
            super(dPhi, seq);
        }
        public RotationPerturbationQuaternion(double[] dPhi, EulerSequence seq, double[] quaternion) {
            super(dPhi, seq);
            this.quaternion = quaternion;
        }
        public RotationPerturbationQuaternion copy() {
            return new RotationPerturbationQuaternion(Arrays.copyOf(this.dPhi, dPhi.length), this.seq,
                    Arrays.copyOf(quaternion, quaternion.length));
        }
    }

    /**
     * calculate R(angle_z, angle_y, angle_x) = R_x(angle_x)*R_y(angle_y)*R_z(angle_z).
     * The internal composition follows passive (right hand), CCW rotations
     * with order of operations ((R(X) * R(Y)) * R(Z)), that is left to right,
     * intrinsic composition.
     *
     * If the context is aeronautical, you may want to use createRotationZYX instead.
     *
     * <pre>
     This method uses a passive (right-hand) rotation system with CCW (=CC) rotations.
       cc rotation about z-axis (yaw):   cc about the y-axis (pitch):    cc about x-axis (roll):    
            | cos φ   -sin φ    0 |          |  cos ψ    0  sin ψ |         |    1       0       0 |  
            | sin φ    cos φ    0 |          |      0    1      0 |         |    0   cos θ  -sin θ |  
            |     0        0    1 |          | -sin ψ    0  cos ψ |         |    0   sin θ   cos θ |  
      </pre>
     * @param angleX angle of rotation about x-axis (roll) in units of radians.
     * @param angleY angle of rotation about y-axis (pitch) in units of radians.
     * @param angleZ angle of rotation about z-axis (yaw) in units of radians.
     * @return 
     */
    public static double[][] createRotationXYZ(double angleX, double angleY, double angleZ) {
        double[][] rX = Rotation.createRotationRoll(angleX);
        double[][] rY = Rotation.createRotationPitch(angleY);
        double[][] rZ = Rotation.createRotationYaw(angleZ);

        return MatrixUtil.multiply(MatrixUtil.multiply(rX, rY), rZ);
    }

    public static double[][] createRotationXYZ(double angleX, double angleY, double angleZ, boolean passive) {
        if (passive) {
            return createRotationXYZ(angleX, angleY, angleZ);
        }
        return createRotationXYZ(-angleX, -angleY, -angleZ);
    }

    public static double[][] createRotationZYX(double angleX, double angleY, double angleZ, boolean passive) {
        if (passive) {
            return createRotationZYX(angleX, angleY, angleZ);
        }
        return createRotationZYX(-angleX, -angleY, -angleZ);
    }

    /**
     * create matrix for rotation about the X-axis, a.k.a. roll in passive system (right-hand system, CCW).
       <pre>
       about x-axis (roll):       
       |    1       0       0 |  
       |    0   cos θ  -sin θ |  
       |    0   sin θ   cos θ | 
       </pre>
     * @param angle angle of rotation about x-axis (roll) in units of radians.
     * @return 
     */
    public static double[][] createRotationRoll(double angle) {
        
        double[][] rot = MatrixUtil.zeros(3, 3);
        
        double c = Math.cos(angle);
        double s = Math.sin(angle);
        
        rot[0][0] = 1;
        rot[1][1] = c;
        rot[1][2] = -s;
        rot[2][1] = s;
        rot[2][2] = c;
        
        return rot;
    }

    public static double[][] createRotationRoll(double angle, boolean passive) {
        if (passive) {
            return createRotationRoll(angle);
        }
        return createRotationRoll(-angle);
    }
    
    /**
     * create matrix for rotation about the X-axis, a.k.a. roll.
     * This method uses a passive (right-hand) rotation system with CCW (=CC) rotations.
       <pre>
       about x-axis (roll):       
       |    1       0       0 |  
       |    0   cos θ  -sin θ |  
       |    0   sin θ   cos θ | 
       </pre>
     * @param angle angle of rotation about x-axis (roll) in units of radians.
     * @param out holds values for rotation matrix for roll. 
     */
    public static void createRotationRoll(double angle, double[][] out) {

        if (out.length != 3 || out[0].length != 3) {
            throw new IllegalArgumentException("out must be 3x3");
        }
        
        
        double c = Math.cos(angle);
        double s = Math.sin(angle);
        
        out[0][0] = 1;
        out[0][1] = 0;
        out[0][2] = 0;
        out[1][0] = 0;
        out[1][1] = c;
        out[1][2] = -s;
        out[2][0] = 0;
        out[2][1] = s;
        out[2][2] = c;        
    }        
    
    /**
    create matrix for rotation about the Y-axis, a.k.a. pitch in passive system (right-hand, CCW).
      <pre>
      about the y-axis (pitch):
      |  cos ψ    0  sin ψ |
      |      0    1      0 |
      | -sin ψ    0  cos ψ |
      </pre>
     * @param angle angle of rotation about y-axis (pitch) in units of radians.  Euler notation
     *              uses the right-hand rule for rotations.
     * @return 
    */
    public static double[][] createRotationPitch(double angle) {
        
        double[][] rot = MatrixUtil.zeros(3, 3);
        
        double c = Math.cos(angle);
        double s = Math.sin(angle);
        
        rot[0][0] = c;
        rot[0][2] = s;
        rot[1][1] = 1;
        rot[2][0] = -s;
        rot[2][2] = c;
        
        return rot;
    }

    public static double[][] createRotationPitch(double angle, boolean passive) {
        if (passive) {
            return createRotationPitch(angle);
        }
        return createRotationPitch(-angle);
    }
    
    /**
    create matrix for rotation about the Y-axis, a.k.a. pitch.
     This method uses a passive (right-hand) rotation system with CCW (=CC) rotations.
      <pre>
      about the y-axis (pitch):
      |  cos ψ    0  sin ψ |
      |      0    1      0 |
      | -sin ψ    0  cos ψ |
      </pre>
     * @param angle angle of rotation about y-axis (pitch) in units of radians.
     * @param out holds values for rotation matrix for pitch 
    */
    public static void createRotationPitch(double angle,
                                           double[][] out) {

        if (out.length != 3 || out[0].length != 3) {
            throw new IllegalArgumentException("out must be 3x3");
        }
        
        double c = Math.cos(angle);
        double s = Math.sin(angle);
        
        out[0][0] = c;
        out[0][1] = 0;
        out[0][2] = s;
        out[1][0] = 0;
        out[1][1] = 1;
        out[1][2] = 0;
        out[2][0] = -s;
        out[2][1] = 0;
        out[2][2] = c;        
    }
    
    /**
    create matrix for rotation about the Z-axis, a.k.a. yaw in passive system (right-hand system, CCW).
      <pre>
        about z-axis (yaw):          
            | cos φ   -sin φ    0 | 
            | sin φ    cos φ    0 | 
            |     0        0    1 | 
      </pre>
     * @param angle angle of rotation about z-axis (yaw) in units of radians.
     * @return 
    */
    public static double[][] createRotationYaw(double angle) {
        
        double[][] rot = MatrixUtil.zeros(3, 3);
        
        double c = Math.cos(angle);
        double s = Math.sin(angle);
        
        rot[0][0] = c;
        rot[0][1] = -s;
        rot[1][0] = s;
        rot[1][1] = c;
        rot[2][2] = 1;
        
        return rot;
    }

    public static double[][] createRotationYaw(double angle, boolean passive) {
        if (passive) {
            return createRotationYaw(angle);
        }
        return createRotationYaw(-angle);
    }
    
    /**
    create matrix for rotation about the Z-axis, a.k.a. yaw.
     This method uses a passive (right-hand) rotation system with CCW (=CC) rotations.
      <pre>
        about z-axis (yaw):          
            | cos φ   -sin φ    0 | 
            | sin φ    cos φ    0 | 
            |     0        0    1 | 
      </pre>
     * @param angle angle of rotation about z-axis (yaw) in units of radians.
     * @param out results for a rotation matrix for yaw.  size given to method 
     * must be 3X3.
    */
    public static void createRotationYaw(double angle, double[][] out) {

        if (out.length != 3 || out[0].length != 3) {
            throw new IllegalArgumentException("out must be 3x3");
        }        
        
        double c = Math.cos(angle);
        double s = Math.sin(angle);
        
        out[0][0] = c;
        out[0][1] = -s;
        out[0][2] = 0;
        out[1][0] = s;
        out[1][1] = c;
        out[1][2] = 0;
        out[2][0] = 0;
        out[2][1] = 0;
        out[2][2] = 1;
        
    }

    /**
     * calculate the rotation needed to transform direction v1 to direction v2 using
     * a quaternion.
     * @param v1
     * @param v2
     * @return
     */
    public static double[][] rotationBetweenTwoDirections0(double[] v1, double[] v2) {
        double[] n1 = MatrixUtil.normalizeL2(v1);
        double[] n2 = MatrixUtil.normalizeL2(v2);
        double angle = -Math.acos(n1[0]*n2[0] + n1[1]*n2[1] + n1[2]*n2[2]);
        double[] axis = MatrixUtil.crossProduct(n1, n2);
        //TODO: check the method createQuaternionUnitLengthBarfootFromEulerXYZ with results for quatB here
        double[] quatB = Rotation.createQuaternionHamiltonFromAngleAxis(angle, axis);
        double[][] r = Rotation.createRotation4FromQuaternion(quatB);
        r = MatrixUtil.copySubMatrix(r, 0, 2, 0, 2);
        return r;
    }

    /**
     * calculate the rotation needed to transform direction v1 to direction v2 using
     * Rodrigues formula.
     * @param v1
     * @param v2
     * @return
     */
    public static double[][] rotationBetweenTwoDirections1(double[] v1, double[] v2) {
        double[] n1 = MatrixUtil.normalizeL2(v1);
        double[] n2 = MatrixUtil.normalizeL2(v2);
        double angle = -Math.acos(n1[0]*n2[0] + n1[1]*n2[1] + n1[2]*n2[2]);
        double[] axis = MatrixUtil.crossProduct(n1, n2);
        double[][] r = Rotation.createRotationRodriguesFormula(axis, false);
        return r;
    }

    public static double createAngleAxisFromRotationVector(double[] rV, double[] outRAxis) {
        double theta = MatrixUtil.lPSum(rV, 2);
        for (int i = 0; i < rV.length; ++i) {
            outRAxis[i] = rV[i]/theta;
        }
        return theta;
    }

    public static double[] createRotationVectorFromAngleAxis(double[] axis, double angle) {
        double[] rVec = Arrays.copyOf(axis, axis.length);
        MatrixUtil.multiply(rVec, angle);
        return rVec;
    }

    /**
     * convert the euler angles to angle axis representation.
     * NOTE that the resulting vector is the angle axis times the angle, so you can calculate:
     <pre>
     angle = ||result||
     and
     axis = result//angle
     </pre>
     *
     * @param euler euler angles for [x, y, z]
     * @return
     */
    public static double[] createRotationVectorFromEulerAnglesXYZ(double[] euler) {
        double[][] r = createRotationXYZ(euler[0], euler[1], euler[2]);
        return Rotation.extractRotationVectorRodrigues(r);
    }

    public static double[][] createRotationFromEulerAngles(double[] eulerXYZ, EulerSequence seq) {
        if (seq.equals(EulerSequence.XYZ_ACTIVE)) {
            // intrinsic and active
            return createRotationXYZ(-eulerXYZ[0], -eulerXYZ[1], -eulerXYZ[2]);
        } else {
            // ZYX active
            return createRotationZYX(-eulerXYZ[0], -eulerXYZ[1], -eulerXYZ[2]);
        }
    }

    public static double[] createRotationVectorFromEulerAngles(double[] eulerXYZ, EulerSequence seq) {
        double[][] r = createRotationFromEulerAngles(eulerXYZ, seq);
        return Rotation.extractRotationVectorRodrigues(r);
    }

    /**
     * convert the eulerXYZ angles to angle axis representation.
     * NOTE that the resulting vector is the angle axis times the angle, so you can calculate:
     <pre>
     angle = ||result||
     and
     axis = result//angle
     </pre>
     *
     * @param eulerXYZ eulerXYZ angles for [x, y, z]
     * @return
     */
    public static double createAngleAxisFromEulerAnglesXYZ(double[] eulerXYZ, double[] outAxis) {
        double[][] r = createRotationXYZ(eulerXYZ[0], eulerXYZ[1], eulerXYZ[2]);
        double[] rotVec = Rotation.extractRotationVectorRodrigues(r);
        return Rotation.createAngleAxisFromRotationVector(rotVec, outAxis);
    }

    public static double[][] createRotationRodriguesFormula(double[] axis, double angle, boolean passive) {
        if (axis.length != 3) {
            throw new IllegalArgumentException("axis length must be 3");
        }
        double[] rotVec = createRotationVectorFromAngleAxis(axis, angle);
        return createRotationRodriguesFormula(rotVec, passive);
    }
    
    /**
     * given rotVec as an array of rotations about x, y, and z, calculate the
     * rotation matrix.  
     * essentially, excepting a small angle correction:
     *     R^⊤ = cosθ*I + sinθ*[rotVec]_× + (1−cosθ)*rotVec*rotVec^⊤
     *
     * NOTE: if computing the partial derivative of Rotation elsewhere, 
     * can use d(R(ω)*rotVec)/d(ω^T) = -[rotVec]_x (see Equation (2.35) of Szeliski 2010).
     * Also note that a + sign in front of the sine(θ) term is used for the
     * "passive" system of rotations.
     * <pre>
     * references:
     *    Dmitry Berenson https://github.com/robEllenberg/comps-plugins/blob/master/python/rodrigues.py
     *    Szeliski 2010 draft "Computer Vision: Algorithms and Applications"
     *    Rodriguez’s formula (Ayache 1989)
     * </pre>
     * @param rotVec [1x3] array of rotVec of rotations about x, y, and z
     * @param passive if true uses passive right-hand system of CCW rotation, else uses active left-hand CW rotations
     * @return rotation matrix [3X3]
     */
    public static double[][] createRotationRodriguesFormula(double[] rotVec, boolean passive) {
        if (rotVec.length != 3) {
            throw new IllegalArgumentException("rotVec length must be 3");
        }

        /*
        https://github.com/robEllenberg/comps-plugins/blob/master/python/rodrigues.py

        Copyright (c) 2010 Carnegie Mellon University and Intel Corporation
        #   Author: Dmitry Berenson <dberenso@cs.cmu.edu>
        #
        #   Redistribution and use in source and binary forms, with or without
        #   modification, are permitted provided that the following conditions are met:
        #
        #     * Redistributions of source code must retain the above copyright
        #       notice, this list of conditions and the following disclaimer.
        #     * Redistributions in binary form must reproduce the above copyright
        #       notice, this list of conditions and the following disclaimer in the
        #       documentation and/or other materials provided with the distribution.
        #     * Neither the name of Intel Corporation nor Carnegie Mellon University,
        #       nor the names of their contributors, may be used to endorse or
        #       promote products derived from this software without specific prior
        #       written permission.
        #
        #   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
        #   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
        #   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
        #   ARE DISCLAIMED. IN NO EVENT SHALL INTEL CORPORATION OR CARNEGIE MELLON
        #   UNIVERSITY BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
        #   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
        #   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
        #   OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
        #   WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
        #   OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
        #   ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

        # -*- coding: utf-8 -*-
        '''Rodrigues formula
        Input: 1x3 array of rotations about x, y, and z
        Output: 3x3 rotation matrix'''
        from numpy import array,mat,sin,cos,dot,eye
        from numpy.linalg import norm

        def rodrigues(r):
            def S(n):
                Sn = array([[0,-n[2],n[1]],[n[2],0,-n[0]],[-n[1],n[0],0]])
                return Sn
            theta = norm(r)
            if theta > 1e-30:
                n = r/theta
                Sn = S(n)
                R = eye(3) + sin(theta)*Sn + (1-cos(theta))*dot(Sn,Sn)
            else:
                Sr = S(r)
                theta2 = theta**2
                R = eye(3) + (1-theta2/6.)*Sr + (.5-theta2/24.)*dot(Sr,Sr)
            return mat(R)
        */

        double theta = MatrixUtil.lPSum(rotVec, 2);

        double[][] tmp1, tmp2;

        if (theta > 1e-30) {

            double[] n = Arrays.copyOf(rotVec, rotVec.length);
            MatrixUtil.multiply(n, 1./theta);
            // [3X3]
            double[][] sn = MatrixUtil.skewSymmetric(n);

            //R = eye(3) + sin(theta)*Sn + (1-cos(theta))*dot(Sn,Sn)
            tmp1 = MatrixUtil.copy(sn);
            //[3X3]
            double sa = Math.sin(theta);
            if (!passive) {
                sa *= -1;
            }
            MatrixUtil.multiply(tmp1, sa);

            //[3X3]
            tmp2 = MatrixUtil.multiply(sn, sn);
            MatrixUtil.multiply(tmp2, 1. - Math.cos(theta));

        } else {

            double[][] sr = MatrixUtil.skewSymmetric(rotVec);
            double theta2 = theta*theta;

            //R = eye(3) + (1-theta2/6.)*Sr + (.5-theta2/24.)*dot(Sr,Sr)
            tmp1 = MatrixUtil.copy(sr);
            MatrixUtil.multiply(tmp1, (1. - theta2)/6.);

            tmp2 = MatrixUtil.multiply(sr, sr);
            MatrixUtil.multiply(tmp2, 0.5 - theta2/24.);

        }

        double[][] r = MatrixUtil.createIdentityMatrix(3);

        int j;
        for (int i = 0; i < tmp1.length; ++i) {
            for (j = 0; j < tmp1[i].length; ++j) {
                r[i][j] += (tmp1[i][j] + tmp2[i][j]);
            }
        }

        return r;
    }

    public static class RodriguesRotation {
        /**
         * [3X3]
         */
        public double[][] r;

        /**
         * [9X3] or [3X9]
         * 9X3 returned for vector, else 3X9 for matrix
         */
        public double[][] dRdR;

        /**
         * [3X1]
         */
        public double[] rotVec;
    }

    /**
     * calculate the Rodrigues rotation vector from the given rotation matrix.
     * The method is ported from github repositories holding the Bouguet Matlab Toolbox code, rodrigues.m.
     *
     *      <pre>
     *      The Bouguet toolbox webpage is currently at http://robots.stanford.edu/cs223b04/JeanYvesCalib/
     *      and states that the source code is freely available.
     *      The github repositories with the forked Bouguet Matlab code do not have license
     *      information.
     *
     *      https://github.com/fragofer/TOOLBOX_calib
     *      and
     *      https://github.com/hunt0r/Bouguet_cam_cal_toolbox
     *
     *      rodrigues.m includes the comment: Copyright (c) March 1993 -- Pietro Perona, CalTech, before a brief
     *      changelist by Bouguet.
     *      </pre>
     *      Note that this is an ambiguous task.
     * @param in [3X3] rotation matrix
     * @return
     */
    public static RodriguesRotation extractRotationVectorRodriguesBouguet(double[][] in,
                                                                          boolean passive) throws NotConvergedException {

        //[m,n] = size(in);
        int m = in.length;
        int n = in[0].length;

        if (m != 3 || n != 3) {
            throw new IllegalArgumentException("in dimensions must be 3 X 3");
        }

        if (!passive) {
            in = MatrixUtil.transpose(in);
        }

        final double eps = 2.2204e-16;

        //%bigeps = 10e+4*eps;
        //bigeps = 10e+20*eps;
        double bigeps = 10e20 * eps;

        //norm(in' * in - eye(3)) < bigeps)...
        //                & (abs(det(in)-1) < bigeps))
        double check1 = MatrixUtil.spectralNorm(
                MatrixUtil.pointwiseSubtract(MatrixUtil.createATransposedTimesA(in),
                        MatrixUtil.createIdentityMatrix(3)));
        double check2 = Math.abs(MatrixUtil.determinant(in) - 1);
        if ((check1 >= bigeps) || check2 >= bigeps){
            throw new IllegalArgumentException("in does not appear to be a rotation matrix");
        }

        //R = in;
        double[][] R = in;

        //% project the rotation matrix to SO(3);
        //[U,S,V] = svd(R);
        //R = U*V';
        R = Rotation.orthonormalizeUsingSVD(R);

        //tr = (trace(R)-1)/2;
        double tr = (MatrixUtil.trace(R) - 1)/2.;
        //dtrdR = [1 0 0 0 1 0 0 0 1]/2; // [1X9]
        double[] dtrdR = new double[]{1/2., 0, 0, 0, 1/2., 0, 0, 0, 1/2.};
        //theta = real(acos(tr));
        double theta = Math.acos(tr);
        if (Double.isNaN(theta)) {
            throw new IllegalArgumentException("the rotation 'in' is not proper?");
        }

        double[] out = null;
        double[][] dout = null;

        //if sin(theta) >= 1e-4,
        if (Math.sin(theta) >= 1e-4) {
            //dthetadtr = -1 / sqrt(1 - tr ^ 2);
            double dthetadtr = -1./Math.sqrt(1. - tr*tr);

            //dthetadR = dthetadtr * dtrdR; [1X9]
            double[] dthetadR = Arrays.copyOf(dtrdR, dtrdR.length);
            MatrixUtil.multiply(dthetadR, dthetadtr);

            //%var1 = [vth; theta];
            //vth = 1 / (2 * sin(theta));
            double vth = 1./(2.*Math.sin(theta));
            //dvthdtheta = -vth * cos(theta) / sin(theta);
            double dvthdtheta = -vth * Math.cos(theta) / Math.sin(theta);
            //dvar1dtheta = [dvthdtheta; 1]; //[2X1]
            double[] dvar1dtheta = new double[]{dvthdtheta, 1};

            //dvar1dR = dvar1dtheta * dthetadR;  //[2X1][1X9] = [2X9]
            double[][] dvar1dR = MatrixUtil.outerProduct(dvar1dtheta, dthetadR);

            //om1 = [R(3, 2) - R(2, 3), R(1, 3) - R(3, 1), R(2, 1) - R(1, 2)]'; // [3X1]
            double[] om1 = new double[]{R[2][1] - R[1][2], R[0][2] - R[2][0], R[1][0] - R[0][1]};

            //          0 1   2  3 4 5 6  7 8
            //dom1dR = [0 0   0  0 0 1 0 -1 0; //[3 X 9]
            //          0 0  -1  0 0 0 1  0 0;
            //          0 1   0 -1 0 0 0  0 0];
            double[][] dom1dR = MatrixUtil.zeros(3, 9);
            dom1dR[0][5] = 1;
            dom1dR[0][7] = -1;
            dom1dR[1][2] = -1;
            dom1dR[1][6] = 1;
            dom1dR[2][1] = 1;
            dom1dR[2][3] = -1;

            //%var = [om1; vth; theta];
            //dvardR = [dom1dR; dvar1dR]; //[3X9] ; [2X9] => [5X9]
            double[][] dvardR = new double[5][];
            dvardR[0] = Arrays.copyOf(dom1dR[0], dom1dR[0].length);
            dvardR[1] = Arrays.copyOf(dom1dR[1], dom1dR[1].length);
            dvardR[2] = Arrays.copyOf(dom1dR[2], dom1dR[2].length);
            dvardR[3] = Arrays.copyOf(dvar1dR[0], dvar1dR[0].length);
            dvardR[4] = Arrays.copyOf(dvar1dR[1], dvar1dR[1].length);

            int i;
            //%var2 = [om; theta];
            //om = vth * om1; // vth*[3X1] = [3X1]
            double[] om = Arrays.copyOf(om1, om1.length);
            MatrixUtil.multiply(om, vth);
            //domdvar = [vth * eye(3) om1 zeros(3, 1)]; // [3X3] | [3X1] | [3X1] = [3 X 5]
            double[][] tmp = MatrixUtil.createIdentityMatrix(3);
            MatrixUtil.multiply(tmp, vth);
            double[][] domdvar = MatrixUtil.zeros(3, 5);
            for (i = 0; i < 3; ++i) {
                System.arraycopy(tmp[i], 0, domdvar[i], 0, tmp[i].length);
                domdvar[i][tmp[i].length] = om1[i];
            }
            //dthetadvar = [0 0 0 0 1]; // [1X5]
            double[] dthetadvar = new double[]{0, 0, 0, 0, 1};
            //dvar2dvar = [domdvar; dthetadvar];  // [3X5] ; [1X5] => [4 X 5]
            double[][] dvar2dvar = MatrixUtil.zeros(4, 5);
            for (i = 0; i < 3; ++i) {
                System.arraycopy(domdvar[i], 0, dvar2dvar[i], 0, domdvar[i].length);
            }
            System.arraycopy(dthetadvar, 0, dvar2dvar[3], 0, dthetadvar.length);

            //out = om * theta;
            out = Arrays.copyOf(om, om.length);
            MatrixUtil.multiply(out, theta);

            //domegadvar2 = [theta * eye(3) om]; // [3X3 | [3X1] ==> [3 X 4]
            double[][] domegadvar2 = MatrixUtil.zeros(3, 4);
            tmp = MatrixUtil.createIdentityMatrix(3);
            MatrixUtil.multiply(tmp, theta);
            for (i = 0; i < 3; ++i) {
                System.arraycopy(tmp[i], 0, domegadvar2[i], 0, tmp[i].length);
                domegadvar2[i][3] = om[i];
            }

            //dout = domegadvar2 * dvar2dvar * dvardR; // [3X4] [4X5] [5X9] = [3X9]
            dout = MatrixUtil.multiply(MatrixUtil.multiply(domegadvar2, dvar2dvar), dvardR);

        } else {
            if (tr > 0) {
                //out = [0 0 0]';
                out = new double[3];

                //         0 1     2    3 4   5   6    7    8
                // dout = [0 0     0    0 0  1/2   0  -1/2  0;
                //         0 0  -1/2    0 0    0 1/2     0  0;
                //         0 1/2   0 -1/2 0    0   0     0  0];

                dout = MatrixUtil.zeros(3, 9);
                dout[0][5] = 0.5;
                dout[0][7] = -0.5;
                dout[1][2] = -0.5;
                dout[1][6] = 0.5;
                dout[2][1] = 0.5;
                dout[2][3] = -0.5;

            } else {

                //% Solution by Mike Burl on Feb 27, 2007
                //% This is a better way to determine the signs of the
                //% entries of the rotation vector using a hash table on all
                //% the combinations of signs of a pairs of products (in the
                //% rotation matrix)

                //% Define hashvec and Smat
                //           0   1   2  3   4  5  6   7  8   9   10
                //hashvec = [0; -1; -3; -9; 9; 3; 1; 13; 5; -7; -11]; // [11 X 1]
                //            0       1       2       3      4      5      6      7      8
                //Smat = [1,1,1; 1,0,-1; 0,1,-1; 1,-1,0; 1,1,0; 0,1,1; 1,0,1; 1,1,1; 1,1,-1; //[11X 3]
                //1,-1,-1; 1,-1,1];
                //      9       10
                double[] hashvec = new double[]{0, -1, -3, -9, 9, 3, 1, 13, 5, -7, -11};
                double[][] Smat = new double[11][];
                Smat[0] = new double[]{1, 1, 1};
                Smat[1] = new double[]{1, 0, -1};
                Smat[2] = new double[]{0, 1, -1};
                Smat[3] = new double[]{1, -1, 0};
                Smat[4] = new double[]{1, 1, 0};
                Smat[5] = new double[]{0, 1, 1};
                Smat[6] = new double[]{1, 0, 1};
                Smat[7] = new double[]{1, 1, 1};
                Smat[8] = new double[]{1, 1, -1};
                Smat[9] = new double[]{1, -1, -1};
                Smat[10] = new double[]{1, -1, 1};

                // M = (R+eye(3,3))/2;
                double[][] M = MatrixUtil.pointwiseAdd(R, MatrixUtil.createIdentityMatrix(3));
                MatrixUtil.multiply(M, 0.5);

                //uabs = sqrt(M(1,1));
                //vabs = sqrt(M(2,2));
                //wabs = sqrt(M(3,3));
                double uabs = Math.sqrt(M[0][0]);
                double vabs = Math.sqrt(M[1][1]);
                double wabs = Math.sqrt(M[2][2]);

                //mvec = ([M(1,2), M(2,3), M(1,3)] + [M(2,1), M(3,2), M(3,1)])/2;
                double[] mvec0 = new double[]{M[0][1], M[1][2], M[0][2]};
                double[] mvec1 = new double[]{M[1][0], M[2][1], M[2][0]};
                double[] mvec = MatrixUtil.add(mvec0, mvec1);
                MatrixUtil.multiply(mvec, 0.5);

                //syn  = ((mvec > eps) - (mvec < -eps)); % robust sign() function
                double[] syn = new double[mvec.length];
                double t0, t1;
                int i;
                for (i = 0; i < 3; ++i) {
                    t0 = (mvec[i] > eps) ? 1 : 0;
                    t1 = (mvec[i] < -eps) ? 1 : 0;
                    syn[i] = t0 - t1;
                }
                //hash = syn * [9; 3; 1];  [1X3][3X1]=[1X1]
                double hash = MatrixUtil.innerProduct(syn, new double[]{9, 3, 1});

                //idx = find(hash == hashvec);
                // should not need to apply an offset of 1 as they are consistent use of indexes
                int idx = -1;
                for (i = 0; i < hashvec.length; ++i) {
                    if (hash == hashvec[i]) {
                        idx = i;
                    }
                }
                if (idx == -1) {
                    throw new IllegalStateException("ERROR: no solution found");
                }
                //svec = Smat(idx,:)';
                double[] svec = Smat[idx];

                //out = theta * [uabs; vabs; wabs] .* svec; // [3X1] .* [3X1] = [3X1]
                out = MatrixUtil.pointwiseMultiplication(new double[]{theta*uabs, theta*vabs, theta*wabs}, svec);
            }
        }

        RodriguesRotation rRot = new RodriguesRotation();
        rRot.rotVec = out;
        rRot.r = MatrixUtil.copy(in);
        rRot.dRdR = dout;

        return rRot;
    }

    /**
     * calculate the rotation matrix given the Rodrigues rotation vector.
     * The method is ported from github repositories holding the Bouguet Matlab Toolbox code, rodrigues.m.
     *
     * this is using passive transformations.
     *
     *      <pre>
     *      The Bouguet toolbox webpage is currently at http://robots.stanford.edu/cs223b04/JeanYvesCalib/
     *      and states that the source code is freely available.
     *      The github repositories with the forked Bouguet Matlab code do not have license
     *      information.
     *
     *      https://github.com/fragofer/TOOLBOX_calib
     *      and
     *      https://github.com/hunt0r/Bouguet_cam_cal_toolbox
     *
     *      rodrigues.m includes the comment: Copyright (c) March 1993 -- Pietro Perona, CalTech, before a brief
     *      changelist by Bouguet.
     *
     *      </pre>
     * @param rotVec [3X1] rotation vector
     * @return
     */
    public static RodriguesRotation createRotationRodriguesBouguet(double[] rotVec, boolean passive) {

        if (rotVec.length != 3) {
            throw new IllegalArgumentException("rotVec length must be 3");
        }

        if (!passive) {
            rotVec = Arrays.copyOf(rotVec, rotVec.length);
            MatrixUtil.multiply(rotVec, -1);
        }

        final double eps = 2.2204e-16;

        //[m,n] = size(rotVec);
        int m = rotVec.length;
        //int n = 1;
        //%bigeps = 10e+4*eps;
        //bigeps = 10e+20*eps;

        double bigeps = 10e20 * eps;

        double[][] R;
        double[][] dRdin;
        //theta = norm(rotVec);
        double theta = MatrixUtil.lPSum(rotVec, 2);
        //if theta < eps
        if (theta < eps) {
            //R = eye(3);
            R = MatrixUtil.createIdentityMatrix(3);
            //dRdin = [0 0 0;
            //0 0 1;
            //0 -1 0;
            //0 0 -1;
            //0 0 0;
            //1 0 0;
            //0 1 0;
            //-1 0 0;
            //0 0 0];
            dRdin = new double[9][];
            dRdin[0] = new double[]{0, 0, 0};
            dRdin[1] = new double[]{0, 0, 1};
            dRdin[2] = new double[]{0, -1, 0};
            dRdin[3] = new double[]{0, 0, -1};
            dRdin[4] = new double[]{0, 0, 0};
            dRdin[5] = new double[]{1, 0, 0};
            dRdin[6] = new double[]{0, 1, 0};
            dRdin[7] = new double[]{-1, 0, 0};
            dRdin[8] = new double[]{0, 0, 0};

            //out = R;
            //dout = dRdin;
            RodriguesRotation rRot = new RodriguesRotation();
            rRot.r = R;
            rRot.dRdR = dRdin;
            rRot.rotVec = Arrays.copyOf(rotVec, rotVec.length);

            return rRot;
        }

        //%m3 = [rotVec,theta]

        //dm3din = [eye(3);rotVec'/theta];
        double[][] dm3din = new double[4][]; // [4X3]
        dm3din[0] = new double[]{1, 0, 0};
        dm3din[1] = new double[]{0, 1, 0};
        dm3din[2] = new double[]{0, 0, 1};
        dm3din[3] = Arrays.copyOf(rotVec, rotVec.length);
        MatrixUtil.multiply(dm3din[3], 1./theta);

        //omega = rotVec/theta;
        double[] omega = Arrays.copyOf(dm3din[3], dm3din[3].length); // [3X1]

        //%m2 = [omega;theta]

        double invTheta = 1./theta;
        double invTheta2 = invTheta*invTheta;
        //dm2dm3 = [eye(3)/theta -rotVec/theta^2; zeros(1,3) 1];// [3X3] | [3X1] ;
        double[][] dm2dm3 = new double[4][];
        dm2dm3[0] = new double[]{1*invTheta, 0, 0, rotVec[0]*-invTheta2};
        dm2dm3[1] = new double[]{0, 1*invTheta, 0, rotVec[1]*-invTheta2};
        dm2dm3[2] = new double[]{0, 0, 1*invTheta, rotVec[2]*-invTheta2};
        dm2dm3[3] = new double[]{0, 0, 0, 1};

        //alpha = cos(theta);
        //beta = sin(theta);
        //gamma = 1-cos(theta);
        //omegav=[[0 -omega(3) omega(2)];[omega(3) 0 -omega(1)];[-omega(2) omega(1) 0 ]];
        //A = omega*omega';
        double alpha = Math.cos(theta);
        double beta = Math.sin(theta);
        double gamma = 1. - alpha;
        double[][] omegav = MatrixUtil.skewSymmetric(omega);
        double[][] A = MatrixUtil.outerProduct(omega, omega); //[3X3]

        //%m1 = [alpha;beta;gamma;omegav;A];

        //dm1dm2 = zeros(21,4);
        //dm1dm2(1,4) = -sin(theta);
        //dm1dm2(2,4) = cos(theta);
        //dm1dm2(3,4) = sin(theta);
        //dm1dm2(4:12,1:3) = [0 0 0  0 0 1 0 -1 0;
        //                    0 0 -1 0 0 0 1 0 0;
        //                    0 1 0 -1 0 0 0 0 0]';
        double[][] dm1dm2 = MatrixUtil.zeros(21, 4);
        dm1dm2[0][3] = -beta;
        dm1dm2[1][3] = alpha;
        dm1dm2[2][3] = beta;
        System.arraycopy(new double[]{0,0,0}, 0, dm1dm2[3], 0, 3);
        System.arraycopy(new double[]{0,0,1}, 0, dm1dm2[4], 0, 3);
        System.arraycopy(new double[]{0,-1,0}, 0, dm1dm2[5], 0, 3);
        System.arraycopy(new double[]{0,0,-1}, 0, dm1dm2[6], 0, 3);
        System.arraycopy(new double[]{0,0,0}, 0, dm1dm2[7], 0, 3);
        System.arraycopy(new double[]{1,0,0}, 0, dm1dm2[8], 0, 3);
        System.arraycopy(new double[]{0,1,0}, 0, dm1dm2[9], 0, 3);
        System.arraycopy(new double[]{-1,0,0}, 0, dm1dm2[10], 0, 3);
        System.arraycopy(new double[]{0,0,0}, 0, dm1dm2[11], 0, 3);

        //w1 = omega(1);
        //w2 = omega(2);
        //w3 = omega(3);
        double w1 = omega[0];
        double w2 = omega[1];
        double w3 = omega[2];
        //dm1dm2(13:21,1) = [2*w1;w2;w3;w2;0;0;w3;0;0];
        //dm1dm2(13: 21,2) = [0;w1;0;w1;2*w2;w3;0;w3;0];
        //dm1dm2(13:21,3) = [0;0;w1;0;0;w2;w1;w2;2*w3];
        int i;
        double[] c0 = new double[]{2*w1,w2,w3,w2,0,0,w3,0,0};
        double[] c1 = new double[]{0,w1,0,w1,2*w2,w3,0,w3,0};
        double[] c2 = new double[]{0,0,w1,0,0,w2,w1,w2,2*w3};
        for (i = 12; i <= 20; ++i) {
            dm1dm2[i][0] = c0[i - 12];
            dm1dm2[i][1] = c1[i - 12];
            dm1dm2[i][2] = c2[i - 12];
        }

        //R = eye(3)*alpha + omegav*beta + A*gamma; // [3X3] + [3X3] + [3X3]
        R = MatrixUtil.createIdentityMatrix(3);
        MatrixUtil.multiply(R, alpha);
        double[][] t1 = MatrixUtil.copy(omegav);
        MatrixUtil.multiply(t1, beta);
        double[][] t2 = MatrixUtil.copy(A);
        MatrixUtil.multiply(t2, gamma);
        R = MatrixUtil.pointwiseAdd(R, t1);
        R = MatrixUtil.pointwiseAdd(R, t2);

        //dRdm1 = zeros(9,21);
        //dRdm1([1 5 9],1) = ones(3,1);
        //dRdm1(:,2) = omegav(:);
        //dRdm1(:,4:12) = beta*eye(9);
        //dRdm1(:,3) = A(:);
        //dRdm1(:,13:21) = gamma*eye(9);
        double[][] dRdm1 = MatrixUtil.zeros(9, 21);
        dRdm1[0][0] = 1;
        dRdm1[4][0] = 1;
        dRdm1[8][0] = 1;
        double[] omegavStack = MatrixUtil.stack(omegav);
        double[][] iB = MatrixUtil.createIdentityMatrix(9);
        MatrixUtil.multiply(iB, beta);
        double[] AStack = MatrixUtil.stack(A);
        double[][] iG = MatrixUtil.createIdentityMatrix(9);
        MatrixUtil.multiply(iG, gamma);
        for (i = 0; i < dRdm1.length; ++i) {
            dRdm1[i][1] = omegavStack[i];
            dRdm1[i][2] = AStack[i];
            System.arraycopy(iB[i], 0, dRdm1[i], 3, iB[i].length);
            System.arraycopy(iG[i], 0, dRdm1[i], 12, iG[i].length);
        }

        //dRdin = dRdm1 * dm1dm2 * dm2dm3 * dm3din;
        //        [9X21]  [21X4]  [4X4]     [4X3] = [9X3]
        dRdin = MatrixUtil.multiply(dRdm1, dm1dm2);
        dRdin = MatrixUtil.multiply(dRdin, dm2dm3);
        dRdin = MatrixUtil.multiply(dRdin, dm3din);

        RodriguesRotation rRot = new RodriguesRotation();
        rRot.r = R;
        rRot.dRdR = dRdin;
        rRot.rotVec = Arrays.copyOf(rotVec, rotVec.length);

        return rRot;
    }
    
    /**
     * determine the rotation between measurements x1 and x2 when both datasets
     * have the same center, that is, there is no translation between them,
     * only rotation.
     * (see Golub & van Loan "Matrix Computations" 11.12.4,
     * Szeliski 2010, Sect 6.1.5).
     * @param x1 a set of measurements having same center as x2, that is,
     * there is no translation between them, only rotation.
     * the expected format is nData X nDimensions.
     * @param x2 another set of measurements having same center as x1, that is,
     * there is no translation between them, only rotation.
     * the expected format is nData X nDimensions.
     * @return difference in rotation between x2 and x1.  this will be ~ the identity matrix for no difference.
     * @throws no.uib.cipr.matrix.NotConvergedException
    */
    public static double[][] procrustesAlgorithmForRotation(double[][] x1, double[][] x2) 
        throws NotConvergedException {
        int m = x1.length;
        int p = x1[0].length;
        if (x2.length != m || x2[0].length != p) {
            throw new IllegalArgumentException("x1 and x2 must have same sizes");
        }
        // minimize || x1 - x2*Q ||_F
        //    subject to Q^T * Q = I_P
        double[][] c = MatrixUtil.multiply(MatrixUtil.transpose(x2), x1);
        MatrixUtil.SVDProducts svdC = MatrixUtil.performSVD(c);

        double[][] q = MatrixUtil.multiply(svdC.u, svdC.vT);
        return q;
    }
    
    /**
     * extract the euler rotation angles from the given rotation matrix assumed
     * to be a result of R_yxz = R_z(theta_z) * R_x(theta_x) * R_y(theta_y)
     * (aka 2-1-3 angle set?).
     * Note, this way of extracting of the angles is 
     * ambiguous (there are more than one angle combination sets that will
       result in the same matrix).
     <pre>
     the method is from equation (37) from 
     lecture notes of Gordon Wetzstein at Stanford University,
     EE 267 Virtual Reality, "Course Notes: 6-DOF Pose Tracking with the VRduino",
     https://stanford.edu/class/ee267/notes/ee267_notes_tracking.pdf

     It uses the passive, right-hand, CCW transformations system.

     </pre>
     * @param r  rotation matrix assumed
     * to be a result of R_yxz = R_z(theta_z)  * R_x(theta_x) * R_y(theta_y).
     * @return theta_x, theta_y, theta_z
     * which, respectively have ranges [0, 2*pi], [0, pi], and [0, 2*pi]
     */
    public static double[] extractRotationAxisFromZXY(double[][] r) {
        if (r.length != 3 || r[0].length != 3) {
            throw new IllegalArgumentException("r must be 3x3");
        }
        double thetaX = Math.asin(r[2][1]);
        double thetaY = Math.atan2(-r[2][0], r[2][2]);
        double thetaZ = Math.atan2(-r[0][1], r[1][1]);
        return new double[]{thetaX, thetaY, thetaZ};
    }
    
    /**
     * extract euler rotation angles from a rotation matrix which has been built following
     * the convention of R_xyz =  R(theta_Z) * R(theta_Y) * R(theta_X)
     * This method uses a passive (right-hand) rotation system with CCW (=CC) rotations.
     * @param r ZYX rotation matrix
     * @return array of theta_x, theta_y, theta_z
     * which, respectively have ranges [0, 2*pi], [0, pi], and [0, 2*pi]
     */
    public static double[] extractThetaFromZYX(double[][] r) {
        if (r.length != 3 || r[0].length != 3) {
            throw new IllegalArgumentException("r must be 3x3");
        }
        double[] out = new double[3];
        extractThetaFromZYX(r, out);
        return out;
    }

    /**
     *
     * @param r ZYX rotation matrix
     * @param passive
     @return array of theta_x, theta_y, theta_z
     which, respectively have ranges [0, 2*pi], [0, pi], and [0, 2*pi]

     */
    public static double[] extractThetaFromZYX(double[][] r, boolean passive) {
        double[] thetas = extractThetaFromZYX(r);
        if (!passive) MatrixUtil.multiply(thetas, -1);
        return thetas;
    }
    
    /**
     * extract euler rotation angles from a rotation matrix which has been built following
     * the convention of R(theta_Z) * R(theta_Y) * R(theta_X).
     * This method uses a passive (right-hand) rotation system with CCW (=CC) rotations.
     * @param r
     * @param out output variable to hold theta_x, theta_y, and theta_z,
     * which, respectively have ranges [0, 2*pi], [0, pi], and [0, 2*pi]
     */
    public static void extractThetaFromZYX(double[][] r, double[] out) {
        if (r.length != 3 || r[0].length != 3) {
            throw new IllegalArgumentException("r must be 3x3");
        }
        /*
        using euler angles
          cc rotation about z-axis (yaw):   cc about the y-axis (pitch):    cc about x-axis (roll):    
            | cos φ   -sin φ    0 |          |  cos ψ    0  sin ψ |         |    1       0       0 |  
            | sin φ    cos φ    0 |          |      0    1      0 |         |    0   cos θ  -sin θ |  
            |     0        0    1 |          | -sin ψ    0  cos ψ |         |    0   sin θ   cos θ |  

        = | (cos φ * cos ψ)   (-sin φ)   (cos φ * sin ψ) |  * | 1       0       0 |
          | (sin φ * cos ψ)   ( cos φ)   (sin φ * sin ψ) |    | 0   cos θ   sin θ |
          | (-sin ψ)          (   0  )   (   cos ψ )     |    | 0  -sin θ   cos θ |

        = | (cos Z * cos Y)   (-sin Z * cos X + cos Z * sin Y * sin X)   ( sin Z * sin X + cos Z * sin Y * cos X) |
          | (sin Z * cos Y)   ( cos Z * cos X + sin Z * sin Y * sin X)   (-cos Z * sin X + sin Z * sin Y * cos X)  |
          | (-sin Y)          ( cos Y * sin X )                          (cos Y * cos X)                          |

        r20 = -sin ψ  ==> ψ = theta_y = -Math.asin(r20)
        r21/r22 = ( cos ψ * sin θ )/(cos ψ * cos θ) = (sin θ)/(cos θ) = tan(θ) ==> θ = theta_x = Math.atan2(r21, r22)
        r10/r00 = (sin φ * cos ψ) / (cos φ * cos ψ) = Math.atan2(r10, r00)
        */
        
        //        θ      ψ       φ
        double thetaX, thetaY, thetaZ;
        
        //ψ
        double d = r[2][1]*r[2][1] + r[2][2]*r[2][2];
        if (d == 0) {
            thetaY = -Math.asin(r[2][0]);
        } else {
            // Y = atan( -r[2][0] / sqrt(r[2][1]*r[2][1] + r[2][2]*r[2][2]) )
            //         (  +sin Y  / sqrt( (cosY)^2 * (sinX)^2 + (cosY)^2 * (cosX)^2 )
            //         (  +sin Y  / sqrt( (cosY)^2 * 1 )
            thetaY = Math.atan2(-r[2][0], Math.sqrt(d));
        }
        
        //θ
        thetaZ = Double.NEGATIVE_INFINITY;
        if (r[2][2] != 0) {
            thetaX = Math.atan2(r[2][1], r[2][2]);
        } else {
            // cos ψ==0 or/and cos θ==0
            if (r[2][1] != 0) {
                // then cos ψ != 0  and cos θ==0
                thetaX = Math.asin(r[2][1]/Math.cos(thetaY));
            } else {
                // else cos ψ == 0 and possibly cos θ==0
                // need thetaZ solved
                if (r[0][0] != 0) {
                    thetaZ = Math.atan2(r[1][0], r[0][0]);
                    //CALC θ(thetaX), knowing φ(thetaZ) and ψ(thetaY)
                    double cPsi = Math.cos(thetaY);
                    if (cPsi != 0) {
                        // can use r[2][1] or r[2][2] for simplest:
                        thetaX = Math.asin(r[2][1]/cPsi);
                    } else {
                        // can use r[0][1], r[0][2], r[1][1], or r[1][2] or combination
                        /*r[0][1] = (-sφ * cos θ + cφ * sψ * sin θ)
                          r[0][2] =  (sφ * sin θ + cφ * sψ * cos θ)
                          r[1][1] = ( cφ * cos θ + sφ * sψ * sin θ)
                          r[1][2] = (-cφ * sin θ + sφ * sψ * cos θ)
                        
                        looking for ways to factor one or more of the 4 equations for the 1 unknown θ
                        rewritten:
                        r[0][1] = cos θ * -sφ      +  sin θ * cφ * sψ 
                        r[0][2] = cos θ * cφ * sψ  +  sin θ * sφ
                        r[1][1] = cos θ * cφ       +  sin θ * sφ * sψ 
                        r[1][2] = cos θ * sφ * sψ  +  sin θ * -cφ
                        add all   : cos θ * (-sφ + cφ * sψ + cφ + sφ * sψ) + sin θ * (cφ * sψ + sφ + sφ * sψ + -cφ)
                        rewritten : cos θ * (cφ - sφ + (cφ * sψ) + (sφ * sψ)) + sin θ * (-cφ + sφ + (cφ * sψ) + (sφ * sψ))
                                  : cos θ * (cφ - sφ + (cφ * sψ) + (sφ * sψ)) + sin θ * (-cφ + sφ + (cφ * sψ) + (sφ * sψ))
                        let a1 = (cφ * sψ) + (sφ * sψ)
                        let a0 = cφ - sφ
                                  : cos θ * (a0 + a1) + sin θ * (-a0 + a1)
                                  : a0*(cos θ - sin θ) + a1*(cos θ + sin θ)
                                  : a0*cosθ - a0*sinθ + a1*cosθ + a1*sinθ
                        
                        squared:  (a0*cosθ - a0*sinθ + a1*cosθ + a1*sinθ)*(a0*cosθ - a0*sinθ + a1*cosθ + a1*sinθ)
                                  = a0^2*(cosθ)^2 - a0^2*cosθ*sinθ + a0*a1*(cosθ)^2 + a0*a1*cosθ*sinθ
                                    - a0^2*cosθ*sinθ + a0^2*(sinθ)^2 - a0*a1*cosθ*sinθ - a0*a1*(sinθ)^2
                                    + aθ*a1*(cosθ)^2 - a0*a1*cosθ*sinθ + a1^2*(cosθ)^2 + a1^2*cosθ*sinθ
                                    + a0*a1*cosθ*sinθ - a0*a1*(sinθ)^2 + a1^2*cosθ*sinθ + a1^2*(sinθ)^2
                                  = a0^2*(cosθ)^2 + a0^2*(sinθ)^2
                                    + a1^2*(cosθ)^2 + a1^2*(sinθ)^2
                                    - a0^2*cosθ*sinθ - a0^2*cosθ*sinθ
                                    + a0*a1*(cosθ)^2 + aθ*a1*(cosθ)^2 - a0*a1*(sinθ)^2 - a0*a1*(sinθ)^2
                                    + a0*a1*cosθ*sinθ - a0*a1*cosθ*sinθ + a0*a1*cosθ*sinθ - a0*a1*cosθ*sinθ 
                                    + a1^2*cosθ*sinθ + a1^2*cosθ*sinθ
                                  = a0^2 + a1^2
                                    - 2*a0^2*( cosθ*sinθ )               <==== sin(2x) = 2 sin(x) cos(x)
                                    + 2*a0*a1*( (cosθ)^2 - (sinθ)^2 )    <==== cos(2x) = cos2(x) – sin2(x)
                                    + 2*a1^2*( cosθ*sinθ )
                                  = a0^2 + a1^2
                                     - a0^2*sin(2θ) + a1^2*sin(2θ) + 2*a0*a1*cos(2θ)
                                  = a0^2 + a1^2 - sin(2θ)*( a0^2 - a1^2 ) + cos(2θ) * ( 2*a0*a1 )
                        */
                        throw new UnsupportedOperationException("There are "
                                + "0's in the rotation matrix, so factoring of "
                                + "more than one exponential variable is needed."
                                + "This case is not yet implemented.");
                    }
                } else {
                    /*
                    need φ, missing θ, knowing ψ
                    r[0][1] : (-sin φ * cos θ + cos φ * sψ * sin θ)
                    r[0][2] : ( sin φ * sin θ + cos φ * sψ * cos θ)
                    r[1][2] : (-cos φ * sin θ + sin φ * sψ * cos θ)
                    ==> 2 unknowns and 3 equations.
                    */
                    throw new UnsupportedOperationException("There are "
                                + "0's in the rotation matrix, so factoring of "
                                + "more than one exponential variable is needed."
                                + "This case is not yet implemented.");
                }
            }
        }
        
        //φ
        if (thetaZ == Double.NEGATIVE_INFINITY) {
            if (r[0][0] != 0) {
                thetaZ = Math.atan2(r[1][0], r[0][0]);
            } else {
                /* if r[0][0]==0, then so are r[1][0], r[2][1], r[2][2].
                 have θ and ψ
                r[0][1] : (-sin φ * cθ + cos φ * sψ * sθ)
                r[0][2] : ( sin φ * sθ + cos φ * sψ * cθ)
                r[1][2] : (-cos φ * sθ + sin φ * sψ * cθ)
                ==> 1 unknown and 3 equations
                */
                throw new UnsupportedOperationException("There are "
                                + "0's in the rotation matrix, so factoring of "
                                + "more than one exponential variable is needed."
                                + "This case is not yet implemented.");
            }
        }
        
        out[0] = thetaX;
        out[1] = thetaY;
        out[2] = thetaZ;
    }
    
    /**
     * extract euler rotation angles from a rotation matrix which has been built following
     * the convention of R_X(theta_X) * R_Y(theta_Y) * R_Z(theta_Z).
     * This method uses a passive (right-hand) rotation system with CCW (=CC) rotations.
     * @param r
     * @return euler [theta_X, theta_Y, theta_Z] angles extracted from the rotation matrix under assumption
     * that the rotation matrix was constructed with multiplication order
     * R_X(theta_X) * R_Y(theta_Y) * R_Z(theta_Z).
     * theta_x, theta_y, and theta_z, respectively have ranges [0, 2*pi], [0, pi], and [0, 2*pi]
     */
    public static double[] extractThetaFromXYZ(double[][] r) {
        if (r.length != 3 || r[0].length != 3) {
            throw new IllegalArgumentException("r must be 3x3");
        }
        double[] out = new double[3];
        extractThetaFromXYZ(r, out);
        return out;
    }

    /**
     * extract the euler angles from the XYZ rotation matrix
     * @param r
     * @param passive
     * @return Euler angles in order X, Y, Z
     */
    public static double[] extractThetaFromXYZ(double[][] r, boolean passive) {
        if (r.length != 3 || r[0].length != 3) {
            throw new IllegalArgumentException("r must be 3x3");
        }
        double[] thetas = extractThetaFromXYZ(r);
        if (!passive) MatrixUtil.multiply(thetas, -1);
        return thetas;
    }


    
    /**
     * extract euler rotation angles from a rotation matrix which has been built following
     * the convention of R_X(theta_X) * R_Y(theta_Y) * R_Z(theta_Z).
     * This method uses a passive (right-hand) rotation system with CCW (=CC) rotations.
     * @param r
     * @param out output array of size 3 used to place euler angles extracted 
     * from the rotation matrix under assumption
     * that the rotation matrix was constructed with multiplication order
     * R_X(theta_X) * R_Y(theta_Y) * R_Z(theta_Z).
     */
    public static void extractThetaFromXYZ(double[][] r, double[] out) {
        if (r.length != 3 || r[0].length != 3) {
            throw new IllegalArgumentException("r must be 3x3");
        }
        /*
        using euler angles and passive (RH) system
          cc rotation about z-axis (yaw):   cc about the y-axis (pitch):    cc about x-axis (roll):    
            | cos φ   -sin φ    0 |          |  cos ψ    0  sin ψ |         |    1       0       0 |  
            | sin φ    cos φ    0 |          |      0    1      0 |         |    0   cos θ  -sin θ |  
            |     0        0    1 |          | -sin ψ    0  cos ψ |         |    0   sin θ   cos θ |  

        from Appendix A of 
        "Euler Angles and Quaternions and Transformations", NASA Mission Planning
         and Analysis Division, 1977, Shuttle Program.
        
        R_X(theta_X)*R_Y(theta_Y)*R_Z(theta_Z) 
        = |  cosY*cosZ                      -cosY*sinZ                    sinY       |
          |  sinX*sinY*cosZ + cosX*sinZ     -sinX*sinY*sinZ + cosX*cosZ   -sinX*cosY |
          |  -cosX*sinY*cosZ + sinX*sinZ     cosX*sinY*sinZ + sinX*cosZ   cosX*cosY  |
        
        q1 = -sin(X/2)*sin(Y/2)*sin(Z/2) + cos(X/2)*cos(Y/2)*cos(Z/2)
        q2 =  sin(X/2)*cos(Y/2)*cos(Z/2) + sin(Y/2)*sin(Z/2)*cos(X/2)
        q3 = -sin(X/2)*sin(Z/2)*cos(Y/2) + sin(Y/2)*cos(X/2)*cos(Z/2)
        q4 =  sin(X/2)*sin(Y/2)*cos(Z/2) + sin(Z/2)*cos(X/2)*cos(Y/2)
        
        X = theta1 = math.atan2(-r[1][2], r[2][2])
        Y = theta2 = math.atan2( r[0][2], sqrt(1 - r[0][2]*r[0][2]) )
        Z = theta3 = math.atan2(-r[0][1], r[0][0])        
        */
        //        θ      ψ       φ
        double thetaX, thetaY, thetaZ;
        
        if (r[2][2] != 0) {
            thetaX = Math.atan2(-r[1][2], r[2][2]);
        } else {
            throw new UnsupportedOperationException("There are "
                                + "0's in the rotation matrix, so factoring of "
                                + "more than one exponential variable is needed."
                                + "This case is not yet implemented.");
        }
        
        double d = 1. - r[0][2]*r[0][2];
        if (d != 0) {
            thetaY = Math.atan2( r[0][2], Math.sqrt(1. - r[0][2]*r[0][2]) );
        } else {
            throw new UnsupportedOperationException("There are "
                                + "0's in the rotation matrix, so factoring of "
                                + "more than one exponential variable is needed."
                                + "This case is not yet implemented.");
        }
    
        if (r[0][0] != 0) {
            thetaZ = Math.atan2(-r[0][1], r[0][0]);
        } else {
            throw new UnsupportedOperationException("There are "
                                + "0's in the rotation matrix, so factoring of "
                                + "more than one exponential variable is needed."
                                + "This case is not yet implemented.");
        }
    
        out[0] = thetaX;
        out[1] = thetaY;
        out[2] = thetaZ;
        
    }

    /**
     * another method to extract the Rodrigues vector (angle and axis)
     * from the given
     * rotation matrix.  it's an ambiguous task.
     *
     * The rotations are by default passive transformations.

     <pre>
     the method is from  
     lecture notes of Pradit Mittrapiyanuruk at Perdue University,
     ECE 661 Robot Vision Laboratory,
     https://engineering.purdue.edu/kak/computervision/ECE661_Fall2012/homework/hw5_LM_handout.pdf
     </pre>
     * @param r rotation matrix
     * @return Rodriques vector.  the axis and angle representation from this
     * can be constructed as angle of rotation = r_vec/||r_vec|| and angle = ||r_vec||
     * where r_vec is the rodrigues vector.
     */
    public static double[] extractRotationVectorRodrigues(double[][] r) {
        if (r.length != 3 || r[0].length != 3) {
            throw new IllegalArgumentException("r must be 3x3");
        }

        double det = MatrixUtil.determinant(r);
        // numerical precision
        det = Math.round(det * 1E11)/1E11;
        if (Math.abs(det - 1) > 1E-7) {
            throw new IllegalArgumentException("expecting det(r) = 1 for proper rotation matrix");
        }

        /*
        https://github.com/robEllenberg/comps-plugins/blob/master/python/rodrigues.py
        compare to:
        def rodrigues(r):
            def S(n):
                Sn = array([[0,-n[2],n[1]],[n[2],0,-n[0]],[-n[1],n[0],0]])
                return Sn
            theta = norm(r)
            if theta > 1e-30:
                n = r/theta
                Sn = S(n)
                R = eye(3) + sin(theta)*Sn + (1-cos(theta))*dot(Sn,Sn)
            else:
                Sr = S(r)
                theta2 = theta**2
                R = eye(3) + (1-theta2/6.)*Sr + (.5-theta2/24.)*dot(Sr,Sr)
            return mat(R)
        */

        double tol = 1.e-7;
        
        double traceR = MatrixUtil.trace(r);
        System.out.printf("trace(R)=%.3e\n", traceR);
        
        int i;
        
        if ((traceR + 1) < tol ) {
            // trace(R) == -1, and theta = Math.PI
            /*
            http://www2.ece.ohio-state.edu/~zhang/RoboticsClass/docs/LN3_RotationalMotion.pdf
            
            w is the unit vector representing the axis of rotation.
            w is one of the following 3:
            
                (1./Math.sqrt(2.*(1.+r[2][2]))) * [r[0][2], r[1][2], 1 + r[2][2]]
            or  (1./Math.sqrt(2.*(1.+r[1][1]))) * [r[0][1], 1 + r[1][1], r[2][1]]
            or  (1./Math.sqrt(2.*(1.+r[0][0]))) * [1 + r[0][0], r[1][0], r[2][0]]
            */
            double t1 = (1./Math.sqrt(2.*(1.+r[2][2])));
            double[] w = new double[3];
            w[0] = t1 * r[0][2];
            w[1] = t1 * r[1][2];
            w[2] = t1 * (1. + r[2][2]);
            return w;
        }
        
        /*
        theta=acos((trace(R)-1)/2);
        w=(theta/(2*sin(theta)))*[R(3,2)-R(2,3); R(1,3)-R(3,1); R(2,1)-R(1,2)];
        wx=w(1);
        wy=w(2);
        wz=w(3);
        */
        
        // theta in range [0, Math.PI)
        // Math.acos argument must be >= 0 and <= 1
        // precision errors possibly result in > 1 so rounding here assuming machine precision 1E-11
        double arg = Math.round(0.5*(traceR - 1.)*1E11)/1E11;
        double theta = Math.acos(arg);
        
        // but // http://www2.ece.ohio-state.edu/~zhang/RoboticsClass/docs/LN3_RotationalMotion.pdf
        // use t1 without theta factor:
        double t1 = (theta == 0) ? 1 : 0.5*theta/Math.sin(theta);
        double[] w = new double[3];
        w[0] = t1*(r[2][1] - r[1][2]);
        w[1] = t1*(r[0][2] - r[2][0]);
        w[2] = t1*(r[1][0] - r[0][1]);
        
        return w;
    }
    
    /**
    calculate the distance metric between quaternions, using vector inner product
    <pre>
    eqn (19) of J Math Imaging Vis (2009) 35: 155–164 DOI 10.1007/s10851-009-0161-2
    Metrics for 3D Rotations: Comparison and Analysis Du Q. Huynh
    </pre>
     * @param q1 a quaternion
     * @param q2 another quaternion
     * @return distance metric between 2 quaternions.  the range of 
     * values will be [0,π/2] (radians).
     */
    public static double distanceBetweenQuaternions(double[] q1, double[] q2) {
        double p = Math.abs(MatrixUtil.innerProduct(q1, q2));
        double d = Math.acos(p);
        return d;
    }

    /**
     * estimate the
     * <pre>
     *     reference:
     *     Huynh 2009, J Math Imaging Vis, 35, 155-164, eqn (21)
     * </pre>
     * @param r1 a 3X3 rotation matrix
     * @param r2 a 3X3 rotation matrix
     * @param useFrobenius if true, uses Frobenius norm internally, else uses the spectral norm.
     * @return return the distance defined as ∥ I − r1 * r2 ∥_F if useFrobenius is true.
     * In that case the distance is within the range [0, 2sqrt(2)].
     * If useFrobenius is false, the result is spectralNorm(I − r1 * r2) which results in a distance in the range [0, 2].
     */
    public static double distanceUsingRigidBodyDisplacements(double[][] r1, double[][] r2, boolean useFrobenius) throws NotConvergedException {

        if (r1.length != 3 || r1[0].length != 3) {
            throw new IllegalArgumentException("r1 must be 3X3");
        }
        if (r2.length != 3 || r2[0].length != 3) {
            throw new IllegalArgumentException("r2 must be 3X3");
        }
        //Φ5(R1,R2)= ∥I−R1R2∥F,
        double[][] t = MatrixUtil.multiply(r1, r2);
        t = MatrixUtil.pointwiseSubtract(MatrixUtil.createIdentityMatrix(3), t);
        if (useFrobenius) {
            return MatrixUtil.frobeniusNorm(t);
        }
        return MatrixUtil.spectralNorm(t);
    }
    
    /*
    calculate the distance measure between 2 Euclidean transformations of the
    given quaternions, using vector inner product.
    <pre>
    eqn (20) of J Math Imaging Vis (2009) 35: 155–164 DOI 10.1007/s10851-009-0161-2
    Metrics for 3D Rotations: Comparison and Analysis Du Q. Huynh
    </pre>
     * @param q1 a quaternion
     * @param q2 another quaternion
     * @return the distance measure between two Euclidean transformations.  the range of 
     * values will be [0, 1] (radians).
     */
    /*public static double distanceBetweenQuaternionEuclideanTransformations(double[] q1, double[] q2) {
        double p = Math.abs(MatrixUtil.innerProduct(q1, q2));
        return 1 - p;
    }*/
    
    /**
     * given a quaternion of Barfoot format [vector scalar], form the left-hand compound operator
     * (symbol is superscript +)
     * <pre>
     * from Barfoot, Forbes, & Furgale 2010, "Pose estimation using linearized 
     * rotations and quaternion algebra", Acta Astronautica (2010), doi:10.1016/j.actaastro.2010.06.049.
     * 
     * eqn (2):
     * given q = 4X1 column vector of [eps eta] where eta is the scalar,
     * and "1" is a 3X3 identity matrix, also written as I_3.
     * and [eps]_x is the skew-symetric matrix for vector eps.
     * 
     * The left-hand compound operator is a 4X4 matrix:
     * 
     * q^+ = [ eta*I_3-[eps]_x   eps ]  // | [3X3]  [3X1] |
     *       [ -eps^T           eta ]   // | [1X3]  [1X1] |
     * 
     * multiplication of quaternions, u and v, which is typically written 
     * as u ⊗ v may be written equivalently as either
         u^+ v or v^⨁ u,
     * </pre>
     * @param quaternion 4X1 column vector of [eps eta] where eta is the scalar
     * @return 
     */
    public static double[][] quaternionLefthandCompoundOperator(double[] quaternion) {
        if (quaternion.length != 4) {
            throw new IllegalArgumentException("quaternion must be length 4");
        }
        double[] eps = Arrays.copyOfRange(quaternion, 0, 3);
        assert(eps.length == 3);
        double eta = quaternion[3];
        
        double[][] minusSkewEps = MatrixUtil.skewSymmetric(eps);
        MatrixUtil.multiply(minusSkewEps, -1);
        
        double[][] i3Eta = MatrixUtil.createIdentityMatrix(3);
        MatrixUtil.multiply(i3Eta, eta);
        
        // 3X3
        double[][] block00 = MatrixUtil.pointwiseAdd(i3Eta, minusSkewEps);
        
        double[][] lhc = MatrixUtil.zeros(4, 4);
        int i;
        for (i = 0; i < 3; ++i) {
            System.arraycopy(block00[i], 0, lhc[i], 0, 3);
            lhc[i][3] = eps[i];
        }
        for (i = 0; i < 3; ++i) {
            lhc[3][i] = -eps[i];
        }
        lhc[3][3] = eta;
        
        return lhc;
    }
    
    /**
     * given a quaternion, form the right-hand compound operator:
     * (symbol is superscript ⨁)
     * <pre>
     * from Barfoot, Forbes, & Furgale 2010, "Pose estimation using linearized 
     * rotations and quaternion algebra", Acta Astronautica (2010), doi:10.1016/j.actaastro.2010.06.049.
     * 
     * eqn (2):
     * given q = 4X1 column vector of [eps eta] where eta is the scalar,
     * and "1" is a 3X3 identity matrix, also written as I_3.
     * and [q]_x is the skew-symetric matrix for q.
     * 
     * The right-hand compound operator is a 4X4 matrix:
     * 
     * q^⨁ = [ eta*I_3+[q]_x   eps ]
     *        [ -eps^T          eta ]
     * 
     * multiplication of quaternions, u and v, which is typically written 
     * as u ⊗ v may be written equivalently as either
         u^+ v or v^⨁ u,
     * </pre>
     * @param quaternion 4X1 column vector of [eps eta] where eta is the scalar
     * @return [4X4]
     */
    public static double[][] quaternionRighthandCompoundOperator(double[] quaternion) {
        if (quaternion.length != 4) {
            throw new IllegalArgumentException("quaternion must be length 4");
        }
        double[] eps = Arrays.copyOfRange(quaternion, 0, 3);
        assert(eps.length == 3);
        double eta = quaternion[3];
        
        double[][] skewEps = MatrixUtil.skewSymmetric(eps);
        
        double[][] i3Eta = MatrixUtil.createIdentityMatrix(3);
        MatrixUtil.multiply(i3Eta, eta);
        
        // 3X3
        double[][] block00 = MatrixUtil.pointwiseAdd(i3Eta, skewEps);
        
        double[][] lhc = MatrixUtil.zeros(4, 4);
        int i;
        for (i = 0; i < 3; ++i) {
            System.arraycopy(block00[i], 0, lhc[i], 0, 3);
            lhc[i][3] = eps[i];
        }
        for (i = 0; i < 3; ++i) {
            lhc[3][i] = -eps[i];
        }
        lhc[3][3] = eta;
        
        return lhc;
    }
    
    /**
     * given a quaternion, return the inverse
     * (symbol is superscript -1).
     * <pre>
     * from Barfoot, Forbes, & Furgale 2010, "Pose estimation using linearized 
     * rotations and quaternion algebra", Acta Astronautica (2010), doi:10.1016/j.actaastro.2010.06.049.
     * 
     * eqn (5):
     * given q = 4X1 column vector of [eps eta] where eta is the scalar,
     * 
     * the inverse is the column vector:
     * 
     * q^-1 = [ -eps  eta ]
     * 
     * </pre>
     * <em>NOTE that if the quaternion is a unit-length quaternion, the conjugate
     * operation is usually called the inverse operation.</em>
     * @param quaternion 4X1 column vector of [eps eta] where eta is the scalar
     * @return 
     */
    public static double[] quaternionConjugateOperator(double[] quaternion) {
        if (quaternion.length != 4) {
            throw new IllegalArgumentException("quaternion must be length 4");
        }
        double[] inv = Arrays.copyOf(quaternion, 4);
        
        int i;
        for (i = 0; i < 3; ++i) {
            inv[i] *= -1;
        }
        
        return inv;
    }

    /**
     * rotate a point x to a new frame
     * <pre>
     * from Barfoot, Forbes, & Furgale 2010, "Pose estimation using linearized
     * rotations and quaternion algebra", Acta Astronautica (2010), doi:10.1016/j.actaastro.2010.06.049.
     *
     * eqn (11):
     *
     * v3 = [ p ]
     *      [ 0 ]
     *
     *  rotated = q^+ * q^(-1)^⨁ * v3 = R*v3
     *
     *  where R is the canonical 3X3 rotation matrix.
      *
     * @param q
     * @param x
     * @return
     */
    public static double[] rotateAPointByQuaternionBarfoot(double[] q, double[] x) {
        if (x.length != 3) {
            throw new IllegalArgumentException("x length must be 3");
        }
        if (q.length != 3) {
            throw new IllegalArgumentException("q length must be 4");
        }
        // (eqn 11) of Barfoot
        double[][] r = createRotation4FromQuaternion(q);
        double[] v = createQuaternionBarfootFromAPoint(x);
        return MatrixUtil.multiplyMatrixByColumnVector(r, v);
    }

    /**
     * multiply quaternions using active transformations
     * @param q1
     * @param q2
     * @return
     */
    public static double[] multiplyQuaternionsBarfoot(double[] q1, double[] q2) {
        if (q1.length != 3) {
            throw new IllegalArgumentException("q1 length must be 4");
        }
        if (q2.length != 3) {
            throw new IllegalArgumentException("q2 length must be 4");
        }
        double[][] lh1 = quaternionLefthandCompoundOperator(q1);
        return MatrixUtil.multiplyMatrixByColumnVector(lh1, q2);
    }

    public static double[] createQuaternionBarfootFromAPoint(double[] x) {
        if (x.length != 3) {
            throw new IllegalArgumentException("x length must be 3");
        }
        return Arrays.copyOf(x, 4);
    }

    public static double[] inverseQuaternionBarfoot(double[] q) {
        // -1*vector portion
        double[] inv = new double[4];
        for (int i =0; i < 3; ++i) {
            inv[i] = -q[i];
        }
        inv[3] = q[3];
        return q;
    }
    
    /**
     * return the identity quaternion
     * <pre>
     * from Barfoot, Forbes, & Furgale 2010, "Pose estimation using linearized 
     * rotations and quaternion algebra", Acta Astronautica (2010), doi:10.1016/j.actaastro.2010.06.049,
     * near eqn (8): 
     * 
     * creates identity element: [0 0 0 1]
     * 
     * </pre>
     * @return 4X1 column vector of [eps eta] where eta is the scalar, specifically
     * [0, 0, 0, 1]
     */
    public static double[] createIdentityQuaternion() {
        
        double[] i4 = new double[4];
        i4[3] = 1;
        
        return i4;
    }
    
     /**
     * create a [4x4] rotation matrix from the given quaternion.
     <pre>
       Note, the first [3X3] block is transposed compared to results 
           from Rotation.createRotationZYX(eulerAngles)
       where 
       double[] qHamilton = Rotation.createQuaternionHamiltonFromAngleAxis(eulerAngles);
       double[] quaternion = Rotation.createQuaternionBarfootFromHamilton(qHamilton);

      Note that the Hamilton systems uses passive (right-hand) transformations and the Barfoot system uses active
      (left-hand) transformations, so this method uses active transformations.
      
     This method's results are consistent with Shuster 1993.
     * The method follows Barfoot, Forbes, & Furgale 2010, "Pose estimation using linearized 
     * rotations and quaternion algebra", Acta Astronautica (2010), doi:10.1016/j.actaastro.2010.06.049.
     *
     * eqn (12):
     * 
     *  R = q^+ * q^(-1)^⨁ = q^(-1)^⨁ * q^+ = q^(⨁)^T * q^+
     *    
     *    = [ C   0 ]
     *      [ 0^T 1 ]
     * where C is the canonical 3X3 rotation matrix.
     * </pre>
     * @param quaternion rotation as a [4X1] column vector of [eps eta] where eta is the scalar
      *                   (i.e. Barfoot format)
     * @return a 4x4 rotation matrix whose 3X3 block at [0:2, 0:2] is the rotation matrix.
     */
    public static double[][] createRotation4FromQuaternion(double[] quaternion) {
        if (quaternion.length != 4) {
            throw new IllegalArgumentException("quaternion must be length 4");
        }
        double[][] lh = quaternionLefthandCompoundOperator(quaternion);
        double[][] invRH = quaternionRighthandCompoundOperator(quaternionConjugateOperator(quaternion));
        return MatrixUtil.multiply(lh, invRH);
    }
    
     /**
     * rotate point p by the given quaternion.
     * <pre>
     * from Barfoot, Forbes, & Furgale 2010, "Pose estimation using linearized 
     * rotations and quaternion algebra", Acta Astronautica (2010), doi:10.1016/j.actaastro.2010.06.049.
     * 
     * eqn (11):
     * 
     * v3 = [ p ]
     *      [ 0 ]
     * 
     *  rotated = q^+ * q^(-1)^⨁ * v3 = R*v3
     * 
     * </pre>
     * @param quaternion 4X1 column vector of [eps eta] where eta is the scalar
     * @param p the vector as a 3 element array followed by a 0, = [4X1].
     * @return 
     */
    public static double[] rotateVectorByQuaternion4(double[] quaternion, double[] p) {
        if (quaternion.length != 4) {
            throw new IllegalArgumentException("quaternion must be length 4");
        }
        if (p.length != 4) {
            throw new IllegalArgumentException("p must be length 4");
        }
        double[][] r = createRotation4FromQuaternion(quaternion);

        double[] result = MatrixUtil.multiplyMatrixByColumnVector(r, p);
        
        return result;
    }

    /**
     * calculate S_theta which is the matrix relating angular velocity to 
     * rotation angle rates.
     *
     The method uses intrinsic, active transformations.
     *
     * TODO: consider overloading for more rotation sequences.
     * <pre>
     * from Barfoot, Forbes, & Furgale 2010, "Pose estimation using linearized 
     * rotations and quaternion algebra", Acta Astronautica (2010), doi:10.1016/j.actaastro.2010.06.049.
     * 
     * eqn (21):
     * calc C = 3X3 rotation matrix (often written as R)
     * given array of euler rotation angles alpha, beta, gamma
     * 
     * s_theta column 0 = C_gamma(eulerAngles[2]) * C_beta(eulerAngles[1]) * [1, 0, 0]^T
     *         column 1 = C_gamma(eulerAngles[2]) * [0, 1, 0]^T
     *         column 2 = [0, 0, 1]^T
     * 
             C_gamma                        C_beta                            C_alpha
      cc rotation about z-axis (yaw):   cc about the y-axis (pitch):    cc about x-axis (roll):    
            | cos φ   -sin φ    0 |          |  cos ψ    0  sin ψ |         |    1       0       0 |  
            | sin φ    cos φ    0 |          |      0    1      0 |         |    0   cos θ  -sin θ |  
            |     0        0    1 |          | -sin ψ    0  cos ψ |         |    0   sin θ   cos θ |  
     *  eulerAngles[0] = angleX angle of rotation about x-axis (roll) in units of radians.
     *             can use createRotationRoll(eulerAngles[0])
     *  eulerAngles[1] = angleY angle of rotation about y-axis (pitch) in units of radians.
     *             can use createRotationPitch(eulerAngles[1])
     *  eulerAngles[2] = angleZ angle of rotation about z-axis (yaw) in units of radians.
     *             can use createRotationYaw(eulerAngles[2])
     * 
     * see also, pp 479-480, eqn (285) of Shuster 1993, "A Survey of AttitudeRepresentations"
     * http://www.ladispe.polito.it/corsi/Meccatronica/02JHCOR/2011-12/Slides/Shuster_Pub_1993h_J_Repsurv_scan.pdf
     * though the sign conventions of the sine terms are different
     * 
     * </pre>
     * This method uses active (left-hand system) transformations.
     *
     * @param eulerAngles euler X,Y,Z angles as representation of rotation matrix
     * @param seq Euler sequence in use
     * @return [3X3]
     */
    public static double[][] sTheta(double[] eulerAngles, EulerSequence seq) {
        if (eulerAngles.length != 3) {
            throw new IllegalArgumentException("eulerAngles length must be 3");
        }
        double[][] sTheta = MatrixUtil.zeros(3, 3);
        sTheta(eulerAngles, sTheta, seq);
        
        return sTheta;
    }
    
    /**
     * calculate S_theta which is the matrix relating angular velocity to 
     * rotation angle rates.
     * The method uses intrinsic, active transformations.
     *
     * TODO: consider overloading for more rotation sequences.
     *
     * <pre>
     * from Barfoot, Forbes, & Furgale 2010, "Pose estimation using linearized 
     * rotations and quaternion algebra", Acta Astronautica (2010), doi:10.1016/j.actaastro.2010.06.049.
       eqn (21)

     the book "State Estimation for Robotics", 2nd edition, by Barfoot
     page 271, eqn (7.51), (7.52):
     * 
     * </pre>
     *
     * @param eulerAngles euler angles as representation of rotation matrix
     * @param seq euler sequence in use
     * @param output 
     */
    public static void sTheta(double[] eulerAngles, double[][] output, EulerSequence seq) {
        if (eulerAngles.length != 3) {
            throw new IllegalArgumentException("eulerAngles length must be 3");
        }

        /*
        the book "State Estimation for Robotics", 2nd edition, by Barfoot
        page 271, eqn 7.52 for 1-2-3 sequence (active)  as naming for active euler sequence 3-2-1
        S(theta2, theta3)
            [ cos(th2) * cos(th3)   sin(th3)   0 ]
            [ -cos(th2)*sin(th3)    cos(th3)   0 ]
            [ sin(th2)              0          1 ]

        which is [yaw_passive(-th3) * pitch_passive(-th2) * I_col0   yaw_passive(-th3) * I_col2  I_col3 ]

        passive = A(a)*B(b)*C(c),
        active = (A(a)*B(b)*C(c))^-1
               = (C(c)^-1) * (B(b)^-1) * (A(a)^-1)
               = (C(c)^T) * (B(b)^T) * (A(a)^T)
        and transpose of a rotation matrix == untransposed with -1*angle
        so  active = C(-c)*B(-b)*A(-a)

        Barfoot "1-2-3" refers to passive Euler rotation sequence, then converts it to active
        by reversing the order and making the angles negative.

        furthermore in eqn (21) Barfoot et al., the columns for the identity matrix are in the passive order "1-2-3"
        This is eqn (7.51) in Barfoot book.

         */

        /*
        column 0:
         C_gamma(eulerAngles[2]) * C_beta(eulerAngles[1]) * [1, 0, 0]^T
        
        = | cos φ   -sin φ    0 |    *  |  cos ψ    0  sin ψ | * [1, 0, 0]^T
          | sin φ    cos φ    0 |    *  |      0    1      0 |
          |     0        0    1 |    *  | -sin ψ    0  cos ψ |
        = |  cos t2 * cos t1    -sin t2   cos t2 * sin t1 | * [1]
          |  sin t2 * cos t1     cos t2   sin t2 * sin t1 |   [0]
          | -sin t1              0        cos t1          |   [0]
        = | cos t2 * cos t1 |
          | sin t2 * cos t1 |
          | -sin t1         |

        column 1:
        | cos φ   -sin φ    0 |  * [0] = [-sin t2]
        | sin φ    cos φ    0 |    [1]   [cos t2 ]
        |     0        0    1 |    [0]   [0 ]

        column 3: [0 0 1]^T
        */
        
        double[] i0 = new double[]{1, 0, 0};
        double[] i1 = new double[]{0, 1, 0};
        double[] i2 = new double[]{0, 0, 1};

        // Barfoot uses euler active transformations.  we use the positive angles given, but later, if extracting
        //   angles, they should be considered negative in a passive context.

        double[][] cBeta = Rotation.createRotationPitch(eulerAngles[1]);

        double[] col0, col1, col2;
        if (seq.equals(EulerSequence.XYZ_ACTIVE)) {
            // XYZ active = passive w/ negative angles
            // eqn (58) of Barfoot et al. for a 3-2-1 sequence, but he means 1-2-3 w.r.t. euler notation...
            double[][] cAlpha = Rotation.createRotationRoll(eulerAngles[0]);
            col0 = MatrixUtil.multiplyMatrixByColumnVector(MatrixUtil.multiply(cAlpha, cBeta), i2);
            col1 = MatrixUtil.multiplyMatrixByColumnVector(cAlpha, i1);
            col2 = i0;
        } else  {
            // ZYX active = passive w/ negative angles = XYZ passive
            // eqn (21) of Barfoot et al. for a 1-2-3 sequence
            double[][] cGamma = Rotation.createRotationYaw(eulerAngles[2]);
            col0 = MatrixUtil.multiplyMatrixByColumnVector(MatrixUtil.multiply(cGamma, cBeta), i0);
            col1 = MatrixUtil.multiplyMatrixByColumnVector(cGamma, i1);
            col2 = i2;
        }
        
        int row;
        for (row = 0; row < 3; ++row) {
            output[row][0] = col0[row];
            output[row][1] = col1[row];
            output[row][2] = col2[row];
        }
    }

    /**
     * calculate inv S_theta which is the matrix relating angular velocity to
     * rotation angle rates.
     * This method fails for eulerAngles[1]=0, so consider using the singularity safe methods instead.
     * The method uses intrinsic, active transformations.
     *
     * <pre>
     the book "State Estimation for Robotics", 2nd edition, by Barfoot
     eqn (7.53):
     *
     * </pre>
     *
     * @param eulerAngles euler angles as representation of rotation matrix
     * @result inverse of s(eulerAngles)
     */
    private static double[][] invSTheta(double[] eulerAngles) {
        if (eulerAngles.length != 3) {
            throw new IllegalArgumentException("eulerAngles length must be 3");
        }
        if (eulerAngles[1] == 0) {
            throw new IllegalArgumentException("eulerAngles[1] cannot be 0");
        }

        double sec2 = 1./Math.cos(eulerAngles[1]);
        double cos3 = Math.cos(eulerAngles[2]);
        double sin3 = Math.sin(eulerAngles[2]);
        double tan2 = Math.tan(eulerAngles[1]);

        double[][] sInv = new double[][] {
                {sec2 * cos3, -sec2*sin3, 0},
                {sin3, cos3, 0},
                {-tan2*cos3, tan2*sin3, 1}
        };

        return sInv;
    }
    
    /**
     * calculate  R_zyx = (R_z(theta_x)*R_y(theta_y))*R_x(theta_z)
     * where theta_x = eulerZYX[0], theta_y=eulerZYX[1], theta_z=eulerZYX[2].
     * This method uses a passive (right-hand) rotation system with CCW (=CC) rotations.
     * This method uses "intrinsic" composition.
     *
     * @return R_zyx intrinsic
     */
    public static double[][] createRotationZYX(double eulerX, double eulerY, double eulerZ) {

        double[][] rX = Rotation.createRotationRoll(eulerX);
        double[][] rY = Rotation.createRotationPitch(eulerY);
        double[][] rZ = Rotation.createRotationYaw(eulerZ);

        return MatrixUtil.multiply(MatrixUtil.multiply(rZ, rY), rX);
    }

    /**
     * calculate  R_zyx = (R_z(theta_x)*R_y(theta_y))*R_x(theta_z)
     * where theta_x = thetas[0], theta_y=thetas[1], theta_z=thetas[2].
     * This method uses a passive (right-hand) rotation system with CCW (=CC) rotations.
     * This method uses "intrinsic" composition.
     *
     * @param aa auxiliary arrays to use as space for intermediate calculations
     * @out gets populated with R_zyx intrinsic
     */
    public static void createRotationZYX(double eulerX, double eulerY, double eulerZ,
                                         AuxiliaryArrays aa, double[][] out) {
        if (out.length != 3 || out[0].length != 3) {
            throw new IllegalArgumentException("out must be 3 x 3");
        }

        //ZYX = transpose of xyz(-1*thetas)
        
        double[][] rZ = aa.a3X3;
        Rotation.createRotationYaw(eulerZ, rZ);
        
        double[][] rY = aa.b3X3;
        Rotation.createRotationPitch(eulerY, rY);
        
        double[][] rX = aa.c3X3;
        Rotation.createRotationRoll(eulerX, rX);

        MatrixUtil.multiply(MatrixUtil.multiply(rZ, rY), rX, out);
    }
    
    /**
     * calculate  the Hamilton unit-length quaternion which would be extracted from the
     * rotation matrix created by R_xyz = R_z(theta_z)*R_y(theta_y)*R_z(theta_z).
     *
     * This method uses a passive (right-hand) rotation system with CCW (=CC) rotations.
     *
     * from https://en.wikipedia.org/wiki/Davenport_chained_rotations
     * "Any extrinsic rotation is equivalent to an intrinsic rotation by the 
     * same angles but with inverted order of elemental rotations, and vice 
     * versa. For instance, the intrinsic rotations x-y’-z″ by angles α, β, γ 
     * are equivalent to the extrinsic rotations z-y-x by angles γ, β, α."
     * <pre>
     * from Barfoot, Forbes, & Furgale 2010, "Pose estimation using linearized 
     * rotations and quaternion algebra", Acta Astronautica (2010), doi:10.1016/j.actaastro.2010.06.049.
     * 
     * eqn (18)
     * 
     * where  
       cc rotation about z-axis (yaw):   cc about the y-axis (pitch):    cc about x-axis (roll):    
            | cos φ   -sin φ    0 |          |  cos ψ    0  sin ψ |         |    1       0       0 |  
            | sin φ    cos φ    0 |          |      0    1      0 |         |    0   cos θ  -sin θ |  
            |     0        0    1 |          | -sin ψ    0  cos ψ |         |    0   sin θ   cos θ |  
            * 

        = | (cos φ * cos ψ)   (-sin φ * cos θ + cos φ * sin ψ * sin θ)   (sin φ * sin θ + cos φ * sin ψ * cos θ)   |
          | (sin φ * cos ψ)   ( cos φ * cos θ + sin φ * sin ψ * sin θ)   (-cos φ * sin θ + sin φ * sin ψ * cos θ)  |
          | (-sin ψ)          ( cos ψ * sin θ )                          (cos ψ * cos θ)                           |
          
       =  | (cosZ * cosY)   (-sinZ * cosX + cosZ * sinY * sinX)   (sinZ * sinX + cosZ * sinY * cosX)   |
          | (sinZ * cosY)   ( cosZ * cosX + sinZ * sinY * sinX)   (-cosZ * sinX + sinZ * sinY * cosX)  |
          | (-sinY)          ( cosY * sinX )                       (cosY * cosX)                       |
           
     given angle theta and axis n:
         q = [scalar, vector] = 
             [cos(theta/2), nx*sin(theta/2), ny*sin(theta/2), nz*sin(theta/2)];
          
     given euler angles:
         q1 =  sin(z/2)*sin(y/2)*sin(x/2) + cos(z/2)*cos(y/2)*cos(x/2)
         q2 = -sin(z/2)*sin(y/2)*cos(x/2) + sin(x/2)*cos(z/2)*cos(y/2)
         q3 =  sin(z/2)*sin(x/2)*cos(y/2) + sin(y/2)*cos(z/2)*cos(x/2)
         q4 =  sin(z/2)*cos(y/2)*cos(x/2) - sin(y/2)*sin(x/2)*cos(z/2)

     This is what is used in eqn (2.40) of "An Invitation to 3-D Vision",
     Ma, Soatto, Kosecka, and Sastry (MASKS).

     It is an intrinsic 'XYZ' format.
     * </pre>
     * @param angle
     * @param axis
     * @return the quaternion in format scalar as first term.
     */
    public static double[] createQuaternionHamiltonFromAngleAxis(double angle, double[] axis) {
        if (axis.length != 3) {
            throw new IllegalArgumentException("axis must be length 3");
        }
        double d = MatrixUtil.lPSum(axis, 2);
        double nx = axis[0]/d;
        double ny = axis[1]/d;
        double nz = axis[2]/d;
        double ca = Math.cos(angle/2);
        double sa = Math.sin(angle/2);

        double[] q = new double[]{ca, nx*sa, ny*sa, nz*sa};
        
        return q;
    }
    
    /**
     * convert the quaternion from [scalar vector] to [vector scalar].
     * Barfoot uses intrinsic, active transformations with scalar last format.
     * Hamilton uses intrinsic with scalar first format.
     * @param qHamilton quaternion of format [scalar  vector]
     * @return convert the quaternion from [scalar vector] to [vector scalar]
     */
    public static double[] createQuaternionBarfootFromHamilton(double[] qHamilton) {
        if (qHamilton.length != 4) {
            throw new IllegalArgumentException("qHamilton length must be 4");
        }
        double[] out = Arrays.copyOf(qHamilton, 4);
        convertQuaternionHamiltonToBarfoot(out);
        return out;
    }

    /**
     * convert the quaternion from [scalar vector] to [vector scalar].
     * Barfoot uses intrinsic, active transformations with scalar last format.
     * Hamilton uses intrinsic with scalar first format.
     * @param qHamilton quaternion in Hamilton format that is converted to Barfoot format in place (modifying array)
     */
    public static void convertQuaternionHamiltonToBarfoot(double[] qHamilton) {
        double scalar = qHamilton[0];
        for (int i = 1; i < qHamilton.length; ++i) {
            qHamilton[i-1] = qHamilton[i];
        }
        qHamilton[3] = scalar;
    }
    /**
     * convert the quaternion from [scalar vector] to [vector scalar].
     * Barfoot uses intrinsic, active transformations with scalar last format.
     * Hamilton uses intrinsic with scalar first format.
     * @param qBarfoot quaternion in Barfoot format to be converted to Hamilton format in place (modifies array)
     */
    public static void convertQuaternionBarfootToHamilton(double[] qBarfoot) {
        double scalar = qBarfoot[3];
        for (int i = 2; i >= 0; --i) {
            qBarfoot[i+1] = qBarfoot[i];
        }
        qBarfoot[0] = scalar;
    }

    /**
     * convert the quaternion from [scalar vector] to [vector scalar].
     * Barfoot uses intrinsic, active transformations with scalar last format.
     * Hamilton uses intrinsic with scalar first format.
     * @param qBarfoot quaternion of format [scalar  vector]
     * @return convert the quaternion from [scalar vector] to [vector scalar]
     */
    public static double[] createQuaternionHamiltonFromBarfoot(double[] qBarfoot) {
        if (qBarfoot.length != 4) {
            throw new IllegalArgumentException("qHamilton length must be 4");
        }
        double[] out = Arrays.copyOf(qBarfoot, 4);
        convertQuaternionBarfootToHamilton(out);
        return out;
    }

    /*
    create a quaternion 4X1 column vector of [vector scalar] from the given
    angle and axis representation of rotation.
    Barfoot et al. place the scalar as the last item in the quaternion.
    Barfoot uses intrinsic, active transformations.
    <pre>
     * from Barfoot, Forbes, & Furgale 2010, "Pose estimation using linearized 
     * rotations and quaternion algebra", Acta Astronautica (2010), doi:10.1016/j.actaastro.2010.06.049.
     * 
     * eqn (32)
    
    satisfies unit-length constraint q^T*q = 1  (size [1X4]*[4X1]=[1X1].
    q^T*q = q_1^2 + q_2^2 + q_3^2 + q_4^2.
    </pre>
    @param unitLengthAxis axis of rotation normalized to unit length
    @param angle rotation about axis in radians
    */
    public static double[] createQuaternionUnitLengthBarfoot(double angle, double[] unitLengthAxis) {
        if (unitLengthAxis.length != 3) {
            throw new IllegalArgumentException("unitLengthAxis.length must be 3");
        }

        double sumSq = MatrixUtil.lPSum(unitLengthAxis, 2);
        sumSq *= sumSq;

        //TODO: revisit tolerance here:
        if (Math.abs(sumSq - 1.) > 1E-5) {
            // can provide same normalization for all, by using rotation vector as an intermediary
            double[] rotVec = Rotation.createRotationVectorFromAngleAxis(unitLengthAxis, angle);
            double[] _axis = new double[3];
            angle = Rotation.createAngleAxisFromRotationVector(rotVec, _axis);
            unitLengthAxis = _axis;
            sumSq = MatrixUtil.lPSum(unitLengthAxis, 2);
            sumSq *= sumSq;
        }

        // follows from sum of squares of each element in unitLengthAxis should be == 1
        //cos(theta)^2 + sin(theta)^2 = 1
        // scalar term = cos(theta)^2
        // vector portion = |sin(theta)|^2
        /// u = q // ||q|| vector portion
        // q2 = q0 + q = cos(theta) + u*sin(theta) ...

        double cPhi2 = Math.cos(angle/2.);
        double sPhi2 = Math.sin(angle/2.);
        double[] out = Arrays.copyOf(unitLengthAxis, 4);
        for (int i = 0; i < 3; ++i) {
            out[i] *= sPhi2;
        }
        out[3] = cPhi2;
        return out;
    }

    /*
    create a quaternion 4X1 column vector of [vector scalar] from the given
    angle and axis representation of rotation.
    Barfoot et al. place the scalar as the last item in the quaternion.
    Barfoot uses intrinsic, active transformations.
    <pre>
     * from Barfoot, Forbes, & Furgale 2010, "Pose estimation using linearized
     * rotations and quaternion algebra", Acta Astronautica (2010), doi:10.1016/j.actaastro.2010.06.049.
     *
     * eqn (32) w/ conversion from eulerAngles to axis angle representation

    satisfies unit-length constraint q^T*q = 1  (size [1X4]*[4X1]=[1X1].
    q^T*q = q_1^2 + q_2^2 + q_3^2 + q_4^2.
    </pre>
    @param eulerAngles eulerAngles
    */
    public static double[] createQuaternionUnitLengthBarfootFromEulerXYZ(double[] eulerAngles) {
        if (eulerAngles.length != 3) {
            throw new IllegalArgumentException("unitLengthAxis.length must be 3");
        }
        double[] rotVec = createRotationVectorFromEulerAnglesXYZ(eulerAngles);
        double[] axis = new double[3];
        double angle = createAngleAxisFromRotationVector(rotVec, axis);
        return createQuaternionUnitLengthBarfoot(angle, axis);
    }

    public static double[] createQuaternionUnitLengthBarfootFromEuler(double[] eulerAnglesXYZ, EulerSequence seq) {
        if (eulerAnglesXYZ.length != 3) {
            throw new IllegalArgumentException("unitLengthAxis.length must be 3");
        }
        double[] rotVec = createRotationVectorFromEulerAngles(eulerAnglesXYZ, seq);
        double[] axis = new double[3];
        double angle = createAngleAxisFromRotationVector(rotVec, axis);
        return createQuaternionUnitLengthBarfoot(angle, axis);
    }

    public static double[] createQuaternionUnitLengthBarfootFromRotationVector(double[] rotVec) {
        if (rotVec.length != 3) {
            throw new IllegalArgumentException("rotVec.length must be 3");
        }
        double[] axis = new double[3];
        double angle = createAngleAxisFromRotationVector(rotVec, axis);
        return createQuaternionUnitLengthBarfoot(angle, axis);
    }

    /*
    <pre>
     * from Barfoot, Forbes, & Furgale 2010, "Pose estimation using linearized 
     * rotations and quaternion algebra", Acta Astronautica (2010), doi:10.1016/j.actaastro.2010.06.049.
     * 
     * eqn (35)
    
    </pre>
    @param principalAxis the axis of rotation, that is x, y, or z.
    (rotation about x is roll, y is pitch, and z is yaw).
    */
    public static double[] quaternionPrincipalAxisRotation(double angle, int principalAxis) {
        if (principalAxis < 0 || principalAxis > 2) {
            throw new IllegalArgumentException("principalAxis must be 0, 1, or 2");
        }
        double sa = Math.sin(angle/2.);
        double ca = Math.cos(angle/2.);
        double[] q = new double[4];
        // two 0's and a 1 for principal axis:
        q[principalAxis] = sa;
        q[3] = ca;
        return q;
    }
    
    /**
     * create a rotation quaternion from Euler x, y, z angles.
     * note that the rotation order of multiplication operations is z, y, x.
     * 
     * <pre>
     * from Barfoot, Forbes, & Furgale 2010, "Pose estimation using linearized 
     * rotations and quaternion algebra", Acta Astronautica (2010), doi:10.1016/j.actaastro.2010.06.049.
     * 
     * eqn (37)
     *     q(euler) = q_z(euler_z)^+ * q_y(euler_y)^+ + q_x(euler_x) 
     * 
     *     and principal axis rotations in eqn (35)
     * </pre>
     * @param rotVec
     * @return [4X1] quaternion rotation vector.
     */
    public static double[] createQuaternionBarfootFromRotationVector(double[] rotVec,
                                                                     EulerSequence seq) {
        if (rotVec.length != 3) {
            throw new IllegalArgumentException("rotVec.length must be 3");
        }
        // length is 4
        double[] q0 = quaternionPrincipalAxisRotation(-rotVec[0], 0);
        double[] q1 = quaternionPrincipalAxisRotation(-rotVec[1], 1);
        double[] q2 = quaternionPrincipalAxisRotation(-rotVec[2], 2);

        /*
        Barfoot "1-2-3" and "alpha-beta-gamma" as names, refer to the passive Euler rotation sequence XYZ,
        then Barfoot impl. converts the sequence to active by reversing the order and making the angles negative.
        passive A(a)*B(b)*C(c)
        active A(-a)*B(-b)*C(-c)
            = A(a)^T * B(b)^T * C(c)^T
            = (C(c) * B(b) * C(c))^T
         */
        double[][] q1m, q2m;
        double[] q3v;
        if (seq.equals(EulerSequence.XYZ_ACTIVE)) {
            // active XYZ == passive XYZ w/ negative angles
            // [4X4]
            //TODO: these might need to be swapped
            q1m = quaternionRighthandCompoundOperator(q2);
            q2m = quaternionRighthandCompoundOperator(q1);
            q3v = q0;
        } else {
            // active ZYX == passive ZYX w/ negative angles
            // [4X4]
            q1m = quaternionRighthandCompoundOperator(q0);
            q2m = quaternionRighthandCompoundOperator(q1);
            q3v = q2;
        }

        // q2m times q1m times q0 = [4X1]
        // [4X4]
        double[][] t1 = MatrixUtil.multiply(q2m, q1m);
        // 4X1
        double[] q = MatrixUtil.multiplyMatrixByColumnVector(t1, q3v);
        
        return q;
    }
    
    /**
     * create the rotation vector dPhi = S(eulerAngles) * dTheta
     * <pre>
     * from Barfoot, Forbes, & Furgale 2010, "Pose estimation using linearized 
     * rotations and quaternion algebra", Acta Astronautica (2010), doi:10.1016/j.actaastro.2010.06.049.
     * 
     * text under eqn (26)
     *    dPhi= S(eulerAngles) * dTheta
     * 
     * </pre>
     * @param eulerAngles euler rotation angles
     * @param dTheta perturbation to apply to rotation
     * @param seq Euler sequence in use
     * @return dPhi= S(eulerAngles) * dTheta.  length is 3.
     */
    public static double[] createRotationVectorBarfoot(double[] eulerAngles, double[] dTheta, EulerSequence seq) {
        
        double[][] sTheta = sTheta(eulerAngles, seq);
        
        // length 3
        double[] dPhi = MatrixUtil.multiplyMatrixByColumnVector(sTheta, dTheta);
        
        return dPhi;
    }

    /**
     * Apply a perturbation dTheta to the rotation matrix represented by euler angles theta0
     * and return the result as a quaternion in Barfoot format or as a rotation matrix.  The result also contains
     * dPhi, the rotation vector, which can be used as an updatable state structure for future
     * use that avoids the need to extract euler angles from a rotation quaternion or matrix.
     *
     * The scalar is the last term in the quaternion.  The rotation transformations are active and intrinsic.
     <pre>
     from Barfoot, Forbes, & Furgale 2010, "Pose estimation using linearized
     rotations and quaternion algebra", Acta Astronautica (2010), doi:10.1016/j.actaastro.2010.06.049.

     rotation vector:
         dPhi = S(theta0) * dTheta

     eqn (29c)
          C(thetaResult) ~ C(theta0 + dTheta) ~ (I_3 - [dPhi]_x) * C(theta0)

     eqn (49)
        q(thetaResult) = q(theta0 + perturbation) = q(dPhi)^⨁ * q(theta0)

     "This update approach allows us to store and update the rotation as a
     unit-length quaternion, thereby avoiding singularities and the need to
     restore the constraint afterwards (i.e., constraint restoration is
     built in)."

     </pre>
     * @param theta0 euler rotation angles used to define the rotation matrix that will
     *              be perturbed by dTheta.  The order of the angles in the array is X,Y,Z even
     *               if the seq used is not.
     *               In Barfoot paper eqn 26, this is theta with a bar over it.
     * @param dTheta the perturbation to apply to the rotation matrix as euler rotation angles as
     *               X, Y, Z ordered values.
     *
     * @param seq Euler sequence
     * @param returnQuaternion if true, calculates the rotation quaternion and returns a datastructure holding it,
     *                         else if false, calculates the rotation matrix and returns a datastructure holding it.
     * @return a data structure holding the resulting rotation matrix or quaternion, the EulerSequence, and the
     * updatable rotation vector dPhi.  This returned data structure can be used in the overloaded
     * applySingularitySafeRotationPerturbation.
     */
    public static RotationPerturbation applySingularitySafeRotationPerturbation(double[] theta0, double[] dTheta,
        EulerSequence seq, boolean returnQuaternion) {

        //TODO: add consistency checks like allowed range values

        if (theta0.length != 3) {
            throw new IllegalArgumentException("theta0 length must be 3");
        }
        if (dTheta.length != 3) {
            throw new IllegalArgumentException("theta0 must be length 3");
        }

        //let dPhi = S(theta0)*dTheta

        //length 3
        double[] dPhi = createRotationVectorBarfoot(theta0, dTheta, seq);

        if (returnQuaternion) {

            double[] _qDPhi = createQuaternionBarfootFromRotationVector(dPhi, seq);

            // eqn (48c), compare to _qDPhi
            double[] qDPhi = Arrays.copyOf(dPhi, 4);
            MatrixUtil.multiply(qDPhi, 0.5);
            qDPhi[3] = 1;

            double[][] qLH = quaternionLefthandCompoundOperator(qDPhi);

            double[] qTheta = createQuaternionUnitLengthBarfootFromEuler(theta0, seq);

            double[] q2 = MatrixUtil.multiplyMatrixByColumnVector(qLH, qTheta);

            RotationPerturbationQuaternion result = new RotationPerturbationQuaternion(dPhi, seq, q2);

            return result;
        }

        // infinitesimally small rotation matrix for the perturbation
        double[][] rPerturb = MatrixUtil.skewSymmetric(dPhi);
        for (int i = 0; i < rPerturb.length; ++i) {
            rPerturb[i][i] = 1. - rPerturb[i][i];
        }

        // for active transformations, we need -1*euler angles to use with euler rotation matrices that are intrinsic

        //from theta0 create r0
        double[][] r0;
        if (seq.equals(EulerSequence.ZYX_ACTIVE)) {
            r0 = createRotationZYX(-theta0[0], -theta0[1], -theta0[2]);
        } else {
            r0 = createRotationXYZ(-theta0[0], -theta0[1], -theta0[2]);
        }

        double[][] r2 = MatrixUtil.multiply(rPerturb, r0);

        RotationPerturbationMatrix result = new RotationPerturbationMatrix(dPhi, seq, r2);

        return result;
    }

    /**
     * Apply a perturbation dTheta to the given rotation matrix or quaternion and return the updated perturbed data
     * structure.
     *
     <pre>
     from Barfoot, Forbes, & Furgale 2010, "Pose estimation using linearized
     rotations and quaternion algebra", Acta Astronautica (2010), doi:10.1016/j.actaastro.2010.06.049.

     rotation vector:
         dPhi = S(theta0) * dTheta

     eqn (29c)
          C(thetaResult) ~ C(theta0 + dTheta) ~ (I_3 - [dPhi]_x) * C(theta0)

     eqn (49)
        q(thetaResult) = q(theta0 + perturbation) = q(dPhi)^⨁ * q(theta0)

     "This update approach allows us to store and update the rotation as a
     unit-length quaternion, thereby avoiding singularities and the need to
     restore the constraint afterwards (i.e., constraint restoration is
     built in)."

     </pre>
     * @param state data structure holding the current state rotation matrix or quaternion, the EulerSequence, and the
     * current state updatable rotation vector dPhi.
     * @param dTheta the perturbation to apply to the rotation matrix as euler rotation angles.
     *               The order of the angles in the array is X,Y,Z even
     *      *               if the seq used is not.
     * @return a data structure holding the resulting rotation matrix or quaternion, the EulerSequence, and the
     * updatable rotation vector dPhi.
     */
    public static RotationPerturbation applySingularitySafeRotationPerturbation(RotationPerturbation state, double[] dTheta) {
        if (dTheta.length != 3) {
            throw new IllegalArgumentException("theta0 must be length 3");
        }

        /*
        let dPhi = S(theta0)*dTheta
        C(theta)
                ~ C(theta0 + dTheta)
                ~ (I_3 - [dPhi]_x) * C(theta0)
         */

        //length 3
        double[] dPhi = state.dPhi;

        if (state instanceof RotationPerturbationQuaternion) {

            double[] _qDPhi = createQuaternionBarfootFromRotationVector(dPhi, state.seq);

            // eqn (48c), compare to _qDPhi
            double[] qDPhi = Arrays.copyOf(dPhi, 4);
            MatrixUtil.multiply(qDPhi, 0.5);
            qDPhi[3] = 1;

            double[][] qLH = quaternionLefthandCompoundOperator(qDPhi);

            double[] qTheta = ((RotationPerturbationQuaternion) state).quaternion;

            double[] q2 = MatrixUtil.multiplyMatrixByColumnVector(qLH, qTheta);

            RotationPerturbationQuaternion result = new RotationPerturbationQuaternion(dPhi, state.seq, q2);

            return result;
        }

        // infinitesimally small rotation matrix for the perturbation
        double[][] rPerturb = MatrixUtil.skewSymmetric(dPhi);
        for (int i = 0; i < rPerturb.length; ++i) {
            rPerturb[i][i] = 1. - rPerturb[i][i];
        }

        double[][] r0 = ((RotationPerturbationMatrix)state).rotation;

        double[][] r2 = MatrixUtil.multiply(rPerturb, r0);

        RotationPerturbationMatrix result = new RotationPerturbationMatrix(dPhi, state.seq, r2);

        return result;
    }

    /**
     * first term in RHS of eqn (30b) in Barfoot et al. paper.
     * @param theta0 euler angles representing current state.
     *  The order of the angles in the array is X,Y,Z even if the seq used is not.
     * @param dTheta euler angle perturbations.  The order of the angles in the array is X,Y,Z even if the seq used is not.
     * @param seq
     * @return
     */
    protected static double[][] createRotationInfSmallPerturbation(double[] theta0, double[] dTheta, EulerSequence seq) {

        //dPhi = S(theta) * dTheta
        double[] dPhi = createRotationVectorBarfoot(theta0, dTheta, seq);

        // infinitesimally small rotation matrix:
        double[][] qPhiX = MatrixUtil.skewSymmetric(dPhi);
        for (int i = 0; i < qPhiX.length; ++i) {
            qPhiX[i][i] = 1. - qPhiX[i][i];
        }
        return qPhiX;
    }

    /**
     * Apply the perturbation dTheta to rotation matrix theta and
     * return the perturbed rotation.
     * The rotation transformations are active and intrinsic.
     * <pre>
     * from Barfoot, Forbes, & Furgale 2010, "Pose estimation using linearized
     * rotations and quaternion algebra", Acta Astronautica (2010), doi:10.1016/j.actaastro.2010.06.049.
     *
     * eqn (26)
     *    C(theta2) ~ C(theta + dTheta) ~ (I_3 - [dPhi]_x) * C(theta)
     *    where C is for rotation matrix.
     *
     *    (I_3 - [dPhi]_x) is effectively a rotation matrix, but is infinitesimally small.
     *
     * "This update approach allows us to store and update the rotation as a
     * unit-length quaternion, thereby avoiding singularities and the need to
     * restore the constraint afterwards (i.e., constraint restoration is
     * built in)."
     *
     * </pre>
     * NOTE: this method is not efficient, but is kept for use in tests.
     *
     * @param rTheta the current rotation matrix.  note that it should have been formed using intrinsic active ZYX.
     * @param dTheta the perturbation to apply to the rotation matrix as euler rotation angles
     * must be small (cosine(dTheta[i]) ~ 1, sine(dTheta[i] ~ 0 or sine(dTheta[i])/i ~ 1).
     *               The order of the angles in the array is X,Y,Z even if the seq used is not.
     * @return resulting quaternion from perturbation applied to quaternion
     * formed from theta euler angles.
     */
    protected static double[][] _applySingularitySafeRotationPerturbation(double[][] rTheta, double[] dTheta,
                                                                      EulerSequence seq) {
        if (dTheta.length != 3) {
            throw new IllegalArgumentException("dTheta length must be 3");
        }
        if (rTheta.length != 3 || rTheta[0].length != 3) {
            throw new IllegalArgumentException("rTheta must be 3 X 3");
        }

        // eqn (26)
        //C(theta2) ~ C(theta + dTheta) ~ (I_3 - [dPhi]_x) * C(theta)

        double[] theta0;
        if (seq.equals(EulerSequence.ZYX_ACTIVE)) {
            theta0 = extractThetaFromZYX(rTheta, false);
        } else {
            theta0 = extractThetaFromXYZ(rTheta, false);
        }

        double[][] qPhiX = createRotationInfSmallPerturbation(theta0, dTheta, seq);

        return MatrixUtil.multiply(qPhiX, rTheta);
    }

    static double[] _extractRotationVectorFromQuaternionBarfoot(double[] q) {
        if (q.length != 4) {
            throw new IllegalArgumentException("q must be length 4");
        }
        // extracting from eqn (7.14) of Barfoot book
        // scalar term = cos(angle/2);  acos(q[3])= angle/2
        double[] axis = new double[3];
        double angle = extractAngleAxisFromBarfootQuaternion(q, axis);
        if (angle == 0) {
            return new double[]{0, 0, 0};
        }
        return createRotationVectorFromAngleAxis(axis, angle);
    }
    private static double extractAngleAxisFromBarfootQuaternion(double[] q, double[] outAxis) {
        if (q.length != 4) {
            throw new IllegalArgumentException("q must be length 4");
        }
        if (outAxis.length != 3) {
            throw new IllegalArgumentException("outAxis must be length 3");
        }
        // extracting from eqn (7.14) of Barfoot book
        // scalar term = cos(angle/2);  acos(q[3])= angle/2
        double angle;
        angle = 2 * Math.acos(q[3]);
        // wikipedia states it is more mathematically stable to use:
        //angle = 2 * Math.atan2(MatrixUtil.lPSum(Arrays.copyOf(q, 3), 2), q[3]);
        System.arraycopy(q, 0, outAxis, 0, 3);
        if (angle == 0) {
            return 0;
        }
        double s = Math.sin(angle/2);
        for (int i = 0; i < 3; ++i) {
            outAxis[i] /= s;
        }
        return angle;
    }
        
    public static class AuxiliaryArrays {
        final double[][] a3X3;
        final double[][] b3X3;
        final double[][] c3X3;
        final double[][] d3X3;
        public AuxiliaryArrays() {
            a3X3 = MatrixUtil.zeros(3, 3);
            b3X3 = MatrixUtil.zeros(3, 3);
            c3X3 = MatrixUtil.zeros(3, 3);
            d3X3 = MatrixUtil.zeros(3, 3);
        }
    }

    public static double[][] orthonormalizeUsingSVD(double[][] r) throws NotConvergedException {
        SVD svd = SVD.factorize(new DenseMatrix(r));
        DenseMatrix r2 = (DenseMatrix) svd.getU().mult(svd.getU(), svd.getVt());
        return MatrixUtil.convertToRowMajor(r2);
    }
    public static double[][] orthonormalizeUsingSkewCayley(double[][] r) throws NotConvergedException {
        r = cay(r);
        double[][] skew = MatrixUtil.pointwiseSubtract(
                r, MatrixUtil.transpose(r));
        MatrixUtil.multiply(skew, 0.5);
        r = cay(skew);
        return r;
    }

    /**
     * orthogonalize matrix R using skew parameters via Cayley's formula.
     * (I - A)*(I + A)^-1
     * @param r  a rotation matrix (i.e. a skew symmetric matrix A^T = A^-1)
     * @return
     * @throws NotConvergedException
     */
    public static double[][] cay(double[][] r) throws NotConvergedException {
        if (r.length != r[0].length) {
            throw new IllegalArgumentException("r must be a square rotation matrix");
        }
        double[][] identity = MatrixUtil.createIdentityMatrix(r.length);
        double[][] t1 = MatrixUtil.pseudoinverseFullColumnRank(
                MatrixUtil.pointwiseAdd(identity, r));
        double[][] t2 = MatrixUtil.pointwiseSubtract(identity, r);
        return MatrixUtil.multiply(t1, t2);
    }
}
