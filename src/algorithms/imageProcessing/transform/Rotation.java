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
 * This Rotation.java class uses "active" rotations of vectors counterclockwise in a right-handed 
 * coordinate system (y counterclockwise from x) by pre-multiplication 
 * (R on the left). If any one of these is changed (such as rotating axes 
 * instead of vectors, a "passive" transformation), then the inverse of the
 * example matrix should be used, which coincides with its transpose.
 * e.g. (A*B*C)^-1 = (C^-1) * (B^-1) * (A^-1).
 * 
 * from: "Tutorial: Consistent Representations of and Conversions Between 3D Rotations"
   by Rowenhorst, Rollett, Rohrer, Groeber, Jackson, Konijnenberg, and De Graef.
   Modelling Simulation Mater. Sci. Eng.
 * A rotation can be viewed as operating on the object, which is the 
 * active interpretation, or operating on the reference frame, which is the 
 * passive interpretation. 
 * An active rotation transforms object coordinates to new coordinates in the 
 * same reference frame; for the passive interpretation, the initial and final 
 * reference frames are different.
 * 
 * 
 * <b><ul>RIGHT HAND SYSTEM</ul></b>
 * The equations below use a right hand system.
 * The right hand system is consistent with methods in physics, engineering,
 * and computer science in general.  The <b>Hamilton quaternion</b> is consistent
 * with the right hand system and has format [scalar vector].
 * The NASA 1977 publication and Szeliski 2010
 *    define a quaternion as (qw, qx, qy, qz) where qw is a scalar and
 * [qx, qy, qz] is a vector.   That is the Hamilton Quaternion format.
 * 
 * In contrast, Barfoot et al. define a quaternion as (qx, qy, qz, qw).
 * One might need to transform some properties to quaternions and then modify 
 * the placement of the scalar term in order
 * to compare results with other "Hamilton" "Right hand" system results.
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
 
cc rotation about z-axis (yaw):   cc about the y-axis (pitch):    cc about x-axis (roll):    
            | cos φ   -sin φ    0 |          |  cos ψ    0  sin ψ |         |    1       0       0 |  
            | sin φ    cos φ    0 |          |      0    1      0 |         |    0   cos θ  -sin θ |  
            |     0        0    1 |          | -sin ψ    0  cos ψ |         |    0   sin θ   cos θ | 
 
two ways that the skew matrix are expressed are transposed from one another:
 
 [v]_x = |  0  -z   y |   this one is used in MatrixUtil.skewSymmetric()
         |  z   0  -x |
         | -y   x   0 |
         
 [v]_x = |  0   z  -y | 
         | -z   0   x |
         |  y  -x   0 |
 
 * Eigen, ROS, and Google Ceres use Hamilton convention.
*  Also  Wolfram Mathematica, Matlab’s aerospace(!) and robotics toolbox, 
*  Boost, GNU Octave, NASA’s SPICE.

Note that Shuster 1993 use rotation matrices that are transposed from the standard
 in createRotationZYX():
    
        about z-axis (yaw):           about the y-axis (pitch):    about x-axis (roll): 
            | cos φ    sin φ    0 |    |  cos ψ    0 -sin ψ |      |    1       0       0 |  
            |-sin φ    cos φ    0 |    |      0    1      0 |      |    0   cos θ   sin θ |  
            |     0        0    1 |    |  sin ψ    0  cos ψ |      |    0  -sin θ   cos θ | 
            * 
    The Shuster transposed matrices are not used in this java class.
    
 <pre>
 Some references in the "right hand system":
 
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
      
   terms used describing axes of rotation, attitute or orientation:
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
      where [nˆ]_x is the skewsymmetric matrix of nˆ

   can rotate v by 180 degrees:
      nˆ cross (nˆ cross v) = ([nˆ]_x)^2 * v = -v_perpendicular
   which shows that can also write
      v_parallel = (I + ([nˆ]_x)^2 * v)
      
 * </pre>
 *
 * TODO: add a method to extract the quaternion rotation from a 4X4 rotation matrix.
 * see Barfoot et al.
 * see http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/
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
     * calculate R(angle_z, angle_y, angle_x) = R_x(angle_x)*R_y(angle_y)*R_z(angle_z).
     * NOTE: this is not the normal convention.  You might want to use
     * createRotationZYX().
     * <pre>
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
        
        double[][] rotX = Rotation.createRollRotationMatrix(angleX);
        double[][] rotY = Rotation.createPitchRotationMatrix(angleY);
        double[][] rotZ = Rotation.createYawRotationMatrix(angleZ);
        
        double[][] r = MatrixUtil.multiply(rotX, rotY);
        r = MatrixUtil.multiply(r, rotZ);
        
        return r;
    }
    
    /**
     * create matrix for rotation about the X-axis, a.k.a. roll.
       <pre>
       about x-axis (roll):       
       |    1       0       0 |  
       |    0   cos θ  -sin θ |  
       |    0   sin θ   cos θ | 
       </pre>
     * @param angle angle of rotation about x-axis (roll) in units of radians.
     * @return 
     */
    public static double[][] createRollRotationMatrix(double angle) {
        
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
    
    /**
     * create matrix for rotation about the X-axis, a.k.a. roll.
       <pre>
       about x-axis (roll):       
       |    1       0       0 |  
       |    0   cos θ  -sin θ |  
       |    0   sin θ   cos θ | 
       </pre>
     * @param angle angle of rotation about x-axis (roll) in units of radians.
     * @param out holds values for rotation matrix for roll. 
     */
    public static void createRollRotationMatrix(double angle, double[][] out) {

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
    create matrix for rotation about the Y-axis, a.k.a. pitch.
      <pre>
      about the y-axis (pitch):
      |  cos ψ    0  sin ψ |
      |      0    1      0 |
      | -sin ψ    0  cos ψ |
      </pre>
     * @param angle angle of rotation about y-axis (pitch) in units of radians.
     * @return 
    */
    public static double[][] createPitchRotationMatrix(double angle) {
        
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
    
    /**
    create matrix for rotation about the Y-axis, a.k.a. pitch.
      <pre>
      about the y-axis (pitch):
      |  cos ψ    0  sin ψ |
      |      0    1      0 |
      | -sin ψ    0  cos ψ |
      </pre>
     * @param angle angle of rotation about y-axis (pitch) in units of radians.
     * @param out holds values for rotation matrix for pitch 
    */
    public static void createPitchRotationMatrix(double angle, 
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
    create matrix for rotation about the Z-axis, a.k.a. yaw.
      <pre>
        about z-axis (yaw):          
            | cos φ   -sin φ    0 | 
            | sin φ    cos φ    0 | 
            |     0        0    1 | 
      </pre>
     * @param angle angle of rotation about z-axis (yaw) in units of radians.
     * @return 
    */
    public static double[][] createYawRotationMatrix(double angle) {
        
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
    
    /**
    create matrix for rotation about the Z-axis, a.k.a. yaw.
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
    public static void createYawRotationMatrix(double angle, double[][] out) {

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
        //TODO: check the method createUnitLengthQuaternionBarfoot with results for quatB here
        double[] quatH = Rotation.createHamiltonQuaternionZYX(angle, axis);
        double[] quatB = Rotation.convertHamiltonToBarfootQuaternion(quatH);
        double[][] r = Rotation.createRotationMatrixFromQuaternion4(quatB);
        r = MatrixUtil.copySubMatrix(r, 0, 2, 0, 2);
        return r;
    }

    /**
     * calculate the rotation needed to tranform direction v1 to direction v2 using
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
        double[][] r = Rotation.createRodriguesFormulaRotationMatrix(axis);
        return r;
    }
    
    /**
     * given axis as an array of rotations about x, y, and z, calculate the 
     * rotation matrix.  
     * essentially, excepting a small angle correction:
     *     R^⊤ = cosθ*I + sinθ*[v]_× + (1−cosθ)*v*v^⊤
     * 
     * NOTE: if computing the partial derivative of Rotation elsewhere, 
     * can use d(R(ω)*v)/d(ω^T) = -[v]_x (see Equation (2.35) of Szeliski 2010).
     * Also note that a + sign in front of the sinθ term is used for the
     * passive system of rotations and these are used in
     * dynamics, mathematics, and computer science, etc.  A - sign 
     * in front of the sinθ term is used in aerospace engineering and medicine
     * as the later use the "left hand rule" for coordinate axes.
     * <pre>
     * references:
     *    Dmitry Berenson https://github.com/robEllenberg/comps-plugins/blob/master/python/rodrigues.py
     *    Szeliski 2010 draft "Computer Vision: Algorithms and Applications"
     *    Rodriguez’s formula (Ayache 1989)
     * </pre>
     * @param axis [1x3] array of axis of rotations about x, y, and z
     * @return rotation matrix [3X3]
     */
    public static double[][] createRodriguesFormulaRotationMatrix(double[] axis) {
        if (axis.length != 3) {
            throw new IllegalArgumentException("axis length must be 3");
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
        
        double theta = MatrixUtil.lPSum(axis, 2);
        
        double[][] tmp1, tmp2;
        
        if (theta > 1e-30) {
            
            double[] n = Arrays.copyOf(axis, axis.length);
            MatrixUtil.multiply(n, 1./theta);
            // [3X3]
            double[][] sn = MatrixUtil.skewSymmetric(n);
            
            //R = eye(3) + sin(theta)*Sn + (1-cos(theta))*dot(Sn,Sn)
            tmp1 = MatrixUtil.copy(sn);
            //[3X3]
            MatrixUtil.multiply(tmp1, Math.sin(theta));
            
            //[3X3]
            tmp2 = MatrixUtil.multiply(sn, sn);
            MatrixUtil.multiply(tmp2, 1. - Math.cos(theta));

        } else {
            
            double[][] sr = MatrixUtil.skewSymmetric(axis);
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

    /**
     * given v as an array of rotations about x, y, and z, calculate the 
     * transpose of the rotation matrix.
     * It is the same as the rodrigues rotation matrix without transposition, 
     * except that the sign in front of the sine term is negatvie.
     * R^⊤ = cosθ*I − sinθ*[v]_× + (1−cosθ)*v*v^⊤
     *     
     * <pre>
     * references:
     *    Dmitry Berenson https://github.com/robEllenberg/comps-plugins/blob/master/python/rodrigues.py
     *    Szeliski 2010 draft "Computer Vision: Algorithms and Applications"
     *    Rodriguez’s formula (Ayache 1989)
     * 
     *    Metrics for 3D Rotations: Comparison and Analysis.
     *    Huynh 2009.
     *    J Math Imaging Vis (2009) 35: 155–164 DOI 10.1007/s10851-009-0161-2
     *    https://www.cs.cmu.edu/~cga/dynopt/readings/Rmetric.pdf
     *    (note that Huynh uses Hamilton quaternions.)
     * </pre>
     * @param v [1x3] array of rotations about x, y, and z
     * @return rotation matrix [3X3]
     */
    public static double[][] createRodriguesFormulaRotationMatrixTranspose(double[] v) {
        if (v.length != 3) {
            throw new IllegalArgumentException("v length must be 3");
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
        
        double theta = MatrixUtil.lPSum(v, 2);
        
        double[][] tmp1, tmp2;
        
        if (theta > 1e-30) {
            
            double[] n = Arrays.copyOf(v, v.length);
            MatrixUtil.multiply(n, 1./theta);
            // [3X3]
            double[][] sn = MatrixUtil.skewSymmetric(n);
            
            //R = eye(3) + sin(theta)*Sn + (1-cos(theta))*dot(Sn,Sn)
            tmp1 = MatrixUtil.copy(sn);
            //[3X3]
            MatrixUtil.multiply(tmp1, Math.sin(theta));
            
            //[3X3]
            tmp2 = MatrixUtil.multiply(sn, sn);
            MatrixUtil.multiply(tmp2, 1. - Math.cos(theta));
            
            
        } else {
            
            double[][] sr = MatrixUtil.skewSymmetric(v);
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
                r[i][j] += (-tmp1[i][j] + tmp2[i][j]);
            }
        }
        
        return r;
    }

    public static class RodriguesRotation
    {
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
        public double[] om;
    }

    /**
     * calculate the Rodrigues rotation vector from the given rotation matrix.
     * The method is ported from github repositories holding the Bouguet Matlab Toolbox code, rodrigues.m.
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
    public static RodriguesRotation extractRodriguesRotationVectorBouguet(final double[][] in) throws NotConvergedException {

        //[m,n] = size(in);
        int m = in.length;
        int n = in[0].length;

        if (m != 3 || n != 3) {
            throw new IllegalArgumentException("in dimensions must be 3 X 3");
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
        rRot.om = out;
        rRot.r = MatrixUtil.copy(in);
        rRot.dRdR = dout;

        return rRot;
    }

    /**
     * calculate the rotation matrix given the Rodrigues rotation vector.
     * The method is ported from github repositories holding the Bouguet Matlab Toolbox code, rodrigues.m.
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
     * @param in [3X1] rotation vector
     * @return
     */
    public static RodriguesRotation createRodriguesRotationMatrixBouguet(final double[] in) {

        if (in.length != 3) {
            throw new IllegalArgumentException("in length must be 3");
        }

        final double eps = 2.2204e-16;

        //[m,n] = size(in);
        int m = in.length;
        //int n = 1;
        //%bigeps = 10e+4*eps;
        //bigeps = 10e+20*eps;

        double bigeps = 10e20 * eps;

        double[][] R;
        double[][] dRdin;
        //theta = norm(in);
        double theta = MatrixUtil.lPSum(in, 2);
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
            rRot.om = Arrays.copyOf(in, in.length);

            return rRot;
        }

        //%m3 = [in,theta]

        //dm3din = [eye(3);in'/theta];
        double[][] dm3din = new double[4][]; // [4X3]
        dm3din[0] = new double[]{1, 0, 0};
        dm3din[1] = new double[]{0, 1, 0};
        dm3din[2] = new double[]{0, 0, 1};
        dm3din[3] = Arrays.copyOf(in, in.length);
        MatrixUtil.multiply(dm3din[3], 1./theta);

        //omega = in/theta;
        double[] omega = Arrays.copyOf(dm3din[3], dm3din[3].length); // [3X1]

        //%m2 = [omega;theta]

        double invTheta = 1./theta;
        double invTheta2 = invTheta*invTheta;
        //dm2dm3 = [eye(3)/theta -in/theta^2; zeros(1,3) 1];// [3X3] | [3X1] ;
        double[][] dm2dm3 = new double[4][];
        dm2dm3[0] = new double[]{1*invTheta, 0, 0, in[0]*-invTheta2};
        dm2dm3[1] = new double[]{0, 1*invTheta, 0, in[1]*-invTheta2};
        dm2dm3[2] = new double[]{0, 0, 1*invTheta, in[2]*-invTheta2};
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
        rRot.om = Arrays.copyOf(in, in.length);

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
        MatrixUtil.SVDProducts svdC = null;
        try {
            svdC = MatrixUtil.performSVD(c);
        } catch (NotConvergedException e) {
            int t = 2;
        }
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
     * the convention of R_xyz =  R(theta_Z) * R(theta_Y) * R(theta_X).
     * @param r
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
     * extract euler rotation angles from a rotation matrix which has been built following
     * the convention of R(theta_Z) * R(theta_Y) * R(theta_X).
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
     * @param r
     * @return euler angles extracted from the rotation matrix under assumption
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
     * extract euler rotation angles from a rotation matrix which has been built following
     * the convention of R_X(theta_X) * R_Y(theta_Y) * R_Z(theta_Z).
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
        using euler angles
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
    public static double[] extractRodriguesRotationVector(double[][] r) {
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
     * convert the euler rotation angles to the Rodrigues vector (angle and axis)
     * with an expectation that the rotation matrix in between is created by
     * the multiplication order R_z*R_y_R_x.

     * @param eulerAngles euler rotation angles as double[]{x, y, z}
     * @return Rodriques vector.  the axis and angle representation from this
     * can be constructed as angle of rotation = r_vec/||r_vec|| and angle = ||r_vec||
     * where r_vec is the rodrigues vector.
     */
    public static double[] convertEulerAnglesToRodriguesVectorForZYX(double[] eulerAngles) {
        if (eulerAngles.length != 3) {
            throw new IllegalArgumentException("eulerAngles.length must be 3");
        }
        
        double[][] r = createRotationZYX(eulerAngles);
        
        double[] w = extractRodriguesRotationVector(r);
        
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
     * @return return the distance defined as ∥ I − r1 * r2 ∥_F if useFrobenious is true.
     * In that case the distance is within the range [0, 2sqrt(2)].
     * If useFrobenious is false, the result is spectralNorm(I − r1 * r2) which results in a distance in the range [0, 2].
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
    calculate the distance measure between 2 euclidean transformations of the
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
     * given a quaternion, return the inverse
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
       double[] qHamilton = Rotation.createHamiltonQuaternionZYX(eulerAngles);
       double[] quaternion = Rotation.convertHamiltonToBarfootQuaternion(qHamilton);
      
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
    public static double[][] createRotationMatrixFromQuaternion4(double[] quaternion) {
        if (quaternion.length != 4) {
            throw new IllegalArgumentException("quaternion must be length 4");
        }        
        double[][] lhc = quaternionLefthandCompoundOperator(quaternion);
        double[][] rhcT = quaternionRighthandCompoundOperator(quaternion);
        rhcT = MatrixUtil.transpose(rhcT);
        double[][] r = MatrixUtil.multiply(rhcT, lhc);
                
        return r;
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
     *  where R is the canonical 3X3 rotation matrix.
     * </pre>
     * @param quaternion 4X1 column vector of [eps eta] where eta is the scalar.
     * @param v3 the point as 3 element array.
     * @return 
     */
    public static double[] rotateAPointByQuaternion(double[] quaternion, double[] v3) {
        if (quaternion.length != 4) {
            throw new IllegalArgumentException("quaternion must be length 4");
        }
        if (v3.length != 3) {
            throw new IllegalArgumentException("v3 must be length 3");
        }
        double[] v = Arrays.copyOf(v3, 4);
        
        return rotateVectorByQuaternion4(quaternion, v);
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
        double[][] r = createRotationMatrixFromQuaternion4(quaternion);
        
        double[] result = MatrixUtil.multiplyMatrixByColumnVector(r, p);
        
        return result;
    }
    
     /**
     * apply a perturbation in rotation angles to the rotation matrix created by
     * the euler angles theta where the rotation matrix created is R_z*R_y*R_x.
     * To the first order, this is a constraint-sensitive approach, i.e.
     * r*r^T = I for the perturbed matrix to first order as long as it was
     * true for the given matrix r.
     * <pre>
     * from Barfoot, Forbes, & Furgale 2010, "Pose estimation using linearized 
     * rotations and quaternion algebra", Acta Astronautica (2010), doi:10.1016/j.actaastro.2010.06.049.
     * 
     * eqn (25).
     * 
     * "This update approach allows us to store and update the rotation as a 
     * rotation matrix, thereby avoiding singularities and the need to restore 
     * the constraint afterwards (i.e., constraint restoration is built in)."
     * 
     * </pre>
     * @param euler euler angles
     * @param dEuler perturbation to apply to the euler angles
     * @return resulting perturbed rotation matrix
     */
    public static double[][] applySingularitySafeEulerAnglesPerturbationZYX(double[] euler, double[] dEuler) {
        if (euler.length != 3) {
            throw new IllegalArgumentException("euler must be length 3");
        }
        if (dEuler.length != 3) {
            throw new IllegalArgumentException("dEuler must be length 3");
        }

        //[3X3]
        double[][] s = sTheta(euler);
        double[][] sds = MatrixUtil.skewSymmetric(MatrixUtil.multiplyMatrixByColumnVector(s, dEuler));
        int i;
        for (i = 0; i < sds.length; ++i) {
            sds[i][i] = 1. - sds[i][i];
        }

        double[][] r = createRotationZYX(euler);
        
        // ==== eqn (31) =======
        double[][] result = MatrixUtil.multiply(sds, r);
        
        return result;
    }
    
    /**
     * apply a perturbation in the rotation vector to a rotation matrix.
     * <pre>
     * from Barfoot, Forbes, & Furgale 2010, "Pose estimation using linearized 
     * rotations and quaternion algebra", Acta Astronautica (2010), doi:10.1016/j.actaastro.2010.06.049.
     * 
     * eqn (31).
     * "This update approach allows us to store and update the rotation as a 
     * rotation matrix, thereby avoiding singularities and the need to restore 
     * the constraint afterwards (i.e., constraint restoration is built in).
     * 
     * </pre>
     * @return 
     */
    public static double[][] applySingularitySafeRotationVectorPerturbation(
        double[] dRVector, double[][] rotation) {
        if (dRVector.length != 3) {
            throw new IllegalArgumentException("dRVector length must be 3");
        }
        if (rotation.length != 3 || rotation[0].length != 3) {
            throw new IllegalArgumentException("rotation must be [3X3]");
        }

        RodriguesRotation dRR = Rotation.createRodriguesRotationMatrixBouguet(dRVector);

        // this result should be orthonormal
        double[][] result = MatrixUtil.multiply(dRR.r, rotation);

        //NOTE: can compare to non-orthonormal results
        // eqn (26)
        //double[][] result2 = MatrixUtil.multiply(MatrixUtil.skewSymmetric(dRVector), rotation);

        return result;
    }

    /**
     * calculate S_theta which is the matrix relating angular velocity to 
     * rotation angle rates.  
     * <pre>
     * from Barfoot, Forbes, & Furgale 2010, "Pose estimation using linearized 
     * rotations and quaternion algebra", Acta Astronautica (2010), doi:10.1016/j.actaastro.2010.06.049.
     * 
     * eqn (21):
     * calc C = 3X3 rotation matrix (often written as R)
     * given array of euler rotation angles alpha, beta, gamma
     * 
     * s_theta column 0 = C_gamma(theta[2]) * C_beta(theta[1]) * [1, 0, 0]^T
     *         column 1 = C_gamma(theta[2]) * [0, 1, 0]^T 
     *         column 2 = [0, 0, 1]^T
     * 
             C_gamma                        C_beta                            C_alpha
      cc rotation about z-axis (yaw):   cc about the y-axis (pitch):    cc about x-axis (roll):    
            | cos φ   -sin φ    0 |          |  cos ψ    0  sin ψ |         |    1       0       0 |  
            | sin φ    cos φ    0 |          |      0    1      0 |         |    0   cos θ  -sin θ |  
            |     0        0    1 |          | -sin ψ    0  cos ψ |         |    0   sin θ   cos θ |  
     *  theta[0] = angleX angle of rotation about x-axis (roll) in units of radians.
     *             can use createRollRotationMatrix(theta[0])
     *  theta[1] = angleY angle of rotation about y-axis (pitch) in units of radians.
     *             can use createPitchRotationMatrix(theta[1])
     *  theta[2] = angleZ angle of rotation about z-axis (yaw) in units of radians.
     *             can use createYawRotationMatrix(theta[2])
     * 
     * see also, pp 479-480, eqn (285) of Shuster 1993, "A Survey of AttitudeRepresentations"
     * http://www.ladispe.polito.it/corsi/Meccatronica/02JHCOR/2011-12/Slides/Shuster_Pub_1993h_J_Repsurv_scan.pdf
     * though the sign conventions of the sine terms are different
     * 
     * </pre>
     * @param theta euler angles as representation of rotation matrix
     * @return [3X3]
     */
    public static double[][] sTheta(double[] theta) {
        if (theta.length != 3) {
            throw new IllegalArgumentException("theta length must be 3");
        }
        double[][] sTheta = MatrixUtil.zeros(3, 3);
        sTheta(theta, sTheta);
        
        return sTheta;
    }
    
    /**
     * calculate S_theta which is the matrix relating angular velocity to 
     * rotation angle rates.  
     * <pre>
     * from Barfoot, Forbes, & Furgale 2010, "Pose estimation using linearized 
     * rotations and quaternion algebra", Acta Astronautica (2010), doi:10.1016/j.actaastro.2010.06.049.
     * 
     * eqn (21):
     * calc C = 3X3 rotation matrix (often written as R)
     * given array of euler rotation angles  alpha, beta,gamma matrices
     * 
     * s_theta column 0 = C_gamma(theta[2]) * C_beta(theta[1]) * [1, 0, 0]^T
     *         column 1 = C_gamma(theta[2]) * [0, 1, 0]^T 
     *         column 2 = [0, 0, 1]^T
     * 
             C_gamma                        C_beta                            C_alpha
          cc rotation about z-axis (yaw):   cc about the y-axis (pitch):    cc about x-axis (roll):    
            | cos φ   -sin φ    0 |          |  cos ψ    0  sin ψ |         |    1       0       0 |  
            | sin φ    cos φ    0 |          |      0    1      0 |         |    0   cos θ  -sin θ |  
            |     0        0    1 |          | -sin ψ    0  cos ψ |         |    0   sin θ   cos θ |  
     *  theta[0] = angleX angle of rotation about x-axis (roll) in units of radians.
     *             can use createRollRotationMatrix(theta[0])
     *  theta[1] = angleY angle of rotation about y-axis (pitch) in units of radians.
     *             can use createPitchRotationMatrix(theta[1])
     *  theta[2] = angleZ angle of rotation about z-axis (yaw) in units of radians.
     *             can use createYawRotationMatrix(theta[2])
     * 
     * see also, pp 479-480, eqn (285) of Shuster 1993, "A Survey of AttitudeRepresentations"
     * http://www.ladispe.polito.it/corsi/Meccatronica/02JHCOR/2011-12/Slides/Shuster_Pub_1993h_J_Repsurv_scan.pdf
     * though the matrices are transposed from these.
     * 
     * </pre>
     * @param theta euler angles as representation of rotation matrix
     * @param output 
     */
    public static void sTheta(double[] theta, double[][] output) {
        if (theta.length != 3) {
            throw new IllegalArgumentException("theta length must be 3");
        }
        
        /*
         C_gamma(theta[2]) * C_beta(theta[1]) * [1, 0, 0]^T
        
        = | cos φ   -sin φ    0 |    *  |  cos ψ    0  sin ψ | * [1, 0, 0]^T
          | sin φ    cos φ    0 |    *  |      0    1      0 |
          |     0        0    1 |    *  | -sin ψ    0  cos ψ |
        = |  cos t2 * cos t1    -sin t2   cos t2 * sin t1 | * [1, 0, 0]^T
          |  sin t2 * cos t1     cos t2   sin t2 * sin t1 |
          | -sin t1              0        cos t1          |
        = | cos t2 * cos t1 |
          | sin t2 * cos t1 |
          | -sin t1         |
        */
        
        double[] i0 = new double[]{1, 0, 0};
        double[] i1 = new double[]{0, 1, 0};
        double[] i2 = new double[]{0, 0, 1};
        double[][] cBeta = Rotation.createPitchRotationMatrix(theta[1]);
        double[][] cGamma = Rotation.createYawRotationMatrix(theta[2]);
        
        double[] col0 = MatrixUtil.multiplyMatrixByColumnVector(
            MatrixUtil.multiply(cGamma, cBeta), i0);
        double[] col1 = MatrixUtil.multiplyMatrixByColumnVector(
            cGamma, i1);
        double[] col2 = i2;
        
        int i;
        for (i = 0; i < 3; ++i) {
            output[0][i] = col0[i];
            output[1][i] = col1[i];
            output[2][i] = col2[i];
        }        
    }
    
    /**
     * calculate  R_xyz = R_z(theta_z)*R_y(theta_y)*R_x(theta_x)
     * ,
     * that is, given an array of rotation angles, return the rotation matrix
     * as the rotations for z, y, x multiplied in that order.
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
     * </pre>
     * @param thetas euler angles
     * @return 
     */
    public static double[][] createRotationZYX(double[] thetas) {
        if (thetas.length != 3) {
            throw new IllegalArgumentException("thetas must be length 3");
        }
        
        double[][] rZ = Rotation.createYawRotationMatrix(thetas[2]);
        
        double[][] rX = Rotation.createRollRotationMatrix(thetas[0]);
        
        double[][] rY = Rotation.createPitchRotationMatrix(thetas[1]);
        
        double[][] r = MatrixUtil.multiply(rZ, MatrixUtil.multiply(rY, rX));
        
        return r;
    }
    
    /**
     * given an array of euler rotation angles, return the rotation matrix
     * as the rotations for z, y, x multiplied in that order.
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
     * </pre>
     * @param thetas euler rotation angles
     * @param aa auxiliary arrays used for internal calculations.  they're
     * meant to reduce object creation and are created by the invoking code.
     * @param out the output rotation matrix values.
     */
    public static void createRotationZYX(double[] thetas, AuxiliaryArrays aa, double[][] out) {
        if (thetas.length != 3) {
            throw new IllegalArgumentException("thetas must be length 3");
        }
        
        double[][] rZ = aa.a3X3;
        Rotation.createYawRotationMatrix(thetas[2], rZ);
        
        double[][] rY = aa.b3X3;
        Rotation.createPitchRotationMatrix(thetas[1], rY);
        
        double[][] rX = aa.c3X3;
        Rotation.createRollRotationMatrix(thetas[0], rX);
        
        double[][] rZY = aa.d3X3;
        MatrixUtil.multiply(rZ, rY, rZY);
        MatrixUtil.multiply(rZY, rX, out);        
    }
    
    /**
     * calculate  the Hamilton quaternion which would be extracted from the
     * rotation matrix created by R_xyz = R_z(theta_z)*R_y(theta_y)*R_z(theta_z).
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
         q = [scalar, vector] = [
             cos(theta/2), nx*sin(theta/2), ny*sin(theta/2), nz*sin(theta/2)];
          
     given euler angles:
         q1 =  sin(z/2)*sin(y/2)*sin(x/2) + cos(z/2)*cos(y/2)*cos(x/2)
         q2 = -sin(z/2)*sin(y/2)*cos(x/2) + sin(x/2)*cos(z/2)*cos(y/2)
         q3 =  sin(z/2)*sin(x/2)*cos(y/2) + sin(y/2)*cos(z/2)*cos(x/2)
         q4 =  sin(z/2)*cos(y/2)*cos(x/2) - sin(y/2)*sin(x/2)*cos(z/2)
     
     * </pre>
     * @param thetas euler rotation angles
     * @return 
     */
    public static double[] createHamiltonQuaternionZYX(double[] thetas) {
        if (thetas.length != 3) {
            throw new IllegalArgumentException("thetas must be length 3");
        }
        double x2 = thetas[0]/2;
        double y2 = thetas[1]/2;
        double z2 = thetas[2]/2;
        
        double q1 =  Math.sin(z2)*Math.sin(y2)*Math.sin(x2) + 
            Math.cos(z2)*Math.cos(y2)*Math.cos(x2);
        double q2 = -Math.sin(z2)*Math.sin(y2)*Math.cos(x2) + Math.sin(x2)*Math.cos(z2)*Math.cos(y2);
        double q3 =  Math.sin(z2)*Math.sin(x2)*Math.cos(y2) + Math.sin(y2)*Math.cos(z2)*Math.cos(x2);
        double q4 =  Math.sin(z2)*Math.cos(y2)*Math.cos(x2) - Math.sin(y2)*Math.sin(x2)*Math.cos(z2);
        
        return new double[]{q1, q2, q3, q4};
    }
    
    /**
     * calculate  the Hamilton quaternion which would be extracted from the
     * rotation matrix created by R_xyz = R_z(theta_z)*R_y(theta_y)*R_z(theta_z).
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
     * </pre>
     * @param angle
     * @param axis
     * @return 
     */
    public static double[] createHamiltonQuaternionZYX(double angle, double[] axis) {
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
     * convert the quaternion from [scalar vector] to [vector scalar] in-place.
     * @param qHamilton quaternion of format [scalar  vector]
     */
    public static void convertHamiltonToBarfootQuaternionInPlace(double[] qHamilton) {
        if (qHamilton.length != 4) {
            throw new IllegalArgumentException("qHamilton length must be 4");
        }
       
        double swap = qHamilton[0];
        for (int i = 1; i < 4; ++i) {
            qHamilton[i-1] = qHamilton[i];
        }
        qHamilton[3] = swap;
    }
    
    /**
     * convert the quaternion from [scalar vector] to [vector scalar].
     * @param qHamilton quaternion of format [scalar  vector]
     * @return convert the quaternion from [scalar vector] to [vector scalar]
     */
    public static double[] convertHamiltonToBarfootQuaternion(double[] qHamilton) {
        if (qHamilton.length != 4) {
            throw new IllegalArgumentException("qHamilton length must be 4");
        }
        double[] out = new double[4];
        double swap = qHamilton[0];
        for (int i = 1; i < 4; ++i) {
            out[i-1] = qHamilton[i];
        }
        out[3] = swap;
        return out;
    }
        
    /**
     * given an array of euler rotation angles, return the rotation matrix
     * as the rotations for x, y, z multiplied in that order.
     * 
     TODO: consider whether this should be using the Shuster equations
     as R_xyz = R_x(theta_x)*R_y(theta_y)*R_z(theta_z) which is transposed
     from the Z*Y*X operation.
      <pre>
      
      where  
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
         
      </pre>
      
     @param thetas euler rotation angles
     @return 
    */
    public static double[][] createRotationXYZ(double[] thetas) {
        if (thetas.length != 3) {
            throw new IllegalArgumentException("thetas must be length 3");
        }
        
        double[][] rZ = Rotation.createYawRotationMatrix(thetas[2]);
        
        double[][] rX = Rotation.createRollRotationMatrix(thetas[0]);
        
        double[][] rY = Rotation.createPitchRotationMatrix(thetas[1]);
        
        double[][] r = MatrixUtil.multiply(rX, MatrixUtil.multiply(rY, rZ));
        
        return r;
    }
    
    /*
    create a quaternion 4X1 column vector of [vector scalar] from the given
    angle and axis representation of rotation.
    Barfoot et al. place the scalar as the last item in the quaternion.
    <pre>
     * from Barfoot, Forbes, & Furgale 2010, "Pose estimation using linearized 
     * rotations and quaternion algebra", Acta Astronautica (2010), doi:10.1016/j.actaastro.2010.06.049.
     * 
     * eqn (32)
    
    satisfies unit-length constraint q^T*q = 1  (size [1X4]*[4X1]=[1X1].
    q^T*q = q_1^2 + q_2^2 + q_3^2 + q_4^2.
    </pre>
    @param unitLengthAxis axis of rotation normalized to unit length
    @param angle angle of rotation to apply
    */
    public static double[] createUnitLengthQuaternionBarfoot(double[] unitLengthAxis, double angle) {
        if (unitLengthAxis.length != 3) {
            throw new IllegalArgumentException("unitLengthAxis.length must be 3");
        }
        double cPhi2 = Math.cos(angle)/2.;
        double sPhi2 = Math.sin(angle)/2.;
        double[] out = new double[]{
            unitLengthAxis[0]*sPhi2, unitLengthAxis[1]*sPhi2, 
            unitLengthAxis[2]*sPhi2, cPhi2
        };
        return out;
    }
    
    /*
    create a quaternion 4X1 column vector of [scalar vector] from the given
    angle and axis representation of rotation.
    Hamilton and most other "right hand system" notation place the scalar as the
    first item in the quaternion.
    <pre>
     
    satisfies unit-length constraint q^T*q = 1  (size [1X4]*[4X1]=[1X1].
    q^T*q = q_1^2 + q_2^2 + q_3^2 + q_4^2.
    </pre>
    @param unitLengthAxis axis of rotation normalized to unit length
    @param angle angle of rotation to apply
    */
    public static double[] createUnitLengthQuaternionHamilton(double[] unitLengthAxis, double angle) {
        if (unitLengthAxis.length != 3) {
            throw new IllegalArgumentException("unitLengthAxis.length must be 3");
        }
        double cPhi2 = Math.cos(angle)/2.;
        double sPhi2 = Math.sin(angle)/2.;
        double[] out = new double[]{
            unitLengthAxis[0]*sPhi2, unitLengthAxis[1]*sPhi2, 
            unitLengthAxis[2]*sPhi2, cPhi2
        };
        return out;
    }
    
    /**
    <pre>
     * from Barfoot, Forbes, & Furgale 2010, "Pose estimation using linearized 
     * rotations and quaternion algebra", Acta Astronautica (2010), doi:10.1016/j.actaastro.2010.06.049.
     * 
     * eqn (14)
     C(axis, angle) = cos(angle)*I + (1-cos(angle))*axis*axis^T - sin(angle)*[axis]_x
     *
    
    satisfies unit-length constraint q^T*q = 1  (size [1X4]*[4X1]=[1X1].
    q^T*q = q_1^2 + q_2^2 + q_3^2 + q_4^2.
    </pre>
     * @param unitLengthAxis    
     * @param angle    
     * @return     
    */
    public static double[][] createRotationFromUnitLengthAngleAxis(double[] unitLengthAxis, double angle) {
        
        double ca = Math.cos(angle);
        double oneMinusCa = 1. - ca;
        double sa = Math.sin(angle);
        double[][] eye = MatrixUtil.createIdentityMatrix(3);
        double[][] skewA = MatrixUtil.skewSymmetric(unitLengthAxis);
        double[][] aaT = MatrixUtil.outerProduct(unitLengthAxis, unitLengthAxis);
        
        double[][] c = MatrixUtil.zeros(3, 3);
        
        int i, j;
        for (i = 0; i < 3; ++i) {
            for (j = 0; j < 3; ++j) {
                c[i][j] = ca*eye[i][j] + oneMinusCa*aaT[i][j] - sa*skewA[i][j];
            }
        }
        return c;
    }
    
    /**
     * apply perturbation dTheta to quaternion formed from euler rotation angles
     * theta.
     * To the first order, this is a constraint-sensitive approach.
     * <pre>
     * from Barfoot, Forbes, & Furgale 2010, "Pose estimation using linearized 
     * rotations and quaternion algebra", Acta Astronautica (2010), doi:10.1016/j.actaastro.2010.06.049.
     * 
     * eqn (45) and (44).
     * 
     * </pre>
     * @param theta euler rotation angles
     * @param dTheta rotation perturbation
     * @return the resulting quaternion formed from theta and perturbed by dTheta.  [4X1]
    */
    public static double[] applyRotationPerturbationToQuaternion(
        double[] theta, double[] dTheta) {
    
        if (theta.length != 3) {
            throw new IllegalArgumentException("theta must be length 3");
        }
        if (dTheta.length != 3) {
            throw new IllegalArgumentException("dTheta must be length 3");
        }
                 
        //[4X4]
        double[][] deltaQ = createDeltaQ(theta, dTheta);
        
        // this is eqn (37) of Barfoot et al.  [$X1]
        double[] q = createQuaternionZYXFromEuler(theta);
        
        double[] q2 = MatrixUtil.multiplyMatrixByColumnVector(deltaQ, q);
        
        return q2;
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
        double sa = Math.sin(angle);
        double ca = Math.cos(angle);
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
     * @param eulerXYZ
     * @return [4X1] quaternion rotation vector.
     */
    public static double[] createQuaternionZYXFromEuler(double[] eulerXYZ) {
        if (eulerXYZ.length != 3) {
            throw new IllegalArgumentException("eulerXYZ.length must be 3");
        }
        // length is 4
        double[] q0 = quaternionPrincipalAxisRotation(eulerXYZ[0], 0);
        double[] q1 = quaternionPrincipalAxisRotation(eulerXYZ[1], 1);
        double[] q2 = quaternionPrincipalAxisRotation(eulerXYZ[2], 2);
        
        // [4X4]
        double[][] q1m = quaternionRighthandCompoundOperator(q1);
        double[][] q2m = quaternionRighthandCompoundOperator(q2);
        
        // q2m times q1m times q0 = [4X1]
        // [4X4]
        double[][] t1 = MatrixUtil.multiply(q2m, q1m);
        // 4X1
        double[] q = MatrixUtil.multiplyMatrixByColumnVector(t1, q0);
        
        return q;
    }
 
    /**
     * create the rotation matrix from quaternions formed by the euler rotation
     * angles theta.
     * <pre>
     * from Barfoot, Forbes, & Furgale 2010, "Pose estimation using linearized 
     * rotations and quaternion algebra", Acta Astronautica (2010), doi:10.1016/j.actaastro.2010.06.049.
     * 
     * eqn (37)
     * 
     * </pre>
     * @param theta euler rotation angles
     * @return rotation matrix [3X3].
     */
    public static double[][] createRotationFromQuaternion(double[] theta) {
        if (theta.length != 3) {
            throw new IllegalArgumentException("theta.length must be 3");
        }
        
        double[] q0 = quaternionPrincipalAxisRotation(theta[0], 0);
        double[] q1 = quaternionPrincipalAxisRotation(theta[1], 1);
        double[] q2 = quaternionPrincipalAxisRotation(theta[2], 2);
       
        double[][] q0m = quaternionRighthandCompoundOperator(q0);
        double[][] q1m = quaternionRighthandCompoundOperator(q1);
        double[][] q2m = quaternionRighthandCompoundOperator(q2);
        
        double[][] q = MatrixUtil.zeros(3, 3);
        MatrixUtil.multiply(q2m, q1m, q);
        MatrixUtil.multiply(q, q0m, q);
        
        return q;
    }
    
    /**
     * create the rotation vector dPhi = S(theta) * dTheta
     * <pre>
     * from Barfoot, Forbes, & Furgale 2010, "Pose estimation using linearized 
     * rotations and quaternion algebra", Acta Astronautica (2010), doi:10.1016/j.actaastro.2010.06.049.
     * 
     * text under eqn (26)
     *    dPhi= S(theta) * dTheta
     * 
     * </pre>
     * @param theta euler rotation angles
     * @param dTheta perturbation to apply to rotation
     * @return dPhi= S(theta) * dTheta.  length is 3.
     */
    public static double[] createRotationVector(double[] theta, double[] dTheta) {
        
        double[][] sTheta = sTheta(dTheta);
        
        // length 3
        double[] dPhi = MatrixUtil.multiplyMatrixByColumnVector(sTheta, dTheta);
        
        return dPhi;
    }
    
    /**
     * create perturbation for a quaternion from a perturbation to euler
     * rotation angles.
     * <pre>
     * from Barfoot, Forbes, & Furgale 2010, "Pose estimation using linearized 
     * rotations and quaternion algebra", Acta Astronautica (2010), doi:10.1016/j.actaastro.2010.06.049.
     * 
     * part of eqn (45) and text under (26).
     *   | dPhi |^⨁
     *   |   1  |
     *      where dPhi = S(theta) * dTheta 
     * </pre>
     * @param theta euler rotation angles
     * @param dTheta rotation perturbation
     * @return perturbation for a quaternion.  [4X4].
    */
    public static double[][] createDeltaQ(double[] theta, double[] dTheta) {
        
        if (theta.length != 3) {
            throw new IllegalArgumentException("theta must be length 3");
        }
        if (dTheta.length != 3) {
            throw new IllegalArgumentException("dTheta must be length 3");
        }
                        
        double[] dPhi = createRotationVector(theta, dTheta);
        
        double[] dPhiQ = new double[4];
        
        for (int i = 0; i < 3; ++i) {
            dPhiQ[i] = dPhi[i]/2.;
        }
        dPhiQ[3] = 1;
        
        double[][] dQ = quaternionRighthandCompoundOperator(dPhiQ);
        
        return dQ;
    }
    
    /**
     * create a unit-length quaternion update of the quaternion formed
     * from the euler rotation angles theta by the perturbation dTheta.
     * To the first order, this is a constraint-sensitive approach.
     * <pre>
     * from Barfoot, Forbes, & Furgale 2010, "Pose estimation using linearized 
     * rotations and quaternion algebra", Acta Astronautica (2010), doi:10.1016/j.actaastro.2010.06.049.
     * 
     * eqn (49)
     *    q(theta + perturbation) = q(dPhi)^⨁ * q(theta)
     *        where dPhi = S(theta) * dTheta
     * 
     * "This update approach allows us to store and update the rotation as a 
     * unit-length quaternion, thereby avoiding singularities and the need to 
     * restore the constraint afterwards (i.e., constraint restoration is 
     * built in)."
     * 
     * </pre>
     * @param theta euler rotation angles
     * @param dTheta perturbation to apply to rotation
     * @return resulting quaternion from perturbation applied to quaternion 
     * formed from theta euler angles.
     */
    public static double[] applySingularitySafeRotationPerturbationQuaternion(
        double[] theta, double[] dTheta) {
        if (theta.length != 3) {
            throw new IllegalArgumentException("theta length must be 3");
        }
        if (dTheta.length != 3) {
            throw new IllegalArgumentException("theta must be length 3");
        }
           
        double[] dPhi = createRotationVector(theta, dTheta);
        double[] qDPhi = createHamiltonQuaternionZYX(dPhi);
        
        double[] qDPhiBarfoot = convertHamiltonToBarfootQuaternion(qDPhi);
                
        // ====== eqn (49) =====
        double[][] qLH = quaternionLefthandCompoundOperator(
            qDPhiBarfoot);
        
        double[] q = createHamiltonQuaternionZYX(theta);
        double[] qBarfoot = convertHamiltonToBarfootQuaternion(q);
        
        double[] result = MatrixUtil.multiplyMatrixByColumnVector(
              qLH, qBarfoot);
        
        return result;
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
