package algorithms.imageProcessing.transform;

import algorithms.matrix.MatrixUtil;
import java.util.Arrays;
import no.uib.cipr.matrix.NotConvergedException;

/**
 * a utility class holding euler rotations and rodrigues formula.
 * 
 * @author nichole
 */
public class Rotation {

    /*    
    NOTE that Schuster 1993, "A Survey of Attitude Representations"
    Journal of Astronautical Sciences, Vol 41, No. 4, Oct-Dec 1993, pp 439-517
    
    uses different notation, each matrix is transposed compared to what is used
    below commonly (see http://planning.cs.uiuc.edu/node102.html
    
        about z-axis (yaw):           about x-axis (roll):       about the y-axis (pitch):
            | cos φ    sin φ    0 |    |    1       0       0 |  |  cos ψ    0 -sin ψ |
            |-sin φ    cos φ    0 |    |    0   cos θ   sin θ |  |      0    1      0 |
            |     0        0    1 |    |    0  -sin θ   cos θ |  |  sin ψ    0  cos ψ |
    */
    
    /*
    Note that the "left-hand rule", which is opposite of the convention used in
    engineering and physics, is used regarding the z-axis in
    aeromedical communities.
    
    Chapter 4, "Human Response to Acceleration" by Banks, Brinkly, Allnut, and Harding
    in "Fundamentals of Aerospace Medicine"
    
    For example, the Advisory Group for Aerospace
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
    the positive direction rightward, and the z-axis dimension
    has the positive dimension upward.    
    */
    
    /*
     terms used describing axes of rotation:
        heading, attitude, bank,
         pitch, yaw, roll,
         pan, tilt, roll
    */
    
    /**
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
    public static double[][] calculateRotationXYZ(double angleX, double angleY, double angleZ) {
        
        double[][] rotX = createEulerRollRotationMatrix(angleX);
        double[][] rotY = createEulerPitchRotationMatrix(angleY);
        double[][] rotZ = createEulerYawRotationMatrix(angleZ);
        
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
    public static double[][] createEulerRollRotationMatrix(double angle) {
        
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
    public static void createEulerRollRotationMatrix(double angle, double[][] out) {

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
    public static double[][] createEulerPitchRotationMatrix(double angle) {
        
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
    public static void createEulerPitchRotationMatrix(double angle, 
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
    public static double[][] createEulerYawRotationMatrix(double angle) {
        
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
    public static void createEulerYawRotationMatrix(double angle, double[][] out) {

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
     * given v as an array of rotations about x, y, and z, calculate the 
     * rotation matrix.
     * NOTE: if computing the partial derivative of Rotation elsewhere, 
     * can use d(R(ω)*v)/d(ω^T) = -[v]_x (see Equation (2.35) of Szeliski 2010).
     * <pre>
     * references:
     *    Dmitry Berenson https://github.com/robEllenberg/comps-plugins/blob/master/python/rodrigues.py
     *    Szeliski 2010 draft "Computer Vision: Algorithms and Applications"
     *    Rodriguez’s formula (Ayache 1989)
     * </pre>
     * @param v
     * @return 
     */
    public static double[][] createRodriguesFormulaRotationMatrix(double[] v) {
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
            double[][] sn = MatrixUtil.skewSymmetric(n);
            
            //R = eye(3) + sin(theta)*Sn + (1-cos(theta))*dot(Sn,Sn)
            tmp1 = MatrixUtil.copy(sn);
            MatrixUtil.multiply(tmp1, Math.sin(theta));
            
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
                r[i][j] += (tmp1[i][j] + tmp2[i][j]);
            }
        }
        
        return r;
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
     * @return
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
     * extract the rotation angles from the given euler rotation matrix assumed
     * to be a result of R_z(theta_z) * R_x(theta_x) * R_y(theta_y).
     * NOTE: Note, however, that this way of extracting of the Euler angles is 
     * ambiguous. Even though whatever angles you extract this way will result 
     * in the correct rotation matrix, if the latter was generated from a 
     * set of Euler angles in the first place, you are not guaranteed to get 
     * exactly those back.
     <pre>
     the method is from equation (37) from 
     lecture notes of Gordon Wetzstein at Stanford University,
     EE 267 Virtual Reality, "Course Notes: 6-DOF Pose Tracking with the VRduino",
     https://stanford.edu/class/ee267/notes/ee267_notes_tracking.pdf
     </pre>
     * @param r euler rotation matrix assumed
     * to be a result of R_z(theta_z) * R_x(theta_x) * R_y(theta_y)
     * @return theta_x, theta_y, theta_z
     */
    public static double[] extractRotationFromZXY(double[][] r) {
        if (r.length != 3 || r[0].length != 3) {
            throw new IllegalArgumentException("r must be 3x3");
        }
        double thetaX = Math.asin(r[2][1]);
        double thetaY = Math.atan2(-r[2][0], r[2][2]);
        double thetaZ = Math.atan2(-r[0][1], r[1][1]);
        return new double[]{thetaX, thetaY, thetaZ};
    }
    
    /**
     * extract euler angles from a rotation matrix which has been built following
     * the convention of R(theta_Z) * R(theta_Y) * R(theta_X).
     * @param r
     * @return 
     */
    public static double[] extractRotationFromZYX(double[][] r) {
        if (r.length != 3 || r[0].length != 3) {
            throw new IllegalArgumentException("r must be 3x3");
        }
        /*
          cc rotation about z-axis (yaw):   cc about the y-axis (pitch):    cc about x-axis (roll):    
            | cos φ   -sin φ    0 |          |  cos ψ    0  sin ψ |         |    1       0       0 |  
            | sin φ    cos φ    0 |          |      0    1      0 |         |    0   cos θ  -sin θ |  
            |     0        0    1 |          | -sin ψ    0  cos ψ |         |    0   sin θ   cos θ |  

        = | (cos φ * cos ψ)   (-sin φ)   (cos φ * sin ψ) |  * | 1       0       0 |
          | (sin φ * cos ψ)   ( cos φ)   (sin φ * sin ψ) |    | 0   cos θ   sin θ |
          | (-sin ψ)          (   0  )   (   cos ψ )     |    | 0  -sin θ   cos θ |

        = | (cos φ * cos ψ)   (-sin φ * cos θ + cos φ * sin ψ * sin θ)   (sin φ * sin θ + cos φ * sin ψ * cos θ) |
          | (sin φ * cos ψ)   ( cos φ * cos θ + sin φ * sin ψ * sin θ)   (-cos φ * sin θ + sin φ * sin ψ * cos θ)  |
          | (-sin ψ)          ( cos ψ * sin θ )                          (cos ψ * cos θ)                          |

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
                        
        return new double[]{thetaX, thetaY, thetaZ};
    }
            
    /**
     * another method to extract the Rodrigues vector from the given
     * euler rotation matrix (it's an ambiguous task).

     <pre>
     the method is from  
     lecture notes of Pradit Mittrapiyanuruk at Perdue University,
     ECE 661 Robot Vision Laboratory,
     https://engineering.purdue.edu/kak/computervision/ECE661_Fall2012/homework/hw5_LM_handout.pdf
     </pre>
     * @param r euler rotation matrix
     * @return Rodriques vector.  the axis and angle representation from this
     * can be constructed as angle of rotation = r_vec/||r_vec|| and angle = ||r_vec||
     * where r_vec is the rodrigues vector.
     */
    public static double[] extractRotation2(double[][] r) {
        if (r.length != 3 || r[0].length != 3) {
            throw new IllegalArgumentException("r must be 3x3");
        }
        /*
        theta=acos((trace(R)-1)/2);
        w=(theta/(2*sin(theta)))*[R(3,2)-R(2,3); R(1,3)-R(3,1); R(2,1)-R(1,2)];
        wx=w(1);
        wy=w(2);
        wz=w(3);
        */
        int i;
        double[] theta = new double[3];
        double[] w = new double[3];
        for (i = 0; i  < 3; ++i) {
            theta[i] = Math.acos((r[i][i] - 1.)/2.);
            w[i] = theta[i]/(2.*Math.sin(theta[i]));
        }
        w[0] *= r[2][1] - r[1][2];
        w[1] *= r[0][2] - r[2][0];
        w[2] *= r[1][0] - r[0][1];
        return w;
    }
    
    /**
     update the rotation matrix by the perturbation of dx where dx is the
     change in the 3 parameters, theta. 
     <pre>from Danping Zou lecture notes, Shanghai Jiao Tong University,
     EE382-Visual localization & Perception, “Lecture 08- Nonlinear least square & RANSAC”
     http://drone.sjtu.edu.cn/dpzou/teaching/course/lecture07-08-nonlinear_least_square_ransac.pdf
      
     parameter perturbations for a Lie group such as rotation are:
         R * (I - [delta x]_x) where [delta x]_x is the skew-symetric matrix of delta_x
     </pre>
     * @param r rotation matrix of size 3X3.
     * @param dx parameter perturbations for theta.  length is 3.
     * @return the rotation matrix perturbed by dx.
     */
    public static double[][] updateRotation(double[][] r, double[] dx) {
                
        double[][] iMinusDX = MatrixUtil.createIdentityMatrix(3);
        double[][] dxM = MatrixUtil.skewSymmetric(dx);
        double[][] factor = MatrixUtil.multiply(iMinusDX, dxM);
        r = MatrixUtil.multiply(r, factor);
        
        return r;
    }
    
    /**
     * given a quaternion, form the left-hand compound operator:
     * <pre>
     * from Barfoot, Forbes, & Furgale 2010, "Pose estimation using linearized 
     * rotations and quaternion algebra", Acta Astronautica (2010), doi:10.1016/j.actaastro.2010.06.049.
     * 
     * eqn (2):
     * given q = 4X1 column vector of [eps eta] where eta is the scalar,
     * and "1" is a 3X3 identity matrix, also written as I_3.
     * and [q]_x is the skew-symetric matrix for q.
     * 
     * The left-hand compound operator is a 4X4 matrix:
     * 
     * q^+ = [ eta*I_3-[q]_x   eps ]
     *       [ -eps^T          eta ]
     * 
     * multiplication of quaternions, u and v, which is typically written 
     * as u ⊗ v may be written equivalently as either
         u^+ v or v^⨁ u,
     * </pre>
     * @param quaternion
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
        double[][] block00 = MatrixUtil.elementwiseAdd(i3Eta, minusSkewEps);
        
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
     * @param quaternion
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
        double[][] block00 = MatrixUtil.elementwiseAdd(i3Eta, skewEps);
        
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
     * @param quaternion
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
     * @return 
     */
    public static double[] createIdentityQuaternion() {
        
        double[] i4 = new double[4];
        i4[3] = 1;
                
        return i4;
    }
    
     /**
     * create a rotation matrix from the given quaternion
     * <pre>
     * from Barfoot, Forbes, & Furgale 2010, "Pose estimation using linearized 
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
     * @param quaternion the quaternion 4 element array
     * @return a 4x4 rotation matrix
     */
    public static double[][] createRotationMatrixFromQuaternion(double[] quaternion) {
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
     * @param quaternion the quaternion 4 element array
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
     * @param quaternion the quaternion 4 element array
     * @param v the vector as a 3 element array followed by a 0.
     * @return 
     */
    public static double[] rotateVectorByQuaternion4(double[] quaternion, double[] v) {
        if (quaternion.length != 4) {
            throw new IllegalArgumentException("quaternion must be length 4");
        }
        if (v.length != 4) {
            throw new IllegalArgumentException("v3 must be length 4");
        }
        double[][] r = createRotationMatrixFromQuaternion(quaternion);
        
        double[] result = MatrixUtil.multiplyMatrixByColumnVector(r, v);
        
        return result;
    }
    
     /**
     * apply infinitesimal rotation expressed as euler angles ([1X3]) 
     * to a rotation matrix ([3X3]).
     * To the first order, this is a constraint-sensitive approach, i.e.
     * r*r^T = I for the perturbed matrix to first order as long as it was
     * true for the given matrix r.
     * <pre>
     * from Barfoot, Forbes, & Furgale 2010, "Pose estimation using linearized 
     * rotations and quaternion algebra", Acta Astronautica (2010), doi:10.1016/j.actaastro.2010.06.049.
     * 
     * eqn (31) and (30c).
     * 
     * "This update approach allows us to store and update the rotation as a 
     * rotation matrix, thereby avoiding singularities and the need to restore 
     * the constraint afterwards (i.e., constraint restoration is built in)."
     * 
     * </pre>
     * @return 
     */
    public static double[][] applyRotationPerturbationEuler(
        double[][] r, double[] theta, double[] dTheta) {
        if (r.length != 3 || r[0].length != 3) {
            throw new IllegalArgumentException("r size must be [3X3]");
        }
        if (theta.length != 3) {
            throw new IllegalArgumentException("theta must be length 3");
        }
        if (dTheta.length != 3) {
            throw new IllegalArgumentException("dTheta must be length 3");
        }
        
        // ==== eqn (30c) ======
        //3X3        
        double[] dPhi = createRotationVector(theta, dTheta);
        
        double[][] skew = MatrixUtil.skewSymmetric(dPhi);
        
        int i;
        for (i = 0; i < skew.length; ++i) {
            skew[i][i] = 1. - skew[i][i];
        }
        
        // ==== eqn (31) =======
        double[][] result = MatrixUtil.multiply(skew, r);
        
        return result;
    }
    
    /**
     * apply infinitesimal rotation expressed as euler angles ([1X3]) 
     * to a rotation matrix ([3X3]).
     * To the first order, this is a constraint-sensitive approach, i.e.
     * r*r^T = I for the perturbed matrix to first order as long as it was
     * true for the given matrix r.
     * <pre>
     * from Barfoot, Forbes, & Furgale 2010, "Pose estimation using linearized 
     * rotations and quaternion algebra", Acta Astronautica (2010), doi:10.1016/j.actaastro.2010.06.049.
     * 
     * eqn (31) and (30c).
     * 
     * "This update approach allows us to store and update the rotation as a 
     * rotation matrix, thereby avoiding singularities and the need to restore 
     * the constraint afterwards (i.e., constraint restoration is built in)."
     * 
     * </pre>
     */
    public static void applySingularitySafeRotationPerturbationEuler(
        double[][] r, double[] theta, double[] dTheta, double[][] out) {
        if (r.length != 3 || r[0].length != 3) {
            throw new IllegalArgumentException("r size must be [3X3]");
        }
        if (theta.length != 3) {
            throw new IllegalArgumentException("theta must be length 3");
        }
        if (dTheta.length != 3) {
            throw new IllegalArgumentException("dTheta must be length 3");
        }
        if (out.length != 3 || out[0].length != 3) {
            throw new IllegalArgumentException("out size must be [3X3]");
        }
        
        // ==== eqn (30c) ======
        //3X3        
        double[] dPhi = createRotationVector(theta, dTheta);
        
        double[][] skew = MatrixUtil.skewSymmetric(dPhi);
        
        // I - [dPhi]_x
        int i;
        for (i = 0; i < skew.length; ++i) {
            skew[i][i] = 1. - skew[i][i];
        }
        
        // ==== eqn (31) =======
        MatrixUtil.multiply(skew, r, out);
        
    }
    
    /**
     * apply infinitesimal rotation expressed as euler angles ([1X3]) 
     * to a rotation matrix ([3X3]).
     * To the first order, this is a constraint-sensitive approach.
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
    public static double[][] applySingularitySafeRotationPerturbationRVec(
        double[][] r, double[] dPhiRotationVector, double[] theta) {
        if (r.length != 3 || r[0].length != 3) {
            throw new IllegalArgumentException("r size must be [3X3]");
        }
        if (theta.length != 3) {
            throw new IllegalArgumentException("theta must be length 3");
        }
               
        double[][] r1 = Rotation.calculateRotationZYX(dPhiRotationVector);
        double[][] result = MatrixUtil.multiply(r1, r);
        
        return result;
    }
    
    /**
     * calculate S_theta which is the matrix relating angular velocity to 
     * Euler-angle rates.  
     * <pre>
     * from Barfoot, Forbes, & Furgale 2010, "Pose estimation using linearized 
     * rotations and quaternion algebra", Acta Astronautica (2010), doi:10.1016/j.actaastro.2010.06.049.
     * 
     * eqn (21):
     * given C = 3X3 rotation matrix (oftern written as R)
     * given array of rotation angles for euler alpha, beta,gamma matrices
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
     *             can use createEulerRollRotationMatrix(theta[0])
     *  theta[1] = angleY angle of rotation about y-axis (pitch) in units of radians.
     *             can use createEulerPitchRotationMatrix(theta[1])
     *  theta[2] = angleZ angle of rotation about z-axis (yaw) in units of radians.
     *             can use createEulerYawRotationMatrix(theta[2])
     * 
     * see also, pp 479-480, eqn (285) of Shuster 1993, "A Survey of AttitudeRepresentations"
     * http://www.ladispe.polito.it/corsi/Meccatronica/02JHCOR/2011-12/Slides/Shuster_Pub_1993h_J_Repsurv_scan.pdf
     * though the sign conventions of the sine terms are different
     * 
     * </pre>
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
     * Euler-angle rates.  
     * <pre>
     * from Barfoot, Forbes, & Furgale 2010, "Pose estimation using linearized 
     * rotations and quaternion algebra", Acta Astronautica (2010), doi:10.1016/j.actaastro.2010.06.049.
     * 
     * eqn (21):
     * given C = 3X3 rotation matrix (oftern written as R)
     * given array of rotation angles for euler alpha, beta,gamma matrices
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
     *             can use createEulerRollRotationMatrix(theta[0])
     *  theta[1] = angleY angle of rotation about y-axis (pitch) in units of radians.
     *             can use createEulerPitchRotationMatrix(theta[1])
     *  theta[2] = angleZ angle of rotation about z-axis (yaw) in units of radians.
     *             can use createEulerYawRotationMatrix(theta[2])
     * 
     * see also, pp 479-480, eqn (285) of Shuster 1993, "A Survey of AttitudeRepresentations"
     * http://www.ladispe.polito.it/corsi/Meccatronica/02JHCOR/2011-12/Slides/Shuster_Pub_1993h_J_Repsurv_scan.pdf
     * though the sign conventions of the sine terms are different
     * 
     * </pre>
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
        double[][] cBeta = createEulerPitchRotationMatrix(theta[1]);
        double[][] cGamma = createEulerYawRotationMatrix(theta[2]);
        
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
     * given an array of euler rotation angles, return thxe rotation matrix
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
     * @return 
     */
    public static double[][] calculateRotationZYX(double[] thetas) {
        if (thetas.length != 3) {
            throw new IllegalArgumentException("thetas must be length 3");
        }
        
        double[][] rZ = Rotation.createEulerYawRotationMatrix(thetas[2]);
        
        double[][] rX = Rotation.createEulerRollRotationMatrix(thetas[0]);
        
        double[][] rY = Rotation.createEulerPitchRotationMatrix(thetas[1]);
        
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
    public static void calculateRotationZYX(double[] thetas, AuxiliaryArrays aa, double[][] out) {
        if (thetas.length != 3) {
            throw new IllegalArgumentException("thetas must be length 3");
        }
        
        double[][] rZ = aa.a3X3;
        Rotation.createEulerYawRotationMatrix(thetas[2], rZ);
        
        double[][] rY = aa.b3X3;
        Rotation.createEulerPitchRotationMatrix(thetas[1], rY);
        
        double[][] rX = aa.c3X3;
        Rotation.createEulerRollRotationMatrix(thetas[0], rX);
        
        double[][] rZY = aa.d3X3;
        MatrixUtil.multiply(rZ, rY, rZY);
        MatrixUtil.multiply(rZY, rX, out);        
    }
    
    /**
     * given an array of euler rotation angles, return the rotation matrix
     * as the rotations for x, y, z multiplied in that order.
     * <pre>
     * 
     * where  
       cc rotation about z-axis (yaw):   cc about the y-axis (pitch):    cc about x-axis (roll):    
            | cos φ   -sin φ    0 |          |  cos ψ    0  sin ψ |         |    1       0       0 |  
            | sin φ    cos φ    0 |          |      0    1      0 |         |    0   cos θ  -sin θ |  
            |     0        0    1 |          | -sin ψ    0  cos ψ |         |    0   sin θ   cos θ |  
     * </pre>
     * @return 
     */
    public static double[][] calculateRotationXYZ(double[] thetas) {
        if (thetas.length != 3) {
            throw new IllegalArgumentException("thetas must be length 3");
        }
        
        double[][] rZ = Rotation.createEulerYawRotationMatrix(thetas[2]);
        
        double[][] rX = Rotation.createEulerRollRotationMatrix(thetas[0]);
        
        double[][] rY = Rotation.createEulerPitchRotationMatrix(thetas[1]);
        
        double[][] r = MatrixUtil.multiply(rX, MatrixUtil.multiply(rY, rZ));
        
        return r;
    }
    
    /*
    <pre>
     * from Barfoot, Forbes, & Furgale 2010, "Pose estimation using linearized 
     * rotations and quaternion algebra", Acta Astronautica (2010), doi:10.1016/j.actaastro.2010.06.049.
     * 
     * eqn (32)
    
    satisfies unit-length constraint q^T*q = 1  (size [1X4]*[4X1]=[1X1].
    q^T*q = q_1^2 + q_2^2 + q_3^2 + q_4^2.
    </pre>
    */
    public static double[] createUnitLengthQuaternion(double[] unitLengthAxis, double angle) {
        double cPhi2 = Math.cos(angle)/2.;
        double sPhi2 = Math.sin(angle)/2.;
        double[] out = new double[]{
            unitLengthAxis[0]*sPhi2, unitLengthAxis[1]*sPhi2, 
            unitLengthAxis[2]*sPhi2, cPhi2
        };
        return out;
    }
    
    /*
    <pre>
     * from Barfoot, Forbes, & Furgale 2010, "Pose estimation using linearized 
     * rotations and quaternion algebra", Acta Astronautica (2010), doi:10.1016/j.actaastro.2010.06.049.
     * 
     * eqn (14)
    
    satisfies unit-length constraint q^T*q = 1  (size [1X4]*[4X1]=[1X1].
    q^T*q = q_1^2 + q_2^2 + q_3^2 + q_4^2.
    </pre>
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
     * apply infinitesimal rotation expressed as euler angles ([1X3]) 
     * to a rotation matrix ([3X3]).
     * To the first order, this is a constraint-sensitive approach.
     * <pre>
     * from Barfoot, Forbes, & Furgale 2010, "Pose estimation using linearized 
     * rotations and quaternion algebra", Acta Astronautica (2010), doi:10.1016/j.actaastro.2010.06.049.
     * 
     * eqn (44).
     * 
     * </pre>
     * @return 
     */
    public static double[] applyRotationPerturbationToQuaternion(
        double[] q, double[] theta, double[] dTheta) {
        if (q.length != 4) {
            throw new IllegalArgumentException("q.length must be 4");
        }
        if (theta.length != 3) {
            throw new IllegalArgumentException("theta must be length 3");
        }
        if (dTheta.length != 3) {
            throw new IllegalArgumentException("dTheta must be length 3");
        }
                
        double[] dPhi = createRotationVector(theta, dTheta);
        
        double[][] skew = MatrixUtil.skewSymmetric(dPhi);
        
        int i;
        for (i = 0; i < skew.length; ++i) {
            skew[i][i] = 1. - skew[i][i];
        }
        
        double[] result = MatrixUtil.multiplyMatrixByColumnVector(skew, q);
        
        return result;
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
        q[principalAxis] = sa;
        q[3] = ca;
        return q;
    }
 
    /**
     * create the rotation vector dPhi = S(theta) * dTheta
     * <pre>
     * from Barfoot, Forbes, & Furgale 2010, "Pose estimation using linearized 
     * rotations and quaternion algebra", Acta Astronautica (2010), doi:10.1016/j.actaastro.2010.06.049.
     * 
     * eqn (37)
     * 
     * </pre>
     * @return rotation matrix [3X3].
     */
    public static double[][] createRotationFromEulerViaQuaternion(double[] theta) {
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
     * 
     * </pre>
     * @return dPhi= S(theta) * dTheta.  length is 3.
     */
    public static double[] createRotationVector(double[] theta, double[] dTheta) {
        double[][] sTheta = sTheta(dTheta);
        
        // length 3
        double[] dPhi = MatrixUtil.multiplyMatrixByColumnVector(sTheta, dTheta);
        
        return dPhi;
    }
    
    /**
     * apply rotation expressed as euler angles ([1X3]) to a quaternion ([4X1]).
     * To the first order, this is a constraint-sensitive approach.
     * <pre>
     * from Barfoot, Forbes, & Furgale 2010, "Pose estimation using linearized 
     * rotations and quaternion algebra", Acta Astronautica (2010), doi:10.1016/j.actaastro.2010.06.049.
     * 
     * eqn (49), (48a)
     * 
     * "This update approach allows us to store and update the rotation as a 
     * unit-length quaternion, thereby avoiding singularities and the need to 
     * restore the constraint afterwards (i.e., constraint restoration is 
     * built in)."
     * 
     * </pre>
     * @return 
     */
    public static double[] applySingularitySafeRotationPerturbationQuaternion(
        double[] q, double[] theta, double[] dTheta) {
        if (q.length != 4) {
            throw new IllegalArgumentException("q length must be 4");
        }
        if (dTheta.length != 3) {
            throw new IllegalArgumentException("theta must be length 3");
        }
        
        // ============= eqn (48c) =============
        
                
        // length 3
        double[] dPhi = createRotationVector(theta, dTheta);
       
        // length 4
        double[] qdPhi = createIdentityQuaternion();
        int i;
        for (i = 0; i < 3; ++i) {
            qdPhi[i] += 0.5*dPhi[i];
        }
        
        // ====== eqn (49) =====
        double[][] qLH = quaternionLefthandCompoundOperator(qdPhi);
        double[] result = MatrixUtil.multiplyMatrixByColumnVector(qLH, q);
        
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
    
}
