package algorithms.imageProcessing.transform;

/**
 *
 * @author nichole
 */
public class Rectification {
    
    // notes from Serge Belongie lectures from Computer Vision II, CSE 252B, USSD
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
        //    given 2 views of a scene, thegoal is to apply projective transformations
        //    to the images so that all epipolar lines correspond to the horizontal 
        //    scan lines.
        //    so need to find 2 linear transformations, H1 and H2 that map the
        //    epipoles to infinity.
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
        //
        //
        //
        //
        //
        //
        //
        //
        //
        //
        //
        //
        //
        //
        //
        //
        //
        //
        //
        //
        //
        //
        //
        //
        //
        //
        //
        //
        //
        //
    
    
}
