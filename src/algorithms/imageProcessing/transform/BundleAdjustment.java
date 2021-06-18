package algorithms.imageProcessing.transform;

/**
 given intrinsic and extrinsic camera parameters, coordinates for points
 in a world reference frame, and observations of those points in one or more
 camera, return data structures needed by Levenberg-Marquardt algorithm
 in refinement of the intrinsic and extrinsic camera parameters.
 BundleAdjustment calculates partial derivatives of the parameters
 and calculates the re-projection error to form the parameter update steps,
 the gradient vector, and the evaluation of the objective (sum of squares of
 the re-projection error).
 * 
 <pre>
 additional information is present in directory doc as "bundle_adjustment.pdf"
   and "Levenberg-Marquardt_notes.pdf"
   
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
</pre>
 * @author nichole
 */
public class BundleAdjustment {
    
}
