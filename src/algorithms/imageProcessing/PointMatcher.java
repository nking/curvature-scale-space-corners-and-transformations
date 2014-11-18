package algorithms.imageProcessing;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.logging.Logger;

/**
 * class to match the points extracted from two images.
 * 
 * the transformation parameters of translation, rotation and scale are 
 * found given the two sets of points.
 * 
 * <pre>
 * Details of determining the transformation from points:
 * 
 * Estimating scale:
 *     Take 2 pairs of points in both datasets and compute the distance between 
 *     them and then take the ratio:

*      scale = (distance between pair in set 1) / (distance between pair in set 2)
* 
*  Estimating rotation:
*      Take the same 2 pairs and determine the difference in their angles:
*          tan(theta) = delta y / delta x
* 
*      rotation = atan((delta y between pair in set 1)/(delta x between pair in set 1)) 
*                 -
*                 atan((delta y between pair in set 2)/(delta x between pair in set 2))
 *
 * Estimate translation:
 *     Performed on one point in set 1 with its candidate match in set 2:
 *     From the full transformation equation, we can rewrite:
 *         transX = xt0 - xc*scale - 
 *             (((x0-xc)*scale*math.cos(theta)) + ((y0-yc)*scale*math.sin(theta)))
 * 
 *         transY = yt0 - yc*scale - 
 *              ((-(x0-xc)*scale*math.sin(theta)) + ((y0-yc)*scale*math.cos(theta)))
 * 
 *         where (xc, yc) is the center of the first image
 * </pre>
 * 
 * For the contour matching, we were able to choose subsets of the entire
 * sets of contours for better matching characteristics (subsets that had larger
 * sigma values for their peaks could be used).
 * 
 * For the points given here, there isn't a clear indicator for a subset that 
 * could be used preferably so that all nodes might not need to be visited.
 * 
 * The need to use pairs of points in set1 matched against pairs in set 2
 * means that if one tries every combination of pairs, the runtime complexity 
 * is exponential.
 * 
 *     number of ways to make pairs in set 1 times the number of ways to make
 *     pairs in set 2 = 
 *           n_1!            n_2!
         ------------  X  ------------
         2*(n_1 - 2)!     2*(n_2 - 2)!
    
   This is not feasible as soon as the sets get larger than a dozen points.
     
   Alternately, one could try a search algorithm to visit the space of all 
   possible combinations of parameters instead of all combinations of subsets 
   of pairs.
    
        max translation x possible = width of image 2
        max translation y possible = height of image 2
        max scale is set to some reasonable number like 10?
        max rotation is 1 degree changes = 360 (or quadrants?)
    
        then total number of permutations 
            = maxTransX * maxTransY * maxScale * 360
            = 1000 * 500 * 10 * 360 = 1.8 billion
  
        then trying the parameters on all points puts another factor into there
        of nPoints set 1 * nPoints set 2.
        
  Because the cosine and sine terms in the transformation equation due to
  rotation work against one another and don't proceed in a total combined
  single direction with theta, we can't use a gradient descent solver.
  
  We could however use something like the Nelder-Mead Downhill Simplex
  solver with enough initialization points to cover the parameter space
  well.  The solution would not be optimal, but it should result a rough 
  solution that can be improved with another method.  
  This rough first solution could be used to remove some points from the
  dataset and the refinement could proceed on a smaller number of points.
  
       positive Y is down 
       positive X is right
       positive theta starts from Y=0, X>=0 and proceeds CW
                270
           QIII  |  QIV   
                 |
          180--------- 0   +X
                 |   
           QII   |   QI
                 90
                 +Y
           
  theta        cos(theta)  sin(theta)
  0     0        1.0         0    -----
  30    0.5236   0.87        0.5       |
  45    0.785    0.707       0.707     | QI
  60    1.047    0.5         0.866     |
  90    1.57     0.0         1.0  -----
  120   2.09    -0.5         0.866     | QII
  150   2.618   -0.866       0.5       |
  180   3.1416  -1.0         0.0  -----
  210   3.6652  -0.866      -0.5       | QIII
  240                                  |
  270   4.7124   0.0        -1.0  -----
  300   5.236    0.5        -0.866     | QIV
  330                                  |
                                  -----

 In the search through all possible parameter combinations, it looks like
 one should try combinations of scale and rotation,
 then find the best fit of translation of X and Y.
     If rotation and scale are fixed, and transX and transY are to be found,
     one can use either:
         (1) compare the offsets implied by each combination of points in set 1
             and set 2.
             for brown_lowe_2003 image1 and image2, the number of corners is
             n1 = 78
             n2 = 88
             n1 * n2 = 5616
             
             How does one determine the best solution among those translations?  
             One needs an assumed tolerance for the translation, and then to 
             count the number of matches to a point in the 2nd set.
             The best solution has the highest number of matches with the same
             rules for an assumed tolerance.
              
             For the tolerance, the reasons that translation might be different
             as a function of position in the image might be:
               -- due to rounding to a pixel.  these are very small errors.
               -- due to skew.  these are potentially very large and skew is not
                  solved in the transformation with this point matcher.  
                  (NOTE however, that once the points are matched, the true 
                  epipolar geometry can be calculated with the StereoProjection 
                  classes).
               -- due to errors in the location of the corner.  this can be due
                  to the edge detector.  these are small errors.
             
             For the Brown & Lowe 2003 points, 
                 transX=293.1 (stdDev=10.3)
                 transY=14.3 (stdDev=5.9)
             the large standard deviation appears to be due to skew.  the
             skyline is rotated about 13 degrees w.r.t. skyline in second image
             while the features at the bottom remain horizontal in both.
             
             The spread in standard deviation appears to be correlated w/
             the image dimensions, as would be expected with skew being the
             largest reason for a spread in translation.
             That is, the X axis is twice the size of the Y and so are their
             respective standard deviations.
             
             Could make a generous allowance for skew by assuming a maximum
             present such as that in the Brown & Lowe images, that is image
             size times an error due to skew spread over those pixels for
             a maximum skew such as 20 degrees or something.  
             
             If there is no reasonable solution using only scale, rotation, and 
             translation, then a more computationally expensive point matcher 
             that solves for skew too is needed (and hasn't been implemented 
             here yet).  OR, even better, an approach using contours and an
             understanding of occlusion (the later possibly requires shape
             identification) can be made with the contour matcher here.
             The contour matcher approach is currently not commonly possible
             with this project currently, because most edges are not closed 
             curves.
             
       OR
       
         (2) try every possible combination of translation for X and Y which
             would be width of image 2 in pixels times the height of image 2
             in pixels.
             for brown_lowe_2003 image1 and image2,
             image2 width = 517, height = 374
             w * h = 193358
             so this method is 34 times more
             
             The best solution among those translations is found in the same
             way as (1) above.
 * 
 * @author nichole
 */
public final class PointMatcher {
        
    private double solutionRotation = Double.MAX_VALUE;
    
    private double solutionScale = Double.MAX_VALUE;
    
    private double solutionTransX = Double.MAX_VALUE;
    
    private double solutionTransY = Double.MAX_VALUE;
        
    private final Logger log = Logger.getLogger(this.getClass().getName());
  
}
