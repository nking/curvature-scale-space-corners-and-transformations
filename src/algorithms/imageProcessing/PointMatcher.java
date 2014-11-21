package algorithms.imageProcessing;

import algorithms.util.PairIntArray;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.logging.Logger;

/**
 * class to match the points extracted from two images.
 * 
 * the transformation parameters of translation, rotation and scale are 
 * found given the two sets of points.
 * 
 <pre>
   Details of determining the transformation from points.
   
   After points are matched:
       Estimating scale:
           Take 2 pairs of points in both datasets and compute the distance between 
           them and then take the ratio:

           scale = (distance between pair in set 1) 
                   / (distance between pair in set 2)
      
           Note: scale is roughly determined from contour matching too in the
              inflection matcher.
      
       Estimating rotation:
           Take the same 2 pairs and determine the difference in their angles:
               tan(theta) = delta y / delta x
      
           rotation = atan((delta y between pair in set 1)
                          /(delta x between pair in set 1)) 
                      -
                      atan((delta y between pair in set 2)
                          /(delta x between pair in set 2))
      
       Estimate translation:
           Performed on one point in set 1 with its candidate match in set 2:
           From the full transformation equation, we can rewrite:
               transX = xt0 - xc*scale - 
                   (((x0-xc)*scale*math.cos(theta)) + ((y0-yc)*scale*math.sin(theta)))
       
               transY = yt0 - yc*scale - 
                    ((-(x0-xc)*scale*math.sin(theta)) + ((y0-yc)*scale*math.cos(theta)))
       
               where (xc, yc) is the center of the first image
 
 Matching the points:
 
   For the contour matching, we were able to choose subsets of the entire
   sets of contours for better matching characteristics (subsets that had larger
   sigma values for their peaks could be used).
   
   For the points given here, there isn't a clear indicator for a subset that 
   could be used preferably so that all nodes might not need to be visited.
   
   The need to use pairs of points in set1 matched against pairs in set 2
   means that if one tries every combination of pairs, the runtime complexity 
   is exponential.
   
       number of ways to make pairs in set 1 times the number of ways to make
       pairs in set 2 = 
             n_1!            n_2!
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
            = 1000 * 500 * 10 * 360 = 1.8 billion for pure grid search over 
              parameter space
  
        then trying the parameters on all points puts another factor into there
        of nPoints set 1 * nPoints set 2.
        
  Note, because the cosine and sine terms in the transformation equation due to
  rotation work against one another and don't proceed in a total combined
  single direction with theta, we can't use a gradient descent solver.
 
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
  So then, it looks like a grid search over rotation and scale intervals
  followed by fitting for translation to get into the local neighborhood
  of the transformation solution, followed by the use of the 
  Nelder-Mead Downhill Simplex to refine the transformation is a better solution.
  
  runtime complexity of:
  
      grid search: nRotationIntervals * nScaleIntervals * nSet1 * nSet2 
  
      downhill search: not polynomial nor deterministic.  varies by dataset.
 
     In the searches, if rotation and scale are fixed, and transX and transY 
     are to be found, one can use either:

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
             rules for an assumed tolerance and the smallest avg difference
             with predicted positions for matches in set 2.
              
             For the tolerance, the reasons that translation might be different
             as a function of position in the image might be:
               -- due to rounding to a pixel.  these are very small errors.
               -- due to projection effects for different epipolar geometry 
                  (camera nadir perpendicular to different feature, for example).  
                  these are potentially very large and are not
                  solved in the transformation with this point matcher though a
                  tolerance is made for a small amount of it.  
                  (NOTE however, that once the points are matched, the true 
                  epipolar geometry can be calculated with the StereoProjection 
                  classes).
               -- due to errors in the location of the corner.  this can be due
                  to the edge detector.  these are small errors.
               -- due to image warping such as distortions from the shape of the
                  lens.  one would need a point matcher tailored for the 
                  specific geometric projection.
               -- due to a camera not completely perpendicularly aligned with
                  the optical axis. presumably, this is an extreme case and
                  you'd want better data...
             
             For the Brown & Lowe 2003 points, 
                 transX=293.1 (stdDev=10.3)
                 transY=14.3 (stdDev=5.9)
             the large standard deviation appears to be due to projection 
             effects.  the skyline is rotated about 13 degrees w.r.t. skyline 
             in second image while the features at the bottom remain horizontal 
             in both.
             
             The spread in standard deviation appears to be correlated w/
             the image dimensions, as would be expected with projection being 
             the largest reason for a spread in translation.
             That is, the X axis is twice the size of the Y and so are their
             respective standard deviations.
             
             Could make a generous allowance for projection effects by assuming 
             a maximum present such as that in the Brown & Lowe images, that is 
             image size times an error due to a differential spread 
             over those pixels for a maximum difference in rotation such as 
             20 degrees or something.  For Brown & Lowe images, the tolerance 
             would be number of pixels for dimension times 0.02.
             
             If there is no reasonable solution using only scale, rotation, and 
             translation, then a more computationally expensive point matcher 
             that solves for epipolar geometry too or a differential rotation
             and associated variance in other parameters
             is needed (and hasn't been implemented 
             here yet).  **OR, even better, an approach using contours and an
             understanding of occlusion (the later possibly requires shape
             identification) can be made with the contour matcher in this 
             project.**
             The contour matcher approach is currently not **commonly**
             possible with this project, because most edges are not closed 
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
 </pre>

 * @author nichole
 */
public final class PointMatcher {
     
    private final Logger log = Logger.getLogger(this.getClass().getName());
 
    /**
     * calculate the rotation, scale, and translation that can be applied
     * to set1 to match points to set2 where the matches are not already known.  
     * The image1Width and imageHeight are used to create a tolerance in 
     * translation for matches.
     * 
     * NOTE: scale has be >= 1, so if one image has a smaller scale, it has to
     * be the first set given in arguments.
     * @param set1
     * @param set2
     * @param image1Width
     * @param image1Height
     * @return 
     */
    public TransformationPointFit calculateTransformation(PairIntArray set1, 
        PairIntArray set2, int image1Width, int image1Height) {
        
        TransformationPointFit fit = calculateTransformationWithGridSearch(
            set1, set2, image1Width, image1Height);
        
        if (fit == null) {
            return null;
        }
        
        int convergence = (set1.getN() < set2.getN()) ? set1.getN() 
            : set2.getN();
        
        if ((fit.getNumberOfMatchedPoints() == convergence)
            && (fit.getMeanDistFromModel() == 0)) {
            return fit;
        }

        fit = refineTransformationWithDownhillSimplex(fit.getParameters(),
            set1, set2, image1Width, image1Height);

        return fit;
    }

    private boolean fitIsBetter(TransformationPointFit bestFit, 
        TransformationPointFit compareFit) {
        
        if (compareFit == null) {
            return false;
        }
        
        int nMatches = compareFit.getNumberOfMatchedPoints();
        
        if (bestFit == null || (nMatches > bestFit.getNumberOfMatchedPoints())) {
            return true;
        } else if (nMatches == bestFit.getNumberOfMatchedPoints()) {
            if (compareFit.getMeanDistFromModel()
                < bestFit.getMeanDistFromModel()) {
                return true;
            } else if (compareFit.getMeanDistFromModel()
                == bestFit.getMeanDistFromModel()) {
                if (compareFit.getStDevFromMean() < bestFit.getStDevFromMean()) {
                    return true;
                }
            }
        }
        return false;
    }
    
    private TransformationPointFit calculateTransformationWithGridSearch(
        PairIntArray set1, PairIntArray set2, 
        int image1Width, int image1Height) {
        
        int centroidX1 = image1Width >> 1;
        
        int centroidY1 = image1Height >> 1;
        
        // projection effects from different camera positions or orientations
        // are not calculated in this point matcher, 
        // but one can set a
        // tolerance in translation to try to allow for a small amount of it.
        // 
        // as an example: rotation of 13 degrees at top of image, 
        // and 0 at bottom from turning the camera slightly during panoramic
        // photos can result in delta translation of 0.02 per pix in X and 
        // 0.016 per pix in Y
        
        double transXTol = image1Width * 0.02;
        
        double transYTol = image1Height * 0.02;
        
        int convergence = (set1.getN() < set2.getN()) ? set1.getN() 
            : set2.getN();
        
        TransformationPointFit bestFit = null;
        
        for (int rot = 0; rot < 360; rot += 10) {
            for (int scale = 1; scale < 10; scale++) {
                
                TransformationPointFit fit = calculateTranslation(set1, set2, 
                    transXTol, transYTol, 
                    rot*Math.PI/180., scale, centroidX1, centroidY1);

                if (fitIsBetter(bestFit, fit)) {
                
                    bestFit = fit;
                    
                    if ((bestFit.getNumberOfMatchedPoints() == convergence) 
                        && (bestFit.getMeanDistFromModel() == 0)) {
                        return bestFit;
                    }
                }
            }
        }
        
        return bestFit;
    }
   
    /**
     * apply the parameters to set1 and find the matches to points in set2
     * within the given tolerance for translations.
     * 
     * NOTE: scale has be >= 1, so if one image has a smaller scale, it has to
     * be the first set given in arguments.
     * 
     * ALSO NOTE: if you know a better solution exists for translation 
     * parameters that matches fewer points, but has a small avg dist from
     * model and smaller standard deviation from the avg dist from model,
     * then transXTol and transYTol should be smaller.
     * @param edges1
     * @param edges2
     * @param params
     * @param transXTol
     * @param transYTol
     * @param centroidX1
     * @param centroidY1
     * @return 
     */
    TransformationPointFit transform(PairIntArray[] edges1, 
        PairIntArray[] edges2, TransformationParameters params, 
        double transXTol, double transYTol, int centroidX1, int centroidY1) {
        
        if (edges1 == null || edges2 == null) {
            throw new IllegalArgumentException(
            "neither edges1 nor edges2 can be null");
        }
        if (edges1.length != edges2.length) {
            throw new IllegalArgumentException(
            "edges1 must be the same length as edges2");
        }
                
        double scale = params.getScale();
        
        if (scale < 1) {
            // numerical errors in rounding to integer will give wrong solutions
            throw new IllegalStateException("scale cannot be smaller than 1");
        }        
        
        int nMaxMatchable = 0;
        for (int i = 0; i < edges1.length; i++) {
            int n = edges1[i].getN();
            if (edges2[i].getN() > n) {
                n = edges2[i].getN();
            }
            nMaxMatchable += n;
        }
        
        double[] diffFromModel = new double[nMaxMatchable];
        double[] avgDiffModel = new double[1];
        
        int nMatched = 0;
        int avgDiffModelSum = 0;
        for (int i = 0; i < edges1.length; i++) {
            
            PairIntArray set1 = edges1[i];
            PairIntArray set2 = edges2[i];
            
            int n = populateDiffFromModel(set1, set2, params, 
                transXTol, transYTol, centroidX1, centroidY1,
                diffFromModel, nMatched, avgDiffModel);
            
            avgDiffModelSum += (avgDiffModel[0] * n);
            
            nMatched += n;
        }       
        
        double finalAvgDiffModel = avgDiffModelSum / (double)nMatched;
        
        double stDevDiffModel = 0;
        for (int i = 0; i < nMatched; i++) {
            double d = diffFromModel[i] - avgDiffModel[0];
            stDevDiffModel += (d * d);
        }
        
        stDevDiffModel = (Math.sqrt(stDevDiffModel/(nMatched - 1.0f)));
        
        TransformationPointFit fit = new TransformationPointFit(params, 
            nMatched, finalAvgDiffModel, stDevDiffModel);
        
        return fit;
    }

    /**
     * apply the parameters to set1 and find the matches to points in set2
     * within the given tolerance for translations.
     * 
     * NOTE: scale has be >= 1, so if one image has a smaller scale, it has to
     * be the first set given in arguments.
     * 
     * ALSO NOTE: if you know a better solution exists for translation 
     * parameters that matches fewer points, but has a small avg dist from
     * model and smaller standard deviation from the avg dist from model,
     * then transXTol and transYTol should be smaller.
     * @param set1 points matched to set 2
     * @param set2 points matched to set 1
     * @param params
     * @param transXTol
     * @param transYTol
     * @param centroidX1
     * @param centroidY1
     * @return 
     */
    TransformationPointFit transform(PairIntArray set1, PairIntArray set2, 
        TransformationParameters params, double transXTol, double transYTol, 
        int centroidX1, int centroidY1) {
        
        if (set1 == null || set2 == null) {
            throw new IllegalArgumentException(
            "neither set1 nor set2 can be null");
        }
                
        double scale = params.getScale();
        
        if (scale < 1) {
            // numerical errors in rounding to integer will give wrong solutions
            throw new IllegalStateException("scale cannot be smaller than 1");
        }        
        
        int nMaxMatchable = set1.getN();
        if (set2.getN() > nMaxMatchable) {
            nMaxMatchable = set2.getN();
        }
        
        double[] diffFromModel = new double[nMaxMatchable];
        double[] avgDiffModel = new double[1];
        
        int nMatched = populateDiffFromModel(set1, set2, 
            params, transXTol, transYTol, centroidX1, centroidY1,
            diffFromModel, 0, avgDiffModel);
        
        double stDevDiffModel = 0;
        for (int i = 0; i < nMatched; i++) {
            double d = diffFromModel[i] - avgDiffModel[0];
            stDevDiffModel += (d * d);
        }
        
        stDevDiffModel = (Math.sqrt(stDevDiffModel/(nMatched - 1.0f)));
        
        TransformationPointFit fit = new TransformationPointFit(params, 
            nMatched, avgDiffModel[0], stDevDiffModel);
        
        return fit;
    }
     
    /**
     * given the scale, rotation and set 1's reference frame centroids,
     * calculate the translation between set1 and set2 assuming that not all
     * points will match.  transXTol and transYTol allow a tolerance when
     * matching the predicted position of a point in set2.
     * 
     * NOTE: scale has be >= 1, so if one image has a smaller scale, it has to
     * be the first set given in arguments.
     * 
     * ALSO NOTE: if you know a better solution exists for translation 
     * parameters that matches fewer points, but has a small avg dist from
     * model and smaller standard deviation from the avg dist from model,
     * then transXTol and transYTol should be smaller.
     * 
     * @param set1 set of points from image 1 to match to image2.
     * @param set2 set of points from image 2 to be matched with image 1
     * @param transXTol tolerance in x for finding a match for translation.
     * For example, transXTol = image1.getWidth() * 0.02; 
     * @param transYTol tolerance in y for finding a match for translation.
     * For example, transYTol = image1.getHeight() * 0.02;
     * @param rotation given in radians with value between 0 and 2*pi, exclusive
     * @param scale
     * @param centroidX1 the x coordinate of the center of image 1 from which
     * set 1 point are from.
     * @param centroidY1 the y coordinate of the center of image 1 from which
     * set 1 point are from.
     * @return 
     */
    public TransformationPointFit calculateTranslation(PairIntArray set1, 
        PairIntArray set2, double transXTol, double transYTol, double rotation, 
        double scale, int centroidX1, int centroidY1) {
        
        if (scale < 1) {
            // numerical errors in rounding to integer will give wrong solutions
            throw new IllegalStateException("scale cannot be smaller than 1");
        }
        
        int nMax = set1.getN() * set2.getN();
        
        Map<Integer, Integer> combinationsTried = new 
            HashMap<Integer, Integer>();
        
        int maxFound = Integer.MIN_VALUE;
        TransformationPointFit bestFit = null;
        
        double scaleTimesCosine = scale * Math.cos(rotation);
        double scaleTimesSine = scale * Math.sin(rotation);
        
        for (int i = 0; i < set1.getN(); i++) {
            
            int x = set1.getX(i);
            int y = set1.getY(i);
            
            double termsX = centroidX1*scale + ( 
                ((x - centroidX1) * scaleTimesCosine) +
                ((y - centroidY1) * scaleTimesSine));
            
            double termsY = centroidY1*scale + ( 
                (-(x - centroidX1) * scaleTimesSine) +
                ((y - centroidY1) * scaleTimesCosine));
                            
            for (int j = 0; j < set2.getN(); j++) {
                
                int x2 = set2.getX(j);
                
                int y2 = set2.getY(j);
                
                int transX = (int)Math.round(x2 - termsX);
                
                int transY = (int)Math.round(y2 - termsY);
                
                // if (transX, transY) hasn't been tried:
                if (!combinationsTried.containsKey(Integer.valueOf(transX))) {
                
                    TransformationParameters params = 
                        new TransformationParameters();
                    params.setRotationInRadians((float) rotation);
                    params.setScale((float) scale);
                    params.setTranslationX(transX);
                    params.setTranslationY(transY);
                    
                    TransformationPointFit fit = transform(set1, set2, params, 
                        transXTol, transYTol, centroidX1, centroidY1);
                    
                    if (fitIsBetter(bestFit, fit)) {
                        bestFit = fit;
                    }
                    
                    combinationsTried.put(Integer.valueOf(transX), 
                        Integer.valueOf(transY));
                }
            }
        }
        
        log.fine("number of combinations tried=" + combinationsTried.size() +
            " n1*n2=" + nMax);
        
        return (bestFit != null) ? bestFit : null;
    }
    
    /**
     * refine the transformation params for set1 to better match set2.
     * @param params
     * @param set1
     * @param set2
     * @param image1Width
     * @param image1Height
     * @return 
     */
    public TransformationPointFit refineTransformationWithDownhillSimplex(
        TransformationParameters params,
        PairIntArray set1, PairIntArray set2, 
        int image1Width, int image1Height) {
        
        int centroidX1 = image1Width >> 1;
        
        int centroidY1 = image1Height >> 1;
        
        // projection effects from different camera positions or orientations
        // are not calculated in this point matcher, 
        // but one can set a
        // tolerance in translation to try to allow for a small amount of it.
        // 
        // as an example: rotation of 13 degrees at top of image, 
        // and 0 at bottom from turning the camera slightly during panoramic
        // photos can result in delta translation of 0.02 per pix in X and 
        // 0.016 per pix in Y
        
        double transXTol = image1Width * 0.02;
        
        double transYTol = image1Height * 0.02;
        
        double r = params.getRotationInRadians();
        double s = params.getScale();

        double[] drs = new double[] {
            -5.0 * Math.PI/180., -1.0 * Math.PI/180., 
            1.0 * Math.PI/180., 5.0 * Math.PI/180.
        };
        if (r == 0) {
             drs = new double[]{0};
        }

        double rMin = r - (10 * Math.PI/180);
        double rMax = r + (10 * Math.PI/180);
        double sMin = s - 0.5;
        double sMax = s + 0.5;
        
        double[] dss = new double[] {
            0.1
        };
        if (s == 1) {
            dss = new double[]{0};
            sMin = 1;
        }
        if (rMin < 0) {
            rMin = 0;
        }
        
        int n = (1 + dss.length) * (1 + drs.length);
        
        // start with simplex for at least 3 points (fitting 2 parameters)
        TransformationPointFit[] fits = new TransformationPointFit[n];
        
        int count = 0;
        for (int i = 0; i <= dss.length; i++) {
            
            double scale = (i == 0) ? s : s + dss[i - 1];
            
            for (int j = 0; j <= drs.length; j++) {
                
                double rotation = (j == 0) ? r : r + drs[j - 1];
                
                fits[count] = calculateTranslation(set1, set2, 
                    transXTol, transYTol, 
                    rotation, scale, centroidX1, centroidY1);
                
                if (fits[count] != null) {
                    count++;
                }
            }
        }
        
        if (count < n) {
            fits = Arrays.copyOf(fits, count);
        }
        
        float reflectionCoeff = 1;   // > 0
        float expansionCoeff = 2;   // > 1
        float contractionCoeff = -0.5f;
        float reductionCoeff = 0.5f;

        TransformationPointFit bestFit = fitWithDownhillSimplex(
            set1, set2, 
            fits, centroidX1, centroidY1, transXTol, transYTol,
            r, s, rMin, rMax, sMin, sMax,
            reflectionCoeff, expansionCoeff, contractionCoeff, reductionCoeff
        );
        
        return bestFit;
    }
  
    /**
     * 
     * @param set1
     * @param set2
     * @param fits
     * @param centroidX1
     * @param centroidY1
     * @param transXTol
     * @param transYTol
     * @param r
     * @param s
     * @param rMin
     * @param rMax
     * @param sMin
     * @param sMax
     * @param alpha parameter used to adjust size of changes for "reflection"
     * @param gamma parameter used to adjust size of changes for "expansion"
     * @param beta parameter used to adjust size of changes for "contraction"
     * @param tau parameter used to adjust size of changes for "reduction"
     * 
     * @return 
     */
    private TransformationPointFit fitWithDownhillSimplex(
        PairIntArray set1, PairIntArray set2, 
        TransformationPointFit[] fits, int centroidX1, int centroidY1,
        double transXTol, double transYTol,
        double r, double s,
        double rMin, double rMax, double sMin, double sMax,
        float alpha, float gamma, float beta, float tau
        ) {
        
        if (alpha < 0) {
            throw new IllegalArgumentException("alpha must be > 0");
        }
        if (gamma < 1) {
            throw new IllegalArgumentException("gamma must be > 1");
        }
        
        int convergence = (set1.getN() < set2.getN()) ? set1.getN() 
            : set2.getN();
      
        if (s == 1) {
            sMin = 1;
        }
        if (rMin < 0) {
            rMin = 0;
        }

        boolean go = true;

        int nMaxIter = 1000;
        int nIter = 0;
        
        int bestFitIdx = 0;
        int secondWorstFitIdx = fits.length - 2;
        int worstFitIdx = fits.length - 1;

        int lastNMatches = Integer.MIN_VALUE;
        double lastAvgDistModel = Double.MAX_VALUE;
        int nIterSameMin = 0;
        
        while (go && (nIter < nMaxIter)) {

            // this has changed to allow smaller avg +- stdv to be better
            // than more matches when close
            sortByDescendingMatches(fits, 0, (fits.length - 1));
            
            if ((lastNMatches == fits[0].getNumberOfMatchedPoints()) &&
                (Math.abs(lastAvgDistModel - fits[0].getMeanDistFromModel())
                < 0.01)) {
                nIterSameMin++;
                if (nIterSameMin >= 5) {
                    break;
                }
            } else {
                nIterSameMin = 0;
            }
            lastNMatches = fits[0].getNumberOfMatchedPoints();
            lastAvgDistModel = fits[0].getMeanDistFromModel();
            

            // determine center for all points excepting the worse fit
            double rSum = 0.0;
            double sSum = 0.0;
            for (int i = 0; i < (fits.length - 1); i++) {
                rSum += fits[i].getRotationInRadians();
                sSum += fits[i].getScale();
            }
            r = rSum / (double)(fits.length - 1);
            s = sSum / (double)(fits.length - 1);

            // "Reflection"
            double rReflect = r + (alpha * 
                (r - fits[worstFitIdx].getRotationInRadians()));
            double sReflect = s +
                (s - alpha * (fits[worstFitIdx].getScale()));
         
            if (sReflect < 1) {
                sReflect = 1;
            }

            TransformationPointFit fitReflected = 
                calculateTranslation(set1, set2, transXTol, transYTol, 
                    rReflect, sReflect,
                    centroidX1, centroidY1);
            
            boolean relectIsWithinBounds = 
                (rReflect >= rMin) && (rReflect <= rMax) 
                && (sReflect >= sMin) && (sReflect <= sMax);

            if (fitIsBetter(fits[secondWorstFitIdx], fitReflected)
                && relectIsWithinBounds
                && !fitIsBetter(fits[bestFitIdx], fitReflected) ) {
                
                fits[worstFitIdx] = fitReflected;
                
            } else {
                
                if (fitIsBetter(fits[bestFitIdx], fitReflected)
                    && relectIsWithinBounds) {
                    
                    // "Expansion"
                    double rExpansion = r + (gamma * 
                        (r - fits[worstFitIdx].getRotationInRadians()));
                    double sExpansion = s + (gamma * 
                        (s - fits[worstFitIdx].getScale()));
                    
                    if (sExpansion < 1) {
                        sExpansion = 1;
                    }
                    
                    TransformationPointFit fitExpansion = 
                        calculateTranslation(set1, set2, transXTol, transYTol, 
                            rExpansion, sExpansion,
                            centroidX1, centroidY1);
                    
                    if (fitIsBetter(fitReflected, fitExpansion)
                        && ((rExpansion >= rMin) && (rExpansion <= rMax)
                        && (sExpansion >= sMin) && (sExpansion <= sMax))) {

                        fits[worstFitIdx] = fitExpansion;
                        
                    } else {
                        
                        fits[worstFitIdx] = fitReflected;
                    }
                    
                } else {
                
                    // we know that the reflection fit is worse than the 2nd worse

                    // "Contraction"
                    double rContraction = r + (beta * 
                        (r - fits[worstFitIdx].getRotationInRadians()));
                    double sContraction = s + (beta * 
                        (s - fits[worstFitIdx].getScale()));
                    
                    TransformationPointFit fitContraction = 
                        calculateTranslation(set1, set2, transXTol, transYTol, 
                            rContraction, sContraction,
                            centroidX1, centroidY1);
                
                    if (sContraction < 1) {
                        sContraction = 1;
                    }

                    if (fitIsBetter(fits[worstFitIdx], fitContraction)
                        && 
                        (rContraction >= rMin) && (rContraction <= rMax) 
                        && (sContraction >= sMin) && (sContraction <= sMax)
                    ) {

                        fits[worstFitIdx] = fitContraction;
                        
                    } else {
                                                
                        // "Reduction"
                        for (int i = 1; i < fits.length; i++) {
                            
                            double rReduction = 
                                fits[bestFitIdx].getRotationInRadians()
                                + (tau * 
                                (fits[i].getRotationInRadians()
                                - fits[bestFitIdx].getRotationInRadians()));
                            double sReduction = 
                                fits[bestFitIdx].getScale()
                                + (tau * 
                                (fits[i].getScale()
                                - fits[bestFitIdx].getScale()));
                            
                            if (sReduction < 1) {
                                sReduction = 1;
                            }
                            
                            //NOTE: there's a possibility of a null fit.
                            //  instead of re-writing the fits array, will 
                            //  assign a fake infinitely bad fit which will 
                            //  fall to the bottom of the list after the next 
                            //  sort.
                            TransformationPointFit fit = calculateTranslation(
                                set1, set2, transXTol, transYTol, 
                                rReduction, sReduction, centroidX1, centroidY1);
                            
                            if (fit != null) {
                                fits[i] = fit;
                            } else {
                                fits[i] = new TransformationPointFit(
                                    new TransformationParameters(),
                                    0, Double.MAX_VALUE, Double.MAX_VALUE);
                            }
                        }
                    }
                }
            }

            log.finest("best fit so far: " + 
                fits[bestFitIdx].getNumberOfMatchedPoints());
            
            nIter++;

            if ((fits[bestFitIdx].getNumberOfMatchedPoints() == convergence) 
                && (fits[bestFitIdx].getMeanDistFromModel() == 0)) {
                go = false;
            } else if ((r > rMax) || (r < rMin)) {
                go = false;
            } else if ((s > sMax) || (s < sMin)) {
                go = false;
            }
        }
        
        // if rotation > 2PI, subtract 2PI
        if (fits[bestFitIdx].getParameters().getRotationInRadians() > 2*Math.PI) {
            float rot = fits[bestFitIdx].getParameters().getRotationInRadians();
            while (rot >= 2*Math.PI) {
                rot -= 2*Math.PI;
            }
            fits[bestFitIdx].getParameters().setRotationInRadians(rot);
        }
        
        return fits[bestFitIdx];
    }
   
    /**
     * refine the transformation params to make a better match of edges1 to
     * edges2 where the points within edges in both sets are not necessarily
     * 1 to 1 matches (that is, the input is not expected to be matched 
     * already).
     * 
     * @param edges1
     * @param edges2
     * @param params
     * @param centroidX1
     * @param centroidY1
     * @return 
     */
    public TransformationParameters refineTransformation(PairIntArray[] edges1, 
        PairIntArray[] edges2, TransformationParameters params,
        int centroidX1, int centroidY1) {
                    
        //TODO: set this empirically from tests
        double convergence = 0;
                
        double r = params.getRotationInRadians();
        double s = params.getScale();
        
        double dR = (1.0 * Math.PI/180.);
        double dS = 0.1;
        
        double rMin = r - (10 * Math.PI/180);
        double rMax = r + (10 * Math.PI/180);
        double sMin = s - 0.5;
        double sMax = s + 0.5;
        
        //TODO: these starting points could be improved.
        
        // start with simplex for at least 3 points (fitting 2 parameters)
        TransformationPointFit[] fits = new TransformationPointFit[9];
        fits[0] = transformEdges(r, s, edges1, edges2, centroidX1, centroidY1);
        fits[1] = transformEdges(r + dR, s, edges1, edges2, centroidX1, 
            centroidY1);
        fits[2] = transformEdges(r - dR, s, edges1, edges2, centroidX1, 
            centroidY1);
        // adding additional points.  need small changes more than changes
        // near the size of the value
        fits[3] = transformEdges(r, s + dS, edges1, edges2, centroidX1, 
            centroidY1);
        fits[4] = transformEdges(r, s - dS, edges1, edges2, centroidX1, 
            centroidY1);
        fits[5] = transformEdges(r + dR, s + dS, edges1, edges2, centroidX1, 
            centroidY1);
        fits[6] = transformEdges(r + dR, s - dS, edges1, edges2, centroidX1, 
            centroidY1);
        fits[7] = transformEdges(r - dR, s + dS, edges1, edges2, centroidX1, 
            centroidY1);
        fits[8] = transformEdges(r - dR, s - dS, edges1, edges2, centroidX1, 
            centroidY1);

        float alpha = 1;   // > 0
        float gamma = 2;   // > 1
        float beta = -0.5f; 
        float tau = 0.5f;
            
        boolean go = true;

        int nMaxIter = 100;
        int nIter = 0;
        
        int bestFitIdx = 0;
        int secondWorstFitIdx = fits.length - 2;
        int worstFitIdx = fits.length - 1;

        int lastNMatches = Integer.MIN_VALUE;
        double lastAvgDistModel = Double.MAX_VALUE;
        int nIterSameMin = 0;
        
        while (go && (nIter < nMaxIter)) {

            sortByDescendingMatches(fits, 0, (fits.length - 1));
            
            if ((lastNMatches == fits[bestFitIdx].getNumberOfMatchedPoints()) &&
                (Math.abs(lastAvgDistModel - fits[bestFitIdx].getMeanDistFromModel())
                < 0.01)) {
                nIterSameMin++;
                if (nIterSameMin >= 5) {
                    break;
                }
            } else {
                nIterSameMin = 0;
            }
            lastNMatches = fits[bestFitIdx].getNumberOfMatchedPoints();
            lastAvgDistModel = fits[bestFitIdx].getMeanDistFromModel();
            
            // determine center for all points excepting the worse fit
            double rSum = 0.0;
            double sSum = 0.0;
            for (int i = 0; i < (fits.length - 1); i++) {
                rSum += fits[i].getRotationInRadians();
                sSum += fits[i].getScale();
            }
            r = rSum / (double)(fits.length - 1);
            s = sSum / (double)(fits.length - 1);

            // "Reflection"
            double rReflect = r + (alpha * 
                (r - fits[worstFitIdx].getRotationInRadians()));
            double sReflect = s +
                (s - alpha * (fits[worstFitIdx].getScale()));
            
            TransformationPointFit fitReflected = transformEdges(
                rReflect, sReflect, 
                edges1, edges2, centroidX1, centroidY1);
            
            boolean relectIsWithinBounds = 
                (rReflect >= rMin) && (rReflect <= rMax) 
                && (sReflect >= sMin) && (sReflect <= sMax);

            if (fitIsBetter(fits[secondWorstFitIdx], fitReflected)
                && relectIsWithinBounds
                && !fitIsBetter(fits[bestFitIdx], fitReflected) ) {
                                
                fits[worstFitIdx] = fitReflected;
                
            } else {
                
                if (fitIsBetter(fits[bestFitIdx], fitReflected)
                    && relectIsWithinBounds) {
                    
                    // "Expansion"
                    double rExpansion = r + (gamma * 
                        (r - fits[worstFitIdx].getRotationInRadians()));
                    double sExpansion = s + (gamma * 
                        (s - fits[worstFitIdx].getScale()));
                    
                    TransformationPointFit fitExpansion = transformEdges(
                        rExpansion, sExpansion, edges1, edges2,
                        centroidX1, centroidY1);
                    
                    if (fitIsBetter(fitReflected, fitExpansion)
                        && ((rExpansion >= rMin) && (rExpansion <= rMax)
                        && (sExpansion >= sMin) && (sExpansion <= sMax))) {

                        fits[worstFitIdx] = fitExpansion;
                        
                    } else {
                        
                        fits[worstFitIdx] = fitReflected;
                    }
                    
                } else {
                
                    // the reflection fit is worse than the 2nd worse

                    // "Contraction"
                    double rContraction = r + (beta * 
                        (r - fits[worstFitIdx].getRotationInRadians()));
                    double sContraction = s + (beta * 
                        (s - fits[worstFitIdx].getScale()));
                    
                    TransformationPointFit fitContraction = transformEdges(
                        rContraction, sContraction, edges1, edges2,
                        centroidX1, centroidY1);
                    
                    if (fitIsBetter(fits[worstFitIdx], fitContraction)
                        && 
                        (rContraction >= rMin) && (rContraction <= rMax) 
                        && (sContraction >= sMin) && (sContraction <= sMax)
                    ) {

                        fits[worstFitIdx] = fitContraction;
                        
                    } else {
                        
                        // "Reduction"
                        for (int i = 1; i < fits.length; i++) {
                            
                            double rReduction = 
                                fits[bestFitIdx].getRotationInRadians()
                                + (tau * 
                                (fits[i].getRotationInRadians()
                                - fits[bestFitIdx].getRotationInRadians()));
                            
                            double sReduction = 
                                fits[bestFitIdx].getScale()
                                + (tau * 
                                (fits[i].getScale()
                                - fits[bestFitIdx].getScale()));
                                                        
                            //NOTE: there's a possibility of a null fit.
                            //  instead of re-writing the fits array, will 
                            //  assign a fake infinitely bad fit which will 
                            //  fall to the bottom of the list after the next 
                            //  sort.
                            TransformationPointFit fit = transformEdges(
                                rReduction, sReduction, 
                                edges1, edges2, centroidX1, centroidY1);
                            
                            if (fit != null) {
                                fits[i] = fit;
                            } else {
                                fits[i] = new TransformationPointFit(
                                    new TransformationParameters(),
                                    0, Double.MAX_VALUE, Double.MAX_VALUE);
                            }
                        }                        
                    }
                }
            }

            log.finest("best fit so far: nMatches=" 
                + fits[bestFitIdx].getNumberOfMatchedPoints() + 
                " diff from model=" + fits[bestFitIdx].getMeanDistFromModel()
                );
            
            nIter++;

            if ((fits[bestFitIdx].getNumberOfMatchedPoints() == convergence) 
                && (fits[bestFitIdx].getMeanDistFromModel() == 0)) {
                go = false;
            } else if ((r > rMax) || (r < rMin)) {
                go = false;
            } else if ((s > sMax) || (s < sMin)) {
                go = false;
            }
        }
        
        // if rotation > 2PI, subtract 2PI
        if (fits[bestFitIdx].getParameters().getRotationInRadians() > 2*Math.PI) 
        {
            float rot = fits[bestFitIdx].getParameters().getRotationInRadians();
            while (rot >= 2*Math.PI) {
                rot -= 2*Math.PI;
            }
            fits[bestFitIdx].getParameters().setRotationInRadians(rot);
        }
        
        return fits[bestFitIdx].getParameters();
    }
    
    private TransformationPointFit transformEdges(double rotInRad, double scl, 
        PairIntArray[] edges1, PairIntArray[] edges2, 
        int centroidX1, int centroidY1) {
        
        if (edges1 == null) {
            throw new IllegalArgumentException("edges1 cannot be null");
        }
        if (edges2 == null) {
            throw new IllegalArgumentException("edges2 cannot be null");
        }
        if ((edges1.length != edges2.length)) {
            throw new IllegalArgumentException(
            "edges1 and edges2 must be the same length");
        }
            
        /*
        Given fixed rotation and scale,
            find best translation X and Y for each edge.
        
        For final solution, could choose the best among the edges
        OR take a weighted average of the translation and re-apply it
        to all edges to return that as the fit for all edges.
        
        For now, using the longer, more accurate approach.
        */
        
        double transXTolerance = (2 * centroidX1) * 0.02;
        double transYTolerance = (2 * centroidY1) * 0.02;
                
        // make a weighted average of translations
        
        double[] transXArray = new double[edges1.length];
        double[] transYArray = new double[edges1.length];
        
        /*
        weights: 
        
        smaller dist should be larger weight, but the total has to be 1
        
        for example:
           d=[0.5, 0.3, 0.2] so the 3rd item should have largest w 
        
        the normalization is ((1/d1) + (1/d2) + (1/d3))/100.
        
        then each weight is (1/d_i)*norm
        */
        double[] w = new double[edges1.length];
        double norm = 0;
        for (int i = 0; i < edges1.length; i++) {
            
            PairIntArray edge1 = edges1[i];
            PairIntArray edge2 = edges2[i];
            
            TransformationPointFit p = calculateTranslation(
                edge1, edge2, transXTolerance, transYTolerance,
                rotInRad, scl, centroidX1, centroidY1);
            
            transXArray[i] = p.getTranslationX();            
            transYArray[i] = p.getTranslationY();
            
            w[i] = p.getMeanDistFromModel();
            if (w[i] != 0) {
                if (w[i] < 0) {
                    w[i] *= -1;
                }
                w[i] = 1./w[i];
            }
            
            norm += w[i];
        }
        
        norm /= 100.;
        
        for (int i = 0; i < edges1.length; i++) {
            w[i] *= norm;
        }
        
        double weightAvgTransX = 0;
        double weightAvgTransY = 0;
        for (int i = 0; i < edges1.length; i++) {
            weightAvgTransX += transXArray[i] * w[i];
            weightAvgTransY += transYArray[i] * w[i];
        }
        
        TransformationParameters params = new TransformationParameters();
        params.setScale((float)scl);
        params.setRotationInRadians((float)rotInRad);
        params.setTranslationX((float)weightAvgTransX);
        params.setTranslationY((float)weightAvgTransY);
        
        /*
        apply the averaged translation to all edges for the final fit
        */
        
        TransformationPointFit fit = transform(edges1, 
            edges2, params, transXTolerance, transYTolerance,
            centroidX1, centroidY1);
        
        return fit;
    }
    
    /**
     * apply the parameters to set1 and find the matches to points in set2
     * within the given tolerance for translations.
     * 
     * NOTE: scale has be >= 1, so if one image has a smaller scale, it has to
     * be the first set given in arguments.
     * 
     * ALSO NOTE: if you know a better solution exists for translation 
     * parameters that matches fewer points, but has a small avg dist from
     * model and smaller standard deviation from the avg dist from model,
     * then transXTol and transYTol should be smaller.
     * @param set1
     * @param set2
     * @param params
     * @param transXTol
     * @param transYTol
     * @param centroidX1
     * @param centroidY1
     * @param diffFromModel output variable holding difference of matched
     *    points from model
     * @param diffFromModelOffset input variable holding the offset in 
     * array diffFromModel for which points are added from this method.
     * @param avgDiffModel output variable holding the average difference from
     * the model.
     * @return number of matched points
     */
    int populateDiffFromModel(PairIntArray set1, PairIntArray set2, 
        TransformationParameters params, double transXTol, double transYTol, 
        int centroidX1, int centroidY1,
        double[] diffFromModel, final int diffFromModelOffset,
        double[] avgDiffModel) {
                
        if (diffFromModel == null) {
            throw new IllegalArgumentException("diffFromModel cannot be null");
        }
        if (diffFromModelOffset < -1) {
            throw new IllegalArgumentException(
            "diffFromModelOffset cannot be < 0");
        }
        if ((diffFromModelOffset + set1.getN()) > diffFromModel.length) {
            throw new IllegalArgumentException(
            "diffFromModelOffset + potentially all matched points is larger"
            + " than array size.");
        }
            
        double scale = params.getScale();
        double rotation = params.getRotationInRadians();
        
        if (scale < 1) {
            // numerical errors in rounding to integer will give wrong solutions
            throw new IllegalStateException("scale cannot be smaller than 1");
        }
        
        int nMatched = 0;
        double scaleTimesCosine = scale * Math.cos(rotation);
        double scaleTimesSine = scale * Math.sin(rotation);
        
        avgDiffModel[0] = 0;
        
        for (int i = 0; i < set1.getN(); i++) {
        
            int x = set1.getX(i);
            int y = set1.getY(i);
            
            double termsX = centroidX1*scale + ( 
                ((x - centroidX1) * scaleTimesCosine) +
                ((y - centroidY1) * scaleTimesSine));
            
            double termsY = centroidY1*scale + ( 
                (-(x - centroidX1) * scaleTimesSine) +
                ((y - centroidY1) * scaleTimesCosine));
            
            int xt = (int)Math.round(termsX + params.getTranslationX());
            int yt = (int)Math.round(termsY + params.getTranslationY());
            
            int lowerX = xt - (int)transXTol;
            int higherX = xt + (int)transXTol;
            int lowerY = yt - (int)transYTol;
            int higherY = yt + (int)transYTol;
            
            for (int j = 0; j < set2.getN(); j++) {
                
                int x2 = set2.getX(j);
                int y2 = set2.getY(j);
                
                if ((x2 < lowerX) || (x2 > higherX) || (y2 < lowerY) 
                    || (y2 > higherY)) {
                    
                    continue;
                }
                
                diffFromModel[diffFromModelOffset + nMatched] = Math.sqrt(
                    Math.pow(xt - x2, 2) + Math.pow(yt - y2, 2)
                );
                
                avgDiffModel[0] += diffFromModel[diffFromModelOffset + nMatched];
                
                nMatched++;
                
                break;
            }
        }
        
        avgDiffModel[0] /= (double)nMatched;
        
        return nMatched;
    }

    /**
     * sort the fits by descending number of matches. 
     * @param fits
     * @param idxLo
     * @param idxHi 
     */
    void sortByDescendingMatches(TransformationPointFit[] fits, int idxLo, 
        int idxHi) {
        
        if (idxLo < idxHi) {

            int idxMid = partition(fits, idxLo, idxHi);

            sortByDescendingMatches(fits, idxLo, idxMid - 1);

            sortByDescendingMatches(fits, idxMid + 1, idxHi);
        }
    }

    private int partition(TransformationPointFit[] fits, int idxLo, int idxHi) {
        
        TransformationPointFit x = fits[idxHi];
        
        int store = idxLo - 1;
        
        for (int i = idxLo; i < idxHi; i++) {
            if (fitIsBetter(x, fits[i])) {
                store++;
                TransformationPointFit swap = fits[store];
                fits[store] = fits[i];
                fits[i] = swap;
            }
        }
        
        store++;
        TransformationPointFit swap = fits[store];
        fits[store] = fits[idxHi];
        fits[idxHi] = swap;
        
        return store;
    }

    /**
     * create an array of values starting with start and ending before stop
     * with separations of size delta.
     * @param start start of intervals, inclusive
     * @param stop end of intervals, exclusive
     * @param delta difference between returned sequential values
     * @return 
     */
    double[] createIntervals(double start, double stop, double delta) {
        
        int n = (int)Math.ceil((stop - start)/delta);
        
        double[] values = new double[n];
        
        int count = 0;
        for (double v = start; v < stop; v += delta) { 
            values[count] = v;
            count ++;
        }
        
        return values;
    }

}
