package algorithms.imageProcessing;

import algorithms.imageProcessing.util.AngleUtil;
import algorithms.util.PairIntArray;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

/**
 *
 * @author nichole
 */
public class MatchedPointsTransformationCalculator {
    
    protected transient Logger log = Logger.getLogger(
        MatchedPointsTransformationCalculator.class.getName());
    
    private boolean debug = false;
    
    public void useDebugMode() {
        debug = true;
    }
    
    /**
     * coordinate transformations from image 1 to image 2 are calculated from
     * matching lists of x, y coordinates, and given "scale" as a starting
     * parameter.  Scale is determined roughly from the contour matcher,
     * so can be used to get a rough first solution.
     * 
     * positive Y is up 
       positive X is right
       positive theta starts from Y=0, X>=0 and proceeds CW
                270
                 |     
                 |
          180--------- 0   +X
                 |   
                 |   
                 90
                 -Y
     * </pre>
     * @param scale
     * @param matchedXY1
     * @param weights1
     * @param matchedXY2
     * @param weights2
     * @param centroidX1
     * @param centroidY1
     * @return 
     */
    public TransformationParameters calulateEuclideanGivenScale(
        double scale, 
        PairIntArray matchedXY1, float[] weights1,
        PairIntArray matchedXY2, float[] weights2,
        double centroidX1, double centroidY1) {
        
        if (matchedXY1 == null) {
            throw new IllegalArgumentException("matchedXY1 cannot be null");
        }
        if (matchedXY2 == null) {
            throw new IllegalArgumentException("matchedXY2 cannot be null");
        }
        if (matchedXY1.getN() != matchedXY2.getN()) {
            throw new IllegalArgumentException(
                "matchedXY1 and matchedXY2 must have same number of points");
        }
        if (weights1 == null) {
            throw new IllegalArgumentException("weights1 cannot be null");
        }
        if (weights2 == null) {
            throw new IllegalArgumentException("weights2 cannot be null");
        }
        if (weights1.length != weights2.length) {
            throw new IllegalArgumentException(
                "weights1 and weights2 must have same number of points");
        }
        if (matchedXY1.getN() < 2) {
            return null;
        }
        
        log.info("start solution for " + matchedXY1.getN() + " points");
        
        /*
        solve for rotation.
        
        Take the same 2 pairs int both imagesand get the difference in their 
        angles:
            tan(theta) = y / x

        For example:
            theta of pair in image1:
                theta = math.atan( (y1-y0)/(x1-x0) )
                      = 0.7853981633974483 radians
                      = 45 degrees

            theta of pair in image2:
                theta = math.atan( (yt1-yt0)/(xt1-xt0) )
                      = 0.3490522203358645
                      = 20.0

            rotation = theta_image1 - theta_image2 = 25 degrees
        */
            
        /*
        discard outside avg +- stdev
        */
        
        AngleUtil angleUtil = new AngleUtil();
        
        double[] thetas = new double[matchedXY1.getN()];
        double[] scales = new double[matchedXY1.getN()];
        double thetaSum = 0;
        double scaleSum = 0;
        for (int i = 0; i < matchedXY1.getN(); i++) {
            int x0im1 = matchedXY1.getX(i);
            int y0im1 = matchedXY1.getY(i);
            int x0im2 = matchedXY2.getX(i);
            int y0im2 = matchedXY2.getY(i);
            int x1im1, y1im1, x1im2, y1im2;
            if ((i + 1) == matchedXY1.getN()) {
                x1im1 = matchedXY1.getX(0);
                y1im1 = matchedXY1.getY(0);
                x1im2 = matchedXY2.getX(0);
                y1im2 = matchedXY2.getY(0);
            } else {
                x1im1 = matchedXY1.getX(i + 1);
                y1im1 = matchedXY1.getY(i + 1);
                x1im2 = matchedXY2.getX(i + 1);
                y1im2 = matchedXY2.getY(i + 1);
            }
            double diffX1 = (x1im1 - x0im1);
            double diffY1 = (y1im1 - y0im1);
            
            double diffX2 = (x1im2 - x0im2);
            double diffY2 = (y1im2 - y0im2);
                        
            double t = angleUtil.subtract(diffX1, diffY1, diffX2, diffY2);
            
            thetas[i] = t;
            
            thetaSum += thetas[i];
            
            double lenim1 = Math.sqrt(Math.pow(diffX1, 2) 
                + Math.pow(diffY1, 2));
            double lenim2 = Math.sqrt(Math.pow(diffX2, 2) 
                + Math.pow(diffY2, 2));
            scales[i] = lenim2/lenim1;
            scaleSum += scales[i];
        }
        
        double avgScale = scaleSum / (double)matchedXY1.getN();
        double avgTheta = thetaSum / (double)matchedXY1.getN();
        double stDevThetaSum = 0;
        double stDevScaleSum = 0;
        for (int i = 0; i < matchedXY1.getN(); i++) {            
            stDevThetaSum += Math.pow(thetas[i] - avgTheta, 2); 
            stDevScaleSum += Math.pow(scales[i] - avgScale, 2); 
        }
        double stDevTheta= Math.sqrt(stDevThetaSum/(matchedXY1.getN() - 1));
        double stDevScale= Math.sqrt(stDevScaleSum/(matchedXY1.getN() - 1));
        
        double rotSum = 0;
        double rCount = 0;
        scaleSum = 0;
        double sCount = 0;
        
        List<Integer> rm = new ArrayList<Integer>();
        for (int i = 0; i < matchedXY1.getN(); i++) {
            double dss = Math.abs(scales[i] - avgScale);
            if (dss > 1.5*stDevScale) {
                rm.add(Integer.valueOf(i));
                continue;
            }
            double dtt = Math.abs(thetas[i] - avgTheta);
            if (dtt > 1.5*stDevTheta) {
                rm.add(Integer.valueOf(i));
                continue;
            }
            
log.info("scl=" + scales[i] + " stDevScale=" + stDevScale
+ " abs(scale-avg)=" + dss);
            
            scaleSum += scales[i];
            sCount++;
           
log.info("rot=" + thetas[i] + " stDevTheta=" + stDevTheta
+ " abs(theta-avg)=" + dtt 
+ " w1=" + weights1[i] + " w2=" + weights2[i]);
            
            rotSum += thetas[i];
            rCount++;
        }
        
        if (!rm.isEmpty()) {
            
            // remove from datasets and re-try
            int nrm = rm.size();
            PairIntArray xy1 = new PairIntArray();
            PairIntArray xy2 = new PairIntArray();
            float[] w1 = new float[weights1.length - nrm];
            float[] w2 = new float[weights2.length - nrm];
            
            int ii = 0;
            for (int i = 0; i < matchedXY1.getN(); i++) {
                if (rm.contains(Integer.valueOf(i))) {
                    continue;
                }
                xy1.add(matchedXY1.getX(i), matchedXY1.getY(i));
                xy2.add(matchedXY2.getX(i), matchedXY2.getY(i));
                w1[ii] = weights1[i];
                w2[ii] = weights2[i];
                ii++;
            }
            
            return calulateEuclideanGivenScale(scale, xy1, w1,
                xy2, w2, centroidX1, centroidY1);
        }
        
        double theRotation = rotSum/rCount;
        double theScale = scaleSum/sCount;
        
        log.info("given scale=" + scale + " found scale=" + theScale);
        log.info("rotation = " + theRotation);

        if (Math.abs(theScale - scale) > scale*0.1) {
            log.warning("the differences in estimated scale and given scale are"
                + " large.  this can happen if the given scale value was " 
                + " determined from contour matching."
                + " the estimate here uses pairs of points which may"
                + " be close to one another.  choosing the scale given to "
                + " the method and continuing.");
        }
        
        theScale = scale;
        
        /*
        estimate translation:
        
        transX = xt0 - xc*scale - (((x0-xc)*scale*math.cos(theta)) 
            + ((y0-yc)*scale*math.sin(theta)))
        
        transY = yt0 - yc*scale - ((-(x0-xc)*scale*math.sin(theta)) 
            + ((y0-yc)*scale*math.cos(theta)))
        */
        double mc = Math.cos(theRotation);
        double ms = Math.sin(theRotation);
        double transXSum = 0;
        double transYSum = 0;
        for (int i = 0; i < matchedXY1.getN(); i++) {
            int xim1 = matchedXY1.getX(i);
            int yim1 = matchedXY1.getY(i);
            int xim2 = matchedXY2.getX(i);
            int yim2 = matchedXY2.getY(i);
            double transX = xim2 - centroidX1*theScale 
                - (((xim1 - centroidX1) * theScale*mc) 
                + ((yim1 - centroidY1) *theScale*ms));
            
            double transY = yim2 - centroidY1*theScale 
                - ((-(xim1 - centroidX1) *theScale*ms) 
                + ((yim1 - centroidY1) *theScale*mc));
            
            transXSum += transX;
            transYSum += transY;
        }
        double theTranslationX = transXSum/(double)matchedXY1.getN();
        double theTranslationY = transYSum/(double)matchedXY1.getN();
        
        TransformationParameters params = new TransformationParameters();
        params.setRotationInRadians((float)theRotation);
        params.setScale((float)theScale);
        params.setTranslationX((float)theTranslationX);
        params.setTranslationY((float)theTranslationY);
        
        if (debug) {
            log.info("params: " + params.toString());
        }
        
        return params;
    }

    /**
     * from a set of transformation parameters params that transform
     * points in reference frame 1 into reference frame 2,
     * transform the point (x1, y1) into reference frame 2.
     * 
     * @param params
     * @param centroidX1 x centroid of reference frame 1
     * @param centroidY1 y centroid of reference frame 1
     * @param x1 the x coordinate of point (x1, y1) to transform
     * @param y1 the y coordinate of point (x1, y1) to transform
     * @return the coordinates of (x1, y1) transformed into reference frame 2
     */
    public double[] applyTransformation(TransformationParameters params, 
        int centroidX1, int centroidY1,
        double x1, double y1) {
        
        double rot = params.getRotationInRadians();
        double scl = params.getScale();
                
        double mc = Math.cos(rot);
        double ms = Math.sin(rot);
                
        double x2 = centroidX1 * scl
            + (((x1 - centroidX1) * scl * mc)
            + ((y1 - centroidY1) * scl * ms));
        x2 += params.getTranslationX();

        double y2 = centroidY1 * scl
            + ((-(x1 - centroidX1) * scl * ms)
            + ((y1 - centroidY1) * scl * mc));
        y2 += params.getTranslationY();
        
        return new double[]{x2, y2};
    }

    /**
     * from a set of transformation parameters params that transform
     * points in reference frame 1 into reference frame 2, create
     * a transformation that can transform points in reference frame 2
     * into reference frame 1.  (x1, y1) and (x2, y2) are the same point
     * in both frames and are used to calculate the translation.
     * @param params transformation parameters to apply to points in reference
     * frame 1 to put them in reference frame 2
     * @param centroidX2 the x centroid of reference frame 2
     * @param centroidY2 the y centroid of reference frame 2
     * @param x1 x coordinate of point (x1, y1) in reference frame 1
     * @param y1 y coordinate of point (x1, y1) in reference frame 1
     * @param x2 x coordinate of point (x2, y2) in reference frame 2, where 
     *    (x2, y2) is the (x1, y1) transformed to reference frame 2
     * @param y2 y coordinate of point (x2, y2) in reference frame 2, where
     *    (x2, y2) is the (x1, y1) transformed to reference frame 2
     * @return transformation parameters that can transform points in reference
     * frame 2 into reference frame 1
     */
    TransformationParameters swapReferenceFrames(TransformationParameters 
        params, int centroidX2, int centroidY2,
        double x1, double y1, double x2, double y2) {
        
        /*
        use (centroidX of edge1, centroidY of edge1) for (x0, y0) 
        to solve for (x1, y1) here and then reverse and invert the 
        rotation and scale to solve for the "reverse" translations.

        xr_0 = xc*scale + (((x0-xc)*scale*math.cos(theta)) 
            + ((y0-yc)*scale*math.sin(theta)))

        xt_0 = xr_0 + transX = x1

        yr_0 = yc*scale + ((-(x0-xc)*scale*math.sin(theta)) 
            + ((y0-yc)*scale*math.cos(theta)))

        yt_0 = yr_0 + transY = y1
        */
     
        double revRot = -1 * params.getRotationInRadians();
        if (revRot < 0) {
            revRot += 2. * Math.PI;
        }
        double revScale = 1. / params.getScale();

        double rmc = Math.cos(revRot);
        double rms = Math.sin(revRot);

        double x1_ = revScale * (centroidX2 
            + (((x2 - centroidX2) * rmc) + ((y2 - centroidY2) * rms)));
        // (x1,y1) transformed to image 1 needs translation to equal (x0, y0)
        double revTransX = (x1 - x1_);

        double y1_ = revScale * (centroidY2
            + ((-(x2 - centroidX2) * rms) + ((y2 - centroidY2) * rmc)));
        double revTransY = (y1 - y1_);

        TransformationParameters paramsRev = new TransformationParameters();
        paramsRev.setScale((float)revScale);
        paramsRev.setRotationInRadians((float)revRot);
        paramsRev.setTranslationX((float)revTransX);
        paramsRev.setTranslationY((float)revTransY);
        
        return paramsRev;
    }
    
    /**
     * from a set of transformation parameters params that transform
     * points in reference frame 1 into reference frame 2, create
     * a transformation that can transform points in reference frame 2
     * into reference frame 1.
     * @param params transformation parameters to apply to points in reference
     * frame 1 to put them in reference frame 2
     * 
     * @return transformation parameters that can transform points in reference
     * frame 2 into reference frame 1
     */
    public TransformationParameters swapReferenceFrames(TransformationParameters 
        params) {
        
        if (params == null) {
            throw new IllegalArgumentException("params cannot be null");
        }
        
        if (params.getOriginX() == 0 && params.getOriginY() == 0) {
            return swapReferenceFramesWRTOrigin(params);
        }
        
        /*
        use (centroidX of edge1, centroidY of edge1) for (x0, y0) 
        to solve for (x1, y1) here and then reverse and invert the 
        rotation and scale to solve for the "reverse" translations.

        xr_0 = xc*scale + (((x0-xc)*scale*math.cos(theta)) 
            + ((y0-yc)*scale*math.sin(theta)))

        xt_0 = xr_0 + transX = x1

        yr_0 = yc*scale + ((-(x0-xc)*scale*math.sin(theta)) 
            + ((y0-yc)*scale*math.cos(theta)))

        yt_0 = yr_0 + transY = y1
        */
        
        float refX2 = (params.getOriginX() * params.getScale()) +
            params.getTranslationX();
        
        float refY2 = (params.getOriginY() * params.getScale()) +
            params.getTranslationY();
        
        double revRot = -1 * params.getRotationInRadians();
        if (revRot < 0) {
            revRot += 2. * Math.PI;
        }
        double revScale = 1. / params.getScale();

        TransformationParameters paramsRev = new TransformationParameters();
        paramsRev.setScale((float)revScale);
        paramsRev.setRotationInRadians((float)revRot);
        paramsRev.setTranslationX(-params.getTranslationX());
        paramsRev.setTranslationY(-params.getTranslationY());
        paramsRev.setOriginX(-refX2);
        paramsRev.setOriginY(-refY2);
        
        return paramsRev;
    }
    
    /**
     * from a set of transformation parameters params that transform
     * points in reference frame 1 into reference frame 2, create
     * a transformation that can transform points in reference frame 2
     * into reference frame 1.  The center of the rotation of frame 1 has
     * to have been point (0, 0).  This method calculates where that 
     * transformed point would be in frame 2 and that becomes the center of
     * the reverse rotation.
     * @param params transformation parameters to apply to points in reference
     * frame 1 to put them in reference frame 2
     * @param centroidX2 the x centroid of reference frame 2
     * @param centroidY2 the y centroid of reference frame 2
     * @param x1 x coordinate of point (x1, y1) in reference frame 1
     * @param y1 y coordinate of point (x1, y1) in reference frame 1
     * @param x2 x coordinate of point (x2, y2) in reference frame 2, where 
     *    (x2, y2) is the (x1, y1) transformed to reference frame 2
     * @param y2 y coordinate of point (x2, y2) in reference frame 2, where
     *    (x2, y2) is the (x1, y1) transformed to reference frame 2
     * @return transformation parameters that can transform points in reference
     * frame 2 into reference frame 1
     */
    TransformationParameters swapReferenceFramesWRTOrigin(TransformationParameters 
        params) {
        
        if (params == null) {
            throw new IllegalArgumentException("params cannot be null");
        }
        
        /*
        use (centroidX of edge1, centroidY of edge1) for (x0, y0) 
        to solve for (x1, y1) here and then reverse and invert the 
        rotation and scale to solve for the "reverse" translations.

        xr_0 = xc*scale + (((x0-xc)*scale*math.cos(theta)) 
            + ((y0-yc)*scale*math.sin(theta)))

        xt_0 = xr_0 + transX = x1

        yr_0 = yc*scale + ((-(x0-xc)*scale*math.sin(theta)) 
            + ((y0-yc)*scale*math.cos(theta)))

        yt_0 = yr_0 + transY = y1
        
        ---------
        simplified for origin:
        xr_0 = 0, yr_0 = 0 and same for set2
        */
     
        double revRot = -1 * params.getRotationInRadians();
        if (revRot < 0) {
            revRot += 2. * Math.PI;
        }
        double revScale = 1. / params.getScale();

        double rmc = Math.cos(revRot);
        double rms = Math.sin(revRot);

        TransformationParameters paramsRev = new TransformationParameters();
        paramsRev.setScale((float)revScale);
        paramsRev.setRotationInRadians((float)revRot);
        paramsRev.setTranslationX(-1*params.getTranslationX());
        paramsRev.setTranslationY(-1*params.getTranslationY());
        
        return paramsRev;
    }
}
