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
     * @param matchedXY2
     * @param centroidX1
     * @param centroidY1
     * @return 
     */
    public TransformationParameters calulateEuclideanGivenScale(
        double scale, PairIntArray matchedXY1, PairIntArray matchedXY2, 
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
+ " abs(theta-avg)=" + dtt);
            
            rotSum += thetas[i];
            rCount++;
        }
        
        if (!rm.isEmpty()) {
            
            PairIntArray xy1 = new PairIntArray();
            PairIntArray xy2 = new PairIntArray();
            
            for (int i = 0; i < matchedXY1.getN(); i++) {
                if (rm.contains(Integer.valueOf(i))) {
                    continue;
                }
                xy1.add(matchedXY1.getX(i), matchedXY1.getY(i));
                xy2.add(matchedXY2.getX(i), matchedXY2.getY(i));
            }
            
            return calulateEuclideanGivenScale(scale, xy1,
                xy2, centroidX1, centroidY1);
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
          
            double trX1 = centroidX1*scale + (((xim1 - centroidX1) * scale*mc) 
                + ((yim1 - centroidY1) *scale*ms));
      
            double trY1 = centroidY1*scale + ((-(xim1 - centroidX1) * scale*ms) 
                + ((yim1 - centroidY1) * scale*mc));
        
            double transX = xim2 - trX1;

            double transY = yim2 - trY1;
            
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
        params.setOriginX((float)centroidX1);
        params.setOriginY((float)centroidY1);
        
        if (debug) {
            log.info("params: " + params.toString());
        }
        
        return params;
    }
    
     /**
     * coordinate transformations from pair 1 to pair 2 are calculated 
     * given "scale".
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
     * @param centroidX1
     * @param centroidY1
     * @return 
     */
    public TransformationParameters calulateEuclideanGivenScale(
        int set1X1, int set1Y1, int set1X2, int set1Y2,
        int set2X1, int set2Y1, int set2X2, int set2Y2,
        double centroidX1, double centroidY1) {
                
        /*
        solve for rotation.
        
        Take the same 2 pairs int both images and get the difference in their 
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
        
        AngleUtil angleUtil = new AngleUtil();
        
        double dx1 = set1X1 - set1X2;
        double dy1 = set1Y1 - set1Y2;
        
        double dx2 = set2X1 - set2X2;
        double dy2 = set2Y1 - set2Y2;
        
        double theta = angleUtil.subtract(dx1, dy1, dx2, dy2);
        
        double sep1 = Math.sqrt((dx1*dx1) + (dy1*dy1));
        
        double sep2 = Math.sqrt((dx2*dx2) + (dy2*dy2));
        
        double scale = sep2/sep1;
        
        /*
        estimate translation:
        
        transX = xt0 - xc*scale - (((x0-xc)*scale*math.cos(theta)) 
            + ((y0-yc)*scale*math.sin(theta)))
        
        transY = yt0 - yc*scale - ((-(x0-xc)*scale*math.sin(theta)) 
            + ((y0-yc)*scale*math.cos(theta)))
        */
        double mc = Math.cos(theta);
        double ms = Math.sin(theta);
        
        double tr1X1 = centroidX1*scale + (((set1X1 - centroidX1) * scale*mc) 
            + ((set1Y1 - centroidY1) *scale*ms));
      
        double tr1Y1 = centroidY1*scale + ((-(set1X1 - centroidX1) *scale*ms) 
            + ((set1Y1 - centroidY1) *scale*mc));
        
        double tr1X2 = centroidX1*scale + (((set1X2 - centroidX1) * scale*mc) 
            + ((set1Y2 - centroidY1) *scale*ms));
      
        double tr1Y2 = centroidY1*scale + ((-(set1X2 - centroidX1) *scale*ms) 
            + ((set1Y2 - centroidY1) *scale*mc));
        
        double transX = 0.5 * ((set2X1 - tr1X1) + (set2X2 - tr1X2));

        double transY = 0.5 * ((set2Y1 - tr1Y1) + (set2Y2 - tr1Y2));
        
        TransformationParameters params = new TransformationParameters();
        params.setRotationInRadians((float)theta);
        params.setScale((float)scale);
        params.setTranslationX((float)transX);
        params.setTranslationY((float)transY);
        params.setOriginX((float)centroidX1);
        params.setOriginY((float)centroidY1);
        
        if (debug) {
            log.info("params: " + params.toString());
        }
        
        return params;
    }
    
    /**
     * from a set of transformation parameters params that transform
     * points in reference frame 1 into reference frame 2, create
     * a transformation that can transform points in reference frame 2
     * into reference frame 1.
     * @param params transformation parameters to apply to points in reference
     * frame 1 to put them in reference frame 2
     * @return transformation parameters that can transform points in reference
     * frame 2 into reference frame 1
     */
    TransformationParameters swapReferenceFrames(TransformationParameters 
        params) {
        
        if (params == null) {
            throw new IllegalArgumentException("params cannot be null");
        }
        
        /*

        xr_0 = xc*scale + (((x0-xc)*scale*math.cos(theta)) + ((y0-yc)*scale*math.sin(theta)))

        xt_0 = xr_0 + transX = x1

        yr_0 = yc*scale + ((-(x0-xc)*scale*math.sin(theta)) + ((y0-yc)*scale*math.cos(theta)))

        yt_0 = yr_0 + transY = y1        
        */
     
        double revRot = -1 * params.getRotationInRadians();
        if (revRot < 0) {
            revRot += 2. * Math.PI;
        }
        double revScale = 1. / params.getScale();

        double transformedXC =
            (params.getOriginX()*params.getScale()) + params.getTranslationX();
        
        double transformedYC =
            (params.getOriginY()*params.getScale()) + params.getTranslationY();
        
        TransformationParameters paramsRev = new TransformationParameters();
        paramsRev.setScale((float)revScale);
        paramsRev.setRotationInRadians((float)revRot);
        paramsRev.setTranslationX(0);
        paramsRev.setTranslationY(0);
        paramsRev.setOriginX((float)transformedXC);
        paramsRev.setOriginY((float)transformedYC);
        
        // transform the new origin, then the new translation is what is needed for it to equal the params origin
        Transformer transformer = new Transformer();
        
        double[] revTransformedXYOrigin = transformer.applyTransformation(paramsRev, 
            transformedXC, transformedYC);
        
        double revTransX = params.getOriginX() - revTransformedXYOrigin[0];
        double revTransY = params.getOriginY() - revTransformedXYOrigin[1];
        
        paramsRev.setTranslationX((float)revTransX);
        paramsRev.setTranslationY((float)revTransY);
        
        return paramsRev;
    }
}
