package algorithms.imageProcessing;

import algorithms.util.PairIntArray;
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
     * positive Y is down 
       positive X is right
       positive theta starts from Y=0, X>=0 and proceeds CW
                270
                 |     
                 |
          180--------- 0   +X
                 |   
                 |   
                 90
                 +Y
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

            rotation = theta_iamge1 - theta_image2 = 25 degrees
        */
                
        //TODO: this could be improved by comparing furthest pairs
        // or all combinations of pairs
        double rotSum = 0;
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
            double xdiff = (x1im1 - x0im1);
            double thetaim1 = (xdiff == 0) ? Math.PI/2. :
                Math.atan((y1im1 - y0im1)/(x1im1 - x0im1));
            xdiff = (x1im2 - x0im2);
            double thetaim2 = (xdiff == 0) ? Math.PI/2. : 
                Math.atan((y1im2 - y0im2)/(x1im2 - x0im2));
            double rot = thetaim2 - thetaim1;
            rotSum += rot;
            
            double lenim1 = Math.sqrt(Math.pow(x1im1 - x0im1, 2) 
                + Math.pow(y1im1 - y0im1, 2));
            double lenim2 = Math.sqrt(Math.pow(x1im2 - x0im2, 2) 
                + Math.pow(y1im2 - y0im2, 2));
            double scl = lenim2/lenim1;
            scaleSum += scl;
        }
        
        double theRotation = rotSum/(double)matchedXY1.getN();
        double theScale = scaleSum/(double)matchedXY1.getN();
        
        log.info("given scale=" + scale + " found scale=" + theScale);
        log.info("rotation = " + theRotation);

        if (Math.abs(theScale - scale) > scale*0.1) {
            log.warning("the differences in estimated scale and given scale are"
                + " large.  this can happen if the given scale value was " 
                + " determined from contour matching."
                + " the estimate here uses pairs of points which may"
                + " be close to one another.  choosing the scale given to "
                + " the method and continuing.");
            theScale = scale;
        }
        
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
    
}
