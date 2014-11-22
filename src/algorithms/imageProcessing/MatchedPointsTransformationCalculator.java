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
        
        //TODO: this could be improved by making pairs out of points with 
        // furthest x and y from each other
        // or by using all combinations weighted by deltax and deltay.
        // could use bipartite max weight algorithm such as 
        // Hungarian or a max flow algorithm.
        
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
            
            double thetaim1 = (diffX1 == 0) ? Math.PI/2. :
                Math.atan(diffY1/diffX1);
            double thetaim2 = (diffX2 == 0) ? Math.PI/2. : 
                Math.atan(diffY2/diffX2);
            
            thetaim1 *= -1;
            thetaim2 *= -1;
            
            // Q1, Q2, Q3, Q4
            int qim1 = 1;
            if ((diffX1 < 0) && (diffY1 < 0)) {
                qim1 = 2;
            } else if ((diffX1 < 0) && (diffY1 >= 0)) {
                qim1 = 3;
            } else if ((diffX1 >= 0) && (diffY1 >= 0)) {
                qim1 = 4;
            }
            int qim2 = 1;
            if ((diffX2 < 0) && (diffY2 < 0)) {
                qim2 = 2;
            } else if ((diffX2 < 0) && (diffY2 >= 0)) {
                qim2 = 3;
            } else if ((diffX2 >= 0) && (diffY2 >= 0)) {
                qim2 = 4;
            }
            
            // interpretation of subtracting angles
            
            /*
            note that they were made with positive as CCW, so 
            a change is made below after the blocks
            */
            
             /*
            defaults fine for:
            1:1 1.406,1.2375
            4:4 -.54, -.987
            3:3 .588 .423
            2:1 -1.55, 1.24
            4:4 -0.785,-.93
            3:3 .6,.327
            ------
            1:1 1.41,.54
            4:3 -.54, 1.48
            3:2 .588,-.48
            2:1 -1.55,.303
            4:4 -0.79, -1.57
            3:2 0.6,-.4
            --
            GOOD TOO: 1:3 1.4, 0.785
            BAD:      4:2 -.54, -1.45 result should be ~-4, not -2.2
                      3:4 0.59, -0.09
                      2:3 -1.55,.8067
                      4:2 -0.785,-1.48
                      3:4 .6, -.0677
             
            */
            
            double t = thetaim1 - thetaim2;
            if ((qim1 == 1) && (qim2 == 2)) {
                t = Math.PI + thetaim1 - thetaim2;
            } else if ((qim1 == 1) && (qim2 == 3)) {
                if (thetaim1 > 45.*Math.PI/180.) {
                    t = Math.PI + thetaim1 - thetaim2;
                } else {
                    t = Math.PI + thetaim2 - thetaim1;
                }
            } else if ((qim1 == 1) && (qim2 == 4)) {
                //t = thetaim1 - thetaim2;
            } else if ((qim1 == 1) && (qim2 == 1)) {
                if (thetaim1 < thetaim2) {
                    t *= -1;
                }
            } else if ((qim1 == 2) && (qim2 == 1)) {
                t = Math.PI + thetaim1 - thetaim2;
            } else if ((qim1 == 2) && (qim2 == 2)) {
                if (thetaim1 < thetaim2) {
                    //t *= -1;
                    t = -thetaim2 + 2 * Math.PI + thetaim1;
                }
            } else if ((qim1 == 2) && (qim2 == 3)) {
                t = 2 * Math.PI +thetaim1 - thetaim2;
            } else if ((qim1 == 2) && (qim2 == 4)) {
                if (thetaim1 < -45.*Math.PI/180.) {
                    t = Math.PI + thetaim1 - thetaim2;
                } else {
                    t = Math.PI - thetaim1 + thetaim2;
                }
            } else if ((qim1 == 3) && (qim2 == 4)) {
                t = Math.PI + thetaim1 - thetaim2;
               
            } else if ((qim1 == 3) && (qim2 == 1)) {
                if (thetaim1 < 45.*Math.PI/180.) {
                    t = Math.PI + thetaim1 - thetaim2;
                } else {
                    t = Math.PI - thetaim1 + thetaim2;
                }
            } else if ((qim1 == 3) && (qim2 == 2)) {
                t = thetaim1 - thetaim2;
            } else if ((qim1 == 3) && (qim2 == 3)) {
                if (thetaim1 < thetaim2) {
                    t = thetaim2 - 2 * Math.PI - thetaim1;
                } 
            } else if ((qim1 == 4) && (qim2 == 1)) {
                t = 2*Math.PI + thetaim1 - thetaim2;
            } else if ((qim1 == 4) && (qim2 == 2)) {
                if (thetaim1 > -45.*Math.PI/180.) {
                    t = Math.PI + thetaim1 - thetaim2; 
                } else if (thetaim2 < -45.*Math.PI/180.) {
                    t = Math.PI + thetaim1 - thetaim2;
                } else {
                    t = Math.PI - thetaim1 + thetaim2;
                }
            } else if ((qim1 == 4) && (qim2 == 3)) {
                t = Math.PI + thetaim1 - thetaim2;
            } else if ((qim1 == 4) && (qim2 == 4)) {
                if (thetaim1 < -45. * Math.PI/180.) {
                    t = thetaim2 - 2 * Math.PI - thetaim1;
                } else if (thetaim1 > 45. * Math.PI/180.) {
                    t = thetaim2 - 2 * Math.PI - thetaim1;
                }
            }
            
            // reverse the direction to CW
            t *= -1;
            
            if (t < 0) {
                t += 2*Math.PI;
            }
            
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
        for (int i = 0; i < matchedXY1.getN(); i++) {
            
            if (Math.abs(scales[i] - avgScale) <= stDevScale) {
                scaleSum += scales[i];
                sCount++;
            }
            if (Math.abs(thetas[i] - avgTheta) <= stDevTheta) {
                rotSum += thetas[i];
                rCount++;
            }
log.info("rot=" + thetas[i] + " stDevTheta=" + stDevTheta
+ " (theta-avg)=" + (thetas[i] - avgTheta) 
+ " w1=" + weights1[i] + " w2=" + weights2[i]);
log.info("scl=" + scales[i] + " stDevScale=" + stDevScale
+ " (scale-avg)=" + (scales[i] - avgScale));
            
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
