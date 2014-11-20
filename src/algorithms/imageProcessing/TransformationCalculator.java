package algorithms.imageProcessing;

import algorithms.util.PairIntArray;
import java.util.logging.Logger;

/**
 *
 * @author nichole
 */
public class TransformationCalculator {
    
    protected transient Logger log = Logger.getLogger(
        TransformationCalculator.class.getName());
    
    private boolean debug = false;
    
    public void useDebugMode() {
        debug = true;
    }
    
    /**
     * coordinate transformations from image 1 to image 2 are calculated from
     * matching lists of x, y coordinates.
     * 
     * Note that if this is used for contours, only one contour at a time
     * should be passed in because the centroids are w.r.t. one contour.
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
     * @param matchedXY1
     * @param matchedXY2
     * @param centroidX1
     * @param centroidY1
     * @param centroidX2
     * @param centroidY2
     * @return 
     */
    public TransformationParameters calulateEuclidean(
        PairIntArray matchedXY1, PairIntArray matchedXY2, 
        double centroidX1, double centroidY1,
        double centroidX2, double centroidY2) {
        
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
        
        int n = matchedXY1.getN();
        
        float[] weights1 = new float[n];
        float[] weights2 = new float[n];
        
        float invN = 1.f/(float)n;
        
        // make even weights
        for (int i = 0; i < n; i++) {
            weights1[i] = invN;
            weights2[i] = invN;
        }
        
        return calulateEuclidean(matchedXY1, weights1, matchedXY2, weights2,
            centroidX1, centroidY1, centroidX2, centroidY2);
    }
    
    /**
     * coordinate transformations from image 1 to image 2 are calculated from
     * matching lists of x, y coordinates.
     *
     * Note that if this is used for contours, only one contour at a time
     * should be passed in because the centroids are w.r.t. one contour.
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
     * @param matchedXY1
     * @param weights1
     * @param matchedXY2
     * @param weights2
     * @param centroidX1
     * @param centroidY1
     * @param centroidX2
     * @param centroidY2
     * @return 
     */
    public TransformationParameters calulateEuclidean(
        PairIntArray matchedXY1, float[] weights1,
        PairIntArray matchedXY2, float[] weights2, 
        double centroidX1, double centroidY1,
        double centroidX2, double centroidY2
        ) {
        
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
        
        if (debug) {
            log.info("centroidX1=" + centroidX1 + " centroidY1=" + centroidY1
                + "\ncentroidX2=" + centroidX2 + " centroidY2=" + centroidY2
            );
        }
        
        /*
        SCALE:
                      (sum set 1:((x_1 - centroid_x)^2 + (x_2 - centroid_x)^2 )
        scale_x = sqrt(-------------------------------------------------------)
                      (sum set 2:((x_1 - centroid_x)^2 + (x_2 - centroid_x)^2 )
        */     
        double sumDiffX1 = 0;
        double sumDiffY1 = 0;
        double sumDiffX2 = 0;
        double sumDiffY2 = 0;
        for (int i = 0; i < matchedXY1.getN(); i++) {
            
            float w1 = weights1[i];
            float w2 = weights2[i];
            
            double x1 = matchedXY1.getX(i) - centroidX1;
            double y1 = matchedXY1.getY(i) - centroidY1;
            
            double x2 = matchedXY2.getX(i) - centroidX2;
            double y2 = matchedXY2.getY(i) - centroidY2;
            
            sumDiffX1 += Math.pow(x1*w1, 2);
            sumDiffY1 += Math.pow(y1*w1, 2);
            
            sumDiffX2 += Math.pow(x2*w2, 2);
            sumDiffY2 += Math.pow(y2*w2, 2);
        }
                
        float scaleX = (float) Math.sqrt(sumDiffX2/sumDiffX1);
        float scaleY = (float) Math.sqrt(sumDiffY2/sumDiffY1);
        float scale = (scaleX + scaleY)/2.f;
      
        double theta2Minus1WeightedSum = 0;
        
        /*
        positive Y is down 
        positive X is right
        positive theta starts from Y=0, X>=0 and proceeds CW
                270
           Q4    |   Q1
                 |
          180--------- 0   +X
                 |   
           Q3    |   Q2
                 90
                 +Y
        */
        
        StringBuilder sb = new StringBuilder();
        
        if (debug) {
            sb.append("\n  x1    y1       x2    y2   len1  len2 theta1 theta2 (x,y) (x,y) \n");
        }
        
        for (int i = 0; i < matchedXY1.getN(); i++) {
            
            float w1 = weights1[i];
            float w2 = weights2[i];
            
            double x1 = scale*(matchedXY1.getX(i) - centroidX1);
            double y1 = scale*(matchedXY1.getY(i) - centroidY1);
            double len1 = Math.sqrt((x1*x1) + (y1*y1));
            
            double x2 = matchedXY2.getX(i) - centroidX2;
            double y2 = matchedXY2.getY(i) - centroidY2;            
            double len2 = Math.sqrt((x2*x2) + (y2*y2));
            
            double aDivH1 = x1/len1;
            double theta1 = Math.acos(aDivH1) * 180./Math.PI;
            
            double aDivH2 = x2/len2;
            double theta2 = Math.acos(aDivH2) * 180./Math.PI;
          
            if ((y1 < 0) && (x1 < 0)) {
                theta1 = 360 - theta1;
            } else if ((y1 < 0) && (x1 >= 0)) {
                theta1 = 360 - theta1;
            } 
            if ((y2 < 0) && (x2 < 0)) {
                theta2 = 360 - theta2;
            } else if ((y2 < 0) && (x2 >= 0)) {
                theta2 = 360 - theta2;
            }
            
            if (debug) {
                
                sb.append(String.format(
                    "%5.1f %5.1f    %5.1f %5.1f %5.1f %5.1f %5.0f %5.0f (%5d,%5d) (%5d,%5d)\n", 
                    x1, y1, x2, y2, len1, len2, theta1, theta2, 
                    matchedXY1.getX(i), matchedXY1.getY(i), 
                    matchedXY2.getX(i), matchedXY2.getY(i)));
            }

            double thetaDiff;
            if (theta2 >= theta1) {
                thetaDiff = theta2 - theta1;
            } else {
                if ((theta2 < 90) && (theta1 > 270)) {
                    thetaDiff = (360 - theta1) + theta2;
                } else if ((theta1 - theta2) < 10) {
                    thetaDiff = theta1 - theta2;
                } else {
                    thetaDiff = (360 - theta1) + theta2;
                }
            }
            
            theta2Minus1WeightedSum += Math.abs((w1+w2)*thetaDiff);
        }
        
        if (debug) {
            log.info(sb.toString());
        }
        
        theta2Minus1WeightedSum /= 2;
        
        double rotationInDegrees = -1*theta2Minus1WeightedSum;
        
        /*
        solve for an x translation and y translation which can be applied
        to a data set from image 1 after scaled and rotated.
        
        xt = centroidX2 + 
             (x-centroidX1)*scale*cos(theta) + (y-centroidY1)*scale*sin(theta)
           = centroidX2 
             + x*scale*cos(theta) + y*scale*sin(theta)
             - centroidX1*scale*cos(theta) - centroidY1*scale*sin(theta)
        
        So, one can define the portion w/o x and y as constant for the image:
           translationX = centroidX2 - centroidX1*scale*cos(theta) 
               - centroidY1*scale*sin(theta)
        
        yt = centroidY2 + 
             -1*(x-centroidX1)*scale*sin(theta) + (y-centroidY1)*scale*cos(theta)
           = centroidY2 + 
             - x*scale*sin(theta) + y*scale*cos(theta)
             + centroidX1*scale*sin(theta) - centroidY1*scale*cos(theta)
        
        So, one can define
           translationY = centroidY2 - centroidX1*scale*sin(theta) 
               - centroidY1*scale*cos(theta)
        */
        float mc = (float)Math.cos(rotationInDegrees*Math.PI/180.);
        float ms = (float)Math.sin(rotationInDegrees*Math.PI/180.);
        
        log.fine("mc=" + mc + " ms=" + ms);
        
        float translationX = (float)(centroidX2 
             - (centroidX1*scale*mc) - (centroidY1*scale*ms));
        float translationY = (float)(centroidY2 
            + (centroidX1*scale*ms) - (centroidY1*scale*mc));
   
        log.fine("(scale*mc + scale*ms)=" + (scale*mc + scale*ms)
            + " (centroidX1*(scale*mc + scale*ms))=" 
            + (centroidX1*(scale*mc + scale*ms)));
        log.fine("(scale*ms - scale*mc)=" + (scale*ms - scale*mc)
            + " (centroidX1*(scale*ms - scale*mc))=" 
            + (centroidX1*(scale*ms - scale*mc)));
        
        log.fine("\nscaleX=" + scaleX + " scaleY=" + scaleY + " (" + scale + ")"
            + "\nrotation=" + rotationInDegrees
            + "\ntranslationX=" + translationX + " translationY=" + translationY 
            );
    
        TransformationParameters params = new TransformationParameters();
        params.setScale(scale);
        params.setTranslationX(translationX);
        params.setTranslationY(translationY);
        params.setRotationInDegrees((float)rotationInDegrees);
        
        return params;
    }
}
