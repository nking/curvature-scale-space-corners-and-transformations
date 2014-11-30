package algorithms.imageProcessing.util;

/**
 *
 * @author nichole
 */
public class AngleUtil {
    
    /**
    calculates the difference of angles between pairs of points in set1 and 
    set 2.  diffX1, diffX2 are the difference
    * in x and y between 2 points in set 1, and diffX2, diffY2 are the 
    * difference for the same matched point pair in set2.

     * @param diffX1
     * @param diffY1
     * @param diffX2
     * @param diffY2
     * @return 
     */
    public double subtract(double diffX1, double diffY1, double diffX2, 
        double diffY2) {
        
        double theta1 = (diffX1 == 0) ? Math.PI/2. : Math.atan(diffY1/diffX1);
        double theta2 = (diffX2 == 0) ? Math.PI/2. : Math.atan(diffY2/diffX2);

        // Q1, Q2, Q3, Q4
        int q1 = 1;
        if ((diffX1 < 0) && (diffY1 < 0)) {
            q1 = 2;
        } else if ((diffX1 < 0) && (diffY1 >= 0)) {
            q1 = 3;
        } else if ((diffX1 >= 0) && (diffY1 >= 0)) {
            q1 = 4;
        }
        int q2 = 1;
        if ((diffX2 < 0) && (diffY2 < 0)) {
            q2 = 2;
        } else if ((diffX2 < 0) && (diffY2 >= 0)) {
            q2 = 3;
        } else if ((diffX2 >= 0) && (diffY2 >= 0)) {
            q2 = 4;
        }
        
        /*
                  +Y
                 270
        QIII      |       QIV
                  |     
                  |
     180-------------------- +X  0, 360
                  |   
                  |      
         QII      |       QI 
                 90
        */
        
        if (q1 == 1) {
            // angle is (-)
            theta1 *= -1;
        } else if (q1 == 2) {
            theta1 = Math.PI - theta1;
        } else if (q1 == 3) {
            // angle is (-)
            theta1 = Math.PI - theta1;
        } else if (q1 == 4) {
            theta1 = 2.*Math.PI - theta1;
        }
        if (diffX1 == 0) {
            if (diffY1 < 0) {
                theta1 = Math.PI/2.;
            } else {
                theta1 = 3.*Math.PI/2.;
            }
        }
        
        if (q2 == 1) {
            // angle is (-)
            theta2 *= -1;
        } else if (q2 == 2) {
            theta2 = Math.PI - theta2;
        } else if (q2 == 3) {
            // angle is (-)
            theta2 = Math.PI - theta2;
        } else if (q2 == 4) {
            theta2 = 2.*Math.PI - theta2;
        }
        if (diffX2 == 0) {
            if (diffY2 < 0) {
                theta2 = Math.PI/2.;
            } else {
                theta2 = 3.*Math.PI/2.;
            }
        }
        
        double t = theta1 - theta2;
        
        while (t < 0) {
            t += 2.*Math.PI;
        }
        while (t > 2.*Math.PI) {
            t -= 2.*Math.PI;
        }
        
        /*
        TODO: correction to a CCW system temporarily
        */
        t = 2.*Math.PI - t;
        
        return t;  
    }
}
