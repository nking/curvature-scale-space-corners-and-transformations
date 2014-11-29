package algorithms.imageProcessing.util;

/**
 *
 * @author nichole
 */
public class AngleUtil {
    
    /**
    calculates the difference of angles as a dependency upon which unit 
    circle quadrant they are in.

    The unit circle quadrants are defined as:
        positive Y is up.
        positive X is right.
        theta is negative in QII and QIV.
        theta is given by Math.atan with respect to the X-axis (Y=0).
                 +Y
        QIII      |       QIV
                  |     
                  |
        -------------------- +X
                  |   
                  |      
         QII      |       QI     
     *           
     * 
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

        /*
        TODO: there are changes to correct for a CCW system and uncorrect below
        that need to be fixed consistenly throughout code.  For now, focussing 
        on tests to assert the correct answers as a total result.
        */
        
        theta1 *= -1;
        theta2 *= -1;
        
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
        
        if (diffX2 == 0) {
            // the result is 90 degrees, but the sign needs to be adjusted to
            // be usable for the difference of q1 and q2 angles
            if ((q2 == 1) || (q2 == 3)) {
                theta2 = Math.PI/2.;
            } else {
                theta2 = -1*Math.PI/2.;
            }
        }
        if (diffX1 == 0) {
            if ((q1 == 1) || (q1 == 3)) {
                theta1 = Math.PI/2.;
            } else {
                theta1 = -1*Math.PI/2.;
            }
        }
            
        double t = theta1 - theta2;
        if (q1 == 1) {
            if (q2 == 1) {
                if (theta1 < theta2) {
                    t *= -1;
                }
            } else if (q2 == 2) {
                t = Math.PI + theta1 - theta2;
            } else if (q2 == 3) {
                if (theta1 < theta2) {
                    t = Math.PI - theta1 + theta2;
                } else {
                    t = Math.PI + theta1 - theta2;
                }
            } else if (q2 == 4) {
                //t = theta1 - theta2;
            }
        } else if (q1 == 2) {
            if (q2 == 1) {
                t = Math.PI + theta1 - theta2;
            } else if (q2 == 2) {
                if (theta1 < theta2) {
                    t = -theta2 + 2 * Math.PI + theta1;
                }
            } else if (q2 == 3) {
                t = 2 * Math.PI + theta1 - theta2;
            } else if (q2 == 4) {
                if (theta1 < -45. * Math.PI / 180.) {
                    t = Math.PI + theta1 - theta2;
                } else if (theta2 > -45. * Math.PI / 180.) {
                    t = Math.PI + theta1 - theta2;
                } else {
                    t = Math.PI + theta1 - theta2;
                }
            }
        } else if (q1 == 3) {
            if (q2 == 1) {
                //if (theta1 < theta2) {
                    t = Math.PI + theta1 - theta2;
                /*} else {
                    t = Math.PI - theta1 + theta2;
                }*/
            } else if (q2 == 2) {
                t = theta1 - theta2;
            } else if (q2 == 3) {
                if (theta1 < theta2) {
                    t = theta1 - theta2;
                }
            } else if (q2 == 4) {                    
                t = Math.PI + theta1 - theta2; 
            }
        } else if (q1 == 4) {
            if (q2 == 1) {
                t = 2 * Math.PI + theta1 - theta2;
            } else if (q2 == 2) {
                /*if (theta1 > -45. * Math.PI / 180.) {
                    t = Math.PI + theta1 - theta2;
                } else if (theta2 < -45. * Math.PI / 180.) {*/
                    t = Math.PI + theta1 - theta2;
                /*} else {
                    t = Math.PI - theta1 + theta2;
                }*/               
            } else if (q2 == 3) {
                t = Math.PI + theta1 - theta2;
                //t = theta1 + theta2;
            } else if (q2 == 4) {
                if (theta1 < theta2) {
                    t = theta1 - theta2;
                } else {
                    
                }
            }
        }

        // reverse the direction to CW
        t *= -1;

        if (t < 0) {
            t += 2*Math.PI;
        }
        
        String str = String.format("%d %d %f %f %f (%f, %f, %f, %f)", 
            q1, q2, theta1, theta2, t, diffX1, diffY1, diffX2, diffY2);

        System.out.println(str.toString());
        
        return t;  
    }
}
