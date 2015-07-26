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
            if (theta1 != 0) {
                theta1 = 2.*Math.PI - theta1;
            }
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
            if (theta2 != 0) {
                theta2 = 2.*Math.PI - theta2;
            }
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
        if (t != 0) {
            t = 2.*Math.PI - t;
        }
        
        return t;  
    }
    
    /**
    calculates the polar theta given x and y w.r.t. origin.  theta increases
    * in value in a counter clockwise direction (CCW).

     * @param x
     * @param y
     * @return 
     */
    public static double polarAngleCCW(double x, double y) {
        
        /*
                  +Y
                 90
        QII       |       QI
                  |     
                  |
     180-------------------- +X  0, 360
                  |   
                  |      
         QIII     |       QIV 
                 270
        */
        
        if (x == 0) {
            if (y >= 0) {
                return Math.PI/2;
            }
            return (3./2.)*Math.PI;
        }
        if (y == 0) {
            if (x > 0) {
                return 0;
            }
            return Math.PI;
        }
        /*
                  +Y
                 90
        QII       |       QI
                  |     
                  |
     180-------------------- +X  0, 360
                  |   
                  |      
         QIII     |       QIV 
                 270
        */
        double theta = Math.atan(y/x);

        // Q1, Q2, Q3, Q4
        int q = 1;
        if ((x < 0) && (y >= 0)) {
            q = 2;
        } else if ((x < 0) && (y < 0)) {
            q = 3;
        } else if ((x >= 0) && (y < 0)) {
            q = 4;
        }
        
        if (q == 2 || q == 3) {
            theta += Math.PI;
        } else if (q == 4) {
            theta = 2*Math.PI + theta;
        }
        
        return theta;  
    }
    
    public static float getAngleDifference(float rotDegrees0, float rotDegrees1) {
         /*
         I  |  0
        ---------
         II | III
        */
        int q0 = 0;
        if (rotDegrees0 >= 270) {
            q0 = 3;
        } else if (rotDegrees0 >= 180) {
            q0 = 2;
        } else if (rotDegrees0 >= 90) {
            q0 = 1;
        }
        int q1 = 0;
        if (rotDegrees1 >= 270) {
            q1 = 3;
        } else if (rotDegrees1 >= 180) {
            q1 = 2;
        } else if (rotDegrees1 >= 90) {
            q1 = 1;
        }

        /*
         I  |  0
        ---------
         II | III
        */
        float angleDiff = -1;
        if (q0 == 0){
            if (q1 == 0) {
                if (rotDegrees0 > rotDegrees1) {
                    angleDiff = rotDegrees0 - rotDegrees1;
                } else {
                    angleDiff = rotDegrees1 - rotDegrees0;
                }
            } else if (q1 == 1) {
                angleDiff = (rotDegrees1 - rotDegrees0);
            } else if (q1 == 2) {
                float diff = rotDegrees1 - rotDegrees0;
                if (diff > 180) {
                    diff = 360 - diff;
                }
                angleDiff = diff;
            } else {
                angleDiff = Math.abs(360 - rotDegrees1 + rotDegrees0);
            }
        } else if (q0 == 1) {
            /*
             I  |  0
             ---------
             II | III
             */
            if (q1 == 0) {
                angleDiff = (rotDegrees1 - rotDegrees0);
            } else if (q1 == 1) {
                if (rotDegrees0 > rotDegrees1) {
                    angleDiff = rotDegrees0 - rotDegrees1;
                } else {
                    angleDiff = rotDegrees1 - rotDegrees0;
                }
            } else if (q1 == 2) {
                angleDiff = (rotDegrees1 - rotDegrees0);
            } else {
                float diff = rotDegrees1 - rotDegrees0;
                if (diff > 180) {
                    diff = 360 - diff;
                }
                angleDiff = diff;
            }
        } else if (q0 == 2) {
            /*
             I  |  0
             ---------
             II | III
             */
            if (q1 == 0) {
                float diff = rotDegrees1 - rotDegrees0;
                if (diff > 180) {
                    diff = 360 - diff;
                }
                angleDiff = diff;
            } else if (q1 == 1) {
                angleDiff = (rotDegrees0 - rotDegrees1);
            } else if (q1 == 2) {
                if (rotDegrees0 > rotDegrees1) {
                    angleDiff = rotDegrees0 - rotDegrees1;
                } else {
                    angleDiff = rotDegrees1 - rotDegrees0;
                }
            } else {
                angleDiff = (rotDegrees1 - rotDegrees0);
            }
        } else if (q0 == 3) {
            /*
             I  |  0
             ---------
             II | III
             */
            if (q1 == 0) {
                angleDiff = (360 - rotDegrees0 + rotDegrees1);
            } else if (q1 == 1) {
                float diff = (rotDegrees0 - rotDegrees1);
                if (diff > 180) {
                    diff = 360 - diff;
                }
                angleDiff = diff;
            } else if (q1 == 2) {
                angleDiff = (rotDegrees0 - rotDegrees1);
            } else {
                if (rotDegrees0 > rotDegrees1) {
                    angleDiff = rotDegrees0 - rotDegrees1;
                } else {
                    angleDiff = rotDegrees1 - rotDegrees0;
                }
            }
        }

        if (angleDiff > 359) {
            angleDiff = angleDiff - 360;
        }

        return angleDiff;
    }

}
