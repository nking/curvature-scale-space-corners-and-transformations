package algorithms.imageProcessing.util;

import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.Image;
import algorithms.misc.MiscMath;
import algorithms.util.PairInt;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

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
    calculates the polar theta in radians given x and y w.r.t. origin.  theta increases
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
    
    public static float getAngleAverageInDegrees(float rotDegrees0, float rotDegrees1) {
       
        double angleAvg = getAngleAverage(rotDegrees0, rotDegrees1, false);
        
        return (float) angleAvg;
    }
    
    public static double getAngleAverageInRadians(double rotation0, double rotation1) {
       
        double angleAvg = getAngleAverage(rotation0, rotation1, true);
        
        return angleAvg;
    }
    
    protected static int getClockwiseQuadrantForDegrees(double rotationInDegrees) {
        /*
          III | IV
          ---------
          II  |  I
        */
        int q = 1;
        if (rotationInDegrees >= 270) {
            q = 4;
        } else if (rotationInDegrees >= 180) {
            q = 3;
        } else if (rotationInDegrees >= 90) {
            q = 2;
        }
        return q;
    }
    
    protected static int getClockwiseQuadrantForRadians(double rotation) {
        /*
          III | IV
          ---------
          II  |  I
        */
        
        if (rotation >= 2*Math.PI) {
            rotation = rotation % (2.*Math.PI);
        } else if (rotation < 0) {
            while (rotation < 0) {
                rotation += 2.*Math.PI;
            }
        }
        
        int q = 1;
        if (rotation >= 3.*Math.PI/2.) {
            q = 4;
        } else if (rotation >= Math.PI) {
            q = 3;
        } else if (rotation >= Math.PI/2.) {
            q = 2;
        }
        return q;
    }
    
    /**
     * given twoPi in degrees or in radians, return the angle average.
     * @param rot0
     * @param rot1
     * @param useRadians
     * @return 
     */
    protected static double getAngleAverage(double rot0, double rot1, boolean useRadians) {
        
         /*
                  270
               III | IV
          180  --------- 0
               II  |  I
                   90
        */
        double angleAvg = calcAngleAddition(rot0, rot1, useRadians)/2.;
      
        return angleAvg;
    }
    
    /**
     * given twoPi in degrees or in radians, return the angle sum corrected to
     * the larger angle frame, e.g. 0 + 350 = 710.
     * @param rot0
     * @param rot1
     * @param useRadians
     * @return 
     */
    public static double calcAngleAddition(double rot0, double rot1, 
        boolean useRadians) {
        
        if (rot0 < 0 || rot1 < 0) {
            throw new IllegalArgumentException(
                "rot0 and rot1 cannot be negative numbers");
        }
        
        /*
                  270
               III | IV
          180  --------- 0
               II  |  I
                   90
        */
         
        /*
        logic below is for rot0 and rot1 between 0 and 360 so temporarily
        reduce angles to those ranges then add back the corrections.
        */
        double rot0Orig = rot0;
        double rot1Orig = rot1;
        
        rot0 = rot0 % 360;
        rot1 = rot1 % 360;
        
        double twoPI;
        int q0, q1;
        if (useRadians) {
            twoPI = 2. * Math.PI;
            q0 = getClockwiseQuadrantForRadians(rot0);
            q1 = getClockwiseQuadrantForRadians(rot1);
        } else {
            twoPI = 360.;
            q0 = getClockwiseQuadrantForDegrees(rot0);
            q1 = getClockwiseQuadrantForDegrees(rot1);
        }
       
        /*
                  270
               III | IV
          180  --------- 0
               II  |  I
                   90
        */
        double angleSum = -1;
        if (q0 == 1) {
            if (q1 == 1 || q1 == 2) {
                angleSum = (rot0 + rot1);
            } else if (q1 == 3) {
                double diff = rot1 - rot0;
                if (diff > (twoPI/2.)) {
                    angleSum = (twoPI + rot0 + rot1);
                } else {
                    angleSum = (rot0 + rot1);
                }
            } else {
                angleSum = (twoPI + rot0 + rot1);
            }
        } else if (q0 == 2) {
            /*
                  270
               III | IV
          180  --------- 0
               II  |  I
                   90
            */
            if (q1 == 1 || q1 == 2 || q1 == 3) {
                angleSum = (rot0 + rot1);
            } else if (q1 == 4) {
                double diff = rot1 - rot0;
                if (diff > (twoPI/2.)) {
                    angleSum = (twoPI + rot0 + rot1);
                } else {
                    angleSum = (rot0 + rot1);
                }
            }
        } else if (q0 == 3) {
            /*
                  270
               III | IV
          180  --------- 0
               II  |  I
                   90
            */
            if (q1 == 1) {
                double diff = rot0 - rot1;
                if (diff > (twoPI/2.)) {
                    angleSum = (twoPI + rot1 + rot0);
                } else {
                    angleSum = (rot0 + rot1);
                }
            } else {
                angleSum = (rot0 + rot1);
            }
        } else if (q0 == 4) {
            /*
                  270
               III | IV
          180  --------- 0
               II  |  I
                   90
            */
            if (q1 == 1) {
                angleSum = (twoPI + rot0 + rot1);
            } else if (q1 == 2) {
                double diff = (rot0 - rot1);
                if (diff > (twoPI/2.)) {
                    angleSum = (twoPI + rot0 + rot1);
                } else {
                    angleSum = (rot0 + rot1);
                }
            } else {
                angleSum = (rot0 + rot1);
            }
        }
        
        // add back any of the cycles in original values
        if (rot0 < rot0Orig) {
            angleSum += (rot0Orig - rot0);
        }
        if (rot1 < rot1Orig) {
            angleSum += (rot1Orig - rot1);
        }

        return angleSum;
    }
    
    protected static int getClockwiseQuadrant(float rotationInDegrees) {
        /*
          III | IV
          ---------
          II  |  I
        */
        
        if (rotationInDegrees >= 360) {
            rotationInDegrees = rotationInDegrees % 360;
        } else if (rotationInDegrees < 0) {
            while (rotationInDegrees < 0) {
                rotationInDegrees += 360;
            }
        }
        
        int q = 1;
        if (rotationInDegrees >= 270) {
            q = 4;
        } else if (rotationInDegrees >= 180) {
            q = 3;
        } else if (rotationInDegrees >= 90) {
            q = 2;
        }
        return q;
    }
}
