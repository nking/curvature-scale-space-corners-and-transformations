package algorithms.imageProcessing.util;

import algorithms.CountingSort;
import algorithms.MergeSort;
import algorithms.MultiArrayMergeSort;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.Image;
import algorithms.misc.MiscMath;
import algorithms.util.PairInt;
import java.util.Arrays;
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
    
    protected static double[] correctForQuadrants(double rot0, double rot1, 
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
        
        if (q0 == 1) {
            if (q1 == 3) {
                double diff = rot1 - rot0;
                if (diff > (twoPI/2.)) {
                    rot0 += twoPI;
                }
            } else if (q1 == 4) {
                rot0 += twoPI;
            }
        } else if (q0 == 2) {
            /*
                  270
               III | IV
          180  --------- 0
               II  |  I
                   90
            */
            if (q1 == 4) {
                double diff = rot1 - rot0;
                if (diff > (twoPI/2.)) {
                    rot0 += twoPI;
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
                    rot1 += twoPI;
                }
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
                rot1 += twoPI;
            } else if (q1 == 2) {
                double diff = (rot0 - rot1);
                if (diff > (twoPI/2.)) {
                    rot1 += twoPI;
                }
            }
        }
        return new double[]{rot0, rot1};
    }
    
    public static float getAngleDifference(float rotDegrees0, float rotDegrees1) {
        
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
        double rot0Orig = rotDegrees0;
        double rot1Orig = rotDegrees1;
        
        double[] rotations = correctForQuadrants(rotDegrees0,rotDegrees1, 
            false);
        
        return (float)(rotations[0] - rotations[1]);
       
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
        
        double[] rotations = correctForQuadrants(rot0, rot1, useRadians);
        
        double angleSum = (rotations[0] + rotations[1]);
        
        // add back any of the cycles in original values
        if (rot0 < rot0Orig) {
            angleSum += (rot0Orig - rot0);
        }
        if (rot1 < rot1Orig) {
            angleSum += (rot1Orig - rot1);
        }

        return angleSum;
    }
    
    /**
     * given angles in degrees or in radians, return the angle sum corrected to
     * the larger angle frame, e.g. 0 + 350 = 710.  It always chooses the
     * smallest difference in quadrant space 
     * between the two angles (adding a cycle if needed) before adding.
     * @param rot0
     * @param rot1
     * @param useRadians
     * @param outputQuadrantCorrected quadrant corrected values of 
     * [rotationInDegrees0, rotationInDegrees1]
     * @return 
     */
    public static double calcAngleAddition(double rot0, double rot1, 
        boolean useRadians, double[] outputQuadrantCorrected) {
        
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
        
        outputQuadrantCorrected[0] = rot0;
        outputQuadrantCorrected[1] = rot1;
        
        // prefer to keep 360 instead of 0 because some of the invokers of this
        // method might prefer a fraction times results to be > 0
        if (rot0 != 360) {
            rot0 = rot0 % 360;
        }
        if (rot1 != 360) {
            rot1 = rot1 % 360;
        }
        
        double twoPI;
        if (useRadians) {
            twoPI = 2. * Math.PI;
        } else {
            twoPI = 360.;
        }
       
        double[] rotations = correctForQuadrants(rot0, rot1, useRadians);
        outputQuadrantCorrected[0] = rotations[0];
        outputQuadrantCorrected[1] = rotations[1];
        rot0 = rotations[0];
        rot1 = rotations[1];
        
        double angleSum = (rot0 + rot1);
        
        // add back any of the cycles in original values if larger than 360
        if ((rot0 > twoPI) && (rot0 < rot0Orig)) {
            double diff = (rot0Orig - rot0);
            angleSum += diff;
            outputQuadrantCorrected[0] += diff;
        }
        if ((rot1 > twoPI) && (rot1 < rot1Orig)) {
            double diff = (rot1Orig - rot1);
            angleSum += diff;
            outputQuadrantCorrected[1] += diff;
        }

        return angleSum;
    }
     
    /**
     * calculate the angular average of rotDegrees0 and rotDegrees1 with
     * corrections for quadrants if necessary.  For example, 0 averaged with
     * 350 is 355.  Any changes in the given angles for quadrant math
     * are seen in the populated outputQuadrantCorrected.
     * @param outputQuadrantCorrected quadrant corrected values of 
     * [rotationInDegrees0, rotationInDegrees1]
     * @return 
    */
    protected static double getAngleAverageInDegrees(
        double rotationInDegrees0, double rotationInDegrees1, 
        double[] outputQuadrantCorrected) {
        
        boolean useRadians = false;
        
        double sum = calcAngleAddition(rotationInDegrees0, 
            rotationInDegrees1, useRadians, outputQuadrantCorrected);
        
        return sum/2.;
    }

    /**
     Runtime complexity is O(N*lg2N) for sort plus O(N) for pair angle averages
     plus O(N) for pair average corrections to total average.
     * Note, may change to use CountingSort for N > 80 and max(angles) ~ 360
     * in the future.
                       
     * @param angles angles in units of degrees.  Note that this array
     * is modified by use here so pass a copy if it should not be.
     * @param lastIndex last index to use in array angles, inclusive
     * @return 
     */
    public static float calculateAverageWithQuadrantCorrections(int[] angles,
        int lastIndex) {
        
        /*
         because these are angles from 0 to 360, need to add them
         while considering their quadrants.

         e.g. (0 + 350 + 340)/3. should equal 350, but if
         the quadrants are not used, the result is 130.

         Can add in pairs, but need to correct the result.

         averaging 4 numbers:
             (a + b + c + d)/4. = (a/4) + (b/4) + (c/4) + (d/4)

         averaging 4 numbers through 3 pair averages:
             1) (a/2) + (b/2)
             2) ((a/2) + (b/2))/2 + (c/2) = (a/4) + (b/4) + (c/2)
             3) ((a/4) + (b/4) + (c/2))/2 + (d/2)
                = (a/8) + (b/8) + (c/4) + (d/2)

         and that needs to be corrected to (a/4) + (b/4) + (c/4) + (d/4)

         for n=5, the successive pair averages is
             1) (a/2) + (b/2)
             2) ((a/2) + (b/2))/2 + (c/2) = (a/4) + (b/4) + (c/2)
             3) ((a/4) + (b/4) + (c/2))/2 + (d/2)
                = (a/8) + (b/8) + (c/4) + (d/2) 
             4) ((a/8) + (b/8) + (c/4) + (d/2))/2 + (e/2)
                = (a/16) + (b/16) + (c/8) + (d/4) + (e/2)
                   0        1         2       3       4
         and that needs to be corrected to (a/5) + (b/5) + (c/5) + (d/5) + (e/5)

         for all but the first number:
            v[i]                 v[i]
         ----------  +  v[i]*x  = ---
         (1<<(n-i))                n
                
         for first variable, can add the correction:
             total += v[0] * ((1<<(n-1)) - n)/(n*(1<<(n-1)))

         then the correction afterwards to add:
             from i=1 to i < n
                 total += v[i] * ((1<<(n-i)) - n)/(n*(1<<(n-i)))
         */
        
        //TODO: because 360 is usually the max value, if N is > 80, can
        // use CountingSort here for better performance.
        // CountingSort is O(maxValue) while MergeSort is O(N*lg_2(N))
       
        MergeSort.sortByDecr(angles, 0, lastIndex);
        
        // visit in order of large values to low to store quadrant corrections
        
        double[] quadrantCorrected = new double[2];
        
        double avg = 0;
        for (int i = 0; i <= lastIndex; ++i) {
            
            int v = angles[i];
            
            if (i == 0) {
                avg = v;
            } else {
                double r1 = avg;
                avg = getAngleAverageInDegrees(r1, v, quadrantCorrected);
                
                if (quadrantCorrected[0] != r1) {
                    // TODO: this branch might not occur...revisit logic
                    angles[i - 1] = (int)Math.round(quadrantCorrected[0]);
                }
                if (quadrantCorrected[1] != v) {
                    angles[i] = (int)Math.round(quadrantCorrected[1]);
                }
            }
        }

        int n = lastIndex + 1;
        for (int i = 0; i <= lastIndex; ++i) {

            int v = angles[i];

            float c0;

            if (i == 0) {
                c0 = 1 << (n - 1);
            } else {
                c0 = 1 << (n - i);
            }

            float add = v * (c0 - n) / (n * c0);

            avg += add;
        }

        return (float)avg;
    }
    
    /**
     Runtime complexity is O(N*lg2N) for sort plus O(N).
     * Note, may change to use CountingSort for N > 80 and max(angles) ~ 360
     * in the future.
                       
     * @param angles angles in units of degrees.  Note that this array
     * is modified by use here so pass a copy if it should not be.
     * @param lastIndex last index to use in array angles, inclusive
     * @return 
     */
    public static float calculateWeightedAverageWithQuadrantCorrections(
        int[] angles, int[] valuesToMakeIntoWeights, int lastIndex) {
      
        //TODO: because 360 is usually the max value, if N is > 80, can
        // use CountingSort here for better performance.
        // CountingSort is O(maxValue) while MergeSort is O(N*lg_2(N))
       
        MultiArrayMergeSort.sortByDecr(angles, valuesToMakeIntoWeights, 0, lastIndex);
        
        // visit in order of large values to low to store quadrant corrections
        
        double[] quadrantCorrected = new double[2];
        
        double sumWeights = 0;
        double avg = 0;
        for (int i = 0; i <= lastIndex; ++i) {
            
            int v = angles[i];
            
            if (i == 0) {
                
                avg = v;
            
            } else {
                
                double r1 = avg;
                
                avg = getAngleAverageInDegrees(r1, v, quadrantCorrected);
                
                if (quadrantCorrected[0] != r1) {
                    // TODO: this branch might not occur...revisit logic
                    angles[i - 1] = (int)Math.round(quadrantCorrected[0]);
                }
                if (quadrantCorrected[1] != v) {
                    angles[i] = (int)Math.round(quadrantCorrected[1]);
                }
            }
            
            sumWeights += valuesToMakeIntoWeights[i];
        }
        
        // redo the mean using the quadrant corrected angles and the weights
        double mean = 0;

        for (int i = 0; i <= lastIndex; ++i) {

            double v = angles[i] * ((double)valuesToMakeIntoWeights[i]/sumWeights);

            mean += v;
        }

        return (float)mean;
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
