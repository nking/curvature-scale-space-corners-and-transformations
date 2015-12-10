package algorithms.imageProcessing.util;

import algorithms.MergeSort;
import algorithms.MultiArrayMergeSort;
import java.util.List;

/**
 *
 * @author nichole
 */
public class AngleUtil {

    /**
    calculates the difference of angles between pairs of points in set1 and 
    set 2 ( angle of (diffX1, diffY1) - angle of (diffX2, diffY2).  
    * diffX1, diffX2 are the difference
    * in x and y between 2 points in set 1, and diffX2, diffY2 are the 
    * difference for the same matched point pair in set2.  The reference
    * frame is polar clockwise.
    * For example, (diffX1, diffY1) being (1, 4) and (diffX2, diffY2) being
    * (3.536, 2.12) leads to angles 284 minus 329 = -45.  Note that to
    * transform (diffX1, diffY1) to (diffX2, diffY2) one would apply 
    * -1*result of this method.
    * 
    * <pre>
    *       clockwise -->
              +Y        V
              270
           III | IV
       180 --------- 0   +X
           II  |  I
               90
    * </pre>
    * Note that the subtraction is from closest angles.
    * For example, if point 1 is in qI and point 2 is in qIV, 360 is added to 
    * point 1's angle, then the equation is (2*PI + point1 angle) - (point 2 angle).

     * @param diffX1
     * @param diffY1
     * @param diffX2
     * @param diffY2
     * @return 
     */
    public double subtract(double diffX1, double diffY1, double diffX2, 
        double diffY2) {
        
        /*  clockwise -->
              +Y        V
              270
           III | IV
       180 --------- 0   +X
           II  |  I
               90
        */
        
        double theta1 = polarAngleCW(diffX1, diffY1);
        
        double theta2 = polarAngleCW(diffX2, diffY2);
        
        boolean useRadians = true;

        double[] quadCorrThetas = correctForQuadrants(theta1, theta2, 
            useRadians);
        
        return (quadCorrThetas[0] - quadCorrThetas[1]);        
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
        
        
        Math.atan2 angles are:
                 90
           135    |    45
                  |
        180 ---------------  0
                  |
          -135    |   -45
                 -90
        so, for d < 0, need 360+d
        */
        
        double theta = Math.atan2(y, x);
        
        if (theta < 0) {
            theta += 2. * Math.PI;
        }
        
        return theta;  
    }
    
    /**
    calculates the polar theta in radians given x and y w.r.t. origin.  theta increases
    * in value in a clockwise direction (CW).

     * @param x
     * @param y
     * @return 
     */
    public static double polarAngleCW(double x, double y) {
        
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
        
        Math.atan2 angles are:
                 90
           135    |    45
                  |
        180 ---------------  0
                  |
          -135    |   -45
                 -90
        so need -1*d, then for d < 0, need 360+d
        */
        
        double theta = -1 * Math.atan2(y, x);
        
        if (theta < 0) {
            theta += 2. * Math.PI;
        }
        
        return theta;  
    }
    
    public static double[] correctForQuadrants(double rot0, double rot1, 
        boolean useRadians) {
        
        if (rot0 < 0 || rot1 < 0) {
            throw new IllegalArgumentException(
                "rot0 and rot1 cannot be negative numbers");
        }
       
        double[] corrected = new double[2];
        correctForQuadrants(rot0, rot1, useRadians, corrected);
        return corrected;
    }
    
    public static void correctForQuadrants(double rot0, double rot1, 
        boolean useRadians, double[] outputCorrectedForQuadrants) {
        
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
        
        double twoPI = 2. * Math.PI;

        // prefer to keep 360 instead of 0 because some of the invokers of this
        // method might prefer a fraction times results to be > 0
        if (useRadians) {
            if (rot0 != twoPI) {
                rot0 = rot0 % twoPI;
            }
            if (rot1 != twoPI) {
                rot1 = rot1 % twoPI;
            }
        } else {
            if (rot0 != 360) {
                rot0 = rot0 % 360;
            }
            if (rot1 != 360) {
                rot1 = rot1 % 360;
            }
        }
        int q0, q1;
        if (useRadians) {            
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
        
        outputCorrectedForQuadrants[0] = rot0;
        outputCorrectedForQuadrants[1] = rot1;
    }
    
    private static void correctForQuadrantsWithoutMod(double rot0, double rot1, 
        boolean useRadians, double[] outputCorrectedForQuadrants) {
        
        if (rot0 < 0 || rot1 < 0) {
            throw new IllegalArgumentException("angles must be non negative");
        }
        
        /*
                  270
               III | IV
          180  --------- 0
               II  |  I
                   90
        */

        double twoPI;
        int q0, q1;
        if (useRadians) {
            twoPI = 2. * Math.PI;
            q0 = (int)(rot0/(0.5 * Math.PI)) + 1;
            q1 = (int)(rot1/(0.5 * Math.PI)) + 1;
        } else {
            twoPI = 360.;
            q0 = (int)(rot0/90.f) + 1;
            q1 = (int)(rot1/90.f) + 1;
        }
        
        /*
                  270
               III | IV
          180  --------- 0
               II  |  I
                   90
        */
        
        //TODO: revisit these and consider only editing rot1, and adding to
        // method comments that only rot1 will be corrected w.r.t rot0
        
        if (q0 == 1) {
            if (q1 == 3) {
                double diff = rot1 - rot0;
                if (diff > Math.PI) {
                    rot0 += twoPI;
                }
            } else if (q1 == 4) {
                rot0 += twoPI;
            } else if (q1 == 5) {
                //double diff = rot1 - Math.PI - rot0;
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
                if (diff > Math.PI) {
                    rot0 += twoPI;
                }
            } else if (q1 == 5) {
                rot1 -= twoPI;
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
                if (diff > Math.PI) {
                    rot1 += twoPI;
                }
            } else if (q1 == 5) {
                rot1 -= twoPI;
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
                if (diff > Math.PI) {
                    rot1 += twoPI;
                }
            }
        } else if (q0 == 5) {
            /*
                  270
               III | IV
          180  --------- 0
               II  |  I  <---V too
                   90
            */
            if (q1 == 1) {
                rot1 += twoPI;
            } else if (q1 == 2) {
                rot0 -= twoPI;
            } else if (q1 == 3) {  
                rot0 -= twoPI;
            }
        }
        
        outputCorrectedForQuadrants[0] = rot0;
        outputCorrectedForQuadrants[1] = rot1;
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
    
    /**
     Runtime complexity is O(N*lg2N) for sort plus O(N).
     * Note, may change to use CountingSort for N > 80 and max(angles) ~ 360
     * in the future.
                       
     * @param angles angles in units of radians or degrees.  Note that this array
     * is modified by use here so pass a copy if it should not be.
     * @param lastIndex last index to use in array angles, inclusive
     * @return 
     */
    public static float calculateWeightedAverageWithQuadrantCorrections(
        Double[] angles, Float[] weights, int lastIndex, boolean useRadians) {
      
        //TODO: because 360 is usually the max value, if N is > 80, can
        // use CountingSort here for better performance.
        // CountingSort is O(maxValue) while MergeSort is O(N*lg_2(N))
        MultiArrayMergeSort.sortByDecr(angles, weights, 0, lastIndex);
      
        double[] corr = new double[2];
        
        for (int i = 0; i < angles.length; ++i) {
            double a0 = angles[i];
            double a1;
            if (i == (angles.length - 1)) {
                a1 = angles[0];
            } else {
                a1 = angles[i + 1];
            }
            correctForQuadrants(a0, a1, useRadians, corr);
            angles[i] = corr[0];
            if (i == (angles.length - 1)) {
                angles[0] = corr[1];
            } else {
                angles[i + 1] = corr[1];
            }
        }
                   
        // redo the mean using the quadrant corrected angles and the weights
        double mean = 0;

        for (int i = 0; i <= lastIndex; ++i) {
            double v = angles[i] * weights[i];
            mean += v;
        }

        return (float)mean;
    }
    
    /**
     Runtime complexity is O(N*lg2N) for sort plus O(N).
     * Note, may change to use CountingSort for N > 80 and max(angles) ~ 360
     * in the future.
              
     */
    public static void correctForQuadrants(List<Double> thetas, 
        boolean useRadians) {
        
        if (thetas == null || thetas.size() < 2) {
            return;
        }
        
        double[] angles = new double[thetas.size()];
        int[] indexes = new int[thetas.size()];
        
        for (int i = 0; i < thetas.size(); ++i) {
            angles[i] = thetas.get(i).doubleValue();
            indexes[i] = i;
        }
        
        MultiArrayMergeSort.sortByDecr(angles, indexes);
      
        double[] corr = new double[2];
        
        for (int i = 0; i < angles.length; ++i) {
            double a0 = angles[i];
            double a1;
            if (i == (angles.length - 1)) {
                a1 = angles[0];
            } else {
                a1 = angles[i + 1];
            }
            correctForQuadrantsWithoutMod(a0, a1, useRadians, corr);
            angles[i] = corr[0];
            if (i == (angles.length - 1)) {
                angles[0] = corr[1];
            } else {
                angles[i + 1] = corr[1];
            }
        }
        
        for (int i = 0; i < angles.length; ++i) {
            int idx = indexes[i];
            double v = angles[i];
            thetas.set(idx, v);
        }        
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
    
    protected static int getClockwiseQuadrant(double xCoord, double yCoord) {
        /*
          III | IV
          ---------
          II  |  I
        */
        
        int q = 1;
        if (xCoord == 0) {
            if (yCoord == 0) {
                //TODO: might be more consistent for origin to be in 4
                q = 1;
            } else if (yCoord > 0) {
                q = 4;
            } else {
                q = 1;
            }
        } else if (yCoord == 0) {
            if (xCoord > 0) {
                q = 1;
            } else {
                q = 3;
            }            
        } else if ((xCoord > 0) && (yCoord > 0)) {
            q = 4;
        } else if ((xCoord > 0) && (yCoord < 0)) {
            q = 1;
        } else if ((xCoord < 0) && (yCoord > 0)) {
            q = 3;
        } else if ((xCoord < 0) && (yCoord < 0)) {
            q = 2;
        }
        
        return q;
    }
}
