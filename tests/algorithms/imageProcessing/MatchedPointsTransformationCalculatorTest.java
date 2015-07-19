package algorithms.imageProcessing;

import algorithms.util.PairIntArray;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class MatchedPointsTransformationCalculatorTest extends TestCase {
    
    public MatchedPointsTransformationCalculatorTest() {
    }

    public void testCalulateEuclideanGivenScale() {
        
        int set1X1 = 1;
        int set1Y1 = 2;
        int set1X2 = 2;
        int set1Y2 = 2;
        
        int set2X1 = 11;
        int set2Y1 = 2;
        int set2X2 = 12;
        int set2Y2 = 2;
        
        double centroidX1 = 0.0;
        double centroidY1 = 0.0;
        MatchedPointsTransformationCalculator tc = new MatchedPointsTransformationCalculator();
        
        TransformationParameters result = tc.calulateEuclideanGivenScale(
            set1X1, set1Y1, set1X2, set1Y2,
            set2X1, set2Y1, set2X2, set2Y2,
            centroidX1, centroidY1);
        
        float diffRotDeg = getAngleDifference(
           result.getRotationInDegrees(), 9);
        float diffScale = Math.abs(result.getScale() - 1);
        float diffTransX = Math.abs(result.getTranslationX() - 10);
        float diffTransY = Math.abs(result.getTranslationY() - 0);

        assertTrue(diffRotDeg == 0);
        assertTrue(diffScale == 0);
        assertTrue(diffTransX == 0);
        assertTrue(diffTransY == 0);
    }

    private float getAngleDifference(float rotDegrees0, float rotDegrees1) {
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
