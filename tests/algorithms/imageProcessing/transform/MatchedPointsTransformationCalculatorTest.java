package algorithms.imageProcessing.transform;

import algorithms.imageProcessing.util.AngleUtil;
import algorithms.util.PairInt;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class MatchedPointsTransformationCalculatorTest extends TestCase {

    public MatchedPointsTransformationCalculatorTest() {
    }
    
    public void testCalulateEuclidean0() {

        int pt1Img1X = 254;
        int pt1Img1Y = 17;
        int pt1Img2X = 203;
        int pt1Img2Y = 251;

        int pt2Img1X = 265;
        int pt2Img1Y = 9;
        int pt2Img2X = 210;
        int pt2Img2Y = 259;

        double centroidX1 = 0.0;
        double centroidY1 = 0.0;
        MatchedPointsTransformationCalculator tc = new MatchedPointsTransformationCalculator();

        TransformationParameters result = tc.calulateEuclidean(
            pt1Img1X, pt1Img1Y, pt2Img1X, pt2Img1Y,
            pt1Img2X, pt1Img2Y, pt2Img2X, pt2Img2Y,
            centroidX1, centroidY1);

        float diffRotDeg = AngleUtil.getAngleDifference(
           result.getRotationInDegrees(), 360 - 90);
        float diffScale = Math.abs(result.getScale() - 1);
        float diffTransX = Math.abs(result.getTranslationX() - 10);
        float diffTransY = Math.abs(result.getTranslationY() - 0);

        assertTrue(diffRotDeg < 10);
        //assertTrue(diffScale == 0);
        //assertTrue(diffTransX == 0);
        //assertTrue(diffTransY == 0);
    }
    
    public void testCalulateEuclidean() {

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

        TransformationParameters result = tc.calulateEuclidean(
            set1X1, set1Y1, set1X2, set1Y2,
            set2X1, set2Y1, set2X2, set2Y2,
            centroidX1, centroidY1);

        float diffRotDeg = AngleUtil.getAngleDifference(
           result.getRotationInDegrees(), 0);
        float diffScale = Math.abs(result.getScale() - 1);
        float diffTransX = Math.abs(result.getTranslationX() - 10);
        float diffTransY = Math.abs(result.getTranslationY() - 0);

        assertTrue(diffRotDeg == 0);
        assertTrue(diffScale == 0);
        assertTrue(diffTransX == 0);
        assertTrue(diffTransY == 0);
    }

    public void testCalulateEuclidean2() {

        int set1X1 = 0; int set1Y1 = 0;
        int set1X2 = 1; int set1Y2 = 3;

        int set2X1 = 0; int set2Y1 = 2;
        int set2X2 = 3; int set2Y2 = 1;

        double centroidX1 = 1.0;
        double centroidY1 = 1.0;
        MatchedPointsTransformationCalculator tc = new MatchedPointsTransformationCalculator();

        TransformationParameters result = tc.calulateEuclidean(
            set1X1, set1Y1, set1X2, set1Y2,
            set2X1, set2Y1, set2X2, set2Y2,
            centroidX1, centroidY1);

        float diffRotDeg = AngleUtil.getAngleDifference(
           result.getRotationInDegrees(), 90);
        float diffScale = Math.abs(result.getScale() - 1);
        float diffTransX = Math.abs(result.getTranslationX() - 0);
        float diffTransY = Math.abs(result.getTranslationY() - 0);

        assertTrue(diffRotDeg < 0.001);
        assertTrue(diffScale < 0.001);
        assertTrue(diffTransX < 0.001);
        assertTrue(diffTransY < 0.001);
    }
    
    public void testCalulateEuclidean3() {

        int set1X1 = 0; int set1Y1 = 0;
        int set1X2 = 1; int set1Y2 = 3;

        int set2X1 = 2; int set2Y1 = 2;
        int set2X2 = 1; int set2Y2 = -1;

        double centroidX1 = 1.0;
        double centroidY1 = 1.0;
        MatchedPointsTransformationCalculator tc = new MatchedPointsTransformationCalculator();

        TransformationParameters result = tc.calulateEuclidean(
            set1X1, set1Y1, set1X2, set1Y2,
            set2X1, set2Y1, set2X2, set2Y2,
            centroidX1, centroidY1);

        float diffRotDeg = AngleUtil.getAngleDifference(
           result.getRotationInDegrees(), 180);
        float diffScale = Math.abs(result.getScale() - 1);
        float diffTransX = Math.abs(result.getTranslationX() - 0);
        float diffTransY = Math.abs(result.getTranslationY() - 0);

        assertTrue(diffRotDeg < 0.001);
        assertTrue(diffScale < 0.001);
        assertTrue(diffTransX < 0.001);
        assertTrue(diffTransY < 0.001);
    }
    
    public void testCalulateEuclidean4() {

        int set1X1 = 0; int set1Y1 = 0;
        int set1X2 = 1; int set1Y2 = 3;

        int set2X1 = 2; int set2Y1 = 0;
        int set2X2 = -1; int set2Y2 = 1;

        double centroidX1 = 1.0;
        double centroidY1 = 1.0;
        MatchedPointsTransformationCalculator tc = new MatchedPointsTransformationCalculator();

        TransformationParameters result = tc.calulateEuclidean(
            set1X1, set1Y1, set1X2, set1Y2,
            set2X1, set2Y1, set2X2, set2Y2,
            centroidX1, centroidY1);

        float diffRotDeg = AngleUtil.getAngleDifference(
           result.getRotationInDegrees(), 270);
        float diffScale = Math.abs(result.getScale() - 1);
        float diffTransX = Math.abs(result.getTranslationX() - 0);
        float diffTransY = Math.abs(result.getTranslationY() - 0);

        assertTrue(diffRotDeg < 0.001);
        assertTrue(diffScale < 0.001);
        assertTrue(diffTransX < 0.001);
        assertTrue(diffTransY < 0.001);
    }
    
    public void testCalulateEuclidean5() {

        int transX = -2;
        int transY = 1;
        
        int set1X1 = 0; int set1Y1 = 0;
        int set1X2 = 1; int set1Y2 = 3;

        int set2X1 = 2 + transX; int set2Y1 = 0 + transY;
        int set2X2 = -1 + transX; int set2Y2 = 1 + transY;

        double centroidX1 = 1.0;
        double centroidY1 = 1.0;
        MatchedPointsTransformationCalculator tc = new MatchedPointsTransformationCalculator();

        TransformationParameters result = tc.calulateEuclidean(
            set1X1, set1Y1, set1X2, set1Y2,
            set2X1, set2Y1, set2X2, set2Y2,
            centroidX1, centroidY1);

        float diffRotDeg = AngleUtil.getAngleDifference(
           result.getRotationInDegrees(), 270);
        float diffScale = Math.abs(result.getScale() - 1);
        float diffTransX = Math.abs(result.getTranslationX() - transX);
        float diffTransY = Math.abs(result.getTranslationY() - transY);

        assertTrue(diffRotDeg < 0.001);
        assertTrue(diffScale < 0.001);
        assertTrue(diffTransX < 0.001);
        assertTrue(diffTransY < 0.001);
    }
    
    public void testCalulateEuclidean6() {

        int transX = 1;
        int transY = -1;
        
        int scale = 2;
        
        int set1X1 = 0; int set1Y1 = 0;
        int set1X2 = 1; int set1Y2 = 3;

        int set2X1 = 2*scale + transX; int set2Y1 = 0*scale + transY;
        int set2X2 = -1*scale + transX; int set2Y2 = 1*scale + transY;

        double centroidX1 = 1.0;
        double centroidY1 = 1.0;
        MatchedPointsTransformationCalculator tc = new MatchedPointsTransformationCalculator();

        TransformationParameters result = tc.calulateEuclidean(
            set1X1, set1Y1, set1X2, set1Y2,
            set2X1, set2Y1, set2X2, set2Y2,
            centroidX1, centroidY1);

        float diffRotDeg = AngleUtil.getAngleDifference(
           result.getRotationInDegrees(), 270);
        float diffScale = Math.abs(result.getScale() - scale);
        float diffTransX = Math.abs(result.getTranslationX() - transX);
        float diffTransY = Math.abs(result.getTranslationY() - transY);

        assertTrue(diffRotDeg < 0.001);
        assertTrue(diffScale < 0.001);
        assertTrue(diffTransX < 0.001);
        assertTrue(diffTransY < 0.001);
    }

    public void testCalulateEuclidean90() {

        PairInt c1Pt1 = new PairInt(209, 271);
        PairInt c2Pt1 = new PairInt(30, 213);
        
        PairInt c1Pt2 = new PairInt(204, 283);
        PairInt c2Pt2 = new PairInt(15, 209);

        double centroidX1 = 0.0;
        double centroidY1 = 0.0;
        MatchedPointsTransformationCalculator tc = new MatchedPointsTransformationCalculator();

        TransformationParameters result = tc.calulateEuclidean(
            c1Pt1.getX(), c1Pt1.getY(), c1Pt2.getX(), c1Pt2.getY(),
            c2Pt1.getX(), c2Pt1.getY(), c2Pt2.getX(), c2Pt2.getY(), 0, 0);

        float diffRotDeg = AngleUtil.getAngleDifference(
           result.getRotationInDegrees(), 360 - 90);
        float diffScale = Math.abs(result.getScale() - 1);
        float diffTransX = Math.abs(result.getTranslationX() - 10);
        float diffTransY = Math.abs(result.getTranslationY() - 0);

        assertTrue(diffRotDeg < 10);
        //assertTrue(diffScale == 0);
        //assertTrue(diffTransX == 0);
        //assertTrue(diffTransY == 0);
    }
    
    public void test2() throws Exception {
        
        PairInt c1Pt1 = new PairInt(47, 202);
        PairInt c1Pt2 = new PairInt(127, 157);
        
        PairInt c2Pt1 = new PairInt(164, 51);
        PairInt c2Pt2 = new PairInt(212, 128);
        
        MatchedPointsTransformationCalculator tc = new MatchedPointsTransformationCalculator();

        TransformationParameters result = tc.calulateEuclidean(
            c1Pt1.getX(), c1Pt1.getY(), c1Pt2.getX(), c1Pt2.getY(),
            c2Pt1.getX(), c2Pt1.getY(), c2Pt2.getX(), c2Pt2.getY(), 0, 0);

        float diffRotDeg = AngleUtil.getAngleDifference(
           result.getRotationInDegrees(), 360 - 90);
        float diffScale = Math.abs(result.getScale() - 1);
        float diffTransX = Math.abs(result.getTranslationX() - 361);
        float diffTransY = Math.abs(result.getTranslationY() - -4.4f);

        assertTrue(diffRotDeg < 10);
        assertTrue(diffScale < 0.05);
        assertTrue(diffTransX < 5);
        assertTrue(diffTransY < 5);
        
    }
}
