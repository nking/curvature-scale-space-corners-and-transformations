package algorithms.imageProcessing;

import algorithms.util.PairFloatArray;
import junit.framework.TestCase;
import org.ejml.simple.SimpleMatrix;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author nichole
 */
public class RANSACSolverTest extends TestCase {
    
    public RANSACSolverTest() {
    }

    public void testRANSAC() throws Exception {
        
        PairFloatArray leftTrueMatches = new PairFloatArray();
        PairFloatArray rightTrueMatches = new PairFloatArray();
        getMertonCollege10TrueMatches(leftTrueMatches, rightTrueMatches);
        
        PairFloatArray leftFalseMatches = new PairFloatArray();
        PairFloatArray rightFalseMatches = new PairFloatArray();
        getMertonCollegeFalseMatch1(leftFalseMatches, rightFalseMatches);
        getMertonCollegeFalseMatch2(leftFalseMatches, rightFalseMatches);
        getMertonCollegeFalseMatch3(leftFalseMatches, rightFalseMatches);
        
        PairFloatArray left = leftTrueMatches.copy();
        PairFloatArray right = rightTrueMatches.copy();
        getMertonCollegeFalseMatch1(left, right);
        getMertonCollegeFalseMatch2(left, right);
        getMertonCollegeFalseMatch3(left, right);
        
        PairFloatArray outputLeft = new PairFloatArray(); 
        PairFloatArray outputRight = new PairFloatArray();
        
        RANSACSolver solver = new RANSACSolver();
        
        StereoProjectionTransformerFit fit = solver.calculateEpipolarProjection(
            left, right, outputLeft, outputRight);
        
        assertNotNull(fit);
        
        assertTrue(outputLeft.getN() == leftTrueMatches.getN());
        assertTrue(outputRight.getN() == rightTrueMatches.getN());
        
        for (int i = 0; i < outputLeft.getN(); ++i) {
            float xL = outputLeft.getX(i);
            float yL = outputLeft.getY(i);
            float xR = outputRight.getX(i);
            float yR = outputRight.getY(i);
            boolean removed = false;
            for (int j = 0; j < left.getN(); ++j) {
                float diffXL = Math.abs(xL - left.getX(j));
                float diffYL = Math.abs(yL - left.getY(j));
                float diffXR = Math.abs(xR - right.getX(j));
                float diffYR = Math.abs(yR - right.getY(j));
                if ((diffXL < 0.1) && (diffYL < 0.1) && (diffXR < 0.1) && (diffYR < 0.1)) {
                    left.removeRange(j, j);
                    right.removeRange(j, j);
                    removed = true;
                    break;
                }
            }
            assertTrue(removed);
        }
        assertTrue(left.getN() == 0);
        assertTrue(right.getN() == 0);
    }
    
    protected void getMertonCollege10TrueMatches(PairFloatArray left, 
        PairFloatArray right) {
        
        /*
        58, 103   32, 100
        486, 46   474, 49
        845, 127  878, 151
        949, 430  998, 471
        541, 428  533, 460
        225, 453  213, 498
        49, 509   21, 571
        373, 239  365, 258
        737, 305  762, 335
        84, 273   60, 298
        */
        
        left.add(58, 103);  right.add(32, 100);
        left.add(486, 46);   right.add(474, 49);
        left.add(845, 127);   right.add(878, 151);
        left.add(949, 430);   right.add(998, 471);
        left.add(541, 428);   right.add(533, 460);
        left.add(225, 453);   right.add(213, 498);
        left.add(49, 509);   right.add(21, 571);
        left.add(373, 239);   right.add(365, 258);
        left.add(737, 305);   right.add(762, 335);
        left.add(84, 273);   right.add(60, 298);
    }
    
    protected void getMertonCollegeFalseMatch1(PairFloatArray left, 
        PairFloatArray right) {
        //765, 487   753, 552
        left.add(765, 487);   right.add(753, 552);
    }
    protected void getMertonCollegeFalseMatch2(PairFloatArray left, 
        PairFloatArray right) {
        //253, 141    256, 229
        left.add(253, 141);   right.add(256, 229);
    }
    protected void getMertonCollegeFalseMatch3(PairFloatArray left, 
        PairFloatArray right) {
        //459, 354  432, 525
        left.add(459, 354);   right.add(432, 525);
    }
    
}
