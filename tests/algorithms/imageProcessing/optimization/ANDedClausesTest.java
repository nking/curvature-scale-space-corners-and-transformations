package algorithms.imageProcessing.optimization;

import java.util.Map;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class ANDedClausesTest extends TestCase {
    
    public ANDedClausesTest() {
    }
    
    public void testEvaluate0() throws Exception {
        
        boolean skyIsRed = false;
        double pixContrast = 8;
        double pixBlueOrRedDiff = 7;
        double pixCIEXDiff = 6;
        double pixCIEYDiff = 5;
        double skyStDevContrast = 2.0;
        double skyStDevBlueOrRedDiff = 3.0;
        double skyStDevCIEX = 1;
        double skyStDevCIEY = 1;
        int red = 127;
        int green = 127;
        int blue = 127;
        
        ColorData data = new ColorData(skyIsRed,
            pixContrast, pixBlueOrRedDiff, 
            pixCIEXDiff, pixCIEYDiff, 
            skyStDevContrast, skyStDevBlueOrRedDiff,
            skyStDevCIEX, skyStDevCIEY, red, green, blue);
        
        //8/2 > 2
        
        ANDedClauses clauses = new ANDedClauses(1, SKYCONDITIONAL.ALL);
        
        clauses.set(0, PARAM.ABSOLUTE_CONTRAST, 
            PARAM.STDEV_CONTRAST, COMPARISON.GREATER_THAN, 2.0f);
        
        assertTrue(clauses.evaluate(data));
        
        clauses = new ANDedClauses(1, SKYCONDITIONAL.BLUE);
        
        clauses.set(0, PARAM.ABSOLUTE_CONTRAST, 
            PARAM.STDEV_CONTRAST, COMPARISON.GREATER_THAN, 2.0f);
        
        // sky is blue and conditionals are true
        assertTrue(clauses.evaluate(data));
        
        
        clauses = new ANDedClauses(1, SKYCONDITIONAL.RED);
        
        clauses.set(0, PARAM.ABSOLUTE_CONTRAST, 
            PARAM.STDEV_CONTRAST, COMPARISON.GREATER_THAN, 2.0f);
        
        // sky is not red
        assertFalse(clauses.evaluate(data));
        
        
        // ----- assert opposite for coefficient 10. -------
        //8/2 > 10
        
        clauses = new ANDedClauses(1, SKYCONDITIONAL.ALL);
        
        clauses.set(0, PARAM.ABSOLUTE_CONTRAST, 
            PARAM.STDEV_CONTRAST, COMPARISON.GREATER_THAN, 10.0f);
        
        assertFalse(clauses.evaluate(data));
        
        clauses = new ANDedClauses(1, SKYCONDITIONAL.BLUE);
        
        clauses.set(0, PARAM.ABSOLUTE_CONTRAST, 
            PARAM.STDEV_CONTRAST, COMPARISON.GREATER_THAN, 10.0f);
        
        // sky is blue and conditionals are not true
        assertFalse(clauses.evaluate(data));
        
        
        clauses = new ANDedClauses(1, SKYCONDITIONAL.RED);
        
        clauses.set(0, PARAM.ABSOLUTE_CONTRAST, 
            PARAM.STDEV_CONTRAST, COMPARISON.GREATER_THAN, 10.0f);
        
        // sky is not red
        assertFalse(clauses.evaluate(data));
    }
    
    public void testEvaluate1() throws Exception {
        
        boolean skyIsRed = false;
        double pixContrast = 8;
        double pixBlueOrRedDiff = 7;
        double pixCIEXDiff = 6;
        double pixCIEYDiff = 5;
        double skyStDevContrast = 2.0;
        double skyStDevBlueOrRedDiff = 3.0;
        double skyStDevCIEX = 1;
        double skyStDevCIEY = 1;
        int red = 127;
        int green = 127;
        int blue = 127;
        
        ColorData data = new ColorData(skyIsRed,
            pixContrast, pixBlueOrRedDiff, 
            pixCIEXDiff, pixCIEYDiff, 
            skyStDevContrast, skyStDevBlueOrRedDiff,
            skyStDevCIEX, skyStDevCIEY, red, green, blue);
        
        //8/2 < 7
        
        ANDedClauses clauses = new ANDedClauses(1, SKYCONDITIONAL.ALL);
        
        clauses.set(0, PARAM.ABSOLUTE_CONTRAST, 
            PARAM.STDEV_CONTRAST, COMPARISON.LESS_THAN, 7.0f);
        
        assertTrue(clauses.evaluate(data));
        
        clauses = new ANDedClauses(1, SKYCONDITIONAL.BLUE);
        
        clauses.set(0, PARAM.ABSOLUTE_CONTRAST, 
            PARAM.STDEV_CONTRAST, COMPARISON.LESS_THAN, 7.0f);
        
        // sky is blue and conditionals are true
        assertTrue(clauses.evaluate(data));
        
        
        clauses = new ANDedClauses(1, SKYCONDITIONAL.RED);
        
        clauses.set(0, PARAM.ABSOLUTE_CONTRAST, 
            PARAM.STDEV_CONTRAST, COMPARISON.LESS_THAN, 7.0f);
        
        // sky is not red
        assertFalse(clauses.evaluate(data));
        
        
        // ----- assert opposite for coefficient 2. -------
        //8/2 < 2
        clauses = new ANDedClauses(1, SKYCONDITIONAL.ALL);
        
        clauses.set(0, PARAM.ABSOLUTE_CONTRAST, 
            PARAM.STDEV_CONTRAST, COMPARISON.LESS_THAN, 2.0f);
        
        assertFalse(clauses.evaluate(data));
        
        clauses = new ANDedClauses(1, SKYCONDITIONAL.BLUE);
        
        clauses.set(0, PARAM.ABSOLUTE_CONTRAST, 
            PARAM.STDEV_CONTRAST, COMPARISON.LESS_THAN, 2.0f);
        
        // sky is blue and conditionals are false
        assertFalse(clauses.evaluate(data));
        
        
        clauses = new ANDedClauses(1, SKYCONDITIONAL.RED);
        
        clauses.set(0, PARAM.ABSOLUTE_CONTRAST, 
            PARAM.STDEV_CONTRAST, COMPARISON.LESS_THAN, 2.0f);
        
        // sky is not red
        assertFalse(clauses.evaluate(data));
    }
    
    public void testEvaluate2() throws Exception {
        
        boolean skyIsRed = false;
        double pixContrast = 8;
        double pixBlueOrRedDiff = 7;
        double pixCIEXDiff = 6;
        double pixCIEYDiff = 5;
        double skyStDevContrast = 2.0;
        double skyStDevBlueOrRedDiff = 3.0;
        double skyStDevCIEX = 1;
        double skyStDevCIEY = 1;
        int red = 127;
        int green = 127;
        int blue = 127;
        
        ColorData data = new ColorData(skyIsRed,
            pixContrast, pixBlueOrRedDiff, 
            pixCIEXDiff, pixCIEYDiff, 
            skyStDevContrast, skyStDevBlueOrRedDiff,
            skyStDevCIEX, skyStDevCIEY, red, green, blue);
        
        //8/2 > 2
        
        ANDedClauses clauses = new ANDedClauses(
            SKYCONDITIONAL.ALL, 
            new PARAM[]{PARAM.ABSOLUTE_CONTRAST}, 
            new PARAM[]{PARAM.STDEV_CONTRAST}, 
            new COMPARISON[]{COMPARISON.GREATER_THAN}, 
            new float[]{2.0f});
        
        assertTrue(clauses.evaluate(data));
                
        clauses = new ANDedClauses(
            SKYCONDITIONAL.BLUE, 
            new PARAM[]{PARAM.ABSOLUTE_CONTRAST}, 
            new PARAM[]{PARAM.STDEV_CONTRAST}, 
            new COMPARISON[]{COMPARISON.GREATER_THAN}, 
            new float[]{2.0f});
        
        // sky is blue and conditionals are true
        assertTrue(clauses.evaluate(data));
        
        
        clauses = new ANDedClauses(
            SKYCONDITIONAL.RED, 
            new PARAM[]{PARAM.ABSOLUTE_CONTRAST}, 
            new PARAM[]{PARAM.STDEV_CONTRAST}, 
            new COMPARISON[]{COMPARISON.GREATER_THAN}, 
            new float[]{2.0f});
        
        // sky is not red
        assertFalse(clauses.evaluate(data));
        
        
        // ----- assert opposite for coefficient 10. -------
        //8/2 > 10
        
        clauses = new ANDedClauses(
            SKYCONDITIONAL.ALL, 
            new PARAM[]{PARAM.ABSOLUTE_CONTRAST}, 
            new PARAM[]{PARAM.STDEV_CONTRAST}, 
            new COMPARISON[]{COMPARISON.GREATER_THAN}, 
            new float[]{10.0f});
        
        assertFalse(clauses.evaluate(data));
        
        clauses = new ANDedClauses(
            SKYCONDITIONAL.BLUE, 
            new PARAM[]{PARAM.ABSOLUTE_CONTRAST}, 
            new PARAM[]{PARAM.STDEV_CONTRAST}, 
            new COMPARISON[]{COMPARISON.GREATER_THAN}, 
            new float[]{10.0f});
        
        // sky is blue and conditionals are not true
        assertFalse(clauses.evaluate(data));
        
        
        clauses = new ANDedClauses(
            SKYCONDITIONAL.RED, 
            new PARAM[]{PARAM.ABSOLUTE_CONTRAST}, 
            new PARAM[]{PARAM.STDEV_CONTRAST}, 
            new COMPARISON[]{COMPARISON.GREATER_THAN}, 
            new float[]{10.0f});
        
        // sky is not red
        assertFalse(clauses.evaluate(data));
    }
    
    public void testEvaluate3() throws Exception {
        
        boolean skyIsRed = false;
        double pixContrast = 8;
        double pixBlueOrRedDiff = 7;
        double pixCIEXDiff = 6;
        double pixCIEYDiff = 5;
        double skyStDevContrast = 2.0;
        double skyStDevBlueOrRedDiff = 3.0;
        double skyStDevCIEX = 1;
        double skyStDevCIEY = 1;
        int red = 127;
        int green = 127;
        int blue = 127;
        
        ColorData data = new ColorData(skyIsRed,
            pixContrast, pixBlueOrRedDiff, 
            pixCIEXDiff, pixCIEYDiff, 
            skyStDevContrast, skyStDevBlueOrRedDiff,
            skyStDevCIEX, skyStDevCIEY, red, green, blue);
        
        //8/2 < 7
        
        ANDedClauses clauses = new ANDedClauses(
            SKYCONDITIONAL.ALL, 
            new PARAM[]{PARAM.ABSOLUTE_CONTRAST}, 
            new PARAM[]{PARAM.STDEV_CONTRAST}, 
            new COMPARISON[]{COMPARISON.LESS_THAN}, 
            new float[]{7.0f});
        
        assertTrue(clauses.evaluate(data));
        
        clauses = new ANDedClauses(
            SKYCONDITIONAL.BLUE, 
            new PARAM[]{PARAM.ABSOLUTE_CONTRAST}, 
            new PARAM[]{PARAM.STDEV_CONTRAST}, 
            new COMPARISON[]{COMPARISON.LESS_THAN}, 
            new float[]{7.0f});
        
        // sky is blue and conditionals are true
        assertTrue(clauses.evaluate(data));
        
        clauses = new ANDedClauses(
            SKYCONDITIONAL.RED, 
            new PARAM[]{PARAM.ABSOLUTE_CONTRAST}, 
            new PARAM[]{PARAM.STDEV_CONTRAST}, 
            new COMPARISON[]{COMPARISON.LESS_THAN}, 
            new float[]{7.0f});
        
        // sky is not red
        assertFalse(clauses.evaluate(data));
        
        
        // ----- assert opposite for coefficient 2. -------
        //8/2 < 2
        
        clauses = new ANDedClauses(
            SKYCONDITIONAL.ALL, 
            new PARAM[]{PARAM.ABSOLUTE_CONTRAST}, 
            new PARAM[]{PARAM.STDEV_CONTRAST}, 
            new COMPARISON[]{COMPARISON.LESS_THAN}, 
            new float[]{2.0f});
        
        assertFalse(clauses.evaluate(data));
        
        clauses = new ANDedClauses(
            SKYCONDITIONAL.BLUE, 
            new PARAM[]{PARAM.ABSOLUTE_CONTRAST}, 
            new PARAM[]{PARAM.STDEV_CONTRAST}, 
            new COMPARISON[]{COMPARISON.LESS_THAN}, 
            new float[]{2.0f});
        
        // sky is blue and conditionals are false
        assertFalse(clauses.evaluate(data));
        
        clauses = new ANDedClauses(
            SKYCONDITIONAL.RED, 
            new PARAM[]{PARAM.ABSOLUTE_CONTRAST}, 
            new PARAM[]{PARAM.STDEV_CONTRAST}, 
            new COMPARISON[]{COMPARISON.LESS_THAN}, 
            new float[]{2.0f});
        
        // sky is not red
        assertFalse(clauses.evaluate(data));
    }
    
    public void testEvaluate10() throws Exception {
        
        boolean skyIsRed = false;
        double pixContrast = 8;
        double pixBlueOrRedDiff = 7;
        double pixCIEXDiff = 6;
        double pixCIEYDiff = 5;
        double skyStDevContrast = 2.0;
        double skyStDevBlueOrRedDiff = 1.5;
        double skyStDevCIEX = 1;
        double skyStDevCIEY = 1;
        int red = 127;
        int green = 127;
        int blue = 127;
        
        ColorData data = new ColorData(skyIsRed,
            pixContrast, pixBlueOrRedDiff, 
            pixCIEXDiff, pixCIEYDiff, 
            skyStDevContrast, skyStDevBlueOrRedDiff,
            skyStDevCIEX, skyStDevCIEY, red, green, blue);
        
        //contr/stdev > 3.5 and clrdiff/stdv > 2.5
        
        ANDedClauses clauses = new ANDedClauses(
            SKYCONDITIONAL.ALL, 
            new PARAM[]{PARAM.ABSOLUTE_CONTRAST, PARAM.ABSOLUTE_DIFF_BLUE_OR_RED}, 
            new PARAM[]{PARAM.STDEV_CONTRAST, PARAM.STDEV_BLUE_OR_RED}, 
            new COMPARISON[]{COMPARISON.GREATER_THAN, COMPARISON.GREATER_THAN}, 
            new float[]{3.5f, 2.5f});
        
        assertTrue(clauses.evaluate(data));
        
        //contr/stdev > 3.5 and clrdiff/stdv > 5.5
        
        clauses = new ANDedClauses(
            SKYCONDITIONAL.ALL, 
            new PARAM[]{PARAM.ABSOLUTE_CONTRAST, PARAM.ABSOLUTE_DIFF_BLUE_OR_RED}, 
            new PARAM[]{PARAM.STDEV_CONTRAST, PARAM.STDEV_BLUE_OR_RED}, 
            new COMPARISON[]{COMPARISON.GREATER_THAN, COMPARISON.GREATER_THAN}, 
            new float[]{3.5f, 5.5f});
        
        assertFalse(clauses.evaluate(data));
    }
    
    public void testEvaluate11() throws Exception {
        
        boolean skyIsRed = false;
        double pixContrast = 8;
        double pixBlueOrRedDiff = 7;
        double pixCIEXDiff = 6;
        double pixCIEYDiff = 5;
        double skyStDevContrast = 2.0;
        double skyStDevBlueOrRedDiff = 1.5;
        double skyStDevCIEX = 1;
        double skyStDevCIEY = 1;
        int red = 127;
        int green = 127;
        int blue = 127;
        
        ColorData data = new ColorData(skyIsRed, pixContrast, pixBlueOrRedDiff, 
            pixCIEXDiff, pixCIEYDiff, skyStDevContrast, skyStDevBlueOrRedDiff,
            skyStDevCIEX, skyStDevCIEY, red, green, blue);
        
        //contr/stdev > 3.5
        
        CustomCoeff cc = new CustomCoeff() {
            @Override
            public double evaluate(ColorData data, float[] c, 
                Map<Integer, Float> cc) {
                
                Float coeff1 = cc.get(Integer.valueOf(1));
                Float coeff2 = cc.get(Integer.valueOf(2));
                
                float a = coeff1.floatValue() - 
                    (float)data.getParameter(PARAM.CONTRAST)
                    * coeff2.floatValue();
                
                return 2.5f;
            }
        };
        
        ANDedClauses clauses = new ANDedClauses(
            SKYCONDITIONAL.ALL, 
            new PARAM[]{PARAM.ABSOLUTE_CONTRAST}, 
            new PARAM[]{PARAM.STDEV_CONTRAST}, 
            new COMPARISON[]{COMPARISON.GREATER_THAN}, 
            new float[]{3.5f});
        
        clauses.setACustomCoefficient(0, cc);
        clauses.setCustomCoefficientVariable(1, Float.valueOf(2.5f));
        clauses.setCustomCoefficientVariable(2, Float.valueOf(-2.f));
        
        assertTrue(clauses.evaluate(data));
        
        
        // ---- same except make custom coefficient result in a high value ----
        
        data = new ColorData(skyIsRed, pixContrast, pixBlueOrRedDiff, 
            pixCIEXDiff, pixCIEYDiff, skyStDevContrast, skyStDevBlueOrRedDiff,
            skyStDevCIEX, skyStDevCIEY, red, green, blue);
        
        //contr/stdev > 3.5
        
        cc = new CustomCoeff() {
            @Override
            public double evaluate(ColorData data, float[] c, 
                Map<Integer, Float> cc) {
                
                return 100;
            }
        };
        
        clauses = new ANDedClauses(
            SKYCONDITIONAL.ALL, 
            new PARAM[]{PARAM.ABSOLUTE_CONTRAST}, 
            new PARAM[]{PARAM.STDEV_CONTRAST}, 
            new COMPARISON[]{COMPARISON.GREATER_THAN}, 
            new float[]{3.5f});
        
        clauses.setACustomCoefficient(0, cc);
        
        assertFalse(clauses.evaluate(data));
    }
}
