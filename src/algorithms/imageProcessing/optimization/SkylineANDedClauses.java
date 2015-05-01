package algorithms.imageProcessing.optimization;

import algorithms.util.PairFloat;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Map.Entry;

/**
 *
 * @author nichole
 */
public class SkylineANDedClauses {
    
    private float[] allSkiesCoeff = new float[] {
        //00   01
        0.0f, 10.0f,
        //02   03     04    05     06     07
        0.0f, 0.1f,  1.5f,  0.5f, -2.0f,  2.5f,
        //08    09     10    11       12     13    14    15
        0.005f, 5.0f, 5.0f,  0.001f, 1.5f, 0.001f, 1.5f, 1.0f
    };
    
    private float[] redSkiesCoeff = new float[] {
        //00  01     02     03
        0.0f, 15.0f, 0.03f, 15.0f,
        //04  05     06     07
        0.0f, 15.0f, 0.03f, 15.0f
    };
        
    private float[] blueSkiesCoeff = new float[] {
        //00  01     02    03    04     05   06
        0.0f, 0.05f, 1.5f, 3.0f, 0.37f, 199, 199,
        //07   08     09   10    11    12
        //0.0f, 0.05f, 1.5f, 0.9f, 0.9f, 0.0f,
        0.0f, 0.05f, 1.5f, 0.6f, 0.9f, 0.0f,
        
        //13   14      15      16     17     18     19   20
        0.05f, 0.005f, 0.005f, 0.08f, 0.03f, 0.03f, 199, 199
    };
    
    static protected final float[][] allLowerLimits = new float[6][];
    static protected final float[][] allUpperLimits = new float[6][];
    
    static protected final float[][] redLowerLimits = new float[2][];
    static protected final float[][] redUpperLimits = new float[2][];
    
    static protected final float[][] blueLowerLimits = new float[7][];
    static protected final float[][] blueUpperLimits = new float[7][];
    
    static protected final Map<Integer, Map<Integer, 
        Map<Integer, Float>>>  allCustomCoeffLowerLimits 
        = new HashMap<Integer, Map<Integer, Map<Integer, Float>>>();
    static protected final Map<Integer, Map<Integer, 
        Map<Integer, Float>>>  allCustomCoeffUpperLimits 
        = new HashMap<Integer, Map<Integer, Map<Integer, Float>>>();
    
    //key = clause number, value = map of coefficient number, and value=coeff
    static protected final Map<Integer, Map<Integer, 
        Map<Integer, Float>>> redCustomCoeffLowerLimits 
        = new HashMap<Integer, Map<Integer, Map<Integer, Float>>>();
    static protected final Map<Integer, Map<Integer, 
        Map<Integer, Float>>> redCustomCoeffUpperLimits 
        = new HashMap<Integer, Map<Integer, Map<Integer, Float>>>();
    
    static protected final Map<Integer, Map<Integer, 
        Map<Integer, Float>>> blueCustomCoeffLowerLimits 
        = new HashMap<Integer, Map<Integer, Map<Integer, Float>>>();
    static protected final Map<Integer, Map<Integer, 
        Map<Integer, Float>>> blueCustomCoeffUpperLimits 
        = new HashMap<Integer, Map<Integer, Map<Integer, Float>>>();
    
    public SkylineANDedClauses() {
        
        //                                00       01           
        allLowerLimits[0] = new float[]{0.00001f, 2.5f};
        allUpperLimits[0] = new float[]{0.2f,     30.f};
        //                                02       03     custom  07
        allLowerLimits[1] = new float[]{0.00001f, 0.01f, -14.5f,  1.1f};
        allUpperLimits[1] = new float[]{0.2f,     0.5f,   1.001f, 5.0f};
        
        // custom coefficient limits are stored by outerClauseIndex 
        //    (same as used by allLowerLimits for example)
        //    then by innerClauseIndex, then by coefficentIndex as the x in PairFloat
        allCustomCoeffLowerLimits.put(Integer.valueOf(1), 
            new HashMap<Integer, Map<Integer, Float>>());
        allCustomCoeffUpperLimits.put(Integer.valueOf(1), 
            new HashMap<Integer, Map<Integer, Float>>());
        allCustomCoeffLowerLimits.get(Integer.valueOf(1))
            .put(Integer.valueOf(2), new HashMap<Integer, Float>());
        allCustomCoeffUpperLimits.get(Integer.valueOf(1))
            .put(Integer.valueOf(2), new HashMap<Integer, Float>());
        
        allCustomCoeffLowerLimits.get(Integer.valueOf(1)).get(Integer.valueOf(2))
            .put(Integer.valueOf(4), 0.2f);
        allCustomCoeffUpperLimits.get(Integer.valueOf(1)).get(Integer.valueOf(2))
            .put(Integer.valueOf(4), 10.0f);
        
        allCustomCoeffLowerLimits.get(Integer.valueOf(1)).get(Integer.valueOf(2))
            .put(Integer.valueOf(5), 0.05f);
        allCustomCoeffUpperLimits.get(Integer.valueOf(1)).get(Integer.valueOf(2))
            .put(Integer.valueOf(5), 0.95f);
        
        allCustomCoeffLowerLimits.get(Integer.valueOf(1)).get(Integer.valueOf(2))
            .put(Integer.valueOf(6), -0.5f);
        allCustomCoeffUpperLimits.get(Integer.valueOf(1)).get(Integer.valueOf(2))
            .put(Integer.valueOf(6), -10.0f);
        
        //                                08       09     10    11         15
        allLowerLimits[2] = new float[]{0.00001f, 1.1f, 1.1f,  0.00001f,   0.2f};
        allUpperLimits[2] = new float[]{0.5f,     10.f, 10.0f, 0.05f,     100.f};
        
        //                                08       09     10     12   15
        allLowerLimits[3] = new float[]{0.00001f, 1.1f, 1.1f,   1.1f, 0.2f};
        allUpperLimits[3] = new float[]{0.5f,     10.f, 10.0f,  5.0f, 100.f};
        
         //                                08       09     10      13      15
        allLowerLimits[4] = new float[]{0.00001f, 1.1f, 1.1f,  0.00001f, 0.2f};
        allUpperLimits[4] = new float[]{0.5f,     10.f, 10.0f,  0.05f,   100.f};
        
        //                                08       09     10    14    15
        allLowerLimits[5] = new float[]{0.00001f, 1.1f, 1.1f,  1.1f, 0.2f};
        allUpperLimits[5] = new float[]{0.5f,     10.f, 10.0f,  5.0f, 100.f};
        
        
        //                                00       01    02        03
        redLowerLimits[0] = new float[]{0.00001f, 0.01f, 0.00001f, 1.5f};
        redUpperLimits[0] = new float[]{0.2f,     30.f,  0.2f,     30.0f};
        
        // clauses 0, 
        redCustomCoeffLowerLimits.put(Integer.valueOf(0), 
            new HashMap<Integer, Map<Integer, Float>>());
        redCustomCoeffUpperLimits.put(Integer.valueOf(0), 
            new HashMap<Integer, Map<Integer, Float>>());
        
        redCustomCoeffLowerLimits.get(Integer.valueOf(0))
            .put(Integer.valueOf(1), new HashMap<Integer, Float>());
        redCustomCoeffUpperLimits.get(Integer.valueOf(0))
            .put(Integer.valueOf(1), new HashMap<Integer, Float>());
        redCustomCoeffLowerLimits.get(Integer.valueOf(0))
            .put(Integer.valueOf(3), new HashMap<Integer, Float>());
        redCustomCoeffUpperLimits.get(Integer.valueOf(0))
            .put(Integer.valueOf(3), new HashMap<Integer, Float>());
        
        redCustomCoeffLowerLimits.get(Integer.valueOf(0)).get(Integer.valueOf(1))
            .put(Integer.valueOf(1), 4.0f);
        redCustomCoeffUpperLimits.get(Integer.valueOf(0)).get(Integer.valueOf(1))
            .put(Integer.valueOf(1), 30.0f);
        
        redCustomCoeffLowerLimits.get(Integer.valueOf(0)).get(Integer.valueOf(3))
            .put(Integer.valueOf(3), 4.0f);
        
        redCustomCoeffUpperLimits.get(Integer.valueOf(0)).get(Integer.valueOf(3))
            .put(Integer.valueOf(3), 30.0f);
        
        //                                04       05    06        07
        redLowerLimits[1] = new float[]{0.00001f, 0.01f, 0.00001f, 1.5f};
        redUpperLimits[1] = new float[]{0.2f,     30.f,  0.2f,     30.0f};
        
        redCustomCoeffLowerLimits.put(Integer.valueOf(1), 
            new HashMap<Integer, Map<Integer, Float>>());
        redCustomCoeffUpperLimits.put(Integer.valueOf(1), 
            new HashMap<Integer, Map<Integer, Float>>());
        
        redCustomCoeffLowerLimits.get(Integer.valueOf(1))
            .put(Integer.valueOf(1), new HashMap<Integer, Float>());
        redCustomCoeffUpperLimits.get(Integer.valueOf(1))
            .put(Integer.valueOf(1), new HashMap<Integer, Float>());
        
        redCustomCoeffLowerLimits.get(Integer.valueOf(1))
            .put(Integer.valueOf(3), new HashMap<Integer, Float>());
        redCustomCoeffUpperLimits.get(Integer.valueOf(1))
            .put(Integer.valueOf(3), new HashMap<Integer, Float>());
        
        redCustomCoeffLowerLimits.get(Integer.valueOf(1)).get(Integer.valueOf(1))
            .put(Integer.valueOf(5), 4.0f);
        redCustomCoeffUpperLimits.get(Integer.valueOf(1)).get(Integer.valueOf(1))
            .put(Integer.valueOf(5), 30.0f);
        
        redCustomCoeffLowerLimits.get(Integer.valueOf(1)).get(Integer.valueOf(3))
            .put(Integer.valueOf(7), 4.0f);
        redCustomCoeffUpperLimits.get(Integer.valueOf(1)).get(Integer.valueOf(3))
            .put(Integer.valueOf(7), 30.0f);
        
        //                                00       01       02    03     04    05   06
        blueLowerLimits[0] = new float[]{0.00001f, 0.0000f, 1.5f, 1.5f,  0.3f, 100, 100};
        blueUpperLimits[0] = new float[]{0.2f,     0.5f,    30.f, 30.0f, 0.4f, 250, 250};
        //                                07       08       09    10     11     12
        blueLowerLimits[1] = new float[]{0.00001f, 0.0000f, 1.5f, 0.75f, 0.75f, 0.00001f};
        blueUpperLimits[1] = new float[]{0.2f,     0.5f,    30.f, 2.0f,  2.0f,  3.0f};
        
        //                                13       14       15       16   
        blueLowerLimits[2] = new float[]{0.00001f, 0.0000f, 0.0000f, 0.01f};
        blueUpperLimits[2] = new float[]{0.2f,     0.5f,    0.5f,    0.30f};
        
        //                                13       14       15        17  
        blueLowerLimits[3] = new float[]{0.00001f, 0.0000f, 0.0000f, 0.01f};
        blueUpperLimits[3] = new float[]{0.2f,     0.5f,    0.5f,    0.30f};
        
        //                                13       14       15         18  
        blueLowerLimits[4] = new float[]{0.00001f, 0.0000f, 0.0000f, 0.01f};
        blueUpperLimits[4] = new float[]{0.2f,     0.5f,    0.5f,    0.30f};
        
        //                                13       14       15        19 
        blueLowerLimits[5] = new float[]{0.00001f, 0.0000f, 0.0000f, 100};
        blueUpperLimits[5] = new float[]{0.2f,     0.5f,    0.5f,    250};
        
        //                                13       14       15       20
        blueLowerLimits[6] = new float[]{0.00001f, 0.0000f, 0.0000f, 100};
        blueUpperLimits[6] = new float[]{0.2f,     0.5f,    0.5f,    250};
  
    }
    
    public ANDedClauses[] getAllClauses() {
        
        ANDedClauses[] a0 = getForAllSkies();
        ANDedClauses[] a1 = getForBlueSkies();
        ANDedClauses[] a2 = getForRedSkies();
        
        int n = a0.length + a1.length + a2.length;
    
        ANDedClauses[] a = new ANDedClauses[n];
        
        System.arraycopy(a0, 0, a, 0, a0.length);
        System.arraycopy(a1, 0, a, a0.length, a1.length);
        System.arraycopy(a2, 0, a, (a0.length + a1.length), a2.length);
        
        return a;
    }
    
    public ANDedClauses[] getAllAndRedClauses() {
        
        ANDedClauses[] a0 = getForAllSkies();
        ANDedClauses[] a2 = getForRedSkies();
        
        int n = a0.length + a2.length;
    
        ANDedClauses[] a = new ANDedClauses[n];
        
        System.arraycopy(a0, 0, a, 0, a0.length);
        System.arraycopy(a2, 0, a, a0.length, a2.length);
        
        return a;
    }
    
    public ANDedClauses[] getAllAndBlueClauses() {
        
        ANDedClauses[] a0 = getForAllSkies();
        ANDedClauses[] a2 = getForBlueSkies();
        
        int n = a0.length + a2.length;
    
        ANDedClauses[] a = new ANDedClauses[n];
        
        System.arraycopy(a0, 0, a, 0, a0.length);
        System.arraycopy(a2, 0, a, a0.length, a2.length);
        
        return a;
    }
    
    public ANDedClauses[] getForAllSkies() {
    
        /*
        c0:
            00         (skyStDevContrast > 0.)
            01         && ((Math.abs(contrastV)/skyStDevContrast) > 10.)
        */
        ANDedClauses c0 = new ANDedClauses(2, SKYCONDITIONAL.ALL);
        c0.set(0, PARAM.STDEV_CONTRAST, PARAM.INT_ONE, 
            COMPARISON.GREATER_THAN, allSkiesCoeff[0]);
        c0.set(1, PARAM.ABSOLUTE_CONTRAST, PARAM.STDEV_CONTRAST, 
            COMPARISON.GREATER_THAN, allSkiesCoeff[1]);
        
        /*
        c1      
            02         (skyStDevContrast > 0.)
            03         && ((Math.abs(contrastV) > 0.1)
            04,05,06   && ((Math.abs(contrastV)/skyStDevContrast) > (1.5 + (Math.abs(contrastV)-0.5)*(-2.0))))  <-------
            07         && ((Math.abs(colorDiffV)/skyStDevColorDiff) > 2.5)
        */
        ANDedClauses c1 = new ANDedClauses(4, SKYCONDITIONAL.ALL);
        c1.set(0, PARAM.STDEV_CONTRAST, PARAM.INT_ONE, 
            COMPARISON.GREATER_THAN, allSkiesCoeff[2]);
        c1.set(1, PARAM.ABSOLUTE_CONTRAST, PARAM.INT_ONE, 
            COMPARISON.GREATER_THAN, allSkiesCoeff[3]);
        
        c1.set(2, PARAM.ABSOLUTE_CONTRAST, PARAM.STDEV_CONTRAST, 
            COMPARISON.GREATER_THAN, allSkiesCoeff[4]);
        c1.setACustomCoefficient(2, new CustomCoeff00()); 
        c1.setCustomCoefficientVariable(4, Float.valueOf(allSkiesCoeff[4]));
        c1.setCustomCoefficientVariable(5, Float.valueOf(allSkiesCoeff[5]));
        c1.setCustomCoefficientVariable(6, Float.valueOf(allSkiesCoeff[6]));
        
        c1.set(3, PARAM.ABSOLUTE_DIFF_BLUE_OR_RED, PARAM.STDEV_BLUE_OR_RED, 
            COMPARISON.GREATER_THAN, allSkiesCoeff[7]);
       
        /* 
        c3,c4,c5,c6        
                    08         (skyStDevContrast > 0.005)
                    09         && ((Math.abs(contrastV)/skyStDevContrast) > 5.)
                    10         && ((Math.abs(colorDiffV)/skyStDevColorDiff) > 5.)
                               &&
                               // if cieXY diffs are zero and within stdev, these are sky,
                               // so test for opposite for boundary pixel
                    11         (((diffCIEX > 0.001)
                    12             || ((diffCIEX/localSky.getStdDevCIEX()) > 1.5)
                    13             || (diffCIEY > 0.001)
                    14             || ((diffCIEY/localSky.getStdDevCIEY()) > 1.5)))
                    15         && (skyStDevColorDiff > 1.)
        */
        ANDedClauses c2 = new ANDedClauses(5, SKYCONDITIONAL.ALL);
        c2.set(0, PARAM.STDEV_CONTRAST, PARAM.INT_ONE, 
            COMPARISON.GREATER_THAN, allSkiesCoeff[8]);
        c2.set(1, PARAM.ABSOLUTE_CONTRAST, PARAM.STDEV_CONTRAST, 
            COMPARISON.GREATER_THAN, allSkiesCoeff[9]);
        c2.set(2, PARAM.ABSOLUTE_DIFF_BLUE_OR_RED, PARAM.STDEV_BLUE_OR_RED, 
            COMPARISON.GREATER_THAN, allSkiesCoeff[10]);
    c2.set(3, PARAM.DIFF_CIEX, PARAM.INT_ONE, 
            COMPARISON.GREATER_THAN, allSkiesCoeff[11]);
        c2.set(4, PARAM.STDEV_BLUE_OR_RED, PARAM.INT_ONE, 
            COMPARISON.GREATER_THAN, allSkiesCoeff[15]);
        
        ANDedClauses c3 = new ANDedClauses(5, SKYCONDITIONAL.ALL);
        c3.set(0, PARAM.STDEV_CONTRAST, PARAM.INT_ONE, 
            COMPARISON.GREATER_THAN, allSkiesCoeff[8]);
        c3.set(1, PARAM.ABSOLUTE_CONTRAST, PARAM.STDEV_CONTRAST, 
            COMPARISON.GREATER_THAN, allSkiesCoeff[9]);
        c3.set(2, PARAM.ABSOLUTE_DIFF_BLUE_OR_RED, PARAM.STDEV_BLUE_OR_RED, 
            COMPARISON.GREATER_THAN, allSkiesCoeff[10]);
    c3.set(3, PARAM.DIFF_CIEX, PARAM.STDEV_CIEX, 
            COMPARISON.GREATER_THAN, allSkiesCoeff[12]);
        c3.set(4, PARAM.STDEV_BLUE_OR_RED, PARAM.INT_ONE, 
            COMPARISON.GREATER_THAN, allSkiesCoeff[15]);
        
        ANDedClauses c4 = new ANDedClauses(5, SKYCONDITIONAL.ALL);
        c4.set(0, PARAM.STDEV_CONTRAST, PARAM.INT_ONE, 
            COMPARISON.GREATER_THAN, allSkiesCoeff[8]);
        c4.set(1, PARAM.ABSOLUTE_CONTRAST, PARAM.STDEV_CONTRAST, 
            COMPARISON.GREATER_THAN, allSkiesCoeff[9]);
        c4.set(2, PARAM.ABSOLUTE_DIFF_BLUE_OR_RED, PARAM.STDEV_BLUE_OR_RED, 
            COMPARISON.GREATER_THAN, allSkiesCoeff[10]);
    c4.set(3, PARAM.DIFF_CIEY, PARAM.INT_ONE, 
            COMPARISON.GREATER_THAN, allSkiesCoeff[13]);
        c4.set(4, PARAM.STDEV_BLUE_OR_RED, PARAM.INT_ONE, 
            COMPARISON.GREATER_THAN, allSkiesCoeff[15]);
       
        ANDedClauses c5 = new ANDedClauses(5, SKYCONDITIONAL.ALL);
        c5.set(0, PARAM.STDEV_CONTRAST, PARAM.INT_ONE, 
            COMPARISON.GREATER_THAN, allSkiesCoeff[8]);
        c5.set(1, PARAM.ABSOLUTE_CONTRAST, PARAM.STDEV_CONTRAST, 
            COMPARISON.GREATER_THAN, allSkiesCoeff[9]);
        c5.set(2, PARAM.ABSOLUTE_DIFF_BLUE_OR_RED, PARAM.STDEV_BLUE_OR_RED, 
            COMPARISON.GREATER_THAN, allSkiesCoeff[10]);
    c5.set(3, PARAM.DIFF_CIEY, PARAM.STDEV_CIEY, 
            COMPARISON.GREATER_THAN, allSkiesCoeff[14]);
        c5.set(4, PARAM.STDEV_BLUE_OR_RED, PARAM.INT_ONE, 
            COMPARISON.GREATER_THAN, allSkiesCoeff[15]);
      
        return new ANDedClauses[]{c0, c1, c2, c3, c4, c5};
    }
    
    public float[] getFittableCoefficientsForAllSkies() {
        return allSkiesCoeff;
    }
    
    public ANDedClauses[] getForRedSkies() {
    
        /*
        00         (skyStDevContrast > 0.)
        01         && ((Math.abs(colorDiffV)/skyStDevColorDiff) > 15.*diffCIEX) <---
        02         && (diffCIEX > 0.03)
        03         && ((diffCIEX/localSky.getStdDevCIEX()) > 15.*diffCIEX) <---
        */
        ANDedClauses c0 = new ANDedClauses(4, SKYCONDITIONAL.RED);
        c0.set(0, PARAM.STDEV_CONTRAST, PARAM.INT_ONE, 
            COMPARISON.GREATER_THAN, redSkiesCoeff[0]);
        c0.set(1, PARAM.ABSOLUTE_DIFF_BLUE_OR_RED, 
            PARAM.STDEV_BLUE_OR_RED, COMPARISON.GREATER_THAN, redSkiesCoeff[1]);
        
        c0.setACustomCoefficient(1, new CustomCoeff01());
        c0.setCustomCoefficientVariable(1, Float.valueOf(redSkiesCoeff[1]));
        
        c0.set(2, PARAM.DIFF_CIEX, 
            PARAM.INT_ONE, COMPARISON.GREATER_THAN, redSkiesCoeff[2]);
        c0.set(3, PARAM.DIFF_CIEX, 
            PARAM.STDEV_CIEX, COMPARISON.GREATER_THAN, redSkiesCoeff[3]);
        c0.setACustomCoefficient(3, new CustomCoeff02());
        c0.setCustomCoefficientVariable(3, Float.valueOf(redSkiesCoeff[3]));
       
        /*
        04         (skyStDevContrast > 0.)
        05         && ((Math.abs(colorDiffV)/skyStDevColorDiff) > 15.*diffCIEY) <---
        06         && (diffCIEY > 0.03)
        07         && ((diffCIEY/localSky.getStdDevCIEY()) > 15.*diffCIEY) <---
        */
        ANDedClauses c1 = new ANDedClauses(4, SKYCONDITIONAL.RED);
        c1.set(0, PARAM.STDEV_CONTRAST, PARAM.INT_ONE, 
            COMPARISON.GREATER_THAN, redSkiesCoeff[4]);
        c1.set(1, PARAM.ABSOLUTE_DIFF_BLUE_OR_RED, 
            PARAM.STDEV_BLUE_OR_RED, COMPARISON.GREATER_THAN, redSkiesCoeff[5]);
        
        c1.setACustomCoefficient(1, new CustomCoeff03());
        c1.setCustomCoefficientVariable(5, Float.valueOf(redSkiesCoeff[5]));
        
        c1.set(2, PARAM.DIFF_CIEY, 
            PARAM.INT_ONE, COMPARISON.GREATER_THAN, redSkiesCoeff[6]);
        c1.set(3, PARAM.DIFF_CIEY, 
            PARAM.STDEV_CIEY, COMPARISON.GREATER_THAN, redSkiesCoeff[7]);
        
        c1.setACustomCoefficient(3, new CustomCoeff04());
        c1.setCustomCoefficientVariable(7, Float.valueOf(redSkiesCoeff[7]));
        
        return new ANDedClauses[]{c0, c1};
    }
    
    public float[] getFittableCoefficientsForRedSkies() {
        return redSkiesCoeff;
    }
    
    public ANDedClauses[] getForBlueSkies() {
        
        /*
        00         (skyStDevContrast > 0.0)
        01         && (contrastV < 0.05)
        02         && ((Math.abs(contrastV)/skyStDevContrast) > 1.5)
        03         && ((Math.abs(colorDiffV)/skyStDevColorDiff) > 3.0)
        04         && ((bPercentV < 0.37)
        05         && (bV > 199)
        06         && (gV > 199))
        */
        ANDedClauses c0 = new ANDedClauses(7, SKYCONDITIONAL.BLUE);
        c0.set(0, PARAM.STDEV_CONTRAST, PARAM.INT_ONE, 
            COMPARISON.GREATER_THAN, blueSkiesCoeff[0]);
        c0.set(1, PARAM.CONTRAST, 
            PARAM.INT_ONE, COMPARISON.LESS_THAN, blueSkiesCoeff[1]);
        
        c0.set(2, PARAM.ABSOLUTE_CONTRAST, 
            PARAM.STDEV_CONTRAST, COMPARISON.GREATER_THAN, blueSkiesCoeff[2]);
        
        c0.set(3, PARAM.ABSOLUTE_DIFF_BLUE_OR_RED, 
            PARAM.STDEV_BLUE_OR_RED, COMPARISON.GREATER_THAN, blueSkiesCoeff[3]);
        
        c0.set(4, PARAM.B_DIV_TOT, 
            PARAM.INT_ONE, COMPARISON.LESS_THAN, blueSkiesCoeff[4]);
        
        c0.set(5, PARAM.BLUE, 
            PARAM.INT_ONE, COMPARISON.GREATER_THAN, blueSkiesCoeff[5]);
        
        c0.set(6, PARAM.GREEN, 
            PARAM.INT_ONE, COMPARISON.GREATER_THAN, blueSkiesCoeff[6]);
       
        /*
        07         (skyStDevContrast > 0.0)
        08         && (Math.abs(contrastV) > 0.05)
        09         && ((Math.abs(colorDiffV)/skyStDevColorDiff) > 1.5)
        10         //&& ((diffCIEX/localSky.getStdDevCIEX()) > 0.9) // ?
                   ((diffCIEX/localSky.getStdDevCIEX()) > 0.6)
        11         && ((diffCIEY/localSky.getStdDevCIEY()) > 0.9) // ?
        12         && (skyStDevColorDiff > 0.)
        */
        ANDedClauses c1 = new ANDedClauses(6, SKYCONDITIONAL.BLUE);
        c1.set(0, PARAM.STDEV_CONTRAST, PARAM.INT_ONE, 
            COMPARISON.GREATER_THAN, blueSkiesCoeff[7]);
        
        c1.set(1, PARAM.ABSOLUTE_CONTRAST, 
            PARAM.INT_ONE, COMPARISON.GREATER_THAN, blueSkiesCoeff[8]);
        
        c1.set(2, PARAM.ABSOLUTE_DIFF_BLUE_OR_RED, 
            PARAM.STDEV_BLUE_OR_RED, COMPARISON.GREATER_THAN, blueSkiesCoeff[9]);
        
        c1.set(3, PARAM.DIFF_CIEX, 
            PARAM.STDEV_CIEX, COMPARISON.GREATER_THAN, blueSkiesCoeff[10]);
        
        c1.set(4, PARAM.DIFF_CIEY, 
            PARAM.STDEV_CIEY, COMPARISON.GREATER_THAN, blueSkiesCoeff[11]);
        
        c1.set(5, PARAM.STDEV_BLUE_OR_RED, 
            PARAM.INT_ONE, COMPARISON.GREATER_THAN, blueSkiesCoeff[12]);
               
        /*
        13         (contrastV < 0.05)
        14         && (diffCIEX < 0.005)
        15         && (diffCIEY < 0.005)
                   &&
        16         ((Math.abs(0.33 - rPercentV) > 0.08)
        17         || (Math.abs(0.33 - gPercentV) > 0.03)
        18         || (Math.abs(0.33 - bPercentV) > 0.03)
        19         || (gV > 199)
        20         || (bV > 199))
        */
        ANDedClauses c2 = new ANDedClauses(4, SKYCONDITIONAL.BLUE);
        c2.set(0, PARAM.CONTRAST, PARAM.INT_ONE, 
            COMPARISON.LESS_THAN, blueSkiesCoeff[13]);
        
        c2.set(1, PARAM.DIFF_CIEX, PARAM.INT_ONE, 
            COMPARISON.LESS_THAN, blueSkiesCoeff[14]);
        
        c2.set(2, PARAM.DIFF_CIEY, PARAM.INT_ONE, 
            COMPARISON.LESS_THAN, blueSkiesCoeff[15]);
        
        c2.set(3, PARAM.DIFF_R_DIV_TOT_ONE_THIRD, PARAM.INT_ONE, 
            COMPARISON.GREATER_THAN, blueSkiesCoeff[16]);
       
        
        ANDedClauses c3 = new ANDedClauses(4, SKYCONDITIONAL.BLUE);
        c3.set(0, PARAM.CONTRAST, PARAM.INT_ONE, 
            COMPARISON.LESS_THAN, blueSkiesCoeff[13]);
        
        c3.set(1, PARAM.DIFF_CIEX, PARAM.INT_ONE, 
            COMPARISON.LESS_THAN, blueSkiesCoeff[14]);
        
        c3.set(2, PARAM.DIFF_CIEY, PARAM.INT_ONE, 
            COMPARISON.LESS_THAN, blueSkiesCoeff[15]);
        
        c3.set(3, PARAM.DIFF_G_DIV_TOT_ONE_THIRD, PARAM.INT_ONE, 
            COMPARISON.GREATER_THAN, blueSkiesCoeff[17]);
        
        
        ANDedClauses c4 = new ANDedClauses(4, SKYCONDITIONAL.BLUE);
        c4.set(0, PARAM.CONTRAST, PARAM.INT_ONE, 
            COMPARISON.LESS_THAN, blueSkiesCoeff[13]);
        
        c4.set(1, PARAM.DIFF_CIEX, PARAM.INT_ONE, 
            COMPARISON.LESS_THAN, blueSkiesCoeff[14]);
        
        c4.set(2, PARAM.DIFF_CIEY, PARAM.INT_ONE, 
            COMPARISON.LESS_THAN, blueSkiesCoeff[15]);
        
        c4.set(3, PARAM.DIFF_B_DIV_TOT_ONE_THIRD, PARAM.INT_ONE, 
            COMPARISON.GREATER_THAN, blueSkiesCoeff[18]);
        
   
        ANDedClauses c5 = new ANDedClauses(4, SKYCONDITIONAL.BLUE);
        c5.set(0, PARAM.CONTRAST, PARAM.INT_ONE, 
            COMPARISON.LESS_THAN, blueSkiesCoeff[13]);
        
        c5.set(1, PARAM.DIFF_CIEX, PARAM.INT_ONE, 
            COMPARISON.LESS_THAN, blueSkiesCoeff[14]);
        
        c5.set(2, PARAM.DIFF_CIEY, PARAM.INT_ONE, 
            COMPARISON.LESS_THAN, blueSkiesCoeff[15]);
        
        c5.set(3, PARAM.GREEN, PARAM.INT_ONE, 
            COMPARISON.GREATER_THAN, blueSkiesCoeff[19]);
        
        
        ANDedClauses c6 = new ANDedClauses(4, SKYCONDITIONAL.BLUE);
        c6.set(0, PARAM.CONTRAST, PARAM.INT_ONE, 
            COMPARISON.LESS_THAN, blueSkiesCoeff[13]);
        
        c6.set(1, PARAM.DIFF_CIEX, PARAM.INT_ONE, 
            COMPARISON.LESS_THAN, blueSkiesCoeff[14]);
        
        c6.set(2, PARAM.DIFF_CIEY, PARAM.INT_ONE, 
            COMPARISON.LESS_THAN, blueSkiesCoeff[15]);
        
        c6.set(3, PARAM.BLUE, PARAM.INT_ONE, 
            COMPARISON.GREATER_THAN, blueSkiesCoeff[20]);
        
        return new ANDedClauses[]{c0, c1, c2, c3, c4, c5, c6};
    }
    
    public float[] getFittableCoefficientsForBlueSkies() {
        return blueSkiesCoeff;
    }
   
    public float[][] getAllCoeffLowerLimits() {
        
        float[][] a0 = allLowerLimits;
        float[][] a1 = blueLowerLimits;
        float[][] a2 = redLowerLimits;
        
        int n = a0.length + a1.length + a2.length;
    
        int count = 0;
        float[][] a = new float[n][];
        for (int i = 0; i < a0.length; i++) {
            a[count] = Arrays.copyOf(a0[i], a0[i].length);
            count++;
        }
        for (int i = 0; i < a1.length; i++) {
            a[count] = Arrays.copyOf(a1[i], a1[i].length);
            count++;
        }
        for (int i = 0; i < a2.length; i++) {
            a[count] = Arrays.copyOf(a2[i], a2[i].length);
            count++;
        }
        
        return a;
    }
    
    public float[][] getAllCoeffUpperLimits() {
        
        float[][] a0 = allUpperLimits;
        float[][] a1 = blueUpperLimits;
        float[][] a2 = redUpperLimits;
        
        int n = a0.length + a1.length + a2.length;
    
        int count = 0;
        float[][] a = new float[n][];
        for (int i = 0; i < a0.length; i++) {
            a[count] = Arrays.copyOf(a0[i], a0[i].length);
            count++;
        }
        for (int i = 0; i < a1.length; i++) {
            a[count] = Arrays.copyOf(a1[i], a1[i].length);
            count++;
        }
        for (int i = 0; i < a2.length; i++) {
            a[count] = Arrays.copyOf(a2[i], a2[i].length);
            count++;
        }
        
        return a;
    }
    
    public float[][] getAllAndRedCoeffLowerLimits() {
        
        float[][] a0 = allLowerLimits;
        float[][] a1 = redLowerLimits;
        
        int n = a0.length + a1.length;
    
        int count = 0;
        float[][] a = new float[n][];
        for (int i = 0; i < a0.length; i++) {
            a[count] = Arrays.copyOf(a0[i], a0[i].length);
            count++;
        }
        for (int i = 0; i < a1.length; i++) {
            a[count] = Arrays.copyOf(a1[i], a1[i].length);
            count++;
        }
        
        return a;
    }
    
    public float[][] getAllAndRedCoeffUpperLimits() {
        
        float[][] a0 = allUpperLimits;
        float[][] a1 = redUpperLimits;
        
        int n = a0.length + a1.length;
    
        int count = 0;
        float[][] a = new float[n][];
        for (int i = 0; i < a0.length; i++) {
            a[count] = Arrays.copyOf(a0[i], a0[i].length);
            count++;
        }
        for (int i = 0; i < a1.length; i++) {
            a[count] = Arrays.copyOf(a1[i], a1[i].length);
            count++;
        }
        
        return a;
    }
    
    /**
     * custom coefficient limits are stored by outerClauseIndex 
       (same as used by allLowerLimits for example)
       then by innerClauseIndex, then by coefficentIndex as the x in PairFloat
      
     * @return 
     */
    public Map<Integer, Map<Integer, Map<Integer, Float>>> getAllCustomCoeffLowerLimits() {
        
        int nAllCoeff = allLowerLimits.length;
        int nBlueCoeff = blueLowerLimits.length;
        
        // combine the maps, but any blue clause indexes will need nAllCoeff added
        
        Map<Integer, Map<Integer, Map<Integer, Float>>> combined 
            = new HashMap<Integer, Map<Integer, Map<Integer, Float>>>();
        
        if (!allCustomCoeffLowerLimits.isEmpty()) {
            
            Iterator<Entry<Integer, Map<Integer, Map<Integer, Float>>>> iter 
                = allCustomCoeffLowerLimits.entrySet().iterator();
            
            while (iter.hasNext()) {
                
                Entry<Integer, Map<Integer, Map<Integer, Float>>> entry = iter.next();
                
                Integer outerClauseIndex = entry.getKey();
                
                Map<Integer, Map<Integer, Float>> value = entry.getValue();
                
                Iterator<Entry<Integer, Map<Integer, Float>>> iter2 = value.entrySet().iterator();
                
                Map<Integer, Map<Integer, Float>> cMap = combined.get(outerClauseIndex);
                if (cMap == null) {
                    cMap = new HashMap<Integer, Map<Integer, Float>>();
                    combined.put(outerClauseIndex, cMap);
                }
                
                while (iter2.hasNext()) {
                    Entry<Integer, Map<Integer, Float>> entry2 = iter2.next();
                    
                    Integer innerClauseIndex = entry2.getKey();
                    
                    Map<Integer, Float> coeffIdxAndValue = entry2.getValue();
                    
                    Map<Integer, Float> cMap2 = cMap.get(innerClauseIndex);
                    
                    if (cMap2 == null) {
                        cMap2 = new HashMap<Integer, Float>();
                        cMap.put(innerClauseIndex, cMap2);
                    }
                    
                    cMap2.putAll(coeffIdxAndValue);
                }
            }
        }
        
        if (!blueCustomCoeffLowerLimits.isEmpty()) {
            
            Iterator<Entry<Integer, Map<Integer, Map<Integer, Float>>>> iter 
                = blueCustomCoeffLowerLimits.entrySet().iterator();
            
            while (iter.hasNext()) {
                
                Entry<Integer, Map<Integer, Map<Integer, Float>>> entry = iter.next();
                
                Integer outerClauseIndex = entry.getKey();
                Integer outerClauseIndexForList = Integer.valueOf(
                    outerClauseIndex.intValue() + nAllCoeff);
                
                Map<Integer, Map<Integer, Float>> value = entry.getValue();
                
                Map<Integer, Map<Integer, Float>> cMap = combined.get(outerClauseIndexForList);
                if (cMap == null) {
                    cMap = new HashMap<Integer, Map<Integer, Float>>();
                    combined.put(outerClauseIndexForList, cMap);
                }
                
                Iterator<Entry<Integer, Map<Integer, Float>>> iter2 = value.entrySet().iterator();

                while (iter2.hasNext()) {
                    
                    Entry<Integer, Map<Integer, Float>> entry2 = iter2.next();
                    
                    Integer innerClauseIndex = entry2.getKey();
                    
                    Map<Integer, Float> coeffIdxAndValue = entry2.getValue();
                    
                    Map<Integer, Float> cMap2 = cMap.get(innerClauseIndex);
                    
                    if (cMap2 == null) {
                        cMap2 = new HashMap<Integer, Float>();
                        cMap.put(innerClauseIndex, cMap2);
                    }

                    cMap2.putAll(coeffIdxAndValue);
                }
            }
        }
        
        if (!redCustomCoeffLowerLimits.isEmpty()) {
            
            Iterator<Entry<Integer, Map<Integer, Map<Integer, Float>>>> iter 
                = redCustomCoeffLowerLimits.entrySet().iterator();
            
            while (iter.hasNext()) {
                
                Entry<Integer, Map<Integer, Map<Integer, Float>>> entry = iter.next();
                
                Integer outerClauseIndex = entry.getKey();
                Integer outerClauseIndexForList = Integer.valueOf(
                    outerClauseIndex.intValue() + nAllCoeff + nBlueCoeff);
                
                Map<Integer, Map<Integer, Float>> value = entry.getValue();
                
                Map<Integer, Map<Integer, Float>> cMap = combined.get(outerClauseIndexForList);
                if (cMap == null) {
                    cMap = new HashMap<Integer, Map<Integer, Float>>();
                    combined.put(outerClauseIndexForList, cMap);
                }
                
                Iterator<Entry<Integer, Map<Integer, Float>>> iter2 = value.entrySet().iterator();

                while (iter2.hasNext()) {
                    
                    Entry<Integer, Map<Integer, Float>> entry2 = iter2.next();
                    
                    Integer innerClauseIndex = entry2.getKey();
                    
                    Map<Integer, Float> coeffIdxAndValue = entry2.getValue();
                    
                    Map<Integer, Float> cMap2 = cMap.get(innerClauseIndex);
                    
                    if (cMap2 == null) {
                        cMap2 = new HashMap<Integer, Float>();
                        cMap.put(innerClauseIndex, cMap2);
                    }

                    cMap2.putAll(coeffIdxAndValue);
                }
            }
        }
        
        return combined;
    }
    
    /**
     * custom coefficient limits are stored by outerClauseIndex 
       (same as used by allLowerLimits for example)
       then by innerClauseIndex, then by coefficentIndex as the x in PairFloat
      
     * @return 
     */
    public Map<Integer, Map<Integer, Map<Integer, Float>>> getAllCustomCoeffUpperLimits() {
        
        int nAllCoeff = allUpperLimits.length;
        int nBlueCoeff = blueUpperLimits.length;
        
        // combine the maps, but any blue clause indexes will need nAllCoeff added
        
        Map<Integer, Map<Integer, Map<Integer, Float>>> combined 
            = new HashMap<Integer, Map<Integer, Map<Integer, Float>>>();
        
        if (!allCustomCoeffUpperLimits.isEmpty()) {
            
            Iterator<Entry<Integer, Map<Integer, Map<Integer, Float>>>> iter 
                = allCustomCoeffUpperLimits.entrySet().iterator();
            
            while (iter.hasNext()) {
                
                Entry<Integer, Map<Integer, Map<Integer, Float>>> entry = iter.next();
                
                Integer outerClauseIndex = entry.getKey();
                
                Map<Integer, Map<Integer, Float>> value = entry.getValue();
                
                Iterator<Entry<Integer, Map<Integer, Float>>> iter2 = value.entrySet().iterator();
                
                Map<Integer, Map<Integer, Float>> cMap = combined.get(outerClauseIndex);
                if (cMap == null) {
                    cMap = new HashMap<Integer, Map<Integer, Float>>();
                    combined.put(outerClauseIndex, cMap);
                }
                
                while (iter2.hasNext()) {
                    Entry<Integer, Map<Integer, Float>> entry2 = iter2.next();
                    
                    Integer innerClauseIndex = entry2.getKey();
                    
                    Map<Integer, Float> coeffIdxAndValue = entry2.getValue();
                    
                    Map<Integer, Float> cMap2 = cMap.get(innerClauseIndex);
                    
                    if (cMap2 == null) {
                        cMap2 = new HashMap<Integer, Float>();
                        cMap.put(innerClauseIndex, cMap2);
                    }
                    
                    cMap2.putAll(coeffIdxAndValue);
                }
            }
        }
        
        if (!blueCustomCoeffUpperLimits.isEmpty()) {
            
            Iterator<Entry<Integer, Map<Integer, Map<Integer, Float>>>> iter 
                = blueCustomCoeffUpperLimits.entrySet().iterator();
            
            while (iter.hasNext()) {
                
                Entry<Integer, Map<Integer, Map<Integer, Float>>> entry = iter.next();
                
                Integer outerClauseIndex = entry.getKey();
                Integer outerClauseIndexForList = Integer.valueOf(
                    outerClauseIndex.intValue() + nAllCoeff);
                
                Map<Integer, Map<Integer, Float>> value = entry.getValue();
                
                Map<Integer, Map<Integer, Float>> cMap = combined.get(outerClauseIndexForList);
                if (cMap == null) {
                    cMap = new HashMap<Integer, Map<Integer, Float>>();
                    combined.put(outerClauseIndexForList, cMap);
                }
                
                Iterator<Entry<Integer, Map<Integer, Float>>> iter2 = value.entrySet().iterator();

                while (iter2.hasNext()) {
                    
                    Entry<Integer, Map<Integer, Float>> entry2 = iter2.next();
                    
                    Integer innerClauseIndex = entry2.getKey();
                    
                    Map<Integer, Float> coeffIdxAndValue = entry2.getValue();
                    
                    Map<Integer, Float> cMap2 = cMap.get(innerClauseIndex);
                    
                    if (cMap2 == null) {
                        cMap2 = new HashMap<Integer, Float>();
                        cMap.put(innerClauseIndex, cMap2);
                    }

                    cMap2.putAll(coeffIdxAndValue);
                }
            }
        }
        
        if (!redCustomCoeffUpperLimits.isEmpty()) {
            
            Iterator<Entry<Integer, Map<Integer, Map<Integer, Float>>>> iter 
                = redCustomCoeffUpperLimits.entrySet().iterator();
            
            while (iter.hasNext()) {
                
                Entry<Integer, Map<Integer, Map<Integer, Float>>> entry = iter.next();
                
                Integer outerClauseIndex = entry.getKey();
                Integer outerClauseIndexForList = Integer.valueOf(
                    outerClauseIndex.intValue() + nAllCoeff + nBlueCoeff);
                
                Map<Integer, Map<Integer, Float>> value = entry.getValue();
                
                Map<Integer, Map<Integer, Float>> cMap = combined.get(outerClauseIndexForList);
                if (cMap == null) {
                    cMap = new HashMap<Integer, Map<Integer, Float>>();
                    combined.put(outerClauseIndexForList, cMap);
                }
                
                Iterator<Entry<Integer, Map<Integer, Float>>> iter2 = value.entrySet().iterator();

                while (iter2.hasNext()) {
                    
                    Entry<Integer, Map<Integer, Float>> entry2 = iter2.next();
                    
                    Integer innerClauseIndex = entry2.getKey();
                    
                    Map<Integer, Float> coeffIdxAndValue = entry2.getValue();
                    
                    Map<Integer, Float> cMap2 = cMap.get(innerClauseIndex);
                    
                    if (cMap2 == null) {
                        cMap2 = new HashMap<Integer, Float>();
                        cMap.put(innerClauseIndex, cMap2);
                    }

                    cMap2.putAll(coeffIdxAndValue);
                }
            }
        }
        
        return combined;
    }
    
    /**
     * custom coefficient limits are stored by outerClauseIndex 
       (same as used by allLowerLimits for example)
       then by innerClauseIndex, then by coefficentIndex as the x in PairFloat
      
     * @return 
     */
    public Map<Integer, Map<Integer, Map<Integer, Float>>> getAllAndRedCustomCoeffLowerLimits() {
        
        int nAllCoeff = allLowerLimits.length;
        
        // combine the maps, but any blue clause indexes will need nAllCoeff added
        
        Map<Integer, Map<Integer, Map<Integer, Float>>> combined 
            = new HashMap<Integer, Map<Integer, Map<Integer, Float>>>();
        
        if (!allCustomCoeffLowerLimits.isEmpty()) {
            
            Iterator<Entry<Integer, Map<Integer, Map<Integer, Float>>>> iter 
                = allCustomCoeffLowerLimits.entrySet().iterator();
            
            while (iter.hasNext()) {
                
                Entry<Integer, Map<Integer, Map<Integer, Float>>> entry = iter.next();
                
                Integer outerClauseIndex = entry.getKey();
                
                Map<Integer, Map<Integer, Float>> value = entry.getValue();
                
                Iterator<Entry<Integer, Map<Integer, Float>>> iter2 = value.entrySet().iterator();
                
                Map<Integer, Map<Integer, Float>> cMap = combined.get(outerClauseIndex);
                if (cMap == null) {
                    cMap = new HashMap<Integer, Map<Integer, Float>>();
                    combined.put(outerClauseIndex, cMap);
                }
                
                while (iter2.hasNext()) {
                    Entry<Integer, Map<Integer, Float>> entry2 = iter2.next();
                    
                    Integer innerClauseIndex = entry2.getKey();
                    
                    Map<Integer, Float> coeffIdxAndValue = entry2.getValue();
                    
                    Map<Integer, Float> cMap2 = cMap.get(innerClauseIndex);
                    
                    if (cMap2 == null) {
                        cMap2 = new HashMap<Integer, Float>();
                        cMap.put(innerClauseIndex, cMap2);
                    }
                    
                    cMap2.putAll(coeffIdxAndValue);
                }
            }
        }
        
        if (!redCustomCoeffLowerLimits.isEmpty()) {
            
            Iterator<Entry<Integer, Map<Integer, Map<Integer, Float>>>> iter 
                = redCustomCoeffLowerLimits.entrySet().iterator();
            
            while (iter.hasNext()) {
                
                Entry<Integer, Map<Integer, Map<Integer, Float>>> entry = iter.next();
                
                Integer outerClauseIndex = entry.getKey();
                Integer outerClauseIndexForList = Integer.valueOf(
                    outerClauseIndex.intValue() + nAllCoeff);
                
                Map<Integer, Map<Integer, Float>> value = entry.getValue();
                
                Map<Integer, Map<Integer, Float>> cMap = combined.get(outerClauseIndexForList);
                if (cMap == null) {
                    cMap = new HashMap<Integer, Map<Integer, Float>>();
                    combined.put(outerClauseIndexForList, cMap);
                }
                
                Iterator<Entry<Integer, Map<Integer, Float>>> iter2 = value.entrySet().iterator();

                while (iter2.hasNext()) {
                    
                    Entry<Integer, Map<Integer, Float>> entry2 = iter2.next();
                    
                    Integer innerClauseIndex = entry2.getKey();
                    
                    Map<Integer, Float> coeffIdxAndValue = entry2.getValue();
                    
                    Map<Integer, Float> cMap2 = cMap.get(innerClauseIndex);
                    
                    if (cMap2 == null) {
                        cMap2 = new HashMap<Integer, Float>();
                        cMap.put(innerClauseIndex, cMap2);
                    }

                    cMap2.putAll(coeffIdxAndValue);
                }
            }
        }
        
        return combined;
    }
    
    /**
     * custom coefficient limits are stored by outerClauseIndex 
       (same as used by allLowerLimits for example)
       then by innerClauseIndex, then by coefficentIndex as the x in PairFloat
      
     * @return 
     */
    public Map<Integer, Map<Integer, Map<Integer, Float>>> getAllAndRedCustomCoeffUpperLimits() {
        
        int nAllCoeff = allUpperLimits.length;
        
        // combine the maps, but any blue clause indexes will need nAllCoeff added
        
        Map<Integer, Map<Integer, Map<Integer, Float>>> combined 
            = new HashMap<Integer, Map<Integer, Map<Integer, Float>>>();
        
        if (!allCustomCoeffUpperLimits.isEmpty()) {
            
            Iterator<Entry<Integer, Map<Integer, Map<Integer, Float>>>> iter 
                = allCustomCoeffUpperLimits.entrySet().iterator();
            
            while (iter.hasNext()) {
                
                Entry<Integer, Map<Integer, Map<Integer, Float>>> entry = iter.next();
                
                Integer outerClauseIndex = entry.getKey();
                
                Map<Integer, Map<Integer, Float>> value = entry.getValue();
                
                Iterator<Entry<Integer, Map<Integer, Float>>> iter2 = value.entrySet().iterator();
                
                Map<Integer, Map<Integer, Float>> cMap = combined.get(outerClauseIndex);
                if (cMap == null) {
                    cMap = new HashMap<Integer, Map<Integer, Float>>();
                    combined.put(outerClauseIndex, cMap);
                }
                
                while (iter2.hasNext()) {
                    Entry<Integer, Map<Integer, Float>> entry2 = iter2.next();
                    
                    Integer innerClauseIndex = entry2.getKey();
                    
                    Map<Integer, Float> coeffIdxAndValue = entry2.getValue();
                    
                    Map<Integer, Float> cMap2 = cMap.get(innerClauseIndex);
                    
                    if (cMap2 == null) {
                        cMap2 = new HashMap<Integer, Float>();
                        cMap.put(innerClauseIndex, cMap2);
                    }
                    
                    cMap2.putAll(coeffIdxAndValue);
                }
            }
        }
        
        if (!redCustomCoeffUpperLimits.isEmpty()) {
            
            Iterator<Entry<Integer, Map<Integer, Map<Integer, Float>>>> iter 
                = redCustomCoeffUpperLimits.entrySet().iterator();
            
            while (iter.hasNext()) {
                
                Entry<Integer, Map<Integer, Map<Integer, Float>>> entry = iter.next();
                
                Integer outerClauseIndex = entry.getKey();
                Integer outerClauseIndexForList = Integer.valueOf(
                    outerClauseIndex.intValue() + nAllCoeff);
                
                Map<Integer, Map<Integer, Float>> value = entry.getValue();
                
                Map<Integer, Map<Integer, Float>> cMap = combined.get(outerClauseIndexForList);
                if (cMap == null) {
                    cMap = new HashMap<Integer, Map<Integer, Float>>();
                    combined.put(outerClauseIndexForList, cMap);
                }
                
                Iterator<Entry<Integer, Map<Integer, Float>>> iter2 = value.entrySet().iterator();

                while (iter2.hasNext()) {
                    
                    Entry<Integer, Map<Integer, Float>> entry2 = iter2.next();
                    
                    Integer innerClauseIndex = entry2.getKey();
                    
                    Map<Integer, Float> coeffIdxAndValue = entry2.getValue();
                    
                    Map<Integer, Float> cMap2 = cMap.get(innerClauseIndex);
                    
                    if (cMap2 == null) {
                        cMap2 = new HashMap<Integer, Float>();
                        cMap.put(innerClauseIndex, cMap2);
                    }

                    cMap2.putAll(coeffIdxAndValue);
                }
            }
        }
        
        return combined;
    }
    
    public float[][] getAllAndBlueCoeffLowerLimits() {
        
        float[][] a0 = allLowerLimits;
        float[][] a1 = blueLowerLimits;
        
        int n = a0.length + a1.length;
    
        int count = 0;
        float[][] a = new float[n][];
        for (int i = 0; i < a0.length; i++) {
            a[count] = Arrays.copyOf(a0[i], a0[i].length);
            count++;
        }
        for (int i = 0; i < a1.length; i++) {
            a[count] = Arrays.copyOf(a1[i], a1[i].length);
            count++;
        }
        
        return a;
    }
    
    public float[][] getAllAndBlueCoeffUpperLimits() {
        
        float[][] a0 = allUpperLimits;
        float[][] a1 = blueUpperLimits;
        
        int n = a0.length + a1.length;
    
        int count = 0;
        float[][] a = new float[n][];
        for (int i = 0; i < a0.length; i++) {
            a[count] = Arrays.copyOf(a0[i], a0[i].length);
            count++;
        }
        for (int i = 0; i < a1.length; i++) {
            a[count] = Arrays.copyOf(a1[i], a1[i].length);
            count++;
        }
        
        return a;
    }
    
    /**
     * custom coefficient limits are stored by outerClauseIndex 
       (same as used by allLowerLimits for example)
       then by innerClauseIndex, then by coefficentIndex as the x in PairFloat
      
     * @return 
     */
    public Map<Integer, Map<Integer, Map<Integer, Float>>> getAllAndBlueCustomCoeffLowerLimits() {
        
        int nAllCoeff = allLowerLimits.length;
        //int nBlueCoeff = blueLowerLimits.length;
        
        // combine the maps, but any blue clause indexes will need nAllCoeff added
        Map<Integer, Map<Integer, Map<Integer, Float>>> combined 
            = new HashMap<Integer, Map<Integer, Map<Integer, Float>>>();
        
        if (!allCustomCoeffLowerLimits.isEmpty()) {
            
            Iterator<Entry<Integer, Map<Integer, Map<Integer, Float>>>> iter 
                = allCustomCoeffLowerLimits.entrySet().iterator();
            
            while (iter.hasNext()) {
                
                Entry<Integer, Map<Integer, Map<Integer, Float>>> entry = iter.next();
                
                Integer outerClauseIndex = entry.getKey();
                
                Map<Integer, Map<Integer, Float>> value = entry.getValue();
                
                Iterator<Entry<Integer, Map<Integer, Float>>> iter2 = value.entrySet().iterator();
                
                Map<Integer, Map<Integer, Float>> cMap = combined.get(outerClauseIndex);
                if (cMap == null) {
                    cMap = new HashMap<Integer, Map<Integer, Float>>();
                    combined.put(outerClauseIndex, cMap);
                }
                
                while (iter2.hasNext()) {
                    Entry<Integer, Map<Integer, Float>> entry2 = iter2.next();
                    
                    Integer innerClauseIndex = entry2.getKey();
                    
                    Map<Integer, Float> coeffIdxAndValue = entry2.getValue();
                    
                    Map<Integer, Float> cMap2 = cMap.get(innerClauseIndex);
                    
                    if (cMap2 == null) {
                        cMap2 = new HashMap<Integer, Float>();
                        cMap.put(innerClauseIndex, cMap2);
                    }
                    
                    cMap2.putAll(coeffIdxAndValue);
                }
            }
        }
        
        if (!blueCustomCoeffLowerLimits.isEmpty()) {
            
            Iterator<Entry<Integer, Map<Integer, Map<Integer, Float>>>> iter 
                = blueCustomCoeffLowerLimits.entrySet().iterator();
            
            while (iter.hasNext()) {
                
                Entry<Integer, Map<Integer, Map<Integer, Float>>> entry = iter.next();
                
                Integer outerClauseIndex = entry.getKey();
                Integer outerClauseIndexForList = Integer.valueOf(
                    outerClauseIndex.intValue() + nAllCoeff);
                
                Map<Integer, Map<Integer, Float>> value = entry.getValue();
                
                Map<Integer, Map<Integer, Float>> cMap = combined.get(outerClauseIndexForList);
                if (cMap == null) {
                    cMap = new HashMap<Integer, Map<Integer, Float>>();
                    combined.put(outerClauseIndexForList, cMap);
                }
                
                Iterator<Entry<Integer, Map<Integer, Float>>> iter2 = value.entrySet().iterator();

                while (iter2.hasNext()) {
                    
                    Entry<Integer, Map<Integer, Float>> entry2 = iter2.next();
                    
                    Integer innerClauseIndex = entry2.getKey();
                    
                    Map<Integer, Float> coeffIdxAndValue = entry2.getValue();
                    
                    Map<Integer, Float> cMap2 = cMap.get(innerClauseIndex);
                    
                    if (cMap2 == null) {
                        cMap2 = new HashMap<Integer, Float>();
                        cMap.put(innerClauseIndex, cMap2);
                    }

                    cMap2.putAll(coeffIdxAndValue);
                }
            }
        }
        
        return combined;
    }
    
    /**
     * custom coefficient limits are stored by outerClauseIndex 
       (same as used by allLowerLimits for example)
       then by innerClauseIndex, then by coefficentIndex as the x in PairFloat
      
     * @return 
     */
    public Map<Integer, Map<Integer, Map<Integer, Float>>> getAllAndBlueCustomCoeffUpperLimits() {
        
        int nAllCoeff = allLowerLimits.length;
        //int nBlueCoeff = blueLowerLimits.length;
        
        // combine the maps, but any blue clause indexes will need nAllCoeff added
        Map<Integer, Map<Integer, Map<Integer, Float>>> combined 
            = new HashMap<Integer, Map<Integer, Map<Integer, Float>>>();
        
        if (!allCustomCoeffUpperLimits.isEmpty()) {
            
            Iterator<Entry<Integer, Map<Integer, Map<Integer, Float>>>> iter 
                = allCustomCoeffUpperLimits.entrySet().iterator();
            
            while (iter.hasNext()) {
                
                Entry<Integer, Map<Integer, Map<Integer, Float>>> entry = iter.next();
                
                Integer outerClauseIndex = entry.getKey();
                
                Map<Integer, Map<Integer, Float>> value = entry.getValue();
                
                Iterator<Entry<Integer, Map<Integer, Float>>> iter2 = value.entrySet().iterator();
                
                Map<Integer, Map<Integer, Float>> cMap = combined.get(outerClauseIndex);
                if (cMap == null) {
                    cMap = new HashMap<Integer, Map<Integer, Float>>();
                    combined.put(outerClauseIndex, cMap);
                }
                
                while (iter2.hasNext()) {
                    Entry<Integer, Map<Integer, Float>> entry2 = iter2.next();
                    
                    Integer innerClauseIndex = entry2.getKey();
                    
                    Map<Integer, Float> coeffIdxAndValue = entry2.getValue();
                    
                    Map<Integer, Float> cMap2 = cMap.get(innerClauseIndex);
                    
                    if (cMap2 == null) {
                        cMap2 = new HashMap<Integer, Float>();
                        cMap.put(innerClauseIndex, cMap2);
                    }
                    
                    cMap2.putAll(coeffIdxAndValue);
                }
            }
        }
        
        if (!blueCustomCoeffUpperLimits.isEmpty()) {
            
            Iterator<Entry<Integer, Map<Integer, Map<Integer, Float>>>> iter 
                = blueCustomCoeffUpperLimits.entrySet().iterator();
            
            while (iter.hasNext()) {
                
                Entry<Integer, Map<Integer, Map<Integer, Float>>> entry = iter.next();
                
                Integer outerClauseIndex = entry.getKey();
                Integer outerClauseIndexForList = Integer.valueOf(
                    outerClauseIndex.intValue() + nAllCoeff);
                
                Map<Integer, Map<Integer, Float>> value = entry.getValue();
                
                Map<Integer, Map<Integer, Float>> cMap = combined.get(outerClauseIndexForList);
                if (cMap == null) {
                    cMap = new HashMap<Integer, Map<Integer, Float>>();
                    combined.put(outerClauseIndexForList, cMap);
                }
                
                Iterator<Entry<Integer, Map<Integer, Float>>> iter2 = value.entrySet().iterator();

                while (iter2.hasNext()) {
                    
                    Entry<Integer, Map<Integer, Float>> entry2 = iter2.next();
                    
                    Integer innerClauseIndex = entry2.getKey();
                    
                    Map<Integer, Float> coeffIdxAndValue = entry2.getValue();
                    
                    Map<Integer, Float> cMap2 = cMap.get(innerClauseIndex);
                    
                    if (cMap2 == null) {
                        cMap2 = new HashMap<Integer, Float>();
                        cMap.put(innerClauseIndex, cMap2);
                    }

                    cMap2.putAll(coeffIdxAndValue);
                }
            }
        }
        
        return combined;
    } 
    
    /*
    =================
    ALL:
    ==============
    } else if (
        00         (skyStDevContrast > 0.)
        01         && ((Math.abs(contrastV)/skyStDevContrast) > 10.)
        } else if (
        02         (skyStDevContrast > 0.)
        03         && ((Math.abs(contrastV) > 0.1)
        04,05,06   && ((Math.abs(contrastV)/skyStDevContrast) > (1.5 + (Math.abs(contrastV)-0.5)*(-2.0))))
        07         && ((Math.abs(colorDiffV)/skyStDevColorDiff) > 2.5)
        } else if (
        08         (skyStDevContrast > 0.005)
        09         && ((Math.abs(contrastV)/skyStDevContrast) > 5.)
        10         && ((Math.abs(colorDiffV)/skyStDevColorDiff) > 5.)
                   &&
                   // if cieXY diffs are zero and within stdev, these are sky,
                   // so test for opposite for boundary pixel
        11         (((diffCIEX > 0.001)
        12             || ((diffCIEX/localSky.getStdDevCIEX()) > 1.5)
        13             || (diffCIEY > 0.001)
        14             || ((diffCIEY/localSky.getStdDevCIEY()) > 1.5)))
        15         && (skyStDevColorDiff > 1.)
    
    =================
    RED:
    =================
    00         (skyStDevContrast > 0.)
    01         && ((Math.abs(colorDiffV)/skyStDevColorDiff) > 15.*diffCIEX)
    02         && (diffCIEX > 0.03)
    03         && ((diffCIEX/localSky.getStdDevCIEX()) > 15.*diffCIEX)
        ) {
    } else if (
    04         (skyStDevContrast > 0.)
    05         && ((Math.abs(colorDiffV)/skyStDevColorDiff) > 15.*diffCIEY)
    06         && (diffCIEY > 0.03)
    07         && ((diffCIEY/localSky.getStdDevCIEY()) > 15.*diffCIEY)
        ) {
    } else if (skyStDevContrast == 0.) {
        if (contrastV >= 0.) {
            doNotAddToStack = true; <====== this is not handled.  if it's necessary, hard code it back in
        }
    
    =================
    BLUE:
    =================
    00         (skyStDevContrast > 0.0)
    01         && (contrastV < 0.05)
    02         && ((Math.abs(contrastV)/skyStDevContrast) > 1.5)
    03         && ((Math.abs(colorDiffV)/skyStDevColorDiff) > 3.0)
    04         && ((bPercentV < 0.37)
    05         && (bV > 199)
    06         && (gV > 199))
    } else if (
    07         (skyStDevContrast > 0.0)
    08         && (Math.abs(contrastV) > 0.05)
    09         && ((Math.abs(colorDiffV)/skyStDevColorDiff) > 1.5)
    10         && ((diffCIEX/localSky.getStdDevCIEX()) > 0.9) // ?
    11         && ((diffCIEY/localSky.getStdDevCIEY()) > 0.9) // ?
    12         && (skyStDevColorDiff > 0.)
    } else if (
    13         (contrastV < 0.05)
    14         && (diffCIEX < 0.005)
    15         && (diffCIEY < 0.005)
               &&
    16         ((Math.abs(0.33 - rPercentV) > 0.08)
    17         || (Math.abs(0.33 - gPercentV) > 0.03)
    18         || (Math.abs(0.33 - bPercentV) > 0.03)
    19         || (gV > 199)
    20         || (bV > 199))
    
    */
}
