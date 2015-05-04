package algorithms.imageProcessing.optimization;

/**
 *
 * @author nichole
 */
public class SkylineANDedClauses {
    
    private ANDedClauses[] generalClauses = null;
    
    private ANDedClauses[] redOnlyClauses = null;
    
    private ANDedClauses[] blueOnlyClauses = null;
   
    private ANDedClauses[] generalClausesLowerLimits = null;
    private ANDedClauses[] generalClausesUpperLimits = null;
    
    private ANDedClauses[] redOnlyClausesLowerLimits = null;
    private ANDedClauses[] redOnlyClausesUpperLimits = null;
    
    private ANDedClauses[] blueOnlyClausesLowerLimits = null;
    private ANDedClauses[] blueOnlyClausesUpperLimits = null;
    
    public SkylineANDedClauses() {
        
        initGeneralClauses();
        
        initRedOnlyClauses();
        
        initBlueOnlyClauses();
       
    }
    
    public ANDedClauses[] getAllClauses() {
        
        return copyAndAppend(generalClauses, blueOnlyClauses, redOnlyClauses);
    }
    
    public ANDedClauses[] getAllClausesLowerLimits() {
        
        return copyAndAppend(generalClausesLowerLimits, 
            blueOnlyClausesLowerLimits, redOnlyClausesLowerLimits);
    }
    
    public ANDedClauses[] getAllClausesUpperLimits() {
        
        return copyAndAppend(generalClausesUpperLimits, 
            blueOnlyClausesUpperLimits, redOnlyClausesUpperLimits);
    }
    
    public ANDedClauses[] getGeneralAndRedClauses() {
        
        return copyAndAppend(generalClauses, redOnlyClauses);
    }
    
    public ANDedClauses[] getGeneralAndRedClausesLowerLimits() {
        
        return copyAndAppend(generalClausesLowerLimits, redOnlyClausesLowerLimits);
    }
    
    public ANDedClauses[] getGeneralAndRedClausesUpperLimits() {
        
        return copyAndAppend(generalClausesUpperLimits, redOnlyClausesUpperLimits);
    }
    
    public ANDedClauses[] getGeneralAndBlueClauses() {
        
        return copyAndAppend(generalClauses, blueOnlyClauses);
    }
    
    public ANDedClauses[] getGeneralAndBlueClausesLowerLimits() {
        
        return copyAndAppend(generalClausesLowerLimits, blueOnlyClausesLowerLimits);
    }
    
    public ANDedClauses[] getGeneralAndBlueClausesUpperLimits() {
        
        return copyAndAppend(generalClausesUpperLimits, blueOnlyClausesUpperLimits);
    }
    
    private ANDedClauses[] copyAndAppend(ANDedClauses[] a0, ANDedClauses[] a1) {
        
        int n = a0.length + a1.length;
    
        ANDedClauses[] a = new ANDedClauses[n];
        for (int i = 0; i < a0.length; i++) {
            a[i] = a0[i].copy();
        }
        
        for (int i = 0; i < a1.length; i++) {
            a[i + a0.length] = a1[i].copy();
        }
        
        return a;
    }
    
    private ANDedClauses[] copyAndAppend(ANDedClauses[] a0, ANDedClauses[] a1,
        ANDedClauses[] a2) {
        
        int n = a0.length + a1.length + a2.length;
    
        ANDedClauses[] a = new ANDedClauses[n];
        
        for (int i = 0; i < a0.length; i++) {
            a[i] = a0[i].copy();
        }
        
        for (int i = 0; i < a1.length; i++) {
            a[i + a0.length] = a1[i].copy();
        }
        
        for (int i = 0; i < a2.length; i++) {
            a[i + a0.length + a1.length] = a2[i].copy();
        }
        
        return a;
    }
    
    private void initGeneralClauses() {
        
        float customCoeffPlaceHolder = Float.MAX_VALUE;
        
        /*
        c0:
            00         (skyStDevContrast > 0.)
            01         && ((Math.abs(contrastV)/skyStDevContrast) > 10.)
        */
        ANDedClauses c0 = new ANDedClauses(2, SKYCONDITIONAL.ALL);
        c0.set(0, PARAM.STDEV_CONTRAST, PARAM.INT_ONE, 
            COMPARISON.GREATER_THAN, 0.0f);
        c0.set(1, PARAM.ABSOLUTE_CONTRAST, PARAM.STDEV_CONTRAST, 
            COMPARISON.GREATER_THAN, 10.0f);
        
        ANDedClauses c0Lower = c0.copy();
        ANDedClauses c0Upper = c0.copy();
        c0Lower.coefficients[0] = 0.00001f;
        c0Upper.coefficients[0] = 0.2f;
        c0Lower.coefficients[1] = 2.5f;
        c0Upper.coefficients[1] = 30.0f;
        
        /*
        c1      
            02         (skyStDevContrast > 0.)
            03         && ((Math.abs(contrastV) > 0.1)
            04,05,06   && ((Math.abs(contrastV)/skyStDevContrast) > (1.5 + (Math.abs(contrastV)-0.5)*(-2.0))))  <-------
            07         && ((Math.abs(colorDiffV)/skyStDevColorDiff) > 2.5)
        */
        ANDedClauses c1 = new ANDedClauses(4, SKYCONDITIONAL.ALL);
        c1.set(0, PARAM.STDEV_CONTRAST, PARAM.INT_ONE, 
            COMPARISON.GREATER_THAN, 0.0f);
        c1.set(1, PARAM.ABSOLUTE_CONTRAST, PARAM.INT_ONE, 
            COMPARISON.GREATER_THAN, 0.1f);
        
        // this is replaced by the custom coefficient
        c1.set(2, PARAM.ABSOLUTE_CONTRAST, PARAM.STDEV_CONTRAST, 
            COMPARISON.GREATER_THAN, customCoeffPlaceHolder);
        c1.setACustomCoefficient(2, new CustomCoeff00()); 
        c1.setCustomCoefficientVariable(4, Float.valueOf(1.5f));
        c1.setCustomCoefficientVariable(5, Float.valueOf(0.5f));
        c1.setCustomCoefficientVariable(6, Float.valueOf(-2.0f));
        
        c1.set(3, PARAM.ABSOLUTE_DIFF_BLUE_OR_RED, PARAM.STDEV_BLUE_OR_RED, 
            COMPARISON.GREATER_THAN, 2.5f);
    
        ANDedClauses c1Lower = c1.copy();
        ANDedClauses c1Upper = c1.copy();
        c1Lower.coefficients[0] = 0.00001f;
        c1Upper.coefficients[0] = 0.2f;
        c1Lower.coefficients[1] = 0.01f;
        c1Upper.coefficients[1] = 0.5f;
        c1Lower.customCoefficientVariables.put(Integer.valueOf(4), 
            Float.valueOf(0.2f));
        c1Upper.customCoefficientVariables.put(Integer.valueOf(4), 
            Float.valueOf(10.0f));
        c1Lower.customCoefficientVariables.put(Integer.valueOf(5), 
            Float.valueOf(0.05f));
        c1Upper.customCoefficientVariables.put(Integer.valueOf(5), 
            Float.valueOf(0.95f));
        c1Lower.customCoefficientVariables.put(Integer.valueOf(6), 
            Float.valueOf(-10.0f));
        c1Upper.customCoefficientVariables.put(Integer.valueOf(6), 
            Float.valueOf(-0.5f));
        
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
            COMPARISON.GREATER_THAN, 0.005f);
        c2.set(1, PARAM.ABSOLUTE_CONTRAST, PARAM.STDEV_CONTRAST, 
            COMPARISON.GREATER_THAN, 5.0f);
        c2.set(2, PARAM.ABSOLUTE_DIFF_BLUE_OR_RED, PARAM.STDEV_BLUE_OR_RED, 
            COMPARISON.GREATER_THAN, 5.0f);
        c2.set(3, PARAM.DIFF_CIEX, PARAM.INT_ONE, 
            COMPARISON.GREATER_THAN, 0.001f);
        c2.set(4, PARAM.STDEV_BLUE_OR_RED, PARAM.INT_ONE, 
            COMPARISON.GREATER_THAN, 1.0f);
       
        ANDedClauses c2Lower = c2.copy();
        ANDedClauses c2Upper = c2.copy();
        c2Lower.coefficients[0] = 0.00001f;
        c2Upper.coefficients[0] = 0.5f;
        c2Lower.coefficients[1] = 1.1f;
        c2Upper.coefficients[1] = 10.0f;
        c2Lower.coefficients[2] = 0.00001f;
        c2Upper.coefficients[2] = 0.05f;
        c2Lower.coefficients[3] = 0.2f;
        c2Upper.coefficients[3] = 100.0f;
        
        ANDedClauses c3 = new ANDedClauses(5, SKYCONDITIONAL.ALL);
        c3.set(0, PARAM.STDEV_CONTRAST, PARAM.INT_ONE, 
            COMPARISON.GREATER_THAN, 0.005f);
        c3.set(1, PARAM.ABSOLUTE_CONTRAST, PARAM.STDEV_CONTRAST, 
            COMPARISON.GREATER_THAN, 5.0f);
        c3.set(2, PARAM.ABSOLUTE_DIFF_BLUE_OR_RED, PARAM.STDEV_BLUE_OR_RED, 
            COMPARISON.GREATER_THAN, 5.0f);
        c3.set(3, PARAM.DIFF_CIEX, PARAM.STDEV_CIEX, 
            COMPARISON.GREATER_THAN, 1.5f);
        c3.set(4, PARAM.STDEV_BLUE_OR_RED, PARAM.INT_ONE, 
            COMPARISON.GREATER_THAN, 1.0f);
        
        ANDedClauses c3Lower = c3.copy();
        ANDedClauses c3Upper = c3.copy();
        c3Lower.coefficients[0] = 0.00001f;
        c3Upper.coefficients[0] = 0.5f;
        c3Lower.coefficients[1] = 1.1f;
        c3Upper.coefficients[1] = 10.0f;
        c3Lower.coefficients[2] = 1.1f;
        c3Upper.coefficients[2] = 10.0f;
        c3Lower.coefficients[3] = 1.1f;
        c3Upper.coefficients[3] = 5.0f;
        c3Lower.coefficients[4] = 0.2f;
        c3Upper.coefficients[4] = 100.0f;
        
        ANDedClauses c4 = new ANDedClauses(5, SKYCONDITIONAL.ALL);
        c4.set(0, PARAM.STDEV_CONTRAST, PARAM.INT_ONE, 
            COMPARISON.GREATER_THAN, 0.005f);
        c4.set(1, PARAM.ABSOLUTE_CONTRAST, PARAM.STDEV_CONTRAST, 
            COMPARISON.GREATER_THAN, 5.0f);
        c4.set(2, PARAM.ABSOLUTE_DIFF_BLUE_OR_RED, PARAM.STDEV_BLUE_OR_RED, 
            COMPARISON.GREATER_THAN, 5.0f);
        c4.set(3, PARAM.DIFF_CIEY, PARAM.INT_ONE, 
            COMPARISON.GREATER_THAN, 0.001f);
        c4.set(4, PARAM.STDEV_BLUE_OR_RED, PARAM.INT_ONE, 
            COMPARISON.GREATER_THAN, 1.0f);
       
        ANDedClauses c4Lower = c4.copy();
        ANDedClauses c4Upper = c4.copy();
        c4Lower.coefficients[0] = 0.00001f;
        c4Upper.coefficients[0] = 0.5f;
        c4Lower.coefficients[1] = 1.1f;
        c4Upper.coefficients[1] = 10.0f;
        c4Lower.coefficients[2] = 1.1f;
        c4Upper.coefficients[2] = 10.0f;
        c4Lower.coefficients[3] = 0.00001f;
        c4Upper.coefficients[3] = 0.05f;
        c4Lower.coefficients[4] = 0.2f;
        c4Upper.coefficients[4] = 100.0f;
        
        ANDedClauses c5 = new ANDedClauses(5, SKYCONDITIONAL.ALL);
        c5.set(0, PARAM.STDEV_CONTRAST, PARAM.INT_ONE, 
            COMPARISON.GREATER_THAN, 0.005f);
        c5.set(1, PARAM.ABSOLUTE_CONTRAST, PARAM.STDEV_CONTRAST, 
            COMPARISON.GREATER_THAN, 5.0f);
        c5.set(2, PARAM.ABSOLUTE_DIFF_BLUE_OR_RED, PARAM.STDEV_BLUE_OR_RED, 
            COMPARISON.GREATER_THAN, 5.0f);
        c5.set(3, PARAM.DIFF_CIEY, PARAM.STDEV_CIEY, 
            COMPARISON.GREATER_THAN, 1.5f);
        c5.set(4, PARAM.STDEV_BLUE_OR_RED, PARAM.INT_ONE, 
            COMPARISON.GREATER_THAN, 1.0f);
        
        ANDedClauses c5Lower = c5.copy();
        ANDedClauses c5Upper = c5.copy();
        c5Lower.coefficients[0] = 0.00001f;
        c5Upper.coefficients[0] = 0.5f;
        c5Lower.coefficients[1] = 1.1f;
        c5Upper.coefficients[1] = 10.0f;
        c5Lower.coefficients[2] = 1.1f;
        c5Upper.coefficients[2] = 10.0f;
        c5Lower.coefficients[3] = 1.1f;
        c5Upper.coefficients[3] = 5.0f;
        c5Lower.coefficients[4] = 0.2f;
        c5Upper.coefficients[4] = 100.0f;
        
        generalClausesLowerLimits = new ANDedClauses[]{
            c0Lower, c1Lower, c2Lower, c3Lower, c4Lower, c5Lower};
        
        generalClausesUpperLimits = new ANDedClauses[]{
            c0Upper, c1Upper, c2Upper, c3Upper, c4Upper, c5Upper};
    
        generalClauses = new ANDedClauses[]{c0, c1, c2, c3, c4, c5};
    }
    
    private void initRedOnlyClauses() {
    
        float customCoeffPlaceHolder = Float.MAX_VALUE;
        
        /*
        00         (skyStDevContrast > 0.)
        01         && ((Math.abs(colorDiffV)/skyStDevColorDiff) > 15.*diffCIEX) <---
        02         && (diffCIEX > 0.03)
        03         && ((diffCIEX/localSky.getStdDevCIEX()) > 15.*diffCIEX) <---
        */
        ANDedClauses c0 = new ANDedClauses(4, SKYCONDITIONAL.RED);
        c0.set(0, PARAM.STDEV_CONTRAST, PARAM.INT_ONE, 
            COMPARISON.GREATER_THAN, 0.0f);
        
        c0.set(1, PARAM.ABSOLUTE_DIFF_BLUE_OR_RED, 
            PARAM.STDEV_BLUE_OR_RED, COMPARISON.GREATER_THAN, customCoeffPlaceHolder);
        c0.setACustomCoefficient(1, new CustomCoeff01());
        c0.setCustomCoefficientVariable(1, Float.valueOf(15.0f));
        
        c0.set(2, PARAM.DIFF_CIEX, 
            PARAM.INT_ONE, COMPARISON.GREATER_THAN, 0.03f);
        
        c0.set(3, PARAM.DIFF_CIEX, 
            PARAM.STDEV_CIEX, COMPARISON.GREATER_THAN, customCoeffPlaceHolder);
        c0.setACustomCoefficient(3, new CustomCoeff02());
        c0.setCustomCoefficientVariable(3, Float.valueOf(15.0f));
        
        ANDedClauses c0Lower = c0.copy();
        ANDedClauses c0Upper = c0.copy();
        c0Lower.coefficients[0] = 0.00001f;
        c0Upper.coefficients[0] = 0.2f;
        c0Lower.coefficients[1] = 0.01f;
        c0Upper.coefficients[1] = 0.5f;
        c0Lower.coefficients[2] = 0.00001f;
        c0Upper.coefficients[2] = 0.2f;
        c0Lower.customCoefficientVariables.put(Integer.valueOf(1), 
            Float.valueOf(4.0f));
        c0Upper.customCoefficientVariables.put(Integer.valueOf(1), 
            Float.valueOf(30.0f));
        c0Lower.customCoefficientVariables.put(Integer.valueOf(3), 
            Float.valueOf(4.0f));
        c0Upper.customCoefficientVariables.put(Integer.valueOf(3), 
            Float.valueOf(30.0f));
        
        /*
        04         (skyStDevContrast > 0.)
        05         && ((Math.abs(colorDiffV)/skyStDevColorDiff) > 15.*diffCIEY) <---
        06         && (diffCIEY > 0.03)
        07         && ((diffCIEY/localSky.getStdDevCIEY()) > 15.*diffCIEY) <---
        */
        ANDedClauses c1 = new ANDedClauses(4, SKYCONDITIONAL.RED);
        c1.set(0, PARAM.STDEV_CONTRAST, PARAM.INT_ONE, 
            COMPARISON.GREATER_THAN, 0.0f);
        
        c1.set(1, PARAM.ABSOLUTE_DIFF_BLUE_OR_RED, 
            PARAM.STDEV_BLUE_OR_RED, COMPARISON.GREATER_THAN, customCoeffPlaceHolder);
        c1.setACustomCoefficient(1, new CustomCoeff03());
        c1.setCustomCoefficientVariable(5, Float.valueOf(15.f));
        
        c1.set(2, PARAM.DIFF_CIEY, 
            PARAM.INT_ONE, COMPARISON.GREATER_THAN, 0.03f);
        
        c1.set(3, PARAM.DIFF_CIEY, 
            PARAM.STDEV_CIEY, COMPARISON.GREATER_THAN, customCoeffPlaceHolder);
        c1.setACustomCoefficient(3, new CustomCoeff04());
        c1.setCustomCoefficientVariable(7, Float.valueOf(15.0f));
        
        ANDedClauses c1Lower = c1.copy();
        ANDedClauses c1Upper = c1.copy();
        c1Lower.coefficients[0] = 0.00001f;
        c1Upper.coefficients[0] = 0.2f;
        c1Lower.coefficients[1] = 0.01f;
        c1Upper.coefficients[1] = 30.0f;
        c1Lower.coefficients[2] = 0.00001f;
        c1Upper.coefficients[2] = 0.2f;
        c1Lower.customCoefficientVariables.put(Integer.valueOf(5), 
            Float.valueOf(4.0f));
        c1Upper.customCoefficientVariables.put(Integer.valueOf(5), 
            Float.valueOf(30.0f));
        c1Lower.customCoefficientVariables.put(Integer.valueOf(7), 
            Float.valueOf(4.0f));
        c1Upper.customCoefficientVariables.put(Integer.valueOf(7), 
            Float.valueOf(30.0f));
        
        redOnlyClausesLowerLimits = new ANDedClauses[]{c0Lower, c1Lower};
        
        redOnlyClausesUpperLimits = new ANDedClauses[]{c0Upper, c1Upper};
    
        redOnlyClauses = new ANDedClauses[]{c0, c1};
   
    }
    
    private void initBlueOnlyClauses() {
        
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
            COMPARISON.GREATER_THAN, 0.0f);
        c0.set(1, PARAM.CONTRAST, 
            PARAM.INT_ONE, COMPARISON.LESS_THAN, 0.05f);
        c0.set(2, PARAM.ABSOLUTE_CONTRAST, 
            PARAM.STDEV_CONTRAST, COMPARISON.GREATER_THAN, 1.5f);
        c0.set(3, PARAM.ABSOLUTE_DIFF_BLUE_OR_RED, 
            PARAM.STDEV_BLUE_OR_RED, COMPARISON.GREATER_THAN, 3.0f);
        c0.set(4, PARAM.B_DIV_TOT, 
            PARAM.INT_ONE, COMPARISON.LESS_THAN, 0.37f);
        
        c0.set(5, PARAM.BLUE, 
            PARAM.INT_ONE, COMPARISON.GREATER_THAN, 199);
        c0.set(6, PARAM.GREEN, 
            PARAM.INT_ONE, COMPARISON.GREATER_THAN, 199);
    
        ANDedClauses c0Lower = c0.copy();
        ANDedClauses c0Upper = c0.copy();
        c0Lower.coefficients[0] = 0.00001f;
        c0Upper.coefficients[0] = 0.2f;
        c0Lower.coefficients[1] = 0.0000f;
        c0Upper.coefficients[1] = 0.5f;
        c0Lower.coefficients[2] = 1.5f;
        c0Upper.coefficients[2] = 30.0f;
        c0Lower.coefficients[3] = 1.5f;
        c0Upper.coefficients[3] = 30.0f;
        c0Lower.coefficients[4] = 0.3f;
        c0Upper.coefficients[4] = 0.4f;
        c0Lower.coefficients[5] = 100.f;
        c0Upper.coefficients[5] = 250.f;
        c0Lower.coefficients[6] = 100.f;
        c0Upper.coefficients[6] = 250.f;
        
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
            COMPARISON.GREATER_THAN, 0.0f);
        c1.set(1, PARAM.ABSOLUTE_CONTRAST, 
            PARAM.INT_ONE, COMPARISON.GREATER_THAN, 0.05f);
        c1.set(2, PARAM.ABSOLUTE_DIFF_BLUE_OR_RED, 
            PARAM.STDEV_BLUE_OR_RED, COMPARISON.GREATER_THAN, 1.5f);
        c1.set(3, PARAM.DIFF_CIEX, 
            PARAM.STDEV_CIEX, COMPARISON.GREATER_THAN, 0.6f);
        
        c1.set(4, PARAM.DIFF_CIEY, 
            PARAM.STDEV_CIEY, COMPARISON.GREATER_THAN, 0.9f);
        
        c1.set(5, PARAM.STDEV_BLUE_OR_RED, 
            PARAM.INT_ONE, COMPARISON.GREATER_THAN, 0.0f);
               
        ANDedClauses c1Lower = c1.copy();
        ANDedClauses c1Upper = c1.copy();
        c1Lower.coefficients[0] = 0.00001f;
        c1Upper.coefficients[0] = 0.2f;
        c1Lower.coefficients[1] = 0.0000f;
        c1Upper.coefficients[1] = 0.5f;
        c1Lower.coefficients[2] = 1.5f;
        c1Upper.coefficients[2] = 30.0f;
        c1Lower.coefficients[3] = 0.75f;
        c1Upper.coefficients[3] = 2.0f;
        c1Lower.coefficients[4] = 0.75f;
        c1Upper.coefficients[4] = 2.0f;
        c1Lower.coefficients[5] = 0.00001f;
        c1Upper.coefficients[5] = 3.0f;
        
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
            COMPARISON.LESS_THAN, 0.05f);
        
        c2.set(1, PARAM.DIFF_CIEX, PARAM.INT_ONE, 
            COMPARISON.LESS_THAN, 0.005f);
        
        c2.set(2, PARAM.DIFF_CIEY, PARAM.INT_ONE, 
            COMPARISON.LESS_THAN, 0.005f);
        
        c2.set(3, PARAM.DIFF_R_DIV_TOT_ONE_THIRD, PARAM.INT_ONE, 
            COMPARISON.GREATER_THAN, 0.08f);
       
        ANDedClauses c2Lower = c2.copy();
        ANDedClauses c2Upper = c2.copy();
        c2Lower.coefficients[0] = 0.00001f;
        c2Upper.coefficients[0] = 0.2f;
        c2Lower.coefficients[1] = 0.0000f;
        c2Upper.coefficients[1] = 0.5f;
        c2Lower.coefficients[2] = 0.0000f;
        c2Upper.coefficients[2] = 0.5f;
        c2Lower.coefficients[3] = 0.01f;
        c2Upper.coefficients[3] = 0.3f;
        
        ANDedClauses c3 = new ANDedClauses(4, SKYCONDITIONAL.BLUE);
        c3.set(0, PARAM.CONTRAST, PARAM.INT_ONE, 
            COMPARISON.LESS_THAN, 0.05f);
        
        c3.set(1, PARAM.DIFF_CIEX, PARAM.INT_ONE, 
            COMPARISON.LESS_THAN, 0.005f);
        
        c3.set(2, PARAM.DIFF_CIEY, PARAM.INT_ONE, 
            COMPARISON.LESS_THAN, 0.005f);
        
        c3.set(3, PARAM.DIFF_G_DIV_TOT_ONE_THIRD, PARAM.INT_ONE, 
            COMPARISON.GREATER_THAN, 0.03f);
      
        ANDedClauses c3Lower = c3.copy();
        ANDedClauses c3Upper = c3.copy();
        c3Lower.coefficients[0] = 0.00001f;
        c3Upper.coefficients[0] = 0.2f;
        c3Lower.coefficients[1] = 0.0000f;
        c3Upper.coefficients[1] = 0.5f;
        c3Lower.coefficients[2] = 0.0000f;
        c3Upper.coefficients[2] = 0.5f;
        c3Lower.coefficients[3] = 0.01f;
        c3Upper.coefficients[3] = 0.3f;
        
        ANDedClauses c4 = new ANDedClauses(4, SKYCONDITIONAL.BLUE);
        c4.set(0, PARAM.CONTRAST, PARAM.INT_ONE, 
            COMPARISON.LESS_THAN, 0.05f);
        
        c4.set(1, PARAM.DIFF_CIEX, PARAM.INT_ONE, 
            COMPARISON.LESS_THAN, 0.005f);
        
        c4.set(2, PARAM.DIFF_CIEY, PARAM.INT_ONE, 
            COMPARISON.LESS_THAN, 0.005f);
        
        c4.set(3, PARAM.DIFF_B_DIV_TOT_ONE_THIRD, PARAM.INT_ONE, 
            COMPARISON.GREATER_THAN, 0.03f);
        
        ANDedClauses c4Lower = c4.copy();
        ANDedClauses c4Upper = c4.copy();
        c4Lower.coefficients[0] = 0.00001f;
        c4Upper.coefficients[0] = 0.2f;
        c4Lower.coefficients[1] = 0.0000f;
        c4Upper.coefficients[1] = 0.5f;
        c4Lower.coefficients[2] = 0.0000f;
        c4Upper.coefficients[2] = 0.5f;
        c4Lower.coefficients[3] = 0.01f;
        c4Upper.coefficients[3] = 0.3f;
        
        ANDedClauses c5 = new ANDedClauses(4, SKYCONDITIONAL.BLUE);
        c5.set(0, PARAM.CONTRAST, PARAM.INT_ONE, 
            COMPARISON.LESS_THAN, 0.05f);
        
        c5.set(1, PARAM.DIFF_CIEX, PARAM.INT_ONE, 
            COMPARISON.LESS_THAN, 0.005f);
        
        c5.set(2, PARAM.DIFF_CIEY, PARAM.INT_ONE, 
            COMPARISON.LESS_THAN, 0.005f);
        
        c5.set(3, PARAM.GREEN, PARAM.INT_ONE, 
            COMPARISON.GREATER_THAN, 199);
        
        ANDedClauses c5Lower = c5.copy();
        ANDedClauses c5Upper = c5.copy();
        c5Lower.coefficients[0] = 0.00001f;
        c5Upper.coefficients[0] = 0.2f;
        c5Lower.coefficients[1] = 0.0000f;
        c5Upper.coefficients[1] = 0.5f;
        c5Lower.coefficients[2] = 0.0000f;
        c5Upper.coefficients[2] = 0.5f;
        c5Lower.coefficients[3] = 100;
        c5Upper.coefficients[3] = 250;
        
        ANDedClauses c6 = new ANDedClauses(4, SKYCONDITIONAL.BLUE);
        c6.set(0, PARAM.CONTRAST, PARAM.INT_ONE, 
            COMPARISON.LESS_THAN, 0.05f);
        
        c6.set(1, PARAM.DIFF_CIEX, PARAM.INT_ONE, 
            COMPARISON.LESS_THAN, 0.005f);
        
        c6.set(2, PARAM.DIFF_CIEY, PARAM.INT_ONE, 
            COMPARISON.LESS_THAN, 0.005f);
        
        c6.set(3, PARAM.BLUE, PARAM.INT_ONE, 
            COMPARISON.GREATER_THAN, 199);
        
        ANDedClauses c6Lower = c6.copy();
        ANDedClauses c6Upper = c6.copy();
        c6Lower.coefficients[0] = 0.00001f;
        c6Upper.coefficients[0] = 0.2f;
        c6Lower.coefficients[1] = 0.0000f;
        c6Upper.coefficients[1] = 0.5f;
        c6Lower.coefficients[2] = 0.0000f;
        c6Upper.coefficients[2] = 0.5f;
        c6Lower.coefficients[3] = 100;
        c6Upper.coefficients[3] = 250;
        
        blueOnlyClausesLowerLimits = new ANDedClauses[]{c0Lower, c1Lower, 
            c2Lower, c3Lower, c4Lower, c5Lower, c6Lower};
        
        blueOnlyClausesUpperLimits = new ANDedClauses[]{c0Upper, c1Upper, 
            c2Upper, c3Upper, c4Upper, c5Upper, c6Upper};
        
        blueOnlyClauses = new ANDedClauses[]{c0, c1, c2, c3, c4, c5, c6};
   
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
