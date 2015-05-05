package algorithms.imageProcessing.optimization;

/**
 *
 * @author nichole
 */
public class SkylineANDedClauses2 {
    
    private ANDedClauses[] generalClauses = null;
    
    private ANDedClauses[] redOnlyClauses = null;
    
    private ANDedClauses[] blueOnlyClauses = null;
   
    private ANDedClauses[] generalClausesLowerLimits = null;
    private ANDedClauses[] generalClausesUpperLimits = null;
    
    private ANDedClauses[] redOnlyClausesLowerLimits = null;
    private ANDedClauses[] redOnlyClausesUpperLimits = null;
    
    private ANDedClauses[] blueOnlyClausesLowerLimits = null;
    private ANDedClauses[] blueOnlyClausesUpperLimits = null;
    
    public SkylineANDedClauses2() {
        
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
        00         (skyStDevContrast > 0.)
        01         && ((Math.abs(contrastV)/skyStDevContrast) > 10.)   RANGE 0.1 to 10.
        */        
        ANDedClauses c0 = new ANDedClauses(2, SKYCONDITIONAL.ALL);
        c0.set(0, PARAM.STDEV_CONTRAST, PARAM.INT_ONE, 
            COMPARISON.GREATER_THAN, 0.04f);
        c0.set(1, PARAM.ABSOLUTE_CONTRAST, PARAM.STDEV_CONTRAST, 
            COMPARISON.GREATER_THAN, 7.03f);
        
        ANDedClauses c0Lower = c0.copy();
        ANDedClauses c0Upper = c0.copy();
        c0Lower.coefficients[0] = 0.00001f;
        c0Upper.coefficients[0] = 0.2f;
        c0Lower.coefficients[1] = 0.1f;
        c0Upper.coefficients[1] = 10.0f;
        
        /*
        } else if (// 1.5 + (0.006911656819283962 - 0.72999996)*(-1.8))
        02         (skyStDevContrast > 0.)
        03         && ((Math.abs(contrastV) > 0.1)                     RANGE 0.01 to 1.0
        04,05,06   && ((Math.abs(contrastV)/skyStDevContrast) > (1.5 + (Math.abs(contrastV)-0.5)*(-2.0))))
                                                                 0.5                        0.1   -3
                                                                 3.5                        1.0   +1
        */
        
        ANDedClauses c1 = new ANDedClauses(3, SKYCONDITIONAL.ALL);
        c1.set(0, PARAM.STDEV_CONTRAST, PARAM.INT_ONE, 
            COMPARISON.GREATER_THAN, 0.12f);
        c1.set(1, PARAM.ABSOLUTE_CONTRAST, PARAM.INT_ONE, 
            COMPARISON.GREATER_THAN, 0.406f);
        
        // this is replaced by the custom coefficient
        c1.set(2, PARAM.ABSOLUTE_CONTRAST, PARAM.STDEV_CONTRAST, 
            COMPARISON.GREATER_THAN, customCoeffPlaceHolder);
        c1.setACustomCoefficient(2, new CustomCoeff00()); 
        c1.setCustomCoefficientVariable(4, Float.valueOf(1.55f));
        c1.setCustomCoefficientVariable(5, Float.valueOf(0.55f));
        c1.setCustomCoefficientVariable(6, Float.valueOf(-1.8f));
        
        ANDedClauses c1Lower = c1.copy();
        ANDedClauses c1Upper = c1.copy();
        c1Lower.coefficients[0] = 0.00001f;
        c1Upper.coefficients[0] = 0.2f;
        c1Lower.coefficients[1] = 0.01f;
        c1Upper.coefficients[1] = 1.0f;
        c1Lower.customCoefficientVariables.put(Integer.valueOf(4), 
            Float.valueOf(0.5f));
        c1Upper.customCoefficientVariables.put(Integer.valueOf(4), 
            Float.valueOf(4.0f));
        c1Lower.customCoefficientVariables.put(Integer.valueOf(5), 
            Float.valueOf(0.1f));
        c1Upper.customCoefficientVariables.put(Integer.valueOf(5), 
            Float.valueOf(1.0f));
        c1Lower.customCoefficientVariables.put(Integer.valueOf(6), 
            Float.valueOf(-3.0f));
        c1Upper.customCoefficientVariables.put(Integer.valueOf(6), 
            Float.valueOf(1.0f));
        
        //} else if (
        //07         (skyStDevColorDiff > 0.0)
        //08         && ((Math.abs(colorDiffV)/skyStDevColorDiff) > 5.)  RANGE 0.5 to 5.
        ANDedClauses c2 = new ANDedClauses(2, SKYCONDITIONAL.ALL);
        c2.set(0, PARAM.STDEV_BLUE_OR_RED, PARAM.INT_ONE, 
            COMPARISON.GREATER_THAN, 0.0f);
        c2.set(1, PARAM.ABSOLUTE_DIFF_BLUE_OR_RED, PARAM.STDEV_BLUE_OR_RED, 
            COMPARISON.GREATER_THAN, 3.65f);
       
        ANDedClauses c2Lower = c2.copy();
        ANDedClauses c2Upper = c2.copy();
        c2Lower.coefficients[0] = 0.00001f;
        c2Upper.coefficients[0] = 0.5f;
        c2Lower.coefficients[1] = 0.5f;
        c2Upper.coefficients[1] = 5.0f;
        
        /*
        //} else if (
        //08         (diffCIEX > 0.005)
        //09         ((diffCIEX/localSky.getStdDevCIEX()) > 1.5)         RANGE 0.5 to 5.
        
        ANDedClauses c3 = new ANDedClauses(2, SKYCONDITIONAL.ALL);
        c3.set(0, PARAM.DIFF_CIEX, PARAM.INT_ONE, 
            COMPARISON.GREATER_THAN, 0.4f);
        c3.set(1, PARAM.DIFF_CIEX, PARAM.STDEV_CIEX, 
            COMPARISON.GREATER_THAN, 2.75f);
        
        ANDedClauses c3Lower = c3.copy();
        ANDedClauses c3Upper = c3.copy();
        c3Lower.coefficients[0] = 0.00001f;
        c3Upper.coefficients[0] = 0.5f;
        c3Lower.coefficients[1] = 0.5f;
        c3Upper.coefficients[1] = 5.0f;
        
        //} else if (
        //10         (diffCIEX > 0.005)
        //11         ((diffCIEX/localSky.getStdDevCIEX()) > 1.5)         RANGE 0.5 to 5.
        ANDedClauses c4 = new ANDedClauses(2, SKYCONDITIONAL.ALL);
        c4.set(0, PARAM.DIFF_CIEY, PARAM.INT_ONE, 
            COMPARISON.GREATER_THAN, 1.E-5f);
        c4.set(1, PARAM.DIFF_CIEY, PARAM.STDEV_CIEY, 
            COMPARISON.GREATER_THAN, 0.5f);
       
        ANDedClauses c4Lower = c4.copy();
        ANDedClauses c4Upper = c4.copy();
        c4Lower.coefficients[0] = 0.00001f;
        c4Upper.coefficients[0] = 0.5f;
        c4Lower.coefficients[1] = 0.5f;
        c4Upper.coefficients[1] = 5.0f;
        */
        
        generalClausesLowerLimits = new ANDedClauses[]{
            c0Lower, c1Lower, c2Lower/*, c3Lower, c4Lower*/
        };
        
        generalClausesUpperLimits = new ANDedClauses[]{
            c0Upper, c1Upper, c2Upper/*, c3Upper, c4Upper*/
        };
    
        generalClauses = new ANDedClauses[]{c0, c1, c2/*, c3, c4*/
        };
    }
    
    private void initRedOnlyClauses() {        
        
        redOnlyClausesLowerLimits = new ANDedClauses[]{};
        
        redOnlyClausesUpperLimits = new ANDedClauses[]{};
    
        redOnlyClauses = new ANDedClauses[]{};
   
    }
    
    private void initBlueOnlyClauses() {
        
        blueOnlyClausesLowerLimits = new ANDedClauses[]{};
        
        blueOnlyClausesUpperLimits = new ANDedClauses[]{};
        
        blueOnlyClauses = new ANDedClauses[]{};
   
    }
    
    /*
    =================
    ALL:
    ==============
    } else if (
        00         (skyStDevContrast > 0.)
        01         && ((Math.abs(contrastV)/skyStDevContrast) > 10.)   RANGE 0.1 to 10.
        } else if (
        02         (skyStDevContrast > 0.)
        03         && ((Math.abs(contrastV) > 0.1)                     RANGE 0.01 to 1.0
        04,05,06   && ((Math.abs(contrastV)/skyStDevContrast) > (1.5 + (Math.abs(contrastV)-0.5)*(-2.0))))
                                                                 0.5                        0.1   -3
                                                                 3.5                        1.0   +1
        } else if (
        07         (skyStDevColorDiff > 0.0)
        08         && ((Math.abs(colorDiffV)/skyStDevColorDiff) > 5.)  RANGE 0.5 to 5.
        } else if (
        08         (diffCIEX > 0.005)
        09         ((diffCIEX/localSky.getStdDevCIEX()) > 1.5)         RANGE 0.5 to 5.
        } else if (
        10         (diffCIEY > 0.005)
        11         ((diffCIEY/localSky.getStdDevCIEY()) > 1.5)         RANGE 0.5 to 5.
    
    If choose starter points as low, mid, high, there would need to be
    3^12 starter points for a proper range search = 531,441 starter points.
    
    If only had 6 coefficients, a range search would be complete with 729 starter points.
    
    The range search using starter points is then local search thereafter. The top
    10 fits' coefficients are retained and the others discarded.
    */
}
