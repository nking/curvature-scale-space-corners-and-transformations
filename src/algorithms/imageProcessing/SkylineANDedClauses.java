package algorithms.imageProcessing;

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
        0.0f, 0.05f, 1.5f, 0.9f, 0.9f, 0.0f,
        //13   14      15      16     17     18     19   20
        0.05f, 0.005f, 0.005f, 0.08f, 0.03f, 0.03f, 199, 199
    };
     
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
    
    public ANDedClauses[] getForAllSkies() {
    
        /*
        c0:
            00         (skyStDevContrast > 0.)
            01         && ((Math.abs(contrastV)/skyStDevContrast) > 10.)
        */
        ANDedClauses c0 = new ANDedClauses(2);
        c0.set(0, SKYCONDITIONAL.ALL, PARAM.STDEV_CONTRAST, PARAM.INT_ONE, 
            COMPARISON.GREATER_THAN, allSkiesCoeff[0]);
        c0.set(1, SKYCONDITIONAL.ALL, PARAM.ABSOLUTE_CONTRAST, PARAM.STDEV_CONTRAST, 
            COMPARISON.GREATER_THAN, allSkiesCoeff[1]);
        
        /*
        c1      
            02         (skyStDevContrast > 0.)
            03         && ((Math.abs(contrastV) > 0.1)
            04,05,06   && ((Math.abs(contrastV)/skyStDevContrast) > (1.5 + (Math.abs(contrastV)-0.5)*(-2.0))))  <-------
            07         && ((Math.abs(colorDiffV)/skyStDevColorDiff) > 2.5)
        */
        ANDedClauses c1 = new ANDedClauses(4);
        c1.set(0, SKYCONDITIONAL.ALL, PARAM.STDEV_CONTRAST, PARAM.INT_ONE, 
            COMPARISON.GREATER_THAN, allSkiesCoeff[2]);
        c1.set(1, SKYCONDITIONAL.ALL, PARAM.ABSOLUTE_CONTRAST, PARAM.STDEV_CONTRAST, 
            COMPARISON.GREATER_THAN, allSkiesCoeff[3]);
        
        c1.set(2, SKYCONDITIONAL.ALL, PARAM.ABSOLUTE_CONTRAST, PARAM.STDEV_CONTRAST, 
            COMPARISON.GREATER_THAN, allSkiesCoeff[4]);
        c1.setACustomCoefficient(2, new CustomCoeff00()); 
        c1.setCustomCoefficientVariable(4, Float.valueOf(allSkiesCoeff[4]));
        c1.setCustomCoefficientVariable(5, Float.valueOf(allSkiesCoeff[5]));
        c1.setCustomCoefficientVariable(6, Float.valueOf(allSkiesCoeff[6]));
        
        c1.set(3, SKYCONDITIONAL.ALL, PARAM.ABSOLUTE_DIFF_BLUE_OR_RED, PARAM.STDEV_BLUE_OR_RED, 
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
        ANDedClauses c2 = new ANDedClauses(5);
        c2.set(0, SKYCONDITIONAL.ALL, PARAM.STDEV_CONTRAST, PARAM.INT_ONE, 
            COMPARISON.GREATER_THAN, allSkiesCoeff[8]);
        c2.set(1, SKYCONDITIONAL.ALL, PARAM.ABSOLUTE_CONTRAST, PARAM.STDEV_CONTRAST, 
            COMPARISON.GREATER_THAN, allSkiesCoeff[9]);
        c2.set(2, SKYCONDITIONAL.ALL, PARAM.ABSOLUTE_DIFF_BLUE_OR_RED, PARAM.STDEV_BLUE_OR_RED, 
            COMPARISON.GREATER_THAN, allSkiesCoeff[10]);
    c2.set(3, SKYCONDITIONAL.ALL, PARAM.DIFF_CIEX, PARAM.INT_ONE, 
            COMPARISON.GREATER_THAN, allSkiesCoeff[11]);
        c2.set(4, SKYCONDITIONAL.ALL, PARAM.STDEV_BLUE_OR_RED, PARAM.INT_ONE, 
            COMPARISON.GREATER_THAN, allSkiesCoeff[15]);
        
        ANDedClauses c3 = new ANDedClauses(5);
        c3.set(0, SKYCONDITIONAL.ALL, PARAM.STDEV_CONTRAST, PARAM.INT_ONE, 
            COMPARISON.GREATER_THAN, allSkiesCoeff[8]);
        c3.set(1, SKYCONDITIONAL.ALL, PARAM.ABSOLUTE_CONTRAST, PARAM.STDEV_CONTRAST, 
            COMPARISON.GREATER_THAN, allSkiesCoeff[9]);
        c3.set(2, SKYCONDITIONAL.ALL, PARAM.ABSOLUTE_DIFF_BLUE_OR_RED, PARAM.STDEV_BLUE_OR_RED, 
            COMPARISON.GREATER_THAN, allSkiesCoeff[10]);
    c3.set(3, SKYCONDITIONAL.ALL, PARAM.DIFF_CIEX, PARAM.STDEV_CIEX, 
            COMPARISON.GREATER_THAN, allSkiesCoeff[12]);
        c3.set(4, SKYCONDITIONAL.ALL, PARAM.STDEV_BLUE_OR_RED, PARAM.INT_ONE, 
            COMPARISON.GREATER_THAN, allSkiesCoeff[15]);
        
        ANDedClauses c4 = new ANDedClauses(5);
        c4.set(0, SKYCONDITIONAL.ALL, PARAM.STDEV_CONTRAST, PARAM.INT_ONE, 
            COMPARISON.GREATER_THAN, allSkiesCoeff[8]);
        c4.set(1, SKYCONDITIONAL.ALL, PARAM.ABSOLUTE_CONTRAST, PARAM.STDEV_CONTRAST, 
            COMPARISON.GREATER_THAN, allSkiesCoeff[9]);
        c4.set(2, SKYCONDITIONAL.ALL, PARAM.ABSOLUTE_DIFF_BLUE_OR_RED, PARAM.STDEV_BLUE_OR_RED, 
            COMPARISON.GREATER_THAN, allSkiesCoeff[10]);
    c4.set(3, SKYCONDITIONAL.ALL, PARAM.DIFF_CIEY, PARAM.INT_ONE, 
            COMPARISON.GREATER_THAN, allSkiesCoeff[13]);
        c4.set(4, SKYCONDITIONAL.ALL, PARAM.STDEV_BLUE_OR_RED, PARAM.INT_ONE, 
            COMPARISON.GREATER_THAN, allSkiesCoeff[15]);
        
        ANDedClauses c5 = new ANDedClauses(5);
        c5.set(0, SKYCONDITIONAL.ALL, PARAM.STDEV_CONTRAST, PARAM.INT_ONE, 
            COMPARISON.GREATER_THAN, allSkiesCoeff[8]);
        c5.set(1, SKYCONDITIONAL.ALL, PARAM.ABSOLUTE_CONTRAST, PARAM.STDEV_CONTRAST, 
            COMPARISON.GREATER_THAN, allSkiesCoeff[9]);
        c5.set(2, SKYCONDITIONAL.ALL, PARAM.ABSOLUTE_DIFF_BLUE_OR_RED, PARAM.STDEV_BLUE_OR_RED, 
            COMPARISON.GREATER_THAN, allSkiesCoeff[10]);
    c5.set(3, SKYCONDITIONAL.ALL, PARAM.DIFF_CIEY, PARAM.STDEV_CIEY, 
            COMPARISON.GREATER_THAN, allSkiesCoeff[14]);
        c5.set(4, SKYCONDITIONAL.ALL, PARAM.STDEV_BLUE_OR_RED, PARAM.INT_ONE, 
            COMPARISON.GREATER_THAN, allSkiesCoeff[15]);
        
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
        ANDedClauses c0 = new ANDedClauses(4);
        c0.set(0, SKYCONDITIONAL.RED, PARAM.STDEV_CONTRAST, PARAM.INT_ONE, 
            COMPARISON.GREATER_THAN, redSkiesCoeff[0]);
        c0.set(1, SKYCONDITIONAL.RED, PARAM.ABSOLUTE_DIFF_BLUE_OR_RED, 
            PARAM.STDEV_BLUE_OR_RED, COMPARISON.GREATER_THAN, redSkiesCoeff[1]);
        
        c0.setACustomCoefficient(1, new CustomCoeff01());
        c0.setCustomCoefficientVariable(1, Float.valueOf(redSkiesCoeff[1]));
        
        c0.set(2, SKYCONDITIONAL.RED, PARAM.DIFF_CIEX, 
            PARAM.INT_ONE, COMPARISON.GREATER_THAN, redSkiesCoeff[2]);
        c0.set(3, SKYCONDITIONAL.RED, PARAM.DIFF_CIEX, 
            PARAM.STDEV_CIEX, COMPARISON.GREATER_THAN, redSkiesCoeff[3]);
        c0.setACustomCoefficient(3, new CustomCoeff02());
        c0.setCustomCoefficientVariable(3, Float.valueOf(redSkiesCoeff[3]));
       
        /*
        04         (skyStDevContrast > 0.)
        05         && ((Math.abs(colorDiffV)/skyStDevColorDiff) > 15.*diffCIEY) <---
        06         && (diffCIEY > 0.03)
        07         && ((diffCIEY/localSky.getStdDevCIEY()) > 15.*diffCIEY) <---
        */
        ANDedClauses c1 = new ANDedClauses(4);
        c1.set(0, SKYCONDITIONAL.RED, PARAM.STDEV_CONTRAST, PARAM.INT_ONE, 
            COMPARISON.GREATER_THAN, redSkiesCoeff[4]);
        c1.set(1, SKYCONDITIONAL.RED, PARAM.ABSOLUTE_DIFF_BLUE_OR_RED, 
            PARAM.STDEV_BLUE_OR_RED, COMPARISON.GREATER_THAN, redSkiesCoeff[5]);
        
        c1.setACustomCoefficient(1, new CustomCoeff03());
        c1.setCustomCoefficientVariable(5, Float.valueOf(redSkiesCoeff[5]));
        
        c1.set(2, SKYCONDITIONAL.RED, PARAM.DIFF_CIEY, 
            PARAM.INT_ONE, COMPARISON.GREATER_THAN, redSkiesCoeff[6]);
        c1.set(3, SKYCONDITIONAL.RED, PARAM.DIFF_CIEY, 
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
        ANDedClauses c0 = new ANDedClauses(7);
        c0.set(0, SKYCONDITIONAL.BLUE, PARAM.STDEV_CONTRAST, PARAM.INT_ONE, 
            COMPARISON.GREATER_THAN, blueSkiesCoeff[0]);
        c0.set(1, SKYCONDITIONAL.BLUE, PARAM.CONTRAST, 
            PARAM.INT_ONE, COMPARISON.LESS_THAN, blueSkiesCoeff[1]);
        
        c0.set(2, SKYCONDITIONAL.BLUE, PARAM.ABSOLUTE_CONTRAST, 
            PARAM.STDEV_CONTRAST, COMPARISON.GREATER_THAN, blueSkiesCoeff[2]);
        
        c0.set(3, SKYCONDITIONAL.BLUE, PARAM.ABSOLUTE_DIFF_BLUE_OR_RED, 
            PARAM.STDEV_BLUE_OR_RED, COMPARISON.GREATER_THAN, blueSkiesCoeff[3]);
        
        c0.set(4, SKYCONDITIONAL.BLUE, PARAM.B_DIV_TOT, 
            PARAM.INT_ONE, COMPARISON.LESS_THAN, blueSkiesCoeff[4]);
        
        c0.set(5, SKYCONDITIONAL.BLUE, PARAM.BLUE, 
            PARAM.INT_ONE, COMPARISON.GREATER_THAN, blueSkiesCoeff[5]);
        
        c0.set(6, SKYCONDITIONAL.BLUE, PARAM.GREEN, 
            PARAM.INT_ONE, COMPARISON.GREATER_THAN, blueSkiesCoeff[6]);
       
        /*
        07         (skyStDevContrast > 0.0)
        08         && (Math.abs(contrastV) > 0.05)
        09         && ((Math.abs(colorDiffV)/skyStDevColorDiff) > 1.5)
        10         && ((diffCIEX/localSky.getStdDevCIEX()) > 0.9) // ?
        11         && ((diffCIEY/localSky.getStdDevCIEY()) > 0.9) // ?
        12         && (skyStDevColorDiff > 0.)
        */
        ANDedClauses c1 = new ANDedClauses(6);
        c1.set(0, SKYCONDITIONAL.BLUE, PARAM.STDEV_CONTRAST, PARAM.INT_ONE, 
            COMPARISON.GREATER_THAN, blueSkiesCoeff[7]);
        
        c1.set(1, SKYCONDITIONAL.BLUE, PARAM.ABSOLUTE_CONTRAST, 
            PARAM.INT_ONE, COMPARISON.GREATER_THAN, blueSkiesCoeff[8]);
        
        c1.set(2, SKYCONDITIONAL.BLUE, PARAM.ABSOLUTE_DIFF_BLUE_OR_RED, 
            PARAM.STDEV_BLUE_OR_RED, COMPARISON.GREATER_THAN, blueSkiesCoeff[9]);
        
        c1.set(3, SKYCONDITIONAL.BLUE, PARAM.DIFF_CIEX, 
            PARAM.STDEV_CIEX, COMPARISON.GREATER_THAN, blueSkiesCoeff[10]);
        
        c1.set(4, SKYCONDITIONAL.BLUE, PARAM.DIFF_CIEY, 
            PARAM.STDEV_CIEY, COMPARISON.GREATER_THAN, blueSkiesCoeff[11]);
        
        c1.set(5, SKYCONDITIONAL.BLUE, PARAM.STDEV_BLUE_OR_RED, 
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
        ANDedClauses c2 = new ANDedClauses(4);
        c2.set(0, SKYCONDITIONAL.BLUE, PARAM.CONTRAST, PARAM.INT_ONE, 
            COMPARISON.LESS_THAN, blueSkiesCoeff[13]);
        
        c2.set(1, SKYCONDITIONAL.BLUE, PARAM.DIFF_CIEX, PARAM.INT_ONE, 
            COMPARISON.LESS_THAN, blueSkiesCoeff[14]);
        
        c2.set(2, SKYCONDITIONAL.BLUE, PARAM.DIFF_CIEY, PARAM.INT_ONE, 
            COMPARISON.LESS_THAN, blueSkiesCoeff[15]);
        
        c2.set(3, SKYCONDITIONAL.BLUE, PARAM.DIFF_R_DIV_TOT_ONE_THIRD, PARAM.INT_ONE, 
            COMPARISON.GREATER_THAN, blueSkiesCoeff[16]);
       
        
        ANDedClauses c3 = new ANDedClauses(4);
        c3.set(0, SKYCONDITIONAL.BLUE, PARAM.CONTRAST, PARAM.INT_ONE, 
            COMPARISON.LESS_THAN, blueSkiesCoeff[13]);
        
        c3.set(1, SKYCONDITIONAL.BLUE, PARAM.DIFF_CIEX, PARAM.INT_ONE, 
            COMPARISON.LESS_THAN, blueSkiesCoeff[14]);
        
        c3.set(2, SKYCONDITIONAL.BLUE, PARAM.DIFF_CIEY, PARAM.INT_ONE, 
            COMPARISON.LESS_THAN, blueSkiesCoeff[15]);
        
        c3.set(3, SKYCONDITIONAL.BLUE, PARAM.DIFF_G_DIV_TOT_ONE_THIRD, PARAM.INT_ONE, 
            COMPARISON.GREATER_THAN, blueSkiesCoeff[17]);
        
        
        ANDedClauses c4 = new ANDedClauses(4);
        c4.set(0, SKYCONDITIONAL.BLUE, PARAM.CONTRAST, PARAM.INT_ONE, 
            COMPARISON.LESS_THAN, blueSkiesCoeff[13]);
        
        c4.set(1, SKYCONDITIONAL.BLUE, PARAM.DIFF_CIEX, PARAM.INT_ONE, 
            COMPARISON.LESS_THAN, blueSkiesCoeff[14]);
        
        c4.set(2, SKYCONDITIONAL.BLUE, PARAM.DIFF_CIEY, PARAM.INT_ONE, 
            COMPARISON.LESS_THAN, blueSkiesCoeff[15]);
        
        c4.set(3, SKYCONDITIONAL.BLUE, PARAM.DIFF_B_DIV_TOT_ONE_THIRD, PARAM.INT_ONE, 
            COMPARISON.GREATER_THAN, blueSkiesCoeff[18]);
        
   
        ANDedClauses c5 = new ANDedClauses(4);
        c5.set(0, SKYCONDITIONAL.BLUE, PARAM.CONTRAST, PARAM.INT_ONE, 
            COMPARISON.LESS_THAN, blueSkiesCoeff[13]);
        
        c5.set(1, SKYCONDITIONAL.BLUE, PARAM.DIFF_CIEX, PARAM.INT_ONE, 
            COMPARISON.LESS_THAN, blueSkiesCoeff[14]);
        
        c5.set(2, SKYCONDITIONAL.BLUE, PARAM.DIFF_CIEY, PARAM.INT_ONE, 
            COMPARISON.LESS_THAN, blueSkiesCoeff[15]);
        
        c5.set(3, SKYCONDITIONAL.BLUE, PARAM.GREEN, PARAM.INT_ONE, 
            COMPARISON.GREATER_THAN, blueSkiesCoeff[19]);
        
        
        ANDedClauses c6 = new ANDedClauses(4);
        c6.set(0, SKYCONDITIONAL.BLUE, PARAM.CONTRAST, PARAM.INT_ONE, 
            COMPARISON.LESS_THAN, blueSkiesCoeff[13]);
        
        c6.set(1, SKYCONDITIONAL.BLUE, PARAM.DIFF_CIEX, PARAM.INT_ONE, 
            COMPARISON.LESS_THAN, blueSkiesCoeff[14]);
        
        c6.set(2, SKYCONDITIONAL.BLUE, PARAM.DIFF_CIEY, PARAM.INT_ONE, 
            COMPARISON.LESS_THAN, blueSkiesCoeff[15]);
        
        c6.set(3, SKYCONDITIONAL.BLUE, PARAM.BLUE, PARAM.INT_ONE, 
            COMPARISON.GREATER_THAN, blueSkiesCoeff[20]);
        
        return new ANDedClauses[]{c0, c1, c2, c3, c4, c5, c6};
    }
    
    public float[] getFittableCoefficientsForBlueSkies() {
        return blueSkiesCoeff;
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
