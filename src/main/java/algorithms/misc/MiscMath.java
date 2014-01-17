package algorithms.misc;

import java.security.SecureRandom;

/**
 *
 * @author nichole
 */
public class MiscMath {

    public static int findPowerOf10_2(float a) {

        if (a == 0) {
            return 0;
        }

        int power = 0;

        if (a <= 1.0f) {
            while (a < 1.0) {
                a *=  10.0;
                power--;
            }
        } else {
            // precision errors in multiplication here are trouble for non base2 numbers such as powers of 10
            while (a >= 1.0) {
                a /= 10.0;
                power++;
            }
            power--;
        }

        return power;
    }

    public static int findPowerOf10(float a) {

        if (a == 0) {
            return 0;
        }
        if (a < 0.f) {
            a *= -1.0f;
        }
        double b = Math.log10(a);
        if (b > 0) {
            return (int)b;
        } else {
            if (b >= 1) {
                return (int)Math.round(b);
            } else if (b > -1) {
                // fractions between -1 and +1
                 return findPowerOf10_2(a);
            } else {
                return (int)Math.round(b);
            }
        }
    }

    /**
     * find max but ignore values such as FLOAT.MAX_VALUE, infinity, and NAN
     * @param a
     * @return
     */
    public static float findMax(float[] a) {
        float max = Float.MIN_VALUE;
        for (int i = 0; i < a.length; i++) {
            if ((a[i] > max) && !Float.isInfinite(a[i]) && !Float.isNaN(a[i]) && (a[i] < Float.MAX_VALUE)) {
                max = a[i];
            }
        }
        return max;
    }

    /**
     * find max but ignore values such as FLOAT.MAX_VALUE, infinity, and NAN
     * @param a
     * @return
     */
    public static int findMax(int[] a) {
        int max = Integer.MIN_VALUE;
        for (int i = 0; i < a.length; i++) {
            if ((a[i] > max) && !Float.isInfinite(a[i]) && !Float.isNaN(a[i]) && (a[i] < Float.MAX_VALUE)) {
                max = a[i];
            }
        }
        return max;
    }

    public static float findMin(float[] a) {
        float min = Float.MAX_VALUE;
        for (int i = 0; i < a.length; i++) {
            if ((a[i] < min) && !Float.isInfinite(a[i])) {
                min = a[i];
            }
        }
        return min;
    }

    /**
     * find max but ignore values such as FLOAT.MAX_VALUE, infinity, and NAN
     * @param a
     * @return
     */
    public static int findYMaxIndex(float[] a) {
        if (a == null) {
            return -1;
        }
        float max = Float.MIN_VALUE;
        int index = -1;
        for (int i = 0; i < a.length; i++) {
            if ((a[i] > max) && !Float.isInfinite(a[i]) && !Float.isNaN(a[i]) && (a[i] < Float.MAX_VALUE)) {
                max = a[i];
                index = i;
            }
        }
        return index;
    }

    public static int findYMaxIndex(int[] a) {
        if (a == null) {
            return -1;
        }
        int max = Integer.MIN_VALUE;
        int index = -1;
        for (int i = 0; i < a.length; i++) {
            if ((a[i] > max) && (a[i] < Integer.MAX_VALUE)) {
                max = a[i];
                index = i;
            }
        }
        return index;
    }

    public static float[] calculateOuterRoundedMinAndMax(float[] a) {

        // find the powers of 10 for the data min and max
        float xmin = Float.MAX_VALUE;
        float xmax = Float.MIN_VALUE;
        for (int i = 0; i < a.length; i++) {
            if (a[i] > xmax) {
                xmax = a[i];
            }
            if (a[i] < xmin) {
                xmin = a[i];
            }
        }

        xmax = MiscMath.roundUpByLargestPower(xmax);

        xmin = MiscMath.roundDownByLargestPower(xmin);

        // xmax > 1 and xmin is between 0 and 1, round xmin down
        if ((xmax > 1) && (xmin > 0) && (xmin < 1.0)) {
            xmin = 0;
        }

        return new float[]{xmin, xmax};
    }

    /**
     * Round number down to in the largest power.
     * For example,
     *     roundDownByLargestPower(3.1) returns 3.0
     *     roundDownByLargestPower(-3.1) returns -4.0
     *     roundDownByLargestPower(31.1) returns 30.0
     *
     * @param f
     * @return
     */
    public static float roundDownByLargestPower(float f) {

        if (f == 0) {
            return 0;
        }

        int pow = findPowerOf10(f);
        double pow10 = Math.pow(10, pow);

        double r;
        if (f > 0) {
            double m = f % pow10;
            r = f - m;
        } else {
            int m = (int)(f/pow10);
            r = pow10 * (m - 1);
        }

        return (float)r;
    }

    /**
     Round number up to the next digit in the largest power.

        MiscMath.roundUpByLargestPower(31.1f) == 40.0f;

        MiscMath.roundUpByLargestPower(0.11f) == 1.0f;

        MiscMath.roundUpByLargestPower(-0.011f) == -0.02f;

        MiscMath.roundUpByLargestPower(310.1f) == 400.0f;

        MiscMath.roundUpByLargestPower(-3.1) == -4.0f;

        MiscMath.roundUpByLargestPower(3.1) == 4.0f;

    @param f
    @return
    */
    public static float roundUpByLargestPower(float f) {

        int pow = findPowerOf10(f);

        double pow10 = Math.pow(10, pow);

        int d = (int)(f/pow10);
        int m;
        if (f == pow10) {
            m = d;
        } else if (f > 0) {
            // residual ?
            if (f > (d*pow10)) {
                m = d + 1;
            } else {
                m = d;
            }
        } else {
            // decimals
            m = d - 1;
        }

        float r = (float) (m * pow10);

        return r;
    }
    
    /*
     http://en.wikipedia.org/wiki/Poisson_distribution
         algorithm poisson random number (Knuth):
     
     @param sr an instance of secure random which is a strong random number generator
     @param lambda a value less than 4 roughly looks like the tested dataset histograms
     */
    public static int poissonRandom(SecureRandom sr, int lambda) {

        double L = Math.exp(-lambda);

        int k = 0;

        double p = 1;

        do {
            k = k + 1;

            double u = sr.nextDouble();

            p = p * u;

        } while (p > L);

        return k - 1;
    }

    //public static double naturalLogOfNegative(double negativeNumber) {
        /*
         * Using 2 rules to get the natural log of a negative number:
         * 
         * (1) ln(a*b) = ln(a) + ln(b)
         * 
         * (2) and the Taylor series expansion:
         *     ln(1+x)= x - (x^2)/2 + (x^3)/3 - ...
         *     
         *     For x = -1.005, ln(1.0 + -1.001) is approx -3.9120... curve is nearly vertical here, so not a good way to proceed
         *
         * Let z = a*b
         * ln(a*b) = ln(a) + ln(1 - 1.001)
         * ln(z) = ln(z/(-1*(1 - 1.001))) + ln(1 - 1.001)  
         *    
         *      ln(z) = ln(z/0.001) + ln(-0.001)  (with x = -1.001)
         *            = ln(z/0.001) - 3.9120
         */
        /*if (negativeNumber == 0) {
            return 0;
        } if (negativeNumber > 0) {
            return Math.log(negativeNumber);
        }
        
        double onePlusX = 1 + -0.0001;
                
        double a = negativeNumber/onePlusX;
        
        double lnz = Math.log(a) - 3.9120;.
        
        return lnz;
    }*/
    
    public static double taylor(double number) {
        return taylor(number, 100);
    }
    /**
     * use the taylor series to approximate the natural log of a number.
     * Note that internally, the variable x which is number - 1.
     * 
     * ln(1 + x) is computed, so number = 1 + x and therefore x = number - 1;
     * 
     * ln(1+x) = x - (x^2)/2 + (x^3)/3 - ...
     * 
     * @param number
     * @param n
     * @return
     */
    protected static double taylor(double number, int n) {
        double x = number - 1.0;
        
        double sum = 0;
        for (int i = 1; i <= n; i++) {
            double f = Math.pow(x, i)/(float)i;
            if (i % 2 == 0) {
                sum -= f;
                System.out.print("     MINUS i=" + i);
            } else {
                sum += f;
                System.out.print("     PLUS i=" + i);
            }            
        }
        return sum;
    }
}
