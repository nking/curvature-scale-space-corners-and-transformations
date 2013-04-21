package algorithms.misc;

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

    public static int findYMinIndex(float[] a) {
        float min = Float.MAX_VALUE;
        int index = -1;
        for (int i = 0; i < a.length; i++) {
            if ((a[i] < min) && !Float.isInfinite(a[i])) {
                min = a[i];
                index = i;
            }
        }
        return index;
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

    /**
     * find max but ignore values such as FLOAT.MAX_VALUE, infinity, and NAN
     * @param a
     * @return
     */
    public static int findYMaxIndex(float[] a, int lookAfterIndex) {
        if (a == null) {
            return -1;
        }
        float max = Float.MIN_VALUE;
        int index = -1;
        for (int i = (lookAfterIndex + 1); i < a.length; i++) {
            if ((a[i] > max) && !Float.isInfinite(a[i]) && !Float.isNaN(a[i]) && (a[i] < Float.MAX_VALUE)) {
                max = a[i];
                index = i;
            }
        }
        return index;
    }

    public static float[] calculateOuterRoundedMinAndMax(float[] a) {

        // find the powers of 10 for the data min and max
        float xmin = a[0];
        float xmax = a[0];
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

    public static float[] findMeanAndStDev(float[] values) {

        float mean = 0;

        for (int i = 0; i < values.length; i++) {
            mean += values[i];
        }
        mean /= values.length;

        float stdev = 0;
        for (int i = 0; i < values.length; i++) {
            float a = (values[i] - mean);
            stdev += ( a*a );
        }
        stdev = (float) Math.sqrt(stdev/(values.length - 1));

        return new float[]{mean, stdev};
    }

    public static float[] findMeanAndStDev(int[][] values) {

        float mean = 0;

        int count = 0;

        for (int i = 0; i < values.length; i++) {
            for (int j = 0; j < values.length; j++) {
                mean += values[i][j];
                count++;
            }
        }
        mean /= (float)count;

        float stdev = 0;
        for (int i = 0; i < values.length; i++) {

            for (int j = 0; j < values.length; j++) {

                float a = (values[i][j] - mean);
                stdev += ( a*a );
            }
        }
        stdev = (float) Math.sqrt(stdev/(float)(count - 1));

        return new float[]{mean, stdev};
    }

}
