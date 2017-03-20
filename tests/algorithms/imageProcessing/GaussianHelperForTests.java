package algorithms.imageProcessing;

import algorithms.util.VeryLargeNumber;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author nichole
 */
public class GaussianHelperForTests {
    
    private static Logger log = Logger.getLogger(GaussianHelperForTests.class.getName());
    
    public static float expectedFWHM(float sigma) {
        return (float)(sigma * 2.f * Math.sqrt(2*Math.log(2)));
    }
    
    public static float measureFWHM(int[] x, int[] y) {
        
        int yMaxIdx = -1;
        float yMax = Float.NEGATIVE_INFINITY;
        for (int i = 0; i < y.length; i++) {
            if (y[i] > yMax) {
                yMax = y[i];
                yMaxIdx = i;
            }
        }

        int yHalfMax0Idx = -1;
        int yHalfMax1Idx = -1;
        for (int i = 0; i < yMaxIdx; i++) {
            if (y[i] < (yMax / 2.)) {
                yHalfMax0Idx = i;
            }
        }
        for (int i = (y.length - 1); i > yMaxIdx; i--) {
            if (y[i] < (yMax / 2.)) {
                yHalfMax1Idx = i;
            }
        }
        
        if ((yHalfMax1Idx == -1) || (yHalfMax0Idx == -1)) {
            return -1;
        }
        
        return x[yHalfMax1Idx] - x[yHalfMax0Idx];
    }
    
    public static float measureFWHM(int[] x, float[] y) {
        
        int yMaxIdx = -1;
        float yMax = Float.NEGATIVE_INFINITY;
        for (int i = 0; i < y.length; i++) {
            if (y[i] > yMax) {
                yMax = y[i];
                yMaxIdx = i;
            }
        }

        int yHalfMax0Idx = -1;
        int yHalfMax1Idx = -1;
        for (int i = 0; i < yMaxIdx; i++) {
            if (y[i] < (yMax / 2.)) {
                yHalfMax0Idx = i;
            }
        }
        for (int i = (y.length - 1); i > yMaxIdx; i--) {
            if (y[i] < (yMax / 2.)) {
                yHalfMax1Idx = i;
            }
        }
        
        if ((yHalfMax1Idx == -1) || (yHalfMax0Idx == -1)) {
            return -1;
        }
        
        return x[yHalfMax1Idx] - x[yHalfMax0Idx];
    }
    
    public static float measureFWHM(float[] x, float[] y) {
        
        if (x.length == 1) {
            return 1;
        }
        
        int yMaxIdx = -1;
        float yMax = Float.NEGATIVE_INFINITY;
        for (int i = 0; i < y.length; i++) {
            if (y[i] > yMax) {
                yMax = y[i];
                yMaxIdx = i;
            }
        }

        float halfMax = yMax/2.f;
        
        int yHalfMax0Idx = -1;
        int yHalfMax1Idx = -1;
        for (int i = 0; i < yMaxIdx; i++) {
            if (y[i] < halfMax) {
                yHalfMax0Idx = i;
            } else {
                break;
            }
        }
        for (int i = (y.length - 1); i > yMaxIdx; i--) {
            if (y[i] < halfMax) {
                yHalfMax1Idx = i;
            } else {
                break;
            }
        }
        
        if ((yHalfMax1Idx == -1) || (yHalfMax0Idx == -1)) {
            return -1;
        }
        
        return x[yHalfMax1Idx] - x[yHalfMax0Idx];
    }
    
    public static float measureFWHMRoundToSmaller(float[] x, float[] y) {
        
        int yMaxIdx = -1;
        float yMax = Float.NEGATIVE_INFINITY;
        for (int i = 0; i < y.length; i++) {
            if (y[i] > yMax) {
                yMax = y[i];
                yMaxIdx = i;
            }
        }

        float halfMax = yMax/2.f;
        
        int yHalfMax0Idx = -1;
        int yHalfMax1Idx = -1;
        for (int i = 0; i < yMaxIdx; i++) {
            if (y[i] > halfMax) {
                yHalfMax0Idx = i;
                break;
            }
        }
        for (int i = (y.length - 1); i > yMaxIdx; i--) {
            if (y[i] > halfMax) {
                yHalfMax1Idx = i;
                break;
            }
        }
        
        if ((yHalfMax1Idx == -1) || (yHalfMax0Idx == -1)) {
            return -1;
        }
        
        return x[yHalfMax1Idx] - x[yHalfMax0Idx];
    }
    
    public static float measureFWHMRoundToSmaller(int[] x, int[] y) {
        
        int yMaxIdx = -1;
        float yMax = Float.NEGATIVE_INFINITY;
        for (int i = 0; i < y.length; i++) {
            if (y[i] > yMax) {
                yMax = y[i];
                yMaxIdx = i;
            }
        }

        float halfMax = yMax/2.f;
        
        int yHalfMax0Idx = -1;
        int yHalfMax1Idx = -1;
        for (int i = 0; i < yMaxIdx; i++) {
            if (y[i] > halfMax) {
                yHalfMax0Idx = i;
                break;
            }
        }
        for (int i = (y.length - 1); i > yMaxIdx; i--) {
            if (y[i] > halfMax) {
                yHalfMax1Idx = i;
                break;
            }
        }
        if (yHalfMax1Idx == -1) {
            yHalfMax1Idx = yMaxIdx + 1;
        }
        if (yHalfMax0Idx == -1) {
            yHalfMax0Idx = yMaxIdx + 1;
        }
        
        return x[yHalfMax1Idx] - x[yHalfMax0Idx];
    }
    
    public static double areaUnderTheCurve(float[] x, float[] y) {
        
        // using trapezoidal rule
        
        if (x.length == 0) {
            return 0;
        } else if (x.length == 1) {
            // this isn't correct, but is fine for the one case it's used for
            return (y[0]);
        }
        
        double sum = 0;
        
        for (int i = 0; i < (x.length - 1); i++) {
            float yTerm = y[i + 1] + y[i];
            float xLen = x[i + 1] - x[i];
            if (xLen < 0) {
                xLen *= -1;
            }
            sum += (yTerm * xLen);
        }
        
        /*
         float xLen = x[i + 1] - x[i];
            
            /*
                 @          @
             @                  @
             |   |          |   |
            -------         -------
            
             |   |          |   |
             @                  @       
                 @          @          
            
            sum += (y[i] * xLen);
            sum += ((y[i + 1] - y[i]) * xLen * 0.5f);
        */
        
        sum *= 0.5;
        
        return sum;
    }
    
    public static double areaUnderTheCurve(int[] x, int[] y) {
        
        // using trapezoidal rule
        
        if (x.length == 0) {
            return 0;
        } else if (x.length == 1) {
            // this isn't correct, but is fine for the one case it's used for
            return (y[0]);
        }
        
        double sum = 0;
        
        for (int i = 0; i < (x.length - 1); i++) {
            float yTerm = y[i + 1] + y[i];
            float xLen = x[i + 1] - x[i];
            if (xLen < 0) {
                xLen *= -1;
            }
            sum += (yTerm * xLen);
        }
        
        /*
         float xLen = x[i + 1] - x[i];
            
            /*
                 @          @
             @                  @
             |   |          |   |
            -------         -------
            
             |   |          |   |
             @                  @       
                 @          @          
            
            sum += (y[i] * xLen);
            sum += ((y[i + 1] - y[i]) * xLen * 0.5f);
        */
        
        sum *= 0.5;
        
        return sum;
    }
   
    /**
     * @param d
     * @return 
     */
    public static byte[] convertIntToBigEndian(int d) {
        
        if (d == 0) {
            return new byte[]{0};
        }
        if (d == -1) {
            return new byte[]{-1};
        }
        
        int nBytes = 0;
        
        long b = d;
        while ((b != 0) && (b != -1)) {
            nBytes++;
            b >>= 8;
        }
        
        byte[] bytes = new byte[nBytes];
        
        for (int i = 0; i < nBytes; i++) {
            int shift = i * 8;
            int a = (d >> shift) & 255;
            bytes[nBytes - i - 1] = (byte)a;
        }

        return bytes;
    }
    
    /**
     * print Pascal's triangle.  this method handles levels that produce
     * extremely large numbers.  for levels >= 64, it divides a level's numbers
     * by the maximum value to print numbers that should fit within a
     * float.
     * 
     * @param seed
     * @param seedN
     * @param finalN
     * @param printAFactor 
     * @throws java.lang.CloneNotSupportedException 
     */
    public static void printPascalsTriangle2(int[] seed, int seedN, int finalN,
        boolean printAFactor) throws CloneNotSupportedException {
        
        int count = 1;
        
        int finalLen = seed.length + (finalN - seedN);
        
        BigDecimal[] in = new BigDecimal[finalLen];
        BigDecimal[] out = new BigDecimal[finalLen];
        for (int i = 0; i < seed.length; i++) {
            int number = seed[i];
            byte[] numberBytes = convertIntToBigEndian(number);
            in[i] = new BigDecimal(new BigInteger(numberBytes));
            out[i] = BigDecimal.ZERO;
        }
        
        while ((seedN + count) <= finalN) {
            
            out[0] = new BigDecimal(in[0].toString());
            
            BigDecimal max = new BigDecimal(new BigInteger(
                convertIntToBigEndian(Integer.MIN_VALUE)));
            int lastInIdx = (seed.length - 1) + (count - 1);
            int lastOutIdx = lastInIdx + 1;
            for (int i = 0; i < lastInIdx; i++) {
             
                out[i + 1] = new BigDecimal(in[i].toString());
                out[i + 1] = out[i + 1].add(in[i + 1]);
                
                int comp = out[i + 1].compareTo(max);
                if ((out[i + 1].signum() != -1) && (comp > 0)) {
                    max = new BigDecimal(out[i + 1].toString());
                } else if ((out[i + 1].signum() == -1) && (comp < 0)) {
                    // -123456   -123         comp=-1
                    // -123456   1234         comp=-1
                    // -123456   -12345678    comp=1
                    max = new BigDecimal(out[i + 1].toString());
                    max = max.negate();
                }
            }
                        
            out[lastOutIdx] = new BigDecimal(in[lastInIdx].toString());
            
            if (!printAFactor) {
                log.log(Level.INFO, "n={0}:\n{1}", 
                    new Object[]{Integer.valueOf(seedN + count).toString(), 
                        Arrays.toString(out)});
            } else {
                StringBuilder sb = new StringBuilder();
                int indent = 0;
                sb.append("n=").append(seedN + count).append(" has ")
                    .append(lastOutIdx).append(" items:\n");
                for (int ii = 0; ii < indent; ii++) {
                    sb.append(" ");
                }
                
                for (int i = 0; i <= lastOutIdx; i++) {
                    if ((seedN + count) < 32) {
                        sb.append("(float)(").append(out[i]).append("l*a), ");
                    } else {
                        // divide by max to keep numbers small
                        out[i] = out[i].divide(max, 5, BigDecimal.ROUND_DOWN);
                        String divStr = out[i].toString();
                        sb.append("(float)(").append(divStr).append(")*a, ");
                    }
                    
                    if (sb.length() >= (50)) {
                        log.log(Level.INFO, "{0}\n", sb.toString());
                        sb.delete(0, sb.length());
                        for (int ii = 0; ii < indent; ii++) {
                            sb.append(" ");
                        }
                    }
                }
                log.info(sb.toString());
                if ((seedN + count) < 64) {
                    log.log(Level.INFO, " ==> |max|={0}", max);
                }
            }
            
            System.arraycopy(out, 0, in, 0, out.length);
            
            count++;
        }
    }
    
    /**
     * print Pascal's triangle.  this method handles levels that produce
     * extremely large numbers.  for levels >= 64, it divides a level's numbers
     * by the maximum value to print numbers that should fit within a
     * float.
     * 
     * @param seed
     * @param seedN
     * @param finalN
     * @param printAFactor 
     * @throws java.lang.CloneNotSupportedException 
     */
    public static void printPascalsTriangle(int[] seed, int seedN, int finalN,
        boolean printAFactor) throws CloneNotSupportedException {
        
        int count = 1;
        
        int finalLen = seed.length + (finalN - seedN);
        
        VeryLargeNumber[] in = new VeryLargeNumber[finalLen];
        VeryLargeNumber[] out = new VeryLargeNumber[finalLen];
        for (int i = 0; i < seed.length; i++) {
            in[i] = new VeryLargeNumber(seed[i]);
            out[i] = new VeryLargeNumber(0);
        }
     
        while ((seedN + count) <= finalN) {
            
            out[0].resetTo(in[0]);
                   
            VeryLargeNumber max = new VeryLargeNumber(Integer.MIN_VALUE);
            
            int lastInIdx = (seed.length - 1) + (count - 1);
            int lastOutIdx = lastInIdx + 1;
            for (int i = 0; i < lastInIdx; i++) {
             
                out[i + 1] = new VeryLargeNumber(0);
                out[i + 1].add(in[i]);
                out[i + 1].add(in[i + 1]);
                
                int comp = out[i + 1].compareTo(max);
                if (out[i + 1].isPositive() && (comp > 0)) {
                    max = out[i + 1].clone();
                } else if (!out[i + 1].isPositive() && (comp < 0)) {
                    // -123456   -123         comp=-1
                    // -123456   1234         comp=-1
                    // -123456   -12345678    comp=1
                    max = out[i + 1].clone();
                    max.reversePolarity();
                }
            }
                        
            out[lastOutIdx] = new VeryLargeNumber(0);
            out[lastOutIdx].resetTo(in[lastInIdx]);
            
            if (!printAFactor) {
                log.log(Level.INFO, "n={0}:\n{1}", 
                    new Object[]{Integer.valueOf(seedN + count).toString(), 
                        Arrays.toString(out)});
            } else {
                StringBuilder sb = new StringBuilder();
                int indent = 0;
                sb.append("n=").append(seedN + count).append(" has ")
                    .append(lastOutIdx).append(" items:\n");
                for (int ii = 0; ii < indent; ii++) {
                    sb.append(" ");
                }
                
                for (int i = 0; i <= lastOutIdx; i++) {
                    if ((seedN + count) < 32) {
                        sb.append("(float)(").append(out[i]).append("l*a), ");
                    } else {
                        // divide by max to keep numbers small
                        String divStr = out[i].divideByAndPrint(max);
                        sb.append("(float)(").append(divStr).append(")*a, ");
                    }
                    
                    if (sb.length() >= (50)) {
                        log.log(Level.INFO, "{0}\n", sb.toString());
                        sb.delete(0, sb.length());
                        for (int ii = 0; ii < indent; ii++) {
                            sb.append(" ");
                        }
                    }
                }
                log.info(sb.toString());
                if ((seedN + count) < 64) {
                    log.log(Level.INFO, " ==> |max|={0}", max);
                }
            }
            
            System.arraycopy(out, 0, in, 0, out.length);
            
            count++;
        }
    }
    
    public static GreyscaleImage getDeltaDiracImage() {
        
        int w = 7;
        int h = w;
        int xc = w >> 1;
        int yc = xc;
        
        GreyscaleImage img = new GreyscaleImage(w, h);
        
        img.setValue(xc, yc, 255);
                
        return img;
    }
    
    public static GreyscaleImage getCircle() {
        
        double dTheta = 10.0;
        
        int n = (int)(360.f/dTheta);
        
        float r = (float)(((float)n) * 10.f/(2. * Math.PI));//10.0f;
        
        float xc = r + 5;
        float yc = xc;
        
        int w = (int)(xc * 2.f);
        int h = w;
        
        GreyscaleImage img = new GreyscaleImage(w, h);
        
        float expectedCurvature = (1.f/r);
                         
        int pointFactor = 10;
        
        float rDiffMax = 5.0f;
        
        for (int i = 0; i < n*pointFactor; i++) {
            /*
            (x-xc)^2 + (y-yc)^2 = r
            x = xc + r*cos(theta)
            y = yc + r*sin(theta)
            */
            double thetaRadians = (Math.PI*(i*dTheta)/180.)/pointFactor;
            
            double cos = Math.cos(thetaRadians);
            double sin = Math.sin(thetaRadians);
            
            int x = (int)(xc + (r * cos));
            int y = (int)(yc + (r * sin));
            
            img.setValue(x, y, 127);
            
            for (float dr = 0.5f; dr < rDiffMax; dr += 0.25f) {
                img.setValue((int)(xc + ((r - dr) * cos)),
                    (int)(yc + ((r - dr) * sin)), 127);
            }            
        }
        
        return img;
    }
    
    /*
    binomial filter for the second derivative:
   
      n=0? extrapolated            1 -2 1                                    
      n=1? extrapolated          1  -1 -1  1                sigma=sqrt(1)/2 = 0.5
      n=2? extrapolated       1   0  -2  0    1             sigma=sqrt(2)/2 = 0.707
      n=3                   1   1   -2  -2   1   1          sigma=sqrt(3)/2 = 0.866
      n=4                 1   2   -1  -4  -1   2   1        <=== sigma = 1   sigma=sqrt(n)/2
      n=5               1   3   1   -5  -5   1   3    1 
      n=6             1   4   4   -4  -10  -4   4   4    1
      n=7 extr      1   5   8   0  -14   -14    0   8    5    1
      n=8 extr    1   6  13   8   -14  -28   -14   8   13   6   1     <=== sigma=sqrt(2)
      n=9 extr  1  7   19  21   -6  -42  -42    -6   21   19   7    1  <=== sigma=3/2
    
      n=10                 1   8  26  40   15  -48  -84  -48  15  40  26   8   1
      n=11               1   9  34  66   55  -33 -132  -132 -33  55  66  34   9   1
      n=12             1  10  43  100 121  22  -165  -264  -165  22  121  100   43   10   1
      n=13           1  11  53  143 221 143 -143  -429 -429 -143  143  221  143   53   11   1
      n=14         1 12  64  196  364 364  0   -572 -858  -572  0   364 364  196   64  12  1
      n=15       1 13  76  260  560 728  364 -572 -1430  -1430 -572  364  728  560 260 76 13 1
      n=16   * 1 14  89  336 820 1288 1092 -208 -2002 -2860 -2002 -208 1092 1288 820 336 89 14 1 <== sigma=2
    */
    
}
