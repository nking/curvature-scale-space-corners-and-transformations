package algorithms.imageProcessing;

import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.PairIntArrayWithColor;
import java.util.Set;

/**
 *
 * @author nichole
 */
public class Kernel1DHelper {
   
    /**
     * convolve point[xyIdx] with the kernel g.  if calcX is true, the 
     * point xPoints[xyIdx] is used, else yPoints[xyIdx].
     * minPointsValue is used for 'mu', the location parameter in the
     * Gaussian during the convolution.  Note that for extremely large
     * range of point values, minPointsValue should probably the median
     * of the points.
     * @param curve
     * @param cIndex
     * @param g
     * @param calcX
     * @return 
     */
    public double convolvePointWithKernel(PairIntArray curve, int cIndex,
        float[] g, final boolean calcX) {

        int h = g.length >> 1;

        double sum = 0;

        int curveLength = curve.getN();
        
        boolean isClosedCurved = (curve instanceof PairIntArrayWithColor)
            && ((PairIntArrayWithColor) curve).isClosedCurve() &&
            (curve.getN() > 2);
            
        /*
        correct for the edges.
        Can try:
            (1) wrap around - If this is a closed curve, use the points on the
                 other end to continue matching the kernel point for point.
                 ** this is what will use for closed curves **
            (2) expand/pad - expand to fill the area large enough to match
                the entire filter by:
                    -- replicating the last point
                    -- or reflecting the previous points around this boundary
                       ** this is the one will use for open curves **
                    -- or pad with zeroes.  (this one does NOT work well  it's
                       the same as 'crop the area')
            (3) crop the area.  this involves re-normalization of the percentage
                 of the kernel which was used and can introduce large errors
                 thru division by a very small number for example, so don't use
                 it if possible.
        */
                
        for (int gIdx = 0; gIdx < g.length; gIdx++) {

            float gg = g[gIdx];

            if (gg == 0) {
                continue;
            }
            
            int x = gIdx - h;

            int curveIdx = cIndex + x;

            if (curveIdx < 0) {
                
                if (isClosedCurved) {
                    // wrap around
                    /*
                    n-4 n-3 n-2 n-1  0   1   2   3   .   .   .  n-1
                     _   _   _   _   @   @   @   @   @   @   @   @
                    */
                    // for the rare case when the kernel is so much larger than
                    //   the edge, will need to iterate
                    while (curveIdx < 0) {
                        curveIdx = curve.getN() + curveIdx;
                    }
                } else {
                    //TODO: revisit this for range of kernel sizes vs edge sizes
                    // replicate
                    curveIdx = -1*curveIdx - 1;
                    if (curveIdx > (curveLength - 1)) {
                        curveIdx = curveLength - 1;
                    }
                }
            } else if (curveIdx >= (curveLength)) {
                if (isClosedCurved) {
                    // wrap around
                    /*
                    0   1   2   3   .   .   .  n-1  0   1   2   3   4
                    @   @   @   @   @   @   @   @   _   _   _   _   _
                    */
                    // for the rare case when the kernel is so much larger than
                    //   the edge, will need to iterate
                    while (curveIdx >= (curveLength)) {
                        curveIdx = curveIdx - curve.getN();
                    }
                } else {
                    //TODO: revisit this for range of kernel sizes vs edge sizes
                    int diff = curveIdx - curveLength;
                    curveIdx = curveLength - diff - 1;
                    if (curveIdx < 0) {
                        curveIdx = 0;
                    }
                }
            }
            
            float point;
            
            if (calcX) {
                point = curve.getX(curveIdx);
            } else {
                point = curve.getY(curveIdx);
            }

            sum += (point * gg);
        }

        return sum;
    }

    /**
     * convolve curve[cIndex] with the kernel g.
     * @param curve
     * @param cIndex
     * @param g
     * @return 
     */
    public float convolvePointWithKernel(float[] curve, int cIndex,
        float[] g) {

        int h = g.length >> 1;

        float sum = 0;

        int curveLength = curve.length;
        
        /*
        correct for the edges.
        Can try:
            (1) wrap around - If this is a closed curve, use the points on the
                 other end to continue matching the kernel point for point.
                 ** this is what will use for closed curves **
            (2) expand/pad - expand to fill the area large enough to match
                the entire filter by:
                    -- replicating the last point
                    -- or reflecting the previous points around this boundary
                       ** this is the one will use for open curves **
                    -- or pad with zeroes.  (this one does NOT work well  it's
                       the same as 'crop the area')
            (3) crop the area.  this involves re-normalization of the percentage
                 of the kernel which was used and can introduce large errors
                 thru division by a very small number for example, so don't use
                 it if possible.
        */
                
        for (int gIdx = 0; gIdx < g.length; gIdx++) {

            float gg = g[gIdx];

            if (gg == 0) {
                continue;
            }
            
            int x = gIdx - h;

            int curveIdx = cIndex + x;

            if (curveIdx < 0) {
                //TODO: revisit this for range of kernel sizes vs edge sizes
                // replicate
                curveIdx = -1*curveIdx - 1;
                if (curveIdx > (curveLength - 1)) {
                    curveIdx = curveLength - 1;
                }
            } else if (curveIdx >= (curveLength)) {
                //TODO: revisit this for range of kernel sizes vs edge sizes
                int diff = curveIdx - curveLength;
                curveIdx = curveLength - diff - 1;
                if (curveIdx < 0) {
                    curveIdx = 0;
                }
            }
            
            float point = curve[curveIdx];
            
            sum += (point * gg);
        }

        return sum;
    }
    
    /**
     * convolve point[xyIdx] with the kernel g along a column if calcX is true,
     * else along a row if calcX is false.
     * @param img
     * @param col
     * @param row
     * @param g
     * @param calcX convolve along column if true, else row
     * @return 
     */
    public double convolvePointWithKernel(final GreyscaleImage img, int col, 
        int row, float[] g, final boolean calcX) {

        int h = g.length >> 1;

        double sum = 0;

        int len = calcX ? img.getWidth() : img.getHeight();
                
        for (int gIdx = 0; gIdx < g.length; gIdx++) {

            float gg = g[gIdx];

            if (gg == 0) {
                continue;
            }
            
            int idx = gIdx - h;

            int cIdx = calcX ? (col + idx) : (row + idx);

            if (cIdx < 0) {
                // replicate
                cIdx = -1*cIdx - 1;
                if (cIdx > (len - 1)) {
                    cIdx = len - 1;
                }
            } else if (cIdx >= (len)) {
                //TODO: revisit this for range of kernel sizes vs edge sizes
                int diff = cIdx - len;
                cIdx = len - diff - 1;
                if (cIdx < 0) {
                    cIdx = 0;
                }
            }
            
            float point;
            
            if (calcX) {
                // keep row constant
                point = img.getValue(cIdx, row);
            } else {
                // keep col constant
                point = img.getValue(col, cIdx);
            }

            sum += (point * gg);
        }

        return sum;
    }

/**
     * convolve point[xyIdx] with the kernel g along a column if calcX is true,
     * else along a row if calcX is false.
     * @param img
     * @param col
     * @param row
     * @param g
     * @param calcX convolve along column if true, else row
     * @return 
     */
    public double convolvePointWithKernel(final int[][] img, int col, 
        int row, float[] g, final boolean calcX) {

        int h = g.length >> 1;

        double sum = 0;

        int len = calcX ? img.length : img[0].length;
                
        for (int gIdx = 0; gIdx < g.length; gIdx++) {

            float gg = g[gIdx];

            if (gg == 0) {
                continue;
            }
            
            int idx = gIdx - h;

            int cIdx = calcX ? (col + idx) : (row + idx);

            if (cIdx < 0) {
                // replicate
                cIdx = -1*cIdx - 1;
                if (cIdx > (len - 1)) {
                    cIdx = len - 1;
                }
            } else if (cIdx >= (len)) {
                //TODO: revisit this for range of kernel sizes vs edge sizes
                int diff = cIdx - len;
                cIdx = len - diff - 1;
                if (cIdx < 0) {
                    cIdx = 0;
                }
            }
            
            float point;
            
            if (calcX) {
                // keep row constant
                point = img[cIdx][row];
            } else {
                // keep col constant
                point = img[col][cIdx];
            }

            sum += (point * gg);
        }

        return sum;
    }
    
    /**
     * convolve point[xyIdx] with the kernel g along a column if calcX is true,
     * else along a row if calcX is false.
     * @param imageValues values from image stored in
     * an array that uses same indexing as in Image.java.
     * @param imageWidth
     * @param imageHeight
     * @param col
     * @param row
     * @param g
     * @param calcX convolve along column if true, else row
     * @return 
     */
    public double convolvePointWithKernel(
        final int[] imageValues, int imageWidth,
        int imageHeight,
        int col, int row, float[] g, final boolean calcX) {

        int h = g.length >> 1;

        double sum = 0;
        
        int len = calcX ? imageWidth : imageHeight;
                
        for (int gIdx = 0; gIdx < g.length; gIdx++) {

            float gg = g[gIdx];

            if (gg == 0) {
                continue;
            }
            
            int idx = gIdx - h;

            int cIdx = calcX ? (col + idx) : (row + idx);

            if (cIdx < 0) {
                // replicate
                cIdx = -1*cIdx - 1;
                if (cIdx > (len - 1)) {
                    cIdx = len - 1;
                }
            } else if (cIdx >= (len)) {
                //TODO: revisit this for range of kernel sizes vs edge sizes
                int diff = cIdx - len;
                cIdx = len - diff - 1;
                if (cIdx < 0) {
                    cIdx = 0;
                }
            }
            
            double point;
            
            if (calcX) {
                int pixIdx = (row * imageWidth) + cIdx;
                // keep row constant
                point = imageValues[pixIdx];
            } else {
                int pixIdx = (cIdx * imageWidth) + col;
                // keep col constant
                point = imageValues[pixIdx];
            }

            sum += (point * gg);
        }

        return sum;
    }
    
    /**
     * convolve point[xyIdx] with the kernel g along a column if calcX is true,
     * else along a row if calcX is false.
     * @param img
     * @param col
     * @param row
     * @param g
     * @param calcX convolve along column if true, else row
     * @return 
     */
    public double convolvePointWithKernel(final double[][] img, int col, 
        int row, float[] g, final boolean calcX) {

        int h = g.length >> 1;

        double sum = 0;

        int len = calcX ? img.length : img[0].length;
                
        for (int gIdx = 0; gIdx < g.length; gIdx++) {

            float gg = g[gIdx];

            if (gg == 0) {
                continue;
            }
            
            int idx = gIdx - h;

            int cIdx = calcX ? (col + idx) : (row + idx);

            if (cIdx < 0) {
                // replicate
                cIdx = -1*cIdx - 1;
                if (cIdx > (len - 1)) {
                    cIdx = len - 1;
                }
            } else if (cIdx >= (len)) {
                //TODO: revisit this for range of kernel sizes vs edge sizes
                int diff = cIdx - len;
                cIdx = len - diff - 1;
                if (cIdx < 0) {
                    cIdx = 0;
                }
            }
            
            double point;
            
            if (calcX) {
                // keep row constant
                point = img[cIdx][row];
            } else {
                // keep col constant
                point = img[col][cIdx];
            }

            sum += (point * gg);
        }

        return sum;
    }
    
    /**
     * convolve point[xyIdx] with the kernel g along a column if calcX is true.
     * @param img
     * @param col
     * @param row
     * @param g
     * @param calcX
     * @return double[]{rSum, gSum, bSum]
     */
    public double[] convolvePointWithKernel(final Image img, int col, 
        int row, float[] g, final boolean calcX) {

        int h = g.length >> 1;

        double[] sum = new double[3];        

        int len = calcX ? img.getWidth() : img.getHeight();
        
        for (int color = 0; color < 3; color++) {
                
            for (int gIdx = 0; gIdx < g.length; gIdx++) {

                float gg = g[gIdx];

                if (gg == 0) {
                    continue;
                }

                int idx = gIdx - h;

                int cIdx = calcX ? (col + idx) : (row + idx);

                if (cIdx < 0) {
                    // replicate
                    cIdx = -1*cIdx - 1;
                    if (cIdx > (len - 1)) {
                        cIdx = len - 1;
                    }
                } else if (cIdx >= (len)) {
                    //TODO: revisit this for range of kernel sizes vs edge sizes
                    int diff = cIdx - len;
                    cIdx = len - diff - 1;
                    if (cIdx < 0) {
                        cIdx = 0;
                    }
                }

                float point;

                if (calcX) {
                    // keep row constant
                    if (color == 0) {
                        point = img.getR(cIdx, row);
                    } else if (color == 1) {
                        point = img.getG(cIdx, row);
                    } else {
                        point = img.getB(cIdx, row);
                    }
                } else {
                    // keep col constant
                    if (color == 0) {
                        point = img.getR(col, cIdx);
                    } else if (color == 1) {
                        point = img.getG(col, cIdx);
                    } else {
                        point = img.getB(col, cIdx);
                    }
                }

                sum[color] += (point * gg);
            }
        }
        
        return sum;
    }
    
    /**
     * convolve point[xyIdx] with the kernel g along a column if calcX is true,
     * else along a row if calcX is false.
     * @param points
     * @param col
     * @param row
     * @param g
     * @param calcX convolve along column if true, else row
     * @return 
     */
    public double convolvePointWithKernel(final Set<PairInt> points, int col, 
        int row, float[] g, final boolean calcX) {

        int h = g.length >> 1;

        double sum = 0;
                
        for (int gIdx = 0; gIdx < g.length; gIdx++) {

            float gg = g[gIdx];

            if (gg == 0) {
                continue;
            }
            
            int idx = gIdx - h;

            int cIdx = calcX ? (col + idx) : (row + idx);

            float point = 0;
            
            if (calcX) {
                // keep row constant
                if (points.contains(new PairInt(cIdx, row))) {
                    point = 1;
                }
            } else {
                // keep col constant
                if (points.contains(new PairInt(col, cIdx))) {
                    point = 1;
                }
            }

            sum += (point * gg);
        }

        return sum;
    }
}
