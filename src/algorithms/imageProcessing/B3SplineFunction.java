package algorithms.imageProcessing;

import java.util.Arrays;

/**
 * @author nichole
 */
public class B3SplineFunction {

    /**
     * <pre>
     * An interpolation function for B-Spline, 3rd order.
     * The implementation follows pseudocode in
     * http://www.multiresolution.com/svbook.pdf
     * "Handbook of Astronomical Data Analysis" by 
     * Jean-Luc Starck and Fionn Murtagh
     * 
     * The runtime complexity is O(N_pixels) and internally uses 2 1-D operations.
     *
     * Handling boundaries:
     * "mirror" :      c(k + N) = c(N −k)
     * "periodicity" : (c(k + N) = c(N))
     * "continuity"  : (c(k + N) = c(k))
     * </pre>
     * @param input
     * @return
    */
    public GreyscaleImage calculate(GreyscaleImage input) {

        int w = input.getWidth();
        int h = input.getHeight();

        GreyscaleImage output = input.copyImage();

        // use separability, that is 1D operation on columns, then rows

        for (int row = 0; row < h; ++row) {
            for (int col = 0; col < w; ++col) {

                // choosing "continuity" for boundary corrections

                int vSum = 0;
                for (int dx = -2; dx <= 2; ++dx) {

                    int xi = col + dx;
                    if ((xi < 0) || (xi > (w - 1))) {
                        xi = col;
                    }

                    int v = input.getValue(xi, row);

                    /*
                     (1/12)*(|x−2|^3 − 4*|x−1|^3 + 6*|x|^3 − 4*|x+1|^3 + |x+2|^3)
                     1/16, 1/4, 3/8, 1/4, 1/16
                     1/16*(1, 4, 6, 4, 1)
                     */
                    switch (dx) {
                        // -2 and +2
                        case -1:
                        case 1:
                            v <<= 2;
                            break;
                        case 0:
                            v <<= 1;
                            v *= 3;
                            break;
                        // case -2 and +2 are factor 1
                        default:
                            break;
                    }

                    vSum += v;
                }

                vSum >>= 4;

                output.setValue(col, row, vSum);
            }
        }

        GreyscaleImage input2 = output.copyImage();

        for (int col = 0; col < w; ++col) {
            for (int row = 0; row < h; ++row) {

                // choosing "continuity" for boundary corrections

                int vSum = 0;
                for (int dy = -2; dy <= 2; ++dy) {

                    int yi = row + dy;
                    if ((yi < 0) || (yi > (h - 1))) {
                        yi = row;
                    }

                    int v = input2.getValue(col, yi);

                    /*
                     (1/12)*(|x−2|^3 − 4*|x−1|^3 + 6*|x|^3 − 4*|x+1|^3 + |x+2|^3)
                     1/16, 1/4, 3/8, 1/4, 1/16
                     1/16*(1, 4, 6, 4, 1)
                     */
                    switch (dy) {
                        // -2 and +2
                        case -1:
                        case 1:
                            v <<= 2;
                            break;
                        case 0:
                            v <<= 1;
                            v *= 3;
                            break;
                        // case -2 and +2 are factor 1
                        default:
                            break;
                    }

                    vSum += v;
                }

                vSum >>= 4;

                output.setValue(col, row, vSum);
            }
        }

        return output;
    }

    protected int interpolate1D(int x, int y, GreyscaleImage img,
        boolean calcForX) {

        if (calcForX) {
            return interpolate1DX(x, y, img);
        } else {
            return interpolate1DY(x, y, img);
        }
    }

    /**
     * interpolate values around (x,y) along x in img using a B3 spline.
     *
     * @param x
     * @param y
     * @param img
     * @return
     */
    public int interpolate1DX(int x, int y, GreyscaleImage img) {

        int w = img.getWidth();
        int h = img.getHeight();

        /*
        (1/12)*(|x−2|^3 − 4*|x−1|^3 + 6*|x|^3 − 4*|x+1|^3 + |x+2|^3)
        1/16, 1/4, 3/8, 1/4, 1/16
        1/16*(1, 4, 6, 4, 1)
        */
        int vSum = 0;
        for (int dx = -2; dx <= 2; ++dx) {

            int xi = x + dx;
            if ((xi < 0) || (xi > (w - 1))) {
                xi = x;
            }

            int v = img.getValue(xi, y);

            switch(dx) {
                // -2 and +2
                case -1:
                case 1:
                    v <<= 2;
                    break;
                case 0:
                    v <<= 1;
                    v *= 3;
                    break;
                // case -2 and +2 are factor 1
                default:
                    break;
            }

            vSum += v;
        }

        vSum >>= 4;

        return vSum;
    }

    /**
     * interpolate values around (x,y) along y in img using a B3 spline.
     *
     * @param x
     * @param y
     * @param img
     * @return
     */
    protected int interpolate1DY(int x, int y, GreyscaleImage img) {

        int w = img.getWidth();
        int h = img.getHeight();

        /*
        (1/12)*(|x−2|^3 − 4*|x−1|^3 + 6*|x|^3 − 4*|x+1|^3 + |x+2|^3)
        1/16, 1/4, 3/8, 1/4, 1/16
        
        1/16*(1, 4, 6, 4, 1)
        */
        int vSum = 0;
        for (int dy = -2; dy <= 2; ++dy) {

            int yi = y + dy;
            if ((yi < 0) || (yi > (h - 1))) {
                yi = y;
            }

            int v = img.getValue(x, yi);

            switch(dy) {
                // -2 and +2
                case -1:
                case 1:
                    v <<= 2;
                    break;
                case 0:
                    v <<= 1;
                    v *= 3;
                    break;
                // case -2 and +2 are factor 1
                default:
                    break;
            }

            vSum += v;
        }
        
        vSum >>= 4;

        return vSum;
    }

    /**
     * interpolate values around (x,y) along y in img using a B3 spline. 
     *
     * @param x
     * @param y
     * @param img
     * @return
     */
    protected int interpolate2D(int x, int y, GreyscaleImage img) {

        int w = img.getWidth();
        int h = img.getHeight();

        /*        
         1/256  1/64   3/128    1/64    1/256
         1/64   1/16   3/32     1/16    1/64
         3/128  3/32   9/64     3/32    3/128
         1/64   1/16   3/32     1/16    1/64
         1/256  1/64   3/128    1/64    1/256

         (1/256) *   |   1    4   3*2    4    1 |
                     |   4   16   3*8   16    4 |
                     | 3*2  3*8   9*4  3*8  3*2 |
                     |   4   16   3*8   16    4 |
                     |   1    4   3*2    4    1 |
        */

        int vSum = 0;

        for (int dy = -2; dy <= 2; ++dy) {
            int yi = y + dy;
            if ((yi < 0) || (yi > (h - 1))) {
                yi = y;
            }
            for (int dx = -2; dx <= 2; ++dx) {
                int xi = x + dx;
                if ((xi < 0) || (xi > (w - 1))) {
                    xi = x;
                }
                int v = img.getValue(xi, yi);
                switch(dx) {
                    case -2:
                    case 2:
                        switch(dy) {
                            case -1:
                            case 1:
                                v <<= 2;
                                break;
                            case 0:
                                v <<= 1;
                                v *= 3;
                                break;
                            // rows -2 and +2 are factors of 1
                            default:
                                break;
                        }
                        break;
                    case -1:
                    case 1:
                        switch(dy) {
                            // -2, 2 are 4
                            case -2:
                            case 2:
                                v <<= 2;
                                break;
                            // rows -1, 1 are 16
                            case -1:
                            case 1:
                                v <<= 4;
                                break;
                            case 0:
                                v <<= 3;
                                v *= 3;
                                break;
                            default:
                                break;
                        }
                        break;
                    case 0:
                        switch(dy) {
                            // -2, 2 are 6
                            case -2:
                            case 2:
                                v <<= 1;
                                v *= 3;
                                break;
                            // rows -1, 1 are 24
                            case -1:
                            case 1:
                                v <<= 3;
                                v *= 3;
                                break;
                            case 0:
                                v <<= 2;
                                v *= 9;
                                break;
                            default:
                                break;
                        }
                        break;
                    default:
                        break;
                } // end switch(dx)
                vSum += v;
            }
        }

        vSum >>= 8;

        return vSum;
    }

    /**
     * calculate the B3 Spline for every pixel in the image.  The runtime
     * complexity is roughly linear, but is 2.5 times larger than the
     * calculate(img) which uses two 1-D splines for each pixel.
     * @param img
     * @return 
     */
    protected GreyscaleImage calculate2D(GreyscaleImage img) {

        int w = img.getWidth();
        int h = img.getHeight();

        /*
         1/256  1/64   3/128    1/64    1/256
         1/64   1/16   3/32     1/16    1/64
         3/128  3/32   9/64     3/32    3/128
         1/64   1/16   3/32     1/16    1/64
         1/256  1/64   3/128    1/64    1/256

         (1/256) *   |   1    4   3*2    4    1 |
                     |   4   16   3*8   16    4 |
                     | 3*2  3*8   9*4  3*8  3*2 |
                     |   4   16   3*8   16    4 |
                     |   1    4   3*2    4    1 |
        */

        GreyscaleImage output = img.createWithDimensions();

        for (int row = 0; row < h; ++row) {
            for (int col = 0; col < w; ++col) {

                int v = interpolate2D(col, row, img);

                output.setValue(col, row, v);
            }
        }

        return output;
    }

    /**
     * <pre>
     * An interpolation function for B-Spline, 3rd order.
     * The implementation follows pseudocode in
     * http://www.multiresolution.com/svbook.pdf
     * "Handbook of Astronomical Data Analysis" by 
     * Jean-Luc Starck and Fionn Murtagh
     * 
     * The runtime complexity is O(N_pixels) and internally uses 2 1-D operations.
     *
     * Handling boundaries:
     * "mirror" :      c(k + N) = c(N −k)
     * "periodicity" : (c(k + N) = c(N))
     * "continuity"  : (c(k + N) = c(k))
     * </pre>
     * 
     * Note that the method depends upon logic for transforming pixel coordinates
     * x and y into single array indexes that is present in GreyscaleImage,
     * so if that ever changes, the same changes need to be made here.
     * (TODO: refactor for pixel coordinate transformation).
     * 
     * @param input
     * @param imgWidth
     * @param imgHeight
     * @return
    */
    public double[] calculate(double[] input, int imgWidth, int imgHeight) {
        
        int w = imgWidth;
        int h = imgHeight;
        
        double[] output = Arrays.copyOf(input, input.length);

        // use separability, that is 1D operation on columns, then rows

        for (int row = 0; row < h; ++row) {
            for (int col = 0; col < w; ++col) {

                // choosing "continuity" for boundary corrections

                double vSum = 0;
                for (int dx = -2; dx <= 2; ++dx) {

                    int xi = col + dx;
                    if ((xi < 0) || (xi > (w - 1))) {
                        xi = col;
                    }
                    
                    int pixIdx = (row * imgWidth) + xi;

                    double v = input[pixIdx];

                    /*
                     (1/12)*(|x−2|^3 − 4*|x−1|^3 + 6*|x|^3 − 4*|x+1|^3 + |x+2|^3)
                     1/16, 1/4, 3/8, 1/4, 1/16
                     1/16*(1, 4, 6, 4, 1)
                     */
                    switch (dx) {
                        // -2 and +2
                        case -1:
                        case 1:
                            v *= 4;
                            break;
                        case 0:
                            v *= 6;
                            break;
                        // case -2 and +2 are factor 1
                        default:
                            break;
                    }

                    vSum += v;
                }

                vSum /= 16.;

                int pixIdx = (row * imgWidth) + col;
                
                output[pixIdx] = vSum;
            }
        }

        double[] input2 = Arrays.copyOf(output, output.length);

        for (int col = 0; col < w; ++col) {
            for (int row = 0; row < h; ++row) {

                // choosing "continuity" for boundary corrections

                double vSum = 0;
                for (int dy = -2; dy <= 2; ++dy) {

                    int yi = row + dy;
                    if ((yi < 0) || (yi > (h - 1))) {
                        yi = row;
                    }

                    int pixIdx = (yi * imgWidth) + col;

                    double v = input2[pixIdx];
                    
                    /*
                     (1/12)*(|x−2|^3 − 4*|x−1|^3 + 6*|x|^3 − 4*|x+1|^3 + |x+2|^3)
                     1/16, 1/4, 3/8, 1/4, 1/16
                     1/16*(1, 4, 6, 4, 1)
                     */
                    switch (dy) {
                        // -2 and +2
                        case -1:
                        case 1:
                            v *= 4;
                            break;
                        case 0:
                            v *= 6;
                            break;
                        // case -2 and +2 are factor 1
                        default:
                            break;
                    }

                    vSum += v;
                }

                vSum /= 16.;

                int pixIdx = (row * imgWidth) + col;
                output[pixIdx] = vSum;
            }
        }

        return output;
    }
}
