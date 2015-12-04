package algorithms.imageProcessing;

/**
 * @author nichole
 */
public class TriangleFunction {

    /**
     * <pre>
     * An interpolation function called the triangle function
     * that uses a base 2 spacing. The implementation follows pseudocode in
     * http://www.multiresolution.com/svbook.pdf
     *
     * The runtime complexity is O(N_pixels).
     *
     * c_(j + 1,k) = (1/4)*c_(j,k-(2^j)) + (1/2)*c_(j,k) + (1/4)*c_(j,k+(2^j))
     *
     * Handling boundaries:
     * "mirror" :      c(k + N) = c(N −k)
     * "periodicity" : (c(k + N) = c(N))
     * "continuity"  : (c(k + N) = c(k))
     * </pre>
     *
     * @param input
     * @param j level associated with input image. The output is calculated
     * using 2^j as spacing for interpolation points.
     * @return a sampling of input, interpolated over spacings 2^j.
     */
    public GreyscaleImage calculateNextLevel(GreyscaleImage input, int j) {

        return addOrSubtract(input, j, true);
    }
    
    /**
     * <pre>
     * An interpolation function called the triangle function
     * that uses a base 2 spacing to subtract to transformed levels. 
     * The implementation follows pseudocode in
     * http://www.multiresolution.com/svbook.pdf
     *
     * The runtime complexity is O(N_pixels).
     *
     * w_(j+1,k) = c_(j,k) − c_(j+1,k)
     *           = (-1/4)*c_(j,k-(2^j)) + (1/2)*c_(j,k) - (1/4)*c_(j,k+(2^j))
     *
     * Handling boundaries:
     * "mirror" :      c(k + N) = c(N −k)
     * "periodicity" : (c(k + N) = c(N))
     * "continuity"  : (c(k + N) = c(k))
     * </pre>
     *
     * @param input
     * @param j level associated with input image. The output is calculated
     * using 2^j as spacing for interpolation points.
     * @return a sampling of input, interpolated over spacings 2^j.
     */
    public GreyscaleImage subtractLevels(GreyscaleImage input, int j) {

        return addOrSubtract(input, j, false);
    }

    private GreyscaleImage addOrSubtract(GreyscaleImage input, int j, boolean add) {

        int w = input.getWidth();
        int h = input.getHeight();

        int s = 1 << j;

        GreyscaleImage output = input.copyImage();

        // use separability, that is 1D operation on columns, then rows
        
        for (int row = 0; row < h; ++row) {
            for (int col = 0; col < w; ++col) {
                int x0 = col - s;
                int x2 = col + 2;

                // choosing "continuity" for boundary corrections
                if (x0 < 0) {
                    x0 = col;
                }
                if (x2 > (w - 1)) {
                    x2 = col;
                }

                // add:
                //    c_(j + 1,k) = (1/4)*c_(j,k-(2^j)) + (1/2)*c_(j,k) + (1/4)*c_(j,k+(2^j))
                // subtract:
                //    w_(j+1,k) = (-1/4)*c_(j,k-(2^j)) + (1/2)*c_(j,k) - (1/4)*c_(j,k+(2^j))
                double v0 = 0.25 * input.getValue(x0, row);
                double v1 = 0.5 * input.getValue(col, row);
                double v2 = 0.25 * input.getValue(x2, row);
                double vSum;
                if (add) {
                    vSum = v0 + v1 + v2;
                } else {
                    vSum = -1*v0 + v1 - v2;
                }
                int v = (int) Math.round(vSum);

                output.setValue(col, row, v);
            }
        }

        GreyscaleImage input2 = output.copyImage();

        for (int col = 0; col < w; ++col) {
            for (int row = 0; row < h; ++row) {
                int y0 = row - s;
                int y2 = row + 2;

                // choosing "continuity" for boundary corrections
                if (y0 < 0) {
                    y0 = row;
                }
                if (y2 > (h - 1)) {
                    y2 = row;
                }

                // add:
                //    c_(j + 1,k) = (1/4)*c_(j,k-(2^j)) + (1/2)*c_(j,k) + (1/4)*c_(j,k+(2^j))
                // subtract:
                //    w_(j+1,k) = (-1/4)*c_(j,k-(2^j)) + (1/2)*c_(j,k) - (1/4)*c_(j,k+(2^j))
                double v0 = 0.25 * input2.getValue(col, y0);
                double v1 = 0.5 * input2.getValue(col, row);
                double v2 = 0.25 * input2.getValue(col, y2);
                double vSum;
                if (add) {
                    vSum = v0 + v1 + v2;
                } else {
                    vSum = -1*v0 + v1 - v2;
                }
                int v = (int) Math.round(vSum);

                output.setValue(col, row, v);
            }
        }

        return output;
    }

}
