package algorithms.imageProcessing;

/**
 * @author nichole
 */
public class B3SplineFunction {

    /**
     * <pre>
     * An interpolation function for B-Spline, 3rd order.
     * The implementation follows pseudocode in
     * http://www.multiresolution.com/svbook.pdf
     *
     * The runtime complexity is O(N_pixels).
     *
     * Handling boundaries:
     * "mirror" :      c(k + N) = c(N −k)
     * "periodicity" : (c(k + N) = c(N))
     * "continuity"  : (c(k + N) = c(k))
     * </pre>
    */
    public GreyscaleImage calculate(GreyscaleImage input) {

        int w = input.getWidth();
        int h = input.getHeight();

        GreyscaleImage output = input.copyImage();

        // use separability, that is 1D operation on columns, then rows
        
        for (int row = 0; row < h; ++row) {
            for (int col = 0; col < w; ++col) {
                // choosing "continuity" for boundary corrections
                int x0 = col - 2;
                if (x0 < 0) {
                    x0 = col;
                }
                int x1 = col - 1;
                if (x1 < 0) {
                    x1 = col;
                }
                int x2 = col;
                int x3 = col + 1;
                if (x3 > (w - 1)) {
                    x3 = col;
                }
                int x4 = col + 2;
                if (x4 > (w - 1)) {
                    x4 = col;
                }
                /*
                (1/12)*(|x−2|^3 − 4*|x−1|^3 + 6*|x|^3 − 4*|x+1|^3 + |x+2|^3) 
                1/16, 1/4, 3/8, 1/4, 1/16
                */
                double v0 = input.getValue(x0, row);
                double v1 = input.getValue(x1, row);
                double v2 = input.getValue(x2, row);
                double v3 = input.getValue(x3, row);
                double v4 = input.getValue(x4, row);
                v0 *= (1./16.);
                v1 *= (1./4.);
                v2 *= (3./8.);
                v3 *= (1./4.);
                v4 *= (1./16.);
                int v = (int) Math.round(v0 + v1 + v2 + v3 + v4);

                output.setValue(col, row, v);
            }
        }

        GreyscaleImage input2 = output.copyImage();

        for (int col = 0; col < w; ++col) {
            for (int row = 0; row < h; ++row) {
                // choosing "continuity" for boundary corrections
                int y0 = row - 2;
                if (y0 < 0) {
                    y0 = row;
                }
                int y1 = row - 1;
                if (y1 < 0) {
                    y1 = row;
                }
                int y2 = row;
                int y3 = row + 1;
                if (y3 > (h - 1)) {
                    y3 = row;
                }
                int y4 = row + 2;
                if (y4 > (h - 1)) {
                    y4 = row;
                }
                /*
                (1/12)*(|x−2|^3 − 4*|x−1|^3 + 6*|x|^3 − 4*|x+1|^3 + |x+2|^3) 
                1/16, 1/4, 3/8, 1/4, 1/16
                */
                double v0 = input2.getValue(col, y0);
                double v1 = input2.getValue(col, y1);
                double v2 = input2.getValue(col, y2);
                double v3 = input2.getValue(col, y3);
                double v4 = input2.getValue(col, y4);
                v0 *= (1./16.);
                v1 *= (1./4.);
                v2 *= (3./8.);
                v3 *= (1./4.);
                v4 *= (1./16.);
                int v = (int) Math.round(v0 + v1 + v2 + v3 + v4);

                output.setValue(col, row, v);
            }
        }

        return output;
    }

}
