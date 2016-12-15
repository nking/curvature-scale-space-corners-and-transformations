package algorithms.imageProcessing;

/**
 * 
 * @author nichole
 */
public class SummedColumnTable {
    
    /**
     * sum along columns, that is a[i][*]
     * 
     * @param a
     * @return 
     */
    public float[][] create(float[][] a) {

        int w = a.length;
        int h = a[0].length;
        
        ImageProcessor imp = new ImageProcessor();

        float[][] out = imp.copy(a);
        //applyAbsoluteValue
        for (int x = 0; x < w; ++x) {
            for (int y = 0; y < h; ++y) {
                float v = out[x][y];
                if (v < 0) {
                    out[x][y] *= -1;
                }
            }
        }
        
        // sum along columns, that is a[i][*]
        
        for (int i = 0; i < w; ++i) {
            for (int j = 0; j < h; ++j) {
                if (j > 0) {
                    out[i][j] += out[i][j - 1];
                }
            }
        }

        return out;
    }
    
    /**
     * extract the sum of a window along a columnn bound by given start and stop 
     * coordinates and return that value and the number of pixels in the
     * window in the output variable, output.
     * for example, extracting columns 3 through 5 of row 2 is a[2][3:5]
     * @param input
     * @param start coordinate for x start of window, inclusive
     * @param stop coordinate for x stop of window, inclusive
     * @param row coordinate for row, that is first dimension of input array
     * @param output one dimensional array of size 2 in which the
     * sum of the window will be returned and the number of pixels in the 
     * window.  float[]{sum, nPixels}
     */
    public void extractWindowInColumn(float[][] input, 
        int start, int stop, int row, float output[]) {
        
        if (start == 0) {
            output[0] = input[row][stop];
        } else {
            output[0] = input[row][stop] - input[row][start - 1];
        }
        output[1] = stop - start + 1;        
    }
    
}
