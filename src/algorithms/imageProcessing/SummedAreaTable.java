package algorithms.imageProcessing;

/**
 * An algorithm to make cumulative sums at a pixel adding to it
 * the value from the pixel
 * below and the value from left of it.
 * The extraction of a window of any size throughout the image
 * then only takes 4 steps at most.
 * The runtime complexity for building the summed area table
 * is O(N) and the extraction of the sum of any size window
 * centered on a point is O(1).
 * The space complexity is O(N).
 * 
 * From https://en.wikipedia.org/wiki/Summed_area_table
 * The algorithm was introduced to computer graphics in 1984 by Frank Crow 
 * for use with mipmaps. In computer vision it was popularized by 
 * Lewis[1] and then given the name "integral image" and prominently 
 * used within the Viola–Jones object detection framework in 2001.
 * 
 * @author nichole
 */
public class SummedAreaTable {
    
    public GreyscaleImage createAbsoluteSummedAreaTable(GreyscaleImage img) {

        int w = img.getWidth();
        int h = img.getHeight();

        GreyscaleImage out = img.copyToFullRangeIntImage();
        applyAbsoluteValue(out);
        
        for (int x = 0; x < w; ++x) {
            for (int y = 0; y < h; ++y) {
                if (x > 0 && y > 0) {
                    int v = out.getValue(x - 1, y) + out.getValue(x, y - 1) 
                        - out.getValue(x - 1, y - 1);
                    out.setValue(x, y, out.getValue(x, y) + v);
                } else if (x > 0) {
                    int v = out.getValue(x - 1, y);
                    out.setValue(x, y, out.getValue(x, y) + v);
                } else if (y > 0) {
                    int v = out.getValue(x, y - 1);
                    out.setValue(x, y, out.getValue(x, y) + v);
                }
            }
        }

        return out;
    }
    
     /**
     * @param imgS
     * @param d
     * @return 
     */
    public GreyscaleImage applyMeanOfWindowFromSummedAreaTable(
        GreyscaleImage imgS, int d) {
        
        int w = imgS.getWidth();
        int h = imgS.getHeight();
        
        GreyscaleImage img2 = imgS.createFullRangeIntWithDimensions();
        
        int[] sumAndN = new int[2];
        
        // extract the summed area of each dxd window centered on x,y
        // and divide by number of pixels
        for (int x = 0; x < w; ++x) {
            for (int y = 0; y < h; ++y) {
                extractWindowFromSummedAreaTable(imgS, x, y, d, sumAndN);
                int v = sumAndN[0]/sumAndN[1];
                img2.setValue(x, y, v);
            }
        }

        return img2;
    }
    
    /**
     * extract the sum of a window centered at (x,y) of x dimension d and y
     * dimension d and return that value and the number of pixels in the
     * aperture in the output variable, output.
     * NOTE GreyscaleImage, x, and y are in column major format
     * @param imgS
     * @param x coordinate for x center of window
     * @param y coordinate for y center of window
     * @param d diameter of window in x and y
     * @param output one dimensional array of size 2 in which the
     * sum of the window will be returned and the number of pixels in the 
     * window.  int[]{sum, nPixels}
     */
    public void extractWindowFromSummedAreaTable(GreyscaleImage imgS, 
        int x, int y, int d, int output[]) {
        
        if (output == null || output.length != 2) {
            throw new IllegalArgumentException(
                "output must be initialized to size 2");
        }
        
        if (d < 0) {
            throw new IllegalArgumentException(
                "d must be a non-negative number");
        }
        
        int w = imgS.getWidth();
        int h = imgS.getHeight();
        
        if (x < 0 || y < 0 || (x > (w - 1)) || (y > (h - 1))) {
            throw new IllegalArgumentException("x or y is out of bounds of "
                + "image. x=" + x + " y=" + y + " w=" + w + " h=" + h);
        }
        
        final int r = (d >> 1);
        
        // extract the summed area of dxd window centered on x,y
        if (r > 0) {
            if (x > r && x < (w-r) && (y > r) && (y < (h-r))) {
                int nPix = d * d;
                int s1 = imgS.getValue(x+r, y+r) - imgS.getValue(x-r, y+r)
                    - imgS.getValue(x+r, y-r) + imgS.getValue(x-r, y-r);
                output[0] = s1;
                output[1] = nPix;
                return;
            }
        }
                
        // handling borders separately
        
        int startX = x - r - 1;
        int stopX = x + r;
        int startY = y - r - 1;
        int stopY = y + r;
        
        if (stopX > (w - 1)) {
            stopX = w - 1;
        }
        if (stopY > (h - 1)) {
            stopY = h - 1;
        }
        
        //System.out.println("x=" + x + " y=" + y + " r=" + r
        //    + " startX=" + startX +
        //    " stopX=" + stopX + " startY=" + startY + " stopY=" + stopY);
       
        // when r == 0, bounds need another edit or immediate return
        /*
         2            2           2           2           2       *
         1            1 *         1           1    *      1
         0 *          0           0    *      0           0
           0  1  2      0  1  2     0  1  2     0  1  2     0  1  2
        */
        if (r == 0) {
            if (stopX == 0) {
                if (stopY == 0) {
                    output[1] = 1;
                    output[0] = imgS.getValue(stopX, stopY);
                    return;
                }
                startY = stopY - 1;
                if (startY == 0) {
                    output[1] = 1;
                    output[0] = imgS.getValue(stopX, stopY) 
                        - imgS.getValue(stopX, startY);
                    return;
                }
            } else {
                // stopX > 0
                startX = stopX - 1;
                if (stopY == 0) {
                    output[1] = 1;
                    output[0] = imgS.getValue(stopX, stopY) 
                        - imgS.getValue(startX, stopY);
                    return;
                }
                startY = stopY - 1;
            }
            //System.out.println(" --> startX=" + startX +
            //    " stopX=" + stopX + " startY=" + startY + " stopY=" + stopY);            
        }
        
        if (startX >= 0 && startY >= 0) {
            int nPix = (r == 0) ? 1 : (stopX - startX) * (stopY - startY);
            int s1 = imgS.getValue(stopX, stopY) - imgS.getValue(startX, stopY)
                - imgS.getValue(stopX, startY) + imgS.getValue(startX, startY);
            output[0] = s1;
            output[1] = nPix;
            return;
        } else if (startX >= 0) {
            // startY is < 0
            int nPix = (r == 0) ? 1 : (stopX - startX) * (stopY + 1);
            int s1 = imgS.getValue(stopX, stopY) - imgS.getValue(startX, stopY);
            output[0] = s1;
            output[1] = nPix;
            return;
        } else if (startY >= 0) {
            // startX < 0
            int nPix = (r == 0) ? 1 : (stopX + 1) * (stopY - startY);
            int s1 = imgS.getValue(stopX, stopY)
                - imgS.getValue(stopX, startY);
            output[0] = s1;
            output[1] = nPix;
            return;
        } else {
            // startX < 0 && startY < 0
            int nPix = (r == 0) ? 1 : (stopX + 1) * (stopY + 1);
            int s1 = imgS.getValue(stopX, stopY);
            output[0] = s1;
            output[1] = nPix;
            return;
        }
               
    }
    
    private void applyAbsoluteValue(GreyscaleImage img) {
        int w = img.getWidth();
        int h = img.getHeight();
        for (int x = 0; x < w; ++x) {
            for (int y = 0; y < h; ++y) {
                int v = img.getValue(x, y);
                img.setValue(x, y, Math.abs(v));
            }
        }
    }
}
