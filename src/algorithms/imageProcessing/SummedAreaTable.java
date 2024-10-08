package algorithms.imageProcessing;

import java.util.PriorityQueue;

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
 * Lewis (Fast template matching. 1995 Lewis, J.P.)
 * and then given the name "integral image" and prominently 
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
    
    public double[][] createAbsoluteSummedAreaTable(double[][] img) {

        int w = img.length;
        int h = img[0].length;
        
        ImageProcessor imp = new ImageProcessor();

        double[][] out = imp.copy(img);
        //applyAbsoluteValue
        for (int x = 0; x < w; ++x) {
            for (int y = 0; y < h; ++y) {
                double v = out[x][y];
                if (v < 0) {
                    out[x][y] *= -1;
                }
            }
        }
        
        for (int x = 0; x < w; ++x) {
            for (int y = 0; y < h; ++y) {
                if (x > 0 && y > 0) {
                    double v = out[x - 1][y] + out[x][y - 1] - out[x - 1][y - 1];
                    out[x][y] += v;
                } else if (x > 0) {
                    double v = out[x - 1][y];
                    out[x][y] += v;
                } else if (y > 0) {
                    double v = out[x][y - 1];
                    out[x][y] += v;
                }
            }
        }

        return out;
    }
    
     /**
     * @param imgS
     * @param d diameter of window in x and y
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
                int nPix = sumAndN[1];
                int v = sumAndN[0]/nPix;
                img2.setValue(x, y, v);
            }
        }

        return img2;
    }
    
    /**
     * @param imgS
     * @param d diameter of window in x and y
     * @return 
     */
    public double[][] applyMeanOfWindowFromSummedAreaTable(double[][] imgS, 
        int d) {
        
        int w = imgS.length;
        int h = imgS[0].length;
        
        double[][] img2 = new double[w][];
        for (int i = 0; i < w; ++i) {
            img2[i] = new double[h];
        }
            
        double[] sumAndN = new double[2];
        
        // extract the summed area of each dxd window centered on x,y
        // and divide by number of pixels
        for (int x = 0; x < w; ++x) {
            for (int y = 0; y < h; ++y) {
                extractWindowFromSummedAreaTable(imgS, x, y, d, sumAndN);
                double nPix = sumAndN[1];
                double v = sumAndN[0]/nPix;
                img2[x][y] = v;
            }
        }

        return img2;
    }
    
    /**
     * extract the sum of a window bound by given start and stop coordinates
     * for x and y and return that value and the number of pixels in the
     * window in the output variable, output.
     * NOTE GreyscaleImage, x, and y are in column major format
     * @param imgS
     * @param startX coordinate for x start of window
     * @param stopX coordinate for x stop of window
     * @param startY coordinate for y start of window
     * @param stopY coordinate for y stop of window
     * @param output one dimensional array of size 2 in which the
     * sum of the window will be returned and the number of pixels in the 
     * window.  int[]{sum, nPixels}
     */
    public void extractWindowFromSummedAreaTable(GreyscaleImage imgS, 
        int startX, int stopX, int startY, int stopY, int output[]) {
        
        int w = imgS.getWidth();
        int h = imgS.getHeight();
        int v;
                
        if (startX == 0 && startY == 0) {
            if (stopX == startX && stopY == startY) {
                output[1] = 1;
                v = imgS.getValue(stopX, stopY);
                output[0] = v;
            } else if (stopX > startX && stopY > startY) {
                output[1] = (stopX + 1) * (stopY + 1);
                v = imgS.getValue(stopX, stopY);
                output[0] = v;
            } else if (stopX > startX) {
                //startY==0 && stopY=0
                output[1] = (stopX + 1);
                v = imgS.getValue(stopX, stopY);
                output[0] = v;
            } else if (stopY > startY) {
                output[1] = (stopY + 1);
                v = imgS.getValue(stopX, stopY);
                output[0] = v;
            }
        } else if (startX > 0 && startY > 0) {
            output[1] = ((stopX - startX) + 1) * ((stopY - startY) + 1);

            v = imgS.getValue(stopX, stopY)
                - imgS.getValue(startX - 1, stopY)
                - imgS.getValue(stopX, startY - 1)
                + imgS.getValue(startX - 1, startY - 1);
                
            output[0] = v;
            
        } else if (startX > 0) {
            //startY == 0
            if (stopX == startX && stopY == startY) {
                output[1] = 1;
                v = imgS.getValue(stopX, stopY);
                
                output[0] = v;
            } else {
                output[1] = ((stopX - startX) + 1) * ((stopY - startY) + 1);

                v = imgS.getValue(stopX, stopY)
                    - imgS.getValue(startX - 1, stopY);
                
                output[0] = v;
            }       
        } else if (startY > 0) {
            //startX == 0
            if (stopX == startX && stopY == startY) {
                output[1] = 1;
                v = imgS.getValue(stopX, stopY);
                
                output[0] = v;
            } else {
                output[1] = ((stopX - startX) + 1) * ((stopY - startY) + 1);

                v = imgS.getValue(stopX, stopY)
                    - imgS.getValue(stopX, startY - 1);
                
                output[0] = v;
            }   
        }
    }
    
    /**
     * extract the sum of a window bound by given start and stop coordinates
     * for x and y and return that value and the number of pixels in the
     * window in the output variable, output.
     * NOTE GreyscaleImage, x, and y are in column major format
     * @param imgS
     * @param startX coordinate for x start of window
     * @param stopX coordinate for x stop of window
     * @param startY coordinate for y start of window
     * @param stopY coordinate for y stop of window
     * @param output one dimensional array of size 2 in which the
     * sum of the window will be returned and the number of pixels in the 
     * window.  int[]{sum, nPixels}
     */
    public void extractWindowFromSummedAreaTable(double[][] imgS, 
        int startX, int stopX, int startY, int stopY, double[] output) {
      
        if (startX == 0 && startY == 0) {
            if (stopX == startX && stopY == startY) {
                output[1] = 1;
                output[0] = imgS[stopX][stopY];
            } else if (stopX > startX && stopY > startY) {
                output[1] = (stopX + 1) * (stopY + 1);
                output[0] = imgS[stopX][stopY];
            } else if (stopX > startX) {
                //startY==0 && stopY=0
                output[1] = (stopX + 1);
                output[0] = imgS[stopX][stopY];
            } else if (stopY > startY) {
                output[1] = (stopY + 1);
                output[0] = imgS[stopX][stopY];
            }
        } else if (startX > 0 && startY > 0) {
            output[1] = ((stopX - startX) + 1) * ((stopY - startY) + 1);

            output[0] = imgS[stopX][stopY]
                - imgS[startX - 1][stopY]
                - imgS[stopX][startY - 1]
                + imgS[startX - 1][startY - 1];
                
        } else if (startX > 0) {
            //startY == 0
            if (stopX == startX && stopY == startY) {
                output[1] = 1;
                output[0] = imgS[stopX][stopY];
            } else {
                output[1] = ((stopX - startX) + 1) * ((stopY - startY) + 1);

                output[0] = imgS[stopX][stopY]
                    - imgS[startX - 1][stopY];
            }       
        } else if (startY > 0) {
            //startX == 0
            if (stopX == startX && stopY == startY) {
                output[1] = 1;
                output[0] = imgS[stopX][stopY];
            } else {
                output[1] = ((stopX - startX) + 1) * ((stopY - startY) + 1);

                output[0] = imgS[stopX][stopY]
                    - imgS[stopX][startY - 1];
            }   
        }
    }
    
    /**
     * extract the sum of a window bound by given start and stop coordinates
     * for x and y and return that value and the number of pixels in the
     * window in the output variable, output.
     * NOTE GreyscaleImage, x, and y are in column major format
     * @param imgS
     * @param startX coordinate for x start of window
     * @param stopX coordinate for x stop of window
     * @param startY coordinate for y start of window
     * @param stopY coordinate for y stop of window
     * @param output one dimensional array of size 2 in which the
     * sum of the window will be returned and the number of pixels in the 
     * window.  float[]{sum, nPixels}
     */
    public void extractWindowFromSummedAreaTable(float[][] imgS, 
        int startX, int stopX, int startY, int stopY, float output[]) {
        
        int w = imgS.length;
        int h = imgS[0].length;
        
        if (startX == 0 && startY == 0) {
            if (stopX == startX && stopY == startY) {
                output[1] = 1;
                output[0] = imgS[stopX][stopY];
            } else if (stopX > startX && stopY > startY) {
                output[1] = (stopX + 1) * (stopY + 1);
                output[0] = imgS[stopX][stopY];
            } else if (stopX > startX) {
                //startY==0 && stopY=0
                output[1] = (stopX + 1);
                output[0] = imgS[stopX][stopY];
            } else if (stopY > startY) {
                output[1] = (stopY + 1);
                output[0] = imgS[stopX][stopY];
            }
        } else if (startX > 0 && startY > 0) {
            output[1] = ((stopX - startX) + 1) * ((stopY - startY) + 1);

            output[0] = imgS[stopX][stopY]
                - imgS[startX - 1][stopY]
                - imgS[stopX][startY - 1]
                + imgS[startX - 1][startY - 1];
                
        } else if (startX > 0) {
            //startY == 0
            if (stopX == startX && stopY == startY) {
                output[1] = 1;
                output[0] = imgS[stopX][stopY];
            } else {
                output[1] = ((stopX - startX) + 1) * ((stopY - startY) + 1);

                output[0] = imgS[stopX][stopY]
                    - imgS[startX - 1][stopY];
            }       
        } else if (startY > 0) {
            //startX == 0
            if (stopX == startX && stopY == startY) {
                output[1] = 1;
                output[0] = imgS[stopX][stopY];
            } else {
                output[1] = ((stopX - startX) + 1) * ((stopY - startY) + 1);

                output[0] = imgS[stopX][stopY]
                    - imgS[stopX][startY - 1];
            }   
        }
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
        
        int startX = x - r;
        int stopX = x + r;
        int startY = y - r;
        int stopY = y + r;
        if ((r & 1) == 0) {
            startX--;
            startY--;
        }
        if (startX < 0) {
            startX = 0;
        }
        if (startY < 0) {
            startY = 0;
        }
        if (startX >= w) {
            startX = w - 1;
        }
        if (startY >= h) {
            startY = h - 1;
        }
        if (stopX < 0) {
            stopX = 0;
        }
        if (stopY < 0) {
            stopY = 0;
        }
        if (stopX >= w) {
            stopX = w - 1;
        }
        if (stopY >= h) {
            stopY = h - 1;
        }
              
        extractWindowFromSummedAreaTable(imgS, startX, stopX, startY,  
            stopY, output);
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
    public void extractWindowFromSummedAreaTable(double[][] imgS, 
        int x, int y, int d, double output[]) {
        
        if (output == null || output.length != 2) {
            throw new IllegalArgumentException(
                "output must be initialized to size 2");
        }
        
        if (d < 0) {
            throw new IllegalArgumentException(
                "d must be a non-negative number");
        }
        
        int w = imgS.length;
        int h = imgS[0].length;
        
        if (x < 0 || y < 0 || (x > (w - 1)) || (y > (h - 1))) {
            throw new IllegalArgumentException("x or y is out of bounds of "
                + "image. x=" + x + " y=" + y + " w=" + w + " h=" + h);
        }
        
        final int r = (d >> 1);
        
        int startX = x - r;
        int stopX = x + r;
        int startY = y - r;
        int stopY = y + r;
        
        if (startX < 0) {
            startX = 0;
        }
        if (startY < 0) {
            startY = 0;
        }
        if (startX >= w) {
            startX = w - 1;
        }
        if (startY >= h) {
            startY = h - 1;
        }
        if (stopX < 0) {
            stopX = 0;
        }
        if (stopY < 0) {
            stopY = 0;
        }
        if (stopX >= w) {
            stopX = w - 1;
        }
        if (stopY >= h) {
            stopY = h - 1;
        }
              
        extractWindowFromSummedAreaTable(imgS, startX, stopX, startY,  
            stopY, output);
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
