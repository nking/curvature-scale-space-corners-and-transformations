package algorithms.imageProcessing.features;

import algorithms.MultiArrayMergeSort;
import algorithms.QuickSort;
import algorithms.imageProcessing.FixedSizeSortedVector;
import algorithms.imageProcessing.Gaussian1D;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.misc.MiscMath;
import algorithms.misc.StatsInSlidingWindow;
import algorithms.util.PairInt;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.TFloatList;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * An implementation of "ORB: an efficient alternative to SIFT or SURF"
 * (paper reference is Rublee, Rabaud, Konolige, and Bradski, 2011).
 * http://www.vision.cs.chubu.ac.jp/CV-R/pdf/Rublee_iccv2011.pdf
 * 
 * The implementation below is adapted from the scipy implementation which has
 * the following copyright:
   
 https://github.com/scikit-image/scikit-image/blob/master/LICENSE.txt

-- begin scipy, skimage copyright ---  
Unless otherwise specified by LICENSE.txt files in individual
directories, all code is

Copyright (C) 2011, the scikit-image team
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

 1. Redistributions of source code must retain the above copyright
    notice, this list of conditions and the following disclaimer.
 2. Redistributions in binary form must reproduce the above copyright
    notice, this list of conditions and the following disclaimer in
    the documentation and/or other materials provided with the
    distribution.
 3. Neither the name of skimage nor the names of its contributors may be
    used to endorse or promote products derived from this software without
    specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
-- end scipy, skimage copyright ---  

Note, in some places, scipy functions have been
replaced with existing functions in this project in this implementation below.

 Oriented FAST and rotated BRIEF feature detector and binary descriptor
    extractor.
   
NOTE, have chosen to keep the results (which are available publicly via getters)
in Lists in which each list index is a scale in the pyramidal decimation 
to allow the user to recover scale information if wanted.
To use the way that scipy does, all results as single lists of variables 
are available via getAll*getters.

   
NOT READY FOR USE.  NOT YET TESTED.
 <pre>
 Example Use:
    int nKeyPoints = 200;
    ORB orb = new ORB(nKeyPoints);
    orb.detectAndExtract(image);
    List Descriptors descList = orb.getDescriptors();
       or
    Descriptors desc = orb.getAllDescriptors();
    
    int[][] matches = ORB.matchDescriptors(desc1.descriptors, desc2.descriptors);
 </pre>
 */
public class ORB {
     
    // these two could be made static across all instances
    private int[][] OFAST_MASK = null;
    private int[] OFAST_UMAX = null;
    
    private float downscale = 1.2f;
    private int nScales = 8;
    private int fastN = 9;
    private float fastThreshold = 0.08f;
    private float harrisK = 0.04f;
    private final int nKeypoints;
    
    private int[][] POS0 = null;
    private int[][] POS1 = null;
     
    private List<TIntList> keypointsList = null;
    private List<TDoubleList> orientationsList = null;
    private List<TFloatList> harrisResponses = null;
    private List<TFloatList> scalesList = null;
        
    //`True`` or ``False`` representing
    //the outcome of the intensity comparison for i-th keypoint on j-th
    //decision pixel-pair. It is ``Q == np.sum(mask)``.
    // NOTE: this output format may need to be changed
    private List<Descriptors> descriptorsList = null;
    
    public ORB(int nKeypoints) {
        
        initMasks();
        
        this.nKeypoints = 200;
    }
    
    private void initMasks() {
        
        OFAST_MASK = new int[31][31];
        for (int i = 0; i < 31; ++i) {
            OFAST_MASK[i] = new int[31];
        }
        
        OFAST_UMAX = new int[]{15, 15, 15, 15, 14, 14, 14, 13, 13, 12, 11, 10, 
            9, 8, 6, 3};
        
        for (int i = -15; i < 16; ++i) {
            int absI = Math.abs(i);
            for (int j = -OFAST_UMAX[absI]; j < (OFAST_UMAX[absI] + 1); ++j) {
                OFAST_MASK[15 + j][15 + i] = 1;
            }
        }
    }
      
    /**
     * NOT READY FOR USE.  NOT YET TESTED.
     * 
     * Detect oriented FAST keypoints and extract rBRIEF descriptors.
        Note that this is faster than first calling `detect` and then
        `extract`.
      
      code is adapted from
      https://github.com/scikit-image/scikit-image/blob/401c1fd9c7db4b50ae9c4e0a9f4fd7ef1262ea3c/skimage/feature/orb.py
      with some functions replaced by code in this project. 
      
     * @param image 
     */
    public void detectAndExtract(ImageExt image) {
    
        List<TwoDFloatArray> pyramid = buildPyramid(image);

        keypointsList = new ArrayList<TIntList>();
        orientationsList = new ArrayList<TDoubleList>();
        harrisResponses = new ArrayList<TFloatList>();        
        scalesList = new ArrayList<TFloatList>();
        descriptorsList = new ArrayList<Descriptors>();

        int nKeypointsTotal = 0;
        
        for (int octave = 0; octave < pyramid.size(); ++octave) {
                        
            float[][] octaveImage = pyramid.get(octave).a;
            
            Resp r = detectOctave(octaveImage);
            
            if (r == null) {
                keypointsList.add(new TIntArrayList());
                orientationsList.add(new TDoubleArrayList());
                harrisResponses.add(new TFloatArrayList());
                scalesList.add(new TFloatArrayList());
                descriptorsList.add(new Descriptors());
                continue;
            } else if (r.keypoints.isEmpty()) {
                keypointsList.add(r.keypoints);
                orientationsList.add(r.orientations);
                harrisResponses.add(new TFloatArrayList());
                
                // empty scales and descriptors:
                scalesList.add(new TFloatArrayList());
                descriptorsList.add(new Descriptors());
                continue;
            }
            
            // result contains descriptors and mask.  
            // also, modified by mask are the keypoints and orientations
            Descriptors desc = extractOctave(octaveImage, r.keypoints,
                r.orientations, r.responses);
                    
            // masked keypoints * self.downscale ** octave
            //    multiplies the keypoints by a scalar
            //    to put the coordinates into the reference frame of image
            float scale = (float)Math.pow(this.downscale, octave);
            for (int i = 0; i < keypointsList.size(); ++i) {
                int v = Math.round(scale * r.keypoints.get(i));
                r.keypoints.set(i, v);
            }
            
            keypointsList.add(r.keypoints);
            orientationsList.add(r.orientations);
            harrisResponses.add(r.responses);
            descriptorsList.add(desc);
            
            // scales is same size as keypoints
            TFloatList scales = new TFloatArrayList(r.keypoints.size());
            for (int i = 0; i < r.keypoints.size(); ++i) {
                scales.add(scale);
            }
            scalesList.add(scales);

            nKeypointsTotal += r.keypoints.size();
        }
        
        if (nKeypointsTotal > this.nKeypoints) {
            
            // prune the lists to nKeypoints by the harris corner response as a score
            
            // once thru to count first
            int n = 0;
            for (int idx1 = 0; idx1 < harrisResponses.size(); ++idx1) {
                TFloatList hResp = harrisResponses.get(idx1);
                n += hResp.size();
            }
            
            float[] scores = new float[n];
            PairInt[] indexes = new PairInt[n];
            // outer list index = idx1, inner list index = idx2
            int count = 0;
            for (int idx1 = 0; idx1 < harrisResponses.size(); ++idx1) {
                TFloatList hResp = harrisResponses.get(idx1);
                for (int idx2 = 0; idx2 < hResp.size(); ++idx2) {
                    indexes[count] = new PairInt(idx1, idx2);
                    scores[count] = hResp.get(idx2);
                    count++;
                }
            }
            
            QuickSort.sortBy1stArg(scores, indexes);
            
            TIntObjectMap<TIntList> keep = new TIntObjectHashMap<TIntList>();
            
            // visit largest scores to smallest
            count = 0;
            int idx = n-1;
            while (count < nKeypoints) {
                
                PairInt m = indexes[idx];
                
                TIntList list = keep.get(m.getX());
                if (list == null) {
                    list = new TIntArrayList();
                    keep.put(m.getX(), list);
                }
                list.add(m.getY());
                
                idx--;
                count++;
            }
            
            List<TIntList> keypointsList2 = new ArrayList<TIntList>();
            List<TDoubleList> orientationsList2 = new ArrayList<TDoubleList>();
            List<TFloatList> harrisResponses2 = new ArrayList<TFloatList>();        
            List<TFloatList> scalesList2 = new ArrayList<TFloatList>();
            List<Descriptors> descriptorsList2 = new ArrayList<Descriptors>();
            
            for (int idx1 = 0; idx1 < keypointsList.size(); ++idx1) {
                TIntList keepList = keep.get(idx1);
                if (keepList == null) {
                    continue;
                }
                int n2 = keepList.size();
                TIntList kp = new TIntArrayList(n2 * 2);
                TDoubleList or = new TDoubleArrayList(n2);
                TFloatList hr = new TFloatArrayList(n2);
                TFloatList s = new TFloatArrayList(n2);
                TIntList m = new TIntArrayList(n2);
                //holds values 1 or 0.  size is [orientations.size][POS0.length]
                int[][] d = new int[n2][POS0.length];
                int dCount = 0;
                for (int i = 0; i < keepList.size(); ++i) {
                    int idx2 = keepList.get(i);
                    
                    kp.add(this.keypointsList.get(idx1).get(2*idx2));
                    kp.add(this.keypointsList.get(idx1).get(2*idx2 + 1));
                    
                    or.add(this.orientationsList.get(idx1).get(idx2));
                    hr.add(this.harrisResponses.get(idx1).get(idx2));
                    s.add(this.scalesList.get(idx1).get(idx2));
                    m.add(this.descriptorsList.get(idx1).mask.get(idx2));
                    
                    int[] d0 = this.descriptorsList.get(idx1).descriptors[idx2];
                    d[dCount] = Arrays.copyOf(d0, d0.length);
                    dCount++;
                }
                
                keypointsList2.add(kp);
                orientationsList2.add(or);
                harrisResponses2.add(hr);
                scalesList2.add(s);
                Descriptors desc = new Descriptors();
                desc.mask = m;
                desc.descriptors = d;
                descriptorsList2.add(desc);
            }
            
            keypointsList = keypointsList2;
            orientationsList = orientationsList2;
            harrisResponses = harrisResponses2;
            scalesList = scalesList2;
            descriptorsList = descriptorsList2;
        }
    }
    
    /**
     * https://github.com/scikit-image/scikit-image/blob/master/skimage/feature/orb.py
     * @param image
     * @return 
     */
    private List<TwoDFloatArray> buildPyramid(ImageExt image) {
        
        /*
        NOTE TO SELF: 
            can probably replace this with 
            MedianTransform.multiscalePyramidalMedianTransform2
        then normalize to floats
        */
        
        float[][] gsImgF = prepareGrayscaleInput2D(image);
        
        return pyramidGaussian(gsImgF, nScales - 1, downscale);
    }
    
    /**
     * create a greyscale image from the color image
     * and return it as two dimensional float array, scaled to
     * range [0.0, 1.0] and placed in format img[row][col].
     * @param image 
     */
    private float[][] prepareGrayscaleInput2D(ImageExt image) {
        
        GreyscaleImage img = image.copyToGreyscale();
        
        int nRows = img.getHeight();
        int nCols = img.getWidth();
        float[][] out = new float[nRows][nCols];
        for (int j = 0; j < nRows; ++j) {
            out[j] = new float[nCols];
        }
        
        for (int j = 0; j < nRows; ++j) {
            for (int i = 0; i < nCols; ++i) {
                out[j][i] = (float)img.getValue(i, j)/255.f;
            }
        }
        
        return out;
    }

    protected void debugPrint(String label, float[][] a) {
        
        StringBuilder sb = new StringBuilder(label);
        sb.append("\n");
        for (int i = 0; i < a.length; ++i) {
            sb.append(Arrays.toString(a[i])).append("\n");
        }
        
        System.out.println(sb.toString());
    }
    
    private void debugPrint(String label, int[][] a) {
        
        StringBuilder sb = new StringBuilder(label);
        sb.append("\n");
        for (int i = 0; i < a.length; ++i) {
            sb.append(Arrays.toString(a[i])).append("\n");
        }
        
        System.out.println(sb.toString());
    }

    private void applyShift(float[][] imageMax, int minDistance, int nRows,
        int nCols) {
        
        for (int i = 0; i < nRows; ++i) {
            System.arraycopy(imageMax[i], 0, imageMax[i], minDistance, 
                nCols - minDistance);
            for (int j = 0; j < minDistance; ++j) {
                imageMax[i][j] = 0;
            }
            for (int j = (nCols - minDistance); j < nCols; ++j) {
                imageMax[i][j] = 0;
            }
        }
        //debugPrint("shift 0 imageMax=", imageMax);
        for (int j = 0; j < nCols; ++j) {
            for (int i = (nRows - minDistance); i >= minDistance; --i) {
                imageMax[i][j] = imageMax[i - minDistance][j];
            }
            for (int i = 0; i < minDistance; ++i) {
                imageMax[i][j] = 0;
            }
            for (int i = (nRows - minDistance); i < nRows; ++i) {
                imageMax[i][j] = 0;
            }
        }
        //debugPrint("shifted imageMax=", imageMax);        
    }

    private float[][] multiply(float[][] a, float[][] b) {
        
        float[][] c = copy(a);
        
        for (int i = 0; i < c.length; ++i) {
            for (int j = 0; j < c[0].length; ++j) {
                c[i][j] *= b[i][j];
            }
        }
        
        return c;
    }
    
    private float[][] add(float[][] a, float[][] b) {
        
        float[][] c = copy(a);
        
        for (int i = 0; i < c.length; ++i) {
            for (int j = 0; j < c[0].length; ++j) {
                c[i][j] += b[i][j];
            }
        }
        
        return c;
    }
    
    private float[][] subtract(float[][] a, float[][] b) {
        
        float[][] c = copy(a);
        
        for (int i = 0; i < c.length; ++i) {
            for (int j = 0; j < c[0].length; ++j) {
                c[i][j] -= b[i][j];
            }
        }
        
        return c;
    }

    protected static class TwoDFloatArray {
        float[][] a;
        public TwoDFloatArray(float[][] b) {
            a = b;
        }
    }
    
    private float[][] copy(float[][] a) {
        
        int n1 = a.length;
        int n2 = a[0].length;
        
        float[][] out = new float[n1][n2];
        for (int i = 0; i < n1; ++i) {
            out[i] = Arrays.copyOf(a[i], a[i].length);
        }
        
        return out;
    }
    
    /**
     Yield images of the Gaussian pyramid formed by the input image.
     Recursively applies the `pyramid_reduce` function to the image, and yields
     the downscaled images.
     Note that the first image of the pyramid will be the original, unscaled
     image. The total number of images is `max_layer + 1`. In case all layers
     are computed, the last image is either a one-pixel image or the image where
     the reduction does not change its shape.
     
     https://github.com/scikit-image/scikit-image/blob/master/skimage/transform/pyramids.py
     
     @param gsImgF 
     */
    private List<TwoDFloatArray> pyramidGaussian(float[][] gsImgF, int maxLayer, 
        float downscale) {
        
        /*
        -------------------------------------
        sigma = None, 
            Default is `2 * downscale / 6.0
        order = 1, 
        mode = 'reflect', 
            How to handle the borders
        cval = 0
        -------------------------------------
        */
        
        if (downscale <= 1) {
            throw new IllegalArgumentException("downscale must be greater than 1");
        }
        
        // if downscale=default=1.2, sigma=0.4.
        // NOTE: should probably consider using SIGMA.ZEROPOINTFIVE or SIGMA.ZEROPOINTSEVENONE
        //   to use binomial
        
        float sigma = 2.f * downscale / 6.0f;
        int order = 1;
        int mode = 0;//reflect = 0
        int cVal = 0;
        
        int layer = 0;
        int nRows = gsImgF.length;
        int nCols = gsImgF[0].length;
        
        List<TwoDFloatArray> pyramidList = new ArrayList<TwoDFloatArray>();
        
        float[][] prevLayerImg = copy(gsImgF);
        int prevNRows = nRows;
        int prevNCols = nCols;
                
        // build downsampled images until maxLayer is reached or downscale
        // process does not change image size
        while (layer != maxLayer) {
            
            pyramidList.add(new TwoDFloatArray(prevLayerImg));
            
            layer++;
            
            float[][] layerImg = pyramidReduce(prevLayerImg, downscale, 
                sigma, order);
            
            prevNRows = nRows;
            prevNCols = nCols;
            prevLayerImg = layerImg;
            nRows = prevLayerImg.length;
            nCols = prevLayerImg[0].length;
            
            if (prevNCols == nCols && prevNRows == nRows) {
                break;
            }
        }
        
        return pyramidList;
    }
    
    /**
     * Smooth and then downsample image.
     * (NOTE, should replace the invoker of this and this method with
     * MedianTransform.multiscalePyramidalMedianTransform2())
     * 
     * @param img
     * @param downscale
     * @param sigma - Sigma for Gaussian filter. 
        Default is `2 * downscale / 6.0` which
        corresponds to a filter mask twice the size of the scale factor that
        covers more than 99% of the Gaussian distribution.
     * @param order
     *  Order of splines used in interpolation of downsampling
     * @return 
     */
    private float[][] pyramidReduce(float[][] img, 
        float downscale, float sigma, int order) {
        
        if (downscale <= 1) {
            throw new IllegalArgumentException("downscale must be greater than 1");
        }
        
        int nRows = img.length;
        int nCols = img[0].length;
        
        // downscale is > 1
        int outRows = (int)Math.ceil((float)nRows / downscale);
        int outCols = (int)Math.ceil((float)nCols / downscale);
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        float[][] smoothed = copy(img);
        
        float[] kernel = Gaussian1D.getKernel(sigma);
        
        imageProcessor.applyKernelTwo1Ds(smoothed, kernel);
        
        int totalSize = outRows * outCols;
        
        float[][] out = new float[outRows][outCols];
        for (int i = 0; i < outRows; ++i) {
            out[i] = new float[outCols];
        }
        
        if (totalSize == 0) {
            return out;
        }
        
        // unravel input smoothed into row-major, C-style ordered array
        // then reshape into 2-d out size array
        int count = 0;
        float[][] b = new float[outRows][outCols];
        for (int i = 0; i < nRows; ++i) {
            b[i] = new float[outCols];
            for (int j = 0; j < nCols; ++j) {
                b[i][j] = smoothed[i][j];
                count++;
                if (count == totalSize) {
                    break;
                }
            }
            if (count == totalSize) {
                break;
            }
        }
       
        return b;
    }
    
    private class Resp {
        TDoubleList orientations;
        TIntList keypoints;
        TFloatList responses;
    }
    
    public static class Descriptors {
        
        //mask of length orientations.size() containing a 1 or 0
        // indicating if pixels are within the image (``True``) or in the
        // border region of the image (``False``).
        TIntList mask;
        
        //holds values 1 or 0.  size is [orientations.size][POS0.length]
        int[][] descriptors;
    }
     
    /**
     * adapted from 
     * https://github.com/scikit-image/scikit-image/blob/master/skimage/feature/orb.py
     * 
     * @param octaveImage
     * @return 
     */
    private Resp detectOctave(float[][] octaveImage) {
        
        float[][] fastResponse = cornerFast(octaveImage, fastN, fastThreshold);
    
        int nRows = fastResponse.length;
        int nCols = fastResponse[0].length;
        
        // list of format [row, col, ...] of filtered maxima ordered by intensity
        TIntList coords = cornerPeaks(fastResponse, 1);
        if (coords.isEmpty()) {
            return null;
        }
        
        TIntList keypoints2 = maskCoordinates(coords, nRows, nCols, 16);
    
        // size is keyPoints2.size/2
        double[] orientations2 = cornerOrientations(octaveImage, keypoints2);
        assert(orientations2.length == keypoints2.size()/2);
        
        // size is same a octaveImage
        float[][] harrisResponse = cornerHarris(octaveImage);

        TFloatList responses = new TFloatArrayList(orientations2.length);
        for (int i = 0; i < keypoints2.size(); i += 2) {
            int x = keypoints2.get(i);
            int y = keypoints2.get(i + 1);
            responses.add(harrisResponse[x][y]);
        }
        
        Resp r2 = new Resp();
        r2.keypoints = keypoints2;
        r2.responses = responses;
        r2.orientations = new TDoubleArrayList();
        for (int i = 0; i < orientations2.length; ++i) {
            r2.orientations.add(orientations2[i]);
        }
        
        return r2;
    }
    
    /**
     Extract FAST corners for a given image.
     
     code is adapted from :
     https://github.com/scikit-image/scikit-image/blob/master/skimage/feature/corner.py
     
     which has the copyright above in the class level documentation.
     
    <pre>
     References
    ----------
    .. [1] Edward Rosten and Tom Drummond
           "Machine Learning for high-speed corner detection",
           http://www.edwardrosten.com/work/rosten_2006_machine.pdf
    .. [2] Wikipedia, "Features from accelerated segment test",
           https://en.wikipedia.org/wiki/Features_from_accelerated_segment_test 
     </pre>
     
     @param img
     *    Input image as 2D array in row-major, C-style ordered array
     @param n
     *    Minimum number of consecutive pixels out of 16 pixels on the circle
          that should all be either brighter or darker w.r.t testpixel.
          A point c on the circle is darker w.r.t test pixel p if
          `Ic lessThan Ip - threshold` and brighter if `Ic greaterThan Ip + threshold`. 
          Also stands for the n in `FAST-n` corner detector.
          (the scipy default is n=12).
     * @param threshold 
     *    Threshold used in deciding whether the pixels on the circle are
          brighter, darker or similar w.r.t. the test pixel. Decrease the
          threshold when more corners are desired and vice-versa.
          (the scipy default is threshold=0.15)
     @return
        FAST corner response image.
     */
    protected float[][] cornerFast(final float[][] img, final int n, final float threshold) {
        
        // NOTE: below, removed code for n > 12 because internal code invocation
        // is for fastN which is '9'
        assert(n == fastN);
        
        int nRows = img.length;
        int nCols = img[0].length;

        int i, j, k;

        int speed_sum_b, speed_sum_d;
        float curr_pixel;
        float lower_threshold, upper_threshold;
        
        float[][] cornerResponse = new float[nRows][nCols];
        for (i = 0; i < nRows; ++i) {
            cornerResponse[i] = new float[nCols];
        }
        
        int[] rp = new int[]{0, 1, 2, 3, 3, 3, 2, 1, 0, -1, -2, -3, -3, -3, -2, -1};
        int[] cp = new int[]{3, 3, 2, 1, 0, -1, -2, -3, -3, -3, -2, -1, 0, 1, 2, 3};
        char[] bins = new char[16];
        float[] circleIntensities = new float[16];
        float currResponse;
        
        for (i = 3; i < (nRows - 3); ++i) {
            for (j = 3; j < (nCols - 3); ++j) {
                curr_pixel = img[i][j];
                lower_threshold = curr_pixel - threshold;
                upper_threshold = curr_pixel + threshold;

                for (k = 0; k < 16; ++k) {
                    circleIntensities[k] = img[i + rp[k]][j + cp[k]];
                    if (circleIntensities[k] > upper_threshold) {
                        //# Brighter pixel
                        bins[k] = 'b';
                    } else if (circleIntensities[k] < lower_threshold) {
                        //# Darker pixel
                        bins[k] = 'd';
                    } else {
                        //# Similar pixel
                        bins[k] = 's';
                    }
                }
                //# High speed test for n >= 12 removed because using n = fastN = 8
                

                //# Test for bright pixels
                currResponse = _corner_fast_response(
                    curr_pixel, circleIntensities, bins, 'b', n);

                //# Test for dark pixels
                if (currResponse == 0) {
                    currResponse = _corner_fast_response(
                        curr_pixel, circleIntensities, bins, 'd', n);
                }

                cornerResponse[i][j] = currResponse;
            }
        }

        return cornerResponse;
    }
   
    /**
     * https://github.com/scikit-image/scikit-image/blob/master/skimage/feature/corner_cy.pyx
   
     * @param curr_pixel
     * @param circleIntensities
     * @param bins
     * @param state
     * @param n
     * @return 
     */
    private float _corner_fast_response(float curr_pixel, 
        float[] circleIntensities, char[] bins, char state, final int n) {
    
        int consecutiveCount = 0;
        float currResponse;
        int l, m;
        for (l = 0; l < (15 + n); ++l) {
            if (bins[l % 16] == state) {
                consecutiveCount += 1;
                if (consecutiveCount == n) {
                    currResponse = 0;
                    for (m = 0; m < 16; ++m) {
                        currResponse += Math.abs(circleIntensities[m] - curr_pixel);
                    }
                    return currResponse;
                }
            } else {
                consecutiveCount = 0;
            }
        }
        
        return 0;
    }
    
    private TIntList maskCoordinates(TIntList orderedPeaks, 
        int nRows, int nCols, int maskRadius) {
       
        TIntSet peaks = new TIntHashSet();
        
        /* making single pixel index out of coordinates:
        (row * width) + col
        pixIdxs[count] = (j * nRows) + i;
        */
        for (int i = 0; i < orderedPeaks.size(); i += 2) {
            int ii = orderedPeaks.get(i);
            int jj = orderedPeaks.get(i + 1);
            int pixIdx = (jj * nRows) + ii;
            peaks.add(pixIdx);
        }
        
        for (int i = 0; i < orderedPeaks.size(); i += 2) {
            int ii = orderedPeaks.get(i);
            int jj = orderedPeaks.get(i + 1);
            int pixIdx = (jj * nRows) + ii;
            if (!peaks.contains(pixIdx)) {
                continue;
            }
            for (int k0 = ii - maskRadius; k0 <= (ii + maskRadius); ++k0) {
                if ((k0 < 0) || (k0 > (nRows - 1))) {
                    continue;
                }
                for (int k1 = jj - maskRadius; k1 <= (jj + maskRadius); ++k1) {
                    if ((k1 < 0) || (k1 > (nCols - 1)) || (k0 == ii && k1 == jj)) {
                        continue;
                    }
                    int pixIdx2 = (k1 * nRows) + k0;
                    if (peaks.contains(pixIdx2)) {
                        peaks.remove(pixIdx2);
                    }
                }
            }
        }

        TIntList peakRowCols2 = new TIntArrayList(peaks.size() * 2);
        for (int i = 0; i < orderedPeaks.size(); i += 2) {
            int ii = orderedPeaks.get(i);
            int jj = orderedPeaks.get(i + 1);
            int pixIdx = (jj * nRows) + ii;
            if (peaks.contains(pixIdx)) {
                peakRowCols2.add(ii);
                peakRowCols2.add(jj);
            }
        }
        
        return peakRowCols2;
    }
    
    /**
     * 
     * @param input
     * @param nRows
     * @param nCols
     * @param border
     * @return mask of length input.size/2 containing a 1 or 0
     *     indicating if pixels are within the image (``True``) or in the
           border region of the image (``False``).
     */
    private TIntList maskCoordinatesIfBorder(TIntList input, 
        int nRows, int nCols, int border) {
       
        TIntSet set = new TIntHashSet();
        
        /* making single pixel index out of coordinates:
        (row * width) + col
        pixIdxs[count] = (j * nRows) + i;
        */
        for (int i = 0; i < input.size(); i += 2) {
            int ii = input.get(i);
            int jj = input.get(i + 1);
            int pixIdx = (jj * nRows) + ii;
            set.add(pixIdx);
        }
        
        for (int i = 0; i < input.size(); i += 2) {
            int ii = input.get(i);
            int jj = input.get(i + 1);
            int pixIdx = (jj * nRows) + ii;
            if (!set.contains(pixIdx)) {
                continue;
            }
            if ((ii < border) || (jj < border) || (ii > (nRows - border - 1)) ||
                (jj > (nCols - border - 1))) {
                set.remove(pixIdx);
            }
        }
        
        /*
        Mask indicating if pixels are within the image (``True``) or in the
            border region of the image (``False``).
        */

        TIntList mask = new TIntArrayList(input.size()/2);
        for (int i = 0; i < input.size(); i += 2) {
            int ii = input.get(i);
            int jj = input.get(i + 1);
            int pixIdx = (jj * nRows) + ii;
            if (set.contains(pixIdx)) {
                mask.add(1);
            } else {
                mask.add(0);
            }
        }
        
        return mask;
    }
    
    /**
     * https://github.com/scikit-image/scikit-image/blob/d19b60add22b818298c7aefa65f40e7c1467ef4d/skimage/feature/corner.py
     * 
      Find corners in corner measure response image.
      This differs from `skimage.feature.peak_local_max` in that it suppresses
      multiple connected peaks with the same accumulator value.
     * 
     * @param img
     * @param minDistance
     * @return a list of coordinates in format [row, col, ...] 
     * sorted by decreasing pixel intensity.
     * @return 
     */
    protected TIntList cornerPeaks(float[][] img, int minDistance) {

        //threshold_abs=None, 
        float thresholdRel = 0.1f;
        boolean excludeBorder = true; 
        boolean indices = true;
        int numPeaks = Integer.MAX_VALUE;
        //footprint=None, labels=None   
        
        int nRows = img.length;
        int nCols = img[0].length;
        
        // these have been sorted by decreasing intensity        
        TIntList peakRowCols = peakLocalMax(img, minDistance, thresholdRel);
        
        System.out.println("peakRowCols in cornerPeaks=" + peakRowCols.toString()
            + "\nsize=" + peakRowCols.size());
        
        TIntList peakRowCols2 = maskCoordinates(peakRowCols, nRows, nCols, 
            minDistance);
        
        return peakRowCols2;
    }

    /**
     * adapted from
     https://github.com/scikit-image/scikit-image/blob/92a38515ac7222aab5e606f9de46caf5f503a7bd/skimage/feature/peak.py

     Find peaks in an image as coordinate list or boolean mask.
     Peaks are the local maxima in a region of `2 * min_distance + 1`
     (i.e. peaks are separated by at least `min_distance`).
     If peaks are flat (i.e. multiple adjacent pixels have identical
     intensities), the coordinates of all such pixels are returned.
     If both `threshold_abs` and `threshold_rel` are provided, the maximum
     of the two is chosen as the minimum intensity threshold of peaks.
    
     * @param img
     * @param minDistance 
        Minimum number of pixels separating peaks in a region of `2 *
        min_distance + 1` (i.e. peaks are separated by at least
        `min_distance`).
        To find the maximum number of peaks, use `min_distance=1`.
     @return [row, column, ...] coordinates of peaks, sorted by decreasing
     * intensity.
     */
    protected TIntList peakLocalMax(float[][] img, int minDistance,
        float thresholdRel) {
         
        int excludeBorder = minDistance;
        int numPeaks = Integer.MAX_VALUE; 
        //int numPeaksPerLabel = Integer.MAX_VALUE;
        
        /*
        The peak local maximum function returns the coordinates of local peaks
        (maxima) in an image. A maximum filter is used for finding local maxima.
        This operation dilates the original image. After comparison of the dilated
        and original image, this function returns the coordinates or a mask of the
        peaks where the dilated image equals the original image.
        */
        
        int nRows = img.length;
        int nCols = img[0].length;
        
        //# Non maximum filter
        int size = 2 * minDistance + 1;
        float[][] imageMax = maximumFilter(img, size);
        assert(nRows == imageMax.length);
        assert(nCols == imageMax[0].length);
        //mask = image == image_max
        
        //debugPrint("before shift imageMax=", imageMax);
        
        // a fudge to match results of scipy which must store same windows at
        // locations shifted by minDistance or so in x and y from the
        // beginning of the sliding window
        applyShift(imageMax, minDistance, nRows, nCols);
        
        //TODO: should be able to simplify the mask here

        // 1's where same, else 0's
        int[][] mask = new int[nRows][nCols];
        for (int i = 0; i < nRows; ++i) {
            mask[i] = new int[nCols];
            for (int j = 0; j < nCols; ++j) {
                if (img[i][j] == imageMax[i][j]) {
                    mask[i][j] = 1;
                }
            }
        }
        
        //debugPrint("0 mask=", mask);
        
        /*
        // exclude border
        for (int i = 0; i < nRows; ++i) {
            if ((i < excludeBorder) || (i > (nRows - 1 - excludeBorder))){
                Arrays.fill(mask[i], 0);
            } else {
                Arrays.fill(mask[i], 0, excludeBorder, 0);
                Arrays.fill(mask[i], nCols - excludeBorder, nCols, 0);
            }
        }
        */
        
        // find top peak candidates above a threshold.
        // TODO: should this use mask so excluding borders?
        float thresholdAbs = MiscMath.findMin(img);
        float thresholdMax = thresholdRel * MiscMath.findMax(img);        
        thresholdAbs = Math.max(thresholdAbs, thresholdMax);
        
        // mask &= image > 0.1
        for (int i = 0; i < nRows; ++i) {
            for (int j = 0; j < nCols; ++j) {
                if (imageMax[i][j] > thresholdAbs) {
                    mask[i][j] &= 1;
                } else {
                    mask[i][j] = 0;
                }
            }
        }
        
        //debugPrint("mask &= image > " + thresholdAbs, mask);
        
        // Select highest intensities (num_peaks)
        // expected output is [row index, col index, ...]
        
        //TODO: should num_peaks be this.nKeypoints?  re-read paper...
        TIntList pixIndexes = new TIntArrayList();
        if (numPeaks == Integer.MAX_VALUE) {
            // find non-zero pixels in mask
            float[] values = new float[nRows * nCols];
            int[] pixIdxs = new int[values.length];
            int count = 0;
            for (int i = 0; i < mask.length; ++i) {
                for (int j = 0; j < mask[i].length; ++j) {
                    if (mask[i][j] > 0.f) {
                        values[count] = img[i][j];
                        //(row * width) + col
                        pixIdxs[count] = (j * nRows) + i;
                        count++;
                    }
                }
            }
            values = Arrays.copyOf(values, count);
            pixIdxs = Arrays.copyOf(pixIdxs, count);
            MultiArrayMergeSort.sortByDecr(values, pixIdxs);
            for (int i = 0; i < values.length; ++i) {
                int pixIdx = pixIdxs[i];
                int jj = pixIdx/nRows;
                int ii = pixIdx - (jj * nRows);
                pixIndexes.add(ii);
                pixIndexes.add(jj);
            }
        } else {
            //need to sort to keep top numPeaks
            FixedSizeSortedVector<Pix> vec = new
                FixedSizeSortedVector<Pix>(numPeaks, Pix.class);
            for (int i = 0; i < mask.length; ++i) {
                for (int j = 0; j < mask[i].length; ++j) {
                    if (mask[i][j] > 0.f) {
                        Pix pix = new Pix(i, j, Float.valueOf(img[i][j]));
                        vec.add(pix);
                    }
                }
            }
            for (int i = 0; i < vec.getNumberOfItems(); ++i) {
                Pix pix = vec.getArray()[i];
                pixIndexes.add(pix.i);
                pixIndexes.add(pix.j);
            }
        }
        
        return pixIndexes;
    }
    
    private class Pix implements Comparable<Pix> {

        public final int i;
        public final int j;
        public final Float value;
        public Pix(int i, int j, Float v) {
            this.i = i;
            this.j = j;
            this.value = v;
        }
        @Override
        public int compareTo(Pix other) {
            // changed for a descending sort
            return other.value.compareTo(this.value);
        }
        
    }
    
    /**
     * @author nichole
     * @param img
     * @param size
     * @return 
     */
    private float[][] maximumFilter(float[][] img, int size) {
        
        int nRows = img.length;
        int nCols = img[0].length;
        
        // return_value = out
        float[][] out = new float[nRows][nCols];
        for (int i = 0; i < nRows; ++i) {
            out[i] = new float[nCols];
        }
        
        // have adapted median window algorithm for this:
        StatsInSlidingWindow maxWindow = new StatsInSlidingWindow();
        maxWindow.calculateMaximum(img, out, size, size);
        
        return out;
    }

    /**
        Compute the orientation of corners.
        The orientation of corners is computed using the first order central moment
        i.e. the center of mass approach. The corner orientation is the angle of
        the vector from the corner coordinate to the intensity centroid in the
        local neighborhood around the corner calculated using first order central
        moment.
        adapted from
        from https://github.com/scikit-image/scikit-image/blob/master/skimage/feature/corner_cy.pyx
        
        References
        ----------
        .. [1] Ethan Rublee, Vincent Rabaud, Kurt Konolige and Gary Bradski
              "ORB : An efficient alternative to SIFT and SURF"
              http://www.vision.cs.chubu.ac.jp/CV-R/pdf/Rublee_iccv2011.pdf
        .. [2] Paul L. Rosin, "Measuring Corner Properties"
              http://users.cs.cf.ac.uk/Paul.Rosin/corner2.pdf
     * @param octaveImage
           Input grayscale image.            
     * @param corners
     *     Corner coordinates as ``(row, col, ...)``.
     * @localParam OFAST_MASK
     *     Mask defining the local neighborhood of the corner used for the
           calculation of the central moment.
     * @return 
     *    orientations : (N, 1) array
               Orientations of corners in the range [-pi, pi].
     */
    protected double[] cornerOrientations(float[][] octaveImage, 
        TIntList corners) {
        
        //same as mask, same 0's and 1's:
        //cdef unsigned char[:, ::1] cmask = np.ascontiguousarray(mask != 0, dtype=np.uint8)
        
        int i, r, c, r0, c0;
        int mRows = OFAST_MASK.length;
        int mCols = OFAST_MASK[0].length;
        int mRows2 = (mRows - 1) / 2;
        int mCols2 = (mCols - 1) / 2;
        
        //cdef double[:, :] cimage = np.pad(image, (mrows2, mcols2), mode='constant',
        //    constant_values=0)
        int nRows2 = octaveImage.length + (2 * mRows2);
        int nCols2 = octaveImage[0].length + (2 * mCols2);
        double[][] cImage = new double[nRows2][nCols2];
        for (i = 0; i < nRows2; ++i) {
            cImage[i] = new double[nCols2];
            if ((i >= mCols2) && (i < (cImage[i].length - mCols2))) {
                float[] src = octaveImage[i - mCols2];
                double[] dest = cImage[i];
                for (int ii = 0; ii < src.length; ++ii) {
                    dest[ii + mCols2] = src[ii];
                }
            }
        }
        
        // number of corner coord pairs
        int nCorners = corners.size()/2;
        
        double[] orientations = new double[nCorners];
        
        double curr_pixel;
        double m01, m10, m01_tmp;
          
        for (i = 0; i < corners.size(); i += 2) {
            r0 = corners.get(i);
            c0 = corners.get(i + 1);

            m01 = 0;
            m10 = 0;

            for (r = 0; r < mRows; ++r) {
                m01_tmp = 0;
                for (c = 0; c < mCols; ++c) {
                    if (OFAST_MASK[r][c] > 0) {
                        curr_pixel = cImage[r0 + r][c0 + c];
                        m10 += curr_pixel * (c - mCols2);
                        m01_tmp += curr_pixel;
                    }
                }
                m01 += m01_tmp * (r - mRows2);
            }

            //arc tangent of y/x, in the interval [-pi,+pi] radians
            orientations[i/2] = Math.atan2(m01, m10);
        }
        
        return orientations;
    }

    /**
      Compute Harris corner measure response image.
        This corner detector uses information from the auto-correlation matrix A::
            A = [(imx**2)   (imx*imy)] = [Axx Axy]
                [(imx*imy)   (imy**2)]   [Axy Ayy]
        Where imx and imy are first derivatives, averaged with a gaussian filter.
        The corner measure is then defined as::
            det(A) - k * trace(A)**2
      
      adapted from
      from https://github.com/scikit-image/scikit-image/blob/master/skimage/feature/corner.py
     
    @param image
    @return the harris corner response image of same size as image, 
    *   and composed of
    *   response = detA - k * traceA ** 2 built from the 2nd derivatives of
    *   image intensity.
    */    
    protected float[][] cornerHarris(float[][] image) {
        
        // method = 'k'.  k is Sensitivity factor to separate corners from edges, 
        // Small values of k result in detection of sharp corners.
        float k = this.harrisK;
        
        // Standard deviation used for the Gaussian kernel, which is used as
        // weighting function for the auto-correlation matrix.
        float sigma = 1;
        
        TwoDFloatArray[] tensorComponents = structureTensor(image, sigma);
        
        float[][] axxyy = multiply(tensorComponents[0].a, 
            tensorComponents[2].a);
        
        float[][] axyxy = multiply(tensorComponents[1].a, 
            tensorComponents[1].a);
        
        float[][] detA = subtract(axxyy, axyxy);
        
        float[][] traceA = add(tensorComponents[0].a, 
            tensorComponents[2].a);
        
        //response = detA - k * traceA ** 2
        float[][] response = copy(detA);
        for (int i = 0; i < detA.length; ++i) {
            for (int j = 0; j < detA[i].length; ++j) {
                float v = k * (traceA[i][j] * traceA[i][j]);
                response[i][j] -= v;
            }
        }
        
        return response;
    }

    /**
     Compute structure tensor using sum of squared differences.
     The structure tensor A is defined as::
         A = [Axx Axy]
             [Axy Ayy]
     which is approximated by the weighted sum of squared differences in a local
     window around each pixel in the image.
    
     adapted from
     https://github.com/scikit-image/scikit-image/blob/master/skimage/feature/corner.py
     and replaced with existing local project functions.
     
     * @param image
     * @param sigma
     * @return [Axx, Axy, Ayy]
     * Axx : ndarray
          Element of the structure tensor for each pixel in the input image.
       Axy : ndarray
          Element of the structure tensor for each pixel in the input image.
       Ayy : ndarray
          Element of the structure tensor for each pixel in the input image.
     */
    protected TwoDFloatArray[] structureTensor(float[][] image, float sigma) {
   
        // --- create Sobel derivatives ----
        ImageProcessor imageProcessor = new ImageProcessor();
        
        // switch X and Y sobel operations to match scipy
        
        float[][] gX = copy(image);
        imageProcessor.applySobelY(gX);
        
        float[][] gY = copy(image);
        imageProcessor.applySobelX(gY);
        
        //debugPrint("gX", gX);
        //debugPrint("gY", gY);
       
        // --- create structure tensors ----
        float[] kernel = Gaussian1D.getKernel(sigma);
        float[][] axx = multiply(gX, gX);
        imageProcessor.applyKernelTwo1Ds(axx, kernel);
        
        float[][] axy = multiply(gX, gY);
        imageProcessor.applyKernelTwo1Ds(axy, kernel);
        
        float[][] ayy = multiply(gY, gY);
        imageProcessor.applyKernelTwo1Ds(ayy, kernel);

        //TODO: might need to apply a normalization factor
        // to these for downstream use
        
        TwoDFloatArray[] tensorComponents = new TwoDFloatArray[3];
        tensorComponents[0] = new TwoDFloatArray(axx);
        tensorComponents[1] = new TwoDFloatArray(axy);
        tensorComponents[2] = new TwoDFloatArray(ayy);
        
        return tensorComponents;
    }
    
    /**
     * filter the keypoints, orientations, and responses to remove those close 
     * to the border and then create descriptors for the remaining.
     * 
     * @param octaveImage
     * @param keypoints
     * @param orientations
     * @return the encapsulated descriptors and mask
     */
    protected Descriptors extractOctave(float[][] octaveImage,
        TIntList keypoints, TDoubleList orientations, TFloatList responses) {
        
        assert(orientations.size() == keypoints.size()/2);
        
        int nRows = octaveImage.length;
        int nCols = octaveImage[0].length;
        
        //mask of length orientations.size() containing a 1 or 0
        // indicating if pixels are within the image (``True``) or in the
        // border region of the image (``False``).
        TIntList mask = maskCoordinatesIfBorder(keypoints, nRows, nCols, 16);
        
        /*
        i    0    1    2
        idx2 0,1  2,3  4,5
        
        looping backwards to keep indexes correct
        */
        for (int i = (mask.size() - 1); i > -1; --i) {
            if (mask.get(i) == 0) {
                int idx2 = 2*i;
                keypoints.remove(idx2 + 1);
                keypoints.remove(idx2);
                orientations.remove(idx2 + 1);
                orientations.remove(idx2);
                responses.remove(idx2 + 1);
                responses.remove(idx2);
            }
        }
       
        //holds values 1 or 0.  size is [orientations.size][POS0.length]
        int[][] descriptors = orbLoop(octaveImage, keypoints, orientations);
        
        Descriptors desc = new Descriptors();
        desc.descriptors = descriptors;
        desc.mask = mask;
        
        return desc;        
    }
    
    /**
     * create descriptors for the given keypoints.
     * 
     * adapted from
     * https://github.com/scikit-image/scikit-image/blob/master/skimage/feature/orb_cy.pyx
     * 
     * @param octaveImage
     * @param keypoints
     * @param orientations
     * @return 
     *   holds values 1 or 0.  size is [orientations.size][POS0.length]
     */
    protected int[][] orbLoop(float[][] octaveImage, TIntList keypoints,
        TDoubleList orientations) {
        
        assert(orientations.size() == keypoints.size()/2);
        
        if (POS0 == null) {
            POS0 = ORBDescriptorPositions.POS0;
        }
        if (POS1 == null) {
            POS1 = ORBDescriptorPositions.POS1;
        }
        
        int nRows = octaveImage.length;
        int nCols = octaveImage[0].length;
        
        int nKP = orientations.size();
        
        // holds values 1 or 0.  size is [orientations.size][POS0.length]
        int[][] descriptors = new int[nKP][POS0.length];
        for (int i = 0; i < descriptors.length; ++i) {
            descriptors[i] = new int[POS0.length];
        }
        
        double pr0, pc0, pr1, pc1;
        int spr0, spc0, spr1, spc1;
        
        for (int i = 0; i < descriptors.length; ++i) {
            double angle = orientations.get(i);
            double sinA = Math.sin(angle);
            double cosA = Math.cos(angle);
            
            int idx = 2*i;
            int kr = keypoints.get(idx);
            int kc = keypoints.get(idx + 1);

            for (int j = 0; j < descriptors[i].length; ++j) {
                pr0 = POS0[j][0];
                pc0 = POS0[j][1];
                pr1 = POS1[j][0];
                pc1 = POS1[j][1];

                spr0 = (int)Math.round(sinA * pr0 + cosA * pc0);
                spc0 = (int)Math.round(cosA * pr0 - sinA * pc0);
                spr1 = (int)Math.round(sinA * pr1 + cosA * pc1);
                spc1 = (int)Math.round(cosA * pr1 - sinA * pc1);

                int x0 = kr + spr0;
                int y0 = kc + spc0;
                int x1 = kr + spr1;
                int y1 = kc + spc1;
                if (octaveImage[x0][y0] < octaveImage[x1][y1]) {
                    descriptors[i][j] = 1;
                }
            }
        }
        
        return descriptors;
    }
    
    /**
     * get a list of each octave's descriptors.  NOTE that the list
     * is not copied so do not modify.
     * The coordinates of the descriptors can be found in keyPointsList, but
     * with twice the spacing because that stores row and col in same list.
     * @return 
     */
    public List<Descriptors> getDescriptorsList() {
        return descriptorsList;
    }
        
    /**
     * get a list of each octave's descriptors as a combined descriptor.
     * The coordinates of the descriptors can be found in getAllKeypoints, but
     * with twice the spacing because that stores row and col in same list.
     * @return 
     */
    public Descriptors getAllDescriptors() {
        
        int n = 0;
        for (Descriptors dl : descriptorsList) {
            n += dl.descriptors.length;
        }
        
        int[][] combinedD = new int[n][POS0.length];
        for (int i = 0; i < n; ++i) {
            combinedD[i] = new int[POS0.length];
        }
        
        TIntList combinedM = new TIntArrayList(n);
        
        int count = 0;
        for (int i = 0; i < descriptorsList.size(); ++i) {
            
            Descriptors dl = descriptorsList.get(i);
            combinedM.addAll(dl.mask);
            int[][] d = dl.descriptors;
            
            for (int j = 0; j < d.length; ++j) {
                System.arraycopy(d[j], 0, combinedD[count], 0, d[j].length);
                count++;
            }            
        }
        
        Descriptors combined = new Descriptors();
        combined.descriptors = combinedD;
        combined.mask = combinedM;
        
        return combined;
    }
    
    /**
     * get a list of each octave's keypoints as a combined list.
     * The list contains coordinates which have already been scaled to the
     * full image reference frame and are present in format [row0, col0,
     * row1, col1, ...].  The corresponding points in the other lists are
     * at index/2.
     * @return 
     */
    public TIntList getAllKeyPoints() {
        
        int n = 0;
        for (TIntList ks : keypointsList) {
            n += ks.size();
        }
        
        TIntList combined = new TIntArrayList(n);
        
        for (int i = 0; i < keypointsList.size(); ++i) {
            
            TIntList ks = keypointsList.get(i);
            
            combined.addAll(ks);            
        }
        
        return combined;
    }
    
    /**
     * get a list of each octave's orientations as a combined list.
     * The corresponding coordinates can be obtained with getAllKeyPoints();
     * @return 
     */
    public TDoubleList getAllOrientations() {
        
        int n = 0;
        for (TDoubleList ks : orientationsList) {
            n += ks.size();
        }
        
        TDoubleList combined = new TDoubleArrayList(n);
        
        for (int i = 0; i < orientationsList.size(); ++i) {
            
            TDoubleList ks = orientationsList.get(i);
            
            combined.addAll(ks);            
        }
        
        return combined;
    }
    
    /**
     * get a list of each octave's scales as a combined list.
     * The corresponding coordinates can be obtained with getAllKeyPoints();
     * @return 
     */
    public TFloatList getAllScales() {
        
        int n = 0;
        for (TFloatList ks : scalesList) {
            n += ks.size();
        }
        
        TFloatList combined = new TFloatArrayList(n);
        
        for (int i = 0; i < scalesList.size(); ++i) {
            
            TFloatList ks = scalesList.get(i);
            
            combined.addAll(ks);            
        }
        
        return combined;
    }
    
    /**
     * get a list of each octave's harris responses as a combined list.
     * The corresponding coordinates can be obtained with getAllKeyPoints();
     * @return 
     */
    public TFloatList getAllHarrisResponses() {
        
        int n = 0;
        for (TFloatList ks : harrisResponses) {
            n += ks.size();
        }
        
        TFloatList combined = new TFloatArrayList(n);
        
        for (int i = 0; i < harrisResponses.size(); ++i) {
            
            TFloatList ks = harrisResponses.get(i);
            
            combined.addAll(ks);            
        }
        
        return combined;
    }
    
    /**
     * get a list of each octave's keypoints in the reference frame
     * of the original full image size.  NOTE that the list
     * is not copied so do not modify.
     * The format of a keypoint list is [row0, col0, row1, col1, ...] using
     * row-major, C-style notation, that is first dimension is row, second
     * is col w.r.t. image they're extracted from (other classes in this project
     * use col-major notation).
     * @return 
     */
    public List<TIntList> getKeyPointList() {
        return keypointsList;
    }
    
    /**
     * get a list of each octave's orientations.  NOTE that the list
     * is not copied so do not modify.
     * The coordinates of the descriptors can be found in keyPointsList, but
     * with twice the spacing because that stores row and col in same list.
     * @return 
     */
    public List<TDoubleList> getOrientationsList() {
        return orientationsList;
    }
    
    /**
     * get a list of each octave's scales.  NOTE that the list
     * is not copied so do not modify.
     * The coordinates of the descriptors can be found in keyPointsList, but
     * with twice the spacing because that stores row and col in same list.
     * @return 
     */
    public List<TFloatList> getScalesList() {
        return scalesList;
    }
    
    /**
     * get a list of each octave's harris responses.  NOTE that the list
     * is not copied so do not modify.
     * The coordinates of the descriptors can be found in keyPointsList, but
     * with twice the spacing because that stores row and col in same list.
     * @return 
     */
    public List<TFloatList> getHarrisResponseList() {
        return harrisResponses;
    }
    
    /**
     * greedy euclidean SSD matching of d1 to d2, with unique mappings for
     * all indexes.
     * 
     * @param d1
     * @param d2
     * @return matches - two dimensional int array of indexes in d1 and
     * d2 which are matched.
     */
    public static int[][] matchDescriptors(int[][] d1, int[][] d2) {
        
        // 2nd dimension is the descriptor length, same as POS0.length
        assert(d1[0].length == d2[0].length);
        
        int n1 = d1.length;
        int n2 = d2.length;
        
        //[n1][n2]
        float[][] cost = calcDescriptorCostMatrix(d1, d2);
    
        // greedy or optimal match can be performed here.
        
        // NOTE: some matching problems might benefit from using the spatial
        //   information at the same time.  for those, will consider adding 
        //   an evaluation term for these descriptors to a specialization of
        //   PartialShapeMatcher.java
        
        // for the greedy match, separating the index information from the cost
        // and then sorting by cost
        int nTot = d1.length * d2.length;
        
        PairInt[] indexes = new PairInt[nTot];
        float[] costs = new float[nTot];
        int count = 0;
        for (int i = 0; i < n1; ++i) {
            for (int j = 0; j < n2; ++j) {
                indexes[count] = new PairInt(i, j);
                costs[count] = cost[i][j];
                count++;
            }
        }
        assert(count == nTot);
        
        QuickSort.sortBy1stArg(costs, indexes);
        
        TIntSet set1 = new TIntHashSet();
        TIntSet set2 = new TIntHashSet();
        
        List<PairInt> matches = new ArrayList<PairInt>();
        
        // visit lowest costs (== differences) first
        for (int i = 0; i < nTot; ++i) {
            PairInt index12 = indexes[i];
            int idx1 = index12.getX();
            int idx2 = index12.getY();
            if (set1.contains(idx1) || set2.contains(idx2)) {
                continue;
            }
            matches.add(index12);
            set1.add(idx1);
            set2.add(idx2);
        }
        
        int[][] results = new int[matches.size()][2];
        for (int i = 0; i < matches.size(); ++i) {
            results[i][0] = matches.get(i).getX();
            results[i][1] = matches.get(i).getY();
        }
        
        return results;
    }
   
    /**
     * calculate a cost matrix composed of the SSD of each descriptor in d1 to d2.
     * 
     * @param d1 two dimensional array with first being keypoint indexes and
     * second dimension being descriptor.
     * @param d2  two dimensional array with first being keypoint indexes and
     * second dimension being descriptor.
     * @return matches two dimensional int array of indexes in d1 and
     * d2 which are matched.
     */
    public static float[][] calcDescriptorCostMatrix(int[][] d1, int[][] d2) {
        
        // 2nd dimension is the descriptor length, same as POS0.length
        assert(d1[0].length == d2[0].length);
        
        int n1 = d1.length;
        int n2 = d2.length;
        
        float[][] cost = new float[n1][n2];
        for (int i = 0; i < n1; ++i) {
            cost[i] = new float[n2];
            for (int j = 0; j < n2; ++j) {
                cost[i][j] = MiscMath.calculateSSD(d1[i], d2[j], Integer.MIN_VALUE);
            }
        }
        
        return cost;
    }
}
