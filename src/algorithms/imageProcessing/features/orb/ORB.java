package algorithms.imageProcessing.features.orb;

import algorithms.QuickSort;
import algorithms.imageProcessing.ATrousWaveletTransform;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.MedianTransform;
import algorithms.imageProcessing.SIGMA;
import algorithms.imageProcessing.StructureTensor;
import algorithms.imageProcessing.util.MatrixUtil;
import algorithms.util.PairInt;
import algorithms.util.TwoDFloatArray;
import algorithms.util.VeryLongBitString;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.TFloatList;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TFloatIntMap;
import gnu.trove.map.hash.TFloatIntHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.logging.Logger;

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

Have also edited alot of the ordering of logic, but the core content is still adapted
from the scipy code.
 
Note: to include color information, consider combining the results of the
a greyscale ORB with the results of a "C" from LCH colorspace and/or another
* color transformation that uses a chromaticity transformation 
* to distinguish colors.  
* For the later, would need to make a separate ORB because the
* differences in the "H", that is polar theta, would need
* to use a threshold of 20 degrees or so and a result of binary
* difference, and would need a correction for wrap around the
* circular coordinate space, so the keypoint extraction would
* need to be adjusted for that.  The main reason for such 
* consideration is to be able to distinguish between colors whose
* combinations are the same in greyscale color.

NOTE, have chosen to keep the results (which are available publicly via getters)
in Lists in which each list index is a scale in the pyramidal decimation
to allow the user to recover scale information if wanted.
To use the way that scipy does, all results as single lists of variables
are available via getAll*getters.

Still testing the class, there may be bugs present.
 <pre>
 Example Use:
    int nKeyPoints = 200;
    ORB orb = new ORB(nKeyPoints);
    orb.overrideToAlsoCreate1stDerivKeypoints();
    orb.detectAndExtract(image);

    // to get the list of keypoint coordinates in row-major, but separated:
    List TIntList keypoints0 = orb.getAllKeyPoints0(); // rows being first dimension
    List TIntList keypoints1 = orb.getAllKeyPoints1(); // cols being second dimension

    // to get the descriptors:
    List Descriptors descList = orb.getDescriptors();
       or
    Descriptors desc = orb.getAllDescriptors();
    
 * </pre>
 */
public class ORB {
   
    protected boolean useSingleScale = false;

    // these could be made static across all instances, but needs guards for synchronous initialization
    protected int[][] OFAST_MASK = null;
    protected int[] OFAST_UMAX = null;
    protected int[][] POS0 = null;
    protected int[][] POS1 = null;

    protected int fastN = 9;
    // fastThreshold should probably be <= 0.01 for low res images
    protected float fastThreshold = 0.08f;
    protected final float harrisK = 0.04f;
    protected final int nKeypoints;

    /**
    the first dimension of the pyramids is scale
    and the value is the image as a 2 dimensional row major format float array.
    */
    protected TwoDFloatArray[] pyramidImages = null;
    protected float[] scales = null;
    protected StructureTensor[] tensorComponents = null;
    protected TwoDFloatArray[] harrisResponseImages = null;

    /**
     * although there may be multiple color images in the pyramid, there is
     * only one set of coordinates derived for each scale.
     */
    protected List<TIntList> keypoints0List = null;
    protected List<TIntList> keypoints1List = null;
    
    protected List<TDoubleList> orientationsList = null;
    protected List<TFloatList> harrisResponses = null;
    
    //`True`` or ``False`` representing
    //the outcome of the intensity comparison for i-th keypoint on j-th
    //decision pixel-pair. It is ``Q == np.sum(mask)``.
    private List<Descriptors> descriptorsList = null;
    
    protected static double twoPI = 2. * Math.PI;

    protected float curvatureThresh = 0.05f;

    private Logger log = Logger.getLogger(this.getClass().getName());
    
    /**
     * a map holding the scale factor of a pyramid image as key
     * and the octave number of pyramid images.
     * it's meant to be useful for delayed creating of hsv pyramid
     * to read the scalesList and remove an hsv pyramid image if not present in
     * scales.
     */
    private TFloatIntMap octaveScaleMap = new TFloatIntHashMap();

    // pyramid images will be no smaller than this
    private final static int defaultDecimationLimit = 32;
    protected final int decimationLimit = defaultDecimationLimit;

    public static enum DescriptorChoice {
        NONE, GREYSCALE
    }
    private DescriptorChoice descrChoice = DescriptorChoice.GREYSCALE;

    protected boolean doCreate1stDerivKeypoints = false;

    protected boolean doCreateCurvatureKeyPoints = false;

    protected int nPyramidB = 3;

    protected boolean useSmallestPyramid = false;

    protected boolean doCreate1stATrousKeypoints = true;
    
    protected final GreyscaleImage img;
    
    /**
     * Still testing the class, there may be bugs present.
     * @param image
     * @param nKeypoints
     */
    public ORB(GreyscaleImage image, int nKeypoints) {

        initMasks();

        this.img = image;
                
        this.nKeypoints = nKeypoints;
    }
   
    /**
     * @param descrChoice the descrChoice to set
     */
    public void setDescrChoice(DescriptorChoice descrChoice) {
        this.descrChoice = descrChoice;
    }
    
    public void overrideToUseSingleScale() {
        useSingleScale = true;
    }

    protected void overrideFastN(int nFast) {
        this.fastN = nFast;
    }
    public void overrideFastThreshold(float threshold) {
        this.fastThreshold = threshold;
    }

    public void overrideToUseSmallestPyramid() {
        
        if (useSingleScale) {
            throw new IllegalArgumentException("useSingleScale=true so cannot "
                + " set this parameter");
        }
        
        useSmallestPyramid = true;
    }

    /**
     * override the number of divisions of scale to add in
     * between the powers of 2 scalings.  The default number
     * is 3;
     * @param n
     */
    public void overridePyamidalExtraN(int n) {
        
        if (useSingleScale) {
            throw new IllegalArgumentException("useSingleScale=true so cannot "
                + " set this parameter");
        }
        
        this.nPyramidB = n;
    }

    /**
     * set option to not create descriptors.
     */
    public void overrideToNotCreateDescriptors() {
        descrChoice = DescriptorChoice.NONE;
    }
    
    public void overrideToAlsoCreate1stDerivKeypoints() {
        doCreate1stDerivKeypoints = true;
    }
    
    public void overrideToNotCreateATrousKeypoints() {
        doCreate1stATrousKeypoints = false;
    }

    public void overrideToCreateCurvaturePoints() {
        doCreateCurvatureKeyPoints = true;
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
      NOTE: still testing, there may be bugs.

      Detect oriented FAST keypoints and extract rBRIEF descriptors.
      Note that this is faster than first calling `detect` and then
      `extract`.

      code is adapted from
      https://github.com/scikit-image/scikit-image/blob/401c1fd9c7db4b50ae9c4e0a9f4fd7ef1262ea3c/skimage/feature/orb.py
      with some functions replaced by code in this project.

     */
    public void detectAndExtract() {

        buildPyramid();
        
        buildTensors();
        
        buildHarrisImages();
        
        // NOTE that the keypoints are all in the coord frame of the largest image
        extractKeypoints();

        int nKeyPointsTotal = countKeypoints();

        log.info("nKeypointsTotal=" + nKeyPointsTotal +
            " this.nKeypoints=" + this.nKeypoints);
                
        extractDescriptors();
    }
  
    private void buildPyramid() {
       
        List<TwoDFloatArray> pyr = buildPyramid(img);

        int nScales = pyr.size();

        pyramidImages = new TwoDFloatArray[nScales];
        
        scales = new float[nScales];
        
        float prevScl = 1.f;
        scales[0] = prevScl;
        
        for (int i = 0; i < nScales; ++i) {
            
            pyramidImages[i] = pyr.get(i);
            
            float scale = prevScl;
            
            if (i > 0) {
                float x0 = (float)pyramidImages[i - 1].a.length/
                    (float)pyramidImages[i].a.length;
                float y0 = (float)pyramidImages[i - 1].a[0].length/
                    (float)pyramidImages[i].a[0].length;
                scale = prevScl * (x0 + y0)/2.f;
                prevScl = scale;
         
                scales[i] = scale;
            }
        }
    }
  
    /**
     * @param image in col major format
     * @return float array of pyramid images in row major format
     */
    private List<TwoDFloatArray> buildPyramid(GreyscaleImage img) {

        //TODO: need a finer grained decimation as an option
        //   can build one with gaussian and down-sampling
        // of use the ATrous b3 method and down-sampling
        //   see the upsampling code in image processing to reverse...

        ImageProcessor imageProcessor = new ImageProcessor();

        List<GreyscaleImage> output;
        
        int dl = useSingleScale ?
            Math.min(img.getWidth(), img.getHeight()) - 1 :
            decimationLimit;
        
        if (useSmallestPyramid) {
            output = new ArrayList<GreyscaleImage>();
            MedianTransform mt = new MedianTransform();
            mt.multiscalePyramidalMedianTransform2(img, output, dl);
        } else {
           output = imageProcessor.buildPyramid2(img, dl, nPyramidB);
        }

        List<TwoDFloatArray> output2 = new ArrayList<TwoDFloatArray>();
        for (int i = 0; i < output.size(); ++i) {
            float[][] rowMajorImg = imageProcessor.copyToRowMajor(
                output.get(i));
            MatrixUtil.multiply(rowMajorImg, 1.f/255.f);
            TwoDFloatArray f = new TwoDFloatArray(rowMajorImg);
            output2.add(f);
        }

        return output2;
    }
    
    protected void buildTensors() {
    
        // Standard deviation used for the Gaussian kernel, which is used as
        // weighting function for the auto-correlation matrix.
        float sigma = SIGMA.getValue(SIGMA.ZEROPOINTFIVE);
        
        tensorComponents = new StructureTensor[scales.length];
        
        for (int i = 0; i < scales.length; ++i) {
            
            TwoDFloatArray octaveImage = pyramidImages[i];
            
            tensorComponents[i] = new StructureTensor(octaveImage.a,
                sigma, doCreateCurvatureKeyPoints);
        }
        
    }

    private void buildHarrisImages() {
        
        harrisResponseImages = new TwoDFloatArray[scales.length];
        
        for (int i = 0; i < scales.length; ++i) {
            
            TwoDFloatArray octaveImage = pyramidImages[i];
            
            StructureTensor tensors = tensorComponents[i];
            
            float[][] detA = tensors.getDeterminant();

            float[][] traceA = tensors.getTrace();

            // size is same a octaveImage
            float[][] harrisResponse = cornerHarris(octaveImage.a, detA, traceA);

            harrisResponseImages[i] = new TwoDFloatArray(harrisResponse);
        }
    }
    
    protected void debugPrint(String label, float[][] a) {

        StringBuilder sb = new StringBuilder(label);
        sb.append("\n");
        for (int i = 0; i < a.length; ++i) {
            sb.append(Arrays.toString(a[i])).append("\n");
        }

        log.info(sb.toString());
    }

    private void debugPrint(String label, int[][] a) {

        StringBuilder sb = new StringBuilder(label);
        sb.append("\n");
        for (int i = 0; i < a.length; ++i) {
            sb.append(Arrays.toString(a[i])).append("\n");
        }

        log.info(sb.toString());
    }

    protected float[][] multiply(float[][] a, float[][] b) {

        float[][] c = copy(a);

        for (int i = 0; i < c.length; ++i) {
            for (int j = 0; j < c[0].length; ++j) {
                c[i][j] *= b[i][j];
            }
        }

        return c;
    }

    protected float[][] add(float[][] a, float[][] b) {

        float[][] c = copy(a);

        for (int i = 0; i < c.length; ++i) {
            for (int j = 0; j < c[0].length; ++j) {
                c[i][j] += b[i][j];
            }
        }

        return c;
    }

    protected float[][] subtract(float[][] a, float[][] b) {

        float[][] c = copy(a);

        for (int i = 0; i < c.length; ++i) {
            for (int j = 0; j < c[0].length; ++j) {
                c[i][j] -= b[i][j];
            }
        }

        return c;
    }

    /**
     extracts keypoints, orientations, and harris response images
     and stores results in the class level variables.
     NOTE: the octaves are also stored at class level in this method too
     once any keypoints are detected.
     returned is the list of octaves which were used.
    */
    private void extractKeypoints() {

        keypoints0List = new ArrayList<TIntList>();
        keypoints1List = new ArrayList<TIntList>();
        orientationsList = new ArrayList<TDoubleList>();
        harrisResponses = new ArrayList<TFloatList>();
     
        // atrous keypoints if any
        TIntList atk0 = new TIntArrayList();
        TIntList atk1 = new TIntArrayList();
        if (doCreate1stATrousKeypoints) {
            extractATrousKeypoints(atk0, atk1);
        }
        
        Resp r = extractKeypoints(0);

        if ((r == null) || (r.keypoints0 == null) || r.keypoints0.isEmpty()) {
            return;
        }
       
        if (!atk0.isEmpty()) {
            r.keypoints0.addAll(atk0);
            r.keypoints1.addAll(atk1);
        }

        //mask of length orientations.size() containing a 1 or 0
        // indicating if pixels are within the image (``True``) or in the
        // border region of the image (``False``).
        TIntList mask = filterNearBorder(pyramidImages[0].a, 
            r.keypoints0, r.keypoints1);
        
        // order by response and filter
        
        float[][] harrisResponse = harrisResponseImages[0].a;

        int n = r.keypoints0.size();

        float[] scores = new float[n];
        PairInt[] points = new PairInt[n];
        int count = 0;
        for (int ii = 0; ii < n; ++ii) {
            int x = r.keypoints1.get(ii);
            int y = r.keypoints0.get(ii);
            points[count] = new PairInt(x, y);
            //NOTE: because scale[0] is the full size frame, these coords
            // do not need to be rescaled
            scores[count] = harrisResponse[y][x];
            count++;
        }
      
        QuickSort.sortBy1stArg(scores, points);
        
        TIntList kpc0s = new TIntArrayList(this.nKeypoints);
        TIntList kpc1s = new TIntArrayList(this.nKeypoints);
        
        for (int i = (points.length - 1); i > -1; --i) {
            PairInt p = points[i];
            int x = p.getX();
            int y = p.getY();
            kpc1s.add(x);
            kpc0s.add(y);
        }

        // then filter
        maskCoordinates(kpc0s, kpc1s, harrisResponse.length,
            harrisResponse[0].length, 4);//8);

        // at this point kpc0s, kpc1s are in coord frame of largest image

        // filter for number of points if requested
        if (kpc0s.size() > this.nKeypoints) {
            while (kpc0s.size() > this.nKeypoints) {
                int idx = kpc0s.size() - 1;
                kpc0s.removeAt(idx);
                kpc1s.removeAt(idx);
            }
        }
        
        // populate the lists for each octave

        for (int octave = 0; octave < scales.length; ++octave) {

            float scale = scales[octave];
            harrisResponse = harrisResponseImages[octave].a;

            // make a keypoint list for each scale, using set first
            // to remove redundant points.
            // TODO: may need to apply a local max filter if points cluster in
            //    smaller image
            Set<PairInt> set = new HashSet<PairInt>();
            TIntList kp0s = new TIntArrayList(set.size());
            TIntList kp1s = new TIntArrayList(set.size());
            for (int i = 0; i < kpc0s.size(); ++i) {
                // reduce to octave coord system
                int c0 = Math.round((float)kpc0s.get(i)/scale);
                int c1 = Math.round((float)kpc1s.get(i)/scale);
                if (c0 > (harrisResponse.length - 1)) {
                    c0 = harrisResponse.length - 1;
                }
                if (c1 > (harrisResponse[0].length - 1)) {
                    c1 = harrisResponse[0].length - 1;
                }
                PairInt p = new PairInt(c0, c1);
                if (set.contains(p)) {
                    continue;
                }
                set.add(p);
                kp0s.add(c0);
                kp1s.add(c1);
            }

            // then filter for best within a pixel range
            maskCoordinates(kp0s, kp1s, harrisResponse.length,
                harrisResponse[0].length, 4);//8);

            TFloatList responses = new TFloatArrayList(kp0s.size());

            for (int i = 0; i < kp0s.size(); ++i) {
                int c0 = kp0s.get(i);
                int c1 = kp1s.get(i);
                float v = harrisResponse[c0][c1];
                responses.add(v);
            }

            float[][] octaveImage = pyramidImages[octave].a;

            TDoubleList orientations = calculateOrientations(octaveImage,
                kp0s, kp1s);

            // transform keypoints to full size coordinate reference frame
            for (int i = 0; i < kp0s.size(); ++i) {
                int v = Math.round(scale * kp0s.get(i));
                kp0s.set(i, v);

                v = Math.round(scale * kp1s.get(i));
                kp1s.set(i, v);
            }

            // kp0s and kp1s are in coord frame of largest image

            keypoints0List.add(kp0s);
            keypoints1List.add(kp1s);
            orientationsList.add(orientations);
            harrisResponses.add(responses);

        } // end loop over octave

        assert(pyramidImages.length == keypoints0List.size());
        assert(pyramidImages.length == keypoints1List.size());
        assert(pyramidImages.length == orientationsList.size());
        assert(pyramidImages.length == this.harrisResponses.size());
        assert(pyramidImages.length == scales.length);        
    }

    private void extractDescriptors() {

        if (descrChoice.equals(DescriptorChoice.NONE)) {
            return;
        }
        
        int nScales = scales.length;
        
        descriptorsList = new ArrayList<Descriptors>(nScales);
    
        for (int i = 0; i < nScales; ++i) {

            float[][] octaveImage = pyramidImages[i].a;

            TIntList kp0 = this.keypoints0List.get(i);
            TIntList kp1 = this.keypoints1List.get(i);
            TDoubleList or = this.orientationsList.get(i);

            // result contains descriptors and mask.
            // also, modified by mask are the keypoints and orientations
            Descriptors desc = extractOctaveDescriptor(octaveImage, kp0, kp1, or,
                scales[i]);

            descriptorsList.add(desc);
        }
    }

    private int countKeypoints() {
        int n = 0;
        for (TIntList kp0 : keypoints0List) {
            n += kp0.size();
        }
        return n;
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

    protected class Resp {
        TIntList keypoints0;
        TIntList keypoints1;
    }

    public static class Descriptors {
        //NOTE: consider packing 4 descriptors into one or
        // doing as Image.java does with sensing 64 bit and 32 bit to make
        // long or int bit vectors
        public VeryLongBitString[] descriptors;
    }

    /**
     * adapted from
     * https://github.com/scikit-image/scikit-image/blob/master/skimage/feature/orb.py
     *
     * @param octave
     * @return
     */
    private Resp extractKeypoints(int octave) {

        float[][] octaveImage = pyramidImages[octave].a;
        
        float scale = scales[octave];
        
        StructureTensor tensors = tensorComponents[octave];
        
        float[][] fastResponse = cornerFast(octaveImage, fastN, fastThreshold);

        int nRows = fastResponse.length;
        int nCols = fastResponse[0].length;

        TIntList keypoints0 = new TIntArrayList();
        TIntList keypoints1 = new TIntArrayList();
 
        // list of format [row, col, ...] of filtered maxima ordered by intensity
        cornerPeaks(fastResponse, 1, keypoints0, keypoints1);
        if (keypoints0.isEmpty()) {
            return null;
        }

        maskCoordinates(keypoints0, keypoints1, nRows, nCols, 8);//16);

        if (doCreate1stDerivKeypoints) {

            ImageProcessor imageProcessor = new ImageProcessor();

            float hLimit = 0.01f;//0.09f;

            imageProcessor.createFirstDerivKeyPoints(
                tensors, keypoints0, keypoints1, hLimit);
            
            /*try {
                float factor = 255.f;
                Image img2 = new Image(nCols, nRows);
                for (int i = 0; i < nRows; ++i) {
                    for (int j = 0; j < nCols; ++j) {
                        int v = Math.round(factor * octaveImage[i][j]);
                        if (v > 255) {
                            v = 255;
                        } else if (v < 0) {
                            v = 0;
                        }
                        img2.setRGB(j, i, v, v, v);
                    }
                }
                for (int i = 0; i < keypoints0.size(); ++i) {
                    int x = keypoints1.get(i);
                    int y = keypoints0.get(i);
                    img2.setRGB(x, y, 255, 0, 0);
                }
                System.out.println("nRows=" + nRows + " nCols=" + nCols);
                //algorithms.imageProcessing.ImageDisplayer.displayImage("first deriv", img2);
                MiscDebug.writeImage(img2, "_fr_" 
                    + MiscDebug.getCurrentTimeFormatted());
                int z = 1;
            } catch(Exception e) {
                System.out.println(e.getMessage());
            }*/
        }

        if (doCreateCurvatureKeyPoints) {

            ImageProcessor imageProcessor = new ImageProcessor();

            // usually, only create these points for points on an edge.
            // wanting the min and max of curvature,
            // and then those maxima that are 2 or 3 times stronger than
            // one of the adjacent minima.
            // with a single edge, the peak curvature should be larger than
            // 2 times that of the preceding or proceeding minima.

            //default thresh is 0.01f
            imageProcessor.createCurvatureKeyPoints(tensors, 
                keypoints0, keypoints1, curvatureThresh);

            /*try {
                float factor = 255.f;
                Image img2 = new Image(nCols, nRows);
                for (int i = 0; i < nRows; ++i) {
                    for (int j = 0; j < nCols; ++j) {
                        int v = Math.round(factor * octaveImage[i][j]);
                        if (v > 255) {
                            v = 255;
                        } else if (v < 0) {
                            v = 0;
                        }
                        img2.setRGB(j, i, v, v, v);
                    }
                }
                for (int i = 0; i < keypoints0.size(); ++i) {
                    int x = keypoints1.get(i);
                    int y = keypoints0.get(i);
                    img2.setRGB(x, y, 255, 0, 0);
                }
                System.out.println("nRows=" + nRows + " nCols=" + nCols);
                //algorithms.imageProcessing.ImageDisplayer.displayImage("curvature", img2);
                MiscDebug.writeImage(img2, "_curve_" 
                    + MiscDebug.getCurrentTimeFormatted());
                int z = 1;
            } catch(Exception e) {
                System.out.println(e.getMessage());
            }*/
        }
        
        // size is same a octaveImage
        float[][] harrisResponse = harrisResponseImages[octave].a;

        // remove redundant keypoints
        Set<PairInt> exists = new HashSet<PairInt>();
        for (int i = (keypoints0.size() - 1); i > -1; --i) {
            PairInt p = new PairInt(keypoints0.get(i),
               keypoints1.get(i));
            if (exists.contains(p)) {
                keypoints0.removeAt(i);
                keypoints1.removeAt(i);
            } else {
                exists.add(p);
            }
        }

        // --- harris corners from response image ----
        ImageProcessor imageProcessor = new ImageProcessor();
        imageProcessor.peakLocalMax(harrisResponse, 1, 0.1f,
            keypoints0, keypoints1);

        Resp r2 = new Resp();
        r2.keypoints0 = keypoints0;
        r2.keypoints1 = keypoints1;
        
        
        /*{//DEBUG
            Image dbg = new Image(harrisResponse[0].length,
                harrisResponse.length);
            float[] tmp = new float[dbg.getNPixels()];
            tmp = new float[dbg.getNPixels()];
            for (int i = 0; i < dbg.getNPixels(); ++i) {
                int x = dbg.getCol(i);
                int y = dbg.getRow(i);
                tmp[i] = harrisResponse[y][x];
            }
            MiscMath.rescale(tmp, 0, 255);
            for (int i = 0; i < dbg.getNPixels(); ++i) {
                if (tmp[i] > 0) {
                    dbg.setRGB(i, 255, 255, 255);
                }
            }
            PairIntArray tmp2 = new PairIntArray(keypoints0.size());
            for (int i = 0; i < keypoints0.size(); ++i) {
                int x = keypoints1.get(i);
                int y = keypoints0.get(i);
                tmp2.add(x, y);
            }
            ImageIOHelper.addCurveToImage(tmp2, dbg, 2, 255, 0, 0);
            MiscDebug.writeImage(dbg, "_hr_" + 
                MiscDebug.getCurrentTimeFormatted());
            System.out.println("nKP=" + keypoints0.size() + " for octave=" + 
                octave);
        }*/
        
        
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
          Input image as 2D array in row-major, C-style ordered array
     @param n
          Minimum number of consecutive pixels out of 16 pixels on the circle
          that should all be either brighter or darker w.r.t testpixel.
          A point c on the circle is darker w.r.t test pixel p if
          `Ic lessThan Ip - threshold` and brighter if `Ic greaterThan Ip + threshold`.
          Also stands for the n in `FAST-n` corner detector.
          (the scipy default is n=12).
     @param threshold
          Threshold used in deciding whether the pixels on the circle are
          brighter, darker or similar w.r.t. the test pixel. Decrease the
          threshold when more corners are desired and vice-versa.
          (the scipy default is threshold=0.15)
     @return
        FAST corner response image.
     */
    protected float[][] cornerFast(final float[][] img, final int n, final float threshold) {

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
                //# High speed test for n >= 12
                if (n >= 12) {
                    speed_sum_b = 0;
                    speed_sum_d = 0;
                    for (k = 0; k < 16; k += 4) {
                        if (bins[k] == 'b') {
                            speed_sum_b += 1;
                        } else if (bins[k] == 'd') {
                            speed_sum_d += 1;
                        }
                    }
                    if ((speed_sum_d < 3) && (speed_sum_b < 3)) {
                        continue;
                    }
                }

                //# Test for bright pixels
                currResponse = _corner_fast_response(curr_pixel,
                    circleIntensities, bins, 'b', n);

                //# Test for dark pixels
                if (currResponse == 0) {
                    currResponse = _corner_fast_response(curr_pixel,
                        circleIntensities, bins, 'd', n);
                }

                if (currResponse > 0) {
                    cornerResponse[i][j] = currResponse;
                }
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
    private float _corner_fast_response(final float curr_pixel,
        float[] circleIntensities, char[] bins, final char state, final int n) {

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

    private void maskCoordinates(TIntList keypoints0, TIntList keypoints1,
        int nRows, int nCols, int maskRadius) {

        assert(keypoints0.size() == keypoints1.size());

        TIntSet peaks = new TIntHashSet();

        /* making single pixel index out of coordinates:
        (row * width) + col
        pixIdxs[count] = (j * nRows) + i;
        */
        for (int i = 0; i < keypoints0.size(); ++i) {
            int ii = keypoints0.get(i);
            int jj = keypoints1.get(i);
            int pixIdx = (jj * nRows) + ii;
            peaks.add(pixIdx);
        }

        for (int i = 0; i < keypoints0.size(); ++i) {
            int ii = keypoints0.get(i);
            int jj = keypoints1.get(i);
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

        TIntList keypoints0_2 = new TIntArrayList(peaks.size());
        TIntList keypoints1_2 = new TIntArrayList(peaks.size());
        for (int i = 0; i < keypoints0.size(); ++i) {
            int ii = keypoints0.get(i);
            int jj = keypoints1.get(i);
            int pixIdx = (jj * nRows) + ii;
            if (peaks.contains(pixIdx)) {
                keypoints0_2.add(ii);
                keypoints1_2.add(jj);
            }
        }

        keypoints0.clear();
        keypoints0.addAll(keypoints0_2);

        keypoints1.clear();
        keypoints1.addAll(keypoints1_2);
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
    private TIntList maskCoordinatesIfBorder(TIntList coords0, TIntList coords1,
        int nRows, int nCols, int border) {

        TIntSet set = new TIntHashSet();

        /* making single pixel index out of coordinates:
        (row * width) + col
        pixIdxs[count] = (j * nRows) + i;
        */
        for (int i = 0; i < coords0.size(); ++i) {
            int ii = coords0.get(i);
            int jj = coords1.get(i);
            int pixIdx = (jj * nRows) + ii;
            set.add(pixIdx);
        }

        int nBefore = set.size();

        for (int i = 0; i < coords0.size(); ++i) {
            int ii = coords0.get(i);
            int jj = coords1.get(i);
            int pixIdx = (jj * nRows) + ii;
            if (!set.contains(pixIdx)) {
                continue;
            }
            if ((ii < border) || (jj < border) || (ii > (nRows - border - 1)) ||
                (jj > (nCols - border - 1))) {
                set.remove(pixIdx);
            }
        }

        //System.out.println("nBefore border rm=" + nBefore +
        //    " nAfter=" + set.size());

        /*
        Mask indicating if pixels are within the image (``True``) or in the
            border region of the image (``False``).
        */

        TIntList mask = new TIntArrayList(coords0.size());
        for (int i = 0; i < coords0.size(); ++i) {
            int ii = coords0.get(i);
            int jj = coords1.get(i);
            int pixIdx = (jj * nRows) + ii;
            if (set.contains(pixIdx)) {
                mask.add(1);
            } else {
                mask.add(0);
            }
        }
        assert(coords0.size() == mask.size());
        assert(coords0.size() == coords1.size());

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
     * @param outputKeypoints0 the output row coordinates of keypoints
     * @param outputKeypoints1 the output col coordinates of keypoints
     */
    protected void cornerPeaks(float[][] img, int minDistance,
        TIntList outputKeypoints0, TIntList outputKeypoints1) {

        //threshold_abs=None,
        float thresholdRel = 0.1f;
        boolean excludeBorder = true;
        boolean indices = true;
        int numPeaks = Integer.MAX_VALUE;
        //footprint=None, labels=None

        int nRows = img.length;
        int nCols = img[0].length;

        ImageProcessor imageProcessor = new ImageProcessor();
        // these results have been sorted by decreasing intensity
        imageProcessor.peakLocalMax(img, minDistance, thresholdRel,
            outputKeypoints0, outputKeypoints1);

        //System.out.println("keypoints in cornerPeaks="
        //    + "rows=" + outputKeypoints0.toString()
        //    + "cols=" + outputKeypoints1.toString()
        //    + "\nsize=" + outputKeypoints0.size());

        maskCoordinates(outputKeypoints0, outputKeypoints1, nRows, nCols,
            minDistance);
    }

    public static class Pix implements Comparable<Pix> {

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
        Compute the orientation of corners.
        The orientation of corners is computed using the first order central moment
        i.e. the center of mass approach. The corner orientation is the angle of
        the vector from the corner coordinate to the intensity centroid in the
        local neighborhood around the corner calculated using first order central
        moment.
        adapted from
        from https://github.com/scikit-image/scikit-image/blob/master/skimage/feature/corner_cy.pyx

        <pre>
            +90
        135  |  +45
             |
       180-------- 0
             |
       -135  |   -45
            -90
        </pre>


        References
        ----------
        .. [1] Ethan Rublee, Vincent Rabaud, Kurt Konolige and Gary Bradski
              "ORB : An efficient alternative to SIFT and SURF"
              http://www.vision.cs.chubu.ac.jp/CV-R/pdf/Rublee_iccv2011.pdf
        .. [2] Paul L. Rosin, "Measuring Corner Properties"
              http://users.cs.cf.ac.uk/Paul.Rosin/corner2.pdf
     * @param octaveImage
           Input grayscale image.
     * @param keypoints0
     * @param keypoints1
     * @localParam OFAST_MASK
     *     Mask defining the local neighborhood of the corner used for the
           calculation of the central moment.
     * @return
     *    orientations : (N, 1) array
               Orientations of corners in the range [-pi, pi].
     */
    public TDoubleList calculateOrientations(float[][] octaveImage,
        TIntList keypoints0, TIntList keypoints1) {

        //same as mask, same 0's and 1's:
        //cdef unsigned char[:, ::1] cmask = np.ascontiguousarray(mask != 0, dtype=np.uint8)

        int i, r, c, r0, c0;
        int nMaskRows = OFAST_MASK.length;
        int nMaskCols = OFAST_MASK[0].length;
        int nMaskRows2 = (nMaskRows - 1) / 2;
        int nMaskCols2 = (nMaskCols - 1) / 2;

        //cdef double[:, :] cimage = np.pad(image, (mrows2, mcols2), mode='constant',
        //    constant_values=0)
        int nRows2 = octaveImage.length + (2 * nMaskRows2);
        int nCols2 = octaveImage[0].length + (2 * nMaskCols2);
        double[][] cImage = new double[nRows2][nCols2];
        for (i = 0; i < nRows2; ++i) {
            cImage[i] = new double[nCols2];
            if ((i >= nMaskRows2) && (i < (nRows2 - nMaskRows2 - 1))) {
                float[] src = octaveImage[i - nMaskRows2];
                double[] dest = cImage[i];
                for (int ii = 0; ii < src.length; ++ii) {
                    dest[ii + nMaskCols2] = src[ii];
                }
            }
        }

        // number of corner coord pairs
        int nCorners = keypoints0.size();

        TDoubleList orientations = new TDoubleArrayList(nCorners);

        double curr_pixel;
        double m01, m10, m01_tmp;

        for (i = 0; i < nCorners; ++i) {
            r0 = keypoints0.get(i);
            c0 = keypoints1.get(i);

            m01 = 0;
            m10 = 0;

            for (r = 0; r < nMaskRows; ++r) {
                m01_tmp = 0;
                for (c = 0; c < nMaskCols; ++c) {
                    if (OFAST_MASK[r][c] > 0) {
                        curr_pixel = cImage[r0 + r][c0 + c];
                        m10 += curr_pixel * (c - nMaskCols2);
                        m01_tmp += curr_pixel;
                    }
                }
                m01 += m01_tmp * (r - nMaskRows2);
            }

            //arc tangent of y/x, in the interval [-pi,+pi] radians
            orientations.add(Math.atan2(m01, m10));

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
    protected float[][] cornerHarris(float[][] image, float[][] detA,
        float[][] traceA) {

        // method = 'k'.  k is Sensitivity factor to separate corners from edges,
        // Small values of k result in detection of sharp corners.
        float k = this.harrisK;

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
     * filter the keypoints, orientations, and responses to remove those close
     * to the border and then create descriptors for the remaining.
     *
     * @param octaveImage
     * @param keypoints0
     * @param keypoints1
     * @return the encapsulated descriptors and mask
     */
    protected TIntList filterNearBorder(float[][] octaveImage,
        TIntList keypoints0, TIntList keypoints1) {

        int nRows = octaveImage.length;
        int nCols = octaveImage[0].length;

        TIntList mask = maskCoordinatesIfBorder(keypoints0, keypoints1,
            nRows, nCols, 5);//16);

        for (int i = (mask.size() - 1); i > -1; --i) {
            if (mask.get(i) == 0) {
                keypoints0.removeAt(i);
                keypoints1.removeAt(i);
            }
        }

        return mask;
    }

    /**
     * filter the keypoints, orientations, and responses to remove those close
     * to the border and then create descriptors for the remaining.
     *
     * @param octaveImage
     * @param keypoints0
     * @param keypoints1
     * @param orientations
     * @param scale
     * @return the encapsulated descriptors and mask
     */
    protected Descriptors extractOctaveDescriptor(float[][] octaveImage,
        TIntList keypoints0, TIntList keypoints1,
        TDoubleList orientations, float scale) {

        assert(orientations.size() == keypoints0.size());
        assert(orientations.size() == keypoints1.size());

        Descriptors desc = new Descriptors();

        if (descrChoice.equals(DescriptorChoice.NONE)) {
            desc.descriptors = new VeryLongBitString[0];
        } else {
            desc.descriptors = orbLoop(octaveImage, keypoints0, keypoints1,
                orientations, scale);
        }

        return desc;
    }

    /**
     * create descriptors for the given keypoints.
     *
     * adapted from
     * https://github.com/scikit-image/scikit-image/blob/master/skimage/feature/orb_cy.pyx
     *
     * @param octaveImage
     * @param keypoints0
     * @param keypoints1
     * @param orientations
     * @return
     * array of bit vectors of which only 256 bits are used
     * length is [orientations.size]
     */
    private VeryLongBitString[] orbLoop(float[][] octaveImage, TIntList keypoints0,
        TIntList keypoints1, TDoubleList orientations, float scale) {
      
        assert(orientations.size() == keypoints0.size());

        if (POS0 == null) {
            POS0 = ORBDescriptorPositions.POS0;
        }
        if (POS1 == null) {
            POS1 = ORBDescriptorPositions.POS1;
        }

        int nKP = orientations.size();

        log.fine("nKP=" + nKP);

        // each item is a 256 length bit string
        VeryLongBitString[] descriptors = new VeryLongBitString[nKP];

        double pr0, pc0, pr1, pc1;
        int spr0, spc0, spr1, spc1;

        for (int i = 0; i < descriptors.length; ++i) {

            descriptors[i] = new VeryLongBitString(POS1.length);

            double angle = orientations.get(i);
            double sinA = Math.sin(angle);
            double cosA = Math.cos(angle);

            int kr = keypoints0.get(i);
            int kc = keypoints1.get(i);

            // put kr and kc into pyramid image reference frame.
            kr = (int)(kr/scale);
            kc = (int)(kc/scale);

            for (int j = 0; j < POS0.length; ++j) {
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

                if (x0 < 0 || y0 < 0 || x1 < 0 || y1 < 0 ||
                    (x0 > (octaveImage.length - 1)) ||
                    (x1 > (octaveImage.length - 1)) ||
                    (y0 > (octaveImage[0].length - 1)) ||
                    (y1 > (octaveImage[0].length - 1))
                    ) {
                    continue;
                }
                if (octaveImage[x0][y0] < octaveImage[x1][y1]) {
                    descriptors[i].setBit(j);
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

        return combineDescriptors(descriptorsList);
    }

    /**
     * get a list of each octave's descriptors as a combined descriptor.
     * The coordinates of the descriptors can be found in getAllKeypoints, but
     * with twice the spacing because that stores row and col in same list.
     * @return
     */
    public static Descriptors combineDescriptors(List<Descriptors> list) {

        int n = 0;
        for (Descriptors dl : list) {
            n += dl.descriptors.length;
        }

        int nPos0 = ORBDescriptorPositions.POS0.length;

        VeryLongBitString[] combinedD = new VeryLongBitString[n];

        int count = 0;
        for (int i = 0; i < list.size(); ++i) {

            VeryLongBitString[] d = list.get(i).descriptors;

            //TODO: may want to revisit this and use the instance copy methods individually
            System.arraycopy(d, 0, combinedD, count, d.length);
            count += d.length;
        }

        Descriptors combined = new Descriptors();
        combined.descriptors = combinedD;

        return combined;
    }

    /**
     * get a list of each octave's keypoint rows as a combined list.
     * The list contains coordinates which have already been scaled to the
     * full image reference frame.
     * @return
     */
    public TIntList getAllKeyPoints0() {

        int n = 0;
        for (TIntList ks : keypoints0List) {
            n += ks.size();
        }

        TIntList combined = new TIntArrayList(n);

        for (int i = 0; i < keypoints0List.size(); ++i) {

            TIntList ks = keypoints0List.get(i);

            combined.addAll(ks);
        }

        return combined;
    }

    /**
     * get a list of each octave's keypoint rows as a combined list.
     * The list contains coordinates which have already been scaled to the
     * full image reference frame and are in column major format.
     * @return
     */
    public List<PairInt> getAllKeyPoints() {

        int n = 0;
        for (TIntList ks : keypoints0List) {
            n += ks.size();
        }

        List<PairInt> combined = new ArrayList<PairInt>(n);

        for (int i = 0; i < keypoints0List.size(); ++i) {

            TIntList kp0 = keypoints0List.get(i);
            TIntList kp1 = keypoints1List.get(i);

            for (int j = 0; j < kp0.size(); ++j) {
                int y = kp0.get(j);
                int x = kp1.get(j);
                PairInt p = new PairInt(x, y);
                combined.add(p);
            }
        }

        return combined;
    }

    /**
     * get a list of each octave's keypoint cols as a combined list.
     * The list contains coordinates which have already been scaled to the
     * full image reference frame.
     * @return
     */
    public TIntList getAllKeyPoints1() {

        int n = 0;
        for (TIntList ks : keypoints1List) {
            n += ks.size();
        }

        TIntList combined = new TIntArrayList(n);

        for (int i = 0; i < keypoints1List.size(); ++i) {

            TIntList ks = keypoints1List.get(i);

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
     * get a list of an octave's keypoints in column major format,
     * that is pairint(x=column, y=row).
     * @return
     */
    public List<PairInt> getKeyPointListColMaj(int octave) {
        List<PairInt> out = new ArrayList<PairInt>();
        for (int i = 0; i < keypoints0List.get(octave).size(); ++i) {
            int x = keypoints1List.get(octave).get(i);
            int y = keypoints0List.get(octave).get(i);
            out.add(new PairInt(x, y));
        }
        return out;
    }
    
    /**
     * get a list of an octave's keypoints in pixel indexes,
     * that is pixIdx = (y * imageWidth) + x.
     * @return
     */
    public TIntList getKeyPointListPix(int octave) {
        int w = this.pyramidImages[octave].a[0].length;
        
        TIntList out = new TIntArrayList();
        for (int i = 0; i < keypoints0List.get(octave).size(); ++i) {
            int x = keypoints1List.get(octave).get(i);
            int y = keypoints0List.get(octave).get(i);
            int pixIdx = (y * w) + x;
            out.add(pixIdx);
        }
        return out;
    }
    
    /**
     * get a list of each octave's keypoint rows in the reference frame
     * of the original full image size.  NOTE that the list
     * is not copied so do not modify.
     * @return
     */
    public List<TIntList> getKeyPoint0List() {
        return keypoints0List;
    }

    /**
     * get a list of each octave's keypoint cols in the reference frame
     * of the original full image size.  NOTE that the list
     * is not copied so do not modify.
     * @return
     */
    public List<TIntList> getKeyPoint1List() {
        return keypoints1List;
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
     * @return
     */
    public float[] getScales() {
        return scales;
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
     * calculate a cost matrix composed of the sum of XOR of each descriptor in d1 to d2.
     *
     * @param d1 array of bit vectors wherein 256 bits of the integer are used
     * @param d2 array of bit vectors wherein 256 bits of the integer are used
     * @return matches two dimensional int array of indexes in d1 and
     * d2 which are matched.
     */
    public static int[][] calcDescriptorCostMatrix(
        VeryLongBitString[] d1, VeryLongBitString[] d2) {

        int n1 = d1.length;
        int n2 = d2.length;

        int[][] cost = new int[n1][n2];
        for (int i = 0; i < n1; ++i) {
            cost[i] = new int[n2];
            for (int j = 0; j < n2; ++j) {
                cost[i][j] = (int)d1[i].nBitsDifferent(d2[j]);
            }
        }

        return cost;
    }

    /**
     * calculate a cost matrix composed of the sum of XOR of each descriptor in d1 to d2.
     *
     * @param desc1 two dimensional array with first being keypoint indexes and
     * second dimension being descriptor.
     * @param desc2  two dimensional array with first being keypoint indexes and
     * second dimension being descriptor.
     * @return matches two dimensional int array of indexes in d1 and
     * d2 which are matched.
     */
    public static int[][] calcDescriptorCostMatrix(
        Descriptors[] desc1, Descriptors[] desc2) {

        assert(desc1.length == desc2.length);

        int nd = desc1.length;

        int n1 = desc1[0].descriptors.length;
        int n2 = desc2[0].descriptors.length;

        // d1 contains multiple descriptors for same points, such as
        // descriptors for H, S, and V

        int[][] cost = new int[n1][n2];
        for (int i = 0; i < n1; ++i) {
            cost[i] = new int[n2];
        }

        for (int k = 0; k < nd; ++k) {
            VeryLongBitString[] d1 = desc1[k].descriptors;
            VeryLongBitString[] d2 = desc2[k].descriptors;
            assert(d1.length == n1);
            assert(d2.length == n2);
            for (int i = 0; i < n1; ++i) {
                for (int j = 0; j < n2; ++j) {
                    //cost[i][j] += (nbits * nbits); could use sq diff...
                    cost[i][j] += (int)d1[i].nBitsDifferent(d2[j]);
                }
            }
        }

        return cost;
    }

    /**
     * create col major image from row major input
     * @param a
     * @return
     */
    public static Image convertToImage(TwoDFloatArray a) {
        int n1 = a.a.length;
        int n2 = a.a[0].length;
        Image img = new Image(n2, n1);
        for (int i = 0; i < n1; ++i) {
            for (int j = 0; j < n2; ++j) {
                float v = 255.f * a.a[i][j];
                int vInt = Math.round(v);
                if (vInt > 255) {
                    vInt = 255;
                }
                img.setRGB(j, i, vInt, vInt, vInt);
            }
        }
        return img;
    }
    
    private void extractATrousKeypoints(TIntList outputKP0, TIntList outputKP1) {
        
        List<GreyscaleImage> transformed = new ArrayList<GreyscaleImage>();
        List<GreyscaleImage> coeffs = new ArrayList<GreyscaleImage>();

        ATrousWaveletTransform at = new ATrousWaveletTransform();
        at.calculateWithB3SplineScalingFunction(img, transformed, coeffs, 2);

        if (coeffs.size() >= 2) {

            ImageProcessor imageProcessor = new ImageProcessor();
            float[][] rowMajorImg = imageProcessor.copyToRowMajor(
                coeffs.get(1));
            MatrixUtil.multiply(rowMajorImg, 1.f / 255.f);

            // lowering this results in more keypoints
            float ff = 0.08f;

            float[][] fastResponse = cornerFast(rowMajorImg, fastN, ff);

            float f = 0.1f;
            imageProcessor.peakLocalMax(fastResponse, 1, f, outputKP0, outputKP1);
            /*
            float factor = 255.f;
            Image img2 = img.createWithDimensions().copyToColorGreyscale();
            for (int i = 0; i < outputKP0.size(); ++i) {
                int x = outputKP1.get(i);
                int y = outputKP0.get(i);
                img2.setRGB(x, y, 255, 0, 0);
            }
            algorithms.misc.MiscDebug.writeImage(img2, "_atrous_" 
                + algorithms.misc.MiscDebug.getCurrentTimeFormatted());
            */
            log.fine("  atrous nKeypoints= " + outputKP0.size());
        }
    }
    
    /**
     * create col major image from row major input
     * @param a
     * @return
     */
    public static GreyscaleImage convertToImageGS(TwoDFloatArray a) {
        int n1 = a.a.length;
        int n2 = a.a[0].length;
        GreyscaleImage img = new GreyscaleImage(n2, n1);
        for (int i = 0; i < n1; ++i) {
            for (int j = 0; j < n2; ++j) {
                float v = 255.f * a.a[i][j];
                int vInt = Math.round(v);
                if (vInt > 255) {
                    vInt = 255;
                }
                img.setValue(j, i, vInt);
            }
        }
        return img;
    }

     /**
     * @return the pyramidImages
     */
    public TwoDFloatArray[] getPyramidImages() {
        return pyramidImages;
    }

    /**
     * @return the descrChoice
     */
    public DescriptorChoice getDescrChoice() {
        return descrChoice;
    }

}
