package algorithms.imageProcessing.features;

import algorithms.QuickSort;
import algorithms.imageProcessing.FixedSizeSortedVector;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.MedianTransform;
import algorithms.imageProcessing.StructureTensor;
import algorithms.imageProcessing.transform.MatchedPointsTransformationCalculator;
import algorithms.imageProcessing.transform.TransformationParameters;
import algorithms.imageProcessing.transform.Transformer;
import algorithms.imageProcessing.util.MatrixUtil;
import algorithms.misc.Misc;
import algorithms.misc.MiscDebug;
import algorithms.search.NearestNeighbor2D;
import algorithms.util.CorrespondencePlotter;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.QuadInt;
import algorithms.util.TwoDFloatArray;
import algorithms.util.TwoDIntArray;
import algorithms.util.VeryLongBitString;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.TFloatList;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TFloatIntMap;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.hash.TFloatIntHashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.awt.Color;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.logging.Level;
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

    // to use brute force, greedy best matching to make a correspondence list:
    int[][] matches = ORB.matchDescriptors(desc1.descriptors, desc2.descriptors);
 </pre>
 */
public class ORB {

    private static Descriptors[] getDescriptors(ORB orb, int i) {

        Descriptors[] d = null;

        if (orb.descrChoice.equals(ORB.DescriptorChoice.ALT)) {
            d = new Descriptors[]{
                orb.getDescriptorsListAlt().get(i)};
        } else if (orb.descrChoice.equals(ORB.DescriptorChoice.HSV)) {
            d = new Descriptors[]{
                orb.getDescriptorsH().get(i),
                orb.getDescriptorsS().get(i),
                orb.getDescriptorsV().get(i)};
        } else if (orb.descrChoice.equals(ORB.DescriptorChoice.GREYSCALE)) {
            d = new Descriptors[]{
                orb.getDescriptorsList().get(i)};
        }

        return d;
    }

    /*
    TODO:
       -- considering adding alternative pyramid building methods
          such as Laplacian pyramids or
          half-octave or quarter-octave pyramids (Lowe 2004; Triggs 2004)
    */

    protected boolean repeatScale1At2 = true;

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

    private List<TwoDFloatArray> pyramidImages = null;
    protected List<TwoDFloatArray> harrisResponseImages = null;

    // if HSV descriptors are requested, these are built too
    private List<TwoDFloatArray> pyramidImagesH = null;
    private List<TwoDFloatArray> pyramidImagesS = null;
    private List<TwoDFloatArray> pyramidImagesV = null;

    // alternate colorspace for descriptors created later
    private List<TwoDFloatArray> pyramidImagesAlt = null;

    protected List<TIntList> keypoints0List = null;
    protected List<TIntList> keypoints1List = null;
    protected List<TDoubleList> orientationsList = null;
    protected List<TFloatList> harrisResponses = null;
    protected List<TFloatList> scalesList = null;

    //`True`` or ``False`` representing
    //the outcome of the intensity comparison for i-th keypoint on j-th
    //decision pixel-pair. It is ``Q == np.sum(mask)``.
    private List<Descriptors> descriptorsList = null;

    private List<Descriptors> descriptorsListAlt = null;

    private List<Descriptors> descriptorsListH = null;
    private List<Descriptors> descriptorsListS = null;
    private List<Descriptors> descriptorsListV = null;

    /**
     * if method createMasks is invoked, these are populated.
     * a set bit indicates a bit outside of the segmented
     * cell of the keypoint.
     */
    private List<Descriptors> descriptorsMaskList = null;

    protected static double twoPI = 2. * Math.PI;

    protected float curvatureThresh = 0.05f;

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

    /**
     * @return the pyramidImages
     */
    public List<TwoDFloatArray> getPyramidImages() {
        return pyramidImages;
    }

    /**
     * @return the pyramidImagesAlt
     */
    public List<TwoDFloatArray> getPyramidImagesAlt() {
        return pyramidImagesAlt;
    }

    /**
     * @return the pyramidImagesH
     */
    public List<TwoDFloatArray> getPyramidImagesH() {
        return pyramidImagesH;
    }

    /**
     * @return the pyramidImagesS
     */
    public List<TwoDFloatArray> getPyramidImagesS() {
        return pyramidImagesS;
    }

    /**
     * @return the pyramidImagesV
     */
    public List<TwoDFloatArray> getPyramidImagesV() {
        return pyramidImagesV;
    }

    protected static enum DescriptorChoice {
        NONE, GREYSCALE, HSV, ALT
    }
    protected DescriptorChoice descrChoice = DescriptorChoice.GREYSCALE;

    public static enum DescriptorDithers {
        NONE, FIFTEEN, TWENTY, FORTY_FIVE, NINETY, ONE_HUNDRED_EIGHTY;
    }
    protected DescriptorDithers descrDithers = DescriptorDithers.NONE;

    protected boolean doCreate1stDerivKeypoints = false;

    protected boolean doCreateCurvatureKeyPoints = false;

    protected int nPyramidB = 3;

    protected boolean useSmallestPyramid = false;

    /**
     * Still testing the class, there may be bugs present.
     * @param nKeypoints
     */
    public ORB(int nKeypoints) {

        initMasks();

        this.nKeypoints = nKeypoints;
    }

    protected void overrideFastN(int nFast) {
        this.fastN = nFast;
    }
    public void overrideFastThreshold(float threshold) {
        this.fastThreshold = threshold;
    }

    public void overrideToUseSmallestPyramid() {
        useSmallestPyramid = true;
    }

    /**
     * override the number of divisions of scale to add in
     * between the powers of 2 scalings.  The default number
     * is 3;
     * @param n
     */
    public void overridePyamidalExtraN(int n) {
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

    public void overrideToCreateHSVDescriptors() {
        descrChoice = DescriptorChoice.HSV;
    }

    /**
     * when creating descriptors, also create descriptors at degrees of
     * dithers rotation about the calculated orientation.
     * The default setting is no offsets from the calculated orientation.
     *
     * @param dithers
     */
    public void overrideToCreateOffsetsToDescriptors(DescriptorDithers
        dithers) {
        if (dithers == null) {
            return;
        }
        descrDithers = dithers;
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

      @param image
     */
    public void detectAndExtract(Image image) {

        // extract keypoints and store results and images at class level.
        // NOTE that the keypoints are all in the coord frame of the largest image
        extractKeypoints(image);

        int nKeyPointsTotal = countKeypoints();

        System.out.println("nKeypointsTotal=" + nKeyPointsTotal +
            " this.nKeypoints=" + this.nKeypoints);

        // if descrDithers is set to other than none, this increases the
        // members of instance variables such as keypoint lists too
        replicateListsForOrientationOffsets();

        assert(pyramidImages.size() == keypoints0List.size());
        assert(pyramidImages.size() == keypoints1List.size());
        assert(pyramidImages.size() == orientationsList.size());
        assert(pyramidImages.size() == this.harrisResponses.size());
        assert(pyramidImages.size() == scalesList.size());
        if (descrChoice.equals(DescriptorChoice.HSV)) {
            assert(pyramidImages.size() == pyramidImagesH.size());
            assert(pyramidImages.size() == pyramidImagesS.size());
            assert(pyramidImages.size() == pyramidImagesV.size());
        }

        boolean useDefaultSize = true;
        extractDescriptors(useDefaultSize);
    }

    void createDescriptorsHSV(Image image) {
        createDescriptorsHSV(image, true);
    }
    void createSmallDescriptorsHSV(Image image) {
        createDescriptorsHSV(image, false);
    }

    void createSmallDescriptorsLABTheta(Image image) {
        createDescriptorsLABTheta(image, false);
    }
    void createDescriptorsLABTheta(Image image) {
        createDescriptorsLABTheta(image, true);
    }

    private void createDescriptorsHSV(Image image, boolean useDefaultSize) {

        if (descriptorsListH != null) {
            throw new IllegalStateException("hsv descriptors are already built");
        }

        this.descrChoice = DescriptorChoice.HSV;

        pyramidImagesH = new ArrayList<TwoDFloatArray>();
        pyramidImagesS = new ArrayList<TwoDFloatArray>();
        pyramidImagesV = new ArrayList<TwoDFloatArray>();
        buildHSVPyramid(image, pyramidImagesH, pyramidImagesS,
            pyramidImagesV);

        // -- may need to remove some of the images if keypoints filtering
        //    has removed the data for a scale.
        if (pyramidImages.size() < pyramidImagesH.size()) {
            // use scalesList and octaveScaleMap to remove images

            List<TwoDFloatArray> h = new ArrayList<TwoDFloatArray>();
            List<TwoDFloatArray> s = new ArrayList<TwoDFloatArray>();
            List<TwoDFloatArray> v = new ArrayList<TwoDFloatArray>();

            for (int i = 0; i < scalesList.size(); ++i) {
                float scale = scalesList.get(i).get(0);
                assert(octaveScaleMap.containsKey(scale));
                int octave = octaveScaleMap.get(scale);
                h.add(pyramidImagesH.get(octave));
                s.add(pyramidImagesS.get(octave));
                v.add(pyramidImagesV.get(octave));
            }

            assert(pyramidImages.size() == h.size());
            assert(pyramidImages.size() == s.size());
            assert(pyramidImages.size() == v.size());

            pyramidImagesH.clear();
            pyramidImagesS.clear();
            pyramidImagesV.clear();

            pyramidImagesH.addAll(h);
            pyramidImagesS.addAll(s);
            pyramidImagesV.addAll(v);
        }

        assert(pyramidImages.size() == keypoints0List.size());
        assert(pyramidImages.size() == keypoints1List.size());
        assert(pyramidImages.size() == pyramidImagesH.size());
        assert(pyramidImages.size() == pyramidImagesS.size());
        assert(pyramidImages.size() == pyramidImagesV.size());

        extractDescriptors(useDefaultSize);
    }

    private void createDescriptorsLABTheta(Image image,
        boolean useDefaultSize) {

        if (descriptorsListAlt != null) {
            throw new IllegalStateException("alt descriptors are already built");
        }

        this.descrChoice = DescriptorChoice.ALT;

        pyramidImagesAlt = new ArrayList<TwoDFloatArray>();

        buildLABThetaPyramid(image, pyramidImagesAlt);

        // -- may need to remove some of the images if keypoints filtering
        //    has removed the data for a scale.
        if (pyramidImages.size() < pyramidImagesAlt.size()) {
            // use scalesList and octaveScaleMap to remove images

            List<TwoDFloatArray> alt = new ArrayList<TwoDFloatArray>();

            for (int i = 0; i < scalesList.size(); ++i) {
                float scale = scalesList.get(i).get(0);
                assert(octaveScaleMap.containsKey(scale));
                int octave = octaveScaleMap.get(scale);
                alt.add(pyramidImagesAlt.get(octave));
            }

            assert(pyramidImages.size() == alt.size());

            pyramidImagesAlt.clear();

            pyramidImagesAlt.addAll(alt);
        }

        assert(pyramidImages.size() == keypoints0List.size());
        assert(pyramidImages.size() == keypoints1List.size());
        assert(pyramidImages.size() == pyramidImagesAlt.size());

        extractDescriptors(useDefaultSize);
    }

    /**
     * @param image
     * @return
     */
    private List<TwoDFloatArray> buildPyramid(Image image) {

        GreyscaleImage img = image.copyToGreyscale2();

        return buildPyramid(img);
    }

    /**
     *
     * @param image
     * @param pyramidH empty list in which to put resulting H pyramid
     * @param pyramidS empty list in which to put resulting S pyramid
     * @param pyramidV empty list in which to put resulting V pyramid
     */
    private void buildHSVPyramid(Image image, List<TwoDFloatArray> pyramidH,
        List<TwoDFloatArray> pyramidS, List<TwoDFloatArray> pyramidV) {

        GreyscaleImage imageR = image.copyRedToGreyscale();
        GreyscaleImage imageG = image.copyGreenToGreyscale();
        GreyscaleImage imageB = image.copyBlueToGreyscale();

        List<GreyscaleImage> outputR;
        List<GreyscaleImage> outputG;
        List<GreyscaleImage> outputB;

        if (useSmallestPyramid) {

            outputR = new ArrayList<GreyscaleImage>();
            outputG = new ArrayList<GreyscaleImage>();
            outputB = new ArrayList<GreyscaleImage>();

            MedianTransform mt = new MedianTransform();
            mt.multiscalePyramidalMedianTransform2(imageR, outputR, decimationLimit);
            mt.multiscalePyramidalMedianTransform2(imageG, outputG, decimationLimit);
            mt.multiscalePyramidalMedianTransform2(imageB, outputB, decimationLimit);

        } else {

            ImageProcessor imageProcessor = new ImageProcessor();

            outputR = imageProcessor.buildPyramid2(
                imageR, decimationLimit, nPyramidB);

            outputG = imageProcessor.buildPyramid2(
                imageG, decimationLimit, nPyramidB);

            outputB = imageProcessor.buildPyramid2(
                imageB, decimationLimit, nPyramidB);
        }

        float[] hsv = new float[3];

        for (int i = 0; i < outputR.size(); ++i) {
            GreyscaleImage imgR = outputR.get(i);
            GreyscaleImage imgG = outputG.get(i);
            GreyscaleImage imgB = outputB.get(i);

            int nRows = imgR.getHeight();
            int nCols = imgR.getWidth();
            float[][] outH = new float[nRows][nCols];
            float[][] outS = new float[nRows][nCols];
            float[][] outV = new float[nRows][nCols];
            for (int j = 0; j < nRows; ++j) {
                outH[j] = new float[nCols];
                outS[j] = new float[nCols];
                outV[j] = new float[nCols];
                for (int k = 0; k < nCols; ++k) {
                    int r = imgR.getValue(k, j);
                    int g = imgG.getValue(k, j);
                    int b = imgB.getValue(k, j);
                    Color.RGBtoHSB(r, g, b, hsv);
                    outH[j][k] = hsv[0];
                    outS[j][k] = hsv[1];
                    outV[j][k] = hsv[2];
                }
            }

            TwoDFloatArray h = new TwoDFloatArray(outH);
            TwoDFloatArray s = new TwoDFloatArray(outS);
            TwoDFloatArray v = new TwoDFloatArray(outV);

            pyramidH.add(h);
            pyramidS.add(s);
            pyramidV.add(v);
        }
    }

    /**
     *
     * @param image
     * @param pyramidAlt empty list in which to put resulting H pyramid
     */
    private void buildLABThetaPyramid(Image image,
        List<TwoDFloatArray> pyramidAlt) {

        ImageProcessor imageProcessor = new ImageProcessor();

        GreyscaleImage imageAlt = imageProcessor.
            createCIELABTheta(image, 255);

        List<GreyscaleImage> outputAlt;

        if (useSmallestPyramid) {

            outputAlt = new ArrayList<GreyscaleImage>();

            MedianTransform mt = new MedianTransform();
            mt.multiscalePyramidalMedianTransform2(imageAlt,
                outputAlt, decimationLimit);

        } else {

            outputAlt = imageProcessor.buildPyramid2(
                imageAlt, decimationLimit, nPyramidB);
        }

        for (int i = 0; i < outputAlt.size(); ++i) {
            GreyscaleImage imgAlt = outputAlt.get(i);

            int nRows = imgAlt.getHeight();
            int nCols = imgAlt.getWidth();

            float[][] rowMajorImg = imageProcessor.copyToRowMajor(
                imgAlt);
            MatrixUtil.multiply(rowMajorImg, 1.f/255.f);
            TwoDFloatArray f = new TwoDFloatArray(rowMajorImg);
            pyramidAlt.add(f);
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

        if (useSmallestPyramid) {
            output = new ArrayList<GreyscaleImage>();
            MedianTransform mt = new MedianTransform();
            mt.multiscalePyramidalMedianTransform2(img, output, decimationLimit);
        } else {
           output = imageProcessor.buildPyramid2(
               img, decimationLimit, nPyramidB);
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

    /**
     * create a greyscale image from the color image
     * and return it as two dimensional float array, scaled to
     * range [0.0, 1.0] and placed in format img[row][col].
     * @param image
     */
    private float[][] prepareGrayscaleInput2D(ImageExt image) {

        GreyscaleImage img = image.copyToGreyscale2();

        ImageProcessor imageProcessor = new ImageProcessor();
        return imageProcessor.multiply(img, 1.f/255.f);
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
    private void extractKeypoints(Image image) {

        List<TwoDFloatArray> pyramid = buildPyramid(image);

        assert(pyramid.get(0).a.length == image.getHeight());
        assert(pyramid.get(0).a[0].length == image.getWidth());

        keypoints0List = new ArrayList<TIntList>();
        keypoints1List = new ArrayList<TIntList>();
        orientationsList = new ArrayList<TDoubleList>();
        harrisResponses = new ArrayList<TFloatList>();
        scalesList = new ArrayList<TFloatList>();
        pyramidImages = new ArrayList<TwoDFloatArray>();
        harrisResponseImages = new ArrayList<TwoDFloatArray>();

        List<TwoDFloatArray> pyramidH = null;
        List<TwoDFloatArray> pyramidS = null;
        List<TwoDFloatArray> pyramidV = null;
        if (descrChoice.equals(DescriptorChoice.HSV)) {
            pyramidImagesH = new ArrayList<TwoDFloatArray>();
            pyramidImagesS = new ArrayList<TwoDFloatArray>();
            pyramidImagesV = new ArrayList<TwoDFloatArray>();
            pyramidH = new ArrayList<TwoDFloatArray>();
            pyramidS = new ArrayList<TwoDFloatArray>();
            pyramidV = new ArrayList<TwoDFloatArray>();
            buildHSVPyramid(image, pyramidH, pyramidS, pyramidV);
            assert(pyramidH.size() == pyramid.size());
            assert(pyramidH.get(0).a.length == image.getHeight());
            assert(pyramidH.get(0).a[0].length == image.getWidth());
        }

        float prevScl = 1;

        TIntList octavesUsed = new TIntArrayList();
        TFloatList octavesScales = new TFloatArrayList();
        List<TwoDFloatArray> hrList = new ArrayList<TwoDFloatArray>();
        Set<PairInt> kp01CombinedSet = new HashSet<PairInt>();
        for (int octave = 0; octave < pyramid.size(); ++octave) {

            //System.out.println("octave=" + octave);

            float[][] octaveImage = pyramid.get(octave).a;

            float scale = prevScl;
            {
                if (octave > 0) {
                    float x0 = (float)pyramid.get(octave - 1).a.length/
                        (float)pyramid.get(octave).a.length;
                    float y0 = (float)pyramid.get(octave - 1).a[0].length/
                        (float)pyramid.get(octave).a[0].length;
                    scale = prevScl * (x0 + y0)/2.f;
                    System.out.println("scl=" + scale + " prevScl=" + prevScl);
                    prevScl = scale;
                }
            }

            if (octaveImage.length < decimationLimit ||
                octaveImage[0].length < decimationLimit) {
                continue;
            }

            Resp r = extractFastKeypointsAndHarrisImage(octaveImage, scale);

            if ((r == null) || (r.keypoints0 == null) || r.keypoints0.isEmpty()) {
                continue;
            }

            System.out.println("  octave " + octave + " nKeypoints="
                + r.keypoints0.size());

            //mask of length orientations.size() containing a 1 or 0
            // indicating if pixels are within the image (``True``) or in the
            // border region of the image (``False``).
            TIntList mask = filterNearBorder(octaveImage, r.keypoints0,
                r.keypoints1);

            // mult by scale and store locally:
            for (int k = 0; k < r.keypoints0.size(); ++k) {
                int c0 = Math.round(scale * r.keypoints0.get(k));
                int c1 = Math.round(scale * r.keypoints1.get(k));
                PairInt p = new PairInt(c0, c1);
                kp01CombinedSet.add(p);
            }

            hrList.add(new TwoDFloatArray(r.harrisResponse));
            octavesUsed.add(octave);
            octavesScales.add(scale);
            octaveScaleMap.put(scale, octave);
        }

        TIntList kpc0s = new TIntArrayList(this.nKeypoints);
        TIntList kpc1s = new TIntArrayList(this.nKeypoints);
        // order by response and filter
        {
            float[][] harrisResponse = hrList.get(0).a;

            // kp01CombinedSet are in coord frame of largest image in pyramid, [0]
            int n = kp01CombinedSet.size();

            float[] scores = new float[n];
            PairInt[] points = new PairInt[n];
            int count = 0;
            for (PairInt p : kp01CombinedSet) {
                points[count] = p;
                //NOTE: because scale[0] is the full size frame, these coords
                // do not need to be rescaled
                scores[count] = harrisResponse[p.getX()][p.getY()];
                count++;
            }

            QuickSort.sortBy1stArg(scores, points);

            kp01CombinedSet.clear();

            for (int i = (points.length - 1); i > -1; --i) {
                PairInt p = points[i];
                int x = p.getX();
                int y = p.getY();
                kpc0s.add(x);
                kpc1s.add(y);
            }

            // then filter
            maskCoordinates(kpc0s, kpc1s, harrisResponse.length,
                harrisResponse[0].length, 4);//8);
        }

        // at this point kpc0s, kpc1s are in coord frame of largest image

        // filter for number of points if requested
        if (kpc0s.size() > this.nKeypoints) {
            while (kpc0s.size() > this.nKeypoints) {
                int idx = kpc0s.size() - 1;
                kpc0s.removeAt(idx);
                kpc1s.removeAt(idx);
            }
        }

        for (int ij = 0; ij < octavesUsed.size(); ++ij) {

            int octave = octavesUsed.get(ij);
            float scale = octavesScales.get(ij);
            float[][] harrisResponse = hrList.get(ij).a;

            // make a keypoint list for each scale, using set first
            // to remove redundant points.
            // TODO: may need to apply a local max filter if points cluster in
            //    smaller image
            Set<PairInt> set = new HashSet<PairInt>();
            TIntList kp0s = new TIntArrayList(set.size());
            TIntList kp1s = new TIntArrayList(set.size());
            for (int i = 0; i < kpc0s.size(); ++i) {
                // reduce to octave coord system
                int c0 = Math.round(kpc0s.get(i)/scale);
                int c1 = Math.round(kpc1s.get(i)/scale);
                if (c0 > (harrisResponse.length - 1)) {
                    c0 = harrisResponse.length - 1;
                }
                if (c1 > (harrisResponse[0].length - 1)) {
                    c1 = harrisResponse[0].length - 1;
                }
                PairInt p = new PairInt(c0, c1);
                if (!set.contains(p)) {
                    set.add(p);
                }
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

            float[][] octaveImage = pyramid.get(octave).a;

            TDoubleList orientations = cornerOrientations(octaveImage,
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

            pyramidImages.add(pyramid.get(octave));
            harrisResponseImages.add(new TwoDFloatArray(harrisResponse));

            TFloatList scales = new TFloatArrayList(kp0s.size());
            for (int i = 0; i < kp0s.size(); ++i) {
                scales.add(scale);
            }
            scalesList.add(scales);

            if (pyramidH != null) {
                pyramidImagesH.add(pyramidH.get(octave));
                pyramidImagesS.add(pyramidS.get(octave));
                pyramidImagesV.add(pyramidV.get(octave));
            }
        }

        assert(pyramidImages.size() == keypoints0List.size());
        assert(pyramidImages.size() == keypoints1List.size());
        assert(pyramidImages.size() == orientationsList.size());
        assert(pyramidImages.size() == this.harrisResponses.size());
        assert(pyramidImages.size() == scalesList.size());
        if (descrChoice.equals(DescriptorChoice.HSV)) {
            assert(pyramidImages.size() == pyramidImagesH.size());
            assert(pyramidImages.size() == pyramidImagesS.size());
            assert(pyramidImages.size() == pyramidImagesV.size());
        }
    }

    /**
     * if descrDithers is not none, replicates lists such as keypoints and
     * then copies and modifies orientation for offsets setting.
     */
    private void replicateListsForOrientationOffsets() {

        if (descrDithers.equals(DescriptorDithers.NONE)) {
             return;
        }

        for (int i = 0; i < pyramidImages.size(); ++i) {

            TFloatList scales = new TFloatArrayList(this.scalesList.get(i));
            TIntList kp0 = new TIntArrayList(this.keypoints0List.get(i));
            TIntList kp1 = new TIntArrayList(this.keypoints1List.get(i));
            TDoubleList or = new TDoubleArrayList(this.orientationsList.get(i));
            TFloatList harrisResp = new TFloatArrayList(this.harrisResponses.get(i));

            int nBefore = kp0.size();
            int nExpected = kp0.size();

            double[] rotations = null;
            if (descrDithers.equals(DescriptorDithers.FIFTEEN)) {
                rotations = new double[]{
                    Math.PI/12., Math.PI/6., 3.*Math.PI/12.,
                    4.*Math.PI/12., 5.*Math.PI/12., 6.*Math.PI/12.,
                    7.*Math.PI/12., 8.*Math.PI/12., 9.*Math.PI/12.,
                    10.*Math.PI/12., 11.*Math.PI/12., Math.PI,
                    13.*Math.PI/12., 14.*Math.PI/12., 15.*Math.PI/12.,
                    16.*Math.PI/12., 17.*Math.PI/12., 18.*Math.PI/12.,
                    19.*Math.PI/12., 20.*Math.PI/12., 21.*Math.PI/12.,
                    22.*Math.PI/12., 23.*Math.PI/12.};
                nExpected *= 24;
            } else if (descrDithers.equals(DescriptorDithers.TWENTY)) {
                rotations = new double[]{
                    Math.PI/9., 2.*Math.PI/9., 3.*Math.PI/9.,
                    4.*Math.PI/9., 5.*Math.PI/9., 6.*Math.PI/9.,
                    7.*Math.PI/9., 8.*Math.PI/9., Math.PI,
                    10.*Math.PI/9., 11.*Math.PI/9., 12.*Math.PI/9.,
                    13.*Math.PI/9., 14.*Math.PI/9., 15.*Math.PI/9.,
                    16.*Math.PI/9., 17.*Math.PI/9.
                };
                nExpected *= 18;
            } else if (descrDithers.equals(DescriptorDithers.FORTY_FIVE)) {
                rotations = new double[]{Math.PI/4., Math.PI/2.,
                    3.*Math.PI/4., Math.PI, 5.*Math.PI/4.,
                    6.*Math.PI/4., 7.*Math.PI/4.};
                nExpected *= 8;
            } else if (descrDithers.equals(DescriptorDithers.NINETY)) {
                rotations = new double[]{Math.PI/2., Math.PI, 3.*Math.PI/2.};
                nExpected *= 4;
            } else if (descrDithers.equals(DescriptorDithers.ONE_HUNDRED_EIGHTY)) {
                rotations = new double[]{Math.PI};
                nExpected *= 2;
            }

            for (double rotation : rotations) {

                scalesList.get(i).addAll(scales);
                keypoints0List.get(i).addAll(kp0);
                keypoints1List.get(i).addAll(kp1);
                harrisResponses.get(i).addAll(harrisResp);

                TDoubleList orientation = new TDoubleArrayList(or);

                /* orientation needs to stay in this range:
                    +90
                135  |  +45
                     |
               180---.---- 0
                     |
               -135  |   -45
                    -90
                */
                for (int ii = 0; ii < orientation.size(); ++ii) {
                    double d = orientation.get(ii);
                    d += rotation;
                    if (d > Math.PI) {
                        d -= twoPI;
                    } else if (d <= -twoPI) {
                        d += twoPI;
                    }
                    assert(d > -180. && d <= 180.);
                    orientation.set(ii, d);
                }

                orientationsList.get(i).addAll(orientation);
            }

            assert(nExpected == scalesList.get(i).size());
            assert(nExpected == keypoints0List.get(i).size());
            assert(nExpected == keypoints1List.get(i).size());
            assert(nExpected == harrisResponses.get(i).size());
            assert(nExpected == orientationsList.get(i).size());
        }
    }

    private void extractDescriptors(boolean useDefaultSize) {

        if (descrChoice.equals(DescriptorChoice.NONE)) {
            return;
        }

        if (descrChoice.equals(DescriptorChoice.GREYSCALE)) {
            descriptorsList = extractGreyscaleDescriptors(pyramidImages, useDefaultSize);
            return;
        }

        if (descrChoice.equals(DescriptorChoice.HSV)) {
            descriptorsListH = extractGreyscaleDescriptors(pyramidImagesH, useDefaultSize);
            descriptorsListS = extractGreyscaleDescriptors(pyramidImagesS, useDefaultSize);
            descriptorsListV = extractGreyscaleDescriptors(pyramidImagesV, useDefaultSize);
            return;
        }

        if (descrChoice.equals(DescriptorChoice.ALT)) {
            descriptorsListAlt = extractGreyscaleDescriptors(
                pyramidImagesAlt, useDefaultSize);
            return;
        }
    }

    protected List<Descriptors> extractGreyscaleDescriptors(List<TwoDFloatArray>
        pyramid) {

        boolean useDefaultSize = true;
        return extractGreyscaleDescriptors(pyramid, useDefaultSize);
    }

    private List<Descriptors> extractGreyscaleDescriptors(List<TwoDFloatArray>
        pyramid, boolean useDefaultSize) {

        assert(pyramid.size() == this.keypoints0List.size());
        assert(pyramid.size() == this.keypoints1List.size());
        assert(pyramid.size() == this.scalesList.size());
        assert(pyramid.size() == this.orientationsList.size());

        List<Descriptors> output = new ArrayList<Descriptors>();

        for (int i = 0; i < pyramid.size(); ++i) {

            float[][] octaveImage = pyramid.get(i).a;

            TFloatList scales = this.scalesList.get(i);
            TIntList kp0 = this.keypoints0List.get(i);
            TIntList kp1 = this.keypoints1List.get(i);
            TDoubleList or = this.orientationsList.get(i);

            // result contains descriptors and mask.
            // also, modified by mask are the keypoints and orientations
            Descriptors desc = extractOctave(octaveImage, kp0, kp1, or,
                useDefaultSize, scales.get(0));

            output.add(desc);
        }

        return output;
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

    private class Resp {
        TIntList keypoints0;
        TIntList keypoints1;
        float[][] harrisResponse;
    }

    public static class Descriptors {
        //NOTE: consider packing 4 descriptors into one or
        // doing as Image.java does with sensing 64 bit and 32 bit to make
        // long or int bit vectors
        VeryLongBitString[] descriptors;
    }

    /**
     * adapted from
     * https://github.com/scikit-image/scikit-image/blob/master/skimage/feature/orb.py
     *
     * @param octaveImage
     * @return
     */
    private Resp extractFastKeypointsAndHarrisImage(float[][] octaveImage,
        float scale) {

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

        // Standard deviation used for the Gaussian kernel, which is used as
        // weighting function for the auto-correlation matrix.
        float sigma = 1;

        StructureTensor tensorComponents = new StructureTensor(octaveImage,
            sigma, doCreateCurvatureKeyPoints);

        if (doCreate1stDerivKeypoints) {

            ImageProcessor imageProcessor = new ImageProcessor();

            float hLimit = 0.01f;//0.09f;

            imageProcessor.createFirstDerivKeyPoints(
                tensorComponents, keypoints0, keypoints1, hLimit);

            /*
            try {
                float factor = 255.f;
                Image img2 = new Image(nRows, nCols);
                for (int i = 0; i < nRows; ++i) {
                    for (int j = 0; j < nCols; ++j) {
                        int v = Math.round(factor * octaveImage[i][j]);
                        if (v > 255) {
                            v = 255;
                        }
                        img2.setRGB(i, j, v, v, v);
                    }
                }
                for (int i = 0; i < keypoints0.size(); ++i) {
                    int y = keypoints1.get(i);
                    int x = keypoints0.get(i);
                    img2.setRGB(x, y, 255, 0, 0);
                }
                System.out.println("nRows=" + nRows + " nCols=" + nCols);
                algorithms.imageProcessing.ImageDisplayer.displayImage("first deriv", img2);
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
            imageProcessor.createCurvatureKeyPoints(
                tensorComponents, keypoints0, keypoints1,
                curvatureThresh);

            /*try {
                float factor = 255.f;
                Image img2 = new Image(nRows, nCols);
                for (int i = 0; i < nRows; ++i) {
                    for (int j = 0; j < nCols; ++j) {
                        int v = Math.round(factor * octaveImage[i][j]);
                        if (v > 255) {
                            v = 255;
                        }
                        img2.setRGB(i, j, v, v, v);
                    }
                }
                for (int i = 0; i < keypoints0.size(); ++i) {
                    int y = keypoints1.get(i);
                    int x = keypoints0.get(i);
                    img2.setRGB(x, y, 255, 0, 0);
                }
                System.out.println("nRows=" + nRows + " nCols=" + nCols);
                algorithms.imageProcessing.ImageDisplayer.displayImage("curvature", img2);
                int z = 1;
            } catch(Exception e) {
                System.out.println(e.getMessage());
            }*/
        }

        float[][] detA = tensorComponents.getDeterminant();

        float[][] traceA = tensorComponents.getTrace();

        // size is same a octaveImage
        float[][] harrisResponse = cornerHarris(octaveImage, detA, traceA);

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
        r2.harrisResponse = harrisResponse;

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
    protected TDoubleList cornerOrientations(float[][] octaveImage,
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

        for (i = 0; i < keypoints0.size(); ++i) {
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
     * @param responses
     * @return the encapsulated descriptors and mask
     */
    protected Descriptors extractOctave(float[][] octaveImage,
        TIntList keypoints0, TIntList keypoints1,
        TDoubleList orientations, float scale) {

        boolean useDefaultSize = true;

        return extractOctave(octaveImage, keypoints0, keypoints1, orientations,
            useDefaultSize, scale);
    }

    /**
     * filter the keypoints, orientations, and responses to remove those close
     * to the border and then create descriptors for the remaining.
     *
     * @param octaveImage
     * @param keypoints0
     * @param keypoints1
     * @param orientations
     * @param responses
     * @return the encapsulated descriptors and mask
     */
    protected Descriptors extractOctave(float[][] octaveImage,
        TIntList keypoints0, TIntList keypoints1,
        TDoubleList orientations, boolean useDefaultSize,
        float scale) {

        if (POS0 == null) {
            if (useDefaultSize) {
                POS0 = ORBDescriptorPositions.POS0;
            } else {
                POS0 = ORBSmallDescriptorPositions.POS3;
            }
        }
        if (POS1 == null) {
            if (useDefaultSize) {
                POS1 = ORBDescriptorPositions.POS1;
            } else {
                POS1 = ORBSmallDescriptorPositions.POS4;
            }
        }

        assert(orientations.size() == keypoints0.size());
        assert(orientations.size() == keypoints1.size());

        VeryLongBitString[] descriptors = null;

        if (descrChoice.equals(DescriptorChoice.NONE)) {
            descriptors = new VeryLongBitString[0];
        } else {
            descriptors = orbLoop(octaveImage, keypoints0, keypoints1,
                orientations, useDefaultSize, scale);
        }

        Descriptors desc = new Descriptors();
        desc.descriptors = descriptors;

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
    protected VeryLongBitString[] orbLoop(float[][] octaveImage, TIntList keypoints0,
        TIntList keypoints1, TDoubleList orientations, float scale) {

        boolean useDefaultSize = true;

        return orbLoop(octaveImage, keypoints0, keypoints1, orientations,
            useDefaultSize, scale);
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
        TIntList keypoints1, TDoubleList orientations, boolean useDefaultSize,
        float scale) {

        assert(orientations.size() == keypoints0.size());

        if (POS0 == null) {
            if (useDefaultSize) {
                POS0 = ORBDescriptorPositions.POS0;
            } else {
                POS0 = ORBSmallDescriptorPositions.POS3;
            }
        }
        if (POS1 == null) {
            if (useDefaultSize) {
                POS1 = ORBDescriptorPositions.POS1;
            } else {
                POS1 = ORBSmallDescriptorPositions.POS4;
            }
        }

        int nKP = orientations.size();

        System.out.println("nKP=" + nKP);

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
     * create masks for the descriptors based upon the given segmentation
     * pointIndexMap.
     * descriptor pixels that are not in the same segmented cell as the keypoint
     * for the descriptor are set bits in the descriptor mask.
     * @param pointIndexMap map of point segmented cell indexes.  key of map
     * is the x,y coordinates (col major data structure compatible)
     * and value is the index of the segmented cell.  if a point is not in the
     * map, it is not in a segmented cell.  points in the same cell will have
     * the same value.
     */
    public void createDescriptorMasks(TObjectIntMap<PairInt> pointIndexMap) {

        if (descrChoice.equals(DescriptorChoice.GREYSCALE)) {
            if (descriptorsList == null) {
                throw new IllegalStateException("descriptors must be created first");
            }
        } else if (descrChoice.equals(DescriptorChoice.HSV)) {
            if (descriptorsListH == null) {
                throw new IllegalStateException("descriptors must be created first");
            }
        } else if (descrChoice.equals(DescriptorChoice.ALT)) {
            if (descriptorsListAlt == null) {
                throw new IllegalStateException("descriptors must be created first");
            }
        }

        if (descriptorsMaskList != null) {
            throw new IllegalStateException("a mask for the descriptors has "
                + " already been created");
        }

        assert(!descrChoice.equals(DescriptorChoice.NONE));
        assert(POS0 != null);
        assert(POS1 != null);
        assert(orientationsList.size() == keypoints0List.size());
        assert(keypoints0List.size() == keypoints1List.size());
        assert(pyramidImages.size() == keypoints0List.size());

        int ns = keypoints0List.size();
        descriptorsMaskList = new ArrayList<Descriptors>();

        for (int i = 0; i < ns; ++i) {
            assert(!scalesList.get(i).isEmpty());
            float scale = scalesList.get(i).get(0);
            //TwoDFloatArray octaveImage = pyramidImages.get(i);
            TDoubleList or = orientationsList.get(i);
            TIntList kp0 = keypoints0List.get(i);
            TIntList kp1 = keypoints1List.get(i);
            Descriptors mask = createMask(kp0, kp1, or,
                pointIndexMap, scale);

            descriptorsMaskList.add(mask);
        }
    }

    /**
     *
     * @param octaveImage
     * @param keypoints0 keypoints coordinate for row.  note that they are
     * assumed to be already scaled up to full frame size, so scale is used
     * to reduce them to reference frame of octave image.
     * @param keypoints1 keypoints coordinates for col.  see scale notes in
     * keypoints0
     * @param orientations
     * @param pointIndexMap
     * @param scale
     * @return
     */
    protected Descriptors createMask(TIntList keypoints0,
        TIntList keypoints1, TDoubleList orientations,
        TObjectIntMap<PairInt> pointIndexMap, float scale) {

        assert(orientations.size() == keypoints0.size());
        assert(POS0 != null);
        assert(POS1 != null);

        int nKP = orientations.size();

        System.out.println("nKP=" + nKP);

        VeryLongBitString[] descriptors = new VeryLongBitString[nKP];
        Descriptors desc = new Descriptors();
        desc.descriptors = descriptors;

        int r = 256;
        if (this.descriptorsList != null) {
            r = (int)this.descriptorsList.get(0).descriptors[0]
                .getInstantiatedBitSize();
        } else if (this.descriptorsListH != null) {
            r = (int)this.descriptorsListH.get(0).descriptors[0]
                .getInstantiatedBitSize();
        }
        assert(r == POS0.length);

        boolean useDefaultSize = (r < 256);

        double pr0, pc0, pr1, pc1;
        int spr0, spc0, spr1, spc1;

        for (int i = 0; i < descriptors.length; ++i) {

            descriptors[i] = new VeryLongBitString(r);

            double angle = orientations.get(i);
            double sinA = Math.sin(angle);
            double cosA = Math.cos(angle);

            int kr = keypoints0.get(i);
            int kc = keypoints1.get(i);

            PairInt p = new PairInt(Math.round(kc), Math.round(kr));
            int label = pointIndexMap.get(p);
            assert(pointIndexMap.containsKey(p));

            // put kr and kc into pyrmid image reference frame.
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

                // x0,y0 is row major notation and pointIndexMap is col major.
                // scale up the coordinates to full frame for use with map.
                PairInt p0 = new PairInt(Math.round(y0*scale), Math.round(x0*scale));
                PairInt p1 = new PairInt(Math.round(y1*scale), Math.round(x1*scale));

                if (!pointIndexMap.containsKey(p0)
                    || pointIndexMap.get(p0) != label
                    || !pointIndexMap.containsKey(p1)
                    || pointIndexMap.get(p1) != label) {

                    descriptors[i].setBit(j);
                }
            }
        }

        return desc;
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

    public List<Descriptors> getDescriptorsListAlt() {
        return descriptorsListAlt;
    }

    /**
     * get a list of each octave's descriptors mask.  NOTE that the list
     * is not copied so do not modify.
     * The coordinates of the descriptors mask can be found in keyPointsList, but
     * with twice the spacing because that stores row and col in same list.
     * @return
     */
    public List<Descriptors> getDescriptorsMaskList() {
        return descriptorsMaskList;
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

    public List<Descriptors> getDescriptorsH() {
        return descriptorsListH;
    }
    public List<Descriptors> getDescriptorsS() {
        return descriptorsListS;
    }
    public List<Descriptors> getDescriptorsV() {
        return descriptorsListV;
    }

    /**
     * get a list of each octave's descriptors as a combined descriptor.
     * The coordinates of the descriptors can be found in getAllKeypoints, but
     * with twice the spacing because that stores row and col in same list.
     * @return
     */
    public Descriptors[] getAllDescriptorsHSV() {

        Descriptors descrH = combineDescriptors(descriptorsListH);
        Descriptors descrS = combineDescriptors(descriptorsListS);
        Descriptors descrV = combineDescriptors(descriptorsListV);

        return new Descriptors[]{descrH, descrS, descrV};
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
     * remove items at indexes in list.
     * Note that a final step of removing all data at a scale size
     * is performed if there are no keypoints left for that scale.
     * @param rmIndexesList
     */
    public void removeAtIndexes(List<TIntList> rmIndexesList) {

        int ns = keypoints0List.size();

        if (rmIndexesList.size() != ns) {
            throw new IllegalArgumentException("indexes list size must "
                + "be the same size as internal lists");
        }

        for (int i = 0; i < ns; ++i) {
            TIntList kp0 = keypoints0List.get(i);
            TIntList kp1 = keypoints1List.get(i);
            TDoubleList or = orientationsList.get(i);
            TFloatList s = scalesList.get(i);
            TFloatList r = harrisResponses.get(i);
            Descriptors d = null;
            Descriptors dAlt = null;
            Descriptors dH = null;
            Descriptors dS = null;
            Descriptors dV = null;
            if (descrChoice.equals(DescriptorChoice.HSV)) {
                dH = descriptorsListH.get(i);
                dS = descriptorsListS.get(i);
                dV = descriptorsListV.get(i);
            } else if (descrChoice.equals(DescriptorChoice.GREYSCALE)) {
                d = descriptorsList.get(i);
            } else if (descrChoice.equals(DescriptorChoice.ALT)) {
                dAlt = descriptorsListAlt.get(i);
            }
            Descriptors m = null;
            if (descriptorsMaskList != null) {
                m = descriptorsMaskList.get(i);
            }

            int np = kp0.size();
            TIntList rm = rmIndexesList.get(i);

            if (!rm.isEmpty()) {
                int nb = np - rm.size();

                Descriptors d2 = null;
                Descriptors dAlt2 = null;
                Descriptors dH2 = null;
                Descriptors dS2 = null;
                Descriptors dV2 = null;

                if (d != null) {
                    d2 = new Descriptors();
                    d2.descriptors = new VeryLongBitString[nb];
                }
                if (dAlt != null) {
                    dAlt2 = new Descriptors();
                    dAlt2.descriptors = new VeryLongBitString[nb];
                }
                if (dH != null) {
                    dH2 = new Descriptors();
                    dH2.descriptors = new VeryLongBitString[nb];
                    dS2 = new Descriptors();
                    dS2.descriptors = new VeryLongBitString[nb];
                    dV2 = new Descriptors();
                    dV2.descriptors = new VeryLongBitString[nb];
                }

                Descriptors m2 = null;
                if (m != null) {
                    m2 = new Descriptors();
                    m2.descriptors = new VeryLongBitString[nb];
                }

                TIntSet rmSet = new TIntHashSet(rm);
                for (int j = (rm.size() - 1); j > -1; --j) {
                    int idx = rm.get(j);
                    kp0.removeAt(idx);
                    kp1.removeAt(idx);
                    or.removeAt(idx);
                    s.removeAt(idx);
                    r.removeAt(idx);
                }
                int count = 0;
                for (int j = 0; j < np; ++j) {
                    if (rmSet.contains(j)) {
                        continue;
                    }
                    if (d2 != null) {
                        d2.descriptors[count] = d.descriptors[j];
                    }
                    if (dAlt2 != null) {
                        dAlt2.descriptors[count] = dAlt.descriptors[j];
                    }
                    if (dH2 != null) {
                        dH2.descriptors[count] = dH.descriptors[j];
                        dS2.descriptors[count] = dS.descriptors[j];
                        dV2.descriptors[count] = dV.descriptors[j];
                    }
                    if (m2 != null) {
                        m2.descriptors[count] = m.descriptors[j];
                    }

                    count++;
                }
                assert(count == nb);
                if (d2 != null) {
                    d.descriptors = d2.descriptors;
                }
                if (dAlt2 != null) {
                    dAlt.descriptors = dAlt2.descriptors;
                }
                if (dH2 != null) {
                    dH.descriptors = dH2.descriptors;
                    dS.descriptors = dS2.descriptors;
                    dV.descriptors = dV2.descriptors;
                }
                if (m2 != null) {
                    m.descriptors = m2.descriptors;
                }
            }
        }

        assert(keypoints0List.size() == ns);
        assert(keypoints1List.size() == ns);
        assert(orientationsList.size() == ns);
        assert(harrisResponseImages.size() == ns);
        assert(harrisResponses.size() == ns);
        assert(scalesList.size() == ns);
        if (descrChoice.equals(DescriptorChoice.HSV)) {
            assert(descriptorsListH.size() == ns);
            assert(descriptorsListS.size() == ns);
            assert(descriptorsListV.size() == ns);
            assert(pyramidImagesH.size() == ns);
            assert(pyramidImagesS.size() == ns);
            assert(pyramidImagesV.size() == ns);
        } else if (descrChoice.equals(DescriptorChoice.GREYSCALE)) {
            assert(descriptorsList.size() == ns);
        }  else if (descrChoice.equals(DescriptorChoice.ALT)) {
            assert(descriptorsListAlt.size() == ns);
        }
        if (descriptorsMaskList != null) {
            assert(descriptorsMaskList.size() == ns);
        }

        ns = keypoints0List.size();
        for (int i = (ns - 1); i > -1; --i) {

            if (!keypoints0List.get(i).isEmpty()) {
                continue;
            }
            assert(keypoints1List.get(i).isEmpty());

            keypoints0List.remove(i);
            keypoints1List.remove(i);
            orientationsList.remove(i);
            scalesList.remove(i);
            pyramidImages.remove(i);
            harrisResponseImages.remove(i);
            harrisResponses.remove(i);

            if (descrChoice.equals(DescriptorChoice.HSV)) {
                descriptorsListH.remove(i);
                descriptorsListS.remove(i);
                descriptorsListV.remove(i);
                pyramidImagesH.remove(i);
                pyramidImagesS.remove(i);
                pyramidImagesV.remove(i);
            } else if (descrChoice.equals(DescriptorChoice.GREYSCALE)) {
                descriptorsList.remove(i);
            } else if (descrChoice.equals(DescriptorChoice.ALT)) {
                descriptorsListAlt.remove(i);
            }
            if (descriptorsMaskList != null) {
                descriptorsMaskList.remove(i);
            }
        }

        ns = keypoints0List.size();
        assert(keypoints0List.size() == ns);
        assert(keypoints1List.size() == ns);
        assert(orientationsList.size() == ns);
        assert(harrisResponseImages.size() == ns);
        assert(harrisResponses.size() == ns);
        assert(scalesList.size() == ns);
        if (descrChoice.equals(DescriptorChoice.HSV)) {
            assert(descriptorsListH.size() == ns);
            assert(descriptorsListS.size() == ns);
            assert(descriptorsListV.size() == ns);
            assert(pyramidImagesH.size() == ns);
            assert(pyramidImagesS.size() == ns);
            assert(pyramidImagesV.size() == ns);
        } else if (descrChoice.equals(DescriptorChoice.GREYSCALE)) {
            assert(descriptorsList.size() == ns);
        } else if (descrChoice.equals(DescriptorChoice.ALT)) {
            assert(descriptorsListAlt.size() == ns);
        }
        if (descriptorsMaskList != null) {
            assert(descriptorsMaskList.size() == ns);
        }
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

    public static int estimateNumberOfDefaultScales(int imageWidth, int imageHeight) {

        // duplicating the default loop logic to calculate number
        // of images in a pyramid for this image size.

        int nBetween = 3;

        int imgDimen = Math.min(imageWidth, imageHeight);

        int nr = (int)(Math.log(imgDimen)/Math.log(2));
        int s = 1;
        int winL = 2*s + 1;

        int w = imageWidth;
        int h = imageHeight;

        int ns = 1;

        for (int j = 0; j < (nr - 1); ++j) {
            if ((w < winL) || (h < winL)) {
                break;
            }
            if ((w <= defaultDecimationLimit) &&
                (h <= defaultDecimationLimit)) {
                break;
            }
            w /= 2;
            h /= 2;
            ns++;
        }

        int stop = ns;

        ns += 4;

        float f = 1.f/(nBetween + 1);

        int start = 1;
        if (nBetween > 4) {
            start = 0;
        }

        for (int i = start; i < stop - 1; ++i) {
            for (int j = 0; j < (nBetween + 1); ++j) {
                ns++;
            }
        }

        return ns;
    }

    /**
     * match descriptors using euclidean transformation evaluation from pairs in
     * feasible combinations of best matches.
     * This method is useful when a geometrical model is necessary to find the
     * object in an image where the background and foreground have changed and the
     * lighting may have changed, for example.
     * (In other words, the number of true positives is very low compared to the
     * total number of keypoints, so greedy or optimal matching and RANSAC are not
     * feasible solutions for this case).
     * Fpr known simpler conditions such as stereo-projections, RANSAC with
     * an epipolar or euclidean evaluator can be used instead with the top 45
     * results and that method will
     * be offered in this class one day.
     * @param scales1
     * @param scales2
     * @param descH1 H descriptors itemized by scales1.
     *     each descriptor is for a coordinate in similar location in
     *     lists keypointsX1, keypointsy1.
     * @param descS1 S descriptors itemized by scales1.
     *     each descriptor is for a coordinate in similar location in
     *     lists keypointsX1, keypointsy1.
     * @param descV1 V descriptors itemized by scales1.
     *     each descriptor is for a coordinate in similar location in
     *     lists keypointsX1, keypointsy1.
     * @param descH2 H descriptors itemized by scales2.
     *     each descriptor is for a coordinate in similar location in
     *     lists keypointsX2, keypointsy2.
     * @param descS2 S descriptors itemized by scales2.
     *     each descriptor is for a coordinate in similar location in
     *     lists keypointsX2, keypointsy2.
     * @param descV2 V descriptors itemized by scales2.
     *     each descriptor is for a coordinate in similar location in
     *     lists keypointsX2, keypointsy2.
     * @param keypointsX1 x coordinates of keypoints itemized by scales1
     * @param keypointsY1 y coordinates of keypoints itemized by scales1
     * @param keypointsX2 x coordinates of keypoints itemized by scales2
     * @param keypointsY2 y coordinates of keypoints itemized by scales2
     * @param scaleFactor the factor that will be the maximum dimension
     * of the object in the search image, for example, 2.

     * @return
     */
    public static CorrespondenceList matchDescriptors2(
        List<TwoDFloatArray> pyr1,
        List<TwoDFloatArray> pyr2,
        TFloatList scales1, TFloatList scales2,
        List<Descriptors> descH1, List<Descriptors> descS1,
        List<Descriptors> descV1,
        List<Descriptors> descH2, List<Descriptors> descS2,
        List<Descriptors> descV2,
        List<TIntList> keypointsX1, List<TIntList> keypointsY1,
        List<TIntList> keypointsX2, List<TIntList> keypointsY2,
        float scaleFactor) {

        List<CorrespondenceList> corList =
            matchDescriptors2(pyr1, pyr2,
            scales1, scales2,
            descH1, descS1, descV1, descH2, descS2, descV2,
            keypointsX1, keypointsY1, keypointsX2, keypointsY2,
            scaleFactor, 0, true);

        if (corList == null || corList.isEmpty()) {
            return null;
        }

        return corList.get(0);
    }

    /**
     * match descriptors using euclidean transformation evaluation from pairs in
     * feasible combinations of best matches.
     * This method is useful when a geometrical model is necessary to find the
     * object in an image where the background and foreground have changed and the
     * lighting may have changed, for example.
     * (In other words, the number of true positives is very low compared to the
     * total number of keypoints, so greedy or optimal matching and RANSAC are not
     * feasible solutions for this case).
     * Fpr known simpler conditions such as stereo-projections, RANSAC with
     * an epipolar or euclidean evaluator can be used instead with the top 45
     * results and that method will
     * be offered in this class one day.
     * NOTE: that one of the internal cost parameters is an estimate of the
     * maximum number of matchable points of any given scale compared to another,
     * so reduce the number of keypoints to truly matchable for best results.
     * @param scales1
     * @param scales2
     * @param descH1 H descriptors itemized by scales1.
     *     each descriptor is for a coordinate in similar location in
     *     lists keypointsX1, keypointsy1.
     * @param descS1 S descriptors itemized by scales1.
     *     each descriptor is for a coordinate in similar location in
     *     lists keypointsX1, keypointsy1.
     * @param descV1 V descriptors itemized by scales1.
     *     each descriptor is for a coordinate in similar location in
     *     lists keypointsX1, keypointsy1.
     * @param descH2 H descriptors itemized by scales2.
     *     each descriptor is for a coordinate in similar location in
     *     lists keypointsX2, keypointsy2.
     * @param descS2 S descriptors itemized by scales2.
     *     each descriptor is for a coordinate in similar location in
     *     lists keypointsX2, keypointsy2.
     * @param descV2 V descriptors itemized by scales2.
     *     each descriptor is for a coordinate in similar location in
     *     lists keypointsX2, keypointsy2.
     * @param keypointsX1 x coordinates of keypoints itemized by scales1
     * @param keypointsY1 y coordinates of keypoints itemized by scales1
     * @param keypointsX2 x coordinates of keypoints itemized by scales2
     * @param keypointsY2 y coordinates of keypoints itemized by scales2
     * @param scaleFactor the factor that will be the maximum dimension
     * of the object in the search image, for example, 2.
       @param sizeScaleFraction (suggested default is 0.12).
     * the expected largest difference in fraction
     * of image size between the template object and search object size.
     * For example, if the scale is 70% then that is represented by
     * a pyramidal image scale of approx 1.4, so is approx 0.1 difference from
     * an existing scale.  This number is used as a tolerance for error in
     * the descriptor.  For example, the same keypoint in two same images,
     * but one image scaled by 0.7 produces a best difference in descriptors
     * of about 20 bits (which is approx 0.08 * 256 bits) per band.
     * This number is used to calculate a tolerance for equivalent costs so
     * set it as accurately as possible.
     * @return
     */
    public static List<CorrespondenceList> matchDescriptors2(
        List<TwoDFloatArray> pyr1,
        List<TwoDFloatArray> pyr2,
        TFloatList scales1, TFloatList scales2,
        List<Descriptors> descH1, List<Descriptors> descS1,
        List<Descriptors> descV1,
        List<Descriptors> descH2, List<Descriptors> descS2,
        List<Descriptors> descV2,
        List<TIntList> keypointsX1, List<TIntList> keypointsY1,
        List<TIntList> keypointsX2, List<TIntList> keypointsY2,
        float scaleFactor, float sizeScaleFraction) {

        return matchDescriptors2(
            pyr1, pyr2,
            scales1, scales2,
            descH1, descS1, descV1, descH2, descS2, descV2,
            keypointsX1, keypointsY1, keypointsX2, keypointsY2,
            scaleFactor, sizeScaleFraction, false);
    }

    private static List<CorrespondenceList> matchDescriptors2(
        List<TwoDFloatArray> pyr1,
        List<TwoDFloatArray> pyr2,
        TFloatList scales1, TFloatList scales2,
        List<Descriptors> descH1, List<Descriptors> descS1,
        List<Descriptors> descV1,
        List<Descriptors> descH2, List<Descriptors> descS2,
        List<Descriptors> descV2,
        List<TIntList> keypointsX1, List<TIntList> keypointsY1,
        List<TIntList> keypointsX2, List<TIntList> keypointsY2,
        float scaleFactor, float sizeScaleFraction,
        boolean returnSingleAnswer) {

        if (scales1.size() != descH1.size() ||
            scales1.size() != descS1.size() ||
            scales1.size() != descV1.size() ||
            scales1.size() != keypointsX1.size() ||
            scales1.size() != keypointsY1.size()
            ) {
            throw new IllegalArgumentException("lists for datasets 1"
                + " must all be same lengths as scales1");
        }
        if (scales2.size() != descH2.size() ||
            scales2.size() != descS2.size() ||
            scales2.size() != descV2.size() ||
            scales2.size() != keypointsX2.size() ||
            scales2.size() != keypointsY2.size()
            ) {
            throw new IllegalArgumentException("lists for datasets 2"
                + " must all be same lengths as scales1");
        }

        //TODO: may need to revise this or allow it as a method argument:
        int pixTolerance = 10;

        int nBands = 3;
        int bitTolerance;
        if (returnSingleAnswer) {
            bitTolerance = 0;
        } else {
            bitTolerance = Math.round(sizeScaleFraction * nBands * 256);
        }
        int topLimit = Math.round(
            //0.17f
            //0.3f
            0.99f
            * nBands * 256) + bitTolerance;
        int minP1Diff = 5;

        MatchedPointsTransformationCalculator tc = new
            MatchedPointsTransformationCalculator();

        Transformer transformer = new Transformer();

        // distance portion of costs gets transformed to this reference frame
        float minScale1 = scales1.min();
        float diag1 = calculateDiagonal(keypointsX1, keypointsY1,
            scales1.indexOf(minScale1));

        // these 3 are used to normalize costs
        final double maxCost = nBands * 256;
        final double maxDist = diag1;
        final double dbgBitTol = bitTolerance / maxCost;
        // a rough estimate of maximum number of matchable points in any
        //     scale dataset comparison
        final int nMaxMatchable = calculateNMaxMatchable(keypointsX1, keypointsX2);
        System.out.println("nMaxMatchable=" + nMaxMatchable);

        int nMax1 = maxSize(keypointsX1);
        int nMax2 = maxSize(keypointsX2);
        int nMax = nMax1 * nMax2;

        if (diag1 <= 16) {
            minP1Diff = 2;
        }

        // --- best cost data ----
        double minCostTotal = Double.MAX_VALUE;
        double minCost1 = Double.MAX_VALUE;
        double minCost2 = Double.MAX_VALUE;
        double minCost3 = Double.MAX_VALUE;
        float minCostTScale = Float.MAX_VALUE;

        //runtime complexity of this vector depends upon the number of items
        // it is currently holding, so can set the capacity high and fill vector only
        // with items within bitTolerance of best, but too high might affect jvm
        // performance.
        // (note, can optimize this for very large results by occassionally ejecting
        // all values with cost > best + bitTolerance.)
        // TODO: a safe size is to set capacity to the number of unique
        // transformation parameter sets, but since that isn't known
        // until later without refactoring here, will make an assumption for now,
        // that size 100 is generous for number of top solutions.
        FixedSizeSortedVector<CObject> vec;
        if (returnSingleAnswer) {
            vec = new FixedSizeSortedVector<CObject>(1, CObject.class);
        } else {
            vec = new FixedSizeSortedVector<CObject>(10, CObject.class);
        }

        //CorrespondenceList minCostCor = null;
        //PairIntArray minCostTr2 = null;
        double[] minCostI = new double[nMax];
        double[] minDistI = new double[nMax];

        // temporary storage of corresp coords until object construction
        int[] m1x = new int[nMax];
        int[] m1y = new int[nMax];
        int[] m2x = new int[nMax];
        int[] m2y = new int[nMax];
        int mCount = 0;

        for (int i = 0; i < scales1.size(); ++i) {
            float pScale1 = scales1.get(i);
            Descriptors dH1 = descH1.get(i);
            Descriptors dS1 = descS1.get(i);
            Descriptors dV1 = descV1.get(i);
            TIntList kpX1 = keypointsX1.get(i);
            TIntList kpY1 = keypointsY1.get(i);
            int n1 = kpX1.size();
            if (n1 != dH1.descriptors.length) {
                throw new IllegalArgumentException("number of descriptors in "
                    + " d1 bitstrings must be same as keypoints1 length");
            }

            float factorToMinScale = pScale1 / minScale1;

            int minX = kpX1.min();
            int maxX = kpX1.max();
            int minY = kpY1.min();
            int maxY = kpY1.max();

            int objDimension = Math.max(maxX - minX, maxY - minY);
            int limit = Math.round(scaleFactor * objDimension);
            int limitSq = limit * limit;

            NearestNeighbor2D nn = new NearestNeighbor2D(
                makeSet(kpX1, kpY1), maxX + limit, maxY + limit);

            TObjectIntMap<PairInt> p1IndexMap = createIndexMap(kpX1, kpY1);

            for (int j = 0; j < scales2.size(); ++j) {
                Descriptors dH2 = descH2.get(j);
                Descriptors dS2 = descS2.get(j);
                Descriptors dV2 = descV2.get(j);
                TIntList kpX2 = keypointsX2.get(j);
                TIntList kpY2 = keypointsY2.get(j);
                int n2 = kpX2.size();
                if (n2 != dH2.descriptors.length) {
                    throw new IllegalArgumentException("number of descriptors in "
                        + " d2 bitstrings must be same as keypoints2 length");
                }

                //[n1][n2]
                int[][] cost = calcDescriptorCostMatrix(
                    new Descriptors[]{dH1, dS1, dV1},
                    new Descriptors[]{dH2, dS2, dV2});

                int nTot = n1 * n2;

                // storing points1 in pairintarray to transform
                // storing points2 in set to create a nearest neighbors

                PairIntArray a1 = new PairIntArray(nTot/2);
                PairIntArray a2 = new PairIntArray(nTot/2);

                TIntList a2Indexes = new TIntArrayList(nTot/2);

                Set<PairInt> s1 = new HashSet<PairInt>(nTot/2);
                Set<PairInt> s2 = new HashSet<PairInt>(nTot/2);

                FixedSizeSortedVector<CObject> vecD;
                FixedSizeSortedVector<CObject> vecF;
                if (returnSingleAnswer) {
                    vecD = new FixedSizeSortedVector<CObject>(1, CObject.class);
                    vecF = new FixedSizeSortedVector<CObject>(1, CObject.class);
                } else {
                    vecD = new FixedSizeSortedVector<CObject>(5, CObject.class);
                    vecF = new FixedSizeSortedVector<CObject>(5, CObject.class);
                }
                double[] minCostI_D = new double[nMax];
                double[] minDistI_D = new double[nMax];
                double[] minCostI_F = new double[nMax];
                double[] minDistI_F = new double[nMax];

                double minCostTotal_D = Double.MAX_VALUE;
                double minCost1_D = Double.MAX_VALUE;
                double minCost2_D = Double.MAX_VALUE;
                double minCost3_D = Double.MAX_VALUE;
                float minCostTScale_D = Float.MAX_VALUE;
                double minCostTotal_F = Double.MAX_VALUE;
                double minCost1_F = Double.MAX_VALUE;
                double minCost2_F = Double.MAX_VALUE;
                double minCost3_F = Double.MAX_VALUE;
                float minCostTScale_F = Float.MAX_VALUE;

                // may need to iterate to reduce pairs.size to szLimit
                float topLimit2 = topLimit;
                while (true) {
                    int count = 0;
                    for (int ii = 0; ii < n1; ++ii) {
                        for (int jj = 0; jj < n2; ++jj) {
                            int c = cost[ii][jj];
                            if (c > topLimit2) {
                                continue;
                            }
                            count++;
                        }
                    }
                    System.out.println(count);
                    //TODO: may need to revise this limit
                    if (count < 200) {
                         break;
                    }
                    topLimit2 -= (0.025f * nBands * 256);
                    if (topLimit2 == 0) {
                        topLimit2 = bitTolerance;
                        if (topLimit2 == 0) {
                            topLimit2 = nBands * 256 * 0.025f;
                        }
                        break;
                    }
                }

                List<QuadInt> pairs = new ArrayList<QuadInt>(nTot/2);
                TIntList costs = new TIntArrayList(nTot);
                for (int ii = 0; ii < n1; ++ii) {
                    PairInt p1 = new PairInt(kpX1.get(ii), kpY1.get(ii));
                    for (int jj = 0; jj < n2; ++jj) {
                        int c = cost[ii][jj];

                        if (c > topLimit2) {
                            continue;
                        }

                        PairInt p2 = new PairInt(kpX2.get(jj), kpY2.get(jj));

                        if (!s1.contains(p1)) {
                            a1.add(p1.getX(), p1.getY());
                            s1.add(p1);
                        }
                        if (!s2.contains(p2)) {
                            a2.add(p2.getX(), p2.getY());
                            a2Indexes.add(jj);
                            s2.add(p2);
                        }
                        pairs.add(new QuadInt(p1, p2));
                        costs.add(c);
                    }
                }

                System.out.println("i=" + i + " j=" + j + " nPairs=" + pairs.size());

                // --- calculate transformations in pairs and evaluate ----
                for (int ii = 0; ii < pairs.size(); ++ii) {

                    QuadInt pair1 = pairs.get(ii);

                    // image 1 point:
                    int t1X = pair1.getA();
                    int t1Y = pair1.getB();
                    // image 2 point:
                    int s1X = pair1.getC();
                    int s1Y = pair1.getD();

                    // choose all combinations of 2nd point within distance
                    // limit of point s1.
                    for (int jj = 0; jj < pairs.size(); ++jj) {

                        if (ii == jj) {
                            continue;
                        }

                        QuadInt pair2 = pairs.get(jj);

                        // image 1 point:
                        int t2X = pair2.getA();
                        int t2Y = pair2.getB();
                        // image 2 point:
                        int s2X = pair2.getC();
                        int s2Y = pair2.getD();

                        if ((t1X == t2X && t1Y == t2Y)
                            || (s1X == s2X && s1Y == s2Y)) {
                            continue;
                        }

                        int diffX = s1X - s2X;
                        int diffY = s1Y - s2Y;
                        int distSq = diffX * diffX + diffY * diffY;
                        if (distSq > limitSq) {
                            continue;
                        }
                        if ((distSq < minP1Diff*minP1Diff) ||
                            ((t1X - t2X)*(t1X - t2X) +
                             (t1Y - t2Y)*(t1Y - t2Y)
                            < minP1Diff)) {

                            continue;
                        }

                        // transform dataset 2 into frame 1
                        TransformationParameters params = tc.calulateEuclidean(
                            s1X, s1Y,
                            s2X, s2Y,
                            t1X, t1Y,
                            t2X, t2Y,
                            0, 0);

                        float tScale = params.getScale();

                        if (Math.abs(tScale - 1.0) > 0.20) {
                             continue;
                        }

                        if (Math.abs(tScale - 1)
                            > Math.abs(minCostTScale - 1)) {
                            continue;
                        }
                        //TODO: this may need to be revised.
                        // wanting to approach scale of '1' from whichever
                        // direction it is from.
                        // The current tests pass arguments such that the
                        // first images in pyramid1 and pyramid2 are the largerst
                        // for both.
                        // so if img1 is smaller than img2,
                        // tScale will be less than 1 and the approach will
                        // be towards 1, and skip when reach it.
                        // else if img2 is larger than img2, the
                        // tScale will be larger than 1, so decresing
                        // scales until reach "1"....
                        //might need to note the direction
                        // of approaching scale "1" and skip when reaches it

                        if (minCostTScale < Float.MAX_VALUE &&
                            (Math.abs(tScale - 1) >
                            (Math.abs(minCostTScale - 1)))) {
                            continue;
                        }

                        mCount = 0;

                        PairIntArray tr2 =
                            transformer.applyTransformation(params, a2);

                        double sum1 = 0;
                        double sum2 = 0;
                        double sum3 = 0;
                        double sum = 0;

                        double sum2_F = 0;
                        double sum_F = 0;

                        double sum2_D = 0;
                        double sum_D = 0;

                        for (int k = 0; k < tr2.getN(); ++k) {
                            int x2Tr = tr2.getX(k);
                            int y2Tr = tr2.getY(k);
                            int idx2 = a2Indexes.get(k);

                            Set<PairInt> nearest = null;
                            if ((x2Tr >= 0) && (y2Tr >= 0)
                                && (x2Tr <= (maxX + pixTolerance))
                                && (y2Tr <= (maxY + pixTolerance))) {
                                nearest = nn.findClosest(x2Tr, y2Tr, pixTolerance);
                            }

                            int minC = Integer.MAX_VALUE;
                            PairInt minCP1 = null;
                            int minIdx1 = 0;
                            if (nearest != null && !nearest.isEmpty()) {
                                for (PairInt p1 : nearest) {
                                    int idx1 = p1IndexMap.get(p1);
                                    int c = cost[idx1][idx2];
                                    if (c < minC) {
                                        minC = c;
                                        minCP1 = p1;
                                        minIdx1 = idx1;
                                    }
                                }
                            }

                            if (minCP1 != null) {
                                double scoreNorm = (nBands*256 - minC)/maxCost;
                                double costNorm = 1. - scoreNorm;
                                sum1 += costNorm;

                                // distances get multiplied by factorToMinScale
                                // to put them into reference frame of largest
                                // set 1 image (== minScale1 frame)

                                // TODO: applying tScale correctly needs knowledge
                                // of the real scale factor brtween object in the
                                // template and search images,
                                // so may need to calculate solutions for 2 vectors.
                                // one soln uses the tScale applied as a factor
                                // and the other soln uses the tScale applied as a
                                // divisor. the decider isn't completely clear yet,
                                //   but tentatively looks like the vec which
                                //   has the smaller minCost1 (== descriptor cost)
                                //   for the top item.
                                //
                                double dist = distance(x2Tr, y2Tr, minCP1);
                                double distNorm = dist *
                                    factorToMinScale / maxDist;
                                sum2 += distNorm;
                                sum2_F += (dist * tScale * factorToMinScale / maxDist);
                                sum2_D += (dist * factorToMinScale / (tScale * maxDist));

                                m2x[mCount] = kpX2.get(idx2);
                                m2y[mCount] = kpY2.get(idx2);
                                m1x[mCount] = minCP1.getX();
                                m1y[mCount] = minCP1.getY();
                                minCostI[mCount] = costNorm;
                                minDistI[mCount] = distNorm;
                                mCount++;

                            } else {
                                sum1 += 1;
                                sum2 += 1;
                                sum2_F += 1;
                                sum2_D += 1;
                            }
                        }

                        double cf = mCount;
                        if (cf > nMaxMatchable) {
                            cf = nMaxMatchable;
                        }
                        cf /= nMaxMatchable;
                        sum3 = 1. - cf;

                        sum = sum1 + sum2 + sum3;

                        sum_D = sum1 + sum2_D + sum3;

                        sum_F = sum1 + sum2_F + sum3;

                        if ((minCostTotal_D == Double.MAX_VALUE) ||
                            (sum_D <= (minCostTotal_D + bitTolerance))
                        ) {

                            if (sum_D < minCostTotal_D) {
                                minCostTotal_D = sum_D;
                                minCost1_D = sum1;
                                minCost2_D = sum2_D;
                                minCost3_D = sum3;
                                minCostTScale_D = tScale;
                            }

                            CorrespondenceList corr
                                = new CorrespondenceList(params.getScale(),
                                Math.round(params.getRotationInDegrees()),
                                Math.round(params.getTranslationX()),
                                Math.round(params.getTranslationY()),
                                0, 0, 0,
                                new ArrayList<PairInt>(), new ArrayList<PairInt>());

                            for (int mi = 0; mi < mCount; ++mi) {
                                corr.addMatch(
                                    new PairInt(m1x[mi], m1y[mi]),
                                    new PairInt(m2x[mi], m2y[mi]),
                                    (minCostI[mi] + minCostI[mi])
                                );
                            }

                            CObject cObj = new CObject(sum_D, corr, tr2);
                            vecD.add(cObj);
                        }

                        if ((minCostTotal_F == Double.MAX_VALUE) ||
                            (sum_F <= (minCostTotal_F + bitTolerance))
                        ) {

                            if (sum_F < minCostTotal_F) {
                                minCostTotal_F = sum_F;
                                minCost1_F = sum1;
                                minCost2_F = sum2_F;
                                minCost3_F = sum3;
                                minCostTScale_F = tScale;
                            }

                            CorrespondenceList corr
                                = new CorrespondenceList(params.getScale(),
                                Math.round(params.getRotationInDegrees()),
                                Math.round(params.getTranslationX()),
                                Math.round(params.getTranslationY()),
                                0, 0, 0,
                                new ArrayList<PairInt>(), new ArrayList<PairInt>());

                            for (int mi = 0; mi < mCount; ++mi) {
                                corr.addMatch(
                                    new PairInt(m1x[mi], m1y[mi]),
                                    new PairInt(m2x[mi], m2y[mi]),
                                    (minCostI[mi] + minCostI[mi])
                                );
                            }

                            CObject cObj = new CObject(sum_F, corr, tr2);
                            vecF.add(cObj);
                        }
                    }
                }

                //TODO: seems obvious now that the correct solution will be
                //      when the transformation scale is "1"...then the descriptor
                //      apertures will be the same size.
                //      making changes now.
                //      NOTE that masked descriptors could still be useful because
                //      they are more precise, but they might not be necessary for best cases,
                //      excepting those such as the android test 01 in which the object
                //      is small in the search image.

                System.out.println(
                    String.format(
                "i=%d j=%d minCost=%.2f c1=%.2f c2=%.2f c3=%.2f  c1Tol=%.2f  tS=%.2f",
                        i, j, (float) minCostTotal, (float) minCost1,
                        (float) minCost2, (float) minCost3, (float)dbgBitTol,
                        minCostTScale));
                System.out.println(
                    String.format(
                        "minCostD=%.2f c1=%.2f c2=%.2f c3=%.2f tS=%.2f",
                        (float) minCostTotal_D, (float) minCost1_D,
                        (float) minCost2_D, (float) minCost3_D,
                        minCostTScale_D));
                System.out.println(
                    String.format(
                        "minCostF=%.2f c1=%.2f c2=%.2f c3=%.2f tS=%.2f",
                        (float) minCostTotal_F, (float) minCost1_F,
                        (float) minCost2_F, (float) minCost3_F,
                        minCostTScale_F));

                if (vecD.getNumberOfItems() == 0) {
                    System.out.println("no matches for i=" + i + " j=" + j);
                    continue;
                }

                // if any of the top costs from descriptors is 0,
                // that vector should be chosen,
                // else, choose the one with smallest mincost3
                float[] c0 = new float[]{(float) minCostTotal,
                    (float) minCostTotal_D, (float) minCostTotal_F};
                float[] c1 = new float[]{(float) minCost1,
                    (float) minCost1_D, (float) minCost1_F};
                float[] c2 = new float[]{(float) minCost2,
                    (float) minCost2_D, (float) minCost2_F};
                float[] c3 = new float[]{(float) minCost3,
                    (float) minCost3_D, (float) minCost3_F};
                int[] indexes = new int[]{0, 1, 2};
                QuickSort.sortBy1stThen2ndThen3rd(c1, c3, c0, indexes);
                int vecIdx = -1;
                if (Math.abs(c1[0] - 0) < 0.001f) {
                    if (indexes[0] == 0) {
                        vecIdx = 0;
                    } else if (indexes[0] == 1) {
                        vecIdx = 1;
                    } else {
                        vecIdx = 2;
                    }
                } else {
                    c0 = new float[]{(float) minCostTotal,
                        (float) minCostTotal_D, (float) minCostTotal_F};
                    c1 = new float[]{(float) minCost1,
                        (float) minCost1_D, (float) minCost1_F};
                    c2 = new float[]{(float) minCost2,
                        (float) minCost2_D, (float) minCost2_F};
                    c3 = new float[]{(float) minCost3,
                        (float) minCost3_D, (float) minCost3_F};
                    indexes = new int[]{0, 1, 2};
                    // TODO: might need to sort by c3 first here and use bitTolerance w/ c0
                    QuickSort.sortBy1stThen2ndThen3rd(c0, c3, c1, indexes);
                    for (int ia = 0; ia < vec.getNumberOfItems(); ++ia) {
                        if (indexes[1] == 1) {
                            vecIdx = 1;
                        } else {
                            vecIdx = 2;
                        }
                    }
                }
                System.out.println("vecIdx=" + vecIdx);
                if (vecIdx == 1) {
                    if (minCostTotal_D < minCostTotal) {
                        System.out.println("Choosing scale as divisor");
                        vec = vecD;
                        minCostTotal = minCostTotal_D;
                        minCost1 = minCost1_D;
                        minCost2 = minCost2_D;
                        minCost3 = minCost3_D;
                        minCostTScale = minCostTScale_D;
                    }
                } else if (vecIdx != 0) { // vecIdx == 2
                    if (minCostTotal_F < minCostTotal) {
                        System.out.println("choosing scale as factor");
                        vec = vecF;
                        minCostTotal = minCostTotal_F;
                        minCost1 = minCost1_F;
                        minCost2 = minCost2_F;
                        minCost3 = minCost3_F;
                        minCostTScale = minCostTScale_F;
                    }
                }
            }// end loop over image j
        }

        if (vec.getNumberOfItems() == 0) {
            return null;
        }

        List<CorrespondenceList> topResults =
            new ArrayList<CorrespondenceList>();

        for (int i = 0; i < vec.getNumberOfItems(); ++i) {
            CObject a = vec.getArray()[i];
            if (a.cost > (minCostTotal + bitTolerance)) {
                break;
            }
            topResults.add(a.cCor);
        }

        return topResults;
    }

    //NOT READY FOR USE.
    public static List<CorrespondenceList> match0(
        ORB orb1, ORB orb2,
        Set<PairInt> labeledPoints1,
        List<Set<PairInt>> labeledPoints2) {

        /*
        uses the descriptors given and then makes masks
        for them using the labeled points.

        -- visits each octave pair
           -- calculates cost of descriptors
           -- uses the segmentation to calculate every
              permutation of 2 pairs of points.
              -- filter out high cost pairs.
           -- filters out 2 pair combinations with transformation scales not near 1
           -- keeps only the top 10 percent cost of items
              from the 2 pair list.
           -- evaluates the transformation using the transformed
              keypoints cost difference, distance from nearest
              neighbor and number of matches
           -- keeps the best of each j
           -- further compares bestJs with SSDs of intersecting
               transformed point sets of the matching keypoints
           -- top of those best is the returned result
        */

        if (!orb1.descrChoice.equals(orb2.descrChoice)) {
            throw new IllegalStateException(
            "orbs must contain same kind of descirptors");
        }
        int nBands = 3;
        if (orb1.descrChoice.equals(ORB.DescriptorChoice.HSV)) {
            if (orb1.getDescriptorsH() == null ||
                orb2.getDescriptorsH() == null) {
                throw new IllegalStateException("hsv descriptors must be created first");
            }
        } else if (orb1.descrChoice.equals(ORB.DescriptorChoice.ALT)) {
            if (orb1.getDescriptorsListAlt() == null ||
                orb2.getDescriptorsListAlt() == null) {
                throw new IllegalStateException(
                    "alt descriptors must be created first");
            }
            nBands = 1;
        } else if (orb1.descrChoice.equals(ORB.DescriptorChoice.GREYSCALE)) {
            if (orb1.getDescriptorsList() == null ||
                orb2.getDescriptorsList() == null) {
                throw new IllegalStateException(
                    "descriptors must be created first");
            }
            nBands = 1;
        }
       
        boolean useMasks = false;
        if (useMasks) {// initialize the masks, but discard the maps
            TObjectIntMap<PairInt> pointLabels1 = new TObjectIntHashMap<PairInt>();
            Set<PairInt> set = labeledPoints1;
            for (PairInt p : set) {
                pointLabels1.put(p, 0);
            }
            TObjectIntMap<PairInt> pointLabels2 = new TObjectIntHashMap<PairInt>();
            for (int i = 0; i < labeledPoints2.size(); ++i) {
                set = labeledPoints2.get(i);
                for (PairInt p : set) {
                    pointLabels2.put(p, i);
                }
            }
            orb1.createDescriptorMasks(pointLabels1);
            orb2.createDescriptorMasks(pointLabels2);
        }

        //TODO: may need to revise this or allow it as a method argument:
        int pixTolerance = 10;

        MatchedPointsTransformationCalculator tc = new
            MatchedPointsTransformationCalculator();

        Transformer transformer = new Transformer();

        TFloatList scales1 = extractScales(orb1.getScalesList());
        TFloatList scales2 = extractScales(orb2.getScalesList());

        // a rough estimate of maximum number of matchable points in any
        //     scale dataset comparison
        final int nMaxMatchable =
            Math.round(0.5f * calculateNMaxMatchable(
                orb1.getKeyPoint1List(), orb2.getKeyPoint1List()));
        //TODO: allow a factor to be passed in
        System.out.println("nMaxMatchable=" + nMaxMatchable);

        int nMax1 = maxSize(orb1.getKeyPoint1List());
        int nMax2 = maxSize(orb2.getKeyPoint1List());
        int nMax = nMax1 * nMax2;

        // --- best cost data ----
        double minCostTotal = Double.MAX_VALUE;
        double minCost1 = Double.MAX_VALUE;
        double minCost2 = Double.MAX_VALUE;
        double minCost3 = Double.MAX_VALUE;
        float minCostTScale = Float.MAX_VALUE;

        //runtime complexity of this vector depends upon the number of items
        // it is currently holding, so can set the capacity high and fill vector only
        // with items within bitTolerance of best, but too high might affect jvm
        // performance.
        // (note, can optimize this for very large results by occassionally ejecting
        // all values with cost > best + bitTolerance.)
        // TODO: a safe size is to set capacity to the number of unique
        // transformation parameter sets, but since that isn't known
        // until later without refactoring here, will make an assumption for now,
        // that size 100 is generous for number of top solutions.
        FixedSizeSortedVector<CObject3> minVec =
            new FixedSizeSortedVector<CObject3>(1, CObject3.class);

        for (int i = 0; i < scales1.size(); ++i) {
        //for (int i = 0; i < 1; ++i) {

            float scale1 = scales1.get(i);

            // coords are in ref frame of scale=1 of their pyramids
            TIntList kpX1 = orb1.getKeyPoint1List().get(i);
            TIntList kpY1 = orb1.getKeyPoint0List().get(i);
            int n1 = kpX1.size();

            TwoDFloatArray octaveImg1 = orb1.getPyramidImages().get(i);

            float diag1 = (float)Math.sqrt(octaveImg1.a.length *
                octaveImg1.a[0].length);
            final double maxDist = diag1;

            // create data structures in scaled reference frame
            TObjectIntMap<PairInt> p1KPIndexMap = new TObjectIntHashMap<PairInt>();
            TIntList kpX1_2 = new TIntArrayList(n1);
            TIntList kpY1_2 = new TIntArrayList(n1);
            for (int i3 = 0; i3 < n1; ++i3) {
                int x = Math.round((float)kpX1.get(i3)/scale1);
                int y = Math.round((float)kpY1.get(i3)/scale1);
                kpX1_2.add(x);
                kpY1_2.add(y);
                p1KPIndexMap.put(new PairInt(x, y), i3);
            }

            List<TIntList> pointIndexLists1 = new ArrayList<TIntList>();
            int ns = 1;
            for (int i3 = 0; i3 < ns; ++i3) {
                pointIndexLists1.add(new TIntArrayList());
            }

            TObjectIntMap<PairInt> pointLabels1 = new TObjectIntHashMap<PairInt>();
            Set<PairInt> set = labeledPoints1;
            Set<PairInt> setScaled = new HashSet<PairInt>();
            TIntList list = pointIndexLists1.get(0);
            assert(list != null);
            for (PairInt p : set) {
                int x = Math.round((float)p.getX()/scale1);
                int y = Math.round((float)p.getY()/scale1);
                PairInt p2 = new PairInt(x, y);
                pointLabels1.put(p2, 0);

                int idx = p1KPIndexMap.get(p2);
                list.add(idx);
                setScaled.add(p2);
            }
            Set<PairInt> shape = new HashSet<PairInt>(setScaled);

            int objDimension = Math.max(kpX1_2.max() - kpX1_2.min(),
                kpY1_2.max() - kpY1_2.min());
            int limit = Math.round(1.15f * objDimension);
            int limitSq = limit * limit;

            PairIntArray a1 = new PairIntArray(kpX1_2.size());
            TIntList a1Indexes = new TIntArrayList(kpX1_2.size());
            for (int ii = 0; ii < kpX1.size(); ++ii) {
                int x = kpX1.get(ii);
                int y = kpY1.get(ii);
                a1.add(x, y);
                a1Indexes.add(ii);
            }

            for (int j = 0; j < scales2.size(); ++j) {
            //for (int j = 0; j < 1; ++j) {

                float scale2 = scales2.get(j);

                // coords are in ref frame of scale=1 of their pyramids
                TIntList kpX2 = orb2.getKeyPoint1List().get(j);
                TIntList kpY2 = orb2.getKeyPoint0List().get(j);
                int n2 = kpX2.size();

                // create data structures in scaled reference frame
                TObjectIntMap<PairInt> p2KPIndexMap = new TObjectIntHashMap<PairInt>();
                TObjectIntMap<PairInt> p2KPIndexMap_2 = new TObjectIntHashMap<PairInt>();
                TIntList kpX2_2 = new TIntArrayList(n2);
                TIntList kpY2_2 = new TIntArrayList(n2);
                for (int j3 = 0; j3 < n2; ++j3) {
                    int x = Math.round((float) kpX2.get(j3) / scale2);
                    int y = Math.round((float) kpY2.get(j3) / scale2);
                    kpX2_2.add(x);
                    kpY2_2.add(y);
                    p2KPIndexMap_2.put(new PairInt(x, y), j3);
                    p2KPIndexMap.put(new PairInt(kpX2.get(j3), kpY2.get(j3)), j3);
                }

                List<TIntList> pointIndexLists2 = new ArrayList<TIntList>();
                int ns2 = labeledPoints2.size();
                for (int j3 = 0; j3 < ns2; ++j3) {
                    pointIndexLists2.add(new TIntArrayList());
                }
                TObjectIntMap<PairInt> pointLabels2 = new TObjectIntHashMap<PairInt>();
                for (int j3 = 0; j3 < ns2; ++j3) {
                    Set<PairInt> set2 = labeledPoints2.get(j3);
                    TIntList list2 = pointIndexLists2.get(j3);
                    assert (list2 != null);
                    for (PairInt p : set2) {
                        int x = Math.round((float) p.getX() / scale2);
                        int y = Math.round((float) p.getY() / scale2);
                        PairInt p2 = new PairInt(x, y);
                        pointLabels2.put(p2, j3);

                        int idx = p2KPIndexMap_2.get(p2);
                        list2.add(idx);
                    }
                }

                TwoDFloatArray octaveImg2 = orb2.getPyramidImages().get(j);

debugPrint(octaveImg1, octaveImg2, kpX1_2, kpY1_2, kpX2_2, kpY2_2, i, j);

                int maxX2 = orb2.getPyramidImages().get(0).a[0].length;
                int maxY2 = orb2.getPyramidImages().get(0).a.length;
                int maxX2_2 = octaveImg2.a[0].length;
                int maxY2_2 = octaveImg2.a.length;

                NearestNeighbor2D nn2 = new NearestNeighbor2D(
                    makeSet(kpX2, kpY2), maxX2 + limit, maxY2 + limit);

                int nTot = n1 * n2;

                //use descriptors with params here to reduce paramsList
                int[][] cost = null;
                if (useMasks) {
                    Descriptors[] desc1 = getDescriptors(orb1, i);
                    Descriptors[] desc2 = getDescriptors(orb2, j);

                    cost = calcMaskedDescriptorCostMatrixes(desc1, desc2,
                        orb1.getDescriptorsMaskList().get(i),
                        orb2.getDescriptorsMaskList().get(j))[1].a;
                } else {
                    Descriptors[] desc1 = getDescriptors(orb1, i);
                    Descriptors[] desc2 = getDescriptors(orb2, j);
                    cost = calcDescriptorCostMatrix(desc1, desc2);
                }

//combinations of pairs with same labels

                // storing them all to reduce nesting
                // quadint is idx1, idx2, idx3, idx4

                //TODO: can use the cost to more quickly filter the
                // pairs at creation time
                List<QuadInt> pairIndexes =
                    createPairLabelIndexes(cost, nBands,
                    pointIndexLists1, kpX1_2, kpY1_2,
                    pointIndexLists2, kpX2_2, kpY2_2);

                System.out.println("i=" + i + " j=" + j + " nPairs="
                    + pairIndexes.size());

                FixedSizeSortedVector<CObject4> vecP =
                    new FixedSizeSortedVector<CObject4>(
                    100,
                    //Math.round(0.1f * pairIndexes.size()),
                    //Math.round(0.01f * pairIndexes.size()),
                    CObject4.class);

                for (int ipi = 0; ipi < pairIndexes.size(); ++ipi) {

                    QuadInt q = pairIndexes.get(ipi);

                    int t1X = kpX1_2.get(q.getA());
                    int t1Y = kpY1_2.get(q.getA());
                    int t2X = kpX1_2.get(q.getB());
                    int t2Y = kpY1_2.get(q.getB());

                    int s1X = kpX2_2.get(q.getC());
                    int s1Y = kpY2_2.get(q.getC());
                    int s2X = kpX2_2.get(q.getD());
                    int s2Y = kpY2_2.get(q.getD());

                    // transform dataset 1 into frame 2
                    TransformationParameters params = tc.calulateEuclidean(
                        t1X, t1Y,
                        t2X, t2Y,
                        s1X, s1Y,
                        s2X, s2Y,
                        0, 0);

                    float tScale = params.getScale();

                    if (Math.abs(tScale - 1.0) > 0.15) {
                        continue;
                    }

                    int idx1_1 = p1KPIndexMap.get(new PairInt(t1X, t1Y));
                    int idx1_2 = p1KPIndexMap.get(new PairInt(t2X, t2Y));
                    int idx2_1 = p2KPIndexMap_2.get(new PairInt(s1X, s1Y));
                    int idx2_2 = p2KPIndexMap_2.get(new PairInt(s2X, s2Y));

                    int sum = cost[idx1_1][idx2_1] + cost[idx1_2][idx2_2];

                    CObject4 cObj = new CObject4(sum, params, q);

                    boolean added = vecP.add(cObj);

                }

                System.out.println("for i=" + i + " j=" + j
                    + " filtered nPairs=" + vecP.getNumberOfItems());

                double minCostJTotal = Double.MAX_VALUE;
                double minCostJ1 = Double.MAX_VALUE;
                double minCostJ2 = Double.MAX_VALUE;
                double minCostJ3 = Double.MAX_VALUE;
                float minCostJTScale = Float.MAX_VALUE;

                FixedSizeSortedVector<CObject3> vecJ =
                    new FixedSizeSortedVector<CObject3>(
                    1, CObject3.class);

                // --- evaluate cost of all keypoints, transformed

                for (int ipi = 0; ipi < vecP.getNumberOfItems(); ++ipi) {

                    CObject4 c = vecP.getArray()[ipi];

                    TransformationParameters params = c.params;

                    float tScale = params.getScale();

                    QuadInt q = c.q;

                    int t1X = kpX1_2.get(q.getA());
                    int t1Y = kpY1_2.get(q.getA());
                    int t2X = kpX1_2.get(q.getB());
                    int t2Y = kpY1_2.get(q.getB());

                    int s1X = kpX2_2.get(q.getC());
                    int s1Y = kpY2_2.get(q.getC());
                    int s2X = kpX2_2.get(q.getD());
                    int s2Y = kpY2_2.get(q.getD());

                    // ----- transform keypoints and sum the distance differences ----
                    PairIntArray tr1 = transformer
                        .applyTransformation(params, a1);
                    // trim to image dimensions
                    tr1 = trimToImageBounds(octaveImg2, tr1);

                    if (tr1.getN() == 0) {
                        continue;
                    }

                    //the matched kpx1,kpy1 kpx2,kpy2 coordinate pairs
                    int[] mp1 = new int[kpX1.size()];
                    int[] mp2 = new int[kpX1.size()];

                    double[] distAndCount =
                        sumKeypointDescAndDist(cost, 3,
                        a1Indexes, tr1, kpX1, kpY1,
                        nn2, p2KPIndexMap,
                        maxX2, maxY2, pixTolerance, maxDist,
                        mp1, mp2
                    );

                    double sumDesc = distAndCount[0];
                    double sumDist = distAndCount[1];
                    int np = (int)distAndCount[2];
                    int count = np;
                    
                    if (count < 2) {
                        continue;
                    }
                    if (count < mp1.length) {
                        mp1 = Arrays.copyOf(mp1, count);
                        mp2 = Arrays.copyOf(mp2, count);
                    }

                    if (count > nMaxMatchable) {
                        count = nMaxMatchable;
                    }

                    double cf = count;
                    if (cf > nMaxMatchable) {
                        cf = nMaxMatchable;
                    }
                    cf /= nMaxMatchable;
                    double sum3 = 1. - cf;
                    //double sum3 = nMaxMatchable - count;

                    sumDesc /= (double)count;
                    sumDist /= (double)count;
                    
                    // adding the count component for each
                    // descriptor
                    sum3 *= 2;
                    
                    double sum = sumDesc + sumDist + sum3;

                    //TODO: this will be streamlined when change to use
                    // minCost vars top value instead of a vector
                    PairInt[] m1 = new PairInt[np];
                    PairInt[] m2 = new PairInt[mp1.length];
                    for (int j3 = 0; j3 < m1.length; ++j3) {
                        int idx1 = mp1[j3];
                        int idx2 = mp2[j3];
                        assert(idx1 < kpX1.size() && idx1 > -1);
                        assert(idx2 < kpX2.size() && idx2 > -1);
                        m1[j3] = new PairInt(kpX1.get(idx1), kpY1.get(idx1));
                        m2[j3] = new PairInt(kpX2.get(idx2), kpY2.get(idx2));
                        assert(labeledPoints1.contains(m1[j3]));
                        assert(p1KPIndexMap.get(
                            new PairInt(kpX1_2.get(idx1), kpY1_2.get(idx1)))
                            == idx1);
                        assert(p2KPIndexMap_2.get(
                            new PairInt(kpX2_2.get(idx2), kpY2_2.get(idx2)))
                            == idx2);
                    }
                    CObject2 cObj2 = new CObject2(ipi, sum, sumDesc, sumDist,
                       sum3, m1, m2);
                    CObject3 cObj = new CObject3(cObj2, sum, 0, params);
                    cObj.keypointCount = count;
                    boolean added = vecJ.add(cObj);
                    if (added) {
                        minCostJTotal = sum;
                        minCostJ1 = sumDesc;
                        minCostJ2 = sumDist;
                        minCostJ3 = sum3;
                        minCostJTScale = tScale;

System.out.println(String.format(
    "i=%d j=%d ipi=%d ts=%.2f  c=%.2f c1=%.2f c2=%.2f c3=%.2f count=%d",
    i, j, ipi, tScale, 
    (float) sum, (float) sumDesc, (float) sumDist, (float) sum3,
    count));

                        if (true) {
                            CorrespondencePlotter plotter = new CorrespondencePlotter(
                                convertToImage(orb1.getPyramidImages().get(i)),
                                convertToImage(orb2.getPyramidImages().get(j)));
                            for (int ii = 0; ii < cObj.m1.length; ++ii) {
                                PairInt p1 = cObj.m1[ii];
                                PairInt p2 = cObj.m2[ii];
                                int x1 = Math.round((float) p1.getX() / scale1);
                                int y1 = Math.round((float) p1.getY() / scale1);
                                int x2 = Math.round((float) p2.getX() / scale2);
                                int y2 = Math.round((float) p2.getY() / scale2);
                                plotter.drawLineInAlternatingColors(
                                    x1, y1, x2, y2, 0);
                            }
                            String str = Integer.toString(i);
                            while (str.length() < 3) {
                                str = "0" + str;
                            }
                            String str2 = Integer.toString(j);
                            while (str2.length() < 3) {
                                str2 = "0" + str2;
                            }
                            str = str + "_" + str2;
                            try {
                                plotter.writeImage(
                                    "_indiv_masked_corres2_" + str + "_"
                                    + ipi);
                            } catch (IOException ex) {
                                Logger.getLogger(ORB.class.getName()).log(Level.SEVERE, null, ex);
                            }
                        }
                    }

                }// end loop over paramsList

                if (vecJ.getNumberOfItems() == 0) {
                    continue;
                }

                if (false) { //DEBUG
                    for (int k = 0; k < vecJ.getNumberOfItems(); ++k) {
                        CObject3 cobj = vecJ.getArray()[k];
                        CorrespondencePlotter plotter = new CorrespondencePlotter(
                            convertToImage(orb1.getPyramidImages().get(i)),
                            convertToImage(orb2.getPyramidImages().get(j)));
                        for (int ii = 0; ii < cobj.m1.length; ++ii) {
                            PairInt p1 = cobj.m1[ii];
                            PairInt p2 = cobj.m2[ii];
                            int x1 = Math.round((float) p1.getX() / scale1);
                            int y1 = Math.round((float) p1.getY() / scale1);
                            int x2 = Math.round((float) p2.getX() / scale2);
                            int y2 = Math.round((float) p2.getY() / scale2);
                            plotter.drawLineInAlternatingColors(
                                x1, y1, x2, y2, 0);
                        }
                        String str = Integer.toString(i);
                        while (str.length() < 3) {
                            str = "0" + str;
                        }
                        String str2 = Integer.toString(j);
                        while (str2.length() < 3) {
                            str2 = "0" + str2;
                        }
                        str = str + "_" + str2;

                        try {
                            plotter.writeImage("_mindiv_masked_corres3_" + str + "_"
                                + MiscDebug.getCurrentTimeFormatted());
                        } catch (IOException ex) {
                            Logger.getLogger(ORB.class.getName()).log(Level.SEVERE, null, ex);
                        }
                    
                        System.out.println(String.format(
                        "* %d %d ts=%.2f  c=%.2f c1=%.2f c2=%.2f c3=%.2f",
                        i, j, cobj.params.getScale(),
                        (float) cobj.cost, 
                        (float) cobj.costDesc, 
                        (float) cobj.costDist, 
                        (float) cobj.costCount));
                    }
                }
                
                if (vecJ.getNumberOfItems() == 0) {
                    System.out.println("no matches for i=" + i + " j=" + j);
                    continue;
                }

                // if expand capacity of minVec, add up to capacity here
                minVec.add(vecJ.getArray()[0]);
               
            }// end loop over image j
        }

        if (minVec.getNumberOfItems() == 0) {
            return null;
        }

        List<CorrespondenceList> topResults =
            new ArrayList<CorrespondenceList>();

        for (int i = 0; i < minVec.getNumberOfItems(); ++i) {
            CObject3 a = minVec.getArray()[i];
            if (a.cost > minCostTotal) {
                break;
            }
            CorrespondenceList cor = new CorrespondenceList(a.params,
                a.m1, a.m2);
            topResults.add(cor);
        }

        return topResults;
    }

    /**
     * NOTE: preliminary results show that this matches the right pattern as
     * a subset of the object, but needs to be followed by a slightly larger
     * aggregated search by segmentation cells using partial shape matcher
     * for example.  This was started in ShapeFinder, but needs to be
     * adjusted for a search given seed cells and possibly improved for the
     * other TODO items).
     * @param keypoints1
     * @param keypoints2
     * @param mT
     * @param mS
     * @param nn
     * @param minMaxXY2
     * @param limit
     * @param tIndexes
     * @param idx1P2CostMap
     * @param indexes
     * @param costs
     * @return
     */
    private static List<CorrespondenceList> completeUsingCombinations(
        List<PairInt> keypoints1, List<PairInt> keypoints2,
        PairIntArray mT, PairIntArray mS,
        NearestNeighbor2D nn, int[] minMaxXY2, int limit,
        TIntList tIndexes,
        TIntObjectMap<TObjectIntMap<PairInt>> idx1P2CostMap,
        PairInt[] indexes, int[] costs, int bitTolerance,
        int nBands
    ) {

        int nTop = mT.getN();

        System.out.println("have " + nTop + " sets of points for "
            + " n of k=2 combinations");

        // need to make pairs of combinations from mT,mS
        //  to calcuate euclidean transformations and evaluate them.
        // -- can reduce the number of combinations by imposing a
        //    distance limit on separation of feasible pairs

        int limitSq = limit * limit;

        MatchedPointsTransformationCalculator tc = new
            MatchedPointsTransformationCalculator();

        Transformer transformer = new Transformer();

        /* can try to use the tolerance in 2 different ways:
           (1) to adjust the weight of cost.
               weight = 1 - (tolerance/(256*nBands)).
               This approach might favor regions of many false keypoints
               such as highly textured regions.
           (2) keep the minCost solution, but also keep all solutions
               that are within cost + tolerance of it.
               This doesn't discard the true solution.
               ---> these will need additional information to distinguish
                    between solutions.
        procedding with (2).
        */

        // this fixed size sorted vector is faster for shorter arrays.
        // TODO: consider ways to robustly set the size from the cost
        // statistics to ensure the vector will always contain the
        // correct solution even if not in top position.
        int nt = mT.getN();
        FixedSizeSortedVector<CObject> vec = new
            FixedSizeSortedVector<CObject>(nt, CObject.class);

        double minCost = Double.MAX_VALUE;
        //CorrespondenceList minCostCor = null;
        //PairIntArray minCostTrT = null;
        double[] minCostI = new double[nTop];
        double[] minDistI = new double[nTop];

        // temporary storage of corresp coords until object construction
        int[] m1x = new int[nTop];
        int[] m1y = new int[nTop];
        int[] m2x = new int[nTop];
        int[] m2y = new int[nTop];
        int mCount = 0;

        for (int i = 0; i < nTop; ++i) {
            int t1X = mT.getX(i);
            int t1Y = mT.getY(i);
            int s1X = mS.getX(i);
            int s1Y = mS.getY(i);

            // choose all combinations of 2nd point within distance
            // limit of point s1.
            for (int j = (i + 1); j < mS.getN(); ++j) {
                int t2X = mT.getX(j);
                int t2Y = mT.getY(j);
                int s2X = mS.getX(j);
                int s2Y = mS.getY(j);

                if ((t1X == t2X && t1Y == t2Y) ||
                    (s1X == s2X && s1Y == s2Y)) {
                    continue;
                }

                int diffX = s1X - s2X;
                int diffY = s1Y - s2Y;
                int distSq = diffX * diffX + diffY * diffY;
                if (distSq > limitSq) {
                    continue;
                }

                // -- calculate euclid transformation
                // -- evaluate the fit
                TransformationParameters params = tc.calulateEuclidean(
                    t1X, t1Y,
                    t2X, t2Y,
                    s1X, s1Y,
                    s2X, s2Y,
                    0, 0);

                float scale = params.getScale();

                mCount = 0;

                // template object transformed
                PairIntArray trT =
                    transformer.applyTransformation(params, mT);

                /*
                two components to the evaluation and both need normalizations
                so that their contributions to total result are
                equally weighted.

                (1) descriptors:
                    -- score is sum of each matched (3*256 - cost)
                    -- the normalization is the maximum possible score,
                       so will use the number of template points.
                       --> norm = nTemplate * 3 * 256
                    -- normalized score = (3*256 - cost)/norm
                   ==> normalized cost = 1 - ((3*256 - cost)/norm)
                (2) spatial distances from transformed points:
                   -- sum of distances within limit
                      and replacement of distance by limit if no matching
                      nearest neighbor is found.
                   -- divide each distance by the transformation scale
                      to compare same values
                   -- divide the total sum by the total max possible
                      --> norm = nTemplate * limit / scale

                Then the total cost is (1) + (2) and the min cost
                among all of these combinations is the resulting
                correspondence list
                */

                double maxCost = nBands * 256;
                double maxDist = limit/scale;

                double sum1 = 0;
                double sum2 = 0;
                double sum = 0;

                for (int k = 0; k < trT.getN(); ++k) {
                    int xTr = trT.getX(k);
                    int yTr = trT.getY(k);

                    int idx1 = tIndexes.get(k);

                    Set<PairInt> nearest = null;
                    if ((xTr >= 0) && (yTr >= 0) &&
                        (xTr <= (minMaxXY2[1] + limit)) &&
                        (yTr <= (minMaxXY2[3] + limit))) {
                        nearest = nn.findClosest(xTr, yTr, limit);
                    }

                    int minC = Integer.MAX_VALUE;
                    PairInt minCP2 = null;

                    if (nearest != null && !nearest.isEmpty()) {
                        TObjectIntMap<PairInt> cMap = idx1P2CostMap.get(idx1);
                        for (PairInt p2 : nearest) {
                            if (!cMap.containsKey(p2)) {
                                continue;
                            }
                            int c = cMap.get(p2);
                            if (c < minC) {
                                minC = c;
                                minCP2 = p2;
                            }
                        }
                    }
                    if (minCP2 != null) {
                        double scoreNorm = (3*256 - minC)/maxCost;
                        double costNorm = 1. - scoreNorm;
                        sum1 += costNorm;

                        double dist = distance(xTr, yTr, minCP2);
                        double distNorm = dist/maxDist;
                        sum2 += distNorm;

                        m1x[mCount] = keypoints1.get(idx1).getX();
                        m1y[mCount] = keypoints1.get(idx1).getY();
                        m2x[mCount] = minCP2.getX();
                        m2y[mCount] = minCP2.getY();
                        minCostI[mCount] = costNorm;
                        minDistI[mCount] = distNorm;
                        mCount++;

                    } else {
                        sum1 += 1;
                        sum2 += 1;
                    }
                }
                sum = sum1 + sum2;

                if ((minCost == Double.MAX_VALUE) ||
                    (sum < (minCost + bitTolerance))) {

                    if (sum < minCost) {
                        minCost = sum;
                    }

                    List<PairInt> m1 = new ArrayList<PairInt>();
                    List<PairInt> m2 = new ArrayList<PairInt>();
                    CorrespondenceList corr =
                        new CorrespondenceList(
                        params.getScale(),
                        Math.round(params.getRotationInDegrees()),
                        Math.round(params.getTranslationX()),
                        Math.round(params.getTranslationY()),
                        0, 0, 0, m1, m2);

                    for (int mi = 0; mi < mCount; ++mi) {
                        m1.add(new PairInt(m1x[mi], m1y[mi]));
                        m2.add(new PairInt(m2x[mi], m2y[mi]));
                    }

                    CObject cObj = new CObject(sum, corr, trT);

                    vec.add(cObj);
                }
            }
        }

        if (vec.getNumberOfItems() == 0) {
            return null;
        }

        List<CorrespondenceList> topResults =
            new ArrayList<CorrespondenceList>();

        for (int i = 0; i < vec.getNumberOfItems(); ++i) {
            CObject a = vec.getArray()[i];
            if (a.cost > (minCost + bitTolerance)) {
                break;
            }
            topResults.add(a.cCor);
        }

        return topResults;
    }

    private static class CObject implements Comparable<CObject> {
        final double cost;
        final CorrespondenceList cCor;
        final PairIntArray transformedTemplate;
        public CObject(double cost, CorrespondenceList cL,
            PairIntArray templTr) {
            this.cost = cost;
            this.cCor = cL;
            this.transformedTemplate = templTr;
        }

        @Override
        public int compareTo(CObject other) {
            if (cost < other.cost) {
                return -1;
            } else if (cost > other.cost) {
                return 1;
            } else {
                int n1 = cCor.getPoints1().size();
                int n2 = other.cCor.getPoints1().size();
                if (n1 > n2) {
                    return -1;
                } else if (n1 < n2) {
                    return 1;
                }
            }
            return 0;
        }
    }

    private static class CObject4 implements Comparable<CObject4>{
        final double cost;
        final TransformationParameters params;
        final QuadInt q;
        public CObject4(double sum, TransformationParameters params,
            QuadInt q) {
            this.cost = sum;
            this.q = q;
            this.params = params;
        }

        @Override
        public int compareTo(CObject4 other) {
            if (cost < other.cost) {
                return -1;
            } else if (cost > other.cost) {
                return 1;
            }
            return 0;
        }
    }

    private static class CObject3 implements Comparable<CObject3>{
        final double cost;
        final double costDesc;
        final double costDist;
        final double costCount;
        final int index;
        final PairInt[] m1;
        final PairInt[] m2;
        final double sumPatch;
        final TransformationParameters params;
        QuadInt q;
        int keypointCount;
        public CObject3(CObject2 cObject2, double sum, double sumPatch,
            TransformationParameters params) {
            this.sumPatch = sumPatch;
            this.cost = sum;
            this.costDesc = cObject2.costDesc;
            this.costDist = cObject2.costDist;
            this.costCount = cObject2.costCount;
            this.index = cObject2.index;
            this.m1 = cObject2.m1;
            this.m2 = cObject2.m2;
            this.params = params;
        }

        @Override
        public int compareTo(CObject3 other) {
            if (cost < other.cost) {
                return -1;
            } else if (cost > other.cost) {
                return 1;
            }
            return 0;
        }
    }

    private static class CObject2 implements Comparable<CObject2> {
        final double cost;
        final double costDesc;
        final double costDist;
        final double costCount;
        final int index;
        final PairInt[] m1;
        final PairInt[] m2;
        public CObject2(int index, double cost, double costDesc, double costDist,
            double costCount, PairInt[] matched1, PairInt[] matched2) {
            this.cost = cost;
            this.index = index;
            this.m1 = matched1;
            this.m2 = matched2;
            this.costDesc = costDesc;
            this.costDist = costDist;
            this.costCount = costCount;
        }

        @Override
        public int compareTo(CObject2 other) {
            if (cost < other.cost) {
                return -1;
            } else if (cost > other.cost) {
                return 1;
            }
            return 0;
        }
    }

    /**
     * greedy matching of d1 to d2 by min cost, with unique mappings for
     * all indexes.
     *
     * @param d1
     * @param d2
     * @return matches - two dimensional int array of indexes in d1 and
     * d2 which are matched.
     */
    public static int[][] matchDescriptors(VeryLongBitString[] d1,
        VeryLongBitString[] d2,
        List<PairInt> keypoints1, List<PairInt> keypoints2) {

        int n1 = d1.length;
        int n2 = d2.length;

        //[n1][n2]
        int[][] cost = calcDescriptorCostMatrix(d1, d2);

        int[][] matches = greedyMatch(keypoints1, keypoints2, cost);
        // greedy or optimal match can be performed here.

        // NOTE: some matching problems might benefit from using the spatial
        //   information at the same time.  for those, will consider adding
        //   an evaluation term for these descriptors to a specialization of
        //   PartialShapeMatcher.java

        return matches;
    }

    /**
     * greedy matching of d1 to d2 by min difference, with unique mappings for
     * all indexes.
     * NOTE that if 2 descriptors match equally well, either one
     * might get the assignment.
     * Consider using instead, matchDescriptors2 which matches
     * by descriptor and relative spatial location.
     *
     * @param d1
     * @param d2
     * @param keypoints2
     * @param keypoints1
     * @return matches - two dimensional int array of indexes in d1 and
     * d2 which are matched.
     */
    public static int[][] matchDescriptors(Descriptors[] d1, Descriptors[] d2,
        List<PairInt> keypoints1, List<PairInt> keypoints2) {

        if (d1.length != d2.length) {
            throw new IllegalArgumentException("d1 and d2 must"
                + " be same length");
        }

        int n1 = d1[0].descriptors.length;
        int n2 = d2[0].descriptors.length;

        if (n1 != keypoints1.size()) {
            throw new IllegalArgumentException("number of descriptors in "
                + " d1 bitstrings must be same as keypoints1 length");
        }
        if (n2 != keypoints2.size()) {
            throw new IllegalArgumentException("number of descriptors in "
                + " d2 bitstrings must be same as keypoints2 length");
        }

        //[n1][n2]
        int[][] cost = calcDescriptorCostMatrix(d1, d2);

        int[][] matches = greedyMatch(keypoints1, keypoints2, cost);
        // greedy or optimal match can be performed here.

        // NOTE: some matching problems might benefit from using the spatial
        //   information at the same time.  for those, will consider adding
        //   an evaluation term for these descriptors to a specialization of
        //   PartialShapeMatcher.java

        return matches;
    }

    private static int[][] greedyMatch(List<PairInt> keypoints1,
        List<PairInt> keypoints2, int[][] cost) {

        int n1 = keypoints1.size();
        int n2 = keypoints2.size();

        // for the greedy match, separating the index information from the cost
        // and then sorting by cost
        int nTot = n1 * n2;

        PairInt[] indexes = new PairInt[nTot];
        int[] costs = new int[nTot];
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

        Set<PairInt> set1 = new HashSet<PairInt>();
        Set<PairInt> set2 = new HashSet<PairInt>();

        List<PairInt> matches = new ArrayList<PairInt>();

        // visit lowest costs (== differences) first
        for (int i = 0; i < nTot; ++i) {
            PairInt index12 = indexes[i];
            int idx1 = index12.getX();
            int idx2 = index12.getY();
            PairInt p1 = keypoints1.get(idx1);
            PairInt p2 = keypoints2.get(idx2);
            if (set1.contains(p1) || set2.contains(p2)) {
                continue;
            }

            //System.out.println("p1=" + p1 + " " + " p2=" + p2 + " cost=" + costs[i]);

            matches.add(index12);
            set1.add(p1);
            set2.add(p2);
        }

        int[][] results = new int[matches.size()][2];
        for (int i = 0; i < matches.size(); ++i) {
            results[i][0] = matches.get(i).getX();
            results[i][1] = matches.get(i).getY();
        }

        return results;
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
     * calculate a cost matrix composed of the sum of XOR of each descriptor
     * in d1 to d2.
     * Three matrixes are currently returned.
     * The first is the cost calculated from non-masked bits only.
     * The second is the first then multiplied by the area ratio of
     * masked bits to non-masked bits.  The second is the the data
     * which is present, averaged over the entire descriptor aperture.
     * The third is the cost calculated from non-masked bits, then
     * all masked bits are added as set bits.
     * The second matrix is probably most useful while the first and
     * third are min and max bounds of the actual costs.
     *
     * @param desc1 two dimensional array with first being keypoint indexes and
     * second dimension being descriptor.
     * @param desc2  two dimensional array with first being keypoint indexes and
     * second dimension being descriptor.
     * @param descMask1 bit mask for desc1 descriptors.  set bits are
     * pixels outside of the segmentation cell for the keypoint of the
     * descriptor.
     * @param descMask2 bit mask for desc2 descriptors.  set bits are
     * pixels outside of the segmentation cell for the keypoint of the
     * descriptor.
     * @return Three cost matrixes:
     * The first is the cost calculated from non-masked bits only.
     * The second is the first then multiplied by the area ratio of
     * masked bits to non-masked bits.  The second is the the data
     * which is present, averaged over the entire descriptor aperture.
     * The third is the cost calculated from non-masked bits, then
     * all masked bits are added as set bits.
     * The second matrix is probably most useful while the first and
     * third are min and max bounds of the actual costs.
     */
    public static TwoDIntArray[] calcMaskedDescriptorCostMatrixes(
        Descriptors[] desc1, Descriptors[] desc2,
        Descriptors descMask1, Descriptors descMask2) {

        if (desc1.length != desc2.length) {
            throw new IllegalArgumentException("desc1 and desc2 lengths"
                + " must be same");
        }

        int nBands = desc1.length;

        int n1 = desc1[0].descriptors.length;
        int n2 = desc2[0].descriptors.length;
        if (n1 != descMask1.descriptors.length) {
            throw new IllegalArgumentException("desc1 descriptors and mask"
                + " must be same lengths");
        }
        if (n2 != descMask2.descriptors.length) {
            throw new IllegalArgumentException("desc2 descriptors and mask"
                + " must be same lengths");
        }

        // d1 contains multiple descriptors for same points, such as
        // descriptors for H, S, and V

        int[][] cost0 = new int[n1][n2];
        int[][] cost1 = new int[n1][n2];
        int[][] cost2 = new int[n1][n2];
        for (int i = 0; i < n1; ++i) {
            cost0[i] = new int[n2];
            cost1[i] = new int[n2];
            cost2[i] = new int[n2];
        }

        int maxCost = nBands * 256;

        float dLength = desc1[0].descriptors[0].getInstantiatedBitSize();
        float fracDiv, nNonMasked;
        long nSetBits, nMaskedBits;
        VeryLongBitString costIJ, d12, mComb;
        for (int k = 0; k < nBands; ++k) {

            VeryLongBitString[] d1 = desc1[k].descriptors;
            VeryLongBitString[] d2 = desc2[k].descriptors;
            assert(d1.length == n1);
            assert(d2.length == n2);

            for (int i = 0; i < n1; ++i) {

                VeryLongBitString m1 = descMask1.descriptors[i];

                for (int j = 0; j < n2; ++j) {

                    VeryLongBitString m2 = descMask2.descriptors[j];

                    // the bits which are different:
                    d12 = d1[i].xor(d2[j]);

                    // combine the set bits of the masks
                    mComb = m1.or(m2);

                    // any place where m1 or m2 is '1' should be set to 0
                    costIJ = VeryLongBitString.subtract(d12, mComb);

                    // the remaining set bits in costIJ are the minimum cost
                    //     of the descriptors
                    nSetBits = costIJ.getNSetBits();
                    assert(nSetBits <= 256);

                    cost0[i][j] += nSetBits;
                    assert(cost0[i][j] <= maxCost);

                    // cost1 averages costIJ over set non-masked bits
                    //     and scales to full aperature by that
                    nMaskedBits = mComb.getNSetBits();
                    nNonMasked = dLength - nMaskedBits;
                    fracDiv = (nSetBits/nNonMasked) * nMaskedBits;
                    if (nNonMasked > 0) {

                        int c = Math.round(nSetBits + fracDiv);
                        assert(c <= 256);

                        cost1[i][j] += c;

                        assert(cost1[i][j] <= maxCost);
                    }

                    // cost 2 adds all masked bits as worse case scenario
                    cost2[i][j] += (nSetBits + nMaskedBits);

                    assert((nSetBits + nMaskedBits) <= 256);

                    assert(cost1[i][j] <= maxCost);
                }
            }
        }

        TwoDIntArray[] costs = new TwoDIntArray[] {
            new TwoDIntArray(cost0),
            new TwoDIntArray(cost1),
            new TwoDIntArray(cost2)
        };

        return costs;
    }

    private static void populateCorrespondence(PairIntArray left,
        PairIntArray right, PairInt[] kpIdxs,
        List<PairInt> keypoints1, List<PairInt> keypoints2) {

        for (PairInt p : kpIdxs) {
            PairInt p1 = keypoints1.get(p.getX());
            PairInt p2 = keypoints2.get(p.getY());
            left.add(p1.getX(), p1.getY());
            right.add(p2.getX(), p2.getY());
        }
    }

    private static double distance(int x, int y, PairInt b) {

        int diffX = x - b.getX();
        int diffY = y - b.getY();

        double dist = Math.sqrt(diffX * diffX + diffY * diffY);
        return dist;
    }

    private static TObjectIntMap<PairInt> createIndexMap(
        List<Set<PairInt>> segmentedCells) {

        TObjectIntMap<PairInt> map = new TObjectIntHashMap<PairInt>();

        for (int i = 0; i < segmentedCells.size(); ++i) {
            for (PairInt p : segmentedCells.get(i)) {
                map.put(p, i);
            }
        }

        return map;
    }

    private static CorrespondenceList setToBestUnique(
        CorrespondenceList minCostCor,
        double[] minCostI, double[] minDistI) {

        if (minCostCor == null) {
            return null;
        }

        /*
        make indexes array and sort by incr cost

        then uniquely assign pairs from lowest costs
        */
        int n = minCostCor.getPoints1().size();
        int[] idxs = new int[n];
        float[] totCost = new float[n];
        for (int i = 0; i < n; ++i) {
            idxs[i] = i;
            totCost[i] = (float)(minCostI[i] + minDistI[i]);
        }
        QuickSort.sortBy1stArg(totCost, idxs);

        Set<PairInt> set1 = new HashSet<PairInt>();
        Set<PairInt> set2 = new HashSet<PairInt>();
        List<PairInt> m1 = new ArrayList<PairInt>();
        List<PairInt> m2 = new ArrayList<PairInt>();

        for (int i = 0; i < n; ++i) {
            int idx = idxs[i];
            PairInt p1 = minCostCor.getPoints1().get(idx);
            PairInt p2 = minCostCor.getPoints2().get(idx);
            if (set1.contains(p1) || set2.contains(p2)) {
                continue;
            }
            m1.add(p1);
            m2.add(p2);
            set1.add(p1);
            set2.add(p2);
        }
        minCostCor.getPoints1().clear();
        minCostCor.getPoints2().clear();

        minCostCor.getPoints1().addAll(m1);
        minCostCor.getPoints2().addAll(m2);

        return minCostCor;
    }

    private static int distance(PairInt p1, PairInt p2) {
        int diffX = p1.getX() - p2.getX();
        int diffY = p1.getY() - p2.getY();
        return (int)Math.sqrt(diffX * diffX + diffY * diffY);
    }

    private static TObjectIntMap<PairInt> createIndexMap(
        TIntList xList, TIntList yList) {

        TObjectIntMap<PairInt> map = new TObjectIntHashMap<PairInt>();

        for (int i = 0; i < xList.size(); ++i) {
            int x = xList.get(i);
            int y = yList.get(i);
            map.put(new PairInt(x, y), i);
        }

        return map;
    }

    private static int maxSize(List<TIntList> a) {

        int maxSz = Integer.MIN_VALUE;
        for (TIntList b : a) {
            int sz = b.size();
            if (sz > maxSz) {
                maxSz = sz;
            }
        }
        return maxSz;
    }

    private static float calculateDiagonal(List<TIntList> keypointsX1,
        List<TIntList> keypointsY1, int idx) {

        TIntList x1 = keypointsX1.get(idx);
        TIntList y1 = keypointsY1.get(idx);

        int maxX = x1.max();
        int maxY = y1.max();

        return (float)Math.sqrt(maxX * maxX + maxY * maxY);
    }

    private static float calculateDiagonal2(List<TwoDFloatArray> pyramidImages,
        int idx) {

        int w = pyramidImages.get(idx).a.length;
        int h = pyramidImages.get(idx).a[0].length;

        double diag = Math.sqrt(w * w + h * h);

        return (float)diag;
    }

    private static Set<PairInt> makeSet(TIntList kpX1, TIntList kpY1) {

        Set<PairInt> set = new HashSet<PairInt>();
        for (int i = 0; i < kpX1.size(); ++i) {
            PairInt p = new PairInt(kpX1.get(i), kpY1.get(i));
            set.add(p);
        }

        return set;
    }

    private static int calculateNMaxMatchable(List<TIntList> keypointsX1,
        List<TIntList> keypointsX2) {

        int nMaxM = Integer.MIN_VALUE;

        for (int i = 0; i < keypointsX1.size(); ++i) {
            int n1 = keypointsX1.get(i).size();
            for (int j = 0; j < keypointsX2.size(); ++j) {
                int n2 = keypointsX2.get(j).size();
                int min = Math.min(n1, n2);
                if (min > nMaxM) {
                    nMaxM = min;
                }
            }
        }

        return nMaxM;
    }

    private static void debugPlot(int i, int j,
        FixedSizeSortedVector<CObject> vecD,
        FixedSizeSortedVector<CObject> vecF,
        TwoDFloatArray pyr1, TwoDFloatArray pyr2, float s1, float s2) {

        Image img1 = convertToImage(pyr1);
        Image img2 = convertToImage(pyr2);

        try {
            for (int i0 = 0; i0 < 2; ++i0) {
                CorrespondenceList cor = null;
                if (i0 == 0) {
                    cor = vecD.getArray()[0].cCor;
                } else {
                    cor = vecF.getArray()[0].cCor;
                }
                Image img1Cp = img1.copyImage();
                Image img2Cp = img2.copyImage();
                CorrespondencePlotter plotter = new CorrespondencePlotter(
                    img1Cp, img2Cp);
                for (int ii = 0; ii < cor.getPoints1().size(); ++ii) {
                    PairInt p1 = cor.getPoints1().get(ii);
                    PairInt p2 = cor.getPoints2().get(ii);
                    int x1 = Math.round(p1.getX()/s1);
                    int y1 = Math.round(p1.getY()/s1);
                    int x2 = Math.round(p2.getX()/s2);
                    int y2 = Math.round(p2.getY()/s2);
                    plotter.drawLineInAlternatingColors(x1, y1, x2, y2, 0);
                }
                String strI = Integer.toString(i);
                while (strI.length() < 3) {
                    strI = "0" + strI;
                }
                String strJ = Integer.toString(j);
                while (strJ.length() < 3) {
                    strJ = "0" + strJ;
                }
                String str = strI + "_" + strJ + "_";
                if (i0 == 0) {
                    str = str + "factor";
                } else {
                    str = str + "divisor";
                }
                plotter.writeImage("_MATCH_" + str);
            }
        } catch(Exception e) {}
    }

    private static void debugPlot(int i, int j,
        FixedSizeSortedVector<CObject> vec,
        TwoDFloatArray pyr1, TwoDFloatArray pyr2, float s1, float s2) {

        Image img1 = convertToImage(pyr1);
        Image img2 = convertToImage(pyr2);

        try {
            CorrespondenceList cor = vec.getArray()[0].cCor;

            Image img1Cp = img1.copyImage();
            Image img2Cp = img2.copyImage();
            CorrespondencePlotter plotter = new CorrespondencePlotter(
                img1Cp, img2Cp);
            for (int ii = 0; ii < cor.getPoints1().size(); ++ii) {
                PairInt p1 = cor.getPoints1().get(ii);
                PairInt p2 = cor.getPoints2().get(ii);
                int x1 = Math.round(p1.getX()/s1);
                int y1 = Math.round(p1.getY()/s1);
                int x2 = Math.round(p2.getX()/s2);
                int y2 = Math.round(p2.getY()/s2);
                plotter.drawLineInAlternatingColors(x1, y1, x2, y2, 0);
            }
            String strI = Integer.toString(i);
            while (strI.length() < 3) {
                strI = "0" + strI;
            }
            String strJ = Integer.toString(j);
            while (strJ.length() < 3) {
                strJ = "0" + strJ;
            }
            String str = strI + "_" + strJ + "_";
            plotter.writeImage("_MATCH_" + str);
        } catch(Exception e) {}
    }

    private static void debugPlot2(int i, int j,
        FixedSizeSortedVector<CObject3> vec,
        TwoDFloatArray pyr1, TwoDFloatArray pyr2, float s1, float s2) {

        Image img1 = convertToImage(pyr1);
        Image img2 = convertToImage(pyr2);

        try {
            PairInt[] m1 = vec.getArray()[0].m1;
            PairInt[] m2 = vec.getArray()[0].m2;

            Image img1Cp = img1.copyImage();
            Image img2Cp = img2.copyImage();
            CorrespondencePlotter plotter = new CorrespondencePlotter(
                img1Cp, img2Cp);
            for (int ii = 0; ii < m1.length; ++ii) {
                PairInt p1 = m1[ii];
                PairInt p2 = m2[ii];
                int x1 = Math.round(p1.getX()/s1);
                int y1 = Math.round(p1.getY()/s1);
                int x2 = Math.round(p2.getX()/s2);
                int y2 = Math.round(p2.getY()/s2);
                plotter.drawLineInAlternatingColors(x1, y1, x2, y2, 0);
            }
            String strI = Integer.toString(i);
            while (strI.length() < 3) {
                strI = "0" + strI;
            }
            String strJ = Integer.toString(j);
            while (strJ.length() < 3) {
                strJ = "0" + strJ;
            }
            String str = strI + "_" + strJ + "_";
            plotter.writeImage("_MATCH_" + str);
        } catch(Exception e) {}
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

    private static void debugPrint(List<QuadInt> pairs, int i, int j,
        TwoDFloatArray pyr1, TwoDFloatArray pyr2, float s1, float s2) {

        Image img1 = convertToImage(pyr1);
        Image img2 = convertToImage(pyr2);

        try {

            for (int ii = 0; ii < pairs.size(); ++ii) {
                QuadInt q = pairs.get(ii);
                int x1 = Math.round(q.getA()/s1);
                int y1 = Math.round(q.getB()/s1);
                int x2 = Math.round(q.getC()/s2);
                int y2 = Math.round(q.getD()/s2);
                ImageIOHelper.addPointToImage(x1, y1, img1, 1, 255, 0, 0);
                ImageIOHelper.addPointToImage(x2, y2, img2, 1, 255, 0, 0);
            }
            String strI = Integer.toString(i);
            while (strI.length() < 3) {
                strI = "0" + strI;
            }
            String strJ = Integer.toString(j);
            while (strJ.length() < 3) {
                strJ = "0" + strJ;
            }
            String str = "_pairs_" + strI + "_" + strJ + "_";
            MiscDebug.writeImage(img1, str + "_" + strI);
            MiscDebug.writeImage(img2, str + "_" + strJ);
        } catch(Exception e) {}
    }

    private static List<TIntList> createLabeledLists(
        TIntList keypointsX, TIntList keypointsY,
        List<Set<PairInt>> labeledPoints,
        TObjectIntMap<PairInt> pointLabels) {

        int ns = labeledPoints.size();
        List<TIntList> output = new ArrayList<TIntList>();
        for (int i = 0; i < ns; ++i) {
            output.add(new TIntArrayList());
        }

        for (int i = 0; i < keypointsX.size(); ++i) {
            int x = keypointsX.get(i);
            int y = keypointsY.get(i);
            int label = pointLabels.get(new PairInt(x, y));
            TIntList list = output.get(label);
            assert(list != null);
            list.add(i);
        }

        return output;
    }

    private static List<QuadInt> createPairLabelIndexes(
        int[][] cost, int nBands,
        List<TIntList> pointIndexLists1, TIntList kpX1, TIntList kpY1,
        List<TIntList> pointIndexLists2, TIntList kpX2, TIntList kpY2) {

        int costLimit = Math.round(
            (float)(nBands * 256) * 0.65f);

        int minP1Diff = 3;

        Set<QuadInt> exists = new HashSet<QuadInt>();

        // pairs of idx from set1 and idx from set 2
        Set<PairInt> skip = new HashSet<PairInt>();

        List<QuadInt> pairIndexes = new ArrayList<QuadInt>();

        List<PairInt> pair2Indexes = calculatePairIndexes(
            pointIndexLists2, kpX2, kpY2, minP1Diff);

        for (int ii = 0; ii < pointIndexLists1.size(); ++ii) {
            TIntList kpIndexes1 = pointIndexLists1.get(ii);
            if (kpIndexes1.size() < 2) {
                continue;
            }
            // draw 2 from kpIndexes1
            for (int ii1 = 0; ii1 < kpIndexes1.size(); ++ii1) {
                int idx1 = kpIndexes1.get(ii1);
                int t1X = kpX1.get(idx1);
                int t1Y = kpY1.get(idx1);
                boolean skipIdx1 = false;
                for (int ii2 = 0; ii2 < kpIndexes1.size(); ++ii2) {
                    if (ii1 == ii2) {
                        continue;
                    }
                    int idx2 = kpIndexes1.get(ii2);
                    int t2X = kpX1.get(idx2);
                    int t2Y = kpY1.get(idx2);

                    if ((t1X == t2X && t1Y == t2Y)) {
                        continue;
                    }

                    int diffX = t1X - t2X;
                    int diffY = t1Y - t2Y;
                    int distSq = diffX * diffX + diffY * diffY;
                    //if (distSq > limitSq) {
                    //    continue;
                    //}
                    if ((distSq < minP1Diff * minP1Diff)) {
                        continue;
                    }

                    for (PairInt p2Index : pair2Indexes) {
                        int idx3 = p2Index.getX();
                        int idx4 = p2Index.getY();

                        PairInt p13 = new PairInt(idx1, idx3);
                        if (skip.contains(p13)) {
                            skipIdx1 = true;
                            break;
                        }
                        PairInt p24 = new PairInt(idx2, idx4);
                        if (skip.contains(p24)) {
                            continue;
                        }
                        int c13 = cost[idx1][idx3];
                        // if idx1 and idx3 cost is above limit, skip
                        if (c13 > costLimit) {
                            skip.add(p13);
                            skipIdx1 = true;
                            break;
                        }
                        int c24 = cost[idx2][idx4];
                        if (c24 > costLimit) {
                            skip.add(p24);
                            continue;
                        }

                        QuadInt q = new QuadInt(idx1, idx2, idx3, idx4);
                        QuadInt qChk = new QuadInt(idx2, idx1, idx4, idx3);

                        if (exists.contains(q) || exists.contains(qChk)) {
                            continue;
                        }

                        /*
                        int s1X = kpX2.get(idx3);
                        int s1Y = kpY2.get(idx3);
                        int s2X = kpX2.get(idx4);
                        int s2Y = kpY2.get(idx4);

                        int diffX2 = s1X - s2X;
                        int diffY2 = s1Y - s2Y;
                        int distSq2 = diffX2 * diffX2 + diffY2 * diffY2;
                        //if (distSq2 > limitSq) {
                        //    continue;
                        //}
                        if ((distSq2 < minP1Diff * minP1Diff)) {
                            continue;
                        }
                        */

                        pairIndexes.add(q);
                        exists.add(q);
                    }
                    if(skipIdx1) {
                        break;
                    }
                }
            }
        }

        return pairIndexes;
    }

    private static class TransformationSort implements Comparator<TransformationParameters> {

        @Override
        public int compare(TransformationParameters o1,
            TransformationParameters o2) {

            float s1 = Float.parseFloat(String.format("%.3f",
                Math.abs(1 - o1.getScale())));
            float s2 = Float.parseFloat(String.format("%.3f",
                Math.abs(1 - o2.getScale())));
            float r1 = Float.parseFloat(String.format("%.3f",
                o1.getRotationInDegrees()));
            float r2 = Float.parseFloat(String.format("%.3f",
                o2.getRotationInDegrees()));
            float t1x = o1.getTranslationX();
            float t2x = o2.getTranslationX();
            float t1y = o1.getTranslationY();
            float t2y = o2.getTranslationY();

            int comp = Float.compare(s1, s2);
            if (comp != 0) {
                return comp;
            }
            comp = Float.compare(r1, r2);
            if (comp != 0) {
                return comp;
            }
            comp = Float.compare(t1x, t2x);
            if (comp != 0) {
                return comp;
            }
            comp = Float.compare(t1y, t2y);
            return comp;
        }

    }

    private static void printSumPatchDifference(TwoDFloatArray octaveImg1,
        TwoDFloatArray octaveImg2,
        PairIntArray a1, PairIntArray tr1,
        Set<PairInt> set2, TransformationParameters params) {

        int ns2 = set2.size();

        Image img1 = convertToImage(octaveImg1);
        Image img2 = convertToImage(octaveImg2);

        //MiscDebug.writeImage(img1, "_" + MiscDebug.getCurrentTimeFormatted());
        //MiscDebug.writeImage(img2, "_" + MiscDebug.getCurrentTimeFormatted());

        ImageIOHelper.addCurveToImage(set2, img2, 2, 0, 0, 255);

        ImageIOHelper.addCurveToImage(a1, img1, 2, 255, 0, 0);

        for (int k = 0; k < tr1.getN(); ++k) {
            int x2 = tr1.getX(k);
            int y2 = tr1.getY(k);
            if (!set2.contains(new PairInt(x2, y2))) {
                continue;
            }
            if (y2 > (octaveImg2.a.length - 1)
                || x2 > (octaveImg2.a[0].length - 1)) {
                continue;
            }
            int x1 = a1.getX(k);
            int y1 = a1.getY(k);
            if (y1 > (octaveImg1.a.length - 1)
                || x1 > (octaveImg1.a[0].length - 1)) {
                continue;
            }

            ImageIOHelper.addPointToImage(x2, y2, img2, 1, 0, 255, 0);
        }
        MiscDebug.writeImage(img1, "_1_" + MiscDebug.getCurrentTimeFormatted());
        MiscDebug.writeImage(img2, "_2_" + MiscDebug.getCurrentTimeFormatted());
    }

    //sumDiff, n
    private static double[] sumPatchDifference(TwoDFloatArray octaveImg1,
        TwoDFloatArray octaveImg2,
        PairIntArray a1, PairIntArray tr1,
        Set<PairInt> set2, TransformationParameters params) {

        int ns2 = set2.size();

        // scaled intersection of points
        TDoubleList v1 = new TDoubleArrayList();
        TDoubleList v2 = new TDoubleArrayList();

        for (int k = 0; k < tr1.getN(); ++k) {
            int x2 = tr1.getX(k);
            int y2 = tr1.getY(k);
            if (!set2.contains(new PairInt(x2, y2))) {
                continue;
            }

            if (y2 > (octaveImg2.a.length - 1)
                || x2 > (octaveImg2.a[0].length - 1)) {
                continue;
            }
            int x1 = a1.getX(k);
            int y1 = a1.getY(k);
            if (y1 > (octaveImg1.a.length - 1)
                || x1 > (octaveImg1.a[0].length - 1)) {
                continue;
            }
            v1.add(octaveImg1.a[y1][x1]);
            v2.add(octaveImg2.a[y2][x2]);
        }

        // count should have at least the 2 points of transformation
if (v1.size() == 0) {
    printSumPatchDifference(octaveImg1, octaveImg2,
    a1, tr1, set2, params);
    int z = 1;
}
        if (v1.size() == 0) {
            return null;
        }

        double n = v1.size();

        double sumDiff = 0;

        for (int i = 0; i < v1.size(); ++i) {
            double a = v1.get(i);
            double b = v2.get(i);
            
            // since this is a polar theta 0 to 255 image, need to
            //account for wrap around
            if (a > b) {
                // add a phase to next value if it's closer to current with addition
                if ((a - b) > Math.abs(a - (b + 255))) {
                    b += 255;
                }
            } else if (b > a) {
                // add a phase to next value if it's closer to current with addition
                if ((b - a) > Math.abs(b - (a + 255))) {
                    a += 255;
                }
            }
            double diff = Math.abs(a - b);
            sumDiff += (diff * diff);
            //sumDiff += diff;
        }

        sumDiff /= n;

        //System.out.println("sumDiff=" + sumDiff);
        return new double[]{sumDiff, n};
    }

    private static double[] sumPatchDifference(TwoDFloatArray octaveImg1,
        TwoDFloatArray octaveImg2, Set<PairInt> set1,
        Set<PairInt> set2,
        TransformationParameters params, float scale1, float scale2,
        PairInt p2_1, PairInt p2_2, int nMaxPatchPixels) {

        //TODO: this needs a better comparison method for the intersection
        // of these regions.
        // the current pixel by pixel needs to be enlarged to include neighbors
        // for each comparison...

        int ns2 = set2.size();
        if (ns2 > nMaxPatchPixels) {
            //TODO: consider delaying the trim until intersection of transformed
            set2 = reduceSet(set2, p2_1, p2_2, nMaxPatchPixels);
        }

        Transformer transformer = new Transformer();

        PairIntArray list2 = Misc.convertWithoutOrder(set2);
        PairIntArray tr2
            = transformer.applyTransformation(params, list2);
        assert (tr2.getN() == list2.getN());

        // scaled intersection of points
        double mean1 = 0;
        double mean2 = 0;
        TDoubleList v1 = new TDoubleArrayList();
        TDoubleList v2 = new TDoubleArrayList();

        /*
        will either subtract mean then divide by stdev
        or will subtract mean and divide by (sqrt(2)/mean)
        */

        for (int k = 0; k < tr2.getN(); ++k) {

            int x1Tr = Math.round(tr2.getX(k) / scale1);
            int y1Tr = Math.round(tr2.getY(k) / scale1);

            if (!set1.contains(new PairInt(tr2.getX(k),
                tr2.getY(k)))) {
                continue;
            }

            if (y1Tr > (octaveImg1.a.length - 1)
                || x1Tr > (octaveImg1.a[0].length - 1)) {
                continue;
            }
            int x2 = Math.round(list2.getX(k) / scale2);
            int y2 = Math.round(list2.getY(k) / scale2);
            if (y2 > (octaveImg2.a.length - 1)
                || x2 > (octaveImg2.a[0].length - 1)) {
                continue;
            }
            v1.add(octaveImg1.a[y1Tr][x1Tr]);
            v2.add(octaveImg2.a[y2][x2]);
            mean1 += octaveImg1.a[y1Tr][x1Tr];
            mean2 += octaveImg2.a[y2][x2];
        }

        // count should have at least the 2 points of transformation
        assert (v1.size() > 0);
        double n = v1.size();

        mean1 /= n;
        mean2 /= n;

        // for the segmentation patches, many of the pixels will be
        // similar values within the set, so divifing by standard deviation
        // to leave the most significantly different pixels at high values
        // will tend to make patches look more similar.
        // trying with only mean subtraction first to  reduce grey illumination
        // effects

        /*
        // calculate standard deviation
        double sumStDv1 = 0;
        double sumStDv2 = 0;
        for (int i = 0; i < v1.size(); ++i) {
            double d1 = v1.get(i) - mean1;
            double d2 = v2.get(i) - mean2;
            sumStDv1 = (d1 * d1);
            sumStDv2 = (d2 * d2);
        }
        sumStDv1 = Math.sqrt(sumStDv1/(v1.size() - 1.0f));
        sumStDv2 = Math.sqrt(sumStDv2/(v1.size() - 1.0f));
        */
        double sumDiff = 0;

        for (int i = 0; i < v1.size(); ++i) {
            double a = (v1.get(i) - mean1);
            double b = (v2.get(i) - mean2);
            //double a = (v1.get(i) - mean1) /sumStDv1;
            //double b = (v2.get(i) - mean2) /sumStDv2;
            //double a = v1.get(i);
            //double b = v2.get(i);
            double diff = Math.abs(a - b);
            sumDiff += (diff * diff);
            //sumDiff += diff;
        }

        sumDiff /= n;

        //TODO:
        // then divide by max possible value to normalize
        //   this should be calculated.  will use a
        //   very rough guesstimate for now.
        //    127*n/n
        //sumDiff /= 127.;
        //System.out.println("sumDiff=" + sumDiff);
        return new double[]{sumDiff, n};
    }

    private static double[] sumPatchDifference(TwoDFloatArray octaveImg1,
        TwoDFloatArray octaveImg2,
        PairIntArray a1, PairIntArray tr1, Set<PairInt> set2,
        TransformationParameters params, float scale1, float scale2,
        int nMaxPatchPixels) {

        // scaled intersection of points
        double mean1 = 0;
        double mean2 = 0;
        TDoubleList v1 = new TDoubleArrayList();
        TDoubleList v2 = new TDoubleArrayList();

        /*
        will either subtract mean then divide by stdev
        or will subtract mean and divide by (sqrt(2)/mean)
        */

        for (int k = 0; k < a1.getN(); ++k) {
            int x2Tr = Math.round(tr1.getX(k)/ scale2);
            int y2Tr = Math.round(tr1.getY(k)/ scale2);
            if (!set2.contains(new PairInt(tr1.getX(k),
                tr1.getY(k)))) {
                continue;
            }
            int x1 = Math.round(a1.getX(k) / scale1);
            int y1 = Math.round(a1.getY(k) / scale1);
            if (y1 > (octaveImg1.a.length - 1)
                || x1 > (octaveImg1.a[0].length - 1)) {
                continue;
            }
            if (x2Tr > (octaveImg2.a.length - 1)
                || y2Tr > (octaveImg2.a[0].length - 1)) {
                continue;
            }
            v1.add(octaveImg1.a[y1][x1]);
            v2.add(octaveImg2.a[y2Tr][x2Tr]);
            mean1 += octaveImg1.a[y1][x1];
            mean2 += octaveImg2.a[y2Tr][x2Tr];
        }

        // count should have at least the 2 points of transformation
        if (v1.isEmpty()) {
            return null;
        }

        double n = v1.size();

        mean1 /= n;
        mean2 /= n;

        // for the segmentation patches, many of the pixels will be
        // similar values within the set, so divifing by standard deviation
        // to leave the most significantly different pixels at high values
        // will tend to make patches look more similar.
        // trying with only mean subtraction first to  reduce grey illumination
        // effects

        /*
        // calculate standard deviation
        double sumStDv1 = 0;
        double sumStDv2 = 0;
        for (int i = 0; i < v1.size(); ++i) {
            double d1 = v1.get(i) - mean1;
            double d2 = v2.get(i) - mean2;
            sumStDv1 = (d1 * d1);
            sumStDv2 = (d2 * d2);
        }
        sumStDv1 = Math.sqrt(sumStDv1/(v1.size() - 1.0f));
        sumStDv2 = Math.sqrt(sumStDv2/(v1.size() - 1.0f));
        */
        double sumDiff = 0;

        for (int i = 0; i < v1.size(); ++i) {
            //double a = (v1.get(i) - mean1) /sumStDv1;
            //double b = (v2.get(i) - mean2) /sumStDv2;
            double a = v1.get(i);
            double b = v2.get(i);
            double diff = Math.abs(a - b);
            //sumDiff += (diff * diff);
            sumDiff += diff;
        }

        sumDiff /= n;

        //TODO:
        // then divide by max possible value to normalize
        //   this should be calculated.  will use a
        //   very rough guesstimate for now.
        //    127*n/n
        //sumDiff /= 127.;
        //System.out.println("sumDiff=" + sumDiff);
        return new double[]{sumDiff, n};
    }

    private static void printSumPatchDifference(TwoDFloatArray octaveImg1,
        TwoDFloatArray octaveImg2,
        PairIntArray a1, PairIntArray tr1, Set<PairInt> set2,
        TransformationParameters params, float scale1, float scale2,
        int nMaxPatchPixels) {

        //Image img1 = convertToImage(octaveImg1);
        Image img2 = convertToImage(octaveImg2);

        ImageIOHelper.addCurveToImage(set2, img2, 2, 255, 0, 0);

        for (int k = 0; k < a1.getN(); ++k) {
            int x2Tr = Math.round(tr1.getX(k)/ scale2);
            int y2Tr = Math.round(tr1.getY(k)/ scale2);
            if (!set2.contains(new PairInt(tr1.getX(k),
                tr1.getY(k)))) {
                continue;
            }
            int x1 = Math.round(a1.getX(k) / scale1);
            int y1 = Math.round(a1.getY(k) / scale1);
            if (y1 > (octaveImg1.a.length - 1)
                || x1 > (octaveImg1.a[0].length - 1)) {
                continue;
            }
            if (x2Tr > (octaveImg2.a.length - 1)
                || y2Tr > (octaveImg2.a[0].length - 1)) {
                continue;
            }

            ImageIOHelper.addPointToImage(x2Tr, y2Tr, img2, 1, 0, 255, 0);
        }
        MiscDebug.writeImage(img2, "_" +
            MiscDebug.getCurrentTimeFormatted());
    }

    private static double[] sumKeypointDistanceDifference(
        TIntList a2Indexes,
        PairIntArray tr2,
        TIntList kpX2, TIntList kpY2,
        NearestNeighbor2D nn, TransformationParameters params,
        int maxX, int maxY, int pixTolerance,
        double maxDist, int[] m1x, int[] m1y, int[] m2x, int[] m2y) {

        double sum2 = 0;
        int mCount = 0;

        for (int k = 0; k < tr2.getN(); ++k) {
            int x2Tr = tr2.getX(k);
            int y2Tr = tr2.getY(k);
            int idx2 = a2Indexes.get(k);

            Set<PairInt> nearest = null;
            if ((x2Tr >= 0) && (y2Tr >= 0)
                && (x2Tr <= (maxX + pixTolerance))
                && (y2Tr <= (maxY + pixTolerance))) {
                nearest = nn.findClosest(x2Tr, y2Tr, pixTolerance);
            }
            double minDist = Double.MAX_VALUE;
            PairInt minDistP1 = null;
            if (nearest != null && !nearest.isEmpty()) {
                for (PairInt p11 : nearest) {
                    double dist = distance(x2Tr, y2Tr, p11);
                    if (dist < minDist) {
                        minDist = dist;
                        minDistP1 = p11;
                    }
                }
            }
            if (minDistP1 != null) {
                double dist = minDist;
                double distNorm = dist / maxDist;
                sum2 += distNorm;
                m2x[mCount] = kpX2.get(idx2);
                m2y[mCount] = kpY2.get(idx2);
                m1x[mCount] = minDistP1.getX();
                m1y[mCount] = minDistP1.getY();
                mCount++;
            } else {
                sum2 += 1;
            }
        } // end loop over trnsformed set 2

        return new double[]{sum2, mCount};
    }

    private static void debugPrint(TwoDFloatArray pyr1,
        TwoDFloatArray pyr2, TIntList kpX1, TIntList kpY1,
        TIntList kpX2, TIntList kpY2,
        float scale1, float scale2, int img1Idx, int img2Idx) {

        Image img1 = convertToImage(pyr1);
        Image img2 = convertToImage(pyr2);

        try {

            for (int i = 0; i < kpX1.size(); ++i) {
                int x1 = (int)(kpX1.get(i)/scale1);
                int y1 = (int)(kpY1.get(i)/scale1);
                ImageIOHelper.addPointToImage(x1, y1, img1, 1, 255, 0, 0);
            }
            for (int i = 0; i < kpX2.size(); ++i) {
                int x2 = (int)(kpX2.get(i)/scale2);
                int y2 = (int)(kpY2.get(i)/scale2);
                ImageIOHelper.addPointToImage(x2, y2, img2, 1, 255, 0, 0);
            }
            String strI = Integer.toString(img1Idx);
            while (strI.length() < 3) {
                strI = "0" + strI;
            }
            String strJ = Integer.toString(img2Idx);
            while (strJ.length() < 3) {
                strJ = "0" + strJ;
            }
            String str = "_kp_" + strI + "_" + strJ + "_";
            MiscDebug.writeImage(img1, str + "_i");
            MiscDebug.writeImage(img2, str + "_j");
        } catch(Exception e) {}

    }

    private static void debugPrint(TwoDFloatArray pyr1,
        TwoDFloatArray pyr2, TIntList kpX1, TIntList kpY1,
        TIntList kpX2, TIntList kpY2, int img1Idx, int img2Idx) {

        Image img1 = convertToImage(pyr1);
        Image img2 = convertToImage(pyr2);

        try {

            for (int i = 0; i < kpX1.size(); ++i) {
                int x1 = kpX1.get(i);
                int y1 = kpY1.get(i);
                ImageIOHelper.addPointToImage(x1, y1, img1, 1, 255, 0, 0);
            }
            for (int i = 0; i < kpX2.size(); ++i) {
                int x2 = kpX2.get(i);
                int y2 = kpY2.get(i);
                ImageIOHelper.addPointToImage(x2, y2, img2, 1, 255, 0, 0);
            }
            String strI = Integer.toString(img1Idx);
            while (strI.length() < 3) {
                strI = "0" + strI;
            }
            String strJ = Integer.toString(img2Idx);
            while (strJ.length() < 3) {
                strJ = "0" + strJ;
            }
            String str = "_kp_" + strI + "_" + strJ + "_";
            MiscDebug.writeImage(img1, str + "_i");
            MiscDebug.writeImage(img2, str + "_j");
        } catch(Exception e) {}

    }

    private static Set<PairInt> reduceSet(Set<PairInt> set2,
        PairInt p2_1, PairInt p2_2, int nMaxPatchPixels) {

        float[] minDist = new float[set2.size()];
        PairInt[] points = new PairInt[minDist.length];

        int count = 0;
        for (PairInt p : set2) {
            double d1 = distance(p, p2_1);
            double d2 = distance(p, p2_2);
            minDist[count] = (float)Math.min(d1, d2);
            points[count] = p;
            count++;
        }

        QuickSort.sortBy1stArg(minDist, points);

        Set<PairInt> output = new HashSet<PairInt>();
        for (int i = 0; i < nMaxPatchPixels; ++i) {
            PairInt p = points[i];
            output.add(p);
        }

        return output;
    }

    private static TFloatList extractScales(List<TFloatList> scalesList) {

        TFloatList scales = new TFloatArrayList();

        for (int i = 0; i < scalesList.size(); ++i) {
            scales.add(scalesList.get(i).get(0));
        }

        return scales;
    }

    private static double[] sumKeypointDescAndDist(
        int[][] cost, int nBands,
        TIntList a1Indexes, PairIntArray tr1, TIntList kpX1, TIntList kpY1,
        NearestNeighbor2D nn2, TObjectIntMap<PairInt> p2KPIndexMap,
        TransformationParameters params, int maxX2, int maxY2,
        int pixTolerance, double maxDist,
        PairInt[] m1, PairInt[] m2) {

        double sumDesc = 0;
        double sumDist = 0;
        int count = 0;

        double maxDesc = nBands * 256.;

        for (int k = 0; k < tr1.getN(); ++k) {
            int x1Tr = tr1.getX(k);
            int y1Tr = tr1.getY(k);
            int idx1 = a1Indexes.get(k);

            Set<PairInt> nearest = null;
            if ((x1Tr >= 0) && (y1Tr >= 0)
                && (x1Tr <= (maxX2 + pixTolerance))
                && (y1Tr <= (maxY2 + pixTolerance))) {
                nearest = nn2.findClosest(x1Tr, y1Tr, pixTolerance);
            }
            int minC = Integer.MAX_VALUE;
            PairInt minCP2 = null;
            int minIdx2 = 0;
            if (nearest != null && !nearest.isEmpty()) {
                for (PairInt p2 : nearest) {
                    int idx2 = p2KPIndexMap.get(p2);
                    int c = cost[idx1][idx2];
                    if (c < minC) {
                        minC = c;
                        minCP2 = p2;
                        minIdx2 = idx2;
                    }
                }
            }

            if (minCP2 != null) {
                double scoreNorm = (nBands * 256 - minC) / maxDesc;
                double costNorm = 1. - scoreNorm;
                sumDesc += costNorm;

                double dist = distance(x1Tr, y1Tr, minCP2);
                double distNorm = dist / maxDist;
                sumDist += distNorm;

                m1[count] = new PairInt(kpX1.get(idx1), kpY1.get(idx1));
                m2[count] = minCP2;
                count++;

            } else {
                sumDesc += 1;
                sumDist += 1;
            }
        }

        return new double[]{sumDesc, sumDist, count};
    }

    private static double[] sumKeypointDescAndDist(
        int[][] cost, int nBands,
        TIntList a1Indexes, PairIntArray tr1, TIntList kpX1, TIntList kpY1,
        NearestNeighbor2D nn2, TObjectIntMap<PairInt> p2KPIndexMap,
        int maxX2, int maxY2, int pixTolerance, double maxDist,
        int[] m1, int[] m2) {

        double sumDesc = 0;
        double sumDist = 0;
        int count = 0;

        double maxDesc = nBands * 256.;

        //best first match, after nearest neighbors
        // TODO: consider optimal bipartite matching when have an
        //       implementation of multi-level-buckets
        float[] costA = new float[tr1.getN()];
        float[] costDesc = new float[tr1.getN()];
        float[] costDist = new float[tr1.getN()];
        int[] indexes = new int[tr1.getN()];
            
        for (int k = 0; k < tr1.getN(); ++k) {
            int x1Tr = tr1.getX(k);
            int y1Tr = tr1.getY(k);
            int idx1 = a1Indexes.get(k);

            Set<PairInt> nearest = null;
            if ((x1Tr >= 0) && (y1Tr >= 0)
                && (x1Tr <= (maxX2 + pixTolerance))
                && (y1Tr <= (maxY2 + pixTolerance))) {
                nearest = nn2.findClosest(x1Tr, y1Tr, pixTolerance);
            }
            int minC = Integer.MAX_VALUE;
            PairInt minCP2 = null;
            int minIdx2 = 0;
            if (nearest != null && !nearest.isEmpty()) {
                for (PairInt p2 : nearest) {
                    int idx2 = p2KPIndexMap.get(p2);
                    int c = cost[idx1][idx2];
                    if (c < minC) {
                        minC = c;
                        minCP2 = p2;
                        minIdx2 = idx2;
                    }
                }
            }

            if (minCP2 != null) {
                double scoreNorm = (nBands * 256 - minC) / maxDesc;
                double costNorm = 1. - scoreNorm;
                sumDesc += costNorm;

                double dist = distance(x1Tr, y1Tr, minCP2);
                double distNorm = dist / maxDist;
                sumDist += distNorm;

                m1[count] = idx1;
                m2[count] = minIdx2;
                costA[count] = (float)(costNorm + distNorm);
                costDesc[count] = (float)costNorm;
                costDist[count] = (float)distNorm;
                indexes[count] = count;
                count++;

            } else {
                sumDesc += 1;
                sumDist += 1;
            }
        }
        
        if (count > 1) {
            costA = Arrays.copyOf(costA, count);
            indexes = Arrays.copyOf(indexes, count);
            QuickSort.sortBy1stArg(costA, indexes);
            TIntSet set1 = new TIntHashSet();
            TIntSet set2 = new TIntHashSet();
            List<PairInt> matched = new ArrayList<PairInt>();
            TIntList idxs = new TIntArrayList();
            for (int i = 0; i < count; ++i) {
                int idx = indexes[i];
                int idx1 = m1[idx];
                int idx2 = m2[idx];
                if (set1.contains(idx1) || set2.contains(idx2)) {
                    continue;
                }
                idxs.add(idx);
                matched.add(new PairInt(idx1, idx2));
                set1.add(idx1);
                set2.add(idx2);
            }
            int nRedundant = count - matched.size();
            if (nRedundant > 0) {
                sumDesc = 0;
                sumDist = 0;
                for (int i = 0; i < matched.size(); ++i) {
                    m1[i] = matched.get(i).getX();
                    m2[i] = matched.get(i).getY();
                    int idx = idxs.get(i);
                    sumDesc += costDesc[idx];
                    sumDist += costDist[idx];
                }
                sumDesc += (tr1.getN() - matched.size());
                sumDist += (tr1.getN() - matched.size());
                count = matched.size();
            }
        }

        return new double[]{sumDesc, sumDist, count};
    }

    private static PairIntArray aggregateLabels(
        List<Set<PairInt>> labeledPoints,
        TIntList kpX, TIntList kpY,
        TObjectIntMap<PairInt> pointLabels) {

        Set<PairInt> out = new HashSet<PairInt>();

        for (int i = 0; i < kpX.size(); ++i) {
            PairInt p = new PairInt(kpX.get(i), kpY.get(i));
            int label = pointLabels.get(p);
            for (PairInt p2 : labeledPoints.get(label)) {
                out.add(p2);
            }
        }

        return Misc.convertWithoutOrder(out);
    }

    private static Set<PairInt> extractLabeledRegions(
        PairInt[] points, TObjectIntMap<PairInt> pointLabels,
        List<Set<PairInt>> labeledPoints) {

        Set<PairInt> out = new HashSet<PairInt>();

        for (int i = 0; i < points.length; ++i) {
            PairInt p = points[i];
            int label = pointLabels.get(p);
            out.addAll(labeledPoints.get(label));
        }

        return out;
    }

    private static int count(Set<PairInt> set, float scale) {

        Set<PairInt> out = new HashSet<PairInt>();
        for (PairInt p : set) {
            int x = Math.round(p.getX()/scale);
            int y = Math.round(p.getY()/scale);
            out.add(new PairInt(x, y));
        }

        return out.size();
    }

    private static List<PairInt> calculatePairIndexes(
        List<TIntList> pointIndexLists2, TIntList kpX2,
        TIntList kpY2, int minPDiff) {

        List<PairInt> pairIndexes = new ArrayList<PairInt>();

        Set<PairInt> exists = new HashSet<PairInt>();

        // draw 2 pairs from other dataset
        for (int jj = 0; jj < pointIndexLists2.size(); ++jj) {
            TIntList kpIndexes2 = pointIndexLists2.get(jj);
            if (kpIndexes2.size() < 2) {
                continue;
            }

            // draw 2 from kpIndexes2
            for (int jj1 = 0; jj1 < kpIndexes2.size(); ++jj1) {
                int idx3 = kpIndexes2.get(jj1);
                int s1X = kpX2.get(idx3);
                int s1Y = kpY2.get(idx3);
                for (int jj2 = 0; jj2 < kpIndexes2.size(); ++jj2) {
                    if (jj1 == jj2) {
                        continue;
                    }
                    int idx4 = kpIndexes2.get(jj2);
                    int s2X = kpX2.get(idx4);
                    int s2Y = kpY2.get(idx4);
                    if ((s1X == s2X && s1Y == s2Y)) {
                        continue;
                    }

                    PairInt q = new PairInt(idx3, idx4);

                    if (exists.contains(q)) {
                        continue;
                    }

                    int diffX2 = s1X - s2X;
                    int diffY2 = s1Y - s2Y;
                    int distSq2 = diffX2 * diffX2 + diffY2 * diffY2;
                    //if (distSq2 > limitSq) {
                    //    continue;
                    //}
                    if ((distSq2 < minPDiff * minPDiff)) {
                        continue;
                    }

                    pairIndexes.add(q);
                    exists.add(q);
                }
            }
        }

        return pairIndexes;
    }

    private static PairIntArray trimToImageBounds(
        TwoDFloatArray octaveImg, PairIntArray a) {

        int n0 = octaveImg.a.length;
        int n1 = octaveImg.a[0].length;

        PairIntArray b = new PairIntArray(a.getN());

        for (int i = 0; i < a.getN(); ++i) {
            int x = a.getX(i);
            int y = a.getY(i);
            if (x < 0 || x > (n1 - 1)) {
                continue;
            } else if(y < 0 || y > (n0 - 1)) {
                continue;
            }
            b.add(x, y);
        }

        return b;
    }

}
