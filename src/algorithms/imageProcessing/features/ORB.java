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
import algorithms.compGeometry.FurthestPair;
import algorithms.imageProcessing.ColorHistogram;
import algorithms.imageProcessing.SIGMA;
import algorithms.imageProcessing.matching.ORBMatcher;
import algorithms.imageProcessing.matching.PartialShapeMatcher;
import algorithms.imageProcessing.matching.PartialShapeMatcher.Result;
import algorithms.misc.Misc;
import algorithms.misc.MiscDebug;
import algorithms.search.NearestNeighbor2D;
import algorithms.util.CorrespondencePlotter;
import algorithms.util.OneDIntArray;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.QuadInt;
import algorithms.util.TwoDFloatArray;
import algorithms.util.TwoDIntArray;
import algorithms.util.VeryLongBitString;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.TFloatList;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TFloatIntMap;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.hash.TFloatIntHashMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import gnu.trove.set.TDoubleSet;
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


    private static TIntList calculateOrientations360_0(ORB orb, 
        Set<PairInt> labeledPoints) {
    
        List<Set<PairInt>> list = new ArrayList<Set<PairInt>>();
        list.add(labeledPoints);
    
        List<TIntList> orList = calculateOrientations360(orb, 
            list);
    
        return orList.get(0);
    }
    
    private static List<TIntList> calculateOrientations360(ORB orb, 
        List<Set<PairInt>> labeledPoints) {
        
        TIntList kp0s = new TIntArrayList();
        TIntList kp1s = new TIntArrayList();
        
        for (int i = 0; i < labeledPoints.size(); ++i) {
            for (PairInt p : labeledPoints.get(i)) {
                kp0s.add(p.getY());
                kp1s.add(p.getX());
            }
        }
        
        float[][] octaveImage = orb.getPyramidImages().get(0).a;
        TDoubleList or = orb.cornerOrientations(octaveImage,
            kp0s, kp1s);
        
        // convert to degrees
        TIntList orD = new TIntArrayList(or.size());
        for (int i = 0; i < or.size(); ++i) {
            double d = or.get(i) * 180./Math.PI;
            if (d < 0) {
                d += 360;
            } else if (d > 359) {
                d -= 360;
            }
            orD.add((int)Math.round(d));
        }
        
        List<TIntList> orList = new ArrayList<TIntList>(labeledPoints.size());
        int count = 0;
        for (int i = 0; i < labeledPoints.size(); ++i) {
            int n = labeledPoints.get(i).size();
            TIntList list = new TIntArrayList(n);
            orList.add(list);
            list.addAll(orD.subList(count, count + n));
            count += n;
        }

        return orList;        
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

    /**
     * @return the descrChoice
     */
    public DescriptorChoice getDescrChoice() {
        return descrChoice;
    }

    /**
     * @param descrChoice the descrChoice to set
     */
    public void setDescrChoice(DescriptorChoice descrChoice) {
        this.descrChoice = descrChoice;
    }

    public static enum DescriptorChoice {
        NONE, GREYSCALE, HSV, ALT
    }
    private DescriptorChoice descrChoice = DescriptorChoice.GREYSCALE;

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
        public VeryLongBitString[] descriptors;
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




}
