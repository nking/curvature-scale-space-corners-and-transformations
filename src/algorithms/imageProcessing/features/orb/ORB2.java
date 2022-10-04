package algorithms.imageProcessing.features.orb;

import algorithms.QuickSort;
import algorithms.VeryLongBitString;
import algorithms.imageProcessing.*;
import algorithms.matrix.MatrixUtil;
import algorithms.misc.MiscDebug;
import algorithms.util.PairInt;
import algorithms.util.TwoDFloatArray;
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
import algorithms.imageProcessing.features.orb.ORB.Descriptors;

import java.util.*;
import java.util.logging.Logger;

/**
 * An implementation of "ORB: an efficient alternative to SIFT or SURF"
 * (paper reference is Rublee, Rabaud, Konolige, and Bradski, 2011).
 * http://www.vision.cs.chubu.ac.jp/CV-R/pdf/Rublee_iccv2011.pdf
 *
 * excerpts from the paper:
 * ORB is a very fast binary descriptor based on BRIEF, 
 * which is rotation invariant and resistant to noise.
 * It builds on the well-known FAST keypoint detector and the recently-developed 
 * BRIEF descriptor [6]; for this reason we call it 
 * ORB (Oriented FAST and Rotated BRIEF).
 * 
 * FAST and its variants are the method of choice for finding keypoints in 
 * real-time systems that match visual features, for example, Parallel Tracking 
 * and Mapping. It is efficient and finds reasonable corner keypoints, although 
 * it must be augmented with pyramid schemes for scale, and in our case, 
 * a Harris corner filter [11] to reject edges and provide a reasonable score.
 * 
 * ORB orientation operator is built using the centroid technique of Rosin.
 * 
   Descriptor BRIEF:
   *   is a recent feature descriptor (a bit string description of an image 
   *   patch) that uses simple binary 
   *   tests between pixels in a smoothed image patch. Its performance is similar
   *   to SIFT in many respects, including robustness to lighting, blur, and 
   *   perspective distortion. However, it is very sensitive to in-plane 
   *   rotation.
 
 *     BRIEF grew out of research that uses binary tests to train a set of 
 *     classification trees [4]. Once trained on a set of 500 or so typical 
 *     keypoints, the trees can be used to return a signature for any arbitrary 
 *     keypoint [5]. In a similar manner, we look for the tests least sensitive 
 *     to orientation...NOTE that PCA and SIFT are computationally intensive.
 * 
 *     Visual vocabulary methods [21, 27] use offline clustering to find exemplars 
 *     that are uncorrelated and can be used in matching. These techniques might 
 *     also be useful in finding uncorrelated binary tests.
 * 
 * FAST detector: 
 *    the intensity threshold between the center pixel and those in a circular 
 *    ring about the center (e.g. a circular radius of 9).
 *    FAST does not produce a measure of cornerness, and it has large responses 
 *    along edges. We employ a Harris corner measure [11] to order the FAST 
 *    keypoints. For a target number N of keypoints, we first set the threshold 
 *    low enough to get more than N keypoints, then order them according to the 
 *    Harris measure, and pick the top N points.
 *    FAST is used on each image in a scaled image pyramid.
 * 
 * Orientation:
 *     corner orientation, the intensity centroid [22]. The intensity centroid 
 *     assumes that a corner’s intensity is offset from its center, and this vector 
 *     may be used to impute an orientation. Rosin defines the moments of a patch as:
 * 
 *    m_p_q = moments = summation_over_x_and_y( (x^p) * (y^q) * I(x,y) 
 *       where I(x,y) is the pixel intensity.
 * 
 *    C = the centroid = ( m_1_0/m_0_0, m_0_1/m_0_0)
 *    
 *    orientation = theta = atan2( m_0_1, m_1_0)
 * 
 *     
 * ---------
 * The implementation below is adapted from the scipy implementation which has
 * the following copyright:* 
* 
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
public class ORB2 {

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
    protected StructureTensorR[] tensorComponents = null;
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

    // combinedKeypointList holds row number, column number pairs
    private List<PairInt> combinedKeypointList = null;
    private TDoubleList combinedOrientationList = null;
    private TFloatList combinedHarrisResponseList = null;
    private Descriptors combinedDescriptors = null;

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

    protected boolean extractOnlyHarrisKeypoints = true;

    protected int nPyramidB = 3;

    protected final GreyscaleImage img;

    /**
     * Still testing the class, there may be bugs present.
     * @param image
     * @param nKeypoints
     */
    public ORB2(GreyscaleImage image, int nKeypoints) {

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
    /**
     * default is 0.08f.
     * For a low resolution image, may want to use 0.01f.
     * @param threshold
     */
    public void overrideFastThreshold(float threshold) {
        this.fastThreshold = threshold;
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

    private void initMasks() {

        OFAST_MASK = new int[31][31];
        int i;
        for (i = 0; i < 31; ++i) {
            OFAST_MASK[i] = new int[31];
        }

        OFAST_UMAX = new int[]{15, 15, 15, 15, 14, 14, 14, 13, 13, 12, 11, 10,
            9, 8, 6, 3};

        int absI;
        int j;
        for (i = -15; i < 16; ++i) {
            absI = Math.abs(i);
            for (j = -OFAST_UMAX[absI]; j < (OFAST_UMAX[absI] + 1); ++j) {
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

        extractDescriptors();

        combineAndFilterTotal();

        int nKeyPointsTotal = countKeypoints();

        //log.info("nKeypointsTotal=" + nKeyPointsTotal +
        //    " number of allowed Keypoints=" + this.nKeypoints);
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
     * @param img in col major format
     * @return float array of pyramid images in row major format
     */
    private List<TwoDFloatArray> buildPyramid(GreyscaleImage img) {

        //TODO: need a finer grained decimation as an option
        //   can build one with gaussian and down-sampling
        // of use the ATrous b3 method and down-sampling
        //   see the upsampling code in image processing to reverse...

        ImageProcessor imageProcessor = new ImageProcessor();

        int dl = useSingleScale ?
            Math.min(img.getWidth(), img.getHeight()) - 1 :
            decimationLimit;

        List<GreyscaleImage> output = imageProcessor.buildPyramid2(img, dl);

        List<TwoDFloatArray> output2 = new ArrayList<TwoDFloatArray>();
        float[][] rowMajorImg;
        TwoDFloatArray f;

        for (int i = 0; i < output.size(); ++i) {
            rowMajorImg = imageProcessor.copyToRowMajor(output.get(i));
            MatrixUtil.multiply(rowMajorImg, 1.f/255.f);
            f = new TwoDFloatArray(rowMajorImg);
            output2.add(f);
        }

        return output2;
    }

    protected void buildTensors() {

        // Standard deviation used for the Gaussian kernel, which is used as
        // weighting function for the auto-correlation matrix.
        float sigma = SIGMA.getValue(SIGMA.ZEROPOINTFIVE);

        tensorComponents = new StructureTensorR[scales.length];

        TwoDFloatArray octaveImage;

        for (int i = 0; i < scales.length; ++i) {

            octaveImage = pyramidImages[i];

            // octaveImage.a is in row-major format
            tensorComponents[i] = new StructureTensorR(
                    MatrixUtil.convertToDouble(octaveImage.a),
                sigma, false);
        }

    }

    private void buildHarrisImages() {

        harrisResponseImages = new TwoDFloatArray[scales.length];

        TwoDFloatArray octaveImage;

        for (int i = 0; i < scales.length; ++i) {

            octaveImage = pyramidImages[i];

            StructureTensorR tensors = tensorComponents[i];
            double[][] detA = tensors.getDeterminant();
            double[][] traceA = tensors.getTrace();

            double[][] harrisResponse = cornerHarris(
                    MatrixUtil.convertToDouble(octaveImage.a), detA, traceA);

            // size is same a octaveImage
            //float[][] harrisResponse = cornerHarris(octaveImage.a, detA, traceA);

            harrisResponseImages[i] = new TwoDFloatArray(
                    MatrixUtil.convertToFloat(harrisResponse));
        }
    }

    protected void debugPrint(String label, float[][] a) {

        String str = MiscDebug.getPrintRowMajor(a, label);

        log.info(str);
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

        Resp r;
        TIntList mask;
        float[][] harrisResponse;
        int n, count, i, c0, c1;
        float[] scores;
        PairInt[] points;
        TIntList kpc0s;
        TIntList kpc1s;
        float scale;
        int row, col;

        // largest pyramid harris response image and tensor image size
        int rH = harrisResponseImages[0].a.length;
        int cH = harrisResponseImages[0].a[0].length;

        Set<PairInt> set;
        TIntList kp0s;
        TIntList kp1s;
        float[][] octaveImage;
        TDoubleList orientations;
        TFloatList responses;
        int rmi = this.nKeypoints;
        int len;

        for (int octave = 0; octave < scales.length; ++octave) {
            // order by response and filter

            scale = scales[octave];

            // keypoints lists are in row major format:
            //  row numbers are in keypoints0, column numbers are in keypoints1
            r = extractHarrisCornerKeypoints(octave);

            //mask of length orientations.size() containing a 1 or 0
            // indicating if pixels are within the image (``True``) or in the
            // border region of the image (``False``).
            mask = filterNearBorder(pyramidImages[octave].a, r.keypoints0, r.keypoints1);

            harrisResponse = harrisResponseImages[octave].a;

            n = r.keypoints0.size();

            assert(harrisResponse.length == pyramidImages[octave].a.length);
            assert(harrisResponse.length == tensorComponents[octave].getDXSquared().length);
            assert(harrisResponse[0].length == pyramidImages[octave].a[0].length);
            assert(harrisResponse[0].length == tensorComponents[octave].getDXSquared()[0].length);

            scores = new float[n];
            points = new PairInt[n];
            count = 0;

            for (i = 0; i < n; ++i) {
                // these are in octave reference frame
                col = r.keypoints1.get(i); // col numbers
                row = r.keypoints0.get(i); // row numbers
                points[count] = new PairInt(row, col);
                scores[count] = harrisResponse[row][col];
                count++;
            }

            if (n > 0) {
                // the harris points are negative values.
                //    strongest response is the smallest value (i.e. -2 is stronger than 0)
                QuickSort.sortBy1stArg(scores, points);
            }

            kpc0s = new TIntArrayList(this.nKeypoints);
            kpc1s = new TIntArrayList(this.nKeypoints);

            for (i = 0; i < points.length; i++) {
                PairInt p = points[i];
                row = p.getX();
                col = p.getY();
                kpc1s.add(col);
                kpc0s.add(row);
            }

            // then filter
            maskCoordinates(kpc0s, kpc1s, harrisResponse.length,
                harrisResponse[0].length, 4);//8);

            // filter to remove redundant points when rescaled
            set = new HashSet<PairInt>();
            kp0s = new TIntArrayList(set.size());
            kp1s = new TIntArrayList(set.size());
            for (i = 0; i < kpc0s.size(); ++i) {
                // coords are in octave reference frame so expand to largest pyramid frame
                //    for redundancy and boundary check
                row = kpc0s.get(i);
                col = kpc1s.get(i);
                c0 = Math.round((float)row*scale);
                c1 = Math.round((float)col*scale);
                // if larger than largest pyramid image height
                if (c0 >= rH) {
                    c0 = rH - 1;
                }
                // if larger than largest pyramid image width
                if (c1 >= cH) {
                    c1 = cH - 1;
                }
                PairInt p = new PairInt(c0, c1);
                if (set.contains(p)) {
                    continue;
                }
                set.add(p);
                kp0s.add(row); // these are in ref frame of octave image
                kp1s.add(col);
            }
            // let the VM reclaim the memory if it chooses:
            set = null;
            kpc0s = null;
            kpc1s = null;

            // filter for number of points if requested
            if (kp0s.size() > this.nKeypoints) {
                len = kp0s.size() - this.nKeypoints;
                kp0s.remove(rmi, len);
                kp1s.remove(rmi, len);
            }
            // coords are in octave reference frame so expand to largest pyramid frame
            octaveImage = pyramidImages[octave].a;
            orientations = calculateOrientations(octaveImage, kp0s, kp1s);

            responses = new TFloatArrayList(kp0s.size());
            for (i = 0; i < kp0s.size(); ++i) {
                c0 = kp0s.get(i);
                c1 = kp1s.get(i);
                responses.add(harrisResponse[c0][c1]);
            }

            // transform keypoints to full size coordinate reference frame.
            for (i = 0; i < kp0s.size(); ++i) {
                c0 = Math.round(scale * kp0s.get(i));
                if (c0 >= rH) {
                    c0 = rH - 1;
                }
                kp0s.set(i, c0);

                c1 = Math.round(scale * kp1s.get(i));
                if (c1 >= cH) {
                    c1 = cH - 1;
                }
                kp1s.set(i, c1);
            }

            // kp0s and kp1s are in coord frame of largest image

            keypoints0List.add(kp0s);
            keypoints1List.add(kp1s);
            orientationsList.add(orientations);
            harrisResponses.add(responses);

            assert(orientations.size() == kp0s.size());
            assert(responses.size() == kp0s.size());
            assert(kp1s.size() == kp0s.size());

        } // end loop over octave

        // note, at this point, each octave has the limit of this.nKeypoints,
        //    but the total may be as large as scales.length*this.nKeypoints

        assert(pyramidImages.length == keypoints0List.size());
        assert(pyramidImages.length == keypoints1List.size());
        assert(pyramidImages.length == orientationsList.size());
        assert(pyramidImages.length == this.harrisResponses.size());
        assert(pyramidImages.length == scales.length);
    }

    private void combineAndFilterTotal() {
        //combined into one list and sort by harris response then take first nKeypoints

        combinedKeypointList = new ArrayList<PairInt>();
        combinedOrientationList = new TDoubleArrayList();
        combinedHarrisResponseList = new TFloatArrayList();

        if (descrChoice.ordinal() != DescriptorChoice.NONE.ordinal()) {
            combinedDescriptors = new Descriptors();
        }

        final int nTotal = countKeypoints();

        int octave, i;
        int count = 0;
        TFloatList h;

        if (nTotal < this.nKeypoints) {
            //fill the arrays and return

            VeryLongBitString[] combinedD = new VeryLongBitString[nTotal];
            TDoubleList orList;
            TIntList kp0List, kp1List;
            Descriptors d;

            for (octave = 0; octave < scales.length; ++octave) {
                kp0List = keypoints0List.get(octave);
                kp1List = keypoints1List.get(octave);
                orList = orientationsList.get(octave);
                h = harrisResponses.get(octave);
                for (i = 0; i < kp0List.size(); ++i) {
                    combinedKeypointList.add(new PairInt(kp0List.get(i), kp1List.get(i)));
                    combinedHarrisResponseList.add(h.get(i));
                    combinedOrientationList.add(orList.get(i));
                }
            }
            count = 0;
            if (descrChoice.ordinal() != DescriptorChoice.NONE.ordinal()) {
                for (octave = 0; octave < scales.length; ++octave) {
                    d = descriptorsList.get(octave);
                    System.arraycopy(d.descriptors, 0, combinedD, count, d.descriptors.length);
                    count += d.descriptors.length;
                }
                combinedDescriptors.descriptors = combinedD;
                assert (combinedDescriptors.descriptors.length == nTotal);
            }
            assert(combinedKeypointList.size() == nTotal);
            assert (combinedOrientationList.size() == nTotal);
            assert (combinedHarrisResponseList.size() == nTotal);

            return;
        }

        PairInt[] idxs00 = new PairInt[nTotal];
        int[] idxs0 = new int[nTotal];

        count = 0;
        for (octave = 0; octave < scales.length; ++octave) {
            h = harrisResponses.get(octave);
            for (i = 0; i < h.size(); ++i) {
                idxs00[count] = new PairInt(octave, i);
                idxs0[count] = count;
                count++;
                combinedHarrisResponseList.add(h.get(i));
            }
        }

        QuickSort.sortBy1stArg(combinedHarrisResponseList, idxs0);

        VeryLongBitString[] combinedD = new VeryLongBitString[this.nKeypoints];
        int idxOctave, idxList, idx0;
        PairInt p;

        for (i = 0; i < this.nKeypoints; ++i) {
            idx0 = idxs0[i];
            idxOctave = idxs00[idx0].getX();
            idxList = idxs00[idx0].getY();

            p = new PairInt(keypoints0List.get(idxOctave).get(idxList),
                keypoints1List.get(idxOctave).get(idxList));
            combinedKeypointList.add(p);

            combinedOrientationList.add(this.orientationsList.get(idxOctave).get(idxList));
        }
        if (descrChoice.ordinal() != DescriptorChoice.NONE.ordinal()) {
            count = 0;
            for (i = 0; i < this.nKeypoints; ++i) {
                idx0 = idxs0[i];
                idxOctave = idxs00[idx0].getX();
                idxList = idxs00[idx0].getY();
                // this is the descriptor for one keypoint
                combinedD[count] = descriptorsList.get(idxOctave).descriptors[idxList];
                count++;
            }
            combinedDescriptors.descriptors = combinedD;
            assert(combinedDescriptors.descriptors.length == this.nKeypoints);
        }

        combinedHarrisResponseList.remove(this.nKeypoints, nTotal - this.nKeypoints);

        assert(combinedKeypointList.size() == this.nKeypoints);
        assert(combinedOrientationList.size() == this.nKeypoints);
        assert(combinedHarrisResponseList.size() == this.nKeypoints);
    }

    private void extractDescriptors() {

        if (descrChoice.equals(DescriptorChoice.NONE)) {
            return;
        }

        int nScales = scales.length;

        descriptorsList = new ArrayList<Descriptors>(nScales);

        TIntList kp0;
        TIntList kp1;
        TDoubleList or;
        for (int i = 0; i < nScales; ++i) {

            float[][] octaveImage = pyramidImages[i].a;

            kp0 = this.keypoints0List.get(i); // row number
            kp1 = this.keypoints1List.get(i); // col number
            or = this.orientationsList.get(i);

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

    private Resp extractHarrisCornerKeypoints(int octave) {

        float[][] octaveImage = pyramidImages[octave].a;

        float scale = scales[octave];

        float sigmaBlur = 1;
        //StructureTensorR tensors = tensorComponents[octave];
        //double[][] detA = tensors.getDeterminant();
        //double[][] traceA = tensors.getTrace();

        TIntList keypoints0 = new TIntArrayList();
        TIntList keypoints1 = new TIntArrayList();

        ImageProcessor imageProcessor = new ImageProcessor();

        // size is same a octaveImage
        float[][] harrisResponse = harrisResponseImages[octave].a;

        imageProcessor.peakLocalMax(harrisResponse,
                1, 0.01f,
                true, keypoints0, keypoints1);

        Resp r2 = new Resp();
        r2.keypoints0 = keypoints0;
        r2.keypoints1 = keypoints1;

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

        for (int i = 0; i < keypoints0.size(); ++i) {
            int ii = keypoints0.get(i); // row number
            int jj = keypoints1.get(i); // column number
            int pixIdx = (ii * nCols) + jj;
            peaks.add(pixIdx);
        }

        for (int i = 0; i < keypoints0.size(); ++i) {
            int ii = keypoints0.get(i); // row number
            int jj = keypoints1.get(i); // column number
            int pixIdx = (ii * nCols) + jj;
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
                    int pixIdx2 = (k0 * nCols) + k1;
                    if (peaks.contains(pixIdx2)) {
                        peaks.remove(pixIdx2);
                    }
                }
            }
        }

        TIntList keypoints0_2 = new TIntArrayList(peaks.size());
        TIntList keypoints1_2 = new TIntArrayList(peaks.size());
        for (int i = 0; i < keypoints0.size(); ++i) {
            int ii = keypoints0.get(i); // row number
            int jj = keypoints1.get(i); // column number
            int pixIdx = (ii * nCols) + jj;
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
        TODO: have changed the use of keypoints0 and keypoints1 to hold the
        row numbers and column numbers, respectively:
            pixIdxs[count] = (rowNumber * nCols) + colNumber;
        */

        for (int i = 0; i < coords0.size(); ++i) {
            int ii = coords0.get(i); // row number
            int jj = coords1.get(i); // column number
            int pixIdx = (ii * nCols) + jj;
            set.add(pixIdx);
        }

        int nBefore = set.size();

        for (int i = 0; i < coords0.size(); ++i) {
            int ii = coords0.get(i); // row number
            int jj = coords1.get(i); // column number
            int pixIdx = (ii * nCols) + jj;
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
            int ii = coords0.get(i); // row number
            int jj = coords1.get(i); // column number
            int pixIdx = (ii * nCols) + jj;
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
     * https://github.com/scikit-image/scikit-image/blob/d19b60add22b818298c7aefa65fharris0e7c1467ef4d/skimage/feature/corner.py
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
        boolean indices = true;
        int numPeaks = Integer.MAX_VALUE;
        //footprint=None, labels=None

        int nRows = img.length;
        int nCols = img[0].length;
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        // these results have been sorted by decreasing intensity
        imageProcessor.peakLocalMax(img, minDistance, thresholdRel,
            true, outputKeypoints0, outputKeypoints1);

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
     * @param keypoints0 keypoint y coords in reference frame of this octave image
     * @param keypoints1 keypoint x coords in reference frame of this octave iamge
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
                for (int ii = 0; ii < src.length; ++ii) {
                    cImage[i][ii + nMaskCols2] = src[ii];
                }
            }
        }

        // number of corner coord pairs
        int nCorners = keypoints0.size();

        TDoubleList orientations = new TDoubleArrayList(nCorners);

        // eqns (1), (2), (3) of Rublee et al.
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
        This corner detector uses information from the auto-correlation matrix A
        (a.k.a. he structure tensor):
            A = [(imx**2)   (imx*imy)] = [Axx Axy]
                [(imx*imy)   (imy**2)]   [Axy Ayy]
        Where imx and imy are first derivatives, averaged with a gaussian filter.
        The corner measure is then defined as::
            det(A) - k * trace(A)**2

      adapted from
      from https://github.com/scikit-image/scikit-image/blob/master/skimage/feature/corner.py

     NOTE: made the change for case where entire detA is 0 to use +k*Trace^2;

    @param image
     * @param detA
     * @param traceA
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
        // unless detA is all 0's, then will use only +k * traceA**2
        float[][] response = new float[detA.length][];
        boolean detAAll0s = true;
        for (int i = 0; i < detA.length; ++i) {
            response[i] = new float[detA[0].length];
            for (int j = 0; j < detA[i].length; ++j) {
                if (detA[i][j] != 0.f) {
                    detAAll0s = false;
                    break;
                }
            }
        }
        
        for (int i = 0; i < detA.length; ++i) {
            for (int j = 0; j < detA[i].length; ++j) {
                float v = k * (traceA[i][j] * traceA[i][j]);
                if (detAAll0s) {
                    response[i][j] = v;
                } else {
                    response[i][j] = detA[i][j] - v;
                }
            }
        }

        return response;
    }

    /**
     Compute Harris corner measure response image
     ([Harris and Stephens, 1988).
     This corner detector uses information from the auto-correlation matrix A
     (a.k.a. the structure tensor):
     <pre>
     A = [(imx**2)   (imx*imy)] = [Axx Axy]
        [(imx*imy)   (imy**2)]   [Axy Ayy]
     </pre>
     Where imx and imy are first derivatives, averaged with a gaussian filter.
     The corner measure is then defined as::
     MASKS eqn (11.1)
         det(A) + k * trace(A)**2

     originally adapted
     from https://github.com/scikit-image/scikit-image/blob/master/skimage/feature/corner.py

     @param image
     @return the harris corner response image of same size as image,
      *   and composed of
      *   response = detA + k * traceA ** 2 built from the 2nd derivatives of
      *   image intensity.
     */
    public static double[][] cornerHarris(double[][] image) {

        float sigma = 1.f;

        StructureTensorR tensor = new StructureTensorR(image, sigma, false);
        double[][] detA = tensor.getDeterminant();
        double[][] traceA = tensor.getTrace();
        return cornerHarris(image, detA, traceA);
    }

    public static double[][] cornerHarris(double[][] image, double[][] detA, double[][] traceA) {

        float sigma = 1.f;

        StructureTensorR tensor = new StructureTensorR(image, sigma, false);

        // method = 'k'.  k is Sensitivity factor to separate corners from edges,
        // Small values of k result in detection of sharp corners.
        double k = 0.04;

        int i;
        int j;
        boolean detAAll0s = true;
        //response = detA - k * traceA ** 2
        double[][] response = new double[detA.length][];
        for (i = 0; i < detA.length; ++i) {
            response[i] = new double[detA[0].length];
            for (j = 0; j < detA[i].length; ++j) {
                if (detA[i][j] != 0.f) {
                    detAAll0s = false;
                    break;
                }
            }
        }
        double v;
        for (i = 0; i < detA.length; ++i) {
            for (j = 0; j < detA[i].length; ++j) {
                v = k * (traceA[i][j] * traceA[i][j]);
                if (detAAll0s) {
                    response[i][j] = v;
                } else {
                    response[i][j] = detA[i][j] - v;
                }
            }
        }

        return response;
    }

    /**
     * filter the keypoints, orientations, and responses to remove those close
     * to the border and then create descriptors for the remaining.
     *
     * @param octaveImage
     * @param keypoints0 holds row numbers
     * @param keypoints1 holds column numbers
     * @return the encapsulated descriptors and mask
     */
    protected TIntList filterNearBorder(float[][] octaveImage,
        TIntList keypoints0, TIntList keypoints1) {

        int nRows = octaveImage.length;
        int nCols = octaveImage[0].length;

        int border = 5;
        if (Math.min(nRows, nCols)/20. < border) {
            border = (int)(Math.min(nRows, nCols)/20.);
        }

        TIntList mask = maskCoordinatesIfBorder(keypoints0, keypoints1,
            nRows, nCols, border);//16);

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
     * @param keypoints0 keypoint y coordinates in the 
     * reference frame of the largest pyramid image
     * @param keypoints1 keypoint x coordinates in the 
     * reference frame of the largest pyramid image
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
     * @param keypoints0 keypoint row numbers in reference frame of the largest
     * pyramid image
     * @param keypoints1 keypoint column numbers in reference frame of the largest
     * pyramid image
     * @param orientations
     * @param scale the scale for this octave image.
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
        
        int nKP = keypoints0.size();

        log.fine("nKP=" + nKP);

        // each item is a 256 length bit string
        VeryLongBitString[] descriptors = new VeryLongBitString[nKP];
        
        if (nKP == 0) {
            return descriptors;
        }
        
        //System.out.printf("scale=%.2f, nKP=%d\n", scale, nKP);

        double pr0, pc0, pr1, pc1, angle, sinA, cosA;
        int spr0, spc0, spr1, spc1, kr, kc, i, j, row0, col0, row1, col1;
      
        for (i = 0; i < nKP; ++i) {

            descriptors[i] = new VeryLongBitString(POS1.length);

            angle = orientations.get(i);
            sinA = Math.sin(angle);
            cosA = Math.cos(angle);

            kr = keypoints0.get(i); // row number
            kc = keypoints1.get(i); // col number

            // put kr and kc into this octave's pyramid image reference frame.
            kr = (int)(kr/scale);
            kc = (int)(kc/scale);

            for (j = 0; j < POS0.length; ++j) {
                pr0 = POS0[j][0];
                pc0 = POS0[j][1];
                pr1 = POS1[j][0];
                pc1 = POS1[j][1];

                spr0 = (int)Math.round(sinA * pr0 + cosA * pc0);
                spc0 = (int)Math.round(cosA * pr0 - sinA * pc0);
                spr1 = (int)Math.round(sinA * pr1 + cosA * pc1);
                spc1 = (int)Math.round(cosA * pr1 - sinA * pc1);

                row0 = kr + spr0;
                col0 = kc + spc0;
                row1 = kr + spr1;
                col1 = kc + spc1;

                if (row0 < 0 || col0 < 0 || row1 < 0 || col1 < 0 ||
                    (row0 > (octaveImage.length - 1)) ||
                    (row1 > (octaveImage.length - 1)) ||
                    (col0 > (octaveImage[0].length - 1)) ||
                    (col1 > (octaveImage[0].length - 1))
                    ) {
                    continue;
                }
                // eqn (4) of Rublee et al.
                if (octaveImage[row0][col0] < octaveImage[row1][col1]) {
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
        return combinedDescriptors;
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
     * combine key-point lists k1 and k2, removing redundancies, then
     * create a combined Descriptors for the combined key-point's descriptors,
     * and a combined orientations list for the same combined key-points.
     * @param d1 descriptors from dataset 1
     * @param k1 key-points from dataset 1
     * @param o1 orientations from dataset 1
     * @param d2 descriptors from dataset 2
     * @param k2 key-points from dataset 2
     * @param o2 orientations from dataset 2
     * 
     * @return a 3-dimensional array holding the combined Descriptors in the
     * first element, the combined key-points in the second element, and the
     * combined orientations in the third element.
     * i.e. new Object[]{combinedDescriptors, combinedKeyPoints} with
     */
    public static Object[] combine(Descriptors d1, List<PairInt> k1, TDoubleList o1,
        Descriptors d2, List<PairInt> k2, TDoubleList o2) {
        
        if (d1.descriptors.length != k1.size()) {
            throw new IllegalArgumentException("d1.descriptors.length and k1.size must be the same");
        }
        if (d2.descriptors.length != k2.size()) {
            throw new IllegalArgumentException("d2.descriptors.length and k2.size must be the same");
        }
        if (k1.size() != o1.size()) {
            throw new IllegalArgumentException("k1.size and o1.size must be the same");
        }
        if (k2.size() != o2.size()) {
            throw new IllegalArgumentException("k2.size and o2.size must be the same");
        }
        
        int n = k1.size() + k2.size();

        Set<PairInt> k12 = new HashSet<PairInt>();
                
        List<PairInt> kCombined = new ArrayList<PairInt>();
        Descriptors dCombined = new Descriptors();
        dCombined.descriptors = new VeryLongBitString[n];
        TDoubleList oCombined = new TDoubleArrayList();
        
        PairInt p;
        int idx, j;
        for (idx = 0; idx < k1.size(); ++idx) {
            p = k1.get(idx);
            if (!k12.contains(p)) {
                kCombined.add(k1.get(idx));
                oCombined.add(o1.get(idx));
                j = k12.size();
                dCombined.descriptors[j] = d1.descriptors[idx].copy();
                k12.add(k1.get(idx));
            }
        }    
        for (idx = 0; idx < k2.size(); ++idx) {
            p = k2.get(idx);
            if (!k12.contains(p)) {
                kCombined.add(p);
                oCombined.add(o2.get(idx));
                j = k12.size();
                dCombined.descriptors[j] = d2.descriptors[idx].copy();
                k12.add(p);
            }
        }
        
        if (k12.size() < n) {
            dCombined.descriptors = Arrays.copyOf(dCombined.descriptors, k12.size());
        }
        
        assert(dCombined.descriptors.length == kCombined.size());
        assert(kCombined.size() == oCombined.size());
        assert(kCombined.size() == k12.size());
        
        return new Object[]{dCombined, kCombined, oCombined};
    }

    /**
     * get a list of each octave's keypoint rows as a combined list.
     * The list contains coordinates which have already been scaled to the
     * full image reference frame.  These are the row coordinates of
     * key-points.  (The column coordinates are in keypoints1).
     * @return
     */
    public TIntList getAllKeyPoints0() {

        TIntList combined = new TIntArrayList(combinedKeypointList.size());
        PairInt p;
        for (int i = 0; i < combinedKeypointList.size(); ++i) {
            p = combinedKeypointList.get(i);
            combined.add(p.getX());
        }

        return combined;
    }
    
    /**
     * get a list of each octave's keypoint cols as a combined list.
     * The list contains coordinates which have already been scaled to the
     * full image reference frame.
     * These are the column coordinates of
     * key-points.  (The row coordinates are in keypoints0).
     * @return
     */
    public TIntList getAllKeyPoints1() {

        TIntList combined = new TIntArrayList(combinedKeypointList.size());
        PairInt p;
        for (int i = 0; i < combinedKeypointList.size(); ++i) {
            p = combinedKeypointList.get(i);
            combined.add(p.getY());
        }

        return combined;
    }


    /**
     * get a list of each octave's keypoint rows as a combined list.
     * The list contains coordinates which have already been scaled to the
     * full image reference frame and are in row major format (= row, col).
     * @return
     */
    public List<PairInt> getAllKeyPointsRC() {
        return combinedKeypointList;
    }


    /**
     * get a list of each octave's orientations as a combined list.
     * The corresponding coordinates can be obtained with getAllKeyPoints();
     * @return
     */
    public TDoubleList getAllOrientations() {
        return combinedOrientationList;
    }

    /**
     * get a list of each octave's harris response strengths as a combined list.
     * The corresponding coordinates can be obtained with getAllKeyPoints();
     * @return
     */
    public TFloatList getAllHarrisResponses() {
        return combinedHarrisResponseList;
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
