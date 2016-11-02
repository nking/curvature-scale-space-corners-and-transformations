package algorithms.imageProcessing.features;

import algorithms.QuickSort;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.StructureTensor;
import algorithms.imageProcessing.transform.MatchedPointsTransformationCalculator;
import algorithms.imageProcessing.transform.TransformationParameters;
import algorithms.imageProcessing.transform.Transformer;
import algorithms.misc.MiscMath;
import algorithms.search.NearestNeighbor2D;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.TwoDFloatArray;
import algorithms.util.VeryLongBitString;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.TFloatList;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.awt.Color;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

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

    /*
    TODO:
       -- considering adding alternative pyramid building methods
          such as Laplacian pyramids or
          half-octave or quarter-octave pyramids (Lowe 2004; Triggs 2004)
    */

final static PairInt tp2 = new PairInt(38, 72);
final static PairInt tp1 = new PairInt(33, 61);

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

    protected List<TIntList> keypoints0List = null;
    protected List<TIntList> keypoints1List = null;
    protected List<TDoubleList> orientationsList = null;
    protected List<TFloatList> harrisResponses = null;
    protected List<TFloatList> scalesList = null;

    //`True`` or ``False`` representing
    //the outcome of the intensity comparison for i-th keypoint on j-th
    //decision pixel-pair. It is ``Q == np.sum(mask)``.
    private List<Descriptors> descriptorsList = null;
    
    private List<Descriptors> descriptorsListH = null;
    private List<Descriptors> descriptorsListS = null;
    private List<Descriptors> descriptorsListV = null;
    
    protected static double twoPI = 2. * Math.PI;

    protected float curvatureThresh = 0.05f;

    // pyramid images will be no smaller than this
    protected final int decimationLimit = 32;
    
    /**
     * @return the pyramidImages
     */
    public List<TwoDFloatArray> getPyramidImages() {
        return pyramidImages;
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
        NONE, GREYSCALE, HSV
    }
    protected DescriptorChoice descrChoice = DescriptorChoice.GREYSCALE;
    
    public static enum DescriptorDithers {
        NONE, FIFTEEN, TWENTY, FORTY_FIVE, NINETY, ONE_HUNDRED_EIGHTY;
    }
    protected DescriptorDithers descrDithers = DescriptorDithers.NONE;
    
    protected boolean doCreate1stDerivKeypoints = false;
    
    protected boolean doCreateCurvatureKeyPoints = false;

    protected final int nPyramidB = 3;
    
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

        // extract keypoints and store results and images at class level
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
        
        if (nKeyPointsTotal > this.nKeypoints) {
            
            // prune the lists to nKeypoints by the harris corner response as a score

            // once thru to count first
            int n = 0;
            for (int idx1 = 0; idx1 < harrisResponses.size(); ++idx1) {
                TFloatList hResp = harrisResponses.get(idx1);
                n += hResp.size();
                assert(keypoints0List.get(idx1).size() == hResp.size());
                assert(keypoints1List.get(idx1).size() == hResp.size());
                assert(orientationsList.get(idx1).size() == hResp.size());
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
 {
float scale = scalesList.get(idx1).get(idx2);
int t1XS = Math.round((float)tp1.getX()/scale);
int t1YS = Math.round((float)tp1.getY()/scale);
int x = keypoints1List.get(idx1).get(idx2);
int y = keypoints0List.get(idx1).get(idx2);
if (Math.abs(t1XS - x) < 3 && Math.abs(t1YS - y) < 3) {
    System.out.println("for scale=" + scale + " resp=" + hResp.get(idx2));
}
 }                   
                }
            }

            QuickSort.sortBy1stArg(scores, indexes);

            System.out.println("--> strongest response=" + 
                scores[scores.length - 1] 
                + " and response at last "
                + "requested keypoint=" 
                + scores[scores.length - nKeypoints - 1]);
            
            List<TIntList> keypoints0List2 = new ArrayList<TIntList>();
            List<TIntList> keypoints1List2 = new ArrayList<TIntList>();
            List<TDoubleList> orientationsList2 = new ArrayList<TDoubleList>();
            List<TFloatList> harrisResponses2 = new ArrayList<TFloatList>();
            List<TFloatList> scalesList2 = new ArrayList<TFloatList>();
            for (int ii = 0; ii < keypoints0List.size(); ++ii) {
                keypoints0List2.add(new TIntArrayList());
                keypoints1List2.add(new TIntArrayList());
                orientationsList2.add(new TDoubleArrayList());
                harrisResponses2.add(new TFloatArrayList());
                scalesList2.add(new TFloatArrayList());
            }
            
            // visit largest scores to smallest
            count = 0;
            int idx = n - 1;
            while (count < nKeypoints) {

                PairInt m = indexes[idx];
                int idx1 = m.getX();
                int idx2 = m.getY();

                keypoints0List2.get(idx1).add
                    (keypoints0List.get(idx1).get(idx2));
                keypoints1List2.get(idx1).add
                    (keypoints1List.get(idx1).get(idx2));
                orientationsList2.get(idx1).add
                    (orientationsList.get(idx1).get(idx2));
                harrisResponses2.get(idx1).add
                    (harrisResponses.get(idx1).get(idx2));
                scalesList2.get(idx1).add
                    (scalesList.get(idx1).get(idx2));
                
                idx--;
                count++;
            }

            keypoints0List = keypoints0List2;
            keypoints1List = keypoints1List2;
            orientationsList = orientationsList2;
            harrisResponses = harrisResponses2;
            scalesList = scalesList2;
        
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

        extractDescriptors();
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
     
        ImageProcessor imageProcessor = new ImageProcessor();
        
        List<GreyscaleImage> outputR = imageProcessor.buildPyramid2(
            imageR, decimationLimit, nPyramidB);
        
        List<GreyscaleImage> outputG = imageProcessor.buildPyramid2(
            imageG, decimationLimit, nPyramidB);
        
        List<GreyscaleImage> outputB = imageProcessor.buildPyramid2(
            imageB, decimationLimit, nPyramidB);
        
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
     * @param image
     * @return
     */
    private List<TwoDFloatArray> buildPyramid(GreyscaleImage img) {

        //TODO: need a finer grained decimation as an option
        //   can build one with gaussian and down-sampling
        // of use the ATrous b3 method and down-sampling
        //   see the upsampling code in image processing to reverse...
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        List<GreyscaleImage> output = imageProcessor.buildPyramid2(
            img, decimationLimit, nPyramidB);

        List<TwoDFloatArray> output2 = new ArrayList<TwoDFloatArray>();
        for (int i = 0; i < output.size(); ++i) {
            float[][] gsImgF = imageProcessor.multiply(output.get(i), 1.f/255.f);
            TwoDFloatArray f = new TwoDFloatArray(gsImgF);
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
        }
        
        float prevScl = 1;
        
        TIntList octavesUsed = new TIntArrayList();

        for (int octave = 0; octave < pyramid.size(); ++octave) {

            //System.out.println("octave=" + octave);

            float[][] octaveImage = pyramid.get(octave).a;

            //multiply the keypoints by a scalar
            //   to put the coordinates into the reference frame of image
            //float scale = (float)Math.pow(this.downscale, octave);
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

            if (repeatScale1At2 && (octave == 1)) {
                //NOTE: an experimental look at making sure strong responses
                // in scale=1 image are also sampled at octave=1.
                // If this corrects a problem, then the selection of keypoints
                // from fast response needs to be improved for lower resolution
                // images.
                
                int n0 = keypoints0List.get(0).size();
                Set<PairInt> points0 = new HashSet<PairInt>(n0);
                PairInt[] pointsA0 = new PairInt[n0];
                float[] r0 = new float[n0];
                for (int i = 0; i < n0; ++i) {
                    int row = keypoints0List.get(0).get(i);
                    int col = keypoints1List.get(0).get(i);
                    PairInt p = new PairInt(row, col);
                    points0.add(p);
                    pointsA0[i] = p;
                    r0[i] = this.harrisResponses.get(0).get(i);
                }
                QuickSort.sortBy1stArg(r0, pointsA0);
                System.out.println("strongest response=" +
                    r0[r0.length - 1] + " median response=" +
                    r0[r0.length/2]
                    + " and fast response thresh=" + this.fastThreshold
                );
                int t1XS = Math.round((float)tp1.getX()/scale);
                int t1YS = Math.round((float)tp1.getY()/scale);
                for (int i = 0; i < n0; ++i) {
                    int x = pointsA0[i].getY();
                    int y = pointsA0[i].getX();
                    if (Math.abs(x - t1XS) < 3 && Math.abs(y - t1YS) < 3) {
                        System.out.println(" oct 0 tie response=" + r0[i]);
                    }
                }
                int z = 1;
            }
            
            //mask of length orientations.size() containing a 1 or 0
            // indicating if pixels are within the image (``True``) or in the
            // border region of the image (``False``).
            TIntList mask = filterNearBorder(octaveImage, r.keypoints0,
                r.keypoints1);
            
            TFloatList responses = new TFloatArrayList(r.keypoints0.size());
            for (int i = 0; i < r.keypoints0.size(); ++i) {
                int x = r.keypoints0.get(i);
                int y = r.keypoints1.get(i);
                float v = r.harrisResponse[x][y];
                responses.add(v);
            }
        
            TDoubleList orientations = cornerOrientations(octaveImage,
                r.keypoints0, r.keypoints1);

            // transform keypoints to full size coordinate reference frame  
            for (int i = 0; i < r.keypoints0.size(); ++i) {
                int v = Math.round(scale * r.keypoints0.get(i));
                r.keypoints0.set(i, v);
                v = Math.round(scale * r.keypoints1.get(i));
                r.keypoints1.set(i, v);
            }

            octavesUsed.add(octave);

            keypoints0List.add(r.keypoints0);
            keypoints1List.add(r.keypoints1);
            orientationsList.add(orientations);
            harrisResponses.add(responses);

            pyramidImages.add(pyramid.get(octave));
            harrisResponseImages.add(new TwoDFloatArray(r.harrisResponse));

            TFloatList scales = new TFloatArrayList(r.keypoints0.size());
            for (int i = 0; i < r.keypoints0.size(); ++i) {
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

    private void extractDescriptors() {
        
        if (descrChoice.equals(DescriptorChoice.NONE)) {
            return;
        }

        if (!descrChoice.equals(DescriptorChoice.HSV)) {
            descriptorsList = extractGreyscaleDescriptors(pyramidImages);
            return;
        }
        
        descriptorsListH = extractGreyscaleDescriptors(pyramidImagesH);
        descriptorsListS = extractGreyscaleDescriptors(pyramidImagesS);
        descriptorsListV = extractGreyscaleDescriptors(pyramidImagesV);
    }

    protected List<Descriptors> extractGreyscaleDescriptors(List<TwoDFloatArray> 
        pyramid) {

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
            Descriptors desc = extractOctave(octaveImage, kp0, kp1, or);

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

int t2XS = Math.round((float)tp2.getX()/scale);
int t2YS = Math.round((float)tp2.getY()/scale);
System.out.println(
String.format("srch img tie (%d,%d) fr=%.4f", t2XS, t2YS, 
    fastResponse[t2YS][t2XS])
);
int t1XS = Math.round((float)tp1.getX()/scale);
int t1YS = Math.round((float)tp1.getY()/scale);
System.out.println(
String.format("template img tie (%d,%d) fr=%.4f (scale=%.1f)", t1XS, t1YS, 
    fastResponse[t1YS][t1XS], scale)
);

        // list of format [row, col, ...] of filtered maxima ordered by intensity
        cornerPeaks(fastResponse, 1, keypoints0, keypoints1);
        if (keypoints0.isEmpty()) {
            return null;
        }

        maskCoordinates(keypoints0, keypoints1, nRows, nCols, 8);//16);

if (keypoints0.contains(t2YS) && keypoints1.contains(t2XS)
    && (keypoints1.get( keypoints0.indexOf(t2YS)) == t2XS)) {
    System.out.println("found srch img tie in fast response corner");
}
if (keypoints0.contains(t1YS) && keypoints1.contains(t1XS)
    && (keypoints1.get( keypoints0.indexOf(t1YS)) == t1XS)) {
    System.out.println("found template img tie in fast response corner");
}
        
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
        
if (keypoints0.contains(t2YS) && keypoints1.contains(t2XS)
    && (keypoints1.get( keypoints0.indexOf(t2YS)) == t2XS)) {
    System.out.println("found srch img tie in in or before first deriv");
}
if (keypoints0.contains(t1YS) && keypoints1.contains(t1XS)
    && (keypoints1.get( keypoints0.indexOf(t1YS)) == t1XS)) {
    System.out.println("found template img tie in in or before first deriv");
}

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
 
if (keypoints0.contains(t2YS) && keypoints1.contains(t2XS)
    && (keypoints1.get( keypoints0.indexOf(t2YS)) == t2XS)) {
    System.out.println("found srch img tie in in or before curvature");
}
if (keypoints0.contains(t1YS) && keypoints1.contains(t1XS)
    && (keypoints1.get( keypoints0.indexOf(t1YS)) == t1XS)) {
    System.out.println("found srch img tie in in or before curvature");
}

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
                continue;
            }
            exists.add(p);
        }
        
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
        
        if (POS0 == null) {
            POS0 = ORBDescriptorPositions.POS0;
        }
        if (POS1 == null) {
            POS1 = ORBDescriptorPositions.POS1;
        }

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
        TDoubleList orientations) {
        
        if (POS0 == null) {
            POS0 = ORBDescriptorPositions.POS0;
        }
        if (POS1 == null) {
            POS1 = ORBDescriptorPositions.POS1;
        }

        assert(orientations.size() == keypoints0.size());
        assert(orientations.size() == keypoints1.size());

        VeryLongBitString[] descriptors = null;

        if (descrChoice.equals(DescriptorChoice.NONE)) {
            descriptors = new VeryLongBitString[0];
        } else {
            descriptors = orbLoop(octaveImage, keypoints0, keypoints1,
                orientations);
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
        TIntList keypoints1, TDoubleList orientations) {

        assert(orientations.size() == keypoints0.size());

        if (POS0 == null) {
            POS0 = ORBDescriptorPositions.POS0;
        }
        if (POS1 == null) {
            POS1 = ORBDescriptorPositions.POS1;
        }

        int nKP = orientations.size();
        
        System.out.println("nKP=" + nKP);

        // holds values 1 or 0.  size is [orientations.size][POS0.length]
        VeryLongBitString[] descriptors = new VeryLongBitString[nKP];

        double pr0, pc0, pr1, pc1;
        int spr0, spc0, spr1, spc1;

        for (int i = 0; i < descriptors.length; ++i) {
            
            descriptors[i] = new VeryLongBitString(256);
            
            double angle = orientations.get(i);
            double sinA = Math.sin(angle);
            double cosA = Math.cos(angle);

            int kr = keypoints0.get(i);
            int kc = keypoints1.get(i);

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
     * full image reference frame.
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
    
    /**
     * match descriptors using euclidean transformation evaluation from pairs in
     * feasible combinations of best matches.
     * Thie method is useful when a geometrical model is necessary to find the
     * object in an image where the background and foreground have changed and the
     * lighting may have changed, for example.  
     * Fpr known simpler conditions such as stereo-projections, RANSAC with
     * an epipolar or euclidean evaluator can be used instead with the top 45
     * results and that method will
     * be offered in this class one day.
     * @param d1
     * @param d2
     * @param keypoints1
     * @param keypoints2
     * @param scaleFactor the factor that will be the maximum dimension
     * of the object in the search image, for example, 2.
     * @return 
     */
    public static CorrespondenceList matchDescriptors2(
        Descriptors[] d1, Descriptors[] d2,
        List<PairInt> keypoints1,
        List<PairInt> keypoints2, float scaleFactor) {
        
        assert(d1.length == d2.length);
        
        int nBands = d1.length;
        
        int n1 = d1[0].descriptors.length;
        int n2 = d2[0].descriptors.length;
        
        assert(n1 == keypoints1.size());
        assert(n2 == keypoints2.size());

        //[n1][n2]
        int[][] cost = calcDescriptorCostMatrix(d1, d2);

        int nTot = n1 * n2;
        PairInt[] indexes = new PairInt[nTot];
        int[] costs = new int[nTot];
        int count = 0;
        for (int i = 0; i < n1; ++i) {
PairInt p1 = keypoints1.get(i);
boolean t1 = distance(p1, tp1) < 3;
if (t1)
System.out.println("found tie in template");
            for (int j = 0; j < n2; ++j) {
                indexes[count] = new PairInt(i, j);
                costs[count] = cost[i][j];
                count++;                
PairInt p2 = keypoints2.get(j);
boolean t2 = distance(p2, tp2) < 3;
if (t1 && t2) {
System.out.println("tie c=" + costs[count - 1]);
} else if (t2) {
    System.out.println("found tie in search image");
}
            }
        }
        assert(count == nTot);

        QuickSort.sortBy1stArg(costs, indexes);
        
        // storing the indexes of matches with cost < 127 and
        // fewer in number than about 45
        
        count = 0;
        Set<PairInt> set1 = new HashSet<PairInt>();
        Set<PairInt> set2 = new HashSet<PairInt>();
        PairIntArray mT = new PairIntArray();
        PairIntArray mS = new PairIntArray();
        TIntList tIndexes = new TIntArrayList();
        
        // below, need to look up cost using idx1 and p2: idx1, p2 -> cost
        TIntObjectMap<TObjectIntMap<PairInt>> idx1P2CostMap = 
            new TIntObjectHashMap<TObjectIntMap<PairInt>>();
        
        // for scale 75%, 0.25*3*256 = 192
        int diffLimit = Math.round(0.25f * nBands * 256);
        
        // visit lowest costs (== differences) first
        for (int i = 0; i < nTot; ++i) {
            if (costs[i] > diffLimit 
                //|| count > 45
                ) {
                break;
            }
            PairInt index12 = indexes[i];
            int idx1 = index12.getX();
            int idx2 = index12.getY();
            PairInt p1 = keypoints1.get(idx1);
            PairInt p2 = keypoints2.get(idx2);
            if (set1.contains(p1) //|| set2.contains(p2)
                ) {
                continue;
            }
        System.out.println("p1=" + p1 + " " + " p2=" + p2 + " cost=" + costs[i]);
            mT.add(p1.getX(), p1.getY());
            mS.add(p2.getX(), p2.getY());
            set1.add(p1);
            set2.add(p2);
            tIndexes.add(idx1);
            
            TObjectIntMap<PairInt> cMap = idx1P2CostMap.get(idx1);
            if (cMap == null) {
                cMap = new TObjectIntHashMap<PairInt>();
                idx1P2CostMap.put(idx1, cMap);
            }
            if (!cMap.containsKey(p2)) {
                // only store the smaller cost, that's reached first
                cMap.put(p2, costs[i]);
            }
            
            count++;
        }
        
        int[] minMaxXY = MiscMath.findMinMaxXY(keypoints1);
        int objDimension = Math.max(minMaxXY[1] - minMaxXY[0],
            minMaxXY[3] - minMaxXY[2]);
        int limit = Math.round(scaleFactor * objDimension);
        int limitSq = limit * limit;
       
        int[] minMaxXY2 = MiscMath.findMinMaxXY(set2);
        
        NearestNeighbor2D nn = new NearestNeighbor2D(
           set2, minMaxXY2[1] + limit, minMaxXY2[3] + limit);
       
        return completeUsingCombinations(
            keypoints1, keypoints2, mT, mS, nn,
            minMaxXY2, limit,
            tIndexes, idx1P2CostMap, indexes, costs
        );       
    }

    /**
     * match descriptors using euclidean transformation evaluation from pairs in
     * feasible combinations of best matches.
     * Thie method is useful when a geometrical model is necessary to find the
     * object in an image where the background and foreground have changed and the
     * lighting may have changed, for example.  
     * Fpr known simpler conditions such as stereo-projections, RANSAC with
     * an epipolar or euclidean evaluator can be used instead with the top 45
     * results and that method will
     * be offered in this class one day.
     * @param d1
     * @param d2
     * @param keypoints1
     * @param keypoints2
     * @param segmentedCells the segmented cells from image 2, that is,
     * keypoints2 are contained within this list of point sets.
     * @param sLimit maximum number of top keypoints to use per segmented cell 
     * @param scaleFactor the factor that will be the maximum dimension
     * of the object in the search image, for example, 2.
     * @return 
     */
    public static CorrespondenceList matchDescriptors2(
        Descriptors[] d1, Descriptors[] d2,
        List<PairInt> keypoints1, List<PairInt> keypoints2,
        List<Set<PairInt>> segmentedCells, int sLimit,
        float scaleFactor
        ) {
        
        assert(d1.length == d2.length);
        
        int nBands = d1.length;
        
        int n1 = d1[0].descriptors.length;
        int n2 = d2[0].descriptors.length;
        
        assert(n1 == keypoints1.size());
        assert(n2 == keypoints2.size());

        //[n1][n2]
        int[][] cost = calcDescriptorCostMatrix(d1, d2);

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
        
        // storing the indexes of matches with cost < 127 and
        // fewer in number than about 45
        // and limited to sLimit per segmented cell
        
        // key = segment point, value = list index
        TObjectIntMap<PairInt> pointIndexMap = createIndexMap(segmentedCells);
        
        // key = segmented cell index, value = number of points in comb subset
        TIntIntMap segIndexCount = new TIntIntHashMap();
        
        count = 0;
        Set<PairInt> set1 = new HashSet<PairInt>();
        Set<PairInt> set2 = new HashSet<PairInt>();
        PairIntArray mT = new PairIntArray();
        PairIntArray mS = new PairIntArray();
        TIntList tIndexes = new TIntArrayList();
        
        // below, need to look up cost using idx1 and p2: idx1, p2 -> cost
        TIntObjectMap<TObjectIntMap<PairInt>> idx1P2CostMap = 
            new TIntObjectHashMap<TObjectIntMap<PairInt>>();
        
        // for scale 75%, 0.25*3*256 = 192
        int diffLimit = Math.round(0.25f * nBands * 256);
        
        // visit lowest costs (== differences) first
        for (int i = 0; i < nTot; ++i) {
            if (costs[i] > diffLimit 
                //|| count > 45
                ) {
                break;
            }
            PairInt index12 = indexes[i];
            int idx1 = index12.getX();
            int idx2 = index12.getY();
            PairInt p1 = keypoints1.get(idx1);
            PairInt p2 = keypoints2.get(idx2);
            if (set1.contains(p1) || set2.contains(p2)) {
                continue;
            }
            int lIdx2 = pointIndexMap.get(p2);
            if (segIndexCount.containsKey(lIdx2) && segIndexCount.get(lIdx2) 
                > sLimit) {
                continue;
            }
            //System.out.println("p1=" + p1 + " " + " p2=" + p2 + " cost=" + costs[i]);
            mT.add(p1.getX(), p1.getY());
            mS.add(p2.getX(), p2.getY());
            set1.add(p1);
            set2.add(p2);
            tIndexes.add(idx1);
            
            if (segIndexCount.containsKey(lIdx2)) {
                int c = segIndexCount.get(lIdx2);
                segIndexCount.put(lIdx2, c + 1);
            } else {
                segIndexCount.put(lIdx2, 1);
            }
            
            TObjectIntMap<PairInt> cMap = idx1P2CostMap.get(idx1);
            if (cMap == null) {
                cMap = new TObjectIntHashMap<PairInt>();
                idx1P2CostMap.put(idx1, cMap);
            }
            if (!cMap.containsKey(p2)) {
                // only store the smaller cost, that's reached first
                cMap.put(p2, costs[i]);
            }
            
            count++;
        }
                
        int[] minMaxXY = MiscMath.findMinMaxXY(keypoints1);
        int objDimension = Math.max(minMaxXY[1] - minMaxXY[0],
            minMaxXY[3] - minMaxXY[2]);
        int limit = Math.round(scaleFactor * objDimension);
        int limitSq = limit * limit;
       
        int[] minMaxXY2 = MiscMath.findMinMaxXY(set2);
        
        NearestNeighbor2D nn = new NearestNeighbor2D(
           set2, minMaxXY2[1] + limit, minMaxXY2[3] + limit);
        
        return completeUsingCombinations(
            keypoints1, keypoints2, mT, mS, nn,
            minMaxXY2, limit,
            tIndexes, idx1P2CostMap, indexes, costs
        );            
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
    private static CorrespondenceList completeUsingCombinations(
        List<PairInt> keypoints1, List<PairInt> keypoints2,
        PairIntArray mT, PairIntArray mS, 
        NearestNeighbor2D nn, int[] minMaxXY2, int limit,
        TIntList tIndexes,        
        TIntObjectMap<TObjectIntMap<PairInt>> idx1P2CostMap,
        PairInt[] indexes, int[] costs
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
        
        double minCost = Double.MAX_VALUE;
        CorrespondenceList minCostCor = null;
        PairIntArray minCostTrT = null;
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
                
                double maxCost = 3 * 256;
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
                if (sum < minCost) {
                    
                    minCost = sum;
                    
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
                    
                    minCostCor = corr;
                    
                    minCostTrT = trT;
                }
            }
        }
        
        //minCostCor = setToBestUnique(
        //    minCostCor, minCostI, minDistI);
        
        return minCostCor; 
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

        assert(d1.length == d2.length);
        
        int n1 = d1[0].descriptors.length;
        int n2 = d2[0].descriptors.length;
        
        assert(n1 == keypoints1.size());
        assert(n2 == keypoints2.size());

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

        //TODO: need to add use of auto-correlation too
        
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
}
