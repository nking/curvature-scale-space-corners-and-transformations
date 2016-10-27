package algorithms.imageProcessing.features;

import algorithms.QuickSort;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.MedianTransform;
import algorithms.imageProcessing.StructureTensor;
import algorithms.util.PairInt;
import algorithms.util.TwoDFloatArray;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.TFloatList;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.TObjectIntMap;
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
 
This version is adapted to work with segmented cells and is using
* experimental descriptor masks of 0 for pixels not in the 
* segmented cell of the keypoint.
* 
 */
public class SegmentedORB {

    // these could be made static across all instances, but needs guards for synchronous initialization
    protected int[][] OFAST_MASK = null;
    protected int[] OFAST_UMAX = null;
    protected int[][] POS0 = null;
    protected int[][] POS1 = null;

    protected int fastN = 9;
    protected float fastThreshold = 0.08f;
    protected final float harrisK = 0.04f;
    protected final int nKeypoints;

    //NOTE: may want to edit the segmented cells in the future so that
    //   a single width boundary of pixels between cells is shared by
    //   the cells and that a keypoint on it has 2 different descriptors,
    //   one for its presence in each cell.
    protected final List<Set<PairInt>> segmentedCells;
    protected final TObjectIntMap<PairInt> pointSegmentedIndexMap;
 
    //list index is scale.  map key=index, value=set of keypoints in column major format
    protected List<TIntObjectMap<List<PairInt>>> keypointsList = null;
    
    private List<TIntObjectMap<Descriptors>> descriptorsList = null;
    private List<TIntObjectMap<Descriptors>> descriptorsListH = null;
    private List<TIntObjectMap<Descriptors>> descriptorsListS = null;
    private List<TIntObjectMap<Descriptors>> descriptorsListV = null;
    
    
    protected static double twoPI = 2. * Math.PI;

    protected float curvatureThresh = 0.05f;

    protected static enum DescriptorChoice {
        NONE, GREYSCALE, HSV
    }
    protected DescriptorChoice descrChoice = DescriptorChoice.GREYSCALE;
    
    protected boolean doCreate1stDerivKeypoints = false;
    
    protected boolean doCreateCurvatureKeyPoints = false;

    private final Image image;
    
    /**
     * Still testing the class, there may be bugs present.
     * @param nKeypoints
     */
    public SegmentedORB(int nKeypoints, Image img, List<Set<PairInt>> segmentedCells) {

        initMasks();

        this.nKeypoints = nKeypoints;
    
        this.image = img;
      
        this.segmentedCells = segmentedCells;
        
        this.pointSegmentedIndexMap = new TObjectIntHashMap<PairInt>();
        
        for (int i = 0; i < segmentedCells.size(); ++i) {
            Set<PairInt> set = segmentedCells.get(i);
            for (PairInt p : set) {
                pointSegmentedIndexMap.put(p, i);
            }
        }
        
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

        List<TwoDFloatArray> pyramid = buildPyramid(image);
     
        // these are filtered to only contain points in segmented cells.
        List<TIntList> keypoints0List = new ArrayList<TIntList>();
        List<TIntList> keypoints1List = new ArrayList<TIntList>();
        List<TDoubleList> orientationsList = new ArrayList<TDoubleList>();
        List<TFloatList> harrisResponses = new ArrayList<TFloatList>();
        List<TFloatList> scalesList = new ArrayList<TFloatList>();
        
        TIntList octavesUsed = extractKeypoints(pyramid, keypoints0List, 
            keypoints1List, orientationsList, harrisResponses, scalesList);
        
        int nKeyPointsTotal = 0;
        for (int i = 0; i < keypoints1List.size(); ++i) {
            nKeyPointsTotal += keypoints1List.get(i).size();
        }
        
        System.out.println("nKeyPointsTotal=" + nKeyPointsTotal +
            " this.nKeypoints=" + this.nKeypoints);
    
        this.keypointsList = new ArrayList<TIntObjectMap<List<PairInt>>>();
        
        // if descrDithers is set to other than none, this increases the
        // members of instance variables such as keypoint lists too
        extractDescriptors(image, pyramid, octavesUsed,
            keypoints0List, keypoints1List, orientationsList, 
            harrisResponses, scalesList);  
        
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
                }
            }

            QuickSort.sortBy1stArg(scores, indexes);

            //TODO: this could be improved.  
            //  it's using the response strengths to rank
            //  the results and keep the number of results requested.
            //  then the segmentation map keypoints and descriptors
            //  are re-written from those.
           
            TIntObjectMap<TIntList> keep = new TIntObjectHashMap<TIntList>();

            // visit largest scores to smallest
            count = 0;
            int idx = n - 1;
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

            List<TIntList> keypoints0List2 = new ArrayList<TIntList>();
            List<TIntList> keypoints1List2 = new ArrayList<TIntList>();
            List<TDoubleList> orientationsList2 = new ArrayList<TDoubleList>();
            List<TFloatList> harrisResponses2 = new ArrayList<TFloatList>();
            List<TFloatList> scalesList2 = new ArrayList<TFloatList>();
             
            for (int idx1 = 0; idx1 < keypoints0List.size(); ++idx1) {
                TIntList keepList = keep.get(idx1);
                if (keepList == null) {
                    continue;
                }
                int n2 = keepList.size();
                TIntList kp0 = new TIntArrayList(n2);
                TIntList kp1 = new TIntArrayList(n2);
                TDoubleList or = new TDoubleArrayList(n2);
                TFloatList hr = new TFloatArrayList(n2);
                TFloatList s = new TFloatArrayList(n2);
                TIntList m = new TIntArrayList(n2);
                
                keypoints0List2.add(kp0);
                keypoints1List2.add(kp1);
                orientationsList2.add(or);
                harrisResponses2.add(hr);
                scalesList2.add(s);
            }

            keypoints0List = keypoints0List2;
            keypoints1List = keypoints1List2;
            orientationsList = orientationsList2;
            harrisResponses = harrisResponses2;
            scalesList = scalesList2;
            
            // --- rewrite the segmentation variables ----
            for (int ii = 0; ii < keypoints0List.size(); ++ii) {
                
                Set<PairInt> sPoints = extractColMajPoints(keypoints0List.get(ii),
                    keypoints1List.get(ii));
              
                TIntObjectMap<List<PairInt>> kpToFilter = keypointsList.get(ii);
                      
                TIntObjectIterator<List<PairInt>> iter = kpToFilter.iterator();
                while (iter.hasNext()) {
                    
                    int groupIdx = iter.key();
                    
                    TIntList keep2 = new TIntArrayList();
                    
                    List<PairInt> points = iter.value();
                    for (int k = (points.size() - 1); k > -1; --k) {
                        if (sPoints.contains(points.get(k))) {
                            keep2.add(k);
                        } else {
                            points.remove(k);
                        }
                    }
                    
                    // --- edit descriptors to keep 'keep2' ---
                    if (descrChoice.equals(DescriptorChoice.HSV)) {
                        
                        Descriptors descH = descriptorsListH.get(ii).get(groupIdx);
                        Descriptors descS = descriptorsListS.get(ii).get(groupIdx);
                        Descriptors descV = descriptorsListV.get(ii).get(groupIdx);
                       
                        int[] dH = new int[keep2.size()];
                        int[] dS = new int[dH.length];
                        int[] dV = new int[dH.length];
                        int count2 = 0;
                        for (int k = 0; k < keep2.size(); ++k) {
                            int idx2 = keep2.get(k);
                            dH[count2] = descH.descriptors[idx2];
                            dS[count2] = descS.descriptors[idx2];
                            dV[count2] = descV.descriptors[idx2];
                            count2++;
                        }
                        descH.descriptors = dH;
                        descS.descriptors = dS;
                        descV.descriptors = dV;
                        
                    } else if (!descrChoice.equals(DescriptorChoice.NONE)) {
                        
                        Descriptors desc = descriptorsList.get(ii).get(groupIdx);
                        
                        int[] d = new int[keep2.size()];
                        int count2 = 0;
                        for (int k = 0; k < keep2.size(); ++k) {
                            int idx2 = keep2.get(k);
                            d[count2] = desc.descriptors[idx2];
                            count2++;
                        }
                        desc.descriptors = d;
                    }
                }                
            }
        }
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
     
        int decimationLimit = 8;

        List<GreyscaleImage> outputR = new ArrayList<GreyscaleImage>();
        MedianTransform mt = new MedianTransform();
        mt.multiscalePyramidalMedianTransform2(imageR, outputR, decimationLimit);

        List<GreyscaleImage> outputG = new ArrayList<GreyscaleImage>();
        mt.multiscalePyramidalMedianTransform2(imageG, outputG, decimationLimit);
        
        List<GreyscaleImage> outputB = new ArrayList<GreyscaleImage>();
        mt.multiscalePyramidalMedianTransform2(imageB, outputB, decimationLimit);
        
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

        int decimationLimit = 8;

        ImageProcessor imageProcessor = new ImageProcessor();
        
        List<GreyscaleImage> output = new ArrayList<GreyscaleImage>();

        MedianTransform mt = new MedianTransform();
        mt.multiscalePyramidalMedianTransform2(img, output, decimationLimit);

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

    private TIntList extractKeypoints(List<TwoDFloatArray> pyramid,
        List<TIntList> keypoints0List, List<TIntList> keypoints1List,
        List<TDoubleList> orientationsList, 
        List<TFloatList> harrisResponses, List<TFloatList> scalesList) {
   
        float prevScl = 1;
        
        TIntList octavesUsed = new TIntArrayList();

        for (int octave = 0; octave < pyramid.size(); ++octave) {

            System.out.println("octave=" + octave);

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
            
            if (octaveImage.length < 8 || octaveImage[0].length < 8) {
                continue;
            }

            // note, method filters to keep only those in segmented cells
            Resp r = detectOctave(octaveImage);

            if ((r == null) || (r.keypoints0 == null) || r.keypoints0.isEmpty()) {
                continue;
            }
            
            octavesUsed.add(octave);

            System.out.println("  octave " + octave + " nKeypoints="
                + r.keypoints0.size());
        
            //mask of length orientations.size() containing a 1 or 0
            // indicating if pixels are within the image (``True``) or in the
            // border region of the image (``False``).
            TIntList mask = filterNearBorder(octaveImage, r.keypoints0,
                r.keypoints1, r.orientations, r.responses);
            
            for (int i = 0; i < r.keypoints0.size(); ++i) {
                int v = Math.round(scale * r.keypoints0.get(i));
                r.keypoints0.set(i, v);
                v = Math.round(scale * r.keypoints1.get(i));
                r.keypoints1.set(i, v);
            }

            keypoints0List.add(r.keypoints0);
            keypoints1List.add(r.keypoints1);
            orientationsList.add(r.orientations);
            harrisResponses.add(r.responses);

            TFloatList scales = new TFloatArrayList(r.keypoints0.size());
            for (int i = 0; i < r.keypoints0.size(); ++i) {
                scales.add(scale);
            }
            scalesList.add(scales);
        }

        return octavesUsed;        
    }
    
    private void extractDescriptors(Image image, List<TwoDFloatArray> pyramid,
        TIntList octavesUsed,
        List<TIntList> keypoints0List, List<TIntList> keypoints1List,
        List<TDoubleList> orientationsList, 
        List<TFloatList> harrisResponses, List<TFloatList> scalesList
        ) {
                
        if (!descrChoice.equals(DescriptorChoice.HSV)) {
            descriptorsList = new ArrayList<TIntObjectMap<Descriptors>>();
            extractGreyscaleDescriptors(pyramid, octavesUsed, 
                keypoints0List, keypoints1List,
                orientationsList, harrisResponses, scalesList,                
                descriptorsList, true);
            return;
        }
        
        List<TwoDFloatArray> pyramidH = new ArrayList<TwoDFloatArray>();
        List<TwoDFloatArray> pyramidS = new ArrayList<TwoDFloatArray>();
        List<TwoDFloatArray> pyramidV = new ArrayList<TwoDFloatArray>();
        
        buildHSVPyramid(image, pyramidH, pyramidS, pyramidV);
                
        descriptorsListH = new ArrayList<TIntObjectMap<Descriptors>>();
        extractGreyscaleDescriptors(pyramidH, octavesUsed, 
            keypoints0List, keypoints1List,
            orientationsList, harrisResponses, scalesList,                
            descriptorsListH, true);
        descriptorsListS = new ArrayList<TIntObjectMap<Descriptors>>();
        extractGreyscaleDescriptors(pyramidS, octavesUsed, 
            keypoints0List, keypoints1List,
            orientationsList, harrisResponses, scalesList,                
            descriptorsListS, false);
        descriptorsListV = new ArrayList<TIntObjectMap<Descriptors>>();
        extractGreyscaleDescriptors(pyramidV, octavesUsed, 
            keypoints0List, keypoints1List,
            orientationsList, harrisResponses, scalesList,                
            descriptorsListV, false);
    }

    protected void extractGreyscaleDescriptors(
        List<TwoDFloatArray> pyramid, TIntList octavesUsed,
        List<TIntList> keypoints0List, List<TIntList> keypoints1List,
        List<TDoubleList> orientationsList, 
        List<TFloatList> harrisResponses, List<TFloatList> scalesList,
        List<TIntObjectMap<Descriptors>> descMap,
        boolean populateKeypointList) {
        
        TIntSet octavesUsedSet = new TIntHashSet(octavesUsed);
        
        int listIdx = 0;

        for (int octave = 0; octave < pyramid.size(); ++octave) {

            if (!octavesUsedSet.contains(octave)) {
                continue;
            }
            
            System.out.println("octave=" + octave);

            float[][] octaveImage = pyramid.get(octave).a; 

            TFloatList scales = scalesList.get(listIdx); 
            TIntList kp0 = keypoints0List.get(listIdx);
            TIntList kp1 = keypoints1List.get(listIdx);
            TDoubleList or = orientationsList.get(listIdx);
            TFloatList harrisResp = harrisResponses.get(listIdx);

         
            TIntObjectMap<Descriptors> indexDescMap = new
                TIntObjectHashMap<Descriptors>();
    
            TIntObjectMap<List<PairInt>> keypointMap = null;
            
            if (populateKeypointList) {
                keypointMap = new TIntObjectHashMap<>();
            }
            
            // ---- further group these by segmented cell -------
            TIntObjectMap<TIntList> indexesList = groupByCell(kp0, kp1);
            TIntObjectIterator<TIntList> iter = indexesList.iterator();
            for (int j = 0; j < indexesList.size(); ++j) {
                
                iter.advance();
                
                // segmented cell index:
                int groupIdx = iter.key();
                
                // indexes in above lists:
                TIntList indexes = iter.value();
                
                if (indexes.isEmpty()) {
                    continue;
                }
                assert(groupIdx == pointSegmentedIndexMap.get(
                    new PairInt(
                        kp1.get(indexes.get(0)), kp0.get(indexes.get(0))
                    )));
                
                // result contains descriptors and mask.
                // also, modified by mask are the keypoints and orientations
                Descriptors desc = extractOctave(octaveImage, kp0, kp1, or, 
                    indexes, groupIdx);

                indexDescMap.put(groupIdx, desc);
                
                if (populateKeypointList) {
                    List<PairInt> kpList = new ArrayList<PairInt>();
                    for (int ii = 0; ii < indexes.size(); ++ii) {
                        int y = kp0.get(ii);
                        int x = kp1.get(ii);
                        PairInt p = new PairInt(x, y);
                        kpList.add(p);
                    }
                    keypointMap.put(groupIdx, kpList);
                }
            }
           
            descMap.add(indexDescMap);
            
            if (populateKeypointList) {
                this.keypointsList.add(keypointMap);
            }
            
            // octave list index:
            listIdx++;
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

    private class Resp {
        TDoubleList orientations;
        TIntList keypoints0;
        TIntList keypoints1;
        TFloatList responses;
    }

    public static class Descriptors {
        //NOTE: consider packing 4 descriptors into one or
        // doing as Image.java does with sensing 64 bit and 32 bit to make
        // long or int bit vectors
        int[] descriptors;
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

        TIntList keypoints0 = new TIntArrayList();
        TIntList keypoints1 = new TIntArrayList();

        // list of format [row, col, ...] of filtered maxima ordered by intensity
        cornerPeaks(fastResponse, 1, keypoints0, keypoints1);
        if (keypoints0.isEmpty()) {
            return null;
        }
        
        filterForSegmentation(keypoints0, keypoints1);

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
        
            filterForSegmentation(keypoints0, keypoints1);
            
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
                algorithms.imageProcessing.ImageDisplayer.displayImage("curvature", img2);
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
    
            filterForSegmentation(keypoints0, keypoints1);
            
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

        TFloatList responses = new TFloatArrayList(keypoints0.size());
        for (int i = 0; i < keypoints0.size(); ++i) {
            int x = keypoints0.get(i);
            int y = keypoints1.get(i);
            float v = harrisResponse[x][y];
            responses.add(v);
        }
        
        double[] orientations2 = cornerOrientations(octaveImage,
            keypoints0, keypoints1);
        assert(orientations2.length == keypoints0.size());
        assert(orientations2.length == keypoints1.size());

        
        Resp r2 = new Resp();
        r2.keypoints0 = keypoints0;
        r2.keypoints1 = keypoints1;
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

        System.out.println("nBefore border rm=" + nBefore +
            " nAfter=" + set.size());

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
    protected double[] cornerOrientations(float[][] octaveImage,
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

        double[] orientations = new double[nCorners];

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
            orientations[i] = Math.atan2(m01, m10);
            
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
     * @param orientations
     * @param responses
     * @return the encapsulated descriptors and mask
     */
    protected TIntList filterNearBorder(float[][] octaveImage,
        TIntList keypoints0, TIntList keypoints1,
        TDoubleList orientations, TFloatList responses) {
        
        if (POS0 == null) {
            POS0 = ORBDescriptorPositions.POS0;
        }
        if (POS1 == null) {
            POS1 = ORBDescriptorPositions.POS1;
        }

        assert(orientations.size() == keypoints0.size());
        assert(orientations.size() == keypoints1.size());
        assert(orientations.size() == responses.size());

        int nRows = octaveImage.length;
        int nCols = octaveImage[0].length;

        //mask of length orientations.size() containing a 1 or 0
        // indicating if pixels are within the image (``True``) or in the
        // border region of the image (``False``).
        TIntList mask = maskCoordinatesIfBorder(keypoints0, keypoints1,
            nRows, nCols, 5);//16);

        assert(orientations.size() == keypoints0.size());
        assert(orientations.size() == keypoints1.size());
        assert(orientations.size() == responses.size());
        assert(orientations.size() == mask.size());

        /*
        i    0    1    2
        idx2 0,1  2,3  4,5

        looping backwards to keep indexes correct
        */
        for (int i = (mask.size() - 1); i > -1; --i) {
            if (mask.get(i) == 0) {
                keypoints0.removeAt(i);
                keypoints1.removeAt(i);
                orientations.removeAt(i);
                responses.removeAt(i);
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
     * @param indexes
     * @param groupIdx
     * @return the encapsulated descriptors and mask
     */
    protected Descriptors extractOctave(float[][] octaveImage,
        TIntList keypoints0, TIntList keypoints1,
        TDoubleList orientations,
        TIntList indexes, int groupIdx) {
        
        if (POS0 == null) {
            POS0 = ORBDescriptorPositions.POS0;
        }
        if (POS1 == null) {
            POS1 = ORBDescriptorPositions.POS1;
        }

        assert(orientations.size() == keypoints0.size());
        assert(orientations.size() == keypoints1.size());

        int[] descriptors = null;

        if (descrChoice.equals(DescriptorChoice.NONE)) {
            descriptors = new int[0];
        } else {      
            descriptors = orbLoop(octaveImage, keypoints0, keypoints1,
                orientations, indexes, groupIdx);
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
    protected int[] orbLoop(float[][] octaveImage, TIntList keypoints0,
        TIntList keypoints1, TDoubleList orientations,
        TIntList indexes, int groupIdx) {

        assert(orientations.size() == keypoints0.size());

        if (POS0 == null) {
            POS0 = ORBDescriptorPositions.POS0;
        }
        if (POS1 == null) {
            POS1 = ORBDescriptorPositions.POS1;
        }

        //int nKP = orientations.size();
        int nKP = indexes.size();
        
        System.out.println("nKP=" + nKP);

        // holds values 1 or 0.  size is [orientations.size] 
        int[] descriptors = new int[nKP];

        double pr0, pc0, pr1, pc1;
        int spr0, spc0, spr1, spc1;

        for (int i = 0; i < indexes.size(); ++i) {
            int idx = indexes.get(i);
            double angle = orientations.get(idx);
            double sinA = Math.sin(angle);
            double cosA = Math.cos(angle);

            int kr = keypoints0.get(idx);
            int kc = keypoints1.get(idx);
            
            int descr = 0;

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
                PairInt p0 = new PairInt(y0, x0);
                PairInt p1 = new PairInt(y1, x1);
                
                float v0, v1;
                if (pointSegmentedIndexMap.containsKey(p0) && 
                    pointSegmentedIndexMap.get(p0) == groupIdx) {
                    v0 = octaveImage[x0][y0];
                } else {
                    v0 = 0;
                }
                if (pointSegmentedIndexMap.containsKey(p1) && 
                    pointSegmentedIndexMap.get(p1) == groupIdx) {
                    v1 = octaveImage[x1][y1];
                } else {
                    v1 = 0;
                }
                
                if (v0 < v1) {
                    descr |= (1 << j);
                }
            }
            descriptors[i] = descr;
        }

        return descriptors;
    }

    public List<TIntObjectMap<List<PairInt>>> getJeypointsList() {
        return keypointsList;
    }
    
    public List<TIntObjectMap<Descriptors>> getDescriptorsList() {
        return descriptorsList;
    }
    public List<TIntObjectMap<Descriptors>> getDescriptorsListH() {
        return descriptorsListH;
    }
    public List<TIntObjectMap<Descriptors>> getDescriptorsListS() {
        return descriptorsListS;
    }
    public List<TIntObjectMap<Descriptors>> getDescriptorsListV() {
        return descriptorsListV;
    }
   
    public TIntObjectMap<Descriptors> getAllDescriptorsPerCell() {
        
        if (descriptorsList == null) {
            return null;
        }
        
        return combineByScale(descriptorsList);
    }
    public TIntObjectMap<Descriptors> getAllDescriptorsHPerCell() {
        
        if (descriptorsListH == null) {
            return null;
        }
        
        return combineByScale(descriptorsListH);
    }
    public TIntObjectMap<Descriptors> getAllDescriptorsSPerCell() {
        
        if (descriptorsListS == null) {
            return null;
        }
        
        return combineByScale(descriptorsListS);
    }
    public TIntObjectMap<Descriptors> getAllDescriptorsVPerCell() {
        
        if (descriptorsListV == null) {
            return null;
        }
        
        return combineByScale(descriptorsListV);
    }
    
    private TIntObjectMap<Descriptors> combineByScale(
        List<TIntObjectMap<Descriptors>> list) {
        
        if (list == null) {
            return null;
        }
        
        TIntObjectMap<Descriptors> out = new TIntObjectHashMap<Descriptors>();
    
        TIntSet groupIndexes = new TIntHashSet();
        for (int i = 0; i < list.size(); ++i) {
            TIntObjectMap<Descriptors> map = list.get(i);
            groupIndexes.addAll(map.keySet());
        }
        
        TIntIterator iter = groupIndexes.iterator();
        while (iter.hasNext()) {
            
            int groupIdx = iter.next();
            
            int n = 0;
            for (int i = 0; i < list.size(); ++i) {
                Descriptors desc = list.get(i).get(groupIdx);
                if (desc != null && desc.descriptors != null) {
                    n += desc.descriptors.length;
                }
            }
            
            if (n == 0) {
                continue;
            }
            
            int[] d = new int[n];
            
            n = 0;
            for (int i = 0; i < list.size(); ++i) {
                
                Descriptors desc = list.get(i).get(groupIdx);
        
                if (desc != null && desc.descriptors != null) {
                    
                    int nLen = desc.descriptors.length;
                    
                    System.arraycopy(desc.descriptors, 0, d, n, nLen);
                    
                    n += nLen;
                }
            }
            
            Descriptors desc = new Descriptors();
            desc.descriptors = d;
            
            out.put(groupIdx, desc);
        }
        
        return out;
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
       
        int[] combinedD = new int[n];

        int count = 0;
        for (int i = 0; i < list.size(); ++i) {

            int[] d = list.get(i).descriptors;
            
            System.arraycopy(d, 0, combinedD, count, d.length);
            count += d.length;
        }

        Descriptors combined = new Descriptors();
        combined.descriptors = combinedD;

        return combined;
    }

    private void filterForSegmentation(TIntList keypoints0, 
        TIntList keypoints1) {
        
        TIntList rm = new TIntArrayList();
        for (int i = 0; i < keypoints0.size(); ++i) {
            PairInt p = new PairInt(keypoints1.get(i), keypoints0.get(i));
            if (!pointSegmentedIndexMap.containsKey(p)) {
                rm.add(i);
            }
        }
        
        for (int i = (rm.size() - 1); i > -1; --i) {
            int idx = rm.get(i);
            keypoints0.removeAt(idx);
            keypoints1.removeAt(idx);
        }
    }
    
    private TIntObjectMap<TIntList> groupByCell(TIntList kp0, TIntList kp1) {
     
        //NOTE: kp0,kp1 is row major, but resturned points are col major
        // also note that keypoints have already been filtered to contain only 
        // those present in segmented cells
        
        TIntObjectMap<TIntList> output = new TIntObjectHashMap<TIntList>();
    
        for (int i = 0; i < kp0.size(); ++i) {
            int y = kp0.get(i);
            int x = kp1.get(i);
            PairInt p = new PairInt(x, y);
        
            assert(pointSegmentedIndexMap.containsKey(p));
            int idx = pointSegmentedIndexMap.get(p);
            
            TIntList kpList = output.get(idx);
            if (kpList == null) {
                kpList = new TIntArrayList();
                output.put(idx, kpList);
            }
            kpList.add(i);
        }
        
        return output;
    }

    private Set<PairInt> extractColMajPoints(TIntList kp0, TIntList kp1) {

        Set<PairInt> out = new HashSet<PairInt>();
        
        for (int i = 0; i < kp0.size(); ++i) {
            
            PairInt p = new PairInt(kp1.get(i), kp0.get(i));
                            
            out.add(p);
        }
        
        return out;
    }

}
