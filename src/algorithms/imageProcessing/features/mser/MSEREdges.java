package algorithms.imageProcessing.features.mser;

import algorithms.QuickSort;
import algorithms.compGeometry.PerimeterFinder2;
import algorithms.imageProcessing.CannyEdgeColorAdaptive;
import algorithms.imageProcessing.ConnectedPointsFinder;
import algorithms.imageProcessing.DFSContiguousValueFinder;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.GroupPixelHSV;
import algorithms.imageProcessing.GroupPixelHSV2;
import algorithms.imageProcessing.Heap;
import algorithms.imageProcessing.HeapNode;
import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.ImageSegmentation;
import algorithms.imageProcessing.MiscellaneousCurveHelper;
import algorithms.imageProcessing.SummedAreaTable;
import algorithms.imageProcessing.features.HCPT;
import algorithms.imageProcessing.features.HGS;
import algorithms.imageProcessing.features.mser.Canonicalizer.RegionGeometry;
import algorithms.imageProcessing.features.mser.MSER.Threshold;
import algorithms.imageProcessing.util.AngleUtil;
import algorithms.misc.Misc;
import algorithms.misc.MiscDebug;
import algorithms.misc.MiscMath;
import algorithms.util.OneDIntArray;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.VeryLongBitString;
import gnu.trove.iterator.TIntIntIterator;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import thirdparty.edu.princeton.cs.algs4.Interval;
import thirdparty.edu.princeton.cs.algs4.Interval2D;
import thirdparty.edu.princeton.cs.algs4.QuadTree;

/**
 * class to explore the boundaries of the accumulated points in MSER regions
 * to make edges that are complete contours (the level sets approach has a
 * much larger pixel association than the few pixels of canny edges for
 * example).
 *
 * It currently uses the level sets found in the MSER =made from
 * greyscale and "H" of LCH color space images (caveat, configured
 * to extract MSER regions from images from COTS of past several years
 * binned to size near 256 X 256).
 * The boundaries of the accumulated points in the regions are extracted
 * and the high scoring points are kept where score is a combination
 * of sobel color contrast and sobel greyscale intensity values).
 *
 * The class in not ready for use yet.
 *
 * NOTE that the default mode may change to the low contrast settings.
 *
 * @author nichole
 */
public class MSEREdges {

    /*
    given a color image
      -- makes greyscale image and a polar theta of cie luv
      -- creates positive and negative mser regions for both,
            filtering by variation for the specific mser type
    extractEdges:
        -- uses PerimeterFinder2 extract the outer boundaries from each mser
           region points, and the same for the contiguous embedded regions.
    mergeAndExtractEdges
        -- merges adjacent regions by similarity of color (and maybe texture in future)
        -- extracts edges using PerimeterFinder2 just as above

    */

    private enum STATE {
        INITIALIZED, REGIONS_EXTRACTED, MERGED, EDGES_EXTRACTED
    }

    private final ImageExt clrImg;
    private final GreyscaleImage gsImg;
    private final GreyscaleImage ptImg;
    // ptImg shifted by +60 degrees in pixel values to look for wrap around affects
    private final GreyscaleImage ptImgShifted;
    private List<Region> regions = null;

    //a re-extraction of regions on the gs positive image, but using the
    // MSER parameters used for the negative image.
    // this is not useful in edges, but is useful for other methods.
    private List<Region> sensitiveGS0 = null;

    // the original regions for gs positive, negative, then pt positive
    // and negative.  regions should be prefered for most
    // uses.
    private List<List<Region>> origGsPtRegions = null;

    // NOTE: the edges indexes do not correspond to the regions indexes
    // list of sets of pixel indexes of boundaries of labeled sets
    private List<TIntSet> edgeList = null;

    // list of sets of pixel indexes
    private List<TIntSet> labeledSets = null;

    private boolean debug = false;

    private boolean useLowerContrastLimits = false;

    private GreyscaleImage sobelScores = null;

    private GreyscaleImage cannyEdges = null;

    private long ts = 0;

    private STATE state = null;

    public MSEREdges(ImageExt img) {

        ImageProcessor imageProcesor = new ImageProcessor();

        this.gsImg = img.copyToGreyscale2();
        this.ptImg = imageProcesor.createCIELUVTheta(img, 255);
        this.clrImg = img.copyToImageExt();

        //to correct for wrap around effects in ptImg, making a copy
        //   shifted by -20 degrees in pixel values.
        //   any false regions found due to the range 0-255 being wrap
        //   around from 255 to 0, but perceived as different can be
        //   found by comparing these two image's results, though
        //   there may be faster ways to achieve this.
        this.ptImgShifted = copyAndShift(ptImg, 60);

        state = STATE.INITIALIZED;
    }

    public void setToDebug() {
        debug = true;
        ts = MiscDebug.getCurrentTimeFormatted();
    }

    /**
     * use this to change the limits to keep lower contrast edges.
     * NOTE that this increases the number of edges that are noise.
     */
    public void setToLowerContrast() {
        useLowerContrastLimits = true;
    }

    public void extractEdges() {

        //INITIALIZED, REGIONS_EXTRACTED, MERGED, EDGES_EXTRACTED
        if (!state.equals(STATE.INITIALIZED)) {
            throw new IllegalStateException("can only perform extraction of "
                + "edges once");
        }

        if (debug) {
            MiscDebug.writeImage(ptImg, "_" + ts + "_polartheta");
        }

        extractRegions();

        extractBoundaries();

        assert(labeledSets != null);
        assert(edgeList != null);

        if (debug) {
            printEdges();
        }
    }

    /**
     * merge regions and then extract edges.  this method is most useful when
     * setToLowerContrast() has been invoked before.
     */
    public void mergeAndExtractEdges() {

        //INITIALIZED, REGIONS_EXTRACTED, MERGED, EDGES_EXTRACTED
        if (!state.equals(STATE.INITIALIZED)) {
            throw new IllegalStateException("can only perform extraction of "
                + "edges once");
        }

        extractEdges();

        long ts0 = System.currentTimeMillis();

        mergeRegions3();

        long ts1 = System.currentTimeMillis();

        System.out.format("%.3f sec for merge\n",
            (((float)ts1 - ts0)/1000.f));

        if (debug) {
            printEdges();
        }
    }

    private void extractRegions() {

        //INITIALIZED, REGIONS_EXTRACTED, MERGED, EDGES_EXTRACTED
        if (state.equals(STATE.EDGES_EXTRACTED)) {
            throw new IllegalStateException("can only perform extraction of "
                + "edges once");
        }

        List<List<Region>> gsRegions = extractMSERRegions(gsImg, Threshold.LESS_SENSITIVE);

        List<List<Region>> ptRegions = extractMSERRegions(ptImg, Threshold.DEFAULT);

        List<List<Region>> ptShiftedRegions = extractMSERRegions(
            ptImgShifted, Threshold.DEFAULT);

        //_debugOrigRegions(ptShiftedRegions.get(0), "_shifted_0");
        //_debugOrigRegions(ptShiftedRegions.get(1), "_shifted_1");

        regions = new ArrayList<Region>();

        origGsPtRegions = new ArrayList<List<Region>>();

        int w = ptImg.getWidth();
        int h = ptImg.getHeight();

        // for minArea 0.001, a gsRegions limit for var of 0. is used
        //   but for 0.0001, limit should be near 0.001

        for (int type = 0; type < 2; ++type) {
            List<Region> list = gsRegions.get(type);

            for (int i = (list.size() - 1); i > -1; --i) {
                Region r = list.get(i);
                if ((type == 1) && r.getVariation() > 2.) {
                    list.remove(i);
                } else if ((type == 0)
                    //&& r.getVariation() == 0.0) {
                    && r.getVariation() < 2.) {
                    list.remove(i);
                }
            }

            // copy list into origGsPtRegions
            List<Region> cpList = new ArrayList<Region>();
            origGsPtRegions.add(cpList);
            for (Region r : list) {
                cpList.add(r.copy());
                regions.add(r);
            }
        }

        for (int type = 0; type < 2; ++type) {
            List<Region> list = ptRegions.get(type);
            GreyscaleImage negImg = null;
            if (type == 1) {
                negImg = gsImg.copyImage();
                for (int i = 0; i < gsImg.getNPixels(); ++i) {
                    int v = ~gsImg.getValue(i);
                    if (v < 0) {
                        v += 256;
                    }
                    negImg.setValue(i, v);
                }
            }

            boolean hadWrapAroundArtifacts = false;
            for (int i = (list.size() - 1); i > -1; --i) {
                Region r = list.get(i);

                if ((type == 1)
                    //&& r.getVariation() > 0.001) {
                    && r.getVariation() > 2) {
                    list.remove(i);
                } else if ((type == 0)
                    //&& r.getVariation() == 0.0) {
                    && r.getVariation() < 2) {
                    list.remove(i);
                } else {
                    // check for this Region being an artifact of wrap around
                    if (r.level_ < 25) {
                        list.remove(i);
                        hadWrapAroundArtifacts = true;
                    } else {
                        int avgLevel;
                        if (type == 0) {
                            avgLevel = calcAvg(ptImg, r.getAcc(w));
                        } else {
                            avgLevel = calcAvg(negImg, r.getAcc(w));
                        }
                        //System.out.format(" %d add x,y=%d,%d level=%d  avgLevel=%d\n",
                        //    type, (int)(r.moments_[0]/r.area_),
                        //    (int)(r.moments_[1]/r.area_), r.level_,
                        //    avgLevel);
                        if (avgLevel > 240) {//230?
                            list.remove(i);
                        }
                    }
                }
            }
            if (hadWrapAroundArtifacts) {
                // excluded level < 25.
                //   so shift by +60 and any region found w level approx 60
                //   is a real region possibly excluded near level=0
                // if find a region near level=60, edit the level
                //   and add it to list
                List<Region> list2 = ptShiftedRegions.get(type);
                for (int i = (list2.size() - 1); i > -1; --i) {
                    Region r = list2.get(i);

                    if ((type == 1)
                        //&& r.getVariation() > 0.001) {
                        && r.getVariation() > 2) {
                        list2.remove(i);
                    } else if ((type == 0)
                        //&& r.getVariation() == 0.0) {
                        && r.getVariation() < 2) {
                        list2.remove(i);
                    } else {
                        //TODO: consider checking whether this already exists in
                        //   the list
                        if (Math.abs(r.level_ - 60) < 20) {
                            r.level_ += 60;
                            list.add(r);
                            //System.out.format("  add shifted x,y=%d,%d level=%d\n",
                            // (int)(r.moments_[0]/r.area_),
                            // (int)(r.moments_[1]/r.area_), r.level_);
                        }
                        //TODO: consider adding other regions in list2
                        //   as long as level is > 25 and < 230
                    }
                }
            }

            // copy list into origGsPtRegions
            List<Region> cpList = new ArrayList<Region>();
            origGsPtRegions.add(cpList);
            for (Region r : list) {
                regions.add(r);
                cpList.add(r.copy());
            }
        }

        /*
        adding better configuration for greyscale positive here.
        TODO: consider fixing the process above to include these.
        */
        regions.addAll(_extractSensitiveGS0());
        //_debugOrigRegions(regions,"_GS0__");

        if (debug) {
            int[] xyCen = new int[2];
            Image imCp;
            int n;
            List<Region> list;
            for (int type = 0; type < 2; ++type) {
                imCp = gsImg.copyToColorGreyscale();
                list = origGsPtRegions.get(type);
                n = list.size();
                for (int i = 0; i < n; ++i) {
                    Region r = list.get(i);
                    int[] clr = ImageIOHelper.getNextRGB(i);
                    r.drawEllipse(imCp, 0, clr[0], clr[1], clr[2]);
                    r.calculateXYCentroid(xyCen, imCp.getWidth(), imCp.getHeight());
                    ImageIOHelper.addPointToImage(xyCen[0], xyCen[1], imCp,
                        1, 255, 0, 0);
                }
                MiscDebug.writeImage(imCp, "_" + ts + "_regions_gs_"+ type);
            }
            for (int type = 0; type < 2; ++type) {
                imCp = ptImg.copyToColorGreyscale();
                list = origGsPtRegions.get(type + 2);
                n = list.size();
                for (int i = 0; i < n; ++i) {
                    Region r = list.get(i);
                    int[] clr = ImageIOHelper.getNextRGB(i);
                    r.drawEllipse(imCp, 0, clr[0], clr[1], clr[2]);
                    r.calculateXYCentroid(xyCen, imCp.getWidth(), imCp.getHeight());
                    ImageIOHelper.addPointToImage(xyCen[0], xyCen[1], imCp,
                        1, 255, 0, 0);
                    //System.out.println(type + " xy=" + xyCen[0] + "," + xyCen[1]
                    //    + " variation=" + r.getVariation());
                }
                MiscDebug.writeImage(imCp, "_" + ts + "_regions_pt_"+ type);
            }
        }

        state = STATE.REGIONS_EXTRACTED;
    }

    private List<List<Region>> extractMSERRegions(GreyscaleImage img,
        Threshold threshold) {

        int delta = 2;
        double minArea = 0.001;//.0001
        double maxArea = 0.99;//0.1;
        double maxVariation = 0.9;//0.5;
        double minDiversity = 0.75;//0.1;//0.5
        if (threshold.equals(Threshold.LESS_SENSITIVE)) {
            maxVariation = 2.0;
            minArea = 0.0075;
            if (useLowerContrastLimits) {
                minDiversity = 0.1;//.4
                delta = 3;//5;
            }
        }

        int[] a = MSER.readIntoArray(img);

        MSER mser = new MSER();

        List<List<Region>> regions = mser.findRegions(a, img.getWidth(),
            img.getHeight(), delta, minArea, maxArea, maxVariation,
            minDiversity);

        return regions;
    }
    
    private TIntSet extractNonEdgePixels(TIntSet edgePixels) {
        int w = clrImg.getWidth();
        int h = clrImg.getHeight();
        int n = clrImg.getNPixels();
        
        TIntSet nonEdgePixels = new TIntHashSet(n - edgePixels.size());
        for (int pixIdx = 0; pixIdx < n; ++pixIdx) {
            if (!edgePixels.contains(pixIdx)) {
                nonEdgePixels.add(pixIdx);
            }
        }
        
        return nonEdgePixels;
    }
    
    /**
     * extract the contiguous conneted pixels that are not in edgePixels
     * @param edgePixels
     * @param minGroupSize
     * @return 
     */
    private List<TIntSet> extractContiguousBetweenEdges(TIntSet 
        edgePixels) {
        
        int w = clrImg.getWidth();
        int h = clrImg.getHeight();
        int n = clrImg.getNPixels();
        
        TIntSet nonEdgePixels = extractNonEdgePixels(edgePixels);
    
        if (debug) {
            Image tmp = clrImg.createWithDimensions();
            ImageIOHelper.addCurveToImage(
                nonEdgePixels, tmp, 0, 255, 255, 0);
            MiscDebug.writeImage(tmp, "_" + ts + "_non_edge_");
        }
        
        ConnectedPointsFinder finder = new ConnectedPointsFinder(w, h);
        finder.setMinimumNumberInCluster(1);
        if (debug) {
            finder.setDebug(debug);
        }
        finder.findConnectedPointGroups(nonEdgePixels);
        
        List<TIntSet> output = new ArrayList<TIntSet>();
        
        for (int i = 0; i < finder.getNumberOfGroups(); ++i) {
            TIntSet group = finder.getXY(i);
            output.add(group);
        }
        
        return output;
    }

    private void extractBoundaries() {

        //INITIALIZED, REGIONS_EXTRACTED, MERGED, EDGES_EXTRACTED
        if (state.equals(STATE.EDGES_EXTRACTED)) {
            throw new IllegalStateException("can only perform extraction of "
                + "edges once");
        }
        
        /*
        TODO: consider a version of this which:
           -- uses the non blurred canny edges and same mLimit
              to find a small but precise subset of regions.
              -- might be able to add it's unmatched in without
                 checking gap size.
           -- then, because there are some edges missing like the
              legs of the gingerbread man in android_statues_01,
              could make one more pass through the remaining
              regions list and instead of matching to the canny edges,
              add regions which match the current set of
              allEdgePixels for a high intersection limit.
              -- will need a fraction matching score for this filter.
         Goal is to make the results of this method return a thinner
             set of edges that are still completely connected around
             objects and to do so in a faster runtime than the
             current method.
        */

        long ts0 = System.currentTimeMillis();

        TIntSet boundaries = combineBoundaries();

        long ts1 = System.currentTimeMillis();

        populateEdgeLists(boundaries);
        
        if (debug) {
            Image tmp = clrImg.copyImage();
            ImageIOHelper.addAlternatingColorPointSetsToImage2(
                labeledSets, 0, 0, 0, tmp);
            MiscDebug.writeImage(tmp, "_" + ts + "_labeled_boundaries_");
        }
        
        long ts1_2 = System.currentTimeMillis();
        
        // populate this.edgeList and this.labeledSets
        thinTheBoundaries(2);

        long ts2 = System.currentTimeMillis();

        System.out.format(
            "%.3f sec for boundary extr, "
                + " %.3f sec for reassigning and labels\n",
            ((float)(ts1 - ts0)/1000.f), ((float)(ts2 - ts1_2)/1000.f)
        );

        state = STATE.EDGES_EXTRACTED;
    }

    private void printEdges() {

        Image im = clrImg.copyImage();

        int[] clr = new int[]{255, 0, 0};

        for (int i = 0; i < edgeList.size(); ++i) {
            //int[] clr = ImageIOHelper.getNextRGB(i);
            ImageIOHelper.addCurveToImage(edgeList.get(i), im, 0, clr[0],
                clr[1], clr[2]);
        }
        MiscDebug.writeImage(im, "_" + ts + "_edges_");
        im = clrImg.copyImage();
        ImageIOHelper.addAlternatingColorCurvesToImage3(
            labeledSets, im, 0);
        MiscDebug.writeImage(im, "_" + ts + "_sets_");
    }

    private class OneDFloatArray {
        float[] a;
        public OneDFloatArray(float[] a) {
            this.a = a;
        }
    }

    /**
     * NOT READY FOR USE
     *
     * moderate merging of the labeled regions is performed
     * to remove noisey edges.  Note that the color filters
     * may need to be revised with more testing.
     *
     */
    private void mergeRegions3() {

        //TODO: this needs many edits after
        //   have finished improvements in the canny edges and the
        //   boundary extraction

        //if (true) {return;}

        /*
        TODO: refactoring this to use color and texture
        (need to add use of gradients within labeled regions to
        improve the merging.
        inspired by Alpert, Galun, Basri et al. (2007) though will
        probably use a different distance and update after merge)
        */

        //INITIALIZED, REGIONS_EXTRACTED, MERGED, EDGES_EXTRACTED
        if (!state.equals(STATE.EDGES_EXTRACTED)) {
            throw new IllegalStateException("error in algorithm.  expecting"
                + " edges were extracted.");
        }

        assert(labeledSets != null);
        assert(edgeList != null);
        assert(labeledSets.size() == edgeList.size());

        // hsv difference upper limit

        float[] hsvUL = new float[]{
             0.06f
            , 0.06f, 0.06f, 0.06f, 0.06f,0.02f
            };
        float hsvUL_green = 0.09f;
        float[] hcptLL = new float[]{
            0.7f
            ,  0.6f,  0.5f, 0.65f, 0.35f, 0.8f
            };
        float[] hgsLL = new float[]{
            0.55f
            , 0.55f, 0.5f, 0.45f, 0.8f, 0.35f
            };

        if (sobelScores == null) {
            sobelScores = createSobelScores();
        }

        HGS hgs = new HGS(sobelScores, 1, 6, 12);
        HCPT hcpt = new HCPT(ptImg, 1, 6, 12);

        TIntObjectMap<GroupPixelHSV2> clrs = new TIntObjectHashMap<GroupPixelHSV2>();

        int w = gsImg.getWidth();
        int h = gsImg.getHeight();

        for (int label = 0; label < labeledSets.size(); ++label) {
            TIntSet set = labeledSets.get(label);
            GroupPixelHSV2 hsv = new GroupPixelHSV2();
            hsv.calculateColors(set, clrImg);
            clrs.put(label, hsv);
        }

        ImageProcessor imageProcessor = new ImageProcessor();

        TIntObjectMap<TIntSet> mapOfSets = new TIntObjectHashMap<TIntSet>();
        TIntObjectMap<TIntSet> mapOfBorders = new TIntObjectHashMap<TIntSet>();
        int[] sizes = new int[labeledSets.size()];
        int[] idxs = new int[labeledSets.size()];
        for (int i = 0; i < labeledSets.size(); ++i) {
            TIntSet set = labeledSets.get(i);
            sizes[i] = set.size();
            idxs[i] = i;
            mapOfSets.put(i, set);

            set = edgeList.get(i);
            mapOfBorders.put(i, set);
        }
        assert(mapOfBorders.size() == mapOfSets.size());

        TIntIntMap pointIndexMap = createPointIndexMap(mapOfSets);

        QuickSort.sortBy1stArg(sizes, idxs);

        PerimeterFinder2 finder2 = new PerimeterFinder2();

        TIntObjectMap<VeryLongBitString> adjMap = imageProcessor
            .createAdjacencyMap(pointIndexMap, mapOfSets, w, h);

        MiscellaneousCurveHelper ch = new MiscellaneousCurveHelper();

        int nIter = 0;
        int nIterMax = 10;
        int nMerged = 0;
        int hcptIdx = 0;

        do {
            nMerged = 0;
            for (int i = (idxs.length - 1); i > -1; --i) {

                int idx1 = idxs[i];

                TIntSet set1 = mapOfSets.get(idx1);

                if (set1 == null || set1.isEmpty()) {
                    // merged already
                    continue;
                }

                TIntSet border1 = mapOfBorders.get(idx1);

                //as suggested by Alpert, Galun, Basri et al. (2007),
                //removing the boundaries from grdients
                set1 = new TIntHashSet(set1);
                set1.removeAll(border1);

                int[] hcpt1H = getRegionHistogram(hcpt, set1);
                int[] hgs1H = getRegionHistogram(hgs, set1);
                GroupPixelHSV2 hsv1 = getColors(clrs, set1, idx1);

                // for white and black, the colorspace filters need
                // specialization
                boolean isW1 = isWhite(hsv1);

                boolean isB1 = isBlack(hsv1);

                int n1 = set1.size();

                int[] xyCen1 = ch.calculateRoundedXYCentroids(border1,
                    clrImg.getWidth());

                VeryLongBitString nbrsBS = adjMap.get(idx1);

                if (nbrsBS == null) {
                    continue;
                }

                int[] idxs2 = nbrsBS.getSetBits();

                // find best merge
                // if none meet limits, this is not reset:
                float minCost = Float.MAX_VALUE;
                int minCostIdx2 = -1;

                for (int idx2 : idxs2) {

                    TIntSet set2 = mapOfSets.get(idx2);

                    // skip over set2 smaller than n1 because they've already
                    //   been compared
                    if (set2 == null || set2.size() < n1) {
                        continue;
                    }
                    TIntSet border2 = mapOfBorders.get(idx2);

                    int[] xyCen2 = ch.calculateRoundedXYCentroids(border2,
                        clrImg.getWidth());

                    //as suggested by Alpert, Galun, Basri et al. (2007),
                    //removing the boundaries from grdients
                    TIntHashSet set3 = new TIntHashSet(set2);
                    set3.removeAll(border2);

                    GroupPixelHSV2 hsv2 = getColors(clrs, set3, idx2);

                    float cost = hsv1.calculateDifference(hsv2);
                    float[] hsvDiffs = hsv1.calculateDifferences(hsv2);

                    boolean isW2 = isWhite(hsv2);

                    boolean isB2 = isW2 ? false : isBlack(hsv2);

                    if ((isB1 && isB2) || (isW1 && isW2)) {
                        // do not use hue
                        cost = hsvDiffs[2];
                    }

                    boolean bothGreen = hsv1.isGreen() && hsv2.isGreen();

                    double hsvLimit = bothGreen ? hsvUL_green :
                        hsvUL[hcptIdx];

                    if (!((isB1 && isB2) || (isW1 && isW2)) && (cost > hsvLimit)) {
                        /*System.out.format("skip (%d,%d) (%d,%d) "
                            + "hsvd=%.3f "
                            + " n=%d,%d isWh=%b,%b isBl=%b,%b "
                            + "\n    hsv1=%.3f,%.3f,%.3f"
                            + "\n    hsv2=%.3f,%.3f,%.3f\n",
                            xyCen1[0], xyCen1[1], xyCen2[0], xyCen2[1],
                            cost, set1.size(), set2.size(),
                            isW1, isW2, isB1, isB2,
                            hsv1.getAvgH(), hsv1.getAvgS(),
                            hsv1.getAvgV(),
                            hsv2.getAvgH(), hsv2.getAvgS(),
                            hsv2.getAvgV()
                        );*/
                        continue;
                    }

                    int[] hcpt2H = getRegionHistogram(hcpt, set3);
                    int[] hgs2H = getRegionHistogram(hgs, set3);

                    float hcptInter = hcpt.intersection(hcpt1H, hcpt2H);
                    float hgsInter = hgs.intersection(hgs1H, hgs2H);

                    boolean simW = isW1 && isW2 &&
                        (hcptIdx == (hcptLL.length - 1)) &&
                        ((cost < 0.065 && hcptInter > 0.475
                          && hgsInter > 0.4)
                        || (cost < 0.02 && hcptInter > 0.4)
                        || (cost < 0.01 && hcptInter > 0.15
                        && hgsInter > 0.35)
                        );

                    boolean simB = isB1 && isB2 &&
                        (hcptIdx == (hcptLL.length - 1)) &&
                        (cost < 0.065);

                    boolean greenException = bothGreen &&
                        ((hgsInter >= hgsLL[hcptIdx]) ||
                        (cost < 0.05 && hcptInter >= 0.6
                        && hgsInter >= 0.3));

                    /*
                    System.out.format("m %d %d (%d,%d) (%d,%d) hsvd=%.3f ptInter=%.3f "
                        + " gradInter=%.3f n=%d,%d wh=%b,%b->%b "
                        + " bl=%b,%b->%b grE=%b->%b\n"
                        + "\n    hsv1=%.3f,%.3f,%.3f"
                        + "\n    hsv2=%.3f,%.3f,%.3f\n",
                        idx1, idx2, xyCen1[0], xyCen1[1], xyCen2[0], xyCen2[1],
                        cost, hcptInter, hgsInter,
                        set1.size(), set2.size(),
                        isW1, isW2, simW,
                        isB1, isB2, simB, bothGreen, greenException,
                        hsv1.getAvgH(), hsv1.getAvgS(),
                        hsv1.getAvgV(),
                        hsv2.getAvgH(), hsv2.getAvgS(),
                        hsv2.getAvgV()
                    );*/
                    //System.out.println("gs hists=\n    " +
                    //    Arrays.toString(hgs1H) + "\n    " +
                    //    Arrays.toString(hgs2H));

                    if ((hcptInter < hcptLL[hcptIdx] ||
                        hgsInter < hgsLL[hcptIdx])
                        && !simW && !simB && !greenException
                        ) {
                        //System.out.format("     %.3f %.3f\n",
                        //    hcptLL[hcptIdx], hgsLL[hcptIdx]
                        //);
                        continue;
                    }

                    hcptInter = 1.f - hcptInter;
                    hgsInter = 1.f - hgsInter;

                    cost *= cost;
                    cost += (hcptInter * hcptInter + hgsInter * hgsInter);
                    cost = (float)Math.sqrt(cost/3.f);

                    if (cost < minCost) {
                        minCost = cost;
                        minCostIdx2 = idx2;
                    }
                }

                if (minCostIdx2 > -1) {

                    // merging contents of idx1 into minCostIdx2
                    nMerged++;

                    //System.out.println("    merging " + minCostIdx2);

                    clrs.get(minCostIdx2).add(set1, clrImg);
                    clrs.remove(idx1);

                    TIntIterator iter = mapOfSets.get(idx1).iterator();
                    while (iter.hasNext()) {
                        int pixIdx = iter.next();
                        pointIndexMap.put(pixIdx, minCostIdx2);
                    }

                    // update the adjacency maps to point to minCostIdx2
                    for (int idx2 : idxs2) {
                        VeryLongBitString nbrsBS3 = adjMap.get(idx2);
                        nbrsBS3.clearBit(idx1);
                        nbrsBS3.setBit(minCostIdx2);
                    }

                    VeryLongBitString union = nbrsBS.or(adjMap.get(minCostIdx2));
                    union.clearBit(minCostIdx2);
                    union.clearBit(idx1);
                    adjMap.put(minCostIdx2, union);
                    adjMap.remove(idx1);

                    mapOfSets.get(minCostIdx2).addAll(mapOfSets.get(idx1));
                    //TODO: consider updating borders by clipping out intersection
                    mapOfBorders.get(minCostIdx2).addAll(mapOfBorders.get(idx1));
                    mapOfSets.remove(idx1);
                    mapOfBorders.remove(idx1);

                    assert(mapOfSets.get(idx1) == null);
                    assert(mapOfBorders.get(idx1) == null);
                }
            }
            if (debug) {
                //System.out.println("nMerged=" + nMerged + " nIter=" +
                //    nIter);
            }
            nIter++;

            if (nMerged == 0 && hcptIdx < (hcptLL.length - 1)) {
                hcptIdx++;
                nIter = 0;
                nMerged = 1;
            }
        } while (nIter < nIterMax && nMerged > 0);

        // reset the instance vars
        edgeList.clear();
        labeledSets.clear();;

        TIntObjectIterator<TIntSet> iter2 = mapOfSets.iterator();

        for (int i = 0; i < mapOfSets.size(); ++i) {

            iter2.advance();

            int idx = iter2.key();

            TIntSet set1 = iter2.value();

            if (set1.isEmpty()) {
                continue;
            }

            labeledSets.add(set1);

            TIntSet embedded = new TIntHashSet();
            TIntSet outerBorder = new TIntHashSet();
            finder2.extractBorder2(set1, embedded, outerBorder, w);
            edgeList.add(outerBorder);
        }

        if (debug) {
            Image imgCp = clrImg.copyImage();
            ImageIOHelper.addAlternatingColorCurvesToImage3(edgeList,
                imgCp, 0);
            MiscDebug.writeImage(imgCp, "_" + ts + "_MERGED_");
        }
    }

    /**
     * experimental method for a specific use case.  if this method is retained
     * and used, it could be made more efficient.
     * @return
     */
    public List<Region> _extractSensitiveGS0() {

        GreyscaleImage invImg = gsImg.copyImage();
        // Invert the pixel values
        for (int i = 0; i < invImg.getNPixels(); ++i) {
            int v = invImg.getValue(i);
            v = ~v;
            if (v < 0) {
                v += 256;
            }
            invImg.setValue(i, v);
        }

        List<List<Region>> gsRegions = extractMSERRegions(invImg,
            Threshold.LESS_SENSITIVE);

        List<Region> regions2 = new ArrayList<Region>();

        for (int type = 1; type < 2; ++type) {
            List<Region> list = gsRegions.get(type);
            for (int i = (list.size() - 1); i > -1; --i) {
                Region r = list.get(i);
                if ((type == 1) && r.getVariation() > 2.) {
                    list.remove(i);
                } else if ((type == 0)
                    //&& r.getVariation() == 0.0) {
                    && r.getVariation() < 2.) {
                    list.remove(i);
                } else {
                    regions2.add(r);
                }
            }
        }

        return regions2;
    }

    public List<TIntList> getEmbeddedGS0Levels() {
        if (sensitiveGS0 == null) {
            sensitiveGS0 = _extractSensitiveGS0();
        }

        return getEmbeddedLevels(sensitiveGS0);
    }

    private List<TIntList> getEmbeddedLevels(List<Region> rs) {

        //NOTE: if end up using this heavily, should consider a layered level
        //   set structure from the start

        // looking for level sets of concentric ellipses
        int[] sizes = new int[rs.size()];
        for (int i = 0; i < rs.size(); ++i) {
            sizes[i] = rs.get(i).area_;
        }
        QuickSort.sortBy1stArg(sizes, rs);

        //_debugOrigRegions(rs, "_rs_");

        List<EllipseHelper> ehs = new ArrayList<EllipseHelper>();
        List<GroupPixelHSV> hsvs = new ArrayList<GroupPixelHSV>();

        int[] xyCen = new int[2];
        int w = clrImg.getWidth();
        int h = clrImg.getHeight();

        QuadTree<Integer, Integer> centroidQT = new QuadTree<Integer, Integer>();

        // store the region centroids in quadtree
        for (int i = 0; i < rs.size(); ++i) {
            Region r = rs.get(i);
            r.calculateXYCentroid(xyCen, w, h);
            double[] coeffs = r.calcParamTransCoeff();
            EllipseHelper eh = new EllipseHelper(xyCen[0], xyCen[1], coeffs);
            ehs.add(eh);

            GroupPixelHSV gHSV = new GroupPixelHSV();
            gHSV.calculateColors(r.getAcc(), clrImg);
            hsvs.add(gHSV);

            centroidQT.insert(xyCen[0], xyCen[1], Integer.valueOf(i));
        }

        TIntList skip = new TIntArrayList();

        /*
        concentric regions due to illumintion gradient tend to have
        similar eccentricities, orientation, and direction away from
        the largest region center.
        They also tend to have 3 or 4 of the ellipse semi major and
        minor endpoints within the largest region's bounds,
        that is the direction of the illumination gradient is in one
        direction.
        */

        List<TIntList> concentric = new ArrayList<TIntList>();

        // visit from largest region to smallest
        for (int i = (rs.size() - 1); i > -1; --i) {
            if (skip.contains(i)) { continue; }
            EllipseHelper eh = ehs.get(i);

            if (hsvs.get(i).getStdDevV() > 0.075) {
                continue;
            }

            int[] minMaxXY = eh.getMinMaxXY();

            Interval<Integer> intX = new Interval<Integer>(minMaxXY[0], minMaxXY[1]);
            Interval<Integer> intY = new Interval<Integer>(minMaxXY[2], minMaxXY[3]);
            Interval2D<Integer> rect = new Interval2D<Integer>(intX, intY);

            List<Integer> indexes = centroidQT.query2D(rect);
            if (indexes == null || indexes.size() < 2) {
                continue;
            }

            float orientation = (float)(
                eh.getOrientation() * 180./Math.PI);
            double ecc = eh.getEccentricity();

            System.out.format(
                "OUTER (%d,%d) or=%.3f ecc=%.3f  stdvV=%.3f n=%d\n",
                eh.getXYCenter()[0], eh.getXYCenter()[1],
                (float)eh.getOrientation(),
                (float)eh.getEccentricity(), hsvs.get(i).getStdDevV(),
                rs.get(i).accX.size()
            );

            TIntList idx2s = new TIntArrayList();
            TDoubleList dirs = new TDoubleArrayList();
            for (int j = 0; j < indexes.size(); ++j) {
                int idx2 = indexes.get(j).intValue();
                if (idx2 == i || skip.contains(idx2)) { continue; }

                if (hsvs.get(idx2).getStdDevV() > 0.075) {
                    continue;
                }

                EllipseHelper eh2 = ehs.get(idx2);
                float orientation2 = (float)(
                    eh2.getOrientation() * 180./Math.PI);

                float diffOr = AngleUtil.getAngleDifference(orientation,
                    orientation2);

                if (Math.abs(diffOr) > 5.75) {
                    continue;
                }
                if (Math.abs(ecc - eh2.getEccentricity()) > 0.07) {
                    continue;
                }

                PairInt[] ep = eh2.getSemiAxesEndoints();
                int nTrue = 0;
                for (int k = 0; k < 4; ++k) {
                    if (eh.isWithin(ep[k].getX(), ep[k].getY())) {
                        nTrue++;
                    }
                }

                boolean t0 = eh.isWithin(ep[0].getX(), ep[0].getY());
                boolean t1 = eh.isWithin(ep[1].getX(), ep[1].getY());
                boolean t2 = eh.isWithin(ep[2].getX(), ep[2].getY());
                boolean t3 = eh.isWithin(ep[3].getX(), ep[3].getY());


                double atan2 = Math.atan2(
                    (double)(eh2.getXYCenter()[1] -  eh.getXYCenter()[1]),
                    (double)(eh2.getXYCenter()[0] -  eh.getXYCenter()[0]))
                    * 180./Math.PI;
                if (atan2 < 0) {
                    atan2 += 360;
                }

                //System.out.println("     " + Arrays.toString(ep));

                System.out.format(
                    "     (%d,%d) or=%.3f ecc=%.3f atan2=%.3f "
                        + "isWithin=%b,%b,%b,%b  stdvV=%.3f n=%d\n",
                    eh2.getXYCenter()[0], eh2.getXYCenter()[1],
                    (float)eh2.getOrientation(),
                    (float)eh2.getEccentricity(), (float)atan2,
                    t0, t1, t2, t3, hsvs.get(idx2).getStdDevV(),
                    rs.get(idx2).accX.size()
                );

                if (nTrue < 2) {
                    continue;
                }

   //cannot include if there is alot of variation in the level set

                // can't tell which directions are consistent w/ a gradient
                // yet, so collect all
                idx2s.add(idx2);
                dirs.add(atan2);
            }

            if (dirs.size() < 2) {
                continue;
            }

            // gradient will be apparent as 2 or more regions having
            // similar atan2 (within approx 10 degrees)
            // key = bin number where bin size is 10 degrees, value = number
            //    of items in dirs which have that direction key bin
            float binSz = 10.f;
            TIntIntMap dirFreq = new TIntIntHashMap();
            for (int k = 0; k < dirs.size(); ++k) {
                double direction = dirs.get(k);
                int bin = (int)(direction/binSz);
                if (dirFreq.containsKey(bin)) {
                    int c = dirFreq.get(bin);
                    dirFreq.put(bin, c + 1);
                } else {
                    dirFreq.put(bin, 1);
                }
            }
            int maxFreq = Integer.MIN_VALUE;
            int maxFreqBin = -1;
            TIntIntIterator iter = dirFreq.iterator();
            for (int k = 0; k < dirFreq.size(); ++k) {
                iter.advance();
                int c = iter.value();
                if (c > maxFreq) {
                    maxFreq = c;
                    maxFreqBin = iter.key();
                }
            }
            if (maxFreq < 2) {
                continue;
            }
            // combine those ideally within 10 degrees of each other,
            //  so roughly 7.07 from bin center
            float dirCenter = ((float)maxFreqBin)*binSz + (binSz/2.f);
            double dist = Math.sqrt(2) * binSz/2.;

            System.out.println("   adding");

            //List<TIntList> concentric = new ArrayList<TIntList>();
            TIntList concList = new TIntArrayList();
            concentric.add(concList);
            concList.add(i);
            skip.add(i);
            for (int k = 0; k < dirs.size(); ++k) {
                double direction = dirs.get(k);
                double diff = Math.abs(direction - dirCenter);
                if (diff > dist) {
                    continue;
                }
                int idx2 = idx2s.get(k);
                concList.add(idx2);
                skip.add(idx2);
            }

            /*{//DEBUG
                List<Region> tmp = new ArrayList<Region>();
                for (int k = 0; k < concList.size(); ++k) {
                    tmp.add(rs.get(concList.get(k)));
                }
                _debugOrigRegions(tmp, "_rs_" + concentric.size());
            }*/
        }

        return concentric;
    }

    public void _debugOrigRegions(int idx, String lbl) {
        List<Region> list = origGsPtRegions.get(idx);
        lbl = lbl + "_" + idx + "_";
        _debugOrigRegions(list, lbl);
    }

    public void _debugOrigRegions(List<Region> list, String lbl) {
        int[] xyCen = new int[2];
        Image imCp;
        System.out.println("printing " + list.size());
        for (int j = 0; j < list.size(); ++j) {
            imCp = ptImg.copyToColorGreyscale();
            Region r = list.get(j);
            for (int i = 0; i < r.accX.size(); ++i) {
                int x = r.accX.get(i);
                int y = r.accY.get(i);
                ImageIOHelper.addPointToImage(x, y, imCp, 0, 0, 255, 0);
            }
            r.drawEllipse(imCp, 0, 255, 0, 0);
            r.calculateXYCentroid(xyCen, imCp.getWidth(), imCp.getHeight());
            ImageIOHelper.addPointToImage(xyCen[0], xyCen[1], imCp,
                1, 255, 0, 0);
            //System.out.println(type + " xy=" + xyCen[0] + "," + xyCen[1]
            //    + " variation=" + r.getVariation());
            MiscDebug.writeImage(imCp, "_" + ts + "_orig_regions_" + lbl
                + "_" + "_" + j + "_");
        }
    }

    public List<TIntSet> getEdges() {

        //INITIALIZED, REGIONS_EXTRACTED, MERGED, EDGES_EXTRACTED
        if (!state.equals(STATE.EDGES_EXTRACTED)) {
            throw new IllegalStateException("must use one of the extraction"
                + " methods first");
        }

        return edgeList;
    }

    /**
     * get labeled sets of points separated by the edges.
     * Note that the method returns null unless mergeAndExtractEdges was
     * already used.
     * @return
     */
    public List<TIntSet> getLabeledSets() {

        //INITIALIZED, REGIONS_EXTRACTED, MERGED, EDGES_EXTRACTED
        if (!state.equals(STATE.EDGES_EXTRACTED)) {
            throw new IllegalStateException("must use one of the extraction"
                + " methods first");
        }

        return labeledSets;
    }

    /**
     * @return the clrImg
     */
    public ImageExt getClrImg() {
        return clrImg;
    }

    /**
     * @return the gsImg
     */
    public GreyscaleImage getGsImg() {
        return gsImg;
    }

    /**
     * @return the ptImg
     */
    public GreyscaleImage getPtImg() {
        return ptImg;
    }

    /**
     * @return the regions
     */
    public List<Region> getRegions() {
        return regions;
    }

    /**
     * @return the origGsPtRegions
     */
    public List<List<Region>> getOrigGsPtRegions() {
        return origGsPtRegions;
    }

    /**
     * @return the useLowerContrastLimits
     */
    public boolean isUsingLowerContrastLimits() {
        return useLowerContrastLimits;
    }

    private GreyscaleImage createSobelScores() {

        int w = gsImg.getWidth();
        int h = gsImg.getHeight();

        ImageProcessor imageProcessor = new ImageProcessor();

        /*
        float[] sobelScores = imageProcessor.createSobelColorScores(
            gsImg, ptImg, 20);

        GreyscaleImage scaled = MiscMath.rescaleAndCreateImage(sobelScores,
            w, h);
        */

        CannyEdgeColorAdaptive canny2 = new CannyEdgeColorAdaptive();
        canny2.overrideToNotUseLineThinner();
        canny2.applyFilter(clrImg);
        GreyscaleImage scaled = canny2.getFilterProducts().getGradientXY();
        scaled.multiply(255/scaled.max());
        //MiscDebug.writeImage(scaled, debugLabel
        //    + "_lc_edges_");

        /*
        GreyscaleImage sobelLCH =
            imageProcessor.createSobelLCCombined(img);
        */

        this.cannyEdges = scaled;
                
        // smearing values over a 3 pixel window to avoid the potential
        //   1 pixel displacement of an edge from the level set boundaries
        SummedAreaTable sumTable = new SummedAreaTable();

        GreyscaleImage imgM = sumTable.createAbsoluteSummedAreaTable(scaled);
        imgM = sumTable.applyMeanOfWindowFromSummedAreaTable(imgM, 3);

        return imgM;
    }

    private TIntSet combineBoundaries() {

        if (sobelScores == null) {
            sobelScores = createSobelScores();
            //MiscDebug.writeImage(sobelScores,
            //    "_" + ts + "_canny_blurred_");
            //MiscDebug.writeImage(cannyEdges, "_" + ts + "_canny_");
        }

        /*
        using the level set region boundaries as contours to complete
        the edges from the canny edges of L and C from LCH.

        -- first applying a filter that is the fraction of region
           boundary points which have a canny edge pixel to filter the
           mser regions list.
        -- then finding the contiguous sets of all matched region points
        and the sets of all unmatched points from those same regions
        (where matched means that the region pixel coordinates
        have canny edge pixel values > 0).
        -- then finding the paths in the unmatched points that complete
        gaps between the matched segments.
        */

        // below this removes:
        float mLimit = 0.4001f;

        // maximum number of contiguous pixels that are not already canny
        //   edge points
        int maxGapSize = 16;

        PerimeterFinder2 finder = new PerimeterFinder2();

        TIntSet allEdgePoints = new TIntHashSet();
        TIntSet unmatchedPoints = new TIntHashSet();
        TIntObjectMap<TIntSet> unmatchedRMap = new TIntObjectHashMap<TIntSet>();
        
        TIntSet rmvdImgBorders = new TIntHashSet();
        
        double tt0 = System.currentTimeMillis();
        
        for (int rListIdx = 0; rListIdx < regions.size(); ++rListIdx) {
            Region r = regions.get(rListIdx);
            TIntSet points = r.getAcc(clrImg.getWidth());
            TIntSet embedded = new TIntHashSet();
            TIntSet outerBorder = new TIntHashSet();
            finder.extractBorder2(points, embedded, outerBorder, clrImg.getWidth());

            TIntSet border2 = removeImageBorder(outerBorder, rmvdImgBorders,
                clrImg.getWidth(), clrImg.getHeight());

            TIntSet matched = new TIntHashSet();
            TIntSet unmatched = new TIntHashSet();

            double[] scoreAndMatch = calcAvgScore(border2, sobelScores,
                matched, unmatched, clrImg.getWidth());

            double matchFraction = scoreAndMatch[1]/(double)border2.size();

            /*
            System.out.format(" rIdx=%d score=%.3f "
                + " border2.n=%d nmf=%.3f  m.n=%d um.n=%d\n",
                rListIdx, (float)scoreAndMatch[0],
                border2.size(),
                (float)matchFraction, matched.size(), unmatched.size());
            */

            boolean doNotAdd = (matchFraction < mLimit) ||
                ((int)scoreAndMatch[1] == 0);

//Image tmpImg = sobelScores.copyToColorGreyscale();
//ImageIOHelper.addCurveToImage(border2, tmpImg, 0, 255, 0, 0);
//MiscDebug.writeImage(tmpImg, "_r_" + rListIdx);

            if (doNotAdd) {
                continue;
            }

//Image tmpImg = sobelScores.copyToColorGreyscale();
//ImageIOHelper.addCurveToImage(border2, tmpImg, 0, 255, 0, 0);
//MiscDebug.writeImage(tmpImg, "_r_" + rListIdx);

            /*
            to find the points in border2 which are not in the canny edges, but
            are the missing contour points of an object:

            -- separating the border2 points which are matched to canny edge
                 points from those which do not have a canny edge pixel
                 (done in steps above)
            -- for the matched points, putting them all in allPoints and
                putting unmatched into unmatchedSet.
            -- when this block is complete,
               -- will find the contiguous matched segments.
               -- will find the unmatchedSet points which are adjacent
                  to one of the matched.
                  these are the endpoints in the unmatched which may be path
                  endpoints between matched segments or they might not be.
               -- will search through unmatched points from the starting
                  endpoints to the other endpoints
                  and keep paths that are shorter than
                  maxGapSize
            */

            allEdgePoints.addAll(matched);
            unmatchedPoints.addAll(unmatched);
            
            unmatchedRMap.put(rListIdx, unmatched);
            
//Image tmpImg = sobelScores.copyToColorGreyscale();
//ImageIOHelper.addCurveToImage(unmatched, tmpImg, 0, 255, 0, 0);
//MiscDebug.writeImage(tmpImg, "_u_" + rListIdx);

        }
        
        double tt1 = System.currentTimeMillis();

        if (debug) {
            System.out.println(((tt1 - tt0)/1000.) + 
                " sec for region filter. " + unmatchedRMap.size() +
                " regions added");
        }
        
        /*
        if (debug) {
            Image tmpImg = clrImg.copyToGreyscale2().copyToColorGreyscale();
            ImageIOHelper.addCurveToImage(allEdgePoints, tmpImg, 0, 255, 0, 0);
            MiscDebug.writeImage(tmpImg, "_" + ts + "_matched_");
            tmpImg = clrImg.copyToGreyscale2().copyToColorGreyscale();
            ImageIOHelper.addCurveToImage(unmatchedPoints, tmpImg, 0, 255, 0, 0);
            MiscDebug.writeImage(tmpImg, "_" + ts + "_unmatched_");
            //tmpImg = clrImg.copyToGreyscale2().copyToColorGreyscale();
            //ImageIOHelper.addCurveToImage(rmvdImgBorders, tmpImg, 0, 255, 0, 0);
            //MiscDebug.writeImage(tmpImg, "_" + ts + "_rmvdBounds_");
        }
        */
        
        //make contiguous connected segments of matched set.
        ConnectedPointsFinder finder2 = new ConnectedPointsFinder(
            clrImg.getWidth(), clrImg.getHeight());
        finder2.setMinimumNumberInCluster(1);
        finder2.findConnectedPointGroups(allEdgePoints);

        final int n2 = finder2.getNumberOfGroups();
        
        System.out.println("number of matched segments=" + n2);

        // key = matched point, value = finder2 index of point
        TIntIntMap mpIdxMap = finder2.createPointIndexMap();
        
        //   sections in unmatchedRMap are single width segments
        //   .  when a segment is adjacent to 2 matched segments (from finder2),
        //      that unmatched segment can be added.
        //      note that multiple regions may have adjacent arcs in
        //      an unmatched area, so to thin to the best path,
        //      could either take the best
        //      unmatched segment that connects those 2 matched segments
        //      for that area, or could let the "assignTheUnassigned"
        //      kmeans-style method handle this thinning.
        //      choosing the later since it is needed for the matched
        //      points already.
        //      NOTE: the unmatched dist between endpoints AND the number
        //            of pixels in the unmatched segment still need to remain under
        //            the maxGapSize, helping to avoid adding the region
        //            portions which are not edges.
       
        TIntObjectMap<List<TIntSet>> unmatchedRSegments = 
            extractontiguousSegments(unmatchedRMap);

        TIntObjectIterator<List<TIntSet>> iter10 = 
            unmatchedRSegments.iterator();
        for (int i = 0; i < unmatchedRSegments.size(); ++i) {
            
            iter10.advance();
            
            int rIdx = iter10.key();
            List<TIntSet> segmentList = iter10.value();
            for (TIntSet segment : segmentList) {
        
                if (segment.size() > maxGapSize) {
                    continue;
                }
                
                // key = pixIdx, value=edge segment indexes that key
                //       pixel is adjacent to
                // NOTE, that for pixels adjacent to an image border pixel,
                //       that is one off from image border, a fake edge
                //       is added to the values and those fake edges are
                //       negative numbers starting at -1
                TIntObjectMap<TIntList> umEPIdxMap
                    = findEndpoints(segment, mpIdxMap, 
                    clrImg.getWidth(), clrImg.getHeight());
                
                assert(!umEPIdxMap.isEmpty());
        
                // for any 2 entries in umEPIdxMap
                //   if the set difference operation results in 2 edge
                //      indexes, can add this segment to allEdgePoints.
                
                TIntSet posEdgeIdxs = new TIntHashSet(3);
                boolean connectsToBorder = false;
               
                TIntObjectIterator<TIntList> iter11 = umEPIdxMap.iterator();
                for (int ii = 0; ii < umEPIdxMap.size(); ++ii) {
                    iter11.advance();
                    TIntList eIdxs = iter11.value();
                    for (int jj = 0; jj < eIdxs.size(); ++jj) {
                        int eIdx = eIdxs.get(jj);
                        if (eIdx < 0) {
                            connectsToBorder = true;
                        } else {
                            posEdgeIdxs.add(eIdx);
                        }
                    }
                }
                
                if (posEdgeIdxs.size() >= 1) {
                    allEdgePoints.addAll(segment);
                }
            }
        }
        
        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;
        
        int w = clrImg.getWidth();
        int h = clrImg.getHeight();
                
        //add back any points in rmvdImgBorders adjacent to allEdgePoints
        TIntSet addRmvd = new TIntHashSet();
        TIntIterator iter = rmvdImgBorders.iterator();
        while (iter.hasNext()) {
            int pixIdx = iter.next();
            int y = pixIdx / w;
            int x = pixIdx - (y * w);

            //System.out.println(" imgB xy=" + x + ", " + y);

            for (int k = 0; k < dxs.length; ++k) {
                int x2 = x + dxs[k];
                int y2 = y + dys[k];
                if (x2 < 0 || y2 < 0 || (x2 >= w) || (y2 >= w)) {
                    continue;
                }
                int pixIdx2 = (y2 * w) + x2;
                if (allEdgePoints.contains(pixIdx2)) {
                    addRmvd.add(pixIdx);
                    break;
                }
            }
        }
        allEdgePoints.addAll(addRmvd);

        if (debug) {
            Image tmp = clrImg.copyToGreyscale2().copyToColorGreyscale();
            ImageIOHelper.addCurveToImage(allEdgePoints, tmp, 0, 255, 0, 0);
            MiscDebug.writeImage(tmp, "_" + ts + "_boundaries0_");
        }

        return allEdgePoints;
    }

    /**
     *
     *
     * @param allPoints
     * @return
     */
    private void closeGapsOf1(TIntSet allPoints) {

        int w = gsImg.getWidth();
        int h = gsImg.getHeight();

        GreyscaleImage img2 = gsImg.copyImage();
        TIntIterator iter = allPoints.iterator();
        while (iter.hasNext()) {
            int pixIdx = iter.next();
            img2.setValue(pixIdx, 1);
        }

        ImageSegmentation imageSegmentation = new ImageSegmentation();
        TIntSet outputAddedGaps = new TIntHashSet();
        img2 = imageSegmentation.fillInCompleteGapsOf1(
            img2, outputAddedGaps, 1);

        // restore gap where the gap is completely surrounded
        imageSegmentation.restoreGapsOf1WhereSurrounded(
            img2, outputAddedGaps, 1);

        TIntSet thinned = new TIntHashSet();
        for (int i = 0; i < img2.getNPixels(); ++i) {
            // img2 edges have pixel value=1
            if (img2.getValue(i) > 0) {
                thinned.add(i);
            }
        }
    }

    /**
     * using the given edges as definitions of separation between
     * contiguous pixels in the image (== labels),
     * perform kmeans on the edge pixels to add them back to
     * the labels which are adjacent and closest in color.
     *
     * NOTE: has the side effect of populating this.edgeList and 
     * this.labeledSets
     */
    private void thinTheBoundaries(int minGroupSize) {

        if (this.labeledSets == null || this.edgeList == null) {
            throw new IllegalStateException("instance vars are null: "
                + " labeledSets, edgeLists");
        }
        
        PerimeterFinder2 finder2 = new PerimeterFinder2();

        int w = clrImg.getWidth();
        int h = clrImg.getHeight();

        float[][] hsvs = new float[labeledSets.size()][];

        int[] labels = new int[clrImg.getNPixels()];
        Arrays.fill(labels, -1);

        for (int label = 0; label < labeledSets.size(); ++label) {

            TIntSet set = labeledSets.get(label);

            TIntIterator iter2 = set.iterator();
            while (iter2.hasNext()) {
                int pixIdx = iter2.next();
                labels[pixIdx] = label;
            }

            GroupPixelHSV2 hsv = new GroupPixelHSV2();
            hsv.calculateColors(set, clrImg);
            hsvs[label] = new float[]{hsv.getAvgH(), hsv.getAvgS(),
                hsv.getAvgV()};
        }

        TIntSet unassignedSet = new TIntHashSet();
        for (int i = 0; i < clrImg.getNPixels(); ++i) {
            if (labels[i] == -1) {
                unassignedSet.add(i);
            }
        }

        if (debug) {
            Image tmp = clrImg.copyImage();
            ImageIOHelper.addAlternatingColorPointSetsToImage2(
                labeledSets, 0, 0, 0, tmp);
            MiscDebug.writeImage(tmp, "_" + ts + "_before_assigned_");
        }

        assignTheUnassigned(labeledSets, labels, hsvs, unassignedSet);

        if (debug) {
            Image tmp = clrImg.copyImage();
            ImageIOHelper.addAlternatingColorPointSetsToImage2(
                labeledSets, 0, 0, 0, tmp);
            MiscDebug.writeImage(tmp, "_" + ts + "_reassigned0_");
        }

        // make successive passes through to re-assign the smallest sets,
        //    pixel by pixel.  the gradual merging of smaller sets helps
        //    preserve some of the low SNR object boundaries
        //    such as the gingerbread man in the test image android_statues_01
        //    downsampled to size near 256 pixels per dimension.
        int[] mszs = new int[]{minGroupSize, 4, 6, 12, 18, 24};

        for (int msz : mszs) {

            // make a pass through results to find any sets that do not have
            //  embedded points and re-submit those if any
            List<TIntSet> contigousSets2 = new ArrayList<TIntSet>();
            List<TIntSet> edgeLists2 = new ArrayList<TIntSet>();
            TIntSet reassign = new TIntHashSet();
            for (int i = 0; i < labeledSets.size(); ++i) {
                TIntSet set = labeledSets.get(i);
                TIntSet embedded = new TIntHashSet();
                TIntSet outerBorder = new TIntHashSet();
                finder2.extractBorder2(set, embedded, outerBorder, w);
                if (set.size() - outerBorder.size() >= msz) {
                    contigousSets2.add(set);
                    edgeLists2.add(outerBorder);
                } else {
                    reassign.addAll(set);
                }
            }
            if (reassign.isEmpty()) {
                labeledSets = contigousSets2;
                edgeList = edgeLists2;
            } else {
                labels = new int[clrImg.getNPixels()];
                Arrays.fill(labels, -1);
                hsvs = new float[contigousSets2.size()][];
                for (int label = 0; label < contigousSets2.size(); ++label) {
                    TIntSet set = contigousSets2.get(label);
                    TIntIterator iter2 = set.iterator();
                    while (iter2.hasNext()) {
                        int pixIdx = iter2.next();
                        labels[pixIdx] = label;
                    }
                    GroupPixelHSV2 hsv = new GroupPixelHSV2();
                    hsv.calculateColors(set, clrImg);
                    hsvs[label] = new float[]{hsv.getAvgH(), hsv.getAvgS(),
                        hsv.getAvgV()};
                }

                assignTheUnassigned(contigousSets2, labels, hsvs, reassign);

                labeledSets.clear();
                edgeList.clear();

                for (int i = 0; i < contigousSets2.size(); ++i) {
                    TIntSet set = contigousSets2.get(i);
                    TIntSet embedded = new TIntHashSet();
                    TIntSet outerBorder = new TIntHashSet();
                    finder2.extractBorder2(set, embedded, outerBorder,
                        clrImg.getWidth());

                    labeledSets.add(set);
                    edgeList.add(outerBorder);
                }
            }
        }

        if (debug) {
            Image tmp = clrImg.copyImage();
            ImageIOHelper.addAlternatingColorPointSetsToImage2(
                labeledSets, 0, 0, 0, tmp);
            MiscDebug.writeImage(tmp, "_" + ts + "thinned0_");
        }
    }

    /**
     * NOTE, this method actually finds the connected groups .geq.
     * minGroupSize and then extracts the boundaries from them.
     * If the set minus the boundaries is not at least minGroupSize,
     * then the set is excluded.
     * So the method uses a minGroupSize that does not include the boundary
     * pixels.
     * @param edgePoints
     * @param minGroupSize
     */
    private void populateEdgeLists(TIntSet edgePixIdxs) {

        // find clusters (contiguous pixels of value 0) between edges
        labeledSets = extractContiguousBetweenEdges(edgePixIdxs);

        edgeList = new ArrayList<TIntSet>();

        PerimeterFinder2 finder2 = new PerimeterFinder2();

        for (int i = 0; i < labeledSets.size(); ++i) {
            TIntSet set = labeledSets.get(i);
            TIntSet embedded = new TIntHashSet();
            TIntSet outerBorder = new TIntHashSet();
            finder2.extractBorder2(set, embedded, outerBorder, gsImg.getWidth());

            edgeList.add(outerBorder);
        }

        assert(labeledSets.size() == edgeList.size());

        System.out.println(labeledSets.size() + " labeled sets");
    }

    private void assignTheUnassigned(List<TIntSet> contiguous,
        int[] labels, float[][] hsvs, TIntSet unassignedSet) {

        //TODO: could refactor this to use a fibonacci heap
        //   with the key being one to maximize the number
        //   of assigned neighbors.
        //   currently, the code uses a one time sort by decreasing
        //   number of assigned neighbors.
        //   fiboncci heap decreasekey is O(1) so updates after each assignment
        //   are fast

        // key = pixel index, value =
        TIntObjectMap<TIntSet> unassignedMap
            = new TIntObjectHashMap<TIntSet>();

        int w = clrImg.getWidth();
        int h = clrImg.getHeight();
        int n = clrImg.getNPixels();

        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;

        TIntIterator iter = unassignedSet.iterator();
        while (iter.hasNext()) {
            int pixIdx = iter.next();
            addNeighborLabelsForPoint(labels, unassignedMap, pixIdx,
                dxs, dys);
        }

        ArrayDeque<Integer> queue0 = populateByNumberOfNeighbors(
            unassignedMap);

        ArrayDeque<Integer> queue1 = new ArrayDeque<Integer>();

        int nIter = 0;

        float[] lab2 = null;

        TIntSet visited = new TIntHashSet();

        while (!queue0.isEmpty() || !queue1.isEmpty()) {

            int pixIdx;
            if (!queue1.isEmpty()) {
                pixIdx = queue1.poll().intValue();
            } else {
                pixIdx = queue0.poll().intValue();
            }

            if (visited.contains(pixIdx)) {
                continue;
            }
            visited.add(pixIdx);

            int y1 = pixIdx/w;
            int x1 = pixIdx - (y1 * w);

            TIntSet adjLabels;
            if (nIter == 0) {
                adjLabels = unassignedMap.get(pixIdx);
                assert (adjLabels != null);
            } else {
                adjLabels = new TIntHashSet();
                addNeighborLabelsForPoint(labels, adjLabels,
                    pixIdx, dxs, dys);
            }

            double minD = Double.MAX_VALUE;
            int minLabel2 = -1;

            float[] lab = new float[3];
            lab[0] = clrImg.getHue(x1, y1);
            lab[1] = clrImg.getSaturation(x1, y1);
            lab[2] = clrImg.getBrightness(x1, y1);

            TIntIterator iter2 = adjLabels.iterator();
            while (iter2.hasNext()) {
                int label2 = iter2.next();
                lab2 = hsvs[label2];
                double diffSum = 0;
                for (int i = 0; i < 3; ++i) {
                    float diff = lab[i] - lab2[i];
                    diffSum += (diff * diff);
                }

                if (diffSum < minD) {
                    minD = diffSum;
                    minLabel2 = label2;
                }
            }

            labels[pixIdx] = minLabel2;

            unassignedMap.remove(pixIdx);

            for (int m = 0; m < dxs.length; ++m) {
                int x2 = x1 + dxs[m];
                int y2 = y1 + dys[m];
                if (x2 < 0 || y2 < 0 || (x2 > (w - 1)) || (y2 > (h - 1))) {
                    continue;
                }
                int pixIdx2 = clrImg.getInternalIndex(x2, y2);
                if (labels[pixIdx2] == -1) {
                    queue1.add(pixIdx2);
                    //assert (!visited.contains(pixIdx2));
                }
            }
            nIter++;
        }

        assert(unassignedMap.isEmpty());

        /*for (int pixIdx = 0; pixIdx < labels.length; ++pixIdx) {
            int label = labels[pixIdx];
            contiguous.get(label).add(pixIdx);
        }*/
        iter = unassignedSet.iterator();
        while (iter.hasNext()) {
            int pixIdx = iter.next();
            int label = labels[pixIdx];
            contiguous.get(label).add(pixIdx);
        }
    }

    private void addNeighborLabelsForPoint(int[] labels,
        TIntObjectMap<TIntSet> unassignedMap,
        int pixIdx, int[] dxs, int[] dys) {

        int w = clrImg.getWidth();
        int h = clrImg.getHeight();

        TIntSet adjLabels = unassignedMap.get(pixIdx);
        if (adjLabels == null) {
            adjLabels = new TIntHashSet();
            unassignedMap.put(pixIdx, adjLabels);
        }

        addNeighborLabelsForPoint(labels, adjLabels, pixIdx, dxs, dys);
    }

    private void addNeighborLabelsForPoint(int[] labels, TIntSet adjLabels,
        int pixIdx, int[] dxs, int[] dys) {

        int w = clrImg.getWidth();
        int h = clrImg.getHeight();

        int j = pixIdx/w;
        int i = pixIdx - (j * w);

        for (int m = 0; m < dxs.length; ++m) {
            int x2 = i + dxs[m];
            int y2 = j + dys[m];
            if (x2 < 0 || y2 < 0 || (x2 > (w - 1)) || (y2 > (h - 1))) {
                continue;
            }
            int pixIdx2 = clrImg.getInternalIndex(x2, y2);
            if (labels[pixIdx2] > -1) {
                adjLabels.add(labels[pixIdx2]);
            }
        }
    }

    private ArrayDeque<Integer> populateByNumberOfNeighbors(
        TIntObjectMap<TIntSet> unassignedMap) {

        int n = unassignedMap.size();

        int[] pixIdxs = new int[n];
        int[] nN = new int[n];

        TIntObjectIterator<TIntSet> iter = unassignedMap.iterator();

        for (int count = 0; count < n; ++count) {
            iter.advance();
            int pixIdx = iter.key();
            pixIdxs[count] = pixIdx;
            nN[count] = iter.value().size();
        }

        QuickSort.sortBy1stArg(nN, pixIdxs);

        ArrayDeque<Integer> queue = new ArrayDeque<Integer>();

        for (int i = (n - 1); i > -1; --i) {

            int nP = nN[i];
            if (nP == 0) {
                break;
            }

            queue.add(Integer.valueOf(pixIdxs[i]));
        }

        return queue;
    }

    private GreyscaleImage copyAndShift(GreyscaleImage polarImage, int valueShift) {

        GreyscaleImage shiftedImg = polarImage.createWithDimensions();
        for (int i = 0; i < polarImage.getNPixels(); ++i) {
            int v = polarImage.getValue(i);
            v += valueShift;
            if (v < 0) {
                v += 255;
            } else if (v > 255) {
                v -= 255;
            }
            shiftedImg.setValue(i, v);
        }

        return shiftedImg;
    }

    private TIntSet removeImageBorder(TIntSet pixIdxs, TIntSet outputRmvd,
        int width, int height) {

        TIntSet set = new TIntHashSet(pixIdxs);

        TIntSet rm = new TIntHashSet();
        TIntIterator iter = set.iterator();
        while (iter.hasNext()) {
            int pixIdx = iter.next();
            int y = pixIdx/width;
            int x = pixIdx - (y * width);
            if (x == 0 || y == 0 || x == (width - 1) || (y == (height - 1))) {
                rm.add(pixIdx);
                outputRmvd.add(pixIdx);
            }
        }

        set.removeAll(rm);

        return set;
    }

    // calc avg score and count the points w/ v>0
    private double[] calcAvgScore(TIntSet pixIdxs, GreyscaleImage sobelScores,
        TIntSet outputMatched, TIntSet outputUnmatched, int imgWidth) {

        double sum = 0;
        TIntIterator iter = pixIdxs.iterator();
        while (iter.hasNext()) {
            int pixIdx = iter.next();
            int v = sobelScores.getValue(pixIdx);
            sum += v;
            if (v > 0) {
                outputMatched.add(pixIdx);
            } else {
                outputUnmatched.add(pixIdx);
            }
        }

        sum /= (double) pixIdxs.size();

        return new double[]{sum, outputMatched.size()};
    }

    private boolean similarExists(Region rToCheck, List<RegionGeometry> listSRG) {

        int w = ptImg.getWidth();
        int h = ptImg.getHeight();

        RegionGeometry rg0 = Canonicalizer.calculateEllipseParams(rToCheck, w, h);
        float orientation0 = (float)(
            rg0.orientation * 180./Math.PI);

        for (int i1 = 0; i1 < listSRG.size(); ++i1) {
            RegionGeometry rg1 = listSRG.get(i1);
            float orientation1 = (float)(
                rg1.orientation * 180./Math.PI);

            float diffOr = AngleUtil.getAngleDifference(orientation0,
                orientation1);

            if ((Math.abs(rg0.xC - rg1.xC) > 5) ||
                (Math.abs(rg0.yC - rg1.yC) > 5) ||
                (Math.abs(rg0.minor - rg1.minor) > 5) ||
                (Math.abs(rg0.major - rg1.major) > 5) ||
                (Math.abs(diffOr) > 5)
                ) {
                continue;
            }
            return true;
        }
        return false;
    }

    private int calcAvg(GreyscaleImage img, TIntSet pixIdxs) {

        double sum = 0;
        TIntIterator iter = pixIdxs.iterator();
        while (iter.hasNext()) {
            int pixIdx = iter.next();
            sum += img.getValue(pixIdx);
        }
        sum /= (double)pixIdxs.size();

        return (int)Math.round(sum);
    }

    private int[] getRegionHistogram(HCPT hcpt0, TIntSet pIdxs) {
        int[] h0 = new int[hcpt0.getNumberOfBins()];
        TIntIterator iter = pIdxs.iterator();
        while (iter.hasNext()) {
            int pixIdx = iter.next();
            int y = pixIdx/gsImg.getWidth();
            int x = pixIdx - (y * gsImg.getWidth());
            hcpt0.extractFeature(x, y, h0);
        }
        return h0;
    }

    private int[] getRegionHistogram(HGS hgs0, TIntSet pIdxs) {
        int[] h0 = new int[hgs0.getNumberOfBins()];
        TIntIterator iter = pIdxs.iterator();
        while (iter.hasNext()) {
            int pixIdx = iter.next();
            int y = pixIdx/gsImg.getWidth();
            int x = pixIdx - (y * gsImg.getWidth());
            hgs0.extractFeature2(x, y, h0);
        }
        return h0;
    }

    private GroupPixelHSV2 getColors(TIntObjectMap<GroupPixelHSV2> clrs,
        TIntSet pIdxs, int idx) {

        GroupPixelHSV2 hsv = clrs.get(idx);
        if (hsv != null) {
            return hsv;
        }

        hsv = new GroupPixelHSV2();
        hsv.calculateColors(pIdxs, clrImg);

        return hsv;
    }

    private TIntIntMap createPointIndexMap(TIntObjectMap<TIntSet> mapOfSets) {

        TIntIntMap pointIndexMap = new TIntIntHashMap();

        TIntObjectIterator<TIntSet> iter = mapOfSets.iterator();

        for (int i = 0; i < mapOfSets.size(); ++i) {

            iter.advance();

            int label = iter.key();
            TIntSet set = iter.value();

            TIntIterator iter2 = set.iterator();
            while (iter2.hasNext()) {
                int pixIdx = iter2.next();
                pointIndexMap.put(pixIdx, label);
            }
        }

        return pointIndexMap;
    }

    private boolean isBlack(GroupPixelHSV2 hsv) {

        if (hsv.getAvgV() >= 0.2) {
            //System.out.println("brightness=" + hsv.getAvgV());
            return false;
        }

        return hsv.isGrey(12);
    }

    private boolean isWhite(GroupPixelHSV2 hsv) {

        if (hsv.getAvgV() < 0.625) {
            //System.out.println("brightness=" + hsv.getAvgV());
            return false;
        }

        return hsv.isGrey(12);
    }

    // the values are sorted
    private TIntObjectMap<TIntList> findEndpoints(
        TIntSet unmatchedPoints,
        TIntIntMap matchedPointsIdxMap, int width, int height) {

        int bIdx = -1;
        
        TIntObjectMap<TIntList> umEPIdxMap =
            new TIntObjectHashMap<TIntList>();

        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;

        TIntIterator iter = unmatchedPoints.iterator();
        while (iter.hasNext()) {
            int pixIdx = iter.next();
            int y = pixIdx/width;
            int x = pixIdx - (y * width);

            // add fake point if this is adjacent to an image border pixel
            if ((x == 1) || (y == 1) || (x == (width - 2)) || (y == (height - 2))) {
                TIntList idxs = umEPIdxMap.get(pixIdx);
                if (idxs == null) {
                    idxs = new TIntArrayList(3);
                    umEPIdxMap.put(pixIdx, idxs);
                }
                idxs.add(bIdx);
                --bIdx;
            }
            
            for (int k = 0; k < dxs.length; ++k) {
                int x2 = x + dxs[k];
                int y2 = y + dys[k];
                if (x2 < 0 || y2 < 0 || (x2 >= width) || (y2 >= height)) {
                    continue;
                }
                int pixIdx2 = (y2 * width) + x2;
                if (matchedPointsIdxMap.containsKey(pixIdx2)) {
                    TIntList values = umEPIdxMap.get(pixIdx);
                    if (values == null) {
                        values = new TIntArrayList(3);
                        umEPIdxMap.put(pixIdx, values);
                    }
                    values.add(matchedPointsIdxMap.get(pixIdx2));
                }
            }
        }
        
        return umEPIdxMap;
    }

    private double distance(int pixIdx0, int pixIdx1, int width) {

        int y0 = pixIdx0/width;
        int x0 = pixIdx0 - (y0 * width);

        int y1 = pixIdx1/width;
        int x1 = pixIdx1 - (y1 * width);

        int diffX = x0 - x1;
        int diffY = y0 - y1;

        double dist = Math.sqrt(diffX * diffX + diffY * diffY);

        return dist;
    }
    
    private TIntObjectMap<List<TIntSet>> extractontiguousSegments(
        TIntObjectMap<TIntSet> unmatchedRMap) {
       
        int n = unmatchedRMap.size();
        
        System.out.println("number of unmatched regions=" + n);
        
        TIntObjectMap<List<TIntSet>> output = new 
            TIntObjectHashMap<List<TIntSet>>();
        
        TIntObjectIterator<TIntSet> iter = unmatchedRMap.iterator();
        for (int i = 0; i < unmatchedRMap.size(); ++i) {
            iter.advance();
            
            int rIdx = iter.key();
            
            TIntSet idxs = iter.value();
            
            ConnectedPointsFinder finder3 = new ConnectedPointsFinder(
                clrImg.getWidth(), clrImg.getHeight());
            finder3.setMinimumNumberInCluster(1);
            finder3.setToUse8Neighbors();
            finder3.findConnectedPointGroups(idxs);
            
            List<TIntSet> list = new ArrayList<TIntSet>();
            output.put(rIdx, list);
            
            int n3 = finder3.getNumberOfGroups();
            for (int ii = 0; ii < n3; ++ii) {
                
                TIntSet segment = finder3.getXY(ii);
                
                list.add(segment);
            }
        }
       
        return output;
    }
}
