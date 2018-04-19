package algorithms.imageProcessing.features.mser;

import algorithms.CountingSort;
import algorithms.QuickSort;
import algorithms.bipartite.MinHeapForRT2012;
import algorithms.compGeometry.PerimeterFinder2;
import algorithms.imageProcessing.CannyEdgeColorAdaptive;
import algorithms.imageProcessing.ColorHistogram;
import algorithms.connected.ConnectedPointsFinder;
import algorithms.imageProcessing.CIEChromaticity;
import algorithms.imageProcessing.EdgeFilterProducts;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.GroupPixelHSV;
import algorithms.imageProcessing.GroupPixelHSV2;
import algorithms.imageProcessing.GroupPixelLUVWideRangeLightness;
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
import algorithms.imageProcessing.features.HOGs;
import algorithms.imageProcessing.features.PatchUtil;
import algorithms.imageProcessing.features.mser.Canonicalizer.RegionGeometry;
import algorithms.imageProcessing.features.mser.MSER.Threshold;
import algorithms.imageProcessing.matching.CMODE;
import algorithms.imageProcessing.util.AngleUtil;
import algorithms.misc.Misc;
import algorithms.misc.MiscDebug;
import algorithms.util.PairInt;
import algorithms.util.PixelHelper;
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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
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
 * Also note that their are fixed parameters tailored for the input images
 * scaled to 256 X 256, especially for the use of merging
 * (see int[] mszs = ...).
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
    
    private EdgeFilterProducts edgeProducts = null;

    private long ts = 0;

    private STATE state = null;

    public MSEREdges(Image img) {

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

        mergeRegions();

        long ts1 = System.currentTimeMillis();

        System.out.format("%.3f sec for merge\n", (((float)ts1 - ts0)/1000.f));
        
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
                negImg = ptImg.copyImage();
                for (int i = 0; i < ptImg.getNPixels(); ++i) {
                    int v = ~ptImg.getValue(i);
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
            
            //System.out.println(nonEdgePixels.size() + " non-edge pixels");
        }
        
        finder.findConnectedPointGroups(nonEdgePixels);
        
        List<TIntSet> output = new ArrayList<TIntSet>();
        
        int ns = 0;
        for (int i = 0; i < finder.getNumberOfGroups(); ++i) {
            TIntSet group = finder.getXY(i);
            output.add(group);
            ns += group.size();
        }
        
        assert(ns == nonEdgePixels.size());
        
        return output;
    }

    private void extractBoundaries() {

        //INITIALIZED, REGIONS_EXTRACTED, MERGED, EDGES_EXTRACTED
        if (state.equals(STATE.EDGES_EXTRACTED)) {
            throw new IllegalStateException("can only perform extraction of "
                + "edges once");
        }
        
        long ts0 = System.currentTimeMillis();

        TIntSet boundaries = combineBoundaries();

        long ts1 = System.currentTimeMillis();

        // populate this.edgeList and this.labeledSets
        populateEdgeLists(boundaries);
        
        if (debug) {
            Image tmp = clrImg.copyImage();
            ImageIOHelper.addAlternatingColorPointSetsToImage2(
                labeledSets, 0, 0, 0, tmp);
            MiscDebug.writeImage(tmp, "_" + ts + "_labeled_boundaries_");
        }
        
        long ts1_2 = System.currentTimeMillis();
        
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

        Image im = clrImg.copyToGreyscale2().copyToColorGreyscale();

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
    
    private void mergeRegions() {
        //mergeRegionsWithHGSRGBHSV();
        mergeRegionsWithHOGsHCPTHGS();
    }

    /**
     * moderate merging of the labeled regions is performed
     * to remove noisy edges.  Note that the color filters
     * may need to be revised with more testing.
     *
     */
    private void mergeRegionsWithHGSRGBHSV() {

        /*
        The use of gradient histograms as texture in this method was inspired by
        HOGs in general and by Alpert, Galun, Basri et al. (2007).
        */

        //INITIALIZED, REGIONS_EXTRACTED, MERGED, EDGES_EXTRACTED
        if (!state.equals(STATE.EDGES_EXTRACTED)) {
            throw new IllegalStateException("error in algorithm.  expecting"
                + " edges were extracted.");
        }

        assert(labeledSets != null);
        assert(edgeList != null);
        assert(labeledSets.size() == edgeList.size());

        if (sobelScores == null) {
            sobelScores = createSobelScores();
        }

        HGS hgs = new HGS(cannyEdges, 1, 6, 12);
        
        TIntObjectMap<GroupPixelHSV2> clrs = new TIntObjectHashMap<GroupPixelHSV2>();

        int w = gsImg.getWidth();
        int h = gsImg.getHeight();

        for (int label = 0; label < labeledSets.size(); ++label) {
            TIntSet set = new TIntHashSet(labeledSets.get(label));
            set.removeAll(edgeList.get(label));
            GroupPixelHSV2 hsv = new GroupPixelHSV2();
            hsv.calculateColors(set, clrImg);
            clrs.put(label, hsv);
        }
      
        ImageProcessor imageProcessor = new ImageProcessor();

        int[] sizes = new int[labeledSets.size()];
        int[] indexes = new int[sizes.length];
        
        TIntIntMap pointIndexMap = new TIntIntHashMap();
        for (int i = 0; i < labeledSets.size(); ++i) {
            TIntSet pixIdxs = labeledSets.get(i);
            TIntIterator iter = pixIdxs.iterator();
            while (iter.hasNext()) {
                int pixIdx = iter.next();
                pointIndexMap.put(pixIdx, i);
            }
            sizes[i] = pixIdxs.size();
            indexes[i] = i;
        }
        
        QuickSort.sortBy1stArg(sizes, indexes);
        
        // -- making adjacency map using edgeList ---
        TIntObjectMap<VeryLongBitString> adjMap = imageProcessor
            .createAdjacencyMap(pointIndexMap, edgeList, w, h);

        ColorHistogram clrHist = new ColorHistogram();
        
        MiscellaneousCurveHelper ch = new MiscellaneousCurveHelper();

        //DEBUG: needed just for logging
        TIntObjectMap<PairInt> centroidsMap = new TIntObjectHashMap<PairInt>();
        for (int i = 0; i < labeledSets.size(); ++i) {
            TIntSet pixIdxs = labeledSets.get(i);
            int[] xyCen = ch.calculateRoundedXYCentroids(pixIdxs, w);
            centroidsMap.put(i, new PairInt(xyCen[0], xyCen[1]));
        }
        
        for (int i = 0; i < labeledSets.size(); ++i) {
            
            int label = indexes[i];
           
            VeryLongBitString adjLabels = adjMap.get(label);
            if (adjLabels == null || adjLabels.getNSetBits() == 0) {
                continue;
            }
            
            PairInt xy1 = centroidsMap.get(label);
            GroupPixelHSV2 clr1 = clrs.get(label);
            
            // subtr edges from sets
            TIntSet set1 = new TIntHashSet(labeledSets.get(label));
            set1.removeAll(edgeList.get(label));
            
            int[] hgs1H = getRegionHistogram(hgs, set1);

            // 0=other, 1=black, 2=white
            int clrMode1 = isBlack(clr1) ? 1 :
                (isWhite(clr1) ? 2 : 0);
            
            int[] adjBits = adjLabels.getSetBits();
            for (int label2 : adjBits) {
                
                if (!clrs.containsKey(label2)) {
                    continue;
                }
                
                PairInt xy2 = centroidsMap.get(label2);
                GroupPixelHSV2 clr2 = clrs.get(label2);
                
                // subtr edges from sets
                TIntSet set2 = new TIntHashSet(labeledSets.get(label2));
                set2.removeAll(edgeList.get(label2));
                
                int[] hgs2H = getRegionHistogram(hgs, set2);
            
                float hsvDiff = clr1.calculateDifference(clr2);
                float rgbDiff = clr1.calculateRGBDifference(clr2);
                float hgsInter = hgs.intersection(hgs1H, hgs2H);

                // 0=other, 1=black, 2=white
                int clrMode2 = isBlack(clr2) ? 1 :
                    (isWhite(clr2) ? 2 : 0);
            
                //System.out.format(
                //    "%s to %s : clr1M=%d clr2M=%d "
                //    + " hsvDiff=%.3f rgbDiff=%.3f hgsI=%.3f\n",
                //    xy1.toString(), xy2.toString(), 
                //    clrMode1, clrMode2, hsvDiff, rgbDiff, hgsInter);
            
                if (clrMode1 == 0 && clrMode2 == 0) {
                    if (hsvDiff > 0.09 || hgsInter < 0.89 || rgbDiff > 0.1) {
                        continue;
                    }
                } else if (clrMode1 == 1 && clrMode2 == 1) {
                    if (hgsInter < 0.8 || rgbDiff > 0.08) {
                        continue;
                    }
                } else if (clrMode1 == 2 && clrMode2 == 2) {
                    //whiteish
                    if (hgsInter < 0.8 || rgbDiff > 0.08) {
                        continue;
                    }
                } else {
                    continue;
                }
                
                //System.out.println("  merging");
                
                // -- merge the adjacent into the current label ---
                
                clr1.add(clr2);
                clrs.remove(label2);

                // not updating the pointIndexMap because not using it here
                
                /*
                label :  label2, label3, label4, ...
                         becomes
                         label
                */
                
                adjLabels.clearBit(label2);
                
                VeryLongBitString adjLabels2 = adjMap.get(label2);
                adjLabels2.clearBit(label);
                int[] setBits2 = adjLabels2.getSetBits();
                for (int label2Adj : setBits2) {
                    VeryLongBitString adjLabel2Adj = adjMap.get(label2Adj);
                    if (adjLabel2Adj == null) {
                        adjLabels2.clearBit(label2Adj);
                        continue;
                    }
                    adjLabel2Adj.clearBit(label2);
                    adjLabel2Adj.setBit(label);
                }
                
                VeryLongBitString union = adjLabels.or(adjLabels2);
                adjMap.put(label, union);
                adjLabels = union;
                adjMap.remove(label2);
                
                labeledSets.get(label).addAll(labeledSets.get(label2));
                edgeList.get(label).addAll(edgeList.get(label2));
                labeledSets.get(label2).clear();
                edgeList.get(label2).clear();
                                
                // subtr edges from sets
                set1 = new TIntHashSet(labeledSets.get(label));
                set1.removeAll(edgeList.get(label));
            
                clrHist.add2To1(hgs1H, hgs2H);
            
                // debug
                int[] xyCen = ch.calculateRoundedXYCentroids(set1, w);
                xy1 = new PairInt(xyCen[0], xyCen[1]);
                centroidsMap.put(label, xy1);
                
            }
        }
        
        // remove empty labeledSets
        for (int i = (labeledSets.size() - 1); i > -1; --i) {
            TIntSet set = labeledSets.get(i);
            if (set.isEmpty()) {
                labeledSets.remove(i);
            }
        }
        
        PerimeterFinder2 finder2 = new PerimeterFinder2();
        
        edgeList.clear();
        
        // redo the edgeList since the merged sets have embedded edges
        for (int i = 0; i < labeledSets.size(); ++i) {
            TIntSet set = labeledSets.get(i);
            TIntSet embedded = new TIntHashSet();
            TIntSet outerBorder = new TIntHashSet();
            finder2.extractBorder2(set, embedded, outerBorder, w);
            edgeList.add(outerBorder);
        }
        assert(edgeList.size() == labeledSets.size());
        
        if (debug) {
            Image imgCp = clrImg.copyImage();
            ImageIOHelper.addAlternatingColorCurvesToImage3(edgeList,
                imgCp, 0);
            MiscDebug.writeImage(imgCp, "_" + ts + "_MERGED_");
        }
    }

    private void mergeRegionsWithHOGsHCPTHGS() {
       
        //INITIALIZED, REGIONS_EXTRACTED, MERGED, EDGES_EXTRACTED
        if (!state.equals(STATE.EDGES_EXTRACTED)) {
            throw new IllegalStateException("error in algorithm.  expecting"
                + " edges were extracted.");
        }

        assert(labeledSets != null);
        assert(edgeList != null);
        assert(labeledSets.size() == edgeList.size());

        if (sobelScores == null) {
            // these are actually CannyEdgeColorAdaptive
            sobelScores = createSobelScores();
        }

        ImageProcessor imageProcessor = new ImageProcessor();
        
        //DEBUG.  replace ptImg with this if keep this method
        GreyscaleImage luvImg = 
            imageProcessor.createCIELUVTheta_WideRangeLightness(clrImg, 255);
        
        int nCellsPerDim = 1;
        int nPixPerCellDim = 6;
        int nBins = 12;
        
        HOGs hogs = new HOGs(gsImg, nCellsPerDim, nPixPerCellDim, nBins);
        HCPT hcpt = new HCPT(luvImg, nCellsPerDim, nPixPerCellDim, nBins);
        HGS hgs = new HGS(gsImg, nCellsPerDim, nPixPerCellDim, nBins);
        
        // integrating the histograms over the regions to be compared for merging
        TIntObjectMap<List<PatchUtil>> clrs = new TIntObjectHashMap<List<PatchUtil>>();

        int w = gsImg.getWidth();
        int h = gsImg.getHeight();

        for (int label = 0; label < labeledSets.size(); ++label) {
            TIntSet set = new TIntHashSet(labeledSets.get(label));
            //set.removeAll(edgeList.get(label));
            
            List<PatchUtil> list = new ArrayList<PatchUtil>();
            PatchUtil p = new PatchUtil(w, h, nBins);
            p.add(set, hogs);
            list.add(p);
            p = new PatchUtil(w, h, nBins);
            p.add(set, hcpt);
            list.add(p);
            p = new PatchUtil(w, h, nBins);
            p.add(set, hgs);
            list.add(p);
            
            clrs.put(label, list);
        }
      
        int[] sizes = new int[labeledSets.size()];
        int[] indexes = new int[sizes.length];
        
        TIntIntMap pointIndexMap = new TIntIntHashMap();
        for (int i = 0; i < labeledSets.size(); ++i) {
            TIntSet pixIdxs = labeledSets.get(i);
            TIntIterator iter = pixIdxs.iterator();
            while (iter.hasNext()) {
                int pixIdx = iter.next();
                pointIndexMap.put(pixIdx, i);
            }
            sizes[i] = pixIdxs.size();
            indexes[i] = i;
        }
        
        QuickSort.sortBy1stArg(sizes, indexes);
        
        // -- making adjacency map using edgeList ---
        TIntObjectMap<VeryLongBitString> adjMap = imageProcessor
            .createAdjacencyMap(pointIndexMap, edgeList, w, h);
        
        MiscellaneousCurveHelper ch = new MiscellaneousCurveHelper();

        //DEBUG: needed just for logging
        TIntObjectMap<PairInt> centroidsMap = new TIntObjectHashMap<PairInt>();
        for (int i = 0; i < labeledSets.size(); ++i) {
            TIntSet pixIdxs = labeledSets.get(i);
            int[] xyCen = ch.calculateRoundedXYCentroids(pixIdxs, w);
            centroidsMap.put(i, new PairInt(xyCen[0], xyCen[1]));
        }
        
        int mCount = 0;
        
        for (int i = 0; i < labeledSets.size(); ++i) {
            
            int label = indexes[i];
           
            VeryLongBitString adjLabels = adjMap.get(label);
            if (adjLabels == null || adjLabels.getNSetBits() == 0) {
                continue;
            }
            
            // subtr edges from sets
            TIntSet set1 = new TIntHashSet(labeledSets.get(label));
            set1.removeAll(edgeList.get(label));
            
            PairInt xy1 = centroidsMap.get(label);
       
            List<PatchUtil> pList1 = clrs.get(label);
            
            Set<PairInt> points1 = pList1.get(0).getPixelSet();
            CMODE cmodeLUV1 = CMODE.determinePolarThetaMode(luvImg, points1);
            CMODE cmode1 = CMODE.determineColorMode(clrImg, points1);

            /*
            GroupPixelHSV2 clr1 = clrs.get(label);
            int[] hgs1H = getRegionHistogram(hgs, set1);
            // 0=other, 1=black, 2=white
            int clrMode1 = isBlack(clr1) ? 1 : (isWhite(clr1) ? 2 : 0);
            */
            
            int[] adjBits = adjLabels.getSetBits();
            for (int label2 : adjBits) {
                
                if (!clrs.containsKey(label2)) {
                    continue;
                }
                
                PairInt xy2 = centroidsMap.get(label2);
                
                // subtr edges from sets
                TIntSet set2 = new TIntHashSet(labeledSets.get(label2));
                //set2.removeAll(edgeList.get(label2));
                
                List<PatchUtil> pList2 = clrs.get(label2);
                
                Set<PairInt> points2 = pList2.get(0).getPixelSet();
                CMODE cmodeLUV2 = CMODE.determinePolarThetaMode(luvImg, points2);
                CMODE cmode2 = CMODE.determineColorMode(clrImg, points2);
            
                //DEBUG
                float inter0 = (float)pList1.get(0).intersection(pList2.get(0));
                float inter1 = (float)pList1.get(1).intersection(pList2.get(1));
                float inter2 = (float)pList1.get(2).intersection(pList2.get(2));
                System.out.format(
                    "(%d,%d) %s %s : (%d,%d) %s %s intersection=%.3f,%.3f,%.3f", 
                    xy1.getX(), xy1.getY(), cmodeLUV1.name(), cmode1.name(),
                    xy2.getX(), xy2.getY(), cmodeLUV2.name(), cmode2.name(),
                    inter0, inter1, inter2
                );
                float limit = 0.6f;
                if (inter0 < limit || inter1 < limit || inter2 < limit) {
                    System.out.format("\n");
                    continue;
                }
                
                /*
                boolean skip = false;
                for (int j = 0; j < pList1.size(); ++j) {
                    if (pList1.get(j).intersection(pList2.get(j)) < 0.5) {
                        // different so do not merge
                        skip = true;
                        break;
                    }
                }
                if (skip) {
                    continue;
                }*/
                
                /*
                GroupPixelHSV2 clr2 = clrs.get(label2);
                
                int[] hgs2H = getRegionHistogram(hgs, set2);
                float hsvDiff = clr1.calculateDifference(clr2);
                float rgbDiff = clr1.calculateRGBDifference(clr2);
                float hgsInter = hgs.intersection(hgs1H, hgs2H);
                // 0=other, 1=black, 2=white
                int clrMode2 = isBlack(clr2) ? 1 : (isWhite(clr2) ? 2 : 0);
            
                if (clrMode1 == 0 && clrMode2 == 0) {
                    if (hsvDiff > 0.09 || hgsInter < 0.89 || rgbDiff > 0.1) {
                        continue;
                    }
                } else if (clrMode1 == 1 && clrMode2 == 1) {
                    if (hgsInter < 0.8 || rgbDiff > 0.08) {
                        continue;
                    }
                } else if (clrMode1 == 2 && clrMode2 == 2) {
                    //whiteish
                    if (hgsInter < 0.8 || rgbDiff > 0.08) {
                        continue;
                    }
                } else {
                    continue;
                }
                clr1.add(clr2);
                clrs.remove(label2);
                */
                
                //System.out.println("  merging");
                
                // -- merge the adjacent into the current label ---
                
                pList1.get(0).add(set2, hogs);
                pList1.get(1).add(set2, hcpt);
                pList1.get(2).add(set2, hgs);
                
                clrs.remove(label2);
                
                // not updating the pointIndexMap because not using it here
                
                /*
                label :  label2, label3, label4, ...
                         becomes
                         label
                */
                
                adjLabels.clearBit(label2);
                
                VeryLongBitString adjLabels2 = adjMap.get(label2);
                adjLabels2.clearBit(label);
                int[] setBits2 = adjLabels2.getSetBits();
                for (int label2Adj : setBits2) {
                    VeryLongBitString adjLabel2Adj = adjMap.get(label2Adj);
                    if (adjLabel2Adj == null) {
                        adjLabels2.clearBit(label2Adj);
                        continue;
                    }
                    adjLabel2Adj.clearBit(label2);
                    adjLabel2Adj.setBit(label);
                }
                
                VeryLongBitString union = adjLabels.or(adjLabels2);
                adjMap.put(label, union);
                adjLabels = union;
                adjMap.remove(label2);
                
                labeledSets.get(label).addAll(labeledSets.get(label2));
                edgeList.get(label).addAll(edgeList.get(label2));
                labeledSets.get(label2).clear();
                edgeList.get(label2).clear();
                                
                // subtr edges from sets
                set1 = new TIntHashSet(labeledSets.get(label));
                //set1.removeAll(edgeList.get(label));
            
                //clrHist.add2To1(hgs1H, hgs2H);
            
                // debug
                if (false && debug) {
                    Image tmp = clrImg.copyImage();
                    ImageIOHelper.addAlternatingColorPointSetsToImage2(
                        labeledSets, 0, 0, 0, tmp);
                    String t = Integer.toString(mCount);
                    while (t.length() < 5) {
                        t = "0" + t;
                    }
                    t = t + "_" + xy1.toString() + "_" + xy2.toString();
                    MiscDebug.writeImage(tmp, "_" + ts + "_merging_" + t);
                }
                
                int[] xyCen = ch.calculateRoundedXYCentroids(set1, w);
                xy1 = new PairInt(xyCen[0], xyCen[1]);
                centroidsMap.put(label, xy1);
                
                System.out.format(" ==> (%d,%d)\n", xy1.getX(), xy1.getY());
                mCount++;
            }
            
        }
        
        // remove empty labeledSets
        for (int i = (labeledSets.size() - 1); i > -1; --i) {
            TIntSet set = labeledSets.get(i);
            if (set.isEmpty()) {
                labeledSets.remove(i);
            }
        }
        
        PerimeterFinder2 finder2 = new PerimeterFinder2();
        
        edgeList.clear();
        
        // redo the edgeList since the merged sets have embedded edges
        for (int i = 0; i < labeledSets.size(); ++i) {
            TIntSet set = labeledSets.get(i);
            TIntSet embedded = new TIntHashSet();
            TIntSet outerBorder = new TIntHashSet();
            finder2.extractBorder2(set, embedded, outerBorder, w);
            edgeList.add(outerBorder);
        }
        assert(edgeList.size() == labeledSets.size());
        
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
    
    /**
     * experimental method for a specific use case.  if this method is retained
     * and used, it could be made more efficient.
     * @return
     */
    public List<Region> _extractSensitivePT0() {

        int delta = 2;
        double minArea = 0.001;//.0001
        double maxArea = 0.99;//0.1;
        double maxVariation = 0.9;//0.5;
        double minDiversity = 0.75;//0.1;//0.5
        
        int[] a = MSER.readIntoArray(ptImg);
        
        int w = ptImg.getWidth();
        int h = ptImg.getHeight();

        List<Region> pt0Regions = new ArrayList<Region>();
        MSER mser8 = new MSER(delta, minArea, maxArea, maxVariation, 
            minDiversity, true);
        mser8.operator(a, w, h, pt0Regions);
        
        a = MSER.readIntoArray(ptImgShifted);
        List<Region> pt0ShiftedRegions = new ArrayList<Region>();
        mser8 = new MSER(delta, minArea, maxArea, maxVariation, 
            minDiversity, true);
        mser8.operator(a, w, h, pt0ShiftedRegions);
        
        //_debugOrigRegions(ptShiftedRegions.get(0), "_shifted_0");
        //_debugOrigRegions(ptShiftedRegions.get(1), "_shifted_1");
            
        boolean hadWrapAroundArtifacts = false;
        for (int i = (pt0Regions.size() - 1); i > -1; --i) {
            
            Region r = pt0Regions.get(i);

            if (r.getVariation() == 0.0) {
                
                pt0Regions.remove(i);
                
            } else {
                // check for this Region being an artifact of wrap around
                if (r.level_ < 25) {
                    pt0Regions.remove(i);
                    hadWrapAroundArtifacts = true;
                } else {
                    int avgLevel = calcAvg(ptImg, r.getAcc(w));
                    
                    //System.out.format(" %d add x,y=%d,%d level=%d  avgLevel=%d\n",
                    //    type, (int)(r.moments_[0]/r.area_),
                    //    (int)(r.moments_[1]/r.area_), r.level_,
                    //    avgLevel);
                    if (avgLevel > 240) {//230?
                        pt0Regions.remove(i);
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
            for (int i = (pt0ShiftedRegions.size() - 1); i > -1; --i) {
                Region r = pt0ShiftedRegions.get(i);

                if (r.getVariation() == 0.0) {
                    pt0ShiftedRegions.remove(i);
                } else {
                    //TODO: consider checking whether this already exists in
                    //   the list
                    if (Math.abs(r.level_ - 60) < 20) {
                        r.level_ += 60;
                        pt0Regions.add(r);
                        //System.out.format("  add shifted x,y=%d,%d level=%d\n",
                        // (int)(r.moments_[0]/r.area_),
                        // (int)(r.moments_[1]/r.area_), r.level_);
                    }
                    //TODO: consider adding other regions in list2
                    //   as long as level is > 25 and < 230
                }
            }
        }
        
        return pt0Regions;
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

        /*
        float[] sobelScores = imageProcessor.createSobelColorScores(
            gsImg, ptImg, 20);

        GreyscaleImage scaled = MiscMath.rescaleAndCreateImage(sobelScores,
            w, h);
        */

        CannyEdgeColorAdaptive canny2 = new CannyEdgeColorAdaptive();
        canny2.overrideToNotUseLineThinner();
        canny2.applyFilter(clrImg);
        this.edgeProducts = canny2.getFilterProducts();
        GreyscaleImage scaled = edgeProducts.getGradientXY();
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
        List<TIntSet> unmatchedRMap = new ArrayList<TIntSet>();
        
        TIntSet rmvdImgBorders = new TIntHashSet();
        
        double tt0 = System.currentTimeMillis();
        
        List<TIntSet> unusedRegionBounds = new ArrayList<TIntSet>();
        
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
                unusedRegionBounds.add(border2);
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
            
            unmatchedRMap.add(unmatched);
            
//Image tmpImg = sobelScores.copyToColorGreyscale();
//ImageIOHelper.addCurveToImage(unmatched, tmpImg, 0, 255, 0, 0);
//MiscDebug.writeImage(tmpImg, "_u_" + rListIdx);

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
        
        double tt1 = System.currentTimeMillis();

        addUnmatchedToEdgePoints(unmatchedPoints, unmatchedRMap,
            allEdgePoints, rmvdImgBorders, maxGapSize);
      
        long tt2 = System.currentTimeMillis();
        
        if (debug) {
            System.out.println(((tt1 - tt0)/1000.) + 
                " sec for region filter " +
                ((tt2 - tt1)/1000.) + " sec for add unmatched "
                );
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

        float[][] luvws = new float[labeledSets.size()][];

        int[] labels = new int[clrImg.getNPixels()];
        Arrays.fill(labels, -1);

        for (int label = 0; label < labeledSets.size(); ++label) {

            TIntSet set = labeledSets.get(label);

            TIntIterator iter2 = set.iterator();
            while (iter2.hasNext()) {
                int pixIdx = iter2.next();
                labels[pixIdx] = label;
            }

            GroupPixelLUVWideRangeLightness luv = new GroupPixelLUVWideRangeLightness();
            luv.calculateColors(set, clrImg);
            luvws[label] = luv.getAvgLUV();
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

        _assignTheUnassigned(labeledSets, labels, luvws, unassignedSet);

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
        //int[] mszs = new int[]{minGroupSize, 4, 6, 12, 18, 24};
        int[] mszs = new int[]{minGroupSize, 4, 6, 12, 18};

        float[][] hsvs;
        
        for (int msz : mszs) {

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
            
            
            if (debug) {
                Image tmp = clrImg.copyImage();
                ImageIOHelper.addAlternatingColorPointSetsToImage2(
                    labeledSets, 0, 0, 0, tmp);
                MiscDebug.writeImage(tmp, "_" + ts + "_dbg_" + msz);
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

        //if (debug) {
        //    System.out.println("populateEdgeLists: nE=" + edgePixIdxs.size() 
        //        + " pixels");
        //}
        
        // find clusters (contiguous pixels of value 0) between edges
        labeledSets = extractContiguousBetweenEdges(edgePixIdxs);

        edgeList = new ArrayList<TIntSet>();

        PerimeterFinder2 finder2 = new PerimeterFinder2();

        int ne = 0;
        for (int i = 0; i < labeledSets.size(); ++i) {
            TIntSet set = labeledSets.get(i);
            TIntSet embedded = new TIntHashSet();
            TIntSet outerBorder = new TIntHashSet();
            finder2.extractBorder2(set, embedded, outerBorder, gsImg.getWidth());

            edgeList.add(outerBorder);
        
            ne += outerBorder.size();
        }
        
        if (debug) {
            System.out.println(ne + " pixels in edgeList");
        }

        assert(labeledSets.size() == edgeList.size());

        System.out.println(labeledSets.size() + " labeled sets");
    }

    private void assignTheUnassigned(List<TIntSet> contiguous,
        int[] labels, float[][] hsvs, TIntSet unassignedSet) {

        int w = clrImg.getWidth();
        int h = clrImg.getHeight();
        int n = clrImg.getNPixels();

        int maxValue = Math.max(w, h);
        int nBits = 1 + (int)Math.ceil(Math.log(maxValue)/Math.log(2));
        if (nBits > 31) {
            nBits = 31;
        }
        
        //System.out.println("assign " + unassignedSet.size() + " out of "
        //    + n + " pixels");
        
        // key = pixel index of unassigned, value = adj pixels that are assigned
        TIntObjectMap<TIntSet> adjAssignedMap = new TIntObjectHashMap<TIntSet>();
        
        // key = pixel index of unassigned, value = adj pixels that are unassigned
        TIntObjectMap<TIntSet> adjUnassignedMap = new TIntObjectHashMap<TIntSet>();

        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;

        TIntIterator iter = unassignedSet.iterator();
        while (iter.hasNext()) {
            int pixIdx = iter.next();
            addNeighborLabelsForPoint(labels, adjAssignedMap, 
                adjUnassignedMap, pixIdx, dxs, dys);
        }
        
        // using a min heap whose priority is to set the nodes
        //    which have the largest number of assigned neighbors.
        //    nAssigned=8 -> key=8-nAssigned = 0.
        MinHeapForRT2012 heap = new MinHeapForRT2012(9, n, nBits);

        // a map of nodes for the unassigned pixels
        TIntObjectMap<HeapNode> unAMap = new TIntObjectHashMap<HeapNode>();
        
        iter = unassignedSet.iterator();
        while (iter.hasNext()) {
            
            int pixIdx = iter.next();
            
            TIntSet neighbors = adjAssignedMap.get(pixIdx);
            assert(neighbors != null);
            
            int nNeigbhors = neighbors.size();
            
            long key = 8 - nNeigbhors;
            HeapNode node = new HeapNode(key);
            node.setData(Integer.valueOf(pixIdx));
            
            unAMap.put(pixIdx, node);
            
            heap.insert(node);
        }
        
        assert(unassignedSet.size() == heap.getNumberOfNodes());
        
        float[] lab2 = null;
        float[] lab = new float[3];
        
        
        //DEBUG
        PixelHelper ph = new PixelHelper();
        int[] xyTmp = new int[2];
        int[] xyTmp2 = new int[2];
        CIEChromaticity cieC = new CIEChromaticity();
        
        
        while (heap.getNumberOfNodes() > 0) {
                        
            HeapNode node = heap.extractMin();
            
            assert(node != null);
            
            int pixIdx = ((Integer)node.getData()).intValue();
            
 Map<PairInt, List<Float>> dbgLUVs2 = null;
 Map<PairInt, Float> _dbgLUVs2 = null;
 ph.toPixelCoords(pixIdx, w, xyTmp);
 List<Float> dbgLUVs1 = null;
 if ((xyTmp[0] == 41 && xyTmp[1] == 6) || (xyTmp[0] == 40 && xyTmp[1] == 5)
     || (xyTmp[0] == 40 && xyTmp[1] == 6)
     ) {
     dbgLUVs2 = new HashMap<PairInt, List<Float>>();
     _dbgLUVs2 = new HashMap<PairInt, Float>();
     dbgLUVs1 = new ArrayList<Float>();
     float[] tmp = cieC.rgbToCIELUV_WideRangeLightness(
         clrImg.getR(pixIdx), clrImg.getG(pixIdx), clrImg.getB(pixIdx));
     for (float t : tmp) {
         dbgLUVs1.add(new Float(t));
     }
     int z = 0;
 }
 
            clrImg.getHSB(pixIdx, lab);

            TIntSet adjAssigned = adjAssignedMap.get(pixIdx);
            if (adjAssigned.isEmpty()) {
                System.out.println("priority=" + node.getKey()
                   + " nUnassigned remaining=" + 
                    adjAssignedMap.size() + " heap.n=" +
                    heap.getNumberOfNodes());
            }
            assert(!adjAssigned.isEmpty());
            
            double minD = Double.MAX_VALUE;
            int minLabel2 = -1;
            int minDBGPIXIDX2 = -1;

            TIntIterator iter2 = adjAssigned.iterator();
            while (iter2.hasNext()) {
                int pixIdx2 = iter2.next();
                int label2 = labels[pixIdx2];
                lab2 = hsvs[label2];
                
                double diffSum = 0;
                for (int i = 0; i < 3; ++i) {
                    float diff = lab[i] - lab2[i];
                    diffSum += (diff * diff);
                }

                if (diffSum < minD) {
                    minD = diffSum;
                    minLabel2 = label2;
                    minDBGPIXIDX2 = pixIdx2;
                }
                
if (dbgLUVs2 != null) {
List<Float> dbg2 = new ArrayList<Float>();
float[] tmp = cieC.rgbToCIELUV_WideRangeLightness(
    clrImg.getR(pixIdx2), clrImg.getG(pixIdx2), clrImg.getB(pixIdx2));
for (float t : tmp) {
    dbg2.add(new Float(t));
}
//diff should be fraction...max min values
ph.toPixelCoords(pixIdx2, w, xyTmp2);
dbgLUVs2.put(new PairInt(xyTmp2[0], xyTmp2[1]), dbg2);
_dbgLUVs2.put(new PairInt(xyTmp2[0], xyTmp2[1]), 
    cieC.calcNormalizedDifferenceLUV_WideRangeLightness(
        clrImg.getR(pixIdx), clrImg.getG(pixIdx), clrImg.getB(pixIdx),
        clrImg.getR(pixIdx2), clrImg.getG(pixIdx2), clrImg.getB(pixIdx2)
    ));
System.out.println("pixIdx2=" + pixIdx2 + " xy=" + Arrays.toString(xyTmp2));
}

            }
            
if (dbgLUVs2 != null) {   
    int z = 0;
}

            labels[pixIdx] = minLabel2;
     
            adjAssignedMap.remove(pixIdx);
            
            unAMap.remove(pixIdx);
            
            // update the adjacent unassigned pixels and their keys in heap.
            // these pixels are not in adjLabels.
            TIntSet adjUnassigned = adjUnassignedMap.get(pixIdx);
            if (adjUnassigned == null) {
                continue;
            }
            adjUnassignedMap.remove(pixIdx);
            
            iter2 = adjUnassigned.iterator();
            while (iter2.hasNext()) {
                
                // this is an unassigned pixel
                int pixIdx2 = iter2.next();
                assert(pixIdx != pixIdx2);
                
                // pixIdx should be in it's adj unassigned and then removed
                TIntSet adj2 = adjUnassignedMap.get(pixIdx2);
                assert(adj2 != null);
                boolean rmvd = adj2.remove(pixIdx);
                assert(rmvd);
                
                // add pixIdx to it's assigned pixels
                adj2 = adjAssignedMap.get(pixIdx2);
                assert(adj2 != null);
                adj2.add(pixIdx);

                HeapNode node2 = unAMap.get(pixIdx2);
                assert(node2 != null);
                long key2 = 8 - adj2.size();
                assert(key2 > -1);
                if (key2 < node2.getKey()) {
                    heap.decreaseKey(node2, key2);
                }
            }
        }

        assert(adjAssignedMap.isEmpty());

        iter = unassignedSet.iterator();
        while (iter.hasNext()) {
            int pixIdx = iter.next();
            int label = labels[pixIdx];
            contiguous.get(label).add(pixIdx);
        }
    }
    
    private void _assignTheUnassigned(List<TIntSet> contiguous,
        int[] labels, float[][] luvws, TIntSet unassignedSet) {

        int w = clrImg.getWidth();
        int h = clrImg.getHeight();
        int n = clrImg.getNPixels();

        int maxValue = Math.max(w, h);
        int nBits = 1 + (int)Math.ceil(Math.log(maxValue)/Math.log(2));
        if (nBits > 31) {
            nBits = 31;
        }
        
        //System.out.println("assign " + unassignedSet.size() + " out of "
        //    + n + " pixels");
        
        // key = pixel index of unassigned, value = adj pixels that are assigned
        TIntObjectMap<TIntSet> adjAssignedMap = new TIntObjectHashMap<TIntSet>();
        
        // key = pixel index of unassigned, value = adj pixels that are unassigned
        TIntObjectMap<TIntSet> adjUnassignedMap = new TIntObjectHashMap<TIntSet>();

        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;

        TIntIterator iter = unassignedSet.iterator();
        while (iter.hasNext()) {
            int pixIdx = iter.next();
            addNeighborLabelsForPoint(labels, adjAssignedMap, 
                adjUnassignedMap, pixIdx, dxs, dys);
        }
        
        // using a min heap whose priority is to set the nodes
        //    which have the largest number of assigned neighbors.
        //    nAssigned=8 -> key=8-nAssigned = 0.
        MinHeapForRT2012 heap = new MinHeapForRT2012(9, n, nBits);

        // a map of nodes for the unassigned pixels
        TIntObjectMap<HeapNode> unAMap = new TIntObjectHashMap<HeapNode>();
        
        iter = unassignedSet.iterator();
        while (iter.hasNext()) {
            
            int pixIdx = iter.next();
            
            TIntSet neighbors = adjAssignedMap.get(pixIdx);
            assert(neighbors != null);
            
            int nNeigbhors = neighbors.size();
            
            long key = 8 - nNeigbhors;
            HeapNode node = new HeapNode(key);
            node.setData(Integer.valueOf(pixIdx));
            
            unAMap.put(pixIdx, node);
            
            heap.insert(node);
        }
        
        assert(unassignedSet.size() == heap.getNumberOfNodes());
        
        float[] luv2 = null;
        float[] luv = new float[3];
        
        
        //DEBUG
        PixelHelper ph = new PixelHelper();
        int[] xyTmp = new int[2];
        int[] xyTmp2 = new int[2];
        CIEChromaticity cieC = new CIEChromaticity();
        
        
        while (heap.getNumberOfNodes() > 0) {
                        
            HeapNode node = heap.extractMin();
            
            assert(node != null);
            
            int pixIdx = ((Integer)node.getData()).intValue();
            
            luv = cieC.rgbToCIELUV_WideRangeLightness(
                clrImg.getR(pixIdx), clrImg.getG(pixIdx), clrImg.getB(pixIdx));

            TIntSet adjAssigned = adjAssignedMap.get(pixIdx);
            if (adjAssigned.isEmpty()) {
                System.out.println("priority=" + node.getKey()
                   + " nUnassigned remaining=" + 
                    adjAssignedMap.size() + " heap.n=" +
                    heap.getNumberOfNodes());
            }
            assert(!adjAssigned.isEmpty());
            
            double minD = Double.MAX_VALUE;
            int minLabel2 = -1;

            TIntIterator iter2 = adjAssigned.iterator();
            while (iter2.hasNext()) {
                int pixIdx2 = iter2.next();
                int label2 = labels[pixIdx2];
                luv2 = luvws[label2];
                
                double diffSum = 0;
                for (int i = 0; i < 3; ++i) {
                    float diff = cieC.calcNormalizedDifferenceLUV_WideRangeLightness(luv, luv2);
                    diffSum += (diff * diff);
                }

                if (diffSum < minD) {
                    minD = diffSum;
                    minLabel2 = label2;
                }
            }

            labels[pixIdx] = minLabel2;
     
            adjAssignedMap.remove(pixIdx);
            
            unAMap.remove(pixIdx);
            
            // update the adjacent unassigned pixels and their keys in heap.
            // these pixels are not in adjLabels.
            TIntSet adjUnassigned = adjUnassignedMap.get(pixIdx);
            if (adjUnassigned == null) {
                continue;
            }
            adjUnassignedMap.remove(pixIdx);
            
            iter2 = adjUnassigned.iterator();
            while (iter2.hasNext()) {
                
                // this is an unassigned pixel
                int pixIdx2 = iter2.next();
                assert(pixIdx != pixIdx2);
                
                // pixIdx should be in it's adj unassigned and then removed
                TIntSet adj2 = adjUnassignedMap.get(pixIdx2);
                assert(adj2 != null);
                boolean rmvd = adj2.remove(pixIdx);
                assert(rmvd);
                
                // add pixIdx to it's assigned pixels
                adj2 = adjAssignedMap.get(pixIdx2);
                assert(adj2 != null);
                adj2.add(pixIdx);

                HeapNode node2 = unAMap.get(pixIdx2);
                assert(node2 != null);
                long key2 = 8 - adj2.size();
                assert(key2 > -1);
                if (key2 < node2.getKey()) {
                    heap.decreaseKey(node2, key2);
                }
            }
        }

        assert(adjAssignedMap.isEmpty());

        iter = unassignedSet.iterator();
        while (iter.hasNext()) {
            int pixIdx = iter.next();
            int label = labels[pixIdx];
            contiguous.get(label).add(pixIdx);
        }
    }

    private void addNeighborLabelsForPoint(int[] labels,
        TIntObjectMap<TIntSet> adjAssignedMap, 
        TIntObjectMap<TIntSet> adjUnassignedMap,
        int pixIdx, int[] dxs, int[] dys) {

        int w = clrImg.getWidth();
        int h = clrImg.getHeight();

        TIntSet adjLabels = adjAssignedMap.get(pixIdx);
        TIntSet adjULabels = adjUnassignedMap.get(pixIdx);
        if (adjLabels == null) {
            assert(adjULabels == null);
            
            adjLabels = new TIntHashSet();
            adjAssignedMap.put(pixIdx, adjLabels);
            
            adjULabels = new TIntHashSet();
            adjUnassignedMap.put(pixIdx, adjULabels);
        }
        
        int j = pixIdx/w;
        int i = pixIdx - (j * w);

        for (int m = 0; m < dxs.length; ++m) {
            int x2 = i + dxs[m];
            int y2 = j + dys[m];
            if (x2 < 0 || y2 < 0 || (x2 > (w - 1)) || (y2 > (h - 1))) {
                continue;
            }
            int pixIdx2 = (y2 * w) + x2;
            if (labels[pixIdx2] == -1) {
                adjULabels.add(pixIdx2);
            } else {
                adjLabels.add(pixIdx2);
            }
        }
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

    private int[] getRegionHistogram(HGS hgs0, TIntSet pIdxs) {
        int[] h0 = new int[hgs0.getNumberOfBins()];
        TIntIterator iter = pIdxs.iterator();
        while (iter.hasNext()) {
            int pixIdx = iter.next();
            int y = pixIdx/gsImg.getWidth();
            int x = pixIdx - (y * gsImg.getWidth());
            hgs0.extractBlock2(x, y, h0);
        }
        return h0;
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
    
    private List<List<TIntSet>> extractontiguousSegments(
        List<TIntSet> unmatchedRMap) {
       
        int n = unmatchedRMap.size();
        
        System.out.println("number of unmatched regions=" + n);
        
        List<List<TIntSet>> output = new 
            ArrayList<List<TIntSet>>();
        
        for (int i = 0; i < unmatchedRMap.size(); ++i) {
                        
            TIntSet idxs = unmatchedRMap.get(i);
            
            ConnectedPointsFinder finder3 = new ConnectedPointsFinder(
                clrImg.getWidth(), clrImg.getHeight());
            finder3.setMinimumNumberInCluster(1);
            finder3.setToUse8Neighbors();
            finder3.findConnectedPointGroups(idxs);
            
            List<TIntSet> list = new ArrayList<TIntSet>();
            output.add(list);
            
            int n3 = finder3.getNumberOfGroups();
            for (int ii = 0; ii < n3; ++ii) {
                
                TIntSet segment = finder3.getXY(ii);
                
                list.add(segment);
            }
        }
       
        return output;
    }
    
    private void addUnmatchedToEdgePoints(TIntSet unmatchedPoints, 
        List<TIntSet> unmatchedRMap, TIntSet allEdgePoints, 
        TIntSet rmvdImgBorders, int maxGapSize) {
        
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
        List<List<TIntSet>> unmatchedRSegments
            = extractontiguousSegments(unmatchedRMap);

        for (List<TIntSet> segmentList : unmatchedRSegments) {

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

                assert (!umEPIdxMap.isEmpty());

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

    }
    
    public EdgeFilterProducts getEdgeFilterProducts() {
        return edgeProducts;
    }
}
