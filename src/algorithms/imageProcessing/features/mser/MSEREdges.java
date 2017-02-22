package algorithms.imageProcessing.features.mser;

import algorithms.QuickSort;
import algorithms.compGeometry.PerimeterFinder2;
import algorithms.imageProcessing.DFSConnectedGroupsFinder;
import algorithms.imageProcessing.DFSContiguousValueFinder;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.GroupPixelHSV;
import algorithms.imageProcessing.GroupPixelRGB0;
import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.ImageProcessor.Colors;
import algorithms.imageProcessing.ImageSegmentation;
import algorithms.imageProcessing.MiscellaneousCurveHelper;
import algorithms.imageProcessing.PostLineThinnerCorrections;
import algorithms.imageProcessing.SummedAreaTable;
import algorithms.imageProcessing.features.mser.Canonicalizer.RegionGeometry;
import algorithms.imageProcessing.features.mser.MSER.Threshold;
import algorithms.imageProcessing.segmentation.ColorSpace;
import algorithms.imageProcessing.segmentation.LabelToColorHelper;
import algorithms.imageProcessing.util.AngleUtil;
import algorithms.misc.Misc;
import algorithms.misc.MiscDebug;
import algorithms.misc.MiscMath;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.QuadInt;
import gnu.trove.iterator.TIntIntIterator;
import gnu.trove.iterator.TIntIterator;
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
import java.awt.Color;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.Stack;
import thirdparty.edu.princeton.cs.algs4.Interval;
import thirdparty.edu.princeton.cs.algs4.Interval2D;
import thirdparty.edu.princeton.cs.algs4.QuadTree;

/**
 * class to explore the boundaries of the accumulated points in MSER regions
 * to make edges that are complete contours (the level sets approach has a
 * much larger pixel association than the few pixels of canny edges for
 * example).
 *
 * It currently uses the level sets found in the MSER regions made from
 * greyscale and "H" of LCH color space images (caveat, configured
 * to extract MSER regions from images from COTS of past several years
 * binned to size near 256 X 256).
 * The boundaries of the accumulated points in the regions are extracted
 * and the high scoring points are kept where score is a combination
 * of sobel color contrast and sobel greyscale intensity values).
 * 
 * The class in not ready for use yet.
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
    private List<Set<PairInt>> edgeList = null;

    private List<Set<PairInt>> labeledSets = null;

    private boolean debug = false;

    private boolean useLowerContrastLimits = false;
    
    private GreyscaleImage sobelScores = null;

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

        this.labeledSets = mergeRegions2();

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
                    && r.getVariation() > 0.001) {
                    list.remove(i);
                } else if ((type == 0) 
                    && r.getVariation() == 0.0) {
                    list.remove(i);
                } else {
                    // check for this Region being an artifact of wrap around
                    if (r.level_ < 25) {
                        list.remove(i);
                        hadWrapAroundArtifacts = true;
                    } else {
                        int avgLevel;
                        if (type == 0) {
                            avgLevel = calcAvg(ptImg, r.getAcc());
                        } else {
                            avgLevel = calcAvg(negImg, r.getAcc());
                        }
                        //System.out.format(" %d add x,y=%d,%d level=%d  acgLevel=%d\n",
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
                        && r.getVariation() > 0.001) {
                        list2.remove(i);
                    } else if ((type == 0)
                        && r.getVariation() == 0.0) {
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

    private List<Set<PairInt>> extractContiguous(GreyscaleImage tImg,
        int value, int minGroupSize) {

        List<Set<PairInt>> out = new ArrayList<Set<PairInt>>();

        DFSContiguousValueFinder dfsFinder = new DFSContiguousValueFinder(tImg);
        dfsFinder.setMinimumNumberInCluster(minGroupSize);
        dfsFinder.findGroups(value);

        int nGroups = dfsFinder.getNumberOfGroups();

        for (int j = 0; j < nGroups; ++j) {
            PairIntArray xy = dfsFinder.getXY(j);
            out.add(Misc.convert(xy));
        }

        return out;
    }

    private void extractBoundaries() {

        //INITIALIZED, REGIONS_EXTRACTED, MERGED, EDGES_EXTRACTED
        if (state.equals(STATE.EDGES_EXTRACTED)) {
            throw new IllegalStateException("can only perform extraction of "
                + "edges once");
        }

        // edgeList and labeledSets are made in this:
        Set<PairInt> thinned = combineAndThinBoundaries();

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
    }

    private class OneDFloatArray {
        float[] a;
        public OneDFloatArray(float[] a) {
            this.a = a;
        }
    }

    private List<Set<PairInt>> reduceToUniquePointToEdge() {

        Set<PairInt> edges = new HashSet<PairInt>();
        for (int i = 0; i < edgeList.size(); ++i) {
            Set<PairInt> set = edgeList.get(i);
            for (PairInt p : set) {
                edges.add(p);
            }
        }

        ImageProcessor imageProcessor = new ImageProcessor();
        imageProcessor.applyThinning(edges, clrImg.getWidth(), clrImg.getHeight(),
            false);

        int w = clrImg.getWidth();
        int h = clrImg.getHeight();

        Set<PairInt> notEdges = new HashSet<PairInt>();
        for (int i = 0; i < w; ++i) {
            for (int j = 0; j < h; ++j) {
                PairInt p = new PairInt(i, j);
                if (!edges.contains(p)) {
                    notEdges.add(p);
                }
            }
        }

        DFSConnectedGroupsFinder finder = new DFSConnectedGroupsFinder();
        finder.setMinimumNumberInCluster(1);
        finder.findConnectedPointGroups(notEdges);

        TIntIntMap pointIndexMap = new TIntIntHashMap();

        List<Set<PairInt>> contigSets = new ArrayList<Set<PairInt>>();
        for (int i = 0; i < finder.getNumberOfGroups(); ++i) {
            Set<PairInt> group = finder.getXY(i);
            contigSets.add(group);
        }

        List<Colors> clrsList = new ArrayList<Colors>();
        for (int i = 0; i < contigSets.size(); ++i) {

            Set<PairInt> group = contigSets.get(i);

            for (PairInt p : group) {
                int pixIdx = clrImg.getInternalIndex(p);
                pointIndexMap.put(pixIdx, i);
            }

            GroupPixelRGB0 gpb = new GroupPixelRGB0();
            gpb.calculateColors(group, clrImg, 0, 0);
            int r = Math.round(gpb.getAvgRed());
            int g = Math.round(gpb.getAvgGreen());
            int b = Math.round(gpb.getAvgBlue());
            float[] hsb = new float[3];
            Color.RGBtoHSB(r, g, b, hsb);
            Colors clr = new Colors(hsb);
            clrsList.add(clr);
        }

        // -- place edges in the contiguous adjacent set, closest to them
        //    in color
        int[] dxs = Misc.dx4;
        int[] dys = Misc.dy4;
        ArrayDeque<PairInt> queue = new ArrayDeque<PairInt>();
        queue.addAll(edges);

        while (!queue.isEmpty()) {
            PairInt p = queue.pop();
            int pixIdx = clrImg.getInternalIndex(p);
            if (pointIndexMap.containsKey(pixIdx)) {
                continue;
            }
            int x = clrImg.getCol(pixIdx);
            int y = clrImg.getRow(pixIdx);

            float[] hsb = new float[3];
            Color.RGBtoHSB(clrImg.getR(pixIdx), clrImg.getG(pixIdx),
                clrImg.getB(pixIdx), hsb);

            int minClrDiff = Integer.MAX_VALUE;
            int minClrDiffIdx = -1;
            for (int k = 0; k < dxs.length; ++k) {
                int x2 = x + dxs[k];
                int y2 = y + dys[k];
                if (x2 < 0 || y2 < 0 || (x2 >= w) || (y2 >= h)) {
                    continue;
                }
                int pixIdx2 = clrImg.getInternalIndex(x2, y2);
                if (!pointIndexMap.containsKey(pixIdx2)) {
                    // an edge pixel adjacent to another
                    continue;
                }
                int gIdx = pointIndexMap.get(pixIdx2);
                Colors clrs2 = clrsList.get(gIdx);
                int diff = 0;
                for (int j = 0; j < hsb.length; ++j) {
                    diff += Math.abs(hsb[j] - clrs2.getColors()[j]);
                }
                if (diff < minClrDiff) {
                    minClrDiff = diff;
                    minClrDiffIdx = pixIdx2;
                }
            }
            if (minClrDiffIdx == -1) {
                queue.add(p);
                continue;
            }
            int gIdx = pointIndexMap.get(minClrDiffIdx);
            pointIndexMap.put(pixIdx, gIdx);
            contigSets.get(gIdx).add(p);

            edges.remove(p);
            notEdges.add(p);
        }
        assert(edges.isEmpty());

        // extract boundaries of contigSets
        edgeList.clear();

        PerimeterFinder2 finder2 = new PerimeterFinder2();
        for (Set<PairInt> set : contigSets) {
            Set<PairInt> embedded = new HashSet<PairInt>();
            Set<PairInt> outerBorder = new HashSet<PairInt>();
            finder2.extractBorder2(set, embedded, outerBorder);
            edgeList.add(outerBorder);
        }

        return contigSets;
    }

    private List<Set<PairInt>> mergeRegions2() {

        //INITIALIZED, REGIONS_EXTRACTED, MERGED, EDGES_EXTRACTED
        if (!state.equals(STATE.EDGES_EXTRACTED)) {
            throw new IllegalStateException("error in algorithm.  expecting"
                + " edges were extracted.");
        }

        List<Set<PairInt>> contigSets = reduceToUniquePointToEdge();

        int w = clrImg.getWidth();
        int h = clrImg.getHeight();

        // -- merge contigSets
        int sizeLimit = 16;//31;
        if (clrImg.getNPixels() < 100) {
            sizeLimit = 5;
        }

        ImageSegmentation imageSegmentation = new ImageSegmentation();
        float hsvLimit = 0.095f;
        int[] labels = imageSegmentation.mergeByColor(clrImg, contigSets,
            ColorSpace.HSV,
            hsvLimit);//0.1f);

        imageSegmentation.mergeSmallSegments(clrImg, labels, sizeLimit,
            ColorSpace.HSV);

        contigSets = LabelToColorHelper
            .extractContiguousLabelPoints(clrImg, labels);

        labels = imageSegmentation.mergeByColor(clrImg, contigSets,
            ColorSpace.CIELUV_NORMALIZED,
            0.01f);

        labels = imageSegmentation.mergeByColor(clrImg, contigSets,
            ColorSpace.POLAR_CIELUV,
            1.f);// in degrees

        sizeLimit = 24;
        if (clrImg.getNPixels() < 100) {
            sizeLimit = 10;
        }
        imageSegmentation.mergeSmallSegments(clrImg, labels, sizeLimit,
            ColorSpace.CIELAB);

        contigSets = LabelToColorHelper
            .extractContiguousLabelPoints(clrImg, labels);

        // - write the final edges (and clear the regions list)
        edgeList.clear();
        PerimeterFinder2 finder2 = new PerimeterFinder2();
        for (Set<PairInt> set : contigSets) {
            Set<PairInt> embedded = new HashSet<PairInt>();
            Set<PairInt> outerBorder = new HashSet<PairInt>();
            finder2.extractBorder2(set, embedded, outerBorder);
            edgeList.add(outerBorder);

            /*
            DFSConnectedGroupsFinder dfsFinder = new DFSConnectedGroupsFinder();
            dfsFinder.setMinimumNumberInCluster(24);
            dfsFinder.findConnectedPointGroups(embedded);
            for (int j = 0; j < dfsFinder.getNumberOfGroups(); ++j) {

                Set<PairInt> eSet = dfsFinder.getXY(j);

                Set<PairInt> embedded2 = new HashSet<PairInt>();
                Set<PairInt> outerBorder2 = new HashSet<PairInt>();
                finder2.extractBorder2(eSet, embedded2, outerBorder2);

                if (outerBorder2.size() > 16) {
                    edgeList.add(outerBorder2);
                }
            }*/
        }
        return contigSets;
    }

    private List<List<Region>> filterOverlapping(List<List<Region>> rlist,
        int w, int h) {

        List<List<Region>> filtered = new ArrayList<List<Region>>();
        for (int i = 0; i < rlist.size(); ++i) {

            List<Region> regions = rlist.get(i);

            List<Region> out = new ArrayList<Region>();
            filtered.add(out);

            int[] xyCen = new int[2];

            List<EllipseHelper> hs = new ArrayList<EllipseHelper>();
            for (int j = 0; j < regions.size(); ++j) {
                Region r = regions.get(j);
                r.calculateXYCentroid(xyCen, w, h);
                double[] coeffs = r.calcParamTransCoeff();
                EllipseHelper eh = new EllipseHelper(xyCen[0], xyCen[1], coeffs);
                hs.add(eh);
            }

            // filter Regions to only keep those without another Region center
            // in them
            for (int j = 0; j < regions.size(); ++j) {
                Region r1 = regions.get(j);
                EllipseHelper eh1 = hs.get(j);

                boolean noneInternal = true;

                for (int k = 0; k < regions.size(); ++k) {
                    if (j == k) {
                        continue;
                    }

                    Region r2 = regions.get(k);
                    EllipseHelper eh2 = hs.get(k);
                    r2.calculateXYCentroid(xyCen, w, h);

                    if (eh1.isWithin(xyCen[0], xyCen[1])
                        && !eh2.surrounds(eh1)) {
                        noneInternal = false;
                        break;
                    }
                }

                if (noneInternal) {
                    out.add(r1);
                }
            }
        }

        return filtered;
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

    public List<Set<PairInt>> getEdges() {

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
    public List<Set<PairInt>> getLabeledSets() {

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
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        float[] sobelScores = imageProcessor.createSobelColorScores(
            gsImg, ptImg, 20);

        int w = gsImg.getWidth();
        int h = gsImg.getHeight();
        GreyscaleImage scaled = MiscMath.rescaleAndCreateImage(sobelScores,
            w, h);

        // smearing values over a 3 pixel window
        SummedAreaTable sumTable = new SummedAreaTable();

        GreyscaleImage imgM = sumTable.createAbsoluteSummedAreaTable(scaled);
        imgM = sumTable.applyMeanOfWindowFromSummedAreaTable(imgM, 3);

        return imgM;
    }

    private Set<PairInt> combineAndThinBoundaries() {

        if (sobelScores == null) {
            sobelScores = createSobelScores();
            MiscDebug.writeImage(sobelScores, "_" + ts + "_sobel_");
        }
        
        // below this removes:
        double limit = 76;//51;//76;//12.75;
        // above this restores
        double limit2 = 5;
        
        //the intersection of overlapping regions, present in many regions
        //   is often a strong edge
        List<Set<PairInt>> rmvd = new ArrayList<Set<PairInt>>();
        
        List<Set<PairInt>> boundaries = new ArrayList<Set<PairInt>>();

        PerimeterFinder2 finder = new PerimeterFinder2();

        TIntObjectMap<TIntList> pointIndexesMap = new TIntObjectHashMap<TIntList>();
        Set<PairInt> allPoints = new HashSet<PairInt>();

        for (int rListIdx = 0; rListIdx < regions.size(); ++rListIdx) {
            Region r = regions.get(rListIdx);
            Set<PairInt> points = r.getAcc();
            Set<PairInt> embedded = new HashSet<PairInt>();
            Set<PairInt> outerBorder = new HashSet<PairInt>();
            finder.extractBorder2(points, embedded, outerBorder);

            double score = calcAvgScore(outerBorder, sobelScores);
            boolean doNotAdd = (score < limit);
            if (doNotAdd) {
                rmvd.add(outerBorder);
            } else {
                boundaries.add(outerBorder);
                allPoints.addAll(outerBorder);
            }
            
            // TODO: revisit this.  when are these vars used?
            for (PairInt p : outerBorder) {
                int pixIdx = clrImg.getInternalIndex(p);
                TIntList bIdxs = pointIndexesMap.get(pixIdx);
                if (bIdxs == null) {
                    bIdxs = new TIntArrayList();
                    pointIndexesMap.put(pixIdx, bIdxs);
                }
                bIdxs.add(rListIdx);
            }

            DFSConnectedGroupsFinder dfsFinder = new DFSConnectedGroupsFinder();
            dfsFinder.setMinimumNumberInCluster(12);
            dfsFinder.findConnectedPointGroups(embedded);
            for (int j = 0; j < dfsFinder.getNumberOfGroups(); ++j) {

                Set<PairInt> eSet = dfsFinder.getXY(j);

                Set<PairInt> embedded2 = new HashSet<PairInt>();
                Set<PairInt> outerBorder2 = new HashSet<PairInt>();
                finder.extractBorder2(eSet, embedded2, outerBorder2);

                if (outerBorder2.size() > 12) {

                    if (doNotAdd) {
                        rmvd.add(outerBorder2);
                    } else {
                        allPoints.addAll(outerBorder2);
                    }
                    
                    for (PairInt p : outerBorder2) {
                        int pixIdx = clrImg.getInternalIndex(p);
                        TIntList bIdxs = pointIndexesMap.get(pixIdx);
                        if (bIdxs == null) {
                            bIdxs = new TIntArrayList();
                            pointIndexesMap.put(pixIdx, bIdxs);
                        }
                        bIdxs.add(rListIdx);
                    }
                }
            }
        }
        
 //TODO: adjust score limits and angle thhresh
        
        if (rmvd.size() > 0) {
            for (Set<PairInt> set : rmvd) {
                for (PairInt p : set) {
                    int score = sobelScores.getValue(p);
                    System.out.println(" ? " + p + " score=" + score);
                    if (score > limit2) {
                        allPoints.add(p);
                    }
                }
            }
        }
        
        if (debug) {
            Image tmp = clrImg.copyImage();
            ImageIOHelper.addCurveToImage(allPoints, tmp, 0, 255, 0, 0);
            MiscDebug.writeImage(tmp, "_" + ts + "_TMP_");
        }
        
        GreyscaleImage img2 = new GreyscaleImage(clrImg.getWidth(),
            clrImg.getHeight());
        for (PairInt p : allPoints) {
            img2.setValue(p.getX(), p.getY(), 1);
        }

        ImageSegmentation imageSegmentation = new ImageSegmentation();
        Set<PairInt> outputAddedGaps = new HashSet<PairInt>();
        img2 = imageSegmentation.fillInGapsOf1(img2, outputAddedGaps, 1);

        // restore gap where the gap is completely surrounded
        imageSegmentation.restoreGapsOf1WhereSurrounded(img2, outputAddedGaps, 1);

        Set<PairInt> thinned = new HashSet<PairInt>();
        for (int i = 0; i < img2.getNPixels(); ++i) {
            // img2 edges have pixel value=1
            if (img2.getValue(i) > 0) {
                thinned.add(new PairInt(img2.getCol(i), img2.getRow(i)));
            }
        }

        if (debug) {
            MiscDebug.writeImage(clrImg, "_" + ts + "_0_");
            Image tmp = clrImg.copyImage();
            for (int i = 0; i < clrImg.getNPixels(); ++i) {
                // img2 has edges w/ value=1
                if (img2.getValue(i) > 0) {
                    tmp.setRGB(i, 255, 0, 0);
                }
            }
            MiscDebug.writeImage(tmp, "_" + ts + "_closing_");
        }
        
        // find clusters (contiguous pixels of value 0) between edges
        List<Set<PairInt>> contigousSets = extractContiguous(img2, 0, 3);//9);
        
        if (debug) {
            Image tmp = clrImg.copyImage();
            ImageIOHelper.addAlternatingColorPointSetsToImage(
                contigousSets, 0, 0, 0, tmp);
            MiscDebug.writeImage(tmp, "_" + ts + "_expanded_");
        }

        assignTheUnassigned(contigousSets);

        if (debug) {
            Image tmp = clrImg.copyImage();
            ImageIOHelper.addAlternatingColorPointSetsToImage(
                contigousSets, 0, 0, 0, tmp);
            MiscDebug.writeImage(tmp, "_" + ts + "_assigned_");
        }
        
        PerimeterFinder2 finder2 = new PerimeterFinder2();
        
        // --- extract the bounds, and if any are empty of internal points,
        //     re-submit those to be reassigned to a cluster which does
        List<Set<PairInt>> extractedBoundaries = new ArrayList<Set<PairInt>>();
        for (int i = (contigousSets.size() - 1); i > -1; --i) {
            Set<PairInt> set = contigousSets.get(i);
            Set<PairInt> embedded = new HashSet<PairInt>();
            Set<PairInt> outerBorder = new HashSet<PairInt>();
            finder2.extractBorder2(set, embedded, outerBorder);
            
            // if outerBorder is same size as set, there are no internal
            // points, so remove the set so can reassign it
            if (set.size() - outerBorder.size() < 2) {
                contigousSets.remove(i);
            } else {
                extractedBoundaries.add(outerBorder);
            }
        }
        
        thinned.clear();
        this.labeledSets = contigousSets;
        edgeList = new ArrayList<Set<PairInt>>();
        
        if (contigousSets.size() < extractedBoundaries.size()) {
            
            //TODO: this needs to be edited to improve the
            // order of assignments
            assignTheUnassigned(contigousSets);
            
            for (Set<PairInt> set : contigousSets) {
                Set<PairInt> embedded = new HashSet<PairInt>();
                Set<PairInt> outerBorder = new HashSet<PairInt>();
                finder2.extractBorder2(set, embedded, outerBorder);

                // for small regions, sometimes outerBorder has no inner
                // points, making the set a 2 pixel thick edge,
                // so for those, need to reassign to the closest
                // set in color            
                edgeList.add(outerBorder);
            }
        } else {
            edgeList = extractedBoundaries;
        }
                
        for (Set<PairInt> set : edgeList) {
            thinned.addAll(set);
        }
        
        if (debug) {
            Image tmp = clrImg.copyImage();
            ImageIOHelper.addCurveToImage(thinned, tmp, 0, 255, 0, 0);
            MiscDebug.writeImage(tmp, "_" + ts + "_thinned_0_");
        }

        return thinned;
    }

    private void assignTheUnassigned(List<Set<PairInt>> contiguous) {

        float[][] hsvs = new float[contiguous.size()][];

        int[] labels = new int[clrImg.getNPixels()];
        Arrays.fill(labels, -1);
        for (int label = 0; label < contiguous.size(); ++label) {
            Set<PairInt> set = contiguous.get(label);
            for (PairInt p : set) {
                int pixIdx = clrImg.getInternalIndex(p);
                labels[pixIdx] = label;
            }

            GroupPixelHSV hsv = new GroupPixelHSV();
            hsv.calculateColors(set, clrImg);
            hsvs[label] = new float[]{hsv.getAvgH(), hsv.getAvgS(),
                hsv.getAvgV()};
        }

        Map<PairInt, TIntSet> unassignedMap = new HashMap<PairInt, TIntSet>();

        int w = clrImg.getWidth();
        int h = clrImg.getHeight();
        int n = clrImg.getNPixels();

        int[] dxs = Misc.dx8;
        int[] dys = Misc.dy8;

        for (int i = 0; i < w; ++i) {
            for (int j = 0; j < h; ++j) {
                int pixIdx = clrImg.getInternalIndex(i, j);
                if (labels[pixIdx] == -1) {
                    addNeighborLabelsForPoint(labels, unassignedMap, i, j, dxs, dys);
                }
            }
        }
        
        ArrayDeque<PairInt> queue0 = populateByNumberOfNeighbors(unassignedMap);

        ArrayDeque<PairInt> queue1 = new ArrayDeque<PairInt>();

        int nIter = 0;

        float[] lab2 = null;

        Set<PairInt> visited = new HashSet<PairInt>();

        while (!queue0.isEmpty() || !queue1.isEmpty()) {

            PairInt p;
            if (!queue1.isEmpty()) {
                p = queue1.poll();
            } else {
                p = queue0.poll();
            }

            if (visited.contains(p)) {
                continue;
            }
            visited.add(p);

            int x1 = p.getX();
            int y1 = p.getY();

            TIntSet adjLabels;
            if (nIter == 0) {
                adjLabels = unassignedMap.get(p);
                assert (adjLabels != null);
            } else {
                adjLabels = new TIntHashSet();
                addNeighborLabelsForPoint(labels, adjLabels, x1, y1, dxs, dys);
            }

            double minD = Double.MAX_VALUE;
            int minLabel2 = -1;

            float[] lab = new float[3];
            lab[0] = clrImg.getHue(x1, y1);
            lab[1] = clrImg.getSaturation(x1, y1);
            lab[2] = clrImg.getBrightness(x1, y1);

            TIntIterator iter = adjLabels.iterator();
            while (iter.hasNext()) {
                int label2 = iter.next();
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

            int pixIdx1 = clrImg.getInternalIndex(p.getX(), p.getY());
            labels[pixIdx1] = minLabel2;

            unassignedMap.remove(p);

            for (int m = 0; m < dxs.length; ++m) {
                int x2 = p.getX() + dxs[m];
                int y2 = p.getY() + dys[m];
                if (x2 < 0 || y2 < 0 || (x2 > (w - 1)) || (y2 > (h - 1))) {
                    continue;
                }
                int pixIdx2 = clrImg.getInternalIndex(x2, y2);
                if (labels[pixIdx2] == -1) {
                    PairInt p2 = new PairInt(x2, y2);
                    queue1.add(p2);
                    //assert (!visited.contains(p2));
                }
            }
            nIter++;
        }

        assert(unassignedMap.isEmpty());

        for (int pixIdx = 0; pixIdx < labels.length; ++pixIdx) {

            int label = labels[pixIdx];
            PairInt p = new PairInt(clrImg.getCol(pixIdx), clrImg.getRow(pixIdx));

            contiguous.get(label).add(p);
        }
    }

    private void addNeighborLabelsForPoint(int[] labels,
        Map<PairInt, TIntSet> unassignedMap,
        int i, int j, int[] dxs, int[] dys) {

        int w = clrImg.getWidth();
        int h = clrImg.getHeight();

        PairInt p = new PairInt(i, j);

        TIntSet adjLabels = unassignedMap.get(p);
        if (adjLabels == null) {
            adjLabels = new TIntHashSet();
            unassignedMap.put(p, adjLabels);
        }

        addNeighborLabelsForPoint(labels, adjLabels, i, j, dxs, dys);
    }

    private void addNeighborLabelsForPoint(int[] labels, TIntSet adjLabels,
        int i, int j, int[] dxs, int[] dys) {

        int w = clrImg.getWidth();
        int h = clrImg.getHeight();

        PairInt p = new PairInt(i, j);

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

    private ArrayDeque<PairInt> populateByNumberOfNeighbors(
        Map<PairInt, TIntSet> unassignedMap) {

        int n = unassignedMap.size();

        PairInt[] points = new PairInt[n];
        int[] nN = new int[n];

        int count = 0;
        for (Map.Entry<PairInt, TIntSet> entry : unassignedMap.entrySet()) {
            points[count] = entry.getKey();
            nN[count] = entry.getValue().size();
            count++;
        }

        QuickSort.sortBy1stArg(nN, points);

        ArrayDeque<PairInt> queue = new ArrayDeque<PairInt>();

        for (int i = (n - 1); i > -1; --i) {

            int nP = nN[i];
            if (nP == 0) {
                break;
            }

            queue.add(points[i]);
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

    private double calcAvgScore(Set<PairInt> points, GreyscaleImage sobelScores) {
   
        double sum = 0;
        for (PairInt p : points) {
            sum += sobelScores.getValue(p);
        }
        
        sum /= (double) points.size();
        
        return sum;
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
    
    private int calcAvg(GreyscaleImage img, Set<PairInt> points) {

        double sum = 0;
        for (PairInt p : points) {
            sum += img.getValue(p);
        }
        sum /= (double)points.size();
        
        return (int)Math.round(sum);
    }    
}
