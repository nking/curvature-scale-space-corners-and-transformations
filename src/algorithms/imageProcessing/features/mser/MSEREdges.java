package algorithms.imageProcessing.features.mser;

import algorithms.QuickSort;
import algorithms.compGeometry.PerimeterFinder2;
import algorithms.imageProcessing.DFSConnectedGroupsFinder;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.GroupPixelHSV;
import algorithms.imageProcessing.GroupPixelRGB0;
import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.ImageProcessor.Colors;
import algorithms.imageProcessing.ImageSegmentation;
import algorithms.imageProcessing.features.mser.MSER.Threshold;
import algorithms.imageProcessing.segmentation.ColorSpace;
import algorithms.imageProcessing.segmentation.LabelToColorHelper;
import algorithms.imageProcessing.util.AngleUtil;
import algorithms.misc.Misc;
import algorithms.misc.MiscDebug;
import algorithms.misc.MiscMath;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import gnu.trove.iterator.TIntIntIterator;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.TFloatList;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.awt.Color;
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
 * class to explore the boundaries of the accumulated points in MSER regions.
 * 
 * The contents may in the future be replaced by an implementation of 
 * http://www.vision.ee.ethz.ch/~rhayko/paper/aapr2009_boundary_detection_SBER_hayko.pdf
 * 
 * but for now is an exploration of the existing MSER and Region class results.
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
    private List<Region> regions = null;
    
    //a re-extraction of regions on the gs positive image, but using the
    // MSER parameters used for the negative image.
    // this is not useful in edges, but is useful for other methods.
    private List<Region> sensitiveGS0 = null;

    // the original regions for gs positive, negative, then pt positive
    // and negative.  regions and filteredRegions should be prefered for most
    // uses.
    private List<List<Region>> origGsPtRegions = null;
    
    // color contrast regions filtered to the smallest that
    // do not have any other ellipse centers within them.
    // (then all positive filtered and negative filtered are compared
    // and kept, but only the largest if intersecting are kept)
    private List<Region> filteredRegions = null;
    
    // NOTE: the edges indexes do not correspond to the regions indexes
    private List<Set<PairInt>> edgeList = null;
    
    private List<Set<PairInt>> labeledSets = null;
    
    private boolean debug = false;
    
    private boolean useLowerContrastLimits = false;
    
    private long ts = 0;
    
    private STATE state = null;
    
    public MSEREdges(ImageExt img) {
        
        ImageProcessor imageProcesor = new ImageProcessor();
        
        this.gsImg = img.copyToGreyscale2();
        this.ptImg = imageProcesor.createCIELUVTheta(img, 255);
        this.clrImg = img.copyToImageExt();
        
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
        
        useFilterOnEdges();
        
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
        
        //NOTE: this method should be cleaned up when have the better pattern
        //   of mser region extraction and filtering
        
        //INITIALIZED, REGIONS_EXTRACTED, MERGED, EDGES_EXTRACTED
        if (state.equals(STATE.EDGES_EXTRACTED)) {
            throw new IllegalStateException("can only perform extraction of "
                + "edges once");
        }
        
        List<List<Region>> gsRegions = extractMSERRegions(gsImg, Threshold.LESS_SENSITIVE);

        List<List<Region>> ptRegions = extractMSERRegions(ptImg, Threshold.DEFAULT);

        regions = new ArrayList<Region>();

        origGsPtRegions = new ArrayList<List<Region>>();
        
        // for minArea 0.001, a gsRegions limit for var of 0. is used
        //   but for 0.0001, limit should be near 0.001

        for (int type = 0; type < 2; ++type) {
            List<Region> list = gsRegions.get(type);
            // temporary change to look at including more sensitive gs0
            if (type==0) {
                list = _extractSensitiveGS0();
            }
            for (int i = (list.size() - 1); i > -1; --i) {
                Region r = list.get(i);
                if ((type == 1) && r.getVariation() > 2.) {
                    list.remove(i);
                } else if ((type == 0)
                    //&& r.getVariation() == 0.0) {
                    && r.getVariation() < 2.) {
                    list.remove(i);
                } else {
                    GroupPixelHSV gHSV = new GroupPixelHSV();
                    gHSV.calculateColors(r.getAcc(), clrImg);
                    if (gHSV.getStdDevV() > 0.3){
                        list.remove(i);
                    }
                }
            }
            if (type == 1) {
                List<TIntList> concList = getEmbeddedLevels(list);
                TIntList rmList = new TIntArrayList();
                for (int k = 0; k < concList.size(); ++k) {
                    for (int i3 = 1; i3 < concList.get(k).size(); ++i3) {
                        rmList.add(concList.get(k).get(i3));
                    }
                }
                rmList.sort();
                for (int k = rmList.size() - 1; k > -1; --k) {
                    int rmIdx = rmList.get(k);
                    list.remove(rmIdx);
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
        
        List<TIntList> concList = getEmbeddedLevels(regions);
        TIntList rmList = new TIntArrayList();
        for (int k = 0; k < concList.size(); ++k) {
            for (int i3 = 1; i3 < concList.get(k).size(); ++i3) {
                rmList.add(concList.get(k).get(i3));
            }
        }
        rmList.sort();
        for (int k = rmList.size() - 1; k > -1; --k) {
            int rmIdx = rmList.get(k);
            regions.remove(rmIdx);
        }
            
        for (int type = 0; type < 2; ++type) {
            List<Region> list = ptRegions.get(type);
            for (int i = (list.size() - 1); i > -1; --i) {
                Region r = list.get(i);
                if ((type == 1) && r.getVariation() > 0.001) {
                    list.remove(i);
                } else if ((type == 0) && r.getVariation() == 0.0) {
                    list.remove(i);
                } else {
                    regions.add(r);
                }
            }
            // copy list into origGsPtRegions
            List<Region> cpList = new ArrayList<Region>();
            origGsPtRegions.add(cpList);
            for (Region r : list) {
                cpList.add(r.copy());
            }
        }

        if (debug) {
            int[] xyCen = new int[2];
            Image imCp;
            for (int type = 0; type < 2; ++type) {
                imCp = gsImg.copyToColorGreyscale();
                int n = gsRegions.get(type).size();
                for (int i = 0; i < n; ++i) {
                    Region r = gsRegions.get(type).get(i);
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
                int n = ptRegions.get(type).size();
                for (int i = 0; i < n; ++i) {
                    Region r = ptRegions.get(type).get(i);
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

        // ---- make a list of filtered color contrast regions while still have
        //      the individual Region lists

        boolean additionalFiltering = true;

        int w = ptImg.getWidth();
        int h = ptImg.getHeight();
        
        if (additionalFiltering) {

            List<List<Region>> filtered = filterOverlapping(ptRegions, w, h);

            // --- compare the two lists of filtered and if there is
            //     an intersection, keep the largest

            List<Region> regions0 = filtered.get(0);

            TIntSet skip0 = new TIntHashSet();
            TIntSet skip1 = new TIntHashSet();

            int[] xyCen = new int[2];

            List<EllipseHelper> hs0 = new ArrayList<EllipseHelper>();
            for (int j = 0; j < regions0.size(); ++j) {
                Region r = regions0.get(j);
                r.calculateXYCentroid(xyCen, w, h);
                double[] coeffs = r.calcParamTransCoeff();
                EllipseHelper eh = new EllipseHelper(xyCen[0], xyCen[1], coeffs);
                hs0.add(eh);
            }

            List<Region> regions1 = filtered.get(1);
            List<EllipseHelper> hs1 = new ArrayList<EllipseHelper>();
            for (int j = 0; j < regions1.size(); ++j) {
                Region r = regions1.get(j);
                r.calculateXYCentroid(xyCen, w, h);
                double[] coeffs = r.calcParamTransCoeff();
                EllipseHelper eh = new EllipseHelper(xyCen[0], xyCen[1], coeffs);
                hs1.add(eh);
            }

            for (int i = 0; i < regions0.size(); ++i) {
                if (skip0.contains(i)) { continue;}

                EllipseHelper eh0 = hs0.get(i);

                TIntSet intersects = new TIntHashSet();

                for (int j = 0; j < regions1.size(); ++j) {
                    if (skip1.contains(j)) { continue;}
                    EllipseHelper eh1 = hs1.get(j);
                    if (eh0.intersects(eh1)) {
                        intersects.add(j);
                    }
                }
                
                if (!intersects.isEmpty()) {
                    // compare eh0 and all in intersects and the largest is kept
                    // while the others are put into skip lists.
                    double maxArea = eh0.getMajorTimesMinor();
                    int maxIdx = -1;
                    TIntIterator iter = intersects.iterator();
                    while (iter.hasNext()) {
                        int j = iter.next();
                        double area = hs1.get(j).getMajorTimesMinor();
                        if (area > maxArea) {
                            maxArea = area;
                            maxIdx = j;
                        }
                    }
                    if (maxIdx == -1) {
                        // current i is largest, so put other in skip set
                        iter = intersects.iterator();
                        while (iter.hasNext()) {
                            skip1.add(iter.next());
                        }
                    } else {
                        // a j is largest, so put i in skip set and all other js
                        skip0.add(i);
                        iter = intersects.iterator();
                        while (iter.hasNext()) {
                            int j = iter.next();
                            if (j == maxIdx) { continue; }
                            skip1.add(j);
                        }
                    }
                }
            }

            // -- write anything not in skip sets to filtered
            filteredRegions = new ArrayList<Region>();
            TIntSet kept0 = new TIntHashSet();
            TIntSet kept1 = new TIntHashSet();
            List<EllipseHelper> keptEHs = new ArrayList<EllipseHelper>();
            for (int i = 0; i < regions0.size(); ++i) {
                if (skip0.contains(i)) { continue;}
                filteredRegions.add(regions0.get(i));
                kept0.add(i);
                keptEHs.add(hs0.get(i));
            }
            for (int i = 0; i < regions1.size(); ++i) {
                if (skip1.contains(i)) { continue;}
                filteredRegions.add(regions1.get(i));
                kept1.add(i);
                keptEHs.add(hs1.get(i));
            }

            // revisit all of filtered regions0 and regions1
            // to pick up any regions not yet
            // included and not intersecting with any in filteredPtRegions
            TIntSet add0 = new TIntHashSet();
            TIntSet add1 = new TIntHashSet();
            for (int i = 0; i < regions0.size(); ++i) {
                if (kept0.contains(i)) { continue;}
                EllipseHelper eh0 = hs0.get(i);
                boolean intersects = false;

                for (int j = 0; j < filteredRegions.size(); ++j) {
                    EllipseHelper eh2 = keptEHs.get(j);
                    if (eh0.intersects(eh2)) {
                        intersects = true;
                        break;
                    }
                }
                if (!intersects) {
                    add0.add(i);
                }
            }
            for (int i = 0; i < regions1.size(); ++i) {
                if (kept1.contains(i)) { continue;}
                EllipseHelper eh1 = hs1.get(i);
                boolean intersects = false;

                for (int j = 0; j < filteredRegions.size(); ++j) {
                    EllipseHelper eh2 = keptEHs.get(j);
                    if (eh1.intersects(eh2)) {
                        intersects = true;
                        break;
                    }
                }
                if (!intersects) {
                    add1.add(i);
                }
            }
            TIntIterator iter = add0.iterator();
            while (iter.hasNext()) {
                int i = iter.next();
                filteredRegions.add(regions0.get(i));
                keptEHs.add(hs0.get(i));
            }
            iter = add1.iterator();
            while (iter.hasNext()) {
                int i = iter.next();
                filteredRegions.add(regions1.get(i));
                keptEHs.add(hs1.get(i));
            }

            if (debug) {

                Image imCp;

                imCp = ptImg.copyToColorGreyscale();
                int n = filteredRegions.size();
                for (int i = 0; i < n; ++i) {
                    Region r = filteredRegions.get(i);
                    int[] clr = ImageIOHelper.getNextRGB(i);
                    r.drawEllipse(imCp, 0, clr[0], clr[1], clr[2]);
                    r.calculateXYCentroid(xyCen, imCp.getWidth(), imCp.getHeight());
                    ImageIOHelper.addPointToImage(xyCen[0], xyCen[1], imCp,
                        1, 255, 0, 0);
                    //System.out.println(type + " xy=" + xyCen[0] + "," + xyCen[1]
                    //    + " variation=" + r.getVariation());
                }
                MiscDebug.writeImage(imCp, "_" + ts + "_regions_pt_filtered_");
            }

            // ---- now add gs regions in that do not intersect with filtered
            List<List<Region>> filteredGS = filterOverlapping(gsRegions, w, h);

            add0 = new TIntHashSet();
            add1 = new TIntHashSet();

            regions0 = filteredGS.get(0);
            regions1 = filteredGS.get(1);
            
            hs0 = new ArrayList<EllipseHelper>();
            for (int j = 0; j < regions0.size(); ++j) {
                Region r = regions0.get(j);
                r.calculateXYCentroid(xyCen, w, h);
                double[] coeffs = r.calcParamTransCoeff();
                EllipseHelper eh = new EllipseHelper(xyCen[0], xyCen[1], coeffs);
                hs0.add(eh);
            }
            hs1 = new ArrayList<EllipseHelper>();
            for (int j = 0; j < regions1.size(); ++j) {
                Region r = regions1.get(j);
                r.calculateXYCentroid(xyCen, w, h);
                double[] coeffs = r.calcParamTransCoeff();
                EllipseHelper eh = new EllipseHelper(xyCen[0], xyCen[1], coeffs);
                hs1.add(eh);
            }

            for (int i = 0; i < regions0.size(); ++i) {
                EllipseHelper eh0 = hs0.get(i);
                boolean intersects = false;
                
                for (int j = 0; j < filteredRegions.size(); ++j) {
                    EllipseHelper eh2 = keptEHs.get(j);
                    if (eh0.intersects(eh2)) {
                        intersects = true;
                        break;
                    }
                }
                if (!intersects) {
                    add0.add(i);
                }
            }
            for (int i = 0; i < regions1.size(); ++i) {
                EllipseHelper eh1 = hs1.get(i);
                boolean intersects = false;
                
                for (int j = 0; j < filteredRegions.size(); ++j) {
                    EllipseHelper eh2 = keptEHs.get(j);
                    if (eh1.intersects(eh2)) {
                        intersects = true;
                        break;
                    }
                }
                if (!intersects) {
                    add1.add(i);
                }
            }
            iter = add0.iterator();
            while (iter.hasNext()) {
                int i = iter.next();
                filteredRegions.add(regions0.get(i));
            }
            iter = add1.iterator();
            while (iter.hasNext()) {
                int i = iter.next();
                filteredRegions.add(regions1.get(i));
            }
            
            if (debug) {

                Image imCp;

                imCp = ptImg.copyToColorGreyscale();
                int n = filteredRegions.size();
                for (int i = 0; i < n; ++i) {
                    Region r = filteredRegions.get(i);
                    int[] clr = ImageIOHelper.getNextRGB(i);
                    r.drawEllipse(imCp, 0, clr[0], clr[1], clr[2]);
                    r.calculateXYCentroid(xyCen, imCp.getWidth(), imCp.getHeight());
                    ImageIOHelper.addPointToImage(xyCen[0], xyCen[1], imCp,
                        1, 255, 0, 0);
                    //System.out.println(type + " xy=" + xyCen[0] + "," + xyCen[1]
                    //    + " variation=" + r.getVariation());
                }
                MiscDebug.writeImage(imCp, "_" + ts + "_regions_pt_filtered2_");
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
    
    /**
     * this extracts the filtered regions and creates edges for them then
     * substitutes them into that area in edgeList.
     */
    private void useFilterOnEdges() {
        
        //INITIALIZED, REGIONS_EXTRACTED, MERGED, EDGES_EXTRACTED
        if (!state.equals(STATE.EDGES_EXTRACTED)) {
            throw new IllegalStateException("edges must be extrcted before"
                + " this method can be used");
        }
        
        Set<PairInt> edges = new HashSet<PairInt>();
        for (int edgeIdx = 0; edgeIdx < edgeList.size(); ++edgeIdx) {
            edges.addAll(edgeList.get(edgeIdx));
        }
        
        PerimeterFinder2 finder = new PerimeterFinder2();
        
        int w = ptImg.getWidth();
        int h = ptImg.getHeight();
        
        Set<PairInt> addEdgePoints = new HashSet<PairInt>();
        for (int i = 0; i < filteredRegions.size(); ++i) {
            Region r = filteredRegions.get(i);
            Set<PairInt> ellipsePoints = new HashSet<PairInt>();
            for (int j = 0; j < r.accX.size(); ++j) {
                ellipsePoints.add(new PairInt(r.accX.get(j), r.accY.get(j)));
            }
            edges.removeAll(ellipsePoints);
            
            Set<PairInt> embedded = new HashSet<PairInt>();
            Set<PairInt> outerBoundary = new HashSet<PairInt>();
            finder.extractBorder2(ellipsePoints, embedded, outerBoundary);
            addEdgePoints.addAll(outerBoundary);
        }
        
        edges.addAll(addEdgePoints);
        
        edgeList.clear();
        
        DFSConnectedGroupsFinder finder2 = new DFSConnectedGroupsFinder();
        finder2.setMinimumNumberInCluster(1);
        finder2.findConnectedPointGroups(edges);
        
        for (int i = 0; i < finder2.getNumberOfGroups(); ++i) {
            edgeList.add(finder2.getXY(i));
        }
    }

    private void extractBoundaries() {
        
        //INITIALIZED, REGIONS_EXTRACTED, MERGED, EDGES_EXTRACTED
        if (state.equals(STATE.EDGES_EXTRACTED)) {
            throw new IllegalStateException("can only perform extraction of "
                + "edges once");
        }
        
        // without a component tree for now
        edgeList = new ArrayList<Set<PairInt>>();
        
        PerimeterFinder2 finder = new PerimeterFinder2();
    
        for (int rListIdx = 0; rListIdx < regions.size(); ++rListIdx) {
            Region r = regions.get(rListIdx);
            //NOTE: accumulated points are a connected group of points
            Set<PairInt> points = new HashSet<PairInt>();
            for (int j = 0; j < r.accX.size(); ++j) {
                int x = r.accX.get(j);
                int y = r.accY.get(j);
                PairInt p2 = new PairInt(x, y);
                points.add(p2);
                // NOTE: here is where could test membership of point in
                //   a disjoint forrest for the region connections
            }
            Set<PairInt> embedded = new HashSet<PairInt>();
            Set<PairInt> outerBorder = new HashSet<PairInt>();
            finder.extractBorder2(points, embedded, outerBorder);

            edgeList.add(outerBorder);

            DFSConnectedGroupsFinder dfsFinder = new DFSConnectedGroupsFinder();
            dfsFinder.setMinimumNumberInCluster(24);
            dfsFinder.findConnectedPointGroups(embedded);
            for (int j = 0; j < dfsFinder.getNumberOfGroups(); ++j) {

                Set<PairInt> eSet = dfsFinder.getXY(j);

                Set<PairInt> embedded2 = new HashSet<PairInt>();
                Set<PairInt> outerBorder2 = new HashSet<PairInt>();
                finder.extractBorder2(eSet, embedded2, outerBorder2);

                if (outerBorder2.size() > 16) {
                    edgeList.add(outerBorder2);
                }
            }
        }
        
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
            
            {//DEBUG
                List<Region> tmp = new ArrayList<Region>();
                for (int k = 0; k < concList.size(); ++k) {
                    tmp.add(rs.get(concList.get(k)));
                }
                _debugOrigRegions(tmp, "_rs_" + concentric.size());
            }
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
     * @return the filteredRegions
     */
    public List<Region> getFilteredRegions() {
        return filteredRegions;
    }

    /**
     * @return the useLowerContrastLimits
     */
    public boolean isUsingLowerContrastLimits() {
        return useLowerContrastLimits;
    }
}
