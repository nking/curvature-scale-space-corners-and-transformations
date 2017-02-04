package algorithms.imageProcessing.features.mser;

import algorithms.compGeometry.PerimeterFinder2;
import algorithms.disjointSets.DisjointSet2Helper;
import algorithms.disjointSets.DisjointSet2Node;
import algorithms.imageProcessing.CIEChromaticity;
import algorithms.imageProcessing.DFSConnectedGroupsFinder;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.features.mser.MSER.Threshold;
import algorithms.misc.Misc;
import algorithms.misc.MiscDebug;
import algorithms.util.OneDIntArray;
import algorithms.util.PairInt;
import algorithms.util.VeryLongBitString;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import java.awt.Color;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

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
    // NOTE: the edges indexes do not correspond to the regions indexes
    private List<Set<PairInt>> edgeList = null;
    
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
        
        if (debug) {
            MiscDebug.writeImage(ptImg, "_" + ts + "_polartheta");
        }
        
        extractRegions();
        
        do {
        } while (mergeRegions() > 0);
        
        extractBoundaries();
        
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
    
        regions = new ArrayList<Region>();
        
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
                } else {
                    regions.add(r);
                }
            }
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
            }
        }
        
        int[] a = MSER.readIntoArray(img);
        
        MSER mser = new MSER();        
        
        List<List<Region>> regions = mser.findRegions(a, img.getWidth(),
            img.getHeight(), delta, minArea, maxArea, maxVariation, 
            minDiversity);
        
        return regions;
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
        
        Image im = gsImg.copyToColorGreyscale();
        
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
    
    //NOT READY FOR USE.  might need another colorspace
    private int mergeRegions() {
                
        // key = pix idx, value = regions idx
        TIntIntMap pRIMap = new TIntIntHashMap();
        
        // key = region idx, value = disjoint forest node
        TIntObjectMap<DisjointSet2Node<Integer>> rMap = new 
            TIntObjectHashMap<DisjointSet2Node<Integer>>();

        // w/ indexes same as regions  
        List<OneDFloatArray> clrs = new ArrayList<OneDFloatArray>();
       
        CIEChromaticity cieC = new CIEChromaticity();
         
        DisjointSet2Helper disjointSetHelper = new DisjointSet2Helper();
        
        float[] clrSum = new float[3];
        
        for (int rIdx = 0; rIdx < regions.size(); ++rIdx) {
            
            Region r = regions.get(rIdx);
            
            Arrays.fill(clrSum, 0);
            
            for (int i = 0; i < r.accX.size(); ++i) {
            
                int pixIdx = clrImg.getInternalIndex(r.accX.get(i), r.accY.get(i));
                pRIMap.put(pixIdx, rIdx);
            
                float[] lab = cieC.rgbToCIELUV(clrImg.getR(pixIdx), 
                    clrImg.getG(pixIdx), clrImg.getB(pixIdx));
                for (int j = 0; j < clrSum.length; ++j) {
                    clrSum[j] += lab[j];
                }
            }
            
            DisjointSet2Node<Integer> node = new 
                DisjointSet2Node<Integer>(Integer.valueOf(rIdx));
            node = disjointSetHelper.makeSet(node);
            assert(disjointSetHelper.findSet(node).getMember().intValue() == rIdx);
            rMap.put(rIdx, node);
        
            float[] clr0 = new float[3];
            for (int j = 0; j < clrSum.length; ++j) {
                clr0[j] = clrSum[j] / (float)r.accX.size();
            }
        
            clrs.add(new OneDFloatArray(clr0));
        }

        // adjacency map of regions indexes
        TIntObjectMap<VeryLongBitString> rAdjMap
            = new TIntObjectHashMap<VeryLongBitString>();
        
        int[] dxs = Misc.dx4;
        int[] dys = Misc.dy4;
        int w = clrImg.getWidth();
        int h = clrImg.getHeight();
        for (int rIdx = 0; rIdx < regions.size(); ++rIdx) {
            Region r = regions.get(rIdx);            
            for (int i = 0; i < r.accX.size(); ++i) {
                int x = r.accX.get(i);
                int y = r.accY.get(i);
                for (int k = 0; k < dxs.length; ++k) {
                    int x2 = x + dxs[k];
                    int y2 = y + dys[k];
                    if (x2 < 0 || y2 < 0 || (x2 >= w) || (y2 >= h)) {
                        continue;
                    }
                    int pixIdx2 = clrImg.getInternalIndex(x2, y2);
                    if (!pRIMap.containsKey(pixIdx2)) {
                        continue;
                    }
                    int rIdx2 = pRIMap.get(pixIdx2);
                    if (rIdx2 == rIdx) {
                        continue;
                    }
                    VeryLongBitString adjIdxs = rAdjMap.get(rIdx);
                    if (adjIdxs == null) {
                        adjIdxs = new VeryLongBitString(regions.size());
                        rAdjMap.put(rIdx, adjIdxs);
                    }
                    adjIdxs.setBit(rIdx2);
                }
            }
        }        
       
        Set<PairInt> compared = new HashSet<PairInt>();
        
        float clrLimit = 3.0f;
        
        for (int rIdx = 0; rIdx < regions.size(); ++rIdx) {
            
            OneDFloatArray clr = clrs.get(rIdx);
            
            VeryLongBitString adjIdxs = rAdjMap.get(rIdx);
            if (adjIdxs == null) {
                continue;
            }
            for (int rIdx2 : adjIdxs.getSetBits()) {
                PairInt c = (rIdx < rIdx2) ? new PairInt(rIdx, rIdx2) : 
                    new PairInt(rIdx2, rIdx);
                if (compared.contains(c)) {
                    continue;
                }
                compared.add(c);
                
                OneDFloatArray clr2 = clrs.get(rIdx2);
                
                float diff0 = clr.a[0] - clr2.a[0];
                float diff1 = clr.a[1] - clr2.a[1];
                float diff2 = clr.a[2] - clr2.a[2];
                double diff = Math.abs(diff0) + Math.abs(diff1) + Math.abs(diff2);
                diff /= 3.;
               
                System.out.format("(%d,%d) (%d,%d) diff=%.3f\n", 
                    regions.get(rIdx).accX.get(0),
                    regions.get(rIdx).accY.get(0),
                    regions.get(rIdx2).accX.get(0),
                    regions.get(rIdx2).accY.get(0), (float)diff
                );
                
                if (diff > clrLimit) {
                    continue;
                }
                
                // merge rIdx2 with rIdx
                DisjointSet2Node<Integer> r2Node = rMap.get(rIdx2);
                DisjointSet2Node<Integer> rNode = rMap.get(rIdx);
                DisjointSet2Node<Integer> mergedNode = 
                    disjointSetHelper.union(rNode, r2Node);
                rMap.put(rIdx, mergedNode);
                rMap.put(rIdx2, mergedNode);
            }
        }

        for (int rIdx = 0; rIdx < regions.size(); ++rIdx) {
            DisjointSet2Node<Integer> rNode = rMap.get(rIdx);
            int repIdx = disjointSetHelper.findSet(rNode).getMember().intValue();
            if (repIdx == rIdx) {
                continue;
            }
            // merge all of rIdx content into repIdx and clear
            Region r = regions.get(rIdx); 
            Region r0 = regions.get(repIdx); 
            for (int i = 0; i < r.accX.size(); ++i) {
                int x = r.accX.get(i);
                int y = r.accY.get(i);
                r0.accumulate(x, y);
            }
            r.accX.clear();
            r.accY.clear();
        }
        
        int nMerged = regions.size();
        for (int rIdx = (regions.size() - 1); rIdx > -1; --rIdx) {
            Region r = regions.get(rIdx); 
            if (r.accX.isEmpty()) {
                regions.remove(rIdx);
            }
        }  
        nMerged -= regions.size();
        System.out.println("merged " + nMerged + " regions");
    
        return nMerged;
    }
    
    public List<Set<PairInt>> getEdges() {
        
        //INITIALIZED, REGIONS_EXTRACTED, MERGED, EDGES_EXTRACTED
        if (!state.equals(STATE.EDGES_EXTRACTED)) {
            throw new IllegalStateException("must use one of the extraction"
                + " methods first");
        }
        
        return edgeList;
    }
}
