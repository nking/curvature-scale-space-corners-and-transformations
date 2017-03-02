package algorithms.imageProcessing;

import algorithms.QuickSort;
import algorithms.imageProcessing.features.mser.MSEREdges;
import algorithms.imageProcessing.features.mser.Region;
import algorithms.misc.Misc;
import algorithms.misc.MiscDebug;
import algorithms.util.OneDIntArray;
import algorithms.util.PairInt;
import algorithms.util.VeryLongBitString;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.list.TIntList;
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
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

/**
 * a class with methods related to daytime sky specific tasks such as finding the sun,
 * rainbow, sky, and skyline.
 * 
 * NOT READY FOR USE.
 * 
 * @author nichole
 */
public class Sky {
    
    private final MSEREdges mserEdges;
    
    private final ImageExt img;
    
    private GreyscaleImage[] lma = null;
    
    private boolean debug = false;
    
    private RainbowFinder rFinder = null;
    
    private SunFinder sFinder = null;
    
    private final float vLLimit0 = 0.2f;
    private final float vLLimit1 = 0.35f;
    private final float hULimit = 0.20f;
    private final float hLLimitPurple = 0.25f;//0.35f;//0.49f;
    private final float hLLimit = 0.49f;
    
    private String debugLabel = "";
    
    private long ts;
    
    public Sky(ImageExt img) {
        
        this.img = img.copyToImageExt();
        
        mserEdges = new MSEREdges(this.img);
        mserEdges.setToDebug();
        mserEdges.setToLowerContrast();
        mserEdges.mergeAndExtractEdges();
        
    }
    
    public SkyObject findSun() {
        if (sFinder == null) {
            sFinder = new SunFinder();
        }
        return sFinder.findSun(mserEdges);
    }
    
    public List<SkyObject> findRainbows() {
        if (rFinder == null) {
            rFinder = new RainbowFinder();
        }
        return rFinder.findRainbows(mserEdges);
    }
    
    public SkyObject findMoonDogs() {
        throw new UnsupportedOperationException("not ready for use");
    }
    
    /**
     * NOT READY FOR USE
     * 
     * NOTE: there are several ways to search for the sky cells in the image,
     * depending upon what information is available.
     * Other methods for use with specific additional information might be
     * made in the future.  For now, this one makes some assumptions
     * about sky color and sky on the border of the image.
     
     Here is an outline in progress of ways to find the sky in the image using
     only the rgb camera image (though w/ transformations to other color spaces
     possibly):
       (1) labeled information
          - to filter out non-sky objects
          - to add context.
            - for example, if know there are buildings and reflection,
              then one knows up and down
              and hence area of image where sky should be.
            - or birds in flight and airplanes in flight or animals present in water
              the are distinctive from their fins, etc.
      (2) information about "up and down" from sensors like gyroscope and compass.
          - that can limit where to look for sky in image
      (3) finding unique atmospheric objects present in the sky.
          - if find the sun or rainbows, then one knows where part of the sky is.
          - some clouds are distinguishable from snow and mountains
      (4) filtering out non sky colors such as green
          then making an assumption of sky being on an image border.
          further search requires a look at effects of sun location, that is,
          blue or red skies in the segmentation regions...
          lt 0.17 gt 0.5 
     (5) could make an assumption about the orientation of the camera place, that is, decreasing y
          pixel coord is direction "up" where sky is found.
      
     @return 
     */
    public GreyscaleImage extractSkyMask(boolean useSun) {
    
        throw new UnsupportedOperationException("not yet implemented");
    }
    
    /**
     * NOT READY FOR USE.
     * 
     * NOTE: this method is possibly over-fitting color terms in conditional
     * logic due to the small number of test images used to approximate the 
     * variables.
     * 
     * Summary: there are several ways to search for the sky cells in the image,
     * depending upon what information is available.
     * Other methods for use with specific additional information might be
     * made in the future.  For now, this one makes some assumptions
     * about sky color and sky on the border of the image.
     
     Here is an outline in progress of ways including this method, 
     to find the sky in the image using
     only the rgb camera image (though w/ transformations to other color spaces
     possibly):
       (1) labeled information
          - to filter out non-sky objects
          - to add context.
            - for example, if know there are buildings and reflection,
              then one knows up and down
              and hence area of image where sky should be.
            - or birds in flight and airplanes in flight or animals present in water
              the are distinctive from their fins, etc.
      (2) information about "up and down" from sensors like gyroscope and compass.
          - that can limit where to look for sky in image
      (3) finding unique atmospheric objects present in the sky.
          - if find the sun or rainbows, then one knows where part of the sky is.
          - some clouds are distinguishable from snow and mountains
      (4) filtering out non sky colors such as green
          then making an assumption of sky being on an image border.
          further search requires a look at effects of sun location, that is,
          blue or red skies in the segmentation regions...
     (5) could make an assumption about the orientation of the camera place, that is, decreasing y
          pixel coord is direction "up" where sky is found.
     (6) polarization can help find boundaries between sky and foreground
         objects or horizon.
     (7) multiple images at same location and pose can help to distinguish
         between moving sky such as clouds and non-sky.
     
     * @return 
     */
    public List<SkyObject> findSky() {
       throw new UnsupportedOperationException("not yet implemented");
    }
    
    public void setToDebug(String dbgLbl) {
        debug = true;
        
        this.debugLabel = dbgLbl;
    
        ts = System.currentTimeMillis();
    }
    
    public void _printGsRegions0() {
        mserEdges._debugOrigRegions(0, "gs");
    }
    public void _printGsRegions1() {
        mserEdges._debugOrigRegions(1, "gs");
    }
    public void _printPtRegions0() {
        mserEdges._debugOrigRegions(2, "pt");
    }
    public void _printPtRegions1() {
        mserEdges._debugOrigRegions(3, "pt");
    }

    private boolean isAVector(BSObj obj, int idx2, 
        List<Set<PairInt>> labeledSets, List<OneDIntArray> xyCenters,
        int[] labels, int width, int height) {
        
        // draw a line from
        // the first cell center to this cell center and asserting
        // that the vector passes through points in the cells
        // in between.
            
        if (obj.getOrderedIdxs().size() == 1) {
            return true;
        }
        
        int idx0 = obj.getOrderedIdxs().get(0);
        int[] xyi = xyCenters.get(idx0).a;
        int[] xyf = xyCenters.get(idx2).a;
        
        List<PairInt> line = new ArrayList<PairInt>();
        BresenhamsLine bLine = new BresenhamsLine();
        bLine.createLinePoints(xyi[0], xyi[1], xyf[0], xyf[1], line);
        
        if (line.isEmpty()) {
            throw new IllegalStateException("error in creating a line "
                + " between points " + Arrays.toString(xyi) + " and " +
                Arrays.toString(xyf));
        }
        
        // NOTE: this method has a flaw if one or both endpoint centers
        //     are not in their point sets.  for example, the centroid
        //     of a "c" shape is not in the "c" points.
        //     so need to make some changes here for that.
        // trimming off the end of bline which is not in the xyf set.
        for (int i = line.size() - 1; i > -1; --i) {
            PairInt p = line.get(i);
            if (labeledSets.get(idx2).contains(p)) {
                break;
            } else {
                line.remove(i);
            }
        }
        // do same trim for xycenter at front of line
        TIntList rm = new TIntArrayList();
        for (int i = 0; i < line.size(); ++i) {
            PairInt p = line.get(i);
            if (labeledSets.get(idx0).contains(p)) {
                break;
            } else {
                rm.add(i);
            }
        }
        for (int i = (rm.size() - 1); i > -1; --i) {
            int rmIdx = rm.get(i);
            line.remove(rmIdx);
        }
        
        TIntList cp = new TIntArrayList(obj.getOrderedIdxs());
        cp.add(idx2);
        
        int cIdx = 0;
        
        TIntSet present = new TIntHashSet();
        
        for (int i = 0; i < line.size(); ++i) {
            
            if (cIdx >= cp.size()) {
                return false;
            }
            
            PairInt p = line.get(i);
            
            int cLabel = cp.get(cIdx);
            
            Set<PairInt> set = labeledSets.get(cLabel);
            
            if (set.contains(p)) {
                present.add(cLabel);
            } else if (!present.contains(cLabel)) {
    
               // int pixIdx = (p.getY() * width) + p.getX();
               // int label = labels[pixIdx];
               // if (!present.contains(label)) {                    
                    return false;
               // } else {
               //     cIdx++;
               // }
            } else {
                cIdx++;
            }
        }
        
        return true;
    }

    /**
     * find any labeledSets with points having y=0 and sky colors and return the
     * index of those sets in labeledSets.
     * @param ptCHs
     * @param labeledSets
     * @return 
     */
    private TIntList findSkyColorsAtTop(List<OneDFloatArray> normPTCHs, 
        List<TIntSet> labeledSets, int imgWidth) {
        
        /*
        ptImg values for histogram bins:
         0:  red = 0 - 18
         1:  orange = 18 - 40
         2:  yellow = 41 - 60ish
         3:  green = 61 - 106
         4:  blue = 107 - 192
         5:  purple = 193 - 255
        */
        
        TIntList indexes = new TIntArrayList();
        
        for (int i = 0; i < labeledSets.size(); ++i) {
            TIntSet pixIdxs = labeledSets.get(i);
            if (hasAPointWithY(0, pixIdxs, imgWidth)) {
                float[] normalizedHist = normPTCHs.get(i).a;
                // all colors are sky colors except green,
                // but this could be improved
                if (normalizedHist[3] < 0.3f) {
                    indexes.add(i);
                }
            }
        }
        
        return indexes;
    }

    private boolean hasAPointWithY(int yCoord, TIntSet pixIdxs, 
        int imgWidth) {
        
        TIntIterator iter = pixIdxs.iterator();
        while (iter.hasNext()) {
            int pixIdx = iter.next();
            int y = pixIdx/imgWidth;
            if (y == yCoord) {
                return true;
            } 
        }
        
        return false;
    }

    private TIntIntMap getLabelSetLs(GreyscaleImage gsImg, 
        List<TIntSet> labeledSets, TIntList indexes) {
        
        TIntIntMap map = new TIntIntHashMap();
        
        TIntIterator iter = indexes.iterator();
        while (iter.hasNext()) {
            int idx = iter.next();
            
            TIntSet pixIdxs = labeledSets.get(idx);
            
            TIntIterator iter2 = pixIdxs.iterator();
            
            int avg = 0;
            while (iter2.hasNext()) {
                int pixIdx = iter2.next();
                avg += gsImg.getValue(pixIdx);
            }
            avg /= pixIdxs.size();
            
            map.put(idx, avg);
        }
        
        return map;
    }

    private TIntObjectMap<PairInt> calcCentroids(List<TIntSet> labeledSets, 
        TIntList indexes, int width) {
        
        TIntObjectMap<PairInt> map = new TIntObjectHashMap<PairInt>();
        
        TIntIterator iter = indexes.iterator();
        while (iter.hasNext()) {
            int idx = iter.next();
            
            TIntSet pixIdxs = labeledSets.get(idx);
            
            TIntIterator iter2 = pixIdxs.iterator();
            
            int xAvg = 0;
            int yAvg = 0;
            while (iter2.hasNext()) {
                int pixIdx = iter2.next();
                int y = pixIdx/width;
                int x = pixIdx - (y * width);
                xAvg += x;
                yAvg += y;
            }
            xAvg /= pixIdxs.size();
            yAvg /= pixIdxs.size();
            
            map.put(idx, new PairInt(xAvg, yAvg));
        }
        
        return map;
    }

    private TIntList getBottomBordering(List<TIntSet> labeledSets, int width,
        int height) {

        TIntList indexes = new TIntArrayList();
        
        for (int i = 0; i < labeledSets.size(); ++i) {
            TIntSet pixIdxs = labeledSets.get(i);
            if (hasAPointWithY(height - 1, pixIdxs, width)) {
                indexes.add(i);
            }
        }
        
        return indexes;
    }

    private Set<PairInt> createPoints(List<TIntSet> labeledSets, 
        TIntList indexes, int width) {
        
        Set<PairInt> set = new HashSet<PairInt>();
        
        TIntIterator iter = indexes.iterator();
        while (iter.hasNext()) {
            int idx = iter.next();
            TIntSet pixIdxs = labeledSets.get(idx);
            TIntIterator iter2 = pixIdxs.iterator();
            while (iter2.hasNext()) {
                int pixIdx = iter2.next();
                int y = pixIdx/width;
                int x = pixIdx - (y * width);
                set.add(new PairInt(x, y));
            }
        }
        
        return set;
    }

    private List<OneDFloatArray> normalize(List<OneDIntArray> chs) {
    
        List<OneDFloatArray> norm = new ArrayList<OneDFloatArray>();
        for (int i = 0; i < chs.size(); ++i) {
            int[] a = chs.get(i).a;
            int tot = 0;
            for (int c : a) {
                tot += c;
            }
            float[] b = new float[a.length];
            for (int j = 0; j < a.length; ++j) {
                b[j] = (float)a[j]/(float)tot;
            }
            norm.add(new OneDFloatArray(b));
        }
        
        return norm;
    }
    
    private class BSObj {
        private VeryLongBitString idxs;
        private TIntList orderedIdxs;
        public BSObj(int nBS) {
            idxs = new VeryLongBitString(nBS);
            orderedIdxs = new TIntArrayList();
        }
        public BSObj() {
        }
        public void add(int idx) {
            idxs.setBit(idx);
            orderedIdxs.add(idx);
        }
        public int latest() {
            return orderedIdxs.get(orderedIdxs.size() - 1);
        }
        public VeryLongBitString getBS() {
            return idxs;
        }
        public TIntList getOrderedIdxs() {
            return orderedIdxs;
        }
        public BSObj copy() {
            BSObj cp = new BSObj();
            cp.idxs = idxs.copy();
            cp.orderedIdxs = new TIntArrayList(orderedIdxs);
            return cp;
        }
    }

    private void addIfGradientIllumination(int[] labels, 
        List<Set<PairInt>> labeledSets, float[][] hsvlch, int startIdx, 
        TIntSet add, float[] avgClrs, float[] stdDvClrs,
        boolean isNearSun) {
        
        MiscellaneousCurveHelper ch = new MiscellaneousCurveHelper();
        
        List<OneDIntArray> xyCenters = new ArrayList<OneDIntArray>();
        for (int i = 0; i < labeledSets.size(); ++i) {
            Set<PairInt> set = labeledSets.get(i);
            int[] xy = ch.calculateRoundedXYCentroids(set);
            xyCenters.add(new OneDIntArray(xy));
        }
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        TIntObjectMap<VeryLongBitString> adjMap = imageProcessor.
            createAdjacencyMap(labels, labeledSets, img.getWidth(),
            img.getHeight());
        
        // 0 = bright blue, 1 = blue, 2 = white, 3 = red or dark
        int avgClrType = 1;
        if (avgClrs[0] < 0.75 && avgClrs[0] > 0.45) {
            if (avgClrs[2] > 0.65 && avgClrs[4] > 40) {
                avgClrType = 0;
            } else if (avgClrs[2] > 0.85) {
                avgClrType = 2;
            } else if (avgClrs[0] < 0.6) {
                // less blue than type 1
                avgClrType = 5;
            } else {
                avgClrType = 1;
            }
        } else if (avgClrs[4] < 20 && stdDvClrs[4] < 20 && 
            stdDvClrs[5] < 20 && avgClrs[2] > 0.75) {
            // whiteish from snow and clouds
            avgClrType = 2;
        } else if (avgClrs[0] < 0.2 && stdDvClrs[0] < 0.5 && stdDvClrs[5] < 20) {
            avgClrType = 3;
        } else if ((avgClrs[0] < 0.1 || avgClrs[0] > 0.85) && stdDvClrs[0]> 0.1) {
            avgClrType = 3;
        } else if (avgClrs[0] < 0.25 && avgClrs[1] > 0.8) {
            avgClrType = 3;
        } else if (avgClrs[0] > 0.75) {
            // purple to red
            avgClrType = 6;
        } 
        
        // find gradients in the 4th item in clrs[] which is the
        //    magnitude lch, that is c
        /*
                 [0.-1]  [0.-1, -2] 
              0  -1         -2
        
                -1
        
                    -2
         if roughly same in h of lch,
            will store the c's of lch until queue is empty
            and for the paths which have increasing or decreasing
            or c, will add those cells into "add"
        */
       
        float eps = isNearSun ? 10 : 0;
        
        List<BSObj> paths = new ArrayList<BSObj>();
        
        int nBS = hsvlch.length;
        
        BSObj strtObj = new BSObj(nBS);
        strtObj.add(startIdx);
       
        // BFS search around items in add
        ArrayDeque<BSObj> queue = new ArrayDeque<BSObj>();
        queue.add(strtObj);
        TIntIterator iter = add.iterator();
        while (iter.hasNext()) {
            int idx = iter.next();
            if (idx == startIdx) { continue; }
            BSObj obj = new BSObj(nBS);
            obj.add(idx);
            queue.add(obj);
        }
                
        Set<VeryLongBitString> visited = new HashSet<VeryLongBitString>();
        
        while (!queue.isEmpty()) {
            
            BSObj obj = queue.pop();
            int idx = obj.latest();
            
            if (visited.contains(obj.getBS())) { continue; }
            visited.add(obj.getBS());
            
            VeryLongBitString ngbhrs = adjMap.get(idx);
            if (ngbhrs == null) {
                continue;
            }
            
            float[] clrs = hsvlch[idx];
            int[] nghbrsIdxs = ngbhrs.getSetBits();
            for (int idx2 : nghbrsIdxs) {
                
                // -- check if path was searched already
                VeryLongBitString tmpBS = obj.getBS().copy();
                tmpBS.setBit(idx2);
                if (visited.contains(tmpBS)) {
                    continue;
                }
                
                // -- skip if this addition to obj.orderedIdxs is not along
                //    a vector.
                //    temporarily, this check will be drawing a line from
                //    the first cell center to this cell center and asserting
                //    that the vector passes through points in the cells
                //    in between.
                if (!isAVector(obj, idx2, labeledSets, xyCenters, labels,
                    img.getWidth(), img.getHeight())) {
                    System.out.println("n0=" + labeledSets.get(
                        obj.getOrderedIdxs().get(0)).size());
                    StringBuilder sb = new StringBuilder("SKIP: ");
                    TIntList tmp = new TIntArrayList(obj.orderedIdxs);
                    tmp.add(idx2);
                    for (int j = 0; j < tmp.size(); ++j) {
                        int _idx = tmp.get(j);
                        int[] xyCen = ch.calculateRoundedXYCentroids(
                            labeledSets.get(_idx));
                        sb.append(Arrays.toString(xyCen)).append(", ");
                    }
                    System.out.println(sb.toString());
                    continue;
                }
                
                float[] clrs2 = hsvlch[idx2];
                
                // diff in h of lch:
                float h_0 = clrs[5];
                float h_1 = clrs2[5];
                if (h_0 > h_1) {
                    // add a phase to next value if it's closer to current with addition
                    if ((h_0 - h_1) > Math.abs(h_0 - (h_1 + 255))) {
                        h_1 += 255;
                    }
                } else if (h_1 > h_0) {
                    // add a phase to next value if it's closer to current with addition
                    if ((h_1 - h_0) > Math.abs(h_1 - (h_0 + 255))) {
                        h_0 += 255;
                    }
                }
                float diffH2 = Math.abs(h_0 - h_1) - eps;
                if (diffH2 < 0) {
                    diffH2 = 0;
                }
                
                float diffC = Math.abs(clrs[4] - clrs2[4]) - eps;
                if (diffC < 0) {
                    diffC = 0;
                }
                
                int[] xyCen = ch.calculateRoundedXYCentroids(
                    labeledSets.get(idx2));
                
                boolean added = false;
                if (avgClrType == 0) {
                    if (diffH2 < 9 && diffC < 30) {
                        added = true;
                    }
                } else if (avgClrType == 1) {
                    if (diffH2 < 5 && diffC < 5) {
                        added = true;
                    }
                } else if (avgClrType == 5) {
                    if (diffH2 < 6 && diffC < 8) {
                        added = true;
                    }
                } else if (avgClrType == 6) {
                    if (diffH2 < 30 && diffC < 20) {
                        added = true;
                    }
                } else if (avgClrType == 3) {
                    if (diffH2 < 9 && diffC < 45) {
                        added = true;
                    }
                }
                
                if (added) {
                    BSObj obj2 = obj.copy();
                    obj2.add(idx2);

                    queue.add(obj2);
                    paths.add(obj2);

                }
                
                System.out.println(debugLabel + 
                    " neighbor " + " xy=" + Arrays.toString(xyCen) + 
                    " hsvlch=" + Arrays.toString(hsvlch[idx2]) 
                    + " diffH2=" + diffH2 + " diffC=" + diffC
                    + " n=" + 
                    labeledSets.get(idx2).size() + 
                    " avgClrType=" + avgClrType + " added=" + added
                );
            }
        }
        
        eps = 1;
        if (avgClrType == 3) {
            eps = 10;
        }
        
        // paths ordered indexes are along vectors, so can now look
        //    for whether a path item has same c's, increasing c's
        //    or decreasing c's
        for (BSObj path : paths) {
            
            TIntList idxs = path.getOrderedIdxs();
            
            float[] deltaCs = new float[idxs.size() - 1];
            for (int i = 0; i < idxs.size() - 1; ++i) {
                deltaCs[i] = hsvlch[idxs.get(i)][4] 
                    - hsvlch[idxs.get(i + 1)][4];
            }
            
            {
                for (int i = 0; i < idxs.size(); ++i) {
            System.out.println(i + " ?: " + hsvlch[idxs.get(i)][4]);
                }
            }
            
            // TODO: improve this.
            // will make an assumption that all are increasing,
            // and sepearately that all are decreasing and see if either are
            // true
            boolean allIncr = true;
            boolean allDecr = true;
            for (int i = 1; i < deltaCs.length;++i) {
                // c[1] > (c[0] - eps)
                if (deltaCs[i] < (deltaCs[i - 1] + eps)) {
                    allIncr = false;
                    break;
                }
            }
            if (!allIncr) {
                for (int i = 1; i < deltaCs.length;++i) {
                    // c[1] < (c[0] - eps)
                    if (deltaCs[i] > (deltaCs[i - 1] + eps)) {
                        allDecr = false;
                        break;
                    }
                }
            }
            if (!allIncr && !allDecr) {
                System.out.println("did not add: " + Arrays.toString(deltaCs));
                continue;
            }
            // add all in path
            for (int i = 0; i < idxs.size(); ++i) {
                int idx = idxs.get(i);
                add.add(idx);
            } 
        }
        
    }
    
    /**
     * NOTE: improvements in segmentation may improve this method for sky and
     * filtering out foreground in the future.
     * 
     * @return 
     */
    public List<SkyObject> findSkyAssumingHorizon() {
        
        /*
        essentially, using level sets for binarization of the image into sky
        and non sky.   
        success depends upon ablilty to evaluate whether 
        pixels are sky or non-sky.
        
        the orientation of the camera is assumed to be up in the images, that
        is, up is towards smaller y.
        
        a bigger picture look at the process would suggest that usually the
        horizon is distinguishable from the sky by color contrast and 
        illumination.
        polarization could help distinguish atmosphere and clouds from the
        foreground at the horizon.  the sunlight source function is not
        polarized (its a mixture of all polarizations), but upon reflection
        or interactions with molecules in clouds and aerosols, the light 
        that the camera receives is more polarized.
        
        for white cloudy skies over snowy mountains, polarization might be
        helpful in differentiating, hence finding the skyline.
        for blue metal buildings upon blue sky backgrounds, polarization difference
        will be large.

        NOTE that polarization data isnt available here in the project currently, 
        nor are multiple images taken at the
        same location and pose.
        */
        
        /*
        starting with the segmentation from MSEREdges,
           -- list 1: all labeled sets touching the y=0 border of the image
                      and having colors not dominated by green.
           -- list 2: all labeled sets touching the y=hieght-1 border of the
                      image.
           -- search for sun and rainbows
        
        general rules used below:
           if list 1 contains only one set
               if no sun nor rainbows, that set is returned 
                  without any further region growing.
               else if has sun or rainbows,
                  return the top border set and the sun or rainbow.
                  it would be difficult to add other sets without
                  other information such as labeling, polarization,
                  multiple images (for cloud or foreground motion), etc.
           if list 1 contains more than one set,
               -- remove the intersection with list 2
               -- if all of the remaining list 1 are predominantly
                  blue histograms,
                     (or alternatively, find brightest 
                      with little to no green (<0.01 blue)
                      as starter and
                     return it and the others resembling its histogram).
               -- if there are only red histograms in list 1
                  return all
               -- if there are blue and red histograms in list 1
                  if all are bright, return all
                  else return brightest
               -- if sun or rainbow are present, return those too.
                  it would be difficult to add other sets without
                  other information such as labeling, polarization,
                  multiple images (for cloud or foreground motion), etc.
        */
        
        GreyscaleImage ptImg = mserEdges.getPtImg();
        
        GreyscaleImage gsImg = mserEdges.getGsImg();
        
        List<TIntSet> labeledSets = mserEdges.getLabeledSets();
       
        if (labeledSets.isEmpty()) {
            return null;
        }
        
        List<OneDIntArray> ptCHs = 
            ColorHistogram.createPTHistograms(ptImg, labeledSets);
        
        List<OneDFloatArray> normPTCHs = normalize(ptCHs);
        
        //NOTE: this may need corrections for some bright white clouds:
        TIntList topSkyIndexes = findSkyColorsAtTop(normPTCHs,
            labeledSets, ptImg.getWidth());
        
        TIntObjectMap<PairInt> topXYs = calcCentroids(labeledSets, topSkyIndexes,
            gsImg.getWidth());
        
        SkyObject sun = findSun();
        
        List<SkyObject> rbs = findRainbows();
        
        if (topSkyIndexes.size() == 1) {
            
            PairInt xy = topXYs.get(topSkyIndexes.get(0));
            
            Set<PairInt> skyPoints = createPoints(labeledSets, topSkyIndexes, 
                ptImg.getWidth());
            
            List<SkyObject> sky = new ArrayList<SkyObject>();
            
            SkyObject obj = new SkyObject();
            obj.points = skyPoints;
            obj.xyCenter = new int[]{xy.getX(), xy.getY()};
            sky.add(obj);
    
            System.out.println(debugLabel + ": hist=" + Arrays.toString(
                normPTCHs.get(0).a));  
            
            if (sun != null) {
                sky.add(sun);
            } else if (rbs != null && !rbs.isEmpty()) {
                sky.addAll(rbs);
            }
            
            return sky;
        }
        
        TIntIntMap topAvgGrey = getLabelSetLs(gsImg, labeledSets, topSkyIndexes);
        
        if (debug) {
            System.out.println(debugLabel);
            TIntIterator iter = topSkyIndexes.iterator();
            while (iter.hasNext()) {
                int idx = iter.next();
                int[] hist = ptCHs.get(idx).a;
                //System.out.println("  top " + Arrays.toString(hist) + ""
                //    + "  inten=" + topAvgGrey.get(idx)
                //    + "  xy=" + topXYs.get(idx));
            }
        }
        
        TIntList bottomBorderIndexes = getBottomBordering(labeledSets,
            ptImg.getWidth(), ptImg.getHeight());
        
        TIntIntMap bottomAvgGrey = getLabelSetLs(gsImg, labeledSets, 
            bottomBorderIndexes);
        
        // DEBUG xy centroids
        TIntObjectMap<PairInt> bottomXYs = calcCentroids(labeledSets, 
            bottomBorderIndexes, gsImg.getWidth());
        
        if (debug) {
            TIntIterator iter = bottomBorderIndexes.iterator();
            while (iter.hasNext()) {
                int idx = iter.next();
                int[] hist = ptCHs.get(idx).a;
                //System.out.println("  bot " + Arrays.toString(hist) + ""
                //    + "  inten=" + bottomAvgGrey.get(idx)
                //    + "  xy=" + bottomXYs.get(idx));
            }
        }
        
        MiscellaneousCurveHelper ch = new MiscellaneousCurveHelper();
        
        // could improve efficiency here:
        TIntList rmvd = new TIntArrayList(topSkyIndexes);
        topSkyIndexes.removeAll(bottomBorderIndexes);
        rmvd.removeAll(topSkyIndexes);
        
        // test for all blue or 
        //    brightest having essentially no green
        //    and return those resembling them
        boolean allAreBlue = true;
        int maxAvgIfBlue = Integer.MIN_VALUE;
        int maxAvgIdx = -1;
        TIntIterator iter = topSkyIndexes.iterator();
        while (iter.hasNext()) {
            int idx = iter.next();
            float[] norm = normPTCHs.get(idx).a;
            System.out.println("  *top " + Arrays.toString(norm) + ""
                + "  inten=" + topAvgGrey.get(idx)
                + "  xy=" + topXYs.get(idx));
            if (norm[0] > 0.1) {
                allAreBlue = false;
            } else if (norm[3] < 0.01 && norm[0] < 0.1) {
                // to try to remove relfection from water, avoiding the
                //    histograms with green and darker for blue skies
                int avg = topAvgGrey.get(idx);
                if (avg > maxAvgIfBlue) {
                    maxAvgIfBlue = avg;
                    maxAvgIdx = idx;
                }
            }
        }
        {
        iter = rmvd.iterator();
        while (iter.hasNext()) {
            int idx = iter.next();
            float[] normHost = normPTCHs.get(idx).a;
            System.out.println("  RMVD " + Arrays.toString(normHost) + ""
                + "  inten=" + topAvgGrey.get(idx)
                + "  xy=" + topXYs.get(idx));
        }
        }
        
        if (allAreBlue) {
            
            if (maxAvgIdx == -1) {
                // possibly an error in algorithm above, espec. 
                //    regarding reflection filter for green
                return null;
            }
            
            //TODO: consider how to correct for sky reflected in water
            //  (see stinson beach test image)
            
            TIntList filtered = new TIntArrayList();
            
            iter = topSkyIndexes.iterator();
            while (iter.hasNext()) {
                int idx = iter.next();
                float[] norm = normPTCHs.get(idx).a;
                // collect the sky w/o green that has similar intensity
                if (norm[3] < 0.01 && norm[0] < 0.1) {
                    int avg = topAvgGrey.get(idx);
                    filtered.add(idx);
                }
            }
            
            Set<PairInt> skyPoints = createPoints(labeledSets, filtered, 
                ptImg.getWidth());
            int[] xyCen = ch.calculateRoundedXYCentroids(skyPoints);
            
            List<SkyObject> sky = new ArrayList<SkyObject>();
            
            SkyObject obj = new SkyObject();
            obj.points = skyPoints;
            obj.xyCenter = xyCen;
            sky.add(obj);
            
            return sky;
        }
        
        /*
        ptImg values for histogram bins:
         0:  red = 0 - 18
         1:  orange = 18 - 40
         2:  yellow = 41 - 60ish
         3:  green = 61 - 106
         4:  blue = 107 - 192
         5:  purple = 193 - 255
        */
        
        // filtering for all red or blue and red
        //    (mostly, trying to remove anything resembling the
        //    removed foreground)
        
        /*
        if list 1 contains more than one set,
           -- remove the intersection with list 2
           -- if all of the remaining list 1 are predominantly
              blue histograms,
                 (or alternatively, find brightest 
                  with little to no green (<0.01 blue)
                  as starter and
                 return it and the others resembling its histogram).
           -- if there are only red histograms in list 1
              return all
           -- if there are blue and red histograms in list 1
              if all are bright, return all
              else return brightest
           -- if sun or rainbow are present, return those too.
              it would be difficult to add other sets without
              other information such as labeling, polarization,
              multiple images (for cloud or foreground motion), etc.
        */
        
        return null;
        //throw new UnsupportedOperationException("not yet implemented");
    }
    
    //TODO: consider adding findSolarEclipse or sun w/ occultation or coronograph...
    // note that the moon can be found with "findSun" since it is illuminated 
    // by sun light.
    
    public static class SkyObject {
        Set<PairInt> points;
        int[] xyCenter;
    }
    
    private class OneDFloatArray {
        float[] a;
        public OneDFloatArray(float[] b) {
            a = b;
        }
    }
}
