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
          lt 0.17 gt 0.5 
     (5) could make an assumption about the orientation of the camera place, that is, decreasing y
          pixel coord is direction "up" where sky is found.
     (6) polarization can help find boundaries between sky and foreground
         objects or horizon.
     (7) multiple images at same location and pose can help to distinguish
         between moving sky such as clouds and non-sky.
     
     * @return 
     */
    public List<SkyObject> findSky() {
        
        /*
        this may change.
        looking for sky colored labeled segmentation cells that are on the border
        of the image.
          - from mserEdges, can see that the polar theta[1] mser regions are usually
            a good way to find the majority of sky and when not, the
            polar theta[0] regions are.
          - that findins is used to decide between two color filters which are
            basically bright blue sky and dark red sky.
          - then, the regions are sorted by a color criteria to find what is hopefully
            sky in a border cell.
          - the color filter learned in the previous steps is then used to filter
            all of the labeled segmentation.
          - then the starter cell and the filtered regions are passed
            to a method to further find the sky or decide its already found 
            completely and return the sky points.
            (logic such as near constancy in color with a gradual change in
            magnitude of luv's u and v for example will be used).
        */
        
        // looking at labeled regions that
        //  have constancy in color but may have a gradient in illumination.
        //  blue - can use the polar theta images or mser edges (ptRegions[1]) 
        //  then the lch, c image, which is the magnitude of the
        //   LUV u and v, can find the gradual change if any.
        //  dark red - same except ptRegions[0]
        //mserEdges._debugOrigRegions(3, "_pt_1_acc_");
        
        List<Set<PairInt>> labeledSets = mserEdges.getLabeledSets();
        
        int[] pointLabels = createPixLabelMap(labeledSets);
        
        int sunLabel = -1;
        // remove sun points
        SkyObject sunObj = findSun();
        if (sunObj != null && !sunObj.points.isEmpty()) {
           
            System.out.println("found sun for " + debugLabel + " at " + 
                Arrays.toString(sunObj.xyCenter));
            
            for (PairInt p : sunObj.points) {
                int pixIdx = img.getInternalIndex(p);
                int label = pointLabels[pixIdx];
                labeledSets.get(label).remove(p);
                sunLabel = label;
            }
        }
        // remove rainbow points
        List<SkyObject> rainbowObjs = findRainbows();
        if (rainbowObjs != null) {
            for (SkyObject obj : rainbowObjs) {
                for (PairInt p : obj.points) {
                    int pixIdx = img.getInternalIndex(p);
                    int label = pointLabels[pixIdx];
                    labeledSets.get(label).remove(p);
                }
            }
        }
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        if (lma == null) {
            lma = imageProcessor.createLCHForLUV(img);
        }
        
        List<Region> regionPt = null;
        // note, idxs not same as other lists
        List<Set<PairInt>> labeledInRegion = null;
        // gather hsv, and lch of each labeled region
        // 6 X n
        float[][] hsvlch = null;
        
        List<Set<PairInt>> filteredLabeledRegion = new
            ArrayList<Set<PairInt>>();
        List<OneDFloatArray> filteredHSVLCH = null; 
        
        int nIter = 0;
        while (filteredLabeledRegion.isEmpty() && nIter < 2) {
            nIter++;
                        
            filteredHSVLCH = new ArrayList<OneDFloatArray>();
            
            if (nIter == 1) {
                regionPt = mserEdges.getOrigGsPtRegions().get(3);
            } else {
                regionPt = mserEdges.getOrigGsPtRegions().get(2);
            }
            
            // note, idxs not same as other lists
            labeledInRegion = createRegionLabeledSubSets(regionPt, pointLabels);

            // gather hsv, and lch of each labeled region
            // 6 X n
            hsvlch = calcHSVLCH(labeledInRegion);

            // filter labeledInRegion1 for sky colors over a generous range in hue
            // NOTE: any changes here should be made in addAdjacentSky too
            // key=index in labeledInRegion, value=set of pairint of the region
            MiscellaneousCurveHelper ch = new MiscellaneousCurveHelper();
            for (int i = 0; i < labeledInRegion.size(); ++i) {
                Set<PairInt> set = labeledInRegion.get(i);
                int[] xyCen = ch.calculateRoundedXYCentroids(set);
                float[] clrs = hsvlch[i];
                if (isInSkyRange(clrs, (nIter == 1))) {
                    if (setDoesBorderImage(set, img)) {
                        System.out.println(debugLabel
                            + " xy=" + Arrays.toString(xyCen)
                            + " hsvlch=" + Arrays.toString(clrs) + " n=" 
                            + set.size());
                        filteredLabeledRegion.add(set);
                        filteredHSVLCH.add(new OneDFloatArray(clrs));
                    }
                }
            }
        }
        
        if (filteredLabeledRegion.isEmpty()) {
            //sky wasn't found.  is the image a tunnel, that is no sky on border
            // of image?
            return null;
        }
         
        float[] avgClrs = new float[filteredHSVLCH.get(0).a.length];
        float[] stdDvClrs = new float[avgClrs.length];
        calcMeanAndStdv(filteredHSVLCH, avgClrs, stdDvClrs);
        
        System.out.println(debugLabel + " mean and stdv sky=" +
            Arrays.toString(avgClrs) + ", " + Arrays.toString(stdDvClrs));
        
        sortForBlue(filteredLabeledRegion, filteredHSVLCH, avgClrs, stdDvClrs);
        
        int startIdx = -1;
        if (sunLabel > -1) {
            //if sunLabel is in a filtered region, select that
            for (int i = 0; i < filteredLabeledRegion.size(); ++i) {
                for (PairInt p : filteredLabeledRegion.get(i)) {
                    int pixIdx = img.getInternalIndex(p);
                    int label = pointLabels[pixIdx];
                    if (sunLabel == label) {
                        startIdx = label;
                        break;
                    }
                }
            }  
        } 
        if (startIdx == -1) {  
            TIntIntMap startIdxMap = new TIntIntHashMap();
            for (PairInt p : filteredLabeledRegion.get(0)) {
                int pixIdx = img.getInternalIndex(p);
                int label = pointLabels[pixIdx];
                if (startIdxMap.containsKey(label)) {
                    startIdxMap.put(label, startIdxMap.get(label) + 1);
                } else {
                    startIdxMap.put(label, 1);
                }
            }
            System.out.println(debugLabel + " sunLabel=" + sunLabel + " strt=" + 
                startIdx);
            assert(startIdxMap.size() == 1);
            startIdx = startIdxMap.keySet().iterator().next();
            //Set<PairInt> starterSet = labeledSets.get(startIdx);
        }
        
        /*
        now that have a starter sky patch and sky filter params
           looking for the gradual change in illumination of sky as
            constancy in lch h and gradual change in lch c.
            add those sets if found.
        */
        TIntSet add = new TIntHashSet();
        add.add(startIdx);
        
        // -- recalc colors, but for all of labeled sets ---
        // 6 X n
        hsvlch = calcHSVLCH(labeledSets);
        
        //findSimilar(labeledSets, hsvlch, startIdx, add);
        
        System.out.println(debugLabel + " starter patch clrs=" + 
            Arrays.toString(hsvlch[startIdx]));
       
        addIfGradientIllumination(pointLabels, labeledSets, hsvlch, startIdx, 
            add, avgClrs, stdDvClrs, false);
        
        // put the "add" sets into output sky0Points
        Set<PairInt> sky0Points = new HashSet<PairInt>();
        TIntIterator iter = add.iterator();
        while (iter.hasNext()) {
            sky0Points.addAll(labeledSets.get(iter.next()));
        }
        
        // start a set for the case where the sun is present but not adjacent
        //    to sky0Points
        Set<PairInt> sky1Points = new HashSet<PairInt>();
        
        //TODO:
        //  if sun was found and is not connected with current found sky, 
        //     search adjacent to it and add results to a wnd list
        //     for the return
        if (sunObj != null && !sunObj.points.isEmpty()) {
            
            TIntSet cSet = findAdjacentToSun(sunObj.points, pointLabels, 
                img.getWidth(), img.getHeight());
            
            boolean isAdjToSky = false;
            iter = cSet.iterator();
            while (iter.hasNext()) {
                for (PairInt p : labeledSets.get(iter.next())) {
                    if (sky0Points.contains(p)) {
                        isAdjToSky = true;
                        break;
                    }
                }
            }
            
            if (!isAdjToSky) {
                addAdjacentSky(sunObj.points, sky0Points, pointLabels, 
                    labeledSets, hsvlch, avgClrs, stdDvClrs, cSet, sky1Points,
                    (nIter == 1));
            }
        }
        
        {// debug
            Image tmpImg = img.copyImage();
            // plotting all filtered segments in green
            // then all edges in red
            // then all added sky cells thus far in purple
            for (Set<PairInt> set : labeledSets) {
                ImageIOHelper.addCurveToImage(set, tmpImg, 0, 0, 255, 0);
            }
            for (Set<PairInt> set : mserEdges.getEdges()) {
                ImageIOHelper.addCurveToImage(set, tmpImg, 0, 255, 0, 0);
            }
            ImageIOHelper.addCurveToImage(sky0Points, tmpImg, 0, 255, 0, 255);
            if (!sky1Points.isEmpty()) {
                ImageIOHelper.addCurveToImage(sky1Points, tmpImg, 0, 255, 0, 255);
            }
            
            MiscDebug.writeImage(tmpImg, "_" + debugLabel + "_CAND_SKY_");
        }
             
        //List<Region> regionPt0 = mserEdges.getOrigGsPtRegions().get(2);
        
        MiscellaneousCurveHelper ch = new MiscellaneousCurveHelper();
        double[] xyCenter = ch.calculateXYCentroids(sky0Points);
        int x = (int)Math.round(xyCenter[0]);
        int y = (int)Math.round(xyCenter[1]);

        List<SkyObject> output = new ArrayList<SkyObject>();
        SkyObject obj = new SkyObject();
        obj.points = sky0Points;
        obj.xyCenter = new int[]{x, y};
        output.add(obj);
        
        if (!sky1Points.isEmpty()) {
            xyCenter = ch.calculateXYCentroids(sky1Points);
            x = (int)Math.round(xyCenter[0]);
            y = (int)Math.round(xyCenter[1]);
            obj = new SkyObject();
            obj.points = sky1Points;
            obj.xyCenter = new int[]{x, y};
            output.add(obj);
        }
        
        return output;
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

    private int[] createPixLabelMap(List<Set<PairInt>> labeledSets) {
        
        int[] labels = new int[img.getNPixels()];
        int count = 0;
        for (int i = 0; i < labeledSets.size(); ++i) {
            Set<PairInt> set = labeledSets.get(i);
            for (PairInt p : set) {
                int pixIdx = img.getInternalIndex(p);
                labels[pixIdx] = i;
                count++;
            }
        }
        System.out.println("count=" + count + " n=" + img.getNPixels());
        assert(count == img.getNPixels());
        
        return labels;
    }

    private List<Set<PairInt>> createRegionLabeledSubSets(
        List<Region> regionPt, int[] pointLabels) {
        
        TIntObjectMap<Set<PairInt>> map = new TIntObjectHashMap<Set<PairInt>>();
        
        for (int i = 0; i < regionPt.size(); ++i) {
            Region r = regionPt.get(i);
            for (int j = 0; j < r.accX.size(); ++j) {
                int x = r.accX.get(j);
                int y = r.accY.get(j);
                int pixIdx = img.getInternalIndex(x, y);
                int label = pointLabels[pixIdx];
                
                Set<PairInt> set = map.get(label);
                if (set == null) {
                    set = new HashSet<PairInt>();
                    map.put(label, set);
                }
                set.add(new PairInt(x, y));
            }
        }
        
        List<Set<PairInt>> out = new ArrayList<Set<PairInt>>();
        TIntObjectIterator<Set<PairInt>> iter = map.iterator();
        for (int j = 0; j < map.size(); ++j) {
            iter.advance();
            Set<PairInt> set = iter.value();
            out.add(set);
        }
        
        return out;
    }

    private float[][] calcHSVLCH(List<Set<PairInt>> regionPointsLabeled) {

        int n = regionPointsLabeled.size();
        
        float[] hsv = new float[3];
        
        float[][] hsvlch = new float[n][];
        for (int i = 0; i < n; ++i) {
            hsvlch[i] = new float[6];
            
            Set<PairInt> set = regionPointsLabeled.get(i);
            for (PairInt p : set) {
                int pixIdx = img.getInternalIndex(p);
                Color.RGBtoHSB(img.getR(pixIdx), img.getG(pixIdx), 
                    img.getB(pixIdx), hsv);
                for (int j = 0; j < 3; ++j) {
                    hsvlch[i][j] += hsv[j];
                }
                hsvlch[i][3] += lma[0].getValue(pixIdx);
                hsvlch[i][4] += lma[1].getValue(pixIdx);
                hsvlch[i][5] += lma[2].getValue(pixIdx);
            }
            float nf = set.size();
            for (int j = 0; j < 6; ++j) {
                hsvlch[i][j] /= nf;
            }
        }
        
        return hsvlch;
    }

    private boolean setDoesBorderImage(Set<PairInt> set, ImageExt img) {
        int lc = img.getWidth() - 1;
        int lr = img.getHeight() - 1;
        for (PairInt p : set) {
            int x = p.getX();
            int y = p.getY();
            if (x == 0 || y == 0 || x == lc || y == lr) {
                return true;
            }
        }
        return false;
    }

    private void sortForBlue(List<Set<PairInt>> filteredLabeledRegion,
        List<OneDFloatArray> filteredHSVLCH, float[] avgClrs,
        float[] stdDvsClrs) {
        
        float hF = 5.0f;
        float sF = 1.0f;
        float vF = 0.5f;
        
        // whiteish from snow and clouds
        if (avgClrs[4] < 20 && stdDvsClrs[4] < 20 && stdDvsClrs[5] < 20) {
            hF = 1.0f;
            sF = 0.5f;
            vF = 1.0f;
        } else if (avgClrs[0] < .2 && stdDvsClrs[0] < 0.5 && stdDvsClrs[5] < 20) {
            // reddish skies
            hF = 1.0f;
            sF = 0.5f;
            vF = 1.0f;
        }
        
        float maxS = Float.MIN_VALUE;
        float maxH = Float.MIN_VALUE;
        float maxV = Float.MIN_VALUE;
        
        for (int i = 0; i < filteredHSVLCH.size(); ++i) {
            
            OneDFloatArray a = filteredHSVLCH.get(i);
            
            if (a.a[0] > maxH) {
                maxH = a.a[0];
            }
            if (a.a[1] > maxS) {
                maxS = a.a[1];
            }
            if (a.a[2] > maxV) {
                maxV = a.a[2];
            }
        }
                
        // trying an objective function of 
        // obj cost = 5*(h-hmax) + (s-smax) + 0.5*(v-vmax).
        float[] costs = new float[filteredHSVLCH.size()];
        for (int i = 0; i < filteredHSVLCH.size(); ++i) {
            float[] a = filteredHSVLCH.get(i).a;
            float cost = 
                hF * (Math.abs(a[0] - maxH)) +
                sF * (Math.abs(a[1] - maxS)) + 
                vF * (Math.abs(a[2] - maxV)) 
                ;
            costs[i] = cost;
        }
        
        // sort filtered lists by costs
        QuickSort.sortBy1stArg(costs, filteredHSVLCH, filteredLabeledRegion);
    }

    private void findSimilar(List<Set<PairInt>> labeledSets,
        float[][] hsvlch, int starterIdx, TIntSet add) {
        
        MiscellaneousCurveHelper ch = new MiscellaneousCurveHelper();
                
        float[] strtClrs = hsvlch[starterIdx];
        
        for (int i = 0; i < hsvlch.length; ++i) {
            if (i == starterIdx) {continue;}
            
            int[] xyCen = ch.calculateRoundedXYCentroids(labeledSets.get(i));
            
            float[] clrs = hsvlch[i];
            
            float diffH = Math.abs(clrs[0] - strtClrs[0]);
            if (diffH > 0.1) {
                //System.out.println("removed by h: " + Arrays.toString(xyCen)
                //    + " " + Arrays.toString(clrs));
                continue;
            }
            
            // h of lch is a wrap around scale from 0 to 255, inclusive,
            // so need to add a phase to check whether closer to other end.
            float h_0 = strtClrs[5];
            float h_1 = clrs[5];
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
            float diff = Math.abs(h_0 - h_1);

            if (diff > 7) {
                //System.out.println("removed by h, lch: diff=" + diff 
                //    + Arrays.toString(xyCen)
                //    + " " + " " 
                //    + Arrays.toString(clrs));
                continue;
            }
            
            float diffV = Math.abs(clrs[2] - strtClrs[2]);
            if (diffV > 0.05) {
                //System.out.println("removed by v: " + Arrays.toString(xyCen)
                //    + " " + Arrays.toString(clrs));
                continue;
            }
            
            add.add(i);
        }
    }

    private void calcMeanAndStdv(List<OneDFloatArray> filteredHSVLCH, 
        float[] avgClrs, float[] stdDvsClrs) {
        
        for (OneDFloatArray a : filteredHSVLCH) {
            float[] clrs = a.a;
            for (int j = 0; j < clrs.length; ++j) {
                avgClrs[j] += clrs[j];
            }
        }
        for (int j = 0; j < avgClrs.length; ++j) {
            avgClrs[j] /= (float)filteredHSVLCH.size();
        }
        
        for (OneDFloatArray a : filteredHSVLCH) {
            float[] clrs = a.a;
            for (int j = 0; j < clrs.length; ++j) {
                float diff = avgClrs[j] - clrs[j];
                stdDvsClrs[j] += (diff * diff);
            }
        }
        for (int j = 0; j < avgClrs.length; ++j) {
            stdDvsClrs[j] = (float)
                (Math.sqrt(stdDvsClrs[j]/(float)filteredHSVLCH.size()));
        }
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

    private TIntSet findAdjacentToSun(Set<PairInt> sunPoints, int[] pointLabels,
        int width, int height) {

        TIntSet adjLabels = new TIntHashSet();
     
        int[] dxs = Misc.dx4;
        int[] dys = Misc.dy4;
        
        for (PairInt p : sunPoints) {
            int x = p.getX();
            int y = p.getY();
            for (int k = 0; k < dxs.length; ++k) {
                int x2 = x + dxs[k];
                int y2 = y + dys[k];
                if (x2 < 0 || y2 < 0 || x2 >= width || y2 >= height) {
                    continue;
                }
                PairInt p2 = new PairInt(x2, y2);
                if (sunPoints.contains(p2)) {
                    continue;
                }
                int pixIdx2 = (y2 * width) + x2;
                adjLabels.add(pointLabels[pixIdx2]);
            }
        }
        
        return adjLabels;
    }

    private void addAdjacentSky(Set<PairInt> sunPoints, Set<PairInt> sky0Points,
        int[] pointLabels, List<Set<PairInt>> labeledSets, float[][] hsvlch, 
        float[] avgClrs, float[] stdDvClrs, TIntSet adjLabels, 
        Set<PairInt> output, boolean type1) {

        MiscellaneousCurveHelper ch = new MiscellaneousCurveHelper();
        
        TIntSet add = new TIntHashSet();
        
        // in the adjacent points, look for those that have sky colors
        TIntIterator iter = adjLabels.iterator();
        while (iter.hasNext()) {
            
            int label = iter.next();
            
            float[] clrs = hsvlch[label];
            
            if (isInSkyRange(clrs, type1)) {
                add.add(label);
            }
        }
        
        if (add.isEmpty()) {
            return;
        }
        
        int n0 = add.size();
        
        addIfGradientIllumination(pointLabels, labeledSets, hsvlch, 
            add.iterator().next(), add, avgClrs, stdDvClrs, true);
        
        int n1 = add.size();
        System.out.println("added to sky=" + (n1 - n0));
                
        iter = add.iterator();
        while (iter.hasNext()) {
            output.addAll(labeledSets.get(iter.next()));
        }
        
    }
    
    private boolean isInSkyRange(float[] clrs, boolean type1) {
        
        if (type1) {
            if ((clrs[2] > 0.8) ||
                ((clrs[2] > vLLimit1) && (clrs[0] < hULimit || 
                clrs[0] > hLLimit))) {

                return true;

                //int[] xyCen = ch.calculateRoundedXYCentroids(
                //    labeledSets.get(label));
                //System.out.println(debugLabel + 
                //    " xy=" + Arrays.toString(xyCen) + 
                //    " hsvlch=" + Arrays.toString(clrs) + " n=" 
                //    + labeledSets.get(label).size());
            }
        } else {
            if ((clrs[2] > 0.8) ||
                ((clrs[2] > vLLimit0) && (clrs[0] < hULimit || 
                clrs[0] > hLLimit))) {

                return true;

                //int[] xyCen = ch.calculateRoundedXYCentroids(
                //    labeledSets.get(label));
                //System.out.println(debugLabel + 
                //    "** xy=" + Arrays.toString(xyCen) + 
                //    " hsvlch=" + Arrays.toString(clrs) + " n=" 
                //    + labeledSets.get(label).size());
            }
        }
        return false;
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
    
    public List<SkyObject> findSkyAssumingHorizon() {
        
        /*
        essentially, using level sets for binarization of the image into sky
        and non sky.   success depends upon ablilty to evaluate whether 
        pixels are sky or non-sky.
        
        the orientation of the camera is assumed to be up in the images, that
        is, up is towards smaller y.
        
        the notes below are details gathered while looking at individual 
        test images in order to design the algorithm.
        
        -- looks as if it is necessary for some cases to find the foreground
           and horizon non-sky pixels as one of the first steps.
           -- (details in progress)
              -- since region growing from bottom of image in general looks like
                 a feasible approach, level set regions should be even better.
              -- re-did regions gs[0] similar to mser sensitivity parameters used
                 for the more sensitive gs[1]:
                 -- looking for whether one region in each provides
                    the right partitions (or more than one region if needed).
                    -- seattle 1_1__29 and 0__42.
                         could curtail overrun by using the intersection with segmentation
                    -- san jose sky: 1_1__29, 0__7, 0__2 
                             nonsky: 1_1__4, 0__11
                    -- venturi_0001 sky: 0__7
                                 nonsky: 0__16, 1_1__1, 1_1__2
                       ==> best partition is single region for non-sky 0__25
                    -- venturi_0010 sky: 
                                 nonsky:
                       ==> best partition is single region for non-sky 0__25
                    -- arches sky: 
                                 nonsky:
                       ==> best partition is single region for non-sky 0__22
                    -- arches-sun sky: 
                                 nonsky:
                       ==> best partition is single region for non-sky 0__36
                    -- stlouis arch sky: 
                                 nonsky:
                       ==> best partition is single region for non-sky 0__23
                       (note: SUN is all empty space in 0__33)
                    --  contrail sky: 
                                 nonsky:
                        ==> best partition is single region for sky 0__7
                    --  klein matterhorn sky: 
                                 nonsky: 0__35, 1_1_19
                    --  patagonia sky: 
                                 nonsky: 1_1__2, 0_18, 0_10, 1_1_23
                    --  rainier snowy field sky: 
                                 nonsky:
                           horizon partly merges with sky within this image
                           so with only this, is difficult to discerne it
                    --  brown and lowe image 1 sky: 
                                            nonsky: 
                        NOTE gs 0 and 1 have several regions which add up to
                            an incomplete partition.
                        BUT, pt1 finds sky completely w/ 1_3__16
                    --  brown and lowe image 2 sky: 
                                            nonsky:
                        NOTE gs 0 and 1 have several regions which add up to
                            an incomplete partition.
                        BUT, pt1 finds sky completely w/ 1_3__18
                    -- stinson sky: 
                                 nonsky:
                       ==> best partition is single region for sky 0__35
                    -- halfdome sky: 
                                 nonsky:
                       ==> best partition is single region for nonsky 0__2 OR 1_1__1
                    -- halfdome 2 sky: 
                                 nonsky:
                       ==> best partition is single region for nonsky 0__24 OR 1_1__6
                    -- halfdome 3 sky: 
                                 nonsky:
                       ==> best partition is single region for nonsky 0__23 
                           OR for sky 1_1__11
                    -- costa-rica   sky: 1_1__13
                                 nonsky: or 0_14
                    -- norwegian    sky: 1_1__5 + 0_5
                                 nonsky:
                    -- stonehenge sky: 
                                 nonsky:
                       ==> best partition is single region for nonsky 0__3
                            OR sky 1_1__4
                    -- new-mexico sky: 
                                 nonsky:
                       ==> best partition is single region for nonsky 0__19
                            OR sky 1_1__0
                    -- arizona sky: overruns horizon 1_1__0
                                 nonsky:
                         NOTE: can find the sun completely in the space of 0___18
                    -- sky with rainbow sky: 
                                     nonsky:
                       ==> best partition is single region for nonsky 0__34
                           OR sky 1_1__3
                    -- sky with rainbow2 sky: 
                                     nonsky:
                       ==> best partition is single region for sky 1_1__1
                    
        
        -- the partition logic when designed, should include a penalty for
           including non-sky pixels...
           might need to find the region which included most possible sky
           pixels while excluding non-sky, then choose among those
           for the best partition.
           either at that prev step or here, need to consider the
           boundary defined by extending the overlap of current region
           boundary with mser edges.
        
        
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

        -- note that the level sets approach to edges 
        results in better boundaries for test with
        blue skies and blue buildings better than other techniques tried
        so far.  in other words, the MSEREdges edges preserve the boundaries
        well, even for low constrast conditions.
        
        (1) so the mseredges are still a good first step in this algorithm.
        (2) foreground/horizon determination
        (3) use of level sets to find good partition(s)
            -- before that need a pass through pixels to characterize pixels
               as possibly sky or not possibly sky.
            -- below, notes about individual tests while designing this.
            -- presumably, a good follow up to the partition is then
               following the edges they coincide with in order to
               finish the boundaries.
            -- need to consider how sun and rainbows when present in
               image affect the process
        -- sometimes the foreground is found really well by the
           regions pt_0 
            (foreground theta is darker than sky theta)
        -- sky is found really well with pt_1 for
           - blue skies in arches, arches_sun, brown and lowe
           - blue and cloudy also in cloudy_sanjose, contrail,
           - *if rainbow had been masked out of bright skies rainbow
            image, a region would have found all of the clouds,
            though not complete sky
           - in no others is the pt_1 useful for finding all sky
         
        -- can see from individual regions in gs_1 that a region can be
           a good partition for the majority of the separation of sky and
           non-sky and the approximation is faster than the
           region growing methods.
        
        NOTE that the regions are the negative only for the final mser
        edges because the positive mser parameters were set for less sensitive
        so no regions are present in gs_0, but that should be re-examined
        for this method if it doesn't use the segmentation edges.
        
       
        */
        
        List<Region> gs0 = mserEdges._extractSensitiveGS0();
        mserEdges._debugOrigRegions(gs0, debugLabel + "_gs_0_");
        
        //mserEdges._debugOrigRegions(1, debugLabel + "_gs_1_");
        //mserEdges._debugOrigRegions(2, debugLabel + "_pt_0_");
        //mserEdges._debugOrigRegions(3, debugLabel + "_pt_1_");
        
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
