package algorithms.imageProcessing;

import algorithms.compGeometry.PerimeterFinder2;
import algorithms.imageProcessing.features.mser.MSEREdges;
import algorithms.util.OneDIntArray;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.VeryLongBitString;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.list.TIntList;
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
    
    //public SkyObject findMoonDogs() {
    //    throw new UnsupportedOperationException("not ready for use");
    //}
    
    /**
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
     */
    
    
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
     
     */
    //public List<SkyObject> findSky() {
    //   throw new UnsupportedOperationException("not yet implemented");
    //}
    
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
    
    /**
     * NOT READY FOR USE
     * uses findSkyAssumingHorizon() to find the sky and then extracts the 
     * outer boundary ordered clockwise, and returns only the points between
     * the image boundaries.
     * 
     * Note, that for some images, more information is needed to determine
     * foreground and clouds.
     * 
     * @return 
     */
    public PairIntArray extractSkyline() {
        
        List<SkyObject> skyObjs = findSkyAssumingHorizon(true);
        if (skyObjs == null || skyObjs.isEmpty()) {
            return null;
        }
        
        PerimeterFinder2 finder = new PerimeterFinder2();
        
        Set<PairInt> sky = skyObjs.get(0).points;
    
        PairIntArray boundary = finder.extractOrderedBorder(sky);
    
        // -- find the smallest x, largest y image boundary point
        //    and the largest x, largest y image boundary point
        //    and extract the subset between them.
        
        int w = mserEdges.getClrImg().getWidth();
        int h = mserEdges.getClrImg().getHeight();
        
        //max x and smallest y
        int minX_0 = Integer.MAX_VALUE;
        int maxY_0 = Integer.MIN_VALUE;
        int idx0 = -1;
        
        int maxX_1 = Integer.MIN_VALUE;
        int maxY_1 = Integer.MIN_VALUE;
        int idx1 = -1;
        
        for (int i = 0; i < boundary.getN(); ++i) {
            int x = boundary.getX(i);
            int y = boundary.getY(i);
            
            if (x > 0 && y > 0 && x < (w - 1) && y < (h - 1)) {
                continue;
            }
            
            if (x < minX_0) {
                minX_0 = x;
                maxY_0 = y;
                idx0 = i;
            } else if (x == minX_0 && y > maxY_0) {
                minX_0 = x;
                maxY_0 = y;
                idx0 = i;
            }
            
            if (x > maxX_1) {
                maxX_1 = x;
                maxY_1 = y;
                idx1 = i;
            } else if (x == maxX_1 && y > maxY_1) {
                maxX_1 = x;
                maxY_1 = y;
                idx1 = i;
            }
        }
        if (idx0 == -1 || idx1 == -1) {
            return boundary;
        }
        
        // TODO: edit here
                
        // trim from idx0 to idx1, so shift idx0 to front to make it easier
        if (idx0 < idx1) {
            
            boundary.rotateLeft(idx0);

            int end = idx1 - idx0;

            PairIntArray skyline = boundary.copyRange(0, end);

            return skyline;
            
        } else {
            
            System.out.println("idx0=" + idx0 + " idx1=" + idx1 
                 + " n=" +  boundary.getN());
            
            return boundary;
        }
    }
    
    /**
     * NOTE: improvements in segmentation may improve this method for sky and
     * filtering out foreground in the future.
     * 
     * @return 
     */
    public List<SkyObject> findSkyAssumingHorizon() {
        return findSkyAssumingHorizon(false);
    }
    
    /**
     * NOTE: improvements in segmentation may improve this method for sky and
     * filtering out foreground in the future.
     * 
     * @return 
     */
    public List<SkyObject> findSkyAssumingHorizon(boolean addOnlyContiguous) {
     
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
        
        TIntIntMap topAvgGrey = getLabelSetLs(gsImg, labeledSets, topSkyIndexes);
        
        if (debug) {
            System.out.println(debugLabel);
            TIntIterator iter = topSkyIndexes.iterator();
            while (iter.hasNext()) {
                int idx = iter.next();
                int[] hist = ptCHs.get(idx).a;
                System.out.println("  top " + Arrays.toString(hist) + ""
                    + "  inten=" + topAvgGrey.get(idx)
                    + "  xy=" + topXYs.get(idx));
            }
        }
        
        TIntList bottomBorderIndexes = getBottomBordering(labeledSets,
            ptImg.getWidth(), ptImg.getHeight());
        
        TIntIntMap bottomAvgGrey = getLabelSetLs(gsImg, labeledSets, 
            bottomBorderIndexes);
        
        // DEBUG xy centroids
        TIntObjectMap<PairInt> bottomXYs = calcCentroids(labeledSets, 
            bottomBorderIndexes, gsImg.getWidth());
        
        if (false && debug) {
            TIntIterator iter = bottomBorderIndexes.iterator();
            while (iter.hasNext()) {
                int idx = iter.next();
                int[] hist = ptCHs.get(idx).a;
                System.out.println("  bot " + Arrays.toString(hist) + ""
                    + "  inten=" + bottomAvgGrey.get(idx)
                    + "  xy=" + bottomXYs.get(idx));
            }
        }
        
        SkyObject sun = findSun();
        
        List<SkyObject> rbs = findRainbows();
        
        // ---- finding the main sky sets at top of image, before further logic
        //      about sky in other sets
        
        TIntList filtered = new TIntArrayList();
        
        ColorHistogram clrHist = new ColorHistogram();
        MiscellaneousCurveHelper ch = new MiscellaneousCurveHelper();
        
        boolean allAreBlue = true;
        int maxAvgIfBlue = Integer.MIN_VALUE;
        int maxAvgIdx = -1;
        
        if (topSkyIndexes.size() == 1) {
            
            int idx = topSkyIndexes.get(0);
            filtered.add(idx);
            
            if (normPTCHs.get(idx).a[0] > 0.1) {
                allAreBlue = false;
            } else {
                maxAvgIdx = idx;
                maxAvgIfBlue = topAvgGrey.get(idx);
            }
            
        } else {
            
            // could improve efficiency here:
            TIntList rmvd = new TIntArrayList(topSkyIndexes);
            //NOTE: assuming all bottom bordering sets are "nonsky"
            topSkyIndexes.removeAll(bottomBorderIndexes);
            rmvd.removeAll(topSkyIndexes);
            
            // test for all blue or 
            //    brightest having essentially no green
            //    and add to filtered, those resembling them
            TIntIterator iter = topSkyIndexes.iterator();
            while (iter.hasNext()) {
                int idx = iter.next();
                float[] norm = normPTCHs.get(idx).a;
                System.out.println("  *top " + Arrays.toString(norm) + ""
                    + "  inten=" + topAvgGrey.get(idx)
                    + "  xy=" + topXYs.get(idx));
                if (norm[0] > 0.1) {
                    allAreBlue = false;
                } else if (norm[3] < 0.03 && norm[0] < 0.1) {
                    // to try to remove relfection from water, avoiding the
                    //    histograms with green and darker for blue skies
                    int avg = topAvgGrey.get(idx);
                    if (avg > maxAvgIfBlue) {
                        
                        // prefer those without any purple if exist
                        boolean doNotSet = maxAvgIdx > -1 && 
                            normPTCHs.get(maxAvgIdx).a[5] == 0.0 &&
                            norm[5] > 0;
                        
                        if (!doNotSet) {
                            maxAvgIfBlue = avg;
                            maxAvgIdx = idx;
                        }
                        
                    } else if (norm[5] == 0.0 && maxAvgIdx > -1 && 
                        normPTCHs.get(maxAvgIdx).a[5] > 0.0) {
                        maxAvgIfBlue = avg;
                        maxAvgIdx = idx;
                    }
                }
            }
            
            if (false && debug) {
                iter = rmvd.iterator();
                while (iter.hasNext()) {
                    int idx = iter.next();
                    float[] normHost = normPTCHs.get(idx).a;
                    System.out.println("  RMVD " + Arrays.toString(normHost) + ""
                        + "  inten=" + topAvgGrey.get(idx)
                        + "  xy=" + topXYs.get(idx));
                }
            
                System.out.println("allAreBlue=" + allAreBlue
                    + " maxAvgIdx=" + maxAvgIdx + 
                    " maxAvgIfBlue=" + maxAvgIfBlue);
            }
            
            boolean blueSkies = (sun != null) && (maxAvgIfBlue > 0);
            
            System.out.println("blueSkies=" + blueSkies);
            
            if (allAreBlue || blueSkies) {

                if (maxAvgIdx == -1) {
                    // possibly an error in algorithm above, espec. 
                    //    regarding reflection filter for green
                    return null;
                }
                
                if (addOnlyContiguous) {

                    filtered.add(maxAvgIdx);
                    
                } else {
                    
                    //TODO: consider how to correct for sky reflected in water
                    //  (see test image for stinson beach)
                    iter = topSkyIndexes.iterator();
                    while (iter.hasNext()) {
                        int idx = iter.next();
                        float[] norm = normPTCHs.get(idx).a;
                        // collect the sky w/o green that has similar intensity
                        if (norm[3] < 0.03 && norm[0] < 0.1) {
                            int avg = topAvgGrey.get(idx);
                            filtered.add(idx);
                        }
                    }
                }
            } else {
                
                //TODO: this section needs to be corrected
               
                if (addOnlyContiguous) {
                    // choose brightest of blue or red if doesn't resemble
                    //     foreground
                    if (maxAvgIdx > -1) {
                        
                        filtered.add(maxAvgIdx);
                        
                    } else {
                        
                        /*
                        ptImg values for histogram bins:
                         0:  red = 0 - 18
                         1:  orange = 18 - 40
                         2:  yellow = 41 - 60ish
                         3:  green = 61 - 106
                         4:  blue = 107 - 192
                         5:  purple = 193 - 255
                         */
                        int maxAvgRedIdx = -1;
                        int maxAvgRed = Integer.MIN_VALUE;
                        iter = topSkyIndexes.iterator();
                        while (iter.hasNext()) {
                            int idx = iter.next();
                            int[] hist = ptCHs.get(idx).a;
                            boolean keep = true;
                            for (int i = 0; i < rmvd.size(); ++i) {
                                int rIdx = rmvd.get(i);
                                int[] rHist = ptCHs.get(rIdx).a;
                                float intersection = clrHist.intersection(hist, rHist);
                                System.out.println(" inter=" + intersection + " "
                                    + Arrays.toString(hist));
                                if (intersection > 0.5) {
                                    keep = false;
                                    break;
                                }
                            }
                            if (keep) {
                                float[] norm = normPTCHs.get(idx).a;
                                if (norm[0] > 0.1 || norm[1] > 0.1) {
                                    int avg = topAvgGrey.get(idx);
                                    if (maxAvgRedIdx == -1) {
                                        maxAvgRedIdx = idx;
                                        maxAvgRed = avg;    
                                    } else if (maxAvgRed < avg) {
                                        maxAvgRedIdx = idx;
                                        maxAvgRed = avg;
                                    }
                                }
                            }
                        }
                        filtered.add(maxAvgRedIdx);
                    }
                    
                } else {
                
                    // filtering for all red or blue and red
                    //    (mostly, trying to remove anything resembling the
                    //    removed foreground)
                    iter = topSkyIndexes.iterator();
                    while (iter.hasNext()) {
                        int idx = iter.next();
                        int[] hist = ptCHs.get(idx).a;
                        boolean keep = true;
                        for (int i = 0; i < rmvd.size(); ++i) {
                            int rIdx = rmvd.get(i);
                            int[] rHist = ptCHs.get(rIdx).a;
                            float intersection = clrHist.intersection(hist, rHist);
                            System.out.println(" inter=" + intersection + " "
                                + Arrays.toString(hist));
                            if (intersection > 0.5) {
                                keep = false;
                                break;
                            }
                        }
                        if (keep) {
                            filtered.add(idx);
                        }
                    }
                }
            }   
        }
        
        // TODO: 
        //    looking at properties of filtered indexes and bottom indexes
        //      and their intersection with remaining sets to see if there
        //      are clear ways to continue adding sky sets        
        

        Set<PairInt> skyPoints = createPoints(labeledSets, filtered, 
            ptImg.getWidth());
        
        int[] xyCen = ch.calculateRoundedXYCentroids(skyPoints);

        List<SkyObject> sky = new ArrayList<SkyObject>();

        SkyObject obj = new SkyObject();
        obj.points = skyPoints;
        obj.xyCenter = xyCen;
        sky.add(obj);

        if (sun != null) {
            sky.add(sun);
        } else if (rbs != null && !rbs.isEmpty()) {
            sky.addAll(rbs);
        }
        return sky;
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
