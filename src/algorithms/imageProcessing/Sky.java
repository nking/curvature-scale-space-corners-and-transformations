package algorithms.imageProcessing;

import algorithms.QuickSort;
import algorithms.imageProcessing.features.mser.MSEREdges;
import algorithms.imageProcessing.features.mser.Region;
import algorithms.misc.MiscDebug;
import algorithms.util.PairInt;
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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * a class with methods related to sky specific tasks such as finding the sun,
 * rainbow, sky, and skyline.
 * 
 * NOT READY FOR USE.
 * still in design phase
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
     * NOTE: there are several ways to search for the sky cells in the image,
     * depending upon what information is available.
     * Other methods for use with specific additional information might be
     * made in the future.  For now, this one makes some assumptions
     * about sky color and sky on the border of the image.
     
     Here is an outline in progress of ways to find the sky in the image using
     only the rgb camera image (though w/ transformations to other color spaces
     possibly):
       (1) labelled information
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
    public GreyscaleImage extractSkyMask() {
        
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
        
        // remove sun points
        SkyObject sunObj = findSun();
        if (sunObj != null && !sunObj.points.isEmpty()) {
            System.out.println("found sun for " + debugLabel + " at " + 
                Arrays.toString(sunObj.xyCenter));
            for (PairInt p : sunObj.points) {
                int pixIdx = img.getInternalIndex(p);
                int label = pointLabels[pixIdx];
                labeledSets.get(label).remove(p);
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
            // key=index in labeledInRegion, value=set of pairint of the region
            MiscellaneousCurveHelper ch = new MiscellaneousCurveHelper();
            for (int i = 0; i < labeledInRegion.size(); ++i) {
                Set<PairInt> set = labeledInRegion.get(i);
                int[] xyCen = ch.calculateRoundedXYCentroids(set);
                float[] clrs = hsvlch[i];
                if (nIter == 1) {
                    if ((clrs[2] > 0.8) ||
                        ((clrs[2] > vLLimit1) && (clrs[0] < hULimit || 
                        clrs[0] > hLLimit))) {
                        if (setDoesBorderImage(set, img)) {
                            System.out.println(debugLabel + 
                                " xy=" + Arrays.toString(xyCen) + 
                                " hsvlch=" + Arrays.toString(clrs) + " n=" + set.size());
                            filteredLabeledRegion.add(set);
                            filteredHSVLCH.add(new OneDFloatArray(clrs));
                        }
                    }
                } else {
                    if ((clrs[2] > 0.8) ||
                        ((clrs[2] > vLLimit0) && (clrs[0] < hULimit || 
                        clrs[0] > hLLimit))) {
                        if (setDoesBorderImage(set, img)) {
                            System.out.println(debugLabel + 
                                " PT0 xy=" + Arrays.toString(xyCen) + 
                                " hsvlch=" + Arrays.toString(clrs) + " n=" + set.size());
                            filteredLabeledRegion.add(set);
                            filteredHSVLCH.add(new OneDFloatArray(clrs));
                        }
                    }
                }
            }
        }
        
        if (filteredLabeledRegion.isEmpty()) {
            //sky wasn't found.  is the image a tunnel, that is no sky on border
            // of image?
            return null;
        }
                
        sortForBlue(filteredLabeledRegion, filteredHSVLCH);
        
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
        System.out.println(debugLabel);
        assert(startIdxMap.size() == 1);
        int startIdx = startIdxMap.keySet().iterator().next();
        Set<PairInt> starterSet = labeledSets.get(startIdx);

        /*
        now that have a starter sky patch and sky filter params
        (1) add labeled sets that are very similar to starter patch
        (2) look for the gradual change in illumination of sky as
            constancy in lch h and gradual change in lch c.
            add those sets if found.
        */
        TIntSet add = new TIntHashSet();
        
        // -- recalc colors, but for all of labeled sets ---
        // 6 X n
        hsvlch = calcHSVLCH(labeledSets);
        
        findSimilar(labeledSets, hsvlch, startIdx, add);
        
        System.out.println("starter patch clrs=" + Arrays.toString(hsvlch[startIdx]));
        
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
            TIntIterator iter = add.iterator();
            while (iter.hasNext()) {
                int label = iter.next();
                ImageIOHelper.addCurveToImage(labeledSets.get(label), tmpImg, 
                    0, 255, 0, 255);
            }
            MiscDebug.writeImage(tmpImg, "_" + debugLabel + "_CAND_SKY_");
        }
        
        /*  a look at sky colors
        filteredHSVLCH = new ArrayList<OneDFloatArray>();
        
        MiscellaneousCurveHelper ch = new MiscellaneousCurveHelper();
        
        TIntList rm = new TIntArrayList();
        if (nIter == 1) {
            for (int i = 0; i < labeledSets.size(); ++i) {
                String str = (i == startIdx) ? "STRT" : "";
                int[] xyCen = ch.calculateRoundedXYCentroids(labeledSets.get(i));
                float[] clrs = hsvlch[i];
                if (!((clrs[2] > vLLimit1) && (clrs[0] < hULimit || 
                    clrs[0] > hLLimit)
                    )) {
                    rm.add(i);
                    
                    System.out.println(debugLabel + 
                        " RM " + str + " xy=" + Arrays.toString(xyCen) + 
                        " hsvlch=" + Arrays.toString(clrs) + " n=" + 
                        labeledSets.get(i).size());
                } else {
                    filteredHSVLCH.add(new OneDFloatArray(clrs));
                    System.out.println(debugLabel + 
                        str + " xy=" + Arrays.toString(xyCen) + 
                        " hsvlch=" + Arrays.toString(clrs) + " n=" + 
                        labeledSets.get(i).size());
                }
            }
        } else {
            for (int i = 0; i < labeledSets.size(); ++i) {
                float[] clrs = hsvlch[i];
                if (!((clrs[2] > vLLimit0) && (clrs[0] < hULimit || 
                    clrs[0] > hLLimit))) {
                    rm.add(i);
                } else {
                    filteredHSVLCH.add(new OneDFloatArray(clrs));
                }
            }
        }
               
        for (int i = rm.size() - 1; i > -1; --i) {
            int rmIdx = rm.get(i);
            labeledSets.remove(rmIdx);
        }
        
        if (debug) {
            Image tmpImg = img.copyImage();
            for (Set<PairInt> set : labeledSets) {
                ImageIOHelper.addCurveToImage(set, tmpImg, 0, 0, 255, 0);
            }
            //ImageIOHelper.addCurveToImage(filteredLabeledRegion.get(0), 
            //    tmpImg, 0, 255, 0, 255);
            MiscDebug.writeImage(tmpImg, "_" + debugLabel + "_CAND_SKY_");
        }
        */
                
        //List<Region> regionPt0 = mserEdges.getOrigGsPtRegions().get(2);
        
        return null;
        //throw new UnsupportedOperationException("not ready for use");
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
        List<OneDFloatArray> filteredHSVLCH) {
        
        /*
        for largely white filtered sets,
             c of lch is small
        */
        float[] avgClrs = new float[filteredHSVLCH.get(0).a.length];
        float[] stdDvsClrs = new float[avgClrs.length];
        calcMeanAndStdv(filteredHSVLCH, avgClrs, stdDvsClrs);
        System.out.println(debugLabel + " mean and stdv sky=" +
            Arrays.toString(avgClrs) + ", " + Arrays.toString(stdDvsClrs));
        
        float hF = 5.0f;
        float sF = 1.0f;
        float vF = 0.5f;
        
        // whiteish from snow and clouds
        if (avgClrs[4] < 20 && stdDvsClrs[4] < 20 && stdDvsClrs[5] < 20) {
            hF = 1.0f;
            sF = 0.5f;
            vF = 1.0f;
        } else if (avgClrs[0] < 20 && stdDvsClrs[0] < 0.5 && stdDvsClrs[5] < 20) {
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
        // obj cost = 5*(h-hmax) + (s-smax) + 0.5*(v-vmax)
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
        
        add.add(starterIdx);
        
        float[] strtClrs = hsvlch[starterIdx];
        
        for (int i = 0; i < hsvlch.length; ++i) {
            if (i == starterIdx) {continue;}
            
            int[] xyCen = ch.calculateRoundedXYCentroids(labeledSets.get(i));
            
            float[] clrs = hsvlch[i];
            
            float diffH = Math.abs(clrs[0] - strtClrs[0]);
            if (diffH > 0.1) {
                System.out.println("removed by h: " + Arrays.toString(xyCen)
                    + " " + Arrays.toString(clrs));
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
                System.out.println("removed by h, lch: diff=" + diff 
                    + Arrays.toString(xyCen)
                    + " " + " " 
                    + Arrays.toString(clrs));
                continue;
            }
            
            float diffV = Math.abs(clrs[2] - strtClrs[2]);
            if (diffV > 0.05) {
                System.out.println("removed by v: " + Arrays.toString(xyCen)
                    + " " + Arrays.toString(clrs));
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
