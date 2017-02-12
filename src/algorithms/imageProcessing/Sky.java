package algorithms.imageProcessing;

import algorithms.imageProcessing.features.mser.MSEREdges;
import algorithms.imageProcessing.features.mser.Region;
import algorithms.util.PairInt;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.map.TIntObjectMap;
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
    
    public GreyscaleImage extractSkyMask() {
        
        /*
        -- sky mask:
       -- need to find starting region then all regions of sky if possible.
          could possibly start w/ a few things.
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
          (5) could make an assumption about the orientation of the camera place, that is, decreasing y
              pixel coord is direction "up" where sky is found.
        */
        
        // looking at labeled regions that
        //  have constancy in color but may have a gradient in illumination.
        //  blue - can use the polar theta images or mser edges (ptRegions[1]) 
        //  then the lch, c image, which is the magnitude of the
        //   LUV u and v, can find the gradual change if any.
        //  dark red - same except ptRegions[0]
        mserEdges._debugOrigRegions(3, "_pt_1_acc_");
        
        List<Set<PairInt>> labeledSets = mserEdges.getLabeledSets();
        
        int[] pointLabelMap = createPixLabelMap(labeledSets);
        
        ImageProcessor imageProcessor = new ImageProcessor();
        
        if (lma == null) {
            lma = imageProcessor.createLCHForLUV(img);
        }
        
        List<Region> regionPt1 = mserEdges.getOrigGsPtRegions().get(3);
        
        // idxs not same as other lists
        List<Set<PairInt>> labeledInRegion1 = createRegionLabeledSubSets(
            regionPt1, pointLabelMap);
        
        // gather h of hsv, c and h of lch of each labeled region
        // 6 X n
        float[][] hsvlch1 = calcHSVLCH(labeledInRegion1);
        
        if (debug) {
            MiscellaneousCurveHelper ch = new MiscellaneousCurveHelper();
            for (int i = 0; i < labeledInRegion1.size(); ++i) {
                Set<PairInt> set = labeledInRegion1.get(i);
                int[] xyCen = ch.calculateRoundedXYCentroids(set);
                float[] clrs = hsvlch1[i];
                System.out.println("xy=" + Arrays.toString(xyCen) + 
                    " hsvlch=" + Arrays.toString(clrs) + " n=" + set.size());
            }
        }
        
        //List<Region> regionPt0 = mserEdges.getOrigGsPtRegions().get(2);
        
        return null;
        //throw new UnsupportedOperationException("not ready for use");
    }
    
    public void setToDebug() {
        debug = true;
    
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
    
    public static class SkyObject {
        Set<PairInt> points;
        int[] xyCenter;
    }
}
