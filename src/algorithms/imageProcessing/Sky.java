package algorithms.imageProcessing;

import algorithms.QuickSort;
import algorithms.imageProcessing.features.mser.MSEREdges;
import algorithms.imageProcessing.features.mser.Region;
import algorithms.misc.MiscDebug;
import algorithms.util.OneDIntArray;
import algorithms.util.PairInt;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntObjectHashMap;
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
        -- sky color hue range
        -- cells on image border w/ color in hue range.
        -- average and std. dev. of those cells.
        -- 
        when find the decider, create method with start region of sky 
        defined.  this is to allow many other means to determine a valid
        starting sky cell.
        
        hsvlch=[0.5953295,  0.39462394, 0.851081, 173.90898, 50.907433, 175.28432] n=7119
        hsvlch=[0.6043746,  0.31119394, 0.7372585, 159.8233, 35.561916, 176.82475] n=19332
        hsvlch=[0.5803546,  0.75092393, 0.9483775, 159.40683, 92.31827, 176.29259] n=6504
        hsvlch=[0.58049035, 0.74868315, 0.9470542, 159.31277, 92.0528, 176.29291] n=6647
        hsvlch=[0.5811613,  0.7009847,  0.5313757, 95.19146, 47.15611, 174.95387] n=7610
        hsvlch=[0.6201314,  0.6317332,  0.68089944, 107.000824, 72.094025, 183.30006] n=6062
        hsvlch=[0.5826857,  0.44832748, 0.9766973, 195.71774, 61.26734, 172.22879] n=7107
        
        contrail:
        hsvlch=[0.56487614, 0.68869406, 0.2976, 56.69307, 21.15481, 169.08089] n=40146
        
        hsvlch=[0.5829515,  0.5649942,  0.47821, 93.125145, 36.189724, 173.60103] n=3484
        
        white clouds:
        hsvlch=[0.57725304, 0.2487225,  0.87804, 201.12279, 30.319664, 167.74161] n=8437
        
        hsvlch=[0.54253477, 0.2853499,  0.77005726, 183.51671, 27.280539, 154.57755] n=10323
        hsvlch=[0.5433781,  0.33766806, 0.7310331, 171.6103, 30.494549, 155.42032] n=8622
        
        hsvlch=[0.6017986,  0.39862663, 0.590903, 122.25383, 36.566723, 176.97502] n=1761
        hsvlch=[0.92348933, 0.34691936, 0.578129, 120.40768, 33.578182, 236.87389] n=7636
        hsvlch=[0.66858685, 0.25035322, 0.42694178, 95.210014, 19.515263, 189.7619] n=819
        
        hsvlch=[0.0940289,  0.5701533,  0.9364302, 197.72493, 73.67476, 35.64526] n=9187
        hsvlch=[0.85543936, 0.42690355, 0.74646324, 146.43176, 58.856434, 224.26616] n=14112

        hsvlch=[0.9579757,  0.14876272, 0.72665846, 173.23404, 16.787233, 224.12766] n=47
        
        hsvlch=[0.85543936, 0.42690355, 0.74646324, 146.43176, 58.856434, 224.26616] n=14112
        
        hsvlch=[0.0940289, 0.5701533, 0.9364302, 197.72493, 73.67476, 35.64526] n=9187
        
        stonehenge:
        hsvlch=[0.96825397, 0.09251101, 0.8901961, 216.0, 12.0, 254.0] n=1
        
        hsvlch=[0.5853093, 0.39359564, 0.7887473, 165.84178, 44.59479, 172.24834] n=5992
        hsvlch=[0.5899351, 0.10313386, 0.249918, 64.911804, 3.5182912, 171.74931] n=7982
        */
        
        // looking at labeled regions that
        //  have constancy in color but may have a gradient in illumination.
        //  blue - can use the polar theta images or mser edges (ptRegions[1]) 
        //  then the lch, c image, which is the magnitude of the
        //   LUV u and v, can find the gradual change if any.
        //  dark red - same except ptRegions[0]
        //mserEdges._debugOrigRegions(3, "_pt_1_acc_");
        
        List<Set<PairInt>> labeledSets = mserEdges.getLabeledSets();
        
        int[] pointLabelMap = createPixLabelMap(labeledSets);
        
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
            labeledInRegion = createRegionLabeledSubSets(regionPt, pointLabelMap);

            // gather hsv, and lch of each labeled region
            // 6 X n
            hsvlch = calcHSVLCH(labeledInRegion);

            // filter labeledInRegion1 for sky colors over a generous range in hue
            // key=index in labeledInRegion, value=set of pairint of the region
            if (debug) {
                MiscellaneousCurveHelper ch = new MiscellaneousCurveHelper();
                for (int i = 0; i < labeledInRegion.size(); ++i) {
                    Set<PairInt> set = labeledInRegion.get(i);
                    int[] xyCen = ch.calculateRoundedXYCentroids(set);
                    float[] clrs = hsvlch[i];
                    if (nIter == 1) {
                        if ((clrs[2] > 0.35) && (clrs[0] < 0.18 || clrs[0] > 0.49)) {
                            if (setDoesBorderImage(set, img)) {
                                System.out.println(debugLabel + 
                                    " xy=" + Arrays.toString(xyCen) + 
                                    " hsvlch=" + Arrays.toString(clrs) + " n=" + set.size());
                                filteredLabeledRegion.add(set);
                                filteredHSVLCH.add(new OneDFloatArray(clrs));
                            }
                        }
                    } else {
                        if ((clrs[2] > 0.20) && (clrs[0] < 0.18 || clrs[0] > 0.49)) {
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
        }
        
        // creating an objective function to sert the segmented cells to find a
        // starting cell which is most likely sky.
        // at this point, we have regions filtered to contain only those on the
        // border of the image having a color within a range of sky colors.
        
        // if nIter==1, the negative polar theta image was used.  blue is in the
        // middle of the theta range, so these pick up the sky regions if they're
        // blue. (NOTE, this logic may need to be edited with more testing).
        if (nIter == 1) {
            sortForBlue(filteredLabeledRegion, filteredHSVLCH);
        }
        
        // plot the first in filtered lists.  it is on the image
        // border, is within sky color range, and is the largest
        // of the filtered segmented cells
        if (debug) {
            Image imgCp = img.copyImage();
            Set<PairInt> s0 = filteredLabeledRegion.get(0);
            ImageIOHelper.addCurveToImage(s0, imgCp, 0, 0, 255, 0);
            for (int i = 0; i < mserEdges.getEdges().size(); ++i) {
                ImageIOHelper.addCurveToImage(mserEdges.getEdges().get(i), 
                    imgCp, 0, 255, 0, 0);
            }
            MiscDebug.writeImage(imgCp, "_" + debugLabel + "_first_");
        }
        
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
        
        float maxS = Float.MIN_VALUE;
        float maxH = Float.MIN_VALUE;
        float maxV = Float.MIN_VALUE;
        
        for (OneDFloatArray a : filteredHSVLCH) {
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
            float cost = 5.f*(Math.abs(a[0] - maxH)) +
                (Math.abs(a[1] - maxS)) + 
                0.5f * (Math.abs(a[2] - maxV));
            costs[i] = cost;
        }
        
        // sort filtered lists by costs
        QuickSort.sortBy1stArg(costs, filteredHSVLCH, filteredLabeledRegion);
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
