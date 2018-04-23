package algorithms.imageProcessing.features;

import algorithms.QuickSort;
import algorithms.compGeometry.PerimeterFinder2;
import algorithms.imageProcessing.ColorHistogram;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.SummedAreaTable;
import algorithms.imageProcessing.TrimmedImage;
import algorithms.imageProcessing.features.mser.Canonicalizer;
import algorithms.imageProcessing.features.mser.Canonicalizer.RegionPoints;
import algorithms.imageProcessing.features.mser.MSER;
import algorithms.imageProcessing.features.mser.MSER.Threshold;
import algorithms.imageProcessing.features.mser.Region;
import algorithms.imageProcessing.matching.CMODE;
import algorithms.imageProcessing.matching.MSERMatcher;
import algorithms.misc.MiscDebug;
import algorithms.misc.MiscMath;
import algorithms.util.PairInt;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * a class that finds a template object in another image where the
 * object may have changed poses, may have different lighting, or
 * may have different foreground or background.
 *
 * NOTE that future options might include ability to choose between
 * 3 or more different light sources and the ability to use an
 * articulated search for small object matches.
 *
 * @author nichole
 */
public class ObjectMatcher {

    private boolean debug = false;

    public void setToDebug() {
        debug = true;
    }

    private void filterCloseToBounds2(List<Region> regions,
        int width, int height, int border) {

        int[] xy = new int[2];
        for (int i = (regions.size() - 1); i > -1; --i) {
            Region r = regions.get(i);
            r.calculateXYCentroid(xy, width, height);
            if (xy[0] < border || xy[1] < border ||
                (xy[0] >= (width - border)) ||
                (xy[1] >= (height - border))) {
                regions.remove(i);
            }
        }
    }

    private TrimmedImage trim(ImageExt img, Set<PairInt> shape,
        int buffer) {

        int[] minMaxXY = MiscMath.findMinMaxXY(shape);

        int x0 = minMaxXY[0] - buffer;
        if (x0 < 0) {
            x0 = 0;
        }
        int x1 = minMaxXY[1] + buffer;
        if (x1 >= img.getWidth()) {
            x1 = img.getWidth() - 1;
        }
        int y0 = minMaxXY[2] - buffer;
        if (y0 < 0) {
            y0 = 0;
        }
        int y1 = minMaxXY[3] + buffer;
        if (y1 >= img.getHeight()) {
            y1 = img.getHeight() - 1;
        }

        TrimmedImage trImg = new TrimmedImage(img, x0, x1, y0, y1);

        return trImg;
    }

    private void mask(Image img, Set<PairInt> shape0) {

        for (int i = 0; i < img.getWidth(); ++i) {
            for (int j = 0; j < img.getHeight(); ++j) {
                PairInt p = new PairInt(i, j);
                if (!shape0.contains(p)) {
                    img.setRGB(i, j, 0, 0, 0);
                }
            }
        }
    }

    private void mask(GreyscaleImage img, Set<PairInt> shape0) {

        for (int i = 0; i < img.getWidth(); ++i) {
            for (int j = 0; j < img.getHeight(); ++j) {
                PairInt p = new PairInt(i, j);
                if (!shape0.contains(p)) {
                    img.setValue(i, j, 0);
                }
            }
        }
    }

    private List<Region> createCombinedMSERRegions(GreyscaleImage gsImg,
        GreyscaleImage luvTheta, CMODE clrMode, CMODE ptMode,
        boolean fewerMSER, String debugLabel) {

        MSER mser = new MSER();

        List<List<Region>> regionsT = new ArrayList<List<Region>>();
        List<List<Region>> regions = new ArrayList<List<Region>>();
        Threshold thrGs;
        Threshold thrPt;
        if (fewerMSER) {
            thrGs = Threshold.LEAST_SENSITIVE;
            thrPt = Threshold.LESS_SENSITIVE;
        } else {
            thrGs = Threshold.SLIGHTLY_LESS_SENSITIVE;
            thrPt = Threshold.DEFAULT;
        }
        
        if (debug) {
            System.out.println(debugLabel + "  clrMode=" + clrMode.name() 
                + " ptMode=" + ptMode.name());
        }
        
        if (clrMode.equals(CMODE.WHITE)) {
            int[] gsA = MSER.readIntoArray(gsImg);
            List<Region> list = mser.findRegionsNeg(gsA,
                gsImg.getWidth(), gsImg.getHeight(), thrGs);
            regions.add(new ArrayList<Region>());
            regions.add(list);
        } else if (clrMode.equals(CMODE.BLACK)) {
            int[] gsA = MSER.readIntoArray(gsImg);
            List<Region> list = mser.findRegionsPos(gsA,
                gsImg.getWidth(), gsImg.getHeight(), thrGs);
            regions.add(list);
            regions.add(new ArrayList<Region>());
        } else {
            regions = mser.findRegions(gsImg, thrGs);
        }
        
        if (ptMode.equals(CMODE.WHITE)) {
            int[] ptA = MSER.readIntoArray(luvTheta);
            List<Region> list = mser.findRegionsNeg(ptA,
                luvTheta.getWidth(), luvTheta.getHeight(), thrPt);
            regionsT.add(new ArrayList<Region>());
            regionsT.add(list);
        } else if (ptMode.equals(CMODE.BLACK)) {
            int[] ptA = MSER.readIntoArray(luvTheta);
            List<Region> list = mser.findRegionsPos(ptA,
                luvTheta.getWidth(), luvTheta.getHeight(), thrPt);
            regionsT.add(list);
            regionsT.add(new ArrayList<Region>());
        } else {
            regionsT = mser.findRegions(luvTheta, thrPt);
        
            int[] xyCen = new int[2];
            
            // filter to remove all w/ variation > 0
            for (int type = 0; type < 2; ++type) {
                List<Region> list = regionsT.get(type);
                for (int i = (list.size() - 1); i > -1; --i) {
                    Region r = list.get(i);
                    if ((type == 1) && r.getVariation() > 0.001) {
                        list.remove(i);
                    } else if ((type == 0) && r.getVariation() == 0.0) {
                        //r.calculateXYCentroid(xyCen, gsImg.getWidth(), 
                        //    gsImg.getHeight());
                       // list.remove(i);
                        
//variation;  First and second moments of the region 
//            (x, y, x^2, xy, y^2)
//0;        1255, 941, 105035, 78724, 59043
                    }
                }
            }
        }
        
        /*
        if (debug){
            long ts = MiscDebug.getCurrentTimeFormatted();
            int[] xyCen = new int[2];
            Image imCp;
            for (int type = 0; type < 2; ++type) {
                imCp = gsImg.copyToColorGreyscale();
                int n9 = regions.get(type).size();
                for (int i = 0; i < n9; ++i) {
                    Region r = regions.get(type).get(i);
                    int[] clr = ImageIOHelper.getNextRGB(i);
                    r.drawEllipse(imCp, 0, clr[0], clr[1], clr[2]);
                    r.calculateXYCentroid(xyCen, imCp.getWidth(), imCp.getHeight());
                    ImageIOHelper.addPointToImage(xyCen[0], xyCen[1], imCp,
                        1, 255, 0, 0);
                    
                    //r.calculateXYCentroid(xyCen, imCp.getWidth(), imCp.getHeight());
                    //System.out.println("r=" + r.toString() + " " + Arrays.toString(xyCen));
                }
                MiscDebug.writeImage(imCp, debugLabel + "__regions_gs_"+ type + "_" + ts);
            }
            
            for (int type = 0; type < 2; ++type) {
                imCp = luvTheta.copyToColorGreyscale();
                int n9 = regionsT.get(type).size();
                for (int i = 0; i < n9; ++i) {
                    Region r = regionsT.get(type).get(i);
                    int[] clr = ImageIOHelper.getNextRGB(i);
                    r.drawEllipse(imCp, 0, clr[0], clr[1], clr[2]);
                    r.calculateXYCentroid(xyCen, imCp.getWidth(), imCp.getHeight());
                    ImageIOHelper.addPointToImage(xyCen[0], xyCen[1], imCp,
                        1, 255, 0, 0);
                    //System.out.println(type + " xy=" + xyCen[0] + "," + xyCen[1] 
                    //    + " variation=" + r.getVariation());
                }
                MiscDebug.writeImage(imCp, debugLabel + "__regions_pt_"+ type + "_" + ts);
            }
        }
        */
 
        List<Region> combined = new ArrayList<Region>();
        
        for (int i = 0; i < 2; ++i) {
            for (Region r : regions.get(i)) {
                combined.add(r);
            }
        }
        
        for (int i = 0; i < 2; ++i) {
            for (Region r : regionsT.get(i)) {
                combined.add(r);
            }
        }

        filterCloseToBounds2(combined, gsImg.getWidth(), gsImg.getHeight(), 10);

        return combined;
    }

    private void applyWindowedMean(List<List<GreyscaleImage>> pyr, int halfDimension) {

        SummedAreaTable sumTable = new SummedAreaTable();
        
        for (int i = 0; i < pyr.size(); ++i) {
            List<GreyscaleImage> imgMs = pyr.get(i);
            for (int j = 0; j < imgMs.size(); ++j) {
                GreyscaleImage imgM = imgMs.get(j);
                imgM = sumTable.createAbsoluteSummedAreaTable(imgM);
                imgM = sumTable.applyMeanOfWindowFromSummedAreaTable(imgM,
                    2 * halfDimension + 1);
                imgMs.set(j, imgM);
            }
        }
    }
    
    private void applyWindowedMean2(List<GreyscaleImage> pyr, int halfDimension) {

        SummedAreaTable sumTable = new SummedAreaTable();
        
        for (int i = 0; i < pyr.size(); ++i) {
            GreyscaleImage imgM = pyr.get(i);
            imgM = sumTable.createAbsoluteSummedAreaTable(imgM);
            imgM = sumTable.applyMeanOfWindowFromSummedAreaTable(imgM,
                2 * halfDimension + 1);
            pyr.set(i, imgM);
        }
    }

    private RegionPoints createARegion(Set<PairInt> points, int w, int h) {
        
        Region r = new Region();
        for (PairInt pl : points) {
            r.accumulate(pl.getX(), pl.getY());
        }

        int[] xyCen = new int[2];
        r.calculateXYCentroid(xyCen, w, h);
        int x = xyCen[0];
        int y = xyCen[1];
        assert (x >= 0 && x < w);
        assert (y >= 0 && y < h);
        double[] m = r.calcParamTransCoeff();

        double angle = Math.atan(m[0] / m[2]);
        if (angle < 0) {
            angle += Math.PI;
        }

        double major = 2. * m[4];
        double minor = 2. * m[5];

        double ecc = Math.sqrt(major * major - minor * minor) / major;
        assert (!Double.isNaN(ecc));

        Canonicalizer.RegionGeometry rg = new Canonicalizer.RegionGeometry();
        rg.eccentricity = ecc;
        rg.major = major;
        rg.minor = minor;
        rg.orientation = angle;
        rg.xC = x;
        rg.yC = y;

        RegionPoints rp = new RegionPoints();
        rp.ellipseParams = rg;
        rp.points = new HashSet<PairInt>(points);
        
        return rp;
    }
    
    private void filterByColorHistograms(ImageExt img0, Set<PairInt> shape0, 
        ImageExt img1, TIntObjectMap<RegionPoints> regions1) {
        
        //filter by color hist of hsv, cielab and by CIECH

        ColorHistogram clrHist = new ColorHistogram();

        float upperLimit = 0.27f;// 0.33f;
        
        // make the template histograms from the first scale only
        int[][] template_ch_HSV = clrHist.histogramHSV(img0, shape0);
        int[][] template_ch_LAB = clrHist.histogramCIELAB(img0, shape0);
        int[] tHist = clrHist.histogramCIECH64(img0, shape0);
        
        TIntObjectIterator<RegionPoints> iter = regions1.iterator();
        
        TIntSet rmSet = new TIntHashSet();
        
        for (int i = 0; i < regions1.size(); ++i) {
            iter.advance();
            
            int rIdx = iter.key();
            RegionPoints r = iter.value();
            
            int[][] ch = clrHist.histogramHSV(img1, r.points);
            float intersection = clrHist.intersection(template_ch_HSV, ch);
            if (intersection < upperLimit) {
                rmSet.add(rIdx);
            } else {
                ch = clrHist.histogramCIELAB(img1, r.points);
                intersection = clrHist.intersection(template_ch_LAB, ch);
                if (intersection < upperLimit) {
                    rmSet.add(rIdx);
                } else {
                    int[] tHist1 = clrHist.histogramCIECH64(img1, r.points);
                    intersection = clrHist.intersection(tHist, tHist1);
                    if (intersection < upperLimit) {
                        rmSet.add(rIdx);
                    }
                }
            }
        }
        
        TIntIterator iter2 = rmSet.iterator();
        while (iter2.hasNext()) {
            int rmIdx = iter2.next();
            regions1.remove(rmIdx);
        }
        
        System.out.println("chist filter removed " + rmSet.size());
    }

    private void replaceWithAccumulatedPoints(TIntObjectMap<RegionPoints> regionPoints) {
    
        PerimeterFinder2 finder = new PerimeterFinder2();
        
        TIntObjectIterator<RegionPoints> iter = regionPoints.iterator();
        for (int i = 0; i < regionPoints.size(); ++i) {
            
            iter.advance();
            int rIdx = iter.key();
            
            RegionPoints rp = iter.value();
            
            //NOTE this may need to be revised.
            //  wanting to trim down the points outside of the ellipse if
            //  they are too far away, such as a line of pixels blended into
            //  what is otherwise a more compact object.
            
            Set<PairInt> ellipse = new HashSet<PairInt>(rp.points);
            
            rp.points.clear();
            
            for (int j = 0; j < rp.accX.size(); ++j) {
                int x = rp.accX.get(j);
                int y = rp.accY.get(j);
                PairInt p2 = new PairInt(x, y);
                if (!ellipse.contains(p2)) {
                    // if far away from center, past major axis, do not add
                    int diffX = x - rp.ellipseParams.xC;
                    int diffY = y - rp.ellipseParams.yC;
                    double d = Math.sqrt(diffX * diffX + diffY * diffY);
                    if (d > 1.2 * rp.ellipseParams.major) {
                        continue;
                    }
                }
                rp.points.add(p2);
            }
            
            // fill in embedded spaces.
            // NOTE: this may need to be reconsidered in special cases.
            Set<PairInt> embedded = finder.findEmbeddedGaps(rp.points);
            if (embedded != null) {
                rp.points.addAll(embedded);
            }
            
        }
    }

    private void filterToLargestPartitions(TIntObjectMap<RegionPoints> 
        regionPoints, ImageExt img, Set<PairInt> shape) {
        
        if (regionPoints.size() < 4) {
            return;
        }
        
        // keep the largest region, which is usually the entire shape
        // and keep the largest 2 regions which sum to a union equal to the whole
        
        // some original regions are present in multiplocity due to having 
        // diferent orientation angles, so need to extract only the unique
        // centers
        Map<PairInt, TIntList> centerRIdxMap = new HashMap<PairInt, TIntList>();
        
        int maxV = Integer.MIN_VALUE;
        PairInt maxVXY = null;
        TIntObjectIterator<RegionPoints> iter = regionPoints.iterator();
        for (int i = 0; i < regionPoints.size(); ++i) {
            iter.advance();
            RegionPoints r = iter.value();
            int rIdx = iter.key();
            PairInt xy = new PairInt(r.ellipseParams.xC, r.ellipseParams.yC);
            TIntList rList = centerRIdxMap.get(xy);
            if (rList == null) {
                rList = new TIntArrayList();
                centerRIdxMap.put(xy, rList);
            }
            rList.add(rIdx);
            int n = r.points.size();
            if (n > maxV) {
                maxV = n;
                maxVXY = xy;
            }
        }
        
        TIntList keep = new TIntArrayList();
        keep.add(centerRIdxMap.get(maxVXY).iterator().next());
       
        // skipping the max value, single partition
        int[] ns = new int[centerRIdxMap.size() - 1];
        int[] rIdxs = new int[ns.length];
        int count = 0;
        
        Set<PairInt> added = new HashSet<PairInt>();
        
        iter = regionPoints.iterator();
        for (int i = 0; i < regionPoints.size(); ++i) {
            iter.advance();
            RegionPoints r = iter.value();
            int rIdx = iter.key();
            PairInt xy = new PairInt(r.ellipseParams.xC, r.ellipseParams.yC);
            if (xy.equals(maxVXY) || added.contains(xy)) {
                continue;
            }
            added.add(xy);
            
            ns[count] = r.points.size();
            rIdxs[count] = rIdx;
            count++;
        }
        assert(ns.length == count);
        
        QuickSort.sortBy1stArg(ns, rIdxs);
         
        int sum = 0;
        for (int i = 0; i < ns.length; ++i) {
            sum += ns[i];
        }
        int half = sum/2;
        
        int mIdx = Arrays.binarySearch(ns, half);
        if (mIdx < 0) {
            //(-(insertion point) - 1)
            mIdx = -mIdx - 2;
        }
        if (mIdx == (ns.length - 1)) {
            mIdx--;
        }
        
        // the sum is actually the sum of non-intersecting points to the
        //   definition of middle of sum at mIdx is an approximate
        //   place to make a partition to choose one from either side
        //   which total to max sum.
        
        int maxSum = Integer.MIN_VALUE;
        int maxIdxA = -1;
        int maxIdxB = -1;
        for (int i = 0; i <= mIdx; ++i) {
            int rIdxA = rIdxs[i];
            Set<PairInt> setA = regionPoints.get(rIdxA).points;
            for (int j = (mIdx + 1); j < ns.length; ++j) {
                int rIdxB = rIdxs[j];
                Set<PairInt> setB = regionPoints.get(rIdxB).points;
                Set<PairInt> aMinusB = new HashSet<PairInt>(setA);
                aMinusB.removeAll(setB);
                
                Set<PairInt> bMinusA = new HashSet<PairInt>(setB);
                bMinusA.removeAll(setA);
                
                int sum2 = aMinusB.size() + bMinusA.size();
                
                if (sum2 > maxSum) {
                    maxSum = sum2;
                    maxIdxA = rIdxA;
                    maxIdxB = rIdxB;
                }
            }
        }
        assert(maxIdxA > -1);
        assert(maxIdxB > -1);
        keep.add(maxIdxA);
        keep.add(maxIdxB);
        
        TIntObjectMap<RegionPoints> regionPoints2 
            = new TIntObjectHashMap<RegionPoints>();
        for (int i = 0; i < keep.size(); ++i) {
            int idx = keep.get(i);
            RegionPoints r = regionPoints.get(idx);
            PairInt xy = new PairInt(r.ellipseParams.xC, r.ellipseParams.yC);
            TIntList idxs = centerRIdxMap.get(xy);
            for (int j = 0; j < idxs.size(); ++j) {
                int idx2 = idxs.get(j);
                regionPoints2.put(idx2, regionPoints.get(idx2));
            }
        }
        regionPoints.clear();
        regionPoints.putAll(regionPoints2);
    }

    public static class Settings {
        private boolean useLargerPyramid0 = false;
        private boolean useLargerPyramid1 = false;

        //TODO: refactor to use an enum to avoid inconsistent state
        private boolean useSmallObjectMethod = false;
        private boolean useShapeFinder = false;

        private boolean findVanishingPoints = false;
        
        private boolean useColorFilter = true;
 
        private String lbl = "";
        
        /**
         * @return the useLargerPyramid0
         */
        public boolean isUseLargerPyramid0() {
            return useLargerPyramid0;
        }
        
        public void setDebugLabel(String dbgLabel) {
            this.lbl = dbgLabel;
        }
        
        public String getDebugLabel() {
            return lbl;
        }

        /**
         if this is set, the default number of pyramid
         image 0 images separated by a factor of 2 in scale is increased
         in number to include finer scales, allowing better descriptor
         * matches at the cost of increased runtime.
         */
        public void setToUseLargerPyramid0() {
            this.useLargerPyramid0 = true;
        }

        /**
         * @return the useLargerPyramid1
         */
        public boolean isUseLargerPyramid1() {
            return useLargerPyramid1;
        }

        public void setToFindVnishingPoints() {
            this.findVanishingPoints = true;
        }

        /**
         if this is set, the default number of pyramid
         image 1 images separated by a factor of 2 in scale is increased
         in number to include finer scales, allowing better descriptor
         * matches at the cost of increased runtime.
         */
        public void setToUseLargerPyramid1() {
            this.useLargerPyramid1 = true;
        }

        /**
         * @return the useSmallObjectMethod
         */
        public boolean isUseSmallObjectMethod() {
            return useSmallObjectMethod;
        }

        public boolean isUseShapeFinder() {
            return useShapeFinder;
        }

        public void setToUseShapeFinderMethod() {
            useShapeFinder = true;
        }

        /**
         * changes the method used to one without descriptors and instead
         * uses a shape matcher and color histograms.
         * This should only be used on images where the object to find
         * is thought to be small (16 or so pixels for example in width
         * and height)...the shape matching needs the segmentation to
         * contain the object shape in one labeled region rather than in
         * more than one over segmented regions.
         */
        public void setToUseSmallObjectMethod() {
            this.useSmallObjectMethod = true;
        }

        private boolean isFindVanishingPoints() {
            return findVanishingPoints;
        }

        public void setToExcludeColorFilter() {
            useColorFilter = false;
        }
        public boolean useColorFilter() {
            return useColorFilter;
        }
    }
    
    /**
     * NOT READY FOR USE
     * 
     * 
     * @param img0
     * @param shape0
     * @param img1
     * @param settings for the method
     * @return
     */
    public List<CorrespondenceList> findObject12(ImageExt img0, 
        Set<PairInt> shape0, ImageExt img1, Settings settings) {

        long ts = 0;
        if (debug) {
            ts = MiscDebug.getCurrentTimeFormatted();
        }

        TrimmedImage img0Trim = trim(img0, shape0, 20);

        Set<PairInt> shape0Trimmed = new HashSet<PairInt>();
        for (PairInt p : shape0) {
            PairInt p2 = new PairInt(p.getX() - img0Trim.getXOffset(),
                p.getY() - img0Trim.getYOffset());
            shape0Trimmed.add(p2);
        }

        ImageExt img0Trimmed = (ImageExt)img0Trim.getTrimmed();

        ImageProcessor imageProcessor = new ImageProcessor();

        mask(img0Trimmed, shape0Trimmed);
        CMODE clrMode = CMODE.determineColorMode(img0Trimmed, shape0Trimmed);
        
        /*
        convert the image to cie luv and then calculate polar angle of u and v
        around 0 in degrees (a.k.a. the "H" of LCH color space, but with
        the 1976 CIE LAB which is LUV).
        If maxV of 360, returns full value image,
        */
        //System.out.println("template clrMode=" + clrMode.name());
        GreyscaleImage luvTheta0;
        GreyscaleImage luvTheta1;
        if (clrMode.equals(CMODE.WHITE)) {
            luvTheta0 = imageProcessor.createCIELUVTheta_WideRangeLightness(img0Trimmed, 255);
            luvTheta1 = imageProcessor.createCIELUVTheta_WideRangeLightness(img1, 255);
        } else {
            luvTheta0 = imageProcessor.createCIELUVTheta(img0Trimmed, 255);
            luvTheta1 = imageProcessor.createCIELUVTheta(img1, 255);
        }
        
        imageProcessor.singlePixelFilter(luvTheta0);
        imageProcessor.singlePixelFilter(luvTheta1);

        mask(luvTheta0, shape0Trimmed);
        CMODE ptMode = CMODE.determinePolarThetaMode(luvTheta0, shape0Trimmed);
        
        // ----- create the cRegions for a masked image pyramid of img 0 ====
        
        GreyscaleImage gsImg0 = img0Trimmed.copyToGreyscale2();
        GreyscaleImage gsImg1 = img1.copyToGreyscale2();
        
        boolean fewerMSER = true;
       
        GreyscaleImage tmp00 = gsImg0.copyImage();
        imageProcessor.enhanceContrast(tmp00, 4);
        GreyscaleImage tmp01 = luvTheta0.copyImage();
        imageProcessor.enhanceContrast(tmp01, 4);
        
        // build combined list of regions
        List<Region> regionsComb0 = createCombinedMSERRegions(
            tmp00, tmp01,
            //gsImg0, luvTheta0, 
            clrMode, ptMode, fewerMSER, settings.getDebugLabel() + "_0_");

        int[] xy = new int[2];
        //remove all regions with centers outside of shape0 points
        for (int i = (regionsComb0.size() - 1); i > -1; --i) {
            Region r = regionsComb0.get(i);
            r.calculateXYCentroid(xy, img0Trimmed.getWidth(), img0Trimmed.getHeight());
            PairInt p = new PairInt(xy[0], xy[1]);
            if (!shape0Trimmed.contains(p)) {
                regionsComb0.remove(i);
            }
        }

        fewerMSER = false;
        
        GreyscaleImage tmp10 = gsImg1.copyImage();
        imageProcessor.enhanceContrast(tmp10, 4);
        GreyscaleImage tmp11 = luvTheta1.copyImage();
        imageProcessor.enhanceContrast(tmp11, 4);
        
        if (debug) {            
            //MiscDebug.writeImage(img0Trimmed, "_shape0_mask_");
            MiscDebug.writeImage(luvTheta0, "_luv_mask_" + ts);
            MiscDebug.writeImage(luvTheta1, "_luv_srch_" + ts);
            
            //MiscDebug.writeImage(tmp00, "_gs_enhanced_0_");
            //MiscDebug.writeImage(tmp01, "_luv_enhanced_0_");
            //MiscDebug.writeImage(tmp10, "_gs_enhanced_1_");
            //MiscDebug.writeImage(tmp11, "_luv_enhanced_1_");
        }
        
        List<Region> regionsComb1 = createCombinedMSERRegions(
            tmp10, tmp11,
            //gsImg1, luvTheta1, 
            clrMode, ptMode, fewerMSER, settings.getDebugLabel() + "_1_");
                
        int critSep = 1;//5;
        Canonicalizer.filterBySpatialProximity(critSep, regionsComb0, 
            img0Trimmed.getWidth(), img0Trimmed.getHeight());
        
        Canonicalizer.filterBySpatialProximity(critSep, regionsComb1, 
            img1.getWidth(), img1.getHeight());
        
        List<List<GreyscaleImage>> pyrRGB0 = imageProcessor.buildColorPyramid(
            img0Trimmed, settings.useLargerPyramid0);
        
        List<GreyscaleImage> pyrPT0 = imageProcessor.buildPyramid(
            luvTheta0, settings.useLargerPyramid0);

        List<List<GreyscaleImage>> pyrRGB1 = imageProcessor.buildColorPyramid(
            img1, settings.useLargerPyramid1);

        List<GreyscaleImage> pyrPT1 = imageProcessor.buildPyramid(
            luvTheta1, settings.useLargerPyramid1);
       
        //// applyWindowedMean(pyrRGB0, 1);
        ////applyWindowedMean(pyrRGB1, 1);
        ////applyWindowedMean2(pyrPT0, 1);
        ////applyWindowedMean2(pyrPT1, 1);
        
        Canonicalizer canonicalizer = new Canonicalizer();

        // ----- create the cRegions for a masked image pyramid of img 0 ====

        //TODO: add filter here for patterns in the MSER regions that
        // are strong, and if present in reference frame1, then
        // anything without it mughr be removable.
        // use of this feature should be a Setting option.

        /*
        if (debug) {
            int[] xyCen = new int[2];
            Image im0Cp, im1Cp;
            im0Cp = img0Trimmed.copyImage();
            int n9 = regionsComb0.size();
            for (int i = 0; i < n9; ++i) {
                Region r = regionsComb0.get(i);
                int[] clr = ImageIOHelper.getNextRGB(i);
                r.drawEllipse(im0Cp, 0, clr[0], clr[1], clr[2]);
                r.calculateXYCentroid(xyCen, im0Cp.getWidth(), im0Cp.getHeight());
                ImageIOHelper.addPointToImage(xyCen[0], xyCen[1], im0Cp,
                    1, 255, 0, 0);
            }
            MiscDebug.writeImage(im0Cp, "_" + settings.getDebugLabel() + 
                "___regions_0_");

            im1Cp = img1.copyImage();
            n9 = regionsComb1.size();
            for (int i = 0; i < n9; ++i) {
                Region r = regionsComb1.get(i);
                int[] clr = ImageIOHelper.getNextRGB(i);
                r.drawEllipse(im1Cp, 0, clr[0], clr[1], clr[2]);
                r.calculateXYCentroid(xyCen, im1Cp.getWidth(), im1Cp.getHeight());
                ImageIOHelper.addPointToImage(xyCen[0], xyCen[1], im1Cp,
                    1, 255, 0, 0);
            //    System.out.println("regIdx1=" + i + " x="+xyCen[0] + " y=" + xyCen[1]);
            }
            MiscDebug.writeImage(im1Cp, "_" + settings.getDebugLabel() 
                + "___regions_1_");
        }
        */
        
        TIntObjectMap<RegionPoints> regionPoints0 =
            canonicalizer.canonicalizeRegions2(regionsComb0, pyrRGB0.get(0).get(1));
   
        TIntObjectMap<RegionPoints> regionPoints1 =
            canonicalizer.canonicalizeRegions2(regionsComb1, pyrRGB1.get(0).get(1));
  
        // filter by color hist of hsv, cielab and CIECH
        if (settings.useColorFilter()) {
            filterByColorHistograms(img0Trimmed, shape0Trimmed, img1, 
                regionPoints1);
        }
        
        //NOTE: not sure this is the best approach, but wanting to keep the 
        //   template shapes as 1 full shape and then the 2 largest 
        //   parts of it to allow a finer fragmented search.
        filterToLargestPartitions(regionPoints0, img0Trimmed, shape0Trimmed);
        
        if (debug) {
            int[] xyCen = new int[2];
            Image im0Cp, im1Cp;
            im0Cp = img0Trimmed.copyImage();
            TIntObjectIterator<RegionPoints> iter = regionPoints0.iterator();
            for (int i = 0; i < regionPoints0.size(); ++i) {
                iter.advance();
                int rIdx = iter.key();
                Region r = regionsComb0.get(rIdx);
                int[] clr = ImageIOHelper.getNextRGB(i);
                r.drawEllipse(im0Cp, 0, clr[0], clr[1], clr[2]);
                r.calculateXYCentroid(xyCen, im0Cp.getWidth(), im0Cp.getHeight());
                ImageIOHelper.addPointToImage(xyCen[0], xyCen[1], im0Cp,
                    1, 255, 0, 0);
            }
            MiscDebug.writeImage(im0Cp, "_" + settings.getDebugLabel() + 
                "_regions_0_filterP_");
        }
        
        /*
        NOTE: tried 2 changes in the region points to see if they improved the
        results.
        (1) modified the ellipse boundaries inward to the nearest bounding edges
            of the accumulated points.
        (2) used just the accumulated points instead of the ellipse filled points.
        */
        
        replaceWithAccumulatedPoints(regionPoints1);
        
        if (debug) {
            int[] xyCen = new int[2];
            Image im1Cp = img1.copyImage();
            TIntObjectIterator<RegionPoints> iter = regionPoints1.iterator();
            for (int i = 0; i < regionPoints1.size(); ++i) {
                iter.advance();
                int rIdx = iter.key();
                Region r = regionsComb1.get(rIdx);
                int[] clr = ImageIOHelper.getNextRGB(i);
                r.drawEllipse(im1Cp, 0, clr[0], clr[1], clr[2]);
                r.calculateXYCentroid(xyCen, im1Cp.getWidth(), im1Cp.getHeight());
                ImageIOHelper.addPointToImage(xyCen[0], xyCen[1], im1Cp,
                    1, 255, 0, 0);
            }
            MiscDebug.writeImage(im1Cp, "_" + settings.getDebugLabel() 
                + "_regions_1_filtered_");
            
            /*
            iter = regionPoints1.iterator();
            for (int i = 0; i < regionPoints1.size(); ++i) {
                iter.advance();
                int rIdx = iter.key();
                RegionPoints rp = iter.value();
                Region r = regionsComb1.get(rIdx);
                int[] clr = ImageIOHelper.getNextRGB(i);

                im1Cp = img1.copyImage();
                r.drawEllipse(im1Cp, 0, 255, 0, 0);
                for (PairInt p : rp.points) {
                    ImageIOHelper.addPointToImage(p.getX(), p.getY(), 
                        im1Cp, 0, 10, 255, 10);
                }
                r.calculateXYCentroid(xyCen, im1Cp.getWidth(), im1Cp.getHeight());
                ImageIOHelper.addPointToImage(xyCen[0], xyCen[1], im1Cp,
                    1, 255, 0, 0);
                MiscDebug.writeImage(im1Cp, "_" + settings.getDebugLabel() 
                    + "_regions_1_acc_" + i + "_");
            }
            */
        }
        
        
        MSERMatcher matcher = new MSERMatcher();

        if (debug) {
            matcher.setToDebug();
        }
        
        List<CorrespondenceList> corList = matcher.matchObject0(
            clrMode, pyrRGB0, pyrPT0, regionPoints0,
            pyrRGB1, pyrPT1, regionPoints1, settings);
        
        if (corList == null) {
            return null;
        }
                
        // apply offsets for having trimmed image 0
        for (int i0 = 0; i0 < corList.size(); ++i0) {
            CorrespondenceList topC = corList.get(i0);
            for (int i = 0; i < topC.getPoints1().size(); ++i) {
                PairInt p = topC.getPoints1().get(i);
                int x = p.getX();
                int y = p.getY();
                p.setX(x + img0Trim.getXOffset());
                p.setY(y + img0Trim.getYOffset());
            }
        }

        return corList;
    }
    
}
