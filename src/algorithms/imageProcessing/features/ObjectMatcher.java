package algorithms.imageProcessing.features;

import algorithms.imageProcessing.CIEChromaticity;
import algorithms.imageProcessing.CannyEdgeFilterAdaptive;
import algorithms.imageProcessing.ColorHistogram;
import algorithms.imageProcessing.EdgeFilterProducts;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.GroupPixelCIELUV;
import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.ImageSegmentation;
import algorithms.imageProcessing.matching.ORBMatcher;
import algorithms.imageProcessing.segmentation.LabelToColorHelper;
import algorithms.misc.MiscDebug;
import algorithms.misc.MiscMath;
import algorithms.util.PairInt;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.TFloatList;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
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
    
    public static class Settings {
        private boolean useLargerPyramid0 = false;
        private boolean useLargerPyramid1 = false;
        
        //TODO: refactor to use an enum to avoid inconsistent state
        private boolean useSmallObjectMethod = false;
        private boolean useShapeFinder = false;
        
        /**
         * @return the useLargerPyramid0
         */
        public boolean isUseLargerPyramid0() {
            return useLargerPyramid0;
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
        
    }
        
    /**
     * given an object in image img0 which is defined by shape0, find the same 
     * object in img1.
     * This method was created to handle a range of lighting changes, poses,
     * and different background and foregrounds, that is, cases where the
     * number of true matching keypoints in img1 is small compared to the
     * total number of keypoints in img1, making the standard best matching
     * descriptors and RANSAC to remove outliers infeasible.
     * 
     * This method in detail creates a segmented img1, removes the labeled regions
     * which have small color histogram intersections with the template,
     * and then compares descriptors scale by scale.
     * Pairs of points are created from same labeled regions and each pair
     * that passes a descriptor cost filter is added to a list of
     * pairs of points.  Those pairs of points are used to calculate 
     * transformations between the keypoints.  Transformations close to "1"
     * are further evaluated for descriptor cost, number of matching points,
     * and the distance between transformed points and nearest matching
     * points.  The lowest cost solution is returned as a correspondence list.
     * Note that the smallest number of pyramid images are used by default, 
     * which locates the object, but precise correspondence is better 
     * achieved by setting the flag to use larger number of pyramid images 
     * to true at the expense of runtime or by refining the results
     * with another method (such as an aggregated partial shape matching
     * as in the unfinished ShapeFinder.java).
     * 
     * Note that if the object in img1 is expected to be smaller than 16 or so
     * pixels wide and high, there will not be enough keypoints, so the user
     * should use setting useSmallObjectMethod.
     * 
     * NOTE: in the future, might provide a method that accepts the labeled
     * regions from the user to allow external segmentation to be used.
     * 
     * For best use,images should be binned down to sizes where a dimension
     * is near 256 pixels or less (findObject2 does that automatically).
     * 
     * @param img0
     * @param shape0
     * @param img1
     * @param settings for the method
     * @return 
     */
    public CorrespondenceList findObject(ImageExt img0, Set<PairInt> shape0,
        ImageExt img1, Settings settings) {
    
        long ts = 0;
        if (debug) {
            ts = MiscDebug.getCurrentTimeFormatted();
        }
        
        ImageExt[] imgs0 = maskImage(img0, shape0);
        
        ImageProcessor imageProcessor = new ImageProcessor();
        ImageSegmentation imageSegmentation = new ImageSegmentation();
        
        int nKeypoints = 200;//200;
        ORB orb0 = extractTemplateORBKeypoints(imgs0[0], shape0,
            nKeypoints, !settings.isUseLargerPyramid0(), ts);

        TFloatList sTempList = new TFloatArrayList(
            orb0.getScalesList().size());
        for (int i = 0; i < orb0.getScalesList().size(); ++i) {
            sTempList.add(orb0.getScalesList().get(i).get(0));
        }

        ORB orb1 = new ORB(500);
        orb1.overrideToNotCreateDescriptors();
        if (!settings.isUseLargerPyramid1()) {
            orb1.overrideToUseSmallestPyramid();
        }
        orb1.detectAndExtract(img1);

        TFloatList sList = new TFloatArrayList(orb1.getScalesList().size());
        for (int i = 0; i < orb1.getScalesList().size(); ++i) {
            sList.add(orb1.getScalesList().get(i).get(0));
        }
        
        int ns = sList.size();
        List<TIntList> rmIndexesList = new ArrayList<TIntList>();
        
        ColorHistogram clrHist = new ColorHistogram();

        // make the template histograms from the first scale only
        int[][] template_ch_HSV = null;
        int[][] template_ch_LAB = null;
        {
            List<TIntList> kp0TempList = orb0.getKeyPoint0List();
            List<TIntList> kp1TempList = orb0.getKeyPoint1List();
            Set<PairInt> points0 = new HashSet<PairInt>();
            for (int i = 0; i < kp0TempList.get(0).size(); ++i) {
                int y = kp0TempList.get(0).get(i);
                int x = kp1TempList.get(0).get(i);
                PairInt p = new PairInt(x, y);
                Set<PairInt> points = imageProcessor.getNeighbors(
                    imgs0[0], p);
                points.add(p);
                points0.addAll(points);
            }
            template_ch_HSV = clrHist.histogramHSV(imgs0[1], points0);
            template_ch_LAB = clrHist.histogramCIELAB(imgs0[1], points0);
        }

        // --- filter out key points at each scale, then associated data ----

        for (int i = 0; i < ns; ++i) {
            TIntList kp0 = orb1.getKeyPoint0List().get(i);
            TIntList kp1 = orb1.getKeyPoint1List().get(i);
            TDoubleList or = orb1.getOrientationsList().get(i);
            TFloatList s = orb1.getScalesList().get(i);

            int np = kp0.size();
            TIntList rm = new TIntArrayList();
            for (int j = 0; j < np; ++j) {
                PairInt p = new PairInt(kp1.get(j), kp0.get(j));
                Set<PairInt> points = imageProcessor.getNeighbors(
                    img1, p);
                points.add(p);
                int[][] ch = clrHist.histogramHSV(img1, points);
                float intersection = clrHist.intersection(
                    template_ch_HSV, ch);

                if (intersection < 0.2) {
                    rm.add(j);
                } else {
                    ch = clrHist.histogramCIELAB(img1, points);
                    intersection = clrHist.intersection(
                        template_ch_LAB, ch);
                    if (intersection < 0.2) {
                        rm.add(j);
                    }
                }
            }
            rmIndexesList.add(rm);
        }
        orb1.removeAtIndexes(rmIndexesList);
        
        if (orb1.getKeyPoint0List().isEmpty()) {
            return null;
        }
  
        float luvDeltaELimit = 10;// between 10 and 20, 25
                    
        rmIndexesList = new ArrayList<TIntList>();

        // ---- calculate the template CIE LUV colors ----
        //      -- also calculate separately, the colors of the 
        //         entire shape to compare results when using the
        //         colors as filters.
        Set<PairInt> set0 = new HashSet<PairInt>();
        for (PairInt p : orb0.getKeyPointListColMaj(0)) {
            Set<PairInt> points = imageProcessor.getNeighbors(
                img0, p);
            set0.add(p);
            set0.addAll(points);
        }
        GroupPixelCIELUV meanLUVKeypoints = new GroupPixelCIELUV(set0, img0);
        meanLUVKeypoints.calculateColors(set0, img0, 0, 0);

        CIEChromaticity cieC = new CIEChromaticity();

        for (int i = 0; i < ns; ++i) {
            TIntList kp0 = orb1.getKeyPoint0List().get(i);
            TIntList kp1 = orb1.getKeyPoint1List().get(i);
            TDoubleList or = orb1.getOrientationsList().get(i);
            TFloatList s = orb1.getScalesList().get(i);

            int np = kp0.size();
            TIntList rm = new TIntArrayList();
            for (int j = 0; j < np; ++j) {

                PairInt p = new PairInt(kp1.get(j), kp0.get(j));
                Set<PairInt> points = imageProcessor.getNeighbors(
                    img1, p);
                points.add(p);

                GroupPixelCIELUV luv1 = new GroupPixelCIELUV(points, img1);
                luv1.calculateColors(points, img1, 0, 0);

                // -- for LUV, a deltaE limit of 9?
                double deltaE = cieC.calcDeltaECIE2000(
                    meanLUVKeypoints.getAvgL(),
                    meanLUVKeypoints.getAvgU(),
                    meanLUVKeypoints.getAvgV(),
                    luv1.getAvgL(),
                    luv1.getAvgU(),
                    luv1.getAvgV());

                if (Math.abs(deltaE) > luvDeltaELimit) {
                    rm.add(j);
                }
            }
            rmIndexesList.add(rm);
            
            orb1.removeAtIndexes(rmIndexesList);
                        
            if (orb1.getKeyPoint0List().isEmpty()) {
                return null;
            }
        }
        
        ImageExt img1Cp = img1.copyToImageExt();
        
        //CannyEdgeFilterAdaptive canny
        //    = new CannyEdgeFilterAdaptive();
        //canny.overrideToUseAdaptiveThreshold();
        //canny.applyFilter(img1.copyToGreyscale2());
        //EdgeFilterProducts edgeProduct = canny.getFilterProducts();
        //int[] labels4 = imageSegmentation.objectSegmentation(
        //    img1Cp, edgeProduct);
        
        int[] labels4 = imageSegmentation.objectSegmentation(img1Cp);
        
        if (debug) {
            ImageExt img11 = img1.copyToImageExt();
            ImageIOHelper.addAlternatingColorLabelsToRegion(
                img11, labels4);
            MiscDebug.writeImage(img11, "_segmented_" + ts);
        }

        List<Set<PairInt>> listOfPointSets2
            = LabelToColorHelper.extractContiguousLabelPoints(img1Cp, labels4);

        boolean useCHist = false;

        boolean changed = false;
        
        changed = imageSegmentation.filterByCIECH(imgs0[0],
            shape0, img1Cp, listOfPointSets2, 0.4f);//0.35f
        /*if (useCHist) {
            // ---- filter segmentation by cie theta histograms
            changed = imageSegmentation.filterByCIETheta(imgs0[0],
                shape0, img1Cp, listOfPointSets2);
        } else {
            // needs a very high limit to not remove objects w/ different illum
            // and results in only a few removals.
            luvDeltaELimit = 40;
            changed = imageSegmentation.filterByLUVDeltaE(imgs0[0],
                shape0, img1Cp, listOfPointSets2, luvDeltaELimit);
        }*/
        
        // ---- remove keypoints if not in segmented cells
        if (changed) {
            TObjectIntMap<PairInt> pointIdxMap
                = new TObjectIntHashMap<PairInt>();
            for (int j = 0; j < listOfPointSets2.size(); ++j) {
                Set<PairInt> set = listOfPointSets2.get(j);
                for (PairInt p : set) {
                    pointIdxMap.put(p, j);
                }
            }

            rmIndexesList = new ArrayList<TIntList>();
            for (int i = 0; i < ns; ++i) {
                TIntList kp0 = orb1.getKeyPoint0List().get(i);
                TIntList kp1 = orb1.getKeyPoint1List().get(i);
                TDoubleList or = orb1.getOrientationsList().get(i);
                TFloatList s = orb1.getScalesList().get(i);

                int np = kp0.size();
                TIntList rm = new TIntArrayList();
                for (int j = 0; j < np; ++j) {
                    PairInt p = new PairInt(kp1.get(j), kp0.get(j));
                    if (!pointIdxMap.containsKey(p)) {
                        rm.add(j);
                    }
                }
                rmIndexesList.add(rm);
                System.out.println("rm at scale " + i + " n=" +
                    rm.size());
            }
            orb1.removeAtIndexes(rmIndexesList);
        }
        
        // TODO: if point sets size < 3 or so, remove it
        
        if (orb1.getKeyPoint0List().isEmpty()) {
            return null;
        }

        if (debug) {   // ----------- segmentation -------
            Set<PairInt> kpSet = new HashSet<PairInt>();
            TIntList kp0 = orb1.getKeyPoint0List().get(0);
            TIntList kp1 = orb1.getKeyPoint1List().get(0);
            for (int i = 0; i < kp0.size(); ++i) {
                int x = kp1.get(i);
                int y = kp0.get(i);
                kpSet.add(new PairInt(x, y));
            }

            ImageExt img11 = img1.createWithDimensions();
            ImageIOHelper.addAlternatingColorPointSetsToImage(
                listOfPointSets2, 0, 0, 1, img11);
            ImageIOHelper.addCurveToImage(kpSet, img11,
                1, 255, 0, 0);
            MiscDebug.writeImage(img11,
                "_filtered_segmentation_1_" + ts);

            // plot the segmentation in black and white
            // and then the sets in alternating color
            img11 = img11.copyToGreyscale2().copyToColorGreyscaleExt();
            TObjectIntMap<PairInt> pointLabels2
                = new TObjectIntHashMap<PairInt>();
            for (int i = 0; i < listOfPointSets2.size(); ++i) {
                int clr = ImageIOHelper.getNextColorRGB(i);
                Set<PairInt> set = listOfPointSets2.get(i);
                for (PairInt p : kpSet) {
                    if (set.contains(p)) {
                        ImageIOHelper.addPointToImage(p.getX(), p.getY(),
                            img11, 1, clr);
                    }
                }
            }
            MiscDebug.writeImage(img11,
                "_filtered_segmentation_2_" + ts);
        }

        sList = new TFloatArrayList(orb1.getScalesList().size());
        for (int i = 0; i < orb1.getScalesList().size(); ++i) {
            sList.add(orb1.getScalesList().get(i).get(0));
        }

        List<CorrespondenceList> corList;

        List<Set<PairInt>> tempListOfPointSets
            = new ArrayList<Set<PairInt>>();
        tempListOfPointSets.add(shape0);

        if (debug) {
            ImageExt img11 = img1.copyToImageExt();
            TIntList kp0 = orb1.getKeyPoint0List().get(0);
            TIntList kp1 = orb1.getKeyPoint1List().get(0);
            for (int i = 0; i < kp1.size(); ++i) {
                int x = kp1.get(i);
                int y = kp0.get(i);
                ImageIOHelper.addPointToImage(x, y, img11, 1, 255, 0, 0);
            }
            MiscDebug.writeImage(img11,"_kp_2_" + ts);
        }

        if (debug) {
            GreyscaleImage theta0 = imageProcessor.createCIELABTheta(
                imgs0[0], 255);
            MiscDebug.writeImage(theta0, "_theta_0_" +  ts);
            GreyscaleImage theta1 = imageProcessor.createCIELABTheta(img1, 
                255);
            MiscDebug.writeImage(theta1, "_theta_1_" +  ts);
        }
        
        orb0.createDescriptorsLABTheta(imgs0[0]);
        orb1.createDescriptorsLABTheta(img1);
        if (settings.isUseSmallObjectMethod()) {
            corList = ORBMatcher.matchSmall(orb0, orb1,
                shape0, listOfPointSets2);
        } else if (settings.isUseShapeFinder()) {
            corList = ORBMatcher.matchAggregatedShape(orb0, orb1,
                shape0, listOfPointSets2);
        } else {
            corList = ORBMatcher.match0(orb0, orb1,
                shape0, listOfPointSets2);
        }

        if (corList == null || corList.isEmpty()) {
            return null;
        }

        return corList.get(0);        
    }
    
    // handles the binning to smaller size automatically
    public CorrespondenceList findObject2(ImageExt img0, Set<PairInt> shape0,
        ImageExt img1, Settings settings) {
    
        int maxDimension = 256;
        
        int w0 = img0.getWidth();
        int h0 = img0.getHeight();
        
        int binFactor0 = (int) Math.ceil(Math.max((float) w0 / maxDimension,
            (float) h0 / maxDimension));
             
        int w1 = img1.getWidth();
        int h1 = img1.getHeight();
        
        int binFactor1 = (int) Math.ceil(Math.max((float) w1 / maxDimension,
            (float) h1 / maxDimension));
        
        ImageProcessor imageProcessor = new ImageProcessor();
        Set<PairInt> shape0_2 = imageProcessor.binPoints(shape0, binFactor0);
        ImageExt img0_2 = imageProcessor.binImage(img0, binFactor0);
        
        ImageExt img1_2 = imageProcessor.binImage(img1, binFactor1);
        
        CorrespondenceList cor = findObject(img0_2, shape0_2, img1_2,
            settings);
        
        // --- TODO: transform results back to full scale
        
        throw new UnsupportedOperationException("not yet implemented");
    }

    private ImageExt[] maskImage(ImageExt img, Set<PairInt> shape) {

        ImageExt imgMasked = img.createWithDimensions();
        
        for (PairInt p : shape) {
            imgMasked.setRGB(p.getX(), p.getY(), img.getRGB(p.getX(), p.getY()));
        }
        
        return new ImageExt[]{img, imgMasked};
    }
      
    private ORB extractTemplateORBKeypoints(ImageExt img,
        Set<PairInt> shape0, int nKeypoints, 
        boolean useSmallPyramid, long debugTs) {

        int[] minMaxXY = MiscMath.findMinMaxXY(shape0);

        int w = img.getWidth();
        int h = img.getHeight();
        
        int buffer = 20;

        int xLL = minMaxXY[0] - buffer;
        if (xLL < 0) {
            xLL = 0;
        }
        int yLL = minMaxXY[2] - buffer;
        if (yLL < 0) {
            yLL = 0;
        }
        int xUR = minMaxXY[1] + buffer;
        if (xUR > (w - 1)) {
            xUR = w - 1;
        }
        int yUR = minMaxXY[3] + buffer;
        if (yUR > (h - 1)) {
            yUR = h - 1;
        }
                
        ORB orb = ORBWrapper.extractKeypointsFromSubImage(
            img, xLL, yLL, xUR, yUR, nKeypoints, useSmallPyramid);
                        
        // trim orb data that is outside of shape
        int ns = orb.getKeyPoint0List().size();
        
        for (int i = 0; i < ns; ++i) {
            TIntList kp0 = orb.getKeyPoint0List().get(i);
            TIntList kp1 = orb.getKeyPoint1List().get(i);
            TDoubleList or = orb.getOrientationsList().get(i);
            TFloatList s = orb.getScalesList().get(i);
            
            int n0 = kp0.size();
            
            TIntList rm = new TIntArrayList();
            for (int j = 0; j < n0; ++j) {
                PairInt p = new PairInt(kp1.get(j), kp0.get(j));
                if (!shape0.contains(p)) {
                    rm.add(j);
                }
            }
            if (!rm.isEmpty()) {
                int nb = n0 - rm.size();
                
                TIntSet rmSet = new TIntHashSet(rm);
                for (int j = (rm.size() - 1); j > -1; --j) {
                    int idx = rm.get(j);
                    kp0.removeAt(idx);
                    kp1.removeAt(idx);
                    or.removeAt(idx);
                    s.removeAt(idx);
                }
                int count = 0;
                for (int j = 0; j < n0; ++j) {
                    if (rmSet.contains(j)) {
                        continue;
                    }
                    count++;
                }
                assert(count == nb);
            }
        }
        
        if (debug) {// DEBUG print each pyramid to see if has matchable points
            // might need to change the ORb response filter to scale by scale level
            for (int i0 = 0; i0 < orb.getKeyPoint0List().size(); ++i0) {
                Image img0Cp = img.copyImage();
                float scale = orb.getScalesList().get(i0).get(0);
                for (int i = 0; i < orb.getKeyPoint0List().get(i0).size(); ++i) {
                    int y = orb.getKeyPoint0List().get(i0).get(i);
                    int x = orb.getKeyPoint1List().get(i0).get(i);
                    ImageIOHelper.addPointToImage(x, y, img0Cp, 
                        1, 255, 0, 0);
                }
                String str = Integer.toString(i0);
                if (str.length() < 2) {
                    str = "0" + str;
                }
                MiscDebug.writeImage(img0Cp, "_template_orb" + str + "_" + debugTs);
            }
        }
        
        return orb;
    }

}
