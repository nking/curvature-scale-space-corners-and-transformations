package algorithms.imageProcessing.features;

import algorithms.compGeometry.FurthestPair;
import algorithms.imageProcessing.CIEChromaticity;
import algorithms.imageProcessing.ColorHistogram;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.GroupPixelCIELUV;
import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.ImageSegmentation;
import algorithms.imageProcessing.MiscellaneousCurveHelper;
import algorithms.imageProcessing.SummedAreaTable;
import algorithms.imageProcessing.TrimmedImage;
import algorithms.imageProcessing.VanishingPoints;
import algorithms.imageProcessing.features.mser.Canonicalizer;
import algorithms.imageProcessing.features.mser.Canonicalizer.CRegion;
import algorithms.imageProcessing.features.mser.Canonicalizer.RegionPoints;
import algorithms.imageProcessing.features.mser.MSER;
import algorithms.imageProcessing.features.mser.Region;
import algorithms.imageProcessing.matching.MSERMatcher;
import algorithms.imageProcessing.matching.ORBMatcher;
import algorithms.imageProcessing.segmentation.LabelToColorHelper;
import algorithms.imageProcessing.util.PairIntWithIndex;
import algorithms.misc.MiscDebug;
import algorithms.misc.MiscMath;
import algorithms.util.OneDIntArray;
import algorithms.util.PairInt;
import com.climbwithyourfeet.clustering.DTClusterFinder;
import gnu.trove.iterator.TIntIntIterator;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.TFloatList;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map.Entry;
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

    private void debugPrint(List<TIntObjectMap<CRegion>> cRegionsList,
        List<List<GreyscaleImage>> pyr, String label) {

        for (int j = 0; j < pyr.size(); ++j) {

            Image img1 = pyr.get(j).get(0).copyToColorGreyscale();

            TIntObjectMap<CRegion> crMap = cRegionsList.get(j);
            TIntObjectIterator<CRegion> iter = crMap.iterator();

            int nExtraDot = 0;

            for (int ii = 0; ii < crMap.size(); ++ii) {
                iter.advance();

                int idx = iter.key();

                CRegion cr = iter.value();

                int[] clr = ImageIOHelper.getNextRGB(ii);

                cr.draw(img1, nExtraDot, clr[0], clr[1], clr[2]);
            }

            MiscDebug.writeImage(img1, label + "_" + j + "_crs_");
        }
    }

    private void debugPrint2(TIntObjectMap<CRegion> cRegions,
        List<GreyscaleImage> rgb, String label) {

        Image img1 = rgb.get(1).copyToColorGreyscale();

        TIntObjectIterator<CRegion> iter = cRegions.iterator();

        int nExtraDot = 0;

        for (int ii = 0; ii < cRegions.size(); ++ii) {
            iter.advance();

            int idx = iter.key();

            CRegion cr = iter.value();

            int[] clr = ImageIOHelper.getNextRGB(ii);

            cr.draw(img1, nExtraDot, clr[0], clr[1], clr[2]);
        }

        MiscDebug.writeImage(img1, label + "_" + "_crs_");
    }
    
    private void debugPrint2(TIntObjectMap<CRegion> cRegions,
        GreyscaleImage gs, String label) {

        Image img1 = gs.copyToColorGreyscale();

        TIntObjectIterator<CRegion> iter = cRegions.iterator();

        int nExtraDot = 0;

        for (int ii = 0; ii < cRegions.size(); ++ii) {
            iter.advance();

            int idx = iter.key();

            CRegion cr = iter.value();

            int[] clr = ImageIOHelper.getNextRGB(ii);
            cr.drawEachPixel(img1, nExtraDot, clr[0], clr[1], clr[2]);
            cr.draw(img1, nExtraDot, clr[0], clr[1], clr[2]);
        }

        MiscDebug.writeImage(img1, label + "_" + "_crs_");
    }

    private double findMinAutoCorrel(TIntObjectMap<CRegion> cRegions) {

        TIntObjectIterator<CRegion> iter = cRegions.iterator();

        double min = Double.MAX_VALUE;

        for (int i = 0; i < cRegions.size(); ++i) {
            iter.advance();
            CRegion cr = iter.value();
            if (cr.autocorrel < min) {
                min = cr.autocorrel;
            }
        }

        return min;
    }

    private void removeSmallAutoCorrel(List<TIntObjectMap<CRegion>>
        cRegionsList, double minAutoCor, float f) {

        float limit = f * (float)minAutoCor;

        TIntList rm = new TIntArrayList();

        TIntObjectMap<CRegion> lgMap = cRegionsList.get(0);

        TIntObjectIterator<CRegion> iter = lgMap.iterator();
        for (int i = 0; i < lgMap.size(); ++i) {
            iter.advance();
            int idx = iter.key();
            CRegion cr = iter.value();
            if (cr.autocorrel < limit) {
                rm.add(idx);
            }
        }

        // === remove indexes rm from each octave
        for (int octave = 0; octave < cRegionsList.size(); ++octave) {
            TIntObjectMap<CRegion> map = cRegionsList.get(octave);
            for (int i = 0; i < rm.size(); ++i) {
                int rmIdx = rm.get(i);
                map.remove(rmIdx);
            }
        }
    }

    private void removeSmallAutoCorrel2(TIntObjectMap<CRegion>
        cRegionsMap, double minAutoCor, float f) {

        float limit = f * (float)minAutoCor;

        TIntList rm = new TIntArrayList();

        TIntObjectIterator<CRegion> iter = cRegionsMap.iterator();
        for (int i = 0; i < cRegionsMap.size(); ++i) {
            iter.advance();
            int idx = iter.key();
            CRegion cr = iter.value();
            if (cr.autocorrel < limit) {
                rm.add(idx);
            }
        }

        for (int i = 0; i < rm.size(); ++i) {
            int rmIdx = rm.get(i);
            cRegionsMap.remove(rmIdx);
        }
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

    private void filterCloseToBounds(List<List<Region>> regions,
        int width, int height, int border) {

        int[] xy = new int[2];
        for (int rIdx = 0; rIdx < 2; ++rIdx) {
            List<Region> list = regions.get(rIdx);
            for (int i = (list.size() - 1); i > -1; --i) {
                Region r = list.get(i);
                r.calculateXYCentroid(xy, width, height);
                if (xy[0] < border || xy[1] < border ||
                    (xy[0] >= (width - border)) ||
                    (xy[1] >= (height - border))) {
                    list.remove(i);
                }
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

    private TIntObjectMap<Set<PairInt>> segment(ImageExt img,
        Set<PairInt> shape) {

        ImageSegmentation imageSegmentation = new ImageSegmentation();

        int[] labels = imageSegmentation.objectSegmentation(img);

        TIntObjectMap<Set<PairInt>> labeledSets =
            LabelToColorHelper.extractLabelPoints(img, labels);

        for (int i = 0; i < labels.length; ++i) {
            int x = img.getCol(i);
            int y = img.getRow(i);

            PairInt p = new PairInt(x, y);
            if (!shape.contains(p)) {
                int label = labels[i];
                labeledSets.remove(label);
            }
        }

        return labeledSets;
    }

    private List<Region> combine(List<List<Region>> regions, int w, int h) {

        int[] xy = new int[2];

        Set<PairInt> centers = new HashSet<PairInt>();

        List<Region> combined = new ArrayList<Region>();
        for (List<Region> list : regions) {
            for (Region r : list) {
                r.calculateXYCentroid(xy, w, h);
                PairInt p = new PairInt((int)Math.round(xy[0]),
                    (int)Math.round(xy[1]));
                if (centers.contains(p)) {
                    continue;
                }
                centers.add(p);
                combined.add(r);
            }
        }

        //TODO: consider a filter for a minimum separation
        return combined;
    }

    private List<Region> createCombinedMSERRegions(GreyscaleImage gsImg,
        GreyscaleImage luvTheta) {

        MSER mser = new MSER();

        List<List<Region>> regions = mser.findRegions(gsImg, true);

        List<List<Region>> regionsT = mser.findRegions(luvTheta);

        List<Region> combined = new ArrayList<Region>();
        for (int i = 0; i < regions.size(); ++i) {
            for (Region r : regions.get(i)) {
                combined.add(r);
            }
        }
        for (int i = 0; i < regionsT.size(); ++i) {
            for (Region r : regionsT.get(i)) {
                combined.add(r);
            }
        }

        filterCloseToBounds2(combined, gsImg.getWidth(), gsImg.getHeight(), 10);

        return combined;
    }

    private void mergeRegionsAndSegmentation(List<Set<PairInt>> labeledSets, 
        TIntObjectMap<RegionPoints> cRegions, GreyscaleImage gsImg, 
        GreyscaleImage luvTheta, boolean replaceWithLabels) {

        TObjectIntMap<PairInt> pointLabelMap = new TObjectIntHashMap<PairInt>();
        for (int i = 0; i < labeledSets.size(); ++i) {
            for (PairInt p : labeledSets.get(i)) {
                pointLabelMap.put(p, i);
            }
        }
        
        MiscellaneousCurveHelper ch = new MiscellaneousCurveHelper();
        
        TIntSet skipRegions = new TIntHashSet();
        
        TIntObjectIterator<RegionPoints> iter = cRegions.iterator();
        for (int i = 0; i < cRegions.size(); ++i) {
            iter.advance();
            
            int rIdx = iter.key();
            RegionPoints cr = iter.value();
            
            int nLabeled = 0;
           
            // for the purpose of removing a labeled region that does not
            // belong in the points, need to store the points keys
            // associated w/ a label for each cRegion in order to remove them
            TIntObjectMap<Set<PairInt>> labelPointKeys 
                = new TIntObjectHashMap<Set<PairInt>>();
            
            // only used if replace is true
            TIntSet rmLabel = new TIntHashSet();
            
            //key=label, nPoints in label
            TIntIntMap labelInRegion = new TIntIntHashMap();
            for (PairInt p : cr.points) {
                
                if (!pointLabelMap.containsKey(p)) {
                    continue;
                }
                
                int label = pointLabelMap.get(p);
                if (labelInRegion.containsKey(label)) {
                    labelInRegion.put(label, labelInRegion.get(label) + 1);
                } else {
                    labelInRegion.put(label, 1);
                }
                
                Set<PairInt> pointKeys = labelPointKeys.get(label);
                if (pointKeys == null) {
                    pointKeys = new HashSet<PairInt>();
                    labelPointKeys.put(label, pointKeys);
                }
                pointKeys.add(p);
                
                nLabeled++;
            }
            
            // calculate percentage of cr covered by label
            //   and percentage of label inside cr
            //   and percentage of cr unassigned.
        
            float nInRegion = cr.points.size();
            
            float fUnassigned = (nInRegion - (float)nLabeled)/nInRegion;
            if (fUnassigned > 0.3333) {
                // missing a labeled region which was removed
                System.out.println("fraction unassigned=" + fUnassigned);
                skipRegions.add(rIdx);
                break;
            }
            
            //Image tmp = gsImg.copyToColorGreyscale();
            
            // for each label within, if any is a large portion of 
            //   cr, and yet a small portion of its labeled region,
            //   exclude this region
            TIntIntIterator iter2 = labelInRegion.iterator();
            for (int j = 0; j < labelInRegion.size(); ++j) {
                iter2.advance();
                int label = iter2.key();
                int count = iter2.value();
                
                int nTotInLabel = labeledSets.get(label).size();
                
                float fOutside = 1.f - ((float)count/(float)nTotInLabel);
                
                float fRegion = (float)count/nInRegion;
                
                PairInt labelXY = ch.calculateXYCentroids2(labeledSets.get(label));
                
                /*
                System.out.println("frcReg=" + fRegion +
                    " frcLblOut=" + fOutside +
                    " frcUnasnd=" + fUnassigned + 
                    " n=" + (int)nInRegion +
                    " cr.x,y=" + cr.xC + "," + cr.yC + 
                    " lbl.xy=" + labelXY + " rIdx=" + rIdx);
                */
                
                //int[] clr = ImageIOHelper.getNextRGB(j);
                //Set<PairInt> set = labeledSets.get(label);
                //ImageIOHelper.addCurveToImage(set, tmp, 0, clr[0], clr[1], clr[2]);
               
                //TODO: this needs revision...it is resolution (and scale) sensitive
                if ((nInRegion < 100 && fOutside > 0.45) || 
                    (nInRegion >= 100 && fOutside > 0.16)) {
                    System.out.println("removing label cen=" + labelXY);
                    // remove the offset points from the cRegion's offsets
                    for (PairInt rm : labelPointKeys.get(label)) {
                        cr.points.remove(rm);
                    }
                    rmLabel.add(label);
                } else {
                   // int[] clr = ImageIOHelper.getNextRGB(j);
                   // Set<PairInt> set = labeledSets.get(label);
                   // ImageIOHelper.addCurveToImage(set, tmp, 0, clr[0], clr[1], clr[2]);
                }
                
                //TODO: if not replacing, consider moving boundary inward to nearest segmentation
                //   bounds and then trimming labels external to the new
                //   moved boundary
                
            }
            
            if (cr.points.isEmpty()) {
                //System.out.println("removing empty region " + rIdx);
                skipRegions.add(rIdx);
                continue;
            }
            
            if (replaceWithLabels) {
                cr.points.clear();
                iter2 = labelInRegion.iterator();
                for (int j = 0; j < labelInRegion.size(); ++j) {
                    iter2.advance();
                    int label = iter2.key();
                    if (rmLabel.contains(label)) {
                        continue;
                    }
                    cr.points.addAll(labeledSets.get(label));
                }
                //TODO: consider re-doing x,y center
            }
            
            /*
            String lbl = Integer.toString(rIdx);
            if (lbl.length() < 4) {
                lbl = "0" + lbl;
            } 
            MiscDebug.writeImage(tmp, "_" + lbl);
            */
            
        } // end loop over regions
        
        System.out.println("nRegions=" + cRegions.size() + " removing=" +
            skipRegions.size());
        
        // remove skipSet
        TIntIterator iter2 = skipRegions.iterator();
        while (iter2.hasNext()) {
            int rmIdx = iter2.next();
            cRegions.remove(rmIdx);
        }        
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

    public static class Settings {
        private boolean useLargerPyramid0 = false;
        private boolean useLargerPyramid1 = false;

        //TODO: refactor to use an enum to avoid inconsistent state
        private boolean useSmallObjectMethod = false;
        private boolean useShapeFinder = false;

        private boolean findVanishingPoints = false;

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

        ORB orb1 = new ORB(500);//ORB(500);
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
        ns = orb1.getScalesList().size();
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
            ImageIOHelper.addCurveToImage(kpSet, img11,
                1, 255, 0, 0);
            MiscDebug.writeImage(img11, "_filtered_1_" + ts);
        }

        float luvDeltaELimit = 40;//10;// between 10 and 30

        boolean filterForCIE = false;
        if (filterForCIE) {
            rmIndexesList = new ArrayList<TIntList>();

            // ---- calculate the template CIE LUV colors ----
            //      -- also calculate separately, the colors of the
            //         entire shape to compare results when using the
            //         colors as filters.
            Set<PairInt> set0 = new HashSet<PairInt>();
            for (PairInt p : orb0.getKeyPointListColMaj(0)) {
                Set<PairInt> points = imageProcessor.getNeighbors(img0, p);
                set0.add(p);
                set0.addAll(points);
            }
            GroupPixelCIELUV meanLUVKeypoints = new GroupPixelCIELUV(set0, img0);
            meanLUVKeypoints.calculateColors(set0, img0, 0, 0);

            CIEChromaticity cieC = new CIEChromaticity();

            ns = orb1.getScalesList().size();
            for (int i = 0; i < ns; ++i) {
                TIntList kp0 = orb1.getKeyPoint0List().get(i);
                TIntList kp1 = orb1.getKeyPoint1List().get(i);
                TDoubleList or = orb1.getOrientationsList().get(i);
                TFloatList s = orb1.getScalesList().get(i);

                int np = kp0.size();
                TIntList rm = new TIntArrayList();
                for (int j = 0; j < np; ++j) {

                    PairInt p = new PairInt(kp1.get(j), kp0.get(j));
                    Set<PairInt> points = imageProcessor.getNeighbors(img1, p);
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
            }

            orb1.removeAtIndexes(rmIndexesList);

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
                ImageIOHelper.addCurveToImage(kpSet, img11,
                    1, 255, 0, 0);
                MiscDebug.writeImage(img11, "_filtered_2_" + ts);
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
            MiscDebug.writeImage(img11, "_segmented_3_" + ts);
        }

        List<Set<PairInt>> listOfPointSets2
            = LabelToColorHelper.extractContiguousLabelPoints(img1Cp, labels4);

        VanishingPoints vp2 = null;
        /*if (settings.isFindVanishingPoints()) {
            // vanishing points are used for attempts to account for
            // foreshortening of object dimensions from euclidean values,
            // and this correction is used at the end of
            // partial shape matching

            //NOTE that the segmentation is good for removing noise and vegetation,
            // etc, but sometimes also has blended in office windows to the
            // building by this point, so may need to use segmentation
            // before filtering

            vp2 = new VanishingPoints();
            //vp2.setToDebug();
            vp2.find(listOfPointSets2, img1.getWidth(),
                img1.getHeight());
            {// DEBUG
                Image img11 = img1.copyImage();
                MiscDebug.writeImage(img11, "_lines_");
            }

            throw new UnsupportedOperationException("not yet implemented");
        }*/

        boolean useCHist = false;

        boolean changed = false;

        changed = imageSegmentation.filterByCIECH(imgs0[0],
            shape0, img1Cp, listOfPointSets2, 0.1f);//0.4f);//0.35f
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

        if (debug) {
            ImageExt img11 = img1.copyToImageExt();
            ImageIOHelper.addAlternatingColorCurvesToImage0(listOfPointSets2,
                img11, 1);
            //ImageIOHelper.addAlternatingColorLabelsToRegion(
            //    img11, labels4);
            MiscDebug.writeImage(img11, "_segmented_2_" + ts);
        }

        // ---- remove keypoints if not in segmented cells
        if (changed && !settings.isUseSmallObjectMethod()
            && !settings.isUseShapeFinder()) {
            TObjectIntMap<PairInt> pointIdxMap
                = new TObjectIntHashMap<PairInt>();
            for (int j = 0; j < listOfPointSets2.size(); ++j) {
                Set<PairInt> set = listOfPointSets2.get(j);
                for (PairInt p : set) {
                    pointIdxMap.put(p, j);
                }
            }

            ns = orb1.getScalesList().size();
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
                //System.out.println("rm at scale " + i + " n=" +
                //    rm.size());
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

        // filter extremes of number density to try to remove
        //   labeled regions that are pure vegetation for example,
        //   when the object is not.  the large number of matches
        //   from features like vegetation points tends to result
        //   in high fraction terms, random good matches to shape,
        //   and roughly equiv matches in hsv descriptors so is often lowest cost
        //TODO: a sigma factor of 2.5 or larger is more safely above
        // effects of resolution and noise,
        // but this is set lower while testing.
        // this factor needs to be revisited one day or the method
        //    replaced with the texture filter which includes other
        //    information such as color too.
        filterBySurfaceDensity(orb0, orb1, listOfPointSets2, 1.25f);

        sList = new TFloatArrayList(orb1.getScalesList().size());
        for (int i = 0; i < orb1.getScalesList().size(); ++i) {
            sList.add(orb1.getScalesList().get(i).get(0));
        }

        List<CorrespondenceList> corList;

        List<Set<PairInt>> tempListOfPointSets
            = new ArrayList<Set<PairInt>>();
        tempListOfPointSets.add(shape0);

        if (debug) {
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
                "_filtered_segmentation_3_" + ts);

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
                "_filtered_segmentation_4_" + ts);
        }

        //TODO: might need to revise intersection limit in
        // these matching methods ... untested changes

        ORBMatcher orbMatcher = new ORBMatcher();
        if (vp2 != null) {
            orbMatcher.setVanishingPointsForSet2(vp2);
        }

        if (settings.isUseSmallObjectMethod()) {
            corList = orbMatcher.matchSmall(orb0, orb1,
                shape0, listOfPointSets2);
        } else if (settings.isUseShapeFinder()) {
            corList = orbMatcher.matchAggregatedShape(orb0, orb1,
                shape0, listOfPointSets2);
        } else {
            orb0.createDescriptorsHSV(imgs0[0]);
            orb1.createDescriptorsHSV(img1);

            boolean followWithShapeSearch = false;

            corList = orbMatcher.match0(orb0, orb1,
                shape0, listOfPointSets2);
        }

        if (corList == null || corList.isEmpty()) {
            return null;
        }

        return corList.get(0);
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
     * This method uses the blob detector MSER and matches descriptors for the
     * regions.
     *
     * @param img0
     * @param shape0
     * @param img1
     * @param settings for the method
     * @return
     */
    public CorrespondenceList findObject10(ImageExt img0, Set<PairInt> shape0,
        ImageExt img1, Settings settings) {

        long ts = 0;
        if (debug) {
            ts = MiscDebug.getCurrentTimeFormatted();
        }

        Canonicalizer canonicalizer = new Canonicalizer();

        GreyscaleImage gsImg0 = img0.copyToGreyscale2();

        // ----- create the cRegions for a masked image pyramid of img 0 ====
        MSER mser = new MSER();

        List<List<Region>> regions0 = mser.findRegions(gsImg0);

        int[] xy = new int[2];
        //remove all regions with centers outside of shape0 points
        for (int rIdx = 0; rIdx < 2; rIdx++) {
            java.util.List<Region> rList = regions0.get(rIdx);
            for (int i = (rList.size() - 1); i > -1; --i) {
                Region r = rList.get(i);
                r.calculateXYCentroid(xy, gsImg0.getWidth(), gsImg0.getHeight());
                PairInt p = new PairInt(xy[0], xy[1]);
                if (!shape0.contains(p)) {
                    rList.remove(i);
                }
            }
        }

        // mask the image, though not necessary
        for (int i = 0; i < gsImg0.getNPixels(); ++i) {
            PairInt p = new PairInt(gsImg0.getCol(i), gsImg0.getRow(i));
            if (!shape0.contains(p)) {
                gsImg0.setValue(i, 0);
            }
        }

        // find regions in dataset 1
        GreyscaleImage gsImg1 = img1.copyToGreyscale2();

        List<List<Region>> regions1 = mser.findRegions(gsImg1);

        filterCloseToBounds(regions1, gsImg1.getWidth(), gsImg1.getHeight(),
            10);

        //TODO: add filter here for patterns in the MSER regions that
        // are strong, and if present in reference frame1, then
        // anything without it should be removed.
        // use of this feature should be a Setting option.

        // filter regions by segmentation
        if (true) {
              List<Set<PairInt>> labeledSets = filterBySegmentation(img0, shape0,
                 img1, regions1, settings);
        }

        ImageProcessor imageProcessor = new ImageProcessor();

        List<List<GreyscaleImage>> pyr0 = imageProcessor
            .buildColorPyramid(img0,
            settings.useLargerPyramid0);

        List<List<GreyscaleImage>> pyr1 = imageProcessor
            .buildColorPyramid(img1,
            settings.useLargerPyramid0);

        // create image0 pyramid and windowing
        SummedAreaTable sumTable = new SummedAreaTable();
        int halfDimension = 1;

        for (int i = 0; i < pyr0.size(); ++i) {
            List<GreyscaleImage> imgMs = pyr0.get(i);
            for (int j = 0; j < imgMs.size(); ++j) {
                GreyscaleImage imgM = imgMs.get(j);
                imgM = sumTable.createAbsoluteSummedAreaTable(imgM);
                imgM = sumTable.applyMeanOfWindowFromSummedAreaTable(imgM,
                    2 * halfDimension + 1);
                imgMs.set(j, imgM);
            }
        }

        for (int i = 0; i < pyr1.size(); ++i) {
            List<GreyscaleImage> imgMs = pyr1.get(i);
            for (int j = 0; j < imgMs.size(); ++j) {
                GreyscaleImage imgM = imgMs.get(j);
                imgM = sumTable.createAbsoluteSummedAreaTable(imgM);
                imgM = sumTable.applyMeanOfWindowFromSummedAreaTable(imgM,
                    2 * halfDimension + 1);
                imgMs.set(j, imgM);
            }
        }

        List<GreyscaleImage> green0 = new ArrayList<GreyscaleImage>();
        for (List<GreyscaleImage> rgb : pyr0) {
            green0.add(rgb.get(1));
        }

        // list item for each octave.  item is for each region

        //regions[0) are found from the image,
        // while regions[1) are found from the inverted image.
        // they're kept separate for now to test matching them separately.
        List<TIntObjectMap<CRegion>> cRegionsList00 =
            canonicalizer.canonicalizeRegions(regions0.get(0), green0);

        List<TIntObjectMap<CRegion>> cRegionsList01 =
            canonicalizer.canonicalizeRegions(regions0.get(1), green0);

        // ---- create the cregion lists for image 1 ------

        if (debug) {
            int[] xyCen = new int[2];
            Image im0Cp, im1Cp;
            im0Cp = img0.copyImage();
            for (int i = 0; i < regions0.get(0).size(); ++i) {
                Region r = regions0.get(0).get(i);
                int[] clr = ImageIOHelper.getNextRGB(i);
                r.drawEllipse(im0Cp, 0, clr[0], clr[1], clr[2]);
                r.calculateXYCentroid(xyCen, im0Cp.getWidth(), im0Cp.getHeight());
                ImageIOHelper.addPointToImage(xyCen[0], xyCen[1], im0Cp,
                    1, 255, 0, 0);
            }
            MiscDebug.writeImage(im0Cp, "regions_0_0");
            im0Cp = img0.copyImage();
            for (int i = 0; i < regions0.get(1).size(); ++i) {
                Region r = regions0.get(1).get(i);
                int[] clr = ImageIOHelper.getNextRGB(i);
                r.drawEllipse(im0Cp, 0, clr[0], clr[1], clr[2]);
                r.calculateXYCentroid(xyCen, im0Cp.getWidth(), im0Cp.getHeight());
                ImageIOHelper.addPointToImage(xyCen[0], xyCen[1], im0Cp,
                    1, 255, 0, 0);
            }
            MiscDebug.writeImage(im0Cp, "regions_0_1");

            im1Cp = img1.copyImage();
            for (int i = 0; i < regions1.get(0).size(); ++i) {
                Region r = regions1.get(0).get(i);
                int[] clr = ImageIOHelper.getNextRGB(i);
                r.drawEllipse(im1Cp, 0, clr[0], clr[1], clr[2]);
                r.calculateXYCentroid(xyCen, im1Cp.getWidth(), im1Cp.getHeight());
                ImageIOHelper.addPointToImage(
                    xyCen[0], xyCen[1], im1Cp,
                    1, 255, 0, 0);
            }
            MiscDebug.writeImage(im1Cp, "regions_1_0");
            im1Cp = img1.copyImage();
            for (int i = 0; i < regions1.get(1).size(); ++i) {
                Region r = regions1.get(1).get(i);
                int[] clr = ImageIOHelper.getNextRGB(i);
                r.drawEllipse(im1Cp, 0, clr[0], clr[1], clr[2]);
                r.calculateXYCentroid(xyCen, im1Cp.getWidth(), im1Cp.getHeight());
                ImageIOHelper.addPointToImage(xyCen[0], xyCen[1], im1Cp,
                    1, 255, 0, 0);
            }
            MiscDebug.writeImage(im1Cp, "regions_1_1");
        }

        // mask the image, though not necessary
        for (int i = 0; i < gsImg0.getNPixels(); ++i) {
            PairInt p = new PairInt(gsImg0.getCol(i), gsImg0.getRow(i));
            if (!shape0.contains(p)) {
                gsImg0.setValue(i, 0);
            }
        }

        List<GreyscaleImage> green1 = new ArrayList<GreyscaleImage>();
        for (List<GreyscaleImage> rgb : pyr1) {
            green1.add(rgb.get(1));
        }

        List<TIntObjectMap<CRegion>> cRegionsList10 =
            canonicalizer.canonicalizeRegions(regions1.get(0), green1);

        List<TIntObjectMap<CRegion>> cRegionsList11 =
            canonicalizer.canonicalizeRegions(regions1.get(1), green1);

        boolean applyAutoCorrFilter = true;
        if (applyAutoCorrFilter) {
            double minAutoCor = findMinAutoCorrel(cRegionsList01.get(0));
            removeSmallAutoCorrel(cRegionsList11, minAutoCor, 0.75f);
            removeSmallAutoCorrel(cRegionsList10, minAutoCor, 0.75f);
        }

        if (debug) {
            debugPrint(cRegionsList00, pyr0, "_cr_0_0_");
            debugPrint(cRegionsList01, pyr0, "_cr_0_1_");
            debugPrint(cRegionsList10, pyr1, "_cr_1_0_");
            debugPrint(cRegionsList11, pyr1, "_cr_1_1_");
        }

        MSERMatcher matcher = new MSERMatcher();

        if (debug) {
            matcher.setToDebug();
        }

        CorrespondenceList corList = matcher.matchObject(
            pyr0, pyr1,
            cRegionsList00, cRegionsList01,
            cRegionsList10, cRegionsList11);

        return corList;
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
     * This method uses the blob detector MSER and matches descriptors for the
     * regions.
     *
     * It is different from findObject10 in several ways:
     * -- instead of providing MSERMatcher with only the
     * filtered canonicalized regions and image pyramids,
     * it also creates segmentation for the shape0
     * and passes both the labeled segmentation sets to the
     * MSERMatcher.
     *
     * @param img0
     * @param shape0
     * @param img1
     * @param settings for the method
     * @return
     */
    public CorrespondenceList findObject11(ImageExt img0, Set<PairInt> shape0,
        ImageExt img1, Settings settings) {

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

        GreyscaleImage luvTheta0 = imageProcessor.createCIELUVTheta(img0Trimmed, 255);
        GreyscaleImage luvTheta1 = imageProcessor.createCIELUVTheta(img1, 255);
        imageProcessor.singlePixelFilter(luvTheta0);
        imageProcessor.singlePixelFilter(luvTheta1);

        GreyscaleImage gsImg0 = img0Trimmed.copyToGreyscale2();
        GreyscaleImage gsImg1 = img1.copyToGreyscale2();

        // ----- create the cRegions for a masked image pyramid of img 0 ====

        // build combined list of regions
        List<Region> regionsComb0 = createCombinedMSERRegions(gsImg0, luvTheta0);

        mask(img0Trimmed, shape0Trimmed);
        mask(luvTheta0, shape0Trimmed);

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

        if (debug) {
            MiscDebug.writeImage(img0Trimmed, "_shape0_mask_");
        }

        List<Region> regionsComb1 = createCombinedMSERRegions(gsImg1, luvTheta1);

        // these are in the trimmed reference frame
        //NOTE: need segmentation for shape0 in order to find partial
        //   object when composed of very different parts, such as the
        //   cupcake cup and frosting, so don't condense to 1 label composed of shape0
        TIntObjectMap<Set<PairInt>> labeledSets0 = segment(img0Trimmed, shape0Trimmed);
        List<Set<PairInt>> labeledSets0List = new ArrayList<Set<PairInt>>();
        TIntObjectIterator<Set<PairInt>> iter = labeledSets0.iterator();
        for (int i = 0; i < labeledSets0.size(); ++i) {
            iter.advance();
            labeledSets0List.add(iter.value());
        }

        // filter regions by segmentation
        List<Set<PairInt>> labeledSets1 = filterBySegmentation2(
            img0Trimmed, shape0Trimmed, img1, regionsComb1, settings);
                
        List<List<GreyscaleImage>> pyrRGB0 = imageProcessor
            .buildColorPyramid(img0Trimmed,
            settings.useLargerPyramid0);
        
        List<GreyscaleImage> pyrPT0 = imageProcessor
            .buildPyramid(luvTheta0, settings.useLargerPyramid0);

        List<List<GreyscaleImage>> pyrRGB1 = imageProcessor
            .buildColorPyramid(img1,
            settings.useLargerPyramid0);

        List<GreyscaleImage> pyrPT1 = imageProcessor
            .buildPyramid(luvTheta1, settings.useLargerPyramid0);
       
        applyWindowedMean(pyrRGB0, 1);
        applyWindowedMean(pyrRGB1, 1);
        applyWindowedMean2(pyrPT0, 1);
        applyWindowedMean2(pyrPT1, 1);
        
        Canonicalizer canonicalizer = new Canonicalizer();

        // ----- create the cRegions for a masked image pyramid of img 0 ====

        //TODO: add filter here for patterns in the MSER regions that
        // are strong, and if present in reference frame1, then
        // anything without it should be removed.
        // use of this feature should be a Setting option.

        if (debug) {
            int[] xyCen = new int[2];
            Image im0Cp, im1Cp;
            im0Cp = img0Trimmed.copyImage();
            int n9 = regionsComb0.size();
            if (n9 > 10) {
            //    n9 = 10;
            }
            for (int i = 0; i < n9; ++i) {
                Region r = regionsComb0.get(i);
                int[] clr = ImageIOHelper.getNextRGB(i);
                r.drawEllipse(im0Cp, 0, clr[0], clr[1], clr[2]);
                r.calculateXYCentroid(xyCen, im0Cp.getWidth(), im0Cp.getHeight());
                ImageIOHelper.addPointToImage(xyCen[0], xyCen[1], im0Cp,
                    1, 255, 0, 0);
            }
            MiscDebug.writeImage(im0Cp, "regions_0_");

            im1Cp = img1.copyImage();
            n9 = regionsComb1.size();
            if (n9 > 10) {
            //    n9 = 10;
            }
            for (int i = 0; i < n9; ++i) {
                Region r = regionsComb1.get(i);
                int[] clr = ImageIOHelper.getNextRGB(i);
                r.drawEllipse(im1Cp, 0, clr[0], clr[1], clr[2]);
                r.calculateXYCentroid(xyCen, im1Cp.getWidth(), im1Cp.getHeight());
                ImageIOHelper.addPointToImage(xyCen[0], xyCen[1], im1Cp,
                    1, 255, 0, 0);
                System.out.println("regIdx1=" + i + " x="+xyCen[0] + " y=" + xyCen[1]);
            }
            MiscDebug.writeImage(im1Cp, "regions_1_");
        }

        //regions[0) are found from the image,
        // while regions[1) are found from the inverted image.
        // they're kept separate for now to test matching them separately.
        TIntObjectMap<RegionPoints> regionPoints0 =
            canonicalizer.canonicalizeRegions2(regionsComb0, pyrRGB0.get(0).get(1));

        TIntObjectMap<RegionPoints> regionPoints1 =
            canonicalizer.canonicalizeRegions2(regionsComb1, pyrRGB1.get(0).get(1));
   
        // these have been filtered to remove excess labeled points where possible.
        // and if last variable is true, the points have been 
        // replaced with filtered segmented label points instead of the ellipse
        mergeRegionsAndSegmentation(labeledSets0List, regionPoints0, gsImg0, luvTheta0,
            true);
        
        mergeRegionsAndSegmentation(labeledSets1, regionPoints1, gsImg1, luvTheta1, true);
        
        TIntObjectMap<CRegion> cRegions0 = canonicalizer.canonicalizeRegions3(
            regionPoints0, pyrRGB0.get(0).get(1));
        
        TIntObjectMap<CRegion> cRegions1 = canonicalizer.canonicalizeRegions3(
            regionPoints1, pyrRGB1.get(0).get(1));
      
        if (debug) {
            debugPrint2(cRegions0, pyrRGB0.get(0), "_cr_0_");
            //debugPrint2(cRegions1, pyrRGB1.get(0), "_cr_1_");
            debugPrint2(cRegions1, luvTheta1, "_cr_1_");
            MiscDebug.writeImage(luvTheta1, "_polar_cieluv_");
        }
        
        MSERMatcher matcher = new MSERMatcher();

        if (debug) {
            matcher.setToDebug();
        }

        List<CorrespondenceList> corList = matcher.matchObject2(
            pyrRGB0, pyrPT0, cRegions0, labeledSets0List,
            pyrRGB1,  pyrPT1, cRegions1, labeledSets1);

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

    private int[] getKeypointsMinMaxXY(ORB orb, int octave) {
        int[] minMaxXY1 = new int[4];
        minMaxXY1[0] = Integer.MAX_VALUE;
        minMaxXY1[2] = Integer.MAX_VALUE;
        minMaxXY1[1] = Integer.MIN_VALUE;
        minMaxXY1[3] = Integer.MIN_VALUE;
        for (int i = 0; i < orb.getKeyPoint0List().get(octave).size(); ++i) {
            int x = orb.getKeyPoint1List().get(octave).get(i);
            int y = orb.getKeyPoint0List().get(octave).get(i);
            if (x < minMaxXY1[0]) {
                minMaxXY1[0] = x;
            }
            if (x > minMaxXY1[1]) {
                minMaxXY1[1] = x;
            }
            if (y < minMaxXY1[2]) {
                minMaxXY1[2] = x;
            }
            if (y > minMaxXY1[3]) {
                minMaxXY1[3] = y;
            }
        }

        return minMaxXY1;
    }

    private void filterBySurfaceDensity(ORB orb1, ORB orb2,
        List<Set<PairInt>> listOfPointSets2, float sigmaFactorAboveMean) {

        Set<PairIntWithIndex> points1 = new HashSet<PairIntWithIndex>();
        int[] minMaxXY1 = getKeypointsMinMaxXY(orb1, 0);

        for (int i = 0; i < orb1.getKeyPoint0List().get(0).size(); ++i) {
            int x = orb1.getKeyPoint1List().get(0).get(i) - minMaxXY1[0];
            int y = orb1.getKeyPoint0List().get(0).get(i) - minMaxXY1[2];
            points1.add(new PairIntWithIndex(x, y, i));
        }

        DTClusterFinder<PairIntWithIndex> cFinder
            = new DTClusterFinder<PairIntWithIndex>(points1,
            minMaxXY1[1] - minMaxXY1[0] + 1,
            minMaxXY1[3] - minMaxXY1[2] + 1);
        cFinder.setToDebug();
        cFinder.calculateCriticalDensity();
        cFinder.findClusters();
        int n = cFinder.getNumberOfClusters();

        if (n == 0) {
            return;
        }

        float cd = cFinder.getCriticalDensity();

        //StringBuilder sb = new StringBuilder();
        //sb.append(String.format("template cd=%.2f n=%d\n", cd, n));
        TFloatList densities1 = new TFloatArrayList();
        TIntList nInGroup = new TIntArrayList();
        for (int i = 0; i < n; ++i) {

            Set<PairIntWithIndex> cluster = cFinder.getCluster(i);

            assert(cluster.size() >= 3);

            double dist = separationOfFurthestPair(cluster);
            float dens = (float)((float)cluster.size()/(dist*dist));
            densities1.add(dens);
            nInGroup.add(cluster.size());

        //    sb.append(String.format("  n[%d]=%d dens=%.2f", i, cluster.size(), dens));
        }
        //System.out.println(sb.toString());

        // average, stdDv, min, max
        float[] densityStats = calcStats(densities1, nInGroup);
        //System.out.printf("  dens avg=%.4f, stdv=%.4f, min=%.2f, max=%.2f\n",
        //    densityStats[0], densityStats[1], densityStats[2],
        //    densityStats[3]);

        double sz1 = Math.max(minMaxXY1[1] - minMaxXY1[0] + 1,
            minMaxXY1[3] - minMaxXY1[2] + 1);

        TObjectIntMap<PairInt> pointLabelMap = new TObjectIntHashMap<PairInt>();
        for (int i = 0; i < listOfPointSets2.size(); ++i) {
            Set<PairInt> set2 = listOfPointSets2.get(i);
            for (PairInt p2 : set2) {
                pointLabelMap.put(p2, i);
            }
        }

        TIntObjectMap<Set<PairInt>> _labelKeypointsMap =
            new TIntObjectHashMap<Set<PairInt>>();
        for (int i = 0; i < orb2.getKeyPoint0List().get(0).size(); ++i) {
            int x = orb2.getKeyPoint1List().get(0).get(i);
            int y = orb2.getKeyPoint0List().get(0).get(i);
            PairInt p0 = new PairInt(x, y);
            if (!pointLabelMap.containsKey(p0)) {
                continue;
            }
            int label2 = pointLabelMap.get(p0);
            Set<PairInt> set2 = _labelKeypointsMap.get(label2);
            if (set2 == null) {
                set2 = new HashSet<PairInt>();
                _labelKeypointsMap.put(label2, set2);
            }
            set2.add(p0);
        }

        TIntObjectMap<OneDIntArray> labelMinMaxXY = new
            TIntObjectHashMap<OneDIntArray>();
        TIntObjectMap<Set<PairIntWithIndex>> labelKeypoints2Map =
            new TIntObjectHashMap<Set<PairIntWithIndex>>();
        TIntObjectIterator<Set<PairInt>> iter = _labelKeypointsMap.iterator();
        for (int i = 0; i < _labelKeypointsMap.size(); ++i) {
            iter.advance();
            int label2 = iter.key();
            Set<PairInt> set = iter.value();

            OneDIntArray minMaxXY = new OneDIntArray(MiscMath.findMinMaxXY(set));
            labelMinMaxXY.put(label2, minMaxXY);

            Set<PairIntWithIndex> set2 = new HashSet<PairIntWithIndex>();
            for (PairInt p : set) {
                int x = p.getX() - minMaxXY.a[0];
                int y = p.getY() - minMaxXY.a[2];
                PairIntWithIndex p2 = new PairIntWithIndex(x, y, i);
                set2.add(p2);
            }
            labelKeypoints2Map.put(label2, set2);
        }

        TIntList rm = new TIntArrayList();
        for (int label2 = 0; label2 < listOfPointSets2.size(); ++label2) {

            Set<PairIntWithIndex> set2 = labelKeypoints2Map.get(label2);
            if (set2 == null || set2.size() < 10) {
                continue;
            }

            OneDIntArray minMaxXY2 = labelMinMaxXY.get(label2);

            double sz2 = Math.max(minMaxXY2.a[1] - minMaxXY2.a[0] + 1,
                minMaxXY2.a[3] - minMaxXY2.a[2] + 1);

            // edit cd to reference frame sz2
            // find clusters in set2
            // add those with significantly higher density to rm

            float cd2;
            int type;
            if (sz1 > sz2 && ((sz1/sz2) > 1.2)) {
                // decrease the critical separation by increasing
                // critical density by factor

                // if the surface density is significantly larger,
                //   the region should be removed

                type = 0;
                cd2 = cd * (float)(sz1/sz2);

            } else if (sz2 > sz1 && ((sz2/sz1) > 1.2)) {
                // increase the critical separation by decreasing
                // critical density by factor

                //NOTE:
                // dataset "2" may have higher resolution, so removing the
                // region because of larger density might not be the correct

                type = 1;
                cd2 = cd / (float)(sz2/sz1);

            } else {
                // these are roughly equiv sized regions, so a direct
                // comparison of densities is needed.

                // if the density of 2 is significantly larger, the region
                // should be removed.

                type = 2;
                cd2 = cd;

            }

            DTClusterFinder<PairIntWithIndex> cFinder2
                = new DTClusterFinder<PairIntWithIndex>(set2,
                minMaxXY2.a[1] - minMaxXY2.a[0] + 1,
                minMaxXY2.a[3] - minMaxXY2.a[2] + 1);
            cFinder2.setCriticalDensity(cd2);
            cFinder2.findClusters();
            int n2 = cFinder2.getNumberOfClusters();

            TFloatList densities2 = new TFloatArrayList();
            TIntList nInGroups2 = new TIntArrayList();

            //sb = new StringBuilder();
            //sb.append(String.format("label2=%d cd2=%.2f n2=%d\n", label2, cd2, n2));
            for (int i = 0; i < n2; ++i) {

                Set<PairIntWithIndex> cluster2 = cFinder2.getCluster(i);

                double dist = separationOfFurthestPair(cluster2);

                if (type == 0 || type == 1) {
                    dist *= (sz2/sz1);
                }
                float dens = (float) ((float) cluster2.size() / (dist * dist));
                densities2.add(dens);
                nInGroups2.add(cluster2.size());

                //sb.append(String.format("  n2[%d]=%d dens=%.2f",
                //    i, cluster2.size(), dens));
            }
            //System.out.println(sb.toString());

            if (n2 > 0) {
                // average, stdDv, min, max
                float[] densityStats2 = calcStats(densities2, nInGroups2);
                //System.out.printf("  dens2 avg=%.4f, stdv=%.4f, min=%.2f, max=%.2f\n",
                //    densityStats2[0], densityStats2[1], densityStats2[2],
                //    densityStats2[3]);

                if (densityStats2[0] > densityStats[0]) {
                    float diff = densityStats2[0] - densityStats[0];
                    float sigmaFactor = diff/densityStats[1];

                    if (sigmaFactor > sigmaFactorAboveMean) {

                        rm.add(label2);

                        if (debug) {
                            MiscellaneousCurveHelper ch = new MiscellaneousCurveHelper();
                            PairInt xyCen = ch.calculateXYCentroids2(
                                listOfPointSets2.get(label2));
                            System.out.format(
                                "removing label2=%d which is %.2f sigma above mean.  coords=%s\n",
                                label2, sigmaFactor, xyCen.toString());
                        }
                    }
                }
            }
        }

        if (rm.isEmpty()) {
            return;
        }

        // remove keypoints in label2
        // and remove label2 from the point set list

        // a list of all points to search each octave2's keypoints for presence
        //    and remove that keypoints data if present
        Set<PairInt> rmAllPoints = new HashSet<PairInt>();
        for (int i = 0; i < rm.size(); ++i) {
            rmAllPoints.addAll(listOfPointSets2.get(rm.get(i)));
        }

        List<TIntList> rmList = new ArrayList<TIntList>();
        for (int octave2 = 0; octave2 < orb2.getKeyPoint1List().size();
            ++ octave2) {

            TIntList kp1 = orb2.getKeyPoint1List().get(octave2);
            TIntList kp0 = orb2.getKeyPoint0List().get(octave2);

            TIntList rm2 = new TIntArrayList();
            rmList.add(rm2);
            for (int i = 0; i < kp1.size(); ++i) {
                PairInt p = new PairInt(kp1.get(i), kp0.get(i));
                if (rmAllPoints.contains(p)) {
                    rm2.add(i);
                }
            }
        }
        orb2.removeAtIndexes(rmList);

        rm.sort();

        for (int i = (rm.size() - 1); i > -1; --i) {
            int rmIdx = rm.get(i);
            listOfPointSets2.remove(rmIdx);
        }

    }

    private double separationOfFurthestPair(Set<PairIntWithIndex> points) {

        if (points == null || points.size() < 2) {
            throw new IllegalArgumentException("points must have at least "
                + " 2 points");
        }

        FurthestPair fp = new FurthestPair();

        Set<PairInt> set = new HashSet<PairInt>(points.size());
        for (PairIntWithIndex p : points) {
            set.add(new PairInt(p.getX(), p.getY()));
        }

        PairInt[] furthest = fp.find(set);

        assert(furthest != null);
        assert(furthest.length > 1);

        int diffX = furthest[0].getX() - furthest[1].getX();
        int diffY = furthest[0].getY() - furthest[1].getY();

        double dist = Math.sqrt(diffX * diffX + diffY * diffY);

        return dist;
    }

    private float[] calcStats(TFloatList densities, TIntList nInGroup) {

        //average, stdDv, min, max

        // using the number in group to form weights
        float nSum = 0;
        for (int i = 0; i < densities.size(); ++i) {
            nSum += nInGroup.get(i);
        }

        TFloatList weights = new TFloatArrayList();
        float totW = 0;
        for (int i = 0; i < densities.size(); ++i) {
            float w = (float)nInGroup.get(i)/nSum;
            //System.out.println(" w=" + w + " n=" + nInGroup.get(i));
            weights.add(w);
            totW += w;
        }
        //System.out.println("totW=" + totW);
        assert(Math.abs(totW - 1.f) < 0.01f);

        double avgDens = 0;
        for (int i = 0; i < densities.size(); ++i) {
            float d = densities.get(i) * weights.get(i);
            avgDens += d;
        }

        double stDev = 0;
        for (int i = 0; i < densities.size(); ++i) {
            double diff = densities.get(i) - avgDens;
            stDev += (diff * diff);
        }
        stDev = Math.sqrt(stDev/((float) densities.size() - 1));

        //average, stdDv, min, max
        return new float[]{(float)avgDens, (float)stDev, densities.min(),
            densities.max()};
    }

    private List<Set<PairInt>> filterBySegmentation2(
        ImageExt img0, Set<PairInt> shape0,
        ImageExt img1, List<Region> regionsList, Settings settings) {

        // creating a fake 2nd item to be able to use existing method
        List<List<Region>> list = new ArrayList<List<Region>>();
        list.add(regionsList);

        return filterBySegmentation(img0, shape0, img1, list, settings);
    }

    private List<Set<PairInt>> filterBySegmentation(
        ImageExt img0, Set<PairInt> shape0,
        ImageExt img1, List<List<Region>> regionsList, Settings settings) {

        int[] xy = new int[2];
        Set<PairInt> keypoints1 = new HashSet<PairInt>();
        for (int i = 0; i < regionsList.size(); ++i) {
            List<Region> regions = regionsList.get(i);
            for (int j = 0; j < regions.size(); ++j) {
                Region r = regions.get(j);
                r.calculateXYCentroid(xy, img1.getWidth(), img1.getHeight());
                keypoints1.add(new PairInt(xy[0], xy[1]));
            }
        }

        int n = keypoints1.size();

        List<Set<PairInt>> setLists = filterBySegmentation(img0, shape0, img1, keypoints1, settings);

        if (keypoints1.size() == n) {
            return setLists;
        }

        for (int i = 0; i < regionsList.size(); ++i) {
            TIntList rm = new TIntArrayList();
            List<Region> regions = regionsList.get(i);
            for (int j = 0; j < regions.size(); ++j) {
                Region r = regions.get(j);
                r.calculateXYCentroid(xy, img1.getWidth(), img1.getHeight());
                PairInt p = new PairInt(xy[0], xy[1]);
                if (!keypoints1.contains(p)) {
                    rm.add(j);
                }
            }
            for (int j = (rm.size() - 1); j > -1; --j) {
                regions.remove(rm.get(j));
            }
        }

        return setLists;
    }

    private List<Set<PairInt>> filterBySegmentation(
        ImageExt img0, Set<PairInt> shape0,
        ImageExt img1, Set<PairInt> keypoints1, Settings settings) {

        //TODO: consider returning the keypoints from the full frames


        //TODO: since the ORBs are not kept, can revise this to
        // work on a single image instead of a pyramid.

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

        ORB orb1 = new ORB(500);//ORB(500);
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
        ns = orb1.getScalesList().size();
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

        // --- filter keypoints1
        Set<PairInt> rm2 = new HashSet<PairInt>();
        for (PairInt p : keypoints1) {
            Set<PairInt> points = imageProcessor.getNeighbors(img1, p);
            points.add(p);
            int[][] ch = clrHist.histogramHSV(img1, points);
            float intersection = clrHist.intersection(
                template_ch_HSV, ch);

            if (intersection < 0.1) {
                rm2.add(p);
            } else {
                ch = clrHist.histogramCIELAB(img1, points);
                intersection = clrHist.intersection(
                    template_ch_LAB, ch);
                if (intersection < 0.2) {
                    rm2.add(p);
                }
            }
        }
        keypoints1.removeAll(rm2);

        if (orb1.getKeyPoint0List().isEmpty()) {
            return null;
        }

        //float luvDeltaELimit = 40;//10;// between 10 and 30

        ImageExt img1Cp = img1.copyToImageExt();

        int[] labels4 = imageSegmentation.objectSegmentation(img1Cp);

        if (debug) {
            ImageExt img11 = img1.copyToImageExt();
            ImageIOHelper.addAlternatingColorLabelsToRegion(
                img11, labels4);
            MiscDebug.writeImage(img11, "_segmented_3_" + ts);
        }

        List<Set<PairInt>> listOfPointSets2
            = LabelToColorHelper.extractContiguousLabelPoints(img1Cp, labels4);

        //boolean changed = false;
        boolean changed = imageSegmentation.filterByCIECH(imgs0[0],
            shape0, img1Cp, listOfPointSets2,
            0.01f);
            //0.1f);//0.4f);//0.35f

        if (debug) {
            ImageExt img11 = img1.copyToImageExt();
            ImageIOHelper.addAlternatingColorCurvesToImage0(listOfPointSets2,
                img11, 1);
            //ImageIOHelper.addAlternatingColorLabelsToRegion(
            //    img11, labels4);
            MiscDebug.writeImage(img11, "_segmented_2_" + ts);
        }

        // ---- remove keypoints if not in segmented cells
        if (changed && !settings.isUseSmallObjectMethod()
            && !settings.isUseShapeFinder()) {
            TObjectIntMap<PairInt> pointIdxMap
                = new TObjectIntHashMap<PairInt>();
            for (int j = 0; j < listOfPointSets2.size(); ++j) {
                Set<PairInt> set = listOfPointSets2.get(j);
                for (PairInt p : set) {
                    pointIdxMap.put(p, j);
                }
            }

            ns = orb1.getScalesList().size();
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
                //System.out.println("rm at scale " + i + " n=" +
                //    rm.size());
            }
            orb1.removeAtIndexes(rmIndexesList);

            // --- filter out keypoints1
            rm2.clear();
            for (PairInt p : keypoints1) {
                if (!pointIdxMap.containsKey(p)) {
                    rm2.add(p);
                }
            }
            keypoints1.removeAll(rm2);
        }

        return listOfPointSets2;
    }
}
