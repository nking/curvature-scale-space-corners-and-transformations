package algorithms.imageProcessing.features;

import algorithms.QuickSort;
import algorithms.compGeometry.PerimeterFinder2;
import algorithms.imageProcessing.ColorHistogram;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageExt;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.ImageSegmentation;
import algorithms.imageProcessing.MiscellaneousCurveHelper;
import algorithms.imageProcessing.SummedAreaTable;
import algorithms.imageProcessing.TrimmedImage;
import algorithms.imageProcessing.features.mser.Canonicalizer;
import algorithms.imageProcessing.features.mser.Canonicalizer.CRegion;
import algorithms.imageProcessing.features.mser.Canonicalizer.RegionGeometry;
import algorithms.imageProcessing.features.mser.Canonicalizer.RegionPoints;
import algorithms.imageProcessing.features.mser.MSER;
import algorithms.imageProcessing.features.mser.MSER.Threshold;
import algorithms.imageProcessing.features.mser.Region;
import algorithms.imageProcessing.matching.MSERMatcher;
import algorithms.imageProcessing.segmentation.LabelToColorHelper;
import algorithms.imageProcessing.util.GroupAverageColors;
import algorithms.imageProcessing.util.PairIntWithIndex;
import algorithms.misc.Misc;
import algorithms.misc.MiscDebug;
import algorithms.misc.MiscMath;
import algorithms.util.PairInt;
import algorithms.util.PixelHelper;
import com.climbwithyourfeet.clustering.ClusterFinder;
import gnu.trove.iterator.TIntIntIterator;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.iterator.TLongIterator;
import gnu.trove.list.TFloatList;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.TLongIntMap;
import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.map.hash.TLongIntHashMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.TLongSet;
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

    private List<Region> combine(List<List<Region>> regions, int w, int h) {

        int[] xy = new int[2];

        Set<PairInt> centers = new HashSet<PairInt>();

        List<Region> combined = new ArrayList<Region>();
        for (List<Region> list : regions) {
            for (Region r : list) {
                r.calculateXYCentroid(xy, w, h);
                PairInt p = new PairInt(Math.round(xy[0]),
                    Math.round(xy[1]));
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
        
        System.out.println(debugLabel + "  clrMode=" + clrMode.name() 
            + " ptMode=" + ptMode.name());
        
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
                }
                MiscDebug.writeImage(imCp, debugLabel + "_regions_gs_"+ type + "_" + ts);
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
                MiscDebug.writeImage(imCp, debugLabel + "_regions_pt_"+ type + "_" + ts);
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
            
            //System.out.println("cr xy=" + cr.ellipseParams.xC + "," + 
            //    cr.ellipseParams.yC + ")");
            
            int nLabeled = 0;
           
            // for the purpose of removing a labeled region that does not
            // belong in the points, need to store the points keys
            // associated w/ a label for each cRegion in order to remove them
            TIntObjectMap<Set<PairInt>> labelPointKeys 
                = new TIntObjectHashMap<Set<PairInt>>();
            
            // only used if replace is true
            TIntSet rmLabel = new TIntHashSet();
            
            Set<PairInt> rmPoints = new HashSet<PairInt>();
            
            //key=label, nPoints in label
            TIntIntMap labelInRegion = new TIntIntHashMap();
            for (PairInt p : cr.points) {
                
                if (!pointLabelMap.containsKey(p)) {
                    rmPoints.add(p);
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
            cr.points.removeAll(rmPoints);
            
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
                    " cr.x,y=" + cr.ellipseParams.xC + "," + cr.ellipseParams.yC + 
                    " lbl.xy=" + labelXY + " rIdx=" + rIdx);
                */
                
                //int[] clr = ImageIOHelper.getNextRGB(j);
                //Set<PairInt> set = labeledSets.get(label);
                //ImageIOHelper.addCurveToImage(set, tmp, 0, clr[0], clr[1], clr[2]);
               
                //TODO: this needs revision...it is resolution (and scale) sensitive
                if (//(nInRegion < 100 && fOutside > 0.45) || 
                    (nInRegion >= 100 && fOutside > 0.16)) {
                    //System.out.println("removing label cen=" + labelXY + 
                    //    " nInRegion=" + nInRegion + " fOutside=" + fOutside + 
                    //    " from region(" + cr.ellipseParams.xC + "," +
                    //    cr.ellipseParams.yC + ")");
                    //// remove the offset points from the cRegion's offsets
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
            
            if (!rmLabel.isEmpty()) {
                TIntIterator iter0 = rmLabel.iterator();
                while (iter0.hasNext()) {
                    int label = iter0.next();
                    Set<PairInt> pts = labelPointKeys.get(label);
                    cr.points.removeAll(pts);
                    labelInRegion.remove(label);
                }
            }
            
            if (replaceWithLabels) {
                cr.points.clear();
                iter2 = labelInRegion.iterator();
                for (int j = 0; j < labelInRegion.size(); ++j) {
                    iter2.advance();
                    int label = iter2.key();
                    cr.points.addAll(labeledSets.get(label));
                }
                // re-calc center and orientation.  tests show this does not
                //  improve final matches
                /*
                Region tmp = new Region();
                for (PairInt pt : cr.points) {
                    tmp.accumulate(pt.getX(), pt.getY());
                }
                cr.ellipseParams = Canonicalizer.calculateEllipseParams(tmp,
                    gsImg.getWidth(), gsImg.getHeight());
                */
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

        //when multiple regions are centered within a spatial limit,
        //  choose one and remove the others
        if (false) {
            int critSep = 15;
            
            System.out.println("before removing near mser, cRegions.n=" + 
                cRegions.size());
            
            // clustering algorithm needs pixel indexes, but would like to
            // track the cRegions indexes too
            PixelHelper ph = new PixelHelper();
            
            // key = pixIdx, value = rIdx
            TLongIntMap pixRIdxMap = new TLongIntHashMap();
             
            iter = cRegions.iterator();
            for (int i = 0; i < cRegions.size(); ++i) {
                iter.advance();
                int rIdx = iter.key();
                RegionPoints cr = iter.value();
                
                long pixIdx = ph.toPixelIndex(
                    cr.ellipseParams.xC, cr.ellipseParams.yC, gsImg.getWidth());
                
                pixRIdxMap.put(pixIdx, rIdx);
            }
            
            ClusterFinder cFinder = new ClusterFinder(pixRIdxMap.keySet(),
                gsImg.getWidth(), gsImg.getHeight());
            cFinder.setThreshholdFactor(1.f);
            cFinder.setMinimumNumberInCluster(2);
            cFinder.setBackgroundSeparation(critSep, critSep);
            cFinder.findClusters();
            List<TLongSet> groupList = cFinder.getGroups();

            //NOTE: may need to revise how to choose best region to keep.
            for (int i = 0; i < groupList.size(); ++i) {
                TLongSet groupPixs = groupList.get(i);
                
                int maxSz = Integer.MIN_VALUE;
                int maxSzIdx = -1;                
                
                TLongIterator iter3 = groupPixs.iterator();
                while (iter3.hasNext()) {
                    long pixIdx = iter3.next();
                    int rIdx = pixRIdxMap.get(pixIdx);
                
                    int sz = calculateObjectSize(cRegions.get(rIdx));
                    if (sz > maxSz) {
                        maxSz = sz;
                        maxSzIdx = rIdx;
                    }
                }
                assert(maxSzIdx > -1);
                iter3 = groupPixs.iterator();
                while (iter3.hasNext()) {
                    long pixIdx = iter3.next();
                    int rIdx = pixRIdxMap.get(pixIdx);
                    if (rIdx == maxSzIdx) {
                        continue;
                    }
                    cRegions.remove(rIdx);
                }
            }
            
            System.out.println("after removing near mser, cRegions.n=" + 
                cRegions.size());
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
    
    private void createAWholeRegion(TIntObjectMap<RegionPoints> regionPoints, 
        Set<PairInt> shape, GreyscaleImage rgb) {
        
        int n = regionPoints.size();
        
        RegionPoints rp = createARegion(shape, rgb.getWidth(), rgb.getHeight());
      
        regionPoints.put(n, rp);
    }

    private void createForLabels(List<Set<PairInt>> labeledSets, 
        TIntObjectMap<CRegion> regions, int imageWidth, int imageHeight) {

        int ns = labeledSets.size();
        TIntSet labels = new TIntHashSet();
        for (int i = 0; i < ns; ++i) {
            labels.add(i);
        }
        
        int idxMax = Integer.MIN_VALUE;
        TIntObjectIterator<CRegion> iter = regions.iterator();
        for (int i = 0; i < regions.size(); ++i) {
            iter.advance();
            int rIdx = iter.key();
            CRegion r = iter.value();
        
            if (rIdx > idxMax) {
                idxMax = rIdx;
            }
            
            if (r.labels.size() == 1) {
                labels.remove(r.labels.iterator().next());
            }
        }
        
        // for labels, make CRegion structures and add to map
        TIntIterator iter2 = labels.iterator();
        while (iter2.hasNext()) {
            int label = iter2.next();
            
            Set<PairInt> set = labeledSets.get(label);
            
            RegionPoints rp = createARegion(set, imageWidth, imageHeight);
            
            RegionGeometry ep = rp.ellipseParams;
            
            Map<PairInt, PairInt> offsetToOrigMap = 
                Canonicalizer.createOffsetToOrigMap(
                ep.xC, ep.yC, Misc.convertWithoutOrder(set), 
                imageWidth, imageHeight, ep.orientation);

            CRegion cRegion = new CRegion();
            cRegion.ellipseParams = ep;
            cRegion.offsetsToOrigCoords = offsetToOrigMap;
               
            idxMax++;
            
            regions.put(idxMax, cRegion);
        }
    }

    private CMODE determineColorMode(ImageExt img, Set<PairInt> set) {

        GroupAverageColors clrs = new GroupAverageColors(img, set);
        
        int limit1 = 150;
        int limit2 = 55;
        if (clrs.getR() >= limit1 && clrs.getG() >= limit1 &&
            clrs.getB() >= limit1) {
            return CMODE.WHITE;
        } else if (clrs.getR() <= limit2 && clrs.getG() <= limit2 &&
            clrs.getB() <= limit2) {
            return CMODE.BLACK;
        } else {
            return CMODE.OTHER;
        }
    }

    private CMODE determinePolarThetaMode(GreyscaleImage luvTheta, 
        Set<PairInt> points) {
    
        double avg = 0;
        for (PairInt p : points) {
            avg += luvTheta.getValue(p);
        }
        avg /= (double)points.size();
        
        int limit1 = 220;
        int limit2 = 25;
        if (avg >= limit1) {
            return CMODE.WHITE;
        } else if (avg <= limit2) {
            return CMODE.BLACK;
        } else {
            return CMODE.OTHER;
        }
    }

    private int calculateObjectSize(RegionPoints region) {
        int[] minMaxXY = MiscMath.findMinMaxXY(region.points);
        int diffX = minMaxXY[1] - minMaxXY[0];
        int diffY = minMaxXY[3] - minMaxXY[2];
        double xy = Math.sqrt(diffX * diffX + diffY * diffY);
        
        return (int)Math.round(xy);
    }

    private void filterByColorHistograms(ImageExt img0, Set<PairInt> shape0, 
        ImageExt img1, TIntObjectMap<RegionPoints> regions1) {
        
        //filter by color hist of hsv, cielab and by CIECH

        ColorHistogram clrHist = new ColorHistogram();

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
            if (intersection < 0.33) {
                rmSet.add(rIdx);
            } else {
                ch = clrHist.histogramCIELAB(img1, r.points);
                intersection = clrHist.intersection(template_ch_LAB, ch);
                if (intersection < 0.33) {
                    rmSet.add(rIdx);
                } else {
                    int[] tHist1 = clrHist.histogramCIECH64(img1, r.points);
                    intersection = clrHist.intersection(tHist, tHist1);
                    if (intersection < 0.33f) {
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

    }

    /**
     * descriptions black, white, or other used in describing the template
     * shape color.  the extremes black and white can be used to limit
     * the regions created.
     */
    private enum CMODE {
        WHITE, BLACK, OTHER
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
    public CorrespondenceList findObject12(ImageExt img0, Set<PairInt> shape0,
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

        mask(img0Trimmed, shape0Trimmed);
        mask(luvTheta0, shape0Trimmed);

        CMODE clrMode = determineColorMode(img0Trimmed, shape0Trimmed);
        
        CMODE ptMode = determinePolarThetaMode(luvTheta0, shape0Trimmed);
        
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
            //MiscDebug.writeImage(luvTheta0, "_luv_mask_");
            
            //MiscDebug.writeImage(tmp00, "_gs_enhanced_0_");
            //MiscDebug.writeImage(tmp01, "_luv_enhanced_0_");
            //MiscDebug.writeImage(tmp10, "_gs_enhanced_1_");
            //MiscDebug.writeImage(tmp11, "_luv_enhanced_1_");
        }
        
        List<Region> regionsComb1 = createCombinedMSERRegions(
            tmp10, tmp11,
            //gsImg1, luvTheta1, 
            clrMode, ptMode, fewerMSER, settings.getDebugLabel() + "_1_");
                
        int critSep = 5;
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
       
       // applyWindowedMean(pyrRGB0, 1);
        //applyWindowedMean(pyrRGB1, 1);
        //applyWindowedMean2(pyrPT0, 1);
        //applyWindowedMean2(pyrPT1, 1);
        
        Canonicalizer canonicalizer = new Canonicalizer();

        // ----- create the cRegions for a masked image pyramid of img 0 ====

        //TODO: add filter here for patterns in the MSER regions that
        // are strong, and if present in reference frame1, then
        // anything without it mughr be removable.
        // use of this feature should be a Setting option.

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
                "_regions_0_");

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
                + "_regions_1_");
        }
        
        TIntObjectMap<RegionPoints> regionPoints0 =
            canonicalizer.canonicalizeRegions2(regionsComb0, pyrRGB0.get(0).get(1));
   
        TIntObjectMap<RegionPoints> regionPoints1 =
            canonicalizer.canonicalizeRegions2(regionsComb1, pyrRGB1.get(0).get(1));
  
        // filter by color hist of hsv, cielab and CIECH
        filterByColorHistograms(img0Trimmed, shape0Trimmed, img1, 
            regionPoints1);

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
        
        List<CorrespondenceList> corList 
            = matcher.matchObject0(
            pyrRGB0, pyrPT0, regionPoints0,
            pyrRGB1,  pyrPT1, regionPoints1,
            settings.getDebugLabel());
        
        if (corList == null) {
            return null;
        }
        
        // apply offsets for having trimmed image 0
        CorrespondenceList topC = corList.get(0);
        for (int i = 0; i < topC.getPoints1().size(); ++i) {
            PairInt p = topC.getPoints1().get(i);
            int x = p.getX();
            int y = p.getY();
            p.setX(x + img0Trim.getXOffset());
            p.setY(y + img0Trim.getYOffset());
        }

        return topC;
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

}
