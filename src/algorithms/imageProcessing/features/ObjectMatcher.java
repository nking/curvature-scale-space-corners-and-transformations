package algorithms.imageProcessing.features;

import algorithms.QuickSort;
import algorithms.compGeometry.FurthestPair;
import algorithms.compGeometry.NearestPoints;
import algorithms.compGeometry.PerimeterFinder2;
import algorithms.imageProcessing.ColorHistogram;
import algorithms.imageProcessing.FixedSizeSortedVector;
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
import algorithms.search.NearestNeighbor2D;
import algorithms.util.OneDIntArray;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
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
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.Stack;

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
            int[] gsA = mser.readIntoArray(gsImg);
            List<Region> list = mser.findRegionsNeg(gsA,
                gsImg.getWidth(), gsImg.getHeight(), thrGs);
            regions.add(new ArrayList<Region>());
            regions.add(list);
        } else if (clrMode.equals(CMODE.BLACK)) {
            int[] gsA = mser.readIntoArray(gsImg);
            List<Region> list = mser.findRegionsPos(gsA,
                gsImg.getWidth(), gsImg.getHeight(), thrGs);
            regions.add(list);
            regions.add(new ArrayList<Region>());
        } else {
            regions = mser.findRegions(gsImg, thrGs);
        }
        
        if (ptMode.equals(CMODE.WHITE)) {
            int[] ptA = mser.readIntoArray(luvTheta);
            List<Region> list = mser.findRegionsNeg(ptA,
                luvTheta.getWidth(), luvTheta.getHeight(), thrPt);
            regionsT.add(new ArrayList<Region>());
            regionsT.add(list);
        } else if (ptMode.equals(CMODE.BLACK)) {
            int[] ptA = mser.readIntoArray(luvTheta);
            List<Region> list = mser.findRegionsPos(ptA,
                luvTheta.getWidth(), luvTheta.getHeight(), thrPt);
            regionsT.add(list);
            regionsT.add(new ArrayList<Region>());
        } else {
            regionsT = mser.findRegions(luvTheta, thrPt);
        
            // filter to remove all w/ variation > 0
            for (int type = 0; type < 2; ++type) {
                List<Region> list = regionsT.get(type);
                for (int i = (list.size() - 1); i > -1; --i) {
                    Region r = list.get(i);
                    if ((type == 1) && r.getVariation() > 0.001) {
                        list.remove(i);
                    } else if ((type == 0) && r.getVariation() == 0.0) {
                        list.remove(i);
                    }
                }
            }
        }
        
        
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
            float critDens = 2.f/15.f;
            
            System.out.println("before removing near mser, cRegions.n=" + 
                cRegions.size());
            
            Set<PairIntWithIndex> points2
                = new HashSet<PairIntWithIndex>();
            
            iter = cRegions.iterator();
            for (int i = 0; i < cRegions.size(); ++i) {
                iter.advance();
                int rIdx = iter.key();
                RegionPoints cr = iter.value();
                PairIntWithIndex pii = new PairIntWithIndex(
                    cr.ellipseParams.xC, cr.ellipseParams.yC, rIdx);
                points2.add(pii);
            }
            
            DTClusterFinder<PairIntWithIndex> cFinder
                = new DTClusterFinder<PairIntWithIndex>(points2,
                gsImg.getWidth() + 1, gsImg.getHeight() + 1);
            cFinder.setMinimumNumberInCluster(2);
            cFinder.setCriticalDensity(critDens);
            cFinder.findClusters();

            //NOTE: may need to revise how to choose best region to keep.
            for (int i = 0; i < cFinder.getNumberOfClusters(); ++i) {
                Set<PairIntWithIndex> set = cFinder.getCluster(i);
                int maxSz = Integer.MIN_VALUE;
                int maxSzIdx = -1;
                
                for (PairIntWithIndex pii : set) {
                    int rIdx = pii.getPixIndex();
                    int sz = calculateObjectSize(cRegions.get(rIdx));
                    if (sz > maxSz) {
                        maxSz = sz;
                        maxSzIdx = rIdx;
                    }
                }
                assert(maxSzIdx > -1);
                for (PairIntWithIndex pii : set) {
                    int rIdx = pii.getPixIndex();
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
                
        float critDens = 2.f/10.f;
        Canonicalizer.filterBySpatialProximity(critDens, regionsComb0, 
            img0Trimmed.getWidth(), img0Trimmed.getHeight());
        
        Canonicalizer.filterBySpatialProximity(critDens, regionsComb1, 
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
                            //System.out.format(
                            //    "removing label2=%d which is %.2f sigma above mean.  coords=%s\n",
                            //    label2, sigmaFactor, xyCen.toString());
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
    
    private TObjectIntMap<PairInt> createPointLabelMap(
        List<Set<PairInt>> labeledSets) {

        TObjectIntMap<PairInt> pointLabelMap = new TObjectIntHashMap<PairInt>();

        for (int label = 0; label < labeledSets.size(); ++label) {
            Set<PairInt> set = labeledSets.get(label);
            for (PairInt p : set) {
                pointLabelMap.put(p, label);
            }
        }

        return pointLabelMap;
    }
    
    private TIntObjectMap<TIntSet> populateLabeLs(
        TIntObjectMap<CRegion> csRegions, TObjectIntMap<PairInt> pointLabelMap) {

        TIntObjectMap<TIntSet> map = new TIntObjectHashMap<TIntSet>();
        
        TIntObjectIterator<CRegion> iter = csRegions.iterator();
        for (int i = 0; i < csRegions.size(); ++i) {
            iter.advance();
            CRegion csr = iter.value();
            int rIdx = iter.key();
            
            for (Entry<PairInt, PairInt> entry : csr.offsetsToOrigCoords.entrySet()) {
                PairInt trXY = entry.getValue();
                int label = pointLabelMap.get(trXY);
                assert(pointLabelMap.containsKey(trXY));
                
                TIntSet set = map.get(rIdx);
                if (set == null) {
                    set = new TIntHashSet();
                    map.put(rIdx, set);
                }
                set.add(label);
            }
        }
        
        return map;
    }
    
    
    
    private void reduceToUnique(TIntObjectMap<CRegion> regions) {

        Map<OneDIntArray, Integer> labelsMap = new HashMap<OneDIntArray, Integer>();
        
        TIntObjectIterator<CRegion> iter = regions.iterator();
        for (int i = 0; i < regions.size(); ++i) {
            iter.advance();
            
            int rIdx = iter.key();
            CRegion r = iter.value();
        
            TIntList a = new TIntArrayList(r.labels);
            a.sort();
            OneDIntArray b = new OneDIntArray(a.toArray(new int[a.size()]));
            
            int n = r.offsetsToOrigCoords.size();
            
            Integer existingRIdx = labelsMap.get(b);
            if (existingRIdx == null) {
                labelsMap.put(b, Integer.valueOf(rIdx));
            } else {
                int count = regions.get(
                    existingRIdx.intValue()).offsetsToOrigCoords.size();
                if (count < n) {
                    labelsMap.put(b, Integer.valueOf(rIdx));
                }
            }
        }

        TIntObjectMap<CRegion> map = new TIntObjectHashMap<CRegion>();
        
        for (Entry<OneDIntArray, Integer> entry : labelsMap.entrySet()) {
            Integer rIndex = entry.getValue();
            int rIdx = rIndex.intValue();
            map.put(rIdx, regions.get(rIdx));
        }
        regions.clear();
        
        iter = map.iterator();
        for (int i = 0; i < map.size(); ++i) {
            iter.advance();
            
            int rIdx = iter.key();
            CRegion r = iter.value();
            
            regions.put(rIdx, r);
        }
    }
    
    // has the side effect of removing unlabeled points also
    private void populateLabels(TIntObjectMap<CRegion> cRegions, 
        TObjectIntMap<PairInt> pointLabelMap) {
        
        TIntObjectIterator<CRegion> iter = cRegions.iterator();
        for (int i = 0; i < cRegions.size(); ++i) {
            iter.advance();
            CRegion cr = iter.value();
            cr.labels.clear();
            Set<PairInt> rm = new HashSet<PairInt>();
            for (Entry<PairInt, PairInt> entry : cr.offsetsToOrigCoords.entrySet()) {
                PairInt xy = entry.getValue();
                if (!pointLabelMap.containsKey(xy)) {
                    rm.add(entry.getKey());
                    //System.out.println("missing for cr[" + iter.key() + "] " + xy);
                    continue;
                }        
                int label = pointLabelMap.get(xy);
                cr.labels.add(label);
            }
            for (PairInt p : rm) {
                cr.offsetsToOrigCoords.remove(p);
            }
        }
    }

}
