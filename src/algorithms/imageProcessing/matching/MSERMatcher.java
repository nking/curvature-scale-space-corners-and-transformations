package algorithms.imageProcessing.matching;

import algorithms.compGeometry.PerimeterFinder2;
import algorithms.imageProcessing.DFSConnectedGroupsFinder;
import algorithms.imageProcessing.FixedSizeSortedVector;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.features.CorrespondenceList;
import algorithms.imageProcessing.features.HCPT;
import algorithms.imageProcessing.features.HOGs;
import algorithms.imageProcessing.features.mser.Canonicalizer;
import algorithms.imageProcessing.features.mser.Canonicalizer.CRegion;
import algorithms.imageProcessing.features.mser.Canonicalizer.RegionGeometry;
import algorithms.imageProcessing.features.mser.Canonicalizer.RegionPoints;
import algorithms.imageProcessing.features.mser.Region;
import algorithms.imageProcessing.transform.TransformationParameters;
import algorithms.imageProcessing.transform.Transformer;
import algorithms.imageProcessing.util.PairIntWithIndex;
import algorithms.misc.Misc;
import algorithms.misc.MiscDebug;
import algorithms.misc.MiscMath;
import algorithms.search.NearestNeighbor2D;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.QuadInt;
import com.climbwithyourfeet.clustering.DTClusterFinder;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.TIntSet;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

/**
 *
 * @author nichole
 */
public class MSERMatcher {

    private boolean debug = false;

    public void setToDebug() {
        debug = true;
    }

    public List<CorrespondenceList> matchObject3(
        List<List<GreyscaleImage>> pyrRGB0, List<GreyscaleImage> pyrPT0,
        TIntObjectMap<CRegion> cRegions0, List<Set<PairInt>> labeledSets0,
        TObjectIntMap<PairInt> pointLabelMap0,
        List<List<GreyscaleImage>> pyrRGB1, List<GreyscaleImage> pyrPT1,
        TIntObjectMap<CRegion> cRegions1, List<Set<PairInt>> labeledSets1,
        TObjectIntMap<PairInt> pointLabelMap1, String dbgLbl) {
                
        if (debug) {
            debugPrint2(cRegions0, pyrRGB0.get(0), "_csr_0_");
            debugPrint2(cRegions1, pyrRGB1.get(0), "_csr_1_");
            System.out.println("cr0.n=" + cRegions0.size() + " cr1.n=" +
                cRegions1.size());
        }
        
        // populated on demand, some are skipped for large size differences
        TIntObjectMap<TIntObjectMap<CRegion>> csr0
            = new TIntObjectHashMap<TIntObjectMap<CRegion>>();
        csr0.put(0, cRegions0);

        TIntObjectMap<TIntObjectMap<CRegion>> csr1
            = new TIntObjectHashMap<TIntObjectMap<CRegion>>();
        csr1.put(0, cRegions1);

        // key = region index, value = Obj w/ cost being hog intersection
        TIntObjectMap<FixedSizeSortedVector<Obj>> rIndexHOGMap
            = new TIntObjectHashMap<FixedSizeSortedVector<Obj>>();

        int n0 = pyrPT0.size();
        int n1 = pyrPT1.size();

        int w0 = pyrPT0.get(0).getWidth();
        int h0 = pyrPT0.get(0).getHeight();
        int w1 = pyrPT1.get(0).getWidth();
        int h1 = pyrPT1.get(0).getHeight();

        TIntObjectMap<HOGs> hogsMap1 = new TIntObjectHashMap<HOGs>();

        for (int imgIdx0 = 0; imgIdx0 < n0; ++imgIdx0) {

            GreyscaleImage gsI0 = combineImages(pyrRGB0.get(imgIdx0));
            GreyscaleImage ptI0 = pyrPT0.get(imgIdx0);

            int w0_i = ptI0.getWidth();
            int h0_i = ptI0.getHeight();
            float scale0 = (((float) w0 / (float) w0_i)
                + ((float) h0 / (float) h0_i)) / 2.f;

            HOGs hogs0 = new HOGs(gsI0, 1, 16);

            HCPT hcpt0 = new HCPT(ptI0, 1, 12, 12);
            
            TIntObjectMap<CRegion> regions0 = getOrCreate(csr0, imgIdx0, gsI0,
                scale0);

            for (int imgIdx1 = 0; imgIdx1 < n1; ++imgIdx1) {

                GreyscaleImage gsI1 = combineImages(pyrRGB1.get(imgIdx1));
                GreyscaleImage ptI1 = pyrPT1.get(imgIdx1);

                int w1_i = ptI1.getWidth();
                int h1_i = ptI1.getHeight();
                float scale1 = (((float) w1 / (float) w1_i)
                    + ((float) h1 / (float) h1_i)) / 2.f;

                TIntObjectMap<CRegion> regions1 = getOrCreate(csr1, imgIdx1,
                    gsI1, scale1);

                HOGs hogs1 = hogsMap1.get(imgIdx1);
                if (hogs1 == null) {
                    hogs1 = new HOGs(gsI1, 1, 16);
                    hogsMap1.put(imgIdx1, hogs1);
                }
                
                HCPT hcpt1 = new HCPT(ptI1, 1, 12, 12);

                TIntObjectIterator<CRegion> iter0 = regions0.iterator();
                for (int i0 = 0; i0 < regions0.size(); ++i0) {
                    iter0.advance();
                    int rIdx0 = iter0.key();
                    CRegion cr0 = iter0.value();

                    int sz0 = calculateObjectSize(cr0.offsetsToOrigCoords.values());

                    //int area0_full = csr0.get(0).get(rIdx0).offsetsToOrigCoords.size();
                    TIntObjectIterator<CRegion> iter1 = regions1.iterator();
                    for (int i1 = 0; i1 < regions1.size(); ++i1) {
                        iter1.advance();
                        int rIdx1 = iter1.key();
                        CRegion cr1 = iter1.value();

                        int sz1 = calculateObjectSize(cr1.offsetsToOrigCoords.values());

                        // size filter
                        if ((sz1 > sz0 && ((sz1 / sz0) > 1.2))
                            || (sz0 > sz1 && ((sz0 / sz1) > 1.2))) {
                            continue;
                        }

                        //double[]{intersectionSSD, f0, f1, count};
                        double[] hogCosts = sumHOGCost2(hogs0, cr0, scale0,
                            hogs1, cr1, scale1
                        );

                        if (hogCosts == null) {
                            continue;
                        }
                        float hogCost = 1.f - (float) hogCosts[0];

                        double fracOfWhole = hogCosts[2];
                        {
                            // temporarily undo the area correction which isn't
                            // quite right over extreme scales
                            double area1 = cr1.offsetsToOrigCoords.size();
                            fracOfWhole = 1. - (hogCosts[3] / area1);
                        }
                        if (fracOfWhole < 0.) {
                            fracOfWhole = 0.;
                        }
                        
                        //double[]{sum, f, err, count}
                        double[] costs2 = sumHCPTCost(hcpt0, cr0, scale0, 
                            hcpt1, cr1, scale1);
                        double hcptCost = 1.f - costs2[0];

                        double cost = (float) Math.sqrt(hogCost * hogCost
                            + fracOfWhole * fracOfWhole
                            + 2. * hcptCost * hcptCost
                        );

                        Obj obj = new Obj();
                        obj.cr0 = cr0;
                        obj.cr1 = cr1;
                        obj.r0Idx = rIdx0;
                        obj.r1Idx = rIdx1;
                        obj.imgIdx0 = imgIdx0;
                        obj.imgIdx1 = imgIdx1;
                        obj.ssd = 1.f - (float) hogCosts[0];
                        obj.nMatched = (int) hogCosts[3];
                        obj.cost = cost;
                        obj.costs = new double[]{hogCost, fracOfWhole, hcptCost};

                        FixedSizeSortedVector<Obj> objVec = rIndexHOGMap.get(rIdx1);
                        if (objVec == null) {
                            objVec = new FixedSizeSortedVector<Obj>(n0, Obj.class);
                            rIndexHOGMap.put(rIdx1, objVec);
                        }

                        boolean added = objVec.add(obj);
                    }
                }
            }
        }
      
        float critDens = 2.f/5.f;
        //filterBySpatialProximity(critDens, rIndexHOGMap, pyrRGB0, pyrRGB1);
        
        System.out.println("r1 map size = " + cRegions1.size()
            + " size filtered = " + rIndexHOGMap.size());

        // re-ordering the best for each rIdx1:
        FixedSizeSortedVector<Obj> tmp
            = new FixedSizeSortedVector<Obj>(
                //rIndexHOGMap.size(), 
                5,
                Obj.class);

        // --- print out phog based rankings -----
        // printing range of hog values for a region1
        TIntObjectIterator<FixedSizeSortedVector<Obj>> iter2
            = rIndexHOGMap.iterator();

        StringBuilder sb = new StringBuilder();

        for (int i3 = 0; i3 < rIndexHOGMap.size(); ++i3) {

            iter2.advance();

            int rIdx = iter2.key();
            FixedSizeSortedVector<Obj> vec = iter2.value();

            int n = vec.getNumberOfItems();
            if (n == 0) {
                continue;
            }

            Obj obj0 = vec.getArray()[0];

            tmp.add(obj0);

            if (debug) {
                int imgIdx0 = obj0.imgIdx0;
                int imgIdx1 = obj0.imgIdx1;

                GreyscaleImage gsI0 = pyrRGB0.get(imgIdx0).get(1);
                GreyscaleImage gsI1 = pyrRGB1.get(imgIdx1).get(1);

                float scale00, scale01;
                {
                    int w0_i = gsI0.getWidth();
                    int h0_i = gsI0.getHeight();
                    scale00 = (((float) w0 / (float) w0_i) + ((float) h0 / (float) h0_i)) / 2.f;

                    int w1_i = gsI1.getWidth();
                    int h1_i = gsI1.getHeight();
                    scale01 = (((float) w1 / (float) w1_i) + ((float) h1 / (float) h1_i)) / 2.f;
                }

                String lbl = "_" + obj0.imgIdx0 + "_" + obj0.imgIdx1 + "_"
                    + obj0.r0Idx + "_" + obj0.r1Idx;

                int or0 = (int) Math.round(
                    obj0.cr0.ellipseParams.orientation * 180. / Math.PI);

                int or1 = (int) Math.round(
                    obj0.cr1.ellipseParams.orientation * 180. / Math.PI);

                String str1 = String.format("angles=(%d,%d ; %d,%d)",
                    or0, or1, obj0.cr0.hogOrientation, obj0.cr1.hogOrientation);

                sb.append(String.format(
                    "1st %s %d (%d,%d) best: %.3f (%d,%d) %s [%.3f,%.3f,%.3f] %s\n",
                    dbgLbl, i3, Math.round(scale01 * obj0.cr1.ellipseParams.xC),
                    Math.round(scale01 * obj0.cr1.ellipseParams.yC),
                    (float) obj0.cost,
                    Math.round(scale00 * obj0.cr0.ellipseParams.xC),
                    Math.round(scale00 * obj0.cr0.ellipseParams.yC), lbl,
                    (float) obj0.costs[0], (float) obj0.costs[1], (float) obj0.costs[2], 
                    str1
                ));
            }
        }
        if (debug) {
            System.out.println(sb.toString());
        }

        StringBuilder sb2 = new StringBuilder();

        for (int i3 = 0; i3 < tmp.getNumberOfItems(); ++i3) {

            Obj obj0 = tmp.getArray()[i3];

            int imgIdx0 = obj0.imgIdx0;
            int imgIdx1 = obj0.imgIdx1;

            GreyscaleImage gsI0 = pyrRGB0.get(imgIdx0).get(1);
            GreyscaleImage gsI1 = pyrRGB1.get(imgIdx1).get(1);

            float scale00, scale01;
            {
                int w0_i = gsI0.getWidth();
                int h0_i = gsI0.getHeight();
                scale00 = (((float) w0 / (float) w0_i) + ((float) h0 / (float) h0_i)) / 2.f;

                int w1_i = gsI1.getWidth();
                int h1_i = gsI1.getHeight();
                scale01 = (((float) w1 / (float) w1_i) + ((float) h1 / (float) h1_i)) / 2.f;
            }

            if (debug) {
                String lbl = "_" + obj0.imgIdx0 + "_" + obj0.imgIdx1 + "_"
                    + obj0.r0Idx + "_" + obj0.r1Idx;

                int or0 = (int) Math.round(
                    obj0.cr0.ellipseParams.orientation * 180. / Math.PI);

                int or1 = (int) Math.round(
                    obj0.cr1.ellipseParams.orientation * 180. / Math.PI);

                String str1 = String.format("angles=(%d,%d ; %d,%d)",
                    or0, or1, obj0.cr0.hogOrientation, obj0.cr1.hogOrientation);

                sb2.append(String.format(
                    "2nd %s %d (%d,%d) best: %.3f (%d,%d) %s [%.3f,%.3f,%.3f] %s\n",
                    dbgLbl, i3, Math.round(scale01 * obj0.cr1.ellipseParams.xC),
                    Math.round(scale01 * obj0.cr1.ellipseParams.yC),
                    (float) obj0.cost,
                    Math.round(scale00 * obj0.cr0.ellipseParams.xC),
                    Math.round(scale00 * obj0.cr0.ellipseParams.yC), lbl,
                    (float) obj0.costs[0], (float) obj0.costs[1], (float) obj0.costs[1],
                    str1
                ));

                /*
                Image im0 = gsI0.copyToColorGreyscale();
                Image im1 = gsI1.copyToColorGreyscale();
                int[] clr = ImageIOHelper.getNextRGB(4);
                obj0.cr0.drawEachPixel(im0, 0, clr[0], clr[1], clr[2]);
                obj0.cr1.drawEachPixel(im1, 0, clr[0], clr[1], clr[2]);
                MiscDebug.writeImage(im0, dbgLbl + "_" + lbl);
                MiscDebug.writeImage(im1, dbgLbl + "_" + lbl);
                */
            }
        }
        if (debug) {
             System.out.println(sb2.toString());
        }

        List<CorrespondenceList> out = new ArrayList<CorrespondenceList>();
        
        for (int i = 0; i < tmp.getNumberOfItems(); ++i) {
            
            List<QuadInt> qs = new ArrayList<QuadInt>();
            
            Obj obj = tmp.getArray()[i];
            
            int imgIdx0 = obj.imgIdx0;
            int imgIdx1 = obj.imgIdx1;
            
            GreyscaleImage gsI0 = pyrRGB0.get(imgIdx0).get(1);
            GreyscaleImage gsI1 = pyrRGB1.get(imgIdx1).get(1);
            
            float scale0, scale1;
            {
                int w0_i = gsI0.getWidth();
                int h0_i = gsI0.getHeight();
                scale0 = (((float)w0/(float)w0_i) + ((float)h0/(float)h0_i))/2.f;
                
                int w1_i = gsI1.getWidth();
                int h1_i = gsI1.getHeight();
                scale1 = (((float)w1/(float)w1_i) + ((float)h1/(float)h1_i))/2.f;                
            }
            
            // NOTE: for now, just mser centers,
            // but should fill out more than this, including centroid of points
            
            int x0 = Math.round(scale0 * obj.cr0.ellipseParams.xC);
            int y0 = Math.round(scale0 * obj.cr0.ellipseParams.yC);
            int x1 = Math.round(scale1 * obj.cr1.ellipseParams.xC);
            int y1 = Math.round(scale1 * obj.cr1.ellipseParams.yC);
            QuadInt q = new QuadInt(x0, y0, x1, y1);
            qs.add(q);
            
            CorrespondenceList cor = new CorrespondenceList(qs);
            out.add(cor);
        }

        return out;
    }
    
    private TIntObjectMap<CRegion> getOrCreate(
        TIntObjectMap<TIntObjectMap<CRegion>> csrs,
        int imgIdx, GreyscaleImage rgb, float scale) {

        TIntObjectMap<CRegion> csrMap = csrs.get(imgIdx);
        if (csrMap != null) {
            return csrMap;
        }
        csrMap = new TIntObjectHashMap<CRegion>();
        csrs.put(imgIdx, csrMap);

        TIntObjectMap<CRegion> csrMap0 = csrs.get(0);

        int w = rgb.getWidth();
        int h = rgb.getHeight();

        TIntObjectIterator<CRegion> iter = csrMap0.iterator();
        for (int i = 0; i < csrMap0.size(); ++i) {
            iter.advance();
            int idx = iter.key();
            CRegion csr = iter.value();

            if (csr.offsetsToOrigCoords.size() < 9) {
                continue;
            }

            // these are in scale of individual octave (not full reference frame)
            Set<PairInt> scaledSet = extractScaledPts(csr, w, h, scale);

            if (scaledSet.size() < 9) {
                continue;
            }

            PairIntArray xy = new PairIntArray();

            Region r = new Region();
            for (PairInt pl : scaledSet) {
                r.accumulate(pl.getX(), pl.getY());
                xy.add(pl.getX(), pl.getY());
            }

            int[] xyCen = new int[2];
            r.calculateXYCentroid(xyCen, w, h);
            int x = xyCen[0];
            int y = xyCen[1];
            assert(x >= 0 && x < w);
            assert(y >= 0 && y < h);
            double[] m = r.calcParamTransCoeff();

            double angle = Math.atan(m[0]/m[2]);
            if (angle < 0) {
                angle += Math.PI;
            }

            double major = 2. * m[4];
            double minor = 2. * m[5];

            double ecc = Math.sqrt(major * major - minor * minor)/major;
            assert(!Double.isNaN(ecc));

            Map<PairInt, PairInt> offsetMap = Canonicalizer.createOffsetToOrigMap(
                x, y, xy, w, h, angle);

            double autocorrel = Canonicalizer.calcAutoCorrel(rgb, x, y, offsetMap);

            CRegion csRegion = new CRegion();
            csRegion.ellipseParams.orientation = angle;
            csRegion.ellipseParams.eccentricity = ecc;
            csRegion.ellipseParams.major = major;
            csRegion.ellipseParams.minor = minor;
            csRegion.ellipseParams.xC = x;
            csRegion.ellipseParams.yC = y;
            csRegion.offsetsToOrigCoords = offsetMap;
            csRegion.autocorrel = Math.sqrt(autocorrel)/255.;
            csRegion.labels.addAll(csr.labels);
            
            csrMap.put(idx, csRegion);
        }

        return csrMap;
    }

    private TIntObjectMap<RegionPoints> getOrCreate2(
        TIntObjectMap<TIntObjectMap<RegionPoints>> csrs,
        int imgIdx, GreyscaleImage rgb, float scale) {

        TIntObjectMap<RegionPoints> csrMap = csrs.get(imgIdx);
        if (csrMap != null) {
            return csrMap;
        }
        csrMap = new TIntObjectHashMap<RegionPoints>();
        csrs.put(imgIdx, csrMap);

        TIntObjectMap<RegionPoints> csrMap0 = csrs.get(0);

        int w = rgb.getWidth();
        int h = rgb.getHeight();

        TIntObjectIterator<RegionPoints> iter = csrMap0.iterator();
        for (int i = 0; i < csrMap0.size(); ++i) {
            iter.advance();
            int idx = iter.key();
            RegionPoints csr = iter.value();

            if (csr.points.size() < 9) {
                continue;
            }

            // these are in scale of individual octave (not full reference frame)
            Set<PairInt> scaledSet = extractScaledPts(csr, w, h, scale);

            if (scaledSet.size() < 4) {
                continue;
            }
                        
            Region r = new Region();
            for (PairInt pl : scaledSet) {
                r.accumulate(pl.getX(), pl.getY());
            }
                        
            RegionPoints rg = new RegionPoints();
            rg.ellipseParams = Canonicalizer.calculateEllipseParams(r,  w, h);
            rg.points = scaledSet;
            
            csrMap.put(idx, rg);
        }

        return csrMap;
    }

    private void match(TransformationParameters params,
        PairIntArray points0, PairIntArray points1,
        PairIntArray outputM0, PairIntArray outputM1,
        TIntObjectMap<PairInt> indexPoint0Map,
        int image1Width, int image1Height, int dMax) {

        Transformer transformer = new Transformer();

        PairIntArray points0Tr = transformer.applyTransformation(params,
            points0);

        NearestNeighbor2D nn1 = new NearestNeighbor2D(
            Misc.convert(points1), image1Width, image1Height);

        Set<PairInt> added1 = new HashSet<PairInt>();

        for (int i = 0; i < points0Tr.getN(); ++i) {
            int x0Tr = points0Tr.getX(i);
            int y0Tr = points0Tr.getY(i);

            Set<PairInt> nearest1 = nn1.findClosest(x0Tr, y0Tr, dMax);
            if (nearest1 == null) {
                continue;
            }

            //TODO: add iteration over nearest1 for closest in color

            PairInt p1 = null;
            for (PairInt p3 : nearest1) {
                if (added1.contains(p3)) {
                    continue;
                }
                p1 = p3;
                break;
            }
            if (p1 == null) {
                continue;
            }

            PairInt p0 = indexPoint0Map.get(i);
            assert(p0 != null);
            outputM0.add(p0);
            outputM1.add(p1);
            added1.add(p1);
        }
    }

    private GreyscaleImage combineImages(List<GreyscaleImage> rgb) {

        GreyscaleImage r = rgb.get(0);
        GreyscaleImage g = rgb.get(1);
        GreyscaleImage b = rgb.get(2);

        if (r.getWidth() != g.getWidth() || r.getWidth() != b.getWidth() ||
            r.getHeight() != g.getHeight() || r.getHeight() != b.getHeight()) {
            throw new IllegalArgumentException("r, g, and b must have same"
                + " width and height");
        }

        int w = r.getWidth();
        int h = r.getHeight();

        GreyscaleImage comb = new GreyscaleImage(w, h, r.getType());

        for (int i = 0; i < r.getNPixels(); ++i) {
            float v0 = r.getValue(i);
            float v1 = g.getValue(i);
            float v2 = b.getValue(i);
            float avg = (v0 + v1 + v2)/3.f;

            comb.setValue(i, Math.round(avg));
        }

        return comb;
    }

    private Set<PairInt> extractScaledPts(CRegion csr, int w, int h,
        float scale) {

        Set<PairInt> scaledSet = new HashSet<PairInt>();
        for (Entry<PairInt, PairInt> entry : csr.offsetsToOrigCoords.entrySet()) {

            PairInt p = entry.getValue();

            int xScaled = Math.round((float) p.getX() / scale);
            int yScaled = Math.round((float) p.getY() / scale);
            if (xScaled == -1) {
                xScaled = 0;
            }
            if (yScaled == -1) {
                yScaled = 0;
            }
            if (xScaled == w) {
                xScaled = w - 1;
            }
            if (yScaled == h) {
                yScaled = h - 1;
            }
            PairInt pOrigScaled = new PairInt(xScaled, yScaled);

            scaledSet.add(pOrigScaled);
        }

        return scaledSet;
    }

    private Set<PairInt> extractScaledPts(RegionPoints csr, int w, int h,
        float scale) {

        Set<PairInt> scaledSet = new HashSet<PairInt>();
        for (PairInt p : csr.points) {

            int xScaled = Math.round((float) p.getX() / scale);
            int yScaled = Math.round((float) p.getY() / scale);
            if (xScaled == -1) {
                xScaled = 0;
            }
            if (yScaled == -1) {
                yScaled = 0;
            }
            if (xScaled == w) {
                xScaled = w - 1;
            }
            if (yScaled == h) {
                yScaled = h - 1;
            }
            PairInt pOrigScaled = new PairInt(xScaled, yScaled);

            scaledSet.add(pOrigScaled);
        }

        return scaledSet;
    }
   
    private TIntIntMap calculateLabelSizes(List<Set<PairInt>> sets) {
    
        TIntIntMap map = new TIntIntHashMap();
        for (int i = 0; i < sets.size(); ++i) {
            int sz = MiscMath.calculateObjectSize(sets.get(i));
            map.put(i, sz);
        }
        return map;
    }

    
    private int calculateObjectSize(Collection<PairInt> values) {

        int[] minMaxXY = MiscMath.findMinMaxXY(values);
        int diffX = minMaxXY[1] - minMaxXY[0];
        int diffY = minMaxXY[3] - minMaxXY[2];
        double xy = Math.sqrt(diffX * diffX + diffY * diffY);
        
        return (int)Math.round(xy);
    }

    private void debugPrint3(List<List<Obj>> list, GreyscaleImage img0, 
        GreyscaleImage img1, float scale0, float scale1, String lbl) {

        for (int i = 0; i < list.size(); ++i) {
            
            List<Obj> objs = list.get(i);
            
            Image im0 = img0.copyToColorGreyscale();
            Image im1 = img1.copyToColorGreyscale();
            
            for (int j = 0; j < objs.size(); ++j) {
                int[] clr = ImageIOHelper.getNextRGB(j);
                Obj obj = objs.get(j);
                obj.cr0.drawEachPixel(im0, 0, clr[0], clr[1], clr[2]);
                obj.cr1.drawEachPixel(im1, 0, clr[0], clr[1], clr[2]);
            }
            
            MiscDebug.writeImage(im0, lbl + "__0__" + i);
            MiscDebug.writeImage(im1, lbl + "__1__" + i);
        }
    }

    //chordDiffSum, nMatches, p.n
    private double[] partialShapeCost(Obj obj, 
        float scale0, float scale1, 
        List<Set<PairInt>> labeledSets0, 
        List<Set<PairInt>> labeledSets1, 
        GreyscaleImage gsI0, GreyscaleImage gsI1) {
        
        // TODO: use caching if keep this method
        Set<PairInt> set0 = new HashSet<PairInt>();
        TIntIterator iter = obj.cr0.labels.iterator();
        while (iter.hasNext()) {
            int label = iter.next();
            Set<PairInt> s0 = labeledSets0.get(label);
            for (PairInt p : s0) {
                int x = Math.round((float)p.getX() / scale0);
                int y = Math.round((float)p.getY() / scale0);
                set0.add(new PairInt(x, y));
            }
        }
        
        set0 = reduceToContiguous(set0);
        
        Set<PairInt> set1 = new HashSet<PairInt>();
        iter = obj.cr1.labels.iterator();
        while (iter.hasNext()) {
            int label = iter.next();
            Set<PairInt> s0 = labeledSets1.get(label);
            for (PairInt p : s0) {
                int x = Math.round((float)p.getX() / scale1);
                int y = Math.round((float)p.getY() / scale1);
                set1.add(new PairInt(x, y));
            }
        }
        
        set1 = reduceToContiguous(set1);
        
        PerimeterFinder2 finder = new PerimeterFinder2();
        PairIntArray p = finder.extractOrderedBorder(set0);
        PairIntArray q = finder.extractOrderedBorder(set1);
    
        int dp = 1;
        if (p.getN() > 500 || q.getN() > 500) {
            int dn = Math.max(p.getN(), q.getN());
            dp += Math.ceil((float)dn/500.f);
        }
        
        PartialShapeMatcher2 matcher = new PartialShapeMatcher2();
        matcher.overrideSamplingDistance(dp);
        matcher.setToUseEuclidean();
        matcher.setToRemoveOutliers();
        PartialShapeMatcher2.Result result = matcher.match(p, q);

        double[] out = new double[] {
            result.chordDiffSum, result.getNumberOfMatches(), p.getN()
        };

        return out;        
    }

    private Set<PairInt> reduceToContiguous(Set<PairInt> set) {
        
        DFSConnectedGroupsFinder finder = new DFSConnectedGroupsFinder();
        finder.setMinimumNumberInCluster(1);
        finder.findConnectedPointGroups(set);
        if (finder.getNumberOfGroups() == 1) {
            return new HashSet<PairInt>(set);
        }
        
        int maxIdx = -1;
        int nMax = Integer.MIN_VALUE;
        for (int i = 0; i < finder.getNumberOfGroups(); ++i) {
            int n = finder.getNumberofGroupMembers(i);
            if (n > nMax) {
                nMax = n;
                maxIdx = i;
            }
        }
        return finder.getXY(maxIdx);
    }

    //double[]{intersectionSSD, f0, f1, count};
    private double[] sumHOGCost2(HOGs hogs0, CRegion cr0, float scale0, 
        HOGs hogs1, CRegion cr1, float scale1) {
        
        //NOTE: the icecream tests show calculating dominant orientation
        // is necessary.
        // that suggests centering and orientation are not precise enough
        //    for strict comparisons.  Note that recalculating the 
        //    center and orientation w/ labeled segmentation only 
        //    at an earlier stage did not improve results.
        
        int orientation0 = hogs0.calculateDominantOrientation(
            cr0.offsetsToOrigCoords.values());
        cr0.hogOrientation = orientation0;
        
        int orientation1 = hogs1.calculateDominantOrientation(
            cr1.offsetsToOrigCoords.values());
        cr1.hogOrientation = orientation1;
        
        Map<PairInt, PairInt> offsetMap1 = cr1.offsetsToOrigCoords;

        double sum = 0;
        int count = 0;
        
        int[] h0 = new int[hogs0.getNumberOfBins()];
        int[] h1 = new int[h0.length];
        
        // key = transformed offsets, value = coords in image ref frame,
        // so, can compare dataset0 and dataset1 points with same
        //  keys
        for (Entry<PairInt, PairInt> entry0 : cr0.offsetsToOrigCoords.entrySet()) {

            PairInt pOffset0 = entry0.getKey();

            PairInt xy1 = offsetMap1.get(pOffset0);

            if (xy1 == null) {
                continue;
            }

            PairInt xy0 = entry0.getValue();

            hogs0.extractFeature(xy0.getX(), xy0.getY(), h0);

            hogs1.extractFeature(xy1.getX(), xy1.getY(), h1);

            float intersection = hogs0.intersection(h0, orientation0, 
                h1, orientation1);
            
            sum += (intersection * intersection);

            count++;
        }
        if (count == 0) {
            return null;
        }

        sum /= (double)count;

        sum = Math.sqrt(sum);
        
        // scaling count/area up to scale=1 reference frames
        double area1 = cr1.offsetsToOrigCoords.size();
        area1 /= scale1;
        double f1 = 1. - ((double) count / area1);
        
        double area0 = cr0.offsetsToOrigCoords.size();
        area0 /= scale0;
        double f0 = 1. - ((double) count / area0);
        
        return new double[]{sum, f0, f1, count};
    }

    private void filterBySpatialProximity(float critDens,
        TIntObjectMap<FixedSizeSortedVector<Obj>> rIndexHOGMap, 
        List<List<GreyscaleImage>> pyrRGB0, List<List<GreyscaleImage>> pyrRGB1) {
                
        System.out.println("before spatial filter rIndexes=" + rIndexHOGMap.size());
         
        int w0 = pyrRGB0.get(0).get(0).getWidth();
        int h0 = pyrRGB0.get(0).get(0).getHeight();
        int w1 = pyrRGB1.get(0).get(0).getWidth();
        int h1 = pyrRGB1.get(0).get(0).getHeight();
        
        Set<PairIntWithIndex> points2
                = new HashSet<PairIntWithIndex>();
        
        TIntObjectIterator<FixedSizeSortedVector<Obj>> iter2
            = rIndexHOGMap.iterator();

        for (int i3 = 0; i3 < rIndexHOGMap.size(); ++i3) {

            iter2.advance();

            int rIdx = iter2.key();
            FixedSizeSortedVector<Obj> vec = iter2.value();

            int n = vec.getNumberOfItems();
            if (n == 0) {
                continue;
            }

            Obj obj0 = vec.getArray()[0];

            int imgIdx0 = obj0.imgIdx0;
            int imgIdx1 = obj0.imgIdx1;

            GreyscaleImage gsI0 = pyrRGB0.get(imgIdx0).get(1);
            GreyscaleImage gsI1 = pyrRGB1.get(imgIdx1).get(1);

            float scale0, scale1;
            {
                int w0_i = gsI0.getWidth();
                int h0_i = gsI0.getHeight();
                scale0 = (((float) w0 / (float) w0_i) + ((float) h0 / (float) h0_i)) / 2.f;

                int w1_i = gsI1.getWidth();
                int h1_i = gsI1.getHeight();
                scale1 = (((float) w1 / (float) w1_i) + ((float) h1 / (float) h1_i)) / 2.f;
            }
            
            int x = Math.round(scale1 * obj0.cr1.ellipseParams.xC);
            int y = Math.round(scale1 * obj0.cr1.ellipseParams.yC);
            
            PairIntWithIndex pii = new PairIntWithIndex(x, y, rIdx);
            points2.add(pii);
        }
        
        DTClusterFinder<PairIntWithIndex> cFinder
            = new DTClusterFinder<PairIntWithIndex>(points2, w1 + 1, h1 + 1);
        cFinder.setMinimumNumberInCluster(2);
        cFinder.setCriticalDensity(critDens);
        cFinder.setThreshholdFactor(1.0f);
        cFinder.findClusters();

        //NOTE: may need to revise how to choose best region to keep.
        for (int i = 0; i < cFinder.getNumberOfClusters(); ++i) {
            Set<PairIntWithIndex> set = cFinder.getCluster(i);
            double minCost = Double.MAX_VALUE;
            int minCostRIdx = -1;

            for (PairIntWithIndex pii : set) {
                int rIdx = pii.getPixIndex();
                double cost = rIndexHOGMap.get(rIdx).getArray()[0].cost;
                if (cost < minCost) {
                    minCost = cost;
                    minCostRIdx = rIdx;
                }
            }
            assert(minCostRIdx > -1);
            for (PairIntWithIndex pii : set) {
                int rIdx = pii.getPixIndex();
                if (rIdx == minCostRIdx) {
                    continue;
                }
                rIndexHOGMap.remove(rIdx);
            }
        }
        System.out.println("after spatial filter rIndexes=" + rIndexHOGMap.size());
    }

    private HOGs getOrCreate(TIntObjectMap<HOGs> hogsMap, GreyscaleImage gs, int idx) {
        
        HOGs hogs = hogsMap.get(idx);
        if (hogs != null) {
            return hogs;
        }
        hogs = new HOGs(gs, 1, 16);
        
        hogsMap.put(idx, hogs);
        
        return hogs;
    }

    private HCPT getOrCreate2(TIntObjectMap<HCPT> hcptMap, GreyscaleImage pt, 
        int idx) {
        
        HCPT hcpt = hcptMap.get(idx);
        if (hcpt != null) {
            return hcpt;
        }
        hcpt = new HCPT(pt, 1, 12, 12);
        
        hcptMap.put(idx, hcpt);
        
        return hcpt;
    }

    private void calculateDominantOrientations(
        TIntObjectMap<RegionPoints> regionPoints, HOGs hogs) {

        TIntObjectIterator<RegionPoints> iter = regionPoints.iterator();
        for (int i = 0; i < regionPoints.size(); ++i) {
            iter.advance();
            
            RegionPoints r = iter.value();
            
            TIntSet orientations = hogs.calculateDominantOrientations(r.points);
        
            r.hogOrientations.addAll(orientations);
        }
    }

    /* NOTE: recalculating this worsens the solutions.
    Using the MSER originally determined ellipse parameters has better
    results for the small number of tests here
    private void recalcOrientationAndTrans(HOGs hogs, 
        TIntObjectMap<CRegion> regions, GreyscaleImage gsImg, float scale) {
        
        MiscellaneousCurveHelper ch = new MiscellaneousCurveHelper();
        
        TIntObjectIterator<CRegion> iter = regions.iterator();
        for (int i = 0; i < regions.size(); ++i) {
            iter.advance();
            int rIdx = iter.key();
            CRegion cr = iter.value();
            
            int orientation = hogs.calculateDominantOrientation(
                cr.offsetsToOrigCoords.values());

            Collection<PairInt> xyp = cr.offsetsToOrigCoords.values();
            
            PairIntArray xy = Misc.convertWithoutOrder(xyp);
            
            PairInt xyCen = ch.calculateXYCentroids2(xyp);
            
            cr.ellipseParams.orientation = Math.PI * orientation/180.;
            cr.ellipseParams.xC = xyCen.getX();
            cr.ellipseParams.yC = xyCen.getY();
            
            Map<PairInt, PairInt> offsetToOrigMap = 
                Canonicalizer.createOffsetToOrigMap(
                xyCen.getX(), xyCen.getY(), 
                xy, gsImg.getWidth(), gsImg.getHeight(), orientation);
        
            cr.offsetsToOrigCoords = offsetToOrigMap;
        }
    }
    */  

    private class Obj implements Comparable<Obj>{
        CRegion cr0;
        CRegion cr1;
        int imgIdx0;
        int imgIdx1;
        double ssd;
        int nMatched;
        double cost = Double.MAX_VALUE;
        double[] costs;
        
        // might not be populatated:
        int r0Idx = -1;
        int r1Idx = -1;

        @Override
        public int compareTo(Obj other) {
            if (cost < other.cost) {
                return -1;
            } else if (cost > other.cost) {
                return 1;
            }
            return 0;
        }
    }

    public static int distance(PairInt p1, PairInt p2) {
        int diffX = p1.getX() - p2.getX();
        int diffY = p1.getY() - p2.getY();
        return (int) Math.sqrt(diffX * diffX + diffY * diffY);
    }

    public static int distance(float x0, float y0, float x1, float y1) {
        float diffX = x0 - x1;
        float diffY = y0 - y1;
        return (int) Math.sqrt(diffX * diffX + diffY * diffY);
    }

    private void debugPrint2(TIntObjectMap<CRegion> cRegions,
        List<GreyscaleImage> rgb, String label) {

        Image img1 = rgb.get(1).copyToColorGreyscale();

        Image img2 = rgb.get(1).copyToColorGreyscale();

        TIntObjectIterator<CRegion> iter = cRegions.iterator();

        int nExtraDot = 0;

        for (int ii = 0; ii < cRegions.size(); ++ii) {
            iter.advance();

            int idx = iter.key();

            CRegion cr = iter.value();

            int[] clr = ImageIOHelper.getNextRGB(ii);

            cr.draw(img1, nExtraDot, clr[0], clr[1], clr[2]);

            cr.drawEachPixel(img2, nExtraDot, clr[0], clr[1], clr[2]);
        }

        MiscDebug.writeImage(img1, label + "_" + "_csrs_");

        MiscDebug.writeImage(img2, label + "_" + "_csrs_pix_");

        System.out.println(cRegions.size() + " labeled regions for " + label);
    }

    public List<CorrespondenceList> matchObject0(
        List<List<GreyscaleImage>> pyrRGB0, List<GreyscaleImage> pyrPT0, 
        TIntObjectMap<Canonicalizer.RegionPoints> regionPoints0, 
        List<List<GreyscaleImage>> pyrRGB1, List<GreyscaleImage> pyrPT1, 
        TIntObjectMap<Canonicalizer.RegionPoints> regionPoints1, 
        String debugLabel) {
        
        TIntObjectMap<HOGs> hogsMap0 = new TIntObjectHashMap<HOGs>();
        TIntObjectMap<HOGs> hogsMap1 = new TIntObjectHashMap<HOGs>();
        TIntObjectMap<HCPT> hcptMap0 = new TIntObjectHashMap<HCPT>();
        TIntObjectMap<HCPT> hcptMap1 = new TIntObjectHashMap<HCPT>();

        // use hogs to calculate the dominant orientations
        calculateDominantOrientations(regionPoints0, 
            getOrCreate(hogsMap0, combineImages(pyrRGB0.get(0)), 0));
        
        calculateDominantOrientations(regionPoints1, 
            getOrCreate(hogsMap1, combineImages(pyrRGB1.get(0)), 0));
        
        Canonicalizer canonicalizer = new Canonicalizer();
        
        // create the CRegion objects which have the rotated points and 
        //    offsets in them
        TIntObjectMap<CRegion> cRegions0 = canonicalizer.canonicalizeRegions4(
            regionPoints0, pyrRGB0.get(0).get(1));
        
        TIntObjectMap<CRegion> cRegions1 = canonicalizer.canonicalizeRegions4(
            regionPoints1, pyrRGB1.get(0).get(1));
                
        // populated on demand, some are skipped for large size differences
        TIntObjectMap<TIntObjectMap<CRegion>> csr0
            = new TIntObjectHashMap<TIntObjectMap<CRegion>>();
        csr0.put(0, cRegions0);

        TIntObjectMap<TIntObjectMap<CRegion>> csr1
            = new TIntObjectHashMap<TIntObjectMap<CRegion>>();
        csr1.put(0, cRegions1);

        // key = region index, value = Obj w/ cost being hog intersection
        TIntObjectMap<FixedSizeSortedVector<Obj>> rIndexHOGMap0
            = new TIntObjectHashMap<FixedSizeSortedVector<Obj>>();

        // key = region index, value = Obj w/ cost being hog intersection
        TIntObjectMap<FixedSizeSortedVector<Obj>> rIndexHOGMap1
            = new TIntObjectHashMap<FixedSizeSortedVector<Obj>>();

        int n0 = pyrPT0.size();
        int n1 = pyrPT1.size();

        int w0 = pyrPT0.get(0).getWidth();
        int h0 = pyrPT0.get(0).getHeight();
        int w1 = pyrPT1.get(0).getWidth();
        int h1 = pyrPT1.get(0).getHeight();
        
        float sizeFactor = 1.2f;

        for (int imgIdx0 = 0; imgIdx0 < n0; ++imgIdx0) {

            GreyscaleImage gsI0 = combineImages(pyrRGB0.get(imgIdx0));
            GreyscaleImage ptI0 = pyrPT0.get(imgIdx0);

            int w0_i = ptI0.getWidth();
            int h0_i = ptI0.getHeight();
            float scale0 = (((float) w0 / (float) w0_i)
                + ((float) h0 / (float) h0_i)) / 2.f;

            HOGs hogs0 = getOrCreate(hogsMap0, gsI0, imgIdx0);

            HCPT hcpt0 = getOrCreate2(hcptMap0, ptI0, imgIdx0);
                
            TIntObjectMap<CRegion> regions0 = getOrCreate(csr0, imgIdx0, gsI0,
                scale0);

            for (int imgIdx1 = 0; imgIdx1 < n1; ++imgIdx1) {

                GreyscaleImage gsI1 = combineImages(pyrRGB1.get(imgIdx1));
                GreyscaleImage ptI1 = pyrPT1.get(imgIdx1);

                int w1_i = ptI1.getWidth();
                int h1_i = ptI1.getHeight();
                float scale1 = (((float) w1 / (float) w1_i)
                    + ((float) h1 / (float) h1_i)) / 2.f;

                TIntObjectMap<CRegion> regions1 = getOrCreate(csr1, imgIdx1,
                    gsI1, scale1);

                HOGs hogs1 = getOrCreate(hogsMap1, gsI1, imgIdx1);

                HCPT hcpt1 = getOrCreate2(hcptMap1, ptI1, imgIdx1);

                TIntObjectIterator<CRegion> iter0 = regions0.iterator();
                for (int i0 = 0; i0 < regions0.size(); ++i0) {
                    iter0.advance();
                    int rIdx0 = iter0.key();
                    CRegion cr0 = iter0.value();

                    int sz0 = calculateObjectSize(cr0.offsetsToOrigCoords.values());

                    //int area0_full = csr0.get(0).get(rIdx0).offsetsToOrigCoords.size();
                    TIntObjectIterator<CRegion> iter1 = regions1.iterator();
                    for (int i1 = 0; i1 < regions1.size(); ++i1) {
                        iter1.advance();
                        
                        // because these regions were made w/ hog orientations,
                        // there may be multiple regions with the same 
                        // original rIdx1 stored as cr.dataIdx, but having
                        // a different rIdx1 here.
                        // so cr.dataIdx is used below for the identity to keep
                        // the best match for the cr.dataIdx (== original rIdx)
                        int rIdx1 = iter1.key();
                        CRegion cr1 = iter1.value();

                        int sz1 = calculateObjectSize(cr1.offsetsToOrigCoords.values());

                        // size filter
                        if ((sz1 > sz0 && ((sz1 / sz0) > sizeFactor))
                            || (sz0 > sz1 && ((sz0 / sz1) > sizeFactor))) {
                            continue;
                        }

                        //double[]{intersectionSSD, f0, f1, count};
                        double[] hogCosts = sumHOGCost2(hogs0, cr0, scale0,
                            hogs1, cr1, scale1
                        );

                        if (hogCosts == null) {
                            continue;
                        }
                        float hogCost = 1.f - (float) hogCosts[0];

                        double fracOfWhole = hogCosts[2];
                        {
                            // temporarily undo the area correction which isn't
                            // quite right over extreme scales
                            double area1 = cr1.offsetsToOrigCoords.size();
                            fracOfWhole = 1. - (hogCosts[3] / area1);
                        }
                        if (fracOfWhole < 0.) {
                            fracOfWhole = 0.;
                        }
                        
                        //double[]{sum, f, err, count}
                        double[] costs2 = sumHCPTCost(hcpt0, cr0, scale0, 
                            hcpt1, cr1, scale1);
                        double hcptCost = 1.f - costs2[0];

                        double cost = (float) Math.sqrt(hogCost * hogCost
                            + fracOfWhole * fracOfWhole
                            + 2. * hcptCost * hcptCost
                        );

                        Obj obj = new Obj();
                        obj.cr0 = cr0;
                        obj.cr1 = cr1;
                        //obj.r1Idx = cr1.dataIdx;
                        obj.r0Idx = rIdx0;
                        obj.r1Idx = rIdx1;
                        obj.imgIdx0 = imgIdx0;
                        obj.imgIdx1 = imgIdx1;
                        obj.ssd = 1.f - (float) hogCosts[0];
                        obj.nMatched = (int) hogCosts[3];
                        obj.cost = cost;
                        obj.costs = new double[]{hogCost, fracOfWhole, hcptCost};

                        //NOTE: may need to consider the best match
                        //  for each rIdx, that is, consider multiple 
                        //  orientations for a region rather than keeping
                        //  the best orienation for a region only
                        
                        // add to r1 map (which is actually cr.dataIdx)
                        FixedSizeSortedVector<Obj> objVec = rIndexHOGMap1.get(
                            cr1.dataIdx);
                        if (objVec == null) {
                            objVec = new FixedSizeSortedVector<Obj>(n0, Obj.class);
                            rIndexHOGMap1.put(cr1.dataIdx, objVec);
                        }
                        boolean added = objVec.add(obj);
                        
                        // add to r0 map
                        objVec = rIndexHOGMap0.get(cr0.dataIdx);
                        if (objVec == null) {
                            // store top 5 for each r0Idx
                            objVec = new FixedSizeSortedVector<Obj>(5, Obj.class);
                            rIndexHOGMap0.put(cr0.dataIdx, objVec);
                        }
                        added = objVec.add(obj);
                    }
                }
            }
        }
      
        float critDens = 1.f/10.f;
        //filterBySpatialProximity(critDens, rIndexHOGMap, pyrRGB0, pyrRGB1);
        
        System.out.println("r1 points size = " + regionPoints1.size()
            + " r1 map size filtered = " + rIndexHOGMap1.size() 
            + " r0 map size filtered = " + rIndexHOGMap0.size());

        // re-ordering the best for each rIdx1:
        FixedSizeSortedVector<Obj> tmp1
            = new FixedSizeSortedVector<Obj>(
                //rIndexHOGMap.size(), 
                5,
                Obj.class);

        // --- print out phog based rankings -----
        // printing range of hog values for a region1
        TIntObjectIterator<FixedSizeSortedVector<Obj>> iter2
            = rIndexHOGMap1.iterator();

        StringBuilder sb = new StringBuilder();

        for (int i3 = 0; i3 < rIndexHOGMap1.size(); ++i3) {

            iter2.advance();

            int rIdx = iter2.key();
            FixedSizeSortedVector<Obj> vec = iter2.value();

            int n = vec.getNumberOfItems();
            if (n == 0) {
                continue;
            }

            Obj obj0 = vec.getArray()[0];

            tmp1.add(obj0);

            if (debug) {
                int imgIdx0 = obj0.imgIdx0;
                int imgIdx1 = obj0.imgIdx1;

                GreyscaleImage gsI0 = pyrRGB0.get(imgIdx0).get(1);
                GreyscaleImage gsI1 = pyrRGB1.get(imgIdx1).get(1);

                float scale00, scale01;
                {
                    int w0_i = gsI0.getWidth();
                    int h0_i = gsI0.getHeight();
                    scale00 = (((float) w0 / (float) w0_i) + ((float) h0 / (float) h0_i)) / 2.f;

                    int w1_i = gsI1.getWidth();
                    int h1_i = gsI1.getHeight();
                    scale01 = (((float) w1 / (float) w1_i) + ((float) h1 / (float) h1_i)) / 2.f;
                }

                String lbl = "_" + obj0.imgIdx0 + "_" + obj0.imgIdx1 + "_"
                    + obj0.r0Idx + "_" + obj0.r1Idx;

                int or0 = (int) Math.round(
                    obj0.cr0.ellipseParams.orientation * 180. / Math.PI);

                int or1 = (int) Math.round(
                    obj0.cr1.ellipseParams.orientation * 180. / Math.PI);

                String str1 = String.format("angles=(%d,%d ; %d,%d)",
                    or0, or1, obj0.cr0.hogOrientation, obj0.cr1.hogOrientation);

                sb.append(String.format(
                    "1] r1 %s %d (%d,%d) best: %.3f (%d,%d) %s [%.3f,%.3f,%.3f] %s\n",
                    debugLabel, i3, 
                    Math.round(scale01 * obj0.cr1.ellipseParams.xC),
                    Math.round(scale01 * obj0.cr1.ellipseParams.yC),
                    (float) obj0.cost,
                    Math.round(scale00 * obj0.cr0.ellipseParams.xC),
                    Math.round(scale00 * obj0.cr0.ellipseParams.yC), lbl,
                    (float) obj0.costs[0], (float) obj0.costs[1], (float) obj0.costs[2], 
                    str1
                ));
            }
        }
        if (debug) {
            System.out.println(sb.toString());
        }

        StringBuilder sb2 = new StringBuilder();

        for (int i3 = 0; i3 < tmp1.getNumberOfItems(); ++i3) {

            Obj obj0 = tmp1.getArray()[i3];

            int imgIdx0 = obj0.imgIdx0;
            int imgIdx1 = obj0.imgIdx1;

            GreyscaleImage gsI0 = pyrRGB0.get(imgIdx0).get(1);
            GreyscaleImage gsI1 = pyrRGB1.get(imgIdx1).get(1);

            float scale00, scale01;
            {
                int w0_i = gsI0.getWidth();
                int h0_i = gsI0.getHeight();
                scale00 = (((float) w0 / (float) w0_i) + ((float) h0 / (float) h0_i)) / 2.f;

                int w1_i = gsI1.getWidth();
                int h1_i = gsI1.getHeight();
                scale01 = (((float) w1 / (float) w1_i) + ((float) h1 / (float) h1_i)) / 2.f;
            }

            if (debug) {
                String lbl = "_" + obj0.imgIdx0 + "_" + obj0.imgIdx1 + "_"
                    + obj0.r0Idx + "_" + obj0.r1Idx;

                int or0 = (int) Math.round(
                    obj0.cr0.ellipseParams.orientation * 180. / Math.PI);

                int or1 = (int) Math.round(
                    obj0.cr1.ellipseParams.orientation * 180. / Math.PI);

                String str1 = String.format("angles=(%d,%d ; %d,%d)",
                    or0, or1, obj0.cr0.hogOrientation, obj0.cr1.hogOrientation);

                sb2.append(String.format(
                    "2] r1 %s %d (%d,%d) best: %.3f (%d,%d) %s [%.3f,%.3f,%.3f] %s\n",
                    debugLabel, i3, 
                    Math.round(scale01 * obj0.cr1.ellipseParams.xC),
                    Math.round(scale01 * obj0.cr1.ellipseParams.yC),
                    (float) obj0.cost,
                    Math.round(scale00 * obj0.cr0.ellipseParams.xC),
                    Math.round(scale00 * obj0.cr0.ellipseParams.yC), lbl,
                    (float) obj0.costs[0], (float) obj0.costs[1], (float) obj0.costs[1],
                    str1
                ));

                /*
                Image im0 = gsI0.copyToColorGreyscale();
                Image im1 = gsI1.copyToColorGreyscale();
                int[] clr = ImageIOHelper.getNextRGB(4);
                obj0.cr0.drawEachPixel(im0, 0, clr[0], clr[1], clr[2]);
                obj0.cr1.drawEachPixel(im1, 0, clr[0], clr[1], clr[2]);
                MiscDebug.writeImage(im0, dbgLbl + "_" + lbl);
                MiscDebug.writeImage(im1, dbgLbl + "_" + lbl);
                */
            }
        }
        if (debug) {
             System.out.println(sb2.toString());
        }
        
        // --- print the top 5 of region0 matches -----
        iter2 = rIndexHOGMap0.iterator();

        for (int i3 = 0; i3 < rIndexHOGMap0.size(); ++i3) {

            iter2.advance();

            int rIdx = iter2.key();
            FixedSizeSortedVector<Obj> vec = iter2.value();

            sb = new StringBuilder();
            
            int n = vec.getNumberOfItems();
            
            for (int j = 0; j < n; ++j) {
                
                Obj obj0 = vec.getArray()[j];

                if (debug) {
                    int imgIdx0 = obj0.imgIdx0;
                    int imgIdx1 = obj0.imgIdx1;

                    GreyscaleImage gsI0 = pyrRGB0.get(imgIdx0).get(1);
                    GreyscaleImage gsI1 = pyrRGB1.get(imgIdx1).get(1);

                    float scale00, scale01;
                    {
                        int w0_i = gsI0.getWidth();
                        int h0_i = gsI0.getHeight();
                        scale00 = (((float) w0 / (float) w0_i) + ((float) h0 / (float) h0_i)) / 2.f;

                        int w1_i = gsI1.getWidth();
                        int h1_i = gsI1.getHeight();
                        scale01 = (((float) w1 / (float) w1_i) + ((float) h1 / (float) h1_i)) / 2.f;
                    }

                    String lbl = "_" + obj0.imgIdx0 + "_" + obj0.imgIdx1 + "_"
                        + obj0.r0Idx + "_" + obj0.r1Idx;

                    int or0 = (int) Math.round(
                        obj0.cr0.ellipseParams.orientation * 180. / Math.PI);

                    int or1 = (int) Math.round(
                        obj0.cr1.ellipseParams.orientation * 180. / Math.PI);

                    String str1 = String.format("angles=(%d,%d ; %d,%d)",
                        or0, or1, obj0.cr0.hogOrientation, 
                        obj0.cr1.hogOrientation);

                    sb.append(String.format(
                        "1] r0 %s %d %d (%d,%d) best: %.3f (%d,%d) %s [%.3f,%.3f,%.3f] %s\n",
                        debugLabel, i3, j, 
                        Math.round(scale01 * obj0.cr1.ellipseParams.xC),
                        Math.round(scale01 * obj0.cr1.ellipseParams.yC),
                        (float) obj0.cost,
                        Math.round(scale00 * obj0.cr0.ellipseParams.xC),
                        Math.round(scale00 * obj0.cr0.ellipseParams.yC), lbl,
                        (float) obj0.costs[0], (float) obj0.costs[1], (float) obj0.costs[2], 
                        str1
                    ));
                }
            }
            if (debug) {
                System.out.println(sb.toString());
            }
        }

        // storing top 5 of r1 matches
        List<CorrespondenceList> out = new ArrayList<CorrespondenceList>();
        
        for (int i = 0; i < tmp1.getNumberOfItems(); ++i) {
            
            List<QuadInt> qs = new ArrayList<QuadInt>();
            
            Obj obj = tmp1.getArray()[i];
            
            int imgIdx0 = obj.imgIdx0;
            int imgIdx1 = obj.imgIdx1;
            
            GreyscaleImage gsI0 = pyrRGB0.get(imgIdx0).get(1);
            GreyscaleImage gsI1 = pyrRGB1.get(imgIdx1).get(1);
            
            float scale0, scale1;
            {
                int w0_i = gsI0.getWidth();
                int h0_i = gsI0.getHeight();
                scale0 = (((float)w0/(float)w0_i) + ((float)h0/(float)h0_i))/2.f;
                
                int w1_i = gsI1.getWidth();
                int h1_i = gsI1.getHeight();
                scale1 = (((float)w1/(float)w1_i) + ((float)h1/(float)h1_i))/2.f;                
            }
            
            // NOTE: for now, just mser centers,
            // but should fill out more than this, including centroid of points
            
            int x0 = Math.round(scale0 * obj.cr0.ellipseParams.xC);
            int y0 = Math.round(scale0 * obj.cr0.ellipseParams.yC);
            int x1 = Math.round(scale1 * obj.cr1.ellipseParams.xC);
            int y1 = Math.round(scale1 * obj.cr1.ellipseParams.yC);
            QuadInt q = new QuadInt(x0, y0, x1, y1);
            qs.add(q);
            
            CorrespondenceList cor = new CorrespondenceList(qs);
            out.add(cor);
        }

        return out;
    }
    
    private double[] sumHCPTCost(HCPT hcpt0, CRegion cr0, float scale0, 
        HCPT hcpt1, CRegion cr1, float scale1) {
        
        Map<PairInt, PairInt> offsetMap1 = cr1.offsetsToOrigCoords;

        double sum = 0;
        int count = 0;
        
        int[] h0 = new int[hcpt0.getNumberOfBins()];
        int[] h1 = new int[h0.length];
        
        // key = transformed offsets, value = coords in image ref frame,
        // so, can compare dataset0 and dataset1 points with same
        //  keys
        for (Entry<PairInt, PairInt> entry0 : cr0.offsetsToOrigCoords.entrySet()) {

            PairInt pOffset0 = entry0.getKey();

            PairInt xy1 = offsetMap1.get(pOffset0);

            if (xy1 == null) {
                continue;
            }

            PairInt xy0 = entry0.getValue();

            hcpt0.extractFeature(xy0.getX(), xy0.getY(), h0);

            hcpt1.extractFeature(xy1.getX(), xy1.getY(), h1);

            float intersection = hcpt0.intersection(h0, h1);
            
            sum += (intersection * intersection);

            count++;
        }
        if (count == 0) {
            return null;
        }

        sum /= (double)count;

        sum = Math.sqrt(sum);
        
        //NOTE: this may need revision.  now assuming that all invoker's 
        // have one object in cRegions0, hence, need to scale fraction
        // of whole so all are in same reference frame
        double area = cr0.offsetsToOrigCoords.size();
        area /= (scale1 * scale1);
        
        double f = 1. - ((double) count / area);

        // TODO: correct this if end up using it.
        //  it's based upon green only
        double err = Math.max(cr0.autocorrel, cr1.autocorrel);

        return new double[]{sum, f, err, count};            
    }
}
