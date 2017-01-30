package algorithms.imageProcessing.matching;

import algorithms.QuickSort;
import algorithms.compGeometry.PerimeterFinder2;
import algorithms.imageProcessing.DFSConnectedGroupsFinder;
import algorithms.imageProcessing.FixedSizeSortedVector;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageIOHelper;
import algorithms.imageProcessing.ImageProcessor;
import algorithms.imageProcessing.MiscellaneousCurveHelper;
import algorithms.imageProcessing.features.CorrespondenceList;
import algorithms.imageProcessing.features.HCPT;
import algorithms.imageProcessing.features.HOGs;
import algorithms.imageProcessing.features.mser.Canonicalizer;
import algorithms.imageProcessing.features.mser.Canonicalizer.CRegion;
import algorithms.imageProcessing.features.mser.Region;
import algorithms.imageProcessing.transform.MatchedPointsTransformationCalculator;
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
import algorithms.util.VeryLongBitString;
import com.climbwithyourfeet.clustering.DTClusterFinder;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.TFloatList;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TIntFloatMap;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.hash.TIntFloatHashMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.awt.Color;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
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
      
        // filter rIndexHogMap to keep best within a spatial proximity,
        //  to help ensure the top results are different regions
        float critDens = 2.f/5.f;
        //filterBySpatialProximity(critDens, rIndexHOGMap, pyrRGB0, pyrRGB1);
        
        // any  indexes still in map passed the size filter.
        // further looking at range of intersection values to see if
        //    can use intersection > 0.5 as a high pass filter.
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

            /*
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
            */
        }
        //System.out.println(sb.toString());

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
        System.out.println(sb2.toString());

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

    private Map<PairInt, CRegion> makePointCRMap(TIntObjectMap<CRegion> crMap) {

        Map<PairInt, CRegion> output = new HashMap<PairInt, CRegion>();

        TIntObjectIterator<CRegion> iter = crMap.iterator();
        for (int i = 0; i < crMap.size(); ++i) {
            iter.advance();
            CRegion cr = iter.value();
            output.put(new PairInt(cr.ellipseParams.xC, cr.ellipseParams.yC), cr);
        }

        return output;
    }

    private double sumCosts(List<GreyscaleImage> rgb0, List<GreyscaleImage> rgb1,
        Map<PairInt, CRegion> cr00Map, Map<PairInt, CRegion> cr01Map,
        Map<PairInt, CRegion> cr10Map, Map<PairInt, CRegion> cr11Map,
        PairIntArray m0, PairIntArray m1) {

        float[] hsv = new float[3];
        TIntList v_0 = new TIntArrayList();
        TIntList v_1 = new TIntArrayList();

        double ssdSumTot = 0;
        double fSumTot = 0;

        int n = m0.getN();
        for (int i = 0; i < n; ++i) {

            PairInt p0 = new PairInt(m0.getX(i), m0.getY(i));
            PairInt p1 = new PairInt(m1.getX(i), m1.getY(i));

            CRegion cr0 = cr01Map.get(p0);
            CRegion cr1 = cr11Map.get(p1);
            if (cr0 == null) {
                cr0 = cr00Map.get(p0);
            }
            if (cr1 == null) {
                cr1 = cr10Map.get(p1);
            }
            if (cr0 == null) {
                int z = 0;
            }
            if (cr1 == null) {
                int z = 0;
            }

            int n0 = cr0.offsetsToOrigCoords.size();
            int n1 = cr1.offsetsToOrigCoords.size();

            int maxMatchable = Math.min(n0, n1);

            v_0.clear();
            v_1.clear();

            long hsDiffSum = 0;
            double ssdSum = 0;
            int ssdCount = 0;

            for (Map.Entry<PairInt, PairInt> entry
                : cr0.offsetsToOrigCoords.entrySet()) {

                PairInt pOffsets = entry.getKey();
                PairInt xy = entry.getValue();

                PairInt xy2 = cr1.offsetsToOrigCoords.get(pOffsets);
                if (xy2 == null) {
                    continue;
                }

                Color.RGBtoHSB(rgb0.get(0).getValue(xy),
                    rgb0.get(1).getValue(xy),
                    rgb0.get(2).getValue(xy), hsv);

                int h_0 = Math.round(hsv[0]*255.f);
                int s_0 = Math.round(hsv[1]*255.f);
                v_0.add(Math.round(hsv[2]*255.f));

                Color.RGBtoHSB(rgb1.get(0).getValue(xy2),
                    rgb1.get(1).getValue(xy2),
                    rgb1.get(2).getValue(xy2), hsv);

                int h_1 = Math.round(hsv[0]*255.f);
                int s_1 = Math.round(hsv[1]*255.f);
                v_1.add(Math.round(hsv[2]*255.f));

                int hDiff = h_0 - h_1;
                int sDiff = s_0 - s_1;
                hsDiffSum += (hDiff * hDiff + sDiff * sDiff);
            }
            if (v_1.isEmpty()) {
                continue;
            }

            // average vs
            long v0Avg = 0;
            long v1Avg = 0;
            for (int jj = 0; jj < v_0.size(); ++jj) {
                v0Avg += v_0.get(jj);
                v1Avg += v_1.get(jj);
            }
            v0Avg /= v_0.size();
            v1Avg /= v_1.size();

            for (int jj = 0; jj < v_0.size(); ++jj) {
                int vDiff = (v_0.get(jj) - (int)v0Avg) -
                    (v_1.get(jj) - (int)v1Avg);
                //int vDiff = v_0.get(jj) - v_1.get(jj);

                ssdSum += (vDiff * vDiff);
            }
            ssdSum += hsDiffSum;
            ssdCount = v_0.size();

            ssdSum /= (double) ssdCount;
            ssdSum = Math.sqrt(ssdSum);
            //ssdSum /= (255. * 3.);
            ssdSum /= (255. * 2.);

            double f = 1. - ((double) ssdCount / (double) maxMatchable);

            // TODO: correct this if end up using it.
            //  it's based upon green only
            double err = Math.max(cr0.autocorrel, cr1.autocorrel);

            ssdSumTot += ssdSum;
            fSumTot += f;
        }

        ssdSumTot /= (double)n;

        fSumTot /= (double)n;

        double costTot = ssdSumTot * ssdSumTot + fSumTot * fSumTot;

        return costTot;
    }

    private PairIntArray reduceByScale(PairIntArray points, float scale, int w,
        int h) {

        Set<PairInt> added = new HashSet<PairInt>();

        PairIntArray out = new PairIntArray(points.getN());
        for (int jj = 0; jj < points.getN(); ++jj) {
            int x = Math.round((float) points.getX(jj) / scale);
            int y = Math.round((float) points.getY(jj) / scale);

            if (x < 0 || y < 0 || x >= w || y >= h) {
                continue;
            }

            PairInt p = new PairInt(x, y);
            if (added.contains(p)) {
                continue;
            }
            out.add(x, y);
            added.add(p);
        }

        return out;
    }

    private void printNUnique(TIntObjectMap<CRegion> cRegions, String label) {

        Set<PairInt> unique = new HashSet<PairInt>();

        TIntObjectIterator<CRegion> iter = cRegions.iterator();
        for (int i = 0; i < cRegions.size(); ++i) {
            iter.advance();
            CRegion cr = iter.value();
            PairInt p = new PairInt(cr.ellipseParams.xC, cr.ellipseParams.yC);
            unique.add(p);
        }

        System.out.println(label + " nUnique=" + unique.size());
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

    private void search(TIntObjectMap<CRegion> csr0,
        GreyscaleImage gs0, GreyscaleImage pt0, float scale0, int imgIdx0,
        TIntObjectMap<CRegion> csr1,
        GreyscaleImage gs1, GreyscaleImage pt1, float scale1, int imgIdx1,
        List<GreyscaleImage> rgb0, List<GreyscaleImage> rgb1, 
        FixedSizeSortedVector<Obj> b) {

        System.out.println("img0=" + imgIdx0 + " img1=" + imgIdx1 + " csr0.n=" +
            csr0.size() + " csr1.n=" + csr1.size());
        
        // coords are in the individual octave reference frames

        float factor = 2.5f;

        float maxArea = Float.MIN_VALUE;
        
        TIntObjectIterator<CRegion> iter0 = csr0.iterator();
        for (int i = 0; i < csr0.size(); ++i) {

            iter0.advance();
            
            int rAIdx = iter0.key();

            CRegion csrA = iter0.value();

            if (csrA.offsetsToOrigCoords.size() < 9) {
                continue;
            }

            double majA = csrA.ellipseParams.major;

            TIntObjectIterator<CRegion> iter1 = csr1.iterator();
            for (int j = 0; j < csr1.size(); ++j) {

                iter1.advance();

                int rBIdx = iter1.key();
                
                CRegion csrB = iter1.value();

                if (csrB.offsetsToOrigCoords.size() < 9) {
                    continue;
                }

                // this method is used in pyramid scaled comparisons,
                // so we expect that similar objects have similar size.
                // NOTE: this filter may need adjustment.
                double majB = csrB.ellipseParams.major;
                if ((majA > majB && (majA/majB) > factor) ||
                    (majB > majA && (majB/majA) > factor)) {

                    /*System.out.format(
                        "RMVD (%d,%d) (%d,%d) im0=%d im1=%d maj0=%.3f, maj1=%.2f\n\n",
                        Math.round(scale0 * csrA.ellipseParams.xC), 
                        Math.round(scale0 * csrA.ellipseParams.yC), 
                        Math.round(scale1 * csrB.ellipseParams.xC), 
                        Math.round(scale1 * csrB.ellipseParams.yC),
                        imgIdx0, imgIdx1, (float)csrA.ellipseParams.major, 
                        (float)csrB.ellipseParams.major);
                     */
                    continue;
                } /*else {
                    System.out.format(
                        "KEPT (%d,%d) (%d,%d) im0=%d im1=%d maj0=%.3f, maj1=%.2f\n\n",
                        csrA.xC, csrA.yC, csrB.xC, csrB.yC,
                        imgIdx0, imgIdx1, (float)csrA.major, (float)csrB.major);
                }*/

                // an assumption until the maxArea is calculated
                float nMaxMatchable = Math.min(csrA.offsetsToOrigCoords.size(),
                    csrB.offsetsToOrigCoords.size());

                //TODO: consider whether to return the matched pixels
                // [ssdSumRGB, ssdSumPTCIELUV, ssdHSV, f, err, ssdCount]
                double[] costs = matchRegion(
                    csrA, rgb0, gs0, pt0, scale0, 
                    csrB, rgb1, gs1, pt1, scale1);

                if (costs == null) {
                    continue;
                }
                                
                //TODO: may need to revise. trying to filter out false
                //   matches to keep the maxArea smaller
                if (costs[0] > 0.5 || costs[1] > 0.5) {
                    continue;
                }
                
                // NOTE, bA will be re-populated once maxArea is found amoun
                // all items that make it to this level.
                // because these sizes are roughly the same, will make an
                // assumption for now that can make the initial list using
                // fraction composed of min areas being compared meanwhile.

                Obj obj = new Obj();
                obj.cr0 = csrA;
                obj.cr1 = csrB;
                obj.r0Idx = rAIdx;
                obj.r1Idx = rBIdx;
                obj.imgIdx0 = imgIdx0;
                obj.imgIdx1 = imgIdx1;
                obj.ssd = costs[0] + costs[1] + costs[2];
                obj.nMatched = (int)costs[5];
                obj.cost = costs[0] + costs[1] + costs[2];// + costs[3];
                boolean addedA = b.add(obj);

                // [ssdSumRGB, ssdSumPTCIELUV, ssdHSV, f, err, ssdCount]
                
                if (addedA) {
                    if (costs[4] > maxArea) {
                        maxArea = (float)costs[5];
                    }
                }

                System.out.format(
                    "(%d,%d) (%d,%d) im0=%d im1=%d c=%.3f (%.3f,%.3f, %.3f) rs=(%d,%d) added=%b\n",
                    Math.round(scale0 * csrA.ellipseParams.xC), 
                    Math.round(scale0 * csrA.ellipseParams.yC), 
                    Math.round(scale1 * csrB.ellipseParams.xC), 
                    Math.round(scale1 * csrB.ellipseParams.yC),
                    imgIdx0, imgIdx1, (float)obj.cost,
                    (float)costs[0], (float)costs[1],  (float)costs[2],
                    rAIdx, rBIdx, addedA);
            }
        }

        // re-calculate the fraction of whole for bA and replace the
        //  entries with new ordered list
        /*
        FixedSizeSortedVector<Obj> bA2 
            = new FixedSizeSortedVector<Obj>(b.getFixedCapacity(), Obj.class);
        
        int n = b.getNumberOfItems(); 
        for (int i = 0; i < n; ++i) {
            Obj obj = b.getArray()[i];
            double ssd2 = obj.cost - obj.ssd;
            float f = 1.f - ((float)obj.nMatched/maxArea);
            obj.cost = ssd2 + f;
            bA2.add(obj);
        }
        
        for (int i = 0; i < n; ++i) {
            b.getArray()[i] = bA2.getArray()[i];
        }
        */
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

    /**
     *
     * @return [ssdSumRGB, ssdSumPTCIELUV, ssdSumHSV, f, err, ssdCount]
     */
    private double[] matchRegion(
        CRegion csr0, List<GreyscaleImage> rgb0, GreyscaleImage gs0, 
        GreyscaleImage pt0, float scale0,
        CRegion csr1, List<GreyscaleImage> rgb1, GreyscaleImage gs1, 
        GreyscaleImage pt1, float scale1) {

        Map<PairInt, PairInt> offsetMap1 = csr1.offsetsToOrigCoords;

        double area1 = offsetMap1.size();
        area1 /= scale1;
        
        double ssdSumGS = 0;
        double ssdSumPT = 0;
        double ssdSumHSV = 0;
        int ssdCount = 0;

        float[] hsv0 = new float[3];
        float[] hsv1 = new float[3];
        
        // key = transformed offsets, value = coords in image ref frame,
        // so, can compare dataset0 and dataset1 points with same
        //  keys
        for (Entry<PairInt, PairInt> entry0 : csr0.offsetsToOrigCoords.entrySet()) {

            PairInt pOffset0 = entry0.getKey();

            PairInt xy1 = offsetMap1.get(pOffset0);

            if (xy1 == null) {
                continue;
            }

            PairInt xy0 = entry0.getValue();

            int vRGB0 = gs0.getValue(xy0);

            int vPT0 = pt0.getValue(xy0);

            int vRGB1 = gs1.getValue(xy1);

            int vPT1 = pt1.getValue(xy1);

            int diffRGB = vRGB0 - vRGB1;

            // NOTE: because these are angles, calculating a correction for
            //   wrap around (note, scale of 360 i 0 to 255).
            //int diffPT = vPT0 - vPT1;
            if (vPT0 > vPT1) {
                // add a phase to next value if it's closer to current with phase
                if ((vPT0 - vPT1) > Math.abs(vPT0 - (vPT1 + 256))) {
                    vPT1 += 256;
                }
            } else if (vPT1 > vPT0) {
                // add a phase to next value if it's closer to current with phase
                if ((vPT1 - vPT0) > Math.abs(vPT1 - (vPT0 + 256.))) {
                    vPT0 += 256.;
                }
            }
            int diffPT = vPT0 - vPT1;

            Color.RGBtoHSB(rgb0.get(0).getValue(xy0), 
                rgb0.get(1).getValue(xy0), 
                rgb0.get(2).getValue(xy0), hsv0);
            Color.RGBtoHSB(rgb1.get(0).getValue(xy1), 
                rgb1.get(1).getValue(xy1), 
                rgb1.get(2).getValue(xy1), hsv1);

            double diffHSV = diffHSV(hsv0, hsv1);
            
            ssdSumGS += (diffRGB * diffRGB);

            ssdSumPT += (diffPT * diffPT);

            ssdSumHSV += (diffHSV * diffHSV);
            
            ssdCount++;
        }
        if (ssdCount == 0) {
            return null;
        }

        ssdSumGS /= (double)ssdCount;

        ssdSumPT /= (double)ssdCount;

        ssdSumHSV /= (double)ssdCount;
        
        ssdSumGS = Math.sqrt(ssdSumGS)/255.;

        ssdSumPT = Math.sqrt(ssdSumPT)/255.;

        ssdSumHSV = Math.sqrt(ssdSumHSV);        

        double f = 1. - ((double) ssdCount / area1);

        // TODO: correct this if end up using it.
        //  it's based upon green only
        double err = Math.max(csr0.autocorrel, csr1.autocorrel);

        return new double[]{ssdSumGS, ssdSumPT, ssdSumHSV, f, err, ssdCount};
    }

    /**
     * given the adjacency maps of the template, dataset 0, and the searchable
     * image, dataset 1, this looks for consistent adjacency and distance
     * between pairings in bestR and returns those.
     *
     * @param bestR
     * @param adjMap0
     * @param adjMap1
     * @param scale0
     * @param scale1
     * @param imgIdx0
     * @param imgIdx1
     * @return
     */
    private List<List<Obj>> filterForAdjacency(
        FixedSizeSortedVector<Obj> bestR,
        TObjectIntMap<PairInt> pointLabelMap0, TObjectIntMap<PairInt> pointLabelMap1,
        TIntObjectMap<TIntSet> regionAdjacencyMap0,
        TIntObjectMap<TIntSet> regionAdjacencyMap1,
        List<Set<PairInt>> labelSets0, List<Set<PairInt>> labelSets1,
        TIntObjectMap<VeryLongBitString> adjMap0, TIntObjectMap<VeryLongBitString> adjMap1,
        NearestNeighbor2D nn0, NearestNeighbor2D nn1,
        TIntObjectMap<PairInt> xyCenLabel0, TIntObjectMap<PairInt> xyCenLabel1,
        TIntIntMap labelSizes0, TIntIntMap labelSizes1,
        GreyscaleImage rgb0, GreyscaleImage rgb1,
        GreyscaleImage pt0, GreyscaleImage pt1,
        float scale0, float scale1, 
        int imgIdx0, int imgIdx1, int distTol) {
       
        if (bestR.getNumberOfItems() == 0) {
            return null;
        }
        
        if (adjMap0.size() == 1) {
            List<Obj> out = new ArrayList<Obj>();
            for (int i = 0; i < bestR.getArray().length; ++i) {
                out.add(bestR.getArray()[i]);
            }
            List<List<Obj>> out2 = new ArrayList<List<Obj>>();
            out2.add(out);
            return out2;
        }
        
        int imgWidth0 = rgb0.getWidth();
        int imgHeight0 = rgb0.getHeight();
        int imgWidth1 = rgb1.getWidth();
        int imgHeight1 = rgb1.getHeight();
                
        // these are in octave coord frames
        //key = unique xyCen0, value = index in xyCen0Objs list
        TObjectIntMap<PairInt> xyCen0Indexes = new TObjectIntHashMap<PairInt>();
        List<List<Obj>> xyCen0Objs = new ArrayList<List<Obj>>();
        List<PairInt> xyCen0s = new ArrayList<PairInt>();
        int n = bestR.getNumberOfItems();
        for (int i = 0; i < n; ++i) {
            Obj obj = bestR.getArray()[i];
            PairInt p = new PairInt(obj.cr0.ellipseParams.xC, obj.cr0.ellipseParams.yC);
            if (!xyCen0Indexes.containsKey(p)) {
                int sz = xyCen0Indexes.size();
                xyCen0s.add(p);
                xyCen0Objs.add(new ArrayList<Obj>());
                xyCen0Indexes.put(p, sz);
            }
            int idx = xyCen0Indexes.get(p);
            xyCen0Objs.get(idx).add(obj);
        }
        assert(xyCen0Indexes.size() == xyCen0Objs.size());
        assert(xyCen0Indexes.size() == xyCen0s.size());

        // key = r0Idx, r1Idx, value=index of b
        TObjectIntMap<PairInt> regionIndexesBMap = new TObjectIntHashMap<PairInt>();
        for (int i = 0; i < n; ++i) {
            Obj obj = bestR.getArray()[i];
            PairInt p = new PairInt(obj.r0Idx, obj.r1Idx);
            regionIndexesBMap.put(p, i);
        }
        
        List<List<Obj>> filtered = new ArrayList<List<Obj>>();

        int distTol2 = Math.round((float)distTol/scale1);
        if (distTol2 < 1) {
            distTol2 = 1;
        }

        System.out.println("xyCen0Objs.size()=" + xyCen0Objs.size()
            + " b.n=" + bestR.getNumberOfItems());

        /*{
            int n1 = labelsNotInRegions1.size();
            System.out.println(" also searching " + Arrays.toString(
                labelsNotInRegions1.toArray(new int[n1])));
        }*/
        
        int dMax = 2;

        float factor = 1.3f;
        
        Set<QuadInt> searchedRegions = new HashSet<QuadInt>();
        
        // List<List<Obj>> xyCen0Objs, index=unique first, item=list of
        //      first mapped to a dataset1 point

        for (int i = 0; i < xyCen0Objs.size(); ++i) {

            List<Obj> listI = xyCen0Objs.get(i);

            for (int ii = 0; ii < listI.size(); ++ii) {

                Obj objI = listI.get(ii);
                assert (objI.imgIdx0 == imgIdx0);
                assert (objI.imgIdx1 == imgIdx1);

                // ---- for the matching regions in objI, look for adjacent
                //      matching regions in the list,
                //      that have an objJ which creates a pair in
                //      reference frame 0 and reference frame 1 which have
                //      the same separation of centers in their frames
                //      and that the adjacent labeled regions have 
                //      similar sizes to one another

                TIntSet rIdxsI_0 = regionAdjacencyMap0.get(objI.r0Idx);
boolean dbg = false;
/*    
System.out.println(" objI (" + objI.cr0.ellipseParams.xC + "," +
objI.cr0.ellipseParams.yC + ") " + 
" (" + objI.cr1.ellipseParams.xC + "," + objI.cr1.ellipseParams.yC
);
//(13,15)  (15,55      // ridx 109 302
//(14,9) and (16,48)  // ridx 45 61
if (imgIdx0==2 && imgIdx1==0 &&
    (Math.abs(objI.cr0.ellipseParams.xC - 12) < 2) &&
    (Math.abs(objI.cr0.ellipseParams.yC - 15) < 2) &&
    (Math.abs(objI.cr1.ellipseParams.xC - 15) < 2) &&
    (Math.abs(objI.cr1.ellipseParams.yC - 55) < 2)
    ) {
    int bAIdx = regionIndexesBMap.get(new PairInt(109, 302));
    int bBIdx = regionIndexesBMap.get(new PairInt(45, 61));
    System.out.println("bA=" + bAIdx + " bBIdx=" + bBIdx);
    // error in the region adjacency maps
    dbg = true;
}
*/
                if (rIdxsI_0 == null) {
                    continue;
                }
                TIntSet rIdxsI_1 = regionAdjacencyMap1.get(objI.r1Idx);
                if (rIdxsI_1 == null) {
                    continue;
                }

                TIntSet bIdxs = new TIntHashSet();

                // try combinations within b that have a region in rIdxsI_0 
                //   and in rIdxsI_1
                TIntIterator iter4 = rIdxsI_0.iterator();
                while (iter4.hasNext()) {
                    int r0 = iter4.next();
                    TIntIterator iter5 = rIdxsI_1.iterator();
                    while (iter5.hasNext()) {
                        int r1 = iter5.next();
if (dbg) System.out.println("    r0=" + r0 + " r1=" + r1);       
                        PairInt p12 = new PairInt(r0, r1);
                        if (regionIndexesBMap.containsKey(p12)) {
                            bIdxs.add(regionIndexesBMap.get(p12));
                        }
                    }
                }
           
if (dbg) {
    int z = 0;
}               
//System.out.println("  (" + objI.cr0.ellipseParams.xC + "," +
//objI.cr0.ellipseParams.yC + ") " + 
//" (" + objI.cr1.ellipseParams.xC + "," + objI.cr1.ellipseParams.yC
//+ " nAdj=" + bIdxs.size());                
                    
                iter4 = bIdxs.iterator();
                while (iter4.hasNext()) {
                    int bIdx = iter4.next();
                    Obj objJ = bestR.getArray()[bIdx];                 
 
if (dbg) {
    System.out.println("       objJ (" + objJ.cr0.ellipseParams.xC + "," +
objJ.cr0.ellipseParams.yC + ") " + 
" (" + objJ.cr1.ellipseParams.xC + "," + objJ.cr1.ellipseParams.yC
);
}                    
                    QuadInt qIJ = new QuadInt(objI.r0Idx, objI.r1Idx,
                        objJ.r0Idx, objJ.r1Idx);
                    if (searchedRegions.contains(qIJ)) {
                        continue;
                    }
                    qIJ = new QuadInt(objJ.r0Idx, objJ.r1Idx, objI.r0Idx, objI.r1Idx);
                    if (searchedRegions.contains(qIJ)) {
                        continue;
                    }
                    searchedRegions.add(qIJ);
                                        
                    // compare separation and size
                    double d0 = distance(objI.cr0.ellipseParams.xC,
                        objI.cr0.ellipseParams.yC,
                        objJ.cr0.ellipseParams.xC, objJ.cr0.ellipseParams.yC);

                    double d1 = distance(objI.cr1.ellipseParams.xC,
                        objI.cr1.ellipseParams.yC,
                        objJ.cr1.ellipseParams.xC, objJ.cr1.ellipseParams.yC);
if (dbg) {
    System.out.println("     d0=" + d0 + " d1=" + d1);
}
                    //NOTE: caveat is that with occlusion, this dould
                    // remove a real match and that is important for
                    // case of matching a single object within images
                    // of all else being unmatchable, different scenes
                    if ((d0 > d1 && (d0 / d1) > factor)
                        || (d1 > d0 && (d1 / d0) > factor)) {
                        continue;
                    }

                    // compare region sizes
                    int szJ0 = calculateObjectSize(
                        objJ.cr0.offsetsToOrigCoords.values());

                    int szJ1 = calculateObjectSize(
                        objJ.cr1.offsetsToOrigCoords.values());
if (dbg) {
    System.out.println(
    "     szJ0=" + szJ0 + " szJ1=" + szJ1
    );
}
                    if ((szJ0 > szJ1 && (szJ0 / szJ1) > factor)
                        || (szJ1 > szJ0 && (szJ1 / szJ0) > factor)) {
                        continue;
                    }

                    List<Obj> out = new ArrayList<Obj>();
                    out.add(objI);
                    out.add(objJ);

                    filtered.add(out);
                }
            }

            /*
            if (dbg) {
                System.out.format(
                    "* RMVD (%d,%d)to(%d,%d) (%d,%d)to(%d,%d) d0=%.1f d1=%.1f\n",
                    Math.round(scale0 * objI.cr0.ellipseParams.xC),
                    Math.round(scale0 * objI.cr0.ellipseParams.yC),
                    xyCen0B.getX(), xyCen0B.getY(),
                    Math.round(scale1 * objI.cr1.ellipseParams.xC),
                    Math.round(scale1 * objI.cr1.ellipseParams.yC),
                    xyCen1B.getX(), xyCen1B.getY(),
                    (float)d0, (float)d1); 

                continue;
            }
            */                        
        }
        
        boolean restrictToAdj = true;

        System.out.println("filtered.size=" + filtered.size());

        // -- merge the consistent lists in filtered ---

        PairIntArray points0 = new PairIntArray();
        populateWithUniquePoints(filtered, points0);

        TObjectIntMap<PairInt> pIdxMap0 = new TObjectIntHashMap<PairInt>();
        for (int i = 0; i < points0.getN(); ++i) {
            PairInt p = new PairInt(points0.getX(i), points0.getY(i));
            pIdxMap0.put(p, i);
        }

        MatchedPointsTransformationCalculator tc =
            new MatchedPointsTransformationCalculator();

        Transformer transformer = new Transformer();

        TIntSet skip = new TIntHashSet();
        for (int i = 0; i < filtered.size(); ++i) {
            if (skip.contains(i)) {
                continue;
            }
            List<Obj> listI = filtered.get(i);
            Set<QuadInt> idxsI = new HashSet<QuadInt>();
            TIntSet adjI_0 = restrictToAdj ? new TIntHashSet() : null;
            TIntSet adjI_1 = restrictToAdj ? new TIntHashSet() : null;
            for (Obj obj : listI) {
                idxsI.add(new QuadInt(obj.cr0.ellipseParams.xC, 
                    obj.cr0.ellipseParams.yC,
                    obj.cr1.ellipseParams.xC, obj.cr1.ellipseParams.yC));
                if (restrictToAdj) {
                    TIntSet a = extractAdjacentLabels(obj.cr0, adjMap0);
                    if (a != null) {
                        adjI_0.addAll(a);
                    }
                    a = extractAdjacentLabels(obj.cr1, adjMap1);
                    if (a != null) {
                        adjI_1.addAll(a);
                    }
                }
            }

            TransformationParameters params = null;
            PairIntArray points0Tr = null;

            for (int j = (i + 1); j < filtered.size(); ++j) {
                if (skip.contains(j)) {
                    continue;
                }
                List<Obj> listJ = filtered.get(j);
        
                if (restrictToAdj) {
                    boolean foundAdjI = false;
                    boolean foundAdjJ = false;
                    for (Obj obj : listJ) {
                        for (int label0 : obj.cr0.labels.toArray(new int[obj.cr0.labels.size()])) {
                            if (adjI_0.contains(label0)) {foundAdjI = true;break;}
                        }
                        if (foundAdjI && foundAdjJ) {break;}
                        for (int label1 : obj.cr1.labels.toArray(new int[obj.cr1.labels.size()])) {
                            if (adjI_1.contains(label1)) {foundAdjJ = true; break;}
                        }
                        if (foundAdjI && foundAdjJ) {break;}
                    }
                    if (!foundAdjI || !foundAdjJ) {continue;}
                }
                
                int nIter = 0;
                for (Obj obj : listJ) {
                    QuadInt q = new QuadInt(obj.cr0.ellipseParams.xC, 
                        obj.cr0.ellipseParams.yC,
                        obj.cr1.ellipseParams.xC, obj.cr1.ellipseParams.yC);
                    if (idxsI.contains(q)) {
                        nIter++;
                    }
                }
                if (nIter == 0) {
                    continue;
                }
                if (params == null) {
                    Obj obj0 = listI.get(0);
                    Obj obj1 = listI.get(1);
                    params = tc.calulateEuclidean(
                        obj0.cr0.ellipseParams.xC, 
                        obj0.cr0.ellipseParams.yC, 
                        obj1.cr0.ellipseParams.xC, obj1.cr0.ellipseParams.yC,
                        obj0.cr1.ellipseParams.xC, obj0.cr1.ellipseParams.yC, 
                        obj1.cr1.ellipseParams.xC, obj1.cr1.ellipseParams.yC,
                        0, 0);
                    if (params == null) {
                        continue;
                    }
                    points0Tr = transformer.applyTransformation(params, points0);
                }
                /*
                if the listJ points transformed are
                   consistent with transformation,
                   merge them.
                */
                boolean isConsistent = true;
                for (Obj obj : listJ) {
                    PairInt p0 = new PairInt(obj.cr0.ellipseParams.xC, 
                        obj.cr0.ellipseParams.yC);
                    PairInt p1 = new PairInt(obj.cr1.ellipseParams.xC, 
                        obj.cr1.ellipseParams.yC);

                    int p0Idx = pIdxMap0.get(p0);
                    int xTr0 = points0Tr.getX(p0Idx);
                    int yTr0 = points0Tr.getY(p0Idx);

                    double dist = distance(xTr0, yTr0, p1.getX(), p1.getY());
                    if (dist > distTol2) {
                        isConsistent = false;
                        break;
                    }
                }
                if (isConsistent) {
                    // merge
                    for (Obj obj : listJ) {
                        QuadInt q = new QuadInt(obj.cr0.ellipseParams.xC, 
                            obj.cr0.ellipseParams.yC,
                            obj.cr1.ellipseParams.xC, obj.cr1.ellipseParams.yC);
                        if (!idxsI.contains(q)) {
                            listI.add(obj);
                        }
                    }
                    skip.add(j);
                }
            }
        }

        if (!skip.isEmpty()) {
            TIntList rm = new TIntArrayList(skip);
            rm.sort();
            for (int i = (rm.size() - 1); i > -1; --i) {
                int rmIdx = rm.get(i);
                filtered.remove(rmIdx);
            }
        }

        TFloatList sumCosts = new TFloatArrayList();

        for (int i = 0; i < filtered.size(); ++i) {

            List<Obj> list = filtered.get(i);

            float sumCost = 0;

            for (int j = 0; j < list.size(); ++j) {
                Obj obj = list.get(j);
                sumCost += (float)obj.cost;
                int x0 = Math.round(scale0 * obj.cr0.ellipseParams.xC);
                int y0 = Math.round(scale0 * obj.cr0.ellipseParams.yC);
            }
            sumCosts.add(sumCost);
        }

        QuickSort.sortBy1stArg(sumCosts, filtered);

        { // DEBUG
            for (int i = 0; i < filtered.size(); ++i) {
                List<Obj> list = filtered.get(i);
                StringBuilder sb = new StringBuilder("filter ").append(i);
                for (int j = 0; j < list.size(); ++j) {
                    Obj obj = list.get(j);
                    int x0 = Math.round(scale0 * obj.cr0.ellipseParams.xC);
                    int y0 = Math.round(scale0 * obj.cr0.ellipseParams.yC);
                    int x1 = Math.round(scale1 * obj.cr1.ellipseParams.xC);
                    int y1 = Math.round(scale1 * obj.cr1.ellipseParams.yC);
                    sb.append(String.format(" (%d,%d):(%d,%d) ", x0, y0, x1, y1));
                }
                for (int j = 0; j < list.size(); ++j) {
                    Obj obj = list.get(j);
                    sb.append(String.format(" c=%.3f ", obj.cost));
                }
                System.out.println(sb.toString());
            }
        }

        return filtered;
    }

    private NearestNeighbor2D createNN(List<Set<PairInt>> labeledSets) {

        Set<PairInt> allPoints = new HashSet<PairInt>();
        for (Set<PairInt> set : labeledSets) {
            allPoints.addAll(set);
        }

        int[] xyMinMax = MiscMath.findMinMaxXY(allPoints);

        NearestNeighbor2D nn = new NearestNeighbor2D(allPoints,
            xyMinMax[1] + 10, xyMinMax[3] + 10);

        return nn;
    }
    
    private List<Set<PairInt>> extractUnorderedBounds(List<Set<PairInt>> list) {
        
        PerimeterFinder2 finder = new PerimeterFinder2();
        
        List<Set<PairInt>> bounds = new ArrayList<Set<PairInt>>();
        for (Set<PairInt> set : list) {
            bounds.add(finder.extractOuterBorder(set));
        }
        
        return bounds;
    }    

    private void populateWithUniquePoints(List<List<Obj>> pairs,
        PairIntArray out0) {

        Set<PairInt> set0 = new HashSet<PairInt>();
        for (int i = 0; i < pairs.size(); ++i) {
            List<Obj> list = pairs.get(i);
            for (Obj obj : list) {
                PairInt p0 = new PairInt(obj.cr0.ellipseParams.xC, 
                    obj.cr0.ellipseParams.yC);
                PairInt p1 = new PairInt(obj.cr1.ellipseParams.xC, 
                    obj.cr1.ellipseParams.yC);
                if (!set0.contains(p0)) {
                    set0.add(p0);
                    out0.add(p0);
                }
            }
        }
    }

    private CRegion createScaledCRegion(Set<PairInt> labeledSet, 
        int imgWidth, int imgHeight, float scale) {
        
        Region r = new Region();
        
        for (PairInt p : labeledSet) {
            int x = Math.round((float)p.getX()/scale);
            int y = Math.round((float)p.getY()/scale);
            if (x == -1) {
                x = 0;
            }
            if (x == imgWidth) {
                x = imgWidth - 1;
            }
            if (y == -1) {
                y = 0;
            }
            if (y == imgHeight) {
                y = imgHeight - 1;
            }
            r.accumulate(x, y);
        }
        
        Canonicalizer c = new Canonicalizer();
        
        CRegion cr = c.canonicalizeRegion(r, imgWidth, imgHeight);
        if (cr == null) {
            return null;
        }
        
        CRegion csr = new CRegion();
        csr.ellipseParams.eccentricity = cr.ellipseParams.eccentricity;
        csr.ellipseParams.major = cr.ellipseParams.major;
        csr.ellipseParams.minor = cr.ellipseParams.minor;
        csr.offsetsToOrigCoords = cr.offsetsToOrigCoords;
        csr.ellipseParams.orientation = cr.ellipseParams.orientation;
        csr.ellipseParams.xC = cr.ellipseParams.xC;
        csr.ellipseParams.yC = cr.ellipseParams.yC;
        
        return csr;
    }

    private TIntObjectMap<PairInt> createCentroidMap(List<Set<PairInt>> sets) {

        MiscellaneousCurveHelper ch = new MiscellaneousCurveHelper();
        
        TIntObjectMap<PairInt> map = new TIntObjectHashMap<PairInt>();
        for (int i = 0; i < sets.size(); ++i) {
            int[] xy = ch.calculateRoundedXYCentroids(sets.get(i));
            PairInt p = new PairInt(xy[0], xy[1]);
            map.put(i, p);
        System.out.println("label " + i + "  xycen=" + p + 
            " n=" + sets.get(i).size());
        }
        
        return map;
    }

    private TIntIntMap calculateLabelSizes(List<Set<PairInt>> sets) {
    
        TIntIntMap map = new TIntIntHashMap();
        for (int i = 0; i < sets.size(); ++i) {
            int sz = MiscMath.calculateObjectSize(sets.get(i));
            map.put(i, sz);
        }
        return map;
    }

    private TIntSet extractAdjacentLabels(CRegion cr, 
        TIntObjectMap<VeryLongBitString> adjMap) {
        
        TIntSet labels = new TIntHashSet();
        TIntIterator iter = cr.labels.iterator();
        while (iter.hasNext()) {
            int label = iter.next();
            VeryLongBitString bs = adjMap.get(label);
            if (bs != null) {
                labels.addAll(bs.getSetBits());
            }
        }
        labels.removeAll(cr.labels);
    
        return labels;
    }

    private int calculateObjectSize(Collection<PairInt> values) {

        int[] minMaxXY = MiscMath.findMinMaxXY(values);
        int diffX = minMaxXY[1] - minMaxXY[0];
        int diffY = minMaxXY[3] - minMaxXY[2];
        double xy = Math.sqrt(diffX * diffX + diffY * diffY);
        
        return (int)Math.round(xy);
    }

    private Set<Obj> extractObjs(TIntSet labels, 
        TIntObjectMap<TIntSet> labelBestRIdxs, FixedSizeSortedVector<Obj> bestR) {

        Set<Obj> set = new HashSet<Obj>();
        
        TIntIterator iter = labels.iterator();
        while (iter.hasNext()) {
            int label = iter.next();
            TIntSet idxs = labelBestRIdxs.get(label);
            if (idxs == null) {
                continue;
            }
            TIntIterator iter2 = idxs.iterator();
            while (iter2.hasNext()) {
                int idx = iter2.next();
                Obj obj = bestR.getArray()[idx];
                set.add(obj);
            }
        }
        
        return set;
    }

    // key = region idx, value = adj regions w/o overlapping labels
    private TIntObjectMap<TIntSet> createRegionAdjacency(
        TIntObjectMap<CRegion> cRegions, TObjectIntMap<PairInt> pointLabelMap, 
        TIntObjectMap<VeryLongBitString> adjMap) {

        // key = region idx, value = adjacent regions w/o overlapping labels
        TIntObjectMap<TIntSet> out = new TIntObjectHashMap<TIntSet>();
        
        // key = label, value = region indexes w/ that label in it
        TIntObjectMap<TIntSet> labelRegions = new TIntObjectHashMap<TIntSet>();
        TIntObjectIterator<CRegion> iter = cRegions.iterator();
        for (int i = 0; i < cRegions.size(); ++i) {
            iter.advance();
            
            int rIdx = iter.key();
            CRegion cr = iter.value();
            TIntSet labels = cr.labels;
                        
            TIntIterator iter2 = labels.iterator();
            while (iter2.hasNext()) {
                int label = iter2.next();
                TIntSet regSet = labelRegions.get(label);
                if (regSet == null) {
                    regSet = new TIntHashSet();
                    labelRegions.put(label, regSet);
                }
                regSet.add(rIdx);
            }
        }
        
        // now have labelRegions w/ key = label, value = region indexes
        
        iter = cRegions.iterator();
        for (int i = 0; i < cRegions.size(); ++i) {
            iter.advance();
            
            int rIdx = iter.key();
            CRegion cr = iter.value();
            TIntSet labels = cr.labels;
                        
            TIntSet adjLabels = new TIntHashSet();
            TIntIterator iter2 = labels.iterator();
            while (iter2.hasNext()) {
                int label = iter2.next();
                VeryLongBitString bs = adjMap.get(label);
                if (bs != null) {
                    adjLabels.addAll(bs.getSetBits());
                }
            }
            adjLabels.removeAll(labels);
            
            // find the adjacent regions.
            // then keep the adjacent regions which do not have an 
            // intersection with cr
            iter2 = adjLabels.iterator();
            while (iter2.hasNext()) {
                int adjLabel = iter2.next();
                
                TIntSet regIdxs = labelRegions.get(adjLabel);
                if (regIdxs == null) {
                    continue;
                }
                TIntIterator iter3 = regIdxs.iterator();
                while (iter3.hasNext()) {
                    int r2Idx = iter3.next();
                    if (r2Idx == rIdx) {
                        continue;
                    }
                    
                    CRegion cr2 = cRegions.get(r2Idx);

                    int n2 = cr2.labels.size();
                    
                    TIntSet cr2Labels = new TIntHashSet(cr2.labels);
                    cr2Labels.removeAll(labels);
                    
                    //NOTE: this might need revision.
                    //  would prefer to only add an association if the intersection
                    //  is null, but can see that it's necessary to add one
                    //  when they are mostly not intersecting too.
                    //    might need to look at the fraction of intersection
                    //      over total
                    
                    if ((n2 - cr2Labels.size()) < 2) {
                        
                        //TODO: may need to improve this.  double checking that
                        // there are adjacenct points in cr and cr2
                        if (!hasALOAdjPoint(cr, cr2)) {
                            continue;
                        }
                        
                        TIntSet set = out.get(rIdx);
                        if (set == null) {
                            set = new TIntHashSet();
                            out.put(rIdx, set);
                        }
                        set.add(r2Idx);                        
                    }
                }
            } 
        }
        
        return out;
    }

    private TIntSet extractLabelsNotInRegions(TIntSet labelsSearched, 
        int nLabels) {
        
        TIntSet out = new TIntHashSet();
        for (int i = 0; i < nLabels; ++i) {
            if (!labelsSearched.contains(i)) {
                out.add(i);
            }
        }
        
        return out;
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

    private boolean hasALOAdjPoint(CRegion cr0, CRegion cr1) {

        CRegion t0, t1;
        if (cr0.offsetsToOrigCoords.size() < cr1.offsetsToOrigCoords.size()) {
            t0 = cr0;
            t1 = cr1;
        } else {
            t0 = cr1;
            t1 = cr0;
        }
        
        Set<PairInt> pts1 = new HashSet<PairInt>(t1.offsetsToOrigCoords.values());
        
        int[] dxs = Misc.dx4;
        int[] dys = Misc.dy4;
        for (Entry<PairInt, PairInt> entry : t0.offsetsToOrigCoords.entrySet()) {
            int x = entry.getValue().getX();
            int y = entry.getValue().getY();
            for (int k = 0; k < dxs.length; ++k) {
                PairInt p2 = new PairInt(x + dxs[k], y + dys[k]);
                if (pts1.contains(p2)) {
                    return true;
                }
            }
        }
        
        return false;
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

    private double diffHSV(float[] hsv0, float[] hsv1) {
        double sum = 0;
        for (int i = 0; i < 3; ++i) {
            sum += Math.abs(hsv0[i] - hsv1[i]);
        }
        return sum/3.;
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

    //double[]{sum, f, err, count}
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
        double area = cr1.offsetsToOrigCoords.size();
        area /= scale1;
        
        double f = 1. - ((double) count / area);

        // TODO: correct this if end up using it.
        //  it's based upon green only
        double err = Math.max(cr0.autocorrel, cr1.autocorrel);

        return new double[]{sum, f, err, count};            
    }

    // clculating deviation from a constant line of difference
    private double calcGreyTerm(float[] rgb0, float[] rgb1) {
        
        float dR = Math.abs(rgb0[0] - rgb1[0]);
        float dG = Math.abs(rgb0[1] - rgb1[1]);
        float dB = Math.abs(rgb0[2] - rgb1[2]);
        
        // deviation from a constant line
        
        float dAvg = (dR + dG + dB)/3.f;
        
        dR -= dAvg;
        dG -= dAvg;
        dB -= dAvg;

        double cost = Math.sqrt((dR*dR + dG*dG + dB*dB)/3.);
        
        // max value would be r,g,b diffs spread over entire range of 255
        // so that is 85.3...85.3/2=42.7,0,42.7 --> 34.8
        cost /= 34.82;
        
        return cost;
    }

    private float[] calcAvgRGB(CRegion cr, List<GreyscaleImage> rgb) {
        
        double r = 0;
        double g = 0;
        double b = 0;
        for (Entry<PairInt, PairInt> entry : cr.offsetsToOrigCoords.entrySet()) {
            PairInt xy = entry.getValue();
            r += rgb.get(0).getValue(xy);
            g += rgb.get(1).getValue(xy);
            b += rgb.get(2).getValue(xy);
        }
        
        r /= (double)cr.offsetsToOrigCoords.size();
        g /= (double)cr.offsetsToOrigCoords.size();
        b /= (double)cr.offsetsToOrigCoords.size();
        
        return new float[]{(float)r, (float)g, (float)b};
    }

    private boolean isGrey(float[] rgb, int deltaLimit) {

        float dRG = Math.abs(rgb[0] - rgb[1]);
        float dRB = Math.abs(rgb[0] - rgb[2]);
        float dBG = Math.abs(rgb[1] - rgb[2]);
        
        if (dRG < deltaLimit && dRB < deltaLimit && dBG < deltaLimit) {
            return true;
        }
        
        return false;
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

}
