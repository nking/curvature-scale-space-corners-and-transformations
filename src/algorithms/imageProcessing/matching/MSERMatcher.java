package algorithms.imageProcessing.matching;

import algorithms.MultiArrayMergeSort;
import algorithms.QuickSort;
import algorithms.imageProcessing.FixedSizeSortedVector;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.features.CorrespondenceList;
import algorithms.imageProcessing.features.mser.Canonicalizer;
import algorithms.imageProcessing.features.mser.Canonicalizer.CRegion;
import algorithms.imageProcessing.transform.MatchedPointsTransformationCalculator;
import algorithms.imageProcessing.transform.TransformationParameters;
import algorithms.imageProcessing.transform.Transformer;
import algorithms.misc.Misc;
import algorithms.search.NearestNeighbor2D;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import algorithms.util.QuadInt;
import algorithms.util.TrioInt;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.list.TFloatList;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TIntFloatMap;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntFloatHashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import java.awt.Color;
import java.util.ArrayList;
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
    
    /**
     * a single object is expected as input in reference frame 0 data,
     * else another method of matching using ORB and RANSAC would be faster
     * for images with nearly all regions in common.
     * 
     * @param pyr0 image pyramid for data0 ordered by descending size
     * @param pyr1 image pyramid for data1 ordered by descending size
     * @param cRegionsList00 list of pyramid scaled canonicalized regions for
     *    the images of pyr0.
     * @param cRegionsList01 list of pyramid scaled canonicalized regions for
     *    the inverted images images of pyr0.
     * @param cRegionsList10 list of pyramid scaled canonicalized regions for
     *    the images of pyr1.
     * @param cRegionsList11 list of pyramid scaled canonicalized regions for
     *    the inverted images images of pyr1.
     * @return 
     */
    public CorrespondenceList matchObject(List<List<GreyscaleImage>> pyr0,
        List<List<GreyscaleImage>> pyr1, 
        List<TIntObjectMap<CRegion>> cRegionsList00,
        List<TIntObjectMap<CRegion>> cRegionsList01,
        List<TIntObjectMap<CRegion>> cRegionsList10,
        List<TIntObjectMap<CRegion>> cRegionsList11) {
       
        int distTol = 10;
        
        int nImg0 = pyr0.size();
        int nImg1 = pyr1.size();
        
        if (nImg0 != cRegionsList00.size() || 
            nImg0 != cRegionsList01.size()) {
            throw new IllegalArgumentException("reference frame 0 lists"
                + " must be the same length");
        }
        
        if (nImg1 != cRegionsList10.size() || 
            nImg1 != cRegionsList11.size()) {
            throw new IllegalArgumentException("reference frame 1 lists"
                + " must be the same length");
        } 
        
        // key=(imgIdx1, imgIdx2), value=map w/ key=idx1, value=FixedSizeSortedVector
        Map<PairInt, TIntObjectMap<FixedSizeSortedVector<Obj>>> bestPerOctave
            = new HashMap<PairInt, TIntObjectMap<FixedSizeSortedVector<Obj>>>();
        
        TIntFloatMap img0Scales = new TIntFloatHashMap();
        TIntFloatMap img1Scales = new TIntFloatHashMap();
       
        //TODO: edit this to handle results for
        //    the cregions inverted list separate from the not-inverted
        //    and combine results.
        
        // collect hsv to normalize v afterwards
        TIntList v_0 = new TIntArrayList();
        TIntList v_1 = new TIntArrayList();
        
        float[] hsv = new float[3];
        
        for (int imgIdx0 = 0; imgIdx0 < nImg0; ++imgIdx0) {
            
            List<GreyscaleImage> rgb0 = pyr0.get(imgIdx0);
            
            TIntObjectMap<CRegion> cMap0 = cRegionsList01.get(imgIdx0);
            
            int w0 = rgb0.get(0).getWidth();
            int h0 = rgb0.get(0).getHeight();
            
            int np0 = cMap0.size();
            
            float scale0 = ((float)pyr0.get(0).get(0).getWidth()/(float)w0) 
                + ((float)pyr0.get(0).get(0).getHeight()/(float)h0);
            scale0 /= 2.f;
            img0Scales.put(imgIdx0, scale0);
            
            for (int imgIdx1 = 0; imgIdx1 < nImg1; ++imgIdx1) {
            
                List<GreyscaleImage> rgb1 = pyr1.get(imgIdx1);
                
                TIntObjectMap<CRegion> cMap1 = cRegionsList11.get(imgIdx1);
            
                if (imgIdx0 == imgIdx1 && imgIdx0 > 0) {
                    // same comparison as 0,0 but lower resolution
               //     continue;
                }
                
                int w1 = rgb1.get(0).getWidth();
                int h1 = rgb1.get(0).getHeight();
            
                int np1 = cMap1.size();
                
                float scale1 = ((float)pyr1.get(0).get(0).getWidth()/(float)w1) 
                    + ((float)pyr1.get(0).get(0).getHeight()/(float)h1);
                scale1 /= 2.f;
                img1Scales.put(imgIdx1, scale1);
                
                PairInt imgKey = new PairInt(imgIdx0, imgIdx1);
                bestPerOctave.put(imgKey, 
                    new TIntObjectHashMap<FixedSizeSortedVector<Obj>>());
                
                TIntObjectIterator<CRegion> iteri0 = cMap0.iterator();
                for (int i0 = 0; i0 < cMap0.size(); ++i0) {
                    iteri0.advance();
                    
                    int idx0 = iteri0.key();

                    FixedSizeSortedVector<Obj> best01 = new 
                        FixedSizeSortedVector<Obj>(np0, Obj.class);
                    bestPerOctave.get(imgKey).put(idx0, best01);
                    
                    CRegion cr0 = iteri0.value();
                    int n0 = cr0.nTrEllipsePixels;

                    TIntObjectIterator<CRegion> iteri1 = cMap1.iterator();
                    for (int i1 = 0; i1 < cMap1.size(); ++i1) {
                        iteri1.advance();

                        CRegion cr1 = iteri1.value();
                        
                        //TODO: revisit this.  size filter which
                        //   will be the wrong thing to use if occlusion
                        //   is present...
                        if (true) {
                            double t1, t2;
                            if (cr0.minor > cr1.minor) {
                                t1 = cr0.minor/cr1.minor;
                            } else {
                                t1 = cr1.minor/cr0.minor;
                            }
                            if (cr0.major > cr1.major) {
                                t2 = cr0.major/cr1.major;
                            } else {
                                t2 = cr1.major/cr0.major;
                            }
                            if (t1 > 1.5 || t2 > 1.5) {
                                continue;
                            }
                        }
                        
                        int n1 = cr1.nTrEllipsePixels;

                        int maxMatchable = Math.min(n0, n1);

                        double hsDiffSum = 0;
                        double ssdSum = 0;
                        int ssdCount = 0;

                        v_0.clear();
                        v_1.clear();
                        
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
                        hsDiffSum += ssdSum;
                        ssdCount = v_0.size();
                       
                        ssdSum /= (double) ssdCount;
                        ssdSum = Math.sqrt(ssdSum);
                        //ssdSum /= (255. * 3.);
                        ssdSum /= (255. * 2.);

                        double f = 1. - ((double) ssdCount / (double) maxMatchable);

                        // TODO: correct this if end up using it.
                        //  it's based upon green only
                        double err = Math.max(cr0.autocorrel, cr1.autocorrel);

                        if (ssdSum <= 3 * err) {
                            
                            Obj obj = new Obj();
                            obj.cr0 = cr0;
                            obj.cr1 = cr1;
                            obj.imgIdx0 = imgIdx0;
                            obj.imgIdx1 = imgIdx1;
                            obj.ssd = ssdSum;
                            obj.nMatched = ssdCount;
                            obj.cost = ssdSum + f;
                            
                            boolean added = best01.add(obj);
                            
                            if (debug) {
                                String str1 = String.format(
                                    "im1=%d im2=%d (%d,%d) (%d,%d) or1=%d or2=%d\n",
                                    imgIdx0, imgIdx1,
                                    cr0.xC, cr0.yC, cr1.xC, cr1.yC,
                                    (int)(cr0.orientation*180./Math.PI),
                                    (int)(cr1.orientation*180./Math.PI));
                                String str2 = String.format(
                                    "   ecc0=%.2f ecc1=%.2f min0=%.2f min1=%.2f maj0=%.2f maj1=%.2f\n",
                                    (float) cr0.eccentricity, (float) cr1.eccentricity,
                                    (float) cr0.minor, (float) cr1.minor,
                                    (float) cr0.major, (float) cr1.major);
                                String str3 = String.format(
                                    "   ssd=%.2f autoc=%.2f,%.2f f=%.3f ssd.n=%d added=%b",
                                    (float) ssdSum,
                                    (float) cr0.autocorrel, (float) cr1.autocorrel,
                                    (float) f, ssdCount, added);
                                System.out.println(str1 + str2 + str3);
                            }
                        }
                    }
                }
            } // end loop over imgIdx1
        } // end loop over imgIdx0
        
        /*
        // once thru to count members
        int count = 0;
        for (Entry<PairInt, TIntObjectMap<FixedSizeSortedVector<Obj>>> entry :
            bestPerOctave.entrySet()) {
                        
            TIntObjectMap<FixedSizeSortedVector<Obj>> map12 = 
                entry.getValue();
            
            TIntObjectIterator<FixedSizeSortedVector<Obj>> iter12 =
                map12.iterator();
            
            for (int ii = 0; ii < map12.size(); ++ii) {
                iter12.advance();
                FixedSizeSortedVector<Obj> vector = iter12.value();
                int n = vector.getNumberOfItems();
                count += n;
            }
        }
        
        // to more easily see matches, sorting by
        // imgidx0, imgidx1, cost
        int[] img0Idxs = new int[count];
        int[] img1Idxs = new int[count];
        float[] costs = new float[count];
        int[] idxs = new int[count];
        Obj[] objs = new Obj[count];
        
        count = 0;
        for (Entry<PairInt, TIntObjectMap<FixedSizeSortedVector<Obj>>> entry :
            bestPerOctave.entrySet()) {
            
            PairInt imgIndexes = entry.getKey();
            
            TIntObjectMap<FixedSizeSortedVector<Obj>> map12 = 
                entry.getValue();
            
            TIntObjectIterator<FixedSizeSortedVector<Obj>> iter12 =
                map12.iterator();
            
            Set<PairInt> uniquePt0s = new HashSet<PairInt>();
            
            for (int ii = 0; ii < map12.size(); ++ii) {
                iter12.advance();
                int idx = iter12.key();
                FixedSizeSortedVector<Obj> vector = iter12.value();
                int n = vector.getNumberOfItems();
                for (int j = 0; j < n; ++j) {
                    Obj obj = vector.getArray()[j];
                    img0Idxs[count] = obj.imgIdx0;
                    img1Idxs[count] = obj.imgIdx1;
                    costs[count] = (float)obj.cost;
                    idxs[count] = count;
                    objs[count] = obj;
                    count++;
                    
                    assert(imgIndexes.getX() == obj.imgIdx0);
                    assert(imgIndexes.getY() == obj.imgIdx1);
                
                    uniquePt0s.add(new PairInt(obj.cr0.xC, obj.cr0.yC));
                }       
            }
            System.out.println("for indexes " + imgIndexes.toString() + ""
                + " there are " + uniquePt0s.size() + " unique 0 coords");
        }
        assert(count == objs.length);
            
        MultiArrayMergeSort.sortBy1stThen2ndThen3rd(img0Idxs, img1Idxs, 
            costs, idxs);

        for (int j = 0; j < count; ++j) {
            int idx = idxs[j];
            Obj obj = objs[idx];
            int x0 = Math.round(
                (float)obj.cr0.xC * img0Scales.get(obj.imgIdx0));
            int y0 = Math.round(
                (float)obj.cr0.yC * img0Scales.get(obj.imgIdx0));
            int x1 = Math.round(
                (float)obj.cr1.xC * img1Scales.get(obj.imgIdx1));
            int y1 = Math.round(
                (float)obj.cr1.yC * img1Scales.get(obj.imgIdx1));

            String str1 = String.format(
                "pyr0=%d pyr1=%d (%d, %d) (%d, %d) c=%.3f\n", 
                obj.imgIdx0, obj.imgIdx1, x0, y0, x1, y1, (float)obj.cost);
        
            String str2 = String.format(
                "   ecc0=%.2f ecc1=%.2f min0=%.2f min1=%.2f maj0=%.2f maj1=%.2f",
                (float) obj.cr0.eccentricity, (float) obj.cr1.eccentricity,
                (float) obj.cr0.minor, (float) obj.cr1.minor,
                (float) obj.cr0.major, (float) obj.cr1.major);
            
            System.out.println(str1 + str2);
        }
        */
        /*
        NOTE: for one of the more difficult tests,
          that is, matching a template gbman from andriod statues 03
          to find within android statues 02 image the same template
          object, the true solution is possible, but
          marginally so.
          in the 16 pairs of octaves, there are 8 or less unique 
             regions from the template object,
             there are at least 2 true matches in 2 ovtave pairs,
             so euclidean pairwise template matches of all unique
             template coords within an octave pair could find the correct 
             solution, but only marginally so.
        
        will create an ObjectMatcher pre-processing method that uses
        segmentation and color filters to reduce some of the false 
        matches.
        NOTE: the keypoints in the full frame ORBs in that method
        could be brought here to improve evaluation.
        */
      
        PairIntArray points0 = new PairIntArray();
        PairIntArray points1 = new PairIntArray();
        populatePointLists(cRegionsList00, cRegionsList01,
            cRegionsList10, cRegionsList11, points0, points1);
        
        int w0 = pyr0.get(0).get(0).getWidth();
        int h0 = pyr0.get(0).get(0).getHeight();
        int w1 = pyr1.get(0).get(0).getWidth();
        int h1 = pyr1.get(0).get(0).getHeight();
        
        double bestCostOverall = Double.MAX_VALUE;
        PairIntArray bestM0Overall = null;
        PairIntArray bestM1Overall = null;
               
        for (Entry<PairInt, TIntObjectMap<FixedSizeSortedVector<Obj>>> entry :
            bestPerOctave.entrySet()) {
           
            int imgIdx0 = entry.getKey().getX();
            int imgIdx1 = entry.getKey().getY();
            
            TIntObjectMap<FixedSizeSortedVector<Obj>> map12 = 
                entry.getValue();
            
            TIntObjectIterator<FixedSizeSortedVector<Obj>> iter12 =
                map12.iterator();
        
            // for this octave pair, store the unique img0 keypoints
            // and then use pairwise euclidean transformations to
            // evaluate and find best solution
            
            // these coordinates are in their respective pyramid scales
            Map<PairInt, List<PairInt>> unique0 = new HashMap<PairInt,
                List<PairInt>>();
            
            for (int ii = 0; ii < map12.size(); ++ii) {
                iter12.advance();
                
                FixedSizeSortedVector<Obj> vector = iter12.value();
                
                int n = vector.getNumberOfItems();
                for (int jj = 0; jj < n; ++jj) {
                    Obj obj = vector.getArray()[jj];
                    PairInt p0 = new PairInt(obj.cr0.xC, obj.cr0.yC);
                    PairInt p1 = new PairInt(obj.cr1.xC, obj.cr1.yC);
           
                    assert(obj.imgIdx0 == imgIdx0);
                    assert(obj.imgIdx1 == imgIdx1);
                    
                    List<PairInt> list = unique0.get(p0);
                    if (list == null) {
                        list = new ArrayList<PairInt>();
                        unique0.put(p0, list);
                    }
                    list.add(p1);
                }
            }
            
            MatchedPointsTransformationCalculator tc = 
                new MatchedPointsTransformationCalculator();
         
            //List<Obj> bestPair = new ArrayList<Obj>();
            double bestCost = Double.MAX_VALUE;
            PairIntArray bestM0 = null;
            PairIntArray bestM1 = null;
            double bestC = Double.MAX_VALUE;
            int bestNMaxMatchable = points0.getN();
            
            List<GreyscaleImage> rgb0 = pyr0.get(imgIdx0);
            List<GreyscaleImage> rgb1 = pyr1.get(imgIdx1);
            
            float scale0 = ((float)w0/(float)rgb0.get(0).getWidth()) 
                + ((float)h0/(float)rgb0.get(0).getHeight());
            scale0 /= 2.f;
                
            float scale1 = ((float)w1/(float)rgb1.get(0).getWidth()) 
                + ((float)h1/(float)rgb1.get(0).getHeight());
            scale1 /= 2.f;
            
            int distTol2 = Math.round((float)distTol/scale1);
            if (distTol2 < 1) {
                distTol2 = 1;
            }
          
            Map<PairInt, CRegion> cr00Map = makePointCRMap(cRegionsList00.get(imgIdx0));
            Map<PairInt, CRegion> cr01Map = makePointCRMap(cRegionsList01.get(imgIdx0));
            Map<PairInt, CRegion> cr10Map = makePointCRMap(cRegionsList10.get(imgIdx1));
            Map<PairInt, CRegion> cr11Map = makePointCRMap(cRegionsList11.get(imgIdx1));
            
            PairIntArray points0Scale0 = new PairIntArray();
            PairIntArray points1Scale1 = new PairIntArray();
            points0Scale0.addAll(cr00Map.keySet());
            points0Scale0.addAll(cr01Map.keySet());
            points1Scale1.addAll(cr10Map.keySet());
            points1Scale1.addAll(cr11Map.keySet());
        
            TIntObjectMap<PairInt> indexPoint0Map = new TIntObjectHashMap<PairInt>();
            for (int i = 0; i < points0Scale0.getN(); ++i) {
                PairInt p = new PairInt(points0Scale0.getX(i), points0Scale0.getY(i));
                indexPoint0Map.put(i, p);
            }
        
            // combinations of unique0
            Set<QuadInt> visited = new HashSet<QuadInt>();
            
            for (Entry<PairInt, List<PairInt>> entry2 : unique0.entrySet()) {
                PairInt p2_0 = entry2.getKey();
                List<PairInt> p2_pairs_1 = entry2.getValue();
                
                for (Entry<PairInt, List<PairInt>> entry3 : unique0.entrySet()) {
                    PairInt p3_0 = entry3.getKey();
                    if (p2_0.equals(p3_0)) {
                        continue;
                    }
                    QuadInt q23 = new QuadInt(p2_0, p3_0);
                    QuadInt q32 = new QuadInt(p3_0, p2_0);
                    if (visited.contains(q23) || visited.contains(q32)) {
                        continue;
                    }
                    visited.add(q23);
                    
                    List<PairInt> p3_pairs_1 = entry3.getValue();
                    
                    for (PairInt p2_1 : p2_pairs_1) {
                        for (PairInt p3_1 : p3_pairs_1) {
                            if (p2_1.equals(p3_1)) {
                                continue;
                            }
                            TransformationParameters params 
                                = tc.calulateEuclidean(
                                p2_0.getX(), p2_0.getY(), 
                                p3_0.getX(), p3_0.getY(), 
                                p2_1.getX(), p2_1.getY(), 
                                p3_1.getX(), p3_1.getY(),
                                0, 0);
 
                            PairIntArray m0 = new PairIntArray(points0Scale0.getN());
                            PairIntArray m1 = new PairIntArray(points0Scale0.getN());
                                                        
                            match(params, 
                                points0Scale0, points1Scale1, 
                                m0, m1, indexPoint0Map, 
                                rgb1.get(0).getWidth(), rgb1.get(0).getHeight(), distTol2);
                         
                            // this has been corrected for area size nnd
                            // number of matches.
                            // it is ssd*ssd + f*f
                            double summedCost = sumCosts(rgb0, rgb1,
                                cr00Map, cr01Map, cr10Map, cr11Map, m0, m1);
                                                    
                            float f2 = 1.f - ((float)m0.getN()/
                                ((float)points0Scale0.getN()/scale0));
                            
                            double c = summedCost + (f2 * f2);
            
                            if (debug) {
                                System.out.format(
                                "___im0=%d im1=%d desc=%.3f c=%.3f n=%d nmaxm=%d\n",
                                imgIdx0, imgIdx1, (float) summedCost,
                                (float)c, m0.getN(), points0Scale0.getN());
                            }
                            
                            if (c < bestCost) {
                                bestM0 = m0;
                                bestM1 = m1;
                                bestCost = c;
                                bestC = summedCost;
                                bestNMaxMatchable = points0Scale0.getN();
                            }
                        }
                    }    
                }
            }
           
            if (bestM0 != null) {
                
                System.out.format(
                    "im0=%d im1=%d c=%.3f d=%.3f n=%d nmaxm=%d\n",
                    imgIdx0, imgIdx1, (float)bestCost,
                    (float)bestC, bestM0.getN(), bestNMaxMatchable);

                if (bestCost < bestCostOverall) {
                    bestCostOverall = bestCost;
                    bestM0Overall = bestM0;
                    bestM1Overall = bestM1;
                    
                    // rewrite m0 and m1 to full scale coords
                    for (int jj = 0; jj < bestM0.getN(); ++jj) {
                        int x0 = bestM0.getX(jj);
                        int y0 = bestM0.getY(jj);
                        x0 = Math.round((float)x0 * scale0);
                        y0 = Math.round((float)y0 * scale0);
                        bestM0.set(jj, x0, y0);
                        
                        int x1 = bestM1.getX(jj);
                        int y1 = bestM1.getY(jj);
                        x1 = Math.round((float)x1 * scale1);
                        y1 = Math.round((float)y1 * scale1);
                        bestM1.set(jj, x1, y1);
                    }
                }
            }
        }
        
        if (bestM0Overall == null) {
            return null;
        }
                    
        int n = bestM0Overall.getN();
        List<QuadInt> qs = new ArrayList<QuadInt>();
        for (int i = 0; i < n; ++i) {
            int x0 = bestM0Overall.getX(i);
            int y0 = bestM0Overall.getY(i);
            int x1 = bestM1Overall.getX(i);
            int y1 = bestM1Overall.getY(i);
            QuadInt q = new QuadInt(x0, y0, x1, y1);
            qs.add(q);
        }
        CorrespondenceList cor = new CorrespondenceList(qs);

        return cor;        
    }

    private void populatePointLists(
        List<TIntObjectMap<CRegion>> cRegionsList00, 
        List<TIntObjectMap<CRegion>> cRegionsList01, 
        List<TIntObjectMap<CRegion>> cRegionsList10, 
        List<TIntObjectMap<CRegion>> cRegionsList11,
        PairIntArray points0, PairIntArray points1) {
    
        populatePointList(cRegionsList00, points0);
        populatePointList(cRegionsList01, points0);
        populatePointList(cRegionsList10, points1);
        populatePointList(cRegionsList11, points1);
    }
    
    private void populatePointList(
        List<TIntObjectMap<CRegion>> cRegionsList, 
        PairIntArray points) {
        
        Set<PairInt> added = new HashSet<PairInt>();
        
        for (TIntObjectMap<CRegion> cRegionMap : cRegionsList) {
            TIntObjectIterator<CRegion> iter = cRegionMap.iterator();
            for (int i = 0; i < cRegionMap.size(); ++i) {
                iter.advance();
                CRegion cr = iter.value();
                PairInt p = new PairInt(cr.xC, cr.yC);
                if (added.contains(p)) {
                    continue;
                }
                added.add(p);
                points.add(p);
            }
        }
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
            output.put(new PairInt(cr.xC, cr.yC), cr);
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
            
            int n0 = cr0.nTrEllipsePixels;
            int n1 = cr1.nTrEllipsePixels;

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
    
    private class Obj implements Comparable<Obj>{
        CRegion cr0;
        CRegion cr1;
        int imgIdx0;
        int imgIdx1;
        double ssd;
        int nMatched;
        double cost;

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
    
}
