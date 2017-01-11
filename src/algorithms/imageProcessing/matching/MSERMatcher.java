package algorithms.imageProcessing.matching;

import algorithms.MultiArrayMergeSort;
import algorithms.QuickSort;
import algorithms.imageProcessing.FixedSizeSortedVector;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.features.CorrespondenceList;
import algorithms.imageProcessing.features.mser.Canonicalizer;
import algorithms.imageProcessing.features.mser.Canonicalizer.CRegion;
import algorithms.util.PairInt;
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
        TIntList h_0 = new TIntArrayList();
        TIntList s_0 = new TIntArrayList();
        TIntList v_0 = new TIntArrayList();
        TIntList h_1 = new TIntArrayList();
        TIntList s_1 = new TIntArrayList();
        TIntList v_1 = new TIntArrayList();
        
        float[] hsv = new float[3];
        
        for (int imgIdx0 = 0; imgIdx0 < nImg0; ++imgIdx0) {
            
            List<GreyscaleImage> rgb0 = pyr0.get(imgIdx0);
            
            TIntObjectMap<CRegion> cMap0 = cRegionsList00.get(imgIdx0);
            
            int w0 = rgb0.get(0).getWidth();
            int h0 = rgb0.get(0).getHeight();
            
            int np0 = cMap0.size();
            
            float scale0 = ((float)pyr0.get(0).get(0).getWidth()/(float)w0) 
                + ((float)pyr0.get(0).get(0).getHeight()/(float)h0);
            scale0 /= 2.f;
            img0Scales.put(imgIdx0, scale0);
            
            for (int imgIdx1 = 0; imgIdx1 < nImg1; ++imgIdx1) {
            
                List<GreyscaleImage> rgb1 = pyr1.get(imgIdx1);
                
                TIntObjectMap<CRegion> cMap1 = cRegionsList10.get(imgIdx1);
            
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

                        double ssdSum = 0;
                        int ssdCount = 0;

                        h_0.clear();
                        s_0.clear();
                        v_0.clear();
                        h_1.clear();
                        s_1.clear();
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
                            
                            h_0.add(Math.round(hsv[0]*255.f));
                            s_0.add(Math.round(hsv[1]*255.f));
                            v_0.add(Math.round(hsv[2]*255.f));
                            
                            Color.RGBtoHSB(rgb1.get(0).getValue(xy2), 
                                rgb1.get(1).getValue(xy2), 
                                rgb1.get(2).getValue(xy2), hsv);
                            
                            h_1.add(Math.round(hsv[0]*255.f));
                            s_1.add(Math.round(hsv[1]*255.f));
                            v_1.add(Math.round(hsv[2]*255.f));
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
                            int hDiff = h_0.get(jj) - h_1.get(jj);
                            int sDiff = s_0.get(jj) - s_1.get(jj);
                            int vDiff = (v_0.get(jj) - (int)v0Avg) - 
                                (v_1.get(jj) - (int)v1Avg);
                            //int vDiff = v_0.get(jj) - v_1.get(jj);
                            
                            ssdSum += (hDiff * hDiff + sDiff * sDiff + vDiff * vDiff);
                        }
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
       
        
        throw new UnsupportedOperationException("not yet implemented");
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
