package algorithms.imageProcessing.matching;

import algorithms.imageProcessing.FixedSizeSortedVector;
import algorithms.imageProcessing.GreyscaleImage;
import algorithms.imageProcessing.features.CorrespondenceList;
import algorithms.imageProcessing.features.mser.Canonicalizer;
import algorithms.imageProcessing.features.mser.Canonicalizer.CRegion;
import algorithms.util.PairInt;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.map.TIntObjectMap;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

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
    public CorrespondenceList matchObject(List<GreyscaleImage> pyr0,
        List<GreyscaleImage> pyr1, 
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
        
        List<FixedSizeSortedVector<Obj>> bestPerOctave 
            = new ArrayList<FixedSizeSortedVector<Obj>>();
        
        for (int imgIdx0 = 0; imgIdx0 < nImg0; ++imgIdx0) {
            GreyscaleImage mImg0 = pyr0.get(imgIdx0);
            TIntObjectMap<CRegion> cMap0 = cRegionsList00.get(imgIdx0);
                
            int np0 = cMap0.size();
            
            for (int imgIdx1 = 0; imgIdx1 < nImg1; ++imgIdx1) {
                GreyscaleImage mImg1 = pyr1.get(imgIdx1);
                TIntObjectMap<CRegion> cMap1 = cRegionsList10.get(imgIdx1);
                
                int np1 = cMap1.size();
                
                FixedSizeSortedVector<Obj> best01 = new 
                    FixedSizeSortedVector<Obj>(np0, Obj.class);
                bestPerOctave.add(best01);
               
                TIntObjectIterator<CRegion> iteri0 = cMap0.iterator();
                for (int i0 = 0; i0 < cMap0.size(); ++i0) {
                    iteri0.advance();

                    CRegion cr0 = iteri0.value();
                    int n0 = cr0.nTrEllipsePixels;

                    TIntObjectIterator<CRegion> iteri1 = cMap1.iterator();
                    for (int i1 = 0; i1 < cMap1.size(); ++i1) {
                        iteri1.advance();

                        CRegion cr1 = iteri1.value();
                        int n1 = cr1.nTrEllipsePixels;

                        int maxMatchable = Math.min(n0, n1);

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
                            int v1 = mImg1.getValue(xy2);

                            int v0 = mImg0.getValue(xy);

                            int diff = v0 - v1;

                            ssdSum += (diff * diff);
                            ssdCount++;
                        }
                        if (ssdCount == 0) {
                            continue;
                        }

                        ssdSum /= (double) ssdCount;
                        ssdSum = Math.sqrt(ssdSum);
                        ssdSum /= 255.;

                        double f = 1. - ((double) ssdCount / (double) maxMatchable);

                        double err = Math.max(cr0.autocorrel, cr1.autocorrel);

                        if (ssdSum <= 1.1 * err) {
                            
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
                                    "im1Idx=%d im2Idx=%d (%d,%d) (%d,%d) ",
                                    imgIdx0, imgIdx1,
                                    cr0.xC, cr0.yC, cr1.xC, cr1.yC);
                                String str2 = String.format(
                                    " ssd=%.2f autoc=%.2f,%.2f f=%.3f ssd.n=%d added=%b",
                                    (float) ssdSum,
                                    (float) cr0.autocorrel, (float) cr1.autocorrel,
                                    (float) f, ssdCount, added);
                                System.out.println(str1 + str2);
                            }
                        }
                    }
                }
            } // end loop over imgIdx1
        } // end loop over imgIdx0
        
        // process List<FixedSizeSortedVector<Obj>> bestPerOctave
    
        /*
        expecting that octave image pair with the 
            largest number of best costs for cregions is the matched
            image scales.
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
