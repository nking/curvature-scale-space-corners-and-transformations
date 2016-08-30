package algorithms.imageProcessing.segmentation;

import algorithms.imageProcessing.ImageExt;
import algorithms.misc.MiscMath;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.logging.Logger;

/**
 * a class to subdivide the data, run
 * NormalizedCuts on the subdivisions,
 * then combine the results.
 * 
 * @author nichole
 */
public class NormalizedCutsWrapper {
    
    private int nPref = 200;
    
    private ColorSpace colorSpace = ColorSpace.RGB;
    private boolean ltRGB = false;
    
    private Logger log = Logger.getLogger(this.getClass().getName());
   
    public NormalizedCutsWrapper() {
    }
    
    public void setToLowThresholdRGB() {
        ltRGB = true;
        colorSpace = ColorSpace.RGB;
    }
   
    public void setColorSpaceToRGB() {
        ltRGB = false;
        colorSpace = ColorSpace.RGB;
    }
    
    public void setColorSpaceToHSV() {
        ltRGB = false;
        colorSpace = ColorSpace.HSV;
    }
    
    public void setColorSpaceToCIELAB() {
        ltRGB = false;
        colorSpace = ColorSpace.CIELAB;
    }
    
    public int[] normalizedCut(ImageExt img, int[] labels) {
        
        int maxLabel = MiscMath.findMax(labels);
        
        if (maxLabel <= nPref) {
            
            return run(img, labels);
        }
        
        float wDivH = (float)img.getWidth()/(float)img.getHeight();
        
        float r = (float)maxLabel/(float)nPref;

        // nX = wDivH * nY
        // nX + nY = (maxLabel/nPref)
        // nY = (maxLabel/nPref)/(wDivH + 1)
        // nX = wDivH * nY
        
        // need an overlap too expressed in pixels or percentage...
        // Ideally, would find the minimum dimension of any
        // contiguous labelled region and then the overlap
        // would need to be a little larger than that.
        
        int nBuffer = 20;
        
        float nYF = r/(wDivH + 1.f);
        int nY = (int)Math.ceil(nYF);
        int nX = (int)Math.ceil(wDivH * nYF);
        
        int binX = img.getWidth()/nX;
        int binY = img.getHeight()/nY;
        
        log.info("nX=" + nX + " nY=" + nY + " binX=" + binX
            + " binY=" + binY);
        
        if (binX < 2*nBuffer) {
            log.warning("WARNING: algorithm may need to be edited"
               + " for this size and number of labels");
            binX = 2*nBuffer;
            nX = img.getWidth()/binX;
            if (nX < 1) {
                nX = 1;
            }
        }
        if (binY < 2*nBuffer) {
            log.warning("WARNING: algorithm may need to be edited"
               + " for this size and number of labels");
            binY = 2*nBuffer;
            nY = img.getHeight()/binY;
            if (nY < 1) {
                nY = 1;
            }
        }        
        
        /* TODO: improve:
        maxLabel=874, nPref=200
        w=320 h=213
        nX=3
        nY=2
        i=0,j=0 has nLabels=4hundred something
            w/ stopX=125,stopY=125
        
        NOTE: the merging needs alot of improvement
        too...perhaps a very low resolution
        pass over the whole image to use the 
        results as a template
        to combine the individual pieces here.
        current results are over segmented and
        not necessarily an improvement over the input.
        */

        int nSub = nX * nY;
        
        int w = img.getWidth();
        int h = img.getHeight();

        int maxLabel2 = -1;
        
        int[] cLabels = new int[labels.length];
       
        for (int i = 0; i < nX; ++i) {
            int startX = i*binX - nBuffer;
            if (startX < 0) {
                startX = 0;
            }
            int stopX = startX + binX + nBuffer - 1;
            if (stopX >= w) {
                stopX = (w - 1);
            }
            for (int j = 0; j < nY; ++j) {
                int startY = j * binY - nBuffer;
                if (startY < 0) {
                    startY = 0;
                }
                int stopY = startY + binY + nBuffer - 1;
                if (stopY >= h) {
                    stopY = (h - 1);
                }
                
                ImageExt img2 = (ImageExt)img.copySubImage(
                    startX, stopX, startY, stopY);
                int nPix = img2.getNPixels();
                
                int[] sLabels = new int[nPix];
                for (int ii = startX; ii < stopX; ++ii) {
                    for (int jj = startY; jj < stopY; ++jj) {
                        int pixIdx = img.getInternalIndex(ii, jj);
                        int pixIdx2 = img2.getInternalIndex(
                            ii - startX, jj - startY);
                        sLabels[pixIdx2] = labels[pixIdx];
                    }
                }
                LabelToColorHelper.condenseLabels(sLabels);
                
if (MiscMath.findMax(sLabels) > 1.2*nPref) {
    int z = 1;
}

                int[] labels2 = run(img2, sLabels);
                               
                int x0 = startX;
                int x1 = stopX;
                int y0 = startY;
                int y1 = stopY;

                // -- merge any results with overlap to left
                // or overlap below in cLabels
                if (i == 0) {
                    maxLabel2 = reassignAboveMaxLabel(labels2, maxLabel2);
                } else {
                    // merge the overlapping buffer region to left
                    maxLabel2 = processLeft(cLabels, w, h,
                        labels2, img2.getWidth(), img2.getHeight(),
                        x0, y0, y1, 
                        nBuffer, maxLabel2);
                }
                if (j == 0) {
                    maxLabel2 = reassignAboveMaxLabel(labels2, maxLabel2);
                } else {
                    // merge the overlapping buffer region below
                    maxLabel2 = processBelow(cLabels, w, h,
                        labels2, img2.getWidth(), img2.getHeight(),
                        y0, x0, x1, 
                        nBuffer, maxLabel2);
                }
                
                log.info("maxLabel2=" + maxLabel2);
                
                // -- convert pixel indexes and place in cLabels
                for (int ii = x0; ii < x1; ++ii) {
                    for (int jj = y0; jj < y1; ++jj) {
                        int pixIdx = img.getInternalIndex(ii, jj);
                        int pixIdx2 = img2.getInternalIndex(
                            ii - startX, jj - startY);
                        cLabels[pixIdx] = sLabels[pixIdx2];
                    }
                }                                    
            }
        }
        
        return cLabels;
    }

    private int[] run(ImageExt img, int[] labels) {
        NormalizedCuts normCuts = new NormalizedCuts();

        if (ltRGB) {
            normCuts.setToLowThresholdRGB();
        } else if (ColorSpace.RGB.equals(colorSpace)) {
            normCuts.setColorSpaceToRGB();
        } else if (ColorSpace.HSV.equals(colorSpace)) {
            normCuts.setColorSpaceToHSV();
        } else if (ColorSpace.CIELAB.equals(colorSpace)) {
            normCuts.setColorSpaceToCIELAB();
        }

        return normCuts.normalizedCut(img, labels);
    }

    private int processLeft(int[] labels1, int w1, int h1, 
        int[] labels2, int w2, int h2,
        final int x0, final int y0, final int y1, 
        int nBuffer, int maxLabel) {
        
        return process(labels1, w1, h1, labels2, w2, h2, 
            true, x0, y0, y1, nBuffer, maxLabel);
    }
    
    private int processBelow(int[] labels1, int w1, int h1, 
        int[] labels2, int w2, int h2,
        final int y0, final int x0, final int x1, 
        int nBuffer, int maxLabel) {
        
        return process(labels1, w1, h1, labels2, w2, h2, 
            false, y0, x0, x1, nBuffer, maxLabel);
    }
    
    private int process(int[] labels1, int w1, int h1, 
        int[] labels2, int w2, int h2, boolean mergeLeft,
        final int a0, final int b0, final int b1, 
        int nBuffer, int maxLabel) {
        
        // key = label from labels2, value=existing label in
        //       labels1
        TIntObjectMap<TIntSet> overlap2 
            = new TIntObjectHashMap<TIntSet>();
        
        // key = label from labels1, value = set of pixels 
        //       in labels2 at that overlapping location
        TIntObjectMap<TIntSet> map1 = 
            new TIntObjectHashMap<TIntSet>();
       
        // key = label from labels2, value = set of pixels 
        //       in labels2 with that label
        TIntObjectMap<TIntSet> map2 = 
            new TIntObjectHashMap<TIntSet>();
        
        // --- handle the overlapping region ----
        if (mergeLeft) {
            /* merge left
                | binX   |    i=1
                |      [ | ]       
                      x0 x0+ x0+
                         nb  nb+nb
            */
            int x0 = a0;
            int y0 = b0;
            int y1 = b1;
            for (int i = x0; i < x0 + 2*nBuffer; ++i) {
                int i2 = i - x0 ;
                for (int j = y0; j < y1; ++j) {
                    int j2 = j - y0;
                    int pixIdx1 = (j * w1) + i;
                    int pixIdx2 = (j2 * w2) + i2;
                    int l1 = labels1[pixIdx1];
                    int l2 = labels2[pixIdx2];

                    TIntSet set0 = overlap2.get(l2);
                    if (set0 == null) {
                        set0 = new TIntHashSet();
                        overlap2.put(l2, set0);
                    }
                    set0.add(l1);

                    TIntSet set2 = map2.get(l2);
                    if (set2 == null) {
                        set2 = new TIntHashSet();
                        map2.put(l2, set2);
                    }
                    set2.add(pixIdx2);

                    TIntSet set1 = map1.get(l1);
                    if (set1 == null) {
                        set1 = new TIntHashSet();
                        map1.put(l1, set1);
                    }
                    set1.add(pixIdx2);
                }
            }
        } else {
            /*  merge below
                .......
                ------- y0
                .......        
            */
            int y0 = a0;
            int x0 = b0;
            int x1 = b1;
            for (int j = y0; j < y0 + 2*nBuffer; ++j) {
                int j2 = j - y0;
                for (int i = x0; i < x1; ++i) {
                    int i2 = i - x0;
                    int pixIdx1 = (j * w1) + i;
                    int pixIdx2 = (j2 * w2) + i2;
                    int l1 = labels1[pixIdx1];
                    int l2 = labels2[pixIdx2];

                    TIntSet set0 = overlap2.get(l2);
                    if (set0 == null) {
                        set0 = new TIntHashSet();
                        overlap2.put(l2, set0);
                    }
                    set0.add(l1);

                    TIntSet set2 = map2.get(l2);
                    if (set2 == null) {
                        set2 = new TIntHashSet();
                        map2.put(l2, set2);
                    }
                    set2.add(pixIdx2);

                    TIntSet set1 = map1.get(l1);
                    if (set1 == null) {
                        set1 = new TIntHashSet();
                        map1.put(l1, set1);
                    }
                    set1.add(pixIdx2);
                }
            }
        }
        
        // --------
        
        // key=original label2, value=new label2 assignment
        TIntIntMap labels2Map = new TIntIntHashMap();
        
        TIntObjectIterator<TIntSet> iter = overlap2.iterator();
        for (int i = 0; i < overlap2.size(); ++i) {
            iter.advance();
            int l2 = iter.key();
            TIntSet l1Set = iter.value();
            int l1Assign = -1;
            if (l1Set.size() > 1) {
                int maxL1Sz = Integer.MIN_VALUE;
                int maxL1 = -1;
                TIntIterator iter2 = l1Set.iterator();
                while (iter2.hasNext()) {
                    int l1T = iter2.next();
                    int sz = map1.get(l1T).size();
                    if (sz > maxL1Sz) {
                        maxL1Sz = sz;
                        maxL1 = l1T;
                    }
                }
                l1Assign = maxL1;
            } else {
                l1Assign = l1Set.iterator().next();
            }
            labels2Map.put(l2, l1Assign);
            if (l1Assign > maxLabel) {
                maxLabel = l1Assign;
            }
        }
        
        // -- iterate over map2, skipping existing in labels2Map
        //    else, assign value of maxValue+1
        iter = map2.iterator();
        for (int i = 0; i < map2.size(); ++i) {
            iter.advance();
            int l2Orig = iter.key();
            if (!labels2Map.containsKey(l2Orig)) {
                maxLabel++;
                labels2Map.put(l2Orig, maxLabel);
            }
        }
        
        // reassign labels2 using labels2Map
        for (int i = 0; i < labels2.length; ++i) {
            int l2 = labels2[i];
            labels2[i] = labels2Map.get(l2);
        }
        
        return maxLabel;
    }

    private int reassignAboveMaxLabel(int[] labels, 
        int maxLabel) {
        
        // key=original label2, value=new label2 assignment
        TIntIntMap labelsMap = new TIntIntHashMap();
                
        for (int i = 0; i < labels.length; ++i) {
            int l2 = labels[i];
            if (!labelsMap.containsKey(l2)) {
                maxLabel++;
                labelsMap.put(l2, maxLabel);
            }
            labels[i] = labelsMap.get(l2);
        }
        
        return maxLabel;
    }
    
}
