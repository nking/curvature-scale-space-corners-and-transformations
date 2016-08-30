package algorithms.imageProcessing.segmentation;

import algorithms.imageProcessing.ImageExt;
import algorithms.misc.MiscMath;

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
        
        if (binX < 2*nBuffer) {
            //TODO: handle this
        }
        if (binY < 2*nBuffer) {
            //TODO: handle this
        }
        
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
            int stopX = startX + binX + nBuffer;
            if (stopX > w) {
                stopX = w;
            }
            for (int j = 0; j < nY; ++j) {
                int startY = j * binY - nBuffer;
                if (startY < 0) {
                    startY = 0;
                }
                int stopY = startY + binY + nBuffer;
                if (stopY > h) {
                    stopY = h;
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
                
                int[] labels2 = run(img2, sLabels);
                
                int x0 = startX;
                int x1 = stopX;
                int y0 = startY;
                int y1 = stopY;
                
                // -- merge any results with overlap to left
                // or overlap below in cLabels
                if (i > 0) {
                    // merge the overlapping buffer region to left
                    maxLabel2 = mergeLeft(cLabels, x0, y0, y1, 
                        nBuffer, maxLabel2);
                    x0 += nBuffer;
                }
                if (j > 0) {
                    // merge the overlapping buffer region below
                    maxLabel2 = mergeBelow(cLabels, y0, x0, x1, 
                        nBuffer, maxLabel2);
                    y0 += nBuffer;
                }
                
                // TODO: increment the values in sLabels (exclude the 
                // buffer regions) by creating map w/ key=label
                // value = pixIdxs, then increment labels past
                // maxLabel2
                
                
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

    private int mergeLeft(int[] labels, int x0, 
        int y0, int y1, int nBuffer, int maxLabel) {
        
        /*
                | binX   |    i=1
                |      [ | ]       
                         x0
        */

        throw new UnsupportedOperationException(
            "Not implemented yet."); 
    }
    
    private int mergeBelow(int[] labels, int y0, 
        int x0, int x1, int nBuffer, int maxLabel) {
        
        /*
        
                .......
                ------- y0
                .......        
        
        */

        throw new UnsupportedOperationException(
            "Not implemented yet."); 
    }
   
}
