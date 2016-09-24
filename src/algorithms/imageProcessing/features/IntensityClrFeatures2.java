package algorithms.imageProcessing.features;

import algorithms.compGeometry.RotatedOffsets;
import algorithms.imageProcessing.Image;
import algorithms.imageProcessing.ImageExt;
import algorithms.misc.Misc;
import algorithms.util.PairInt;
import java.awt.Color;
import java.util.HashMap;
import java.util.Map;

/**
 * a more recent version of IntensityClrFeatures, but
 * designed to use different data.
 * (IntensityClrFeatures may be deprecated soon...the features and matching
 * packages will be refactored after a few more additions).
 * 
 * @author nichole
 */
public class IntensityClrFeatures2 {
    
    protected final int bHalfW;
    
    protected final float[][] xyOffsets;
    
    protected final RotatedOffsets rotatedOffsets;
    
    protected final ImageExt img;
    
    protected final boolean useNormalizedIntensities = false;
    
    protected final boolean useBinnedCellIntensities = true;
    
    /**
    key = pixel coordinates of center of frame;
    value = map with key = rotation (in degrees) of the frame and value = the
        extracted descriptor
    NOTE that the IntensityDescriptor is ClrIntensityDescriptor, but instead
    of r, g, b is populated with h, s, v scaled by a factor to make them 
    integers.
    */
    protected Map<PairInt, Map<Integer, IntensityDescriptor>> 
        intensityHBlocks = new HashMap<PairInt, Map<Integer, IntensityDescriptor>>();
    
    protected Map<PairInt, Map<Integer, IntensityDescriptor>> 
        intensitySBlocks = new HashMap<PairInt, Map<Integer, IntensityDescriptor>>();
    
    protected Map<PairInt, Map<Integer, IntensityDescriptor>> 
        intensityVBlocks = new HashMap<PairInt, Map<Integer, IntensityDescriptor>>();
    
    /**
     * note, RotatedOffsets caches the results of expensive transcendental
     * operations and is used in a shared manner, that is,
     * it is a singleton.  caution is needed in it's use because it 
     * contains unsynchronized unguarded cache variables that could be 
     * corrupted by multi-threaded use.  The design is meant to keep access
     * to it fast.
     * 
     * @param img
     * @param blockHalfWidths the descriptor block size.  A value of 5 is
     * recommended.
     * @param rotatedOffsets 
     */
    public IntensityClrFeatures2(final ImageExt img, final int blockHalfWidths, 
        final RotatedOffsets rotatedOffsets) {
        
        this.bHalfW = blockHalfWidths;
        this.xyOffsets = Misc.createNeighborOffsets(bHalfW);
        this.rotatedOffsets = rotatedOffsets;
        this.img = img;
    }
    
    public IntensityDescriptor[] extractIntensityHSV(
        final int xCenter, final int yCenter, final int rotation) {
    
        checkBounds(img, xCenter, yCenter);
        
        PairInt p = new PairInt(xCenter, yCenter);
        
        Integer rotationKey = Integer.valueOf(rotation);
        
        Map<Integer, IntensityDescriptor> descriptorsH = intensityHBlocks.get(p);
        Map<Integer, IntensityDescriptor> descriptorsS = intensitySBlocks.get(p);
        Map<Integer, IntensityDescriptor> descriptorsV = intensityVBlocks.get(p);
        
        IntensityDescriptor descriptorH = null;
        
        if (descriptorsH != null) {
            
            descriptorH = descriptorsH.get(rotationKey);
            
            if (descriptorH != null) {
                // if H is present, then so are S and V
                IntensityDescriptor descriptorS = descriptorsS.get(rotationKey);
                IntensityDescriptor descriptorV = descriptorsV.get(rotationKey);
                
                return new IntensityDescriptor[]{descriptorH, descriptorS,
                    descriptorV};
            }
        } else {
            descriptorsH = new HashMap<Integer, IntensityDescriptor>();
            intensityHBlocks.put(p, descriptorsH);
            descriptorsS = new HashMap<Integer, IntensityDescriptor>();
            intensitySBlocks.put(p, descriptorsS);
            descriptorsV = new HashMap<Integer, IntensityDescriptor>();
            intensityVBlocks.put(p, descriptorsV);
        }
        
        IntensityDescriptor[] descriptors 
            = extractHSVIntensitiesForCells(xCenter, yCenter, rotation);
       
        if (descriptors == null) {
            return null;
        }
       
        descriptorH = descriptors[0];
        IntensityDescriptor descriptorS = descriptors[1];
        IntensityDescriptor descriptorV = descriptors[2];
        
        if (useNormalizedIntensities && (descriptorH != null)) {
            descriptorH.applyNormalization();
            descriptorS.applyNormalization();
            descriptorV.applyNormalization();
        }
        
        if (descriptorH != null) {
            descriptorsH.put(rotationKey, descriptorH);
            descriptorsS.put(rotationKey, descriptorS);
            descriptorsV.put(rotationKey, descriptorV);
        }
        
        return descriptors;        
    }
    
    /**
     * note, the internal auto-correlation is determined w.r.t. the 2nd
     * set of descriptors which might be helpful to know when comparing
     * results.
     * @param desc1_h
     * @param desc1_s
     * @param desc1_v
     * @param x1
     * @param y1
     * @param useTop1
     * @param desc2_h
     * @param desc2_s
     * @param desc2_v
     * @param x2
     * @param y2
     * @param useTop2
     * @return 
     */
    public static FeatureComparisonStat calculateHalfStats(
        IntensityDescriptor desc1_h, IntensityDescriptor desc1_s, 
        IntensityDescriptor desc1_v, int x1, int y1, boolean useTop1, 
        IntensityDescriptor desc2_h, IntensityDescriptor desc2_s, 
        IntensityDescriptor desc2_v, int x2, int y2, boolean useTop2) {
        
        if (desc1_h == null) {
            throw new IllegalArgumentException("desc1_h cannot be null");
        }
        if (desc1_s == null) {
            throw new IllegalArgumentException("desc1_s cannot be null");
        }
        if (desc1_v == null) {
            throw new IllegalArgumentException("desc1_v cannot be null");
        }
        if (desc2_h == null) {
            throw new IllegalArgumentException("desc2_h cannot be null");
        }
        if (desc2_s == null) {
            throw new IllegalArgumentException("desc2_s cannot be null");
        }
        if (desc2_v == null) {
            throw new IllegalArgumentException("desc2_v cannot be null");
        }
        
        int[] indexes1 = useTop1 ? 
            ((GsIntensityDescriptor)desc1_h).getUpperHalfIndexes() :
            ((GsIntensityDescriptor)desc1_h).getLowerHalfIndexes();
        
        int[] indexes2 = useTop2 ? 
            ((GsIntensityDescriptor)desc2_h).getUpperHalfIndexes() :
            ((GsIntensityDescriptor)desc2_h).getLowerHalfIndexes();
        
        if (indexes1.length != indexes2.length) {
            throw new IllegalArgumentException(
            "indexes1 and indexes2 must be same length");
        }
        
        float sentinel = GsIntensityDescriptor.sentinel;
        
        float[] h1 = ((GsIntensityDescriptor)desc1_h).grey;
        float[] s1 = ((GsIntensityDescriptor)desc1_s).grey;
        float[] v1 = ((GsIntensityDescriptor)desc1_v).grey;
        
        float[] h2 = ((GsIntensityDescriptor)desc2_h).grey;
        float[] s2 = ((GsIntensityDescriptor)desc2_s).grey;
        float[] v2 = ((GsIntensityDescriptor)desc2_v).grey;
        
        int n = h1.length;
        int centralPixIdx1 = desc1_h.getCentralIndex();
        int centralPixIdx2 = desc2_h.getCentralIndex();
                
        float vcH1 = h1[centralPixIdx1];
        if (vcH1 == sentinel) {
            throw new IllegalStateException(
            "ERROR: the central value for the array is somehow sentinel");
        }
        //float vcA1 = a1[centralPixIdx1];
        //float vcB1 = b1[centralPixIdx1];
        
        float vcH2 = h2[centralPixIdx2];
        if (vcH2 == sentinel) {
            throw new IllegalStateException(
            "ERROR: the central value for the array is somehow sentinel");
        }
        float vcS2 = s2[centralPixIdx2];
        float vcV2 = v2[centralPixIdx2];
        
        assert(indexes1.length == indexes2.length);
                        
        int count = 0;
        //TODO: review the math for auto-correlation here since it is estimated 
        // from 2 points instead of 1
        double autoCorrel = 0;
        double deltaSum = 0;
        
        n = indexes1.length;
        
        for (int i = 0; i < n; ++i) {
            
            int idx1 = indexes1[i];
            int idx2 = indexes2[i];
            
            float vH1 = h1[idx1];
            if (vH1 == sentinel) {
                continue;
            }
            float vH2 = h2[idx2];
            if (vH2 == sentinel) {
                continue;
            }
            float vS1 = s1[idx1];
            float vV1 = v1[idx1];
            float vS2 = s2[idx2];
            float vV2 = v2[idx2];
            
            double diffH = vH1 - vH2;
            double diffS = vS1 - vS2;
            double diffV = vV1 - vV2;
            double delta = Math.sqrt(diffH * diffH + diffS * diffS + 
                diffV * diffV);
            
            deltaSum += (delta * delta);
                        
            diffH = vH1 - vcH2;
            diffS = vS1 - vcS2;
            diffV = vV1 - vcV2;
            double deltaC = Math.sqrt(diffH * diffH + diffS * diffS + 
                diffV * diffV);
            
            autoCorrel += (deltaC * deltaC);            
            
            count++;
        }
        
        autoCorrel /= (double)count;
        
        deltaSum /= (double)count;
         
        FeatureComparisonStat stat = new FeatureComparisonStat();
        stat.setImg1Point(new PairInt(x1, y1));
        stat.setImg2Point(new PairInt(x2, y2));
        stat.setSumIntensitySqDiff((float)deltaSum);
        stat.setImg2PointIntensityErr((float)autoCorrel);
        
        return stat;
    }
    
    protected void checkBounds(Image img, int x, int y) {
        if (!isWithinXBounds(img, x)) {
            throw new IllegalArgumentException("x is out of bounds of image");
        }
        if (!isWithinYBounds(img, y)) {
            throw new IllegalArgumentException("y is out of bounds of image");
        }
    }
    
    protected boolean isWithinXBounds(Image img, int x) {
        if (x < 0) {
            return false;
        }
        if (x > (img.getWidth() - 1)) {
            return false;
        }
        return true;
    }
    
    protected boolean isWithinYBounds(Image img, int y) {
        if (y < 0) {
            return false;
        }
        if (y > (img.getHeight() - 1)) {
            return false;
        }
        return true;
    }

    /**
     * extract the H, S, V, intensities from the image in 2X2 cells surrounding
     * (xCenter, yCenter) for 16 cells.
     * @param xCenter
     * @param yCenter
     * @param rotation in degrees
     * @param oInt represents whether to calculate O1, O2, or O3
     * @return
     */
    private IntensityDescriptor[] extractHSVIntensitiesForCells(int xCenter, 
        int yCenter, int rotation) {
        
        float sentinel = GsIntensityDescriptor.sentinel;
        
        float toIntFactor = 255.f;
        
        int cellDim = 2;        
        int nCellsAcross = 6;
        int nColsHalf = nCellsAcross / 2;
        int range0 = cellDim * nColsHalf;
        
        int[] xOffsets = rotatedOffsets.getXOffsets(rotation);
        int[] yOffsets = rotatedOffsets.getYOffsets(rotation);
        
        int w = img.getWidth();
        int h = img.getHeight();

        float[] outputH = new float[nCellsAcross * nCellsAcross];
        float[] outputS = new float[nCellsAcross * nCellsAcross];
        float[] outputV = new float[nCellsAcross * nCellsAcross];

        int n2 = cellDim * cellDim;
        
        //index of center pixel in array of descriptor:
        int centralPixelIndex = (nCellsAcross*nColsHalf) + nColsHalf;
        
        int count = 0;
        int idx = 0;
        for (int dx = -range0; dx < range0; dx += cellDim) {
            for (int dy = -range0; dy < range0; dy += cellDim) {
                boolean withinBounds = true;
                int cCount = 0;
                // ---- sum within the cell ----
                double rV = 0;
                double gV = 0;
                double bV = 0;
                for (int i = 0; i < n2; ++i) {
                    int xOff = xOffsets[idx];
                    int yOff = yOffsets[idx];
                    idx++;
                    int x = xOff + xCenter;
                    int y = yOff + yCenter;
                    if ((x < 0) || (x > (w - 1)) || (y < 0) || (y > (h - 1))) {
                        withinBounds = false;
                        break;
                    }
                    rV += img.getR(x, y);
                    gV += img.getG(x, y);
                    bV += img.getB(x, y);
                    cCount++;
                }
                if (!withinBounds || (cCount == 0)) {
                    if (count == centralPixelIndex) {
                        return null;
                    }
                    outputH[count] = sentinel;
                    outputS[count] = sentinel;
                    outputV[count] = sentinel;
                    count++;
                    continue;
                }
                rV /= (float) cCount;
                gV /= (float) cCount;
                bV /= (float) cCount;
                
                float[] hsb = new float[3];
                Color.RGBtoHSB((int)Math.round(rV), 
                    (int)Math.round(gV), (int)Math.round(bV), hsb);
                
                int hV = Math.round(toIntFactor * hsb[0]);
                int sV = Math.round(toIntFactor * hsb[1]);
                int vV = Math.round(toIntFactor * hsb[2]);
                
                outputH[count] = hV;
                outputS[count] = sV;
                outputV[count] = vV;
                count++;
            }
        }

        IntensityDescriptor descH = new GsIntensityDescriptor(outputH, 
            centralPixelIndex, nCellsAcross);
        IntensityDescriptor descS = new GsIntensityDescriptor(outputS, 
            centralPixelIndex, nCellsAcross);
        IntensityDescriptor descV = new GsIntensityDescriptor(outputV, 
            centralPixelIndex, nCellsAcross);
        
        return new IntensityDescriptor[]{descH, descS, descV};
    }
    
}
