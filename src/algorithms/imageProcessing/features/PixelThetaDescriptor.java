package algorithms.imageProcessing.features;

import algorithms.util.AngleUtil;
import algorithms.misc.MiscMath;
import java.util.Arrays;

/**
 *
 * @author nichole
 */
public class PixelThetaDescriptor extends ThetaDescriptor {
    
    protected final int[] a;
    
    protected float sumSquaredError = Float.NaN;
        
    /**
     * the index within array a that the central pixel
     * value is stored in.
     */
    protected final int centralIndex;
    
    public PixelThetaDescriptor(int[] intensities, int centralPixelIndex) {
        this.a = intensities;
        this.centralIndex = centralPixelIndex;
    }
    
    @Override
    public float calculateDifference(ThetaDescriptor otherDesc) {
                
        return calculateSSD(otherDesc);
    }
    
    @Override
    public float calculateError() {
        return sumSquaredError();
    }
    
    @Override
    public float calculateSSD(IDescriptor otherDesc) {
        
        if (otherDesc == null) {
            throw new IllegalArgumentException("otherDesc cannot be null");
        }
        
        if (!(otherDesc instanceof PixelThetaDescriptor)) {
            throw new IllegalArgumentException(
            "otherDesc has to be type ThetaDescriptor");
        }
        
        PixelThetaDescriptor other = (PixelThetaDescriptor)otherDesc;
        
        if (this.a.length != other.a.length) {
            throw new IllegalArgumentException(
            "this and other arrays must have the same lengths");
        }
         
        float ssd = MiscMath.calculateAngular360SSD(a, other.a, sentinel);
                
        return ssd;
    }

    @Override
    public float sumSquaredError() {
        
        if (!Float.isNaN(sumSquaredError)) {
            return sumSquaredError;
        }
        
        int cIdx = centralIndex;
                
        int vc = a[cIdx];
        
        if (vc == sentinel) {
            // sometimes arrive here because the gradient for center pixel
            // was too low intensity, so workaround is to reassign a nearby value
            // if possible
            int range = Math.min(((a.length - 1) - centralIndex), centralIndex);
            for (int dIdx = 1; dIdx < range; ++dIdx) {
                if (a[cIdx + dIdx] != sentinel) {
                    cIdx = cIdx + dIdx;
                    break;
                } else if (a[cIdx - dIdx] != sentinel) {
                    cIdx = cIdx - dIdx;
                    break;
                }
            }
        }
        
        float sqErr = MiscMath.sumSquaredAngular360Error(a, sentinel, cIdx);
            
        this.sumSquaredError = sqErr;
        
        return sumSquaredError;
    }
    
    public float[] calculateMeanAndStDev(ThetaDescriptor otherDesc) {
        
        if (otherDesc == null) {
            throw new IllegalArgumentException("otherDesc cannot be null");
        }
        
        if (!(otherDesc instanceof PixelThetaDescriptor)) {
            throw new IllegalArgumentException(
            "otherDesc has to be type ThetaDescriptor");
        }
        
        PixelThetaDescriptor other = (PixelThetaDescriptor)otherDesc;
        
        if (this.a.length != other.a.length) {
            throw new IllegalArgumentException(
            "this and other arrays must have the same lengths");
        }
        
        int[] diffs = new int[a.length];
        
        for (int i = 0; i < diffs.length; ++i) {
            
            if (a[i] == sentinel || other.a[i] == sentinel) {
                
                diffs[i] = sentinel;
                
            } else {
                
                diffs[i] = Math.round(AngleUtil.getAngleDifference(a[i], 
                    other.a[i]));
            }
        }
         
        float[] mnAndStDev = MiscMath.getAvgAndStDevIgnoreForSentinel(diffs,
            diffs.length, sentinel);
        
        return mnAndStDev;
    }

    @Override
    public float calculateCosineSimilarity(IDescriptor otherDesc) {
        
        if (otherDesc == null) {
            throw new IllegalArgumentException("otherDesc cannot be null");
        }
        
        if (!(otherDesc instanceof PixelThetaDescriptor)) {
            throw new IllegalArgumentException(
            "otherDesc has to be type ThetaDescriptor");
        }
        
        PixelThetaDescriptor other = (PixelThetaDescriptor)otherDesc;
        
        if (this.a.length != other.a.length) {
            throw new IllegalArgumentException(
            "this and other arrays must have the same lengths");
        }
         
        float cSim = MiscMath.calculateCosineSimilarity(a, other.a, sentinel);
                
        return cSim;
    }
    
    @Override
    public int getCentralIndex() {
        return centralIndex;
    }

    @Override
    public String toString() {
        return Arrays.toString(a);
    }
    
}
