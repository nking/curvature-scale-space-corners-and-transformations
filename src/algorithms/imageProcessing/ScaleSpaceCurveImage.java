package algorithms.imageProcessing;

import java.util.Arrays;
import java.util.logging.Logger;

/**
 * Holds the inflection points of a curve convolved over many sigma until
 * there are no new inflection points.  the y axis of the image is sigma and
 * can be retrieved with getImageSigmas().  the x axis of the image is
 * the scale free parameter t which is the index of the curve that the 
 * zero-crossing occurs normalized by the number of points in the curve.
 * @author nichole
 */
public class ScaleSpaceCurveImage {
    
    private final int nSigmaLevels;

    private final float[] imageSigmas;
    
    private final float[][] scaleSpaceImage;
    
    private final int[][] xCoords;
    
    private final int[][] yCoords;
    
    private int edgeNumber = -1;
    
    private int edgeSize = 0;
    
    private Logger log = Logger.getLogger(this.getClass().getName());
    
    public ScaleSpaceCurveImage(int nSigmaLevels) {

        this.nSigmaLevels = nSigmaLevels;
        
        imageSigmas = new float[nSigmaLevels];
        
        scaleSpaceImage = new float[nSigmaLevels][];
                
        xCoords = new int[nSigmaLevels][];
        
        yCoords = new int[nSigmaLevels][];
    }
    
    public void setEdgeNumber(int number) {
        edgeNumber = number;
    }
    
    public void setEdgeSize(int size) {
        edgeSize = size;
    }
    
    public int getEdgeSize() {
        return edgeSize;
    }
    
    public void setRow(int sigmaIndex, float[] scaleFreeZeroCrossings) {
        if (sigmaIndex < 0 || sigmaIndex > (nSigmaLevels - 1)) {
            throw new IllegalStateException("sigmaIndex is out of bounds");
        }
        
        scaleSpaceImage[sigmaIndex] = scaleFreeZeroCrossings;
    }
    
    public void setXYCoords(int sigmaIndex, int[] theXCoords, int[] theYCoords) {
        
        if (sigmaIndex < 0 || sigmaIndex > (nSigmaLevels - 1)) {
            throw new IllegalStateException("sigmaIndex is out of bounds");
        }
        
        this.xCoords[sigmaIndex] = theXCoords;
        
        this.yCoords[sigmaIndex] = theYCoords;
    }
    
    public void setSigma(int sigmaIndex, float sigma) {
        if (sigmaIndex < 0 || sigmaIndex > (nSigmaLevels - 1)) {
            throw new IllegalStateException("sigmaIndex is out of bounds");
        }
        
        imageSigmas[sigmaIndex] = sigma;
    }
    
    public float[] getImageSigmas() {
        return imageSigmas;
    }
    
    public float[][] getScaleSpaceImage() {
        return scaleSpaceImage;
    }
    
    public int getEdgeNumber() {
        return edgeNumber;
    }
    
    public int getXCoord(int sigmaIndex, int tIndex) {
        if (sigmaIndex < 0 || sigmaIndex > (nSigmaLevels - 1)) {
            throw new IllegalStateException("sigmaIndex is out of bounds");
        }
        if (tIndex < 0 || tIndex > (xCoords[sigmaIndex].length - 1)) {
            throw new IllegalStateException("tIndex is out of bounds");
        }
        return xCoords[sigmaIndex][tIndex];
    }
    
    public int getYCoord(int sigmaIndex, int tIndex) {
        if (sigmaIndex < 0 || sigmaIndex > (nSigmaLevels - 1)) {
            throw new IllegalStateException("sigmaIndex is out of bounds");
        }
        if (tIndex < 0 || tIndex > (yCoords[sigmaIndex].length - 1)) {
            throw new IllegalStateException("tIndex is out of bounds");
        }
        return yCoords[sigmaIndex][tIndex];
    }
   
    public ScaleSpaceCurveImage copy() {
        
        ScaleSpaceCurveImage c = new ScaleSpaceCurveImage(nSigmaLevels);
        
        c.setEdgeNumber(edgeNumber);
        
        c.setEdgeSize(edgeSize);
        
        System.arraycopy(imageSigmas, 0, c.imageSigmas, 0, imageSigmas.length);
        
        for (int j = 0; j < imageSigmas.length; j++) {
            
            int n = scaleSpaceImage[j].length;
            
            c.setRow(j, Arrays.copyOf(scaleSpaceImage[j], n));
            
            c.setXYCoords(j, Arrays.copyOf(xCoords[j], n), 
                Arrays.copyOf(yCoords[j], n));
        }
        
        return c;
    }
}
