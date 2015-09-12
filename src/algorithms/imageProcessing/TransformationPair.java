package algorithms.imageProcessing;

/**
 *
 * @author nichole
 */
public class TransformationPair {
    
    private final int contourIndex1;
    private final int contourIndex2;
    
    private double scale;
    private double shift;
    
    private boolean someScaleAreSmallerThanOne = false;
    
    /**
     * object used to track visits and help choose the next contour.
     * NOTE: this may be null until first search to reduce memory used.
     */
    private NextContour nextContour;
   
    public TransformationPair(int indexForContour1, int indexForContour2) {
        
        this.contourIndex1 = indexForContour1;
        
        this.contourIndex2 = indexForContour2;
        
    }

    /**
     * @return the contourIndex1
     */
    public int getContourIndex1() {
        return contourIndex1;
    }

    /**
     * @return the contourIndex2
     */
    public int getContourIndex2() {
        return contourIndex2;
    }

    /**
     * @return the scale
     */
    public double getScale() {
        return scale;
    }

    /**
     * @param scale the scale to set
     */
    public void setScale(double scale) {
        this.scale = scale;
    }

    /**
     * @return the shift
     */
    public double getShift() {
        return shift;
    }

    /**
     * @param shift the shift to set
     */
    public void setShift(double shift) {
        this.shift = shift;
    }

    /**
     * @return the nextContour
     */
    public NextContour getNextContour() {
        return nextContour;
    }

    /**
     * @param nextContour the nextContour to set
     */
    public void setNextContour(NextContour nextContour) {
        this.nextContour = nextContour;
    }

    public void setSomeScaleAreSmallerThanOne() {
        someScaleAreSmallerThanOne = true;
    }
    
    public boolean scaleIsPossiblyAmbiguous() {
        return someScaleAreSmallerThanOne;
    }

}
