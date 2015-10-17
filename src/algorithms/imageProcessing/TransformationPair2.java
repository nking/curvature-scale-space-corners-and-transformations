package algorithms.imageProcessing;

/**
 *
 * @author nichole
 */
public class TransformationPair2 {
    
    private final int corner1Curve1Index;
    private final int corner1Curve2Index;
    private final int corner2Curve1Index;
    private final int corner2Curve2Index;
    
    private TransformationParameters params;
    
    private boolean someScaleAreSmallerThanOne = false;
    
    /**
     * object used to track visits and help choose the next contour.
     */
    private NextCorner nextCorner;
   
    public TransformationPair2(int indexForCorner1Curve1, 
        int indexForCorner1Curve2, int indexForCorner2Curve1, 
        int indexForCorner2Curve2) {
        
        this.corner1Curve1Index = indexForCorner1Curve1;
        
        this.corner1Curve2Index = indexForCorner1Curve2;
        
        this.corner2Curve1Index = indexForCorner2Curve1;
        
        this.corner2Curve2Index = indexForCorner2Curve2;
    }
    
    /**
     * @return the params
     */
    public TransformationParameters getTransformationParameters() {
        return params;
    }

    /**
     * @param theParameters
     */
    public void setTransformationParameters(TransformationParameters theParameters) {
        this.params = theParameters;
    }

    /**
     * @return the nextContour
     */
    public NextCorner getNextCorner() {
        return nextCorner;
    }

    /**
     * @param theNextCorner the theNextCorner to set
     */
    public void setNextCorner(NextCorner theNextCorner) {
        this.nextCorner = theNextCorner;
    }

    public void setSomeScaleAreSmallerThanOne() {
        someScaleAreSmallerThanOne = true;
    }
    
    public boolean scaleIsPossiblyAmbiguous() {
        return someScaleAreSmallerThanOne;
    }

}
