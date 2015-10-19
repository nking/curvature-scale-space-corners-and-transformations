package algorithms.imageProcessing;

/**
 *
 * @author nichole
 */
public class TransformationPair2 {
    
    private final CornerRegion corner1Curve1;
    private final CornerRegion corner1Curve2;
    private final CornerRegion corner2Curve1;
    private final CornerRegion corner2Curve2;
    
    private TransformationParameters params;
    
    private boolean someScaleAreSmallerThanOne = false;
    
    /**
     * object used to track visits and help choose the next contour.
     */
    private NextCorner nextCorner;
   
    public TransformationPair2(CornerRegion corner1Curve1, 
        CornerRegion corner1Curve2, CornerRegion corner2Curve1, 
        CornerRegion corner2Curve2) {
        
        this.corner1Curve1 = corner1Curve1;
        
        this.corner1Curve2 = corner1Curve2;
        
        this.corner2Curve1 = corner2Curve1;
        
        this.corner2Curve2 = corner2Curve2;
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
