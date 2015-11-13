package algorithms.imageProcessing;

/**
 *
 * @author nichole
 */
public class TransformationPair2<T extends CornerRegion> {
    
    private int cornerListIndex1 = -1;
    private int cornerListIndex2 = -1;
    
    private final T corner1Curve1;
    private final T corner1Curve2;
    private final T corner2Curve1;
    private final T corner2Curve2;
    
    private TransformationParameters params = null;
    
    private boolean someScaleAreSmallerThanOne = false;
    
    private double costAsSSD = 0;
    
    private double costAsDist = 0;
    
    /**
     * object used to track visits and help choose the next contour.
     */
    private NextCorner<T> nextCorner;
   
    public TransformationPair2(T corner1Curve1, T corner1Curve2, T corner2Curve1, 
        T corner2Curve2) {
        
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
    
    public void addToCostAsSSD(double cost) {
        costAsSSD += cost;
    }
    
    public void addToCostAsDistance(double cost) {
        costAsDist += cost;
    }
    
    public double getCostAsSSD() {
        return costAsSSD;
    }
    
    public double getCostAsDistance() {
        return costAsDist;
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
    public NextCorner<T> getNextCorner() {
        return nextCorner;
    }

    /**
     * @param theNextCorner the theNextCorner to set
     */
    public void setNextCorner(NextCorner<T> theNextCorner) {
        this.nextCorner = theNextCorner;
    }

    public void setSomeScaleAreSmallerThanOne() {
        someScaleAreSmallerThanOne = true;
    }
    
    public boolean scaleIsPossiblyAmbiguous() {
        return someScaleAreSmallerThanOne;
    }
    
    public void setCornerListIndex1(int idx) {
        this.cornerListIndex1 = idx;
    }

    public void setCornerListIndex2(int idx) {
        this.cornerListIndex2 = idx;
    }
    
    public int getCornerListIndex1() {
        return this.cornerListIndex1;
    }

    public int getCornerListIndex2() {
        return this.cornerListIndex2;
    }
    
    @Override
    public String toString() {
        
        StringBuilder sb = new StringBuilder();
        sb.append("costAsSSD=").append(Double.toString(costAsSSD))
            .append(" costAsDist=").append(Double.toString(costAsDist)).append(" ");
        
        if (params != null) {
            sb.append(params.toString());
        }
        
        if (nextCorner != null) {
            int n = nextCorner.getMatchedCornerIndexes1().size();
            sb.append(" nMatches=").append(Integer.toString(n)).append(" ");
        }
        
        return sb.toString();
    }

}
