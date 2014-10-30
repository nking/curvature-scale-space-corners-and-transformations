package algorithms.imageProcessing;

/**
 *
 * @author nichole
 */
public class CurvatureScaleSpaceImagePoint {
    
    private int xCoord = -1;
    
    private int yCoord = -1;
    
    private float sigma = Float.MAX_VALUE;
    
    private float scaleFreeLength = Float.MAX_VALUE;

    public CurvatureScaleSpaceImagePoint(float theSigma, 
        float theScaleFreeLength, int theXCoordinate, int theYCoordinate) {
        
        xCoord = theXCoordinate;
        
        yCoord = theYCoordinate;
        
        sigma = theSigma;
        
        scaleFreeLength = theScaleFreeLength;
    }
    
    /**
     * @return the xCoord
     */
    public int getXCoord() {
        return xCoord;
    }
    
    /**
     * @return the yCoord
     */
    public int getYCoord() {
        return yCoord;
    }

    /**
     * @return the sigma
     */
    public float getSigma() {
        return sigma;
    }

    /**
     * @return the scaleFreeLength
     */
    public float getScaleFreeLength() {
        return scaleFreeLength;
    }

    /**
     * @param scaleFreeLength the scaleFreeLength to set
     */
    public void setScaleFreeLength(float scaleFreeLength) {
        this.scaleFreeLength = scaleFreeLength;
    }
    
    @Override
    public String toString() {
        
        StringBuilder sb = new StringBuilder();
        
        sb.append("sigma=").append(sigma)
            .append(" scaleFreeLength=").append(scaleFreeLength)
            .append(" (").append(xCoord).append(",").append(yCoord).append(")");
        
        return sb.toString();
    }

}
