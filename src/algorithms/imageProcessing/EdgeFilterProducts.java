package algorithms.imageProcessing;

import algorithms.util.PairInt;
import java.util.Map;
import java.util.Set;

/**
 * a wrapper to hold the results of an edge image filter, specifically the
 * gradient X, gradient Y, the combined gradients XY, and the theta image.
 * 
 * @author nichole
 */
public class EdgeFilterProducts {
    
    private GreyscaleImage gradientXY = null;
        
    private GreyscaleImage gradientX = null;
        
    private GreyscaleImage gradientY = null;
        
    private GreyscaleImage theta = null;
    
    /**
     * hough lines as a map with key = hough line set of of points and value =
     * theta and radius.
     */
    private Map<Set<PairInt>, PairInt> houghLines = null;

    /**
     * @return the gradientXY
     */
    public GreyscaleImage getGradientXY() {
        return gradientXY;
    }

    /**
     * @param gradientXY the gradientXY to set
     */
    public void setGradientXY(GreyscaleImage gradientXY) {
        this.gradientXY = gradientXY;
    }

    /**
     * @return the gradientX
     */
    public GreyscaleImage getGradientX() {
        return gradientX;
    }

    /**
     * @param gradientX the gradientX to set
     */
    public void setGradientX(GreyscaleImage gradientX) {
        this.gradientX = gradientX;
    }

    /**
     * @return the gradientY
     */
    public GreyscaleImage getGradientY() {
        return gradientY;
    }

    /**
     * @param gradientY the gradientY to set
     */
    public void setGradientY(GreyscaleImage gradientY) {
        this.gradientY = gradientY;
    }

    /**
     * @return the theta
     */
    public GreyscaleImage getTheta() {
        return theta;
    }

    /**
     * @param theta the theta to set
     */
    public void setTheta(GreyscaleImage theta) {
        this.theta = theta;
    }
    
    /**
     * hough lines as a map with key = hough line set of of points and value =
     * theta and radius.
     * @return 
     */
    public Map<Set<PairInt>, PairInt> getHoughLines() {
        return houghLines;
    }
    
    /**
     * hough lines as a map with key = hough line set of of points and value =
     * theta and radius.
     * @param theHoughLines
     */
    public void setHoughLines( Map<Set<PairInt>, PairInt> theHoughLines) {
        this.houghLines = theHoughLines;
    }

}
