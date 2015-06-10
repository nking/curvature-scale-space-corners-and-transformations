package algorithms.imageProcessing;

/**
 *
 * @author nichole
 */
public interface ILineThinner {
    
    /**
     * set the gradientXY image to be an edge guide image as a final stage
     * in applyFilter.  Use of this is discouraged as the changes are few
     * and do not improve the corners or contours.
     * @param gradientXY 
     */
    public void setEdgeGuideImage(final GreyscaleImage gradientXY);
    
    /**
     * apply the filter
     * @param image the input image (it's modified internally)
     */
    public void applyFilter(GreyscaleImage image);

    /**
     * turn on the extra output to help debugging.  this includes plot displays
     * with awt, html output and additional logging
     * 
     * @param setToDebug 
     */
    public void setDebug(boolean setToDebug);

    /**
    if the filter has special handling for line drawings, this applies that.
     */
    public void useLineDrawingMode();
    
}
