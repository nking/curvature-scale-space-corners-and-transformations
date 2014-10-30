package algorithms.imageProcessing;

/**
 * class holding data for a contour for an edge in a scale space image.
 * The class holds the peak of the contour.  It also carries information to 
 * help identify the peak in the original image from which the edge was 
 * extracted.  Note, that might be 2 points (the left and right of the peak).
 * Also note that the peak sigma's real location is relative to an image
 * which has been convolved with a Gaussian and so the exact location is not
 * the same.  The coordinates provided are useful for debugging.
 * 
 * @author nichole
 */
public class CurvatureScaleSpaceContour {
    
    protected final PairFloat peak;
    
    private int edgeNumber = -1;
    
    /**
     * x, y indexes of inflection points for the peak of the contour.
     * the indexes are relative to the coordinates in the edge found by
     * edgeNumber, owned by another class instance.
     * For example, CurvatureScaleSpaceImageMaker has a member called
     * closedCurves which are closed curve edges.  The return of
     * createInflectionContours() are relative to that.
     * Expect that peakIndexes will have 1 or 2 indexes in it only
     * for a single point peak, or a peak defined by a left and right point.
     * CurvatureScaleSpaceImagePoint holds the peakIndex that can be used
     * with the original edge to get the x, y digital image 
     * coordinates when needed.
     */
    private CurvatureScaleSpaceImagePoint[] peakDetailPoints = null;
    
    public CurvatureScaleSpaceContour(float sigma, float t) {
        peak = new PairFloat(sigma, t);
    }
    
    public float getPeakSigma() {
        return peak.getX();
    }
    
    public float getPeakScaleFreeLength() {
        return peak.getY();
    }

    public void setEdgeNumber(int number) {
        this.edgeNumber = number;
    }
    
    public void setPeakDetails(CurvatureScaleSpaceImagePoint[] points) {
        peakDetailPoints = points;
    }
    
    public int getEdgeNumber() {
        return this.edgeNumber;
    }
    
    public CurvatureScaleSpaceImagePoint[] getPeakDetails() {
        return this.peakDetailPoints;
    }

    @Override
    public String toString() {
        
        StringBuilder sb = new StringBuilder();
        
        sb.append("edgeNumber=").append(edgeNumber);
        
        sb.append(" sigma=").append(getPeakSigma()).append(" scaleFreeLength=")
            .append(getPeakScaleFreeLength());
        
        for (int i = 0; i < peakDetailPoints.length; i++) {
            sb.append("[").append(peakDetailPoints[i].toString()).append("] ");
        }
        
        return sb.toString();
    }
    
}
