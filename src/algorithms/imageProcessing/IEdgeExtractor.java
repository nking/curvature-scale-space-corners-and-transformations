package algorithms.imageProcessing;

import algorithms.util.PairIntArray;
import java.util.List;

/**
 * interface for edge extractor. the input should be an image which has
 * already been processed by the CannyEdgeFilter to have lines that are only
 * 1 pixel wide.
 * 
 * @author nichole
 */
public interface IEdgeExtractor {

    /**
     * find the edges and return as a list of points.  The method uses a
     * DFS search through all points in the image with values > 0 to link
     * adjacent sequential points into edges.
     * As a side effect, the method also populates
     * member variables edgeJunctionMap and outputIndexLocatorForJunctionPoints.
     *
     * @return
     */
    List<PairIntArray> findEdges();

    GreyscaleImage getImage();

    void overrideEdgeSizeLowerLimit(int length);
    
    public void removeShorterEdges(boolean doRemove);
    
}
