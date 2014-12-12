package algorithms.compGeometry.convexHull;

/**
 *adapted from 
 * https://code.google.com/p/two-point-correlation/source/browse/src/test/java/algorithms/compGeometry/convexHull/
 * under MIT License (MIT), Nichole King 2013
 */
public class GrahamScanTooFewPointsException extends Exception {

    protected static final long serialVersionUID = 12345678;

    GrahamScanTooFewPointsException(String errorMessage) {
        super(errorMessage);
    }

}
