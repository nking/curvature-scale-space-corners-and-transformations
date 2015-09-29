package algorithms.compGeometry.convexHull;

/**
 *
 * @author nichole
 */
public class GrahamScanTooFewPointsException extends Exception {

    /**
     *
     */
    protected static final long serialVersionUID = 12345678;

    GrahamScanTooFewPointsException(String errorMessage) {
        super(errorMessage);
    }

}
