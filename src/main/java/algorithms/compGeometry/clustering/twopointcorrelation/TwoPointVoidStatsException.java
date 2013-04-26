package algorithms.compGeometry.clustering.twopointcorrelation;

public class TwoPointVoidStatsException extends Exception {

    protected static final long serialVersionUID = 123456789;

    TwoPointVoidStatsException(String errorMessage) {
        super(errorMessage);
    }

    TwoPointVoidStatsException(Exception e) {
        super(e);
    }

}
