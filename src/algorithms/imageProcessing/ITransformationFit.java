package algorithms.imageProcessing;

import java.util.List;

/**
 *
 * @author nichole
 */
public interface ITransformationFit {

    void calculateErrorStatistics();

    /**
     * @return the errors
     */
    List<Double> getErrors();

    List<Integer> getInlierIndexes();

    double getMeanError();

    /**
     * @return the nMatches
     */
    long getNMatches();

    /**
     * @return the nMaxMatchable
     */
    long getNMaxMatchable();

    /**
     * @return the stDevFromAvg
     */
    double getStDevFromMean();

    /**
     * @return the tolerance
     */
    double getTolerance();

    /**
     * @param theNumberOfPossibleMatches
     */
    void setNMaxMatchable(long theNumberOfPossibleMatches);
    
    public boolean isBetter(ITransformationFit other);
}
