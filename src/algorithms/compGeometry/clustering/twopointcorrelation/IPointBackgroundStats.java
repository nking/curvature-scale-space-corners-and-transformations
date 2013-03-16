package algorithms.compGeometry.clustering.twopointcorrelation;

import java.io.IOException;

/**
 * interface for classes which supply an estimate of the background surface density
 *
 * @author nichole
 */
public interface IPointBackgroundStats {

    public String persistIndexer() throws IOException;

    public void setDebug(boolean turnDebugOn);

    /**
     * calculate stats for the largest minima areas.
     *
     * @throws TwoPointVoidStatsException
     */
    public void calc() throws TwoPointVoidStatsException;

    /**
     * use a previously persisted set of two-point background data to run the code.
     *
     * @param persistedTwoPointBackgroundFileName
     * @throws TwoPointVoidStatsException
     * @throws IOException
     */
    public void calc(String persistedTwoPointBackgroundFileName) throws TwoPointVoidStatsException, IOException;

    /**
     * @return the backgroundSurfaceDensity
     */
    public float getBackgroundSurfaceDensity();

    /**
     * @return the backgroundSurfaceDensityError
     */
    public float getBackgroundSurfaceDensityError();

}
