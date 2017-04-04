package algorithms.compGeometry.clustering.twopointcorrelation;

import java.io.IOException;

/**
 * interface for classes which supply an estimate of the background void density
 *
 * @author nichole
 */
public interface IPointBackgroundStats {

    /**
     *
     * @return
     * @throws IOException
     */
    public String persistIndexer() throws IOException;

    /**
     *
     * @param turnDebugOn
     */
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
    public float getBackgroundDensity();

    /**
     * @return the backgroundSurfaceDensityError
     */
    public float getBackgroundDensityError();
    
    /**
     *
     * @return
     */
    public long approximateMemoryUsed();
}