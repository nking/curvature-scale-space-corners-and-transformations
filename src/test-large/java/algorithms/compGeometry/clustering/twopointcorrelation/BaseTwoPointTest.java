package algorithms.compGeometry.clustering.twopointcorrelation;

import algorithms.compGeometry.clustering.twopointcorrelation.RandomClusterAndBackgroundGenerator.CLUSTER_SEPARATION;
import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;
import junit.framework.TestCase;

/**
 * utility class for unit tests which need to create random background points and
 * points in clusters.  It's meant to be extended.
 * TODO: change to use composition instead of inheritance?
 *
 * @author nichole
 */
public class BaseTwoPointTest extends TestCase {

    /**
     *
     */
    protected RandomClusterAndBackgroundGenerator generator = null;

    /**
     *
     * @throws Exception
     */
    @Override
    protected void setUp() throws Exception {

        super.setUp();

        generator = new RandomClusterAndBackgroundGenerator();
    }

    /**
     *
     * @throws Exception
     */
    @Override
    protected void tearDown() throws Exception {
        super.tearDown();
    }

    /**
     *
     * @return
     */
    protected int getExpectedNumberOfClusters() {
        return generator.getExpectedNumberOfClusters();
    }

    /**
     *
     * @return
     */
    protected int getTotalNumberOfPoints() {
        return generator.getTotalNumberOfPoints();
    }

    /**
     *
     */
    public void testNeededForJunitRuntime() {}

    /**
     *
     * @return
     * @throws NoSuchAlgorithmException
     */
    protected AxisIndexer createIndexerWithRandomPoints() throws NoSuchAlgorithmException {
        return generator.createIndexerWithRandomPoints();
    }

    /**
     *
     * @param xmin
     * @param xmax
     * @param ymin
     * @param ymax
     * @return
     * @throws NoSuchAlgorithmException
     */
    protected AxisIndexer createIndexerWithRandomPoints(float xmin,
        float xmax, float ymin, float ymax) throws NoSuchAlgorithmException {

        return generator.createIndexerWithRandomPoints(xmin, xmax, ymin, ymax);
    }

    /**
     *
     * @param sr
     * @param xmin
     * @param xmax
     * @param ymin
     * @param ymax
     * @param nClusters
     * @param nBackgroundPoints
     * @param clusterSeparation
     * @return
     */
    protected AxisIndexer createIndexerWithRandomPoints(SecureRandom sr, float xmin, float xmax, float ymin, float ymax,
        int[] nClusters, int nBackgroundPoints, CLUSTER_SEPARATION clusterSeparation) {

        return generator.createIndexerWithRandomPoints(sr, xmin, xmax, ymin, ymax,
            nClusters, nBackgroundPoints, clusterSeparation);
    }

    /**
     *
     * @param sr
     * @param xmin
     * @param xmax
     * @param ymin
     * @param ymax
     * @param numberOfClusters
     * @param minimumNumberOfPointsPerCluster
     * @param maximumNumberOfPointsPerCluster
     * @param backgroundPointFractionToClusters
     * @return
     */
    protected AxisIndexer createIndexerWithRandomPoints(SecureRandom sr, float xmin, float xmax, float ymin, float ymax,
        int numberOfClusters, int minimumNumberOfPointsPerCluster, int maximumNumberOfPointsPerCluster,
        float backgroundPointFractionToClusters) {

        return generator.createIndexerWithRandomPoints(sr, xmin, xmax, ymin, ymax,
            numberOfClusters, minimumNumberOfPointsPerCluster, maximumNumberOfPointsPerCluster,
            backgroundPointFractionToClusters);
    }

    /**
     *
     * @param sr
     * @param maxRadius
     * @param numberOfPoints
     * @param xc
     * @param yc
     * @param x
     * @param y
     * @param xyStartOffset
     */
    protected void createRandomPointsAroundCenter(SecureRandom sr, float maxRadius,
        int numberOfPoints, float xc, float yc, float[] x, float[] y, int xyStartOffset) {

        generator.createRandomPointsAroundCenter(sr, maxRadius, numberOfPoints,
            xc, yc, x, y, xyStartOffset);
    }

    /**
     *
     * @param sr
     * @param xmin
     * @param xmax
     * @param ymin
     * @param ymax
     * @param numberOfPointsToCreate
     * @param xPoints
     * @param yPoints
     * @param minSeparationBetweenPoints
     */
    protected void createRandomSeparatedPoints(SecureRandom sr, float xmin, float xmax, float ymin, float ymax,
        int numberOfPointsToCreate, float[] xPoints, float[] yPoints, float minSeparationBetweenPoints) {

        generator.createRandomSeparatedPoints(sr, xmin, xmax, ymin, ymax,
            numberOfPointsToCreate, xPoints, yPoints, minSeparationBetweenPoints);
    }

    /**
     *
     * @param x
     * @param y
     * @param nXY
     * @param xp
     * @param yp
     * @param minimumSeparation
     * @return
     */
    protected boolean separationBetweenExistingPointsIsLargerThanMin(float[] x, float[] y, int nXY,
        float xp, float yp, float minimumSeparation) {

        return generator.separationBetweenExistingPointsIsLargerThanMin(x, y, nXY,
            xp, yp, minimumSeparation);
    }

    /**
     *
     * @param sr
     * @param nBackgroundPoints
     * @param xmin
     * @param xmax
     * @param ymin
     * @param ymax
     * @param x
     * @param y
     * @param xyStartOffset
     */
    protected void createRandomPointsInRectangle(SecureRandom sr, int nBackgroundPoints,
        float xmin, float xmax, float ymin, float ymax,
        float[] x, float[] y,  int xyStartOffset) {

        generator.createRandomPointsInRectangle(sr, nBackgroundPoints,
            xmin, xmax, ymin, ymax, x, y, xyStartOffset);
    }

    /**
     *
     * @param numberOfBackgroundPoints
     * @param numberOfClusterPoints
     * @param clusterSeparation
     * @param xmin
     * @param xmax
     * @param ymin
     * @param ymax
     * @param sr
     * @param useRandomForErrors
     */
    protected void createPoints(int numberOfBackgroundPoints, int[] numberOfClusterPoints,
        CLUSTER_SEPARATION clusterSeparation,
        float xmin, float xmax, float ymin, float ymax, SecureRandom sr, boolean useRandomForErrors) {

        generator.createPoints(numberOfBackgroundPoints, numberOfClusterPoints,
            clusterSeparation, xmin, xmax, ymin, ymax, sr, useRandomForErrors);
    }

     /**
     *      |
     *      |
     * -----|.....  <---- angle is w.r.t y=0, x=xc.  increases in CW order
     *      |
     *      |
     *
     * @param xc
     * @param yc
     * @param radius
     * @param angleInDegreesFromYEQ0XGT0  angle in degrees, CW from point y=0, x=xc
     * @return
     */
    public static float[] calculateXAndYFromXcYcAndRadius(float xc, float yc, float radius, double angleInDegreesFromYEQ0XGT0) {

       return RandomClusterAndBackgroundGenerator.calculateXAndYFromXcYcAndRadius(xc, yc, radius, angleInDegreesFromYEQ0XGT0);
    }

    /**
     *   *  |  *
     *      |
     * -----|.....  <---- angle is w.r.t y=0, x=xc.  increases in CCW order
     *  *   |
     *      |  *
     *
     * @param xc
     * @param yc
     * @param radius
     * @param angleInDegreesFromYEQ0XGT0  angle in degrees, CCW from point y=0, x=xc
     * @return
     */
    static float[] calculateXAndYFromXcYcAndRadiusCCW(float xc, float yc, float radius, double angleInDegreesFromYEQ0XGT0) {

        return RandomClusterAndBackgroundGenerator.calculateXAndYFromXcYcAndRadiusCCW(xc, yc, radius, angleInDegreesFromYEQ0XGT0);
    }
}
