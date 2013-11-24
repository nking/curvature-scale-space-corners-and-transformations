package algorithms.compGeometry.clustering.twopointcorrelation;

import algorithms.compGeometry.LinesAndAngles;
import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;
import java.util.logging.Logger;

/**
 * utility class for unit tests which need to create random background points and
 * points in clusters.
 *
 * @author nichole
 */
public class RandomClusterAndBackgroundGenerator {

    float[] x = null;
    float[] y = null;
    float[] xc = null;
    float[] yc = null;
    float[] xErrors = null;
    float[] yErrors = null;

    Logger log = Logger.getLogger(this.getClass().getName());

    static enum CLUSTER_SEPARATION {
        SMALL, MODERATE, LARGE
    }

    protected int getExpectedNumberOfClusters() {
        return (xc == null) ? 0 : xc.length;
    }
    protected int getTotalNumberOfPoints() {
        return (x == null) ? 0 : x.length;
    }

    protected DoubleAxisIndexer createIndexerWithRandomPoints() throws NoSuchAlgorithmException {

        float xmin = 0;
        float xmax = 300;
        float ymin = 0;
        float ymax = 300;

        return createIndexerWithRandomPoints(xmin, xmax, ymin, ymax);
    }

    protected DoubleAxisIndexer createIndexerWithRandomPoints(float xmin,
        float xmax, float ymin, float ymax) throws NoSuchAlgorithmException {

        SecureRandom srr = SecureRandom.getInstance("SHA1PRNG");
        srr.setSeed( System.currentTimeMillis() );
        long seed = srr.nextLong();

        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        sr.setSeed( seed );
        //sr.setSeed(-2384802679227907254l);

        // randomly choose between point sets:
        //  (0) sparsely populated background with clusters
        //  (1) moderately populated background with clusters
        //  (2) densely populated background with clusters
        //  (3) only background points

        // the limit to the number of points here is comparable to the number
        // within the tests in TwoPointVoidTests, and those are kept somewhat
        // small due to runtime with 'useCompleteSampling'

        int setChoice = sr.nextInt(4);

        int nBackgroundPoints = 0;
        int[] nClusters = null;

        CLUSTER_SEPARATION clusterSep = null;

        int sum = 0;
        switch(setChoice) {
            case 0:
                nClusters = new int[]{30, 40, 60};
                sum = 0;
                for (int i = 0; i < nClusters.length; i++) {
                    sum += nClusters[i];
                }
                nBackgroundPoints = (int)0.1*sum;
                clusterSep = CLUSTER_SEPARATION.values()[sr.nextInt(2)];
                break;
            case 1:
                nClusters = new int[]{30, 40, 60};
                sum = 0;
                for (int i = 0; i < nClusters.length; i++) {
                    sum += nClusters[i];
                }
                nBackgroundPoints = sum;
                clusterSep = CLUSTER_SEPARATION.values()[sr.nextInt(2)];
                break;
            case 2:
                nClusters = new int[]{30, 40, 60};
                sum = 0;
                for (int i = 0; i < nClusters.length; i++) {
                    sum += nClusters[i];
                }
                nBackgroundPoints = 10*sum;
                clusterSep = CLUSTER_SEPARATION.values()[sr.nextInt(2)];
                break;
            default:
            //case 3:
                nClusters = new int[0];
                nBackgroundPoints = 150;
                break;
        }

        int nn = (nClusters == null) ? 0 : nClusters.length;
        String ns = (clusterSep == null) ? "" : clusterSep.name();

        log.info("Creating points: " + nn + " clusters, "
            + nBackgroundPoints + " background points, " + ns);

        return createIndexerWithRandomPoints(sr, xmin, xmax, ymin, ymax, nClusters, nBackgroundPoints, clusterSep);
    }

    protected DoubleAxisIndexer createIndexerWithRandomPoints(SecureRandom sr, float xmin, float xmax, float ymin, float ymax,
        int[] nClusters, int nBackgroundPoints, CLUSTER_SEPARATION clusterSeparation) {

        createPoints(nBackgroundPoints, nClusters, clusterSeparation, xmin, xmax, ymin, ymax, sr, false);

        DoubleAxisIndexer indexer = new DoubleAxisIndexer();
        indexer.sortAndIndexXThenY(x, y, xErrors, yErrors, x.length);

        return indexer;
    }

    protected DoubleAxisIndexer createIndexerWithRandomPoints(SecureRandom sr, float xmin, float xmax, float ymin, float ymax,
        int numberOfClusters, int minimumNumberOfPointsPerCluster, int maximumNumberOfPointsPerCluster,
        float backgroundPointFractionToClusters) {

        int dn = maximumNumberOfPointsPerCluster - minimumNumberOfPointsPerCluster;

        int count = 0;
        int[] nClusters = new int[numberOfClusters];
        for (int i = 0; i < numberOfClusters; i++) {
            int add = (dn == 0) ? 0 : sr.nextInt(dn);
            nClusters[i] = + minimumNumberOfPointsPerCluster + add;
            count += nClusters[i];
        }

        int nBackgroundPoints = (int)backgroundPointFractionToClusters*count;

        CLUSTER_SEPARATION clusterSeparation = CLUSTER_SEPARATION.values()[sr.nextInt(2)];


        createPoints(nBackgroundPoints, nClusters, clusterSeparation, xmin, xmax, ymin, ymax, sr, false);

        DoubleAxisIndexer indexer = new DoubleAxisIndexer();
        indexer.sortAndIndexXThenY(x, y, xErrors, yErrors, x.length);

        return indexer;
    }

    protected DoubleAxisIndexer createIndexerWithRandomPoints(SecureRandom sr, float xmin, float xmax, float ymin, float ymax,
        int numberOfClusters, int minimumNumberOfPointsPerCluster, int maximumNumberOfPointsPerCluster,
        float backgroundPointFractionToClusters, CLUSTER_SEPARATION clusterSeparation) {

        int dn = maximumNumberOfPointsPerCluster - minimumNumberOfPointsPerCluster;

        int count = 0;
        int[] nClusters = new int[numberOfClusters];
        for (int i = 0; i < numberOfClusters; i++) {
            int add = (dn == 0) ? 0 : sr.nextInt(dn);
            nClusters[i] = + minimumNumberOfPointsPerCluster + add;
            count += nClusters[i];
        }

        int nBackgroundPoints = (int)backgroundPointFractionToClusters*count;

        createPoints(nBackgroundPoints, nClusters, clusterSeparation, xmin, xmax, ymin, ymax, sr, false);

        DoubleAxisIndexer indexer = new DoubleAxisIndexer();
        indexer.sortAndIndexXThenY(x, y, xErrors, yErrors, x.length);

        return indexer;
    }

    protected void createRandomPointsAroundCenter(SecureRandom sr, float maxRadius,
        int numberOfPoints, float xc0, float yc0, float[] x0, float[] y0, int xyStartOffset) {

        for (int i = 0; i < numberOfPoints; i++) {

            float radius = maxRadius * sr.nextFloat();
            double angle = 360. * sr.nextDouble();

            float[] xy = calculateXAndYFromXcYcAndRadius(xc0, yc0, radius, angle);

            if ((xy[0] > 0) && (xy[1] > 0)) {
                x0[xyStartOffset + i] = xy[0];
                y0[xyStartOffset + i] = xy[1];
            } else {
                i--;
            }
        }
    }

    protected void createRandomSeparatedPoints(SecureRandom sr, float xmin, float xmax, float ymin, float ymax,
        int numberOfPointsToCreate, float[] xPoints, float[] yPoints, float minSeparationBetweenPoints) {

        float xWidth = xmax - xmin;
        float yHeight = ymax - ymin;

        float xx = -1;
        float yy = -1;

        for (int i = 0; i < numberOfPointsToCreate; i++) {

            boolean sepIsLarger = false;

            while (!sepIsLarger) {
                xx = xmin + sr.nextFloat()*xWidth;
                yy = ymin + sr.nextFloat()*yHeight;

                sepIsLarger = separationBetweenExistingPointsIsLargerThanMin(
                    xPoints, yPoints, i, xx, yy, minSeparationBetweenPoints);

            }
            xPoints[i] = xx;
            yPoints[i] = yy;
        }
    }

    protected boolean separationBetweenExistingPointsIsLargerThanMin(float[] x0, float[] y0, int nXY,
        float xp, float yp, float minimumSeparation) {

        if (nXY == 0) {
            return true;
        }

        float eps = minimumSeparation/10.f;

        float minSq = minimumSeparation * minimumSeparation;

        for (int i = 0; i < nXY; i++) {

            double distSq = LinesAndAngles.distSquared(x0[i], y0[i], xp, yp) + eps;

            if (distSq < minSq) {
                return false;
            }
        }
        return true;
    }

    protected void createRandomPointsInRectangle(SecureRandom sr, int nBackgroundPoints,
        float xmin, float xmax, float ymin, float ymax,
        float[] x0, float[] y0,  int xyStartOffset) {

        float xWidth = xmax - xmin;
        float yHeight = ymax - ymin;

        for (int i = 0; i < nBackgroundPoints; i++) {
            x0[xyStartOffset + i] = sr.nextFloat()*xWidth;
            y0[xyStartOffset + i] = sr.nextFloat()*yHeight;
        }
    }

    protected void createPoints(int numberOfBackgroundPoints, int[] numberOfClusterPoints,
        CLUSTER_SEPARATION clusterSeparation,
        float xmin, float xmax, float ymin, float ymax, SecureRandom sr, boolean useRandomForErrors) {

        int nClusters = (numberOfClusterPoints == null) ? 0 : numberOfClusterPoints.length;

        int nTotalPoints = 0;

        for (int i = 0; i < nClusters; i++) {
            nTotalPoints += numberOfClusterPoints[i];
        }

        nTotalPoints += numberOfBackgroundPoints;

        // contains points sequentially for: background, cluster[0],...cluster[n-1]
        x = new float[nTotalPoints];
        y = new float[nTotalPoints];

        // compare these to findings.  they are used to generate the clusters
        xc = new float[nClusters];
        yc = new float[nClusters];

        createRandomPointsInRectangle(sr, numberOfBackgroundPoints, xmin, xmax, ymin, ymax, x, y, 0);

        /*
         *  For n = 3
         *
         *  |----max/n----|             |             |
         *  |             |             |             |
         *         *             *             *
         *      |  .  |       |  .  |       |  .  |
         *            |       |                .  |
         *             d=max/2*n               .  |
         *                                     r=max/4*n
         */
        if (nClusters > 0) {

            float maxClusterRadius = (xmax - xmin) / (4.0f * nClusters);
            float minDistanceBetweenClusterCenters;
            float factor;
            if (clusterSeparation.ordinal() == CLUSTER_SEPARATION.LARGE.ordinal()) {
                factor = 4.0f;
            } else if (clusterSeparation.ordinal() == CLUSTER_SEPARATION.SMALL.ordinal()) {
                factor = 2.0f;
            } else {
                factor = 3.0f;
            }
            maxClusterRadius = 0.8f*(xmax - xmin) / (factor * nClusters);
            minDistanceBetweenClusterCenters = (factor * maxClusterRadius);


            createRandomSeparatedPoints(sr, xmin + maxClusterRadius, xmax - maxClusterRadius,
                ymin + maxClusterRadius, ymax - maxClusterRadius,
                nClusters, xc, yc, minDistanceBetweenClusterCenters);

            int startOffset = numberOfBackgroundPoints;

            for (int i = 0; i < nClusters; i++) {

                int n = numberOfClusterPoints[i];

                createRandomPointsAroundCenter(sr, maxClusterRadius, n, xc[i], yc[i], x, y, startOffset);

                startOffset += n;
            }
        }

        xErrors = new float[nTotalPoints];
        yErrors = new float[nTotalPoints];
        for (int i = 0; i < nTotalPoints; i++) {
            // simulate x error as a percent error of 0.03 for each bin
            xErrors[i] = x[i] * 0.03f;
            yErrors[i] = (float) (Math.sqrt(y[i]));
        }
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

        double dx = radius * Math.cos(angleInDegreesFromYEQ0XGT0 * (Math.PI/180.f));
        double dy = radius * Math.sin(angleInDegreesFromYEQ0XGT0 * (Math.PI/180.f));

        float x = (float) (xc + dx);
        float y = (float) (yc - dy);
        return new float[]{x, y};
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

        return calculateXAndYFromXcYcAndRadius(xc, yc, radius, 360 - angleInDegreesFromYEQ0XGT0);
    }
}
