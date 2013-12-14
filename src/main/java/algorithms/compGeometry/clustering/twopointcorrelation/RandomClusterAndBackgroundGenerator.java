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

    public DoubleAxisIndexer createIndexerWithRandomPoints() throws NoSuchAlgorithmException {

        float xmin = 0;
        float xmax = 300;
        float ymin = 0;
        float ymax = 300;

        return createIndexerWithRandomPoints(xmin, xmax, ymin, ymax);
    }

    public DoubleAxisIndexer createIndexerWithRandomPoints(float xmin,
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

    public DoubleAxisIndexer createIndexerWithRandomPoints(SecureRandom sr, float xmin, float xmax, float ymin, float ymax,
        int[] nClusters, int nBackgroundPoints, CLUSTER_SEPARATION clusterSeparation) {

        createPoints(nBackgroundPoints, nClusters, clusterSeparation, xmin, xmax, ymin, ymax, sr, false);

        DoubleAxisIndexer indexer = new DoubleAxisIndexer();
        indexer.sortAndIndexXThenY(x, y, xErrors, yErrors, x.length);

        return indexer;
    }

    public DoubleAxisIndexer createIndexerWithRandomPoints(SecureRandom sr, float xmin, float xmax, float ymin, float ymax,
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

    public DoubleAxisIndexer createIndexerWithRandomPoints(SecureRandom sr, float xmin, float xmax, float ymin, float ymax,
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
    
    /**
     * create a group of points around a center whose density distribution decreases
     * by distance squared from the center.  the points are randomly chosen to populate
     * that distribution.
     * 
     * @param sr - instance of secure random to use
     * @param maxRadius - maximum radius of group of points
     * @param numberOfPoints - number of points in the group
     * @param xc0 - x coordinate of the group center to create
     * @param yc0 - y coordinate of the group center to create
     * @param x0 - the output array to hold results for the x coordinates of the group
     * @param y0 - the output array to hold results for the y coordinates of the group
     * @param xyStartOffset the offset in the arrays x0 and y0 indexes with which to start adding points.
     */
    protected void createRandomPointsAroundCenterWithDSquared(SecureRandom sr, float maxRadius,
        int numberOfPoints, float xc0, float yc0, float[] x0, float[] y0, int xyStartOffset) {
        
        /* want to represent the increasing density towards the center of the group 
         * using annular radii that result in equal areas (though smaller annular widths).
         *     
         *     The area of an annulus:
         *         area = pi * (r0^2 - r1^2) where r0 is outer annuli and r1 is inner radius of annulus
         *         
         *         try 2 annuli to start relationship:
         *         
         *         pi*(r0^2 - r1^2) = pi*(r1^2 - 0)
         *            (r0^2 - r1^2) = r1^2
         *            r0^2 = 2 * (r1^2)
         *            r0 = sqrt 2 * r1 for them to have equal areas
         *       
         *         then 4 annuli:
         *         |  |    |      |         |
         *         r0 r1   r2     r3        c
         *         (r0^2 - r1^2) = (r1^2 - r2^2) = (r2^2 - r3^2) = r3^2
         *            --> r2 = (sqrt 2) * r3
         *            --> r1^2 =  2*(r2^2) - r3^3 = 4*r3^2 - r3^2 = 3r3^2
         *            
         *            (r0^2 - r1^2) = (r1^2 - r2^2)
         *            r0^2 = (2*r1^2 - r2^2)
         *                 = 2*( 3r3^2 ) - 2*(r3^2)
         *                 = 6*r3^2 - 2*r3^2
         *                 = 4*r3^2
         *
         *         in summary, for 4 annuli,
         *             r0^2 = 4 * r3^2  = maxr^2
         *             r1^2 = 3 * r3^2
         *             r2^2 = 2 * r3^2
         *           
         *         extrapolate to nbins from i=0 at outer edge to i=nbins-1 at center:
         *             (r_i)^2 = (nbins - i) * (r3)^2
         *             
         *             and r3^2 = maxr^2/4
         *             
         *             (r_i)^2 = (nbins - i) * (maxr)^2 / nbins
         *             
         *             ==> r_i = math.sqrt(  (nbins - i)/nbins ) * maxr
         *      
         *     For those equal area annuli, we need n_i to increase with i.
         *     Knowing that we held the area constant over i,
         *     we have n_i = n * fraction / (r_i)^2
         *     
         *     for 4 annuli:
         *         |  |    |      |         |
         *         r0 r1   r2     r3        c
         *         
         *                                 1               1               1            1
         *         n =  n*fraction *  ------------  + ------------ + ------------ + -----------
         *                              (r_0)^2          (r_1)^2        (r_2)^2       (r_3)^2
         *         
         *               1              1               1               1            1
         *         ----------- =   ------------  + ------------ + ------------ + -----------
         *           fraction        (r_0)^2          (r_1)^2        (r_2)^2       (r_3)^2
         *           
         *  Then to create a density distribution that is radially increasing by r^2 towards center of group,
         *     choose the number of bins.
         *     solve the fraction.
         *     create n_i points within that  r_i annulus randomly.
         */ 
                
        int nBins = 5;
        
        double sum = 0;
        for (int i = 0; i < nBins; i++) {
            
            float x = (nBins - i)/(float)nBins;
            float xInner = (nBins - i - 1.f)/(float)nBins;
            
            // r_i = math.sqrt(  (nbins - i)/nbins ) * maxr
            // this is rOuter for current bin
            double rOuter = maxRadius * Math.sqrt(x);
            double rInner = maxRadius * Math.sqrt(xInner);
            double ri = (rOuter + rInner)/2.;
            
            double f = Math.pow(maxRadius/(maxRadius - ri), 2);
                                                
            sum += (numberOfPoints/f);
        }
        float fraction = (float)(numberOfPoints/sum);
        
        int offset = xyStartOffset;
        
        for (int i = 0; i < nBins; i++) {
                        
            // index 0 is outer edge, rOuter = maxRadius
            
            float x = (nBins - i)/(float)nBins;
            float xInner = (nBins - i - 1.f)/(float)nBins;
            
            // r_i = math.sqrt(  (nbins - i)/nbins ) * maxr
            // this is rOuter for current bin
            double rOuter = maxRadius * Math.sqrt(x);
            double rInner = maxRadius * Math.sqrt(xInner);
            double ri = (rOuter + rInner)/2.;
            
            // maxradius = 80,  n=300
            // i=0   r_i=76   n ~ 0            f=large number 1/(maxRadius-r_i)
            // i=x   r_x=40   n ~ 300/(2^2)     
            // i=n-1 r_n-1=0  n close to 300
            
            double f = Math.pow(maxRadius/(maxRadius - ri), 2);
            
            double n = fraction * numberOfPoints/f;
            
            System.out.println("n_" + i + " = " + n + " r_i=" + ri + " f_i=" + f
                + " (maxRadius=" + maxRadius 
                + " rOuter=" + rOuter + " rInner=" + rInner + ")");
            
            int np = (int)n;
            //if (i == 2) {
            for (int ii = 0; ii < np; ii++) {
                float[] xy = calculateRandomXAndYWithinAnnulus(sr, xc0, yc0, (float)rInner, (float)rOuter);
                if ((xy[0] > 0) && (xy[1] > 0)) {
                    x0[offset + ii] = xy[0];
                    y0[offset + ii] = xy[1];
                } else {
                    // redo random point
                    ii--;
                }
            }
            //}
            offset += np;
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

    public DoubleAxisIndexer createIndexerWithRandomPointsAroundCenterWithDSquared(
        SecureRandom sr, int numberOfClusterPoints,
        float xmin, float xmax, float ymin, float ymax, float maximumRadius) {

        // contains points sequentially for: background, cluster[0],...cluster[n-1]
        x = new float[numberOfClusterPoints];
        y = new float[numberOfClusterPoints];

        // compare these to findings.  they are used to generate the clusters
        xc = new float[1];
        yc = new float[1];
        
        /*
         *  draw center of group as anywhere from xmin + maxRadius to xmax-maxRadius and similar for y
         */

        float xcd = (xmax + xmin)/2.f;
        float ycd = (ymax + ymin)/2.f;
        float xdiff = ((xmax - xmin)/2.f) - maximumRadius;
        float ydiff = ((ymax - ymin)/2.f) - maximumRadius;
        float xd = (xdiff*sr.nextFloat());
        float yd = (ydiff*sr.nextFloat());
        float xCenter = (sr.nextBoolean()) ? xcd + xd : xcd - xd;
        float yCenter = (sr.nextBoolean()) ? ycd + yd : ycd - yd;
        
        xc[0] = xCenter;
        yc[0] = yCenter;
        
        createRandomPointsAroundCenterWithDSquared(sr, maximumRadius,
            numberOfClusterPoints, xCenter, yCenter, x, y, 0);

        xErrors = new float[numberOfClusterPoints];
        yErrors = new float[numberOfClusterPoints];
        for (int i = 0; i < numberOfClusterPoints; i++) {
            // simulate x error as a percent error of 0.03 for each bin
            xErrors[i] = x[i] * 0.03f;
            yErrors[i] = (float) (Math.sqrt(y[i]));
        }
        
        DoubleAxisIndexer indexer = new DoubleAxisIndexer();
        indexer.sortAndIndexXThenY(x, y, xErrors, yErrors, x.length);

        return indexer;
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
    
    static float[] calculateRandomXAndYWithinAnnulus(SecureRandom sr, 
        float xc, float yc, float innerRadius, float outerRadius) {

        /* 
        *    -- choose a random number within annular width, then add inner radius to it.
        *    -- then choose a random angle zero through 360.
        *    -- then determine x and y from the distance and the angle
        */
        float radius = innerRadius + (outerRadius - innerRadius) * sr.nextFloat();
        
        double angle = 360. * sr.nextDouble();
        
        float[] xy = calculateXAndYFromXcYcAndRadius(xc, yc, radius, angle);
        
        return xy;
    }

}
