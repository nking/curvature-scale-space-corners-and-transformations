package algorithms.compGeometry.clustering.twopointcorrelation;

import java.util.Arrays;
import java.util.logging.Logger;

import algorithms.compGeometry.LinesAndAngles;

public abstract class AbstractVoidFinder implements IVoidFinder {

    protected float[] allTwoPointSurfaceDensities = null;

    protected float[] allTwoPointSurfaceDensitiesErrors = null;

    protected int nTwoPointSurfaceDensities = 0;

    protected int[] point1 = null;

    protected int[] point2 = null;

    protected ITwoPointIdentity twoPointIdentities = null;

    protected VoidSampling sampling = null;

    protected Logger log = null;

    protected boolean debug = false;

    protected AxisIndexer indexer = null;

    public AbstractVoidFinder() {
    }

    @Override
    public float[] getTwoPointDensities() {
        return allTwoPointSurfaceDensities;
    }

    @Override
    public float[] getTwoPointDensityErrors() {
        return allTwoPointSurfaceDensitiesErrors;
    }

    @Override
    public int getNumberOfTwoPointDensities() {
        return nTwoPointSurfaceDensities;
    }

    public int[] getPoint1() {
        return point1;
    }

    public int[] getPoint2() {
        return point2;
    }

    public void setSampling(VoidSampling sampling) {
        this.sampling = sampling;
    }
    public VoidSampling getSampling() {
        return this.sampling;
    }

    public void setDebug(boolean setDebugToTrue) {
        this.debug = setDebugToTrue;
    }

    protected void initializeVariables() {

        allTwoPointSurfaceDensities = new float[100];
        point1 = new int[100];
        point2 = new int[100];
        twoPointIdentities = TwoPointIdentityFactory.create(this.indexer.getNXY());
    }

    protected abstract void findVoidsImpl();

    public abstract void constructLogger();

    public void findVoids(AxisIndexer indexer) throws TwoPointVoidStatsException {

        if (indexer == null) {
            throw new IllegalArgumentException("indexer cannot be null");
        }

        this.indexer = indexer;

        constructLogger();

        initializeVariables();

        findVoidsImpl();

        condenseArrays();
    }

    /**
     * process the region bounded by the pair if the region is found to be a void
     * bounded only by the 2 given points referenced by the indexes.
     *
     * @param x
     * @param y
     * @param xSortedIndex0 index within indexer.sortedXIndexes which is an index to reference a point in indexer's x and y
     * @param xSortedIndex1 index within indexer.sortedXIndexes which is an index to reference a point in indexer's x and y
     */
    public void processIndexedRegion(int xSortedIndex0, int xSortedIndex1) {

        if (xSortedIndex0 > xSortedIndex1) {
            int tmp = xSortedIndex0;
            xSortedIndex0 = xSortedIndex1;
            xSortedIndex1 = tmp;
        }

        int[] sortedXIndexes = indexer.getSortedXIndexes();

        float[] y = indexer.getY();

        int idx1 = sortedXIndexes[xSortedIndex0];
        int idx2 = sortedXIndexes[xSortedIndex1];

        // x's are already ordered by increasing x
        float y0 = y[idx1];
        float y1 = y[idx2];

        if (y0 > y1) {
            float tmp = y0;
            y0 = y1;
            y1 = tmp;
        }

        float x0 = indexer.getX()[idx1];
        float x1 = indexer.getX()[idx2];

        boolean doProcess = true;

        // O(1) to O(N-2) at worse
        for (int i = (xSortedIndex0 + 1); i < xSortedIndex1; i++) {

            int idx = sortedXIndexes[i];

            float xt = indexer.getX()[idx];

            float yt = y[idx];

            // quickly let points outside of boundaries pass
            if (xt < x0) {
                continue;
            } else if (xt > x1) {
                continue;
            } else if (yt < y0) {
                continue;
            } else if (yt > y1) {
                continue;
            }

            if ( ((xt > x0) && (xt < x1)) || ((yt > y0) && (yt < y1)) ) {

                // else we're in bounds of a rectangle or line whose corners are the 2 points x0,y0  x1,y1
                doProcess = false;

                break;
            }

            log.finest("(" + x0 + "," + y0 +") (" + x1 + "," + y1 + ")  test " + xt + "," + yt + "  passed");
        }

        if (doProcess) {

            // if it has passed, we still need to look for the same x values
            //   for sortedXIndexes just before xSortedIndex0 and just after xSortedIndex1

            int t2Idx = xSortedIndex0 - 1;

            while (t2Idx > -1 && ( indexer.getX()[sortedXIndexes[t2Idx]] == x0 )) {
                // check whether it's within boundaries
                float y2t = y[t2Idx];
                if ((y2t > y0) && (y2t < y1)) {
                    doProcess = false;
                    break;
                }
                t2Idx--;
            }

            if (doProcess) {

                t2Idx = xSortedIndex1 + 1;

                while ((t2Idx < (indexer.getNumberOfPoints() - 1)) && ( indexer.getX()[sortedXIndexes[t2Idx]] == x1 )) {
                    float y2t = y[t2Idx];
                    if ((y2t > y0) && (y2t < y1)) {
                        doProcess = false;
                        break;
                    }
                    t2Idx++;
                }

                if (doProcess) {
                    processIndexedPair(idx1, idx2);
                }
            }
        }
    }

    /**
     * process the pair of points by calculating the linear density and storing it
     * in the instance arrays if not already stored.
     *
     * @param idx0 index w.r.t. indexer.x and indexer.y
     * @param idx1 index w.r.t. indexer.x and indexer.y
     */
    public void processIndexedPair(int idx0, int idx1) {

        // Note: using 1-D instead of 2-D rectangles seems to be a better choice because it is not
        // dependent on rotation of the reference frame

        float d = (float) Math.sqrt(LinesAndAngles.distSquared(
            indexer.getX()[idx0], indexer.getY()[idx0],
            indexer.getX()[idx1], indexer.getY()[idx1]));

        if (d == 0) {
            return;
        }

        float linearDensity = (float) (2.f / Math.sqrt(d));

        log.finest("(" + indexer.getX()[idx0] + "," + indexer.getY()[idx0] + ") ("
            + indexer.getX()[idx1] + "," + indexer.getY()[idx1] + ")  ld=" + linearDensity);

        // expand arrays by 100 if needed
        if ((nTwoPointSurfaceDensities + 2) > allTwoPointSurfaceDensities.length) {
            allTwoPointSurfaceDensities = Arrays.copyOf(allTwoPointSurfaceDensities, nTwoPointSurfaceDensities + 100);
            point1 = Arrays.copyOf(point1, nTwoPointSurfaceDensities + 100);
            point2 = Arrays.copyOf(point2, nTwoPointSurfaceDensities + 100);
        }

        if (twoPointIdentities.storeIfDoesNotContain(idx0, idx1)) {
            allTwoPointSurfaceDensities[nTwoPointSurfaceDensities] = linearDensity;
            point1[nTwoPointSurfaceDensities] = idx0;
            point2[nTwoPointSurfaceDensities] = idx1;
            nTwoPointSurfaceDensities++;
        }
    }


    protected void condenseArrays() throws TwoPointVoidStatsException {

        if (nTwoPointSurfaceDensities == 0) {
            //throw new TwoPointVoidStatsException("No pairs were found isolated within an area");
        }

        // condense arrays
        allTwoPointSurfaceDensities = Arrays.copyOf(allTwoPointSurfaceDensities, nTwoPointSurfaceDensities);
        point1 = Arrays.copyOf(point1, nTwoPointSurfaceDensities);
        point2 = Arrays.copyOf(point2, nTwoPointSurfaceDensities);

        allTwoPointSurfaceDensitiesErrors =
            calulateTwoPointDensityErrors(allTwoPointSurfaceDensities, point1, point2,
            indexer.getX(), indexer.getY(), indexer.getXErrors(), indexer.getYErrors());

        // release twoPointIdentities to free up memory
        /*twoPointIdentities = null;
        MemoryMXBean mb = ManagementFactory.getMemoryMXBean();
        mb.gc();*/

    }

    /**
     * Calculate the 2 point density errors following the chain rule
     *
     *                                | df |^2               | df |^2         df   df
     *      (sigma_f)^2 =  (sigma_x)^2|----|   +  (sigma_y)^2|----|    +  2 * -- * -- * cov_ab
     *                                | dx |                 | dy |           dx   dy
     *
     *      For uncorrelated variables the covariance terms are zero.
     *
     * For two-point density:
     *                 N             dF      -N
     *      F(x,y) = ------   ==>   ---- =  -----
     *                 X             dx      X^2
     *
     *                                  |   -N   |^2
     *     (sigma_f)^2 =  (sigma_x)^2 * |--------|
     *                                  |   X^2  |
     *
     * @param densities - the two point densities
     * @param point1Indexes indexes to xp, yp of one of the 2 points in the two-point densities
     * @param point2Indexes indexes to xp, yp of the other of the 2 points in the two-point densities
     * @param xp x array of points
     * @param yp y array of points
     * @param xpe errors in x
     * @param ype errors in y
     * @return errors for densities
     */
    public float[] calulateTwoPointDensityErrors(float[] densities, int[] point1Indexes, int[] point2Indexes,
        float[] xp, float[] yp, float[] xpe, float[] ype) {

        if (densities == null) {
            throw new IllegalArgumentException("densities cannot be null");
        }
        if (point1Indexes == null || point2Indexes == null) {
            throw new IllegalArgumentException("pointIndexes cannot be null");
        }
        if (xp == null || yp == null) {
            throw new IllegalArgumentException("xp and yp cannot be null");
        }
        if (xpe == null || ype == null) {
            throw new IllegalArgumentException("errors are required");
        }

        float[] densityErrors = new float[densities.length];

        for (int i = 0; i < densities.length; i++) {

            if (Float.isInfinite(densities[i])) {
                continue;
            }

            int pt1 = point1[i];
            int pt2 = point2[i];

            double xe1 = xpe[pt1] ; // xp[pt1];
            double xe2 = xpe[pt2] ; // xp[pt2];
            double ye1 = ype[pt1] ; // yp[pt1];
            double ye2 = ype[pt2] ; // yp[pt2];

            double xsigmasq = ( xe1 * xe1 ) + ( xe2 * xe2 );
            double ysigmasq = ( ye1 * ye1 ) + ( ye2 * ye2 );
            double xysigmasq = xsigmasq + ysigmasq;

            double linearDensity = allTwoPointSurfaceDensities[i];
            double linearDensitySq = Math.pow(linearDensity, 2);

            double err = Math.sqrt(xysigmasq * linearDensitySq);

            densityErrors[i] = (float)err;
        }
        return densityErrors;
    }


    /*
     (non-Javadoc)
     * @see algorithms.compGeometry.clustering.twopointcorrelation.IVoidFinder#approximateMemoryUsed()
     */
    @Override
    public long approximateMemoryUsed() {

        String arch = System.getProperty("sun.arch.data.model");

        boolean is32Bit = ((arch != null) && arch.equals("64")) ? false : true;

        int nbits = (is32Bit) ? 32 : 64;

        int overheadBytes = 16;

        int intBytes = (is32Bit) ? 4 : 8;
        int arrayBytes = 32/8;
        int refBytes = nbits/8;

        long sumBytes = 0;

        // 4 array references
        sumBytes += (4*arrayBytes);

        // 2 word size primitives
        sumBytes += (2*intBytes);

        // 1 object reference
        sumBytes += refBytes;


        int n = (allTwoPointSurfaceDensities == null) ? 0 : allTwoPointSurfaceDensities.length;

        if (twoPointIdentities != null) {
            sumBytes += twoPointIdentities.approximateMemoryUsed();
        }

        if (allTwoPointSurfaceDensities != null) {
            sumBytes +=  (allTwoPointSurfaceDensities.length*intBytes);
        }

        if (allTwoPointSurfaceDensitiesErrors != null) {
            sumBytes += (allTwoPointSurfaceDensities.length*intBytes);
        }

        if (point1 != null) {
            sumBytes += (point1.length*intBytes);
        }

        if (point2 != null) {
            sumBytes += (point2.length*intBytes);
        }

        sumBytes += overheadBytes;

        // amount of padding needed to make it a round 8 bytes
        long padding = (sumBytes % 8);

        sumBytes += padding;

        return sumBytes;
    }

}
