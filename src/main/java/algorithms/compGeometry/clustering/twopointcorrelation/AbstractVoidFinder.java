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
    
    protected DoubleAxisIndexer indexer = null;
 
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
        
    public void findVoids(DoubleAxisIndexer indexer) throws TwoPointVoidStatsException {
                
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
     * this returns the intersection of points within the area defined by the x
     * range and y range where those are given in the reference frame of the
     * indexes sorted by each axis.
     *
     * @param x
     * @param y
     * @param xIndexLo index given w.r.t. indexer.sortedXIndexes
     * @param xIndexHi index given w.r.t. indexer.sortedXIndexes
     * @param yIndexLo index given w.r.t. indexer.sortedYIndexes
     * @param yIndexHi index given w.r.t. indexer.sortedYIndexes
     */
    public void processIndexedRegion(int xIndexLo, int xIndexHi, int yIndexLo, int yIndexHi, boolean useCompleteSampling) {

        // returns indexes w.r.t. indexer.x and indexer.y
        int[] regionIndexes = indexer.findIntersectingRegionIndexesIfOnlyTwo(
            xIndexLo, xIndexHi, yIndexLo, yIndexHi, useCompleteSampling);

        int nPointsInRegion = regionIndexes.length;

        // an area which has 2 points in it and has the smallest nPoints/area, that is the
        // largest area is the value we are looking for.
        // we won't be finding the furthest pair because that 'rectangular area' would presumably
        // contain other points too.
        //    *an exception is made for useCompleteSampling if more than one same y is present
        if ( (nPointsInRegion < 2) || ((nPointsInRegion != 2) && !useCompleteSampling) ) {
            return;
        }

        if (nPointsInRegion == 2) {

            processIndexedPair(regionIndexes[0], regionIndexes[1]);

        } else if (nPointsInRegion == 3) {

            // find the point that doesn't have the same y as the others.
            int uniqueYIndex = -1;
            for (int i = 0; i < regionIndexes.length; i++) {
                boolean foundSameYAsI = false;
                float ypi = indexer.getY()[ regionIndexes[i] ];
                for (int j = 0; j < regionIndexes.length; j++) {
                    if (i == j) {
                        continue;
                    }
                    float ypj = indexer.getY()[ regionIndexes[j] ];
                    if (ypj == ypi) {
                        foundSameYAsI = true;
                        break;
                    }
                }
                if (!foundSameYAsI) {
                    // we found the unique y
                    uniqueYIndex = i;
                    break;
                }
            }
            if (uniqueYIndex == -1) {
                StringBuilder err = new StringBuilder("ERROR: intersecting region contained more than 2 points, but unique y wasn't found");
                for (int i = 0; i < regionIndexes.length; i++) {
                    err.append("\n  (").append( indexer.getX()[ regionIndexes[i] ] )
                        .append(", ").append( indexer.getY()[ regionIndexes[i] ] ).append(")");
                }
                throw new IllegalStateException(err.toString());
            }
            int regionIndex0 = regionIndexes[uniqueYIndex];
            for (int i = 0; i < regionIndexes.length; i++) {
                if (i != uniqueYIndex) {
                    int regionIndex1 = regionIndexes[i];
                    processIndexedPair(regionIndex0, regionIndex1);
                }
            }
        } else if (nPointsInRegion == 4) {
            // this is a rectangle, so write all combinations
            for (int i = 0; i < regionIndexes.length; i++) {
                int regionIndex0 = regionIndexes[i];
                for (int j = (i + 1); j < regionIndexes.length; j++) {
                    int regionIndex1 = regionIndexes[j];
                    processIndexedPair(regionIndex0, regionIndex1);
                }
            }
        }
    }

    /**
     * process the pair of points by calculating the linear density and storing it
     * in the instance arrays if not already stored.
     *
     * @param regionIndex0 index w.r.t. indexer.x and indexer.y
     * @param regionIndex1 index w.r.t. indexer.x and indexer.y
     */
    public void processIndexedPair(int regionIndex0, int regionIndex1) {

        // Note: using 1-D instead of 2-D rectangles seems to be a better choice because it is not
        // dependent on rotation of the reference frame

        float d = (float) Math.sqrt(LinesAndAngles.distSquared(
            indexer.getX()[regionIndex0], indexer.getY()[regionIndex0],
            indexer.getX()[regionIndex1], indexer.getY()[regionIndex1]));

        if (d == 0) {
            return;
        }

        float linearDensity = 2.f / d;

        // expand arrays by 100 if needed
        if ((nTwoPointSurfaceDensities + 2) > allTwoPointSurfaceDensities.length) {
            allTwoPointSurfaceDensities = Arrays.copyOf(allTwoPointSurfaceDensities, nTwoPointSurfaceDensities + 100);
            point1 = Arrays.copyOf(point1, nTwoPointSurfaceDensities + 100);
            point2 = Arrays.copyOf(point2, nTwoPointSurfaceDensities + 100);
        }

        if (twoPointIdentities.storeIfDoesNotContain(regionIndex0, regionIndex1)) {
            allTwoPointSurfaceDensities[nTwoPointSurfaceDensities] = linearDensity;
            point1[nTwoPointSurfaceDensities] = regionIndex0;
            point2[nTwoPointSurfaceDensities] = regionIndex1;
            nTwoPointSurfaceDensities++;
        }
    }

    
    protected void condenseArrays() throws TwoPointVoidStatsException {
        
        if (nTwoPointSurfaceDensities == 0) {
            throw new TwoPointVoidStatsException("No pairs were found isolated within an area");
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
