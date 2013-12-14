package algorithms.compGeometry.clustering.twopointcorrelation;

import algorithms.compGeometry.LinesAndAngles;
import algorithms.compGeometry.convexHull.GrahamScanTooFewPointsException;
import algorithms.compGeometry.convexHull.GrahamScan;
import algorithms.curves.GEVYFit;
import algorithms.misc.MiscMath;
import algorithms.util.Errors;
import algorithms.util.PolygonAndPointPlotter;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.lang.management.ManagementFactory;
import java.lang.management.MemoryMXBean;
import java.lang.management.MemoryUsage;
import java.util.logging.Logger;

/**
  Find clusters in data by looking for regions whose density is
      2 or 3 times the background density (that is 2 or 3 sigma above 'error').
      The default is 3.

  The background density can be determined with the following methods:
      FIT_TWO_POINT_VOIDS -
         clusterFinder.calculateBackgroundVia2PtVoidFit(false);
         returns an estimate of the background density
         by calculating the density of rectangles holding only 2 points,
         fitting a GEV curve to that distribution, and returning the point
         at which 10% of the total area under the curve has occurred.
      OR the background density can be set manually:
          setBackground(float backgroundSurfaceDensity, float standardDeviationOfBackground);

  Use as an API:
      TwoPointCorrelation clusterFinder = new TwoPointCorrelation(x, y, xErrors, yErrors, getTotalNumberOfPoints());
      clusterFinder.calculateBackground();
      clusterFinder.findClusters();
      clusterFinder.calculateHullsOfClusters();
      String plotFilePath = clusterFinder.plotClusters();


      TwoPointCorrelation clusterFinder = new TwoPointCorrelation(x, y, xErrors, yErrors, getTotalNumberOfPoints());
      setBackground(0.03f, 0.003f);
      clusterFinder.findClusters();
      clusterFinder.calculateHullsOfClusters();
      String plotFilePath = clusterFinder.plotClusters();

  Note:  For datasets in which the density of background points is high, if you don't
     have the ability to reduce the data by a key characteristic, you might consider
     the results of this code as seeds for a Voronoi diagram or other code.
     float[] xSeeds = clusterFinder.getXHullCentroids();
     float[] ySeeds = clusterFinder.getYHullCentroids();

  Use from the command line:
      Requires a tab delimited text file with 4 columns: x, y, xErrors, yErrors.

          java -cp bin/classes  algorithms.compGeometry.clustering.twopointcorrelation.TwoPointCorrelation --file /path/to/file/fileName.txt

  @author nichole
 */
public class TwoPointCorrelation {

    protected enum STATE {
        INITIALIZED, BACKGROUND_SET, CLUSTERS_FOUND, CLUSTER_HULLS_CALCULATED
    }

    protected enum BACKGROUND_METHOD {
        FIT_TWO_POINT_VOIDS, USER_SUPPLIED, DESERIALIZED
    }

    protected final DoubleAxisIndexer indexer;
    
    protected Boolean refineSolution = Boolean.FALSE;
    
    protected boolean allowRefinement = false;
    
    protected DoubleAxisIndexer tempRefineSolnIndexer = null;

    private float backgroundSurfaceDensity;
    private float backgroundError;
    private float sigmaFactor = 2.5f;
    // we are looking for points which have density > sigmaFactor*backgroundAverage

    protected int minimumNumberInCluster = 10;

    protected STATE state = null;
    protected BACKGROUND_METHOD bMethod = null;

    protected IGroupFinder groupFinder = null;

    /**
     * for groups already calculated in groupMembership, this is how to find the points
     * of the convex hull for each group.
     * For example, to find the convex hull of group # 0
     *    groupHullIndexes[0] will return a linked list to the indexes of x and y
     *    that form the hull for that group.
     */
    protected SimpleLinkedListNode[] groupHullIndexes = null;

    /**
     * centroid coordinates of the hulls for the groups.  note that the centroids
     * are not derived from all points in the group, only from the hull polygon.
     */
    protected float[] xGroupHullCentroids = null;
    protected float[] yGroupHullCentroids = null;
    protected float[] groupHullSurfaceAreas = null;

    protected boolean persistTheMinimaStats = false;
    protected String indexerFilePath = null;
    protected String minimaStatsFilePath = null;

    // for debugging plots, keeping a handle on TwoPointVoidStats.
    public IPointBackgroundStats backgroundStats = null;

    protected boolean debug = false;

    protected boolean doLogPerformanceMetrics = false;
    
    protected boolean setUseDownhillSimplexHistogramFitting = false;

    protected Logger log = Logger.getLogger(this.getClass().getName());

    /**
     * constructor without errors on xPoints and yPoints.  Note that the
     * errors are estimated internally as rms, shot noise and used throughout
     * the code.
     *
     * @param xPoints
     * @param yPoints
     * @param nXYPoints
     */
    public TwoPointCorrelation(float[] xPoints, float[] yPoints, int nXYPoints) {

        this.indexer = new DoubleAxisIndexer();

        float[] xPointErrors = Errors.populateYErrorsBySqrt(xPoints);
        float[] yPointErrors = Errors.populateYErrorsBySqrt(yPoints);

        indexer.sortAndIndexXThenY(xPoints, yPoints, xPointErrors, yPointErrors, nXYPoints);

        state = STATE.INITIALIZED;
    }

    public TwoPointCorrelation(float[] xPoints, float[] yPoints, float[] xPointErrors, float[] yPointErrors, int nXYPoints) {

        this.indexer = new DoubleAxisIndexer();

        indexer.sortAndIndexXThenY(xPoints, yPoints, xPointErrors, yPointErrors, nXYPoints);

        state = STATE.INITIALIZED;
    }

    public TwoPointCorrelation(String indexerFilePath) throws IOException {

        this.indexer = SerializerUtil.readPersistedPoints(indexerFilePath, true);

        state = STATE.INITIALIZED;
    }

    public TwoPointCorrelation(DoubleAxisIndexer doubleAxisIndexer) throws IOException {

        this.indexer = doubleAxisIndexer;

        state = STATE.INITIALIZED;
    }

    public void setSigmaFactorToTwo() {
        if (debug) {
            log.info("threshhold=2.0");
        }
        sigmaFactor = 2.0f;
    }
    public void setSigmaFactorToTwoPointFive() {
        if (debug) {
            log.info("threshhold=2.5");
        }
        sigmaFactor = 2.5f;
    }
    public void setSigmaFactorToThree() {
        if (debug) {
            log.info("threshhold=3.0");
        }
        sigmaFactor = 3.0f;
    }

    public void setDebug(boolean turnDebugOn) {
        this.debug = turnDebugOn;
    }
    
    /**
     * if letting the code fit the background distribution, this method will use
     * a downhill simplex method wrapped in range searches to fit the distribution,
     * else the default method will be used which is a non-quadratic conjugate gradient
     * solver. 
     */
    /*
    public void setUseDownhillSimplexHistogramFitting() {
        this.setUseDownhillSimplexHistogramFitting = true;
    }*/

    protected void logPerformanceMetrics() {
        this.doLogPerformanceMetrics = true;
    }

    public void persistIndexer(boolean doPersistIndexer) throws IOException {
        indexerFilePath = SerializerUtil.serializeIndexer(indexer);
    }

    public void setPersistMinimaStats(boolean doPersistMinimaStats) {
        persistTheMinimaStats = doPersistMinimaStats;
    }

    public void setMinimumNumberInCluster(int minimumNumberForClusterMembership) {
        this.minimumNumberInCluster = minimumNumberForClusterMembership;
    }
    
    public void setAllowRefinement() {
        allowRefinement = true;
    }

    public void setBackground(float backgroundSurfaceDensity, float standardDeviationOfBackground) {

        this.backgroundSurfaceDensity = backgroundSurfaceDensity;

        this.backgroundError = standardDeviationOfBackground;

        state = STATE.BACKGROUND_SET;

        bMethod = BACKGROUND_METHOD.USER_SUPPLIED;
    }

    /**
     * calculate background using the default method.
     *
     * @see TwoPointVoidStats.calc()
     *
     * @throws TwoPointVoidStatsException
     * @throws IOException
     */
    public void calculateBackground() throws TwoPointVoidStatsException, IOException {

        if ((bMethod == null) || (bMethod.ordinal() != BACKGROUND_METHOD.USER_SUPPLIED.ordinal())) {

            calculateBackgroundVia2PtVoidFit();
        }
    }

    public void reuseStatsForBackgroundCalculation(String minimaFilePath) throws TwoPointVoidStatsException, IOException {

        TwoPointVoidStats minStats = new TwoPointVoidStats(indexer);
        minStats.setDebug(debug);
        if (setUseDownhillSimplexHistogramFitting) {
            minStats.setUseDownhillSimplexHistogramFitting();
        }
        
        minStats.setStandardDeviationFactor(sigmaFactor);
        
        minStats.calc(minimaFilePath);

        if (debug) {
            backgroundStats = minStats;
        }

        this.backgroundSurfaceDensity = minStats.getBackgroundSurfaceDensity();
        this.backgroundError = minStats.getBackgroundSurfaceDensityError();

        if (debug) {
            log.info("background density ="
                + this.backgroundSurfaceDensity + " with error =" + this.backgroundError);
        }

        state = STATE.BACKGROUND_SET;

        bMethod = BACKGROUND_METHOD.DESERIALIZED;
    }

    protected void calculateBackgroundVia2PtVoidFit() throws TwoPointVoidStatsException, IOException {

        if ((bMethod != null) && (bMethod.ordinal() == BACKGROUND_METHOD.USER_SUPPLIED.ordinal())) {
            return;
        }

        TwoPointVoidStats voidStats = null;
        
        if (refineSolution.booleanValue()) {

            voidStats = new TwoPointVoidStats(tempRefineSolnIndexer);
            
            voidStats.setUseCompleteSampling();
           
            voidStats.setInterpretForSparseBackgroundToTrue();
                                
        } else {
            
            voidStats = new TwoPointVoidStats(indexer);
        }
        
        voidStats.setDebug(debug);
        
        voidStats.setStandardDeviationFactor(sigmaFactor);

        if (doLogPerformanceMetrics) {
            voidStats.logPerformanceMetrics();
        }
        
        if (setUseDownhillSimplexHistogramFitting) {
            voidStats.setUseDownhillSimplexHistogramFitting();
        }

        voidStats.calc();

        backgroundStats = voidStats;

        if (persistTheMinimaStats) {
            minimaStatsFilePath = voidStats.persistTwoPointBackground();
        }

        this.backgroundSurfaceDensity = voidStats.getBackgroundSurfaceDensity();
        this.backgroundError = voidStats.getBackgroundSurfaceDensityError();
        
        //backgroundStats.releaseLargeVariables();
        
        if (refineSolution) {
            tempRefineSolnIndexer = null;
        }

        if (debug) {
            log.info("==>background density ="
                + this.backgroundSurfaceDensity + " with error =" + this.backgroundError);
        }

        state = STATE.BACKGROUND_SET;

        bMethod = BACKGROUND_METHOD.FIT_TWO_POINT_VOIDS;
    }

    protected void printPerformanceMetrics(long startTimeMillis, long stopTimeMillis, String methodName, int nPoints) {

        long diffSec = (stopTimeMillis - startTimeMillis)/1000;

        MemoryMXBean mbean = ManagementFactory.getMemoryMXBean();
        MemoryUsage heapUsage = mbean.getHeapMemoryUsage();
        MemoryUsage nonHeapUsage = mbean.getNonHeapMemoryUsage();

        String str = String.format("%35s:  N=%9d  RT(sec)=%8d  instance estimates(bytes)=%9d   heapUsed(bytes)=%9d   memoryPoolsSum(bytes)=%9d",
            methodName,
            nPoints, diffSec, approximateMemoryUsed(),
            heapUsage.getUsed(), nonHeapUsage.getUsed() );

        Logger.getLogger(this.getClass().getSimpleName()).info(str);
    }

    public long approximateMemoryUsed() {

        String arch = System.getProperty("sun.arch.data.model");

        boolean is32Bit = ((arch != null) && arch.equals("64")) ? false : true;

        int nbits = (is32Bit) ? 32 : 64;

        int overheadBytes = 16;

        int intBytes = (is32Bit) ? 4 : 8;
        int arrayBytes = 32/8;
        int refBytes = nbits/8;

        /*
         * enums:  one has 4 items
         *         one has 3 items
         */
        long sumBits = 4*nbits;
        long tmpSumBytes = (sumBits/8) + overheadBytes;
        long padding = (tmpSumBytes % 8);
        long sumBytes = tmpSumBytes + padding;
        sumBits = 3*nbits;
        tmpSumBytes = (sumBits/8) + overheadBytes;
        padding = (tmpSumBytes % 8);
        sumBytes += (tmpSumBytes + padding);

        // a reference to each of the enums
        sumBytes += (2*intBytes);

        sumBytes += indexer.approximateMemoryUsed();

        if (groupHullIndexes != null) {
            sumBytes += (arrayBytes + (groupHullIndexes.length*(overheadBytes + refBytes + intBytes)));
        }
        if (xGroupHullCentroids != null) {
            // float[] on the heap
            sumBytes += 3*(arrayBytes + (xGroupHullCentroids.length*(arrayBytes)));
        }
        if (groupFinder != null) {
            sumBytes += groupFinder.approximateMemoryUsed();
        }
        if (backgroundStats != null) {
            sumBytes += backgroundStats.approximateMemoryUsed();
        }
        if (indexerFilePath != null) {
            // String size on the heap = reference size + content size?
            sumBytes += (intBytes + (indexerFilePath.length()*intBytes));
        }
        if (minimaStatsFilePath != null) {
            sumBytes += (intBytes + (minimaStatsFilePath.length()*intBytes));
        }

        // 8 variables in the stack are each word size
        sumBytes += 8 * intBytes;

        // log is reference size on the heap
        sumBytes += intBytes;

        sumBytes += overheadBytes;

        // amount of padding needed to make it a round 8 bytes
        padding = (sumBytes % 8);

        sumBytes += padding;

        return sumBytes;
    }

    public void findClusters() throws TwoPointVoidStatsException, IOException {

        if (state.ordinal() < STATE.BACKGROUND_SET.ordinal()) {
            calculateBackgroundVia2PtVoidFit();
        }

        findGroups();
    }
    
    public IGroupFinder getGroupFinder() {
        return groupFinder;
    }

    public String plotClusters(float xMin, float xMax, float yMin, float yMax)
        throws FileNotFoundException, IOException, TwoPointVoidStatsException {

        if (state.ordinal() < STATE.CLUSTER_HULLS_CALCULATED.ordinal()) {
            calculateHullsOfClusters();
        }

        TwoPointCorrelationPlotter plotter = new TwoPointCorrelationPlotter(xMin, xMax, yMin, yMax);
        plotter.addPlot(this);

        return plotter.writeFile();
    }

    public String plotClusters() throws FileNotFoundException, IOException, TwoPointVoidStatsException {

        if (state.ordinal() < STATE.CLUSTER_HULLS_CALCULATED.ordinal()) {
            calculateHullsOfClusters();
        }

        float[] xMinMax = MiscMath.calculateOuterRoundedMinAndMax(indexer.getX());
        float[] yMinMax = MiscMath.calculateOuterRoundedMinAndMax(indexer.getY());

        TwoPointCorrelationPlotter plotter = new TwoPointCorrelationPlotter(xMinMax[0], xMinMax[1], yMinMax[0], yMinMax[1]);
        String label = "";
        if ((this.backgroundStats != null) && (this.backgroundStats instanceof TwoPointVoidStats)) {
            GEVYFit bestFit = ((TwoPointVoidStats)this.backgroundStats).bestFit;
            if (bestFit != null) {
                label = String.format("k=%.7f sigma=%.7f mu=%.7f", bestFit.getK(), bestFit.getSigma(), bestFit.getMu());
            }
        }
        plotter.addPlot(this, label);

        return plotter.writeFile3();
    }

    protected void findGroups() throws TwoPointVoidStatsException, IOException {
        
        long startTimeMillis = System.currentTimeMillis();

        groupFinder = new DFSGroupFinder(backgroundSurfaceDensity, sigmaFactor);
        
        groupFinder.setMinimumNumberInCluster(minimumNumberInCluster);
        
        groupFinder.findGroups(indexer);
        
        state = STATE.CLUSTERS_FOUND;

        if (doLogPerformanceMetrics) {

            long stopTimeMillis = System.currentTimeMillis();

            printPerformanceMetrics(startTimeMillis, stopTimeMillis, "findGroups", indexer.getNXY());
        }
        
        if (allowRefinement && !refineSolution && backgroundStats != null && backgroundStats instanceof TwoPointVoidStats) {
            
            TwoPointVoidStats tmp = (TwoPointVoidStats)backgroundStats;
            
            if (tmp.getInterpretForSparseBackground() != null && tmp.getInterpretForSparseBackground().booleanValue()) {
            
                if (tmp.getSampling() != null && tmp.getSampling().ordinal() == VoidSampling.COMPLETE.ordinal()) {
                    
                    tempRefineSolnIndexer = createIndexerMinusGroupPoints();
                    
                    if (tempRefineSolnIndexer != null) {

                        // subtract the groups to create a new indexer
                        refineSolution = Boolean.TRUE;

                        state = STATE.INITIALIZED;

                        findClusters();
                        
                    }
                }
            }
        }
    }
  
    private DoubleAxisIndexer createIndexerMinusGroupPoints() {
        
        if (indexer == null) {
            throw new IllegalStateException("indexer cannot be null");
        }
        if (groupFinder == null) {
            throw new IllegalStateException("groupFinder cannot be null");
        }
        
        int[] pointToGroupIndexes = groupFinder.getPointToGroupIndexes();
        
        int numberInGroups = 0;
        for (int idx : pointToGroupIndexes) {
            if (idx > -1) {
                numberInGroups++;
            }
        }
        
        int n = indexer.getNumberOfPoints() - numberInGroups;
        
        if (n == 0) {
            return null;
        }
        
        float[] tmpx = new float[n];
        float[] tmpy = new float[n];
        float[] tmpxe = new float[n];
        float[] tmpye = new float[n];
        
        int count = 0;
        
        for (int i = 0; i < pointToGroupIndexes.length; i++) {
            
            int groupId = pointToGroupIndexes[i];
            
            if (groupId == -1) {
             
                tmpx[count] = indexer.getX()[i];
                tmpy[count] = indexer.getY()[i];
                tmpxe[count] = indexer.getXErrors()[i];
                tmpye[count] = indexer.getYErrors()[i];
                
                count++;
            }
        }
        
        DoubleAxisIndexer tmpIndexer = new DoubleAxisIndexer();
        tmpIndexer.sortAndIndexXThenY(tmpx, tmpy, tmpxe, tmpye, tmpx.length);
                
        return tmpIndexer;
    }

    protected float calculateFractionOfPointsOutsideOfClusters() {

        if (groupFinder == null) {
            // there are no clusters, so fraction outside is 1.0
            return 1.0f;
        }

        int nGroups = groupFinder.getNumberOfGroups();
        
        boolean[] insideClusters = new boolean[indexer.getNXY()];

        for (int i = 0; i < nGroups; i++) {

            SimpleLinkedListNode groupNode = groupFinder.getGroupMembershipList()[i];

            while ((groupNode != null) && (groupNode.key != -1)) {

                int pointIndex = groupNode.key;

                insideClusters[pointIndex] = true;

                groupNode = groupNode.next;
            }
        }

        // count number outside
        int count = 0;
        for (int i = 0; i < insideClusters.length; i++) {
            if (!insideClusters[i]) {
                count++;
            }
        }

        float frac = (float)count/(float)indexer.getNXY();

        return frac;
    }

    /**
     * calculate convex hulls of clusters and calculate the hull centroids and area
     *
     * @throws IOException
     * @throws TwoPointVoidStatsException
     */
    public void calculateHullsOfClusters() throws IOException, TwoPointVoidStatsException {

        if (state.ordinal() < STATE.CLUSTERS_FOUND.ordinal()) {
            findClusters();
        }

        long startTimeMillis = System.currentTimeMillis();
        
        int nGroups = (groupFinder != null) ? groupFinder.getNumberOfGroups() : 0;

        groupHullIndexes = new SimpleLinkedListNode[nGroups];

        xGroupHullCentroids = new float[nGroups];
        yGroupHullCentroids = new float[nGroups];
        groupHullSurfaceAreas = new float[nGroups];

        // the indexes stored in the instance vars such as pointToGroupIndex are w.r.t. the arrays sorted by y
        float[] x = indexer.getXSortedByY();
        float[] y = indexer.getYSortedByY();

        int nXY = indexer.getNXY();

        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter(indexer.getX()[indexer.getSortedXIndexes()[0]],
            indexer.getX()[indexer.getSortedXIndexes()[nXY - 1]],
            indexer.getY()[indexer.getSortedYIndexes()[0]], indexer.getY()[indexer.getSortedYIndexes()[nXY - 1]]);

        for (int i = 0; i < nGroups; i++) {

            groupHullIndexes[i] = new SimpleLinkedListNode();

            int[] memberIndexes = getGroupFinder().getIndexes(i);
            float[] xMember = getGroupFinder().getX(i, indexer);
            float[] yMember = getGroupFinder().getY(i, indexer);
            
            float[] xhull = null;
            float[] yhull = null;

            if (xMember.length > 2) {
                try {
                    GrahamScan scan = new GrahamScan();
                    scan.computeHull(xMember, yMember);

                    xhull = scan.getXHull();
                    yhull = scan.getYHull();
                } catch (GrahamScanTooFewPointsException e) {
                }
            }

            if (xhull == null) {

                xhull = new float[xMember.length];
                yhull = new float[xMember.length];
                for (int ii = 0; ii < xMember.length; ii++) {
                    xhull[ii] = xMember[ii];
                    yhull[ii] = yMember[ii];
                }
            }

            int hullMatchCount = 0;
            // the hull is a subset of memberIndexes so point to the original coordinates
            //   instead of creating a new instance data structure
            for (int j = 0; j < xhull.length; j++) {

                float xh = xhull[j];
                float yh = yhull[j];

                for (int k = 0; k < xMember.length; k++) {

                    int pointIndex = memberIndexes[k];

                    if ( (x[pointIndex] == xh) && (y[pointIndex] == yh) ) {
                        groupHullIndexes[i].insert(pointIndex);
                        hullMatchCount++;
                        break;
                    }
                }
            }

            if (debug) {
                try {
                    plotter.addPlot(xMember, yMember, xhull, yhull, "");
                    plotter.writeFile();
                } catch (IOException e) {

                }
            }

            float[] ca = LinesAndAngles.calcAreaAndCentroidOfSimplePolygon(xhull, yhull);

            if (ca != null) {
                groupHullSurfaceAreas[i] = ca[0];
                xGroupHullCentroids[i] = ca[1];
                yGroupHullCentroids[i] = ca[2];
            }
        }

        state = STATE.CLUSTER_HULLS_CALCULATED;

        if (doLogPerformanceMetrics) {

            long stopTimeMillis = System.currentTimeMillis();

            printPerformanceMetrics(startTimeMillis, stopTimeMillis, "calculateHullsOfClusters", nGroups);
        }
    }

    public int getNumberOfGroups() {
        return (groupFinder != null) ? groupFinder.getNumberOfGroups() : 0;
    }

    public float[] calculateGroupCentroidUsingAllPointsEquallyWeighted(int groupNumber) {

        if (groupFinder == null) {
            return null;
        }
        
        int nGroups = groupFinder.getNumberOfGroups() ;
        
        if (groupNumber >= nGroups) {
            throw new IllegalArgumentException("groupNumber is larger than existing number of groups");
        }

        float[] xMember = getGroupFinder().getX(groupNumber, indexer);
        float[] yMember = getGroupFinder().getY(groupNumber, indexer);
        
        float xCoordsAvg = 0;
        float yCoordsAvg = 0;
        int count = 0;
        
        for (int i = 0; i < xMember.length; i++) {
            xCoordsAvg += xMember[i];
            yCoordsAvg += yMember[i];
        }

        xCoordsAvg /= (float)count;
        yCoordsAvg /= (float)count;

        return new float[]{xCoordsAvg, yCoordsAvg};
    }

    public float[] getXGroupHull(int groupNumber) {

        if (groupFinder == null) {
            return null;
        }
        
        int nGroups = groupFinder.getNumberOfGroups() ;
        
        if (groupNumber >= nGroups) {
            throw new IllegalArgumentException("groupNumber is larger than existing number of groups");
        }

        float[] x = indexer.getX();

        SimpleLinkedListNode hullNode = groupHullIndexes[groupNumber];

        float[] xhull = new float[hullNode.getKeys().length];
        int count = 0;

        while ((hullNode != null) && (hullNode.key != -1)) {

            int pointIndex = hullNode.key;

            xhull[count] = x[pointIndex];

            hullNode = hullNode.next;
            count++;
        }

        return xhull;
    }

    public float[] getYGroupHull(int groupNumber) {

        if (groupFinder == null) {
            return null;
        }
        
        int nGroups = groupFinder.getNumberOfGroups() ;
        
        if (groupNumber >= nGroups) {
            throw new IllegalArgumentException("groupNumber is larger than existing number of groups");
        }

        // the indexes stored in the instance vars such as pointToGroupIndex are w.r.t. the arrays sorted by y
        float[] y = indexer.getY();

        SimpleLinkedListNode hullNode = groupHullIndexes[groupNumber];

        float[] yhull = new float[hullNode.getKeys().length];
        int count = 0;

        while ((hullNode != null) && (hullNode.key != -1)) {

            int pointIndex = hullNode.key;

            yhull[count] = y[pointIndex];

            hullNode = hullNode.next;
            count++;
        }

        return yhull;
    }

    public float[] getXHullCentroids() {
        
        if (groupFinder == null) {
            return null;
        }
        
        int nGroups = groupFinder.getNumberOfGroups() ;
        
        float[] seeds = new float[nGroups];
        for (int i = 0; i < nGroups; i++) {
            seeds[i] = getXGroupHullCentroid(i);
        }
        return seeds;
    }
    public float[] getYHullCentroids() {
        
        if (groupFinder == null) {
            return null;
        }
        
        int nGroups = groupFinder.getNumberOfGroups() ;
        
        float[] seeds = new float[nGroups];
        for (int i = 0; i < nGroups; i++) {
            seeds[i] = getYGroupHullCentroid(i);
        }
        return seeds;
    }

    public float getXGroupHullCentroid(int groupNumber) {

        if (groupFinder == null) {
            return Float.MIN_VALUE;
        }
        
        int nGroups = groupFinder.getNumberOfGroups();
        
        if (groupNumber >= nGroups) {
            throw new IllegalArgumentException("groupNumber is larger than existing number of groups");
        }

        return xGroupHullCentroids[groupNumber];
    }
    public float getYGroupHullCentroid(int groupNumber) {

        if (groupFinder == null) {
            return Float.MIN_VALUE;
        }
        
        int nGroups = groupFinder.getNumberOfGroups();
        
        if (groupNumber >= nGroups) {
            throw new IllegalArgumentException("groupNumber is larger than existing number of groups");
        }

        return yGroupHullCentroids[groupNumber];
    }

    public float getBackgroundSurfaceDensity() {
        return backgroundSurfaceDensity;
    }

    public float getBackgroundSurfaceDensityError() {
        return backgroundError;
    }

    public float[] getX() {
        return indexer.getX();
    }

    public float[] getY() {
        return indexer.getY();
    }
    public float[] getXErrors() {
        return indexer.getXErrors();
    }
    public float[] getYErrors() {
        return indexer.getYErrors();
    }
    DoubleAxisIndexer getIndexer() {
        return indexer;
    }

    public float[] getXGroup(int groupNumber) {

        return (groupFinder != null) ? groupFinder.getX(groupNumber, indexer) : new float[0];
    }

    public float[] getYGroup(int groupNumber) {

        return (groupFinder != null) ? groupFinder.getY(groupNumber, indexer) : new float[0];
    }

}
