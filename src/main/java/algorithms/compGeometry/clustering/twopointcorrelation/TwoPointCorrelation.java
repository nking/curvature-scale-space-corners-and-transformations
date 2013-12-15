package algorithms.compGeometry.clustering.twopointcorrelation;

import algorithms.compGeometry.LinesAndAngles;
import algorithms.compGeometry.convexHull.GrahamScanTooFewPointsException;
import algorithms.compGeometry.convexHull.GrahamScan;
import algorithms.curves.GEVYFit;
import algorithms.misc.MiscMath;
import algorithms.util.ArrayPair;
import algorithms.util.Errors;
import algorithms.util.PolygonAndPointPlotter;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.lang.management.ManagementFactory;
import java.lang.management.MemoryMXBean;
import java.lang.management.MemoryUsage;
import java.util.logging.Logger;

/**
  Find clusters in data.
  
  The default method expects that there are background points in the data which are not part
  of the groups.  The groups are found by statistics:  they are clustered 2 to 3 times more
  than the background points.  The density of the background is determined, and pairs of 
  points that are closer than 2 to 3 times that background density are in a group.
     
  There is an alternative method that can be invoked on datasets in which there are no
  points outside of groups, that is no background points.
     to use the alternative method:
         useFindMethodForDataWithoutBackgroundPoints()
         
  You should chose the default method or the alternative method prior to using the code
  because they are performed on very different kinds of datasets.
  If for some reason, there is a need to have the code automate the decision, you can use:
         automateTheChoiceOfFindMethod();
  *Note that this isn't recommended simply because the dataset types for the different 
  methods are very different, so not likely to be present in different runs of your datasets.
  
  
  Use as an API:
      TwoPointCorrelation clusterFinder = new TwoPointCorrelation(x, y, xErrors, yErrors, getTotalNumberOfPoints());
      clusterFinder.calculateBackground();
      clusterFinder.findClusters();
      
      // the results are available as group points or as convex hulls surrounding the groups:
      int n = clusterFinder.getNumberOfGroups()
      
      // get the hull for groupId 0
      ArrayPair hull0 = clusterFinder.getGroupHull(0)

      // get the points in groupId 0
      ArrayPair group0 = clusterFinder.getGroup(int groupNumber)
      
      // plot the results:
      String plotFilePath = clusterFinder.plotClusters();

  To set the background density manually:
      TwoPointCorrelation clusterFinder = new TwoPointCorrelation(x, y, xErrors, yErrors, getTotalNumberOfPoints());
      clusterFinder.setBackground(0.03f, 0.003f);
      clusterFinder.findClusters();
      String plotFilePath = clusterFinder.plotClusters();


  Note:  For datasets in which the density of background points is high, if you don't
     have the ability to filter the data by a key characteristic, you might consider
     the results of this code as seeds for a Voronoi diagram or other code.
     
         ArrayPair seeds = clusterFinder.getHullCentroids();


  Use from the command line:
      Requires a tab delimited text file with 4 columns: x, y, xErrors, yErrors.

          java -cp bin/classes  algorithms.compGeometry.clustering.twopointcorrelation.TwoPointCorrelation --file /path/to/file/fileName.txt

  @author nichole
 */
public class TwoPointCorrelation {

    protected enum STATE {
        INITIALIZED, BACKGROUND_SET, CLUSTERS_FOUND
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

    protected boolean persistTheMinimaStats = false;
    protected String indexerFilePath = null;
    protected String minimaStatsFilePath = null;

    // for debugging plots, keeping a handle on TwoPointVoidStats.
    public IPointBackgroundStats backgroundStats = null;

    protected boolean debug = false;

    protected boolean doLogPerformanceMetrics = false;

    protected boolean setUseDownhillSimplexHistogramFitting = false;
    
    protected boolean useFindMethodForSparseBackground = false;
    
    protected boolean automateTheFindMethodChoice = false;

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
     * for datasets where you know that there are no points outside of the groups,
     * that is, there are no background points, set the method choice to
     * that for no background points here.  This has to be set before
     * findClusters() is invoked.
     */
    public void useFindMethodForDataWithoutBackgroundPoints() {
        if (automateTheFindMethodChoice) {
            throw new IllegalStateException("useFindMethodForSparseBackground and automateTheChoiceOfFindMethod cannot both be set");
        }
        this.useFindMethodForSparseBackground = true;
    }
    
    /**
     * this method is not recommended, but is to allow the code to choose between
     * the default method of finding clusters based upon a density above the background
     * density or the alternative method for finding clusters when there are no background
     * points between the groups.
     * The reason that it is not recommended is that the 2 types of datasets that the methods
     * would be used on are so different
     * that you likely know which you have before use of this code.
     */
    public void automateTheChoiceOfFindMethod() {
        if (useFindMethodForSparseBackground) {
            throw new IllegalStateException("useFindMethodForSparseBackground and automateTheChoiceOfFindMethod cannot both be set");
        }
        this.automateTheFindMethodChoice = true;
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
        
        if (useFindMethodForSparseBackground) {
            minStats.setInterpretForSparseBackgroundToTrue();
        } else if (automateTheFindMethodChoice) {
            minStats.automateTheFindMethodChoice();
        }

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

        if (useFindMethodForSparseBackground) {
            voidStats.setInterpretForSparseBackgroundToTrue();
        } else if (automateTheFindMethodChoice) {
            voidStats.automateTheFindMethodChoice();
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

        TwoPointCorrelationPlotter plotter = new TwoPointCorrelationPlotter(xMin, xMax, yMin, yMax);
        plotter.addPlot(this);

        return plotter.writeFile();
    }

    public String plotClusters() throws FileNotFoundException, IOException, TwoPointVoidStatsException {

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

        if ((this.backgroundStats != null) && (this.backgroundStats instanceof TwoPointVoidStats)
            && ((TwoPointVoidStats)backgroundStats).getInterpretForSparseBackground().booleanValue()) {

            // for method without any background points, we use the density of the edge points without a factor
            groupFinder = new DFSGroupFinder(backgroundSurfaceDensity, 1.0f);

        } else {

            groupFinder = new DFSGroupFinder(backgroundSurfaceDensity, sigmaFactor);
        }

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

                        /*
                        float[] xymm = tempRefineSolnIndexer.findXYMinMax();
                        PolygonAndPointPlotter p0 = new PolygonAndPointPlotter();
                        p0.addPlot(indexer.getX(), indexer.getY(),
                            indexer.getXErrors(), indexer.getYErrors(), "original");
                        p0.writeFile3();
                        p0.addPlot(tempRefineSolnIndexer.getX(), tempRefineSolnIndexer.getY(),
                            tempRefineSolnIndexer.getXErrors(), tempRefineSolnIndexer.getYErrors(), "refining...");
                        System.out.println(p0.writeFile3());
                        */

                        findClusters();

                        tempRefineSolnIndexer = null;
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
    
    public float[] calculateAreaAndCentroidOfHull(float[] xHull, float[] yHull) {
        if (xHull == null || yHull == null) {
            throw new IllegalArgumentException("neither xHull nor yHull can be null");
        }
        return LinesAndAngles.calcAreaAndCentroidOfSimplePolygon(xHull, yHull);
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

    public ArrayPair getGroupHull(int groupNumber) {

        if (groupFinder == null) {
            return null;
        }

        int nGroups = groupFinder.getNumberOfGroups() ;

        if (groupNumber >= nGroups) {
            throw new IllegalArgumentException("groupNumber is larger than existing number of groups");
        }

        float[] xg = groupFinder.getX(groupNumber, indexer);
        
        float[] yg = groupFinder.getY(groupNumber, indexer);
        
        try {
            
            GrahamScan scan = new GrahamScan();
            
            scan.computeHull(xg, yg);

            return new ArrayPair(scan.getXHull(), scan.getYHull());

        } catch (GrahamScanTooFewPointsException e) {

            return new ArrayPair(new float[0], new float[0]);
        }
    }

    public ArrayPair getHullCentroids() {

        if (groupFinder == null) {
            return null;
        }

        int nGroups = groupFinder.getNumberOfGroups() ;

        float[] xc = new float[nGroups];
        float[] yc = new float[nGroups];
        
        for (int i = 0; i < nGroups; i++) {
            
            ArrayPair hull = getGroupHull(i);
            
            float[] ca = LinesAndAngles.calcAreaAndCentroidOfSimplePolygon(hull.getX(), hull.getY());
            
            xc[i] = ca[1];
            
            yc[i] = ca[2];
            
        }
        
        return new ArrayPair(xc, yc);
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
    
    public ArrayPair getGroup(int groupNumber) {
        
        if (groupFinder == null) {
            throw new IllegalStateException("groupFinder is null.  Please run findClusters() first.");
        }
                    
        float[] xg = groupFinder.getX(groupNumber, indexer);
        
        float[] yg = groupFinder.getY(groupNumber, indexer);
        
        return new ArrayPair(xg, yg);      
    }

}
