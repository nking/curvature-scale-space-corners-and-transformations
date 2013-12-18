package algorithms.compGeometry.clustering.twopointcorrelation;

import algorithms.misc.HistogramHolder;
import algorithms.misc.MiscMath;
import algorithms.util.ArrayPair;
import algorithms.util.ResourceFinder;
import java.util.logging.Logger;
import junit.framework.TestCase;

/**
 * @author nichole
 */
public class FindClusters3Test extends TestCase {

    public FindClusters3Test(String testName) {
        super(testName);
    }

    @Override
    protected void setUp() throws Exception {
        super.setUp();
    }

    @Override
    protected void tearDown() throws Exception {
        super.tearDown();
    }

    boolean debug = true;

    boolean persistHistogram = false;

    protected Logger log = Logger.getLogger(this.getClass().getSimpleName());

    public void testFindClusters() throws Exception {

        log.info("testFindClusters()");

        DoubleAxisIndexer indexer = CreateClusterDataTest.getWikipediaDBScanExampleData();
        // dug these points out of http://upload.wikimedia.org/wikipedia/commons/0/05/DBSCAN-density-data.svg

        float xmin = MiscMath.findMin(indexer.getX());
        float xmax = MiscMath.findMax(indexer.getX());
        float ymin = MiscMath.findMin(indexer.getY());
        float ymax = MiscMath.findMax(indexer.getY());

        TwoPointCorrelationPlotter plotter = new TwoPointCorrelationPlotter(xmin, xmax, ymin, ymax);

        // generates fake errors when errors not present:
        TwoPointCorrelation twoPtC = new TwoPointCorrelation(indexer);
        twoPtC.setDebug(true);

        // 390 points

        //twoPtC.setBackground(0.2725f, 0.001f);
        //twoPtC.findClusters();

        //twoPtC.useFindMethodForDataWithoutBackgroundPoints();

        twoPtC.calculateBackground();
        twoPtC.findClusters();

        TwoPointVoidStats stats = (TwoPointVoidStats)twoPtC.backgroundStats;

        String plotLabel = "";

        if (twoPtC.backgroundStats != null && twoPtC.backgroundStats instanceof TwoPointVoidStats 
            && ((TwoPointVoidStats)twoPtC.backgroundStats).bestFit != null) {

            // centroid of area defined by the top portion of the fit where y >= ypeak/2
            float[] areaAndXYTopCentroid = TwoPointVoidStats.calculateCentroidOfTop(
                stats.bestFit.getOriginalScaleX(), stats.bestFit.getOriginalScaleYFit(), 0.5f);

            if (areaAndXYTopCentroid != null) {

                plotLabel = String.format(
                    "  (%4d)  k=%.4f  sigma=%.4f  mu=%.4f chst=%.1f",
                    twoPtC.indexer.nXY, ((TwoPointVoidStats)twoPtC.backgroundStats).bestFit.getK(),
                    ((TwoPointVoidStats)twoPtC.backgroundStats).bestFit.getSigma(),
                    ((TwoPointVoidStats)twoPtC.backgroundStats).bestFit.getMu(),
                    areaAndXYTopCentroid[1], ((TwoPointVoidStats)twoPtC.backgroundStats).bestFit.getChiSqStatistic());

            } else {

                plotLabel = String.format("  (%4d)", twoPtC.indexer.nXY);
            }
        }

        plotter.addPlotWithoutHull(twoPtC, plotLabel);
        plotter.writeFile();

        if (persistHistogram) {
            float[] xf = null;
            float[] yf = null;
            HistogramHolder histogram = stats.getStatsHistogram();
            String fileNamePostfix = "wikipedia_dbscan.dat";
            String fileName = CreateClusterDataTest.histogramFileNamePrefix + fileNamePostfix;
            String filePath = ResourceFinder.getAFilePathInTmpData(fileName);
            CreateClusterDataTest.writeHistogram(filePath, histogram);
        }

        assertTrue(twoPtC.getNumberOfGroups() == 2);
        
        ArrayPair centroids = twoPtC.getHullCentroids();
        
        assertTrue(centroids.getX().length == 2);
        assertTrue(centroids.getY().length == 2);
        
        log.info( twoPtC.indexer.nXY + " points ... ");

        log.info("  END ");
    }

}
