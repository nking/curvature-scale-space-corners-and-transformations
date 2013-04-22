package algorithms.compGeometry.clustering.twopointcorrelation;

import algorithms.misc.HistogramHolder;
import algorithms.util.ResourceFinder;
import java.util.logging.Logger;

/**
 * @author nichole
 */
public class FindClusters2Test extends BaseTwoPointTest {

    boolean debug = true;

    protected Logger log = Logger.getLogger(this.getClass().getSimpleName());

    public void test_Find_Clusters_Stats() throws Exception {

        log.info("test_Find_Clusters_Stats()");

        float xmin = 0;
        float xmax = 300;
        float ymin = 0;
        float ymax = 300;

        TwoPointCorrelationPlotter plotter = new TwoPointCorrelationPlotter(xmin, xmax, ymin, ymax);

        String[] filePaths = CreateClusterDataTest.getIndexerFilePaths();

        if (filePaths == null || filePaths.length == 0) {
            return;
        }

        for (int i = 0; i < filePaths.length; i++) {
        //for (int i = 150; i < filePaths.length; i++) {

            String filePath = filePaths[i];

            DoubleAxisIndexer indexer = CreateClusterDataTest.readIndexer(filePath);

            String srchFor = "indexer_random_background_with_";
            int i0 = filePath.indexOf(srchFor);
            String numberOfClusters = filePaths[i].substring(i0 + srchFor.length());
            i0 = numberOfClusters.indexOf("_clusters_");
            numberOfClusters = numberOfClusters.substring(0, i0);

            String strCount;
            i0 = filePath.indexOf("clusters_");
            i0 += "clusters_".length();
            int i1 = filePath.indexOf(".dat");
            strCount = filePath.substring(i0, i1);

            TwoPointCorrelation twoPtC = new TwoPointCorrelation(indexer);
            twoPtC.setDebug(false);

            log.info(" " + i + " (" + twoPtC.indexer.nXY + " points) ... ");

            //twoPtC.setBackground(0.035f, 0.00994f/10);
            twoPtC.findClusters();

            TwoPointVoidStats stats = (TwoPointVoidStats)twoPtC.backgroundStats;

            String plotLabel = "";

            if (stats != null && ((TwoPointVoidStats)twoPtC.backgroundStats).bestFit != null) {

                // centroid of area defined by the top portion of the fit where y >= ypeak/2
                float[] areaAndXYTopCentroid = stats.calculateCentroidOfTop(
                    stats.bestFit.getOriginalScaleX(), stats.bestFit.getOriginalScaleYFit(), 0.5f);

                if (areaAndXYTopCentroid != null) {

                    String samplingAbbrev = "";
                    switch(((TwoPointVoidStats)twoPtC.backgroundStats).getSampling().ordinal()) {
                        case 0:
                            samplingAbbrev = "C";
                            break;
                        case 1:
                            samplingAbbrev = "L_C";
                            break;
                        case 2:
                            samplingAbbrev = "S_C";
                            break;
                        case 3:
                            samplingAbbrev = "S_C_R_S";
                            break;
                    }

                    plotLabel = String.format(
                    "  (%3d %s %4d)  peak=%.4f  xcen=%.4f  chst=%.1f  %s",
                    i, numberOfClusters, twoPtC.indexer.nXY,
                    ((TwoPointVoidStats)twoPtC.backgroundStats).bestFit.getXPeak(),
                    areaAndXYTopCentroid[1],
                    ((TwoPointVoidStats)twoPtC.backgroundStats).bestFit.getChiSqStatistic(),
                    samplingAbbrev
                    );
                }

            } else {

                plotLabel = String.format(
                    "  (%3d %s %4d)",
                    i, numberOfClusters, twoPtC.indexer.nXY);
            }

            twoPtC.calculateHullsOfClusters();

            float[] xf = null;
            float[] yf = null;
            HistogramHolder histogram = stats.getStatsHistogram();
            String fileNamePostfix = "_clusters_" + strCount + ".dat";
            String fileName = CreateClusterDataTest.histogramFileNamePrefix + numberOfClusters + fileNamePostfix;
            filePath = ResourceFinder.getAFilePathInTmpData(fileName);
            CreateClusterDataTest.writeHistogram(filePath, histogram);

            plotter.addPlot(twoPtC, plotLabel);
            plotter.writeFile();

            log.info("  END " + i);
        }
    }

}
