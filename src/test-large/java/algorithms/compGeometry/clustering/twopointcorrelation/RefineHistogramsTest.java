package algorithms.compGeometry.clustering.twopointcorrelation;

import algorithms.curves.GEVChiSquareMinimization;
import algorithms.curves.GEVYFit;
import algorithms.misc.HistogramHolder;
import algorithms.util.PolygonAndPointPlotter;
import algorithms.util.ResourceFinder;
import java.io.File;
import java.io.FileFilter;
import java.io.IOException;
import java.util.logging.Logger;

/**
 * @author nichole
 */
public class RefineHistogramsTest extends BaseTwoPointTest {

    boolean debug = true;

    boolean writeToTmpData = false;

    boolean findClusters = true;

    protected Logger log = Logger.getLogger(this.getClass().getName());

    public void testWriteDensityFiles() throws Exception {

        if (!writeToTmpData) {
            return;
        }

        log.info("test_Find_Clusters_Stats()");

        String[] filePaths = CreateClusterDataTest.getIndexerFilePaths();

        if (filePaths == null || filePaths.length == 0) {
            return;
        }

        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter();

        for (int i = 0; i < filePaths.length; i++) {

            // filtering for files indexer_random_background_with_3_clusters_ 150-179

            String filePath = filePaths[i];

            if (!filePath.contains("indexer_random_background_with_3_clusters_")) {
                continue;
            }
            int d0 = filePath.lastIndexOf(".dat");
            String nFileStr = filePath.substring(d0-3, d0);

            int nFile = Integer.valueOf(nFileStr);
            if ((nFile < 150) || (nFile > 179) ) {
                continue;
            }

            AxisIndexer indexer = CreateClusterDataTest.readIndexer(filePath);

            plotter.addPlot(indexer.getX(), indexer.getY(), indexer.getXErrors(), indexer.getYErrors(), null, null, String.valueOf(i));
            plotter.writeFile();


            // use the various sampling and persist the points.
            TwoPointVoidStats voidStats = null;
            String fileRootName = "stats_2pt_voids_" + nFileStr;
            String fileRootNameIndiv, statsFilePath;

            fileRootNameIndiv = fileRootName + "_" + TwoPointVoidStats.Sampling.LEAST_COMPLETE.name() + "_";
            voidStats = new TwoPointVoidStats(indexer);
            voidStats.setDebug(debug);
            voidStats.setUseLeastCompleteSampling();
            voidStats.calculateTwoPointVoidDensities();
            statsFilePath =
                voidStats.serializeTwoPointBackground(fileRootNameIndiv);
            log.info("writing " + statsFilePath);


            fileRootNameIndiv = fileRootName + "_" + TwoPointVoidStats.Sampling.SEMI_COMPLETE_RANGE_SEARCH.name() + "_";
            voidStats = new TwoPointVoidStats(indexer);
            voidStats.setDebug(debug);
            voidStats.setUseSemiCompleteRangeSampling();
            voidStats.calculateTwoPointVoidDensities();
            statsFilePath =
                voidStats.serializeTwoPointBackground(fileRootNameIndiv);
            log.info("writing " + statsFilePath);

/*
            fileRootNameIndiv = fileRootName + "_" + TwoPointVoidStats.Sampling.SEMI_COMPLETE.name() + "_";
            voidStats = new TwoPointVoidStats(indexer);
            voidStats.setDebug(debug);
            voidStats.setUseSemiCompleteSampling();
            voidStats.calculateTwoPointVoidDensities();
            statsFilePath =
                voidStats.serializeTwoPointBackground(fileRootNameIndiv);
            log.info("writing " + statsFilePath);
*/
            log.info("  END " + i);
        }
    }

    public void testRefineHistograms() throws Exception {

        log.info("testRefineHistograms()");

        String[] filePaths = getDensityFilePaths();

        if (filePaths == null || filePaths.length == 0) {
            return;
        }

        String[] indexerFilePaths = CreateClusterDataTest.getIndexerFilePaths();

        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter();

        //String fltr = TwoPointVoidStats.Sampling.LEAST_COMPLETE.name();
        //String fltr = TwoPointVoidStats.Sampling.SEMI_COMPLETE.name();
        String fltr = TwoPointVoidStats.Sampling.SEMI_COMPLETE_RANGE_SEARCH.name();

        for (int i = 0; i < filePaths.length; i++) {

            // filtering for files indexer_random_background_with_3_clusters_ 150-179

            String filePath = filePaths[i];

            if (!filePath.contains(fltr)) {
                continue;
            }

            String indexerFilePath = getIndexerFilePath(indexerFilePaths, filePath);

            if (indexerFilePath == null) {
                continue;
            }


            TwoPointVoidStats voidStats = new TwoPointVoidStats(indexerFilePath);
            voidStats.readTwoPointBackground(filePath);

            HistogramHolder histogram = voidStats.createHistogramWithHigherPeakResolution();

            GEVChiSquareMinimization chiSqMin
                = new GEVChiSquareMinimization(histogram.getXHist(), histogram.getYHistFloat(),
                histogram.getXErrors(), histogram.getYErrors());

            float k = 1.04f;// changes sharpness of left side slope
            float s = 0.195f;
            float kMin = k/2.f;
            float kMax = k*2.f;
            float sMin = s/2;
            float sMax = s*2.0f;

            float mu = 0.06f;
            float yNorm = 1.0f;
            float yErrSqSum = chiSqMin.calcYErrSquareSum();

            chiSqMin.setDebug(debug);

            // k=0.67109966 sigma=0.16777492 mu=0.016949153 chiSqSum=96.46002 chiSqStatistic=3.7100008
            //  trial and error:
            //   k=1.0402603 sigma=0.19504745 mu=0.06 chiSqSum=86.89037 chiSqStatistic=3.3419375

            //GEVYFit yfit = chiSqMin.fitCurve(kMin, kMax, sMin, sMax, mu, yErrSqSum, GEVChiSquareMinimization.WEIGHTS_DURING_CHISQSUM.ERRORS, yNorm);

            //GEVYFit yfit = chiSqMin.fitCurveKGreaterThanZero(GEVChiSquareMinimization.WEIGHTS_DURING_CHISQSUM.ERRORS, k, s, mu, yNorm);

            //GEVYFit yfit = chiSqMin.fitCurveKGreaterThanZero(GEVChiSquareMinimization.WEIGHTS_DURING_CHISQSUM.ERRORS);

            GEVYFit yfit = chiSqMin.fitCurveKGreaterThanZeroAndMu(GEVChiSquareMinimization.WEIGHTS_DURING_CHISQSUM.ERRORS);

            float[] xf = yfit.getOriginalScaleX();
            float[] yf = yfit.getOriginalScaleYFit();

            plotter.addPlot(histogram.getXHist(), histogram.getYHistFloat(),
                histogram.getXErrors(), histogram.getYErrors(), xf, yf, "");
            plotter.writeFile();

            log.info(yfit.toString());

            if (findClusters) {

                voidStats.finalizeStats(histogram, yfit);

                TwoPointCorrelation clusterFinder = new TwoPointCorrelation(indexerFilePath);

                clusterFinder.setBackground(voidStats.getBackgroundSurfaceDensity(), voidStats.getBackgroundSurfaceDensityError());

                clusterFinder.findClusters();

                clusterFinder.plotClusters();
            }

            log.info("  END " + i);
        }
    }

    public static String[] getDensityFilePaths() throws IOException {

        String tmpDataDirPath = ResourceFinder.findTmpDataDirectory();

        File tmpDataDir = new File(tmpDataDirPath);

        File[] indexerFiles = tmpDataDir.listFiles(new DensityFileFilter());

        String[] filePaths = new String[indexerFiles.length];

        for (int i = 0; i < indexerFiles.length; i++) {
            filePaths[i] = indexerFiles[i].getPath();
        }

        return filePaths;
    }

    protected static class DensityFileFilter implements FileFilter {
        @Override
        public boolean accept(File pathname) {
            if (pathname.getName().startsWith("stats_2pt_voids_")) {
                return true;
            }
            return false;
        }
    }

    protected String getIndexerFilePath(String[] indexerFilePaths, String twoPtDensityFilePath) {
        //"stats_2pt_voids_" + String.valueOf(nFile)

        int i0 = twoPtDensityFilePath.indexOf("voids_");
        i0 += "voids_".length();
        String fileNumber = twoPtDensityFilePath.substring(i0, i0 + 3);
        fileNumber = "_" + fileNumber;

        for (String indexerFilePath : indexerFilePaths) {
            if (indexerFilePath.contains("indexer_random_background_with_3_clusters_")) {
                if (indexerFilePath.contains( fileNumber )) {
                    return indexerFilePath;
                }
            }
        }
        return null;
    }

    public class TwoPointVoidStatsExt extends TwoPointVoidStats {

        public TwoPointVoidStatsExt(AxisIndexer indexedSortedPoints) {
            super(indexedSortedPoints);
        }

        public TwoPointVoidStatsExt(String persistedIndexerFileName) throws IOException {
            super(persistedIndexerFileName);
        }

        @Override
        protected void finalizeStats(HistogramHolder histogram, GEVYFit yfit) throws TwoPointVoidStatsException {

            state = State.HISTOGRAM_FITTED;

            super.finalizeStats(histogram, yfit);
        }
    }
}
