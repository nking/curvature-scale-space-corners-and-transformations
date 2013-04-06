package algorithms.compGeometry.clustering.twopointcorrelation;

import algorithms.misc.HistogramHolder;
import algorithms.util.PolygonAndPointPlotter;
import algorithms.util.ResourceFinder;
import java.io.File;
import java.io.FileFilter;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FilterInputStream;
import java.io.FilterOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.security.SecureRandom;
import java.util.logging.Logger;

/**
 * This persists the cluster data in files for use elsewhere.  the test
 * that generates the data is disabled when not needed.
 *
 * @author nichole
 */
public class CreateClusterDataTest extends BaseTwoPointTest {

    protected static String histogramFileNamePrefix = "histogram_random_background_with_";

    protected static String indexerFileNamePrefix = "indexer_random_background_with_";

    protected Logger log = Logger.getLogger(this.getClass().getSimpleName());

    protected boolean enable = false;

    public void testCreateData() throws Exception {

        if (!enable) {
            return;
        }

        log.info("testCreateData()");

        float xmin = 0;
        float xmax = 300;
        float ymin = 0;
        float ymax = 300;

        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter();

        SecureRandom srr = SecureRandom.getInstance("SHA1PRNG");
        srr.setSeed(System.currentTimeMillis());
        long seed = srr.nextLong();

        seed = 310357278571620991l;

        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        sr.setSeed(seed);
        //sr.setSeed(4066852271294423764l);
        //sr.setSeed(310357278571620991l);
        //sr.setSeed(-1993887065899734688l);
        //sr.setSeed(-6886733535826319879l);
        //sr.setSeed(-3765842324512485314l);
        //sr.setSeed(1152752110035096347l);
        //sr.setSeed(-6221198867223436351l);
        //sr.setSeed(6899554926901724961l);

        log.info("SEED=" + seed);

        // a long running test to calculcate and print the stats of fits
        //  for sparse, moderate, and densely populated backgrounds,
        //  all with the same number of clusters and cluster points, though
        //  randomly distributed.

        int nIterPerBackground = 30;

        int m = nIterPerBackground*6;

        DoubleAxisIndexer indexer = null;

        int count = 0;

        for (int i = 0; i < 6; i++) {

            for (int ii = 0; ii < nIterPerBackground; ii++) {

                String numberOfClusters = "";
                String str = String.valueOf(count);
                while (str.length() < 3) {
                    str = "0" + str;
                }
                String fileNamePostfix = "_clusters_" + str + ".dat";

                switch(i) {
                    case 0:
                        numberOfClusters = "0";
                        indexer = createIndexerWithRandomPoints(sr, xmin, xmax, ymin, ymax,
                            new int[0], 100, null);
                        break;
                    case 1:
                        numberOfClusters = "3";
                        indexer = createIndexerWithRandomPoints(sr, xmin, xmax, ymin, ymax,
                            3, 30, 60, 0.1f);
                        break;
                    case 2:
                        numberOfClusters = "1";
                        indexer = createIndexerWithRandomPoints(sr, xmin, xmax, ymin, ymax,
                            1, 30, 60, 1f);
                        break;
                    case 3:
                        numberOfClusters = "2";
                        indexer = createIndexerWithRandomPoints(sr, xmin, xmax, ymin, ymax,
                            2, 30, 60, 1f);
                        break;
                    case 4:
                        numberOfClusters = "3";
                        indexer = createIndexerWithRandomPoints(sr, xmin, xmax, ymin, ymax,
                            3, 30, 60, 1f);
                        break;
                    case 5:
                        numberOfClusters = "3";
                        indexer = createIndexerWithRandomPoints(sr, xmin, xmax, ymin, ymax,
                            3, 30, 60, 10.0f);
                        break;
                    default:
                        break;
                }

                indexer.sortAndIndexXThenY(x, y, xErrors, yErrors, x.length);

                String fileName = indexerFileNamePrefix + numberOfClusters + fileNamePostfix;
                String filePath = ResourceFinder.getAFilePathInTmpData(fileName);


                writeIndexer(filePath, indexer);


                log.info(" " + count + " " + i + " number of clusters=" + numberOfClusters + " (" + indexer.nXY + " points) ... ");

                TwoPointVoidStatsExt twoPtC = new TwoPointVoidStatsExt(indexer);
                twoPtC.setDebug(true);

                float[] xf = null;
                float[] yf = null;
                HistogramHolder histogram = twoPtC.constructAndFitHistogram();
                if (twoPtC.bestFit != null) {
                    xf = twoPtC.bestFit.getOriginalScaleX();
                    yf = twoPtC.bestFit.getOriginalScaleYFit();
                }

                fileName = histogramFileNamePrefix + numberOfClusters + fileNamePostfix;
                filePath = ResourceFinder.getAFilePathInTmpData(fileName);


                writeHistogram(filePath, histogram);


                plotter.addPlot(histogram.getXHist(), histogram.getYHistFloat(),
                    histogram.getXErrors(), histogram.getYErrors(), xf, yf,
                    String.valueOf(count) + " : " + i + " : " + numberOfClusters + " => " + twoPtC.bestFit.getChiSqStatistic());

                plotter.writeFile();

                count++;
            }
        }

        log.info("SEED=" + seed);
    }

    public static String[] getHistogramFilePaths() throws IOException {

        String tmpDataDirPath = ResourceFinder.findTmpDataDirectory();

        File tmpDataDir = new File(tmpDataDirPath);

        File[] histogramFiles = tmpDataDir.listFiles(new HistogramFileFilter());

        String[] filePaths = new String[histogramFiles.length];

        for (int i = 0; i < histogramFiles.length; i++) {
            filePaths[i] = histogramFiles[i].getPath();
        }

        return filePaths;
    }

    public static String[] getIndexerFilePaths() throws IOException {

        String tmpDataDirPath = ResourceFinder.findTmpDataDirectory();

        File tmpDataDir = new File(tmpDataDirPath);

        File[] indexerFiles = tmpDataDir.listFiles(new IndexerFileFilter());

        String[] filePaths = new String[indexerFiles.length];

        for (int i = 0; i < indexerFiles.length; i++) {
            filePaths[i] = indexerFiles[i].getPath();
        }

        return filePaths;
    }

    protected static class HistogramFileFilter implements FileFilter {
        @Override
        public boolean accept(File pathname) {
            if (pathname.getName().startsWith(histogramFileNamePrefix)) {
                return true;
            }
            return false;
        }
    }

    protected static class IndexerFileFilter implements FileFilter {
        @Override
        public boolean accept(File pathname) {
            if (pathname.getName().startsWith(indexerFileNamePrefix)) {
                return true;
            }
            return false;
        }
    }

    public static DoubleAxisIndexer readIndexer(String filePath) throws IOException {
        return SerializerUtil.readPersistedPoints(filePath, true);
    }

    public static HistogramHolder readHistogram(String filePath) throws IOException {

        HistogramHolder histogram = new HistogramHolder();

        FileInputStream fileInputStream = null;
        FilterInputStream filterInputStream = null;
        ObjectInputStream objectInputStream = null;

        try {

            fileInputStream = new FileInputStream(filePath);

            objectInputStream = new ObjectInputStream(fileInputStream);

            histogram.readExternal(objectInputStream);

            return histogram;

        } finally {
            if (fileInputStream != null) {
                fileInputStream.close();
            }
            if (filterInputStream != null) {
                filterInputStream.close();
            }
            if (objectInputStream != null) {
                objectInputStream.close();
            }
        }
    }

    protected static void writeIndexer(String filePath, DoubleAxisIndexer indexer) throws IOException {
        SerializerUtil.serializeIndexer(indexer, filePath);
    }

    protected static void writeHistogram(String filePath, HistogramHolder histogram) throws IOException {

        FileOutputStream fileOutputStream = null;
        FilterOutputStream filterOutputStream = null;
        ObjectOutputStream objectOutputStream = null;

        try {

            fileOutputStream = new FileOutputStream(filePath);

            objectOutputStream = new ObjectOutputStream(fileOutputStream);

            histogram.writeExternal(objectOutputStream);

        } finally {
            if (fileOutputStream != null) {
                fileOutputStream.close();
            }
            if (filterOutputStream != null) {
                filterOutputStream.close();
            }
            if (objectOutputStream != null) {
                objectOutputStream.close();
            }
        }
    }

    protected class TwoPointVoidStatsExt extends TwoPointVoidStats {

        public TwoPointVoidStatsExt(DoubleAxisIndexer indexer) {
            super(indexer);
        }

        public HistogramHolder constructAndFitHistogram() throws TwoPointVoidStatsException {

            calc();

            return statsHistogram;
        }
    }
}
