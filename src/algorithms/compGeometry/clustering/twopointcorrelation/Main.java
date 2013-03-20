package algorithms.compGeometry.clustering.twopointcorrelation;

import algorithms.util.InputFileReader;
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


  Usage from the command line:
      Requires a tab delimited text file with 4 columns: x, y, xErrors, yErrors.

          java -jar two-point-correlation.jar --file /path/to/file/fileName.txt

  @author nichole
 */
public class Main {

    public static void main(String[] args) {

        if (args == null || args.length < 2) {
            System.out.println("Requires file:  --file /path/to/file/fileName.txt");
            return;
        }

        String filePath = null;
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("--file")) {
                if ( (i+1) < args.length) {
                    filePath = args[i+1];
                }
            }
        }
        if (filePath == null) {
            System.out.println("Requires file:  --file /path/to/file/fileName.txt");
            return;
        }

        try {

            InputFileReader reader = new InputFileReader(filePath);
            reader.read();

            Logger.getLogger(Main.class.getName()).info("begin TwoPointCorrelation");

            TwoPointCorrelation clusterFinder = new TwoPointCorrelation(
                reader.getX(), reader.getY(), reader.getXErrors(), reader.getYErrors(), reader.getX().length);

            clusterFinder.setDebug(false);

            Logger.getLogger(Main.class.getName()).info("calculate background");

            clusterFinder.calculateBackground();

            Logger.getLogger(Main.class.getName()).info("find clusters");

            clusterFinder.findClusters();

            Logger.getLogger(Main.class.getName()).info("calculate hull of clusters");

            clusterFinder.calculateHullsOfClusters();

            String plotFilePath = clusterFinder.plotClusters();

        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
