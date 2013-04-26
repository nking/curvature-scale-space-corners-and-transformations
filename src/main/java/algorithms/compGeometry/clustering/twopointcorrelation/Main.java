package algorithms.compGeometry.clustering.twopointcorrelation;

import algorithms.util.MainRunner;

/**

  Usage from the command line:
      Requires a tab delimited text file with 4 columns: x, y, xErrors, yErrors.

          java -jar two-point-correlation.jar --file /path/to/file/fileName.txt

  @author nichole
 */
public class Main {

    public static void main(String[] args) {

        MainRunner runner = new MainRunner();

        runner.run(args);
    }
}
