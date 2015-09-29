package algorithms.compGeometry.clustering.twopointcorrelation;

import algorithms.util.IInputFileReader;
import algorithms.util.InputFileReader;
import algorithms.util.MainRunner;

/**
<pre>
  Usage from the command line:
      Requires a tab delimited text file with 4 columns: x, y, xErrors, yErrors.

          java -jar two-point-correlation.jar --file /path/to/file/fileName.txt
</pre>
  @author nichole
 */
public class Main {

    /**
     *
     * @param args
     */
    public static void main(String[] args) {

        MainRunner runner = new MainRunner();

        IInputFileReader reader = new InputFileReader();

        runner.run(reader, args);
    }
}
