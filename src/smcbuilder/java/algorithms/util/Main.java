package algorithms.util;

import algorithms.compGeometry.clustering.twopointcorrelation.TwoPointVoidStatsException;
import java.io.FileNotFoundException;
import java.io.IOException;

/**
  Usage from the command line:
      Requires a tab delimited text file with 4 columns: x, y, xErrors, yErrors.

          java -jar two-point-correlation.jar --file /path/to/file/fileName.txt

  @author nichole
 */
public class Main {

    public static void main(String[] args) {

        String dir = System.getProperty("user.dir") + "/photfiles";

        String[] files = new String[]{
            //"smc118.1.phot", //  3 min
            //"smc118.6.phot", //    min
            //"smc118.8.phot", //    min
            //"smc118.7.phot", //    min
            //"smc115.5.phot", //    min
            //"smc116.3.phot", //    min
            //"smc111.4.phot", // 10 min         <==== REDO with trim for xp < 2100
            //"smc110.3.phot", // 10 min         <==== REDO
            //"smc110.2.phot", // 13.5 min
            //"smc114.7.phot", //  9 - 14 min when using _3.  3 min when using _4
            //"smc101.4.phot", //  9 min using _4
            "smc125.4.phot",   //    min using _4 26117 points
            "smc100.1.phot",   //
            "smc100.5.phot",   //
            "smc125.3.phot",   //
            "smc108.7.phot",   //   ~46000 points
            "smc108.8.phot"    //
        };

        try {

            for (int i = 0; i < files.length; i++) {

                String filePath = dir + "/" + files[i];

                String[] args2 = new String[] {
                    "--file",
                    filePath,
                    "--twosigma"
                };

                run(args2);
            }

        } catch (Exception e) {

            System.err.println(e.getMessage());
        }
    }

    public static void run(String[] args) throws IOException, FileNotFoundException, TwoPointVoidStatsException {

        MainRunner runner = new MainRunner();

        runner.run(args);

        String filePath = runner.getArguments().get("filePath");

        //photfiles/smc118.1.phot ==> photfiles/smc118.1_clusters.txt
        String outFilePath = filePath.replace(".phot", "_clusters.txt");
        outFilePath = outFilePath.replace("photfiles", "clusterfiles");

        OutputFileWriter writer = new OutputFileWriter(outFilePath);
        writer.writeFile(runner.getTwoPointCorrelation(), outFilePath);
    }

}
