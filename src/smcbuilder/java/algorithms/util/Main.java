package algorithms.util;

import algorithms.compGeometry.clustering.twopointcorrelation.TwoPointVoidStatsException;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map.Entry;

/**
  Usage from the command line:
      Requires a tab delimited text file with 4 columns: x, y, xErrors, yErrors.

          java -jar two-point-correlation.jar

  @author nichole
 */
public class Main {

    public static void main(String[] args) {

        String dir = System.getProperty("user.dir") + "/photfiles";

        String[] files = new String[]{
            /*"smc118.1.phot", // 4108 points
            "smc118.6.phot",
            "smc118.8.phot",
            "smc118.7.phot",
            "smc115.5.phot",
            "smc116.3.phot",
            "smc111.4.phot",
            "smc110.3.phot",
            "smc110.2.phot",
            "smc114.7.phot",
            "smc101.4.phot",
            "smc125.4.phot",
            "smc100.1.phot",
            "smc100.5.phot",
            "smc125.3.phot",
            "smc108.7.phot",
            "smc108.8.phot"*/
            "smc110.3.phot"
        };

        try {

            for (int i = 0; i < files.length; i++) {

                String filePath = dir + "/" + files[i];

                String[] args2 = new String[] {
                    "--file", filePath,
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

        IInputFileReader reader = new InputPhotFileReader();

        runner.run(reader, args);

        HashMap<String, String> arguments = runner.getArguments();
        String filePath = arguments.get("filePath");

        Iterator<Entry<String,String>> iter = arguments.entrySet().iterator();
        while (iter.hasNext()) {
            Entry<String, String> entry = iter.next();
            System.out.println("key=" + entry.getKey() + " value=" + entry.getValue());
        }

        String outFilePath = filePath.replace(".phot", "_clusters.txt");

        outFilePath = outFilePath.replace("photfiles", "clusterfiles");

        OutputFileWriter writer = new OutputFileWriter(outFilePath);
        writer.writeFile(runner.getTwoPointCorrelation(), outFilePath);
    }

}
