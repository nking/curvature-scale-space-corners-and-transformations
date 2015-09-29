package algorithms.util;

import algorithms.compGeometry.clustering.twopointcorrelation.TwoPointCorrelation;
import algorithms.compGeometry.clustering.twopointcorrelation.TwoPointVoidStatsException;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author nichole
 */
public class OutputFileWriter {

    /**
     *
     */
    protected final Logger log;

    /**
     *
     * @param pathToFile
     */
    public OutputFileWriter(String pathToFile) {

        log = Logger.getLogger(this.getClass().getName());
    }

    /**
     *
     * @param clusterFinder
     * @param filePath
     * @throws IOException
     * @throws FileNotFoundException
     * @throws algorithms.compGeometry.clustering.twopointcorrelation.TwoPointVoidStatsException
     */
    public void writeFile(TwoPointCorrelation clusterFinder, String filePath)
        throws IOException, FileNotFoundException, TwoPointVoidStatsException {

        FileWriter writer = null;

        try {
            File fl = new File(filePath);
            fl.createNewFile();

            writer = new FileWriter(filePath);

            int nGroups = clusterFinder.getNumberOfGroups();

            System.out.println("Begin writing " + nGroups + " to file " + filePath);

            for (int i = 0; i < nGroups; i++) {

                float[] x = clusterFinder.getXGroup(i);
                float[] y = clusterFinder.getYGroup(i);

                for (int ii = 0; ii < x.length; ii++) {
                    String line = String.format("%.4f\t%.4f\t%4d\n", x[ii], y[ii], i);
                    writer.write(line);
                    writer.flush();
                }
            }

            String htmlFilePath = clusterFinder.plotClusters();

            String htmlFile = filePath.replace(".txt", ".html");

            // base directory
            String baseDir = System.getProperty("user.dir");
            String voidStatsInFile = baseDir + "/"+ "twoptvoid_stats.html";
            String voidStatsOutFile = htmlFile.replace("clusters", "twoptvoid_stats");

            Runtime runtime = Runtime.getRuntime();
            try {
                StringBuilder sb = new StringBuilder();

                String cp = "cp " + htmlFilePath + " " + htmlFile;
                Process proc = runtime.exec(cp);
                int exit = proc.waitFor();
                if (exit == 0) {
                    String content = getContent(proc.getInputStream());
                    int index = content.indexOf(')');
                    if (index >= 0) {
                        sb.append(content.substring(0, index + 1) );
                    }
                }
                System.out.println(sb.toString());

                if ( (new File(voidStatsInFile)).exists()) {

                    cp = "cp " + voidStatsInFile + " " + voidStatsOutFile;
                    proc = runtime.exec(cp);
                    exit = proc.waitFor();
                    if (exit == 0) {
                        String content = getContent(proc.getInputStream());
                        int index = content.indexOf(')');
                        if (index >= 0) {
                            sb.append(content.substring(0, index + 1) );
                        }
                    }
                    System.out.println(sb.toString());
                }

            } catch (IOException e) {
                e.printStackTrace();
            } catch (InterruptedException e) {
                e.printStackTrace();
            }

        } finally {
            if (writer != null) {
                writer.close();
            }

            log.log(Level.INFO, "finished writing file: {0}", filePath);
        }
    }

    /**
     *
     * @param input
     * @return
     * @throws IOException
     */
    public static String getContent(InputStream input) throws IOException {
        if (input == null) {
            return null;
        }
        byte[] b = new byte[1024];
        int readBytes = 0;
        StringBuilder result = new StringBuilder();
        while ((readBytes = input.read(b)) >= 0) {
            result.append(new String(b, 0, readBytes));
        }
        return result.toString();
    }
}
