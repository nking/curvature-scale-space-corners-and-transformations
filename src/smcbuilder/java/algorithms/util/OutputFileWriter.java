package algorithms.util;

import algorithms.compGeometry.clustering.twopointcorrelation.TwoPointCorrelation;
import algorithms.compGeometry.clustering.twopointcorrelation.TwoPointVoidStatsException;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.util.logging.Logger;

/**
 *
 * @author nichole
 */
public class OutputFileWriter {

    protected Logger log = null;

    public OutputFileWriter(String pathToFile) {

        log = Logger.getLogger(this.getClass().getName());
    }

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

            // copy html file

            //photfiles/smc118.1_clusters.txt
            String htmlFilePath = clusterFinder.plotClusters();
            String htmlFile = filePath.replace(".txt", ".html");

            // base directory
            String baseDir = getBaseDirectory();
            String voidStatsInFile = baseDir + "/"+ "twoptvoid_stats.html";
            String voidStatsOutFile = htmlFile.replace("clusters", "twoptvoid_stats");

            Runtime runtime = Runtime.getRuntime();
            try {
                StringBuffer sb = new StringBuffer();

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

            log.info("finished writing file: " + filePath);
        }
    }

    public static String getContent(InputStream input) throws
        IOException {
        if (input == null) {
            return null;
        }
        byte[] b = new byte[1024];
        int readBytes = 0;
        StringBuffer result = new StringBuffer();
        while ((readBytes = input.read(b)) >= 0) {
            result.append(new String(b, 0, readBytes));
        }
        return result.toString();
    }

    private String getBaseDirectory() {

        ClassLoader cls = this.getClass().getClassLoader();

        return cls.getResource(".").getPath();
    }

}
