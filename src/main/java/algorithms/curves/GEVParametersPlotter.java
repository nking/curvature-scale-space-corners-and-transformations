package algorithms.curves;

import algorithms.util.ResourceFinder;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;

/**
 * fills in data for gev_plot.html to compare a few
 * curve measurements with GEV curve parameters
 *
 * @author nichole
 */
public class GEVParametersPlotter {

    /**
     *
     */
    protected final StringBuffer plotContent;

    /**
     *
     */
    protected int plotNumber = 0;

    /**
     *
     * @throws FileNotFoundException
     * @throws IOException
     */
    public GEVParametersPlotter() throws FileNotFoundException, IOException {
        plotContent = getTemplateHtmlPlot();
    }

    /**
     *
     * @param k
     * @param sigma
     * @param mu
     * @param slope0
     * @param slope1
     * @param slope2
     */
    public void addPlot(float[] k, float[] sigma, float[] mu, float[] slope0, float[] slope1, float[] slope2) {

        float kMin      = Float.MAX_VALUE;
        float sigmaMin  = Float.MAX_VALUE;
        float muMin     = Float.MAX_VALUE;

        float kMax      = Float.MIN_VALUE;
        float sigmaMax  = Float.MIN_VALUE;
        float muMax     = Float.MIN_VALUE;

        float slope0Min = Float.MAX_VALUE;
        float slope1Min = Float.MAX_VALUE;
        float slope2Min = Float.MAX_VALUE;

        float slope0Max = Float.MIN_VALUE;
        float slope1Max = Float.MIN_VALUE;
        float slope2Max = Float.MIN_VALUE;

        StringBuffer dataSB = new StringBuffer();

        //  ===== add points data =====
        dataSB.append("\n\n").append("var data_").append(plotNumber).append(" = [\n");
        for (int i = 0; i < k.length; i++) {

            if (i > 0) {
                dataSB.append(",\n");
            }

            dataSB.append("    {k:").append(k[i]).append(", sigma:").append(sigma[i])
                .append(", mu:").append(mu[i]).append(", slope0:").append(slope0[i])
                .append(", slope1:").append(slope1[i]).append(", slope2:").append(slope2[i]).append("}");

            if (k[i] < kMin) {
                kMin = k[i];
            }
            if (k[i] > kMax) {
                kMax = k[i];
            }
            if (sigma[i] < sigmaMin) {
                sigmaMin = sigma[i];
            }
            if (sigma[i] > sigmaMax) {
                sigmaMax = sigma[i];
            }
            if (mu[i] < muMin) {
                muMin = mu[i];
            }
            if (mu[i] > muMax) {
                muMax = mu[i];
            }
            if (slope0[i] < slope0Min) {
                slope0Min = slope0[i];
            }
            if (slope0[i] > slope0Max) {
                slope0Max = slope0[i];
            }
            if (slope1[i] < slope1Min) {
                slope1Min = slope1[i];
            }
            if (slope1[i] > slope1Max) {
                slope1Max = slope1[i];
            }
            if (slope2[i] < slope2Min) {
                slope2Min = slope2[i];
            }
            if (slope2[i] > slope2Max) {
                slope2Max = slope2[i];
            }
        }
        dataSB.append("\n];\n");

        String kMinStr = String.format("%.7f", kMin);
        String kMaxStr = String.format("%.7f", kMax);
        String sigmaMinStr = String.format("%.7f", sigmaMin);
        String sigmaMaxStr = String.format("%.7f", sigmaMax);
        String muMinStr = String.format("%.7f", muMin);
        String muMaxStr = String.format("%.7f", muMax);

        String slope0MinStr = String.format("%.7f", slope0Min);
        String slope1MinStr = String.format("%.7f", slope1Min);
        String slope2MinStr = String.format("%.7f", slope2Min);

        String slope0MaxStr = String.format("%.7f", slope0Max);
        String slope1MaxStr = String.format("%.7f", slope1Max);
        String slope2MaxStr = String.format("%.7f", slope2Max);

        dataSB.append("\n").append("var min_maxes_").append(plotNumber).append(" = {\n");
        dataSB.append("    k_min:").append(kMinStr).append(", k_max:").append(kMaxStr)
            .append(", sigma_min:").append(sigmaMinStr).append(", sigma_max:").append(sigmaMaxStr)
            .append(", mu_min:").append(muMinStr).append(", mu_max:").append(muMaxStr)
            .append(", slope0_min:").append(slope0MinStr).append(", slope0_max:").append(slope0MaxStr)
            .append(", slope1_min:").append(slope1MinStr).append(", slope1_max:").append(slope1MaxStr)
            .append(", slope2_min:").append(slope2MinStr).append(", slope2_max:").append(slope2MaxStr)
            .append("};\n");


        // ======= add RENDER statement ==========
        dataSB.append("\nrenderPlot('plot").append(plotNumber).append("', data_").append(plotNumber)
            .append(", min_maxes_").append(plotNumber)
            .append(");\n\n");


        String srchFor = "/* === DO NOT REMOVE THIS == END DATA */";
        int insertOffset = plotContent.indexOf(srchFor);
        if (insertOffset == -1) {
            throw new IllegalStateException("cannot find END DATA marker");
        }
        plotContent.insert(insertOffset, dataSB.toString());
        dataSB = null;

        // ========== add the PLOT DIVS ==============
        StringBuffer plotDivs = new StringBuffer();
        plotDivs.append("<div id='plot").append(plotNumber).append("'              class='plot'></div>\n");
        plotDivs.append("<div id='plot").append(plotNumber).append("_k_slope0'     class='plot'></div>\n");
        plotDivs.append("<div id='plot").append(plotNumber).append("_k_slope1'     class='plot'></div>\n");
        plotDivs.append("<div id='plot").append(plotNumber).append("_k_slope2'     class='plot'></div>\n");
        plotDivs.append("<div id='plot").append(plotNumber).append("_sigma_slope0' class='plot'></div>\n");
        plotDivs.append("<div id='plot").append(plotNumber).append("_sigma_slope1' class='plot'></div>\n");
        plotDivs.append("<div id='plot").append(plotNumber).append("_sigma_slope2' class='plot'></div>\n");
        plotDivs.append("<div id='plot").append(plotNumber).append("_mu_slope0'    class='plot'></div>\n");
        plotDivs.append("<div id='plot").append(plotNumber).append("_mu_slope1'    class='plot'></div>\n");
        plotDivs.append("<div id='plot").append(plotNumber).append("_mu_slope2'    class='plot'></div>\n");

        srchFor = "<!-- === DO NOT REMOVE THIS == END PLOT DIVS -->";
        insertOffset = plotContent.indexOf(srchFor);
        if (insertOffset == -1) {
            throw new IllegalStateException("cannot find END DATA marker");
        }
        plotContent.insert(insertOffset, plotDivs.toString());
        plotDivs = null;

        plotNumber++;
    }

    /**
     *
     * @return
     * @throws FileNotFoundException
     * @throws IOException
     */
    protected StringBuffer getTemplateHtmlPlot() throws FileNotFoundException, IOException {
        return getTemplateHtmlPlot("plot_gev.html");
    }

    /**
     *
     * @param fileName
     * @return
     * @throws FileNotFoundException
     * @throws IOException
     */
    protected StringBuffer getTemplateHtmlPlot(String fileName) throws FileNotFoundException, IOException {

        StringBuffer sb = new StringBuffer();

        Reader reader = null;
        BufferedReader in = null;

        try {
            String path = ResourceFinder.findFileInResources(fileName);

            reader = new FileReader(new File(path));
            in = new BufferedReader(reader);

            String line = in.readLine();

            while (line != null) {
                sb.append(line).append("\n");
                line = in.readLine();
            }

        } catch(IOException e) {

            ClassLoader cls = ResourceFinder.class.getClassLoader();

            InputStream input = cls.getResourceAsStream(fileName);

            if (input == null) {
                throw new IOException("could not find file " + fileName);
            }

            reader = new InputStreamReader(input);
            in = new BufferedReader(reader);

            String line = in.readLine();

            while (line != null) {
                sb.append(line).append("\n");
                line = in.readLine();
            }

        } finally {
            if (in != null) {
                in.close();
            }
            if (reader != null) {
                reader.close();
            }
        }

        return sb;
    }

    /**
     *
     * @throws IOException
     */
    public void writeFile() throws IOException {
        writeToFile(this.plotContent.toString(), "gev.html");
    }

    /**
     *
     * @param fileContent
     * @param fileName
     * @throws IOException
     */
    protected void writeToFile(String fileContent, String fileName) throws IOException {

        String copyFilePath = ResourceFinder.writeToCWD(fileContent, fileName);
    }
}
